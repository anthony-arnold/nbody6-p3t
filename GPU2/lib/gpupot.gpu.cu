// #include <iostream>
#include <cstdio>
// #include <cutil.h>
#ifdef WITH_CUDA5
#  include <helper_cuda.h>
#  define CUDA_SAFE_CALL checkCudaErrors
#else
#  include <cutil.h>
#endif
#include "cuda_pointer.h"
#define NTHREAD 128

#if 0
#define PROFILE
#endif

extern "C" double double_env(const char*);

namespace pot {

#ifdef PROFILE
#include <sys/time.h>
double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#else
double get_wtime(){
	return 0.0;
}
#endif

   double get_eps() {
      static double EPS = double_env("EPS");
      return EPS;
   }

float2 float2_split(double x){
	const int shift = 20;
	float2 ret;
	x *= (1<<shift);
	double xi = (int)x;
	double xf = x - xi;
	ret.x = xi * (1./(1<<shift));
	ret.y = xf * (1./(1<<shift));
	return ret;
}
__device__ float2 float2_accum(float2 acc, float x){
	float tmp = acc.x + x;
	acc.y -= (tmp - acc.x) - x;
	acc.x = tmp;
	return acc;
}

__device__ float2 float2_regularize(float2 acc){
	float tmp = acc.x + acc.y;
	acc.y = acc.y -(tmp - acc.x);
	acc.x = tmp;
	return acc;
}

struct Particle{
	float2 pos[3];
	float mass;
	float pad;

	Particle(double x[3], double m){
		pos[0] = float2_split(x[0]);
		pos[1] = float2_split(x[1]);
		pos[2] = float2_split(x[2]);
		mass = (float)m;
	}
	Particle(int){
		pos[0].x = pos[0].y = pos[1].x = pos[1].y = pos[2].x = pos[2].y = mass = pad = 0.f;
	}
	__device__ Particle() {}
};

__global__ void pot_kernel(int n, double eps2, Particle *ptcl, float2 *phi){
	__shared__ Particle jpbuf[NTHREAD];
	int i = NTHREAD * blockIdx.x + threadIdx.x;
	Particle ip = ptcl[i];
	float2 phii = make_float2(0.f, 0.f);
	for(int j=0; j<n; j+= NTHREAD){
		__syncthreads();
		jpbuf[threadIdx.x] = ptcl[j + threadIdx.x];
		__syncthreads();
#pragma unroll 4
		for(int jj=0; jj<NTHREAD; jj++){
			// if(j+jj == i) continue;
			Particle &jp = jpbuf[jj];
			float dx = (jp.pos[0].x - ip.pos[0].x) + (jp.pos[0].y - ip.pos[0].y);
			float dy = (jp.pos[1].x - ip.pos[1].x) + (jp.pos[1].y - ip.pos[1].y);
			float dz = (jp.pos[2].x - ip.pos[2].x) + (jp.pos[2].y - ip.pos[2].y);
			float r2 = dx*dx + dy*dy + dz*dz + eps2;
			if(r2>0.f) {
                            float pij = jp.mass * rsqrtf(r2);
                            phii = float2_accum(phii, pij);
                        }
		}
		phii = float2_regularize(phii);
	}
	phi[i] = phii;
}

extern "C"  void gpunb_devinit_();

void gpupot(
		int n,
		double m[],
		double x[][3],
		double pot[]){
        gpunb_devinit_();

	double t0 = get_wtime();
#if 0
	float2 *phi_d, *phi_h;
	Particle *ptcl_d, *ptcl_h;
#else
	cudaPointer <float2> phi;
	cudaPointer <Particle> ptcl;
#endif
	int ng = NTHREAD * (n/NTHREAD + (n%NTHREAD ? 1 : 0));

#if 0
	cudaMalloc    ((void **)&phi_d, ng * sizeof(float2));
	cudaMallocHost((void **)&phi_h, ng * sizeof(float2));
	cudaMalloc    ((void **)&ptcl_d, ng * sizeof(Particle));
	cudaMallocHost((void **)&ptcl_h, ng * sizeof(Particle));
#else
	phi.allocate(ng);
	ptcl.allocate(ng);
#endif

	// std::cout << n << " " << ng << std::endl;
	for(int i=0; i<n; i++){
		// ptcl_h[i] = Particle(x[i], m[i]);
		ptcl[i] = Particle(x[i], m[i]);
	}
	for(int i=n; i<ng; i++){
		// ptcl_h[i] = Particle(0);
		ptcl[i] = Particle(0);
	}

	// cudaMemcpy(ptcl_d, ptcl_h, ng * sizeof(Particle), cudaMemcpyHostToDevice);
	ptcl.htod(ng);

	dim3 grid(ng/NTHREAD, 1, 1);
	dim3 threads(NTHREAD, 1, 1);
	int sharedMemSize = NTHREAD * sizeof(Particle);
        double eps2 = get_eps() * get_eps();
	// pot_kernel <<<grid, threads, sharedMemSize >>> (n, ptcl_d, phi_d);
	pot_kernel <<<grid, threads, sharedMemSize >>> (n, eps2, ptcl, phi);

	// cudaMemcpy(phi_h, phi_d, n * sizeof(float2), cudaMemcpyDeviceToHost);
	phi.dtoh(n);
	for(int i=0; i<n; i++){
		// pot[i] = (double)phi_h[i].x + (double)phi_h[i].y;
		pot[i] = (double)phi[i].x + (double)phi[i].y;
	}

#if 0
	cudaFree    (phi_d);
	cudaFreeHost(phi_h);
	cudaFree    (ptcl_d);
	cudaFreeHost(ptcl_h);
#else
	phi.free();
	ptcl.free();
#endif
	double t1 = get_wtime();
#ifdef PROFILE
	fprintf(stderr, "gpupot: %f sec\n", t1 - t0);
#endif
}

} // namespace pot

//#include <fenv.h>
extern "C"{
	void gpupot_(
			int *n,
			double m[],
			double x[][3],
			double pot[]){
            //fedisableexcept(FE_INVALID | FE_DIVBYZERO);
            pot::gpupot(*n, m, x, pot);
            //feenableexcept(FE_INVALID | FE_DIVBYZERO);
	}
}
