#define __USE_GNU
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <unistd.h>
#include <sched.h>
#include <omp.h>
#include "vector3.h"

#define _out_
#define PROFILE

#ifdef PROFILE
#include <sys/time.h>
static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#else
static double get_wtime(){
	return 0.0;
}
#endif

typedef float  v4sf __attribute__((vector_size(16)));
typedef double v2df __attribute__((vector_size(16)));

static v4sf rsqrt_NR(const v4sf x){
	const v4sf y = __builtin_ia32_rsqrtps(x);
	const v4sf c1 = {-0.5f, -0.5f, -0.5f, -0.5f};
	const v4sf c2 = {-3.0f, -3.0f, -3.0f, -3.0f};
	return (c1 * y) * (x*y*y + c2);

}

#if 0
struct float4{
	float x, y, z, w;

	float4(){}
	float4(const double p[3]){
		x = float(p[0]);
		y = float(p[1]);
		z = float(p[2]);
		w = 0.0f;
	}
	operator v4sf() const{
		return *(v4sf *)&x;
	}
} __attribute__((aligned(16))) ;

struct float2{
	float x, y;

	float2(const double d){
		x = float(d);
		y = float(d - double(x));
	}
};

struct Particle{
	float4 posH; // elemen w is for mass
	float4 posL;
	float4 vel;
	float4 acc2;
	float4 jrk6;
	double time, pad; // 6 xmm words;

	Particle(
		const double _pos [3],
		const double _vel [3],
		const double _acc2[3],
		const double _jrk6[3],
		const double _mass,
		const double _time)
		: vel(_vel), acc2(_acc2), jrk6(_jrk6), time(_time)
	{
		posH.x = float2(_pos[0]).x;
		posL.x = float2(_pos[0]).y;
		posH.y = float2(_pos[1]).x;
		posL.y = float2(_pos[1]).y;
		posH.z = float2(_pos[2]).x;
		posL.z = float2(_pos[2]).y;
		posH.w = float(_mass);
	}

    template <int rw, int locality>
	void prefetch() const {
		__builtin_prefetch(&posH, rw , locality);
		__builtin_prefetch(&jrk6, rw , locality);
	}
	bool is_aligned() const {
		unsigned long addr = (unsigned long)(&posH);
		return (0 == addr%32);
	}
	void sstore(Particle *pdst) const {
		const v4sf *src = (const v4sf *)&posH;
		const v4sf xm0 = src[0];
		const v4sf xm1 = src[1];
		const v4sf xm2 = src[2];
		const v4sf xm3 = src[3];
		const v4sf xm4 = src[4];
		const v2df xm5 = *(v2df *)&time;
		v4sf *dst = (v4sf *)pdst;
		__builtin_ia32_movntps((float *)(dst + 0), xm0);
		__builtin_ia32_movntps((float *)(dst + 1), xm1);
		__builtin_ia32_movntps((float *)(dst + 2), xm2);
		__builtin_ia32_movntps((float *)(dst + 3), xm3);
		__builtin_ia32_movntps((float *)(dst + 4), xm4);
		__builtin_ia32_movntps((float *)(dst + 5), (v4sf)xm5);
	}
};
#else
struct Particle{
	v4sf posH; // elemen w is for mass
	v4sf posL;
	v4sf vel;
	v4sf acc2;
	v4sf jrk6;
	double time, pad; // 6 xmm words;

	void make_pos(const double *x, const double *m){
		v2df xy = {x[0], x[1]};
		v2df zw = {x[2], m[0]};
		v4sf xyH = __builtin_ia32_cvtpd2ps(xy);
		v4sf xyL = __builtin_ia32_cvtpd2ps(xy - __builtin_ia32_cvtps2pd(xyH));
		v4sf zwH = __builtin_ia32_cvtpd2ps(zw);
		v4sf zwL = __builtin_ia32_cvtpd2ps(zw - __builtin_ia32_cvtps2pd(zwH));

		posH = __builtin_ia32_movlhps(xyH, zwH);
		posL = __builtin_ia32_movlhps(xyL, zwL);
	}
	static v4sf make_v4sf(const double *x){
		v2df xy = {x[0], x[1]};
		v2df zw = {x[2],     };
		return __builtin_ia32_movlhps(
				__builtin_ia32_cvtpd2ps(xy),
				__builtin_ia32_cvtpd2ps(zw));
	}

	Particle(
		const double _pos [3],
		const double _vel [3],
		const double _acc2[3],
		const double _jrk6[3],
		const double _mass,
		const double _time)
		: time(_time)
	{
		make_pos(_pos, &_mass);
		vel  = make_v4sf(_vel);
		acc2 = make_v4sf(_acc2);
		jrk6 = make_v4sf(_jrk6);
	}

    template <int rw, int locality>
	void prefetch() const {
		__builtin_prefetch(&posH, rw , locality);
		__builtin_prefetch(&jrk6, rw , locality);
	}
	bool is_aligned() const {
		unsigned long addr = (unsigned long)(&posH);
		return (0 == addr%32);
	}
	void sstore(Particle *pdst) const {
		const v4sf xm0 = posH;
		const v4sf xm1 = posL;
		const v4sf xm2 = vel;
		const v4sf xm3 = acc2;
		const v4sf xm4 = jrk6;
		const v2df xm5 = *(v2df *)&time;
		v4sf *dst = (v4sf *)pdst;
		__builtin_ia32_movntps((float *)(dst + 0), xm0);
		__builtin_ia32_movntps((float *)(dst + 1), xm1);
		__builtin_ia32_movntps((float *)(dst + 2), xm2);
		__builtin_ia32_movntps((float *)(dst + 3), xm3);
		__builtin_ia32_movntps((float *)(dst + 4), xm4);
		__builtin_ia32_movntps((float *)(dst + 5), (v4sf)xm5);
	}
	void store(Particle *pdst) const {
		const v4sf xm0 = posH;
		const v4sf xm1 = posL;
		const v4sf xm2 = vel;
		const v4sf xm3 = acc2;
		const v4sf xm4 = jrk6;
		const v2df xm5 = *(v2df *)&time;
		v4sf *dst = (v4sf *)pdst;
		dst[0] = xm0;
		dst[1] = xm1;
		dst[2] = xm2;
		dst[3] = xm3;
		dst[4] = xm4;
		dst[5] = (v4sf)xm5;
	}
};
#endif

struct Predictor{
	v4sf posH; // elemen w is for mass
	v4sf posL;
	v4sf vel;  // 3 xmm words

	Predictor(const int){ // generate dummy
		posH = (v4sf){255.0f, 255.0f, 255.0f, 0.0f};
		posL = vel = (v4sf){0.f, 0.f, 0.f, 0.f};
	}
	Predictor(
			const Particle &p,
			const double ti)
	{
		const float dt = float(ti - p.time);
		const v4sf s0 = {dt, dt, dt, dt};
/*
		const v4sf s1 = s0 + s0;
		const v4sf s2 = s0 * (v4sf){1.5f, 1.5f, 1.5f, 1.5f};

		posH = p.posH;
		posL = v4sf(p.posL) + s0*(v4sf(p.vel) + s0*(v4sf(p.acc2) + s0*(v4sf(p.jrk6))));
		vel = v4sf(p.vel) + s1*(v4sf(p.acc2) + s2*(v4sf(p.jrk6)));
*/
		const v4sf s2 = s0 + s0;
		const v4sf s1 = s0 * (v4sf){ 0.5f, 0.5f, 0.5f, 0.5f };

		posH = p.posH;
		posL = v4sf(p.posL) + s0*(v4sf(p.vel) + s1 * (v4sf(p.acc2) + s2 * v4sf(p.jrk6)));
		vel = v4sf(p.vel) + s0*(v4sf(p.acc2) + s1 * v4sf(p.jrk6));

	}

    template<int rw, int locality>
	void prefetch() const {
		__builtin_prefetch(&posH, rw, locality);
		__builtin_prefetch(&vel,  rw, locality);
	}
	void store(Predictor *dst) const {
		dst->posH = posH;
		dst->posL = posL;
		dst->vel  = vel ;
	}
	void sstore(Predictor *dst) const {
		__builtin_ia32_movntps((float*)dst + 0, posH);;
		__builtin_ia32_movntps((float*)dst + 4, posL);;
		__builtin_ia32_movntps((float*)dst + 8, vel );;
	}
};

struct NBlist{
	// enum{ NB_MAX = 512 + 88 };
	enum{ NB_MAX = 699 + 128}; // 2015 Oct
	int nnb;
	int nb[NB_MAX];

	void print(const int i, FILE *fp = stdout) const{
		fprintf(fp, "%6d%6d :", i, nnb);
		for(int k=0; k<nnb; k++){
			fprintf(fp, " %d", nb[k]);
		}
		fprintf(fp, "\n");
		fflush(fp);
	}

	NBlist(){
		nnb = 0;
		for(int i=0; i<NB_MAX; i++){
			nb[i] = 0;
		}
	}

	void prefetch() const{
		const char *ptr = (const char *)&nnb;
		__builtin_prefetch(ptr, 0, 0);
		asm volatile ("#bar");
		const int nnb = this->nnb;
		// const int nnb = ((NBlist * volatile)(this))->nnb;
		ptr += 64;
		for(int i=15; i<nnb; i+=16, ptr+=64){
			__builtin_prefetch(ptr, 0, 0);
		}
	}
};

static bool       is_open = false;
static Particle  *ptcl = NULL;
static Predictor *pred = NULL;
// static std::vector<NBlist> list;
static NBlist    *list = NULL;
static int G_nmax;

static double time_pred, time_pact, time_grav, time_onep, time_push, time_sort;
static unsigned long long num_inter, num_fcall, num_steps, num_pred_all, num_pred_act;
static int num_threads;

extern "C" void cflush(){
	for(double *itr = (double *)(ptcl); itr < (double *)(ptcl + G_nmax); itr += 8){
		__builtin_ia32_clflush(itr);
	}
	for(double *itr = (double *)(pred); itr < (double *)(pred + G_nmax); itr += 8){
		__builtin_ia32_clflush(itr);
	}
}

static void init_affinity(){
#ifdef CORE_I7_AFFINITY
#pragma omp parallel
	{
		const int th2core[8] = {0, 4, 1, 5, 2, 6, 3, 7};
		const int tid = omp_get_thread_num();
		const int core = th2core[tid];
		cpu_set_t set;
		if(sched_getaffinity(0, sizeof(cpu_set_t), &set)) fprintf(stderr, "error in getting affinity\n");
		CPU_ZERO(&set);
		CPU_SET(core, &set);
		if(sched_setaffinity(0, sizeof(cpu_set_t), &set)) fprintf(stderr, "error in setting affinity\n");
	}
#endif
}

extern "C" int cudaMallocHost(void **, size_t);
extern "C" int cudaFreeHost(void *);

static void gpuirr_open(
		const int nmax,
		const int lmax)
{
	if(is_open){
		fprintf(stderr, "gpuirr: it is already open\n");
		return;
	}

	assert(lmax <= 1 + NBlist::NB_MAX);

	fprintf(stderr, "**************************** \n");
	fprintf(stderr, "Opening GPUIRR lib. SSE ver. \n");
	fprintf(stderr, "  (debug ver. at 2013/11/15) \n");
	fprintf(stderr, " nmax = %d, lmax = %d\n", nmax, lmax);
	fprintf(stderr, "**************************** \n");

	init_affinity();
	// ptcl.resize(nmax);
	// pred.resize(nmax);
	void *ptr;
	assert(0 == posix_memalign(&ptr, 64, nmax * sizeof(Particle )));
	ptcl = (Particle *)ptr;
#if 0
	float *f = (float *)ptr;
	for(int i=0; i<nmax*24; i++){
		f[i] = drand48();
	}
#endif

	assert(0 == posix_memalign(&ptr, 64, (1+nmax) * sizeof(Predictor)));
	pred = (Predictor *)ptr;
	// store dummy predictor to the last of the array
	pred[nmax] = Predictor(0);

	// list.resize(nmax);
	assert(0 == posix_memalign(&ptr, 64, nmax * sizeof(NBlist)));
	list = (NBlist *)ptr;

	G_nmax = nmax;
#pragma omp parallel
	{
		num_threads = omp_get_num_threads();
	}
	time_pred = time_pact =  time_grav = time_onep = time_push = time_sort = 0.0;
	num_inter = num_fcall = num_steps = num_pred_all = num_pred_act = 0;

	is_open = true;
}

static void gpuirr_close(){
	if(!is_open){
		fprintf(stderr, "gpuirr: it is already close\n");
		return;
	}
#if 0
	cudaFreeHost(ptcl); ptcl = NULL;
	cudaFreeHost(pred); pred = NULL;
#else
	free(pred); pred = NULL;
	free(ptcl); ptcl = NULL;
	free(list); list = NULL;
#endif

	const double Gflops = 60.0 * double(num_inter) * 1.e-9 / time_grav;
	const double usec_fcall    = 1.e6 * (time_grav / num_fcall);
	const double nsec_pred_all = 1.e9 * (time_pred / num_pred_all);
	const double nsec_pred_act = 1.e9 * (time_pact / num_pred_act);
	const double nnb_avr = double(num_inter) / double(num_steps);
	const double time_tot = time_pred + time_pact +  time_grav + time_onep + time_push + time_sort;

	fprintf(stderr, "**************************** \n");
	fprintf(stderr, "Closing GPUIRR lib. CPU ver. \n");
	fprintf(stderr, "time pred  : %f sec\n", time_pred);
	fprintf(stderr, "time pact  : %f sec\n", time_pact);
	fprintf(stderr, "time grav  : %f sec\n", time_grav);
	fprintf(stderr, "time onep  : %f sec\n", time_onep);
	fprintf(stderr, "time push  : %f sec\n", time_push);
	fprintf(stderr, "time sort  : %f sec\n", time_sort);
	fprintf(stderr, "-------------------\n");
	fprintf(stderr, "time total : %f sec\n", time_tot );
	fprintf(stderr, "\n");
	fprintf(stderr, "perf grav  : %f Gflops\n", Gflops);
	fprintf(stderr, "perf grav  : %f usec\n", usec_fcall);
	fprintf(stderr, "perf pred  : %f nsec\n", nsec_pred_all);
	fprintf(stderr, "perf pact  : %f nsec\n", nsec_pred_act);
	fprintf(stderr, "<#NB>      : %f \n",     nnb_avr);
	fprintf(stderr, "**************************** \n");

	is_open = false;
}

static void gpuirr_set_jp(
		const int addr,
		const double pos [3],
		const double vel [3],
		const double acc2[3],
		const double jrk6[3],
		const double mass,
		const double time)
{
	// ptcl[addr] = Particle(pos, vel, acc2, jrk6, mass, time);
	// Particle(pos, vel, acc2, jrk6, mass, time).sstore(&ptcl[addr]);
	Particle(pos, vel, acc2, jrk6, mass, time).store(&ptcl[addr]);
}

static void gpuirr_set_list(
		const int addr,
		const int nnb,
		const int nblist[])
{
	assert(nnb <= NBlist::NB_MAX);
	list[addr].nnb = nnb;
#if 0
	for(int k=0; k<nnb; k++){
		list[addr].nb[k] = nblist[k] - 1;
	}
#else
	const int *src = nblist;
	      int *dst = list[addr].nb;
	for(int k=0; k<nnb; k+=4){
		const int i0 = src[k+0] - 1;
		const int i1 = src[k+1] - 1;
		const int i2 = src[k+2] - 1;
		const int i3 = src[k+3] - 1;
		dst[k+0] = i0;
		dst[k+1] = i1;
		dst[k+2] = i2;
		dst[k+3] = i3;
	}
	// fill dummy
	const int nmax = G_nmax;
	const int kmax = 4 + 4 * (1 + (nnb-1)/4); // 4 more dummies for prefetching
	for(int k=nnb; k<kmax; k++){
		dst[k] = nmax;
	}
#endif
}

static void gpuirr_get_list(
		const int addr,
		_out_ int &nnb,
		_out_ int nblist[])
{
	nnb = list[addr].nnb;
	assert(nnb <= NBlist::NB_MAX);
	const int n = nnb;
	for(int i=0; i<n; i++){
		nblist[i] = list[addr].nb[i] + 1;
	}
}

static void gpuirr_pred_all(
		const int    js,
		const int    je,
		const double ti)
{
	const int nmax = G_nmax;
	assert(js >= 0);
	assert(je <= nmax);
	Particle  *ptcl = ::ptcl;
	Predictor *pred = ::pred;
	const double t0 = get_wtime();
#if 1
#pragma omp parallel for
	for(int j=js; j<je; j+=2){
		char *ptr = (char *)&ptcl[j+2];
		__builtin_prefetch(ptr + 0,   0, 0);
		__builtin_prefetch(ptr + 64,  0, 0);
		__builtin_prefetch(ptr + 128, 0, 0);
		// asm volatile ("#bar");
		const Predictor pr0(ptcl[j+0], ti);
		const Predictor pr1(ptcl[j+1], ti);
		pr0.sstore(&pred[j+0]);
		pr1.sstore(&pred[j+1]);
	}
#else
#pragma omp parallel for
	for(int j=js; j<je; j+=4){
		char *ptr = (char *)&ptcl[j+4];
		__builtin_prefetch(ptr + 0,   0, 0);
		__builtin_prefetch(ptr + 64,  0, 0);
		__builtin_prefetch(ptr + 128, 0, 0);
		__builtin_prefetch(ptr + 192,   0, 0);
		__builtin_prefetch(ptr + 256,  0, 0);
		__builtin_prefetch(ptr + 320, 0, 0);
		const Predictor pr0(ptcl[j+0], ti);
		const Predictor pr1(ptcl[j+1], ti);
		const Predictor pr2(ptcl[j+2], ti);
		const Predictor pr3(ptcl[j+3], ti);
		pr0.sstore(&pred[j+0]);
		pr1.sstore(&pred[j+1]);
		pr2.sstore(&pred[j+2]);
		pr3.sstore(&pred[j+3]);
	}
#endif
	// store dummy predictor to the last of the array
	Predictor(0).sstore(&pred[nmax]);;
	const double t1 = get_wtime();
	::time_pred += t1-t0;
	::num_pred_all += (je - js);
}

static void gpuirr_pred_act(
		const int ni,
		const int addr[],
		const double ti)
{
#if 0
	static std::vector<int> slist(1024), ulist(1024);
	const double t0 = get_wtime();
	slist.clear();
	ulist.clear();
	asm("#push_back");
#if 0
	for(int i=0; i<ni; i++){
		const int ii = addr[i]-1;
		slist.push_back(ii);
		for(int k=0; k<list[ii].nnb; k++){
			slist.push_back(list[ii].nb[k]);
		}
	}
#else
	int length = 0;
	for(int i=0; i<ni; i++){
		const int ii = addr[i]-1;
		length += 1 + list[ii].nnb;
	}
	slist.resize(length);
	int *itr = slist.data();
	for(int i=0; i<ni; i++){
		const int ii = addr[i]-1;
		*itr++ = ii;
		const int nnb = list[ii].nnb;
		const int *nb = list[ii].nb;
		for(int k=0; k<nnb; k++){
			*itr++ = nb[k];
		}
	}
#endif
	asm("#sort");
	std::sort(slist.begin(), slist.end());
	ulist.resize( slist.size() );
	asm("#unique_copy");
	std::vector<int>::iterator it
		= std::unique_copy(slist.begin(), slist.end(), ulist.begin());
	ulist.erase(it, ulist.end());

	const int npred = ulist.size();
	const double t1 = get_wtime();

	if(npred < (1<<8)){
#if 1
		for(int i=0; i<npred; i++){
			const int j = ulist[i];
			ptcl[j].prefetch<0,1>();
			pred[j].prefetch<0,1>();
		}
#endif
		for(int i=0; i<npred; i++){
			const int j = ulist[i];
			Predictor(ptcl[j], ti).store(&pred[j]);
		}
	}else{
#pragma omp parallel for
		for(int i=0; i<npred; i++){
			const int j = ulist[i];
			Predictor(ptcl[j], ti).sstore(&pred[j]);
		}
	}
	const double t2 = get_wtime();
	::time_sort += t1-t0;
	::time_pact += t2-t1;
	::num_pred_act += npred;
#else
	struct Bin{
		enum{
			NBIN = 16,
			NBAC = 16,
		};
		std::vector<int> item;
		std::vector<int> backet[NBAC];
		int pad[16]; // avoid false sharing

		Bin(){
			item.reserve(1024);
			for(int i=0; i<NBAC; i++){
				backet[i].reserve(64);
			}
		}
#if 0
		void sort(){
			asm("#sort");
			std::sort(item.begin(), item.end());
			asm("#unique");
			item.erase(
					std::unique(item.begin(), item.end()),
					item.end());
		}
#endif
		void sort(const Particle ptcl[], const Predictor pred[]){
			asm("#unique_sort");
			const int len = item.size();
			for(int i=0; i<len; i++){
				const int ii = item[i];
				const unsigned key = (unsigned(ii) / NBIN) % NBAC;
				if(std::find(backet[key].begin(), backet[key].end(), ii) == backet[key].end()){
					ptcl[ii].prefetch<0,0>();
					pred[ii].prefetch<1, 3>();
					backet[key].push_back(ii);
				}
			}
			asm("#write_bacwrite_back");
			std::vector<int>::iterator it = item.begin();
			for(int k=0; k<NBAC; k++){
				const int len = backet[k].size();
				for(int i=0; i<len; i++){
					*it++ = backet[k][i];
				}
				backet[k].clear();
			}
			item.erase(it, item.end());
		}
		int size(){
			return item.size();
		}
		void predict(const double ti, const Particle ptcl[], Predictor pred[]){
#if 0
			for(std::vector<int>::iterator it = item.begin(); it != item.end(); ++it){
				const int j = *it;
				ptcl[j].prefetch<0,0>();
			}
#endif
			for(std::vector<int>::iterator it = item.begin(); it != item.end(); ++it){
				const int j = *it;
				Predictor(ptcl[j], ti).store(&pred[j]);
			}
		}
	};

	const double t0 = get_wtime();
	const int NBIN = Bin::NBIN;
	static Bin bin[NBIN];
	for(int ibin=0; ibin<NBIN; ibin++){
		bin[ibin].item.clear();
	}
	asm("#push_back");
	for(int i=0; i<ni; i++){
		const int ii    = addr[i]-1;
		// const int inext = addr[i+1]-1;
		// list[inext].prefetch<0,1>();
		const int ibin = unsigned(ii) % NBIN;
		bin[ibin].item.push_back(ii);
		const int nnb = list[ii].nnb;
		const int *nb = list[ii].nb;
		for(int k=0; k<nnb; k++){
			const int iii = nb[k];
			const int ibin = unsigned(iii) % NBIN;
			bin[ibin].item.push_back(iii);
		}
	}
	const double t1 = get_wtime();
#pragma omp parallel for
	for(int ibin=0; ibin<NBIN; ibin++){
		// bin[ibin].sort();
		bin[ibin].sort(::ptcl, ::pred);
	}
	for(int ibin=0; ibin<NBIN; ibin++){
		num_pred_act += bin[ibin].size();
	}
	const double t2 = get_wtime();
#pragma omp parallel for
	for(int ibin=0; ibin<NBIN; ibin++){
		bin[ibin].predict(ti, ::ptcl, ::pred);
	}
	const double t3 = get_wtime();

	::time_push += t1-t0;
	::time_sort += t2-t1;
	::time_pact += t3-t2;
#endif
}

struct DPAccum{
	v2df hi, lo;
	DPAccum() : hi((v2df){0.0, 0.0}), lo((v2df){0.0, 0.0}) {}
	void accumulate(const v4sf data){
		lo += __builtin_ia32_cvtps2pd(data);
		hi += __builtin_ia32_cvtps2pd(
				__builtin_ia32_movhlps(data, data));
	}
	double reduce() const {
		v2df tmpd = lo + hi;
		v2df sum  = __builtin_ia32_haddpd(tmpd, tmpd);
		return __builtin_ia32_vec_ext_v2df(sum , 0);
	}
};

struct SPAccum{
	v4sf accum;

	SPAccum() : accum((v4sf){0.f, 0.f, 0.f, 0.f}) {}
	void accumulate(const v4sf data){
		accum += data;
	}
	double reduce() const {
		v4sf tmpf = __builtin_ia32_haddps(accum, accum);
		v2df tmpd = __builtin_ia32_cvtps2pd(tmpf);
		v2df sum  = __builtin_ia32_haddpd(tmpd, tmpd);
		return __builtin_ia32_vec_ext_v2df(sum , 0);
	}
};

struct Pred4{
	v4sf xH, yH, zH, mass;
	v4sf xL, yL, zL;
	v4sf vx, vy, vz;

	static v4sf bcast(const float f){
		v4sf v = {f, f, f, f};
		return v;
	}
	static v4sf bcast0(const v4sf v){
		return __builtin_ia32_shufps(v, v, 0);
	}
	static v4sf bcast1(const v4sf v){
		return __builtin_ia32_shufps(v, v, 0x55);
	}
	static v4sf bcast2(const v4sf v){
		return __builtin_ia32_shufps(v, v, 0xaa);
	}
	static v4sf bcast3(const v4sf v){
		return __builtin_ia32_shufps(v, v, 0xff);
	}
	static void transpose(v4sf &v0, v4sf &v1, v4sf &v2, v4sf &v3){
		const v4sf t0 = __builtin_ia32_unpcklps(v0, v2);
		const v4sf t1 = __builtin_ia32_unpckhps(v0, v2);
		const v4sf t2 = __builtin_ia32_unpcklps(v1, v3);
		const v4sf t3 = __builtin_ia32_unpckhps(v1, v3);

		v0 = __builtin_ia32_unpcklps(t0, t2);
		v1 = __builtin_ia32_unpckhps(t0, t2);
		v2 = __builtin_ia32_unpcklps(t1, t3);
		v3 = __builtin_ia32_unpckhps(t1, t3);
	}

	// for the i-particle
	Pred4(const Predictor &pi){
		xH   = bcast0(pi.posH);
		yH   = bcast1(pi.posH);
		zH   = bcast2(pi.posH);
		mass = bcast3(pi.posH);
		xL   = bcast0(pi.posL);
		yL   = bcast1(pi.posL);
		zL   = bcast2(pi.posL);
		vx   = bcast0(pi.vel );
		vy   = bcast1(pi.vel );
		vz   = bcast2(pi.vel );
	}
	// for the j-particle
	Pred4(const Predictor &p0, const Predictor &p1, const Predictor &p2, const Predictor &p3){
		xH   = p0.posH;
		yH   = p1.posH;
		zH   = p2.posH;
		mass = p3.posH;
		transpose(xH, yH, zH, mass);
		v4sf dum0;
		xL   = p0.posL;
		yL   = p1.posL;
		zL   = p2.posL;
		dum0 = p3.posL;
		transpose(xL, yL, zL, dum0);
		v4sf dum1;
		vx   = p0.vel;
		vy   = p1.vel;
		vz   = p2.vel;
		dum1 = p3.vel;
		transpose(vx, vy, vz, dum1);
	}
};

static inline void v4sf_print(const v4sf v){
	printf("%16.6e %16.6e %16.6e %16.6e ",
			__builtin_ia32_vec_ext_v4sf(v, 0),
			__builtin_ia32_vec_ext_v4sf(v, 1),
			__builtin_ia32_vec_ext_v4sf(v, 2),
			__builtin_ia32_vec_ext_v4sf(v, 3));
}

static void gpuirr_firr(
		const int addr,
		_out_ double accout[3],
		_out_ double jrkout[3])
{
	// list[addr].prefetch<0,1>();
	const Predictor *pred = ::pred;
	const int *jptr = &(list[addr].nb[0]);
	const Predictor *jp0 = &pred[jptr[0]];
	const Predictor *jp1 = &pred[jptr[1]];
	const Predictor *jp2 = &pred[jptr[2]];
	const Predictor *jp3 = &pred[jptr[3]];
	jp0->prefetch<0, 0>();
	jp1->prefetch<0, 0>();
	jp2->prefetch<0, 0>();
	jp3->prefetch<0, 0>();

	const int nnb = list[addr].nnb;
	const Pred4 ip(pred[addr]);

	SPAccum ax, ay, az;
	SPAccum jx, jy, jz;
#if 0
	for(int k=0; k<nnb; k++){
		const int j = list[addr].nb[k];
		pred[j].prefetch<0,1>();
	}
#endif

	for(int k=0; k<nnb; k+=4){
		const int *jptr = &(list[addr].nb[k]);
#if 0
		pred[jptr[4]].prefetch<0,0>();
		pred[jptr[5]].prefetch<0,0>();
		pred[jptr[6]].prefetch<0,0>();
		pred[jptr[7]].prefetch<0,0>();
		const Predictor &jp0 = pred[jptr[0]];
		const Predictor &jp1 = pred[jptr[1]];
		const Predictor &jp2 = pred[jptr[2]];
		const Predictor &jp3 = pred[jptr[3]];
		const Pred4 jp(jp0, jp1, jp2, jp3);
#else
		const Predictor *jj0 = &pred[jptr[4]];
		const Predictor *jj1 = &pred[jptr[5]];
		const Predictor *jj2 = &pred[jptr[6]];
		const Predictor *jj3 = &pred[jptr[7]];
		jj0->prefetch<0, 0>();
		jj1->prefetch<0, 0>();
		jj2->prefetch<0, 0>();
		jj3->prefetch<0, 0>();
		const Pred4 jp(*jp0, *jp1, *jp2, *jp3);
		jp0 = jj0;
		jp1 = jj1;
		jp2 = jj2;
		jp3 = jj3;
#endif

		const v4sf dx = (jp.xH - ip.xH) + (jp.xL - ip.xL);
		const v4sf dy = (jp.yH - ip.yH) + (jp.yL - ip.yL);
		const v4sf dz = (jp.zH - ip.zH) + (jp.zL - ip.zL);

		const v4sf dvx = jp.vx - ip.vx;
		const v4sf dvy = jp.vy - ip.vy;
		const v4sf dvz = jp.vz - ip.vz;

		const v4sf r2 = dx*dx + dy*dy + dz*dz;
#if 0
		v4sf_print(jp.xH);
		v4sf_print(jp.yH);
		v4sf_print(jp.zH);
		printf("\n");
#endif
		const v4sf rv = dx*dvx + dy*dvy + dz*dvz;
		const v4sf rinv   = rsqrt_NR(r2);
		const v4sf rinv2  = rinv * rinv;
		const v4sf c1     = {-3.0f, -3.0f, -3.0f, -3.0f};
		const v4sf alpha  = c1 * rinv2 * rv;
		const v4sf mrinv3 = jp.mass * rinv * rinv2;

		ax.accumulate(mrinv3 * dx);
		ay.accumulate(mrinv3 * dy);
		az.accumulate(mrinv3 * dz);
#if 0
		jx.accumulate(mrinv3 * (dvx + alpha * dx));
		jy.accumulate(mrinv3 * (dvy + alpha * dy));
		jz.accumulate(mrinv3 * (dvz + alpha * dz));
#else
		jx.accumulate(mrinv3 * dvx); jx.accumulate(alpha * (mrinv3 * dx));
		jy.accumulate(mrinv3 * dvy); jy.accumulate(alpha * (mrinv3 * dy));
		jz.accumulate(mrinv3 * dvz); jz.accumulate(alpha * (mrinv3 * dz));
#endif
	}

	accout[0] = ax.reduce();
	accout[1] = ay.reduce();
	accout[2] = az.reduce();
	jrkout[0] = jx.reduce();
	jrkout[1] = jy.reduce();
	jrkout[2] = jz.reduce();
}

static void gpuirr_firr_vec(
		const int ni,
		const int addr[],
		_out_ double accout[][3],
		_out_ double jrkout[][3])
{
	const double t0 = get_wtime();
#if 1
	int nint = 0;
#pragma omp parallel for reduction(+: nint) schedule(guided)
	for(int i=0; i<ni; i++){
		gpuirr_firr(addr[i]-1, accout[i], jrkout[i]);
		nint += list[addr[i]-1].nnb;
	}
	::num_inter += nint;
#else
#pragma omp parallel
	{
		const int tid = omp_get_thread_num();
		const int nth = omp_get_num_threads();
		for(int i=0; i<ni; i++){
			const int iaddr = addr[i] - 1;
			if(tid == iaddr % nth){
				gpuirr_firr(iaddr, accout[i], jrkout[i]);
			}
			if(0 == tid){
				::num_inter += list[iaddr].nnb;
			}
		}
	}
#endif
	const double t1 = get_wtime();
	::time_grav += t1-t0;
	::num_fcall++;
	::num_steps += ni;

#if 0
	for(int i=0; i<ni; i++){
		if(fabs(accout[i][0]) > 1.0e6){
			fprintf(stdout, "big force on %d, (%e, %e, %e)\n",
					addr[i], accout[i][0], accout[i][1], accout[i][2]);
		}
	}
#endif
}

extern "C"{
	void gpuirr_open_(int *nmax, int *lmax){
		gpuirr_open(*nmax, *lmax);
	}
	void gpuirr_close_(){
		gpuirr_close();
	}
	void gpuirr_set_jp_(
		int    *addr,
		double  pos [3],
		double  vel [3],
		double  acc2[3],
		double  jrk6[3],
		double *mass,
		double *time)
	{
		gpuirr_set_jp((*addr)-1, pos, vel, acc2, jrk6, *mass, *time);
	}
	void gpuirr_set_list_(
		int *addr,
		int *nblist)
	{
		gpuirr_set_list((*addr)-1, *nblist, nblist+1);
	}
	void gpuirr_get_list_(
		int *addr,
		int *nblist)
	{
		gpuirr_get_list((*addr)-1, *nblist, nblist+1);
	}
	void gpuirr_pred_all_(
			int    *js,
			int    *je,
			double *ti)
	{
		gpuirr_pred_all((*js)-1, *(je),  *ti);
	}
	void gpuirr_pred_act_(
			int    *ni,
			int     addr[],
			double *ti)
	{
		gpuirr_pred_act(*ni, addr, *ti);
	}
#if 1
	void gpuirr_firr_(
			int    *addr,
			double  acc[3],
			double  jrk[3])
	{
		const double t0 = get_wtime();
		gpuirr_firr((*addr)-1, acc, jrk);
		const double t1 = get_wtime();
		::time_onep += t1-t0;
	}
	void gpuirr_firr_vec_(
			int   *ni,
			int    addr[],
			double acc [][3],
			double jrk [][3])
	{
		gpuirr_firr_vec(*ni, addr, acc, jrk);
	}
	void gpuirr_debug_pair_(
			const int *_i,
			const int *_j)
	{
		const int i = *_i - 1;
		const int j = *_j - 1;
		const Predictor *pred = ::pred;
		v4sf dr = (pred[j].posH - pred[i].posH)
		        + (pred[j].posL - pred[i].posL);
		union{
			v4sf  v;
			float f[4];
		} m128;
		m128.v = dr;
		fprintf(stdout, "x[%d] - x[%d] = (%e, %e, %e)\n",
				i+1, j+1, m128.f[0], m128.f[1], m128.f[2]);
		fflush(stdout);
	}
# if 0
	void gpuirr_debug_shoot_j_(
			const int *_i,
			int       *jout,
			double    *rmin)
	{
		const int i = *_i - 1;
		const Predictor *pred = ::pred;
		const int nnb = list[i].nnb;
		int jmin = -1;
		float r2min = HUGE;
		for(int k=0; k<nnb; k++){
			const int j = list[i].nb[k];
			v4sf dr = (pred[j].posH - pred[i].posH)
					+ (pred[j].posL - pred[i].posL);
			union{
				v4sf  v;
				float f[4];
			} m128;
			m128.v = dr*dr;
			const float r2 = m128.f[0] + m128.f[1] + m128.f[2];
			if(r2 < r2min){
				r2min = r2;
				jmin = j;
			}
		}
		*jout = jmin + 1;
		*rmin = sqrt((double)r2min);
	}
# endif
#endif
}
