// #include "support_kernels.cu"
#include <stdio.h>
#include "../profiling/bonsai_timing.h"
PROF_MODULE(dev_approximate_gravity);

#include "node_specs.h"

#ifdef WIN32
#define M_PI        3.14159265358979323846264338328
#endif

__forceinline__ __device__ float Wkernel(const float q)
{
  const float sigma = 8.0f/M_PI;

  const float qm = 1.0f - q;
  const float f1 = sigma * (1.0f + (-6.0f)*q*q*qm);
  const float f2 = sigma * 2.0f*qm*qm*qm;

  return fmaxf(0.0f, fminf(f1, f2));
}

__forceinline__ __device__ float interact(
    const float3 ipos,
    const float  h,
    const float  hinv,
    const float3 jpos,
    const float  jmass)
{
  const float3 dr = make_float3(jpos.x - ipos.x, jpos.y - ipos.y, jpos.z - ipos.z);
  const float  r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
  if (r2 >= h*h) return 0.0f;
  const float q  = sqrtf(r2) * hinv;
  const float hinv3 = hinv*hinv*hinv;

  return jmass * Wkernel(q) * hinv3;
}


/***
**** --> prefix calculation via Horn(2005) data-parallel algoritm
***/
#define BTEST(x) (-(int)(x))
template<int DIM2>
__device__ int calc_prefix(int N, int* prefix_in, int tid) {
  int x, y = 0;

  const int DIM = 1 << DIM2;
  
  for (int p = 0; p < N; p += DIM) {
    int *prefix = &prefix_in[p];

    x = prefix[tid -  1]; __syncthreads(); prefix[tid] += x & BTEST(tid >=  1); __syncthreads();
    x = prefix[tid -  2]; __syncthreads(); prefix[tid] += x & BTEST(tid >=  2); __syncthreads();
    x = prefix[tid -  4]; __syncthreads(); prefix[tid] += x & BTEST(tid >=  4); __syncthreads();
    x = prefix[tid -  8]; __syncthreads(); prefix[tid] += x & BTEST(tid >=  8); __syncthreads();
    x = prefix[tid - 16]; __syncthreads(); prefix[tid] += x & BTEST(tid >= 16); __syncthreads();
    if (DIM2 >= 6) {x = prefix[tid - 32]; __syncthreads(); prefix[tid] += x & BTEST(tid >= 32); __syncthreads();}
    if (DIM2 >= 7) {x = prefix[tid - 64]; __syncthreads(); prefix[tid] += x & BTEST(tid >= 64); __syncthreads();}
    if (DIM2 >= 8) {x = prefix[tid -128]; __syncthreads(); prefix[tid] += x & BTEST(tid >=128); __syncthreads();}
    

    prefix[tid] += y;
    __syncthreads();

    y = prefix[DIM-1];
    __syncthreads();
  }

  return y;
} 

template<int DIM2>
__device__ int calc_prefix(int* prefix, int tid, int value) {
  int  x;
  
  const int DIM = 1 << DIM2;

  prefix[tid] = value;
  __syncthreads();

#if 1
  x = prefix[tid -  1]; __syncthreads(); prefix[tid] += x & BTEST(tid >=  1); __syncthreads();
  x = prefix[tid -  2]; __syncthreads(); prefix[tid] += x & BTEST(tid >=  2); __syncthreads();
  x = prefix[tid -  4]; __syncthreads(); prefix[tid] += x & BTEST(tid >=  4); __syncthreads();
  x = prefix[tid -  8]; __syncthreads(); prefix[tid] += x & BTEST(tid >=  8); __syncthreads();
  x = prefix[tid - 16]; __syncthreads(); prefix[tid] += x & BTEST(tid >= 16); __syncthreads();
  if (DIM2 >= 6) {x = prefix[tid - 32]; __syncthreads(); prefix[tid] += x & BTEST(tid >= 32); __syncthreads();}
  if (DIM2 >= 7) {x = prefix[tid - 64]; __syncthreads(); prefix[tid] += x & BTEST(tid >= 64); __syncthreads();}
  if (DIM2 >= 8) {x = prefix[tid -128]; __syncthreads(); prefix[tid] += x & BTEST(tid >=128); __syncthreads();}

  x = prefix[DIM - 1];
  __syncthreads();
  return x;
#else
  
  int offset = 0;
  int tid2 = tid << 1;

#pragma unroll
  for (int d = DIM >> 1; d > 0; d >>= 1) {
    __syncthreads();

    int iflag = BTEST(tid < d);
    int ai = (((tid2 + 1) << offset) - 1) & iflag;
    int bi = (((tid2 + 2) << offset) - 1) & iflag;
    
    prefix[bi] += prefix[ai] & iflag;
    offset++;
  }

  // clear the last element
  if (tid == 0) prefix[DIM - 1] = 0;

  // traverse down the tree building the scan in place
#pragma unroll
  for (int d = 1; d < DIM; d <<= 1) {
    offset--;
    __syncthreads();
    
    int iflag = BTEST(tid < d);
    int ai = (((tid2 + 1) << offset) - 1) & iflag;
    int bi = (((tid2 + 2) << offset) - 1) & iflag;
    
    int t       = prefix[ai];
    if (tid < d) {
      prefix[ai]  = (prefix[bi] & iflag) + (t & BTEST(tid >= d));
      prefix[bi] += t & iflag;
    }
  }
  __syncthreads();

  prefix[tid] += value;
  __syncthreads();
  
  x = prefix[DIM - 1];
  __syncthreads();
  return x;
#endif
}

template<int SHIFT>
__forceinline__ __device__ int ACCS(const int i)
{
  return (i & ((LMEM_STACK_SIZE << SHIFT) - 1))*blockDim.x + threadIdx.x;
}


#define BTEST(x) (-(int)(x))

texture<float4, 1, cudaReadModeElementType> texNodeSize;
texture<float4, 1, cudaReadModeElementType> texNodeCenter;
texture<float4, 1, cudaReadModeElementType> texMultipole;
texture<float4, 1, cudaReadModeElementType> texBody;

template<class T>
 struct ADDOP {
  __device__ static inline T identity()           {return (T)(0);}
  __device__ static inline T apply(T a, T b)      {return (T)(a + b);};
  __device__ static inline T unapply(T a, T b)    {return (T)(a - b);};
  __device__ static inline T mask(bool flag, T b) {return (T)(-(int)(flag) & b);};
};

template<class OP, class T>
// __device__ T inclusive_scan_warp(volatile T *ptr, T mysum,  const unsigned int idx = threadIdx.x) {
__device__ __forceinline__ T inclusive_scan_warp(volatile T *ptr, T mysum,  const unsigned int idx ) {
  const unsigned int lane = idx & 31;

  if (lane >=  1) ptr[idx] = mysum = OP::apply(ptr[idx -  1], mysum);
  if (lane >=  2) ptr[idx] = mysum = OP::apply(ptr[idx -  2], mysum);
  if (lane >=  4) ptr[idx] = mysum = OP::apply(ptr[idx -  4], mysum);
  if (lane >=  8) ptr[idx] = mysum = OP::apply(ptr[idx -  8], mysum);
  if (lane >= 16) ptr[idx] = mysum = OP::apply(ptr[idx - 16], mysum);

  return ptr[idx];
}


__device__ __forceinline__ int inclusive_scan_warp(volatile int *ptr, int mysum, const unsigned int idx) {

  const unsigned int lane = idx & 31;

  if (lane >=  1) ptr[idx] = mysum = ptr[idx -  1]   + mysum;
  if (lane >=  2) ptr[idx] = mysum = ptr[idx -  2]   + mysum;
  if (lane >=  4) ptr[idx] = mysum = ptr[idx -  4]   + mysum;
  if (lane >=  8) ptr[idx] = mysum = ptr[idx -  8]   + mysum;
  if (lane >= 16) ptr[idx] = mysum = ptr[idx -  16]  + mysum;

  return ptr[idx];
}


template<class OP, class T>
__device__ __inline__ T inclusive_scan_block(volatile T *ptr, const T v0, const unsigned int idx) {
  const unsigned int lane   = idx & 31;
  const unsigned int warpid = idx >> 5;

  // step 0: Write the valume from the thread to the memory
  ptr[idx] = v0;
  T mysum = v0;
  __syncthreads();

  // step 1: Intra-warp scan in each warp
//   T val = inclusive_scan_warp<OP, T>(ptr, mysum, idx);
  T val = inclusive_scan_warp(ptr, mysum, idx);
  __syncthreads();

  // step 2: Collect per-warp particle results
  if (lane == 31) ptr[warpid] =  ptr[idx];
  __syncthreads();

  mysum =  ptr[idx];

  // step 3: Use 1st warp to scan per-warp results
  if (warpid == 0) inclusive_scan_warp<OP, T>(ptr,mysum, idx);
  __syncthreads();

  // step 4: Accumulate results from Steps 1 and 3;
  if (warpid > 0) val = OP::apply(ptr[warpid - 1], val);
  __syncthreads();

  // Step 5: Write and return the final result
  ptr[idx] = val;
  __syncthreads();

  return val; //ptr[blockDim.x - 1];
}



template<class OP, class T>
// __device__ T inclusive_scan_block(volatile T *ptr, const unsigned int idx = threadIdx.x) {
__device__ T inclusive_scan_block(volatile T *ptr, const unsigned int idx) {
  const unsigned int lane   = idx & 31;
  const unsigned int warpid = idx >> 5;

   T mysum = ptr[idx];
   __syncthreads();

  // step 1: Intra-warp scan in each warp
  T val = inclusive_scan_warp<OP, T>(ptr, mysum, idx);
  __syncthreads();

  // step 2: Collect per-warp particle results
  if (lane == 31) ptr[warpid] = ptr[idx];
  __syncthreads();

  mysum = ptr[idx];

  // step 3: Use 1st warp to scan per-warp results
  if (warpid == 0) inclusive_scan_warp<OP, T>(ptr,mysum, idx);
  __syncthreads();

  // step 4: Accumulate results from Steps 1 and 3;
  if (warpid > 0) val = OP::apply(ptr[warpid - 1], val);
  __syncthreads();

  // Step 5: Write and return the final result
  ptr[idx] = val;
  __syncthreads();

  return val; //ptr[blockDim.x - 1];
}


template<class OP, class T>
// __device__ T inclusive_scan_array(volatile T *ptr_global, const int N, const unsigned int idx = threadIdx.x) {
__device__ T inclusive_scan_array(volatile T *ptr_global, const int N, const unsigned int idx) {


  T y = OP::identity();
  volatile T *ptr = ptr_global;

  for (int p = 0; p < N; p += blockDim.x) {
    ptr = &ptr_global[p];
    inclusive_scan_block<OP, T>(ptr, idx);
    ptr[idx] = OP::apply(ptr[idx], y);
    __syncthreads();

    y = ptr[blockDim.x - 1];
    __syncthreads();
  }

  return y;

}

/*********** Forces *************/

__device__ float4 add_acc(
        float4 acc,  const float4 pos,
			  const float massj, const float3 posj,
			  const float eps2)
{
#if 1  /* to test performance of a tree-walk */
  const float3 dr = make_float3(posj.x - pos.x, posj.y - pos.y, posj.z - pos.z);

  const float r2     = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z + eps2;
  const float rinv   = rsqrtf(r2);
  const float rinv2  = rinv*rinv;
  const float mrinv  = massj * rinv;
  const float mrinv3 = mrinv * rinv2;

  acc.w -= mrinv;
  acc.x += mrinv3 * dr.x;
  acc.y += mrinv3 * dr.y;
  acc.z += mrinv3 * dr.z;
#endif

  return acc;
}


//Improved Barnes Hut criterium
__device__ bool split_node_grav_impbh(
    const float4 nodeCOM, 
    const float4 groupCenter, 
    const float4 groupSize)
{
  //Compute the distance between the group and the cell
  float3 dr = make_float3(
      fabsf(groupCenter.x - nodeCOM.x) - (groupSize.x),
      fabsf(groupCenter.y - nodeCOM.y) - (groupSize.y),
      fabsf(groupCenter.z - nodeCOM.z) - (groupSize.z)
      );

  dr.x += fabsf(dr.x); dr.x *= 0.5f;
  dr.y += fabsf(dr.y); dr.y *= 0.5f;
  dr.z += fabsf(dr.z); dr.z *= 0.5f;

  //Distance squared, no need to do sqrt since opening criteria has been squared
  const float ds2    = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

  return (ds2 <= fabsf(nodeCOM.w));
}



#define TEXTURES
#define OLDPREFIX
#if 0
#define _ORIG_SHMEM_
#endif


template<int DIM2, int SHIFT>
__device__ float4 approximate_gravity(int DIM2x, int DIM2y,
    int tid, int tx, int ty,
    int body_i, float4 pos_i,
    real4 group_pos,
    float eps2,
    uint2 node_begend,
    real4 *multipole_data,
    real4 *body_pos,
    int *shmem,
    int *lmem,
    int &ngb,
    int &apprCount, int &direCount,
    volatile float4 *boxSizeInfo,
    float4 groupSize,
    volatile float4 *boxCenterInfo,
    float group_eps,
    real4 *body_vel) {

  float4 acc_i = {0.0f, 0.0f, 0.0f, 0.0f};


  /*********** set necessary thread constants **********/

  const int DIMx = 1  << DIM2x;
  const int DIMy = 1  << DIM2y;
  const int DIM  = 1  << DIM2;
  const int offs = ty << DIM2x;

  /*********** shared memory distribution **********/

  //  begin,    end,   size
  // -----------------------
#ifdef _ORIG_SHMEM_
  
  int *approxS = (int*)&shmem  [     0];            //  0*DIM,  2*DIM,  2*DIM
  int *directS = (int*)&approxS[ 2*DIM];            //  2*DIM,  3*DIM,  1*DIM
  int *nodesS = (int*)&directS [   DIM];            //  3*DIM, 12*DIM,  9*DIM
  int *prefix = (int*)&nodesS  [9 *DIM];            // 12*DIM, 14*DIM,  2*DIM
  int *sh_body = &approxS[DIM];
  
  int *prefix0 = &prefix[  0];
  int *prefix1 = &prefix[DIM];
  
  const int NJMAX = DIM*2;
  int    *body_list = (int*   )&nodesS   [  DIM]; //  4*DIM,  6*DIM,  2*DIM
  float  *sh_mass   = (float* )&body_list[NJMAX]; //  6*DIM,  7*DIM,  1*DIM
  float3 *sh_pos    = (float3*)&sh_mass  [  DIM]; //  7*DIM, 10*DIM   3*DIM
  
  int *approxM = approxS;
  int *directM = directS;
  int * nodesM =  nodesS;

#else   /* !_ORIG_SHMEM_ */

  const int stack_sz = (LMEM_STACK_SIZE << SHIFT) << DIM2;  /* stack allocated per thread-block */
  int *approxL = lmem + stack_sz; 

  int *directS = shmem;                              //  0*DIM,  1*DIM,  1*DIM
  int *nodesS  = directS + DIM;                      //  1*DIM, 10*DIM,  9*DIM
  int *prefix  = nodesS  + DIM*9;                    // 10*DIM, 12*DIM,  2*DIM
  
  int *prefix0 = &prefix[  0];
  int *prefix1 = &prefix[DIM];
  
  const int NJMAX = DIM*3;
  int    *body_list = (int*   )&nodesS   [  DIM]; //  2*DIM,   5*DIM,  2*DIM
  float  *sh_mass   = (float* )&body_list[NJMAX]; //  5*DIM,   6*DIM,  1*DIM
  float3 *sh_pos    = (float3*)&sh_mass  [  DIM]; //  6*DIM,   9*DIM   3*DIM
  int    *sh_body   = nodesS + DIM*8;             //  9*DIM,  10*DIM,  1*DIM
  
  int *approxM = approxL;
  int *directM = directS;
  int * nodesM =  nodesS;

#endif /* _ORIG_SHMEM_ */


  float  *node_mon0 = sh_mass;
  float3 *node_mon1 = sh_pos; 
  
  float  *sh_pot = sh_mass;
  float3 *sh_acc = sh_pos;

  /*********** stack **********/

  int *nstack = lmem;

  /*********** begin tree-walk **********/

  int n_approx = 0;
  int n_direct = 0;


  for (int root_node = node_begend.x; root_node < node_begend.y; root_node += DIM) 
  {
    int n_nodes0 = min(node_begend.y - root_node, DIM);
    int n_stack0 = 0;
    int n_stack_pre = 0;

    { nstack[ACCS<SHIFT>(n_stack0)] = root_node + tid;   n_stack0++; }

    /*********** walk each level **********/
    while (n_nodes0 > 0) {


      int n_nodes1 = 0;
      int n_offset = 0;

      int n_stack1 = n_stack0;
      int c_stack0 = n_stack_pre;

      /*********** walk a level **********/
      while(c_stack0 < n_stack0) 
      {

        /***
         **** --> fetch the list of nodes rom LMEM
         ***/
        bool use_node = tid <  n_nodes0;
#if 0
        { prefix[tid] = nstack[ACCS<SHIFT>(c_stack0)];   c_stack0++; }
        __syncthreads();
        int node  = prefix[min(tid, n_nodes0 - 1)];
#else  /* eg: seems to work, but I do not remember if that will *always* work */
        int node;
        { node  = nstack[ACCS<SHIFT>(c_stack0)];   c_stack0++; }
#endif

#if 0
        if(n_nodes0 > 0){       //Work around pre 4.1 compiler bug
          n_nodes0 -= DIM;
        }
#else
        n_nodes0 -= DIM;
#endif

        /***
         **** --> process each of the nodes in the list in parallel
         ***/

#ifndef TEXTURES
        float4 nodeSize = boxSizeInfo[node];                   //Fetch the size of the box. Size.w = child info
        float4 node_pos = boxCenterInfo[node];                 //Fetch the center of the box. center.w = opening info
#else
        float4 nodeSize =  tex1Dfetch(texNodeSize, node);
        float4 node_pos =  tex1Dfetch(texNodeCenter, node);
#endif

        int node_data = __float_as_int(nodeSize.w);

        //Check if a cell has to be opened
#ifndef TEXTURES
        float4 nodeCOM = multipole_data[node*3];
#else
        float4 nodeCOM = tex1Dfetch(texMultipole,node*3);
#endif

        nodeCOM.w      = node_pos.w;
        bool   split   = split_node_grav_impbh(nodeCOM, group_pos, groupSize);


        bool leaf       = node_pos.w <= 0;  //Small AND equal incase of a 1 particle cell       //Check if it is a leaf
        //         split = true;

        uint mask    = BTEST((split && !leaf) && use_node);               // mask = #FFFFFFFF if use_node+split+not_a_leaf==true, otherwise zero
        int child    =    node_data & 0x0FFFFFFF;                         //Index to the first child of the node
        int nchild   = (((node_data & 0xF0000000) >> 28)) & mask;         //The number of children this node has

        /***
         **** --> calculate prefix
         ***/


#ifdef OLDPREFIX
        int n_total = calc_prefix<DIM2>(prefix, tid,  nchild);
        prefix[tid] += n_offset - nchild;
        __syncthreads();
#else
        inclusive_scan_block<ADDOP<int>, int>(prefix, nchild, tid);        // inclusive scan to compute memory offset of each child
        int n_total = prefix[blockDim.x - 1];                              // fetch total number of children, i.e. offset of the last child -1
        __syncthreads();                                                   // thread barrier to make sure that warps completed their jobs
        prefix[tid] += n_offset - nchild;                                  // convert inclusive into exclusive scan for referencing purpose
        __syncthreads();                                                   // thread barrier
#endif

        for (int i = n_offset; i < n_offset + n_total; i += DIM)         //nullify part of the array that will be filled with children
          nodesM[tid + i] = 0;                                          //but do not touch those parts which has already been filled
        __syncthreads();                                                 //Thread barrier to make sure all warps finished writing data

        bool flag = (split && !leaf) && use_node;                        //Flag = use_node + split + not_a_leaf;Use only non_leaf nodes that are to be split
#if 1
        if (flag) nodesM[prefix[tid]] = child;                            //Thread with the node that is about to be split
        __syncthreads();                                                 //writes the first child in the array of nodes

        /*** in the following 8 lines, we calculate indexes of all the children that have to be walked from the index of the first child***/
        if (flag && nodesM[prefix[tid] + 1] == 0) nodesM[prefix[tid] + 1] = child + 1; __syncthreads();
        if (flag && nodesM[prefix[tid] + 2] == 0) nodesM[prefix[tid] + 2] = child + 2; __syncthreads();
        if (flag && nodesM[prefix[tid] + 3] == 0) nodesM[prefix[tid] + 3] = child + 3; __syncthreads();
        if (flag && nodesM[prefix[tid] + 4] == 0) nodesM[prefix[tid] + 4] = child + 4; __syncthreads();
        if (flag && nodesM[prefix[tid] + 5] == 0) nodesM[prefix[tid] + 5] = child + 5; __syncthreads();
        if (flag && nodesM[prefix[tid] + 6] == 0) nodesM[prefix[tid] + 6] = child + 6; __syncthreads();
        if (flag && nodesM[prefix[tid] + 7] == 0) nodesM[prefix[tid] + 7] = child + 7; __syncthreads();
#else
#if 1
        if (flag) nodesM[prefix[tid]] = child;                            //Thread with the node that is about to be split
                                                                          //writes the first child in the array of nodes
#else
        const int maskT = flag ? 0xFFFFFFFF : 0x0;
        const int maskF = ~maskT;
        const int addr = (prefix[tid] & maskT) + (-1 & maskF);;
        nodesM[addr] = (maskF & nodesM[addr]) + (maskT & child);
#endif
        
        if (flag && nodesM[prefix[tid] + 1] == 0) nodesM[prefix[tid] + 1] = child + 1; 
        if (flag && nodesM[prefix[tid] + 2] == 0) nodesM[prefix[tid] + 2] = child + 2;
        if (flag && nodesM[prefix[tid] + 3] == 0) nodesM[prefix[tid] + 3] = child + 3;
        if (flag && nodesM[prefix[tid] + 4] == 0) nodesM[prefix[tid] + 4] = child + 4;
        if (flag && nodesM[prefix[tid] + 5] == 0) nodesM[prefix[tid] + 5] = child + 5;
        if (flag && nodesM[prefix[tid] + 6] == 0) nodesM[prefix[tid] + 6] = child + 6;
        if (flag && nodesM[prefix[tid] + 7] == 0) nodesM[prefix[tid] + 7] = child + 7;
        __syncthreads();
#endif

        n_offset += n_total;    //Increase the offset in the array by the number of newly added nodes


        /***
         **** --> save list of nodes to LMEM
         ***/

        /*** if half of shared memory or more is filled with the the nodes, dump these into slowmem stack ***/
        while(n_offset >= DIM) 
        {
          n_offset -= DIM;
          const int offs1 = ACCS<SHIFT>(n_stack1);
          nstack[offs1] = nodesM[n_offset + tid];   n_stack1++;
          n_nodes1 += DIM;

          if((n_stack1 - c_stack0) >= (LMEM_STACK_SIZE << SHIFT))
          {
            //We overwrote our current stack
            apprCount = -1; 
            return acc_i;	 
          }
        }

        __syncthreads();



        /******************************/
        /******************************/
        /*****     EVALUATION     *****/
        /******************************/
        /******************************/
#if 1
        /***********************************/
        /******       APPROX          ******/
        /***********************************/

#ifdef OLDPREFIX
        n_total = calc_prefix<DIM2>(prefix, tid,  1 - (split || !use_node));
#else
        inclusive_scan_block<ADDOP<int>, int>(prefix, 1 - (split || !use_node), tid);
        n_total = prefix[blockDim.x - 1];
#endif


        // 	n_total = calc_prefix<DIM2>(prefix, tid,  !split && use_node);         // for some unkown reason this does not work right on the GPU
        if (!split && use_node) approxM[n_approx + prefix[tid] - 1] = node;
        __syncthreads();
        n_approx += n_total;

        while (n_approx >= DIM) 
        {
          n_approx -= DIM;
          const int address      = (approxM[n_approx + tid] << 1) + approxM[n_approx + tid];
#ifndef TEXTURES
          const float4 monopole  = multipole_data[address    ];
#if 0
          float4 octopole0 = multipole_data[address + 1];
          float4 octopole1 = multipole_data[address + 2];
#endif
#else
          const float4 monopole  = tex1Dfetch(texMultipole, address);
#if 0
          float4 octopole0 = tex1Dfetch(texMultipole, address + 1);
          float4 octopole1 = tex1Dfetch(texMultipole, address + 2);
#endif
#endif

          node_mon0[tid] = monopole.w;
          node_mon1[tid] = make_float3(monopole.x,  monopole.y,  monopole.z);
          __syncthreads();

#if 0
          const float f_dm   = 0.0f;
          const float f_star = 1.0f
            const float darkMatterMass = f_dm   * octopole1.w;
          /* eg: we need to be careful with the line below to avoid truncation error due to 
             subtraction of two large numbers, monopole.w and darkMatterMass both could be
             very large.
             Instead, we can use octopole1.w to be stellar mass, and DM mass to be 
             monopole.w, then we add the two together to get total mass, but this will
             require more changes to the kernel */
          const float    stellarMass = f_star * (monopole.w - darkMatterMass);
          const float hinv = 1.0f/hi;   /* eg: this can be precomputing to avoid division */
          density += interact(
              make_float3(pos_i.x, pos_i.y, pos_i.z), h, hinv,
              make_float3(monopole.x, monople.y, monopole.z), darkMatterMass + stellarMass);
          /* eg: the interact function still calls sqrtf(f), which invloves 1 div and 1 rsqrtf,
             so ideally we would like to take advantage of rsqrtf in add_acc, and then we only
             do 1 div */
#endif


#if 1
#pragma unroll 16
          for (int i = 0; i < DIMx; i++)
            acc_i = add_acc(acc_i, pos_i, node_mon0[offs + i], node_mon1[offs+i], eps2);
          apprCount += DIMx;
          __syncthreads();
#endif
        }
        __syncthreads();
#endif

#if 1
        /***********************************/
        /******       DIRECT          ******/
        /***********************************/


        flag         = split && leaf && use_node;                                //flag = split + leaf + use_node
        int  jbody   = node_data & BODYMASK;                                     //the first body in the leaf
        int  nbody   = (((node_data & INVBMASK) >> LEAFBIT)+1) & BTEST(flag);    //number of bodies in the leaf masked with the flag

        body_list[tid] = directM[tid];                                            //copy list of bodies from previous pass to body_list
        sh_body  [tid] = jbody;                                                  //store the leafs first body id into shared memory

        // step 1
#ifdef OLDPREFIX
        calc_prefix<DIM2>(prefix0, tid, flag);
#else
        inclusive_scan_block<ADDOP<int>, int>(prefix0, (int)flag, tid);       // inclusive scan on flags to construct array
#endif

        if (flag) prefix1[prefix0[tid] - 1] = tid;                             //with tidś whose leaves have to be opened
        __syncthreads();                                                      //thread barrier, make sure all warps completed the job

        // step 2
#ifdef OLDPREFIX
        int n_bodies  = calc_prefix<DIM2>(prefix0, tid, nbody);
#else
        inclusive_scan_block<ADDOP<int>, int>(prefix0, nbody, tid);        // inclusive scan to compute memory offset for each body
        int n_bodies = prefix0[blockDim.x - 1];                            //Total number of bides extract from the leaves
        __syncthreads();                                                   // thread barrier to make sure that warps completed their jobs
#endif

        directM[tid]  = prefix0[tid];                                       //Store a copy of inclusive scan in direct
        prefix0[tid] -= nbody;                                              //convert inclusive int oexclusive scan
        prefix0[tid] += 1;                                                  //add unity, since later prefix0[tid] == 0 used to check barrier

        int nl_pre = 0;                                                     //Number of leaves that have already been processed

        while (n_bodies > 0) 
        {
          int nb    = min(n_bodies, NJMAX - n_direct);                    //Make sure number of bides to be extracted does not exceed
          //the amount of allocated shared memory

          // step 0                                                      //nullify part of the body_list that will be filled with bodies
          for (int i = n_direct; i < n_direct + nb; i += DIM){           //from the leaves that are being processed
            body_list[i + tid] = 0;
          }
          __syncthreads();

          //step 1:
          if (flag && (directM[tid] <= nb) && (prefix0[tid] > 0))        //make sure that the thread indeed carries a leaf
            body_list[n_direct + prefix0[tid] - 1] = 1;                 //whose bodies will be extracted
          __syncthreads();

          //step 2:
#ifdef OLDPREFIX
          int nl = calc_prefix<DIM2>(nb, &body_list[n_direct], tid);
#else
          int nl = inclusive_scan_array<ADDOP<int>, int>              // inclusive scan to compute number of leaves to process
            (&body_list[n_direct], nb, tid);            // to make sure that there is enough shared memory for bodies
#endif
          nb = directM[prefix1[nl_pre + nl - 1]];                        // number of bodies stored in these leaves

          // step 3:
          for (int i = n_direct; i < n_direct + nb; i += DIM) {          //segmented fill of the body_list
            int j = prefix1[nl_pre + body_list[i + tid] - 1];            // compute the first body in shared j-body array
            body_list[i + tid] = (i + tid - n_direct) -                 //add to the index of the first j-body in a child
              (prefix0[j] - 1) + sh_body[j];         //the index of the first child in body_list array
          }
          __syncthreads();


          /**************************************************
           *  example of what is accomplished in steps 0-4   *
           *       ---------------------------               *
           * step 0: body_list = 000000000000000000000       *
           * step 1: body_list = 100010001000000100100       *
           * step 2: body_list = 111122223333333444555       *
           * step 3: body_list = 012301230123456012012       *
           *         assuming that sh_body[j] = 0            *
           ***************************************************/

          n_bodies     -= nb;                                   //subtract from n_bodies number of bodies that have been extracted
          nl_pre       += nl;                                   //increase the number of leaves that where processed
          directM[tid] -= nb;                                   //subtract the number of extracted bodies in this pass
          prefix0[tid] = max(prefix0[tid] - nb, 0);             //same here, but do not let the number be negative (GT200 bug!?)
          n_direct     += nb;                                  //increase the number of bodies to be procssed

          while(n_direct >= DIM) 
          {
            n_direct -= DIM;


            const float4 posj  = body_pos[body_list[n_direct + tid]];
#if 0
            const float4 posj  = tex1Dfetch(texBody, body_list[n_direct + tid]);
#endif
            sh_mass[tid] = posj.w;
            sh_pos [tid] = make_float3(posj.x, posj.y, posj.z);

            __syncthreads();
#if 1
#pragma unroll 16
            for (int j = 0; j < DIMx; j++)
              acc_i = add_acc(acc_i, pos_i, sh_mass[offs + j], sh_pos[offs + j], eps2);
            direCount += DIMx;
            __syncthreads();
#endif
          }

        }
        directM[tid] = body_list[tid];
        __syncthreads();
#endif
      } //end lvl


      n_nodes1 += n_offset;
      if (n_offset > 0)
      { 
        nstack[ACCS<SHIFT>(n_stack1)] = nodesM[tid];   n_stack1++; 
        if((n_stack1 - c_stack0) >= (LMEM_STACK_SIZE << SHIFT))
        {
          //We overwrote our current stack
          apprCount = -1; 
          return acc_i;	 
        }
      }
      __syncthreads();


      /***
       **** --> copy nodes1 to nodes0: done by reassigning the pointers
       ***/
      n_nodes0    = n_nodes1;

      n_stack_pre = n_stack0;
      n_stack0    = n_stack1;

    }//end while   levels
  }//end for


  if(n_approx > 0)
  {
    if (tid < n_approx) 
    {
      const int address = (approxM[tid] << 1) + approxM[tid];
#ifndef TEXTURES
      float4 monopole  = multipole_data[address    ];
      float4 octopole0 = multipole_data[address + 1];
      float4 octopole1 = multipole_data[address + 2];
#else
      float4 monopole  = tex1Dfetch(texMultipole, address);
      float4 octopole0 = tex1Dfetch(texMultipole, address + 1);
      float4 octopole1 = tex1Dfetch(texMultipole, address + 2);
#endif

      node_mon0[tid] = monopole.w;
      node_mon1[tid] = make_float3(monopole.x,  monopole.y,  monopole.z);

    } else {

      //Set non-active memory locations to zero
      node_mon0[tid] = 0.0f;
      node_mon1[tid] = make_float3(1.0e10f, 1.0e10f, 1.0e10f);

    }
    __syncthreads();
#pragma unroll
    for (int i = 0; i < DIMx; i++)
      acc_i = add_acc(acc_i, pos_i, node_mon0[offs + i], node_mon1[offs+i],eps2);
    apprCount += DIMx;

    __syncthreads();
  } //if n_approx > 0

  if(n_direct > 0)
  {
    if (tid < n_direct) 
    {
      const float4 posj = body_pos[directM[tid]];
#if 0
      const float4 posj  = tex1Dfetch(texBody, direct[tid]);
#endif
      sh_mass[tid] = posj.w;
      sh_pos [tid] = make_float3(posj.x, posj.y, posj.z);
    } else {
      sh_mass[tid] = 0.0f;
      sh_pos [tid] = make_float3(1.0e10f, 1.0e10f, 1.0e10f);
    }

    __syncthreads();
#pragma unroll
    for (int j = 0; j < DIMx; j++) 
      acc_i = add_acc(acc_i, pos_i, sh_mass[offs + j], sh_pos[offs + j], eps2);
    direCount += DIMx;
    __syncthreads();
  }

  /***
   **** --> reduce data between threads
   ***/
  sh_pot[tid] = acc_i.w;
  sh_acc[tid] = make_float3(acc_i.x, acc_i.y, acc_i.z);
  __syncthreads();

  if (ty == 0) 
#pragma unroll
    for (int i = 1; i < DIMy; i++) 
    {
      const int idx = (i << DIM2x) + tx;
      acc_i.w += sh_pot[idx];
      acc_i.x += sh_acc[idx].x;
      acc_i.y += sh_acc[idx].y;
      acc_i.z += sh_acc[idx].z;
    }
  __syncthreads();


  //Sum the interaction counters
  float  *sh_ds2 = (float*)&sh_acc[DIM];
  int    *sh_ngb = (int*  )&sh_ds2[DIM];
  sh_ds2[tid] = direCount;
  sh_ngb[tid] = apprCount;

  __syncthreads();


  if (ty == 0) {
#pragma unroll
    for (int i = 1; i < DIMy; i++){
      int idx = (i << DIM2x) + tx;
      direCount  += sh_ds2[idx];
      apprCount  += sh_ngb[idx];
    }
  }
  __syncthreads();

  return acc_i;
}


  extern "C" __global__ void
__launch_bounds__(NTHREAD)
  dev_approximate_gravity(const int n_active_groups,
      int    n_bodies,
      float eps2,
      uint2 node_begend,
      int    *active_groups,
      real4  *body_pos,
      real4  *multipole_data,
      float4 *acc_out,
      real4  *group_body_pos,
      int    *ngb_out,
      int    *active_inout,
      int2   *interactions,
      float4  *boxSizeInfo,
      float4  *groupSizeInfo,
      float4  *boxCenterInfo,
      float4  *groupCenterInfo,
      real4   *body_vel,
      int     *MEM_BUF) {
    //                                                    int     grpOffset){


    const int blockDim2 = NTHREAD2;
#ifdef _ORIG_SHMEM_
    __shared__ int shmem[15*(1 << blockDim2)];
#else
    __shared__ int shmem[12*(1 << blockDim2)];
#endif
    //    __shared__ int shmem[24*(1 << blockDim2)]; is possible on FERMI
    //    int             lmem[LMEM_STACK_SIZE];



    /*********** check if this block is linked to a leaf **********/

    int bid = gridDim.x * blockIdx.y + blockIdx.x;

    while(true)
    {

      if(threadIdx.x == 0)
      {
        bid         = atomicAdd(&active_inout[n_bodies], 1);
        shmem[0]    = bid;
      }
      __syncthreads();

      bid   = shmem[0];

      if (bid >= n_active_groups) return;


      int tid = threadIdx.y * blockDim.x + threadIdx.x;

      int grpOffset = 0;

      //   volatile int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x + threadIdx.x*LMEM_STACK_SIZE];
      //   int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x + threadIdx.x*LMEM_STACK_SIZE];
#ifdef _ORIG_SHMEM_
      int *lmem = &MEM_BUF[blockIdx.x* LMEM_STACK_SIZE*blockDim.x];
#else
      int *lmem = &MEM_BUF[blockIdx.x*(LMEM_STACK_SIZE*blockDim.x + LMEM_EXTRA_SIZE)];
#endif


      /*********** set necessary thread constants **********/
#ifdef DO_BLOCK_TIMESTEP
      real4 curGroupSize    = groupSizeInfo[active_groups[bid + grpOffset]];
#else
      real4 curGroupSize    = groupSizeInfo[bid + grpOffset];
#endif
      int   groupData       = __float_as_int(curGroupSize.w);
      uint body_i           =   groupData & CRITMASK;
      uint nb_i             = ((groupData & INVCMASK) >> CRITBIT) + 1;

#ifdef DO_BLOCK_TIMESTEP
      real4 group_pos       = groupCenterInfo[active_groups[bid + grpOffset]];
#else
      real4 group_pos       = groupCenterInfo[bid + grpOffset];
#endif
      //   if(tid == 0)
      //   printf("[%f %f %f %f ] \n [%f %f %f %f ] %d %d \n",
      //           curGroupSize.x, curGroupSize.y, curGroupSize.z, curGroupSize.w,
      //           group_pos.x, group_pos.y, group_pos.z, group_pos.w, body_i, nb_i);


      int DIM2x = 0;
      while (((nb_i - 1) >> DIM2x) > 0) DIM2x++;

      DIM2x     = max(DIM2x,4);
      int DIM2y = blockDim2 - DIM2x;

      int tx = tid & ((1 << DIM2x) - 1);
      int ty = tid >> DIM2x;

      body_i += tx%nb_i;

      //float4 pos_i = tex1Dfetch(bodies_pos_ref, body_i);   // texture read: 4 floats


//       float4 pos_i = body_pos[body_i];
      float4 pos_i = group_body_pos[body_i];



      int ngb_i;

      float4 acc_i = {0.0f, 0.0f, 0.0f, 0.0f};

#ifdef INDSOFT
      eps2 = body_vel[body_i].w;
      float group_eps = eps2;

      volatile float *reduc = (float*) &shmem[0];
      reduc[threadIdx.x] = eps2;

      //Find the maximum softening value for the particles in this group
      __syncthreads();
      // do reduction in shared mem
      if(blockDim.x >= 512) if (tid < 256) {reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 256]);} __syncthreads();
      if(blockDim.x >= 256) if (tid < 128) {reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 128]);} __syncthreads();
      if(blockDim.x >= 128) if (tid < 64)  {reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 64]);} __syncthreads();
      if(blockDim.x >= 64) if (tid < 32) { reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 32]);}
      if(blockDim.x >= 32) if (tid < 16) { reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 16]);}

      if(tid < 8)
      {
        reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 8]);
        reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 4]);
        reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 2]);
        reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 1]);
      }
      __syncthreads();

      group_eps = reduc[0];
#else
      float group_eps  = 0;
#endif

      int apprCount = 0;
      int direCount = 0;


      acc_i = approximate_gravity<blockDim2, 0>( DIM2x, DIM2y, tid, tx, ty,
          body_i, pos_i, group_pos,
          eps2, node_begend,
          multipole_data, body_pos,
          shmem, lmem, ngb_i, apprCount, direCount, boxSizeInfo, curGroupSize, boxCenterInfo,
          group_eps, body_vel);
      if(apprCount < 0)
      {

        //Try to get access to the big stack, only one block per time is allowed
        if(threadIdx.x == 0)
        {
          int res = atomicExch(&active_inout[n_bodies+1], 1); //If the old value (res) is 0 we can go otherwise sleep
          int waitCounter  = 0;
          while(res != 0)
          {
            //Sleep
            for(int i=0; i < (1024); i++)
            {
              waitCounter += 1;
            }
            //Test again
            shmem[0] = waitCounter;
            res = atomicExch(&active_inout[n_bodies+1], 1); 
          }
        }

        __syncthreads();

#ifdef _ORIG_SHMEM_
        lmem = &MEM_BUF[gridDim.x*LMEM_STACK_SIZE*blockDim.x];    //Use the extra large buffer
#else
        lmem = &MEM_BUF[gridDim.x*(LMEM_STACK_SIZE*blockDim.x + LMEM_EXTRA_SIZE)];    //Use the extra large buffer
#endif
        apprCount = direCount = 0;
        acc_i = approximate_gravity<blockDim2, 8>( DIM2x, DIM2y, tid, tx, ty,
            body_i, pos_i, group_pos,
            eps2, node_begend,
            multipole_data, body_pos,
            shmem, lmem, ngb_i, apprCount, direCount, boxSizeInfo, curGroupSize, boxCenterInfo,
            group_eps, body_vel);

#ifdef _ORIG_SHMEM_
        lmem = &MEM_BUF[blockIdx.x* LMEM_STACK_SIZE*blockDim.x]; //Back to normal location
#else
        lmem = &MEM_BUF[blockIdx.x*(LMEM_STACK_SIZE*blockDim.x + LMEM_EXTRA_SIZE)];
#endif

        if(threadIdx.x == 0)
        {
          atomicExch(&active_inout[n_bodies+1], 0); //Release the lock
        }
      }//end if apprCount < 0

      if (tid < nb_i) {
        acc_out     [body_i] = acc_i;
        ngb_out     [body_i] = -1;
        active_inout[body_i] = 1;
        interactions[body_i].x = apprCount;
        interactions[body_i].y = direCount ;
      }


    }     //end while
  }
