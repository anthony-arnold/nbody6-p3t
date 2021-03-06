#include "bonsai.h"
//Definitions
#ifndef _SUPPORT_KERNELS_
#define _SUPPORT_KERNELS_

#include "node_specs.h"

#include <stdio.h>



static __device__ int undilate3(uint2 key) {
  int x, value = 0;
  
  key.x = key.x & 0x09249249;
  key.y = key.y & 0x09249249;
  
  // undilate first 10 bits

  x = key.y & 0x3FFFF;
  x = ((x <<  4) + (x << 2) + x) & 0x0E070381;
  x = ((x << 12) + (x << 6) + x) & 0x0FF80001;
  x = ((x << 18) + x) & 0x0FFC0000;
  value = value | (x >> 18);
  
  x = (key.y >> 18) & 0x3FFFF;
  x = ((x <<  4) + (x << 2) + x) & 0x0E070381;
  x = ((x << 12) + (x << 6) + x) & 0x0FF80001;
  x = ((x << 18) + x) & 0x0FFC0000;
  value = value | (x >> 12);
  

  // undilate second 10 bits

  x = key.x & 0x3FFFF;
  x = ((x <<  4) + (x << 2) + x) & 0x0E070381;
  x = ((x << 12) + (x << 6) + x) & 0x0FF80001;
  x = ((x << 18) + x) & 0x0FFC0000;
  value = value | ((x >> 18) << 10);
  
  x = (key.x >> 18) & 0x3FFFF;
  x = ((x <<  4) + (x << 2) + x) & 0x0E070381;
  x = ((x << 12) + (x << 6) + x) & 0x0FF80001;
  x = ((x << 18) + x) & 0x0FFC0000;
  value = value | ((x >> 12) << 10);
  
  return value;
}


static __device__ uint2 dilate3(int value) {
  unsigned int x;
  uint2 key;
  
  // dilate first 10 bits

  x = value & 0x03FF;
  x = ((x << 16) + x) & 0xFF0000FF;
  x = ((x <<  8) + x) & 0x0F00F00F;
  x = ((x <<  4) + x) & 0xC30C30C3;
  x = ((x <<  2) + x) & 0x49249249;
  key.y = x;

  // dilate second 10 bits

  x = (value >> 10) & 0x03FF;
  x = ((x << 16) + x) & 0xFF0000FF;
  x = ((x <<  8) + x) & 0x0F00F00F;
  x = ((x <<  4) + x) & 0xC30C30C3;
  x = ((x <<  2) + x) & 0x49249249;
  key.x = x;

  return key;
} 

//#if 0
//Morton order
static __device__ uint2 get_key_morton(int4 crd) {
  uint2 key, key1;
  key  = dilate3(crd.x);

  key1 = dilate3(crd.y);
  key.x = key.x | (key1.x << 1);
  key.y = key.y | (key1.y << 1);

  key1 = dilate3(crd.z);
  key.x = key.x | (key1.x << 2);
  key.y = key.y | (key1.y << 2);

  return key;
}

//#else

static __device__ uint4 get_key(int4 crd)
{
  const int bits = 30;  //20 to make it same number as morton order
  int i,xi, yi, zi;
  int mask;
  int key;

  mask = crd.y;
  crd.y = crd.z;
  crd.z = mask;
    
  //0= 000, 1=001, 2=011, 3=010, 4=110, 5=111, 6=101, 7=100
  //000=0=0, 001=1=1, 011=3=2, 010=2=3, 110=6=4, 111=7=5, 101=5=6, 100=4=7
  const int C[8] = {0, 1, 7, 6, 3, 2, 4, 5};
    
  int temp;
    
  mask = 1 << (bits - 1);
  key  = 0;

  uint4 key_new;
    
  for(i = 0; i < bits; i++, mask >>= 1)
  {
    xi = (crd.x & mask) ? 1 : 0;
    yi = (crd.y & mask) ? 1 : 0;
    zi = (crd.z & mask) ? 1 : 0;        

    int index = (xi << 2) + (yi << 1) + zi;
      
    if(index == 0)
    {
      temp = crd.z; crd.z = crd.y; crd.y = temp;
    }
    else  if(index == 1 || index == 5)
    {
      temp = crd.x; crd.x = crd.y; crd.y = temp;
    }
    else  if(index == 4 || index == 6)
    {
      crd.x = (crd.x) ^ (-1);
      crd.z = (crd.z) ^ (-1);
    }
    else  if(index == 7 || index == 3)
    {
      temp = (crd.x) ^ (-1);         
      crd.x = (crd.y) ^ (-1);
      crd.y = temp;
    }
    else
    {
      temp = (crd.z) ^ (-1);         
      crd.z = (crd.y) ^ (-1);
      crd.y = temp;          
    }   

    key = (key << 3) + C[index];

    if(i == 19)
    {
      key_new.y = key;
      key = 0;
    }
    if(i == 9)
    {
      key_new.x = key;
      key = 0;
    }
  } //end for

   key_new.z = key;

  return key_new;
}

//#endif


static __device__ uint4 get_mask(int level) {
  int mask_levels = 3*max(MAXLEVELS - level, 0);
  uint4 mask = {0x3FFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,0xFFFFFFFF};

  if (mask_levels > 60)
  {
    mask.z = 0;
    mask.y = 0;
    mask.x = (mask.x >> (mask_levels - 60)) << (mask_levels - 60);
  }
  else if (mask_levels > 30) {
    mask.z = 0;
    mask.y = (mask.y >> (mask_levels - 30)) << (mask_levels - 30);
  } else {
    mask.z = (mask.z >> mask_levels) << mask_levels;
  }

  return mask;
}

static __device__ uint4 get_imask(uint4 mask) {
  return make_uint4(0x3FFFFFFF ^ mask.x, 0xFFFFFFFF ^ mask.y, 0xFFFFFFFF ^ mask.z, 0);
}



static __device__ int4 get_crd(uint2 key) {
  int4 crd;

  crd.x = undilate3(key);
  crd.y = undilate3(make_uint2(key.x >> 1, key.y >> 1));
  crd.z = undilate3(make_uint2(key.x >> 2, key.y >> 2));
  
  return crd;
}

static __device__ inline int cmp_uint2(uint2 a, uint2 b) {
  if      (a.x < b.x) return -1;
  else if (a.x > b.x) return +1;
  else {
    if       (a.y < b.y) return -1;
    else  if (a.y > b.y) return +1;
    return 0;
  }  
}

static __device__ int cmp_uint4(uint4 a, uint4 b) {
  if      (a.x < b.x) return -1;
  else if (a.x > b.x) return +1;
  else {
    if       (a.y < b.y) return -1;
    else  if (a.y > b.y) return +1;
    else {
      if       (a.z < b.z) return -1;
      else  if (a.z > b.z) return +1;
      return 0;
    } //end z    
  }  //end y
} //end x, function


//Binary search of the key within certain bounds (cij.x, cij.y)
static __device__ int find_key(uint4 key, uint2 cij, uint4 *keys) {
  int l = cij.x;
  int r = cij.y - 1;
  while (r - l > 1) {
    int m = (r + l) >> 1;
    int cmp = cmp_uint4(keys[m], key);
    if (cmp == -1) {
      l = m;
    } else { 
      r = m;
    }
  }
  if (cmp_uint4(keys[l], key) >= 0) return l;

  return r;
}



static __device__ float2 ds_accumulate(float2 a, float b){
  float tmp = a.x + b;
  float del = (tmp - a.x) - b;
  a.x = tmp;
  a.y -= del;
  return a;
}
static __device__ float2 ds_regularise(float2 a){
  float tmp = a.x + a.y;
  a.y -= (tmp - a.x);
  a.x = tmp;
  return a;
}


static __device__ __forceinline__ void sh_MinMax(int i, int j, float3 *r_min, float3 *r_max,
                          volatile float3 *sh_rmin, volatile  float3 *sh_rmax)
{
  sh_rmin[i].x  = (*r_min).x = fminf((*r_min).x, sh_rmin[j].x);
  sh_rmin[i].y  = (*r_min).y = fminf((*r_min).y, sh_rmin[j].y);
  sh_rmin[i].z  = (*r_min).z = fminf((*r_min).z, sh_rmin[j].z);
  sh_rmax[i].x  = (*r_max).x = fmaxf((*r_max).x, sh_rmax[j].x);
  sh_rmax[i].y  = (*r_max).y = fmaxf((*r_max).y, sh_rmax[j].y);
  sh_rmax[i].z  = (*r_max).z = fmaxf((*r_max).z, sh_rmax[j].z);
}

static __device__ void MinMaxPos(float4 pos, float4 &rmax, float4 &rmin)
{
      rmin.x  = fminf(pos.x, rmin.x);
      rmin.y  = fminf(pos.y, rmin.y);
      rmin.z  = fminf(pos.z, rmin.z);
      rmax.x  = fmaxf(pos.x, rmax.x); 
      rmax.y  = fmaxf(pos.y, rmax.y); 
      rmax.z  = fmaxf(pos.z, rmax.z); 
}


static __device__ real4 get_pos(uint2 key, float size, float4 corner) {
  real4 pos;
  pos.w = size;
  
  int4 crd = get_crd(key);
  float domain_fac = corner.w;
  pos.x = crd.x*domain_fac + corner.x;
  pos.y = crd.y*domain_fac + corner.y;
  pos.z = crd.z*domain_fac + corner.z;

  return pos;
}


#endif
