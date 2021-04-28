/*
 *
 * Functions to read / write to various file formats
 *
 *
 *
 */

#pragma once



#include "IDType.h"

#ifdef USE_MPI
    #include "BonsaiIO.h"
#endif

/************* data exchange containers for async IO ***************/

struct IOSharedData_t
{
  volatile bool writingFinished;
  volatile float t_current;
  volatile int   nBodies;
  unsigned long long * volatile  IDs;
  real4 * volatile Pos, * volatile Vel;
  IOSharedData_t() : writingFinished(true), nBodies(0), IDs(NULL), Pos(NULL), Vel(NULL) {}
  void malloc(const int n) volatile
  {
    assert(nBodies == 0);
    nBodies = n;
    IDs = (unsigned long long*volatile)::malloc(n*sizeof(unsigned long long));
    Pos = (real4*volatile)::malloc(n*sizeof(real4));
    Vel = (real4*volatile)::malloc(n*sizeof(real4));
  }
  void free() volatile
  {
    assert(nBodies > 0);
    nBodies = 0;
    ::free(IDs);
    ::free(Pos);
    ::free(Vel);
  }
  ~IOSharedData_t()
  {
    if (nBodies > 0)
      free();
  }
};


extern volatile IOSharedData_t ioSharedData;


template<typename T>
static void lHandShake(SharedMemoryBase<T> &header)
{
  header.acquireLock();
  header[0].handshake = false;
  header.releaseLock();

  while (!header[0].handshake)
    usleep(10000);

  header.acquireLock();
  header[0].handshake = false;
  header.releaseLock();
}

static IDType lGetIDType(const long long id)
{
  IDType ID;
  ID.setID(id);
  ID.setType(3);     /* Everything is Dust until told otherwise */
  if(id >= DISKID  && id < BULGEID)
  {
    ID.setType(2);  /* Disk */
    ID.setID(id - DISKID);
  }
  else if(id >= BULGEID && id < DARKMATTERID)
  {
    ID.setType(1);  /* Bulge */
    ID.setID(id - BULGEID);
  }
  else if (id >= DARKMATTERID)
  {
    ID.setType(0);  /* DM */
    ID.setID(id - DARKMATTERID);
  }
  return ID;
};
