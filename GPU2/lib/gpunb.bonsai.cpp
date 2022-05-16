/******************************************************************************
 * Author: Anthony D. Arnold
 * Institution: School of Maths and Physics, University of Queensland
 * Contact: anthony.arnold@uqconnect.edu.au
 * File: gpunb.bonsai.cpp
 * Description: The glue between NBODY6 and the Bonsai code.
 *****************************************************************************/

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>
#include <omp.h>
#include "octree.h"

#include <thrust/system/cuda/execution_policy.h>
#include <thrust/remove.h>

#define PROFILE 1

extern "C" double double_env(const char* name);

// Bonsai static data
bool ENABLE_RUNTIME_LOG;
bool PREPEND_RANK;
int  PREPEND_RANK_PROCID;
int  PREPEND_RANK_NPROCS;

volatile IOSharedData_t ioSharedData;

long long my_dev::base_mem::currentMemUsage;
long long my_dev::base_mem::maxMemUsage;

extern cudaEvent_t startLocalGrav;
extern cudaEvent_t startRemoteGrav;
extern cudaEvent_t endLocalGrav;
extern cudaEvent_t endRemoteGrav;

#if PROFILE
#include <chrono>
namespace {
using clock_type = std::chrono::high_resolution_clock;
using tick_type = clock_type::duration;

tick_type total_time;
int num_ops;

auto now() -> decltype(clock_type::now()) {
    return clock_type::now();
}
}
#endif

namespace {

constexpr MPI_Comm mpi_comm = 0;

class p3tree : public octree {
public:
    p3tree(double theta, double eps)
        : octree(mpi_comm, &context_, nullptr, 0, theta, eps),
          exec_stream_(std::make_shared<my_dev::dev_stream>(0)),
          grav_stream_(std::make_shared<my_dev::dev_stream>(0)),
          copy_stream_(std::make_shared<my_dev::dev_stream>(0))
    {
        context_.create(log_, false);
        context_.createQueue(0);
        context_.setLogPreamble("P3T");

        execStream = exec_stream_.get();
        gravStream = grav_stream_.get();
        copyStream = copy_stream_.get();
    }

    cudaStream_t copy_stream() const {
        return copyStream->s();
    }

    cudaStream_t exec_stream() const {
        return execStream->s();
    }

    void sync_copy() {
        copyStream->sync();
    }
    void sync_exec() {
        execStream->sync();
    }
    void sync_grav() {
        gravStream->sync();
    }
private:
    std::stringstream log_;
    my_dev::context context_;

    std::shared_ptr<my_dev::dev_stream> exec_stream_;
    std::shared_ptr<my_dev::dev_stream> grav_stream_;
    std::shared_ptr<my_dev::dev_stream> copy_stream_;

};

void copy_to_bonsai(
    int ni,
    double mi[],
    double xi[][3],
    double vi[][3],
    double rcut[],
    p3tree& tree,
    bool initial = false)
{
    tree.sync_exec();
#pragma omp parallel for
    for (int i = 0; i < ni; i++) {
        real4 pos;

        pos.x = static_cast<float>(xi[i][0]);
        pos.y = static_cast<float>(xi[i][1]);
        pos.z = static_cast<float>(xi[i][2]);
        pos.w = static_cast<float>(mi[i]);

        tree.localTree.bodies_pos[i] = pos;
        tree.localTree.bodies_Ppos[i] = pos;

        tree.localTree.bodies_h[i] = rcut[i];

        real4 vel;
        vel.x = static_cast<float>(vi[i][0]);
        vel.y = static_cast<float>(vi[i][1]);
        vel.z = static_cast<float>(vi[i][2]);
        vel.w = 0.0f;

        tree.localTree.bodies_vel[i] = vel;
        tree.localTree.bodies_Pvel[i] = vel;

        tree.localTree.bodies_ids[i] = i;

        tree.localTree.bodies_time[i] = make_float2(0.0f, 0.0f);
    }

    tree.localTree.bodies_pos. h2d(false, tree.copy_stream());
    tree.localTree.bodies_Ppos.h2d(false, tree.exec_stream());
    tree.localTree.bodies_h.h2d(false, tree.copy_stream());
    tree.localTree.bodies_Pvel.h2d(false, tree.exec_stream());
    tree.localTree.bodies_vel. h2d();
    tree.localTree.bodies_ids. h2d();
    tree.localTree.bodies_time.h2d();

    tree.sync_copy();
    if (initial) {
        tree.sort_bodies(tree.localTree, true, true);
    }
    tree.sort_bodies(tree.localTree, true);
    tree.build(tree.localTree);
    tree.allocateTreePropMemory(tree.localTree);
    tree.compute_properties(tree.localTree);
}

p3tree* SYSTEM(NULL);
int NMAX;
int NBMAX;
double THETA;
double EPS;
bool FIRST_SORT;
std::vector<int> ID2IDX;

// init
bool devinit = false;
int numCPU, numGPU;
constexpr int MAX_CPU = 38;
constexpr int MAX_GPU = 4;
int devid[MAX_GPU];



void approx_grav(
        int *ni,
        double mi[],
        double xi[][3],
        double vi[][3],
        double acc[][3],
        double jrk[][3],
        double rcut[]) {

   SYSTEM->localTree.setN(*ni);
   SYSTEM->set_t_current(0.0f);
   copy_to_bonsai(*ni, mi, xi, vi, rcut, *SYSTEM, FIRST_SORT);
   FIRST_SORT = false;
   SYSTEM->approximate_gravity(SYSTEM->localTree);

   ID2IDX.resize(*ni);

   SYSTEM->sync_grav();
   SYSTEM->localTree.neighbours.d2h(false, SYSTEM->copy_stream());
   SYSTEM->localTree.bodies_dens.d2h(false, SYSTEM->exec_stream());


   SYSTEM->localTree.bodies_ids.d2h();
#pragma omp parallel for
   for (int i = 0; i < *ni; i++) {
      ID2IDX[SYSTEM->localTree.bodies_ids[i]] = i;
   }

   SYSTEM->localTree.bodies_acc1.d2h();
#pragma omp parallel for schedule(static)
   for (int i = 0; i < *ni; i++) {
      int idx = ID2IDX[i];
      const real4& dev_acc = SYSTEM->localTree.bodies_acc1[idx];
      acc[i][0] = static_cast<double>(dev_acc.x);
      acc[i][1] = static_cast<double>(dev_acc.y);
      acc[i][2] = static_cast<double>(dev_acc.z);
      jrk[i][0] = 0.0;
      jrk[i][1] = 0.0;
      jrk[i][2] = 0.0;
   }

   SYSTEM->sync_exec();

}

}


extern "C" {
    void gpunb_devinit_() {
        devinit = true;
    }
    void gpu_stat_() {
    }
    void gpunb_open_(int* nmax, int* nbmax) {
        NMAX = *nmax;
        NBMAX = *nbmax;
        THETA = double_env("THETA");
        EPS = double_env("EPS");
        ID2IDX.reserve(*nmax);

        if (!SYSTEM) {
            SYSTEM = new ::p3tree(THETA, EPS);

            SYSTEM->load_kernels();

            SYSTEM->localTree.setN(NMAX);
            SYSTEM->localTree.n_nbmax = NBMAX;
            SYSTEM->allocateParticleMemory(SYSTEM->localTree);

            CU_SAFE_CALL(cudaEventCreate(&startLocalGrav));
            CU_SAFE_CALL(cudaEventCreate(&endLocalGrav));
            CU_SAFE_CALL(cudaEventCreate(&startRemoteGrav));
            CU_SAFE_CALL(cudaEventCreate(&endRemoteGrav));

            FIRST_SORT = true;
        }
        else {
            SYSTEM->localTree.setN(NMAX);
            SYSTEM->localTree.n_nbmax = NBMAX;
            SYSTEM->reallocateParticleMemory(SYSTEM->localTree);
        }

#if PROFILE
        total_time = tick_type(0);
        num_ops = 0;
#endif

        std::cerr << "***********************" << std::endl
                  << "Opened NBODY6/Bonsai library" << std::endl
                  << "  NMAX " << NMAX << std::endl
                  << "  NBMAX " << NBMAX << std::endl
                  << "  THETA " << std::fixed << std::setprecision(4) << THETA << std::endl
                  << "  EPS " << std::fixed << std::setprecision(4) << EPS << std::endl
                  << "***********************" << std::endl;
    }

    void gpunb_close_() {
#if PROFILE
        std::cerr << "**************************" << std::endl
                  << "Closed NBODY6/Bonsai library" << std::endl
                  << " time: "
                  << std::chrono::duration<double>(total_time).count() << std::endl
                  << " ops:  " << num_ops << std::endl
                  << " avg:  "
                  << std::chrono::duration<double>(total_time).count() / num_ops
                  << std::endl
                  << "**************************" << std::endl;
#endif
    }

   void gpunb_lf_(
        int *ni,
        double mi[],
        double xi[][3],
        double vi[][3],
        double acc[][3],
        double jrk[][3],
        double phi[],
        double rcut[]) {

        auto start = now();

        ::approx_grav(ni, mi, xi, vi, acc, jrk, rcut);

#pragma omp parallel for schedule(static)
        for (int i = 0; i < *ni; i++) {
           int idx = ID2IDX[i];
           const real4& dev_acc = SYSTEM->localTree.bodies_acc1[idx];
           phi[i] = static_cast<double>(dev_acc.w);
        }

        total_time += now() - start;
        num_ops++;
   }


    void gpunb_regf_(
        int *ni,
        double mi[],
        double xi[][3],
        double vi[][3],
        double acc[][3],
        double jrk[][3],
        int *lmax,
        int *list,
        double rcut[])
    {
        auto start = now();

        ::approx_grav(ni, mi, xi, vi, acc, jrk, rcut);

        int dst_lst_size = *lmax;
        int src_lst_size = SYSTEM->localTree.n_nbmax;

#pragma omp parallel for schedule(static)
        for (int i = 0; i < *ni; i++) {
           int idx = ID2IDX[i];
           int nnb = static_cast<int>(SYSTEM->localTree.bodies_dens[idx].y);
           if (nnb > NBMAX) {
              nnb = -nnb;
           }
           list[i*dst_lst_size] = nnb;
        }
        SYSTEM->sync_copy();

#pragma omp parallel for
        for (int i = 0; i < *ni; i++) {
           int idx = ID2IDX[i];
           int nnb = list[i*dst_lst_size];
           if (nnb < 1) continue;

           int* dst_lst = &list[i*dst_lst_size+1];
           int* src_lst = &SYSTEM->localTree.neighbours[idx*src_lst_size];

           for (int j = 0; j < nnb; j++) {
              int xfm = SYSTEM->localTree.bodies_ids[src_lst[j]];
              auto bound = std::lower_bound(dst_lst, dst_lst + j, xfm);
              std::move_backward(bound, dst_lst + j, dst_lst + j + 1);
              *bound = xfm;
           }
        }
        total_time += now() - start;
        num_ops++;
    }
}
