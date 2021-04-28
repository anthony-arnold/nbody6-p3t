#include <algorithm>
#include <cmath>
#include <cassert>
#include "v4df.h"
#include "v4sf.h"

static const v4sf zero = v4sf(0.0);
static const v4sf one = v4sf(1.0);
static const v4sf nf20 = v4sf(-20.0f);
static const v4sf f70 = v4sf(70.0f);
static const v4sf nf84 = v4sf(-84.0f);
static const v4sf f35 = v4sf(35.0f);
static const v4sf f140 = v4sf(140.0f);
static const v4sf f420 = v4sf(420.0f);

void cnbint(
    int i,
    const double pos[][3],
    const double vel[][3],
    const double mass[],
    int nnb,
    int list[],
    double f[3],
    double fdot[3],
    double rcut){
    const int NBMAX = 512;
    // assert(nnb <= NBMAX);
    const double rin = rcut * 0.1;
    const double dricut = 1.0 / (rcut - rin);

    float xbuf[NBMAX] __attribute__ ((aligned(16)));
    float ybuf[NBMAX] __attribute__ ((aligned(16)));
    float zbuf[NBMAX] __attribute__ ((aligned(16)));
    float vxbuf[NBMAX] __attribute__ ((aligned(16)));
    float vybuf[NBMAX] __attribute__ ((aligned(16)));
    float vzbuf[NBMAX] __attribute__ ((aligned(16)));
    float mbuf[NBMAX] __attribute__ ((aligned(16)));
    assert((unsigned long)xbuf % 16 == 0);

    double xi = pos[i][0];
    double yi = pos[i][1];
    double zi = pos[i][2];
    float vxi = vel[i][0];
    float vyi = vel[i][1];
    float vzi = vel[i][2];
    v4df ax(0.0), ay(0.0), az(0.0);
    v4sf jx(0.0), jy(0.0), jz(0.0);
    for(int koff=0; koff<nnb; koff+=NBMAX){
        int nk = std::min(nnb-koff, NBMAX);
        for(int k=0; k<nk; k++){
            int j = list[k+koff];
            assert(j != i);
#if 1
            int jj = list[k+koff+4];
            __builtin_prefetch(pos[jj]);
            __builtin_prefetch(pos[jj]+2);
            __builtin_prefetch(vel[jj]);
            __builtin_prefetch(vel[jj]+2);
            __builtin_prefetch(&mass[jj]);
#endif
            double xj = pos[j][0];
            double yj = pos[j][1];
            double zj = pos[j][2];
            float vxj = vel[j][0];
            float vyj = vel[j][1];
            float vzj = vel[j][2];
            float mj = mass[j];
            xj -= xi;
            yj -= yi;
            zj -= zi;
            vxj -= vxi;
            vyj -= vyi;
            vzj -= vzi;
            xbuf[k] = xj;
            ybuf[k] = yj;
            zbuf[k] = zj;
            vxbuf[k] = vxj;
            vybuf[k] = vyj;
            vzbuf[k] = vzj;
            mbuf[k] = mj;
        }
        for(int k=nk; k%4; k++){
            xbuf[k] = 16.0f;
            ybuf[k] = 16.0f;
            zbuf[k] = 16.0f;
            vxbuf[k] = 0.0f;
            vybuf[k] = 0.0f;
            vzbuf[k] = 0.0f;
            mbuf[k] = 0.0f;
        }


        for(int k=0; k<nk; k+=4){
            const v4sf dx(xbuf + k);
            const v4sf dy(ybuf + k);
            const v4sf dz(zbuf + k);
            const v4sf dvx(vxbuf + k);
            const v4sf dvy(vybuf + k);
            const v4sf dvz(vzbuf + k);
            const v4sf mj(mbuf + k);

            const v4sf r2 = dx*dx + dy*dy + dz*dz;
            const v4sf rv = dx*dvx + dy*dvy + dz*dvz;

            const v4sf rinv1 = r2.rsqrt();
            const v4sf rinv2 = rinv1 * rinv1;
            const v4sf rinv3 = mj * rinv1 * rinv2;

            const v4sf alpha = v4sf(-3.0f) * rinv2 * rv;

            v4sf x = (r2.sqrt() - rin) * dricut;
            v4sf xv = rv * rinv1 * dricut;

            x = __builtin_ia32_andps(x, __builtin_ia32_cmpps(zero, x, 0x1));
            x = v4sf(__builtin_ia32_andps(x-one, __builtin_ia32_cmpps(x, one, 0x1)))+one;
            xv = __builtin_ia32_andps(xv, __builtin_ia32_cmpps(zero, x, 0x1));
            xv = __builtin_ia32_andps(xv, __builtin_ia32_cmpps(x, one, 0x1));

            const v4sf x2 = x*x;
            const v4sf x3 = x2*x;
            const v4sf x4 = x2*x2;

            const v4sf K = (((nf20*x + f70)*x + nf84)*x + f35)*x4;
            const v4sf KD = ((((-f140*x + f420)*x - f420)*x + f140)*x3)*xv;

            const v4sf kmrinv3 = (one - K) * rinv3;
            const v4sf kdmrinv3 = KD * rinv3;

            const v4sf a1 = kmrinv3 * dx;
            const v4sf a2 = kmrinv3 * dy;
            const v4sf a3 = kmrinv3 * dz;

            const v4sf j1 = kmrinv3 * dvx + alpha*a1 - kdmrinv3*dx;
            const v4sf j2 = kmrinv3 * dvy + alpha*a2 - kdmrinv3*dy;
            const v4sf j3 = kmrinv3 * dvz + alpha*a3 - kdmrinv3*dz;

            ax += v4df(a1);
            ay += v4df(a2);
            az += v4df(a3);
            jx += j1;
            jy += j2;
            jz += j3;
        }
    }
    f[0] = ax.sum();
    f[1] = ay.sum();
    f[2] = az.sum();
    fdot[0] = (v4df(jx)).sum();
    fdot[1] = (v4df(jy)).sum();
    fdot[2] = (v4df(jz)).sum();
    assert(f[0] == f[0]);
    assert(f[1] == f[1]);
    assert(f[2] == f[2]);
    assert(fdot[0] == fdot[0]);
    assert(fdot[1] == fdot[1]);
    assert(fdot[2] == fdot[2]);
}

extern "C" {
    void cnbint_(
        int *i,
        double pos[][3],
        double vel[][3],
        double mass[],
        int *nnb,
        int *nblist,
        double f[3],
        double fdot[3],
        double* rcut){
        cnbint(*i, pos-1, vel-1, mass-1, *nnb-1, nblist, f, fdot, *rcut);
    }
}
