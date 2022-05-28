#include <math.h>
#include <string.h>
#include "outputs.h"

static const double m1 = 9.5e9,
    b1 = 230.0,
    m2 = 6.6e10,
    a2 = 4220.0,
    b2 = 292.0,
    m3 = 2.335e10,
    a3 = 2562.0;

static void bulge(const double ru[3], const double vu[3],
                  double f[3], double fd[3]) {

    const double dis = sqrt(mag2(ru));
    const double rv = dot(ru, vu);
    const double xd1 = sqrt(b1*b1 + dis*dis);
    const double xd3 = xd1*xd1*xd1;
    const double xd5 = xd3*xd1*xd1;

    for (int i = 0; i < 3; i++) {
        f[i] += -m1 * ru[i] / xd3;
    }
    for (int i = 0; i < 3; i++) {
        fd[i] += -m1 * (vu[i] / xd3 - 3.0 * rv * ru[i] / xd5);
    }
}

static void disk(const double ru[3], const double vu[3],
                 double f[3], double fd[3]) {

    const double rv = dot(ru, vu);
    const double xc3 = sqrt(ru[2]*ru[2] + b2*b2);
    const double xc2 = a2 + xc3;
    const double xc1 = sqrt(ru[0]*ru[0] + ru[1]*ru[1] + xc2*xc2);

    const double f2[] = {
        -m2 * ru[0] / pow(xc1, 3.0),
        -m2 * ru[1] / pow(xc1, 3.0),
        -m2 * ru[2] / pow(xc1, 3.0) / xc3 * xc2
    };
    const double fd2[] = {
        -m2 * vu[0] / pow(xc1, 3.0) + 3*m2 * ru[0] / pow(xc1, 5.0) *
        (rv + ru[2] * vu[2] * a2 / xc3),

        -m2 * vu[1] / pow(xc1, 3.0) + 3*m2 * ru[1] / pow(xc1, 5.0) *
        (rv + ru[2] * vu[2] * a2 / xc3),

        -m2 * vu[2] / pow(xc1, 3.0) + m2*xc2/xc3/pow(xc1, 5.0) *
        (3 * ru[2] + rv + ru[2]*ru[2] + vu[2]*a2 / xc3) +
        m2 * vu[2]/pow(xc1, 3.0) * (ru[2]*ru[2] * a2/pow(xc3, 3.0) - xc2/xc3)
    };

    for (int i = 0; i < 3; i++) {
        f[i] += f2[i];
        fd[i] += fd2[i];
    }
}

static void halo(const double ru[3], const double vu[3],
                 double f[3], double fd[3]) {

    const double dis = sqrt(mag2(ru));
    const double rv = dot(ru, vu);
    const double fx = (dis*a3*a3 + dis*dis*a3);
    const double fx2 = fx*fx;
    const double fy = (a3*a3/dis+2.0*a3);

    for (int i = 0; i < 3; i++) {
        f[i] += -m3 * ru[i] / fx;
    }
    for (int i = 0; i < 3; i++) {
        fd[i] += -m3 * vu[i] / fx + m3 * ru[i] * rv / fx2 * fy;
    }
}

void forceir13(const struct frm_hdr_t* hdr,
               const double rg[3], const double vg[3],
               double fp[3], double fd[3])
{
    const double ru[] = {
        rg[0] * hdr->rbar,
        rg[1] * hdr->rbar,
        rg[2] * hdr->rbar
    };
    const double vu[] = {
        vg[0] * hdr->vstar,
        vg[1] * hdr->vstar,
        vg[2] * hdr->vstar
    };
    const double rbar2 = hdr->rbar*hdr->rbar;
    const double zmbar = hdr->zmbar;
    const double vstar = hdr->vstar;

    for (int i = 0; i < 3; i++) {
        fp[i] = 0.0;
        fd[i] = 0.0;
    }
    bulge(ru, vu, fp, fd);
    disk(ru, vu, fp, fd);
    halo(ru, vu, fp, fd);

    for (int i = 0; i < 3; i++) {
        fp[i] = fp[i] / zmbar * rbar2;
    }
    for (int i = 0; i < 3; i++) {
        fd[i] = fd[i] / zmbar * rbar2 / vstar;
    }
}
