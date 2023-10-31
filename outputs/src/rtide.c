#include <math.h>
#include "outputs.h"

double rtide(const struct frm_hdr_t* hdr, const char* galmodel, double mass) {
    /* Find external force. */
    double fx[3], fdx[3];
    galforce(hdr,
             galmodel,
             hdr->rgal,
             hdr->vgal,
             fx,
             fdx);

    /* Find mass and circular velocity corresponding to fx. */
    double r[3];
    vsub(hdr->rgal, hdr->rdens, r);
    double r2 = mag2(r);
    double m = sqrt(mag2(fx)) * r2;
    double vc = sqrt(m / sqrt(r2));
    double rt = cbrt(mass / (2.0 * vc*vc)) * pow(sqrt(mag2(hdr->rgal)), 2.0/3.0);

    return rt;
}
