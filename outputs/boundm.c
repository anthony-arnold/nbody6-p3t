/**
 * A program to determine the tidal radius using the
 * circular velocity of the cluster. Will print the
 * tidal radius for each frame in the output.
 */

#include "outputs.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

double find_rt(const struct frm_t* frame, double mass) {
    /* Find external force. */
    double fx[3], fdx[3];
    forceir13(frame->hdr,
              frame->hdr->rgal,
              frame->hdr->vgal,
              fx,
              fdx);

    /* Find mass and circular velocity corresponding to fx. */
    double r[3];
    vsub(frame->hdr->rgal, frame->hdr->rdens, r);
    double r2 = mag2(r);
    double m = sqrt(mag2(fx)) * r2;
    double vc = sqrt(m / sqrt(r2));
    double rt = cbrt(mass / (2.0 * vc*vc)) * pow(sqrt(mag2(frame->hdr->rgal)), 2.0/3.0);

    return rt;
}

double bound_m(const struct frm_t* frame, double mass) {
    /*
      double rt = rtide(0.5 * frame->hdr->tidal4);
    */
    double rt = find_rt(frame, mass);
    double rt2 = rt*rt;

    double m = 0.0;
    for (int i = 0; i < frame->ntot; i++) {
        double x[3];
        vsub(frame->ptcls[i].x, frame->hdr->rdens, x);

        if (mag2(x) < rt2) {
            m += frame->ptcls[i].m;
        }
    }
    return m;
}

int main(int argc, char* argv[]) {
    int skip = 0;
    if (argc < 2) {
        fprintf(stderr, "Usage: rtide <file>\n");
        return 1;
    }
    if (argc > 2) {
        skip = atoi(argv[2]);
    }

    FILE* fp = opnout(argv[1]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

    double mass = 0;
    int skipc = 0;
    while (!feof(fp)) {
        if (skipc--) {
            skpfrm(fp);
            continue;
        }
        skipc = skip;

        struct frm_t* frame = rdfrm(fp, -1);
        if (ofail()) {
            fprintf(stderr, "%s\n", oerror());
            return 1;
        }
        if (frame) {
            if (mass == 0) {
                for (int i = 0; i < frame->ntot; i++) {
                    mass += frame->ptcls[i].m;
                }
            }

            double m = bound_m(frame, mass);
            printf("%lf    %lf\n",
                   frame->hdr->t * frame->hdr->tscale,
                   m / mass);

            freefrm(frame);
        }
    }

    fclose(fp);
    return 0;
}
