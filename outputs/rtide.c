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

//static const double G = 6.6743e-08;

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
    double rt = cbrt(mass / 2 * vc*vc) * pow(sqrt(mag2(frame->hdr->rgal)), 2.0/3.0);

    return rt;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: rtide <file>\n");
        return 1;
    }

    FILE* fp = opnout(argv[1]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

    while (!feof(fp)) {
        double mass = 0;
        struct frm_t* frame = rdfrm(fp, -1);
        if (ofail()) {
            fprintf(stderr, "%s\n", oerror());
            return 1;
        }
        if (frame) {
            for (int i = 0; i < frame->ntot; i++) {
                mass += frame->ptcls[i].m;
            }

            double rt = find_rt(frame, mass);
            printf("%lf    %lf [%lf]\n",
                   frame->hdr->t * frame->hdr->tscale,
                   rt,
                   frame->hdr->rtide);

            freefrm(frame);
        }
    }

    fclose(fp);
    return 0;
}
