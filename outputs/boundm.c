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

double bound_m(const struct frm_t* frame) {
    double rt = rtide(0.5 * frame->hdr->tidal4);
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
            if (mass == 0) {
                for (int i = 0; i < frame->ntot; i++) {
                    mass += frame->ptcls[i].m;
                }
            }

            double m = bound_m(frame);
            printf("%lf    %lf\n",
                   frame->hdr->t * frame->hdr->tscale,
                   m / mass);

            freefrm(frame);
        }
    }

    fclose(fp);
    return 0;
}
