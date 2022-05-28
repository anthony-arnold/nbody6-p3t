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

double find_rt(const struct frm_t* frame) {
    /* Find external force. */
    double fx[3], fdx[3];
    forceir13(frame->hdr,
              frame->hdr->rgal,
              frame->hdr->vgal,
              fx,
              fdx);

    /* Find mass and circular velocity corresponding to fx. */


    /* Find omega for circular velocity. */

    double omega = 0.5 * frame->hdr->tidal4;

    return rtide(omega);
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

            double rt = find_rt(frame);
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
