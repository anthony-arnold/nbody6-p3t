/**
 * Read each frame and write an output file containing, in each row:
 * x, y, z, m, bound
 */

#include "outputs.h"
#include <stdlib.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

void snapshot(const char* prefix, const char* galmodel, const struct frm_t* frame) {
    /* Create the output file name. */
    char snapname[256];
    double gyr = frame->hdr->t * frame->hdr->tscale;
    snprintf(snapname, 255, "%s_%5.2lf.frame", prefix, gyr);

    double clm = clmass(frame);
    double rt = rtide(frame->hdr, galmodel, clm);
    double rt2 = rt * rt;

    FILE* snap = fopen(snapname, "w");
    for (int i = 0; i < frame->ntot; i++) {
        int bound = 0;
        double x[3], y[3];
        vsub(frame->ptcls[i].x, frame->hdr->rdens, x);
        if (mag2(x) < rt2) {
            bound = 1;
        }

        vmul(frame->ptcls[i].x, frame->hdr->rbar, x);
        vmul(frame->hdr->rgal, frame->hdr->rbar, y);
        vsub(y, x, x);

        fprintf(snap, "%lg %lg %lg %d\n",
                x[0],
                x[1],
                frame->ptcls[i].m * frame->hdr->zmbar,
                bound);
    }
    fclose(snap);
}

int main(int argc, char* argv[]) {
    int skipc = 0;
    int skip = 0;
    if (argc < 3) {
        fprintf(stderr, "Usage: snapshot [NONE|IRRGANG|BOVY] <file> [skip]\n");
        return 1;
    }
    if (argc > 3) {
        skip = atoi(argv[3]);
        if (skip < 0) {
            fprintf(stderr, "bad skip %s\n", argv[2]);
            return 1;
        }
    }

    FILE* fp = opnout(argv[2]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

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
            snapshot(argv[2], argv[1], frame);
            freefrm(frame);
        }
    }

    fclose(fp);
    return 0;

}
