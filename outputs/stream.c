/**
  Output a stream of frame data onto stdout for consumption
  by 3rd part processing (see movie.py).
 */

#include "outputs.h"
#include <stdlib.h>

void stream(struct frm_t* frame) {
    /* Output frame time at start of frame */
    printf("T = %lf\n", frame->hdr->t * frame->hdr->tscale);

    double y[3];
    vmul(frame->hdr->rgal, frame->hdr->rbar, y);

    for (int i = 0; i < frame->ntot; i++) {
        double x[3];
        double x2[3];
        vmul(frame->ptcls[i].x, frame->hdr->rbar, x);
        vsub(y, x, x2);

        printf("%lf %lf %lf %lf %lf\n",
               frame->ptcls[i].m * frame->hdr->zmbar,
               x[0],
               x[1],
               x2[0],
               x2[1]);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: stream <file>\n");
        return 1;
    }

    FILE* fp = opnout(argv[1]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

    struct frm_t* frame = NULL;
    while (!feof(fp)) {
        struct frm_t* tmp = rdfrm2(fp, -1, frame);
        if (ofail()) {
            fprintf(stderr, "%s\n", oerror());
            return 1;
        }
        if (tmp) {
            frame = tmp;
            stream(frame);
        }
    }

    if (frame) {
        freefrm(frame);
    }
    fclose(fp);
    return 0;
}
