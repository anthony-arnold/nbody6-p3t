/**
 * A simple program to output some critical data about each frame.
 */

#include "outputs.h"
#include <stdlib.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char* argv[]) {

    if (argc < 2) {
        fprintf(stderr, "Usage: info <file>\n");
        return 1;
    }

    FILE* fp = opnout(argv[1]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

    printf("T     R*   <M>   V*\n");
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

            printf("%lf   %lf   %lf   %lf",
                   frame->hdr->t * frame->hdr->tscale,
                   frame->hdr->rbar,
                   mass,
                   frame->hdr->vstar);

            freefrm(frame);
        }
    }

    fclose(fp);
    return 0;

}
