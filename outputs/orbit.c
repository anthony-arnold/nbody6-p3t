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

    while (!feof(fp)) {
        struct frm_t* frame = rdfrm(fp, -1);
        if (ofail()) {
            fprintf(stderr, "%s\n", oerror());
            return 1;
        }
        if (frame) {
            printf("%lf   ", frame->hdr->t * frame->hdr->tscale);
            for (int i = 0; i < 3; i++) {
                printf("%lf   ", frame->hdr->rgal[i] * frame->hdr->rbar);
            }
            for (int i = 0; i < 3; i++) {
                printf("%lf   ", frame->hdr->vgal[i] * frame->hdr->vstar);
            }
            printf("\n");

            freefrm(frame);
        }
    }

    fclose(fp);
    return 0;

}
