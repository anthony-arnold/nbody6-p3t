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

void print_bound_masses(FILE* fp, int skip) {
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
            exit(1);
        }
        else if (!frame) {
            break;
        }

        double m, n, rt;
        bound(frame, &m, &n, &rt);
        printf("%lf    %lf    %lf     %lf\n",
               frame->hdr->t * frame->hdr->tscale,
               rt * frame->hdr->rbar,
               m,
               n);

        freefrm(frame);
    }
}

int main(int argc, char* argv[]) {
    int skip = 0;
    if (argc < 2) {
        fprintf(stderr, "Usage: rtide <file> [skip]\n");
        return 1;
    }
    if (argc > 2) {
        skip = atoi(argv[2]);
        if (skip < 0) {
            fprintf(stderr, "bad skip %s\n", argv[2]);
            return 1;
        }
    }

    FILE* fp = opnout(argv[1]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }
    print_bound_masses(fp, skip);

    fclose(fp);
    return 0;
}
