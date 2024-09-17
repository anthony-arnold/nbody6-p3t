/**
  Output the mean velocity at various radii.
 */

#include "outputs.h"
#include <math.h>
#include <stdlib.h>

#define BINS 25
#define MAX 10000.0
#define PI 3.14159

double* d = NULL;
void d_data(int n) {
    double* tmp = realloc(d, n * sizeof(double));
    if (!tmp) {
        fprintf(stderr, "Out of memory\n");
        exit(1);
    }
    d = tmp;
}

void vel(struct frm_t* frame) {
    double r0 = 0.0;
    d_data(frame->ntot);

    double y[3];
    vmul(frame->hdr->rgal, frame->hdr->rbar, y);

#pragma omp parallel for
    for (int i = 0; i < frame->ntot; i++) {
        double x[3];
        vmul(frame->ptcls[i].x, frame->hdr->rbar, x);
        vsub(y, x, x);
        d[i] = sqrt(mag2(x));
    }
    for (int i = 0; i < BINS; i++) {
        double r = (MAX / BINS) * (i + 1);

        double v = 0.0;
        int n = 0;
#pragma omp parallel for reduction(+:v,n)
        for (int j = 0; j < frame->ntot; j++) {
            if (d[j] <= r&& d[j] > r0) {
                double x[3];
                n++;
                vmul(frame->ptcls[i].dx, frame->hdr->vstar, x);
                v += sqrt(mag2(x));
            }
        }
        double area = (r*r - r0*r0) * PI;
        v = (v*1000) / frame->ntot / area;

        printf(" %lg %d", v, n);

        r0 = r;
    }
    printf("\n");
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
            vel(frame);
        }
    }

    if (frame) {
        freefrm(frame);
    }
    fclose(fp);
    return 0;
}
