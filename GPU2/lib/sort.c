#include <stdlib.h>

int comp_dbl(const void* l, const void* r) {
    double a = (*(double*)l) - (*(double*)r);
    if (a < 0) {
        return -1;
    }
    if (a > 0) {
        return 1;
    }
    return 0;
}

void sort_(int* n, double* v) {
    qsort(v, *n, sizeof(double), &comp_dbl);
}
