#include <stdio.h>
#include <stdbool.h>

extern void frmlst(FILE* fp, long* ptrs, int bufsz, int* nptrs);

int numfrm(FILE* fp) {
    int num = 0;
    frmlst(fp, NULL, 0, &num);
    return num;
}
