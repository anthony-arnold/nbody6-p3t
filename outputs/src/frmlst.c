#include <stdio.h>
#include <stdbool.h>
#include "oerrno.h"

extern bool ofail();
extern int frmsz(int ntot, int nk);
extern bool frmmeta(FILE* fp, long* ptr, int* ntot, int* nk);

void frmlst(FILE* fp, long* ptrs, int bufsz, int* nptrs) {
    int num = 0;

    /* Clear error flag */
    _oseterrno(OERR_SUCCESS);

    if (!fp || !nptrs) {
        _oseterrno(OERR_NULL);
        return;
    }

    /* Reset the stream */
    if(0 == fseek(fp, 0, SEEK_SET)) {
        /* Read each frame */
        while (frmmeta(fp, bufsz > num ? ptrs : NULL, NULL, NULL)) {
            num++;
            if (ptrs) {
                ptrs++;
            }
        }
    }
    *nptrs = num;
}
