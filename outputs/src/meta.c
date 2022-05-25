#include <stdbool.h>
#include <stdio.h>
#include "readrec.h"
#include "oerrno.h"

extern bool ofail();

bool _meta(FILE* fp, int* ntot, int* nk) {
    char meta[16];
    int len;

    /* Read the meta data */
    if (_readrec(fp, meta, sizeof(meta), &len) < 0) {
        /* If the buffer isn't big enough, that's a major problem. */
        if (!ofail() && (size_t)len > sizeof(meta)) {
            _oseterrno(OERR_FORMAT);
        }
        return false;
    }

    *ntot = *(int*)&meta[0];
    if (*ntot < 0) {
        _oseterrno(OERR_BAD_META);
        return false;
    }

    *nk = *(int*)&meta[12];
    if (*nk < 0) {
        _oseterrno(OERR_BAD_META);
        return false;
    }

    return true;
}
