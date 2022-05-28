#include <stdio.h>
#include <stdbool.h>

#include "oerrno.h"

bool _data (FILE* fp, size_t n, double* out) {
    if (fread(out, sizeof(double), n, fp) != n) {
        if (feof(fp)) {
            _oseterrno(OERR_FRAME_TRUNC);
        }
        else {
            _oseterrno(OERR_READ);
        }
        return false;
    }
    return true;
}
