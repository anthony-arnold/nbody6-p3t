#include <stdio.h>
#include <stdbool.h>

#include "oerrno.h"

bool _data (FILE* fp, double ptcl[7]) {
    if (fread(ptcl, sizeof(double), 7, fp) != 7) {
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
