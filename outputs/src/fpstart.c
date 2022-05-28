#include <stdio.h>
#include <stdbool.h>
#include "oerrno.h"

bool _fpstart(FILE* fp, long ptr) {
    if (ptr >= 0) {
        if (fseek(fp, ptr, SEEK_SET)) {
            _oseterrno(OERR_READ);
            return false;
        }
    }
    return true;
}
