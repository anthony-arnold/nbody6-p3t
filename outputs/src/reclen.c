#include <stdio.h>
#include <stdint.h>
#include "oerrno.h"

long _reclen(FILE* fp, long expect) {
    int_least32_t len;
    size_t size = fread(&len, 4, 1, fp);

    if (size != 1) {
        return -1;
    }
    if (len < 0) {
        _oseterrno(OERR_PANIC);
        return -1;
    }
    if (expect > 0 && expect != len) {
        _oseterrno(OERR_FORMAT);
        return -1;
    }
    return len;
}
