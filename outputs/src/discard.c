#include <stdio.h>
#include "reclen.h"
#include "oerrno.h"

int _discard(FILE* fp) {
    /* Read the record length. */
    int len = _reclen(fp);
    if (len < 0) {
        return -1;
    }

    /* discard the record. */
    if (fseek(fp, len, SEEK_CUR) != 0) {
        return -1;
    }

    /* Discard the trailing record length. */
    int trail = _reclen(fp);
    if (trail < 0) {
        return -1;
    }
    if (trail != len) {
        /* Unexpected */
        _oseterrno(OERR_FORMAT);
        return -1;
    }

    return len;
}
