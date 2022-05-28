#include <stdio.h>
#include "reclen.h"
#include "oerrno.h"

int _discard(FILE* fp) {
    /* Read the record length. */
    int len = _reclen(fp, 0);
    if (len < 0) {
        return -1;
    }

    /* discard the record. */
    if (fseek(fp, len, SEEK_CUR) != 0) {
        return -1;
    }

    /* Discard the trailing record length. */
    return _reclen(fp, len);
}
