#include <stdio.h>
#include "oerrno.h"
#include "reclen.h"

int _readrec(FILE* fp, char* buff, int buffsize, int *reclen) {
    /* Read the record length. */
    *reclen = _reclen(fp, 0);
    if (*reclen < 0) {
        return -1;
    }

    if (*reclen > buffsize) {
        /* Not enough room. Rewind and report. */
        if (fseek(fp, -4, SEEK_CUR) != 0) {
            _oseterrno(OERR_READ);
        }
        return -1;
    }

    /* Read the record. */
    if (fread(buff, *reclen, 1, fp) != 1) {
        return -1;
    }

    /* Discard the trailing record length. */
    return _reclen(fp, *reclen);
}
