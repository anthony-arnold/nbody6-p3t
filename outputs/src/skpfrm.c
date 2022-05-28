#include <stdio.h>
#include <stdbool.h>
#include "reclen.h"
#include "oerrno.h"

extern bool ofail();

static void maybe_eof(FILE* fp) {
    if (feof(fp)) {
        _oseterrno(OERR_FRAME_TRUNC);
    }
}

static void skprec(FILE* fp) {
    long r = _reclen(fp, -1);
    if (r < 0) {
        maybe_eof(fp);
        return;
    }

    /* skip record */
    if (fseek(fp, r, SEEK_CUR)) {
        _oseterrno(OERR_READ);
        return;
    }

    if (_reclen(fp, r) != r) {
        maybe_eof(fp);
    }
}

void skpfrm(FILE* fp) {
    skprec(fp);
    if (!ofail()) {
        skprec(fp);
    }
}
