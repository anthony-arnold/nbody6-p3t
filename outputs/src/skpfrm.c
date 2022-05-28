#include <stdio.h>
#include <stdbool.h>
#include "reclen.h"
#include "oerrno.h"
#include "discard.h"

extern bool ofail();

static void maybe_eof(FILE* fp) {
    if (feof(fp)) {
        _oseterrno(OERR_FRAME_TRUNC);
    }
    else {
        _oseterrno(OERR_READ);
    }
}

static void skprec(FILE* fp) {
    if (_discard(fp) < 0 && !ofail()) {
        maybe_eof(fp);
    }
}

void skpfrm(FILE* fp) {
    skprec(fp);
    if (!ofail()) {
        skprec(fp);
    }
}
