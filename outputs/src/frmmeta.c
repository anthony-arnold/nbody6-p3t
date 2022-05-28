#include <stdio.h>
#include <stdbool.h>
#include "oerrno.h"
#include "meta.h"
#include "discard.h"
#include "reclen.h"

extern bool ofail();

/**
 * One of the read functions failed in the middle of the frame.
 * Either way it's an error, but which one?
 */
static void maybe_eof(FILE* fp) {
    /* Could be eof */
    if (feof(fp)) {
        _oseterrno(OERR_FRAME_TRUNC);
    }
    /* Generic read failure. */
    else {
        _oseterrno(OERR_READ);
    }
}

static bool discard_data(FILE* fp) {
    /* Discard the data record */
    int datalen = _discard(fp);
    if (datalen < 0) {
        if (!ofail()) {
            /* Error mode not set yet. */
            maybe_eof(fp);
        }
        return false;
    }

    return true;
}

bool frmmeta(FILE* fp, long* ptr, int* ntot, int* nk) {
    int _ntot, _nk;

    /* Clear error flag */
    _oseterrno(OERR_SUCCESS);

    if (!fp) {
        _oseterrno(OERR_NULL);
        return false;
    }

    /* Record the file location. */
    if (ptr) {
        *ptr = (int)ftell(fp);
        if (*ptr < 0) {
            _oseterrno(OERR_READ);
            return false;
        }
    }

    /* Read a record length to test for more data. */
    if (_reclen(fp, 0) < 0) {
        /* Check for read failure.
           EOF indicates no more frames and isn't an error. */
        if (!feof(fp)) {
            _oseterrno(OERR_READ);
        }
        else {
            _oseterrno(OERR_SUCCESS);
        }
        return false;
    }

    /* Rewind the record length so readrec can get it. */
    if (fseek(fp, -4, SEEK_CUR) != 0) {
        _oseterrno(OERR_READ);
        return false;
    }

    /* Read the metadata */
    if (!_meta(fp, &_ntot, &_nk)) {
        if (!ofail()) {
            /* Error mode not set yet. */
            maybe_eof(fp);
        }
        return false;
    }

    /* Return the data if requested */
    if (ntot) {
        *ntot = _ntot;
    }
    if (nk) {
        *nk = _nk;
    }

    return discard_data(fp);
}
