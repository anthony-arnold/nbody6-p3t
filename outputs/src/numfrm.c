#include <stdio.h>
#include <stdbool.h>
#include "oerrno.h"
#include "meta.h"
#include "discard.h"
#include "reclen.h"

extern bool ofail();

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

static int data_length(int ntot, int nk) {
    static const int VAR_SIZE = 8;
    static const int SCALAR_SIZE = 8;
    static const int VEC_SIZE = 3 * 8;
    static const int NAME_SIZE = 4;


    return VAR_SIZE * nk +   /* HEADER */
        SCALAR_SIZE * ntot + /* MASS */
        VEC_SIZE * ntot +    /* POS */
        VEC_SIZE * ntot +    /* VEL */
        NAME_SIZE * ntot;    /* NAME */
}


static bool discard_data(FILE* fp, int ntot, int nk) {
    /* Discard the data record */
    int datalen = _discard(fp);
    if (datalen < 0) {
        if (!ofail()) {
            /* Error mode not set yet. */
            maybe_eof(fp);
        }
        return false;
    }
    if (datalen != data_length(ntot, nk)) {
        /* Unexpected */
        _oseterrno(OERR_FORMAT);
        return false;
    }

    return true;
}

static bool read_frame(FILE* fp) {
    int ntot, nk;

    /* Read a record length to test for more data. */
    if (_reclen(fp) < 0) {
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
    if (!_meta(fp, &ntot, &nk)) {
        if (!ofail()) {
            /* Error mode not set yet. */
            maybe_eof(fp);
        }
        return false;
    }

    return discard_data(fp, ntot, nk);
}

int numfrm(FILE* fp) {
    int num = 0;

    /* Clear error flag */
    _oseterrno(OERR_SUCCESS);

    if (!fp) {
        _oseterrno(OERR_NULL);
        return -1;
    }

    /* Reset the stream */
    if(0 == fseek(fp, 0, SEEK_SET)) {
        /* Read each frame */
        while (read_frame(fp)) {
            num++;
        }
    }

    return num;
}
