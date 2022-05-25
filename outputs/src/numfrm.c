#include <stdio.h>
#include <stdbool.h>
#include "oerrno.h"
#include "reclen.h"

static bool read_chunk(void* buffer, size_t size, FILE* fp) {
    size_t len = fread(buffer, size, 1, fp);
    if (len != 1) {
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

static bool discard(FILE* fp, long size) {
    return fseek(fp, size, SEEK_CUR) == 0;
}

static size_t frame_length(int ntot, int nk) {
    static const size_t VAR_SIZE = 8;
    static const size_t SCALAR_SIZE = 8;
    static const size_t VEC_SIZE = 3 * 8;
    static const size_t NAME_SIZE = 4;


    return VAR_SIZE * nk +   /* HEADER */
        SCALAR_SIZE * ntot + /* MASS */
        VEC_SIZE * ntot +    /* POS */
        VEC_SIZE * ntot +    /* VEL */
        NAME_SIZE * ntot;    /* NAME */
}

static bool read_frame(FILE* fp) {
    char header[16];
    int ntot, nk;

    /* Read the meta data */
    if (!read_chunk(header, 16, fp)) {
        return false;
    }

    ntot = *(int*)&header[0];
    if (ntot < 0) {
        _oseterrno(OERR_BAD_HEADER);
        return false;
    }

    nk = *(int*)&header[4];
    if (nk < 0) {
        _oseterrno(OERR_BAD_HEADER);
        return false;
    }

    /* Read the final record length of the meta data. */
    if (!_reclen(fp)) {
        return false;
    }

    /* Read the opening record length of the data. */
    if (!_reclen(fp)) {
        return false;
    }

    /* Discard the frame */
    if (!discard(fp, frame_length(ntot, nk))) {
        return false;
    }

    /* Read the final record length of the data. */
    if (!_reclen(fp)) {
        return false;
    }
    return true;
}

int numfrm(FILE* fp) {
    int num = 0;
    bool ok = true;

    if (!fp) {
        _oseterrno(OERR_NULL);
        return -1;
    }

    /* Reset the stream */
    if(0 == fseek(fp, 0, SEEK_SET)) {
        /* Read each frame */
        while (ok && _reclen(fp)) {
            ok = read_frame(fp);
        }
    }
    if (!ok) {
        if (feof(fp)) {
            _oseterrno(OERR_FRAME_TRUNC);
        }
        else if (ferror(fp)) {
            _oseterrno(OERR_READ);
        }
    }

    return num;
}
