#include <stdlib.h>
#include "outputs.h"
#include "oerrno.h"

bool frmhdr(FILE* fp, int nk, struct frm_hdr_t** hdr) {
    size_t size = sizeof(double) * nk;
    struct frm_hdr_t* tmp;

    /* Clear error flag */
    _oseterrno(OERR_SUCCESS);

    if (!fp || !hdr) {
        _oseterrno(OERR_NULL);
        return false;
    }
    *hdr = NULL;

    if (nk < 1) {
        _oseterrno(OERR_PARAM);
        return false;
    }

    if (size < sizeof(struct frm_hdr_t)) {
        size = sizeof(struct frm_hdr_t);
    }
    tmp = malloc(size);
    if (!tmp) {
        _oseterrno(OERR_OOM);
        return false;
    }

    if (fread(tmp->as, sizeof(double), nk, fp) != (size_t)nk) {
        if (feof(fp)) {
            _oseterrno(OERR_FRAME_TRUNC);
        }
        else {
            _oseterrno(OERR_READ);
        }
        free(tmp);
        return false;
    }

    *hdr = tmp;
    return true;
}
