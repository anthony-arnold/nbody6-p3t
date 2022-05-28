#include <assert.h>
#include <stdlib.h>
#include "outputs.h"
#include "oerrno.h"
#include "meta.h"
#include "reclen.h"

extern bool _hdr(FILE* fp, int nk, struct frm_hdr_t** hdr);
extern bool _fpstart(FILE* fp, long ptr);

void frmhdr(FILE* fp, long ptr, int* ntot, int* nk,
            struct frm_hdr_t** hdr) {
    int _ntot, _nk;
    long size, fstart;
    /* Clear error flag */
    _oseterrno(OERR_SUCCESS);

    if (!fp || !hdr) {
        _oseterrno(OERR_NULL);
        return;
    }
    *hdr = NULL;

    if (!_fpstart(fp, ptr)) {
        return;
    }

    /* Read the frame meta. */
    if (!_meta(fp, &_ntot, &_nk)) {
        /* If meta fails, then either we're at the end of the file,
           or meta will set the error condition. */
        assert(feof(fp) || ofail());
        return;
    }
    if (ntot) {
        *ntot = _ntot;
    }
    if (nk) {
        *nk = _nk;
    }

    /* Read the opening record length. */
    size = _reclen(fp, 0);
    if (size < 0) {
        return;
    }
    fstart = ftell(fp);

    /* Read the frame header. */
    if (!_hdr(fp, _nk, hdr)) {
        return;
    }
    /* Get to the end of the frame if not already there. */
    if (!ofail()) {
        if (fseek(fp, fstart + size, SEEK_SET)) {
            _oseterrno(OERR_READ);
        }
        else {
            _reclen(fp, size);
        }
    }
}
