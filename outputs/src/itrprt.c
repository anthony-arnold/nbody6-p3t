#include <stdlib.h>

#include "outputs.h"
#include "meta.h"
#include "oerrno.h"
#include "reclen.h"
#include "data.h"

extern bool frmhdr(FILE* fp, int nk, struct frm_hdr_t** hdr);

void itrprt(FILE* fp, long ptr, fn_frmrdr_t cb, void* dat) {
    int ntot, nk;
    struct frm_hdr_t* header;

    /* Clear error flag */
    _oseterrno(OERR_SUCCESS);

    if (!fp) {
        _oseterrno(OERR_NULL);
        return;
    }

    if (ptr >= 0) {
        if (fseek(fp, ptr, SEEK_SET)) {
            _oseterrno(OERR_READ);
            return;
        }
    }

    /* Read the frame meta. */
    if (!_meta(fp, &ntot, &nk)) {
        return;
    }

    /* Read the opening record length. */
    if (_reclen(fp, frmsz(ntot, nk)) < 0) {
        return;
    }

    /* Read the frame header. */
    if (!frmhdr(fp, nk, &header)) {
        return;
    }

    /* --- iterate ---- */
    for (int n = 0; n < ntot; n++) {
        double cmpnt[7];
        if (!_data(fp, cmpnt)) {
            break;
        }

        if (!cb(ntot, nk, header, cmpnt, dat)) {
            break;
        }
    }

    /* Clean up */
    free(header);
}
