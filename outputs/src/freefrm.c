#include <stdlib.h>
#include "outputs.h"

void freefrm(struct frm_t* frm) {
    if (frm) {
        if (frm->hdr) {
            free(frm->hdr);
        }
        frm->hdr = NULL;
        free(frm);
    }
}
