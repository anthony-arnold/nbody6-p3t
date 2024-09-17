#include <stdlib.h>
#include <assert.h>

#include "outputs.h"
#include "oerrno.h"
#include "meta.h"
#include "reclen.h"
#include "data.h"

extern bool _fpstart(FILE* fp, long ptr);
extern bool _hdr(FILE* fp, int nk, struct frm_hdr_t** hdr);

static struct frm_t* alloc(struct frm_t* frame,
                           int ntot,
                           struct frm_hdr_t* header)
{


    /* Allocate enough space for all the data. */
    struct frm_t* tmp = realloc(frame,
                                sizeof(struct frm_t) +
                                sizeof(struct particle_t) * ntot);
    if (tmp == NULL) {
        _oseterrno(OERR_OOM);
        return NULL;
    }
    tmp->hdr = header;
    return tmp;
}

static bool mass(FILE* fp, int ntot, struct frm_t* frame) {
    for (int i = 0; i < ntot; i++) {
        if (!_data(fp, 1, &frame->ptcls[i].m)) {
            return false;
        }
    }
    return true;
}
static bool pos(FILE* fp, int ntot, struct frm_t* frame) {
    for (int i = 0; i < ntot; i++) {
        if (!_data(fp, 3, frame->ptcls[i].x)) {
            return false;
        }
    }
    return true;
}
static bool vel(FILE* fp, int ntot, struct frm_t* frame) {
    for (int i = 0; i < ntot; i++) {
        if (!_data(fp, 3, frame->ptcls[i].dx)) {
            return false;
        }
    }
    return true;
}

static bool readall(FILE* fp, int ntot, struct frm_t* frame) {
    return mass(fp, ntot, frame) && pos(fp, ntot, frame) && vel(fp, ntot, frame);
}

static struct frm_t* getfrm(FILE* fp,
                            int ntot,
                            int nk,
                            struct frm_hdr_t* header,
                            struct frm_t* old) {

    /* Allocate enough space for all the data. */
    struct frm_t* frame = alloc(old, ntot, header);
    if (!frame) {
        free(header);
        return NULL;
    }

    if (!readall(fp, ntot, frame)) {
        freefrm(frame);
        return NULL;
    }

    frame->ntot = ntot;
    frame->nk = nk;
    return frame;
}

static bool fastforward(FILE* fp, long fstart, long size) {
    if (ofail()) {
        return false;
    }

    /* Get to the end of the frame if not already there. */
    if (fseek(fp, fstart + size, SEEK_SET)) {
        _oseterrno(OERR_READ);
        return false;
    }

    _reclen(fp, size);
    return true;
}

struct frm_t* rdfrm(FILE* fp, long ptr) {
    return rdfrm2(fp, ptr, NULL);
}

struct frm_t* rdfrm2(FILE* fp, long ptr, struct frm_t* frame) {

    int ntot, nk;
    long fstart;
    struct frm_hdr_t* header = NULL;
    long size;

    /* Clear error flag */
    _oseterrno(OERR_SUCCESS);

    if (!fp) {
        _oseterrno(OERR_NULL);
        return NULL;
    }

    if (!_fpstart(fp, ptr)) {
        return NULL;
    }

    if (frame) {
        header = frame->hdr;
    }

    /* Read the frame meta. */
    if (!_meta(fp, &ntot, &nk)) {
        /* If meta fails, then either we're at the end of the file,
           or meta will set the error condition. */
        assert(feof(fp) || ofail());
        return NULL;
    }

    /* Read the opening record length. */
    size = _reclen(fp, -1);
    if (size < 0) {
        return NULL;
    }
    fstart = ftell(fp);

    /* Read the frame header. */
    if (!_hdr(fp, nk, &header)) {
        return NULL;
    }
    assert(size == frmsz(ntot, nk, header->kz19));

    /* Allocate enough space for all the data. */
    frame = getfrm(fp, ntot, nk, header, frame);

    /* Get to the end of the frame if not already there. */
    if (!fastforward(fp, fstart, size)) {
        freefrm(frame);
        frame = NULL;
    }
    return frame;
}
