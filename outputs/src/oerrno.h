#ifndef OERRNO_HEADER_INCLUDED
#define OERRNO_HEADER_INCLUDED

typedef enum {
    OERR_SUCCESS,
    OERR_PANIC,
    OERR_FILE,
    OERR_NULL,
    OERR_READ,
    OERR_FRAME_TRUNC,
    OERR_BAD_META,
    OERR_RECLEN,
    OERR_FORMAT,
    OERR_PARAM,
    OERR_OOM,
    OERR_GALMODEL,

    /* MUST GO LAST */
    OERR_NUM_ERRORS
} OError;

extern int _ogeterrno();
extern void _oseterrno(int e);

#endif
