#include "oerrno.h"

static const char* errstrs[OERR_NUM_ERRORS] = {
    "No error",
    "General badness",
    "Failed to open file",
    "Unexpected NULL pointer",
    "Read error",
    "Frame truncated",
    "Bad header"
};

const char* _oerrstr() {
    int e = _ogeterrno();

    // This should never happen.
    if (e < 0 || e >= OERR_NUM_ERRORS) {
        e = OERR_PANIC;
        _oseterrno(e);
    }

    return errstrs[e];
}
