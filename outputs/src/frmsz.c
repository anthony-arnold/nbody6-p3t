#include <stdbool.h>

/**
 * Determine the size, in bytes, of the frame data.
 */
long frmsz(int ntot, int nk, bool kz19) {
    static const long VAR_SIZE = 8;
    static const long SCALAR_SIZE = 8;
    static const long VEC_SIZE = 3 * 8;
    static const long NAME_SIZE = 4;
    static const long KSTAR_SIZE = 4;
    static const long SEV_SIZE = 4;

    long size = VAR_SIZE * nk +   /* HEADER */
        SCALAR_SIZE * ntot + /* MASS */
        VEC_SIZE * ntot +    /* POS */
        VEC_SIZE * ntot +    /* VEL */
        NAME_SIZE * ntot;    /* NAME */

    if (kz19) {
        size += KSTAR_SIZE * ntot + /* KSTAR */
            2 * SEV_SIZE * ntot; /* LSEV & RSEV */
    }
    return size;
}
