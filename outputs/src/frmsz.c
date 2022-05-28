/**
 * Determine the size, in bytes, of the frame data.
 */
long frmsz(int ntot, int nk) {
    static const long VAR_SIZE = 8;
    static const long SCALAR_SIZE = 8;
    static const long VEC_SIZE = 3 * 8;
    static const long NAME_SIZE = 4;

    return VAR_SIZE * nk +   /* HEADER */
        SCALAR_SIZE * ntot + /* MASS */
        VEC_SIZE * ntot +    /* POS */
        VEC_SIZE * ntot +    /* VEL */
        NAME_SIZE * ntot;    /* NAME */
}
