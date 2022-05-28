#ifndef OUTPUTS_HEADER_INCLUDED
#define OUTPUTS_HEADER_INCLUDED

#include <stdio.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif


#define OEXTRN extern

/**
 * \brief Contains frame header data.
 */
struct frm_hdr_t {
    union {
        struct {
            /**
             * The time in NBODY units.
             */
            double t;

            /**
             * The number of pairs.
             */
            double npairs;

            /**
             * The radius scale factor.
             */
            double rbar;

            /**
             * The mass scale factor.
             */
            double zmbar;

            /**
             * The tidal radius.
             */
            double rtide;

            /**
             * ???
             */
            double tidal4;

            /**
             * The cluster density center.
             */
            double rdens[3];

            /**
             * The physical time.
             */
            double tphys;

            /**
             * The time scale.
             */
            double tscale;

            /**
             * The velocity scale factor.
             */
            double vstar;

            double rc;
            double nc;
            double vc;
            double rhom;
            double cmax;
            double rscale;
            double rsmin;
            double dmin1;

            /**
             * The galaxy particle position in NBODY units.
             */
            double rgal[3];

            /**
             * The galaxy particle velocity in NBODY units.
             */
            double vgal[3];

            double empty[6];

            /**
             * The cluster density center velocity.
             */
            double dvr[3];

            double kz19;
        };
        double as[30];
    };
};

/**
 * \brief A callback function to read every particle in a frame.
 */
typedef bool (*fn_frmrdr_t)(int,
                            int,
                            struct frm_hdr_t*,
                            double[7],
                            void*);

/**
 * \brief Get the current error state.
 * \return Returns true if the previous function call failed.
 */
OEXTRN bool ofail();

/**
 * \brief Get the error message for the current error state.
 * \return A human-readable string indicating the error state.
 */
OEXTRN const char* oerror();

/**
 * \brief Open an output file for reading.
 * \param path The path to the file to open.
 * \return The handle to the output file, or NULL if an error occurred.
 */
OEXTRN FILE* opnout(const char* path);

/**
 * \brief Count the number of frames in the output file.
 * \param fp The output file.
 * \return the number of frames in the file, or a negative integer if an error
 * occurred.
 */
OEXTRN int numfrm(FILE* fp);

/**
 * \brief Return a list of offsets into the file at which each frame begins.
 *
 * Pass in a zero-length buffer to get the required size of `ptrs`, or
 * call `numfrm`.
 *
 * \param fp The output file.
 * \param ptrs A buffer of integers which will be populated with offsets.
 * \param bufsz The maximum number of integers that will fit into `ptrs`.
 * \param nptrs The total number of offsets in the file.
 */
OEXTRN void frmlst(FILE* fp, long* ptrs, int bufsz, int* nptrs);

/**
 * \brief Read the frame starting at the current position of the stream,
 * optionally extracting some meta information about the frame.
 *
 * The entire frame will be read and discared. To go back to the start of
 * the frame, fseek to the position given in `ptr`.
 *
 * \param fp The file stream which is currently pointing to the start
 * of a frame.
 * \param ptr If not NULL, will be populated with the current offset
 * into the stream.
 * \param ntot If not NULL, will be populated with the number of
 * particles contained in the frame.
 * \param nk If not NULL, will be populated with the number of
 * header items contained in the frame.
 *
 * \return True if the call was successful.
 */
OEXTRN bool frmmeta(FILE* fp, long* ptr, int* ntot, int* nk);

/**
 * \brief Determine the size, in bytes, of the frame data, not including
 * the metadata.
 *
 * \param ntot The total number of particles.
 * \param nk The total number of header vars.
 *
 * \return The size of the data record, in bytes.
 */
OEXTRN long frmsz(int ntot, int nk);


/**
 * \brief Read each particle in the frame.
 *
 * The `cb` callback will be passed the following, in order:
 * 1. The total number of particles in the frame.
 * 2. The total number of header values.
 * 3. A struct containing the header values.
 * 4. The particle's components, in the order: mass, position*3, velocity*3.
 * 7. Arbitrary data passed to this function.
 *
 * Return false from the callback to stop iterating.
 *
 * \param fp The file handle.
 * \param ptr The position in the stream containing the frame of interest.
 * Pass a negative number to use the current position of the file.
 * \param cb The function to call with particle data.
 * \param dat Any arbitrary data to pass to the callback.
 */
OEXTRN void itrprt(FILE* fp, long ptr, fn_frmrdr_t cb, void* dat);

#ifdef __cplusplus
}
#endif
#endif
