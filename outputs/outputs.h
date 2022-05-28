#ifndef OUTPUTS_HEADER_INCLUDED
#define OUTPUTS_HEADER_INCLUDED

#define OEXTRN extern

#include <stdio.h>
#include <stdbool.h>
#include "outputs/vmath.h"

#ifdef __cplusplus
extern "C" {
#endif


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

struct particle_t {
    double x[3];
    double dx[3];
    double m;
};

struct frm_t {
    int ntot;
    int nk;
    struct frm_hdr_t* hdr;
    struct particle_t ptcls[];
};

/**
 * \brief A callback function to read every particle in a frame.
 */
typedef bool (*fn_frmrdr_t)(int,
                            int,
                            const struct frm_hdr_t*,
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
 * \brief Read a frame header, discarding the particle data.
 *
 * \param fp The file handle.
 * \param ptr The position in the stream containing the frame of interest.
 * Pass a negative number to use the current position of the file.
 * \param ntot If not NULL, will be populated with the number of
 * particles contained in the frame.
 * \param nk If not NULL, will be populated with the number of
 * header items contained in the frame.
 * \param hdr Will contain the frame hader. Must be `free`d.
 *
 * \return True if the call was successful.
 */
OEXTRN void frmhdr(FILE* fp, long ptr, int* ntot, int* nk,
                   struct frm_hdr_t** hdr);

/**
 * \brief Determine the size, in bytes, of the frame data, not including
 * the metadata.
 *
 * \param ntot The total number of particles.
 * \param nk The total number of header vars.
 * \param kz19 Set to true if the output file was generated with KZ(19) on.
 *
 * \return The size of the data record, in bytes.
 */
OEXTRN long frmsz(int ntot, int nk, bool kz19);

/**
 * \brief Read one frame, including all of its data, from the stream
 * at the given location.
 *
 * To read from the stream at its current location, pass a negative integer
 * to `ptr`.
 *
 * \param fp The file pointer.
 * \param ptr The location of the frame to read, or a negative value to
 * read the next frame in the stream.
 * \param frm Will be populated with the frame data. Must be released
 * with a call to `freefrm`.
 */
OEXTRN struct frm_t* rdfrm(FILE* fp, long ptr);

/**
 * \brief Free the memory held by the frame structure.
 *
 * \param frm The frame to free. Can be NULL (no op).
 */
OEXTRN void freefrm(struct frm_t* frm);


/**
 * \brief Skip the current frame.
 *
 * \param fp A file open to the beginning of a frame.
 */
OEXTRN void skpfrm(FILE* fp);

/**
 * \brief Calculate the force and first time derivative from
 * the galaxy particle according to Irrgang et al. 2013.
 *
 * \param hdr A frame header.
 * \param fp A 3-vector to be populated with the force.
 * \param fdp A 3-vector to be populated with the first time derivative of
 * the force.
 */
OEXTRN void forceir13(const struct frm_hdr_t* hdr,
                      const double rg[3], const double vg[3],
                      double fp[3], double fd[3]);


/**
 * \brief Find the tidal radius from circular angular velocity.
 *
 * \param omega Circular angular velocity.
 * \return The tidal radius.
 */
OEXTRN double rtide(double omega);

#ifdef __cplusplus
}
#endif
#endif
