#ifndef OUTPUTS_HEADER_INCLUDED
#define OUTPUTS_HEADER_INCLUDED

#include <stdio.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#define OEXTRN extern

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


#ifdef __cplusplus
}
#endif
#endif
