/**
 * A sample file for doing something silly - creating an audio
 * rendering of the cluster evolution.
 */
#include "outputs.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

const double TWO_PI = 6.28318530718;
const uint16_t CHANNELS = 2;
const uint16_t BITS_SAMPLE = 16;
const uint16_t BLOCK_ALIGN = CHANNELS * (BITS_SAMPLE / 8);
const uint32_t FRAMES_S = 24;
const uint32_t SAMPLES_S = FRAMES_S * (48000 / FRAMES_S);
const uint32_t SAMPLES_FRAME = SAMPLES_S / FRAMES_S;
const uint32_t BYTES_S = SAMPLES_S * BLOCK_ALIGN;
const uint32_t BYTES_FRAME = BLOCK_ALIGN * SAMPLES_FRAME;
const double SAMPLE_DT = 1.0 / SAMPLES_S;

const double MAX_FREQ = 4434.922; // above high C
const double MIN_FREQ = 27.5; // close to low A
const int16_t MAX_AMP = 32767;
const int16_t MIN_AMP = -32768;
const uint32_t FREQ_BINS = 128;

const union {
    uint32_t _TEST_ENDIAN;
    char IS_BIG_ENDIAN[4];
} _CHECK_BE = {0x01020304};

#define IS_BE _CHECK_BE.IS_BIG_ENDIAN[0] == 1
#define BYTE_SWAP_32(i) IS_BE ? (         \
            (((i) & 0xff) << 24)        | \
            ((((i) >> 8) & 0xff) << 16) | \
            ((((i) >> 16)& 0xff) << 8)  | \
            (((i) >> 24) & 0xff)) : (i);

#define BYTE_SWAP_16(i) IS_BE ? ( \
            (((i) & 0xff) << 8) | \
            (((i) >> 8) & 0xff)) : (i);

void writeu32(FILE* fp, uint32_t i) {
    uint32_t be = BYTE_SWAP_32(i);
    fwrite(&be, 4, 1, fp);
}
void writeu16(FILE* fp, uint16_t i) {
    uint16_t be = BYTE_SWAP_16(i);
    fwrite(&be, 4, 1, fp);
}
void write32(FILE* fp, int32_t i) {
    uint32_t be = BYTE_SWAP_32(i);
    fwrite(&be, 4, 1, fp);
}
void write16(FILE* fp, int16_t i) {
    uint16_t be = BYTE_SWAP_16(i);
    fwrite(&be, 4, 1, fp);
}

void write_header(FILE* fp) {
    fwrite("WAVEfmt ", 1, 8, fp);
    write32(fp, 16); // chunk size
    write16(fp, 1); // no compression
    writeu16(fp, CHANNELS); // stereo
    writeu32(fp, SAMPLES_S); // 44100 Hz
    writeu32(fp, BYTES_S); // bytes per second
    writeu16(fp, BLOCK_ALIGN);
    writeu16(fp, BITS_SAMPLE);

    // Write the data chunk ID and size
    // we don't actually know what the size will be, so
    // just write 0 with the intention of overwriting it
    // later.
    fwrite("data", 1, 4, fp);
    write32(fp, 0);
}

FILE* open_wav(const char* name) {
    char* fname = malloc(strlen(name) + 5);
    if (!fname) {
        fprintf(stderr, "OOM\n");
        exit(1);
    }

    sprintf(fname, "%s.wav", name);

    FILE* fp = fopen(fname, "w+");
    if (!fp) {
        fprintf(stderr, "Could not open output file %s\n", fname);
    }
    else {
        write_header(fp);
    }
    free(fname);
    return fp;
}

double omega(double r[3], double v[3]) {
    double r2 = mag2(r);
    double ome[3];
    cross(r, v, ome);
    for(int i = 0; i < 3; i++) {
        ome[i] /= r2;
    }
    return mag2(ome);
}

double freq(struct particle_t* ptcl) {
    return omega(ptcl->x, ptcl->dx) / TWO_PI;
}

int16_t clamp(int16_t imin, int16_t imax, int32_t i) {
    if (i > imax) {
        return imax;
    }
    if (i < imin) {
        return imin;
    }
    return i;
}
struct stream_t {
    struct frm_t* frame;
    FILE* fp;
    int16_t* buffer;
    double* left_domain;
    double* right_domain;
    double mmass;
    double mfreq;
};

void harmonics(double f, double a1, double a2, int16_t* buffer) {
    const double a[] = { a1, a2 };
#pragma omp parallel for
    for (uint32_t i = 0; i < SAMPLES_FRAME; i++) {
        double b = sin(SAMPLE_DT * i * f);
        for(int c = 0; c < 2; c++) {
            int32_t clip = buffer[i*2+c] + a[c] + b;
            buffer[i*2+c] = clamp(MIN_AMP, MAX_AMP, clip);
        }
    }
}

double binned_freq(uint32_t bin) {
    return MIN_FREQ + (MAX_FREQ - MIN_FREQ) / FREQ_BINS * bin;
}

void sample(struct stream_t* stream) {
    memset(stream->buffer, 0, BYTES_FRAME);
    /* Add up each harmonic */
    for (uint32_t i = 0; i < FREQ_BINS; i++) {
        harmonics(binned_freq(i),
                  stream->left_domain[i],
                  stream->right_domain[i],
                  stream->buffer);
    }
}

double max_mass(struct frm_t* frame) {
    double maxm = 0;
    for (int i = 0; i < frame->ntot; i++) {
        double m = frame->ptcls[i].m;
        if (m > maxm) {
            maxm = m;
        }
    }
    return maxm;
}
double max_freq(struct frm_t* frame) {
    double maxf = 0;
    for (int i = 0; i < frame->ntot; i++) {
        double f = freq(&frame->ptcls[i]);
        if (f > maxf) {
            maxf = f;
        }
    }
    return maxf;
}

void teardown_stream(struct stream_t* stream) {
    if (stream->left_domain)
        free(stream->left_domain);
    if (stream->right_domain)
        free(stream->right_domain);
    if (stream->buffer)
        free(stream->buffer);
    if (stream->frame)
        freefrm(stream->frame);
}

void setup_stream(struct stream_t* stream) {
    stream->left_domain = malloc(FREQ_BINS * sizeof(double));
    stream->right_domain = malloc(FREQ_BINS * sizeof(double));
    stream->buffer = malloc(BYTES_FRAME);
    if (!stream->buffer || !stream->left_domain || !stream->right_domain) {
        teardown_stream(stream);
        fprintf(stderr, "OOM\n");
        exit(1);
    }

    stream->mmass = max_mass(stream->frame);
    stream->mfreq = max_freq(stream->frame);
}

void fade(double a, struct particle_t* ptcl, double* left, double* right) {
    // Disipate left or right amplitude depending on the position.
    double posx = ptcl->x[0];

    double la = a, ra = a;
    if (posx > 0) {
        la = -a / log(posx);
    }
    else if (posx < 0) {
        ra = -a / log(fabs(posx));
    }

    *left += la;
    *right += ra;
}

void stream_frame(struct stream_t* stream) {
    memset(stream->buffer, 0, BYTES_FRAME);
    memset(stream->left_domain, 0, FREQ_BINS * sizeof(double));
    memset(stream->right_domain, 0, FREQ_BINS * sizeof(double));

    for (int i = 0; i < stream->frame->ntot; i++) {
        struct particle_t* p = &stream->frame->ptcls[i];
        // normalised amplitude and frequency
        double a = p->m / stream->mmass;
        double f = freq(p) / stream->mfreq;

        // binned frequency
        uint32_t bin = f * (FREQ_BINS - 1);

        fade(a, p, &stream->left_domain[bin], &stream->right_domain[bin]);
    }

    // Sum harmonics
    sample(stream);

    // write the frame
    fwrite(stream->buffer, 1, BYTES_FRAME, stream->fp);
}


void finalise_wav(FILE* fp, int32_t data_size) {
    rewind(fp);
    // skip the 8 bytes of header
    fseek(fp, 8, SEEK_CUR);

    // read the chunk size
    int32_t size;
    fread(&size, 4, 1, fp);
    size = BYTE_SWAP_32(size);

    // skip to the end of the header
    // plus the data tag
    fseek(fp, size + 4, SEEK_CUR);

    write32(fp, data_size);
    fclose(fp);
}

void make_wav(FILE* fin, const char* name) {
    struct stream_t stream = {0};
    stream.fp = open_wav(name);
    if (!stream.fp) {
        exit(1);
    }

    int32_t data_size = 0;
    while(!feof(fin)) {
       stream.frame = rdfrm(fin, -1);
       if (ofail()) {
           teardown_stream(&stream);
           fprintf(stderr, "%s\n", oerror());
           exit(1);
       }
       if (stream.frame == NULL) {
           break;
       }
       if (stream.buffer == NULL) {
           setup_stream(&stream);
       }
       stream_frame(&stream);
       data_size += BYTES_FRAME;

       freefrm(stream.frame);
       stream.frame = NULL;

       printf("%d\n", data_size);
    }

    teardown_stream(&stream);
    finalise_wav(stream.fp, data_size);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: wav <output>\n");
        return 1;
    }
    FILE* fp = opnout(argv[1]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

    make_wav(fp, argv[1]);
    fclose(fp);
    return 0;
}
