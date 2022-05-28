#include "outputs.h"

#define UNUSED(x) (void)(x)

bool list_frame(int ntot, int nk, struct frm_hdr_t* hdr, double* c, void* d) {
    UNUSED(nk);
    UNUSED(c);
    UNUSED(d);

    double t = hdr->t * hdr->tscale;
    printf("%lf %i\n", t, ntot);
    return false;
}

int main(int argc, char* argv[]) {
    int nframes;
    long ptrs[1024];

    if (argc < 2) {
        fprintf(stderr, "Usage: list_frames <file>\n");
        return 1;
    }

    FILE* fp = opnout(argv[1]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

    frmlst(fp, ptrs, 1024, &nframes);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }
    if (nframes > 1024) {
        fprintf(stderr, "Too many frames: %i", nframes);
        return 1;
    }

    for (int i = 0; i < nframes; i++) {
        itrprt(fp, ptrs[i], list_frame, NULL);
        if (ofail()) {
            fprintf(stderr, "%s\n", oerror());
            break;
        }
    }

    fclose(fp);
    return 0;
}
