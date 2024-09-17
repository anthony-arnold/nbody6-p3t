#include <stdlib.h>
#include "outputs.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: list_frames <file>\n");
        return 1;
    }

    FILE* fp = opnout(argv[1]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

    struct frm_hdr_t* hdr = NULL;
    while (!feof(fp)) {
        int ntot;
        frmhdr(fp, -1, &ntot, NULL, &hdr);
        if (ofail()) {
            fprintf(stderr, "%s\n", oerror());
            return 1;
        }
        if (hdr) {
            double t = hdr->t * hdr->tscale;
            printf("%lf %i\n", t, ntot);
        }
    }
    if (hdr) {
        free(hdr);
    }

    fclose(fp);
    return 0;
}
