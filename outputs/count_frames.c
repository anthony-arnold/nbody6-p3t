#include "outputs.h"
#include <stdio.h>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: count_frames <file>\n");
        return 1;
    }

    FILE* fp = opnout(argv[1]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

    int nfrm = numfrm(fp);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }
    printf("%i frames in %s\n", nfrm, argv[1]);

    fclose(fp);
    return 0;
}
