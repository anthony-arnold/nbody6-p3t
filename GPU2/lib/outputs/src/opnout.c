#include <stdio.h>
#include "oerrno.h"

FILE* opnout(const char* path) {
    FILE* fp = fopen(path, "rb");
    if (fp == NULL) {
        _oseterrno(OERR_FILE);
    }
    return fp;
}
