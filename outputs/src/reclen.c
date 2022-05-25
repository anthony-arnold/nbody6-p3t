#include <stdio.h>

int _reclen(FILE* fp) {
    int len;
    size_t size = fread(&len, sizeof(len), 1, fp);

    if (size != 1) {
        return 0;
    }
    return len;
}
