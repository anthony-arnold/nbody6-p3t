/**
 *
 */

int* _oerrno() {
    static _Thread_local int oerrno = 0;
    return &oerrno;
}

int _ogeterrno() {
    return *_oerrno();
}

void _oseterrno(int e) {
    *_oerrno() = e;
}
