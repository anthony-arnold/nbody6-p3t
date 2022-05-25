/**
 *
 */

#ifdef THREAD_SAFE
#define LOCAL _Thread_local
#else
#define LOCAL
#endif

int* _oerrno() {
    static LOCAL int oerrno = 0;
    return &oerrno;
}

int _ogeterrno() {
    return *_oerrno();
}

void _oseterrno(int e) {
    *_oerrno() = e;
}
