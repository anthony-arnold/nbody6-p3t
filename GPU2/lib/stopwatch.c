double start = -1.0;

#include <sys/time.h>
void stopwatch_(double* dt){
    struct timeval tv;
    gettimeofday(&tv, 0);
    double t = tv.tv_sec + 1.e-6 * tv.tv_usec;
    if (start < 0) {
       start = t;
    }
    *dt = t - start;
}
