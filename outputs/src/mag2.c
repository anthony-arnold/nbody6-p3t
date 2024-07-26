#include "outputs/vmath.h"

double mag2(const double v[3]) {
    return dot(v, v);
}

double mag2_2d(const double v[2]) {
    return dot_2d(v, v);
}
