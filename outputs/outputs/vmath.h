#ifndef VMATH_HEADER_INCLUDED
#define VMATH_HEADER_INCLUDED

#ifndef OEXTRN
#define OEXTRN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

OEXTRN double dot(const double a[3], const double b[3]);
OEXTRN void cross(const double l[3], const double r[3], double out[3]);
OEXTRN double mag2(const double v[3]);
OEXTRN void vsub(const double l[3], const double r[3], double out[3]);
OEXTRN void vmul(const double l[3], const double f, double out[3]);


OEXTRN double dot_2d(const double a[2], const double b[2]);
OEXTRN double mag2_2d(const double v[2]);
OEXTRN void vsub_2d(const double l[2], const double r[2], double out[2]);

#ifdef __cplusplus
}
#endif
#endif
