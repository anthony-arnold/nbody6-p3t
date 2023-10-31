#include "outputs.h"
#include <math.h>
#include "outputs/vmath.h"

static const double mb = 5.0e9;
static const double rhobm = -1.8e0;
static const double bulgecut=1900.e0;
static const double m2=6.8E10;
static const double a2=3000.e0;
static const double b2=280.e0;
static const double rhodm=7.1111E-3;
static const double rsdm=16000.0e0;

static double gammln(double xx) {
    static double coeff[6] = {
        76.18009172947146e0,
        -86.50532032941677e0,
        24.01409824083091e0,
        -1.231739572450155e0,
        0.1208650973866179e-2,
        -0.5395239384953e-5
    };
    static double stp = 2.5066282746310005e0;

    double ser = 1.000000000190015e0;
    double x = xx;
    double y = x;
    double tmp = x + 5.5;
    tmp = (x + 0.5) * log(tmp) - tmp;

    for (int j = 0; j < 5; j++) {
        y++;
        ser += coeff[j] / y;
    }
    return tmp + log(stp * ser/x);
}

static double gser(double a, double x) {
    static double eps = 3.0e-7;
    static int itmax = 100;

    double gln = gammln(a);
    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    for (int n = 1; n <= itmax; n++) {
        ap++;
        del = del * x / ap;
        sum += del;
        if (fabs(del) < fabs(sum)*eps) {
            break;
        }
    }
    return sum * exp(-x + a * log(x) - gln);
}

static double gcf(double a, double x) {
    static double fpmin = 1.0e-30;
    static double eps = 3.0e-7;
    static int itmax = 100;
    double gln = gammln(a);
    double b = x + 1 - a;
    double c = 1.0 / fpmin;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i <= itmax; i++) {
        double an = -i * (i - a);
        b += 2;
        d = an * d + b;
        if (fabs(d) < fpmin) {
            d = fpmin;
        }
        c = b + an / c;
        if (fabs(c) < fpmin) {
            c = fpmin;
        }
        d = 1.0 / d;
        double del = d * c;
        h = h * del;
        if (fabs(del - 1) < eps) {
            break;
        }
    }
    return exp(-x + a * log(x) - gln) * h;
}

static double gammp(double a, double x) {
    if (x < a + 1) {
        return gser(a, x);
    }
    else {
        return 1.0 - gcf(a, x);
    }
}

static void bulge(const double dis,
                  const double rv,
                  const double ru[3],
                  const double vu[3],
                  double f1[3],
                  double fd1[3])
{
    double dis3 = dis*dis*dis;
    double dis5 = dis3 * dis * dis;
    double mbulge = mb*gammp(1.5+rhobm/2.0,dis*dis/bulgecut/bulgecut);
    double xf = rv/dis3/dis*pow(dis/bulgecut,3.0-rhobm)*
        exp(-dis*dis/bulgecut/bulgecut);

    for (int i = 0; i < 3; i++) {
        f1[i] += -mbulge*ru[i]/dis3;
    }
    for (int i = 0; i < 3; i++) {
        fd1[i] += -mbulge*vu[i]/dis3+(3.e0*mbulge*rv/dis5-mbulge*xf)*ru[i];
    }
}

static void disc(double rv, const double ru[3], const double vu[3], double f2[3], double fd2[3]) {

    double xc1 =
        sqrt(ru[0]*ru[0] + ru[1]*ru[1] + pow(a2+sqrt(ru[2]*ru[2]+b2*b2), 2));
    double xc2 = a2+sqrt(ru[2]*ru[2]+b2*b2);
    double xc3 = sqrt(ru[2]*ru[2]+b2*b2);

    f2[0] += -m2*ru[0]/pow(xc1,3);
    f2[2] += -m2*ru[1]/pow(xc1,3);
    f2[3] += -m2*ru[2]/pow(xc1,3)/
        sqrt(ru[2]*ru[2]+b2*b2)*(a2+sqrt(ru[2]*ru[2]+b2*b2));

    fd2[0] += -m2*vu[0]/pow(xc1,3)+3*m2*ru[0]/pow(xc1,5)*(rv+
                                                 ru[2]*vu[2]*a2/xc3);
    fd2[1] += -m2*vu[1]/pow(xc1,3)+3*m2*ru[1]/pow(xc1,5)*(rv+
                                                 ru[2]*vu[2]*a2/xc3);
    fd2[2] += -m2*vu[2]/pow(xc1,3)+
        m2*xc2/xc3/pow(xc1,5)*(3*ru[2]*rv+ru[2]*ru[2]*vu[2]*a2/xc3)+
        m2*vu[2]/pow(xc1,3)*(ru[2]*ru[2]*a2/pow(xc3,3)-xc2/xc3);

}

static void halo(const double dis,
                 const double rv,
                 const double ru[3],
                 const double vu[3],
                 double f3[3],
                 double fd3[3])
{
    double dis3 = dis*dis*dis;
    double dis5 = dis3 * dis * dis;
    double mhalo  = 8.0*atan(1.0)*rhodm*pow(rsdm,3.0)*(log(1.0+dis/rsdm)-dis/(dis+rsdm));
    double mhalo2 = 8.0*atan(1.0)*rhodm*pow(rsdm,3.0)/pow(dis+rsdm, 2);

    for (int i = 0; i < 3; i++) {
        f3[i] +=-mhalo*ru[i]/dis3;
    }

    for (int i = 0; i < 3; i++) {
        fd3[i] += -mhalo*vu[i]/dis3+3.0*mhalo*rv*ru[i]/dis5
            -mhalo2*rv*ru[i]/dis3;
    }
}

void forcebv15(const struct frm_hdr_t* hdr,
               const double rg[3],
               const double vg[3],
               double fp[3],
               double fd[3])
{
    double ru[3];
    double vu[3];

    for (int i = 0; i < 3; i++) {
        ru[i] = rg[i] * hdr->rbar;
    }
    for (int i = 0; i < 3; i++) {
        vu[i] = vg[i] * hdr->vstar;
    }

    double dis = sqrt(mag2(ru));
    double rv = dot(ru, vu);

    for (int i = 0; i < 3; i++) {
        fp[i] = 0.0;
        fd[i] = 0.0;
    }
    bulge(dis, rv, ru, vu, fp, fd);
    disc(rv, ru, vu, fp, fd);
    halo(dis, rv, ru, vu, fp, fd);


    for (int i = 0; i < 3; i++) {
        fp[i] = fp[i] / hdr->zmbar * pow(hdr->rbar, 2);
    }

    for (int i = 0; i < 3; i++) {
        fd[i] = fd[i] / hdr->zmbar * pow(hdr->rbar, 2) / hdr->vstar;
    }
}
