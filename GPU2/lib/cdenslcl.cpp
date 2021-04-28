// Get local density of core particle neighbour lists.
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>

namespace {

enum {
    NBMIN = 6
};

double denslcl(
    int i, double rs, double xi[3], double x[][3], double mass[], int* list)
{
    int nnb = list[0];
    if (nnb < NBMIN) {
        return 0.0;
    }

    // Find the six closest neighbours.
    double r2n[NBMIN];
    int in[NBMIN];

    for (int nb = 0; nb < nnb; nb++) {
        int j = list[nb + 1] - 1;

        double rx = xi[0]-x[j][0];
        double ry = xi[1]-x[j][1];
        double rz = xi[2]-x[j][2];

        double r2 = rx*rx + ry*ry + rz*rz;

        double *end = r2n + std::min(static_cast<int>(NBMIN), nb);
        double *pos = std::lower_bound(r2n, end, r2);

        int *npos = in + std::distance(r2n, pos);
        double *last;
        double *dest;
        int *nlast;
        int *ndest;

        if (nb < NBMIN) {
            last = end;
            dest = end + 1;

            nlast = in + nb + 1;
            ndest = nlast + 1;
        }
        else if (pos != end) {
            last = end - 1;
            dest = end;

            ndest = in + NBMIN;
            nlast = ndest - 1;
        }
        else {
            // Not in the closests neighbours.
            continue;
        }

        std::move_backward(pos, last, dest);
        *pos = r2;

        std::move_backward(npos, nlast, ndest);
        *npos = j;
    }

    double xmass = std::accumulate(in, in + NBMIN - 1, 0.0,
                                   [&mass](double tot, int n) {
                                       return tot + mass[n];
                                   });

    double r6_2 = r2n[NBMIN-1];
    double r6 = std::sqrt(r6_2);
    return xmass / (r6 * r6_2);
}

}

extern "C" void cdenslcl_(int* nc,
                          int jlist[],
                          double rho[],
                          double rs[],
                          double x[][3],
                          double mass[],
                          int *lmax,
                          int *list)
{
#pragma omp parallel for
    for (int i = 0; i < *nc; i++) {
        int idx = jlist[i] - 1;
        rho[idx] = denslcl(idx, rs[idx], x[idx], x, mass, &list[idx**lmax]);
    }
}
