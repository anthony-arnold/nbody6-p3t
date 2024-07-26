#include <math.h>
#include <stdlib.h>
#include "outputs.h"
#include "oerrno.h"

#define PI 3.1415927
#define GND  (90-27.128251)/180.0*PI
#define GNR  (90+192.859481)/180.0*PI
#define L0   32.931918/180.0*PI
#define GCDIST 8178.0    // Gravity collaboration

void get_radecl(const struct particle_t* ptcl, double* pra, double* pdec) {
    double xst = -(ptcl->x[0] - GCDIST) / 1000.0;
    double yst = ptcl->x[1] / 1000.0;
    double zst = ptcl->x[2] / 1000.0;

    double dis  = sqrt(xst*xst+yst*yst+zst*zst);
    double b = asin(zst/dis);
    double l = acos(xst/dis/cos(b));
    double l2= asin(yst/dis/cos(b));

    if (l2<0.0) {
        l = 2.0 * PI - l;
    }

    double decst = asin(cos(b)*sin(l-L0)*sin(GND)+sin(b)*cos(GND));
    double x1 = cos(b)*cos(l-L0);
    double x2 = sin(b)*sin(GND)-cos(b)*cos(GND)*sin(l-L0);

    double rast;
    if (x2>=0.0) {
        rast = (GNR-PI/2.0)+atan(x1/x2);
    }
    else if (x1>=0.0) {
        rast = (GNR-PI/2.0)+atan(x1/x2)+PI;
    }
    else {
        rast = (GNR-PI/2.0)+atan(x1/x2)-PI;
    }

    if (rast>2.0*PI) {
        rast=rast-2.0*PI;
    }

    rast  *= 180.0/PI;
    decst *= 180.0/PI;

    *pra = rast;
    *pdec = decst;
}

void center(const struct frm_t* frame, double* pra, double* pdec) {
    /* Clear error flag */
    _oseterrno(OERR_SUCCESS);

    if (!frame || !pra || !pdec) {
        _oseterrno(OERR_NULL);
        return;
    }

    int nstar = frame->ntot;
    double* rast  = malloc(nstar*sizeof(double));
    double* decst = malloc(nstar*sizeof(double));
    double* disij = malloc(nstar*sizeof(double));

    if (!rast || !decst || !disij) {
        _oseterrno(OERR_OOM);
        goto cleanup;
    }

    for (int i = 0; i < nstar; i++) {
        get_radecl(&frame->ptcls[i], &rast[i], &decst[i]);
    }

    double ramean  = 0.0;
    double decmean = 0.0;
    double weight  = 0.0;
    int samples = frame->ntot / 10;

    for (int i = 0; i < samples; i++) {
        for (int j = 0; j  < nstar; j++) {
            double di[2] = { rast[i], decst[i] };
            double dj[2] = { rast[j], decst[j] };
            double dij[2];
            vsub_2d(di, dj, dij);
            disij[j] = sqrt(mag2_2d(dij));
        }
        disij[i] = 1.E10;

        double mind;
        int merk;
        for (int k = 0; k < 10; k++) {
            mind = disij[0];
            merk = 0;
            for (int j = 0; j < nstar; j++) {
                if (disij[j] < mind) {
                    mind = disij[j];
                    merk = j;
                }
            }
            disij[merk] = 1E10;
        }

        double dens = 1.0/(mind*mind*mind);

        ramean  += rast[i]*dens;
        decmean += decst[i]*dens;
        weight  += dens;
    }

    ramean /= weight;
    decmean /= weight;

    *pra = ramean;
    *pdec = decmean;

 cleanup:
    if (rast) {
        free(rast);
    }
    if (decst) {
        free(decst);
    }
    if (disij) {
        free(disij);
    }
}
