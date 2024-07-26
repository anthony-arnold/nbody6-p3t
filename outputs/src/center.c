#include <math.h>
#include "outputs.h"
#include "oerrno.h"

#define PI 3.1415927
#define GND  (90-27.128251)/180.0*PI
#define GNR  (90+192.859481)/180.0*PI
#define L0   32.931918/180.0*PI
#define GCDIST 8178.0    // Gravity collaboration

void center(const struct frm_t* frame, double* pra, double* pdec) {
    /* Clear error flag */
    _oseterrno(OERR_SUCCESS);

    if (!frame || !pra || !pdec) {
        _oseterrno(OERR_NULL);
        return;
    }

    double x = frame->hdr->rdens[0];
    double y = frame->hdr->rdens[1];
    double z = frame->hdr->rdens[2];

    double xst = -(x-GCDIST)/1000.0;
    double yst = y/1000.0;
    double zst = z/1000.0;

    double dis  = sqrt(xst*xst+yst*yst+zst*zst);
    double b = asin(zst/dis);
    double l = acos(xst/dis/cos(b));
    double l2= asin(yst/dis/cos(b));

    if (l2 < 0.0) {
        l=2.0 * PI - l;
    }

    double dec = asin(cos(b)*sin(l-L0)*sin(GND)+sin(b)*cos(GND));
    double x1 = cos(b)*cos(l-L0);
    double x2 = sin(b)*sin(GND)-cos(b)*cos(GND)*sin(l-L0);

    double ra;
    if (x2>=0.0) {
        ra = (GNR-PI/2.0)+atan(x1/x2);
    }
    else if (x1>=0.0) {
        ra = (GNR-PI/2.0)+atan(x1/x2)+PI;
    }
    else {
        ra = (GNR-PI/2.0)+atan(x1/x2)-PI;
    }

    if (ra>2.0*PI) {
        ra -= 2.0*PI;
    }

    ra  *= 180.0/PI;
    dec *= 180.0/PI;

    *pra = ra;
    *pdec = dec;
}
