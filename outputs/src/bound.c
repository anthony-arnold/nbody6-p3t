#include "outputs.h"
#include "oerrno.h"

void bound(const struct frm_t* frame,
           const char* galmodel,
           double* boundm,
           double* boundn,
           double* rt)
{
    if (!frame || !boundm || !boundn || !rt) {
        _oseterrno(OERR_NULL);
        return;
    }

    double clm = clmass(frame);
    *rt = rtide(frame->hdr, galmodel, clm);
    double rt2 = (*rt) * (*rt);

    double m = 0.0;
    int n = 0;
    for (int i = 0; i < frame->ntot; i++) {
        double x[3];
        vsub(frame->ptcls[i].x, frame->hdr->rdens, x);

        if (mag2(x) < rt2) {
            m += frame->ptcls[i].m;
            n++;
        }
    }

    *boundm = m / clm;
    *boundn = (double)n / (double)frame->ntot;
}
