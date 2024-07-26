/**
 * Print some details about the run.
 */


#include "outputs.h"
#include <stdlib.h>

double mass(struct frm_t* frame) {
    double tot = 0.0;
    for (int i = 0; i < frame->ntot; i++) {
        tot += frame->ptcls[i].m;
    }
    return tot * frame->hdr->zmbar;
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: info <file> <galmodel>\n");
        return 1;
    }

    FILE* fp = opnout(argv[1]);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

    struct frm_t* frame = rdfrm(fp, -1);
    if (ofail()) {
        fprintf(stderr, "%s\n", oerror());
        return 1;
    }

    // Get the initial time, mass and half-mass radius
    int n = frame->ntot;
    double m_init = mass(frame);
    double r_init = frame->hdr->rscale * frame->hdr->rbar;
    freefrm(frame);

    // Find the final frame
    double m_final = m_init;
    double r_final = r_init;
    double t_final = 0.0;
    double m_bound = m_init;
    double rafin = 0.0;
    double decfin = 0.0;
    while (!feof(fp)) {
        frame = rdfrm(fp, -1);
        if (frame) {
            m_final = mass(frame);
            bound(frame, argv[2], &m_bound, NULL, NULL);
            r_final = frame->hdr->rscale * frame->hdr->rbar;
            t_final = frame->hdr->tphys;
            center(frame, &rafin, &decfin);
            freefrm(frame);
        }
    }
    fclose(fp);

    // Print the info
    printf("%6i %d %0.4lf %0.4lf %0.4lf %0.4lf %0.4lf %0.4lf %0.6lf %0.6lf\n",
           n,
           m_init == m_final,
           m_init,
           r_init,
           m_bound,
           r_final,
           1.0 - m_final / m_init,
           t_final,
           rafin,
           decfin);
}
