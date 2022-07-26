#include "outputs.h"

double clmass(const struct frm_t* frame) {
    double mass = 0;
    for (int i = 0; i < frame->ntot; i++) {
        mass += frame->ptcls[i].m;
    }
    return mass;
}
