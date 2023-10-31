#include <string.h>
#include <ctype.h>
#include "outputs.h"
#include "oerrno.h"

typedef enum {
    OGM_NONE,
    OGM_IR13,
    OGM_BV15,

    // Save for last.
    OGM_MAX
} OGalModel;

static const char* galmodel_strs[OGM_MAX] = {
    "NONE",
    "IRRGANG",
    "BOVY"
};

int cmp_gal_model(const char* galmodel, const char* to) {
    const char* s, *t;
    for (s = to, t = galmodel; *s && *t; t++, s++) {
        if (toupper(*s) != toupper(*t)) {
            return 0;
        }
    }

    // both end on a nul
    return *s == *t;
}

OGalModel _get_gal_model(const char* galmodel) {
    for (int i = 0; i < OGM_MAX; i++) {
        if (cmp_gal_model(galmodel, galmodel_strs[i])) {
            return (OGalModel)i;
        }
    }
    return OGM_MAX;
}

void galforce(const struct frm_hdr_t* hdr,
              const char* galmodel,
              const double rg[3],
              const double vg[3],
              double fp[3],
              double fd[3])
{
    if (!galmodel) {
        _oseterrno(OERR_PARAM);
        return;
    }

    OGalModel model = _get_gal_model(galmodel);
    if (model == OGM_NONE) {
        memset(fp, 0, sizeof(double) * 3);
        memset(fd, 0, sizeof(double) * 3);
    }
    else if (model == OGM_IR13) {
        forceir13(hdr, rg, vg, fp, fd);
    }
    else if (model == OGM_BV15) {
        forcebv15(hdr, rg, vg, fp, fd);
    }
    else {
        _oseterrno(OERR_GALMODEL);
    }
}
