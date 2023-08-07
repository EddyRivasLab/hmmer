#ifndef h4ASCMX_INCLUDED
#define h4ASCMX_INCLUDED
#include <h4_config.h>

#include "h4_anchorset.h"
#include "h4_refmx.h"

extern int h4_ascmx_fb_compare_std (const H4_REFMX *std, const H4_REFMX *ascu, const H4_REFMX *ascd, const H4_ANCHORSET *anch, float abstol);
extern int h4_ascmx_pp_compare_std (const H4_REFMX *rxd, const H4_REFMX *apu,  const H4_REFMX *apd,  const H4_ANCHORSET *anch, float abstol);
extern int h4_ascmx_pp_compare_path(const H4_PATH  *pi,  const H4_REFMX *apu,  const H4_REFMX *apd,  const H4_ANCHORSET *anch, float abstol);
extern int h4_ascmx_compare_asc    (const H4_REFMX *au1, const H4_REFMX *ad1,  const H4_REFMX *au2,  const H4_REFMX *ad2, const H4_ANCHORSET *anch, float abstol);

extern int h4_ascmx_maxdiff_std    (const H4_REFMX *std, const H4_REFMX *ascu, const H4_REFMX *ascd, const H4_ANCHORSET *anch, float *ret_maxdiff);
extern int h4_ascmx_pp_maxdiff_path(const H4_PATH *pi,   const H4_REFMX *apu,  const H4_REFMX *apd,  const H4_ANCHORSET *anch, float *ret_maxdiff);
extern int h4_ascmx_pp_max_overage (const H4_REFMX *apu, const H4_REFMX *apd,  const H4_ANCHORSET *anch, float *ret_max_overage);

extern int h4_ascmx_pp_Validate    (const H4_REFMX *apu, const H4_REFMX *apd, const H4_ANCHORSET *anch, float abstol, char *errbuf);

#endif // h4ASCMX_INCLUDED
