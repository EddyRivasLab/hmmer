#ifndef h4REFERENCE_MPAS_INCLUDED
#define h4REFERENCE_MPAS_INCLUDED

#include <h4_config.h>

#include "easel.h"
#include "esl_random.h"

#include "h4_anchorhash.h"
#include "h4_anchorset.h"
#include "h4_mode.h"
#include "h4_mpas.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_refmx.h"

extern int h4_reference_MPAS(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo,
                             const H4_REFMX *rxf, const H4_REFMX *rxd, H4_PATH *pi,
                             float **byp_wrk,  H4_ANCHORHASH *ah,
                             H4_REFMX *afu, H4_REFMX *afd, H4_ANCHORSET *anch, float *ret_asc,
                             const H4_MPAS_PARAMS *prm, H4_MPAS_STATS *stats);
extern int  h4_reference_mpas_path2anchors(const H4_PATH *pi, const H4_REFMX *rxd, H4_ANCHORSET *anch);


#endif //h4REFERENCE_MPAS_INCLUDED
