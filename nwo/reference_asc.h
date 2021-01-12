#ifndef h4REFERENCE_ASC_INCLUDED
#define h4REFERENCE_ASC_INCLUDED

#include "h4_config.h"

#include "easel.h"
#include "h4_profile.h"
#include "h4_mode.h"
#include "h4_anchorset.h"
#include "h4_refmx.h"

extern int h4_reference_asc_Forward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo,
                                    const H4_ANCHORSET *anch, H4_REFMX *mxu, H4_REFMX *mxd, float *opt_sc);


#endif // h4REFERENCE_ASC_INCLUDED
