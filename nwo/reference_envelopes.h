#ifndef h4REFERENCE_ENVELOPES_INCLUDED
#define h4REFERENCE_ENVELOPES_INCLUDED
#include "h4_config.h"

#include "easel.h"
#include "h4_profile.h"
#include "h4_mode.h"
#include "h4_anchorset.h"
#include "h4_refmx.h"
#include "h4_envset.h"

extern int h4_reference_Envelopes(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_ANCHORSET *anch,
                                  const H4_REFMX *apu, const H4_REFMX *apd, H4_REFMX *afu, H4_REFMX *afd, H4_ENVSET *env);


#endif //h4REFERENCE_ENVELOPES_INCLUDED
