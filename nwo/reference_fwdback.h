#ifndef h4REFERENCE_FWDBACK_INCLUDED
#define h4REFERENCE_FWDBACK_INCLUDED
#include "h4_config.h"

#include "easel.h"

#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_refmx.h"

extern int h4_ReferenceForward (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx, float *opt_sc);
extern int h4_ReferenceBackward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx, float *opt_sc);

#endif /* p7REFERENCE_FWDBACK_INCLUDED */

