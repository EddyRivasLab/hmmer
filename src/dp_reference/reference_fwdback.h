#ifndef p7REFERENCE_FWDBACK_INCLUDED
#define p7REFERENCE_FWDBACK_INCLUDED

#include <p7_config.h>

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"

extern int p7_ReferenceForward (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, float *opt_sc);
extern int p7_ReferenceBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, float *opt_sc);

#endif /* p7REFERENCE_FWDBACK_INCLUDED */

