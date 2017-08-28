#ifndef p7REFERENCE_ASC_FWDBACK_INCLUDED
#define p7REFERENCE_ASC_FWDBACK_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_anchors.h"

#include "dp_reference/p7_refmx.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif

extern int p7_ReferenceASCForward (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D, 
				   P7_REFMX *afu, P7_REFMX *afd, float *opt_sc);
extern int p7_ReferenceASCBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D, 
				   P7_REFMX *abu, P7_REFMX *abd, float *opt_sc);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7REFERENCE_ASC_FWDBACK_INCLUDED*/
