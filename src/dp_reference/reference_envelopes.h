
#ifndef p7REFERENCE_ENVELOPES_INCLUDED
#define p7REFERENCE_ENVELOPES_INCLUDED
#include "p7_config.h"

#include "easel.h"

#include "base/p7_envelopes.h" 
#include "base/p7_anchors.h"
#include "base/p7_profile.h"

#include "dp_reference/p7_refmx.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
	
extern int p7_reference_Envelopes(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D,
				  const P7_REFMX *apu, const P7_REFMX *apd,
				  P7_REFMX *afu, P7_REFMX *afd, P7_ENVELOPES *env);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7REFERENCE_ENVELOPES_INCLUDED*/
