#ifndef p7REFERENCE_DECODING_INCLUDED
#define p7REFERENCE_DECODING_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"

#include "dp_reference/p7_refmx.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif

	
extern int p7_ReferenceDecoding(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_REFMX *fwd, const P7_REFMX *bck, P7_REFMX *pp);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /* p7REFERENCE_DECODING_INCLUDED */

