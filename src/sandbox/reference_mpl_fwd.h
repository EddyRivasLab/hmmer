#ifndef p7REFERENCE_MPL_FWD_INCLUDED
#define p7REFERENCE_MPL_FWD_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_coords2.h"

#include "dp_reference/p7_refmx.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
extern int p7_ReferenceMPLForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_COORD2 *dom, int ndom, P7_REFMX *rmx, float *opt_sc);
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7REFERENCE_MPL_FWD_INCLUDED*/
