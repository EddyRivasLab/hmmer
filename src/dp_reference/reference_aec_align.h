#ifndef p7REFERENCE_AEC_ALIGN_INCLUDED
#define p7REFERENCE_AEC_ALIGN_INCLUDED

#include "p7_config.h"

#include "base/p7_profile.h"
#include "base/p7_envelopes.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif

extern int p7_reference_AEC_Align(const P7_PROFILE *gm, P7_ENVELOPES *env, const P7_REFMX *apu, const P7_REFMX *apd, P7_REFMX *mx, P7_TRACE *tr);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7REFERENCE_AEC_ALIGN_INCLUDED*/


