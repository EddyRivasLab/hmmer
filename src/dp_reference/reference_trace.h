#ifndef p7REFERENCE_TRACE_INCLUDED
#define p7REFERENCE_TRACE_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "base/p7_profile.h"
#include "base/p7_trace.h"
#include "dp_reference/p7_refmx.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif

extern int p7_reference_trace_Viterbi   (                                      const P7_PROFILE *gm, const P7_REFMX *rmx, P7_TRACE *tr);
extern int p7_reference_trace_Stochastic(ESL_RANDOMNESS *rng, float **wrk_byp, const P7_PROFILE *gm, const P7_REFMX *rmx, P7_TRACE *tr);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif

#endif /*p7REFERENCE_TRACE_INCLUDED*/

