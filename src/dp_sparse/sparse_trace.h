#ifndef p7SPARSE_TRACE_INCLUDED
#define p7SPARSE_TRACE_INCLUDED

#include "p7_config.h"

#include "esl_random.h"

#include "base/p7_profile.h"
#include "base/p7_trace.h"
#include "dp_sparse/p7_sparsemx.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
extern int p7_sparse_trace_Viterbi      (                                      const P7_PROFILE *gm, const P7_SPARSEMX *sx, P7_TRACE *tr);
extern int p7_sparse_trace_Stochastic   (ESL_RANDOMNESS *rng, float **wrk_byp, const P7_PROFILE *gm, const P7_SPARSEMX *sx, P7_TRACE *tr);
extern int p7_sparse_trace_StochasticSeg(ESL_RANDOMNESS *rng, float **wrk_byp, const P7_PROFILE *gm, const P7_SPARSEMX *sx, float p0, float sF, int g, float *dpc, float *xc, P7_TRACE *tr);
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7SPARSE_TRACE_INCLUDED*/

