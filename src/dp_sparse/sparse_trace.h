#ifndef P7_SPARSE_TRACE_INCLUDED
#define P7_SPARSE_TRACE_INCLUDED

#include "p7_config.h"

#include "esl_random.h"

#include "base/p7_profile.h"
#include "base/p7_trace.h"
#include "dp_sparse/p7_sparsemx.h"

extern int p7_sparse_trace_Viterbi      (                                      const P7_PROFILE *gm, const P7_SPARSEMX *sx, P7_TRACE *tr);
extern int p7_sparse_trace_Stochastic   (ESL_RANDOMNESS *rng, float **wrk_byp, const P7_PROFILE *gm, const P7_SPARSEMX *sx, P7_TRACE *tr);
extern int p7_sparse_trace_StochasticSeg(ESL_RANDOMNESS *rng, float **wrk_byp, const P7_PROFILE *gm, const P7_SPARSEMX *sx, float p0, float sF, int g, float *dpc, float *xc, P7_TRACE *tr);

#endif /*P7_SPARSE_TRACE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/ 
