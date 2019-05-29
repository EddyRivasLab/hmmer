#ifndef h4REFERENCE_DP_INCLUDED
#define h4REFERENCE_DP_INCLUDED
#include "h4_config.h"

#include "easel.h"

#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_refmx.h"

extern int h4_reference_Viterbi (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx, H4_PATH *opt_pi, float *opt_sc);
extern int h4_reference_Forward (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx,                  float *opt_sc);
extern int h4_reference_Backward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx,                  float *opt_sc);

extern int h4_reference_Decoding(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_REFMX *fwd, const H4_REFMX *bck, H4_REFMX *pp);

extern int h4_reference_ViterbiTrace   (                                      const H4_PROFILE *hmm, const H4_MODE *mo, const H4_REFMX *rmx, H4_PATH *pi);
extern int h4_reference_StochasticTrace(ESL_RANDOMNESS *rng, float **wrk_byp, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_REFMX *rxf, H4_PATH *pi);


#endif /* h4REFERENCE_DP_INCLUDED */

