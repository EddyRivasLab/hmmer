#ifndef h4SPARSE_DP_INCLUDED
#define h4SPARSE_DP_INCLUDED
#include <h4_config.h>

#include "easel.h"
#include "esl_random.h"

#include "h4_path.h"
#include "h4_profile.h"
#include "h4_sparsemx.h"

extern int h4_sparse_Viterbi (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMASK *sm, H4_SPARSEMX *sx, H4_PATH *opt_pi, float *opt_vsc);
extern int h4_sparse_Forward (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMASK *sm, H4_SPARSEMX *sx,                  float *opt_fsc);
extern int h4_sparse_Backward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMASK *sm, H4_SPARSEMX *sx,                  float *opt_bsc);
extern int h4_sparse_Decoding(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMX *sxf, const H4_SPARSEMX *sxb, H4_SPARSEMX *sxd);
extern int h4_sparse_ViterbiTrace(const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMX *sx, H4_PATH *pi);
extern int h4_sparse_StochasticTrace(ESL_RANDOMNESS *rng, float **wrk_byp, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMX *sxf, H4_PATH *pi);

#endif //h4SPARSE_DP_INCLUDED
