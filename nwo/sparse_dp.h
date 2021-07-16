#ifndef h4SPARSE_DP_INCLUDED
#define h4SPARSE_DP_INCLUDED
#include "h4_config.h"

#include "easel.h"

#include "h4_path.h"
#include "h4_profile.h"
#include "h4_sparsemx.h"

extern int h4_sparse_Viterbi(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMASK *sm, H4_SPARSEMX *sx, H4_PATH *opt_pi, float *opt_vsc);
extern int h4_sparse_Forward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMASK *sm, H4_SPARSEMX *sx,                  float *opt_fsc);


#endif //h4SPARSE_DP_INCLUDED
