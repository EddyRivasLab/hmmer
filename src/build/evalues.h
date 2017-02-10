#ifndef p7EVALUES_INCLUDED
#define p7EVALUES_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "base/p7_bg.h"
#include "base/p7_hmm.h"
#include "base/p7_profile.h"

#include "dp_vector/p7_oprofile.h"

#include "build/p7_builder.h"

extern int p7_Calibrate(P7_HMM *hmm, P7_BUILDER *cfg_b, ESL_RANDOMNESS **byp_rng, P7_BG **byp_bg, P7_PROFILE **byp_gm, P7_OPROFILE **byp_om);
extern int p7_Lambda(P7_HMM *hmm, P7_BG *bg, double *ret_lambda);
extern int p7_MSVMu     (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda,               double *ret_mmu);
extern int p7_ViterbiMu (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda,               double *ret_vmu);
extern int p7_Tau       (ESL_RANDOMNESS *r, P7_OPROFILE *om, P7_BG *bg, int L, int N, double lambda, double tailp, double *ret_tau);

#endif /*p7EVALUES_INCLUDED*/
