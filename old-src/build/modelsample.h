#ifndef p7MODELSAMPLE_INCLUDED
#define p7MODELSAMPLE_INCLUDED

#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"

#include "base/p7_bg.h"
#include "base/p7_anchors.h"
#include "base/p7_hmm.h"
#include "base/p7_prior.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

/* 1. Model sampling routines */
extern int     p7_modelsample            (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,                      P7_HMM **ret_hmm);
extern int     p7_modelsample_Prior      (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, const P7_PRIOR *pri, P7_HMM **ret_hmm);
extern int     p7_modelsample_Ungapped   (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,                      P7_HMM **ret_hmm);
extern int     p7_modelsample_Uniform    (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, 
					  float tmi, float tii, float tmd, float tdd,  P7_HMM **ret_hmm);

/* 2. Models with finite enumerable path # */
extern int     p7_modelsample_Enumerable (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,                      P7_HMM **ret_hmm);
extern int     p7_modelsample_Enumerable2(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,                      P7_HMM **ret_hmm);

/* 3. Models with only a single valid path */
extern int     p7_modelsample_SinglePathed(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
extern int     p7_modelsample_SinglePathedSeq(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
					      P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
					      P7_TRACE **opt_tr, P7_ANCHOR **opt_anch, int *opt_D, float *opt_sc);
extern int     p7_modelsample_SinglePathedASC(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
					      P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
					      P7_TRACE **opt_tr, P7_ANCHOR **opt_anch, int *opt_D, float *opt_sc);

/* 4. Models such that ASC paths == all paths */
extern int     p7_modelsample_AnchoredUni(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
					  P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
					  P7_TRACE **opt_tr, P7_ANCHOR **opt_anch, int *opt_D, float *opt_sc);
extern int     p7_modelsample_AnchoredLocal(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
					    P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
					    P7_TRACE **opt_tr, P7_ANCHOR **opt_anch, int *opt_D, float *opt_sc);
extern int     p7_modelsample_AnchoredMulti(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
					    P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
					    P7_TRACE **opt_tr, P7_ANCHOR **opt_anch, int *opt_D, float *opt_sc);
#endif /*p7MODELSAMPLE_INCLUDED*/

