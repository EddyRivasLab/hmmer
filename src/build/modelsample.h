#ifndef P7_HMM_SAMPLE_INCLUDED
#define P7_HMM_SAMPLE_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"

#include "base/p7_bg.h"
#include "base/p7_coords2.h"
#include "base/p7_hmm.h"
#include "base/p7_prior.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

extern int     p7_hmm_Sample           (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,                      P7_HMM **ret_hmm);
extern int     p7_hmm_SamplePrior      (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, const P7_PRIOR *pri, P7_HMM **ret_hmm);
extern int     p7_hmm_SampleUngapped   (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,                      P7_HMM **ret_hmm);
extern int     p7_hmm_SampleEnumerable (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,                      P7_HMM **ret_hmm);
extern int     p7_hmm_SampleEnumerable2(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc,                      P7_HMM **ret_hmm);
extern int     p7_hmm_SampleUniform    (ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, 
					float tmi, float tii, float tmd, float tdd,  P7_HMM **ret_hmm);
extern int     p7_hmm_SampleSinglePathed(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
extern int     p7_hmm_SampleSinglePathedSeq(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
					    P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
					    P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc);
#endif /*P7_HMM_SAMPLE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

