/* modelsample : sampling random profile HMMs
 */
#ifndef h4MODELSAMPLE_INCLUDED
#define h4MODELSAMPLE_INCLUDED
#include "h4_config.h"

#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "h4_anchorset.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"


/* 1. Profile sampling routines */
extern int h4_modelsample             (ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);
extern int h4_modelsample_zeropeppered(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);


/* 2. Profiles with a readily enumerable number of possible paths/sequences */
extern int h4_modelsample_Enumerable  (ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);
extern int h4_modelsample_Enumerable2 (ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);


/* 3. Profiles with only a single valid path (or profile/seq pairs, or profile/seq/anchorset triplets) */
extern int h4_modelsample_SinglePath   (ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);
extern int h4_modelsample_SinglePathSeq(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, 
                                        H4_PROFILE **ret_hmm, ESL_SQ **ret_sq,
                                        H4_MODE **opt_mo, H4_PATH **opt_pi, H4_ANCHORSET **opt_anch, float *opt_sc);
extern int h4_modelsample_SinglePathASC(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, 
                                        H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                                        H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_sc);


/* 4. Profiles such that ASC paths = all valid paths */
extern int h4_modelsample_AnchoredUni  (ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, 
                                        H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                                        H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_tsc);
extern int h4_modelsample_AnchoredLocal(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M,
                                        H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                                        H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_tsc);
extern int h4_modelsample_AnchoredMulti(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M,
                                        H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                                        H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_tsc);

#endif /*h4MODELSAMPLE_INCLUDED*/
