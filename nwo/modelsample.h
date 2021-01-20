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


extern int h4_modelsample             (ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);
extern int h4_modelsample_zeropeppered(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);
extern int h4_modelsample_Enumerable  (ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);
extern int h4_modelsample_Enumerable2 (ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);

/* Profiles, profile/seq pairs, profile/seq/anchorset triplets with only a single path */
extern int h4_modelsample_SinglePath   (ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);
extern int h4_modelsample_SinglePathSeq(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, 
                                        H4_PROFILE **ret_hmm, ESL_SQ **ret_sq,
                                        H4_MODE **opt_mo, H4_PATH **opt_pi, H4_ANCHORSET **opt_anch, float *opt_sc);
extern int h4_modelsample_SinglePathASC(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, 
                                        H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                                        H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_sc);



#endif /*h4MODELSAMPLE_INCLUDED*/
