/* modelsample : sampling random profile HMMs
 */
#ifndef h4MODELSAMPLE_INCLUDED
#define h4MODELSAMPLE_INCLUDED
#include "h4_config.h"

#include "esl_alphabet.h"
#include "esl_random.h"

#include "h4_profile.h"

extern int h4_modelsample(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm);

#endif /*h4MODELSAMPLE_INCLUDED*/
