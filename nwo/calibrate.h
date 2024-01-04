#ifndef h4CALIBRATE_INCLUDED
#define h4CALIBRATE_INCLUDED
#include <h4_config.h>

#include "esl_random.h"

#include "h4_profile.h"

extern int h4_lambda(const H4_PROFILE *hmm, double *ret_lambda);
extern int h4_ssv_mu(ESL_RANDOMNESS *rng, const H4_PROFILE *hmm, int L, int N, double lambda, double *ret_smu);

#endif //h4CALIBRATE_INCLUDED  
