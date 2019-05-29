/* Emitting (sampling) sequence from a profile HMM
 */
#ifndef h4EMIT_INCLUDED
#define h4EMIT_INCLUDED
#include "h4_config.h"

#include "easel.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"

extern int h4_emit(ESL_RANDOMNESS *rng, const H4_PROFILE *hmm, const H4_MODE *mo, ESL_SQ *sq, H4_PATH *pi);


#endif // h4EMIT_INCLUDED
