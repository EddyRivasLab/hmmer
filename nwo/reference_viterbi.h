/* Reference implementation of Viterbi scoring, alignment.
 */
#ifndef h4REFERENCE_VITERBI_INCLUDED
#define h4REFERENCE_VITERBI_INCLUDED

#include "easel.h"

#include "h4_profile.h"
#include "h4_mode.h"
#include "h4_refmx.h"
#include "h4_path.h"

extern int h4_ReferenceViterbi(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx, H4_PATH *opt_tr, float *opt_sc);

#endif /*h4REFERENCE_VITERBI_INCLUDED*/ 
