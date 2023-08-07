#include <h4_config.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "h4_checkptmx.h"
#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_sparsemask.h"

#include "fbfilter.h"

#ifdef eslENABLE_AVX

int
h4_fwdfilter_avx(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc)
{
  return eslFAIL; 		// SRE: TBW
}

#else  // ! eslENABLE_AVX
/* provide callable functions even when we're `./configure --disable-avx` */
int
h4_fwdfilter_avx(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(hmm); ESL_UNUSED(mo); ESL_UNUSED(cpx); ESL_UNUSED(opt_sc); 
  esl_fatal("AVX support was not enabled at compile time. Can't use h4_fwdfilter_avx().");
  return eslFAIL; // NOTREACHED
}
int
h4_bckfilter_avx(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(hmm); ESL_UNUSED(mo); ESL_UNUSED(cpx); ESL_UNUSED(sm); ESL_UNUSED(sm_thresh); 
  esl_fatal("AVX support was not enabled at compile time. Can't use h4_bckfilter_avx().");
  return eslFAIL; // NOTREACHED
}
#endif // eslENABLE_AVX or not
