#include "h4_config.h"

#include "easel.h"
#include "esl_sq.h"

#include "h4_profile.h"

#ifdef eslENABLE_NEON  // includes are deliberately outside the ifdef

int
h4_ssvfilter_neon(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc)
{
  return eslFAIL;  // SRE: TBW
}


#else // !eslENABLE_NEON:  provide a callable function even when we're `./configure --disable-neon` 
int
h4_ssvfilter_neon(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(hmm); ESL_UNUSED(ret_sc);  
  esl_fatal("ARM NEON support was not enabled at compile time. Can't use h4_ssvfilter_neon().");
  return eslFAIL; // NOTREACHED
}
#endif // eslENABLE_NEON or not
