#include <h4_config.h>

#include "easel.h"
#include "esl_sq.h"

#include "h4_profile.h"

#ifdef eslENABLE_AVX512  // includes are deliberately outside the ifdef

int
h4_ssvfilter_avx512(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc)
{
  return eslFAIL;  // SRE: TBW
}


#else // !eslENABLE_AVX512:  provide a callable function even when we're `./configure --disable-avx512` 
int
h4_ssvfilter_avx512(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(hmm); ESL_UNUSED(ret_sc);  
  esl_fatal("AVX512 support was not enabled at compile time. Can't use h4_ssvfilter_avx512().");
  return eslFAIL; // NOTREACHED
}
#endif // eslENABLE_AVX512 or not
