#include <h4_config.h>


#include "easel.h"
#include "esl_sq.h"

#include "h4_filtermx.h"
#include "h4_mode.h"
#include "h4_profile.h"

/* The includes are deliberately outside the #ifdef.  If AVX512 is
 * enabled at compile-time, we provide an h4_vitfilter_avx512() function
 * that works; otherwise we provide one that issues a fatal error.
 * This allows driver programs (such as vitfilter_benchmark) to have
 * --sse/--avx/--avx512 options to call vitfilter_{sse,avx,avx512}
 * directly, overriding the CPU dispatcher.
 */
#ifdef eslENABLE_AVX512

int
h4_vitfilter_avx512(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc)
{
  return eslFAIL;  // SRE: TBW
}


#else // ! eslENABLE_AVX512
/* provide a callable function even when we're `./configure --disable-avx512` */
int
h4_vitfilter_avx512(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(hmm); ESL_UNUSED(mo); ESL_UNUSED(fx); ESL_UNUSED(ret_sc);  // shut up, compiler, I know what I'm doing
  esl_fatal("AVX512 support was not enabled at compile time. Can't use h4_vitfilter_avx512().");
  return eslFAIL; // NOTREACHED
}
#endif // eslENABLE_AVX512 or not
