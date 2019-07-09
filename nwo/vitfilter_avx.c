#include "h4_config.h"
#ifdef eslENABLE_AVX

#include "easel.h"
#include "esl_sq.h"

#include "h4_filtermx.h"
#include "h4_mode.h"
#include "h4_profile.h"

int
h4_vitfilter_avx(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc)
{
  return eslFAIL;  // SRE: TBW
}


#else // ! eslENABLE_AVX
/* Standard compiler-pleasing mantra for an #ifdef'd-out, empty code file. */
void h4_vitfilter_avx_silence_hack(void) { return; }
#if defined h4VITFILTER_AVX_TESTDRIVE || h4VITFILTER_AVX_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_AVX or not
