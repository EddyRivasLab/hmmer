#include "h4_config.h"
#ifdef eslENABLE_AVX512

#include "easel.h"
#include "esl_sq.h"

#include "h4_filtermx.h"
#include "h4_mode.h"
#include "h4_profile.h"


int
h4_vitfilter_avx512(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc)
{
  return eslFAIL;  // SRE: TBW
}


#else // ! eslENABLE_AVX512
/* Standard compiler-pleasing mantra for an #ifdef'd-out, empty code file. */
void h4_vitfilter_avx_silence_hack(void) { return; }
#if defined h4VITFILTER_AVX512_TESTDRIVE || h4VITFILTER_AVX512_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_AVX512 or not
