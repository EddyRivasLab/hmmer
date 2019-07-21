#include "h4_config.h"


#include "easel.h"
#include "esl_sq.h"

#include "h4_filtermx.h"
#include "h4_mode.h"
#include "h4_profile.h"

#ifdef eslENABLE_VMX

int
h4_vitfilter_vmx(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc)
{
  return eslFAIL;  // SRE: TBW
}


#else // ! eslENABLE_VMX
/* provide a callable function even when we're `./configure --disable-vmx` */
int
h4_vitfilter_vmx(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(hmm); ESL_UNUSED(mo); ESL_UNUSED(fx); ESL_UNUSED(ret_sc);  // shut up, compiler, I know what I'm doing
  esl_fatal("Altivec/VMX support was not enabled at compile time. Can't use h4_vitfilter_vmx().");
  return eslFAIL; // NOTREACHED
}
#endif // eslENABLE_VMX or not
