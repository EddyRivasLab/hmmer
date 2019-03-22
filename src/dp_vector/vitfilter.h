#ifndef p7VITFILTER_INCLUDED
#define p7VITFILTER_INCLUDED
#include "p7_config.h"

#include "easel.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
extern int (*p7_ViterbiFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);

#ifdef eslENABLE_SSE4
extern int p7_ViterbiFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
#endif

#ifdef eslENABLE_AVX
extern int p7_ViterbiFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
#endif

#ifdef eslENABLE_AVX512
extern int p7_ViterbiFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
#endif

#ifdef eslENABLE_NEON
extern int p7_ViterbiFilter_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
#endif

#ifdef eslENABLE_VMX
extern int p7_ViterbiFilter_vmx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
#endif
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7VITFILTER_INCLUDED*/
