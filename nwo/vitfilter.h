#ifndef h4VITFILTER_INCLUDED
#define h4VITFILTER_INCLUDED
#include "h4_config.h"

#include "easel.h"

#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_filtermx.h"

extern int (*h4_vitfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);

#ifdef eslENABLE_SSE4
extern int h4_vitfilter_sse(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);
#endif

#ifdef eslENABLE_AVX
extern int h4_vitfilter_avx(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);
#endif

#ifdef eslENABLE_AVX512
extern int h4_vitfilter_avx512(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);
#endif

#ifdef eslENABLE_NEON
extern int h4_vitfilter_neon(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);
#endif

#ifdef eslENABLE_VMX
extern int h4_vitfilter_vmx(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);
#endif

#endif /*p7VITFILTER_INCLUDED*/
