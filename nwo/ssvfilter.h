#ifndef h4SSVFILTER_INCLUDED
#define h4SSVFILTER_INCLUDED
#include "h4_config.h"

#include "easel.h"

#include "h4_profile.h"

extern int (*h4_ssvfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc);


extern int h4_ssvfilter_sse   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc);
extern int h4_ssvfilter_avx   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc);
extern int h4_ssvfilter_avx512(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc);
extern int h4_ssvfilter_neon  (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc);
extern int h4_ssvfilter_vmx   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc);

#endif // h4SSVFILTER_INCLUDED
