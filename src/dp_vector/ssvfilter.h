#ifndef p7SSVFILTER_INCLUDED
#define p7SSVFILTER_INCLUDED
#include <p7_config.h>

#include "easel.h"

#include "dp_vector/p7_filtermx.h"
#include "dp_vector/p7_oprofile.h"

extern int (*p7_SSVFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);

extern int p7_SSVFilter_sse     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);
extern int p7_SSVFilter_base_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *fx, float *ret_sc);
extern int p7_SSVFilter_avx     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc); 
extern int p7_SSVFilter_avx512  (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc); 
extern int p7_SSVFilter_neon    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc); 
extern int p7_SSVFilter_vmx     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc); 

#endif /*p7SSVFILTER_INCLUDED*/
