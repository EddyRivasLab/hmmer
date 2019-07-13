#ifndef p7VITFILTER_INCLUDED
#define p7VITFILTER_INCLUDED
#include "p7_config.h"

#include "easel.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"

extern int (*p7_ViterbiFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);

extern int p7_ViterbiFilter_sse   (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_avx   (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_neon  (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_vmx   (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);

#endif /*p7VITFILTER_INCLUDED*/
