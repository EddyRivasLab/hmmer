#ifndef p7SSVFILTER_INCLUDED
#define p7SSVFILTER_INCLUDED

#include "p7_config.h"
#include "easel.h"
#include "dp_vector/p7_oprofile.h"

extern int p7_SSVFilter    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);
extern int p7_SSVFilter_sse    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc); //SSE-only version
extern int p7_SSVFilter_avx    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc); //SSE-only version
extern int p7_SSVFilter_avx512    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc); //SSE-only version
#endif /*p7SSVFILTER_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
