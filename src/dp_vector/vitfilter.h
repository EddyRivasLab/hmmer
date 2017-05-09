#ifndef p7VITFILTER_INCLUDED
#define p7VITFILTER_INCLUDED
#include "p7_config.h"

#include "easel.h"

#include "base/p7_hmmwindow.h"
#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"

extern int (*p7_ViterbiFilter)           (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int (*p7_ViterbiFilter_longtarget)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);

#ifdef eslENABLE_SSE
extern int p7_ViterbiFilter_sse           (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_longtarget_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);
#endif

#ifdef eslENABLE_AVX
extern int p7_ViterbiFilter_avx           (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_longtarget_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);
#endif

#ifdef eslENABLE_AVX512
extern int p7_ViterbiFilter_avx512           (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_longtarget_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);
#endif

#ifdef eslENABLE_NEON
extern int p7_ViterbiFilter_neon           (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_longtarget_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);
#endif

#endif /*p7VITFILTER_INCLUDED*/
