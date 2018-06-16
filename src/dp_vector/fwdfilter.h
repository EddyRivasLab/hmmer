#ifndef p7FWDFILTER_INCLUDED
#define p7FWDFILTER_INCLUDED
#include "p7_config.h"

#include "easel.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_sparse/p7_sparsemx.h"

extern int (*p7_ForwardFilter) (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int (*p7_BackwardFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);

#ifdef eslENABLE_SSE4
extern int p7_ForwardFilter_sse (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);
#endif

#ifdef eslENABLE_AVX
extern int p7_ForwardFilter_avx (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);
#endif

#ifdef eslENABLE_AVX512
extern int p7_ForwardFilter_avx512 (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);
#endif

#ifdef eslENABLE_NEON
extern int p7_ForwardFilter_neon (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);
#endif

#ifdef eslENABLE_VMX
extern int p7_ForwardFilter_vmx (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter_vmx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);
#endif
#endif /*p7FWDFILTER_INCLUDED*/

