#ifndef p7FWDFILTER_INCLUDED
#define p7FWDFILTER_INCLUDED

#include "p7_config.h"

#include "easel.h"
#if p7_CPU_ARCH == x86
#include <xmmintrin.h>
#include <emmintrin.h>
#endif
#ifdef HAVE_NEON
#include "esl_neon.h"
#endif

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_sparse/p7_sparsemx.h"

typedef
#if p7_CPU_ARCH == intel 
#if defined HAVE_AVX512
__m512
#elif defined HAVE_AVX2
__m256
#else
__m128
#endif
#endif /* intel */
#if p7_CPU_ARCH == arm || p7_CPU_ARCH == arm64
esl_neon_128f_t
#endif
debug_print;

extern int p7_ForwardFilter (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);
extern int p7_ForwardFilter_sse (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);
extern int p7_ForwardFilter_avx (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);
extern int p7_ForwardFilter_avx512 (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);
extern int p7_ForwardFilter_neon (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);
extern int p7_ForwardFilter_neon64 (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
extern int p7_BackwardFilter_neon64(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);



#endif /*p7FWDFILTER_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
