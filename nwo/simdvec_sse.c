/* SSE-specific initialization of SIMD vector environment
 *
 * simdvec.c itself is not allowed to depend on SIMD headers or
 * CFLAGS, but to turn on FTZ/DAZ CPU flags on x86 platforms, we need
 * to include SSE and SSE3 headers.  If you put this code in
 * simdvec.c, it will fail to compile on the few platforms that need a
 * specific compilation flag in SIMD_CFLAGS to compile SSE3 code
 * (gcc4 is an example).
 *
 * See simdvec.c, simdvec.md for further explanations.
 */
#include <h4_config.h>

#include "simdvec.h"

#ifdef eslENABLE_SSE4
#ifdef HAVE_FLUSH_ZERO_MODE
#include <xmmintrin.h>    // x86 SSE : FTZ
#endif
#ifdef HAVE_DENORMALS_ZERO_MODE
#include <pmmintrin.h>    // x86 SSE3 : DAZ
#endif

/* Function:  h4_simdvec_init_sse()
 * Synopsis:  SSE-specific initialization of SIMD vector environment
 * Incept:    SRE, Thu 05 Mar 2020
 *
 * Purpose:   On x86 platforms, we set two CPU mode flags, FTZ
 *            (flush-to-zero) and DAZ (denormals-are-zero), for
 *            floating point performance reasons. See simdvec.md for
 *            further explanation.
 *
 *            We need to do the same thing for sse, avx, avx512
 *            implementations, but our Makefiles don't provide the
 *            ability to compile one _x86() function for all three
 *            implementations to call.  Moreover we can't just provide
 *            the _sse() version and call it from avx/avx512 because
 *            it's possible to compile for avx/avx512 alone, without
 *            sse. Therefore we provide all three init functions,
 *            doing the same thing; only one will be called by
 *            h4_simdvec_Init().
 */
void
h4_simdvec_init_sse(void)
{
#ifdef HAVE_FLUSH_ZERO_MODE
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif  
#ifdef HAVE_DENORMAL_ZERO_MODE
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
}


#else // ! eslENABLE_SSE4
/* provide a callable function even when we're `./configure --disable-sse` */
void
h4_simdvec_init_sse(void)
{
  esl_fatal("SSE support was not enabled at compile time. Can't use h4_simdvec_init_sse().");
}
#endif // eslENABLE_SSE4 or not

