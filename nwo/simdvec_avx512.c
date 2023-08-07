/* AVX512-specific initialization of SIMD vector environment
 *
 * Essentially identical to simdvec_sse.c; we just need an AVX512 version
 * in case we're compiling an AVX512 implementation without an SSE
 * implementation. See simdvec_sse.c, simdvec.c, and simdvec.md for
 * more explanation.
 */
#include <h4_config.h>

#include "simdvec.h"

#ifdef eslENABLE_AVX512
#ifdef HAVE_FLUSH_ZERO_MODE
#include <xmmintrin.h>    // x86 SSE : FTZ
#endif
#ifdef HAVE_DENORMALS_ZERO_MODE
#include <pmmintrin.h>    // x86 SSE3 : DAZ
#endif

/* Function:  h4_simdvec_init_avx512()
 * Synopsis:  AVX512-specific initialization of SIMD vector environment
 */
void
h4_simdvec_init_avx512(void)
{
#ifdef HAVE_FLUSH_ZERO_MODE
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif  
#ifdef HAVE_DENORMAL_ZERO_MODE
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
}


#else // ! eslENABLE_AVX512
/* provide a callable function even when we're `./configure --disable-avx512` */
void
h4_simdvec_init_avx512(void)
{
  esl_fatal("AVX512 support was not enabled at compile time. Can't use h4_simdvec_init_avx512().");
}
#endif // eslENABLE_AVX512 or not

