/* AVX-specific initialization of SIMD vector environment
 *
 * Essentially identical to simdvec_sse.c; we just need an AVX version
 * in case we're compiling an AVX implementation without an SSE
 * implementation. See simdvec_sse.c, simdvec.c, and simdvec.md for
 * more explanation.
 */
#include <h4_config.h>

#include "simdvec.h"

#ifdef eslENABLE_AVX
#ifdef HAVE_FLUSH_ZERO_MODE
#include <xmmintrin.h>    // x86 SSE : FTZ
#endif
#ifdef HAVE_DENORMALS_ZERO_MODE
#include <pmmintrin.h>    // x86 SSE3 : DAZ
#endif

/* Function:  h4_simdvec_init_avx()
 * Synopsis:  AVX-specific initialization of SIMD vector environment
 */
void
h4_simdvec_init_avx(void)
{
#ifdef HAVE_FLUSH_ZERO_MODE
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif  
#ifdef HAVE_DENORMAL_ZERO_MODE
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
}


#else // ! eslENABLE_AVX
/* provide a callable function even when we're `./configure --disable-avx` */
void
h4_simdvec_init_avx(void)
{
  esl_fatal("AVX support was not enabled at compile time. Can't use h4_simdvec_init_avx().");
}
#endif // eslENABLE_AVX or not

