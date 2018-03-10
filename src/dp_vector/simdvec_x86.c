#include "p7_config.h"
#if defined eslENABLE_SSE || eslENABLE_AVX || eslENABLE_AVX512

#ifdef HAVE_FLUSH_ZERO_MODE
#include <xmmintrin.h>    // x86 SSE 
#endif
#ifdef HAVE_DENORMALS_ZERO_MODE
#include <pmmintrin.h>    // x86 SSE3 
#endif

void
p7_simdvec_x86_Init(void)
{
   /* Turn off denormalized floating-point, if possible.  Vectorized
   * prob-space Fwd/Bck filter underflows by design, and underflows
   * are negligible. See notes in simdvec.md.
   */
#ifdef HAVE_FLUSH_ZERO_MODE
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif  
#ifdef HAVE_DENORMAL_ZERO_MODE
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
}

#endif 
