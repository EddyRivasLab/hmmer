/* Processor-specific initialization of SIMD vector environment.
 * Called by h4_Init() whenever a new thread or process starts.
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "simdvec.h"

#if defined eslENABLE_SSE4 || eslENABLE_AVX || eslENABLE_AVX512
#ifdef HAVE_FLUSH_ZERO_MODE
#include <xmmintrin.h>    // x86 SSE : FTZ
#endif
#ifdef HAVE_DENORMALS_ZERO_MODE
#include <pmmintrin.h>    // x86 SSE3 : DAZ
#endif
#endif

void
h4_simdvec_Init(void)
{
#if defined eslENABLE_SSE4 || eslENABLE_AVX || eslENABLE_AVX512
  /* On x86 platforms, turn off denormalized floating-point, if possible.  
   * Vectorized prob-space Fwd/Bck filter underflows by design, and underflows
   * are negligible. See notes in simdvec.md.
   */
#ifdef HAVE_FLUSH_ZERO_MODE
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif  
#ifdef HAVE_DENORMAL_ZERO_MODE
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
#endif
}



/* Function:  h4_simdvec_width()
 * Synopsis:  Returns SIMD vector width, in bytes
 * Incept:    SRE, Tue Feb 21 08:17:34 2017 [Bear McCreary]
 *
 * Purpose:   Returns the SIMD vector width, in bytes, for the vector
 *            implementation that this process/thread is (or will be)
 *            using. Possible answers are 16 (SSE, NEON, VMX); 32
 *            (AVX); 64 (AVX-512).
 *            
 *            Having this one function allows other code to be vector
 *            ISA independent. For example, <h4_profile> routines only
 *            need to know the vector width <V> to stripe data; they
 *            don't need to know which ISA is being used.
 */
int
h4_simdvec_width(void)
{
  static int V = -1;  // decide once, and store the answer. This is threadsafe.

  return 16; // SRE: temporary, while we're only testing SSE.

  if (V == -1) 
    {
#ifdef eslENABLE_AVX512
      if (esl_cpu_has_avx512()) { V = 64; goto DONE; }
#endif
#ifdef eslENABLE_AVX
      if (esl_cpu_has_avx())    { V = 32;  goto DONE; }
#endif
#ifdef eslENABLE_SSE4
      if (esl_cpu_has_sse4())   { V = 16;  goto DONE; }
#endif
#ifdef eslENABLE_NEON
      V = 16; goto DONE;
#endif
#ifdef eslENABLE_VMX
      V = 16; goto DONE;
#endif
    }

 DONE:
  if (V == -1) esl_fatal("found no vector implementation - this shouldn't happen");
  return V;
}



/* Function:  h4_simdvec_byteify()
 * Synopsis:  Convert log2-odds score to rounded, scaled int8_t
 * Incept:    SRE, Tue 28 May 2019
 *
 * Purpose:   For SSV filter emission scores.
 */
int8_t
h4_simdvec_byteify(float sc)
{
  sc = roundf(h4_SCALE_B * sc);
  if      (sc < -128.) return -128;          // can happen for sc = -inf. Otherwise shouldn't happen.
  else if (sc > 127.)  return 127;
  else                 return (int8_t) sc;
}

/* Function:  h4_simdvec_wordify()
 * Synopsis:  Convert log2-odds score to rounded, scaled int16_t
 * Incept:    SRE, Tue 28 May 2019
 *
 * Purpose:   For Viterbi filter emission and transition scores,
 *            including specials in H4_MODE.
 */
int16_t 
h4_simdvec_wordify(float sc)
{
  sc  = roundf(h4_SCALE_W * sc);
  if      (sc < -32768.) return -32768;
  else if (sc >  32727.) return  32767;
  else                   return (int16_t) sc;
}



