
#include "p7_config.h"

#include <xmmintrin.h>    /* SSE  */
#include <emmintrin.h>    /* SSE2 */
#ifdef _PMMINTRIN_H_INCLUDED
#include <pmmintrin.h>   /* DENORMAL_MODE */
#endif

#include "dp_vector/simdvec.h"

void
p7_simdvec_Init(void)
{
#ifdef HAVE_FLUSH_ZERO_MODE
  /* In order to avoid the performance penalty dealing with sub-normal
   * values in the floating point calculations, set the processor flag
   * so sub-normals are "flushed" immediately to zero.
   */
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif  
}

void
p7_simdvec_InitThread(void)
{
#ifdef HAVE_FLUSH_ZERO_MODE
  /* On OS/X, need to reset this flag for each thread
   * (see TJW notes 05/08/10 for details)
   */
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

#ifdef _PMMINTRIN_H_INCLUDED
  /*
   * FLUSH_ZERO doesn't necessarily work in non-SIMD calculations
   * (yes on 64-bit, maybe not of 32-bit). This ensures that those
   * scalar calculations will agree across architectures.
   * (See TW notes  2012/0106_printf_underflow_bug/00NOTES for details)
   */
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
}

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
