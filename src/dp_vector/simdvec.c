/* Processor-specific initialization of SIMD vector environment.
 * Called by p7_Init() whenever a new thread or process starts.
 * 
 */

#include "p7_config.h"

#include <xmmintrin.h>    /* SSE  */
#include <emmintrin.h>    /* SSE2 */
#ifdef HAVE_PMMINTRIN_H
#include <pmmintrin.h>   /* DENORMAL_MODE */
#endif


#include "dp_vector/simdvec.h"

void
p7_simdvec_Init(void)
{
  /* In order to avoid the performance penalty dealing with denormalized
   * values in the floating point calculations, set the processor flag
   * so denormals are "flushed" immediately to zero.
   * This is the FTZ flag on an Intel CPU.
   */
#ifdef HAVE_FLUSH_ZERO_MODE
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif  

  /*
   * FLUSH_ZERO doesn't necessarily work in non-SIMD calculations
   * (yes on 64-bit, maybe not of 32-bit). This ensures that those
   * scalar calculations will agree across architectures.
   * (See TW notes  2012/0106_printf_underflow_bug/00NOTES for details)
   * This is the DAZ flag on an Intel CPU.
   */
#ifdef HAVE_PMMINTRIN_H
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
}


/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
