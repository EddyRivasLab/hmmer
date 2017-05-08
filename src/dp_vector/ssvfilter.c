/* SSV (Single Segment Viterbi) filter
 * 
 * See ssvfilter.md for notes.
 *
 * This is only the CPU dispatch front end, which dispatches a
 * p7_SSVFilter() call to the appropriate vector code. See
 * ssvfilter_{sse,avx...} for the various available vector
 * implementations.

 * Contents:
 *   1. p7_SSVFilter() API
 *   2. CPU dispatching to vector implementations.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/ssvfilter.h"

static int ssvfilter_dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);



/*****************************************************************
 * 1. SSVFilter API call
 *****************************************************************/

/* Function:  p7_SSVFilter()
 * Synopsis:  The SSV filter, the first step in the acceleration pipeline
 * Incept:    SRE, Mon May  8 07:37:41 2017 [HHGTTG, The Campaign for Real Time]
 *
 * Purpose:   Calculates approximate SSV score for digital sequence <dsq>
 *            of length <L> residues, using vector profile <om>. Return
 *            the SSV score, in nats, in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring sequences,
 *            but will not underflow.
 *            
 *            The model <om> may be in any mode. Only its match
 *            emission scores will be used. The SSV filter inherently
 *            assumes a singlehit local mode, and uses its own special
 *            state transitions scores, not the scores in the profile.
 *
 * Args:      dsq    - digital target sequence 1..L
 *            L      - length of <dsq> in residues
 *            om     - optimized profile
 *            ret_sc - RETURN: SSV score in nats
 *
 * Returns:   <eslOK> on success, and <*ret_sc> is the SSV score.
 * 
 *            <eslERANGE> if score overflows limited range. In this
 *            case, this is a high-scoring hit that passes the filter,
 *            and <*ret_sc> is set to its maximum allowed value.
 *
 * Throws:    (no abnormal error conditions)
 */
int
(*p7_SSVFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc) = 
  ssvfilter_dispatcher;




/*****************************************************************
 * 2. CPU dispatching to vector implementations.
 *****************************************************************/

static int 
ssvfilter_dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc) 
{
#ifdef eslENABLE_AVX512  // Fastest implementations first here.
  if (esl_cpu_has_avx512())
    {
      p7_SSVFilter = p7_SSVFilter_sse;
      return p7_SSVFilter_sse(dsq, L, om, ret_sc);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      p7_SSVFilter = p7_SSVFilter_avx;
      return p7_SSVFilter_avx(dsq, L, om, ret_sc);
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse())
    {
      p7_SSVFilter = p7_SSVFilter_sse;
      return p7_SSVFilter_sse(dsq, L, om, ret_sc);
    }
#endif
  
#ifdef eslENABLE_NEON
  p7_SSVFilter = p7_SSVFilter_neon;
  return p7_SSVFilter_neon(dsq, L, om, ret_sc);
#endif

  //#ifdef eslENABLE_VMX
  //  p7_SSVFilter = p7_SSVFilter_vmx;
  //  return p7_SSVFilter_vmx(dsq, L, om, ret_sc);
  //#endif

  p7_Die("ssvfilter_dispatcher found no vector implementation - that shouldn't happen.");
}

