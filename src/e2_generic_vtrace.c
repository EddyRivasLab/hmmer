/* Forward/Backward algorithms; generic (non-SIMD) versions.
 * 
 * Contents:
 *   1. Forward, Backward implementations.  
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "e2.h"
#include "e2_generic_vtrace.h"
#include "e2_profile.h"
#include "e2_profilesq.h"
#include "e2_trace.h"

/*****************************************************************
 * 1. Viterbi traceback implementation.
 *****************************************************************/

/* Function:  e2_GTrace()
 *
 * Purpose:  Traceback of a Viterbi matrix: retrieval 
 *           of optimum alignment.
 *           
 *           This function is currently implemented as a
 *           reconstruction traceback, rather than using a shadow
 *           matrix. Because H3 uses floating point scores, and we
 *           can't compare floats for equality, we have to compare
 *           floats for near-equality and therefore, formally, we can
 *           only guarantee a near-optimal traceback. However, even in
 *           the unlikely event that a suboptimal is returned, the
 *           score difference from true optimal will be negligible.
 *           
 * Args:     dsq    - digital sequence aligned to, 1..L 
 *           L      - length of <dsq>
 *           gm     - profile
 *           mx     - Viterbi matrix to trace, L x M
 *           tr     - storage for the recovered traceback.
 *           
 * Return:   <eslOK> on success.
 *           <eslFAIL> if even the optimal path has zero probability;
 *           in this case, the trace is set blank (<tr->N = 0>).
 *
 * Note:     Care is taken to evaluate the prev+tsc+emission
 *           calculations in exactly the same order that Viterbi did
 *           them, lest you get numerical problems with
 *           a+b+c = d; d-c != a+b because d,c are nearly equal.
 *           (This bug appeared in dev: xref J1/121.)
 */

int
e2_GTrace(const PSQ *psql, const PSQ *psqr, const E2_PROFILE *gm, E2_GMX *vit, E2_TRACE *tr)
{
 
 
  return eslOK;
}

