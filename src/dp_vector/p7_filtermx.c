/* Implementation of P7_FILTERMX, the one-row DP memory for MSV and
 * Viterbi filters.
 * 
* Note: many of the functions in this file now just dispatch to the 
* appropriate SIMD-versioned function.  Where the SIMD type
* is known in advance, performance can be improved by calling
* the appropriate version directly.
 * Contents:
 *   1. The P7_FILTERMX object
 *   2. Debugging and development routines
 *   3. Copyright and license information
 */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>


#include "easel.h"
#include "esl_random.h"
#include "dp_vector/simdvec.h"
#include "dp_vector/p7_filtermx.h"
#include "hardware/hardware.h"
#include "base/general.h"
/*****************************************************************
 * 1. The P7_FILTERMX object.
 *****************************************************************/


/* Function:  p7_filtermx_Create()
 * Synopsis:  Create a one-row DP matrix for MSV, VF.
 *
 * Purpose:   Allocate a reusable, resizeable one-row <P7_FILTERMX>
 *            suitable for MSV and Viterbi filter calculations on 
 *            query profiles of up to <allocM> consensus positions.
 *            
 *            <allocM> must be $\leq$ 100,000. This is an H3 design
 *            limit.
 *            
 * Args:      allocM - initial allocation size, in profile positions. (>=1, <=100000)
 *
 * Returns:   ptr to new <P7_FILTERMX>
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_FILTERMX *
p7_filtermx_Create(int allocM, SIMD_TYPE simd)
{

  switch(simd){
    case SSE:
      return p7_filtermx_Create_sse(allocM);
      break;
    case AVX:
      return p7_filtermx_Create_avx(allocM);
      break;
    case AVX512:
      return p7_filtermx_Create_avx512(allocM);
      break;
    case NEON:
      return p7_filtermx_Create_neon(allocM);
      break;
    case NEON64:
      return p7_filtermx_Create_neon64(allocM);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_filtermx_Create");  
  }
  
}


/* Function:  p7_filtermx_GrowTo()
 * Synopsis:  Resize filter DP matrix for new profile size.
 *
 * Purpose:   Given an existing filter matrix structure <fx>,
 *            and the dimension <M> of a new profile that 
 *            we're going to use (in consensus positions),
 *            assure that <fx> is large enough for such a 
 *            profile; reallocate and reinitialize as needed.
 *
 *            <p7_filtermx_Reuse(fx); p7_filtermx_GrowTo(fx, M)>
 *            is essentially equivalent to <p7_filtermx_Create(M)>,
 *            while minimizing reallocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. The state of
 *            <fx> is now undefined, and it should not be used.
 */
int
p7_filtermx_GrowTo(P7_FILTERMX *fx, int allocM)
{
 switch(fx->simd){
    case SSE:
      return p7_filtermx_GrowTo_sse(fx, allocM);
      break;
    case AVX:
      return p7_filtermx_GrowTo_avx(fx, allocM);
      break;
    case AVX512:
      return p7_filtermx_GrowTo_avx512(fx, allocM);
      break;
    case NEON:
      return p7_filtermx_GrowTo_neon(fx, allocM);
      break;
    case NEON64:
      return p7_filtermx_GrowTo_neon64(fx, allocM);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_filtermx_GrowTo");  
  }
}


/* Function:  p7_filtermx_Sizeof()
 * Synopsis:  Calculate and return the current size, in bytes.
 *
 * Purpose:   Calculate and return the current allocated size
 *            of the filter matrix <fx>, in bytes.
 *            
 *            Used in diagnostics, benchmarking of memory usage.
 *
 * Returns:   the allocation size.
 */
size_t 
p7_filtermx_Sizeof(const P7_FILTERMX *fx)
{

  switch(fx->simd){
    case SSE:
      return p7_filtermx_Sizeof_sse(fx);
      break;
    case AVX:
      return p7_filtermx_Sizeof_avx(fx);
      break;
    case AVX512:
      return p7_filtermx_Sizeof_avx512(fx);
      break;
    case NEON:
      return p7_filtermx_Sizeof_neon(fx);
      break;
    case NEON64:
      return p7_filtermx_Sizeof_neon64(fx);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_filtermx_Sizeof");  
  }
}


/* Function:  p7_filtermx_MinSizeof()
 * Synopsis:  Calculate minimum size of a filter matrix, in bytes.
 *
 * Purpose:   Calculate and return the minimum allocation size
 *            required for a one-row filter matrix for a query 
 *            profile of <M> consensus positions.
 */
size_t
p7_filtermx_MinSizeof(int M, SIMD_TYPE simd)
{
switch(simd){
    case SSE:
      return p7_filtermx_MinSizeof_sse(M);
      break;
    case AVX:
      return p7_filtermx_MinSizeof_avx(M);
      break;
    case AVX512:
      return p7_filtermx_MinSizeof_avx512(M);
      break;
    case NEON:
      return p7_filtermx_MinSizeof_neon(M);
      break;
    case NEON64:
      return p7_filtermx_MinSizeof_neon64(M);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_filtermx_MinSizeof");  
  }
}


/* Function:  p7_filtermx_Reuse()
 * Synopsis:  Recycle a P7_FILTERMX.
 *
 * Purpose:   Recycle the filter DP matrix <fx>: reinitialize
 *            information having to do with the current calculation,
 *            but leave the allocation so we can reuse it.
 *            
 *            Debugging flags are reset too. If the caller reuses this
 *            matrix on another profile, and wants to dump those
 *            comparisons for diagnostics, it must call
 *            <p7_filtermx_SetDumpMode()> again.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_filtermx_Reuse(P7_FILTERMX *fx)
{
  fx->M         = 0;
  fx->type      = p7F_NONE;
#ifdef p7_DEBUGGING
  fx->do_dumping= FALSE;
  fx->dfp       = NULL;
#endif
  return eslOK;
}


/* Function:  p7_filtermx_Destroy()
 * Synopsis:  Frees a one-row MSV/VF filter DP matrix.
 *
 * Purpose:   Frees the one-row MSV/VF filter DP matrix <fx>.
 *
 * Returns:   (void)
 */
void
p7_filtermx_Destroy(P7_FILTERMX *fx)
{
  switch(fx->simd){
    case SSE:
      p7_filtermx_Destroy_sse(fx);
      break;
    case AVX:
      p7_filtermx_Destroy_avx(fx);
      break;
    case AVX512:
      p7_filtermx_Destroy_avx512(fx);
      break;
    case NEON:
      p7_filtermx_Destroy_neon(fx);
      break;
    case NEON64:
      p7_filtermx_Destroy_neon64(fx);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_filtermx_Destroy");  
  }
  return;
}
  

/*****************************************************************
 * 2. Debugging and development routines
 *****************************************************************/

/* Function:  p7_filtermx_SetDumpMode()
 * Synopsis:  Set the dump mode flag in a P7_FILTERMX.
 *
 * Purpose:   Set the dump mode flag in <fx> to <truefalse>,
 *            and the diagnostic stream to <dfp>.
 *            
 *            To turn on dumping, <truefalse> is <TRUE>, and <dfp> is
 *            an open file pointer (usually <stdout> or <stderr>).
 *            
 *            To turn off dumping, <truefalse> is <FALSE>, and <dfp>
 *            is <NULL>.
 *            
 *            For this to have any effect, <p7_DEBUGGING> compile-time
 *            flag must be up. If it is not, no dumping will occur,
 *            and this call is a no-op.
 *            
 * Args:      fx        - DP matrix object to set for diagnostic dumping
 *            dfp       - open FILE * for diagnostic output, or NULL     
 *            truefalse - TRUE to set, FALSE to unset.
 * 
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_filtermx_SetDumpMode(P7_FILTERMX *fx, FILE *dfp, int truefalse)
{
#ifdef p7_DEBUGGING
  fx->do_dumping = truefalse;
  fx->dfp        = dfp;
#endif
  return eslOK;
}

#ifdef p7_DEBUGGING
/* Function:  p7_filtermx_DumpMFRow()
 * Synopsis:  Dump one row from MSV version of a DP matrix.
 *
 * Purpose:   Dump current row of MSV calculations from DP matrix <fx>
 *            for diagnostics, and include the values of specials
 *            <xE>, etc. The index <rowi> for the current row is used
 *            as a row label. This routine has to be specialized for
 *            the layout of the MSVFilter() row, because it's all
 *            match scores dp[0..q..Q-1], rather than triplets of
 *            M,D,I.
 * 
 *            If <rowi> is 0, print a header first too.
 * 
 *            The output format is coordinated with <p7_refmx_Dump()> to
 *            facilitate comparison to a known answer.
 *            
 *            This also works for an SSV filter row, for SSV implementations
 *            that use a single row of DP memory (like <_longtarget>). 
 *            The Knudsen assembly code SSV does not use any RAM.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. 
 */
int
p7_filtermx_DumpMFRow(const P7_FILTERMX *fx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC)
{
 switch(fx->simd){
    case SSE:
      return p7_filtermx_DumpMFRow_sse(fx, rowi, xE, xN, xJ, xB, xC);
      break;
    case AVX:
      return p7_filtermx_DumpMFRow_avx(fx, rowi, xE, xN, xJ, xB, xC);
      break;
    case AVX512:
      return p7_filtermx_DumpMFRow_avx512(fx, rowi, xE, xN, xJ, xB, xC);
      break;
    case NEON:
      return p7_filtermx_DumpMFRow_neon(fx, rowi, xE, xN, xJ, xB, xC);
      break;
    case NEON64:
      return p7_filtermx_DumpMFRow_neon64(fx, rowi, xE, xN, xJ, xB, xC);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_filtermx_DumpMFRow");  
  }
}



/* Function:  p7_filtermx_DumpVFRow()
 * Synopsis:  Dump current row of ViterbiFilter (int16) filter matrix.
 *
 * Purpose:   Dump current row of ViterbiFilter (int16) filter DP
 *            matrix <fx> for diagnostics, and include the values of
 *            specials <xE>, etc. The index <rowi> for the current row
 *            is used as a row label.
 *
 *            If <rowi> is 0, print a header first too.
 * 
 *            The output format is coordinated with <p7_refmx_Dump()> to
 *            facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_filtermx_DumpVFRow(const P7_FILTERMX *fx, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC)
{
  switch(fx->simd){
    case SSE:
      return p7_filtermx_DumpVFRow_sse(fx, rowi, xE, xN, xJ, xB, xC);
      break;
    case AVX:
      return p7_filtermx_DumpVFRow_avx(fx, rowi, xE, xN, xJ, xB, xC);
      break;
    case AVX512:
      return p7_filtermx_DumpVFRow_avx512(fx, rowi, xE, xN, xJ, xB, xC);
      break;
    case NEON:
      return p7_filtermx_DumpVFRow_neon(fx, rowi, xE, xN, xJ, xB, xC);
      break;
    case NEON64:
      return p7_filtermx_DumpVFRow_neon64(fx, rowi, xE, xN, xJ, xB, xC);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_filtermx_DumpVFRow");  
  }
}
#endif /*p7_DEBUGGING*/



/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
