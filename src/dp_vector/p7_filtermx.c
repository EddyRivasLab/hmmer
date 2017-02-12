/* P7_FILTERMX: one-row of DP memory for MSV and Viterbi filters.
 * 
 *   1. The P7_FILTERMX object
 *   2. Debugging and development routines
 */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>


#include "easel.h"
#include "esl_random.h"

#include "dp_vector/simdvec.h"
#include "dp_vector/p7_filtermx.h"

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
 *            <allocM> must be $\leq$ 100,000. This is a design
 *            limit.
 *            
 *            As a special case, <allocM> may be 0, in which case an
 *            empty structure is created, and <p7_filtermx_Reinit()>
 *            must be called on it to finish the allocation before it
 *            is used. This is a frill, supporting a pattern where we
 *            have to have a <_Reinit()> inside a loop anyway,
 *            following something that gives us the size of the first
 *            (and subsequent) problems.
 *            
 * Args:      allocM - initial allocation size, in profile positions. (>=1, <=100000)
 *
 * Returns:   ptr to new <P7_FILTERMX>
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_FILTERMX *
p7_filtermx_Create(int allocM)
{
  P7_FILTERMX *fx = NULL;
  int status;

  /* Contract checks / argument validation */
  ESL_DASSERT1( (allocM >= 0 && allocM <= 100000) ); // 0 is ok here, gives you a shell.

  ESL_ALLOC(fx, sizeof(P7_FILTERMX));
  fx->M         = 0;
  fx->V         = 0;
  fx->dp        = NULL;
  fx->allocM    = 0;
  fx->type      = p7F_NONE;
#if eslDEBUGLEVEL > 0
  fx->do_dumping= FALSE;
  fx->dfp       = NULL;
#endif 
  
  if (allocM) {
    /* VF needs, in bytes:         3 MDI cells *     Qmax vectors         * bytes/vector; (aligned) */
    fx->dp     = esl_alloc_aligned(p7F_NSCELLS * P7_NVW(allocM, p7_MAXVW) * p7_MAXVW,     p7_MAXVB);
    if (fx->dp == NULL) { status = eslEMEM; goto ERROR; }
    fx->allocM = allocM;
  }
  return fx;

 ERROR:
  p7_filtermx_Destroy(fx);
  return status;
}


/* Function:  p7_filtermx_Reinit()
 * Synopsis:  Resize filter DP matrix for new profile size.
 *
 * Purpose:   Given an existing filter matrix structure <fx>,
 *            and the dimension <M> of a new profile that 
 *            we're going to use (in consensus positions),
 *            assure that <fx> is large enough for such a 
 *            profile; reallocate and reinitialize as needed.
 *
 *            Debugging flags (if any) are also cleared, so
 *            caller may need to do <_SetDumpMode()> again.
 *
 *            <p7_filtermx_Reinit(fx,M)> is essentially equivalent to
 *            <p7_filtermx_Create(M)>, while minimizing reallocation.
 *            It should only be called in a similar context as a
 *            p7_filtermx_Create() call. Existing data may be lost,
 *            because for efficiency, old data are not copied to new
 *            internal allocations.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. <fx> is unchanged
 *            and its data remain valid.
 */
int
p7_filtermx_Reinit(P7_FILTERMX *fx, int allocM)
{
  char *p;
  int   status;

  /* Contract checks / argument validation */
  ESL_DASSERT1( (allocM >= 1 && allocM <= 100000) );

  /* Reinitialization to an empty structure */
  fx->M          = 0;
  fx->V          = 0;
  fx->type       = p7F_NONE;
#if eslDEBUGLEVEL > 0
  fx->do_dumping = FALSE;
  fx->dfp        = NULL;
#endif 

  /* Reallocate if needed */
  if (allocM > fx->allocM) {
    free(fx->dp);
    fx->dp     = esl_alloc_aligned(p7F_NSCELLS * P7_NVW(allocM, p7_MAXVW) * p7_MAXVW,  p7_MAXVB);
    if (fx->dp == NULL) { status = eslEMEM; goto ERROR; }
    fx->allocM = allocM;
  }
  return eslOK;

 ERROR:
  fx->allocM = 0;
  return status;
}

/* Function:  p7_filtermx_GetScore()
 * Synopsis:  Get a score from striped DP memory.
 * Incept:    SRE, Sun Feb 12 10:46:28 2017 [The Chemical Brothers]
 *
 * Purpose:   Return the score at position <k>, state <s> in the striped
 *            vector DP row of <fx>. <k> is a model coordinate from 1..M.
 *            <s> is a state flag, one of <p7F_M | p7F_D | p7F_I>.
 *
 *            Slow, because it accesses the data in scalar form, and
 *            has to do a convoluted coord conversion from <k,s> to
 *            striped <q,r> coords. Only use this in noncritical code
 *            like debugging comparisons and dumps.
 *
 *            Independent of vector ISA, because all it needs to
 *            unstripe is the vector width <V>. Works for MSVFilter,
 *            Vitfilter, and SSV_longtarget matrices.
 *
 *            MSV and SSV matrices only have p7F_M scores. If a D or I
 *            score is asked for, returns 0.
 *
 *            If k==0, returns 0, as a special case that we use in
 *            debugging dumps to sync to the output of
 *            p7_refmx_Dump().
 *
 * Args:      k : position (0).1..M
 *            s : state,   p7F_M | p7F_D | p7F_I
 *
 * Returns:   the score
 */
int
p7_filtermx_GetScore(P7_FILTERMX *fx, int k, enum p7f_scells_e s)
{
  int q, r;

  ESL_DASSERT1(( k >= 0 &&  k <= fx->allocM ));
  ESL_DASSERT1(( fx->V &&  !(fx->V & (fx->V-1)) ));     // V is set, to a power of two.
  ESL_DASSERT1(( fx->M ));                              // matrix has to have some data in it

  if (k == 0) return 0;

  if (fx->type == p7F_MSVFILTER || fx->type == p7F_SSVFILTER)
    {
      if (s != p7F_M) return 0;
      uint8_t *dp = (uint8_t *) fx->dp;
      int      q  = (k-1) % fx->V;
      int      r  = (k-1) / fx->V;
      return dp[q*fx->V + r];
    }
  else if (fx->type == p7F_VITFILTER)
    {
      int16_t *dp = (int16_t *) fx->dp;
      int      q  = (k-1) % fx->V;
      int      r  = (k-1) / fx->V;
      return dp[ (q*P7F_NSCELLS + s) * fx->V + r];
    }
  
}

/* Function:  p7_filtermx_Sizeof()
 * Synopsis:  Calculate and return current size, in bytes.
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
  size_t n = sizeof(P7_FILTERMX);
  n       += p7F_NSCELLS * P7_NVW(fx->allocM, p7_MAXVW) * p7_MAXVW; 
  // neglects any alignment overhead, which may be up to p7_MAXVB bytes.
  return n;
}


/* Function:  p7_filtermx_MinSizeof()
 * Synopsis:  Calculate minimum size of a filter matrix, in bytes.
 *
 * Purpose:   Calculate and return the minimum allocation size
 *            required for a one-row filter matrix for a query 
 *            profile of <M> consensus positions.
 */
size_t
p7_filtermx_MinSizeof(int M)
{
  size_t n = sizeof(P7_FILTERMX);
  n       += p7F_NSCELLS * P7_NVW(fx->allocM, p7_MAXVW) * p7_MAXVW; 
  // neglects any alignment overhead
  return n;
}


/* Function:  p7_filtermx_Destroy()
 * Synopsis:  Frees a one-row MSV/VF filter DP matrix.
 *
 * Purpose:   Frees the one-row MSV/VF filter DP matrix <fx>.
 */
void
p7_filtermx_Destroy(P7_FILTERMX *fx)
{
  if (fx) {
    if (fx->dp) esl_alloc_free(fx->dp);  // not free(); this is an aligned mem ptr
    free(fx);
  }
}
  



/*****************************************************************
 * 2. Debugging and development routines
 *****************************************************************/
#if eslDEBUGLEVEL > 0

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
  fx->do_dumping = truefalse;
  fx->dfp        = dfp;
  return eslOK;
}


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
 */
int
p7_filtermx_DumpMFRow(const P7_FILTERMX *fx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC)
{
  ESL_DASSERT1(( fx->type == p7F_MSVFILTER || fx->type == p7F_SSVFILTER ));
  ESL_DASSERT1(( fx->V &&  !(fx->V & (fx->V-1)) ));     // V is set, to a power of two.
  ESL_DASSERT1(( fx->M ));
  
  /* Header (if we're on the 0th row)  */
  if (rowi == 0)
    {
      fprintf(fx->dfp, "       ");
      for (k = 0; k <= fx->M;  k++) fprintf(fx->dfp, "%3d ", k);
      fprintf(fx->dfp, "%3s %3s %3s %3s %3s\n", "E", "N", "J", "B", "C");
      fprintf(fx->dfp, "       ");
      for (k = 0; k <= fx->M+5;  k++) fprintf(fx->dfp, "%3s ", "---");
      fprintf(fx->dfp, "\n");
    }

  /* Unstripe and print match scores. */
  fprintf(fx->dfp, "%4d M ", rowi);
  for (k = 0; k <= fx->M; k++) fprintf(fx->dfp, "%3d ", p7_filtermx_GetScore(fx, k, p7F_M);

  /* Specials */
  fprintf(fx->dfp, "%3d %3d %3d %3d %3d\n", xE, xN, xJ, xB, xC);

  /* I's are all 0's; print just to facilitate comparison to refmx. */
  fprintf(fx->dfp, "%4d I ", rowi);
  for (k = 0; k <= fx->M; k++) fprintf(fx->dfp, "%3d ", 0);
  fprintf(fx->dfp, "\n");

  /* D's are all 0's too */
  fprintf(fx->dfp, "%4d D ", rowi);
  for (k = 0; k <= fx->M; k++) fprintf(fx->dfp, "%3d ", 0);
  fprintf(fx->dfp, "\n\n");

  return eslOK;
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
  ESL_DASSERT1(( fx->type == p7F_VITFILTER ));
  ESL_DASSERT1(( fx->V &&  !(fx->V & (fx->V-1)) ));     // V is set, to a power of two.
  ESL_DASSERT1(( fx->M ));

  /* Header (if we're on the 0th row) */
  if (rowi == 0)
    {
      fprintf(fx->dfp, "       ");
      for (k = 0; k <= fx->M;  k++) fprintf(fx->dfp, "%6d ", k);
      fprintf(fx->dfp, "%6s %6s %6s %6s %6s\n", "E", "N", "J", "B", "C");
      fprintf(fx->dfp, "       ");
      for (k = 0; k <= fx->M+5;  k++) fprintf(fx->dfp, "%6s ", "------");
      fprintf(fx->dfp, "\n");
    }

  /* Unstripe and print match scores. */
  fprintf(fx->dfp, "%4d M ", rowi);
  for (k = 0; k <= fx->M; k++) fprintf(fx->dfp, "%6d ", p7_filtermx_GetScore(fx, k, p7F_M));
  
  /* Specials */
  fprintf(fx->dfp, "%6d %6d %6d %6d %6d\n", xE, xN, xJ, xB, xC);

  /* Insert scores. */
  fprintf(fx->dfp, "%4d I ", rowi);
  for (k = 0; k <= fx->M; k++) fprintf(fx->dfp, "%6d ", p7_filtermx_GetScore(fx, k, p7F_I));
  fprintf(fx->dfp, "\n");

  /* Delete scores. */
  fprintf(fx->dfp, "%4d D ", rowi);
  for (k = 0; k <= fx->M; k++) fprintf(fx->dfp, "%6d ", p7_filtermx_GetScore(fx, k, p7F_D));
  fprintf(fx->dfp, "\n\n");

  return eslOK;
}
#endif // eslDEBUGLEVEL

