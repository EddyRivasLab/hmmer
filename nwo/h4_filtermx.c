/* H4_FILTERMX: one row of DP memory for Viterbi filter.
 * 
 * Independent of vector ISA. (Do not add any ISA-specific code.)
 *
 *   1. The H4_FILTERMX object
 *   2. Debugging and development routines
 */
#include "h4_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alloc.h"

#include "simdvec.h"
#include "h4_filtermx.h"


/*****************************************************************
 * 1. The H4_FILTERMX object.
 *****************************************************************/

/* Function:  h4_filtermx_Create()
 * Synopsis:  Create a one-row DP matrix for VF.
 *
 * Purpose:   Allocate a reusable, resizeable one-row <H4_FILTERMX>
 *            suitable for Viterbi filter calculations on 
 *            query profiles of up to <allocM> consensus positions.
 *            
 *            <allocM> must be $\leq$ 100,000. This is a design
 *            limit.
 *            
 *            As a special case, <allocM> may be 0, in which case an
 *            empty structure is created, and <h4_filtermx_Reinit()>
 *            must be called on it to finish the allocation before it
 *            is used. This is a frill, supporting a pattern where we
 *            have to have a <_Reinit()> inside a loop anyway,
 *            following something that gives us the size of the first
 *            (and subsequent) problems.
 *            
 * Args:      allocM - initial allocation size, in profile positions. (0, or >=1 && <=100000)
 *
 * Returns:   ptr to new <H4_FILTERMX>
 *
 * Throws:    <NULL> on allocation failure.
 */
H4_FILTERMX *
h4_filtermx_Create(int allocM)
{
  H4_FILTERMX *fx = NULL;
  int          status;

  /* Contract checks / argument validation */
  ESL_DASSERT1( (allocM >= 0 && allocM <= 100000) ); // 0 is ok here; gives you a shell.

  ESL_ALLOC(fx, sizeof(H4_FILTERMX));
  fx->M          = 0;
  fx->Vw         = 0;
  fx->dp         = NULL;
  fx->allocM     = 0;
#if eslDEBUGLEVEL > 0
  fx->dfp        = NULL;
#endif 
  
  if (allocM) {
    /* VF needs, in bytes:        3 MID cells *     maxQw vectors        * bytes/vector; (aligned on) */
    fx->dp     = esl_alloc_aligned(h4F_NCELLS * H4_Q(allocM, h4_VMAX_VF) * h4_VWIDTH,    h4_VALIGN);
    if (fx->dp == NULL) goto ERROR;
    fx->allocM = allocM;
  }
  return fx;

 ERROR:
  h4_filtermx_Destroy(fx);
  return NULL;
}


/* Function:  h4_filtermx_Reinit()
 * Synopsis:  Resize filter DP matrix for new profile size
 * Incept:    SRE, Sun 30 Jun 2019 (The Heavy, Short Change Hero)
 *
 * Purpose:   Given an existing filter matrix structure <fx>,
 *            and the dimension <allocM> of a new profile that 
 *            we're going to use (in consensus positions),
 *            assure that <fx> is large enough; reallocate 
 *            and reinitialize as needed.
 *
 *            <allocM> must be $\leq$ 100,000. This is a design
 *            limit on profile length.
 *           
 *            If <fx> has been set in dump mode for debugging, it's
 *            left that way. If caller wants to toggle dump mode, see
 *            <h4_filtermx_SetDumpMode()>.

 *            <h4_filtermx_Reinit(fx,M)> is essentially equivalent to
 *            <h47_filtermx_Create(M)>, while minimizing reallocation.
 *            It should only be called in a similar context as a
 *            <h4_filtermx_Create()> call. Existing data may be lost,
 *            because for efficiency, old data are not copied to new
 *            internal allocations.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. <fx> is unchanged
 *            and its data remain valid.
 *
 * Xref:      from H3's dp_vector/p7_filtermx.c
 */
int
h4_filtermx_Reinit(H4_FILTERMX *fx, int allocM)
{
  int status;

  /* Contract checks / argument validation */
  ESL_DASSERT1( (allocM >= 1 && allocM <= 100000) );

  /* Reinitialization to an empty structure */
  fx->M          = 0;
  fx->Vw         = 0;
  // Don't change debugging status (do_dumping, dfp).
  // _Reinit() gets called inside VF; caller needs to control dump mode.
  

  /* Reallocate if needed */
  if (allocM > fx->allocM) {
    free(fx->dp);
    /*                             3 MDI cells *   maxQw vectors/row      * bytes/vector; aligned */ 
    fx->dp     = esl_alloc_aligned(h4F_NCELLS * H4_Q(allocM, h4_VMAX_VF) * h4_VWIDTH,    h4_VALIGN);
    if (fx->dp == NULL) { status = eslEMEM; goto ERROR; }
    fx->allocM = allocM;
  }
  return eslOK;

 ERROR:
  fx->allocM = 0;
  return status;
}


/* Function:  h4_filtermx_GetScore()
 * Synopsis:  Get a score from striped DP memory.
 * Incept:    SRE, Sun Feb 12 10:46:28 2017 [The Chemical Brothers]
 *
 * Purpose:   Return the score at position <k>, state <s> in the striped
 *            vector DP row of <fx>. <k> is a model coordinate from (0).1..M.
 *            <s> is a state flag, one of <h4F_M | h4F_D | h4F_I>.
 *
 *            Slow, because it accesses the data in scalar form, and
 *            has to do a convoluted coord conversion from <k,s> to
 *            striped <q,z> coords. Only use this in noncritical code
 *            like debugging comparisons and dumps.
 *
 *            Independent of vector ISA, because all it needs to
 *            unstripe is the vector width <fx->Vw>. 
 *
 *            If k==0, returns 0, as a special case that we use to
 *            dump matrices in the same format as <h4_refmx_Dump()>.
 *
 * Args:      k : position (0).1..M
 *            s : state,   h4F_M | h4F_D | h4F_I
 *
 * Returns:   the score
 */
int
h4_filtermx_GetScore(const H4_FILTERMX *fx, int k, int s)
{
  int Q;
  int q, z;

  ESL_DASSERT1(( k >= 0 &&  k <= fx->allocM ));
  ESL_DASSERT1(( fx->Vw &&  !(fx->Vw & (fx->Vw-1)) ));     // Vw is set to a power of two.
  ESL_DASSERT1(( fx->M ));                                 // matrix has to have some data in it

  if (k == 0) return 0;

  Q = H4_Q(fx->M, fx->Vw);
  q = H4_Q_FROM_K(k,Q);
  z = H4_Z_FROM_K(k,Q);
  return fx->dp[ (q * h4F_NCELLS + s) * fx->Vw + z];
}


/* Function:  h4_filtermx_Reuse()
 * Synopsis:  Reuse memory of an existing H4_FILTERMX
 * Incept:    SRE, Sun 07 Jul 2019
 */
int
h4_filtermx_Reuse(H4_FILTERMX *fx)
{
  fx->M  = 0;
  fx->Vw = 0;
  return eslOK;
}


/* Function:  h4_filtermx_Destroy()
 * Synopsis:  Frees a H4_FILTERMX
 * Incept:    SRE, Sun 30 Jun 2019
 */
void
h4_filtermx_Destroy(H4_FILTERMX *fx)
{
  if (fx) {
    esl_alloc_free(fx->dp);
    free(fx);
  }
}



/*****************************************************************
 * 2. Debugging and development routines
 *****************************************************************/
#if eslDEBUGLEVEL > 0

/* Function:  h4_filtermx_SetDumpMode()
 * Synopsis:  Set an H4_FILTERMX so it gets dumps to a stream.
 *
 * Purpose:   Set the Viterbi filter matrix <fx> so VF dumps
 *            the matrix row by row to stream <dfp>.
 * 
 *            To turn dumping off, pass <NULL> for <dfp>.
 *            
 *            
 * Args:      fx        - DP matrix object to set for diagnostic dumping
 *            dfp       - open FILE * for diagnostic output, or NULL     
 * 
 * Returns:   <eslOK> on success.
 */
int
h4_filtermx_SetDumpMode(H4_FILTERMX *fx, FILE *dfp)
{
  fx->dfp        = dfp;
  return eslOK;
}

/* Function:  h4_filtermx_DumpRow()
 * Synopsis:  Dump current row of ViterbiFilter (int16) filter matrix.
 *
 * Purpose:   Dump current row of ViterbiFilter (int16) filter DP
 *            matrix <fx> for diagnostics, and include the values of
 *            specials <xE>, etc. The index <rowi> for the current row
 *            is used as a row label.
 *
 *            If <rowi> is 0, print a header first too.
 * 
 *            The output format is coordinated with <h4_refmx_Dump()> to
 *            facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
h4_filtermx_DumpRow(const H4_FILTERMX *fx, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC)
{
  int k;

  ESL_DASSERT1(( fx->Vw &&  !(fx->Vw & (fx->Vw-1)) ));     // V is set, to a power of two. bit trickery!
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
  for (k = 0; k <= fx->M; k++) fprintf(fx->dfp, "%6d ", h4_filtermx_GetScore(fx, k, h4F_M));
  
  /* Specials */
  fprintf(fx->dfp, "%6d %6d %6d %6d %6d\n", xE, xN, xJ, xB, xC);

  /* Insert scores. */
  fprintf(fx->dfp, "%4d I ", rowi);
  for (k = 0; k <= fx->M; k++) fprintf(fx->dfp, "%6d ", h4_filtermx_GetScore(fx, k, h4F_I));
  fprintf(fx->dfp, "\n");

  /* Delete scores. */
  fprintf(fx->dfp, "%4d D ", rowi);
  for (k = 0; k <= fx->M; k++) fprintf(fx->dfp, "%6d ", h4_filtermx_GetScore(fx, k, h4F_D));
  fprintf(fx->dfp, "\n\n");

  return eslOK;
}


#endif // eslDEBUGLEVEL > 0
