/* P7_REFMX implementation: dynamic programming matrix for reference 
 * implementation of dual-mode (local/glocal) alignment.
 *
 * See also: reference_fwdback.c
 * 
 * Contents:
 *   1. The <P7_REFMX> object
 *   2. Debugging aids
 *   3. Validation 
 *   4. Copyright and license information.
 */

#include "p7_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_dmatrix.h"	/* ESL_DMATRIX is used for heat map visualization, in p7_refmx_PlotHeatMap() */
#include "esl_vectorops.h"

#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"


/*****************************************************************
 * 1. The <P7_REFMX> object
 *****************************************************************/

/* Function:  p7_refmx_Create()
 * Synopsis:  Create a new <P7_REFMX> DP matrix.
 *
 * Purpose:   Create a new <P7_REFMX> matrix for a model of
 *            length <M> and a target sequence of length
 *            <L>.
 *
 * Args:      M - model length
 *            L - target sequence length
 * 
 * Returns:   ptr to the new <P7_REFMX>.
 *
 * Throws:    <NULL> on any allocation failure.
 */
P7_REFMX *
p7_refmx_Create(int M, int L)
{
  P7_REFMX *rmx = NULL;
  int      r,x;
  int      status;

  ESL_ALLOC(rmx, sizeof(P7_REFMX));
  rmx->dp_mem = NULL;
  rmx->dp     = NULL;

  rmx->allocR = L+1;
  rmx->allocW = (M+1) * p7R_NSCELLS + p7R_NXCELLS;
  rmx->allocN = (int64_t) rmx->allocR * (int64_t) rmx->allocW;

  ESL_ALLOC(rmx->dp_mem, sizeof(float  ) * rmx->allocN);
  ESL_ALLOC(rmx->dp,     sizeof(float *) * rmx->allocR);
  for (r = 0; r < rmx->allocR; r++)
    rmx->dp[r] = rmx->dp_mem + (r * rmx->allocW);

  /* Initialize all k=0 cells to -inf. */
  for (r = 0; r < rmx->allocR; r++)
    for (x = 0; x < p7R_NSCELLS; x++)
      rmx->dp[r][x] = -eslINFINITY;

  rmx->validR = rmx->allocR;
  rmx->M      = 0;
  rmx->L      = 0;
  rmx->type   = p7R_UNSET;
  return rmx;

 ERROR:
  if (rmx) p7_refmx_Destroy(rmx);
  return NULL;
}

/* Function:  p7_refmx_GrowTo()
 * Synopsis:  Efficiently reallocate a <P7_REFMX>.
 *
 * Purpose:   Efficiently reallocate the matrix <rmx> to a new
 *            DP problem size, for model length <M> and target 
 *            sequence length <L>. Reuse the existing allocation
 *            as much as possible, to minimize reallocation calls.
 *
 * Args:      rmx - existing DP matrix
 *            M   - new model length
 *            L   - new target sequence length
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on memory allocation failure.
 */
int
p7_refmx_GrowTo(P7_REFMX *rmx, int M, int L)
{
  int      W        = (M+1) * p7R_NSCELLS + p7R_NXCELLS;
  int      R        = L+1;
  uint64_t N        = (int64_t) R * (int64_t) W;
  int      do_reset = FALSE;
  int      r,x;
  int      status;

  /* are we already big enough? */
  if (W <= rmx->allocW && R <= rmx->validR) return eslOK;

  /* must we reallocate the matrix cells? */
  if (N > rmx->allocN)
    {
      ESL_REALLOC(rmx->dp_mem, sizeof(float) * N);
      rmx->allocN = N;
      do_reset    = TRUE;
    }
  
  /* must we reallocate the row pointers? */
  if (R > rmx->allocR)
    {
      ESL_REALLOC(rmx->dp, sizeof(float *) * R);
      rmx->allocR = R;
      do_reset    = TRUE;
    }

  /* must we widen the rows? */
  if (W > rmx->allocW) do_reset = TRUE;

  /* must we set some more valid row pointers? */
  if (R > rmx->validR) do_reset = TRUE;

  /* resize rows, reset valid row pointers */
  if (do_reset)
    {
      rmx->allocW = W;
      rmx->validR = ESL_MIN(rmx->allocR, (int) ( rmx->allocN / (uint64_t) rmx->allocW));
      for (r = 0; r < rmx->validR; r++)
	{
	  rmx->dp[r] = rmx->dp_mem + (r * rmx->allocW);
	  for (x = 0; x < p7R_NSCELLS; x++)  
	    rmx->dp[r][x] = -eslINFINITY;
	}
    }
  rmx->M = 0;
  rmx->L = 0;
  return eslOK;

 ERROR:
  return status;
}
	
/* Function:  p7_refmx_Zero()
 * Synopsis:  Initialize a decoding matrix to all zeros.
 *
 * Purpose:   We're going to use <rmx> to count stochastic traces, to
 *            create an approximate posterior decoding matrix.  (We do
 *            this in unit testing of posterior decoding, for
 *            example.) First we need to set the whole matrix to
 *            zeros, before we start counting traces into it.  A newly
 *            allocated <rmx> only knows its allocation size, not the
 *            dimensions <M>,<L> of a DP comparison, so we need to
 *            pass that information too, for the comparison we're
 *            about to collect a stochastic trace ensemble from.
 *            
 *            Upon return, the <M> and <L> fields in the <rmx>
 *            structure have been set, the <type> field has been set
 *            to <p7R_DECODING>, and all DP cells are 0.0.
 *            
 * Args:      rmx  - posterior decoding matrix to zero
 *            M    - profile length
 *            L    - sequence length
 *
 * Returns:   <eslOK> on success.
 */
int
p7_refmx_Zero(P7_REFMX *rmx, int M, int L)
{
  int i;

  rmx->type = p7R_DECODING; 
  rmx->M    = M; 
  rmx->L    = L; 
  for (i = 0; i <= L; i++)
    esl_vec_FSet(rmx->dp[i], (rmx->M+1)*p7R_NSCELLS + p7R_NXCELLS, 0.0f);
  return eslOK;
}

int
p7_refmx_Copy(const P7_REFMX *src, P7_REFMX *dst)
{
  int i;
  for (i = 0; i <= src->L; i++)
    esl_vec_FCopy(src->dp[i], (src->M+1)*p7R_NSCELLS + p7R_NXCELLS, dst->dp[i]);
  dst->M    = src->M;
  dst->L    = src->L;
  dst->type = src->type;
  return eslOK;
}

int
p7_refmx_Rescale(P7_REFMX *rmx, float scale)
{
  int i;
  for (i = 0; i <= rmx->L; i++)
    esl_vec_FScale(rmx->dp[i], (rmx->M+1)*p7R_NSCELLS + p7R_NXCELLS, scale);
  return eslOK;

}

/* Function:  p7_refmx_Sizeof()
 * Synopsis:  Returns current allocated size of reference DP matrix, in bytes.
 */
size_t
p7_refmx_Sizeof(const P7_REFMX *rmx)
{
  size_t n = sizeof(P7_REFMX);
  n += rmx->allocN * sizeof(float);     /* dp_mem */
  n += rmx->allocR * sizeof(float *);	/* dp[]   */
  return n;
}

/* Function:  p7_refmx_MinSizeof()
 * Synopsis:  Returns minimum required size of reference DP matrix, in bytes.
 * 
 * Purpose:   Calculate and return the minimal required size of 
 *            a reference DP matrix for a comparison of a 
 *            profile of length <M> to a sequence of length <L>.
 *            
 *            This could overflow, on machines with 32b size_t. 
 */
size_t
p7_refmx_MinSizeof(int M, int L)
{
  size_t n = sizeof(P7_REFMX);
  n += sizeof(float)   * (size_t) (L+1) * (size_t) ( (M+1)*p7R_NSCELLS + p7R_NXCELLS); // dp_mem
  n += sizeof(float *) * (size_t) (L+1);                                               // dp[]
  return n;
}

/* Function:  p7_refmx_Reuse()
 * Synopsis:  Finish using a <P7_REFMX> without deallocation.
 *
 * Purpose:   Caller says it is done using <rmx> for now, but is
 *            soon going to use it again for a new problem;
 *            so reinitialize, saving deallocation/reallocation.
 *            Equiv to <p7_refmx_Destroy(); p7_refmx_Create()> in 
 *            effect, but faster.
 *
 * Args:      rmx - matrix that will be reused.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_refmx_Reuse(P7_REFMX *rmx)
{
  rmx->M    = 0;
  rmx->L    = 0;
  rmx->type = p7R_UNSET;
  return eslOK;
}


/* Function:  p7_refmx_Destroy()
 * Synopsis:  Free a <P7_REFMX>.
 *
 * Purpose:   Free the matrix <rmx>.
 */
void
p7_refmx_Destroy(P7_REFMX *rmx)
{
  if (!rmx) return;
  if (rmx->dp_mem) free(rmx->dp_mem);
  if (rmx->dp)     free(rmx->dp);
  free(rmx);
}


/*****************************************************************
 * 2. Debugging aids
 *****************************************************************/

/* Function:  p7_refmx_Compare()
 * Synopsis:  Compare two DP matrices for equality (within given tolerance)
 *
 * Purpose:   Compare all the values in DP matrices <rx1>, <rx2> for
 *            equality within absolute epsilon <tolerance>, using
 *            <esl_FCompareAbs()> calls. Return <eslOK> if all cell
 *            comparisons succeed; <eslFAIL> if not.
 *            
 *            Absolute difference comparison is preferred over
 *            relative differences. Numerical error accumulation in DP
 *            scales more with the number of terms than their
 *            magnitude. DP cells with values close to zero (and hence
 *            small absolute differences) may reasonably have large
 *            relative differences.
 */
int
p7_refmx_Compare(const P7_REFMX *rx1, const P7_REFMX *rx2, float tolerance)
{
  int i,k,s;
  int killmenow = FALSE;
#ifdef p7_DEBUGGING
  killmenow = TRUE;	/* Setting <killmenow> during debugging helps find where a comparison fails */
#endif /*p7_DEBUGGING*/
  
  if (rx1->M    != rx2->M)    { if (killmenow) abort(); return eslFAIL; }
  if (rx1->L    != rx2->L)    { if (killmenow) abort(); return eslFAIL; }
  if (rx1->type != rx2->type) { if (killmenow) abort(); return eslFAIL; }
  
  for (i = 0; i <= rx1->L; i++)
    {
      for (k = 0; k <= rx1->M; k++)   
	for (s = 0; s < p7R_NSCELLS; s++)
	  if ( esl_FCompareAbs(P7R_MX(rx1,i,k,s), P7R_MX(rx2,i,k,s), tolerance) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
      for (s = 0; s < p7R_NXCELLS; s++)
	if ( esl_FCompareAbs(P7R_XMX(rx1,i,s), P7R_XMX(rx2,i,s), tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
    }
  return eslOK;	
}

/* Function:  p7_refmx_CompareLocal()
 * Synopsis:  Compare two DP matrices (local paths only) for equality 
 *
 * Purpose:   A variant of <p7_refmx_Compare()> that compares only 
 *            cells on local paths. <MG,IG,DG,G> states are excluded.
 *            
 *            This gets used in unit tests that compare Backwards
 *            matrices computed by the ForwardFilter() with the
 *            reference implementation. You can't expect these Bck
 *            matrices to compare completely equal.  In the fwdfilter,
 *            glocal paths are all -inf by construction (they're not
 *            evaluated at all).  In the reference impl, in local
 *            mode, backwards values are finite along glocal paths
 *            (G/MG/IG/DG) because the zero transition that prohibits
 *            the path is the B->G transition, which isn't evaluated
 *            in the recursion until the *end* of glocal paths. 
 */
int
p7_refmx_CompareLocal(const P7_REFMX *rx1, const P7_REFMX *rx2, float tolerance)
{
  int i,k;
  int killmenow = FALSE;
#ifdef p7_DEBUGGING
  killmenow = TRUE;	/* Setting <killmenow> during debugging helps find where a comparison fails */
#endif /*p7_DEBUGGING*/
  
  if (rx1->M    != rx2->M)    { if (killmenow) abort(); return eslFAIL; }
  if (rx1->L    != rx2->L)    { if (killmenow) abort(); return eslFAIL; }
  if (rx1->type != rx2->type) { if (killmenow) abort(); return eslFAIL; }
  
  for (i = 0; i <= rx1->L; i++)
    {
      for (k = 0; k <= rx1->M; k++)   
	{
	  if ( esl_FCompareAbs(P7R_MX(rx1,i,k,p7R_ML), P7R_MX(rx2,i,k,p7R_ML), tolerance) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	  if ( esl_FCompareAbs(P7R_MX(rx1,i,k,p7R_IL), P7R_MX(rx2,i,k,p7R_IL), tolerance) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	  if ( esl_FCompareAbs(P7R_MX(rx1,i,k,p7R_DL), P7R_MX(rx2,i,k,p7R_DL), tolerance) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	}
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_E),  P7R_XMX(rx2,i,p7R_E),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_N),  P7R_XMX(rx2,i,p7R_N),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_J),  P7R_XMX(rx2,i,p7R_J),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_B),  P7R_XMX(rx2,i,p7R_B),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_L),  P7R_XMX(rx2,i,p7R_L),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_C),  P7R_XMX(rx2,i,p7R_C),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_JJ), P7R_XMX(rx2,i,p7R_JJ), tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_CC), P7R_XMX(rx2,i,p7R_CC), tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
    }
  return eslOK;	
}

/* Function:  p7_refmx_CompareDecoding()
 * Synopsis:  Compare exact, approximate posterior decoding matrices.
 *
 * Purpose:   Compare exact sparse decoding matrix <ppe> (calculated by 
 *            <p7_ReferenceDecoding()> to an approximate one <ppa> 
 *            (calculated by sampling lots of stochastic traces).
 *            Make sure that no nonzero value in <ppe> differs by
 *            more than absolute difference <tol> in <ppa>, and make
 *            sure that an exact zero value in <ppe> is also zero
 *            in <ppa>. Return <eslOK> if these tests succeed; <eslFAIL>
 *            if they do not.
 *            
 *            This comparison is used in a unit test of
 *            posterior decoding. See
 *            <reference_fwdback.c::utest_approx_decoding()>.
 *            
 * Args:      ppe  - exact posterior decoding matrix, from p7_ReferenceDecoding
 *            ppa  - approximate decoding mx, from stochastic trace ensemble
 *            tol  - absolute difference to tolerate per cell value
 *
 * Returns:   <eslOK> if all comparisons pass. <eslFAIL> if not.
 */
int
p7_refmx_CompareDecoding(const P7_REFMX *ppe, const P7_REFMX *ppa, float tol)
{
  float *dpe, *dpa;
  int i,k,s;
  int killmenow = FALSE;
#ifdef p7_DEBUGGING
  killmenow = TRUE;
#endif

  if (ppe->M    != ppa->M)       { if (killmenow) abort(); return eslFAIL; }
  if (ppe->L    != ppa->L)       { if (killmenow) abort(); return eslFAIL; }
  if (ppe->type != ppa->type)    { if (killmenow) abort(); return eslFAIL; }
  if (ppe->type != p7R_DECODING) { if (killmenow) abort(); return eslFAIL; }

  for (i = 0; i <= ppe->L; i++)	/* include i=0; DG0's have nonzero pp, as do N/B/L/G */
    {
      for (dpe = ppe->dp[i], dpa = ppa->dp[i], k = 0; k <= ppe->M; k++)
	for (s = 0; s < p7R_NSCELLS; s++, dpe++, dpa++)
	  if ( (*dpe == 0. && *dpa != 0.0) || esl_FCompareAbs(*dpe, *dpa, tol) != eslOK) 
	    { if (killmenow) abort(); return eslFAIL; }
      for (s = 0; s < p7R_NXCELLS; s++, dpe++, dpa++)
	if ( (*dpe == 0. && *dpa != 0.0) || esl_FCompareAbs(*dpe, *dpa, tol) != eslOK) 
	  { if (killmenow) abort(); return eslFAIL; }
    }
  return eslOK;
}

int
p7_refmx_CountTrace(const P7_TRACE *tr, P7_REFMX *rx)
{
  static int sttype[p7T_NSTATETYPES] = { -1, p7R_ML, p7R_MG, p7R_IL, p7R_IG, p7R_DL, p7R_DG, -1, p7R_N, p7R_B, p7R_L, p7R_G, p7R_E, p7R_C, p7R_J, -1 }; /* sttype[] translates trace idx to DP matrix idx*/
  int i = 0;
  int z;

  /* Contract checks || arg validation */
  ESL_DASSERT1( (rx->type == p7R_DECODING) );
  ESL_DASSERT1( (rx->L == tr->L) );
  ESL_DASSERT1( (rx->M == tr->M) );

  for (z = 1; z < tr->N-1; z++)	/* z=0 is S; z=N-1 is T; neither is represented in DP matrix, so avoid */
    {
      if (tr->i[z]) i = tr->i[z];

      if (tr->k[z]) 
	P7R_MX (rx, i, tr->k[z], sttype[(int)tr->st[z]]) += 1.0;
      else {
	P7R_XMX(rx, i,           sttype[(int)tr->st[z]]) += 1.0;
	if (tr->st[z] == p7T_C && tr->i[z]) P7R_XMX(rx, i, p7R_CC) += 1.0;
	if (tr->st[z] == p7T_J && tr->i[z]) P7R_XMX(rx, i, p7R_JJ) += 1.0;
      }	  
    }
  return eslOK;
}


/* Function:  p7_refmx_DecodeSpecial()
 * Synopsis:  Convert special state code to string for debugging output.
 *
 * Purpose:   Given a special state code (such as <p7R_E>), return a
 *            string label, suitable for state label in debugging output.
 *
 * Args:      type  - special state code, such as <p7R_E>
 *
 * Returns:   ptr to a static string representation
 *
 * Throws:    (no abnormal error conditions)
 */
char *
p7_refmx_DecodeSpecial(int type)
{
  switch (type) {
  case p7R_E:  return "E";
  case p7R_N:  return "N";
  case p7R_J:  return "J";
  case p7R_B:  return "B";
  case p7R_L:  return "L";
  case p7R_G:  return "G";
  case p7R_C:  return "C";
  case p7R_JJ: return "JJ";
  case p7R_CC: return "CC";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_REFMX special state code %d\n", type);
  return NULL;
}

/* Function:  p7_refmx_DecodeState()
 * Synopsis:  Convert main state code to string for debugging output.
 *
 * Purpose:   Given a main state code (such as <p7R_MG>), return a
 *            string label, suitable as a state label in debugging
 *            output.
 *
 * Args:      type - main state code, such as <p7R_MG>
 *
 * Returns:   ptr to static string representation.
 */
char *
p7_refmx_DecodeState(int type)
{
  switch (type) {
  case p7R_ML: return "ML";
  case p7R_MG: return "MG";
  case p7R_IL: return "IL";
  case p7R_IG: return "IG";
  case p7R_DL: return "DL";
  case p7R_DG: return "DG";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_REFMX main state code %d\n", type);
  return NULL;
}
    

/* Function:  p7_refmx_Dump()
 * Synopsis:  Dump a <P7_REFMX> for examination.
 *
 * Purpose:   Print the contents of <rmx> to stream <ofp>, for
 *            examination/debugging.
 */
int
p7_refmx_Dump(FILE *ofp, P7_REFMX *rmx)
{
  return p7_refmx_DumpWindow(ofp, rmx, 0, rmx->L, 0, rmx->M);
}

int
p7_refmx_DumpWindow(FILE *ofp, P7_REFMX *rmx, int istart, int iend, int kstart, int kend)
{
  int   width     = 9;
  int   precision = 4;
  int   i,k,x;

  /* Header */
  fprintf(ofp, "       ");
  for (k = kstart; k <= kend;       k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0;      x < p7R_NXCELLS; x++) fprintf(ofp, "%*s ", width, p7_refmx_DecodeSpecial(x));
  fprintf(ofp, "\n");

  fprintf(ofp, "       ");
  for (k = kstart; k <= kend;  k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < p7R_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");

   /* DP matrix data */
  for (i = istart; i <= iend; i++)
    {
      fprintf(ofp, "%3d ML ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_ML]);
      for (x = 0;  x <  p7R_NXCELLS; x++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][ (rmx->M+1) * p7R_NSCELLS + x]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IL ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_IL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DL ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_DL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d MG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_MG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_IG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_DG]);
      fprintf(ofp, "\n\n");
  }
  return eslOK;
}


/* a private hack for making heatmaps via a CSV file, using R, rather than
 * going via PlotHeatMap(). Using R from raw data, we have more flexibility
 * in formatting the plot.
 */
int
p7_refmx_DumpCSV(FILE *fp, P7_REFMX *pp, int istart, int iend, int kstart, int kend)
{
  int   width     = 7;
  int   precision = 5;
  int   i, k;
  float val;

  fprintf(fp, "i,");
  for (k = kend; k >= kstart; k--)
    fprintf(fp, "%-d%s", k, k==kstart ? "\n" : ",");

  for (i = istart; i <= iend; i++)
    {
      fprintf(fp, "%-d,", i);
      for (k = kend; k >= kstart; k--)
	{
	  val = 
	    pp->dp[i][k * p7R_NSCELLS + p7R_ML] + pp->dp[i][k * p7R_NSCELLS + p7R_MG] + 
 	    pp->dp[i][k * p7R_NSCELLS + p7R_IL] + pp->dp[i][k * p7R_NSCELLS + p7R_IG] + 
	    pp->dp[i][k * p7R_NSCELLS + p7R_DL] + pp->dp[i][k * p7R_NSCELLS + p7R_DG];

	  fprintf(fp, "%*.*f%s", width, precision, val, k==kstart ? "\n" : ", ");
	}
    }
  return eslOK;
}


/* Function:  p7_refmx_DumpBestDecoding()
 * Synopsis:  Dump best posterior decoding info on each seq position.
 *
 * Purpose:   Dump best posterior decoding info for each sequence position.
 *            For each row (sequence position) in posterior decoding matrix <rxd>,
 *            output a table with 8 columns to open stream <ofp>.
 *            
 *            <dsq>, <L>, and <gm> are optional, and are only used to
 *            obtain extra annotation for the table. They can be
 *            passed as <NULL>, <0>, <NULL>.
 *            
 *            The first two columns are just position i, and symbol at that 
 *            position. 
 *            
 *            3rd column is the decoded probability that i is
 *            homologous to the model: a sum over all emitting main
 *            states, as in <_PlotDomainInference()>.
 *            
 *            4th, 5th column: best (i,k) decoding, marginal over all
 *            main states at k, including D states; as in heat map
 *            functions. 4th col = best k; 5th col = best probability.
 *            
 *            6th, 7th, 8th column: best (i,k,s) decoding, for main
 *            model emitting states (ML, MG, IL, IG), as in
 *            <_PlotDomainInference()>. 6th = best k; 7th = best s (as
 *            string, i.e. MG); 8th = best probability.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_refmx_DumpBestDecoding(FILE *ofp, ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_REFMX *rxd)
{
  float *dpc;
  int   i, k;
  int   kmax, smax;
  float ppval, ppmax;

  /* Header 
   * <i> <sym> <p_homologous> <best_k> <best_pp> <best_k> <best_s> <best_pp>
   */
  fprintf(ofp, "# %12s %12s %15s\n", "", "(i,k)", "(i,k,s)");
  fprintf(ofp, "# %12s %12s %15s\n", "", "------------", "---------------");
  fprintf(ofp, "# %3s %1s %6s %5s %6s %5s %2s %6s\n",
	  "i", "X", "p(hom)", "bestk", "pp", "bestk", "s", "pp");
  fprintf(ofp, "# %3s %1s %6s %5s %6s %5s %2s %6s\n",
	  "---", "-", "------", "-----", "------", "-----", "--", "------");

  for (i = 1; i <= rxd->L; i++)
    {
      /* <i> <sym> <p_homologous>
       */
      fprintf(ofp, "%-5d %c %6.4f ",
	      i, 
	      ( (dsq && gm) ? gm->abc->sym[ dsq[i] ] : '?'),
	      1.0 - ( P7R_XMX(rxd, i, p7R_N) + P7R_XMX(rxd, i, p7R_JJ) + P7R_XMX(rxd, i, p7R_CC)));
 
      /* Best (i,k) decoding. D states included, as in heat maps.
       * <best_k> <best_pp>
       */
      for ( ppmax = -1., k = 1; k <= rxd->M; k++)
	{
	  ppval = 
	    rxd->dp[i][k * p7R_NSCELLS + p7R_ML] + rxd->dp[i][k * p7R_NSCELLS + p7R_MG] + 
	    rxd->dp[i][k * p7R_NSCELLS + p7R_IL] + rxd->dp[i][k * p7R_NSCELLS + p7R_IG] + 
	    rxd->dp[i][k * p7R_NSCELLS + p7R_DL] + rxd->dp[i][k * p7R_NSCELLS + p7R_DG];

	  if (ppval > ppmax) { ppmax = ppval;  kmax  = k; }
	}
      fprintf(ofp, "%-5d %6.4f ", kmax, ppmax);

      /* Best (i,k,s) decoding. Emitters only, as in PlotDomainInference().
       * <best_k> <best_s> <best_pp>
       */
      dpc = rxd->dp[i] + p7R_NSCELLS; /* on k=1; skip k=0 cells */
      for ( ppmax = -1., k = 1; k <= rxd->M; k++)
	{
	  if (dpc[p7R_ML] > ppmax) { ppmax = dpc[p7R_ML]; kmax = k; smax = p7R_ML; }
	  if (dpc[p7R_MG] > ppmax) { ppmax = dpc[p7R_MG]; kmax = k; smax = p7R_MG; }
	  if (dpc[p7R_IL] > ppmax) { ppmax = dpc[p7R_IL]; kmax = k; smax = p7R_IL; }
	  if (dpc[p7R_IG] > ppmax) { ppmax = dpc[p7R_IG]; kmax = k; smax = p7R_IG; }
	  dpc += p7R_NSCELLS;
	}
      fprintf(ofp, "%-5d %2s %6.4f\n", kmax, p7_refmx_DecodeState(smax), ppmax);
    }

  return eslOK;
}


/* Function:  p7_refmx_PlotHeatMap()
 * Synopsis:  Plot heat map representation of decoding matrix, PostScript.
 *
 * Purpose:   Plot a heat map representation of a window of the decoding
 *            matrix <pp> to open stream <ofp>, in PostScript format.
 *            Window is sequence position <istart>..<iend>, and model node
 *            <kstart>..<kend>. To plot the whole matrix, pass <istart=1>,
 *            <iend=pp->L>, <kstart=1>, <kend=pp->M>.
 *            
 *            One cell (i,k) is plotted for each sequence position i
 *            and model position k. The value in this cell is the sum
 *            of the posterior probabilities of all states s at
 *            (i,k,s): including M,D,I and including both glocal and
 *            local.
 *
 *            The format of the plot is hardcoded in several respects,
 *            just sufficient for some visualization purposes I had.
 *            The plot is always a single PostScript 8x11.5" page,
 *            792x612pt, with 20pt margins. The color scheme is a
 *            10-level scheme modified from ColorBrewer
 *            (colorbrewer2.org), their 9-class Red plus a light blue
 *            zero bin (pp $<$ 0.1). No coords are plotted on the axes
 *            of the matrix. All this is dictated by
 *            <esl_dmatrix_PlotHeatMap()>, and the documentation there
 *            may be more up to date than this.
 *            
 *            The sneakiest thing to watch for is that the sequence
 *            coord (vertical on the page) is numbered with 1 at the
 *            bottom, and L at the top; whereas the model coord axis
 *            1..M runs left to right horizontally. This is so I can
 *            rotate the plot to the right, and get plots that have
 *            the (usually longer) sequence running left to right
 *            across a slide or figure.
 *            
 *            Don't try to plot anything but a posterior probability
 *            matrix, because another thing that's hardcoded here is
 *            the expectation that DP matrix values range from 0 to 1.
 *
 * Args:      ofp    - open stream for writing PostScript output
 *            pp     - decoding matrix to visualize
 *            istart - seq coord start (1..L; perhaps 1)
 *            iend   - seq coord end   (1..L; perhaps pp->L)
 *            kstart - model coord start (1..M; perhaps 1)
 *            kend   - model coord end   (1..M; perhaps pp->M)
 *
 * Returns:   <eslOK > on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_refmx_PlotHeatMap(FILE *ofp, P7_REFMX *pp, int istart, int iend, int kstart, int kend)
{
  ESL_DMATRIX *dmx  = NULL;
  int          nrow = iend - istart + 1;
  int          ncol = kend - kstart + 1;
  int          i,k;
  float        val;
  int          status;

  /* Copy to an <ESL_DMATRIX>, which has heat mapping tools.
   * Marginalize local+glocal paths, MDI.
   */
  if ((dmx = esl_dmatrix_Create(nrow, ncol)) == NULL) { status = eslEMEM; goto ERROR; }
  for (i = istart; i <= iend; i++)
    for (k = kstart; k <= kend; k++)
      {
	val = 
	  pp->dp[i][k * p7R_NSCELLS + p7R_ML] + pp->dp[i][k * p7R_NSCELLS + p7R_MG] + 
	  pp->dp[i][k * p7R_NSCELLS + p7R_IL] + pp->dp[i][k * p7R_NSCELLS + p7R_IG] + 
	  pp->dp[i][k * p7R_NSCELLS + p7R_DL] + pp->dp[i][k * p7R_NSCELLS + p7R_DG];

	dmx->mx[iend-i][k-kstart] = val;
      }

  /* Plot it. */
  if ((status = esl_dmatrix_PlotHeatMap(ofp, dmx, 0.0, 1.0)) != eslOK) goto ERROR;

  esl_dmatrix_Destroy(dmx);
  return eslOK;

 ERROR:
  if (dmx) esl_dmatrix_Destroy(dmx);
  return status;
}

/* Function:  p7_refmx_PlotDomainInference()
 * Synopsis:  Plot graph of domain inference data, in XMGRACE xy format
 *
 * Purpose:   Given posterior decoding matrix <rxd>, plot a 2D graph of
 *            various quantities that could be useful for inferring 
 *            sequence coordinates of a domain location. Output is
 *            in XMGRACE xy format, and is sent to open writable
 *            stream <ofp>.
 *            
 *            Restrict the plot to sequence coords <ia..ib>, which can
 *            range from <1..rxd->L>. To get the full sequence, pass
 *            <ia=1>, <ib=rxd->L>.
 *            
 *            At least four datasets are plotted. Set 0 is
 *            P(homology): the posterior probability that this residue
 *            is 'in' the model (as opposed to being emitted by N,C,
 *            or J states).  Set 1 is P(B), the posterior probability
 *            of the BEGIN state. Set 2 is P(E), the posterior
 *            probability of the END state. Set 3 plots the maximum
 *            posterior probability of any emitting model state (M or
 *            I, glocal or local).
 *            
 *            Optionally, caller may also provide an indexed trace
 *            <tr>, showing where domains have been defined (as
 *            horizontal lines above the plot); each domain is an
 *            xmgrace dataset (so, sets 4..4+ndom-1 are these
 *            horizontal lines).
 */
int
p7_refmx_PlotDomainInference(FILE *ofp, const P7_REFMX *rxd, int ia, int ib, const P7_TRACE *tr)
{
  float *dpc;
  int    i,k,d;
  float  pmax;
  float  tr_height = 1.2;

  for (i = ia; i <= ib; i++)    /* Set #0 is P(homology) */
    fprintf(ofp, "%-6d %.5f\n", i, 1.0 - ( P7R_XMX(rxd, i, p7R_N) + P7R_XMX(rxd, i, p7R_JJ) + P7R_XMX(rxd, i, p7R_CC)));
  fprintf(ofp, "&\n");
  
  for (i = ia-1; i <= ib; i++)    /* Set #1 is P(B), the begin state. Start at ia-1, because P(B) is on start-1 for domain */
    fprintf(ofp, "%-6d %.5f\n", i, P7R_XMX(rxd, i, p7R_B));
  fprintf(ofp, "&\n");

  for (i = ia; i <= ib; i++)      /* Set #2 is P(E), the end state */
    fprintf(ofp, "%-6d %.5f\n", i, P7R_XMX(rxd, i, p7R_E));
  fprintf(ofp, "&\n");

  for (i = ia; i <= ib; i++)	  /* Set #3 is max P(i,k,s) for s=M,I states (G and L) */
    {
      dpc = rxd->dp[i] + p7R_NSCELLS; /* initialize at k=1 by skipping k=0 cells */
      for (pmax = 0., k = 1; k <= rxd->M; k++)
	{
	  pmax = ESL_MAX(pmax,
			 ESL_MAX( ESL_MAX( dpc[p7R_ML], dpc[p7R_MG] ),
				  ESL_MAX( dpc[p7R_IL], dpc[p7R_IG] )));
	  dpc += p7R_NSCELLS;
	}
      fprintf(ofp, "%-6d %.5f\n", i, pmax);
    }
  fprintf(ofp, "&\n");

  if (tr && tr->ndom)  /* Remaining sets are horiz lines, representing individual domains appearing in the optional <tr> */
    {
      for (d = 0; d < tr->ndom; d++)
	{
	  fprintf(ofp, "%-6d %.5f\n", tr->sqfrom[d], tr_height);
	  fprintf(ofp, "%-6d %.5f\n", tr->sqto[d],   tr_height);
	  fprintf(ofp, "&\n");
	}
    }
  return eslOK;
}



/*------------- end, (most) debugging tools ---------------------*/


/*****************************************************************
 * 3. Validation of matrices.
 *****************************************************************/
/* Gets its own section because it's a little fiddly and detailed.
 * See p7_refmx.h for notes on the patterns that get tested.
 */

static inline int
validate_dimensions(P7_REFMX *rmx, char *errbuf)
{
  if ( rmx->M <= 0)               ESL_FAIL(eslFAIL, errbuf, "nonpositive M");
  if ( rmx->L <= 0)               ESL_FAIL(eslFAIL, errbuf, "nonpositive L");
  if ( rmx->L+1    > rmx->validR) ESL_FAIL(eslFAIL, errbuf, "L+1 is larger than validR");
  if ( rmx->validR > rmx->allocR) ESL_FAIL(eslFAIL, errbuf, "validR larger than allocR");
  if ( (rmx->M+1)*p7R_NSCELLS+p7R_NXCELLS < rmx->allocW) ESL_FAIL(eslFAIL, errbuf, "M is too large for allocW");
  return eslOK;
}

static inline int 
validate_mainstate(P7_REFMX *rmx, int i, int k, int s, float val, char *errbuf)
{
  if (P7R_MX(rmx, i, k, s) != val) ESL_FAIL(eslFAIL, errbuf, "expected %f at i=%d, k=%d, %s", val, i, k, p7_refmx_DecodeState(s));
  return eslOK;
}

static inline int
validate_special(P7_REFMX *rmx, int i, int s, float val, char *errbuf)
{
  if (P7R_XMX(rmx, i, s) != val) ESL_FAIL(eslFAIL, errbuf, "expected %f at i=%d, %s", val, i, p7_refmx_DecodeSpecial(s));
  return eslOK;
}

static inline int
validate_no_nan(P7_REFMX *rmx, char *errbuf)
{
  int i,k,s;
  for (i = 0; i <= rmx->L; i++)
    for (k = 0; k <= rmx->M; k++)
      {
	for (s = 0; s < p7R_NSCELLS; s++)
	  if (isnan(P7R_MX(rmx, i, k, s))) ESL_FAIL(eslFAIL, errbuf, "found NaN at i=%d, k=%d, %s", i, k, p7_refmx_DecodeState(s));      
	for (s = 0; s < p7R_NXCELLS; s++)
	  if (isnan(P7R_XMX(rmx, i, s)))   ESL_FAIL(eslFAIL, errbuf, "found NaN at i=%d, %s", i, p7_refmx_DecodeSpecial(s));      
      }
  return eslOK;
}

static inline int
validate_column_zero(P7_REFMX *rmx, float val, char *errbuf)
{
  int i, s, status;
  for (i = 0; i <= rmx->L; i++)
    for (s = 0; s < p7R_NSCELLS; s++)
      if (( status = validate_mainstate(rmx, i, 0, s, val, errbuf)) != eslOK) return status;
  return eslOK;
}


static inline int
validate_forward(P7_REFMX *rmx, char *errbuf)
{
  int i,k,s;
  int status;

  if (( status = validate_column_zero(rmx, -eslINFINITY, errbuf)) != eslOK) return status;

  /* Row i=0 */
  for (k = 1; k <= rmx->M; k++)
    for (s = 0; s <= p7R_NSCELLS; s++)
      if ( (status = validate_mainstate(rmx, 0, k, s, -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_E,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_N,     0.0f,         errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_J,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_C,     -eslINFINITY, errbuf)) != eslOK) return status;

  /* Row i=1 has some additional expected -inf's, different from the remaining rows 2..L */
  for (k = 1; k <= rmx->M; k++) {
    if ( (status = validate_mainstate(rmx, 1, k, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
    if ( (status = validate_mainstate(rmx, 1, k, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
  }

  /* Rows 2..L, plus the row i=1 cells we didn't just check */
  for (i = 1; i <= rmx->L; i++)
    {
      if ((status =    validate_mainstate(rmx, i,      1, p7R_DL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i,      1, p7R_DG, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i, rmx->M, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i, rmx->M, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
    }

  /* CC/JJ are only for decoding; Forward matrix has them as -inf */
  for (i = 0; i <= rmx->L; i++)
    {
      if ( ( status = validate_special(rmx, 0, p7R_JJ,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(rmx, 0, p7R_CC,    -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}

static inline int
validate_backward(P7_REFMX *rmx, char *errbuf)
{
  int i,k,s;
  int status;

  if (( status = validate_column_zero(rmx, -eslINFINITY, errbuf)) != eslOK) return status;  

  /* Row i=0 */
  for (k = 1; k <= rmx->M; k++)
    for (s = 0; s < p7R_NSCELLS; s++)
      if ( (status = validate_mainstate(rmx, 0, k, s, -eslINFINITY, errbuf)) != eslOK) return status;

  /* Rows 1..L */
  for (i = 1; i <= rmx->L; i++)
    {
      if ((status = validate_mainstate(rmx, i, rmx->M, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status = validate_mainstate(rmx, i, rmx->M, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
    }

  /* Row i=L has some additional expected -inf's, different from the remaining rows */
  for (k = 1; k <= rmx->M; k++) {
    if ( (status = validate_mainstate(rmx, rmx->L, k, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
    if ( (status = validate_mainstate(rmx, rmx->L, k, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
  }
  if ( (status = validate_special  (rmx, rmx->L,    p7R_N,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_J,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_B,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_L,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_G,  -eslINFINITY, errbuf)) != eslOK) return status;

  /* CC/JJ are only for decoding; Backward matrix has them as -inf */
  for (i = 0; i <= rmx->L; i++)
    {
      if ( ( status = validate_special(rmx, 0, p7R_JJ,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(rmx, 0, p7R_CC,    -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}

static inline int
validate_decoding(P7_REFMX *rmx, char *errbuf)
{
  int i,k;
  int status;

  /* All of k=0 column is 0.0 */
  if (( status = validate_column_zero(rmx, 0.0, errbuf)) != eslOK) return status;  

  /* Row i=0: main states all 0.0 except DGk for k=1..M-1, which are reachable by unfolded wing-retracted entry */
  for (k = 1; k <= rmx->M; k++)
    {
      if ( (status = validate_mainstate(rmx, 0, k, p7R_ML, 0.0f, errbuf)) != eslOK) return status;
      if ( (status = validate_mainstate(rmx, 0, k, p7R_MG, 0.0f, errbuf)) != eslOK) return status;
      if ( (status = validate_mainstate(rmx, 0, k, p7R_IL, 0.0f, errbuf)) != eslOK) return status;
      if ( (status = validate_mainstate(rmx, 0, k, p7R_IG, 0.0f, errbuf)) != eslOK) return status;
      if ( (status = validate_mainstate(rmx, 0, k, p7R_DL, 0.0f, errbuf)) != eslOK) return status;
    }
  if ( (status = validate_mainstate(rmx, 0, rmx->M, p7R_DG, 0.0f, errbuf)) != eslOK) return status;
  
  /*          and E,J,C are unreached */
  if ( ( status = validate_special(rmx, 0, p7R_E, 0.0f, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_J, 0.0f, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_C, 0.0f, errbuf)) != eslOK) return status;

  /* Rows 1..L */
  for (i = 1; i <= rmx->L; i++)
    {
      if ((status = validate_mainstate(rmx, i,      1, p7R_DL, 0.0f, errbuf)) != eslOK) return status; /* DL1 doesn't exist. DG1 does, in decoding! */
      if ((status = validate_mainstate(rmx, i, rmx->M, p7R_IL, 0.0f, errbuf)) != eslOK) return status;
      if ((status = validate_mainstate(rmx, i, rmx->M, p7R_IG, 0.0f, errbuf)) != eslOK) return status;
    }
  /* on row L, DG1 cannot be reached */
  if ( (status = validate_mainstate(rmx, rmx->L,    1, p7R_DG, 0.0f, errbuf)) != eslOK) return status;

  /* Rows i=1,L have some additional zeros */
  for (k = 1; k < rmx->M; k++) {
    if ((status = validate_mainstate(rmx,      1, k, p7R_IL, 0.0f, errbuf)) != eslOK) return status;
    if ((status = validate_mainstate(rmx,      1, k, p7R_IG, 0.0f, errbuf)) != eslOK) return status;
    if ((status = validate_mainstate(rmx, rmx->L, k, p7R_IL, 0.0f, errbuf)) != eslOK) return status;
    if ((status = validate_mainstate(rmx, rmx->L, k, p7R_IG, 0.0f, errbuf)) != eslOK) return status;
  }
  /* on row 1, JJ and CC cannot be reached */
  if ((status = validate_special  (rmx,      1, p7R_JJ, 0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx,      1, p7R_CC, 0.0f, errbuf)) != eslOK) return status;
  /* on row L, only ...E->C(CC)->T specials can be reached. */
  if ((status = validate_special  (rmx, rmx->L, p7R_N,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx, rmx->L, p7R_J,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx, rmx->L, p7R_B,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx, rmx->L, p7R_L,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx, rmx->L, p7R_G,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx, rmx->L, p7R_JJ, 0.0f, errbuf)) != eslOK) return status;

  return eslOK;
}

static inline int
validate_alignment(P7_REFMX *rmx, char *errbuf)
{
  int i,k,s;
  int status;

  if (( status = validate_column_zero(rmx, -eslINFINITY, errbuf)) != eslOK) return status;  

  /* Row i=0 */
  for (k = 1; k <= rmx->M; k++)
    for (s = 0; s <= p7R_NSCELLS; s++)
      if ( (status = validate_mainstate(rmx, 0, k, s, -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_E,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_J,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_C,     -eslINFINITY, errbuf)) != eslOK) return status;

  /* Row i=1 has some additional expected -inf's, different from the remaining rows 2..L */
  for (k = 1; k <= rmx->M; k++) {
    if ( (status = validate_mainstate(rmx, 1, k, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
    if ( (status = validate_mainstate(rmx, 1, k, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
  }

  /* Rows 2..L, plus the row i=1 cells we didn't just check */
  for (i = 1; i <= rmx->L; i++)
    {
      if ((status =    validate_mainstate(rmx, i,      1, p7R_DL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i,      1, p7R_DG, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i, rmx->M, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i, rmx->M, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
    }

  /* Row i=L has some additional expected -inf's, different from the remaining rows */
  for (k = 1; k <= rmx->M; k++) {
    if ( (status = validate_mainstate(rmx, rmx->L, k, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
    if ( (status = validate_mainstate(rmx, rmx->L, k, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
  }
  if ( (status = validate_special  (rmx, rmx->L,    p7R_N,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_J,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_B,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_L,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_G,  -eslINFINITY, errbuf)) != eslOK) return status;

  /* and throughout, CC/JJ are only for decoding; Alignment matrix has them as -inf */
  for (i = 0; i <= rmx->L; i++)
    {
      if ( ( status = validate_special(rmx, 0, p7R_JJ,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(rmx, 0, p7R_CC,    -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}


/* Function:  p7_refmx_Validate()
 * Synopsis:  Validates a reference DP matrix.
 *
 * Purpose:   Validates the internals of the
 *            DP matrix <rmx>. Returns <eslOK> if
 *            it passes. Returns <eslFAIL> if it fails,
 *            and sets <errbuf> to contain an explanation,
 *            if caller provided an <errbuf>.
 *
 * Args:      rmx    - reference DP matrix to validate.
 *            errbuf - allocated space for err msg, or NULL
 *
 * Returns:   <eslOK> on success; <eslFAIL> on failure, and
 *            <errbuf> (if provided) contains explanation.
 *
 * Xref:      See notes in <p7_refmx.h> for an explanation
 *            of the particular patterns we're checking.
 */
int
p7_refmx_Validate(P7_REFMX *rmx, char *errbuf)
{
  int status;

  if ( (status = validate_dimensions(rmx, errbuf)) != eslOK) return status;
  if ( (status = validate_no_nan    (rmx, errbuf)) != eslOK) return status;
  
  switch (rmx->type) {
  case p7R_FORWARD:   if ( (status = validate_forward    (rmx,               errbuf)) != eslOK) return status;  break;
  case p7R_VITERBI:   if ( (status = validate_forward    (rmx,               errbuf)) != eslOK) return status;  break; // Viterbi has same pattern as Forward
  case p7R_BACKWARD:  if ( (status = validate_backward   (rmx,               errbuf)) != eslOK) return status;  break;
  case p7R_DECODING:  if ( (status = validate_decoding   (rmx,               errbuf)) != eslOK) return status;  break;     
  case p7R_ALIGNMENT: if ( (status = validate_alignment  (rmx,               errbuf)) != eslOK) return status;  break; 
  case p7R_UNSET:     if ( (status = validate_column_zero(rmx, -eslINFINITY, errbuf)) != eslOK) return status;  break;
  default:            ESL_FAIL(eslFAIL, errbuf, "no such reference DP algorithm type %d", rmx->type);
  }
  return eslOK;
}
/*------------------ end, validation ----------------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
