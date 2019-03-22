/* H4_REFMX: DP matrix for reference implementations
 *
 * Contents:
 *   1. The <H4_REFMX> object
 *   2. Debugging aids
 *   3. Validation 
 *
 * See also: 
 *   reference_viterbi.c
 */
#include "h4_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "h4_refmx.h"




/*****************************************************************
 * 1. The <H4_REFMX> object
 *****************************************************************/

/* Function:  h4_refmx_Create()
 * Synopsis:  Create a new <H4_REFMX> DP matrix.
 *
 * Purpose:   Create a new <H4_REFMX> matrix for a model of
 *            length <M> and a target sequence of length
 *            <L>.
 *
 * Args:      M - model length
 *            L - target sequence length
 * 
 * Returns:   ptr to the new <H4_REFMX>.
 *
 * Throws:    <NULL> on any allocation failure.
 */
H4_REFMX *
h4_refmx_Create(int M, int L)
{
  H4_REFMX *rmx = NULL;
  int      r,x;
  int      status;

  ESL_ALLOC(rmx, sizeof(H4_REFMX));
  rmx->dp_mem = NULL;
  rmx->dp     = NULL;

  rmx->allocR = L+1;
  rmx->allocW = (M+1) * h4R_NSCELLS + h4R_NXCELLS;
  rmx->allocN = (int64_t) rmx->allocR * (int64_t) rmx->allocW;

  ESL_ALLOC(rmx->dp_mem, sizeof(float  ) * rmx->allocN);
  ESL_ALLOC(rmx->dp,     sizeof(float *) * rmx->allocR);
  for (r = 0; r < rmx->allocR; r++)
    rmx->dp[r] = rmx->dp_mem + (r * rmx->allocW);

  /* Initialize all k=0 cells to -inf. */
  for (r = 0; r < rmx->allocR; r++)
    for (x = 0; x < h4R_NSCELLS; x++)
      rmx->dp[r][x] = -eslINFINITY;

  rmx->validR = rmx->allocR;
  rmx->M      = 0;
  rmx->L      = 0;
  rmx->type   = h4R_UNSET;
  return rmx;

 ERROR:
  h4_refmx_Destroy(rmx);
  return NULL;
}

/* Function:  h4_refmx_GrowTo()
 * Synopsis:  Efficiently reallocate a <H4_REFMX>.
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
h4_refmx_GrowTo(H4_REFMX *rmx, int M, int L)
{
  int      W        = (M+1) * h4R_NSCELLS + h4R_NXCELLS;
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
	  for (x = 0; x < h4R_NSCELLS; x++)  
	    rmx->dp[r][x] = -eslINFINITY;
	}
    }
  rmx->M = 0;
  rmx->L = 0;
  return eslOK;

 ERROR:
  return status;
}
	


/* Function:  h4_refmx_Reuse()
 * Synopsis:  Finish using a <H4_REFMX> without deallocation.
 *
 * Purpose:   Caller says it is done using <rmx> for now, but is
 *            soon going to use it again for a new problem;
 *            so reinitialize, saving deallocation/reallocation.
 *            Equiv to <h4_refmx_Destroy(); h4_refmx_Create()> in 
 *            effect, but faster.
 *
 * Args:      rmx - matrix that will be reused.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
h4_refmx_Reuse(H4_REFMX *rmx)
{
  rmx->M    = 0;
  rmx->L    = 0;
  rmx->type = h4R_UNSET;
  return eslOK;
}


/* Function:  h4_refmx_Destroy()
 * Synopsis:  Free a <H4_REFMX>.
 *
 * Purpose:   Free the matrix <rmx>.
 */
void
h4_refmx_Destroy(H4_REFMX *rmx)
{
  if (!rmx) return;
  if (rmx->dp_mem) free(rmx->dp_mem);
  if (rmx->dp)     free(rmx->dp);
  free(rmx);
}


/*****************************************************************
 * 2. Debugging and development tools
 *****************************************************************/

/* Function:  h4_refmx_DecodeSpecial()
 * Synopsis:  Convert special state code to string for debugging output.
 *
 * Purpose:   Given a special state code (such as <h4R_E>), return a
 *            string label suitable for state label in debugging output.
 *
 * Args:      type  - special state code, such as <h4R_E>
 *
 * Returns:   ptr to a static string representation
 *
 * Throws:    (no abnormal error conditions)
 */
char *
h4_refmx_DecodeSpecial(int type)
{
  switch (type) {
  case h4R_E:  return "E";
  case h4R_N:  return "N";
  case h4R_J:  return "J";
  case h4R_B:  return "B";
  case h4R_L:  return "L";
  case h4R_G:  return "G";
  case h4R_C:  return "C";
  case h4R_JJ: return "JJ";
  case h4R_CC: return "CC";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such H4_REFMX special state code %d\n", type);
  return NULL;
}

/* Function:  h4_refmx_DecodeState()
 * Synopsis:  Convert main state code to string for debugging output.
 *
 * Purpose:   Given a main state code (such as <h4R_MG>), return a
 *            string label, suitable as a state label in debugging
 *            output.
 *
 * Args:      type - main state code, such as <h4R_MG>
 *
 * Returns:   ptr to static string representation.
 */
char *
h4_refmx_DecodeState(int type)
{
  switch (type) {
  case h4R_ML: return "ML";
  case h4R_MG: return "MG";
  case h4R_IL: return "IL";
  case h4R_IG: return "IG";
  case h4R_DL: return "DL";
  case h4R_DG: return "DG";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such H4_REFMX main state code %d\n", type);
  return NULL;
}
    

/* Function:  h4_refmx_Dump(), h4_refmx_DumpWindow()
 * Synopsis:  Dump a <H4_REFMX> for examination.
 *
 * Purpose:   Print the contents of <rmx> to stream <ofp>, for
 *            examination/debugging.
 *            
 *            <h4_refmx_Dump()> prints the entire matrix.
 *            
 *            <h4_refmx_DumpWindow()> prints a chunk of it
 *            from <istart..iend>, <kstart..kend>.
 */
int
h4_refmx_Dump(FILE *ofp, H4_REFMX *rmx)
{
  return h4_refmx_DumpWindow(ofp, rmx, 0, rmx->L, 0, rmx->M);
}
int
h4_refmx_DumpWindow(FILE *ofp, H4_REFMX *rmx, int istart, int iend, int kstart, int kend)
{
  int   width     = 9;
  int   precision = 4;
  int   i,k,x;

  /* Header */
  fprintf(ofp, "       ");
  for (k = kstart; k <= kend;       k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0;      x < h4R_NXCELLS; x++) fprintf(ofp, "%*s ", width, h4_refmx_DecodeSpecial(x));
  fprintf(ofp, "\n");

  fprintf(ofp, "       ");
  for (k = kstart; k <= kend;  k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < h4R_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");

   /* DP matrix data */
  for (i = istart; i <= iend; i++)
    {
      fprintf(ofp, "%3d ML ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * h4R_NSCELLS + h4R_ML]);
      for (x = 0;  x <  h4R_NXCELLS; x++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][ (rmx->M+1) * h4R_NSCELLS + x]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d MG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * h4R_NSCELLS + h4R_MG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IL ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * h4R_NSCELLS + h4R_IL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * h4R_NSCELLS + h4R_IG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DL ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * h4R_NSCELLS + h4R_DL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * h4R_NSCELLS + h4R_DG]);
      fprintf(ofp, "\n\n");
  }
  return eslOK;
}
