/* P7_REFMX implementation: dynamic programming matrix for reference 
 * implementation of dual-mode (local/glocal) alignment.
 *
 * See also: reference_fwdback.c
 * 
 * Contents:
 *   1. The <P7_REFMX> object
 *   2. Debugging aids
 *   3. Copyright and license information.
 */

#include "p7_config.h"

#include "hmmer.h"
#include "p7_refmx.h"


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
  int      r;
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

  rmx->validR = rmx->allocR;
  rmx->M      = 0;
  rmx->L      = 0;
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
  int      r;
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
	rmx->dp[r] = rmx->dp_mem + (r * rmx->allocW);
    }

  rmx->M = 0;
  rmx->L = 0;
  return eslOK;

 ERROR:
  return status;
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
  rmx->M = 0;
  rmx->L = 0;
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
  case p7R_E: return "E";
  case p7R_N: return "N";
  case p7R_J: return "J";
  case p7R_B: return "B";
  case p7R_L: return "L";
  case p7R_G: return "G";
  case p7R_C: return "C";
  }
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
  fprintf(ofp, "     ");
  for (k = kstart; k <= kend;        k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0;      x < p7R_NXCELLS; x++) fprintf(ofp, "%*s ", width, p7_refmx_DecodeSpecial(x));
  fprintf(ofp, "\n");

  fprintf(ofp, "      ");
  for (k = kstart; k <= kend; k++)   fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < p7R_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");

   /* DP matrix data */
  for (i = istart; i <= iend; i++)
    {
      fprintf(ofp, "%3d ML ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_ML]);
      for (x = 0;  x <  p7R_NXCELLS; x++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][ (rmx->M+1) * p7R_NSCELLS + x]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IL ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_IL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DL ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_DL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d MG ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_MG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IG ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_IG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DG ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_DG]);
      fprintf(ofp, "\n\n");
  }
  return eslOK;
}


/* a private hack for making heatmaps */
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

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
