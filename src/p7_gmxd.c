/* P7_GMXD implementation: dynamic programming matrix for generic dual-mode 
 * (local/glocal) alignment.
 *
 * See also: generic_fwdback_dual.c
 * 
 * Contents:
 *   1. The <P7_GMXD> object
 *   2. Debugging aids
 *   3. Copyright and license information.
 */

#include "p7_config.h"

#include "hmmer.h"
#include "p7_gmxd.h"


/*****************************************************************
 * 1. The <P7_GMXD> object
 *****************************************************************/


/* Function:  p7_gmxd_Create()
 * Synopsis:  Create a new <P7_GMXD> DP matrix.
 *
 * Purpose:   Create a new <P7_GMXD> matrix for a model of
 *            length <M> and a target sequence of length
 *            <L>.
 *
 * Args:      M - model length
 *            L - target sequence length
 * 
 * Returns:   ptr to the new <P7_GMXD>.
 *
 * Throws:    <NULL> on any allocation failure.
 */
P7_GMXD *
p7_gmxd_Create(int M, int L)
{
  P7_GMXD *gxd = NULL;
  int      r;
  int      status;

  ESL_ALLOC(gxd, sizeof(P7_GMXD));
  gxd->dp_mem = NULL;
  gxd->dp     = NULL;

  gxd->allocR = L+1;
  gxd->allocW = (M+1) * p7GD_NSCELLS + p7GD_NXCELLS;
  gxd->allocN = (int64_t) gxd->allocR * (int64_t) gxd->allocW;

  ESL_ALLOC(gxd->dp_mem, sizeof(float  ) * gxd->allocN);
  ESL_ALLOC(gxd->dp,     sizeof(float *) * gxd->allocR);
  for (r = 0; r < gxd->allocR; r++)
    gxd->dp[r] = gxd->dp_mem + (r * gxd->allocW);

  gxd->validR = gxd->allocR;
  gxd->M      = 0;
  gxd->L      = 0;
  return gxd;

 ERROR:
  if (gxd) p7_gmxd_Destroy(gxd);
  return NULL;
}

/* Function:  p7_gmxd_GrowTo()
 * Synopsis:  Efficiently reallocate a <P7_GMXD>.
 *
 * Purpose:   Efficiently reallocate the matrix <gxd> to a new
 *            DP problem size, for model length <M> and target 
 *            sequence length <L>. Reuse the existing allocation
 *            as much as possible, to minimize reallocation calls.
 *
 * Args:      gxd - existing DP matrix
 *            M   - new model length
 *            L   - new target sequence length
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on memory allocation failure.
 */
int
p7_gmxd_GrowTo(P7_GMXD *gxd, int M, int L)
{
  int      W        = (M+1) * p7GD_NSCELLS + p7GD_NXCELLS;
  int      R        = L+1;
  uint64_t N        = (int64_t) R * (int64_t) W;
  int      do_reset = FALSE;
  int      r;
  int      status;

  /* are we already big enough? */
  if (W <= gxd->allocW && R <= gxd->validR) return eslOK;

  /* must we reallocate the matrix cells? */
  if (N > gxd->allocN)
    {
      ESL_REALLOC(gxd->dp_mem, sizeof(float) * N);
      gxd->allocN = N;
      do_reset    = TRUE;
    }
  
  /* must we reallocate the row pointers? */
  if (R > gxd->allocR)
    {
      ESL_REALLOC(gxd->dp, sizeof(float *) * R);
      gxd->allocR = R;
      do_reset    = TRUE;
    }

  /* must we widen the rows? */
  if (W > gxd->allocW) do_reset = TRUE;

  /* must we set some more valid row pointers? */
  if (R > gxd->validR) do_reset = TRUE;

  /* resize rows, reset valid row pointers */
  if (do_reset)
    {
      gxd->allocW = W;
      gxd->validR = ESL_MIN(gxd->allocR, (int) ( gxd->allocN / (uint64_t) gxd->allocW));
      for (r = 0; r < gxd->validR; r++)
	gxd->dp[r] = gxd->dp_mem + (r * gxd->allocW);
    }

  gxd->M = 0;
  gxd->L = 0;
  return eslOK;

 ERROR:
  return status;
}
	



/* Function:  p7_gmxd_Reuse()
 * Synopsis:  Finish using a <P7_GMXD> without deallocation.
 *
 * Purpose:   Caller says it is done using <gxd> for now, but is
 *            soon going to use it again for a new problem;
 *            so reinitialize, saving deallocation/reallocation.
 *            Equiv to <p7_gmxd_Destroy(); p7_gmxd_Create()> in 
 *            effect, but far faster.
 *
 * Args:      gxd - matrix that will be reused.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_gmxd_Reuse(P7_GMXD *gxd)
{
  gxd->M = 0;
  gxd->L = 0;
  return eslOK;
}


/* Function:  p7_gmxd_Destroy()
 * Synopsis:  Free a <P7_GMXD>.
 *
 * Purpose:   Free the matrix <gxd>.
 */
void
p7_gmxd_Destroy(P7_GMXD *gxd)
{
  if (!gxd) return;
  if (gxd->dp_mem) free(gxd->dp_mem);
  if (gxd->dp)     free(gxd->dp);
  free(gxd);
}


/*****************************************************************
 * 2. Debugging aids
 *****************************************************************/


/* Function:  p7_gmxd_DecodeSpecial()
 * Synopsis:  Convert special state code to string for debugging output.
 *
 * Purpose:   Given a special state code (such as <p7GD_E>), return a
 *            string label, suitable for state label in debugging output.
 *
 * Args:      type  - special state code, such as <p7GD_E>
 *
 * Returns:   ptr to a static string representation
 *
 * Throws:    (no abnormal error conditions)
 */
char *
p7_gmxd_DecodeSpecial(int type)
{
  switch (type) {
  case p7GD_E: return "E";
  case p7GD_N: return "N";
  case p7GD_J: return "J";
  case p7GD_B: return "B";
  case p7GD_L: return "L";
  case p7GD_G: return "G";
  case p7GD_C: return "C";
  }
  return NULL;
}

/* Function:  p7_gmxd_Dump()
 * Synopsis:  Dump a <P7_GMXD> for examination.
 *
 * Purpose:   Print the contents of <gxd> to stream <ofp>, for
 *            examination/debugging.
 */
int
p7_gmxd_Dump(FILE *ofp, P7_GMXD *gxd)
{
  return p7_gmxd_DumpWindow(ofp, gxd, 0, gxd->L, 0, gxd->M);
}

int
p7_gmxd_DumpWindow(FILE *ofp, P7_GMXD *gxd, int istart, int iend, int kstart, int kend)
{
  int   width     = 9;
  int   precision = 4;
  int   i,k,x;

  /* Header */
  fprintf(ofp, "     ");
  for (k = kstart; k <= kend;        k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0;      x < p7GD_NXCELLS; x++) fprintf(ofp, "%*s ", width, p7_gmxd_DecodeSpecial(x));
  fprintf(ofp, "\n");

  fprintf(ofp, "      ");
  for (k = kstart; k <= kend; k++)   fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < p7GD_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");

   /* DP matrix data */
  for (i = istart; i <= iend; i++)
    {
      fprintf(ofp, "%3d ML ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, gxd->dp[i][k * p7GD_NSCELLS + p7GD_ML]);
      for (x = 0;  x <  p7GD_NXCELLS; x++)  fprintf(ofp, "%*.*f ", width, precision, gxd->dp[i][ (gxd->M+1) * p7GD_NSCELLS + x]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IL ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, gxd->dp[i][k * p7GD_NSCELLS + p7GD_IL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DL ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, gxd->dp[i][k * p7GD_NSCELLS + p7GD_DL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d MG ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, gxd->dp[i][k * p7GD_NSCELLS + p7GD_MG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IG ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, gxd->dp[i][k * p7GD_NSCELLS + p7GD_IG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DG ", i);
      for (k = kstart; k <= kend;     k++)  fprintf(ofp, "%*.*f ", width, precision, gxd->dp[i][k * p7GD_NSCELLS + p7GD_DG]);
      fprintf(ofp, "\n\n");
  }
  return eslOK;
}


/* a private hack for making heatmaps */
int
p7_gmxd_DumpCSV(FILE *fp, P7_GMXD *pp, int istart, int iend, int kstart, int kend)
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
	    pp->dp[i][k * p7GD_NSCELLS + p7GD_ML] + pp->dp[i][k * p7GD_NSCELLS + p7GD_MG] + 
 	    pp->dp[i][k * p7GD_NSCELLS + p7GD_IL] + pp->dp[i][k * p7GD_NSCELLS + p7GD_IG] + 
	    pp->dp[i][k * p7GD_NSCELLS + p7GD_DL] + pp->dp[i][k * p7GD_NSCELLS + p7GD_DG];

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
