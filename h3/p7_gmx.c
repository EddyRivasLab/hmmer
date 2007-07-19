/* P7_GMX implementation: a generic dynamic programming matrix
 * 
 * SRE, Tue Jan 30 11:14:11 2007 [Einstein's, in St. Louis]
 * SVN $Id$
 */

#include "p7_config.h"
#include "hmmer.h"

/* Function:  p7_gmx_Create()
 * Incept:    SRE, Tue Jan 30 11:20:33 2007 [Einstein's, in St. Louis]
 *
 * Purpose:   Allocate a reusable, resizeable <P7_GMX> for models up to
 *            size <allocM> and sequences up to length <allocL>.
 *            
 *            We've set this up so it should be easy to allocate
 *            aligned memory, though we're not doing this yet.
 *
 * Returns:   a pointer to the new <P7_GMX>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_GMX *
p7_gmx_Create(int allocM, int allocL)
{
  int     status;
  P7_GMX *gx = NULL;
  int     i;

  /* level 1: the structure itself */
  ESL_ALLOC(gx, sizeof(P7_GMX));
  gx->dp  = NULL;
  gx->xmx = NULL;
  gx->xmx_mem = gx->dp_mem = NULL;

  /* level 2: row pointers, 0.1..L */
  ESL_ALLOC(gx->dp,  sizeof(float *) * (allocL+1));
 
  /* level 3: dp cell memory */
  ESL_ALLOC(gx->dp_mem,  sizeof(float)   * (allocL+1) * (allocM+1) * p7G_NSCELLS);
  ESL_ALLOC(gx->xmx_mem, sizeof(float)   * (allocL+1) * p7G_NXCELLS);

  /* Set the row pointers. */
  gx->xmx = gx->xmx_mem;
  for (i = 0; i <= allocL; i++) 
    gx->dp[i] = gx->dp_mem + i * (allocM+1) * p7G_NSCELLS;

  /* Initialize memory that's allocated but unused, only to keep
   * valgrind and friends happy.
   */
  for (i = 0; i <= allocL; i++) 
    {
      gx->dp[i][0      * p7G_NSCELLS + p7G_M] = -eslINFINITY; /* M_0 */
      gx->dp[i][0      * p7G_NSCELLS + p7G_I] = -eslINFINITY; /* I_0 */      
      gx->dp[i][0      * p7G_NSCELLS + p7G_D] = -eslINFINITY; /* D_0 */
      gx->dp[i][1      * p7G_NSCELLS + p7G_D] = -eslINFINITY; /* D_1 */
      gx->dp[i][allocM * p7G_NSCELLS + p7G_I] = -eslINFINITY; /* I_M */
    }

  gx->M      = allocM;
  gx->L      = allocL;
  gx->ncells = (allocM+1)*(allocL+1);
  gx->nrows  = (allocL+1);
  return gx;

 ERROR:
  if (gx != NULL) p7_gmx_Destroy(gx);
  return NULL;
}

/* Function:  p7_gmx_GrowTo()
 * Incept:    SRE, Tue Jan 30 11:31:23 2007 [Olin Library, St. Louis]
 *
 * Purpose:   Assures that a DP matrix <gx> is allocated
 *            for a model of size up to <allocM> and a sequence of
 *            length up to <allocL>; reallocates if necessary.
 *            
 *            This function does not respect the configured
 *            <RAMLIMIT>; it will allocate what it's told to
 *            allocate. 
 *
 * Returns:   <eslOK> on success, and <gx> may be reallocated upon
 *            return; any data that may have been in <gx> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <gx> must be assumed to be invalidated.
 */
int
p7_gmx_GrowTo(P7_GMX *gx, int allocM, int allocL)
{
  int     status;
  void   *p;
  int     i;
  size_t  ncells, nrows;

  ncells = (allocM+1)*(allocL+1);
  nrows  = allocL+1;
  if (ncells <= gx->ncells && nrows <= gx->nrows) { gx->M = allocM; gx->L = allocL; return eslOK; }

  /* must we realloc the 2D matrices? (or can we get away with just
   * jiggering the row pointers, if we are growing in one dimension
   * while shrinking in another?)
   */
  if (ncells > gx->ncells) 
    {
      ESL_RALLOC(gx->dp_mem, p, sizeof(float) * ncells * p7G_NSCELLS);
      gx->ncells = ncells;
    }

  /* must we reallocate the row pointers? */
  if (nrows > gx->nrows)
    {
      ESL_RALLOC(gx->xmx_mem, p, sizeof(float)   * nrows * p7G_NXCELLS);
      ESL_RALLOC(gx->dp,      p, sizeof(float *) * nrows);
      gx->nrows  = nrows;
    }

  /* reset all the row pointers (even though we may not have to: if we
   * increased allocL without reallocating cells or changing M,
   * there's some wasted cycles here - but I'm not sure this case
   * arises, and I'm not sure it's worth thinking hard about. )
   */
  gx->xmx = gx->xmx_mem;
  for (i = 0; i <= allocL; i++) 
    gx->dp[i] = gx->dp_mem + i * (allocM+1) * p7G_NSCELLS;

  gx->M      = allocM;
  gx->L      = allocL;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_gmx_Destroy()
 * Synopsis:  Frees a DP matrix.
 * Incept:    SRE, Tue Jan 30 11:17:36 2007 [Einstein's, in St. Louis]
 *
 * Purpose:   Frees a <P7_GMX>.
 *
 * Returns:   (void)
 */
void
p7_gmx_Destroy(P7_GMX *gx)
{
  if (gx == NULL) return;

  if (gx->dp      != NULL)  free(gx->dp);
  if (gx->dp_mem  != NULL)  free(gx->dp_mem);
  if (gx->xmx_mem != NULL)  free(gx->xmx_mem);
  free(gx);
  return;
}

/* Function:  p7_gmx_Dump()
 * Synopsis:  Dump a DP matrix to a stream, for diagnostics.
 * Incept:    SRE, Fri Jul 13 09:56:04 2007 [Janelia]
 *
 * Purpose:   Dump matrix <gx> to stream <fp> for diagnostics.
 */
int
p7_gmx_Dump(FILE *ofp, P7_GMX *gx)
{
  int i, k, x;

  /* Header */
  fprintf(ofp, "     ");
  for (k = 0; k <= gx->M;  k++) fprintf(ofp, "%8d ", k);
  fprintf(ofp, "%8s ", "E");
  fprintf(ofp, "%8s ", "N");
  fprintf(ofp, "%8s ", "J");
  fprintf(ofp, "%8s ", "B");
  fprintf(ofp, "%8s\n","C");
  for (k = 0; k <= gx->M+5; k++) fprintf(ofp, "%8s ", "--------");
  fprintf(ofp, "\n");
  
  /* DP matrix data */
  for (i = 0; i <= gx->L; i++)
    {
      fprintf(ofp, "%3d M ", i);
      for (k = 0; k <= gx->M;       k++) fprintf(ofp, "%8.4f ", gx->dp[i][k * p7G_NSCELLS + p7G_M]);
      for (x = 0; x <  p7G_NXCELLS; x++) fprintf(ofp, "%8.4f ", gx->xmx[  i * p7G_NXCELLS + x]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d I ", i);
      for (k = 0; k <= gx->M;       k++) fprintf(ofp, "%8.4f ", gx->dp[i][k * p7G_NSCELLS + p7G_I]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d D ", i);
      for (k = 0; k <= gx->M;       k++) fprintf(ofp, "%8.4f ", gx->dp[i][k * p7G_NSCELLS + p7G_D]);
      fprintf(ofp, "\n\n");
    }
  return eslOK;
}




/*****************************************************************
 * @LICENSE@
 *****************************************************************/
