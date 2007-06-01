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

  ESL_ALLOC(gx, sizeof(P7_GMX));
  gx->xmx     = gx->mmx     = gx->imx     = gx->dmx     = NULL;
  gx->xmx_mem = gx->mmx_mem = gx->imx_mem = gx->dmx_mem = NULL;

  ESL_ALLOC(gx->xmx, sizeof(int *) * (allocL+1)); 
  ESL_ALLOC(gx->mmx, sizeof(int *) * (allocL+1)); 
  ESL_ALLOC(gx->imx, sizeof(int *) * (allocL+1)); 
  ESL_ALLOC(gx->dmx, sizeof(int *) * (allocL+1)); 
  
  ESL_ALLOC(gx->xmx_mem, 5 * (allocL+1) * sizeof(int));
  ESL_ALLOC(gx->mmx_mem, (allocM+1) * (allocL+1) * sizeof(int));
  ESL_ALLOC(gx->imx_mem, (allocM+1) * (allocL+1) * sizeof(int));
  ESL_ALLOC(gx->dmx_mem, (allocM+1) * (allocL+1) * sizeof(int));

  for (i = 0; i <= allocL; i++) {
    gx->xmx[i] = gx->xmx_mem + i * 5;
    gx->mmx[i] = gx->mmx_mem + i * (allocM+1);
    gx->imx[i] = gx->imx_mem + i * (allocM+1);
    gx->dmx[i] = gx->dmx_mem + i * (allocM+1);
  }

  gx->M      = 0;
  gx->L      = 0;
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
  if (ncells <= gx->ncells && nrows <= gx->nrows) return eslOK;

  /* must we realloc the 2D matrices? (or can we get away with just
   * jiggering the row pointers, if we are growing in one dimension
   * while shrinking in another?)
   */
  if (ncells > gx->ncells) {
    ESL_RALLOC(gx->mmx_mem, p, sizeof(int) * ncells);
    ESL_RALLOC(gx->imx_mem, p, sizeof(int) * ncells);
    ESL_RALLOC(gx->dmx_mem, p, sizeof(int) * ncells);
  }

  /* must we reallocate the row pointers? */
  if (nrows > gx->nrows)
    {
      ESL_RALLOC(gx->xmx_mem, p, sizeof(int *) * nrows * 5);
      ESL_RALLOC(gx->xmx,     p, sizeof(int *) * nrows);
      ESL_RALLOC(gx->mmx,     p, sizeof(int *) * nrows);
      ESL_RALLOC(gx->imx,     p, sizeof(int *) * nrows);
      ESL_RALLOC(gx->dmx,     p, sizeof(int *) * nrows);
    }

  /* reset all the row pointers (even though we may not have to: if we
   * increased allocL without reallocating cells or changing M,
   * there's some wasted cycles here - but I'm not sure this case
   * arises, and I'm not sure it's worth thinking hard about. )
   */
  for (i = 0; i <= allocL; i++) {
    gx->xmx[i] = gx->xmx_mem + i * 5;
    gx->mmx[i] = gx->mmx_mem + i * (allocM+1);
    gx->imx[i] = gx->imx_mem + i * (allocM+1);
    gx->dmx[i] = gx->dmx_mem + i * (allocM+1);
  }

  gx->M      = 0;
  gx->L      = 0;
  gx->ncells = ncells;
  gx->nrows  = nrows;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_gmx_Destroy()
 * Incept:    SRE, Tue Jan 30 11:17:36 2007 [Einstein's, in St. Louis]
 *
 * Purpose:   Frees a <P7_GMX>.
 *
 * Returns:   (void)
 */
void
p7_gmx_Destroy(P7_GMX *gx)
{
  if (gx != NULL) {
    if (gx->xmx     != NULL) free(gx->xmx);
    if (gx->mmx     != NULL) free(gx->mmx);
    if (gx->imx     != NULL) free(gx->imx);
    if (gx->dmx     != NULL) free(gx->dmx);
    if (gx->xmx_mem != NULL) free(gx->xmx_mem);
    if (gx->mmx_mem != NULL) free(gx->mmx_mem);
    if (gx->imx_mem != NULL) free(gx->imx_mem);
    if (gx->dmx_mem != NULL) free(gx->dmx_mem);
    free(gx);
  }
  return;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
