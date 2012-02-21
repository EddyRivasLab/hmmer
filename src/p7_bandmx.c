/* Creation, manipulation of a P7_BANDMX object -
 * banded DP matrices for banded Forward/Backward, posterior decoding,
 * and optimal accuracy alignment, with a dual-mode local/glocal model.
 *
 * Contents:
 *    1. The P7_BANDMX object.
 *    x. Copyright and license information.
 */
#include "p7_config.h"

#include "easel.h"

#include "hmmer.h"
#include "p7_gbands.h"
#include "p7_bandmx.h"


/*****************************************************************
 * 1. The P7_BANDMX object
 *****************************************************************/

/* Function:  p7_bandmx_Create() 
 * Synopsis:  Creates banded DP matrix object, given bands.
 *
 * Purpose:   Create a new banded DP matrix object, given the bands
 *            specified in <bnd>. Return the pointer to the new
 *            <P7_BANDMX> object.
 *            
 *            The object contains two DP matrices, which are
 *            sufficient for all stages of processing: Forward,
 *            Backward, Decoding, and Alignment.
 *
 *            The <P7_BANDMX> object takes a copy of the <bnd>
 *            pointer. The caller remains responsible for free'ing
 *            <bnd>. If the caller subsequently alters <bnd>
 *            (including calculating new bands, either for the same or
 *            a new sequence), the caller must also reinitialize using
 *            <p7_bandmx_Reinit()>. 
 *
 * Args:      bnd - bands to use
 *
 * Returns:   ptr to new <P7_BANDMX> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_BANDMX *
p7_bandmx_Create(P7_GBANDS *bnd)
{
  P7_BANDMX *bmx = NULL;
  int        status;

  ESL_ALLOC(bmx, sizeof(P7_BANDMX));
  bmx->dp1    = NULL;
  bmx->dp2    = NULL;
  bmx->xmx1   = NULL;
  bmx->xmx2   = NULL;
  bmx->bnd    = bnd;
  bmx->dalloc = 0;
  bmx->xalloc = 0;

  ESL_ALLOC(bmx->dp1,  sizeof(float) * bnd->ncell * p7B_NSCELLS); /* i.e. *6, for [ ML MG IL IG DL DG ] */
  ESL_ALLOC(bmx->dp2,  sizeof(float) * bnd->ncell * p7B_NSCELLS); 
  ESL_ALLOC(bmx->xmx1, sizeof(float) * bnd->nrow  * p7B_NXCELLS); /* i.e. *7, for [ E N J B L G C ]     */
  ESL_ALLOC(bmx->xmx2, sizeof(float) * bnd->nrow  * p7B_NXCELLS); 
  bmx->dalloc = bnd->ncell;
  bmx->xalloc = bnd->nrow;
  return bmx;

 ERROR:
  p7_bandmx_Destroy(bmx);
  return NULL;
}


/* Function:  p7_bandmx_Sizeof()
 * Synopsis:  Returns the allocated size of a <P7_BANDMX>
 *
 * Purpose:   Returns the allocated size of <P7_BANDMX> <bmx>.
 *            
 *            This includes DP matrix space allocated for the
 *            <P7_BANDMX> but does not include the size of the
 *            associated <P7_GBANDS> structure.
 */
size_t
p7_bandmx_Sizeof(P7_BANDMX *bmx)
{
  size_t n = sizeof(P7_BANDMX);
  n += 2 * bmx->dalloc * p7B_NSCELLS * sizeof(float);
  n += 2 * bmx->xalloc * p7B_NXCELLS * sizeof(float);
  return n;
}


/* Function:  p7_bandmx_GrowTo()
 * Synopsis:  Resize a <P7_BANDMX> with new bands.
 *
 * Purpose:   Reallocate an existing <P7_BANDMX> <bmx> as needed,
 *            to make sure it can accommodate the bands in <bnd>.
 *            Essentially equivalent to <p7_bandmx_Destroy(); 
 *            p7_bandmx_Create()> while saving as much <free()/malloc()>
 *            activity as possible.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_bandmx_GrowTo(P7_BANDMX *bmx, P7_GBANDS *bnd)
{
  int status;

  if (bnd->ncell > bmx->dalloc) {
    ESL_REALLOC(bmx->dp1,  sizeof(float) * bnd->ncell * p7B_NSCELLS); 
    ESL_REALLOC(bmx->dp2,  sizeof(float) * bnd->ncell * p7B_NSCELLS); 
    bmx->dalloc = bnd->ncell;
  }
  if (bnd->nrow  > bmx->xalloc) {
    ESL_REALLOC(bmx->xmx1, sizeof(float) * bnd->nrow  * p7B_NXCELLS); 
    ESL_REALLOC(bmx->xmx2, sizeof(float) * bnd->nrow  * p7B_NXCELLS); 
    bmx->xalloc = bnd->nrow;
  }
  bmx->bnd = bnd;
  return eslOK;

 ERROR:
  return status;
}
  

/* Function:  p7_bandmx_Reuse()
 * Synopsis:  Reuse a <P7_BANDMX> instead of free'ing it.
 *
 * Purpose:   Reuse the banded DP matrix <bmx>, instead of free'ing
 *            it. Caller will also need to call <p7_bandmx_GrowTo()>
 *            to reinitialize it for a new set of bands for a new
 *            target sequence.
 */
int
p7_bandmx_Reuse(P7_BANDMX *bmx)
{
  bmx->bnd = NULL;
  return eslOK;
}

/* Function:  p7_bandmx_Destroy()
 * Synopsis:  Frees a <P7_BANDMX>.
 *
 * Purpose:   Free the banded DP matrix <bmx>. 
 * 
 *            The caller remains responsible for freeing the set of
 *            bands (a <P7_GBANDS> object) that it used to create
 *            the DP matrix.
 */
void
p7_bandmx_Destroy(P7_BANDMX *bmx)
{
  if (bmx)
    {
      if (bmx->dp1)  free(bmx->dp1);
      if (bmx->dp2)  free(bmx->dp2);
      if (bmx->xmx1) free(bmx->xmx1);
      if (bmx->xmx2) free(bmx->xmx2);
      /* gxb->bnd is a reference ptr copy; memory remains caller's responsibility */
      free(bmx);
    }
}


/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
