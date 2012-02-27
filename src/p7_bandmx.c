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
 *            If <bnd> is <NULL>, allocate for a generic default size;
 *            caller will call <p7_bandmx_GrowTo()> when it has its
 *            bands ready. This allows the caller to create one
 *            matrix that it will use later in a loop over lots of 
 *            sequences.
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
 *            <p7_bandmx_GrowTo()>. 
 *
 * Args:      bnd - bands to use; or <NULL> if creating when bands are still unknown 
 *
 * Returns:   ptr to new <P7_BANDMX> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_BANDMX *
p7_bandmx_Create(P7_GBANDS *bnd)
{
  P7_BANDMX *bmx           = NULL;
  int64_t    default_ncell = 4096;
  int        default_nrow  = 256;
  int        status;

  ESL_ALLOC(bmx, sizeof(P7_BANDMX));
  bmx->dp1    = NULL;
  bmx->dp2    = NULL;
  bmx->xmx1   = NULL;
  bmx->xmx2   = NULL;
  bmx->bnd    = bnd;
  bmx->dalloc = (bnd ? bnd->ncell : default_ncell);
  bmx->xalloc = (bnd ? bnd->nrow  : default_nrow);

  ESL_ALLOC(bmx->dp1,  sizeof(float) * bmx->dalloc * p7B_NSCELLS); /* i.e. *6, for [ ML MG IL IG DL DG ] */
  ESL_ALLOC(bmx->dp2,  sizeof(float) * bmx->dalloc * p7B_NSCELLS); 
  ESL_ALLOC(bmx->xmx1, sizeof(float) * bmx->xalloc * p7B_NXCELLS); /* i.e. *7, for [ E N J B L G C ]     */
  ESL_ALLOC(bmx->xmx2, sizeof(float) * bmx->xalloc * p7B_NXCELLS); 
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

/****************************************************************
 * x. Debugging and development tools
 ****************************************************************/

char *
p7_bandmx_DecodeSpecial(int type)
{
  switch (type) {
  case p7B_E: return "E";
  case p7B_N: return "N";
  case p7B_J: return "J";
  case p7B_B: return "B";
  case p7B_L: return "L";
  case p7B_G: return "G";
  case p7B_C: return "C";
  default:    break;
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_BANDMX special state type code %d\n", type);
  return NULL;
}

static void
print_val(FILE *ofp, float val, int width, int precision, int which)
{
  if (which == p7B_DECODING) val = log(val);
  fprintf(ofp, "%*.*f ", width, precision, val);
}


/* Function:  p7_bandmx_Dump()
 * Synopsis:  Dumps a banded dual DP matrix for examination, debugging
 *
 * Purpose:   Dump one of the matrices in <bmx> to the stream <ofp> for
 *            examination. Which matrix you want is specified by the
 *            <which> argument, which is one of <p7B_FORWARD>, 
 *            <p7B_BACKWARD>, <p7B_DECODING>, <p7B_ALIGN>. 
 */
int
p7_bandmx_Dump(FILE *ofp, P7_BANDMX *bmx, int which)
{
  int   *bnd_ip = bmx->bnd->imem;  	/* array of ia..ib pairs for 0..nseg-1 segments: [ia0 ib0][ia1 ib1]... */
  int   *bnd_kp = bmx->bnd->kmem;	/* array of ka..kb pairs for each banded row i in segments */
  int    M      = bmx->bnd->M;
  float *dp;  
  float *xp;  
  int    g, i, k, x;
  int    ia, ib;
  int    ka, kb;
  int    width     = 9;
  int    precision = 4;

  /* Header; and initialization of dp, xp */
  switch (which) {
  case p7B_FORWARD:  fputs("Forward banded dual matrix dump.\n",       ofp); dp = bmx->dp1; xp = bmx->xmx1; break;
  case p7B_BACKWARD: fputs("Backward banded dual matrix dump.\n",      ofp); dp = bmx->dp2; xp = bmx->xmx2; break;
  case p7B_DECODING: fputs("Decoding banded dual matrix dump.\n",      ofp); dp = bmx->dp2; xp = bmx->xmx2; break;
  case p7B_ALIGN:    fputs("OA Alignment banded dual matrix dump.\n",  ofp); dp = bmx->dp1; xp = bmx->xmx1; break;
  default:           ESL_EXCEPTION(eslEINVAL, "no such matrix type");
  } 

  fputs("      ", ofp);
  for (k = 0; k <= M;          k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0; x < p7B_NXCELLS; x++) fprintf(ofp, "%*s ", width, p7_bandmx_DecodeSpecial(x));
  fputc('\n', ofp);

  fputs("       ", ofp);
  for (k = 0; k <= M+p7B_NXCELLS; k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fputc('\n', ofp);

  i = 0;
  for (g = 0; g < bmx->bnd->nseg; g++)
    {
      ia = *bnd_ip++;
      ib = *bnd_ip++;
      if (ia > i+1) fputs("...\n\n", ofp);
      for (i = ia; i <= ib; i++)
	{
	  ka = *bnd_kp++;
	  kb = *bnd_kp++;

	  fprintf(ofp, "%3d ML ", i);
	  for (k = 0; k <  ka; k++) fprintf  (ofp, "%*s ", width, ".....");
	  for (     ; k <= kb; k++) print_val(ofp, dp[(k-ka)*p7B_NSCELLS + p7B_ML],  width, precision, which);
	  for (     ; k <= M;  k++) fprintf  (ofp, "%*s ", width, ".....");

	  for (x = 0; x < p7B_NXCELLS; x++)   print_val(ofp, xp[x], width, precision, which);
	  fputc('\n', ofp);

	  fprintf(ofp, "%3d IL ", i);
	  for (k = 0; k <  ka; k++) fprintf  (ofp, "%*s ", width, ".....");
	  for (     ; k <= kb; k++) print_val(ofp, dp[(k-ka)*p7B_NSCELLS + p7B_IL], width, precision, which);
	  for (     ; k <= M;  k++) fprintf  (ofp, "%*s ", width, ".....");
	  fputc('\n', ofp);

	  fprintf(ofp, "%3d DL ", i);
	  for (k = 0; k <  ka; k++) fprintf  (ofp, "%*s ", width, ".....");
	  for (     ; k <= kb; k++) print_val(ofp, dp[(k-ka)*p7B_NSCELLS + p7B_DL], width, precision, which);
	  for (     ; k <= M;  k++) fprintf  (ofp, "%*s ", width, ".....");
	  fputc('\n', ofp);

	  fprintf(ofp, "%3d MG ", i);
	  for (k = 0; k <  ka; k++) fprintf  (ofp, "%*s ", width, ".....");
	  for (     ; k <= kb; k++) print_val(ofp, dp[(k-ka)*p7B_NSCELLS + p7B_MG], width, precision, which);
	  for (     ; k <= M;  k++) fprintf  (ofp, "%*s ", width, ".....");
	  fputc('\n', ofp);

	  fprintf(ofp, "%3d IG ", i);
	  for (k = 0; k <  ka; k++) fprintf  (ofp, "%*s ", width, ".....");
	  for (     ; k <= kb; k++) print_val(ofp, dp[(k-ka)*p7B_NSCELLS + p7B_IG], width, precision, which);
	  for (     ; k <= M;  k++) fprintf  (ofp, "%*s ", width, ".....");
	  fputc('\n', ofp);

	  fprintf(ofp, "%3d DG ", i);
	  for (k = 0; k <  ka; k++) fprintf  (ofp, "%*s ", width, ".....");
	  for (     ; k <= kb; k++) print_val(ofp, dp[(k-ka)*p7B_NSCELLS + p7B_DG], width, precision, which);
	  for (     ; k <= M;  k++) fprintf  (ofp, "%*s ", width, ".....");
	  fputs("\n\n", ofp);

	  dp += p7B_NSCELLS * (kb-ka+1);	/* skip ahead to next dp sparse "row" */
	  xp += p7B_NXCELLS;
	}
    }
  if (i <= bmx->bnd->L) fputs("...\n", ofp);
  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
