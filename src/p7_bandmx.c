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
 *            caller will call <p7_bandmx_Reinit()> when it has its
 *            bands ready. This allows the caller to create one
 *            matrix that it will use later in a loop over lots of 
 *            sequences.
 *
 *            The object contains one DP matrix. The caller will
 *            likely require up for four separate ones: Forward,
 *            Backward, Decoding, and Alignment. We may create a
 *            unified structure with all four matrices at some point,
 *            to streamline, once our logical flow is a little more
 *            fixed.
 *
 *            The <P7_BANDMX> object takes a copy of the <bnd>
 *            pointer. The caller remains responsible for free'ing
 *            <bnd>. If the caller subsequently alters <bnd>
 *            (including calculating new bands, either for the same or
 *            a new sequence), the caller must also reinitialize using
 *            <p7_bandmx_Reinit()>. 
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
  int        default_nx    = 256;
  int        status;

  ESL_ALLOC(bmx, sizeof(P7_BANDMX));
  bmx->dp     = NULL;
  bmx->xmx    = NULL;
  bmx->bnd    = bnd;
  bmx->dalloc = (bnd ? bnd->ncell            : default_ncell);
  bmx->xalloc = (bnd ? bnd->nrow + bnd->nseg : default_nx);

  ESL_ALLOC(bmx->dp,  sizeof(float) * bmx->dalloc * p7B_NSCELLS); /* i.e. *6, for [ ML MG IL IG DL DG ]  */
  ESL_ALLOC(bmx->xmx, sizeof(float) * bmx->xalloc * p7B_NXCELLS); /* i.e. *9, for [ E N J B L G C JJ CC] */
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
  n += bmx->dalloc * p7B_NSCELLS * sizeof(float);
  n += bmx->xalloc * p7B_NXCELLS * sizeof(float);
  return n;
}


/* Function:  p7_bandmx_Reinit()
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
p7_bandmx_Reinit(P7_BANDMX *bmx, P7_GBANDS *bnd)
{
  int status;

  if (bnd->ncell > bmx->dalloc) {
    ESL_REALLOC(bmx->dp,  sizeof(float) * bnd->ncell * p7B_NSCELLS); 
    bmx->dalloc = bnd->ncell;
  }
  if (bnd->nrow + bnd->nseg > bmx->xalloc) {
    ESL_REALLOC(bmx->xmx, sizeof(float) * (bnd->nrow+bnd->nseg) * p7B_NXCELLS); 
    bmx->xalloc = bnd->nrow + bnd->nseg;
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
 *            it. Caller will also need to call <p7_bandmx_Reinit()>
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
      if (bmx->dp)  free(bmx->dp);
      if (bmx->xmx) free(bmx->xmx);
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
  case p7B_E:  return "E";
  case p7B_N:  return "N";
  case p7B_J:  return "J";
  case p7B_B:  return "B";
  case p7B_L:  return "L";
  case p7B_G:  return "G";
  case p7B_C:  return "C";
  case p7B_JJ: return "JJ";
  case p7B_CC: return "CC";
  default:     break;
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_BANDMX special state type code %d\n", type);
  return NULL;
}



/* Function:  p7_bandmx_Dump()
 * Synopsis:  Dumps a banded dual DP matrix for examination, debugging
 *
 * Purpose:   Dump the banded DP matrix <bmx> to the stream <ofp> for
 *            examination. 
 */
int
p7_bandmx_Dump(FILE *ofp, P7_BANDMX *bmx)
{
  return p7_bandmx_DumpWindow(ofp, bmx, 0, bmx->bnd->L, 0, bmx->bnd->M);
}


/* Function:  p7_bandmx_DumpWindow()
 * Synopsis:  Dumps a window of a banded dual DP matrix for debugging.
 */
int
p7_bandmx_DumpWindow(FILE *ofp, P7_BANDMX *bmx, int istart, int iend, int kstart, int kend)
{
  int   *bnd_ip = bmx->bnd->imem;  	/* array of ia..ib pairs for 0..nseg-1 segments: [ia0 ib0][ia1 ib1]... */
  int   *bnd_kp = bmx->bnd->kmem;	/* array of ka..kb pairs for each banded row i in segments */
  float *dp     = bmx->dp;	/* we'll use <dp> to step through banded DP matrix, traversing banded structure */
  float *xp     = bmx->xmx;	/* ditto for <xp> stepping through xmx - remember it has extra ia-1 row for each segment */
  int    g, i, k, x;
  int    ia, ib;
  int    ka, kb;
  int    width     = 9;
  int    precision = 4;
  int    show_interseg = TRUE;

  /* Header line 1: model coords, special state labels */
  fputs("      ", ofp);
  for (k = kstart; k <= kend;  k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0; x < p7B_NXCELLS; x++) fprintf(ofp, "%*s ", width, p7_bandmx_DecodeSpecial(x));
  fputc('\n', ofp);

  /* Header line 2: underlining */
  fputs("       ", ofp);
  for (k = kstart; k <= kend;  k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < p7B_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fputc('\n', ofp);

  i = 0;
  for (g = 0; g < bmx->bnd->nseg; g++)
    {
      ia = *bnd_ip++;
      ib = *bnd_ip++;

      if (ia > i+1 && show_interseg) {
	fputs("...\n\n", ofp); 
	show_interseg = FALSE;
      }

      /* Row ia-1 off segment:
       * Show the specials on unbanded row ia-1 at start of each segment. 
       */
      if (ia-1 >= istart && ia-1 <= iend)
	{
	  fprintf(ofp, "%3d -- ", ia-1);
	  for (k = kstart; k <= kend;  k++) fprintf  (ofp, "%*s ", width, ".....");
	  for (x = 0; x < p7B_NXCELLS; x++) fprintf(ofp, "%*.*f ", width, precision, xp[x]);
	  fputs("\n\n", ofp);
	}
      xp += p7B_NXCELLS;

      for (i = ia; i <= ib; i++)
	{
	  /* Must always iterate through all bands, regardless of where the window is */
	  ka = *bnd_kp++;
	  kb = *bnd_kp++;

	  if (i >= istart && i <= iend) 
	    {
	      fprintf(ofp, "%3d ML ", i);
	      for (k = kstart; k <= ESL_MIN(kend,ka-1); k++) fprintf(ofp, "%*s ",   width, ".....");
	      for (          ; k <= ESL_MIN(kend,  kb); k++) fprintf(ofp, "%*.*f ", width, precision, dp[(k-ka)*p7B_NSCELLS + p7B_ML]);
	      for (          ; k <= kend;               k++) fprintf(ofp, "%*s ",   width, ".....");
	      for (x = 0; x < p7B_NXCELLS; x++) fprintf(ofp, "%*.*f ", width, precision, xp[x]);
	      fputc('\n', ofp);

	      fprintf(ofp, "%3d IL ", i);
	      for (k = kstart; k <= ESL_MIN(kend,ka-1); k++) fprintf(ofp, "%*s ",   width, ".....");
	      for (          ; k <= ESL_MIN(kend,  kb); k++) fprintf(ofp, "%*.*f ", width, precision, dp[(k-ka)*p7B_NSCELLS + p7B_IL]);
	      for (          ; k <= kend;               k++) fprintf(ofp, "%*s ",   width, ".....");
	      fputc('\n', ofp);

	      fprintf(ofp, "%3d DL ", i);
	      for (k = kstart; k <= ESL_MIN(kend,ka-1); k++) fprintf(ofp, "%*s ",   width, ".....");
	      for (          ; k <= ESL_MIN(kend,  kb); k++) fprintf(ofp, "%*.*f ", width, precision, dp[(k-ka)*p7B_NSCELLS + p7B_DL]);
	      for (          ; k <= kend;               k++) fprintf(ofp, "%*s ",   width, ".....");
	      fputc('\n', ofp);

	      fprintf(ofp, "%3d MG ", i);
	      for (k = kstart; k <= ESL_MIN(kend,ka-1); k++) fprintf(ofp, "%*s ",   width, ".....");
	      for (          ; k <= ESL_MIN(kend,  kb); k++) fprintf(ofp, "%*.*f ", width, precision, dp[(k-ka)*p7B_NSCELLS + p7B_MG]);
	      for (          ; k <= kend;               k++) fprintf(ofp, "%*s ",   width, ".....");
	      fputc('\n', ofp);

	      fprintf(ofp, "%3d IG ", i);
	      for (k = kstart; k <= ESL_MIN(kend,ka-1); k++) fprintf(ofp, "%*s ",   width, ".....");
	      for (          ; k <= ESL_MIN(kend,  kb); k++) fprintf(ofp, "%*.*f ", width, precision, dp[(k-ka)*p7B_NSCELLS + p7B_IG]);
	      for (          ; k <= kend;               k++) fprintf(ofp, "%*s ",   width, ".....");
	      fputc('\n', ofp);

	      fprintf(ofp, "%3d DG ", i);
	      for (k = kstart; k <= ESL_MIN(kend,ka-1); k++) fprintf(ofp, "%*s ",   width, ".....");
	      for (          ; k <= ESL_MIN(kend,  kb); k++) fprintf(ofp, "%*.*f ", width, precision, dp[(k-ka)*p7B_NSCELLS + p7B_DG]);
	      for (          ; k <= kend;               k++) fprintf(ofp, "%*s ",   width, ".....");
	      fputs("\n\n", ofp);

	      show_interseg = TRUE;
	    }

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
