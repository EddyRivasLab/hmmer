/* Implementation plugin:
 * Original HMMER2 implementation of dynamic programming algorithms.
 * 
 * Intended for reference and regression testing.
 * 
 * Implementation plugins provide a standard API, including support
 * for two objects, a P7_OPROFILE score profile and a P7_OMX dynamic
 * programming matrix.
 * 
 * Contents:
 *   1. The P7_OPROFILE structure: a score profile.
 *   2. The P7_OMX structure: a dynamic programming matrix.
 *   3. Configuration of/conversion to an optimized profile.
 *   4. Dynamic programming algorithm implementations.
 *   5. MPI communication.
 *   6. Benchmark driver.
 * 
 * SRE, Fri Jul 13 08:05:37 2007 [Janelia] [Enigma]
 * SVN $Id$
 */

#include "p7_config.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_h2.h"


/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create()
 * Synopsis:  Create a profile.
 * Incept:    SRE, Thu Jan 11 15:53:28 2007 [Janelia]
 *
 * Purpose:   Creates a profile of <M> nodes, for digital alphabet <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_OPROFILE *
p7_oprofile_Create(int M, const ESL_ALPHABET *abc)
{
  P7_OPROFILE *om = NULL;
  int          x;
  int          status;

  /* level 0 */
  ESL_ALLOC(om, sizeof(P7_OPROFILE));
  om->tsc   = om->msc = om->isc = NULL;
  om->bsc   = om->esc = NULL;
  
  /* level 1 */
  ESL_ALLOC(om->tsc, sizeof(int *) * p7O_NTRANS); om->tsc[0] = NULL;
  ESL_ALLOC(om->msc, sizeof(int *) * abc->Kp);    om->msc[0] = NULL;
  ESL_ALLOC(om->isc, sizeof(int *) * abc->Kp);    om->isc[0] = NULL;
  ESL_ALLOC(om->bsc, sizeof(int) * (M+1));      
  ESL_ALLOC(om->esc, sizeof(int) * (M+1));      

  /* level 2 */
  ESL_ALLOC(om->tsc[0], sizeof(int) * p7O_NTRANS * M);    
  ESL_ALLOC(om->msc[0], sizeof(int) * abc->Kp * (M+1));
  ESL_ALLOC(om->isc[0], sizeof(int) * abc->Kp * M);
  for (x = 1; x < abc->Kp; x++) {
    om->msc[x] = om->msc[0] + x * (M+1);
    om->isc[x] = om->isc[0] + x * M;
  }
  for (x = 0; x < p7O_NTRANS; x++)
    om->tsc[x] = om->tsc[0] + x * M;

  /* Initialize some pieces of memory that are never used,
   * only there for indexing convenience.
   */
  om->tsc[p7O_MM][0] = p7_IMPOSSIBLE; /* node 0 nonexistent, has no transitions  */
  om->tsc[p7O_MI][0] = p7_IMPOSSIBLE;
  om->tsc[p7O_MD][0] = p7_IMPOSSIBLE;
  om->tsc[p7O_IM][0] = p7_IMPOSSIBLE;
  om->tsc[p7O_II][0] = p7_IMPOSSIBLE;
  om->tsc[p7O_DM][0] = p7_IMPOSSIBLE;
  om->tsc[p7O_DD][0] = p7_IMPOSSIBLE;
  om->tsc[p7O_DM][1] = p7_IMPOSSIBLE; /* delete state D_1 is wing-retracted */
  om->tsc[p7O_DD][1] = p7_IMPOSSIBLE;
  for (x = 0; x < abc->Kp; x++) {     /* no emissions from nonexistent M_0, I_0 */
    om->msc[x][0] = p7_IMPOSSIBLE;
    om->isc[x][0] = p7_IMPOSSIBLE;
  }
  x = esl_abc_XGetGap(abc);	      /* no emission can emit/score gap characters */
  esl_vec_ISet(om->msc[x], M+1, p7_IMPOSSIBLE);
  esl_vec_ISet(om->isc[x], M,   p7_IMPOSSIBLE);
  x = esl_abc_XGetMissing(abc);	      /* no emission can emit/score missing data characters */
  esl_vec_ISet(om->msc[x], M+1, p7_IMPOSSIBLE);
  esl_vec_ISet(om->isc[x], M,   p7_IMPOSSIBLE);
  om->bsc[0] = p7_IMPOSSIBLE;
  om->esc[0] = p7_IMPOSSIBLE;

  /* Set remaining info
   */
  om->mode        = p7_NO_MODE;
  om->M           = M;
  om->abc         = abc;
  return om;

 ERROR:
  p7_oprofile_Destroy(om);
  return NULL;
}

/* Function:  p7_oprofile_Destroy()
 * Synopsis:  Frees a profile.
 * Incept:    SRE, Thu Jan 11 15:54:17 2007 [Janelia]
 *
 * Purpose:   Frees profile <om>.
 */
void
p7_oprofile_Destroy(P7_OPROFILE *om)
{
  if (om != NULL) {
    if (om->tsc   != NULL && om->tsc[0] != NULL) free(om->tsc[0]);
    if (om->msc   != NULL && om->msc[0] != NULL) free(om->msc[0]);
    if (om->isc   != NULL && om->isc[0] != NULL) free(om->isc[0]);
    if (om->tsc   != NULL) free(om->tsc);
    if (om->msc   != NULL) free(om->msc);
    if (om->isc   != NULL) free(om->isc);
    if (om->bsc   != NULL) free(om->bsc);
    if (om->esc   != NULL) free(om->esc);
  }
  free(om);
  return;
}

/*****************************************************************
 * 2. The P7_OMX structure: a dynamic programming matrix.
 *****************************************************************/

/* Function:  p7_omx_Create()
 * Synopsis:  Create a dynamic programming matrix.
 * Incept:    SRE, Tue Jan 30 11:20:33 2007 [Einstein's, in St. Louis]
 *
 * Purpose:   Allocate a reusable, resizeable <P7_OMX> for models up to
 *            size <allocM> and sequences up to length <allocL>.
 *
 * Returns:   a pointer to the new <P7_OMX>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_OMX *
p7_omx_Create(int allocM, int allocL)
{
  int     status;
  P7_OMX *ox = NULL;
  int     i;

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->xmx     = ox->mmx     = ox->imx     = ox->dmx     = NULL;
  ox->xmx_mem = ox->mmx_mem = ox->imx_mem = ox->dmx_mem = NULL;

  ESL_ALLOC(ox->xmx, sizeof(int *) * (allocL+1)); 
  ESL_ALLOC(ox->mmx, sizeof(int *) * (allocL+1)); 
  ESL_ALLOC(ox->imx, sizeof(int *) * (allocL+1)); 
  ESL_ALLOC(ox->dmx, sizeof(int *) * (allocL+1)); 
  
  ESL_ALLOC(ox->xmx_mem, p7X_NXCELLS* (allocL+1) * sizeof(int));
  ESL_ALLOC(ox->mmx_mem, (allocM+1) * (allocL+1) * sizeof(int));
  ESL_ALLOC(ox->imx_mem, (allocM+1) * (allocL+1) * sizeof(int));
  ESL_ALLOC(ox->dmx_mem, (allocM+1) * (allocL+1) * sizeof(int));

  for (i = 0; i <= allocL; i++) {
    ox->xmx[i] = ox->xmx_mem + i * p7X_NXCELLS;
    ox->mmx[i] = ox->mmx_mem + i * (allocM+1);
    ox->imx[i] = ox->imx_mem + i * (allocM+1);
    ox->dmx[i] = ox->dmx_mem + i * (allocM+1);
  }

  ox->M      = allocM;
  ox->L      = allocL;
  ox->ncells = (allocM+1)*(allocL+1);
  ox->nrows  = (allocL+1);
  return ox;

 ERROR:
  if (ox != NULL) p7_omx_Destroy(ox);
  return NULL;
}

/* Function:  p7_omx_GrowTo()
 * Synopsis:  Expand a DP matrix, if needed.
 * Incept:    SRE, Tue Jan 30 11:31:23 2007 [Olin Library, St. Louis]
 *
 * Purpose:   Assures that a DP matrix <ox> is allocated
 *            for a model of size up to <allocM> and a sequence of
 *            length up to <allocL>; reallocates if necessary.
 *            
 *            This function does not respect the configured
 *            <p7_RAMLIMIT>; it will allocate what it's told to
 *            allocate. 
 *
 * Returns:   <eslOK> on success, and <ox> may be reallocated upon
 *            return; any data that may have been in <ox> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <ox> must be assumed to be invalidated.
 */
int
p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL)
{
  int     status;
  void   *p;
  int     i;
  size_t  ncells, nrows;

  ncells = (allocM+1)*(allocL+1);
  nrows  = allocL+1;
  if (ncells <= ox->ncells && nrows <= ox->nrows) return eslOK;

  /* must we realloc the 2D matrices? (or can we get away with just
   * jiggering the row pointers, if we are growing in one dimension
   * while shrinking in another?)
   */
  if (ncells > ox->ncells) {
    ESL_RALLOC(ox->mmx_mem, p, sizeof(int) * ncells);
    ESL_RALLOC(ox->imx_mem, p, sizeof(int) * ncells);
    ESL_RALLOC(ox->dmx_mem, p, sizeof(int) * ncells);
  }

  /* must we reallocate the row pointers? */
  if (nrows > ox->nrows)
    {
      ESL_RALLOC(ox->xmx_mem, p, sizeof(int *) * nrows * p7X_NXCELLS);
      ESL_RALLOC(ox->xmx,     p, sizeof(int *) * nrows);
      ESL_RALLOC(ox->mmx,     p, sizeof(int *) * nrows);
      ESL_RALLOC(ox->imx,     p, sizeof(int *) * nrows);
      ESL_RALLOC(ox->dmx,     p, sizeof(int *) * nrows);
    }

  /* reset all the row pointers (even though we may not have to: if we
   * increased allocL without reallocating cells or changing M,
   * there's some wasted cycles here - but I'm not sure this case
   * arises, and I'm not sure it's worth thinking hard about. )
   */
  for (i = 0; i <= allocL; i++) {
    ox->xmx[i] = ox->xmx_mem + i * p7X_NXCELLS;
    ox->mmx[i] = ox->mmx_mem + i * (allocM+1);
    ox->imx[i] = ox->imx_mem + i * (allocM+1);
    ox->dmx[i] = ox->dmx_mem + i * (allocM+1);
  }

  ox->M      = allocM;
  ox->L      = allocL;
  ox->ncells = ncells;
  ox->nrows  = nrows;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_omx_Destroy()
 * Synopsis:  Free a DP matrix.
 * Incept:    SRE, Tue Jan 30 11:17:36 2007 [Einstein's, in St. Louis]
 *
 * Purpose:   Frees a <P7_OMX>.
 *
 * Returns:   (void)
 */
void
p7_omx_Destroy(P7_OMX *ox)
{
  if (ox != NULL) {
    if (ox->xmx     != NULL) free(ox->xmx);
    if (ox->mmx     != NULL) free(ox->mmx);
    if (ox->imx     != NULL) free(ox->imx);
    if (ox->dmx     != NULL) free(ox->dmx);
    if (ox->xmx_mem != NULL) free(ox->xmx_mem);
    if (ox->mmx_mem != NULL) free(ox->mmx_mem);
    if (ox->imx_mem != NULL) free(ox->imx_mem);
    if (ox->dmx_mem != NULL) free(ox->dmx_mem);
    free(ox);
  }
  return;
}

/* Function:  p7_omx_Dump()
 * Synopsis:  Dump a DP matrix to a stream, for diagnostics.
 * Incept:    SRE, Wed Jul 18 08:35:17 2007 [Janelia]
 *
 * Purpose:   Dump matrix <ox> to stream <fp> for diagnostics.
 */
int
p7_omx_Dump(FILE *ofp, P7_OMX *ox)
{
  int i, k, x;

  /* Header */
  fprintf(ofp, "     ");
  for (k = 0; k <= ox->M;  k++) fprintf(ofp, "%8d ", k);
  fprintf(ofp, "%8s ", "E");
  fprintf(ofp, "%8s ", "N");
  fprintf(ofp, "%8s ", "J");
  fprintf(ofp, "%8s ", "B");
  fprintf(ofp, "%8s\n","C");
  for (k = 0; k <= ox->M+5; k++) fprintf(ofp, "%8s ", "--------");
  fprintf(ofp, "\n");
  
  /* DP matrix data */
  for (i = 0; i <= ox->L; i++)
    {
      fprintf(ofp, "%3d M ", i);
      for (k = 0; k <= ox->M;       k++) fprintf(ofp, "%8.4f ", p7_SILO2Lod(ox->mmx[i][k]));
      for (x = 0; x <  p7G_NXCELLS; x++) fprintf(ofp, "%8.4f ", p7_SILO2Lod(ox->xmx[i][x]));
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d I ", i);
      for (k = 0; k <= ox->M;       k++) fprintf(ofp, "%8.4f ", p7_SILO2Lod(ox->imx[i][k]));
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d D ", i);
      for (k = 0; k <= ox->M;       k++) fprintf(ofp, "%8.4f ", p7_SILO2Lod(ox->dmx[i][k]));
      fprintf(ofp, "\n\n");
    }
  return eslOK;
}


/*****************************************************************
 * 3. Configuration of/conversion to an optimized profile.
 *****************************************************************/

static int convert_lodscore(float x)
{
  if (x == -eslINFINITY) return p7_IMPOSSIBLE;
  else                   return round(x * p7_INTSCALE);
}

/* Function:  p7_oprofile_Convert()
 * Synopsis:  Convert generic profile to optimized profile.
 * Incept:    SRE, Fri Jul 13 08:28:15 2007 [Janelia]
 *
 * Purpose:   Use the <gm> profile scores to fill in the scores of <om>,
 *            rearranging as needed.
 * 
 * Returns:   <eslOK> on success.
 */
int
p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om)
{
  int  M = gm->M;
  int  k, x;

  /* Contract checks */
  if (gm->M           != om->M)       ESL_EXCEPTION(eslEINVAL, "profile sizes don't match");
  if (gm->abc->type != om->abc->type) ESL_EXCEPTION(eslEINVAL, "alphabet types don't match");
  
  /* Transition scores */
  for (k = 1 ; k < M; k++) {
    om->tsc[p7O_MM][k] = convert_lodscore(p7P_TSC(gm, k, p7P_MM));
    om->tsc[p7O_MI][k] = convert_lodscore(p7P_TSC(gm, k, p7P_MI));
    om->tsc[p7O_MD][k] = convert_lodscore(p7P_TSC(gm, k, p7P_MD));
    om->tsc[p7O_IM][k] = convert_lodscore(p7P_TSC(gm, k, p7P_IM));
    om->tsc[p7O_II][k] = convert_lodscore(p7P_TSC(gm, k, p7P_II));
    om->tsc[p7O_DM][k] = convert_lodscore(p7P_TSC(gm, k, p7P_DM));
    om->tsc[p7O_DD][k] = convert_lodscore(p7P_TSC(gm, k, p7P_DD));
  }

  /* Begin scores (watch off-by-one in general profile) */
  for (k = 0; k < M; k++)     om->bsc[k+1] = convert_lodscore(p7P_TSC(gm, k, p7P_BM));

  /* End scores: in H3, these are either 0 or impossible */
  for (k = 1; k < gm->M; k++) om->esc[k] = (p7_profile_IsLocal(gm) ? 0 : p7_IMPOSSIBLE);
  om->esc[M] = 0;

  /* Match, insert scores */
  for (x = 0; x < gm->abc->Kp; x++) {
    for (k = 1; k < M; k++) {
      om->msc[x][k] = convert_lodscore(p7P_MSC(gm, k, x));
      om->isc[x][k] = convert_lodscore(p7P_ISC(gm, k, x));
    }
    om->msc[x][M]   = convert_lodscore(p7P_MSC(gm, M, x));
  }

  /* Specials */
  om->xsc[p7O_E][p7O_LOOP] = convert_lodscore(gm->xsc[p7P_E][p7P_LOOP]);  
  om->xsc[p7O_E][p7O_MOVE] = convert_lodscore(gm->xsc[p7P_E][p7P_MOVE]);
  om->xsc[p7O_N][p7O_LOOP] = convert_lodscore(gm->xsc[p7P_N][p7P_LOOP]);
  om->xsc[p7O_N][p7O_MOVE] = convert_lodscore(gm->xsc[p7P_N][p7P_MOVE]);
  om->xsc[p7O_C][p7O_LOOP] = convert_lodscore(gm->xsc[p7P_C][p7P_LOOP]);
  om->xsc[p7O_C][p7O_MOVE] = convert_lodscore(gm->xsc[p7P_C][p7P_MOVE]);
  om->xsc[p7O_J][p7O_LOOP] = convert_lodscore(gm->xsc[p7P_J][p7P_LOOP]);
  om->xsc[p7O_J][p7O_MOVE] = convert_lodscore(gm->xsc[p7P_J][p7P_MOVE]);

  return eslOK;
}


/*****************************************************************
 * 4. Dynamic programming algorithm implementations.
 *****************************************************************/

/* Function:  p7_Viterbi()
 * Synopsis:  The Viterbi algorithm.
 * Incept:    SRE, Tue Jan 30 10:50:53 2007 [Einstein's, St. Louis]
 * 
 * Purpose:   The standard Viterbi dynamic programming algorithm. 
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <gm>, and DP matrix <gm> allocated for at least <L>
 *            by <gm->M> cells; calculate the maximum scoring path by
 *            Viterbi; return the Viterbi score in <ret_sc>, and the
 *            Viterbi matrix is in <mx>.
 *            
 *            The caller may then retrieve the Viterbi path by calling
 *            <p7_Trace()>.
 *           
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm     - profile. Does not need to contain any
 *                     reference pointers (alphabet, HMM, or null model)
 *            mx     - DP matrix with room for an MxL alignment
 *            ret_sc - RETURN: Viterbi lod score in nats.
 *           
 * Return:   <eslOK> on success.
 */
int
p7_Viterbi(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;

  /* Some convenience */
  xmx = ox->xmx;
  mmx = ox->mmx;
  imx = ox->imx;
  dmx = ox->dmx;

  /* Initialization of the zero row.  */
  xmx[0][p7X_N] = 0;		                                  /* S->N, p=1            */
  xmx[0][p7X_B] = om->xsc[p7O_N][p7O_MOVE];                       /* S->N->B, no N-tail   */
  xmx[0][p7X_E] = xmx[0][p7X_C] = xmx[0][p7X_J] = p7_IMPOSSIBLE;  /* need seq to get here */
  for (k = 0; k <= om->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = p7_IMPOSSIBLE;            /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = impossible for all eight transitions (no node 0)
   *    I_M is wastefully calculated (doesn't exist)
   *    D_M is wastefully calculated (provably can't appear in a Viterbi path)
   */
  for (i = 1; i <= L; i++) {
    mmx[i][0] = imx[i][0] = dmx[i][0] = p7_IMPOSSIBLE;

    for (k = 1; k <= om->M; k++) {
				/* match state */
      mmx[i][k]  = p7_IMPOSSIBLE;
      if ((sc = mmx[i-1][k-1] + om->tsc[p7O_MM][k-1]) > mmx[i][k]) mmx[i][k] = sc;
      if ((sc = imx[i-1][k-1] + om->tsc[p7O_IM][k-1]) > mmx[i][k]) mmx[i][k] = sc;
      if ((sc = xmx[i-1][p7X_B] + om->bsc[k])         > mmx[i][k]) mmx[i][k] = sc;
      if ((sc = dmx[i-1][k-1] + om->tsc[p7O_DM][k-1]) > mmx[i][k]) mmx[i][k] = sc;
      if (om->msc[dsq[i]][k] != p7_IMPOSSIBLE) mmx[i][k] += om->msc[dsq[i]][k];
      else                                     mmx[i][k] =  p7_IMPOSSIBLE;

				/* delete state */
      dmx[i][k] = p7_IMPOSSIBLE;
      if ((sc = mmx[i][k-1] + om->tsc[p7O_MD][k-1]) > dmx[i][k])  dmx[i][k] = sc;
      if ((sc = dmx[i][k-1] + om->tsc[p7O_DD][k-1]) > dmx[i][k])  dmx[i][k] = sc;

				/* insert state */
      if (k < om->M) {
	imx[i][k] = p7_IMPOSSIBLE;
	if ((sc = mmx[i-1][k] + om->tsc[p7O_MI][k]) > imx[i][k])  imx[i][k] = sc;
	if ((sc = imx[i-1][k] + om->tsc[p7O_II][k]) > imx[i][k])  imx[i][k] = sc;

	if (om->isc[dsq[i]][k] != p7_IMPOSSIBLE) 
	  imx[i][k] += om->isc[dsq[i]][k];
	else
	  imx[i][k] = p7_IMPOSSIBLE;   
      }
    }

    /* Now the special states. Order is important here.
     * remember, N, C and J emissions are zero score by definition.
     */
				/* N state */
    xmx[i][p7X_N] = p7_IMPOSSIBLE;
    if ((sc = xmx[i-1][p7X_N] + om->xsc[p7O_N][p7O_LOOP]) > p7_IMPOSSIBLE)
      xmx[i][p7X_N] = sc;

				/* E state */
    xmx[i][p7X_E] = p7_IMPOSSIBLE;
    for (k = 1; k <= om->M; k++) {
      if ((sc =  mmx[i][k] + om->esc[k]) > xmx[i][p7X_E]) xmx[i][p7X_E] = sc;
      if ((sc =  dmx[i][k] + om->esc[k]) > xmx[i][p7X_E]) xmx[i][p7X_E] = sc;
    }
				/* J state */
    xmx[i][p7X_J] = p7_IMPOSSIBLE;
    if ((sc = xmx[i-1][p7X_J] + om->xsc[p7O_J][p7O_LOOP]) > p7_IMPOSSIBLE)
      xmx[i][p7X_J] = sc;
    if ((sc = xmx[i][p7X_E]   + om->xsc[p7O_E][p7O_LOOP]) > xmx[i][p7X_J])
      xmx[i][p7X_J] = sc;
				/* B state */
    xmx[i][p7X_B] = p7_IMPOSSIBLE;
    if ((sc = xmx[i][p7X_N] + om->xsc[p7O_N][p7O_MOVE]) > p7_IMPOSSIBLE)
      xmx[i][p7X_B] = sc;
    if ((sc = xmx[i][p7X_J] + om->xsc[p7O_J][p7O_MOVE]) > xmx[i][p7X_B])
      xmx[i][p7X_B] = sc;

				/* C state */
    xmx[i][p7X_C] = p7_IMPOSSIBLE;
    if ((sc = xmx[i-1][p7X_C] + om->xsc[p7O_C][p7O_LOOP]) > p7_IMPOSSIBLE)
      xmx[i][p7X_C] = sc;
    if ((sc = xmx[i][p7X_E] + om->xsc[p7O_E][p7O_MOVE]) > xmx[i][p7X_C])
      xmx[i][p7X_C] = sc;
  }
				/* T state (not stored) */
  *ret_sc = p7_SILO2Lod(xmx[L][p7X_C] + om->xsc[p7O_C][p7O_MOVE]);
  return eslOK;
}

/* Function:  p7_Forward()
 * Synopsis:  The Forward algorithm.
 * Incept:    SRE, Mon Apr 16 13:57:35 2007 [Janelia]
 *
 * Purpose:   The Forward dynamic programming algorithm. 
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <om>, and DP matrix <ox> allocated for at least <om->M>
 *            by <L> cells; calculate the probability of the sequence
 *            given the model using the Forward algorithm; return the
 *            Forward matrix in <gx>, and the Forward score in <ret_sc>.
 *           
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm     - profile. Does not need to contain any
 *                     reference pointers (alphabet, HMM, or null model)
 *            mx     - DP matrix with room for an MxL alignment
 *            ret_sc - RETURN: Forward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
p7_Forward(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;

  /* Some convenience */
  xmx = ox->xmx;
  mmx = ox->mmx;
  imx = ox->imx;
  dmx = ox->dmx;

  /* Initialization of the zero row.
   * Note that xmx[i][stN] = 0 by definition for all i,
   *    and xmx[i][stT] = xmx[i][stC], so neither stN nor stT need
   *    to be calculated in DP matrices.
   */
  xmx[0][p7X_N] = 0;		                     /* S->N, p=1            */
  xmx[0][p7X_B] = om->xsc[p7O_N][p7O_MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][p7X_E] = xmx[0][p7X_C] = xmx[0][p7X_J] = p7_IMPOSSIBLE;  /* need seq to get here */
  for (k = 0; k <= om->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = p7_IMPOSSIBLE;      /* need seq to get here */
  p7_ILogsumInit();

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = impossible for all eight transitions (no node 0)
   *    I_M is wastefully calculated (doesn't exist)
   */
  for (i = 1; i <= L; i++) {
    mmx[i][0] = imx[i][0] = dmx[i][0] = p7_IMPOSSIBLE;
    
    for (k = 1; k <= om->M; k++)
     {
       /* match state */
       mmx[i][k]  = p7_ILogsum(p7_ILogsum(mmx[i-1][k-1] + om->tsc[p7O_MM][k-1],
					  imx[i-1][k-1] + om->tsc[p7O_IM][k-1]),
			       p7_ILogsum(xmx[i-1][p7X_B] + om->bsc[k],
					  dmx[i-1][k-1] + om->tsc[p7O_DM][k-1]));
       mmx[i][k] += om->msc[dsq[i]][k];
       if (mmx[i][k] < p7_IMPOSSIBLE) mmx[i][k] = p7_IMPOSSIBLE;

       dmx[i][k]  = p7_ILogsum(mmx[i][k-1] + om->tsc[p7O_MD][k-1],
			       dmx[i][k-1] + om->tsc[p7O_DD][k-1]);
       if (dmx[i][k] < p7_IMPOSSIBLE) dmx[i][k] = p7_IMPOSSIBLE;

       if (k < om->M) {
	 imx[i][k]  = p7_ILogsum(mmx[i-1][k] + om->tsc[p7O_MI][k],
				 imx[i-1][k] + om->tsc[p7O_II][k]);
	 imx[i][k] += om->isc[dsq[i]][k];
	 if (imx[i][k] < p7_IMPOSSIBLE) imx[i][k] = p7_IMPOSSIBLE;
       }
     }

   /* Now the special states.
    * remember, C and J emissions are zero score by definition
    */
   xmx[i][p7X_N] = xmx[i-1][p7X_N] + om->xsc[p7O_N][p7O_LOOP];
   if (xmx[i][p7X_N] < p7_IMPOSSIBLE) xmx[i][p7X_N] = p7_IMPOSSIBLE;

   xmx[i][p7X_E] = p7_IMPOSSIBLE;
   for (k = 1; k <= om->M; k++) {
     sc =            p7_ILogsum(    mmx[i][k] + om->esc[k],
				    dmx[i][k] + om->esc[k]);
     xmx[i][p7X_E] = p7_ILogsum(sc, xmx[i][p7X_E]);
   }
   if (xmx[i][p7X_E] < p7_IMPOSSIBLE) xmx[i][p7X_E] = p7_IMPOSSIBLE;

   xmx[i][p7X_J] = p7_ILogsum(xmx[i-1][p7X_J] + om->xsc[p7O_J][p7O_LOOP],
			      xmx[i][p7X_E]   + om->xsc[p7O_E][p7O_LOOP]);
   if (xmx[i][p7X_J] < p7_IMPOSSIBLE) xmx[i][p7X_J] = p7_IMPOSSIBLE;

   xmx[i][p7X_B] = p7_ILogsum(xmx[i][p7X_N] + om->xsc[p7O_N][p7O_MOVE],
			      xmx[i][p7X_J] + om->xsc[p7O_J][p7O_MOVE]);
   if (xmx[i][p7X_B] < p7_IMPOSSIBLE) xmx[i][p7X_B] = p7_IMPOSSIBLE;

   xmx[i][p7X_C] = p7_ILogsum(xmx[i-1][p7X_C] + om->xsc[p7O_C][p7O_LOOP],
			      xmx[i][p7X_E] + om->xsc[p7O_E][p7O_MOVE]);
   if (xmx[i][p7X_C] < p7_IMPOSSIBLE) xmx[i][p7X_C] = p7_IMPOSSIBLE;
  }
  
  *ret_sc = p7_SILO2Lod( xmx[L][p7X_C] + om->xsc[p7O_C][p7O_MOVE]);
  return eslOK;
}

/*****************************************************************
 * 5. MPI communication
 *****************************************************************/
#ifdef HAVE_MPI
/* Function:  p7_oprofile_MPISend()
 * Synopsis:  Send a profile as an MPI message.
 * Incept:    SRE, Fri Apr 20 13:55:47 2007 [Janelia]
 *
 * Purpose:   Sends profile <om> to MPI process <dest> (where
 *            <dest> ranges from 0..<nproc-1>), with MPI tag <tag> 
 *            for MPI communicator <comm>.
 *            
 *            In order to minimize alloc/free cycles in this routine,
 *            caller passes a pointer to a working buffer <*buf> of
 *            size <*nalloc> characters. If necessary (i.e. if <om> is
 *            too big to fit), <*buf> will be reallocated and <*nalloc>
 *            increased to the new size. As a special case, if <*buf>
 *            is <NULL> and <*nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *            
 *            If <om> is NULL, an end-of-data signal is sent, which
 *            <p7_profile_MPIRecv()> knows how to interpret.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Note:      This was tested against a version that simply issues a series
 *            of MPI_Send()'s, rather than Pack()'ing into a buffer
 *            and issuing one MPI_Send(). The packed version seems to
 *            be significantly faster, although benchmarking MPI
 *            programs is difficult, and variance on the results is high.
 *            
 *            To optimize communication still further, one might try to
 *            avoid the MPI_Pack()'s. It might be feasible to change the
 *            allocation of a profile such that it is allocated in one
 *            contiguous chunk of memory. And once one embarks on that, 
 *            the memory layout of the profile should also be optimized
 *            with respect to cache performance during DP alignment.
 */
int
p7_oprofile_MPISend(P7_OPROFILE *om, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   sz, n, position;

  /* First, figure out the size of the profile */
  if (om == NULL) { 
    MPI_Pack_size(1, MPI_INT, comm, &n); 
  } else {
    /* This will look wasteful, but the MPI spec doesn't guarantee that 
     * MPI_Pack_size(x, ...) + MPI_Pack_size(x, ... ) == MPI_Pack_size(2x, ...).
     * Indeed there are some hints in the spec that that's *not* true.
     * So we match our Pack_size calls exactly to our Pack calls.
     */
    n = 0;
    MPI_Pack_size(1,                     MPI_INT,   comm, &sz);   n +=   sz; /* mode      */
    MPI_Pack_size(p7O_NTRANS*om->M,      MPI_INT,   comm, &sz);   n +=   sz; /* tsc[0]    */
    MPI_Pack_size((om->M+1)*om->abc->Kp, MPI_INT,   comm, &sz);   n +=   sz; /* msc[0]    */
    MPI_Pack_size(om->M*om->abc->Kp,     MPI_INT,   comm, &sz);   n +=   sz; /* isc[0]    */
    MPI_Pack_size(p7O_NXTRANS,           MPI_INT,   comm, &sz);   n += p7O_NXSTATES*sz; /* xsc[0..3] */
    MPI_Pack_size(om->M+1,               MPI_INT,   comm, &sz);   n += 2*sz; /* bsc, esc  */
  }
  
  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Pack the profile into the buffer */
  position = 0;
  if (om == NULL) 
    {
      int   eod_code = -1;
      MPI_Pack(&eod_code,                      1, MPI_INT,   *buf, n, &position,  comm);
    } 
  else 
    {    
      MPI_Pack(&(om->M),                       1, MPI_INT,   *buf, n, &position,  comm);
      MPI_Pack(&(om->mode),                    1, MPI_INT,   *buf, n, &position,  comm);
      MPI_Pack(om->tsc[0],      p7O_NTRANS*om->M, MPI_INT,   *buf, n, &position,  comm);
      MPI_Pack(om->msc[0], (om->M+1)*om->abc->Kp, MPI_INT,   *buf, n, &position,  comm);
      MPI_Pack(om->isc[0],     om->M*om->abc->Kp, MPI_INT,   *buf, n, &position,  comm);
      for (x = 0; x < p7O_NXSTATES; x++)
	MPI_Pack(om->xsc[x],         p7O_NXTRANS, MPI_INT,   *buf, n, &position,  comm);
      MPI_Pack(om->bsc,                  om->M+1, MPI_INT,   *buf, n, &position,  comm);
      MPI_Pack(om->esc,                  om->M+1, MPI_INT,   *buf, n, &position,  comm);
    }

  /* Send the packed profile to destination  */
  MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm);
  return eslOK;
  
 ERROR:
  return status;
}


/* Function:  p7_profile_MPIRecv()
 * Synopsis:  Receive a profile as an MPI message.
 * Incept:    SRE, Fri Apr 20 14:19:07 2007 [Janelia]
 *
 * Purpose:   Receive a profile from <source> (where <source> is usually
 *            process 0, the master) with tag <tag> from communicator <comm>,
 *            and return it in <*ret_om>. 
 *            
 *            Caller must also provide the alphabet <abc> and the
 *            background model <bg> for this profile. (Of course, that means
 *            the caller already knows them, by an appropriate
 *            initialization.)
 *            
 *            To minimize alloc/free cycles in this routine, caller
 *            passes a pointer to a buffer <*buf> of size <*nalloc>
 *            characters. These are passed by reference because if
 *            necessary, <buf> will be reallocated and <nalloc>
 *            increased to the new size. As a special case, if <buf>
 *            is <NULL> and <nalloc> is 0, the buffer will be
 *            allocated appropriately, but the caller is still
 *            responsible for free'ing it.
 *
 *            If the packed profile is an end-of-data signal, return
 *            <eslEOD>, and <*ret_om> is <NULL>.
 *            
 * Returns:   <eslOK> on success. <*ret_om> contains the new profile; it
 *            is allocated here, and the caller is responsible for
 *            free'ing it.  <*buf> may have been reallocated to a
 *            larger size, and <*nalloc> may have been increased.
 *            
 */
int
p7_profile_MPIRecv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc, const P7_BG *bg, char **buf, int *nalloc,  P7_OPROFILE **ret_om)
{
  int          status;
  P7_OPROFILE *om    = NULL;
  int          n, x;
  int          position;
  MPI_Status   mpistatus;
  int          M;

  /* Probe first, because we need to know if our buffer is big enough.
   */
  MPI_Probe(source, tag, comm, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Receive the packed profile */
  MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus);

  /* Unpack it - watching out for the EOD signal of M = -1. */
  position = 0;
  MPI_Unpack(*buf, n, &position, &M,                      1, MPI_INT,   comm);
  if (M == -1) { *ret_om = NULL; return eslEOD; }

  om = p7_profile_Create(M, abc);
  MPI_Unpack(*buf, n, &position, &(om->mode),             1, MPI_INT,   comm);
  MPI_Unpack(*buf, n, &position, om->tsc[0],   p7O_NTRANS*M, MPI_INT,   comm);
  MPI_Unpack(*buf, n, &position, om->msc[0],  (M+1)*abc->Kp, MPI_INT,   comm);
  MPI_Unpack(*buf, n, &position, om->isc[0],      M*abc->Kp, MPI_INT,   comm);
  for (x = 0; x < p7O_NXSTATES; x++)
    MPI_Unpack(*buf,n,&position, om->xsc[x],    p7O_NXTRANS, MPI_INT,   comm);  
  MPI_Unpack(*buf, n, &position, om->bsc,               M+1, MPI_INT,   comm);
  MPI_Unpack(*buf, n, &position, om->esc,               M+1, MPI_INT,   comm);
  
  om->abc = abc;
  *ret_om = om;
  return eslOK;

 ERROR:
  if (om  != NULL) p7_oprofile_Destroy(om);
  *ret_om = NULL;
  return status;
}
#endif /*HAVE_MPI*/

/*****************************************************************
 * 6. Benchmark driver.
 *****************************************************************/
#ifdef p7IMPL_H2_BENCHMARK
/* gcc -o benchmark-h2 -g -O2 -I. -L. -I../easel -L../easel -Dp7IMPL_H2_BENCHMARK impl_h2.c -lhmmer -leasel -lm
 * ./benchmark-h2 <hmmfile>
 */
/* 22 Mc/s for Viterbi, 7.3 Mc/s for Forward.
 * (gcc -g -O2, 3.2GHz Xeon, N=50K, L=400, M=72 RRM_1 model)
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-f",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "run Forward algorithm",                          0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual scores",             0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for old H2 implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_RANDOMNESS *r       = esl_randomness_Create(42);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  float           nullsc;
  float           bitscore;
  
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL);

  om = p7_oprofile_Create(hmm->M, abc);
  p7_oprofile_Convert(gm, om);

  ox = p7_omx_Create(om->M, L);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, bg->f, abc->K, L, dsq);

      if (esl_opt_GetBoolean(go, "-f"))  p7_Forward(dsq, L, om, ox, &sc);
      else                               p7_Viterbi(dsq, L, om, ox, &sc);
      p7_bg_NullOne(bg, dsq, L, &nullsc);
      
      if (esl_opt_GetBoolean(go, "-v")) printf("%.4f bits  (%.4f raw)\n", (sc-nullsc)/eslCONST_LOG2, sc); 
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  free(dsq);
  p7_omx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7DP_GENERIC_BENCHMARK*/



/*****************************************************************
 * 7. Unit tests.
 *****************************************************************/
#ifdef p7IMPL_H2_TESTDRIVE
#include <string.h>

#include "esl_getopts.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

/* The "basic" utest is a minimal driver for making a small DNA profile and a small DNA sequence,
 * then running Viterbi and Forward. It's useful for dumping DP matrices and profiles for debugging.
 */
static void
utest_basic(ESL_GETOPTS *go)
{
  char           *query= "# STOCKHOLM 1.0\n\nseq1 GAATTC\nseq2 GAATTC\n//\n";
  int             fmt  = eslMSAFILE_STOCKHOLM;
  char           *targ = "GAATTC";
  ESL_ALPHABET   *abc  = NULL;
  ESL_MSA        *msa  = NULL;
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_PROFILE     *om   = NULL;
  P7_BG          *bg   = NULL;
  P7_DPRIOR      *pri  = NULL;	
  ESL_DSQ        *dsq  = NULL;
  P7_OMX         *ox   = NULL;
  P7_GMX         *gx   = NULL;
  P7_TRACE       *tr   = NULL;
  int             L    = strlen(targ);
  float           vsc, vsc2, fsc, fsc2;

  if ((abc = esl_alphabet_Create(eslDNA))          == NULL)  esl_fatal("failed to create alphabet");
  if ((pri = p7_dprior_CreateNucleic())            == NULL)  esl_fatal("failed to create prior");
  if ((msa = esl_msa_CreateFromString(query, fmt)) == NULL)  esl_fatal("failed to create MSA");
  if (esl_msa_Digitize(abc, msa)                   != eslOK) esl_fatal("failed to digitize MSA");
  if (p7_Fastmodelmaker(msa, 0.5, &hmm, NULL)      != eslOK) esl_fatal("failed to create GAATTC model");
  if (p7_ParameterEstimation(hmm, pri)             != eslOK) esl_fatal("failed to parameterize GAATTC model");
  if ((bg = p7_bg_Create(abc))                     == NULL)  esl_fatal("failed to create DNA null model");
  if ((gm = p7_profile_Create(hmm->M, abc))        == NULL)  esl_fatal("failed to create GAATTC profile");
  if ((om = p7_oprofile_Create(hmm->M, abc))       == NULL)  esl_fatal("failed to create oprofile");
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL)!= eslOK) esl_fatal("failed to config profile");
  if (p7_oprofile_Convert(gm, om)                  != eslOK) esl_fatal("failed to convert profile");
  if (esl_abc_CreateDsq(abc, targ, &dsq)           != eslOK) esl_fatal("failed to create GAATTC digital sequence");
  if ((gx = p7_gmx_Create(gm->M, L))               == NULL)  esl_fatal("failed to create DP matrix, generic form");
  if ((ox = p7_omx_Create(gm->M, L))               == NULL)  esl_fatal("failed to create DP matrix, optimized form");
  if ((tr = p7_trace_Create())                     == NULL)  esl_fatal("trace creation failed");

  p7_Viterbi   (dsq, L, om, ox, &vsc);
  p7_GViterbi  (dsq, L, gm, gx, &vsc2);
  if (esl_FCompare(vsc, vsc2, 1e-5) != eslOK)  esl_fatal("Viterbi scores don't agree (%.4f,%.4f)", vsc, vsc2);

  p7_Forward    (dsq, L, om, ox, &fsc);
  p7_GForward   (dsq, L, gm, gx, &fsc2);
  if (esl_FCompare(fsc, fsc2, 1e-5) != eslOK)  esl_fatal("Forward scores don't agree (%.4f,%.4f)", fsc, fsc2);


  p7_trace_Destroy(tr);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  free(dsq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_msa_Destroy(msa);
  p7_dprior_Destroy(pri);
  esl_alphabet_Destroy(abc);
  return;
}
#endif /*p7IMPL_H2_TESTDRIVE*/

/*****************************************************************
 * 8. Test driver.
 *****************************************************************/
/* gcc -g -Wall -Dp7IMPL_H2_TESTDRIVE -I. -I../easel -L. -L../easel -o impl_h2_utest impl_h2.c -lhmmer -leasel -lm
 */
#ifdef p7IMPL_H2_TESTDRIVE
#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"

#include "p7_config.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the generic implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  char            errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS *r    = NULL;
  ESL_ALPHABET   *abc  = NULL;
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_OPROFILE    *om   = NULL;
  P7_BG          *bg   = NULL;
  int             M    = 100;
  int             L    = 200;
  int             nseq = 20;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  utest_basic(go);
  
  if ((abc = esl_alphabet_Create(eslAMINO))         == NULL)  esl_fatal("failed to create alphabet");
  if (p7_hmm_Sample(r, M, abc, &hmm)                != eslOK) esl_fatal("failed to sample an HMM");
  if ((bg = p7_bg_Create(abc))                      == NULL)  esl_fatal("failed to create null model");
  if ((gm = p7_profile_Create(hmm->M, abc))         == NULL)  esl_fatal("failed to create profile");
  if ((om = p7_oprofile_Create(hmm->M, abc))        == NULL)  esl_fatal("failed to create oprofile");
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL) != eslOK) esl_fatal("failed to config profile");
  if (p7_oprofile_Convert(gm, om)                   != eslOK) esl_fatal("failed to convert profile");
  if (p7_hmm_Validate    (hmm,  errbuf, 0.0001)     != eslOK) esl_fatal("whoops, HMM is bad!: %s", errbuf);
  if (p7_profile_Validate(gm,   errbuf, 0.0001)     != eslOK) esl_fatal("whoops, profile is bad!: %s", errbuf);

  utest_viterbi(r, abc, bg, gm, nseq, L);
  utest_forward(r, abc, bg, gm, nseq, L);
#if 0
  utest_Enumeration(r, abc, 4);	/* can't go much higher than 5; enumeration test is cpu-intensive. */
#endif
  
  p7_oprofile_Destroy(om)
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}

#endif /*p7DP_GENERIC_TESTDRIVE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
