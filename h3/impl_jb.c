/* Implementation plugin:
 * Optimized implementation of HMMER2 from Jeremy Buhler (Washington
 * University, St. Louis); adapted to HMMER3.
 *
 * Intended as a reference implementation of Forward() in scaled
 * integer scores, instead of the default floating-point
 * implementation.
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
 * SRE, Wed Jun 27 07:57:26 2007 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"
#include "impl_jb.h"

/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Wed Jun 27 08:07:01 2007 [Janelia]
 *
 * Purpose:   Create a profile of <M> nodes for digital alphabet <abc>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_OPROFILE *
p7_oprofile_Create(int M, const ESL_ALPHABET *abc)
{
  int          status;
  P7_OPROFILE *om = NULL;
  int          x;

  /* level 0 */
  ESL_ALLOC(om, sizeof(P7_OPROFILE));
  om->tsc = NULL;
  om->rsc = NULL;

  /* level 1 */
  ESL_ALLOC(om->tsc,   sizeof(int)   * M  * p7O_NTRANS);
  ESL_ALLOC(om->rsc,   sizeof(int *) * abc->Kp);                     
  om->rsc[0] = NULL;

  ESL_ALLOC(om->rsc[0],sizeof(int)   * abc->Kp * (M+1) * p7O_NR);
  for (x = 1; x < abc->Kp; x++)
    om->rsc[x] = om->rsc[0] + (x * (M+1) * p7O_NR);

  om->mode  = p7_NO_MODE;
  om->M     = M;
  om->abc   = abc;

  /* Initializations of unused model transitions 
   * tsc[0][0..6] are unused: set to -infty. [7=TBM is used, though]. These values are touched by DP algs.
   */
  for (x = 0; x < p7O_NTRANS; x++) om->tsc[x] = p7_IMPOSSIBLE;
  return om;

 ERROR:
  p7_oprofile_Destroy(om);
  return NULL;
}

/* Function:  p7_oprofile_Destroy()
 * Synopsis:  Free an optimized profile structure.
 * Incept:    SRE, Wed Jun 27 08:07:08 2007 [Janelia]
 *
 * Purpose:   Frees profile <om>.
 */
void
p7_oprofile_Destroy(P7_OPROFILE *om)
{
  if (om != NULL) {
    if (om->rsc   != NULL && om->rsc[0] != NULL) free(om->rsc[0]);
    if (om->tsc   != NULL) free(om->tsc);
    if (om->rsc   != NULL) free(om->rsc);
    free(om);
  }
  return;
}


/*****************************************************************
 * 2. The P7_OMX structure: a dynamic programming matrix.
 *****************************************************************/

/* Function:  p7_omx_Create()
 * Synopsis:  Create an optimized dynamic programming matrix.
 * Incept:    SRE, Wed Jun 27 08:07:29 2007 [Janelia]
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
  int      status;
  P7_OMX  *ox;
  int      i;

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->dp  = NULL;
  ox->xmx = NULL;

  ESL_ALLOC(ox->dp,  sizeof(int *) * (allocL+1));
  ESL_ALLOC(ox->xmx, sizeof(int)   * (allocL+1) * p7X_NXCELLS);
  ox->dp[0] = NULL;
  
  ESL_ALLOC(ox->dp[0], sizeof(int) * (allocL+1) * (allocM+1) * p7X_NSCELLS);
  for (i = 1; i <= allocL; i++)
    ox->dp[i] = ox->dp[0] + i*(allocM+1)*p7X_NSCELLS;
  
  /* Initialize memory that's allocated but unused, only to keep
   * valgrind and friends happy.
   */
  for (i = 0; i <= allocL; i++) 
    {
      ox->dp[i][0      * p7X_NSCELLS + p7X_M] = p7_IMPOSSIBLE; /* M_0 */
      ox->dp[i][0      * p7X_NSCELLS + p7X_I] = p7_IMPOSSIBLE; /* I_0 */      
      ox->dp[i][0      * p7X_NSCELLS + p7X_D] = p7_IMPOSSIBLE; /* D_0 */
      ox->dp[i][1      * p7X_NSCELLS + p7X_D] = p7_IMPOSSIBLE; /* D_1 */
      ox->dp[i][allocM * p7X_NSCELLS + p7X_I] = p7_IMPOSSIBLE; /* I_M */
    }



  ox->nrows  = allocL+1;
  ox->ncells = (allocM+1)*(allocL+1);
  ox->M      = allocM;
  ox->L      = allocL;
  return ox;

 ERROR:
  p7_omx_Destroy(ox);
  return NULL;
}

/* Function:  p7_omx_Destroy()
 * Synopsis:  Frees an optimized DP matrix.
 * Incept:    SRE, Wed Jun 27 08:07:37 2007 [Janelia]
 *
 * Purpose:   Frees optimized DP matrix <ox>.
 *
 * Returns:   (void)
 */
void
p7_omx_Destroy(P7_OMX *ox)
{
  if (ox == NULL) return;
  if (ox->dp != NULL) {
    if (ox->dp[0] != NULL) free(ox->dp[0]);
    free(ox->dp);
  }
  if (ox->xmx != NULL) free(ox->xmx);
  free(ox);
}


/* Function:  p7_omx_GrowTo()
 * Incept:    SRE, Wed Jul 18 09:12:34 2007 [Janelia]
 *
 * Purpose:   Assures that a DP matrix <ox> is allocated
 *            for a model of size up to <allocM> and a sequence of
 *            length up to <allocL>; reallocates if necessary.
 *            
 *            This function does not respect the configured
 *            <RAMLIMIT>; it will allocate what it's told to
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
  if (ncells <= ox->ncells && nrows <= ox->nrows) { ox->M = allocM; ox->L = allocL; return eslOK; }

  /* must we realloc the 2D matrices? (or can we get away with just
   * jiggering the row pointers, if we are growing in one dimension
   * while shrinking in another?)
   */
  if (ncells > ox->ncells) 
    {
      ESL_RALLOC(ox->dp[0], p, sizeof(int) * ncells * p7X_NSCELLS);
      ox->ncells = ncells;
    }

  /* must we reallocate the row pointers? */
  if (nrows > ox->nrows)
    {
      ESL_RALLOC(ox->xmx, p, sizeof(int)   * nrows * p7X_NXCELLS);
      ESL_RALLOC(ox->dp,  p, sizeof(int *) * nrows);
      ox->nrows  = nrows;
    }

  /* reset all the row pointers (even though we may not have to: if we
   * increased allocL without reallocating cells or changing M,
   * there's some wasted cycles here - but I'm not sure this case
   * arises, and I'm not sure it's worth thinking hard about. )
   */
  for (i = 0; i <= allocL; i++) 
    ox->dp[i] = ox->dp[0] + i * (allocM+1) * p7X_NSCELLS;

  ox->M      = allocM;
  ox->L      = allocL;
  return eslOK;

 ERROR:
  return status;
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
      for (k = 0; k <= ox->M;       k++) fprintf(ofp, "%8.4f ", p7_SILO2Lod(ox->dp[i][k * p7G_NSCELLS + p7G_M]));
      for (x = 0; x <  p7G_NXCELLS; x++) fprintf(ofp, "%8.4f ", p7_SILO2Lod(ox->xmx[  i * p7G_NXCELLS + x]));
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d I ", i);
      for (k = 0; k <= ox->M;       k++) fprintf(ofp, "%8.4f ", p7_SILO2Lod(ox->dp[i][k * p7G_NSCELLS + p7G_I]));
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d D ", i);
      for (k = 0; k <= ox->M;       k++) fprintf(ofp, "%8.4f ", p7_SILO2Lod(ox->dp[i][k * p7G_NSCELLS + p7G_D]));
      fprintf(ofp, "\n\n");
    }
  return eslOK;
}




/*****************************************************************
 * Configuration of/conversion to an optimized profile.
 *****************************************************************/

static int convert_lodscore(float x)
{
  if (x <= -eslINFINITY) return p7_IMPOSSIBLE;
  else                   return round(x * p7_INTSCALE);
}

/* Function:  p7_oprofile_Convert()
 * Synopsis:  Convert generic profile to optimized profile.
 * Incept:    SRE, Wed Jun 27 08:09:02 2007 [Janelia]
 *
 * Purpose:   Use the <gm> profile scores to fill in the scores of <om>,
 *            rearranging as needed.
 * 
 * Returns:   <eslOK> on success.
 */
int
p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om)
{
  int  k, x;
  
  /* Transition scores (except begin) */
  for (k = 1 ; k < gm->M; k++) {
    int *tsc = om->tsc + k*p7O_NTRANS;
    
    tsc[p7O_MM] = convert_lodscore(p7P_TSC(gm, k, p7P_MM));
    tsc[p7O_MI] = convert_lodscore(p7P_TSC(gm, k, p7P_MI));
    tsc[p7O_MD] = convert_lodscore(p7P_TSC(gm, k, p7P_MD));
    tsc[p7O_IM] = convert_lodscore(p7P_TSC(gm, k, p7P_IM));
    tsc[p7O_II] = convert_lodscore(p7P_TSC(gm, k, p7P_II));
    tsc[p7O_DM] = convert_lodscore(p7P_TSC(gm, k, p7P_DM));
    tsc[p7O_DD] = convert_lodscore(p7P_TSC(gm, k, p7P_DD));
    tsc[p7O_BM] = convert_lodscore(p7P_TSC(gm, k, p7P_BM));
  }

  /* Begin scores */
  for (k = 0 ; k < gm->M; k++) {
    int *tsc = om->tsc + k*p7O_NTRANS;
    tsc[p7O_BM] = convert_lodscore(p7P_TSC(gm, k, p7P_BM));
  }
  
  /* Match scores */
  for (x = 0; x < gm->abc_r->Kp; x++)
    for (k = 1; k <= gm->M; k++)
      {
	int *rsc = om->rsc[x] + k*p7O_NR;
	rsc[p7O_MSC] = convert_lodscore(p7P_MSC(gm, k, x));
      }

  /* Insert scores */
  for (x = 0; x < gm->abc_r->Kp; x++)
    for (k = 1; k < gm->M; k++)
      {
	int *rsc = om->rsc[x] + k*p7O_NR;
	rsc[p7O_ISC] = convert_lodscore(p7P_ISC(gm, k, x));
      }

  om->xsc[p7O_E][p7O_LOOP] = convert_lodscore(gm->xsc[p7P_E][p7P_LOOP]);  
  om->xsc[p7O_E][p7O_MOVE] = convert_lodscore(gm->xsc[p7P_E][p7P_MOVE]);
  om->xsc[p7O_N][p7O_LOOP] = convert_lodscore(gm->xsc[p7P_N][p7P_LOOP]);
  om->xsc[p7O_N][p7O_MOVE] = convert_lodscore(gm->xsc[p7P_N][p7P_MOVE]);
  om->xsc[p7O_C][p7O_LOOP] = convert_lodscore(gm->xsc[p7P_C][p7P_LOOP]);
  om->xsc[p7O_C][p7O_MOVE] = convert_lodscore(gm->xsc[p7P_C][p7P_MOVE]);
  om->xsc[p7O_J][p7O_LOOP] = convert_lodscore(gm->xsc[p7P_J][p7P_LOOP]);
  om->xsc[p7O_J][p7O_MOVE] = convert_lodscore(gm->xsc[p7P_J][p7P_MOVE]);

  om->mode  = gm->mode;
  om->abc   = gm->abc_r;
  return eslOK;
}

/*****************************************************************
 * 4. Dynamic programming algorithm implementations.
 *****************************************************************/

#define MMX(i,k) (dp[(i)][(k) * p7X_NSCELLS + p7X_M])
#define IMX(i,k) (dp[(i)][(k) * p7X_NSCELLS + p7X_I])
#define DMX(i,k) (dp[(i)][(k) * p7X_NSCELLS + p7X_D])
#define XMX(i,s) (xmx[(i) * p7X_NXCELLS + (s)])

#define TSC(s,k) (tsc[(k) * p7O_NTRANS + (s)])
#define MSC(k)   (rsc[(k) * p7O_NR     + p7O_MSC])
#define ISC(k)   (rsc[(k) * p7O_NR     + p7O_ISC])

/* Function:  p7_Viterbi()
 * Synopsis:  Viterbi algorithm.
 * Incept:    SRE, Wed Jun 27 08:49:53 2007 [Janelia]
 *
 * Purpose:   The Buhler (optimized integer) version of the Viterbi dynamic
 *            programming algorithm.
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <om>, and DP matrix <ox> allocated for at least <L>
 *            by <om->M> cells; calculate the maximum scoring path by
 *            Viterbi; return the Viterbi score in <ret_sc>, and the
 *            Viterbi matrix is in <ox>.
 *            
 *            The caller may then retrieve the Viterbi path by calling
 *            <p7_GTrace()>.
 *           
 *            The Viterbi lod score is returned in nats. The caller
 *            needs to subtract a null model lod score, then convert
 *            to bits.
 *           
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq
 *            om     - profile. Does not need to contain any
 *                     reference pointers (alphabet, HMM, or null model)
 *            ox     - DP matrix with room for an MxL alignment
 *            ret_sc - RETURN: Viterbi lod score in nats
 *           
 * Return:   <eslOK> on success.
 */
int
p7_Viterbi(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  int const *tsc = om->tsc;
  int      **dp  = ox->dp;
  int       *xmx = ox->xmx;
  int        M   = om->M;
  int        i,k;
  int        esc = p7_IsLocal(om->mode) ? 0 : p7_IMPOSSIBLE;

  /* Initialization of the zero row.  */
  XMX(0,p7X_N) = 0;                                                /* S->N, p=1            */
  XMX(0,p7X_B) = om->xsc[p7O_N][p7O_MOVE];                         /* S->N->B, no N-tail   */
  XMX(0,p7X_E) = XMX(0,p7X_C) = XMX(0,p7X_J) = p7_IMPOSSIBLE;      /* need seq to get here */
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = p7_IMPOSSIBLE;               /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   */
  for (i = 1; i <= L; i++) 
    {
      int const * rsc = om->rsc[dsq[i]];
      int sc;
	  
      MMX(i,0) = IMX(i,0) = DMX(i,0) = p7_IMPOSSIBLE;
      XMX(i,p7X_E) = p7_IMPOSSIBLE;
      
      for (k = 1; k < M; k++) 
	{
	  /* match state */
	  sc = ESL_MAX(    MMX(i-1,k-1)   + TSC(p7O_MM,k-1),
			   IMX(i-1,k-1)   + TSC(p7O_IM,k-1));
	  sc = ESL_MAX(sc, DMX(i-1,k-1)   + TSC(p7O_DM,k-1));
	  sc = ESL_MAX(sc, XMX(i-1,p7X_B) + TSC(p7O_BM,k-1)); /* remember, bsc is slid off-by-one */
	  MMX(i,k) = ESL_MAX(sc + MSC(k), p7_IMPOSSIBLE);

	  /* E state update */
	  XMX(i,p7X_E) = ESL_MAX(XMX(i,p7X_E), MMX(i,k) + esc);
	  /* in Viterbi alignments, Dk->E can't win in local mode (and
	   * isn't possible in glocal mode), so don't bother
	   * looking. */
	  
	  /* insert state */
	  sc = ESL_MAX(MMX(i-1,k) + TSC(p7O_MI,k),
		       IMX(i-1,k) + TSC(p7O_II,k));
	  IMX(i,k) = ESL_MAX(sc + ISC(k), p7_IMPOSSIBLE);
	  
	  /* delete state */
	  sc =             MMX(i,k-1) + TSC(p7O_MD,k-1);
	  sc = ESL_MAX(sc, DMX(i,k-1) + TSC(p7O_DD,k-1));
	  DMX(i,k) = ESL_MAX(sc, p7_IMPOSSIBLE);
	}
 
      /* unrolled match state M*/
      sc =             MMX(i-1,M-1)   + TSC(p7O_MM,M-1);
      sc = ESL_MAX(sc, IMX(i-1,M-1)   + TSC(p7O_IM,M-1));
      sc = ESL_MAX(sc, DMX(i-1,M-1)   + TSC(p7O_DM,M-1));
      sc = ESL_MAX(sc, XMX(i-1,p7X_B) + TSC(p7O_BM,M-1));
      MMX(i,M) = ESL_MAX(sc + MSC(M), p7_IMPOSSIBLE);
      
      /* unrolled D state M */
      sc =             MMX(i,M-1) + TSC(p7O_MD,M-1);
      sc = ESL_MAX(sc, DMX(i,M-1) + TSC(p7O_DD,M-1));
      DMX(i,M) = ESL_MAX(sc, p7_IMPOSSIBLE);

      /* E state update: transitions from M_M, D_M score 0 by def'n */
      sc           = ESL_MAX(XMX(i,p7X_E), MMX(i,M));
      XMX(i,p7X_E) = ESL_MAX(sc,           DMX(i,M));

      /* Now the special states. Order is important here.
       * remember, N, C and J emissions are zero score by definition.
       */
      /* J state */
      sc           =             XMX(i-1,p7X_J) + om->xsc[p7O_J][p7O_LOOP];
      sc           = ESL_MAX(sc, XMX(i,  p7X_E) + om->xsc[p7O_E][p7O_LOOP]);
      XMX(i,p7X_J) = ESL_MAX(sc, p7_IMPOSSIBLE);
      
      /* C state */
      sc           =             XMX(i-1,p7X_C) + om->xsc[p7O_C][p7O_LOOP];
      sc           = ESL_MAX(sc, XMX(i,  p7X_E) + om->xsc[p7O_E][p7O_MOVE]);
      XMX(i,p7X_C) = ESL_MAX(sc, p7_IMPOSSIBLE);
      
      /* N state */
      XMX(i,p7X_N) = ESL_MAX(XMX(i-1,p7X_N) + om->xsc[p7O_N][p7O_LOOP], p7_IMPOSSIBLE);
      
      /* B state */
      sc           =             XMX(i,p7X_N) + om->xsc[p7O_N][p7O_MOVE];
      sc           = ESL_MAX(sc, XMX(i,p7X_J) + om->xsc[p7O_J][p7O_MOVE]);
      XMX(i,p7X_B) = ESL_MAX(sc, p7_IMPOSSIBLE);
    }
  
  /* T state (not stored) */
  *ret_sc =  p7_SILO2Lod(XMX(L,p7X_C) + om->xsc[p7O_C][p7O_MOVE]);
  return eslOK;
}


/* Function:  p7_Forward()
 * Synopsis:  Forward algorithm.
 * Incept:    SRE, Tue Jul 17 12:53:59 2007 [Janelia]
 *
 * Purpose:   The Buhler (optimized integer) version of the Forward dynamic
 *            programming algorithm.
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <om>, and DP matrix <ox> allocated for at least <L> by
 *            <om->M> cells; calculate the probability of the sequence
 *            given the model using the Forward algorithm; return the
 *            Forward matrix in <ox>, and the Forward score in <ret_sc>.
 *           
 *            The Viterbi lod score is returned in nats. The caller
 *            needs to subtract a null model lod score, then convert
 *            to bits.
 *           
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq
 *            om     - profile. Does not need to contain any
 *                     reference pointers (alphabet, HMM, or null model)
 *            ox     - DP matrix with room for an MxL alignment
 *            ret_sc - RETURN: Forward lod score in nats
 *           
 * Return:   <eslOK> on success.
 */
int
p7_Forward(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  int const *tsc = om->tsc;
  int      **dp  = ox->dp;
  int       *xmx = ox->xmx;
  int        M   = om->M;
  int        i,k;
  int        esc = p7_IsLocal(om->mode) ? 0 : p7_IMPOSSIBLE;

  /* Initialization of the zero row.  */
  XMX(0,p7X_N) = 0;                                                /* S->N, p=1            */
  XMX(0,p7X_B) = om->xsc[p7O_N][p7O_MOVE];                         /* S->N->B, no N-tail   */
  XMX(0,p7X_E) = XMX(0,p7X_C) = XMX(0,p7X_J) = p7_IMPOSSIBLE;      /* need seq to get here */
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = p7_IMPOSSIBLE;               /* need seq to get here */
  p7_ILogsumInit();

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   */
  for (i = 1; i <= L; i++) 
    {
      int const *rsc = om->rsc[dsq[i]];
      int sc;
	  
      MMX(i,0) = IMX(i,0) = DMX(i,0) = p7_IMPOSSIBLE;
      XMX(i,p7X_E) = p7_IMPOSSIBLE;
      
      for (k = 1; k < M; k++) 
	{
	  /* match state */
	  sc = p7_ILogsum(p7_ILogsum(MMX(i-1,k-1)   + TSC(p7O_MM,k-1), 
				     IMX(i-1,k-1)   + TSC(p7O_IM,k-1)),
			  p7_ILogsum(XMX(i-1,p7G_B) + TSC(p7O_BM,k-1), /* remember, bsc is slid off-by-one */
				     DMX(i-1,k-1)   + TSC(p7O_DM,k-1)));
	  MMX(i,k) = sc + MSC(k);

	  sc = p7_ILogsum(MMX(i-1,k) + TSC(p7O_MI,k),
			  IMX(i-1,k) + TSC(p7O_II,k));
	  IMX(i,k) = sc + ISC(k);

	  /* delete state */
	  DMX(i,k) = p7_ILogsum(MMX(i,k-1) + TSC(p7O_MD,k-1),
				DMX(i,k-1) + TSC(p7O_DD,k-1));

	  /* E state update */
	  XMX(i,p7X_E) = p7_ILogsum(p7_ILogsum(MMX(i,k) + esc,
					       DMX(i,k) + esc),
				               XMX(i,p7X_E));
	}

      /* unrolled match state M_M */
      sc = p7_ILogsum(p7_ILogsum(MMX(i-1,M-1)   + TSC(p7O_MM,M-1), 
				 IMX(i-1,M-1)   + TSC(p7O_IM,M-1)),
		      p7_ILogsum(XMX(i-1,p7G_B) + TSC(p7O_BM,M-1),
				 DMX(i-1,M-1)   + TSC(p7O_DM,M-1)));
      MMX(i,M) = sc + MSC(M);

      /* unrolled delete state D_M */
      DMX(i,M) = p7_ILogsum(MMX(i,M-1) + TSC(p7O_MD,M-1),
			    DMX(i,M-1) + TSC(p7O_DD,M-1));

      /* unrolled E state update */
      XMX(i,p7X_E) = p7_ILogsum(p7_ILogsum(MMX(i,M),
					   DMX(i,M)),
				           XMX(i,p7X_E));
      /* J state */
      XMX(i,p7X_J) = p7_ILogsum(XMX(i-1,p7X_J) + om->xsc[p7O_J][p7O_LOOP],
				XMX(i,  p7X_E) + om->xsc[p7O_E][p7O_LOOP]);
      /* C state */
      XMX(i,p7X_C) = p7_ILogsum(XMX(i-1,p7X_C) + om->xsc[p7O_C][p7O_LOOP],
				XMX(i,  p7X_E) + om->xsc[p7O_E][p7O_MOVE]);
      /* N state */
      XMX(i,p7X_N) = XMX(i-1,p7X_N) + om->xsc[p7O_N][p7O_LOOP];

      /* B state */
      XMX(i,p7X_B) = p7_ILogsum(XMX(i,  p7X_N) + om->xsc[p7O_N][p7O_MOVE],
				XMX(i,  p7X_J) + om->xsc[p7O_J][p7O_MOVE]);
    }

  *ret_sc = p7_SILO2Lod(XMX(L,p7X_C) + om->xsc[p7O_C][p7O_MOVE]);
  return eslOK;
}





/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef p7IMPL_JB_BENCHMARK
/* gcc -o benchmark-jb -g -O2 -I. -L. -I../easel -L../easel -Dp7IMPL_JB_BENCHMARK impl_jb.c -lhmmer -leasel -lm
 * ./benchmark-jb <hmmfile>
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
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-f",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use the Forward algorithm",                      0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual scores",             0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the Buhler integer implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
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

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

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

      if (esl_opt_GetBoolean(go, "-f")) p7_Forward(dsq, L, om, ox, &sc);
      else                              p7_Viterbi(dsq, L, om, ox, &sc);

      p7_bg_NullOne(bg, dsq, L, &nullsc);
      bitscore = (sc - nullsc) / eslCONST_LOG2;
      
      if (esl_opt_GetBoolean(go, "-v")) printf("%.4f bits  (%.4f raw)\n", bitscore, sc); 

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
#endif /*p7IMPL_FP_BENCHMARK*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
