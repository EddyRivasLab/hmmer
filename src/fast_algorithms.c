/************************************************************
 * @LICENSE@
 ************************************************************/

/* fast_algorithms.c
 * SRE, Sun Nov 10 08:54:48 2002 [AA 3080, Denver to StL]
 * CVS $Id$
 * 
 * Optimized routines to replace slower implementations in core_algorithms.c.
 * 
 * The routines in core_algorithms.c are designed for clarity
 * and maintainability, not for speed. Implementations here
 * are designed for speed, not clarity. If you're trying to 
 * understand the code, or optimize for a specific platform,
 * you are probably better off looking at core_algorithms.c.
 * 
 * Which implementation is used is controlled by ifdef's. One and
 * only one of the following must be on:
 * 
 * (default, no flag set)
 *   all implementations provided by core_algorithms.c
 * SPEEDY
 *   enable generic (portable, platform-independent) optimized implementations
 *   P7Viterbi() is replaced.
 * INTEL_VECTORIZED
 *   enable code optimized for icc vectorization (Intel C compiler)
 *   P7Viterbi() is replaced.
 * ALTIVEC
 *   enable Erik Lindahl's Altivec code for Macintosh OSX
 */

#include "structs.h"
#include "config.h"
#include "funcs.h"
#include "squid.h"

#ifdef SPEEDY
/* Function: P7Viterbi() - portably optimized version
 * Incept:   SRE, Sun Nov 10 09:04:21 2002 [Boulder]          
 * 
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           Derived from core_algorithms.c:P7Viterbi().
 *           
 *  Optimized by replacing all array lookups with pointer
 *  arithmetic in the innermost loop of the DP.
 *  
 *  Some of it is dependent on order of evaluation,
 *  so don't rearrange any code without understanding.
 *  
 *  Particularly sneaky is the use of tp to step through hmm->tsc[][]:
 *    Take advantage of the fact that all the memory is alloc'ed
 *    in one chunk, in tsc[0]. We can step through the whole 2D array
 *    in 1D.
 *       |------- k-1 ------|   |------ k ---------|
 *       MM MI MD IM II DM DD   MM MI MD IM II DM DD
 *    
 *    so if we start with a ptr to [k-1][TMM]:
 *      +3 to get to [k-1][TIM] 
 *      +2 to get to [k-1][TDM]
 *      -3 to get to [k-1][TMD]
 *      +4 to get to [k-1][TDD]
 *      +2 to get to [k][TMI]
 *      +3 to get to [k][TII]
 *      -4 to set to [k-1][TMM] for next iteration of k.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float
P7Viterbi(char *dsq, int L, struct plan7_s *hmm, struct p7trace_s **ret_tr)
{
  struct dpmatrix_s *mx;
  struct p7trace_s  *tr;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k,z;
  int   sc;
  int  *mc, *dc, *ic;        /* pointers to rows of mmx, dmx, imx */
  int  *ms, *is;             /* pointers to msc[i], isc[i] */
  int  *mpp, *mpc, *ip;      /* ptrs to mmx[i-1], mmx[i], imx[i-1] */
  int  *bp;		     /* ptr into bsc[] */
  int  *ep;                  /* ptr into esc[] */
  int   xmb;		     /* value of xmx[i-1][XMB] */
  int   xme;                 /* max for xmx[i][XME] */
  int  *dpp, *dpc;           /* ptr into dmx[i-1] (previous) and dmx[i] (current) */
  int  *tpmm, *tpmi, *tpmd, *tpim, *tpii, *tpdm, *tpdd; /* ptrs into tsc */
  
  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 
  mx = AllocPlan7Matrix(L+1, hmm->M, &xmx, &mmx, &imx, &dmx);

  /* Initialization of the zero row.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 1; i <= L; i++) {
    mc = mmx[i];  *mc = -INFTY;  mc++;
    dc = dmx[i];  *dc = -INFTY;  dc++;
    ic = imx[i];  *ic = -INFTY;  ic++;

    ms = hmm->msc[(int) dsq[i]] + 1;  /* start vector at hmm->msc[dsq[i]][1] */
    is = hmm->isc[(int) dsq[i]] + 1;  /* start vector at hmm->isc[dsq[i]][1] */

    mpp = mmx[i-1];
    mpc = mmx[i];
    ip  = imx[i-1];
    bp  = hmm->bsc + 1;
    xmb = xmx[i-1][XMB];
    dpp = dmx[i-1];
    dpc = dmx[i];
    tpmm = hmm->tsc[TMM];
    tpmi = hmm->tsc[TMI] + 1;
    tpmd = hmm->tsc[TMD];
    tpim = hmm->tsc[TIM];
    tpii = hmm->tsc[TII] + 1;
    tpdm = hmm->tsc[TDM];
    tpdd = hmm->tsc[TDD];

    for (k = 1; k <= hmm->M; k++) {
      *mc = *mpp++ + *tpmm++;
      if ((sc = *ip++  + *tpim++) > *mc)  *mc = sc;
      if ((sc = *dpp++ + *tpdm++) > *mc)  *mc = sc;
      if ((sc = xmb  + *bp++)     > *mc)  *mc = sc; 
      *mc += *ms++;
      if (*mc < -INFTY) *mc = -INFTY;
      mc++;

      *dc = *mpc++ + *tpmd++;                         
      if ((sc = *dpc++ + *tpdd++) > *dc) *dc = sc;
      if (*dc < -INFTY) *dc = -INFTY;
      dc++;

      if (k < hmm->M) {
	*ic = *mpp + *tpmi++;
	if ((sc = *ip + *tpii++) > *ic) *ic = sc; 
	*ic += *is++;
	if (*ic < -INFTY) *ic = -INFTY;	
	ic++;
      }
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
				/* N state */
    xmx[i][XMN] = -INFTY;
    if ((sc = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[i][XMN] = sc;

				/* E state */
    xme = -INFTY;
    mpc = mmx[i] + 1;
    ep  = hmm->esc + 1;
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  *mpc++ + *ep++) > xme) xme = sc; 
    xmx[i][XME] = xme;
				/* J state */
    xmx[i][XMJ] = -INFTY;
    if ((sc = xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP]) > -INFTY)
      xmx[i][XMJ] = sc;
    if ((sc = xmx[i][XME]   + hmm->xsc[XTE][LOOP]) > xmx[i][XMJ])
      xmx[i][XMJ] = sc;

				/* B state */
    xmx[i][XMB] = -INFTY;
    if ((sc = xmx[i][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
      xmx[i][XMB] = sc;
    if ((sc = xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[i][XMB])
      xmx[i][XMB] = sc;

				/* C state */
    xmx[i][XMC] = -INFTY;
    if ((sc = xmx[i-1][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[i][XMC] = sc;
    if ((sc = xmx[i][XME] + hmm->xsc[XTE][MOVE]) > xmx[i][XMC])
      xmx[i][XMC] = sc;
  }
				/* T state (not stored) */
  sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

  if (ret_tr != NULL) {
    P7ViterbiTrace(hmm, dsq, L, mx, &tr);
    *ret_tr = tr;
  }

  FreePlan7Matrix(mx);
  return Scorify(sc);		/* the total Viterbi score. */
}
#endif /*SPEEDY*/


#ifdef INTEL_VECTORIZED
/* Function: P7Viterbi() - icc optimized version
 * Incept:   SRE, Fri Nov 15 13:14:33 2002 [St. Louis]
 * 
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           Derived from fast_algorithms.c:P7Viterbi() (-DSPEEDY version)
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float
P7Viterbi(char *dsq, int L, struct plan7_s *hmm, struct p7trace_s **ret_tr)
{
  struct dpmatrix_s *mx;
  struct p7trace_s  *tr;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k,z;
  int   sc;
  int  *mc, *dc, *ic;        /* pointers to rows of mmx, dmx, imx */
  int  *ms, *is;             /* pointers to msc[i], isc[i] */
  int  *mpp, *mpc, *ip;      /* ptrs to mmx[i-1], mmx[i], imx[i-1] */
  int  *bp;		     /* ptr into bsc[] */
  int  *ep;                  /* ptr into esc[] */
  int   xmb;		     /* value of xmx[i-1][XMB] */
  int   xme;                 /* max for xmx[i][XME] */
  int  *dpp, *dpc;           /* ptr into dmx[i-1] (previous) and dmx[i] (current) */
  int  *tpmm, *tpmi, *tpmd, *tpim, *tpii, *tpdm, *tpdd; /* ptrs into tsc */
  int   M;
  int   neginfty;
  
  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 
  mx = AllocPlan7Matrix(L+1, hmm->M, &xmx, &mmx, &imx, &dmx);

  /* Initialization of the zero row.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Initializations that help icc vectorize.
   */
  M        = hmm->M;

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */

  tpmm  = hmm->tsc[TMM];
  tpim  = hmm->tsc[TIM];
  tpdm  = hmm->tsc[TDM];
  tpmd  = hmm->tsc[TMD];
  tpdd  = hmm->tsc[TDD];
  tpmi  = hmm->tsc[TMI];
  tpii  = hmm->tsc[TII];
  bp    = hmm->bsc;
  for (i = 1; i <= L; i++) {
    mc    = mmx[i];    
    dc    = dmx[i];
    ic    = imx[i];
    mpp   = mmx[i-1];
    dpp   = dmx[i-1];
    ip    = imx[i-1];
    xmb   = xmx[i-1][XMB];
    ms    = hmm->msc[(int) dsq[i]];
    is    = hmm->isc[(int) dsq[i]];
    mc[0] = -INFTY;
    dc[0] = -INFTY;
    ic[0] = -INFTY;

    for (k = 1; k <= M; k++) {
      mc[k] = mpp[k-1]   + tpmm[k-1];
      if ((sc = ip[k-1]  + tpim[k-1]) > mc[k])  mc[k] = sc;
      if ((sc = dpp[k-1] + tpdm[k-1]) > mc[k])  mc[k] = sc;
      if ((sc = xmb  + bp[k])         > mc[k])  mc[k] = sc; 
      mc[k] += ms[k];
      if (mc[k] < -INFTY) mc[k] = -INFTY;  

      dc[k] = dc[k-1] + tpdd[k-1];
      if ((sc = mc[k-1] + tpmd[k-1]) > dc[k]) dc[k] = sc;
      if (dc[k] < -INFTY) dc[k] = -INFTY;  

      if (k < M) {
	ic[k] = mpp[k] + tpmi[k];
	if ((sc = ip[k] + tpii[k]) > ic[k]) ic[k] = sc; 
	ic[k] += is[k];
	if (ic[k] < -INFTY) ic[k] = -INFTY; 
      }
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
				/* N state */
    xmx[i][XMN] = -INFTY;
    if ((sc = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[i][XMN] = sc;

				/* E state */
    xme = -INFTY;
    mpc = mmx[i];
    ep  = hmm->esc;
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  mpc[k] + ep[k]) > xme) xme = sc; 
    xmx[i][XME] = xme;
				/* J state */
    xmx[i][XMJ] = -INFTY;
    if ((sc = xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP]) > -INFTY)
      xmx[i][XMJ] = sc;
    if ((sc = xmx[i][XME]   + hmm->xsc[XTE][LOOP]) > xmx[i][XMJ])
      xmx[i][XMJ] = sc;

				/* B state */
    xmx[i][XMB] = -INFTY;
    if ((sc = xmx[i][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
      xmx[i][XMB] = sc;
    if ((sc = xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[i][XMB])
      xmx[i][XMB] = sc;

				/* C state */
    xmx[i][XMC] = -INFTY;
    if ((sc = xmx[i-1][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[i][XMC] = sc;
    if ((sc = xmx[i][XME] + hmm->xsc[XTE][MOVE]) > xmx[i][XMC])
      xmx[i][XMC] = sc;
  }
				/* T state (not stored) */
  sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

  if (ret_tr != NULL) {
    P7ViterbiTrace(hmm, dsq, L, mx, &tr);
    *ret_tr = tr;
  }

  FreePlan7Matrix(mx);
  return Scorify(sc);		/* the total Viterbi score. */
}
#endif /*INTEL_VECTORIZED*/


#ifdef ALTIVEC
void
AllocPlan7Body(struct plan7_s *hmm, int M) 
{
  int k, x;

  hmm->M = M;

  hmm->rf     = MallocOrDie ((M+2) * sizeof(char));
  hmm->cs     = MallocOrDie ((M+2) * sizeof(char));
  hmm->ca     = MallocOrDie ((M+2) * sizeof(char));
  hmm->map    = MallocOrDie ((M+1) * sizeof(int));

  hmm->t      = MallocOrDie (M     *           sizeof(float *));
  hmm->tsc    = MallocOrDie (7     *           sizeof(int *));
  hmm->mat    = MallocOrDie ((M+1) *           sizeof(float *));
  hmm->ins    = MallocOrDie (M     *           sizeof(float *));
  hmm->msc    = MallocOrDie (MAXCODE   *       sizeof(int *));
  hmm->isc    = MallocOrDie (MAXCODE   *       sizeof(int *)); 
  
  hmm->t[0]   = MallocOrDie ((7*M)     *       sizeof(float));
  /* Allocate extra memory so tsc[TMM,TIM,TDM,TMD,TDD] start on the 
   * 16-byte cache boundary, and tsc[TMI,TII] start
   * 12 bytes offset from the boundary. 
   */
  hmm->tsc_mem = MallocOrDie (((7*(M+16)))  *   sizeof(int));
  hmm->mat[0] = MallocOrDie ((MAXABET*(M+1)) * sizeof(float));
  hmm->ins[0] = MallocOrDie ((MAXABET*M) *     sizeof(float));
  /* Allocate extra mem. to make sure all members of msc,isc start
   * on 12-byte offsets from cache boundary.
   */
  hmm->msc_mem = MallocOrDie ((MAXCODE*(M+1+16)) * sizeof(int));
  hmm->isc_mem = MallocOrDie ((MAXCODE*(M+16)) *   sizeof(int));

  /* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   */
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * MAXABET;
    if (k < M) {
      hmm->ins[k] = hmm->ins[0] + k * MAXABET;
      hmm->t[k]   = hmm->t[0]   + k * 7;
    }
  }
  
  /* align tsc pointers */
  hmm->tsc[TMM] = (int *) (((((unsigned long int) hmm->tsc_mem) + 15) & (~0xf)));
  hmm->tsc[TMI] = (int *) (((((unsigned long int) hmm->tsc_mem) + (M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  hmm->tsc[TMD] = (int *) (((((unsigned long int) hmm->tsc_mem) + 2*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TIM] = (int *) (((((unsigned long int) hmm->tsc_mem) + 3*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TII] = (int *) (((((unsigned long int) hmm->tsc_mem) + 4*(M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  hmm->tsc[TDM] = (int *) (((((unsigned long int) hmm->tsc_mem) + 5*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TDD] = (int *) (((((unsigned long int) hmm->tsc_mem) + 6*(M+12)*sizeof(int) + 15) & (~0xf)));

  for (x = 0; x < MAXCODE; x++) {
    hmm->msc[x] = (int *) (((((unsigned long int)hmm->msc_mem) + x*(M+1+12)*sizeof(int) + 15) & (~0xf)) + 12);
    hmm->isc[x] = (int *) (((((unsigned long int)hmm->isc_mem) + x*(M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  }
  /* tsc[0] is used as a boundary condition sometimes [Viterbi()],
   * so set to -inf always.
   */
  for (x = 0; x < 7; x++)
    hmm->tsc[x][0] = -INFTY;

  hmm->begin  = MallocOrDie  ((M+1) * sizeof(float));
  hmm->bsc_mem= MallocOrDie  ((M+1+12) * sizeof(int));
  hmm->end    = MallocOrDie  ((M+1) * sizeof(float));
  hmm->esc_mem= MallocOrDie  ((M+1+12) * sizeof(int));

  hmm->bsc = (int *) (((((unsigned long int) hmm->bsc_mem) + 15) & (~0xf)) + 12);
  hmm->esc = (int *) (((((unsigned long int) hmm->esc_mem) + 15) & (~0xf)) + 12);
  
  return;
}  
#endif /*ALTIVEC*/
