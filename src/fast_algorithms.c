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
 * and maintainability, not for speed. The implementations here
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

