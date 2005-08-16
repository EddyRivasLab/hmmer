#include "structs.h"
#include "funcs.h"

/* Function: Viterbi()
 *
 * Note:     This used to be defined as P7Viterbi() in core_algorithms.c, but I
 *           moved it here and renamed it so that it fit into the new
 *           architecture.  - CRS 20 June 2005
 * 
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           Identical to Forward() except that max's
 *           replace sum's. 
 *           
 *           This is the slower, more understandable version
 *           of P7Viterbi(). The default version in fastfuncs.c
 *           is portably optimized and more difficult to understand;
 *           the ALTIVEC version in altivecfuncs.c is vectorized
 *           with Altivec-specific code, and is pretty opaque.
 *           
 *           This function is not enabled by default; it is only
 *           activated by -DSLOW at compile time.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - reused DP matrix
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float
Viterbi(unsigned char *dsq, int L, struct plan7_s *hmm, cust_dpmatrix_s *mx, 
	struct p7trace_s **ret_tr)
{
  struct p7trace_s  *tr;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;  

  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 
  ResizeDPMatrix(mx, L, hmm->M);
  xmx = mx->xmx;
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;

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
    mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;

    for (k = 1; k <= hmm->M; k++) {
				/* match state */
      mmx[i][k]  = -INFTY;
      if ((sc = mmx[i-1][k-1] + hmm->tsc[TMM][k-1]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = imx[i-1][k-1] + hmm->tsc[TIM][k-1]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = xmx[i-1][XMB] + hmm->bsc[k]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = dmx[i-1][k-1] + hmm->tsc[TDM][k-1]) > mmx[i][k])
	mmx[i][k] = sc;
      if (hmm->msc[dsq[i]][k] != -INFTY) mmx[i][k] += hmm->msc[dsq[i]][k];
      else                                     mmx[i][k] = -INFTY;

				/* delete state */
      dmx[i][k] = -INFTY;
      if ((sc = mmx[i][k-1] + hmm->tsc[TMD][k-1]) > dmx[i][k])
	dmx[i][k] = sc;
      if ((sc = dmx[i][k-1] + hmm->tsc[TDD][k-1]) > dmx[i][k])
	dmx[i][k] = sc;

				/* insert state */
      if (k < hmm->M) {
	imx[i][k] = -INFTY;
	if ((sc = mmx[i-1][k] + hmm->tsc[TMI][k]) > imx[i][k])
	  imx[i][k] = sc;
	if ((sc = imx[i-1][k] + hmm->tsc[TII][k]) > imx[i][k])
	  imx[i][k] = sc;
	if (hmm->isc[dsq[i]][k] != -INFTY) imx[i][k] += hmm->isc[dsq[i]][k];
	else                                    imx[i][k] = -INFTY;   
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
    xmx[i][XME] = -INFTY;
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  mmx[i][k] + hmm->esc[k]) > xmx[i][XME])
	xmx[i][XME] = sc;
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
    ViterbiTrace(hmm, dsq, L, mx, &tr);
    *ret_tr = tr;
  }

  return Scorify(sc);		/* the total Viterbi score. */
}

