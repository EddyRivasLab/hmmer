/* Reference implementation of the DP algorithms.
 * 
 * This implementation is intended for clarity, not speed; it is
 * deliberately unoptimized. It is provided as a reference
 * implementation of the computationally intensive DP algorithms of
 * HMMER, both as documentation and a regression test target for
 * developing optimized code.
 * 
 * SRE, Tue Jan 30 10:49:43 2007 [at Einstein's in St. Louis]
 * SVN $Id$
 */
#include "p7_config.h"

#include <easel.h>

#include "hmmer.h"


/* Function: p7_Viterbi()
 * Incept:   SRE, Tue Jan 30 10:50:53 2007 [Einstein's in St. Louis]
 * 
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           
 * Args:     dsq    - sequence in digitized form, 1..L
 *           L      - length of dsq
 *           gm     - profile. Does not need to contain any
 *                    reference pointers (alphabet, HMM, or null model)
 *           mx     - DP matrix with room for an MxL alignment
 *           tr     - RETURN: traceback; pass NULL if it's not wanted
 *           ret_sc - RETURN: Viterbi score in bits.
 *           
 * Return:   <eslOK> on success.
 */
int
p7_Viterbi(ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_GMX *mx, P7_TRACE *tr, float *ret_sc)
{
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;

  /* Some convenience */
  xmx = mx->xmx;
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;

  /* Initialization of the zero row.  */
  xmx[0][p7_XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][p7_XMB] = gm->xsc[p7_XTN][p7_MOVE];         /* S->N->B, no N-tail   */
  xmx[0][p7_XME] = xmx[p7_0][p7_XMC] = xmx[0][p7_XMJ] = p7_IMPOSSIBLE;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = p7_IMPOSSIBLE;      /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = impossible for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 1; i <= L; i++) {
    mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;

    for (k = 1; k <= gm->M; k++) {
				/* match state */
      mmx[i][k]  = p7_IMPOSSIBLE;
      if ((sc = mmx[i-1][k-1] + gm->tsc[p7_TMM][k-1]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = imx[i-1][k-1] + gm->tsc[p7_TIM][k-1]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = xmx[i-1][p7_XMB] + gm->bsc[k]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = dmx[i-1][k-1] + gm->tsc[p7_TDM][k-1]) > mmx[i][k])
	mmx[i][k] = sc;
      if (gm->msc[dsq[i]][k] != p7_IMPOSSIBLE) mmx[i][k] += gm->msc[dsq[i]][k];
      else                                     mmx[i][k] =  p7_IMPOSSIBLE;

				/* delete state */
      dmx[i][k] = p7_IMPOSSIBLE;
      if ((sc = mmx[i][k-1] + gm->tsc[p7_TMD][k-1]) > dmx[i][k])
	dmx[i][k] = sc;
      if ((sc = dmx[i][k-1] + gm->tsc[p7_TDD][k-1]) > dmx[i][k])
	dmx[i][k] = sc;

				/* insert state */
      if (k < hmm->M) {
	imx[i][k] = p7_IMPOSSIBLE;
	if ((sc = mmx[i-1][k] + gm->tsc[p7_TMI][k]) > imx[i][k])
	  imx[i][k] = sc;
	if ((sc = imx[i-1][k] + gm->tsc[p7_TII][k]) > imx[i][k])
	  imx[i][k] = sc;
	if (gm->isc[dsq[i]][k] != p7_IMPOSSIBLE) imx[i][k] += gm->isc[dsq[i]][k];
	else                                     imx[i][k] = p7_IMPOSSIBLE;   
      }
    }

    /* Now the special states. Order is important here.
     * remember, N, C and J emissions are zero score by definition.
     */
				/* N state */
    xmx[i][p7_XMN] = p7_IMPOSSIBLE;
    if ((sc = xmx[i-1][p7_XMN] + gm->xsc[p7_XTN][p7_LOOP]) > p7_IMPOSSIBLE)
      xmx[i][p7_XMN] = sc;

				/* E state */
    /* We don't need to check D_M->E transition; it provably cannot
     * be used in a Viterbi path in HMMER3 parameterization */
    xmx[i][p7_XME] = p7_IMPOSSIBLE;
    for (k = 1; k <= gm->M; k++)
      if ((sc =  mmx[i][k] + gm->esc[k]) > xmx[i][p7_XME])
	xmx[i][p7_XME] = sc;
				/* J state */
    xmx[i][p7_XMJ] = p7_IMPOSSIBLE;
    if ((sc = xmx[i-1][p7_XMJ] + gm->xsc[p7_XTJ][p7_LOOP]) > p7_IMPOSSIBLE)
      xmx[i][p7_XMJ] = sc;
    if ((sc = xmx[i][p7_XME]   + gm->xsc[p7_XTE][p7_LOOP]) > xmx[i][p7_XMJ])
      xmx[i][p7_XMJ] = sc;

				/* B state */
    xmx[i][p7_XMB] = p7_IMPOSSIBLE;
    if ((sc = xmx[i][XMN] + gm->xsc[p7_XTN][p7_MOVE]) > p7_IMPOSSIBLE)
      xmx[i][p7_XMB] = sc;
    if ((sc = xmx[i][XMJ] + gm->xsc[p7_XTJ][p7_MOVE]) > xmx[i][p7_XMB])
      xmx[i][p7_XMB] = sc;

				/* C state */
    xmx[i][p7_XMC] = p7_IMPOSSIBLE;
    if ((sc = xmx[i-1][p7_XMC] + gm->xsc[p7_XTC][p7_LOOP]) > p7_IMPOSSIBLE)
      xmx[i][p7_XMC] = sc;
    if ((sc = xmx[i][p7_XME] + gm->xsc[p7_XTE][p7_MOVE]) > xmx[i][p7_XMC])
      xmx[i][p7_XMC] = sc;
  }
				/* T state (not stored) */
  sc = xmx[L][p7_XMC] + hmm->xsc[p7_XTC][p7_MOVE];

  if (ret_tr != NULL) {
    P7ViterbiTrace(hmm, dsq, L, mx, &tr);
    *ret_tr = tr;
  }

  return p7_Score2Output(sc);	/* the total Viterbi score, in bits */
}

/* Function: p7_ViterbiTrace()
 * Incept:   SRE, Thu Feb  1 10:25:56 2007 [UA 8018 St. Louis to Dulles]
 * 
 * Purpose:  Traceback of a Viterbi matrix: retrieval 
 *           of optimum alignment.
 *           
 * Args:     dsq    - sequence aligned to (digital form) 1..L 
 *           L      - length of dsq gm    
 *           gm     - profile model; does not need any ref ptrs
 *           mx     - the matrix to trace, L x M
 *           tr     - storage for the recovered traceback.
 *           
 * Return:   <eslOK> on success.
 *           <eslFAIL> if even the optimal path has zero probability;
 *           in this case, the trace is set blank (<tr->N = 0>).
 */
int
p7_ViterbiTrace(ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_GMX *mx, P7_TRACE *tr)
{
  int status;
  int tpos;			/* position in trace */
  char st;			/* last state type we assigned */
  int i;			/* position in seq (1..L) */
  int k;			/* position in model (1..M) */
  int **xmx, **mmx, **imx, **dmx;
  int sc;			/* temp var for pre-emission score */

  if ((status = p7_trace_Reuse(tr)) != eslOK) goto ERROR;
  xmx = mx->xmx;
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;

  /* Initialization of trace
   * We do it back to front; ReverseTrace() is called later.
   */
  if ((status = p7_trace_Append(tr, p7_STT, 0, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_Append(tr, p7_STC, 0, 0)) != eslOK) goto ERROR;
  i    = L;			/* next i (seq pos) to assign    */
  k    = gm->M;			/* next k (model node) to assign */
  st   = p7_STC;		/* last state type we assigned   */

  /* Traceback
   */
  while (st != p7_STS) {
    switch (st) {
    case p7_STC:			/* C comes from C, E */
      if      (xmx[i][p7_XMC] <= p7_IMPOSSIBLE) { status = eslFAIL; goto ERROR; }
      else if (xmx[i][p7_XMC] == xmx[i-1][p7_XMC] + gm->xsc[p7_XTC][p7_LOOP])
	{
	  if ((status = p7_trace_Append(tr, p7_STC, 0, i)) != eslOK) goto ERROR;
	  tr->statetype[tpos] = STC;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;    /* note convention adherence: */
	  tr->pos[tpos-1]     = i--;  /* first C doesn't emit       */
	}
      else if (xmx[i][XMC] == xmx[i][XME] + hmm->xsc[XTE][MOVE])
	{
	  tr->statetype[tpos] = STE;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0; /* E is a nonemitter */
	}
      
      else Die("Traceback failed.");
      break;



    case p7_STM:			/* M connects from i-1,k-1, or B */
      sc = mmx[i][k] - gm->msc[dsq[i]][k];
      if (sc <= p7_IMPOSSIBLE) ESL_XFAIL(eslFAIL, NULL, "impossible trace");
      else if (sc == xmx[i][p7_XMB] + gm->bsc[k+1])
	{
				/* Check for wing unfolding */
	  if (Prob2Score(hmm->begin[k+1], hmm->p1) + 1 * INTSCALE <= hmm->bsc[k+1])
	    while (k > 0)
	      {
		tr->statetype[tpos] = STD;
		tr->nodeidx[tpos]   = k--;
		tr->pos[tpos]       = 0;
		tpos++;
		if (tpos == curralloc) 
		  {				/* grow trace if necessary  */
		    curralloc += N;
		    P7ReallocTrace(tr, curralloc);
		  }
	      }

	  tr->statetype[tpos] = STB;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	}
      else if (sc == mmx[i][k] + hmm->tsc[TMM][k])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == imx[i][k] + hmm->tsc[TIM][k])
	{
	  tr->statetype[tpos] = STI;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == dmx[i][k] + hmm->tsc[TDM][k])
	{
	  tr->statetype[tpos] = STD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	}
      else
	Die("traceback failed");
      break;

    case STD:			/* D connects from M,D */
      if (dmx[i][k+1] <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (dmx[i][k+1] == mmx[i][k] + hmm->tsc[TMD][k])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (dmx[i][k+1] == dmx[i][k] + hmm->tsc[TDD][k]) 
	{
	  tr->statetype[tpos] = STD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	}
      else Die("traceback failed");
      break;

    case STI:			/* I connects from M,I */
      sc = imx[i+1][k] - hmm->isc[dsq[i+1]][k];
      if (sc <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (sc == mmx[i][k] + hmm->tsc[TMI][k])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == imx[i][k] + hmm->tsc[TII][k])
	{
	  tr->statetype[tpos] = STI;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = i--;
	}
      else Die("traceback failed");
      break;

    case STN:			/* N connects from S, N */
      if (i == 0 && xmx[i][XMN] == 0)
	{
	  tr->statetype[tpos] = STS;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	}
      else if (i > 0 && xmx[i+1][XMN] == xmx[i][XMN] + hmm->xsc[XTN][LOOP])
	{
	  tr->statetype[tpos] = STN;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;    /* note convention adherence:  */
	  tr->pos[tpos-1]     = i--;  /* first N doesn't emit        */
	}
      else Die("traceback failed");
      break;

    case STB:			/* B connects from N, J */
      if (xmx[i][XMB] <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (xmx[i][XMB] == xmx[i][XMN] + hmm->xsc[XTN][MOVE])
	{
	  tr->statetype[tpos] = STN;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	}
      else if (xmx[i][XMB] == xmx[i][XMJ] + hmm->xsc[XTJ][MOVE])
	{
	  tr->statetype[tpos] = STJ;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	}

      else Die("traceback failed");
      break;

    case STE:			/* E connects from any M state. k set here */
      if (xmx[i][XME] <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      for (k = hmm->M; k >= 1; k--)
	if (xmx[i][XME] == mmx[i][k] + hmm->esc[k])
	  {
				/* check for wing unfolding */
	    if (Prob2Score(hmm->end[k], 1.) + 1*INTSCALE <=  hmm->esc[k])
	      {
		int dk;		/* need a tmp k while moving thru delete wing */
		for (dk = hmm->M; dk > k; dk--)
		  {
		    tr->statetype[tpos] = STD;
		    tr->nodeidx[tpos]   = dk;
		    tr->pos[tpos]       = 0;
		    tpos++;
		    if (tpos == curralloc) 
		      {				/* grow trace if necessary  */
			curralloc += N;
			P7ReallocTrace(tr, curralloc);
		      }
		  }
	      }

	    tr->statetype[tpos] = STM;
	    tr->nodeidx[tpos]   = k--;
	    tr->pos[tpos]       = i--;
	    break;
	  }
      if (k < 0) Die("traceback failed");
      break;

    case STJ:			/* J connects from E, J */
      if (xmx[i][XMJ] <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (xmx[i][XMJ] == xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP])
	{
	  tr->statetype[tpos] = STJ;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;    /* note convention adherence: */
	  tr->pos[tpos-1]     = i--;  /* first J doesn't emit       */
	}
      else if (xmx[i][XMJ] == xmx[i][XME] + hmm->xsc[XTE][LOOP])
	{
	  tr->statetype[tpos] = STE;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0; /* E is a nonemitter */
	}

      else Die("Traceback failed.");
      break;

    default:
      Die("traceback failed");

    } /* end switch over statetype[tpos-1] */
    
    tpos++;
    if (tpos == curralloc) 
      {				/* grow trace if necessary  */
	curralloc += N;
	P7ReallocTrace(tr, curralloc);
      }

  } /* end traceback, at S state; tpos == tlen now */
  tr->tlen = tpos;
  P7ReverseTrace(tr);
  *ret_tr = tr;
}


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
