#include "config.h"
#include "squidconf.h"

#include "squid.h"
#include "structs.h"
#include "funcs.h"


/* Function: P7ViterbiSize()
 * Date:     SRE, Fri Mar  6 15:13:20 1998 [St. Louis]
 *
 * Note:      This was originally defined as in core_algorithms.c, but it is 
 *            mainly a helper function for ViterbiSpaceOK(), which is 
 *            implementation-dependent.  So, I moved it here because it
 *            should really only be used for the default implementation.
 *	        - CRS 19 July 2005
 *
 * Purpose:  Returns the ballpark predicted memory requirement for a 
 *           Viterbi() alignment, in MB. 
 *           
 *           Currently L must fit in an int (< 2 GB), but we have
 *           to deal with LM > 2 GB - e.g. watch out for overflow, do
 *           the whole calculation in floating point. Bug here detected
 *           in 2.1.1 by David Harper, Sanger Centre.
 *                                    
 * Args:     L  - length of sequence
 *           M  - length of HMM       
 *
 * Returns:  # of MB
 */
int
P7ViterbiSize(int L, int M)
{
  float Mbytes;

  /* We're excessively precise here, but it doesn't cost
   * us anything to be pedantic. The four terms are:
   *   1. the matrix structure itself;
   *   2. the O(NM) main matrix (this dominates!)
   *   3. ptrs into the rows of the matrix
   *   4. storage for 5 special states. (xmx)
   */
  Mbytes =  (float) sizeof(cust_dpmatrix_s); 
  Mbytes += 3. * (float) (L+1) * (float) (M+2) * (float) sizeof(int); 
  Mbytes += 4. * (float) (L+1) * (float) sizeof(int *); 
  Mbytes += 5. * (float) (L+1) * (float) sizeof(int);
  Mbytes /= 1048576.;
  return (int) Mbytes;
}

/* Function: P7ViterbiTrace()
 * Date:     SRE, Sat Aug 23 10:30:11 1997 (St. Louis Lambert Field) 
 * 
 * Note:     This was originally defined in core_algorithms.c, but it is 
 *           really an implementation dependent helper function to Viterbi.  
 *           So I moved it here with the other implementation-dependent
 *           functions. - CRS 15 July 2005
 *
 * Purpose:  Traceback of a Viterbi matrix: i.e. retrieval 
 *           of optimum alignment.
 *           
 * Args:     hmm    - hmm, log odds form, used to make mx
 *           dsq    - sequence aligned to (digital form) 1..N  
 *           N      - length of seq
 *           mx     - the matrix to trace back in, N x hmm->M
 *           ret_tr - RETURN: traceback.
 *           
 * Return:   (void)
 *           ret_tr is allocated here. Free using P7FreeTrace().
 */
void
P7ViterbiTrace(struct plan7_s *hmm, unsigned char *dsq, int N,
	       cust_dpmatrix_s *mx, struct p7trace_s **ret_tr)
{
  struct p7trace_s *tr;
  int curralloc;		/* current allocated length of trace */
  int tpos;			/* position in trace */
  int i;			/* position in seq (1..N) */
  int k;			/* position in model (1..M) */
  int **xmx, **mmx, **imx, **dmx;
  int sc;			/* temp var for pre-emission score */

  /* Overallocate for the trace.
   * S-N-B- ... - E-C-T  : 6 states + N is minimum trace;
   * add N more as buffer.             
   */
  curralloc = N * 2 + 6; 
  P7AllocTrace(curralloc, &tr);

  xmx = mx->xmx;
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;

  /* Initialization of trace
   * We do it back to front; ReverseTrace() is called later.
   */
  tr->statetype[0] = STT;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tr->statetype[1] = STC;
  tr->nodeidx[1]   = 0;
  tr->pos[1]       = 0;
  tpos = 2;
  i    = N;			/* current i (seq pos) we're trying to assign */

  /* Traceback
   */
  while (tr->statetype[tpos-1] != STS) {
    switch (tr->statetype[tpos-1]) {
    case STM:			/* M connects from i-1,k-1, or B */
      sc = mmx[i+1][k+1] - hmm->msc[dsq[i+1]][k+1];
      if (sc <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (sc == xmx[i][XMB] + hmm->bsc[k+1])
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

    case STC:			/* C comes from C, E */
      if (xmx[i][XMC] <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (xmx[i][XMC] == xmx[i-1][XMC] + hmm->xsc[XTC][LOOP])
	{
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



/* Function:  ViterbiSpaceOK()
 * Incept:    SRE, Wed Oct  1 12:53:13 2003 [St. Louis]
 *
 * Note:      This was originally defined as P7ViterbiSpaceOK() in 
 *            core_algorithms.c, but it is implementation-dependent since it
 *	      accesses the customized dpmatrix structure.  Hence, I moved it
 *	      here and renamed it so that it fit into the new architecture.
 *	        - CRS 19 July 2005
 *
 * Purpose:   Returns TRUE if the existing matrix allocation
 *            is already sufficient to hold the requested MxN, or if
 *            the matrix can be expanded in M and/or N without
 *            exceeding RAMLIMIT megabytes. 
 *            
 *            This gets called anytime we're about to do Viterbi().
 *            If it returns FALSE, we switch into the appropriate
 *            small-memory alternative: P7SmallViterbi() or
 *            P7WeeViterbi().
 *            
 *            Checking the DP problem size against P7ViterbiSize()
 *            is not enough, because the DP matrix may be already
 *            allocated in MxN. For example, if we're already
 *            allocated to maxM,maxN of 1447x979, and we try to
 *            do a problem of MxN=12x125000, P7ViterbiSize() may
 *            think that's fine - but when we resize, we can only
 *            grow, so we'll expand to 1447x125000, which is 
 *            likely over the RAMLIMIT. [bug#h26; xref SLT7 p.122]
 *
 * Args:      L  - length of sequence
 *            M  - length of HMM
 *            mx - an allocated model
 *
 * Returns:   TRUE if we can run Viterbi(); FALSE if we need
 *            to use a small memory variant.
 *
 * Xref:      STL7 p.122.
 */
int
ViterbiSpaceOK(int L, int M, cust_dpmatrix_s *mx)
{
  int newM;
  int newN;

  if (M <= mx->maxM && L <= mx->maxN) return TRUE;

  if (M > mx->maxM) newM = M + mx->padM; else newM = mx->maxM;
  if (L > mx->maxN) newN = L + mx->padN; else newN = mx->maxN;

  if (P7ViterbiSize(newN, newM) <= RAMLIMIT)
    return TRUE;
  else
    return FALSE;
}


/* Function: Backward()
 *
 * Note:     This was originally defined as P7Backward() in postprob.c, but it
 *           is implementation dependent, since it accesses the customized 
 *           dpmatrix structure.  So I moved it here and renamed it with a more
 *           generic name that fits into the new architecture.
 *             - CRS 15 July 2005
 * 
 * Purpose:  The Backward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float
Backward(unsigned char *dsq, int L, struct plan7_s *hmm, cust_dpmatrix_s **ret_mx)
{
  cust_dpmatrix_s *mx;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;

  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 
  mx = CreateDPMatrix(L+1, hmm->M, 0, 0);
  xmx = mx->xmx;
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;


  /* Initialization of the L row.
   * Note that xmx[i][stS] = xmx[i][stN] by definition for all i,
   *    so stS need not be calculated in backward DP matrices.
   */
  xmx[L][XMC] = hmm->xsc[XTC][MOVE];                 /* C<-T                          */
  xmx[L][XME] = xmx[L][XMC] + hmm->xsc[XTE][MOVE];   /* E<-C, no C-tail               */
  xmx[L][XMJ] = xmx[L][XMB] = xmx[L][XMN] = -INFTY;  /* need seq to get out from here */
  for (k = hmm->M; k >= 1; k--) {
    mmx[L][k] = xmx[L][XME] + hmm->esc[k];           /* M<-E ...                      */
    mmx[L][k] += hmm->msc[dsq[L]][k];                /* ... + emitted match symbol    */
    imx[L][k] = dmx[L][k] = -INFTY;                  /* need seq to get out from here */
  }

  /* Recursion. Done as a pull.
   * Note slightly wasteful boundary conditions:
   *    M_M precalculated, D_M set to -INFTY,
   *    D_1 wastefully calculated.
   * Scores for transitions to D_M also have to be hacked to -INFTY,
   * as Plan7Logoddsify does not do this for us (I think? - ihh).
   */
  hmm->tsc[TDD][hmm->M-1] = hmm->tsc[TMD][hmm->M-1] = -INFTY;    /* no D_M state -- HACK -- should be in Plan7Logoddsify */
  for (i = L-1; i >= 0; i--)
    {
      /* Do the special states first.
       * remember, C, N and J emissions are zero score by definition
       */
      xmx[i][XMC] = xmx[i+1][XMC] + hmm->xsc[XTC][LOOP];
      
      xmx[i][XMB] = -INFTY;
      /* The following section has been hacked to fit a bug in core_algorithms.c
       * The "correct" code is:
       * for (k = hmm->M; k >= 1; k--)
       * xmx[i][XMB] = ILogsum(xmx[i][XMB], mmx[i+1][k] + hmm->bsc[k];
       *
       * The following code gives the same results as core_algorithms.c:
       */
      xmx[i][XMB] = ILogsum(xmx[i][XMB], mmx[i+1][hmm->M] + hmm->bsc[hmm->M-1]);
      for (k = hmm->M-1; k >= 1; k--)
	xmx[i][XMB] = ILogsum(xmx[i][XMB], mmx[i+1][k] + hmm->bsc[k]);
      
      xmx[i][XMJ] = ILogsum(xmx[i][XMB] + hmm->xsc[XTJ][MOVE],
			    xmx[i+1][XMJ] + hmm->xsc[XTJ][LOOP]);

      xmx[i][XME] = ILogsum(xmx[i][XMC] + hmm->xsc[XTE][MOVE],
			    xmx[i][XMJ] + hmm->xsc[XTE][LOOP]);

      xmx[i][XMN] = ILogsum(xmx[i][XMB] + hmm->xsc[XTN][MOVE],
			    xmx[i+1][XMN] + hmm->xsc[XTN][LOOP]);

      /* Now the main states. Note the boundary conditions at M.
       */

      if (i>0) {
	mmx[i][hmm->M] = xmx[i][XME] + hmm->esc[hmm->M] + hmm->msc[dsq[i]][hmm->M];
	dmx[i][hmm->M] = -INFTY;
	for (k = hmm->M-1; k >= 1; k--)
	  {
	    mmx[i][k]  = ILogsum(ILogsum(xmx[i][XME] + hmm->esc[k],
					 mmx[i+1][k+1] + hmm->tsc[TMM][k]),
				 ILogsum(imx[i+1][k] + hmm->tsc[TMI][k],
					 dmx[i][k+1] + hmm->tsc[TMD][k]));
	    mmx[i][k] += hmm->msc[dsq[i]][k];
	    
	    imx[i][k]  = ILogsum(imx[i+1][k] + hmm->tsc[TII][k],
				 mmx[i+1][k+1] + hmm->tsc[TIM][k]);
	    imx[i][k] += hmm->isc[dsq[i]][k];
	    
	    dmx[i][k]  = ILogsum(dmx[i][k+1] + hmm->tsc[TDD][k],
				 mmx[i+1][k+1] + hmm->tsc[TDM][k]);
	    
	  }
      }
      
    }

  sc = xmx[0][XMN];

  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeDPMatrix(mx);

  return Scorify(sc);		/* the total Backward score. */
}

/* Function: Forward()
 * 
 * Note:     This was originally defined as P7Forward() in core_algorithms.c, 
 *           but it is implementation dependent, since it accesses the 
 *           customized dpmatrix structure.  So I moved it here and renamed it 
 *           with a more generic name that fits into the new architecture.
 *             - CRS 15 July 2005
 *
 * Purpose:  The Forward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float
Forward(unsigned char *dsq, int L, struct plan7_s *hmm, cust_dpmatrix_s **ret_mx)
{
  cust_dpmatrix_s *mx;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;  

  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 
  mx = CreateDPMatrix(L, hmm->M, 0, 0);
  xmx = mx->xmx;
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;

  /* Initialization of the zero row.
   * Note that xmx[i][stN] = 0 by definition for all i,
   *    and xmx[i][stT] = xmx[i][stC], so neither stN nor stT need
   *    to be calculated in DP matrices.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   */
  for (i = 1; i <= L; i++)
    {
      mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;
      for (k = 1; k < hmm->M; k++)
	{
	  mmx[i][k]  = ILogsum(ILogsum(mmx[i-1][k-1] + hmm->tsc[TMM][k-1],
				       imx[i-1][k-1] + hmm->tsc[TIM][k-1]),
			       ILogsum(xmx[i-1][XMB] + hmm->bsc[k],
				       dmx[i-1][k-1] + hmm->tsc[TDM][k-1]));
	  mmx[i][k] += hmm->msc[dsq[i]][k];
	  if (mmx[i][k] < -INFTY) mmx[i][k] = -INFTY;

	  dmx[i][k]  = ILogsum(mmx[i][k-1] + hmm->tsc[TMD][k-1],
			       dmx[i][k-1] + hmm->tsc[TDD][k-1]);
	  if (dmx[i][k] < -INFTY) dmx[i][k] = -INFTY;

	  imx[i][k]  = ILogsum(mmx[i-1][k] + hmm->tsc[TMI][k],
			       imx[i-1][k] + hmm->tsc[TII][k]);
	  imx[i][k] += hmm->isc[dsq[i]][k];
	  if (imx[i][k] < -INFTY) imx[i][k] = -INFTY;
	}
      mmx[i][hmm->M] = ILogsum(ILogsum(mmx[i-1][hmm->M-1] + hmm->tsc[TMM][hmm->M-1],
				       imx[i-1][hmm->M-1] + hmm->tsc[TIM][hmm->M-1]),
			       ILogsum(xmx[i-1][XMB] + hmm->bsc[hmm->M],
				       dmx[i-1][hmm->M-1] + hmm->tsc[TDM][hmm->M-1]));
      mmx[i][hmm->M] += hmm->msc[dsq[i]][hmm->M];
      if (mmx[i][hmm->M] < -INFTY) mmx[i][hmm->M] = -INFTY;

      /* Now the special states.
       * remember, C and J emissions are zero score by definition
       */
      xmx[i][XMN] = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP];

      xmx[i][XME] = -INFTY;
      for (k = 1; k <= hmm->M; k++)
	xmx[i][XME] = ILogsum(xmx[i][XME], mmx[i][k] + hmm->esc[k]);

      xmx[i][XMJ] = ILogsum(xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP],
			    xmx[i][XME]   + hmm->xsc[XTE][LOOP]);

      xmx[i][XMB] = ILogsum(xmx[i][XMN] + hmm->xsc[XTN][MOVE],
			    xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]);

      xmx[i][XMC] = ILogsum(xmx[i-1][XMC] + hmm->xsc[XTC][LOOP],
			    xmx[i][XME] + hmm->xsc[XTE][MOVE]);
    }
			    
  sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeDPMatrix(mx);

  printf("Forward: %d\n", sc);

  return Scorify(sc);		/* the total Forward score. */
}

