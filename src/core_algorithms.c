/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* core_algorithms.c
 * SRE, Mon Nov 11 15:58:52 1996
 * RCS $Id$
 * 
 * Simple and robust "research" implementations of Forward, Backward,
 * and Viterbi for Plan7.
 */

#include "structs.h"
#include "config.h"
#include "funcs.h"
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static float get_wee_midpt(struct plan7_s *hmm, char *dsq, int L, 
			   int k1, enum p7stype t1, int s1,
			   int k3, enum p7stype t3, int s3,
			   int *ret_k2, enum p7stype *ret_t2, int *ret_s2);


/* Function: AllocPlan7Matrix()
 * 
 * Purpose:  Allocate a dynamic programming matrix for standard Forward,
 *           Backward, or Viterbi, with scores kept as scaled log-odds
 *           integers. Keeps 2D arrays compact in RAM in an attempt 
 *           to maximize cache hits. Sets up individual ptrs to the
 *           four matrix components as a convenience.
 *           
 * Args:     rows  - number of rows to allocate; typically L+1
 *           M     - size of model
 *           xmx, mmx, imx, dmx 
 *                 - RETURN: ptrs to four mx components as a convenience
 *                 
 * Return:   mx
 *           mx is allocated here. Caller frees with FreeDPMatrix(mx).
 */

struct dpmatrix_s *
AllocPlan7Matrix(int rows, int M, int ***xmx, int ***mmx, int ***imx, int ***dmx)
{
  struct dpmatrix_s *mx;
  int i;

  mx         = (struct dpmatrix_s *) MallocOrDie (sizeof(struct dpmatrix_s));
  mx->xmx    = (int **) MallocOrDie (sizeof(int *) * rows);
  mx->mmx    = (int **) MallocOrDie (sizeof(int *) * rows);
  mx->imx    = (int **) MallocOrDie (sizeof(int *) * rows);
  mx->dmx    = (int **) MallocOrDie (sizeof(int *) * rows);
  mx->xmx[0] = (int *)  MallocOrDie (sizeof(int) * (rows*5));
  mx->mmx[0] = (int *)  MallocOrDie (sizeof(int) * (rows*(M+2)));
  mx->imx[0] = (int *)  MallocOrDie (sizeof(int) * (rows*(M+2)));
  mx->dmx[0] = (int *)  MallocOrDie (sizeof(int) * (rows*(M+2)));
  for (i = 1; i < rows; i++)
    {
      mx->xmx[i] = mx->xmx[0] + (i*5); 
      mx->mmx[i] = mx->mmx[0] + (i*(M+2));
      mx->imx[i] = mx->imx[0] + (i*(M+2));
      mx->dmx[i] = mx->dmx[0] + (i*(M+2));
    }

  if (xmx != NULL) *xmx = mx->xmx;
  if (mmx != NULL) *mmx = mx->mmx;
  if (imx != NULL) *imx = mx->imx;
  if (dmx != NULL) *dmx = mx->dmx;
  return mx;
}

/* Function: FreePlan7Matrix()
 * 
 * Purpose:  Free a dynamic programming matrix allocated by AllocPlan7Matrix().
 * 
 * Return:   (void)
 */
void
FreePlan7Matrix(struct dpmatrix_s *mx)
{
  free (mx->xmx[0]);
  free (mx->mmx[0]);
  free (mx->imx[0]);
  free (mx->dmx[0]);
  free (mx->xmx);
  free (mx->mmx);
  free (mx->imx);
  free (mx->dmx);
  free (mx);
}

/* Function: P7ViterbiSize()
 * Date:     SRE, Fri Mar  6 15:13:20 1998 [St. Louis]
 *
 * Purpose:  Returns the ballpark predicted memory requirement for a 
 *           P7Viterbi() alignment, in MB.
 *
 * Args:     L  - length of sequence
 *           M  - length of HMM       
 *
 * Returns:  # of MB
 */
int
P7ViterbiSize(int L, int M)
{
  return ((sizeof(struct dpmatrix_s)       + /* matrix structure     */
	   3 * (L+1) * (M+2) * sizeof(int) + /* main matrix is O(NM) */ 
	   4 * (L+1) * sizeof(int *)       + /* ptrs into rows of matrix */
	   5 * (L+1) * sizeof(int))          /* 5 special states     */
	  / 1000000);
}

/* Function: P7SmallViterbiSize()
 * Date:     SRE, Fri Mar  6 15:20:04 1998 [St. Louis]
 *
 * Purpose:  Returns the ballpark predicted memory requirement for
 *           a P7SmallViterbi() alignment, in MB. 
 *           
 *           P7SmallViterbi() is a wrapper, calling both P7ParsingViterbi()
 *           and P7WeeViterbi(). P7ParsingViterbi() typically dominates
 *           the memory requirement, so the value returned
 *           is the P7ParsingViterbi() number.
 *
 * Args:     L - length of sequence
 *           M - length of HMM   
 *               
 * Returns:  # of MB
 */
int
P7SmallViterbiSize(int L, int M)
{
  return ((2 * sizeof(struct dpmatrix_s) +
	   12 * (M+2) * sizeof(int)      + /* 2 matrices w/ 2 rows */ 
	   16 * sizeof(int *)            + /* ptrs into rows of matrix */
	   20 * sizeof(int)              + /* 5 special states     */
	   2  * (L+1) * sizeof(int))       /* traceback indices    */      
	  / 1000000);
}


/* Function: P7WeeViterbiSize()
 * Date:     SRE, Fri Mar  6 15:40:42 1998 [St. Louis]
 *
 * Purpose:  Returns the ballpark predicted memory requirement for
 *           a P7WeeViterbi() alignment, in MB. 
 *
 * Args:     L - length of sequence
 *           M - length of HMM   
 *               
 * Returns:  # of MB 
 */
int
P7WeeViterbiSize(int L, int M)
{
  return ((2 * sizeof(struct dpmatrix_s) +
	   12 * (M+2) * sizeof(int)      + /* 2 matrices w/ 2 rows */ 
	   16 * sizeof(int *)            + /* ptrs into rows of matrix */
	   20 * sizeof(int)              + /* 5 special states      */
	   2  * (L+2) * sizeof(int) +      /* stacks for starts/ends (overkill) */      
	   (L+2) * sizeof(int) +           /* k assignments to seq positions */
	   (L+2) * sizeof(enum p7stype))   /* state assignments to seq positions */
	  / 1000000);
}


/* Function: P7Forward()
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
P7Forward(char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s **ret_mx)
{
  struct dpmatrix_s *mx;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;

  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 
  mx = AllocPlan7Matrix(L+1, hmm->M, &xmx, &mmx, &imx, &dmx);

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
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 1; i <= L; i++)
    {
      mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;
      for (k = 1; k < hmm->M; k++)
	{
	  mmx[i][k]  = ILogsum(ILogsum(mmx[i-1][k-1] + hmm->tsc[k-1][TMM],
				     imx[i-1][k-1] + hmm->tsc[k-1][TIM]),
			      ILogsum(xmx[i-1][XMB] + hmm->bsc[k],
				     dmx[i-1][k-1] + hmm->tsc[k-1][TDM]));
	  mmx[i][k] += hmm->msc[(int) dsq[i]][k];

	  dmx[i][k]  = ILogsum(mmx[i][k-1] + hmm->tsc[k-1][TMD],
			      dmx[i][k-1] + hmm->tsc[k-1][TDD]);
	  imx[i][k]  = ILogsum(mmx[i-1][k] + hmm->tsc[k][TMI],
			      imx[i-1][k] + hmm->tsc[k][TII]);
	  imx[i][k] += hmm->isc[(int) dsq[i]][k];
	}
      mmx[i][hmm->M] = ILogsum(ILogsum(mmx[i-1][hmm->M-1] + hmm->tsc[hmm->M-1][TMM],
				   imx[i-1][hmm->M-1] + hmm->tsc[hmm->M-1][TIM]),
			       ILogsum(xmx[i-1][XMB] + hmm->bsc[hmm->M-1],
				   dmx[i-1][hmm->M-1] + hmm->tsc[hmm->M-1][TDM]));
      mmx[i][hmm->M] += hmm->msc[(int) dsq[i]][hmm->M];

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
  else                FreePlan7Matrix(mx);

  return Scorify(sc);		/* the total Forward score. */
}

      
/* Function: P7Viterbi()
 * 
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           Identical to Forward() except that max's
 *           replace sum's. 
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
  int   i,k;
  int   sc;

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
    mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;

    for (k = 1; k <= hmm->M; k++) {
				/* match state */
      mmx[i][k]  = -INFTY;
      if ((sc = mmx[i-1][k-1] + hmm->tsc[k-1][TMM]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = imx[i-1][k-1] + hmm->tsc[k-1][TIM]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = xmx[i-1][XMB] + hmm->bsc[k]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = dmx[i-1][k-1] + hmm->tsc[k-1][TDM]) > mmx[i][k])
	mmx[i][k] = sc;
      if (hmm->msc[(int) dsq[i]][k] != -INFTY) mmx[i][k] += hmm->msc[(int) dsq[i]][k];
      else                                     mmx[i][k] = -INFTY;

				/* delete state */
      dmx[i][k] = -INFTY;
      if ((sc = mmx[i][k-1] + hmm->tsc[k-1][TMD]) > dmx[i][k])
	dmx[i][k] = sc;
      if ((sc = dmx[i][k-1] + hmm->tsc[k-1][TDD]) > dmx[i][k])
	dmx[i][k] = sc;

				/* insert state */
      if (k < hmm->M) {
	imx[i][k] = -INFTY;
	if ((sc = mmx[i-1][k] + hmm->tsc[k][TMI]) > imx[i][k])
	  imx[i][k] = sc;
	if ((sc = imx[i-1][k] + hmm->tsc[k][TII]) > imx[i][k])
	  imx[i][k] = sc;
	if (hmm->isc[(int)dsq[i]][k] != -INFTY) imx[i][k] += hmm->isc[(int) dsq[i]][k];
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
    P7ViterbiTrace(hmm, dsq, L, mx, &tr);
    *ret_tr = tr;
  }

  FreePlan7Matrix(mx);
  return Scorify(sc);		/* the total Viterbi score. */
}


/* Function: P7ViterbiTrace()
 * Date:     SRE, Sat Aug 23 10:30:11 1997 (St. Louis Lambert Field) 
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
P7ViterbiTrace(struct plan7_s *hmm, char *dsq, int N,
	       struct dpmatrix_s *mx, struct p7trace_s **ret_tr)
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
      sc = mmx[i+1][k+1] - hmm->msc[(int) dsq[i+1]][k+1];
      if (sc == xmx[i][XMB] + hmm->bsc[k+1])
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
      else if (sc == mmx[i][k] + hmm->tsc[k][TMM])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == imx[i][k] + hmm->tsc[k][TIM])
	{
	  tr->statetype[tpos] = STI;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == dmx[i][k] + hmm->tsc[k][TDM])
	{
	  tr->statetype[tpos] = STD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	}
      else Die("traceback failed");
      break;

    case STD:			/* D connects from M,D */
      if (dmx[i][k+1] == mmx[i][k] + hmm->tsc[k][TMD])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (dmx[i][k+1] == dmx[i][k] + hmm->tsc[k][TDD]) 
	{
	  tr->statetype[tpos] = STD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	}
      else Die("traceback failed");
      break;

    case STI:			/* I connects from M,I */
      sc = imx[i+1][k] - hmm->isc[(int) dsq[i+1]][k];
      if (sc == mmx[i][k] + hmm->tsc[k][TMI])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == imx[i][k] + hmm->tsc[k][TII])
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
      if (xmx[i][XMB] == xmx[i][XMN] + hmm->xsc[XTN][MOVE])
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
      if (xmx[i][XMC] == xmx[i-1][XMC] + hmm->xsc[XTC][LOOP])
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
      if (xmx[i][XMJ] == xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP])
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


/* Function: P7SmallViterbi()
 * Date:     SRE, Fri Mar  6 15:29:41 1998 [St. Louis]
 *
 * Purpose:  Wrapper function, for linear memory alignment
 *           with same arguments as P7Viterbi(). 
 *           
 *           Calls P7ParsingViterbi to break the sequence
 *           into fragments. Then, based on size of fragments,
 *           calls either P7Viterbi() or P7WeeViterbi() to 
 *           get traces for them. Finally, assembles all these
 *           traces together to produce an overall optimal
 *           trace for the sequence.
 *           
 *           If the trace isn't needed for some reason,
 *           all we do is call P7ParsingViterbi.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *
 * Returns:  Score of optimal alignment in bits.
 */
float
P7SmallViterbi(char *dsq, int L, struct plan7_s *hmm, struct p7trace_s **ret_tr)
{
  struct p7trace_s *ctr;        /* collapsed trace of optimal parse */
  struct p7trace_s *tr;         /* full trace of optimal alignment */
  struct p7trace_s **tarr;      /* trace array */   
  int   ndom;			/* number of subsequences */
  int   i;			/* counter over domains   */
  int   pos;			/* position in sequence */
  int   tpos;			/* position in trace */
  int   tlen;			/* length of full trace   */
  int   sqlen;			/* length of a subsequence */
  int   totlen;                 /* length of L matched by model (as opposed to N/C/J) */
  float sc;			/* score of optimal alignment */
  int   t2;			/* position in a subtrace */
  
  /* Step 1. Call P7ParsingViterbi to calculate an optimal parse
   *         of the sequence into single-hit subsequences; this parse
   *         is returned in a "collapsed" trace
   */
  sc = P7ParsingViterbi(dsq, L, hmm, &ctr);

  /* If we don't want full trace, we're done */
  if (ret_tr == NULL)
    {
      P7FreeTrace(ctr);
      return sc;
    }
  
  /* Step 2. Call either P7Viterbi or P7WeeViterbi on each subsequence
   *         to recover a full traceback of each, collecting them in 
   *         an array. 
   */
  ndom = ctr->tlen/2 - 1;
  tarr = MallocOrDie(sizeof(struct p7trace_s *) * ndom);
  tlen = totlen = 0;
  for (i = 0; i < ndom; i++)
    {
      sqlen = ctr->pos[i*2+2] - ctr->pos[i*2+1];   /* length of subseq */

      if (P7ViterbiSize(sqlen, hmm->M) > RAMLIMIT)
	P7WeeViterbi(dsq + ctr->pos[i*2+1], sqlen, hmm, &(tarr[i]));
      else
	P7Viterbi(dsq + ctr->pos[i*2+1], sqlen, hmm, &(tarr[i]));

      tlen  += tarr[i]->tlen - 4; /* not counting S->N,...,C->T */
      totlen += sqlen;
    }

  /* Step 3. Compose the subtraces into one big final trace.
   *         This is wasteful because we're going to TraceDecompose()
   *         it again in both hmmsearch and hmmpfam to look at
   *         individual domains; but we do it anyway so the P7SmallViterbi
   *         interface looks exactly like the P7Viterbi interface. Maybe
   *         long traces shouldn't include all the N/J/C states anyway,
   *         since they're unambiguously implied.
   */

  /* Calculate total trace len and alloc;
   * nonemitting SNCT + nonemitting J's + emitting NJC
   */
  tlen += 4 + (ndom-1) + (L-totlen);
  P7AllocTrace(tlen, &tr);
  tr->tlen = tlen;

  /* Add N-terminal trace framework
   */
  tr->statetype[0] = STS;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tr->statetype[1] = STN;
  tr->nodeidx[1]   = 0;
  tr->pos[1]       = 0;
  tpos = 2;
				/* add implied N's */
  for (pos = 1; pos <= ctr->pos[1]; pos++)
    {
      tr->statetype[tpos] = STN;
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = pos;
      tpos++;
    }

  /* Add each subseq trace in, with its appropriate 
   * sequence offset derived from the collapsed trace
   */
  for (i = 0; i < ndom; i++)
    {				/* skip SN, CT framework at ends */
      for (t2 = 2; t2 < tarr[i]->tlen-2; t2++)
	{
	  tr->statetype[tpos] = tarr[i]->statetype[t2];
	  tr->nodeidx[tpos]   = tarr[i]->nodeidx[t2];
	  if (tarr[i]->pos[t2] > 0) 
	    tr->pos[tpos]       = tarr[i]->pos[t2] + ctr->pos[i*2+1];
	  else
	    tr->pos[tpos]       = 0;
	  tpos++;
	}
				/* add nonemitting J or C */
      tr->statetype[tpos] = (i == ndom-1) ? STC : STJ;
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = 0;
      tpos++;
				/* add implied emitting J's */
      if (i != ndom-1)
	for (pos = ctr->pos[i*2+2]+1; pos <= ctr->pos[(i+1)*2+1]; pos++)
	  {
	    tr->statetype[tpos] = STJ;
	    tr->nodeidx[tpos]   = 0;
	    tr->pos[tpos]       = pos;
	    tpos++; 
	  }
    }

				/* add implied C's */
  for (pos = ctr->pos[ndom*2]+1; pos <= L; pos++)
    {
      tr->statetype[tpos] = STC;
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = pos;
      tpos++;
    }
				/* add terminal T */
  tr->statetype[tpos] = STT;
  tr->nodeidx[tpos]   = 0;
  tr->pos[tpos]       = 0;
  tpos++;
	    
  for (i = 0; i < ndom; i++) P7FreeTrace(tarr[i]);
  free(tarr);
  P7FreeTrace(ctr);

  *ret_tr = tr;
  return sc;
}




/* Function: P7ParsingViterbi()
 * Date:     SRE, Wed Mar  4 14:07:31 1998 [St. Louis]
 *
 * Purpose:  The "hmmfs" linear-memory algorithm for finding
 *           the optimal alignment of a very long sequence to
 *           a looping, multihit (e.g. Plan7) model, parsing it into
 *           a series of nonoverlapping subsequences that match  
 *           the model once. Other algorithms (e.g. P7Viterbi()
 *           or P7WeeViterbi()) are applied subsequently to
 *           these subsequences to recover complete alignments.
 *           
 *           The hmmfs algorithm appears briefly in [Durbin98],
 *           but is otherwise unpublished.
 *           
 *           The traceback structure returned is special: a
 *           "collapsed" trace S->B->E->...->B->E->T, where
 *           stateidx is unused and pos is used to indicate the
 *           position of B and E in the sequence. The matched
 *           subsequence is B_pos+1...E_pos. The number of
 *           matches in the trace is (tlen/2)-1.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model (log odds scores ready)
 *           ret_tr - RETURN: a collapsed traceback.      
 *
 * Returns:  Score of the optimal Viterbi alignment, in bits.
 */
float
P7ParsingViterbi(char *dsq, int L, struct plan7_s *hmm, struct p7trace_s **ret_tr)
{
  struct dpmatrix_s *mx;        /* two rows of score matrix */
  struct dpmatrix_s *tmx;       /* two rows of misused score matrix: traceback ptrs */
  struct p7trace_s  *tr;        /* RETURN: collapsed traceback */
  int  **xmx, **mmx, **dmx, **imx; /* convenience ptrs to score matrix */  
  int  **xtr, **mtr, **dtr, **itr; /* convenience ptrs to traceback pointers */
  int   *btr, *etr;             /* O(L) trace ptrs for B, E state pts in seq */    
  int    sc;			/* integer score of optimal alignment  */
  int    i,k,tpos;		/* index for seq, model, trace position */
  int    cur, prv;		/* indices for rolling dp matrix */
  int    curralloc;		/* size of allocation for tr */


  /* Alloc a DP matrix and traceback pointers, two rows each, O(M).
   * Alloc two O(L) arrays to trace back through the sequence thru B and E.
   */
  mx  = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  tmx = AllocPlan7Matrix(2, hmm->M, &xtr, &mtr, &itr, &dtr);
  btr = MallocOrDie(sizeof(int) * (L+1));
  etr = MallocOrDie(sizeof(int) * (L+1));

  /* Initialization of the zero row.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  btr[0]      = 0;
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  etr[0]      = -1; 
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Recursion. Done as a pull. Rolling index trick. Trace ptr propagation trick.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   *    
   * Notes on traceback pointer propagation.
   *  - In the path B->E, we propagate the i that B was aligned to in the optimal
   *    alignment, via mtr, dtr, and itr. 
   *  - When we reach an E, we record the i of the B it started from in etr.
   *  - In a looping path E->J...->B or terminal path E->C...->T, we propagate
   *    the i that E was aligned to in the optimal alignment via xtr[][XMC]
   *    and xtr[][XMJ].
   *  - When we enter B, we record the i of the best previous E, or 0 if there
   *    isn't one, in btr.
   */
  for (i = 1; i <= L; i++) {
    cur = i % 2;
    prv = !cur;
    
    mmx[cur][0] = imx[cur][0] = dmx[cur][0] = -INFTY;

    for (k = 1; k <= hmm->M; k++) {
				/* match state */
      mmx[cur][k] = -INFTY;
      if ((sc = mmx[prv][k-1] + hmm->tsc[k-1][TMM]) > -INFTY)
	{ mmx[cur][k] = sc; mtr[cur][k] = mtr[prv][k-1]; }
      if ((sc = imx[prv][k-1] + hmm->tsc[k-1][TIM]) > mmx[cur][k])
	{ mmx[cur][k] = sc; mtr[cur][k] = itr[prv][k-1]; }
      if ((sc = xmx[prv][XMB] + hmm->bsc[k]) > mmx[cur][k])
	{ mmx[cur][k] = sc; mtr[cur][k] = i-1; }
      if ((sc = dmx[prv][k-1] + hmm->tsc[k-1][TDM]) > mmx[cur][k])
	{ mmx[cur][k] = sc; mtr[cur][k] = dtr[prv][k-1]; }
      if (hmm->msc[(int) dsq[i]][k] != -INFTY)
	mmx[cur][k] += hmm->msc[(int) dsq[i]][k];
      else
	mmx[cur][k] = -INFTY;

				/* delete state */
      dmx[cur][k] = -INFTY;
      if ((sc = mmx[cur][k-1] + hmm->tsc[k-1][TMD]) > -INFTY)
	{ dmx[cur][k] = sc; dtr[cur][k] = mtr[cur][k-1]; }
      if ((sc = dmx[cur][k-1] + hmm->tsc[k-1][TDD]) > dmx[cur][k])
	{ dmx[cur][k] = sc; dtr[cur][k] = dtr[cur][k-1]; }

				/* insert state */
      if (k < hmm->M) {
	imx[cur][k] = -INFTY;
	if ((sc = mmx[prv][k] + hmm->tsc[k][TMI]) > -INFTY)
	  { imx[cur][k] = sc; itr[cur][k] = mtr[prv][k]; }
	if ((sc = imx[prv][k] + hmm->tsc[k][TII]) > imx[cur][k])
	  { imx[cur][k] = sc; itr[cur][k] = itr[prv][k]; }
	if (hmm->isc[(int) dsq[i]][k] != -INFTY)
	  imx[cur][k] += hmm->isc[(int) dsq[i]][k];
	else
	  imx[cur][k] = -INFTY;
      }
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
				/* N state */
    xmx[cur][XMN] = -INFTY;
    if ((sc = xmx[prv][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[cur][XMN] = sc;
				/* E state */
    xmx[cur][XME] = -INFTY;
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  mmx[cur][k] + hmm->esc[k]) > xmx[cur][XME])
	{ xmx[cur][XME] = sc; etr[i] = mtr[cur][k]; }
				/* J state */
    xmx[cur][XMJ] = -INFTY;
    if ((sc = xmx[prv][XMJ] + hmm->xsc[XTJ][LOOP]) > -INFTY)
      { xmx[cur][XMJ] = sc; xtr[cur][XMJ] = xtr[prv][XMJ]; }
    if ((sc = xmx[cur][XME]   + hmm->xsc[XTE][LOOP]) > xmx[cur][XMJ])
      { xmx[cur][XMJ] = sc; xtr[cur][XMJ] = i; }
				/* B state */
    xmx[cur][XMB] = -INFTY;
    if ((sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
      { xmx[cur][XMB] = sc; btr[i] = 0; }
    if ((sc = xmx[cur][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[cur][XMB])
      { xmx[cur][XMB] = sc; btr[i] = xtr[cur][XMJ]; }
				/* C state */
    xmx[cur][XMC] = -INFTY;
    if ((sc = xmx[prv][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      { xmx[cur][XMC] = sc; xtr[cur][XMC] = xtr[prv][XMC]; }
    if ((sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > xmx[cur][XMC])
      { xmx[cur][XMC] = sc; xtr[cur][XMC] = i; }
  }
				/* T state (not stored) */
  sc = xmx[cur][XMC] + hmm->xsc[XTC][MOVE];

  /*****************************************************************
   * Collapsed traceback stage. 
   * xtr[L%2][XMC] contains the position j of the previous E
   * etr[j]        contains the position i of the previous B
   * btr[i]        contains the position j of the previous E, or 0
   * continue until btr[i] = 0.
   *****************************************************************/

  curralloc = 2;		/* minimum: no hits */
  P7AllocTrace(curralloc, &tr);

  /* Init of collapsed trace. Back to front; we ReverseTrace() later.
   */
  tpos = 0;
  tr->statetype[tpos] = STT;
  tr->pos[tpos]       = 0;
  i                   = xtr[L%2][XMC];
  while (i > 0)
    {
      curralloc += 2;
      P7ReallocTrace(tr, curralloc);

      tpos++;
      tr->statetype[tpos] = STE;
      tr->pos[tpos] = i;
      i = etr[i];

      tpos++;
      tr->statetype[tpos] = STB;
      tr->pos[tpos] = i;
      i = btr[i];
    }

  tpos++;
  tr->statetype[tpos] = STS;
  tr->pos[tpos]       = 0;
  tr->tlen = tpos + 1;
  P7ReverseTrace(tr);
  
  FreePlan7Matrix(mx);
  FreePlan7Matrix(tmx);
  free(btr);
  free(etr);

  *ret_tr = tr;
  return Scorify(sc);
}

/* Function: P7WeeViterbi()
 * Date:     SRE, Wed Mar  4 08:24:04 1998 [St. Louis]
 *
 * Purpose:  Hirschberg/Myers/Miller linear memory alignment.
 *           See [Hirschberg75,MyM-88a] for the idea of the algorithm.
 *           Adapted to HMM implementation. 
 *           
 *           Requires that you /know/ that there's only
 *           one hit to the model in the sequence: either
 *           because you're forcing single-hit, or you've
 *           previously called P7ParsingViterbi to parse
 *           the sequence into single-hit segments. The reason
 *           for this is that a cyclic model (a la Plan7)
 *           defeats the nice divide and conquer trick.
 *           (I think some trickery with propagated trace pointers
 *           could get around this but haven't explored it.)
 *           This is implemented by ignoring transitions
 *           to/from J state.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_tr - RETURN: traceback.
 *
 * Returns:  Score of the optimal Viterbi alignment.
 */
float
P7WeeViterbi(char *dsq, int L, struct plan7_s *hmm, struct p7trace_s **ret_tr)
{
  struct p7trace_s *tr;         /* RETURN: traceback */
  int          *kassign;        /* 0..L+1, alignment of seq positions to model nodes */
  enum p7stype *tassign;        /* 0..L+1, alignment of seq positions to state types */
  int          *endlist;        /* stack of end points on sequence to work on */
  int          *startlist;      /* stack of start points on sequence to work on */
  int          lpos;            /* position in endlist, startlist */
  int          k1, k2, k3;	/* start, mid, end in model      */
  enum p7stype t1, t2, t3;	/* start, mid, end in state type */
  int          s1, s2, s3;	/* start, mid, end in sequence   */
  float        sc;		/* score of segment optimal alignment */
  float        ret_sc;		/* optimal score over complete seq */
  int          tlen;		/* length needed for trace */
  int          i, k, tpos;	/* index in sequence, model, trace */


  /* Initialize.
   */
  kassign   = MallocOrDie (sizeof(int) * (L+2));
  tassign   = MallocOrDie (sizeof(enum p7stype) * (L+2));
  endlist   = MallocOrDie (sizeof(int) * (L+2));
  startlist = MallocOrDie (sizeof(int) * (L+2));

  lpos            = 0; 
  startlist[lpos] = 0;
  endlist[lpos]   = L+1;
  kassign[0]      = 0;
  kassign[L+1]    = hmm->M+1;
  tassign[0]      = STN;		/* artificial; should be STS */
  tassign[L+1]    = STC;		/* artificial; should be STT */

  /* Recursive divide-and-conquer alignment.
   */
  while (lpos >= 0)
    {
				/* Pop a segment off the stack */
      s1 = startlist[lpos];
      k1 = kassign[s1];
      t1 = tassign[s1];
      s3 = endlist[lpos];
      k3 = kassign[s3];
      t3 = tassign[s3];
      lpos--;
				/* find optimal midpoint of segment */
      sc = get_wee_midpt(hmm, dsq, L, k1, t1, s1, k3, t3, s3, &k2, &t2, &s2);
      kassign[s2] = k2;
      tassign[s2] = t2;
                               /* score is valid on first pass */
      if (s1 == 0 && s3 == L+1) ret_sc = sc;

				/* push N-terminal segment on stack */
      if (s2 - s1 > 1 && t2 != STN)
	{
	  lpos++;
	  startlist[lpos] = s1;
	  endlist[lpos]   = s2;
	}
				/* push C-terminal segment on stack */
      if (s3 - s2 > 1 && t2 != STC)
	{
          lpos++;
          startlist[lpos] = s2;
          endlist[lpos]   = s3;
	}

      if (t2 == STN)
	{			/* if we see STN midpoint, we know the whole N-term is STN */
	  for (; s2 >= s1; s2--) {
	    kassign[s2] = 0;
	    tassign[s2] = STN;
	  }
	}
      if (t2 == STC)
	{			/* if we see STC midpoint, we know whole C-term is STC */
	  for (; s2 <= s3; s2++) {
	    kassign[s2] = hmm->M+1;
	    tassign[s2] = STC;
	  }
	}
    }

  /*****************************************************************
   * Construct a traceback structure from kassign/tassign by interpolating
   * necessary states. 
   * Trace allocation is as follows. We clearly need L emitting states.
   * We also need nonemitting states as follows:
   * STS,STN,STB,STE,STC,STT = 6
   * STD: count k2-k1-1 in kassign M->M's
   * Also, count N->M's and M->C's (potential wing unfoldings),
   *****************************************************************/ 

  tlen = L + 6;
  for (i = 1; i < L; i++)
    {
      if (tassign[i] == STM && tassign[i+1] == STM)
	tlen += kassign[i+1] - kassign[i] - 1;
      if (tassign[i] == STN && tassign[i+1] == STM)
	tlen += kassign[i+1] - 1;
      if (tassign[i] == STM && tassign[i+1] == STC)
	tlen += hmm->M - kassign[i];
    }
  P7AllocTrace(tlen, &tr);

  tr->statetype[0] = STS;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tr->statetype[1] = STN;
  tr->nodeidx[1]   = 0;
  tr->pos[1]       = 0;
  tpos = 2;
  
  for (i = 1; i <= L; i++)
    {
      switch(tassign[i]) {
      case STM:
				/* check for first match state */
	if (tr->statetype[tpos-1] == STN) {
	  tr->statetype[tpos] = STB;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	  tpos++;
				/* check for wing unfolding */
	  if (Prob2Score(hmm->begin[kassign[i]], hmm->p1) + INTSCALE <= hmm->bsc[kassign[i]])
	    for (k = 1; k < kassign[i]; k++) {
	      tr->statetype[tpos] = STD;
	      tr->nodeidx[tpos]   = k;
	      tr->pos[tpos]       = 0;
	      tpos++;
	    }
	}
				/* do the match state itself */
	tr->statetype[tpos] = STM;
	tr->nodeidx[tpos]   = kassign[i];
	tr->pos[tpos]       = i;
	tpos++;
				/* do any deletes necessary 'til next match */
	if (i < L && tassign[i+1] == STM && kassign[i+1] - kassign[i] > 1)
	  for (k = kassign[i] + 1; k < kassign[i+1]; k++)
	    {
	      tr->statetype[tpos] = STD;
	      tr->nodeidx[tpos]   = k;
	      tr->pos[tpos]       = 0;
	      tpos++;
	    }
				/* check for last match state */
	if (i == L || tassign[i+1] == STC) {
				/* check for wing unfolding */
	  if (Prob2Score(hmm->end[kassign[i-1]], 1.) + INTSCALE <=  hmm->esc[kassign[i-1]])
	    for (k = kassign[i]+1; k <= hmm->M; k++)
	      {
		tr->statetype[tpos] = STD;
		tr->nodeidx[tpos]   = k;
		tr->pos[tpos]       = 0;
		tpos++;
	      }
				/* add on the end state */
	  tr->statetype[tpos] = STE;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0; 
	  tpos++;
				/* and a nonemitting C state */
	  tr->statetype[tpos] = STC;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0; 
	  tpos++;
	}
	break;
	
      case STI:
	tr->statetype[tpos] = STI;
	tr->nodeidx[tpos]   = kassign[i];
	tr->pos[tpos]       = i;
	tpos++;
	break;

      case STN:
	tr->statetype[tpos] = STN;
	tr->nodeidx[tpos]   = 0;
	tr->pos[tpos]       = i;
	tpos++;
	break;

      case STC:
	tr->statetype[tpos] = STC;
	tr->nodeidx[tpos]   = 0;
	tr->pos[tpos]       = i; 
	tpos++;
	break;

      default: Die("Bogus state %s", Statetype(tassign[i]));
      }
    }
				/* terminate the trace */
  tr->statetype[tpos] = STT;
  tr->nodeidx[tpos]   = 0;
  tr->pos[tpos]       = 0; 
  tr->tlen = tpos+1;

  *ret_tr = tr;

  free(kassign);
  free(tassign);
  free(startlist);
  free(endlist);
  return ret_sc;
}


/* Function: Plan7ESTViterbi()
 * 
 * Purpose:  Frameshift-tolerant alignment of protein model to cDNA EST.
 *           
 * 
 */
float
Plan7ESTViterbi(char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s **ret_mx)
{
  struct dpmatrix_s *mx;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;
  int   codon;
  
  /* Allocate a DP matrix with 0..L rows, 0..M+1 columns.
   */ 
  mx = AllocPlan7Matrix(L+1, hmm->M, &xmx, &mmx, &imx, &dmx);

  /* Initialization of the zero row (DNA sequence of length 0)
   * Note that xmx[i][stN] = 0 by definition for all i,
   *    and xmx[i][stT] = xmx[i][stC], so neither stN nor stT need
   *    to be calculated in DP matrices.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Initialization of the first row (DNA sequence of length 1)
   * only N can make this nucleotide.
   */
  xmx[1][XMN] = xmx[0][XMN] + hmm->xsc[XTN][LOOP];
  xmx[1][XMB] = xmx[1][XMN] + hmm->xsc[XTN][MOVE]; 
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need 2 nt to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need 2 nt to get into model */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 2; i <= L; i++) {
    mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;

				/* crude calculation of lookup value for codon */
    if (i > 2) {
      if (dsq[i-2] < 4 && dsq[i-1] < 4 && dsq[i] < 4)
	codon = dsq[i-2] * 16 + dsq[i-1] * 4 + dsq[i];
      else
	codon = 64;		/* ambiguous codon; punt */
    }

    for (k = 1; k <= hmm->M; k++) {
				/* match state */
      if (i > 2) {
	mmx[i][k]  = mmx[i-3][k-1] + hmm->tsc[k-1][TMM];
	if ((sc = imx[i-3][k-1] + hmm->tsc[k-1][TIM]) > mmx[i][k])
	  mmx[i][k] = sc;
	if ((sc = xmx[i-3][XMB] + hmm->bsc[k]) > mmx[i][k])
	  mmx[i][k] = sc;
	if ((sc = dmx[i-3][k-1] + hmm->tsc[k-1][TDM]) > mmx[i][k])
	  mmx[i][k] = sc;
	mmx[i][k] += hmm->dnam[codon][k];
      }
				/* -1 frameshifts into match state */
      if ((sc = mmx[i-2][k-1] + hmm->tsc[k-1][TMM] + hmm->dna2) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = imx[i-2][k-1] + hmm->tsc[k-1][TIM] + hmm->dna2) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = xmx[i-2][XMB] + hmm->bsc[k] + hmm->dna2) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = dmx[i-2][k-1] + hmm->tsc[k-1][TDM] + hmm->dna2) > mmx[i][k])
	mmx[i][k] = sc;
      
				/* +1 frameshifts into match state */
      if (i > 3) {
	if ((sc = mmx[i-4][k-1] + hmm->tsc[k-1][TMM] + hmm->dna4) > mmx[i][k])
	  mmx[i][k] = sc;
	if ((sc = imx[i-4][k-1] + hmm->tsc[k-1][TIM] + hmm->dna4) > mmx[i][k])
	  mmx[i][k] = sc;
	if ((sc = xmx[i-4][XMB] + hmm->bsc[k] + hmm->dna4) > mmx[i][k])
	  mmx[i][k] = sc;
	if ((sc = dmx[i-4][k-1] + hmm->tsc[k-1][TDM] + hmm->dna4) > mmx[i][k])
	  mmx[i][k] = sc;
      }
      				/* delete state */
      dmx[i][k]  = mmx[i][k-1] + hmm->tsc[k-1][TMD];
      if ((sc = dmx[i][k-1] + hmm->tsc[k-1][TDD]) > dmx[i][k])
	dmx[i][k] = sc;

				/* insert state */
      if (i > 2) {
	imx[i][k] = mmx[i-3][k] + hmm->tsc[k][TMI];
	if ((sc = imx[i-3][k] + hmm->tsc[k][TII]) > imx[i][k])
	  imx[i][k] = sc;
	imx[i][k] += hmm->dnai[codon][k];
      }

				/* -1 frameshifts into insert state */
      if ((sc = mmx[i-2][k] + hmm->tsc[k][TMI] + hmm->dna2) > imx[i][k])
	imx[i][k] = sc;
      if ((sc = imx[i-2][k] + hmm->tsc[k][TII] + hmm->dna2) > imx[i][k])
	imx[i][k] = sc;

				/* +1 frameshifts into insert state */
      if (i > 4) {
	if ((sc = mmx[i-4][k] + hmm->tsc[k][TMI] + hmm->dna4) > imx[i][k])
	  imx[i][k] = sc;
	if ((sc = imx[i-4][k] + hmm->tsc[k][TII] + hmm->dna4) > imx[i][k])
	  imx[i][k] = sc;
      }
    }
    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
				/* N state: +1 nucleotide */
    xmx[i][XMN] = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP];
                                /* E state: collect from M's, and last D  */
    xmx[i][XME] = dmx[i][hmm->M];    /* transition prob from last D = 1.0 */
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  mmx[i][k] + hmm->esc[k]) > xmx[i][XME])
        xmx[i][XME] = sc;
                                /* J state: +1 nucleotide */
    xmx[i][XMJ] = xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP];
    if ((sc = xmx[i][XME]   + hmm->xsc[XTE][LOOP]) > xmx[i][XMJ])
      xmx[i][XMJ] = sc;
                                /* B state: collect from N,J */
    xmx[i][XMB] = xmx[i][XMN] + hmm->xsc[XTN][MOVE];
    if ((sc = xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[i][XMB])
      xmx[i][XMB] = sc;
				/* C state: +1 nucleotide */
    xmx[i][XMC] = xmx[i-1][XMC] + hmm->xsc[XTC][LOOP];
    if ((sc = xmx[i][XME] + hmm->xsc[XTE][MOVE]) > xmx[i][XMC])
      xmx[i][XMC] = sc;
  }

  sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

  if (ret_mx != NULL) *ret_mx = mx;
  else                FreePlan7Matrix(mx);

  return Scorify(sc);            /* the total Viterbi score. */
}



/* Function: get_wee_midpt()
 * Date:     SRE, Wed Mar  4 08:27:11 1998 [St. Louis]
 *
 * Purpose:  The heart of the divide and conquer algorithm
 *           for P7WeeViterbi(). This function is called
 *           recursively to find successive optimal midpoints
 *           in the alignment matrix. See P7WeeViterbi() for
 *           further comments on the assumptions of this algorithm.
 *           
 * Args:     hmm   - the model, set up for integer scores
 *           dsq   - the sequence, digitized
 *           L     - length of the sequence
 *           k1    - model node to start with, 0,1..M,M+1
 *           t1    - state type to start with, STM | STI | STN | STC | STJ
 *           s1    - sequence position to start with, 0,1..L,L+1
 *           k3    - model node to end with
 *           t3    - state type to end with
 *           s3    - sequence position to end with
 *          ret_k2 - RETURN: optimal midpoint, node position in model
 *          ret_t2 - RETURN: optimal midpoint, state type 
 *          ret_s2 - RETURN: optimal midpoint, sequence position
 *
 * Returns: score of optimal alignment, in bits. 
 */
static float
get_wee_midpt(struct plan7_s *hmm, char *dsq, int L, 
	      int k1, enum p7stype t1, int s1,
	      int k3, enum p7stype t3, int s3,
	      int *ret_k2, enum p7stype *ret_t2, int *ret_s2)
{
  struct dpmatrix_s *fwd;
  struct dpmatrix_s *bck;
  int        **xmx;             /* convenience ptr into special states */
  int        **mmx;             /* convenience ptr into match states   */
  int        **imx;             /* convenience ptr into insert states  */
  int        **dmx;             /* convenience ptr into delete states  */
  int          k2;
  enum p7stype t2;
  int          s2;
  int          cur, prv, nxt;	/* current, previous, next row index (0 or 1)*/
  int          i,k;		/* indices for seq, model */
  int          sc;		/* integer score */
  int          max;		/* maximum integer score */

 
  /* Choose our midpoint.
   */
  s2 = s1 + (s3-s1) / 2;

  /* Allocate our forward two rows.
   * Initialize row zero.
   */
  fwd = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  cur = s1%2;
  xmx[cur][XMN] = xmx[cur][XMB] = -INFTY;
  xmx[cur][XME] = xmx[cur][XMC] = -INFTY;  
  for (k = k1; k <= k3; k++)
    mmx[cur][k] = imx[cur][k] = dmx[cur][k] = -INFTY;      

  /* Where to put our zero for our start point...
   * (only possible to start on an emitting state; J disallowed)
   */
  switch (t1) {
  case STM: mmx[cur][k1]  = 0; break;
  case STI: imx[cur][k1]  = 0; break;
  case STN: xmx[cur][XMN] = 0; break;
  case STC: xmx[cur][XMC] = 0; break;
  default:  Die("you can't init get_wee_midpt with a %s\n", Statetype(t1));
  }

  /* Still initializing.
   * Deal with pulling horizontal matrix moves in row zero.
   * Some of this we know we only have to do if we started on an STM.
   */ 
  if (t1 == STM)
    {
      for (k = k1+1; k <= k3 && k <= hmm->M; k++)
	{				/* transits into STD */
	  dmx[cur][k] = -INFTY;
	  if ((sc = mmx[cur][k-1] + hmm->tsc[k-1][TMD]) > -INFTY)
	    dmx[cur][k] = sc;
	  if ((sc = dmx[cur][k-1] + hmm->tsc[k-1][TDD]) > dmx[cur][k])
	    dmx[cur][k] = sc;
	}
				/* transit into STE */
      xmx[cur][XME] = -INFTY;
      if ((sc = mmx[cur][k1] + hmm->esc[k1]) > -INFTY)
	xmx[cur][XME] = sc;
    }
				/* transit into STB from STN */
  xmx[cur][XMB] = -INFTY;
  if ((sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
    xmx[cur][XMB] = sc;
				/* transit into STC from STE */
  xmx[cur][XMC] = -INFTY;
  if ((sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > -INFTY)
    xmx[cur][XMC] = sc;
  
  /* Done initializing.
   * Start recursive DP; sweep forward to chosen s2 midpoint. Done as a pull.
   */
  for (i = s1+1; i <= s2; i++) {
    cur = i % 2;
    prv = !cur;

    mmx[cur][k1] = imx[cur][k1] = dmx[cur][k1] = -INFTY;

    /* Insert state in column k1
     */
    if (k1 > 0 && k1 < hmm->M) {
      imx[cur][k1] = -INFTY;
      if ((sc = mmx[prv][k1] + hmm->tsc[k1][TMI]) > -INFTY)
	imx[cur][k1] = sc;
      if ((sc = imx[prv][k1] + hmm->tsc[k1][TII]) > imx[cur][k1])
	imx[cur][k1] = sc;
      if (hmm->isc[(int) dsq[i]][k1] != -INFTY)
	imx[cur][k1] += hmm->isc[(int) dsq[i]][k1];
      else
	imx[cur][k1] = -INFTY;
    }

    /* Main chunk of recursion across model positions
     */
    for (k = k1+1; k <= k3 && k <= hmm->M; k++) {
				/* match state */
      mmx[cur][k]  = -INFTY;
      if ((sc = mmx[prv][k-1] + hmm->tsc[k-1][TMM]) > -INFTY)
	mmx[cur][k] = sc;
      if ((sc = imx[prv][k-1] + hmm->tsc[k-1][TIM]) > mmx[cur][k])
	mmx[cur][k] = sc;
      if ((sc = xmx[prv][XMB] + hmm->bsc[k]) > mmx[cur][k])
	mmx[cur][k] = sc;
      if ((sc = dmx[prv][k-1] + hmm->tsc[k-1][TDM]) > mmx[cur][k])
	mmx[cur][k] = sc;
      if (hmm->msc[(int) dsq[i]][k] != -INFTY)
	mmx[cur][k] += hmm->msc[(int) dsq[i]][k];
      else
	mmx[cur][k] = -INFTY;

				/* delete state */
      dmx[cur][k] = -INFTY;
      if ((sc = mmx[cur][k-1] + hmm->tsc[k-1][TMD]) > -INFTY)
	dmx[cur][k] = sc;
      if ((sc = dmx[cur][k-1] + hmm->tsc[k-1][TDD]) > dmx[cur][k])
	dmx[cur][k] = sc;

				/* insert state */
      if (k < hmm->M) {
	imx[cur][k] = -INFTY;
	if ((sc = mmx[prv][k] + hmm->tsc[k][TMI]) > -INFTY)
	  imx[cur][k] = sc;
	if ((sc = imx[prv][k] + hmm->tsc[k][TII]) > imx[cur][k])
	  imx[cur][k] = sc;
	if (hmm->isc[(int) dsq[i]][k] != -INFTY)
	  imx[cur][k] += hmm->isc[(int) dsq[i]][k];
	else
	  imx[cur][k] = -INFTY;
      }
    }
				/* N state */
    xmx[cur][XMN] = -INFTY;
    if ((sc = xmx[prv][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[cur][XMN] = sc;
				/* E state */
    xmx[cur][XME] = -INFTY;
    for (k = k1; k <= k3 && k <= hmm->M; k++)
      if ((sc =  mmx[cur][k] + hmm->esc[k]) > xmx[cur][XME])
	xmx[cur][XME] = sc;
				/* B state */
    xmx[cur][XMB] = -INFTY;
    if ((sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
      xmx[cur][XMB] = sc;
				/* C state */
    xmx[cur][XMC] = -INFTY;
    if ((sc = xmx[prv][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[cur][XMC] = sc;
    if ((sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > xmx[cur][XMC])
      xmx[cur][XMC] = sc;
  }

  /* Row s2%2 in fwd matrix now contains valid scores from s1 to s2,
   * with J transitions disallowed (no cycles through model). 
   */

  /*****************************************************************
   * Backwards pass.
   *****************************************************************/ 

  /* Allocate our backwards two rows. Init last row.
   */
  bck = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  nxt = s3%2;
  xmx[nxt][XMN] = xmx[nxt][XMB] = -INFTY;
  xmx[nxt][XME] = xmx[nxt][XMC] = -INFTY;  
  for (k = k1; k <= k3 && k <= hmm->M; k++)
    mmx[nxt][k] = imx[nxt][k] = dmx[nxt][k] = -INFTY;      

  /* Where to put the zero for our end point...
   */
  switch (t3) {
  case STM: mmx[nxt][k3]  = 0; break;
  case STI: imx[nxt][k3]  = 0; break;
  case STN: xmx[nxt][XMN] = 0; break;
  case STC: xmx[nxt][XMC] = (s3 == L+1) ? hmm->xsc[XTC][MOVE] : 0; break;
  default:  Die("you can't init get_wee_midpt with a %s\n", Statetype(t3));
  }

  /* Start recursive DP; sweep backwards to chosen s2 midpoint.
   * Done as a pull. M, I scores at current row do /not/ include
   * emission scores.
   */
  for (i = s3-1; i >= s2; i--) {
    cur = i%2;
    nxt = !cur;

    mmx[cur][k3] = dmx[cur][k3] = imx[cur][k3] = -INFTY;

				/* C pulls from C (T is special cased) */
    xmx[cur][XMC] = -INFTY;
    if ((sc = xmx[nxt][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[cur][XMC] = sc;
				/* B pulls from M's */
    xmx[cur][XMB] = -INFTY;
    if (i < L)
      for (k = k1; k <= k3 && k <= hmm->M; k++)
	if ((sc = mmx[nxt][k] + hmm->msc[(int)dsq[i+1]][k] + hmm->bsc[k]) > xmx[cur][XMB])
	  xmx[cur][XMB] = sc;
				/* E pulls from C (J disallowed) */
    xmx[cur][XME] = -INFTY;
    if ((sc = xmx[cur][XMC] + hmm->xsc[XTE][MOVE]) > -INFTY)
      xmx[cur][XME] = sc;
				/* N pulls from B, N */
    xmx[cur][XMN] = -INFTY;
    if ((sc = xmx[cur][XMB] + hmm->xsc[XTN][MOVE]) > -INFTY)
      xmx[cur][XMN] = sc;
    if ((sc = xmx[nxt][XMN] + hmm->xsc[XTN][LOOP]) > xmx[cur][XMN])
      xmx[cur][XMN] = sc;
				/* Special case k=M */
    if (k3-1 == hmm->M) {
      mmx[cur][hmm->M] = xmx[cur][XME];     /* p=1 transit to E, by definition */
      dmx[cur][hmm->M] = -INFTY;            /* by definition; no D state there */
      imx[cur][hmm->M] = -INFTY;	          /* by definition; no I state there */
    }

    /* Deal with right boundary
     */
    if (i < L && k3 < hmm->M) {
      imx[cur][k3] = -INFTY;
      if ((sc = imx[nxt][k3] + hmm->isc[(int)dsq[i+1]][k3] + hmm->tsc[k3][TII]) > -INFTY)
	imx[cur][k3] = sc;
      mmx[cur][k3] = -INFTY;
      if ((sc = imx[nxt][k3] + hmm->isc[(int)dsq[i+1]][k3] + hmm->tsc[k3][TMI]) > -INFTY)
	mmx[cur][k3] = sc;
    }

    /* Main recursion across model
     */
    for (k = (k3 >= hmm->M) ? hmm->M-1 : k3-1; k >= k1 && k >= 1; k--)  {
				/* pull into match state */
      mmx[cur][k] = -INFTY;
      if ((sc = xmx[cur][XME] + hmm->esc[k]) > -INFTY)
	mmx[cur][k] = sc;
      if (i < L && 
	  (sc = mmx[nxt][k+1] + hmm->msc[(int)dsq[i+1]][k+1] + hmm->tsc[k][TMM]) > mmx[cur][k])
	mmx[cur][k] = sc;
      if (i < L && 
	  (sc = imx[nxt][k] + hmm->isc[(int)dsq[i+1]][k] + hmm->tsc[k][TMI]) > mmx[cur][k])
	mmx[cur][k] = sc;
      if ((sc = dmx[cur][k+1] + hmm->tsc[k][TMD]) > mmx[cur][k])
	mmx[cur][k] = sc;
				/* pull into delete state */
      dmx[cur][k] = -INFTY;
      if ((sc = mmx[nxt][k+1] + hmm->tsc[k][TDM]) > -INFTY)
	dmx[cur][k] = sc;
      if (i < L && hmm->msc[(int)dsq[i+1]][k+1] != -INFTY) 
	dmx[cur][k] += hmm->msc[(int)dsq[i+1]][k+1]; 
      else
	dmx[cur][k] = -INFTY;
      if ((sc = dmx[cur][k+1] + hmm->tsc[k][TDD]) > dmx[cur][k])
	dmx[cur][k] = sc;
				/* pull into insert state */
      if (i < L)
	{
	  imx[cur][k] = -INFTY;
	  if ((sc = mmx[nxt][k+1] + hmm->msc[(int)dsq[i+1]][k+1] + hmm->tsc[k][TIM]) > -INFTY)
	    imx[cur][k] = sc;
	  if ((sc = imx[nxt][k] + hmm->isc[(int)dsq[i+1]][k] + hmm->tsc[k][TII]) > imx[cur][k])
	    imx[cur][k] = sc;
	}
    }
  }
   
  /*****************************************************************
   * DP complete; we have both forward and backward passes. Now we
   * look across the s2 row and find the optimal emitting state.
   *****************************************************************/  

  cur = s2%2;
  max = -INFTY;
  for (k = k1; k <= k3; k++)
    {
      if ((sc = fwd->mmx[cur][k] + bck->mmx[cur][k]) > max)
	{ k2 = k; t2 = STM; max = sc; }
      if ((sc = fwd->imx[cur][k] + bck->imx[cur][k]) > max)
	{ k2 = k; t2 = STI; max = sc; }
    }
  if ((sc = fwd->xmx[cur][XMN] + bck->xmx[cur][XMN]) > max)
    { k2 = 0;        t2 = STN; max = sc; }
  if ((sc = fwd->xmx[cur][XMC] + bck->xmx[cur][XMC]) > max)
    { k2 = hmm->M+1; t2 = STC; max = sc; }

  /*****************************************************************
   * Garbage collection, return.
   *****************************************************************/
  
  FreePlan7Matrix(fwd);
  FreePlan7Matrix(bck);
  *ret_k2 = k2;
  *ret_t2 = t2;
  *ret_s2 = s2;
  return Scorify(max);
}
