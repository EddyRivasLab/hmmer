/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* core_algorithms.c
 * SRE, Mon Nov 11 15:58:52 1996
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
  mx->mmx[0] = (int *)  MallocOrDie (sizeof(int) * (rows*(M+1)));
  mx->imx[0] = (int *)  MallocOrDie (sizeof(int) * (rows*(M+1)));
  mx->dmx[0] = (int *)  MallocOrDie (sizeof(int) * (rows*(M+1)));
  for (i = 1; i < rows; i++)
    {
      mx->xmx[i] = mx->xmx[0] + (i*5); 
      mx->mmx[i] = mx->mmx[0] + (i*(M+1));
      mx->imx[i] = mx->imx[0] + (i*(M+1));
      mx->dmx[i] = mx->dmx[0] + (i*(M+1));
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

/* Function: Plan7Forward()
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
Plan7Forward(char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s **ret_mx)
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

      
/* Function: Plan7Viterbi()
 * 
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           Identical to Forward() except that max's
 *           replace sum's. 
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float
Plan7Viterbi(char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s **ret_mx)
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
  for (i = 1; i <= L; i++) {
    mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;

    for (k = 1; k <= hmm->M; k++) {
				/* match state */
      mmx[i][k]  = mmx[i-1][k-1] + hmm->tsc[k-1][TMM];
      if ((sc = imx[i-1][k-1] + hmm->tsc[k-1][TIM]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = xmx[i-1][XMB] + hmm->bsc[k]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = dmx[i-1][k-1] + hmm->tsc[k-1][TDM]) > mmx[i][k])
	mmx[i][k] = sc;
      mmx[i][k] += hmm->msc[(int) dsq[i]][k];

				/* delete state */
      dmx[i][k]  = mmx[i][k-1] + hmm->tsc[k-1][TMD];
      if ((sc = dmx[i][k-1] + hmm->tsc[k-1][TDD]) > dmx[i][k])
	dmx[i][k] = sc;

				/* insert state */
      if (k < hmm->M) {
	imx[i][k] = mmx[i-1][k] + hmm->tsc[k][TMI];
	if ((sc = imx[i-1][k] + hmm->tsc[k][TII]) > imx[i][k])
	  imx[i][k] = sc;
	imx[i][k] += hmm->isc[(int) dsq[i]][k];
      }
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
				/* N state */
    xmx[i][XMN] = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP];

				/* E state */
    xmx[i][XME] = -INFTY;
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  mmx[i][k] + hmm->esc[k]) > xmx[i][XME])
	xmx[i][XME] = sc;
				/* J state */
    xmx[i][XMJ] = xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP];
    if ((sc = xmx[i][XME]   + hmm->xsc[XTE][LOOP]) > xmx[i][XMJ])
      xmx[i][XMJ] = sc;

				/* B state */
    xmx[i][XMB] = xmx[i][XMN] + hmm->xsc[XTN][MOVE];
    if ((sc = xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[i][XMB])
      xmx[i][XMB] = sc;

				/* C state */
    xmx[i][XMC] = xmx[i-1][XMC] + hmm->xsc[XTC][LOOP];
    if ((sc = xmx[i][XME] + hmm->xsc[XTE][MOVE]) > xmx[i][XMC])
      xmx[i][XMC] = sc;
  }
				/* T state (not stored) */
  sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

  if (ret_mx != NULL) *ret_mx = mx;
  else                FreePlan7Matrix(mx);

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
  i    = N;			/* current i for tpos -1*/

  /* Traceback
   */
  while (tr->statetype[tpos-1] != STS) {
    switch (tr->statetype[tpos-1]) {
    case STM:			/* M connects from i-1,k-1, or B */
      sc = mmx[i][k] - hmm->msc[(int) dsq[i]][k];
      if (sc == xmx[i-1][XMB] + hmm->bsc[k])
	{
	  tr->statetype[tpos] = STB;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	  i--;
	}
      else if (sc == mmx[i-1][k-1] + hmm->tsc[k-1][TMM])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = --k;
	  tr->pos[tpos]       = --i;
	}
      else if (sc == imx[i-1][k-1] + hmm->tsc[k-1][TIM])
	{
	  tr->statetype[tpos] = STI;
	  tr->nodeidx[tpos]   = --k;
	  tr->pos[tpos]       = --i;
	}
      else if (sc == dmx[i-1][k-1] + hmm->tsc[k-1][TDM])
	{
	  tr->statetype[tpos] = STD;
	  tr->nodeidx[tpos]   = --k;
	  tr->pos[tpos]       = 0;
	  i--;
	}
      else Die("traceback failed");
      break;

    case STD:			/* D connects from M,D */
      if (dmx[i][k] == mmx[i][k-1] + hmm->tsc[k-1][TMD])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = --k;
	  tr->pos[tpos]       = i;
	}
      else if (dmx[i][k] == dmx[i][k-1] + hmm->tsc[k-1][TDD]) 
	{
	  tr->statetype[tpos] = STD;
	  tr->nodeidx[tpos]   = --k;
	  tr->pos[tpos]       = 0;
	}
      else Die("traceback failed");
      break;

    case STI:			/* I connects from M,I */
      sc = imx[i][k] - hmm->isc[(int) dsq[i]][k];
      if (sc == mmx[i-1][k] + hmm->tsc[k][TMI])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = --i;
	}
      else if (sc == imx[i-1][k] + hmm->tsc[k][TII])
	{
	  tr->statetype[tpos] = STI;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = --i;
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
      else if (i > 0 && xmx[i][XMN] == xmx[i-1][XMN] + hmm->xsc[XTN][LOOP])
	{
	  tr->statetype[tpos] = STN;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos-1]     = --i;
	  tr->pos[tpos]       = 0;
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
	    tr->statetype[tpos] = STM;
	    tr->nodeidx[tpos]   = k;
	    tr->pos[tpos]       = i;
	    break;
	  }
      if (k < 1) Die("traceback failed");
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
  
  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
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
