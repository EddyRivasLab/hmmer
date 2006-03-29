/* sample.c
 * 
 * Sample alignments from the posterior distribution, given an HMM and
 * a sequence.
 * SVN $Id:
 */

#include "config.h"		/* compile-time configuration constants */
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include "squid.h"		/* general sequence analysis library    */

#include "plan7.h"		/* plan 7 profile HMM structure         */
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "easel.h"
#include "esl_random.h"


/* Function: _weighted_choice()
 * Date:     Alex Coventry, Tue Mar 28 3:22:10 2006 (St. Louis)
 *
 * Purpose:  Choose an array index at random, with probability
 *           according to the relative log likelihoods given by the
 *           entries in the array.
 *
 * Args:     randomness - Placeholder for RNG.
 *           N          - Length of array
 *           weights    - Array of weights.  The function monkeys
 *                        with this array.
 *
 * Returns:  The index of the chosen weight.
 */
int
_weighted_choice(ESL_RANDOMNESS *randomness, int N,
		 double *weights) {
  int weight_idx;
  double max_weight, total_weight;

  /* Find the largest weight. */
  max_weight = -INFTY;
  for (weight_idx = 0; weight_idx < N; weight_idx++)
    if (weights[weight_idx] > max_weight) 
      max_weight = weights[weight_idx];

  DNorm(weights, N);
  return esl_rnd_DChoose(randomness, weights, N);  
}
 

/* Function: weighted_choice()
 * Date:     Alex Coventry, Mon Mar 27 4:20:25 2006 (St. Louis)
 *
 * Purpose: Choose an index at random, with probability according to
 *          the relative log likelihoods given by the numerical
 *          entries in the function arguments (which must be doubles.)
 *
 * Args:    The ESL_RANDOMNESS placekeeper, followed by the weights as
 *          doubles, followed by 1e308 (as a sentinel.  C sucks.)
 *
 * Returns: The index of the chosen weight.
 */
int
weighted_choice(ESL_RANDOMNESS *randomness, ...) {
  va_list weight_args;     /* Tracks the variable arguments passed */
  double  weights[100];    /* Array of weights to choose from. */
  int     N;               /* Number of weights */

  /* Take a copy of the weight arguments. */
  va_start(weight_args, randomness);
  N = 0;
  while (N < 100) {
    weights[N] = va_arg(weight_args, double);
    if (weights[N] == 1e308) break; 
    N++;
  }
  va_end(weight_args);
  if (N >= 100) Die("Too many weights.");
  return _weighted_choice(randomness, N, weights);
}

/* Function: P7SampleAlignment()
 * Date:     Alex Coventry, Mon Mar 27 4:09:30 2006 (St. Louis)
 * 
 * Purpose: Given a forward matrix, return a sampled alignment, in
 *          traceback form.  Modeled On
 *          core_algorithms.c:P7ViterbiTrace, subversion revision
 *          number 1508.
 *           
 * Args:     hmm        - hmm, log odds form, used to make mx
 *           dsq        - sequence aligned to (digital form) 1..N  
 *           N          - length of seq
 *           mx         - the matrix to trace back in, N x hmm->M.  Should
 *                        contain forward scores.
 *           ret_tr     - RETURN: traceback.
 *           randomness - Placekeeper for RNG.
 *           
 * Return:   (void)
 *           ret_tr is allocated here. Free using P7FreeTrace().
 */
void
P7SampleAlignment(struct plan7_s *hmm, unsigned char *dsq, int N,
		  struct dpmatrix_s *mx, struct p7trace_s **ret_tr,
		  ESL_RANDOMNESS *randomness)
{
  struct p7trace_s *tr;
  int curralloc;		/* current allocated length of trace */
  int tpos;			/* position in trace */
  int i;			/* position in seq (1..N) */
  int k;			/* position in model (1..M) */
  int **xmx, **mmx, **imx, **dmx;
  int trace_choice;             /* For random choice of traceback */
  double current_null;          /* Null score for the current character. */
  double *match_odds;           /* Keeps the log odds scores for us. */

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

    /* The scores in the forward matrix include the null correction.
       When we consider moves in the forward matrix corresponding to
       moves in the sequence (i.e., those which involve a change in
       the value of "i",) it is important to correct their scores by
       the current null so that they can be fairly compared to moves
       which don't correspond to a move in the sequence space. */
    current_null = sreLOG2(hmm->null[dsq[i]]) + sreLOG2(hmm->p1);
    switch (tr->statetype[tpos-1]) {
    case STM:			/* M connects from i-1,k-1, or B */
      trace_choice = weighted_choice(randomness,
				     xmx[i][XMB] + hmm->bsc[k+1],
				     mmx[i][k] + hmm->tsc[TMM][k] - current_null,
				     imx[i][k] + hmm->tsc[TIM][k] - current_null,
				     dmx[i][k] + hmm->tsc[TDM][k],
				     1e308);
      switch (trace_choice) {
      case 0: /* xmx[i][XMB] + hmm->bsc[k+1] */
	  tr->statetype[tpos] = STB;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	  break;
      case 1: /* mmx[i][k] + hmm->tsc[TMM][k] */
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	  break;
      case 2: /* imx[i][k] + hmm->tsc[TIM][k] */
	  tr->statetype[tpos] = STI;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = i--;
	  break;
      case 3: /* dmx[i][k] + hmm->tsc[TDM][k] */
	  tr->statetype[tpos] = STD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	  break;
      }
      break;

    case STD:			/* D connects from M,D */
      trace_choice = weighted_choice(randomness,
				     mmx[i][k] + hmm->tsc[TMD][k] - current_null,
				     dmx[i][k] + hmm->tsc[TDD][k],
				     1e308);
      switch (trace_choice) {
      case 0: /* mmx[i][k] + hmm->tsc[TMD][k] */
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	  break;
      case 1: /* dmx[i][k] + hmm->tsc[TDD][k] */
	  tr->statetype[tpos] = STD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	  break;
      }
      break;

    case STI:			/* I connects from M,I */
      trace_choice = weighted_choice(randomness,
				     mmx[i][k] + hmm->tsc[TMI][k] - current_null,
				     imx[i][k] + hmm->tsc[TII][k] - current_null,
				     1e308);
      switch (trace_choice) {
      case 0: /* mmx[i][k] + hmm->tsc[TMI][k] */
	tr->statetype[tpos] = STM;
	tr->nodeidx[tpos]   = k--;
	tr->pos[tpos]       = i--;
	break;
      case 1: /* imx[i][k] + hmm->tsc[TII][k] */
	tr->statetype[tpos] = STI;
	tr->nodeidx[tpos]   = k;
	tr->pos[tpos]       = i--;
	break;
      }
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
      trace_choice = weighted_choice(randomness,
				     xmx[i][XMN] + hmm->xsc[XTN][MOVE],
				     xmx[i][XMJ] + hmm->xsc[XTJ][MOVE],
				     1e308);
      switch (trace_choice) {
      case 0:
	tr->statetype[tpos] = STN;
	tr->nodeidx[tpos]   = 0;
	tr->pos[tpos]       = 0;
	break;
      case 1:
	tr->statetype[tpos] = STJ;
	tr->nodeidx[tpos]   = 0;
	tr->pos[tpos]       = 0;
	break;
      }
      break;

    case STE:			/* E connects from any M state. k set here */

      /* Have to do the work of weighted_choice by hand, here, since
	 we don't know how many match states there are, and it would
	 be a bore to write them all out even if we did.  Fill out the
	 weights matrix. */
      match_odds = (double *)MallocOrDie(sizeof(double)*hmm->M+1);
      match_odds[0] = -INFTY;   /* Match-state indexing starts at 1, not 0. */
      for (k = hmm->M; k >= 1; k--)

	/* Not necessary to subtract current_null, here, since *all*
	 * the candidates involve a move in sequence space, so it'll
	 * cancel out in the normalizations. */
	match_odds[k] = mmx[i][k] + hmm->esc[k];
      k = _weighted_choice(randomness, hmm->M, match_odds);
      tr->statetype[tpos] = STM;
      tr->nodeidx[tpos]   = k--;
      tr->pos[tpos]       = i--;
      break;

    case STC:			/* C comes from C, E */
      trace_choice = weighted_choice(randomness,
				     xmx[i-1][XMC] + hmm->xsc[XTC][LOOP] - current_null,
				     xmx[i][XME] + hmm->xsc[XTE][MOVE],
				     1e308);
      switch (trace_choice) {
      case 0:  /* xmx[i-1][XMC] + hmm->xsc[XTC][LOOP] */
	tr->statetype[tpos] = STC;
	tr->nodeidx[tpos]   = 0;
	tr->pos[tpos]       = 0;    /* note convention adherence: */
	tr->pos[tpos-1]     = i--;  /* first C doesn't emit       */
	break;
      case 1: /* xmx[i][XME] + hmm->xsc[XTE][MOVE] */
	tr->statetype[tpos] = STE;
	tr->nodeidx[tpos]   = 0;
	tr->pos[tpos]       = 0; /* E is a nonemitter */
	break;
      }
      break;

    case STJ:			/* J connects from E, J */
      trace_choice = weighted_choice(randomness,
				     xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP] - current_null,
				     xmx[i][XME] + hmm->xsc[XTE][LOOP],
				     1e308);
      switch (trace_choice) {
      case 0: /* xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP] */
	tr->statetype[tpos] = STJ;
	tr->nodeidx[tpos]   = 0;
	tr->pos[tpos]       = 0;    /* note convention adherence: */
	tr->pos[tpos-1]     = i--;  /* first J doesn't emit       */
	break;
      case 1: /* xmx[i][XME] + hmm->xsc[XTE][LOOP] */
	tr->statetype[tpos] = STE;
	tr->nodeidx[tpos]   = 0;
	tr->pos[tpos]       = 0; /* E is a nonemitter */
	break;
      }
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

