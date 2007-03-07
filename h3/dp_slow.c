/* Reference implementation of the DP algorithms.
 * 
 * This implementation is intended for clarity, not speed. It is
 * deliberately not optimized. It is provided as a reference
 * implementation of the computationally intensive DP algorithms of
 * HMMER. It serves both as documentation and as a regression test
 * target for developing optimized code.
 * 
 * Contents:
 *     1. Viterbi and its traceback.
 *     2. Unit tests.
 *     3. Test driver.
 * 
 * SRE, Tue Jan 30 10:49:43 2007 [at Einstein's in St. Louis]
 * SVN $Id$
 */

#include "p7_config.h"
#include <easel.h>
#include "hmmer.h"


/*****************************************************************
 * 1. Viterbi and its traceback.
 *****************************************************************/

/* Function: p7_Viterbi()
 * Incept:   SRE, Tue Jan 30 10:50:53 2007 [Einstein's, St. Louis]
 * 
 * Purpose:  The standard Viterbi dynamic programming algorithm. 
 *
 *           Given a digital sequence <dsq> of length <L>, 
 *           a profile <gm>, and DP matrix <gm> allocated for 
 *           at least <gm->M> by <L> cells;
 *           calculate the maximum scoring path by Viterbi;
 *           if <tr> is provided non-<NULL>, return the path in <tr>;
 *           return the score in <ret_sc> in its internal (SILO) 
 *           form. 
 *           
 *           To convert a SILO score to units of bits, call
 *           p7_SILO2Bitscore().
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
p7_Viterbi(ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_GMX *mx, P7_TRACE *tr, int *ret_sc)
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
  xmx[0][p7_XME] = xmx[0][p7_XMC] = xmx[0][p7_XMJ] = p7_IMPOSSIBLE;  /* need seq to get here */
  for (k = 0; k <= gm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = p7_IMPOSSIBLE;      /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = impossible for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 1; i <= L; i++) {
    mmx[i][0] = imx[i][0] = dmx[i][0] = p7_IMPOSSIBLE;

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
      if (k < gm->M) {
	imx[i][k] = p7_IMPOSSIBLE;
	if ((sc = mmx[i-1][k] + gm->tsc[p7_TMI][k]) > imx[i][k])
	  imx[i][k] = sc;
	if ((sc = imx[i-1][k] + gm->tsc[p7_TII][k]) > imx[i][k])
	  imx[i][k] = sc;

	if (gm->isc[dsq[i]][k] != p7_IMPOSSIBLE) 
	  imx[i][k] += gm->isc[dsq[i]][k];
	else
	  imx[i][k] = p7_IMPOSSIBLE;   
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
    if ((sc = xmx[i][p7_XMN] + gm->xsc[p7_XTN][p7_MOVE]) > p7_IMPOSSIBLE)
      xmx[i][p7_XMB] = sc;
    if ((sc = xmx[i][p7_XMJ] + gm->xsc[p7_XTJ][p7_MOVE]) > xmx[i][p7_XMB])
      xmx[i][p7_XMB] = sc;

				/* C state */
    xmx[i][p7_XMC] = p7_IMPOSSIBLE;
    if ((sc = xmx[i-1][p7_XMC] + gm->xsc[p7_XTC][p7_LOOP]) > p7_IMPOSSIBLE)
      xmx[i][p7_XMC] = sc;
    if ((sc = xmx[i][p7_XME] + gm->xsc[p7_XTE][p7_MOVE]) > xmx[i][p7_XMC])
      xmx[i][p7_XMC] = sc;
  }
				/* T state (not stored) */
  sc = xmx[L][p7_XMC] + gm->xsc[p7_XTC][p7_MOVE];

  if (tr != NULL) p7_ViterbiTrace(dsq, L, gm, mx, tr);
  *ret_sc = sc;
  return eslOK;
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
  int i;			/* position in seq (1..L) */
  int k;			/* position in model (1..M) */
  int **xmx, **mmx, **imx, **dmx;
  int sc;			/* temp var for pre-emission score */

  if ((status = p7_trace_Reuse(tr)) != eslOK) goto ERROR;
  xmx = mx->xmx;
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;

  /* Initialization.
   * (back to front. ReverseTrace() called later.)
   */
  if ((status = p7_trace_Append(tr, p7_STT, 0, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_Append(tr, p7_STC, 0, 0)) != eslOK) goto ERROR;
  i    = L;			/* next position to explain in seq */

  /* Traceback
   */
  while (tr->st[tr->N-1] != p7_STS) {
    switch (tr->st[tr->N-1]) {
    case p7_STC:		/* C(i) comes from C(i-1) or E(i) */
      if   (xmx[i][p7_XMC] <= p7_IMPOSSIBLE)
	ESL_XEXCEPTION(eslFAIL, "impossible C reached at i=%d", i);

      if (xmx[i][p7_XMC] == xmx[i-1][p7_XMC] + gm->xsc[p7_XTC][p7_LOOP]) {
	tr->i[tr->N-1]    = i--;  /* first C doesn't emit: subsequent ones do */
	status = p7_trace_Append(tr, p7_STC, 0, 0);
      } else if (xmx[i][p7_XMC] == xmx[i][p7_XME] + gm->xsc[p7_XTE][p7_MOVE]) 
	status = p7_trace_Append(tr, p7_STE, 0, 0);
      else ESL_XEXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7_STE:		/* E connects from any M state. k set here */
      if (xmx[i][p7_XME] <= p7_IMPOSSIBLE) 
	ESL_XEXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      for (k = gm->M; k >= 1; k--)
	if (xmx[i][p7_XME] == mmx[i][k] + gm->esc[k]) {
	  status = p7_trace_Append(tr, p7_STM, k, i);
	  break;
	}
      if (k < 0) ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      break;

    case p7_STM:			/* M connects from i-1,k-1, or B */
      sc = mmx[i][k] - gm->msc[dsq[i]][k];
      if (sc <= p7_IMPOSSIBLE) ESL_XEXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      if      (sc == xmx[i-1][p7_XMB] + gm->bsc[k])        status = p7_trace_Append(tr, p7_STB, 0,   0);
      else if (sc == mmx[i-1][k-1] + gm->tsc[p7_TMM][k-1]) status = p7_trace_Append(tr, p7_STM, k-1, i-1);
      else if (sc == imx[i-1][k-1] + gm->tsc[p7_TIM][k-1]) status = p7_trace_Append(tr, p7_STI, k-1, i-1);
      else if (sc == dmx[i-1][k-1] + gm->tsc[p7_TDM][k-1]) status = p7_trace_Append(tr, p7_STD, k-1, 0);
      else ESL_XEXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
      if (status != eslOK) goto ERROR;
      k--; 
      i--;
      break;

    case p7_STD:			/* D connects from M,D at i,k-1 */
      if (dmx[i][k] <= p7_IMPOSSIBLE) ESL_XEXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);

      if      (dmx[i][k] == mmx[i][k-1] + gm->tsc[p7_TMD][k-1]) status = p7_trace_Append(tr, p7_STM, k-1, i);
      else if (dmx[i][k] == dmx[i][k-1] + gm->tsc[p7_TDD][k-1]) status = p7_trace_Append(tr, p7_STD, k-1, 0);
      else ESL_XEXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      if (status != eslOK) goto ERROR;
      k--;
      break;

    case p7_STI:			/* I connects from M,I at i-1,k*/
      sc = imx[i][k] - gm->isc[dsq[i]][k];
      if (sc <= p7_IMPOSSIBLE) ESL_XEXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      if      (sc == mmx[i-1][k] + gm->tsc[p7_TMI][k]) status = p7_trace_Append(tr, p7_STM, k, i-1);
      else if (sc == imx[i-1][k] + gm->tsc[p7_TII][k]) status = p7_trace_Append(tr, p7_STI, k, i-1);
      else ESL_XEXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      if (status != eslOK) goto ERROR;
      i--;
      break;

    case p7_STN:			/* N connects from S, N */
      if (xmx[i][p7_XMN] <= p7_IMPOSSIBLE) ESL_XEXCEPTION(eslFAIL, "impossible N reached at i=%d", i);

      if      (i == 0 && xmx[i][p7_XMN] == 0) 
	status = p7_trace_Append(tr, p7_STS, 0, 0);
      else if (i > 0  && xmx[i][p7_XMN] == xmx[i-1][p7_XMN] + gm->xsc[p7_XTN][p7_LOOP]) {
	tr->i[tr->N-1] = i--;
	status = p7_trace_Append(tr, p7_STN, 0, 0);
      } else ESL_XEXCEPTION(eslFAIL, "N at i=%d couldn't be traced", i);
      if (status != eslOK) goto ERROR;
      break;

    case p7_STB:			/* B connects from N, J */
      if (xmx[i][p7_XMB] <= p7_IMPOSSIBLE) ESL_XEXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if (xmx[i][p7_XMB] == xmx[i][p7_XMN] + gm->xsc[p7_XTN][p7_MOVE]) 
	status = p7_trace_Append(tr, p7_STN, 0, 0);
      else if (xmx[i][p7_XMB] == xmx[i][p7_XMJ] + gm->xsc[p7_XTJ][p7_MOVE])
	status = p7_trace_Append(tr, p7_STJ, 0, 0);
      else  ESL_XEXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7_STJ:			/* J connects from E(i) or J(i-1) */
      if (xmx[i][p7_XMJ] <= p7_IMPOSSIBLE) ESL_XEXCEPTION(eslFAIL, "impossible J reached at i=%d", i);

      if (xmx[i][p7_XMJ] == xmx[i-1][p7_XMJ] + gm->xsc[p7_XTJ][p7_LOOP]) {
	tr->i[tr->N-1] = i--;
	status = p7_trace_Append(tr, p7_STJ, 0, 0);
      } else if (xmx[i][p7_XMJ] == xmx[i][p7_XME] + gm->xsc[p7_XTE][p7_LOOP])
	status = p7_trace_Append(tr, p7_STE, 0, 0);
      else  ESL_XEXCEPTION(eslFAIL, "J at i=%d couldn't be traced", i);
      break;

    default: ESL_XEXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    if (status != eslOK) goto ERROR;
  } /* end traceback, at S state */

  if ((status = p7_trace_Reverse(tr)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  return status;
}


/*****************************************************************
 * 2. Unit tests.
 *****************************************************************/
#ifdef p7DP_SLOW_TESTDRIVE

static void
utest_viterbi(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_PROFILE *gm, int nseq, int L)
{
  int       status;
  char      errbuf[eslERRBUFSIZE];
  ESL_DSQ  *dsq = NULL;
  P7_GMX   *mx  = NULL;
  P7_TRACE *tr  = NULL;
  int       idx;
  int       sc1, sc2;
  
  if ((dsq    = malloc(sizeof(ESL_DSQ) *(L+2))) == NULL)  esl_fatal("malloc failed");
  if ((status = p7_trace_Create(L, &tr))        != eslOK) esl_fatal("trace creation failed");
  if ((mx     = p7_gmx_Create(gm->M, L))        == NULL)  esl_fatal("matrix creation failed");

  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rnd_xfIID(r, gm->bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");
      if (p7_Viterbi(dsq, L, gm, mx, tr, &sc1)        != eslOK) esl_fatal("viterbi failed");
      if (p7_trace_Validate(tr, abc, dsq, errbuf)     != eslOK) esl_fatal("trace invalid:\n%s", errbuf);
      if (p7_trace_Score(tr, dsq, gm, &sc2)           != eslOK) esl_fatal("trace score failed");

      if (sc1 != sc2) esl_fatal("Trace score not equal to Viterbi score");
    }

  p7_gmx_Destroy(mx);
  p7_trace_Destroy(tr);
  free(dsq);
  return;
}

#endif /*p7DP_SLOW_TESTDRIVE*/




/*****************************************************************
 * 3. Test driver.
 *****************************************************************/
/* gcc -g -Wall -Dp7DP_SLOW_TESTDRIVE -I. -I../easel -L. -L../easel -o dp_slow_utest dp_slow.c -lhmmer -leasel -lm
 */
#ifdef p7DP_SLOW_TESTDRIVE
#include "easel.h"

#include "p7_config.h"
#include "hmmer.h"

int
main(int argc, char **argv)
{
  char            errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS *r    = NULL;
  ESL_ALPHABET   *abc  = NULL;
  P7_HMM         *hmm  = NULL;
  int             M    = 100;
  int             L    = 200;
  int             nseq = 20;

  if ((r   = esl_randomness_CreateTimeseeded())    == NULL)  esl_fatal("failed to create rng");
  if ((abc = esl_alphabet_Create(eslAMINO))        == NULL)  esl_fatal("failed to create alphabet");

  if (p7_hmm_Sample(r, M, abc, &hmm)               != eslOK) esl_fatal("failed to sample an HMM");
  if ((hmm->bg = p7_bg_Create(abc))                == NULL)  esl_fatal("failed to create null model");
  if ((hmm->gm = p7_profile_Create(hmm->M, abc))   == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig(hmm, hmm->gm, p7_UNILOCAL)  != eslOK) esl_fatal("failed to config profile");
  if (p7_ReconfigLength(hmm->gm, L)                != eslOK) esl_fatal("failed to config profile length");
  if (p7_hmm_Validate    (hmm,     0.0001, errbuf) != eslOK) esl_fatal("whoops, HMM is bad!");
  if (p7_profile_Validate(hmm->gm, 0.0001)         != eslOK) esl_fatal("whoops, profile is bad!");

  utest_viterbi(r, abc, hmm->gm, nseq, L);
  
  p7_profile_Destroy(hmm->gm);
  p7_bg_Destroy(hmm->bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return 0;
}

#endif /*p7DP_SLOW_TESTDRIVE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
