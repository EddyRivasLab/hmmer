/* calibration.c
 * Calibrating parameters necessary for E-value determination.
 * 
 * SRE, Wed Nov 23 19:41:27 2005
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "config.h"
#include "squidconf.h"
#include "squid.h"

#include "structs.h"
#include "funcs.h"
#include "plan7.h"

#include <easel.h>
#include <esl_stats.h>
#include <esl_random.h>
#include <esl_gumbel.h>
#include <esl_exponential.h>
#include <esl_vectorops.h>

/* Function:  P7CalibrateV()
 * Incept:    SRE, Wed Nov 23 19:41:54 2005 [St. Louis]
 *
 * Purpose:   Calibrates the <hmm> for estimating E-values of local
 *            Viterbi scores. 
 *
 *            Assume that local Viterbi scores follow
 *            a Gumbel distribution with a known $\lambda = \log 2 = 0.693$.
 *            Sample <N> iid random sequences of length
 *            <L> using the residue frequencies <fq>. Perform Viterbi
 *            alignments to the model <hmm> and collect observed
 *            scores.  Fit the scores to a Gumbel distribution of
 *            fixed $\lambda=0.693$, and determine a ML estimate of
 *            $\mu$.
 *            
 *            The model is configured into fs mode, with length
 *            modeling to match <L>. sw mode distribution is assumed
 *            to be offset by +1 bit, and scores are assumed to be
 *            length-independent so long as the model is reconfigured
 *            for each target length.
 *            
 *            Typical values of <N> and <L> are N=100, L=100. It does
 *            not take many <N> to get a reasonable estimate of $\mu$,
 *            and <L> only has to be large enough to avoid significant
 *            edge effects.
 *            
 *            Caller provides random number source <r> for iid seq
 *            generation. This way, the caller controls whether the
 *            simulation is reproducible or not.
 *            
 *            Upon successful return, the <vN>, <vL>, <vmu>, and
 *            <vlambda> parameters are set in the model, and the 
 *            <PLAN7_STATS_LV> flag is raised. Previous calibration,
 *            if any, is overwritten.
 *            
 * Args:      r         - source of random numbers
 *            hmm       - the profile HMM
 *            fq        - background residue frequencies to use
 *            N         - number of iid random sequences [100]
 *            L         - length of iid random sequences [100]
 *
 * Returns:   (void).
 */
void
P7CalibrateV(ESL_RANDOMNESS *r, struct plan7_s *hmm, double *fq, int N, int L)
{
  struct dpmatrix_s *mx;
  int     i;
  char   *seq;
  unsigned char *dsq;
  double *sc;			/* Viterbi scores for N random seqs */
  double  mu, lambda;       	/* Easel works in doubles; HMMER in floats. */
  
  sc  = MallocOrDie(sizeof(double) * N);

  if (P7ViterbiSize(L, hmm->M) <= RAMLIMIT)
    mx = CreatePlan7Matrix(L, hmm->M, 0, 0);  /* full matrix  */
  else
    mx = CreatePlan7Matrix(1, hmm->M, 25, 0); /* initially 2 rows, growable */

  /* Configure the model to calibrate Viterbi local scores
   */
  P7Config(hmm, P7_FS_MODE);
  P7ReconfigLength(hmm, L);

  /* Collect scores and optimal alignment lengths.
   * We could parallelize this section (threads and/or MPI) someday;
   * I'm collecting results into arrays sc[] and len[] to facilitate
   * parallelization later.
   * We might also be able to accelerate here using island statistics.
   */
  for (i = 0; i < N; i++)
    {
      seq = esl_rnd_IID(r, Alphabet, fq, Alphabet_size, L);
      dsq = DigitizeSequence(seq, L); /* better: someday generate dsq directly */
      
      /* We only need the score, so we could optimize this by 
       * making a score-only Viterbi that always uses a two-row DP matrix.
       */
      if (P7ViterbiSpaceOK(L, hmm->M, mx))
	sc[i] = P7Viterbi(dsq, L, hmm, mx, NULL);
      else
      	sc[i] = P7SmallViterbi(dsq, L, hmm, mx, NULL);

      sc[i] += P7ScoreCorrection(hmm, L);

      free(dsq);        
      free(seq);
    }

  hmm->vN      = N;		                               /* remember */
  hmm->vL      = L;		                               /* remember */
  hmm->vlambda = 0.693;		                               /* assert   */
  esl_gumbel_FitCompleteLoc(sc, N, hmm->vlambda, &(hmm->vmu)); /* fit      */
  hmm->flags  |= PLAN7_STATS_LV;

  FreePlan7Matrix(mx);
  free(sc);
  return;
}
 
/* Function:  P7CalibrateF()
 * Incept:    SRE, Wed Jan 25 12:25:20 2006 [St. Louis]
 *
 * Purpose:   Calibrates the <hmm> for estimating E-values of local
 *            Forward scores.
 *
 *            We assume that Forward local alignment
 *            scores follow an exponential tail of known lambda=0.693.
 *            We sample <N> iid random sequences of length
 *            <L> using residue frequencies <fq>. Perform Forward
 *            alignments to the model <hmm> and collect observed
 *            Forward scores. Find a threshold <ret_mu> that defines the
 *            the highest <mass>% of the scores; let this be the $\mu$
 *            parameter for an exponential tail of slope lambda=0.693.
 *            
 *            The model is configured into fs mode, with length
 *            modeling to match <L>. sw mode distribution is assumed
 *            to be offset by +1 bit, and scores are assumed to be
 *            length-independent so long as the model is reconfigured
 *            for each target length.
 *            
 *            Typical values of <mass>, <N>, and <L> are mass=0.01 (1%
 *            tail), N=5000, L=100. These values were empirically
 *            tested [xref STL10/90].
 *            
 *            Caller provides random number source <r> for iid sequence
 *            generation, so the caller can control whether a simulation
 *            is reproducible or not.
 *            
 *            Upon successful return, the <fN>, <fL>, <fmass>, <fmu>,
 *            and <flambda> parameters are set in the model, and the
 *            <PLAN7_STATS_LF> flag is raised. Previous calibration,
 *            if any, is overwritten.
 *            
 * Internal notes:            
 *            Some optimizations and cleanup should be done here.
 *            
 *            - we should generate dsq's directly, rather than the
 *              2-step of generating alphabetic seq, then digitizing it.
 *            - we should have a P7Forward() that generates scores only,
 *              to save memory.
 *            - we should reuse the P7Forward() matrix.
 *            - we should abstract and parallelize score collection.
 *            
 * Args:      r         - source of random numbers
 *            hmm       - the profile HMM
 *            fq        - background residue frequencies to use
 *            mass      - the tail mass to fit an exponential tail to (0.01)
 *            N         - number of iid random sequences to sample    (5000)
 *            L         - length of iid random sequences              (100)
 *
 * Returns:   (void).
 */
extern void
P7CalibrateF(ESL_RANDOMNESS *r, struct plan7_s *hmm, double *fq, 
	      double mass, int N, int L)
{
  int            i;
  char          *seq;
  unsigned char *dsq;
  int            rank;
  float         *sc;

  sc = MallocOrDie(sizeof(float) * N);

  /* Configure the model to calibrate Forward local scores
   */
  P7Config(hmm, P7_FS_MODE);
  P7ReconfigLength(hmm, L);

  for (i = 0; i < N; i++)
    {
      seq = esl_rnd_IID(r, Alphabet, fq, Alphabet_size, L);
      dsq = DigitizeSequence(seq, L); /* better: someday generate dsq directly */
      
      sc[i] = P7Forward(dsq, L, hmm, NULL);
      sc[i] += P7ScoreCorrection(hmm, L);
      
      free(dsq);
      free(seq);
    }
  esl_vec_FSortDecreasing(sc, N);
  rank    = N * mass;

  hmm->fmu     = sc[rank]; 
  hmm->flambda = 0.693;
  hmm->fmass   = mass;
  hmm->fN      = N;
  hmm->fL      = L;
  hmm->flags  |= PLAN7_STATS_LF;

  free(sc);
  return;
}


/* Function:  P7PValueV()
 * Incept:    SRE, Mon Feb  6 09:06:19 2006 [St. Louis]
 *
 * Purpose:   Calculates the P-value (statistical significance) 
 *            for an HMM Viterbi alignment bit score <sc>.
 */
double
P7PValueV(struct plan7_s *hmm, double sc)
{
  double pval, pval2;
  double delta;		/* mu offset in single vs. multihit mode */
  
  /* For any score, in any mode, we have an upper bound
   * on the P-value from a probabilistic (Bayesian) argument.
   */
  if      (sc >= sreLOG2(DBL_MAX))       pval = 0.0;
  else if (sc <= -1. * sreLOG2(DBL_MAX)) pval = 1.0;
  else    pval = 1. / (1.+sreEXP2(sc));

  /* For local modes, we can also take advantage of more
   * accurate empirical calibration: Viterbi scores follow
   * a Gumbel distribution.
   */
  if ((hmm->flags & PLAN7_STATS_LV) &&
      (hmm->mode == P7_SW_MODE || hmm->mode == P7_FS_MODE))
    {
      delta = (hmm->mode == P7_SW_MODE) ? 1.0 : 0.0; /* singlehit offset by +1 */
      pval2 = esl_gumbel_surv(sc, hmm->vmu+delta, hmm->vlambda);
      if (pval2 < pval) pval = pval2;
    }

  return pval;
}

/* Function:  P7PValueF()
 * Incept:    SRE, Mon Feb  6 09:23:36 2006 [St. Louis]
 *
 * Purpose:   Calculates the P-value (statistical significance) 
 *            for an HMM Forward alignment bit score <sc>.
 */
double
P7PValueF(struct plan7_s *hmm, double sc)
{
  double pval, pval2;
  double delta;		/* mu offset in single vs. multihit mode */
  
  /* For any score, in any mode, we have an upper bound
   * on the P-value from a probabilistic (Bayesian) argument.
   */
  if      (sc >= sreLOG2(DBL_MAX))       pval = 0.0;
  else if (sc <= -1. * sreLOG2(DBL_MAX)) pval = 1.0;
  else    pval = 1. / (1.+sreEXP2(sc));

  /* For local modes, we can also take advantage of more
   * accurate empirical calibration: Forward scores have an
   * exponential tail.
   */
  if ((hmm->flags & PLAN7_STATS_LF) &&
      (hmm->mode == P7_SW_MODE || hmm->mode == P7_FS_MODE))
    {
      delta = (hmm->mode == P7_SW_MODE) ? 1.0 : 0.0; /* singlehit offset by +1 */
      pval2 = hmm->fmass * esl_exp_surv(sc, hmm->fmu+delta, hmm->flambda);
      if (pval2 < pval) pval = pval2;
    }

  return pval;
}




/************************************************************
 * @LICENSE@
 ************************************************************/
