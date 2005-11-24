/* calibration.c
 * Calibrating parameters necessary for E-value determination.
 * 
 * SRE, Wed Nov 23 19:41:27 2005
 * SVN $Id$
 */

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

/* Function:  P7CalibrateVL()
 * Incept:    SRE, Wed Nov 23 19:41:54 2005 [St. Louis]
 *
 * Purpose:   Calibrate for Viterbi local modes (fs and sw), assuming
 *            that $\lambda = 0.693$. Sample <N> iid random sequences
 *            of length <L> using the residue frequencies in <fq>. Collect
 *            scores and determine ML estimate of Gumbel $\mu$
 *            parameter.  Collect optimal alignment lengths, and
 *            estimate $\kappa$ parameter for edge effect correction.
 *
 *            Caller provides random number source <r> for iid seq
 *            generation. This way, the caller controls whether the
 *            simulation is reproducible or not.
 *            
 *            The <hmm> must be configured in a local alignment
 *            mode (fs or sw) before calling <P7CalibrateVL()>, by
 *            an appropriate call to <P7Config()>. 
 *            
 *            <hmm->kappa> is set to zero during the calibration.
 *            Edge effect is assumed to be negligible during the
 *            simulation ($L' \simeq L - \kappa \simeq L). 
 *            <hmm->kappa> is then set to the mean optimal alignment
 *            length.
 *            
 * Args:      r    - source of random numbers
 8            hmm  - the profile HMM
 *            N    - number of iid random sequences to sample
 *            L    - length of iid random sequences
 *
 * Returns:   (void).
 *            <hmm->mu>, <hmm->lambda>, and <hmm->kappa> are set,
 *            and the <PLAN7_STATS> flag is raised on the model.
 */
void
P7CalibrateVL(ESL_RANDOMNESS *r, struct plan7_s *hmm, double *fq, int N, int L)
{
  struct p7trace_s  *tr;
  struct dpmatrix_s *mx;
  int     i;
  char   *seq;
  unsigned char *dsq;
  double *sc;			/* Viterbi scores for N random seqs */
  double *len;			/* Viterbi alignment lengths for N random seqs */
  double  mu, lambda, kappa;	/* Easel works in doubles; HMMER in floats. */
  
  sc  = malloc(sizeof(double) * N);
  len = malloc(sizeof(int)    * N);

  if (P7ViterbiSize(L, hmm->M) <= RAMLIMIT)
    mx = CreatePlan7Matrix(L, hmm->M, 0, 0);  /* full matrix  */
  else
    mx = CreatePlan7Matrix(1, hmm->M, 25, 0); /* initially 2 rows, growable */

  hmm->kappa = 0;	/* kappa=0 must precede P7ReconfigLength() */
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
      
      if (P7ViterbiSpaceOK(L, hmm->M, mx))
	sc[i] = P7Viterbi(dsq, L, hmm, mx, &tr);
      else
      	sc[i] = P7SmallViterbi(dsq, L, hmm, mx, &tr);

      len[i] = (double) P7TraceCountAnnotated(tr);

      P7FreeTrace(tr);	/* wasteful, all this trace alloc/dealloc, when all  
			   we need is the # of annotated residues */
      free(dsq);        
      free(seq);
    }

  /* Estimate mu, lambda parameters of a complete Gumbel-distributed dataset
   */
  lambda = 0.693;		/* we assert */
  esl_gumbel_FitCompleteLoc(sc, N, lambda, &mu);

  /* Estimate kappa parameter for edge effect correction, as
   * mean optimal alignment length. (With island statistics, we
   * would probably want to regression fit to score.)
   */
  esl_stats_Mean(len, N, &kappa, NULL);

  hmm->lambda = (float) lambda;
  hmm->mu     = (float) mu;
  hmm->kappa  = (float) kappa;
  free(sc);
  free(len);
  return;
}



 
 
/************************************************************
 * @LICENSE@
 ************************************************************/
