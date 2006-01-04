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

/* Function:  P7CalibrateV()
 * Incept:    SRE, Wed Nov 23 19:41:54 2005 [St. Louis]
 *
 * Purpose:   Sample <N> iid random sequences of length <L> using the
 *            residue frequencies <fq>; perform Viterbi alignments to
 *            the model <hmm>; collect scores and optimal alignment
 *            lengths. Fit the scores to a Gumbel distribution with
 *            fixed $\lambda=0.693$, and determine a ML estimate of
 *            $\mu$; return this estimate in <ret_mu>. Calculate the
 *            mean alignment length; return this in <ret_kappa>.
 *            
 *            The model is used exactly as provided; the caller is
 *            responsible for configuring in the desired mode with
 *            <P7Config()>, setting the edge correction parameter
 *            <hmm->kappa>, and configuring the length model with
 *            <P7ReconfigLength()>.
 *            
 *            Typically, to determine an <hmm->mu> Gumbel parameter,
 *            the model would be configured in sw mode, <hmm->kappa>
 *            would be set first by a prior call to <P7CalibrateV()>,
 *            and the length model would be configured to length <L>.
 *            
 *            To determine a new <hmm->kappa> parameter, the model
 *            would be configured in the desired mode, <hmm->kappa>
 *            would be set to 0, and the length model would be 
 *            configured to length <L>.
 *            
 *            Caller provides random number source <r> for iid seq
 *            generation. This way, the caller controls whether the
 *            simulation is reproducible or not.
 *            
 * Args:      r         - source of random numbers
 *            hmm       - the profile HMM
 *            N         - number of iid random sequences to sample
 *            L         - length of iid random sequences
 *            ret_mu    - optRETURN: fitted Gumbel ML mu estimate
 *            ret_kappa - optRETURN: mean alignment length         
 *
 * Returns:   (void).
 */
void
P7CalibrateV(ESL_RANDOMNESS *r, struct plan7_s *hmm, double *fq, int N, int L,
	     float *ret_mu, float *ret_kappa)
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
  len = malloc(sizeof(double) * N);

  if (ViterbiSize(L, hmm->M) <= RAMLIMIT)
    mx = CreateDPMatrix(L, hmm->M, 0, 0);  /* full matrix  */
  else
    mx = CreateDPMatrix(1, hmm->M, 25, 0); /* initially 2 rows, growable */

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
      
      sc[i] = DispatchViterbi(dsq, L, hmm, mx, &tr, TRUE);
      
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


  FreeDPMatrix(mx);
  free(sc);
  free(len);

  if (ret_mu    != NULL) *ret_mu    = (float) mu;
  if (ret_kappa != NULL) *ret_kappa = (float) kappa;
  return;
}
 
/************************************************************
 * @LICENSE@
 ************************************************************/
