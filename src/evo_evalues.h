/* P7_EVOPIPELINE is the standardized pipeline for one profile/sequence
 * comparison, from the fast filters down through domain postprocessing,
 * alignment, and scoring.
 */

#ifndef P7_EVOEVALUES_INCLUDED
#define P7_EVOEVALUES_INCLUDED

#include <p7_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"

#include "evohmmer.h"

/* evo_evalues.c */
extern int p7_EvoCalibrate(P7_RATE *R, P7_HMM *hmm, P7_BUILDER *cfg_b, ESL_RANDOMNESS **byp_rng, P7_BG **byp_bg, P7_PROFILE **byp_gm, P7_OPROFILE **byp_om,
			   int noevo_msv, int noevo_vit, int noevo_fwd, float fixtime, float tol);
extern int p7_EvoMSVMu    (ESL_RANDOMNESS *r, P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg,
			   int L, int N, double lambda,               double *ret_mmu, float fixtime, int noevo, float tol);
extern int p7_EvoViterbiMu(ESL_RANDOMNESS *r, P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg,
			   int L, int N, double lambda,                double *ret_vmu, float fixtime, int noevo, float tol);
extern int p7_EvoTau      (ESL_RANDOMNESS *r, P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, 
			   int L, int N, double lambda, double tailp, double *ret_tau, float fixtime, int noevo, float tol);


#endif /*P7_EVOEVALUES_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
