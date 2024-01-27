/* calibrate.c : determine parameters for E-value statistics
 */
#include <h4_config.h>

#include "easel.h"
#include "esl_dsq.h"
#include "esl_gumbel.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "h4_checkptmx.h"
#include "h4_filtermx.h"
#include "h4_mode.h"
#include "h4_profile.h"

#include "fbfilter.h"
#include "ssvfilter.h"
#include "vitfilter.h"

/* Function:  h4_lambda()
 * Synopsis:  Determine lambda (slope) parameter
 * Incept:    SRE, Tue 02 Jan 2024 [Ave, Zaragoza-Madrid]
 *
 * Purpose:   Determine the $\hat{\lambda}$ parameter to use for profile <hmm>.
 *            This will be used for Viterbi Gumbel distributions (including
 *            SSV and Viterbi acceleration filters) and for Forward exponential
 *            tails.
 *
 * Args:      hmm        : model to calculate lambda for
 *            ret_lambda : RETURN: estimated lambda
 *
 * Returns:   <eslOK> on success and <*ret_lambda> is the result.
 */
int
h4_lambda(const H4_PROFILE *hmm, double *ret_lambda)
{
  *ret_lambda = log(2.0);
  return eslOK;
}


/* Function:  h4_ssv_mu()
 * Synopsis:  Determine the SSV Gumbel mu location parameter for a profile
 * Incept:    SRE, Mon 01 Jan 2024 [Zaragoza]
 *
 * Purpose:   Given profile <hmm>, determine the Gumbel location parameter
 *            $\mu$ for SSV scores using a brief simulation, and return that
 *            parameter in <ret_smu>.
 *
 *            The simulation generates <N> random sequences of length
 *            <L> using the background frequencies set for the profile
 *            (<hmm->f>) and random number generator <rng>.  The $\mu$
 *            parameter is determined by maximum likelihood fitting
 *            the <N> SSV filter scores to a Gumbel distribution with
 *            a given assumed <lambda>.
 *
 * Args:      rng     : random number generator
 *            hmm     : profile (with vectorized scores set)
 *            L       : length of sequences to simulate
 *            N       : number of sequences to simulate
 *            lambda  : known Gumbel lambda parameter
 *            ret_smu : RETURN: \hat{mu} estimated SSV location parameter
 *
 * Returns:   <eslOK> on success, and <ret_smu> contains the maximum likelihood
 *            estimate of the SSV mu parameter.
 *
 * Throws:    <eslEMEM> on allocation failures.
 *            <eslEINVAL> if profile isn't vectorized.
 */
int
h4_ssv_mu(ESL_RANDOMNESS *rng, const H4_PROFILE *hmm, int L, int N, double lambda, double *ret_smu)
{
  ESL_DSQ   *dsq = NULL;
  double    *xv  = NULL;
  H4_MODE   *mo  = NULL;
  float      sc;
  int        i;
  int        status;

  if (!(hmm->flags & h4_HASVECS)) ESL_XEXCEPTION(eslEINVAL, "profile not vectorized");

  if ((mo = h4_mode_Create()) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(xv,  sizeof(double)  * N);
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

  // SSV only uses mo->nullsc. It doesn't use any H4_MODE transitions.
  if ((status = h4_mode_SetLength(mo, L)) != eslOK) goto ERROR;

  for (i = 0; i < N; i++)
    {
      if ((status = esl_rsq_xfIID(rng, hmm->f, hmm->abc->K, L, dsq)) != eslOK) goto ERROR;

      status = h4_ssvfilter(dsq, L, hmm, mo, &sc); 
      if (status != eslOK && status != eslERANGE) goto ERROR;  
      // eslERANGE means we overflowed SSV filter score range, and sc is only a lower bound

      xv[i] = sc;                       // H4 scores are floats, but esl_gumbel takes doubles.
    }
  if ((status = esl_gumbel_FitCompleteLoc(xv, N, lambda, ret_smu)) != eslOK) goto ERROR;
  /* status = OK flows through to cleanup and exit; *ret_smu was just set */
 ERROR:
  h4_mode_Destroy(mo);
  free(xv);
  free(dsq);
  return status;
}


/* Function:  h4_vit_mu()
 * Synopsis:  Determine the local Viterbi Gumbel mu parameter for a profile
 * Incept:    SRE, Mon 08 Jan 2024 [Jay-Z, Numb/Encore]
 *
 * Purpose:   Identical to h4_ssv_mu(), above, except that it fits
 *            Viterbi scores instead of SSV scores. 
 *
 *            The difference between the two mus is small, but can be
 *            up to ~1 bit or so for large, low-info models [J4/126] so
 *            decided to calibrate the two mus separately [J5/8]. 
 *
 * Args:      rng     :  source of random numbers
 *            hmm     :  profile (with vectorized scores set)
 *            L       :  length of sequences to simulate
 *            N	      :  number of sequences to simulate		
 *            lambda  :  known Gumbel lambda parameter to use
 *            ret_vmu :  RETURN: ML estimate of location param mu
 *
 * Returns:   <eslOK> on success, and <ret_mu> contains the ML estimate
 *            of $\mu$.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if profile isn't vectorized.
 */
int
h4_vit_mu(ESL_RANDOMNESS *rng, const H4_PROFILE *hmm, int L, int N, double lambda, double *ret_vmu)
{
  H4_FILTERMX  *fx  = h4_filtermx_Create(hmm->M); 
  ESL_DSQ      *dsq = NULL;
  double       *xv  = NULL;
  H4_MODE      *mo  = NULL;
  float        sc;
  int          i;
  int          status;

  if (!(hmm->flags & h4_HASVECS)) ESL_XEXCEPTION(eslEINVAL, "profile not vectorized");
  if (!fx)                        { status = eslEMEM; goto ERROR; }

  if ((mo = h4_mode_Create()) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(xv,  sizeof(double)  * N);
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

  if ((status = h4_mode_SetLength(mo, L)) != eslOK) goto ERROR;

  for (i = 0; i < N; i++)
    {
      if ((status = esl_rsq_xfIID(rng, hmm->f, hmm->abc->K, L, dsq)) != eslOK) goto ERROR;

      status = h4_vitfilter(dsq, L, hmm, mo, fx, &sc); 
      if (status != eslOK && status != eslERANGE) goto ERROR;  
      if (sc == -eslINFINITY) ESL_XEXCEPTION(eslEINCONCEIVABLE, "you told me the viterbi filter doesn't underflow");
      // eslERANGE (and a finite sc) means we overflowed SSV filter score range, and sc is only a lower bound

      xv[i] = sc;                       // H4 scores are floats, but esl_gumbel takes doubles.
    }
  if ((status = esl_gumbel_FitCompleteLoc(xv, N, lambda, ret_vmu)) != eslOK) goto ERROR;
  /* status = OK flows through to cleanup and exit; *ret_vmu was just set */
 ERROR:
  h4_filtermx_Destroy(fx);
  h4_mode_Destroy(mo);
  free(xv);
  free(dsq);
  return status;
}


/* Function:  h4_fwd_tau()
 * Synopsis:  Determine Forward tau location parameter for a profile
 *
 * Purpose:   Determine $\tau$ location parameter for exponential fit to
 *            Forward score distribution for profile <hmm>.  Done by
 *            fitting scores for <N> i.i.d. random sequences of length
 *            <L> and default composition (in an H4_MODE). Focus
 *            fitting procedure on high-scoring tail of top fraction
 *            of <tailp> scores, using previously determined slope
 *            parameter <lambda>. Return estimated tau in <*ret_tau>.
 *
 *            The $\tau$ parameter is for a complete exponential
 *            distribution Exp($\tau$, $\lambda$), but the fit is only
 *            accurate in the extreme tail of the Forward score
 *            distribution, from about $P(s>x) < 0.001$ or so.
 *
 *            Because we can't afford the time to generate enough
 *            samples N to collect enough points (0.001 * N) in this
 *            high-scoring tail for an actual exponential fit, we use
 *            a heuristic method that works for much smaller N. We fit
 *            all N scores to a complete Gumbel, determine the
 *            location of the base of this Gumbel's high-scoring tail
 *            at a high <tailp> (typically ~0.04), and then
 *            extrapolate that base location back to the base of a
 *            complete exponential (backing it up by
 *            log(tailp)/lambda).
 *
 *            This heuristic works because a delicate choice of
 *            <tailp> will roughly cancel two opposing systematic
 *            errors.  On the one hand, the base location for the top
 *            <tailp> of the *actual* distribution of Forward scores
 *            overestimates <tau> (because the empirical distribution
 *            asymptotically approaches the exponential tail from
 *            above). On the other hand, the Gumbel fit systematically
 *            errs toward a steeper tail (higher fitted lambda), so
 *            the fitted Gumbel distribution *underestimates* the base
 *            location for the top <tailp> scores. The choice of
 *            <tailp> ~ 0.03-0.04 was empirically and carefully tuned
 *            to work for any model. [SRE:J1/134-135].
 *
 *            Typical choices are L=100, N=200, tailp=0.04, yielding
 *            estimates $\hat{\tau}$ with a precision (standard
 *            deviation) of $\pm$ 0.2 bits, corresponding to a $\pm$
 *            15\% error in E-values. [SRE:J1/135].
 *            
 *            
 *
 * Args:      rng     :  source of random numbers
 *            hmm     :  profile (with vectorized scores set)
 *            L       :  length of sequences to simulate
 *            N	      :  number of sequences to simulate		
 *            lambda  :  known lambda parameter to use
 *            tailp   :  tail mass to fit (typically 0.04)
 *            ret_tau :  RETURN: estimate of location param tau
 *
 * Returns:   <eslOK> on success, and <ret_mu> contains the ML estimate
 *            of $\mu$.
 *
 * Throws:    <eslEMEM> on allocation failure, and <*ret_tau> is -eslINFINITY.
 */
int
h4_fwd_tau(ESL_RANDOMNESS *rng, const H4_PROFILE *hmm, int L, int N, double lambda, double tailp, double *ret_tau)
{
  H4_MODE      *mo  = NULL;
  H4_CHECKPTMX *cpx = NULL;
  ESL_DSQ      *dsq = NULL;
  double       *xv  = NULL;
  float         sc;
  double        gmu, glam;
  int           i;
  int           status;

  if ((mo  = h4_mode_Create())                               == NULL) { status = eslEMEM; goto ERROR; }
  if ((cpx = h4_checkptmx_Create(hmm->M, L, ESL_MBYTES(32))) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(xv,  sizeof(double)  * N);
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

  if ((status = h4_mode_SetLength(mo, L)) != eslOK) goto ERROR;

  for (i = 0; i < N; i++)
    {
      if ((status = esl_rsq_xfIID(rng, hmm->f, hmm->abc->K, L, dsq)) != eslOK) goto ERROR;
      if ((status = h4_fwdfilter(dsq, L, hmm, mo, cpx, &sc))         != eslOK) goto ERROR; 
      xv[i] = sc;                       // H4 scores are floats, but esl_gumbel takes doubles.
    }
      
  if ((status = esl_gumbel_FitComplete(xv, N, &gmu, &glam)) != eslOK) goto ERROR;  // this is the heuristic. complete Gumbel fit...
  *ret_tau = esl_gumbel_invcdf(1.0-tailp, gmu, glam) + (log(tailp) / lambda);      // ... then find the location of its tail... then back up to complete exponential base location.
  // flowthrough to cleanup... status == eslOK from the fit above
 ERROR:
  if (status != eslOK) *ret_tau = -eslINFINITY;
  h4_mode_Destroy(mo);
  h4_checkptmx_Destroy(cpx);
  free(xv);
  free(dsq);
  return status;
}


/*****************************************************************
 * x. Example driver
 *****************************************************************/
#ifdef h4CALIBRATE_EXAMPLE

#include <string.h>

#include "esl_getopts.h"
#include "esl_stopwatch.h"

#include "general.h"
#include "h4_hmmfile.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version info",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char banner[] = "example of E-value calibrations";
static char usage[]  = "[options] <hmmfile>";


int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  char           *hmmfile = esl_opt_GetArg(go, 1);
  int             N       = esl_opt_GetInteger(go, "-N");
  int             L       = esl_opt_GetInteger(go, "-L");
  H4_HMMFILE     *hfp     = NULL;
  ESL_ALPHABET   *abc     = NULL;
  H4_PROFILE     *hmm     = NULL;
  double          lambda;
  double          smu;
  int             status;

  status = h4_hmmfile_Open(hmmfile, NULL, &hfp);
  if (status != eslOK) esl_fatal("Error: failed to open %s for reading profile HMM(s)\n%s\n", strcmp(hmmfile, "-") == 0? "<stdin>" : hmmfile, hfp->errmsg);

  status = h4_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT) esl_fatal("Parse failed, bad profile HMM file format in %s:\n   %s",  strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status == eslEOF)     esl_fatal("Empty input? No profile HMM found in %s\n",                strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status != eslOK)      esl_fatal("Unexpected error reading profile HMM from %s (code %d)\n", strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, status);

  esl_stopwatch_Start(w);

  if ( h4_lambda(hmm, &lambda)                 != eslOK) esl_fatal("h4_lambda() failed");
  if ( h4_ssv_mu(rng, hmm, L, N, lambda, &smu) != eslOK) esl_fatal("h4_ssv_mu() failed");

  printf("lambda = %.4f\n", lambda);
  printf("SSV mu = %.2f\n", smu);

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(rng);
  esl_stopwatch_Destroy(w);
  h4_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  h4_profile_Destroy(hmm);
  return 0;
}
#endif // h4CALIBRATE_EXAMPLE
