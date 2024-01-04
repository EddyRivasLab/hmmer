/* calibrate.c : determine parameters for E-value statistics
 */
#include <h4_config.h>

#include "easel.h"
#include "esl_dsq.h"
#include "esl_gumbel.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "h4_mode.h"
#include "h4_profile.h"

#include "ssvfilter.h"

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

  if (!(hmm->flags & h4_HASVECS)) esl_fatal("profile not vectorized");

  if ((mo = h4_mode_Create()) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(xv,  sizeof(double)  * N);
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

  // SSV only uses mo->nullsc. It doesn't use any H4_MODE transitions.
  if ((status = h4_mode_SetLength(mo, L)) != eslOK) goto ERROR;

  for (i = 0; i < N; i++)
    {
      if ((status = esl_rsq_xfIID(rng, hmm->f, hmm->abc->K, L, dsq)) != eslOK) goto ERROR;

      status = h4_ssvfilter(dsq, L, hmm, &sc); 
      if (status != eslOK && status != eslERANGE) goto ERROR;  
      // eslERANGE means we overflowed SSV filter score range, and sc is only a lower bound

      sc -= mo->nullsc;                 // sc is now the log-odds bitscore
      xv[i] = sc;                       // H4 scores are floats, but esl_gumbel takes doubles.
    }
  if ((status = esl_gumbel_FitCompleteLoc(xv, N, lambda, ret_smu)) != eslOK) goto ERROR;

 ERROR:
  h4_mode_Destroy(mo);
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
