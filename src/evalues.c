/* Calculations and simulations relevant to E-value calculations.
 * 
 * SRE, Mon Aug  6 13:00:06 2007
 * SVN $Id$
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_gumbel.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

#include "hmmer.h"

/* Function:  p7_Lambda()
 * Synopsis:  Determines length-corrected local lambda parameter.
 * Incept:    SRE, Wed Aug  8 17:54:55 2007 [Janelia]
 *
 * Purpose:   Determine the effective scale parameter $\hat{\lambda}$ to
 *            use for model <hmm>. This will be applied both to
 *            Viterbi Gumbel distributions and Forward exponential
 *            tails.
 *            
 *            The 'true' $\lambda$ is always $\log 2 = 0.693$. The effective
 *            lambda is corrected for edge effect, using the equation
 *             
 *             \[
 *                \hat{\lambda} = \lambda + \frac{1.44}{MH}
 *             \]
 *             
 *            where $M$ is the model length and $H$ is the model
 *            relative entropy. The model relative entropy is
 *            approximated by the average relative entropy of match
 *            emission distributions.  The 1.44 is an empirically
 *            determined fudge factor [J1/125]. This edge-effect
 *            correction is based largely on \citep{Altschul01},
 *            except for the fudge factor, which we don't understand
 *            and can't theoretically justify.
 *            
 * Args:      hmm        : model to calculate corrected lambda for
 *            bg         : null model (source of background frequencies)
 *            ret_lambda : RETURN: edge-corrected lambda
 *
 * Returns:   <eslOK> on success, and <*ret_lambda> is the result.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_Lambda(P7_HMM *hmm, P7_BG *bg, double *ret_lambda)
{
  double H = p7_MeanMatchRelativeEntropy(hmm, bg);
  
  *ret_lambda = eslCONST_LOG2 + 1.44 / ((double) hmm->M * H);
  return eslOK;
}


/* Function:  p7_Mu()
 * Synopsis:  Determines the local Viterbi Gumbel mu parameter for a model.
 * Incept:    SRE, Mon Aug  6 13:00:57 2007 [Janelia]
 *
 * Purpose:   Given model <gm> configured for local alignment (typically
 *            multihit, but may be unihit), determine the Gumbel
 *            location parameter $\mu$ by brief simulation. The
 *            simulation generates <N> random sequences of length <L>
 *            using background frequencies in the null model <bg> and
 *            the random number generator <r>; scores them with <gm>
 *            and <bg> with the Viterbi algorithm; and fitting the
 *            resulting distribution to a Gumbel of known <lambda>.
 *            
 *            Typical default choices are L=100, N=200, which gives
 *            $\hat{\mu}$ estimates with precision (standard
 *            deviation) of $\pm$ 0.1 bits, corresponding to an error
 *            of $\pm$ 8\% in E-value estimates. [J1/135].
 *            
 *            This function changes the length configuration of both
 *            <gm> and <bg>. The caller must remember to reconfigure
 *            both of their length models appropriately for any
 *            subsequent alignments.
 *            
 * Args:      r      :  source of random numbers
 *            gm     :  score profile
 *            bg     :  null model
 *            L      :  length of sequences to simulate
 *            N	     :  number of sequences to simulate		
 *            lambda :  known Gumbel lambda parameter
 *            ret_mu :  RETURN: ML estimate of location param mu
 *
 * Returns:   <eslOK> on success, and <ret_mu> contains the ML estimate
 *            of $\mu$.
 *
 * Throws:    (no abnormal error conditions)
 * 
 * Note:      The FitCompleteLoc() function is simple, and it's tempting
 *            to inline it here and save the <xv> working memory. However,
 *            realized that the FitCompleteLoc() function is vulnerable
 *            to under/overflow error, and we'll probably fix it
 *            eventually - need to be sure that fix applies here too.
 */
int
p7_Mu(ESL_RANDOMNESS *r, P7_PROFILE *gm, P7_BG *bg, int L, int N, double lambda, double *ret_mu)
{
  P7_GMX  *gx      = p7_gmx_Create(gm->M, L);	 /* DP matrix: L (rows) will grow as needed */
  ESL_DSQ *dsq     = NULL;
  double  *xv      = NULL;
  int      i;
  float    sc, nullsc;
  int      status;

  if (gx == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(xv,  sizeof(double)  * N);
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

  p7_ReconfigLength(gm, L);
  p7_bg_SetLength(bg, L);

  for (i = 0; i < N; i++)
    {
      if ((status = esl_rsq_xfIID(r, bg->f, gm->abc->K, L, dsq)) != eslOK) goto ERROR;
      if ((status = p7_GViterbi(dsq, L, gm, gx, &sc))            != eslOK) goto ERROR;
      if ((status = p7_bg_NullOne(bg, dsq, L, &nullsc))          != eslOK) goto ERROR;   
      xv[i] = (sc - nullsc) / eslCONST_LOG2;
    }

  if ((status = esl_gumbel_FitCompleteLoc(xv, N, lambda, ret_mu))  != eslOK) goto ERROR;
  p7_gmx_Destroy(gx);
  free(xv);
  free(dsq);
  return eslOK;

 ERROR:
  *ret_mu = 0.0;
  if (gx  != NULL) p7_gmx_Destroy(gx);
  if (xv  != NULL) free(xv);
  if (dsq != NULL) free(dsq);
  return status;

}


/* Function:  p7_Tau()
 * Synopsis:  Determine Forward mu by brief simulation.
 * Incept:    SRE, Thu Aug  9 15:08:39 2007 [Janelia]
 *
 * Purpose:   Determine the <mu> parameter for an exponential tail fit
 *            to the Forward score distribution for model <gm>, on
 *            random sequences with the composition of the background
 *            model <bg>. This <mu> parameter is for an exponential
 *            distribution anchored from $P=1.0$, so it's not really a
 *            tail per se; but it's only an accurate fit in the tail
 *            of the Forward score distribution, from about $P=0.001$
 *            or so.
 *            
 *            The determination of <mu> is done by a brief simulation
 *            in which we fit a Gumbel distribution to a small number
 *            of Forward scores of random sequences, and use that to
 *            predict the location of the tail at probability <tailp>.
 *            
 *            The Gumbel is of course inaccurate, but we can use it
 *            here solely as an empirical distribution to determine
 *            the location of a reasonable <mu> more accurately on a
 *            smaller number of samples than we could do with raw
 *            order statistics. 
 *            
 *            Typical choices are L=100, N=200, tailp=0.04, which
 *            typically yield estimates $\hat{\mu}$ with a precision
 *            (standard deviation) of $\pm$ 0.2 bits, corresponding to
 *            a $\pm$ 15\% error in E-values. See [J1/135].
 *            
 *            The use of Gumbel fitting to a small number of $N$
 *            samples and the extrapolation of $\hat{\mu}$ from the
 *            estimated location of the 0.04 tail mass are both
 *            empirical and carefully optimized against several
 *            tradeoffs. Most importantly, around this choice of tail
 *            probability, a systematic error introduced by the use of
 *            the Gumbel fit is being cancelled by systematic error
 *            introduced by the use of a higher tail probability than
 *            the regime in which the exponential tail is a valid
 *            approximation. See [J1/135] for discussion.
 *            
 *            This function changes the length configuration of both
 *            <gm> and <bg>. The caller must remember to reconfigure
 *            both of their length models appropriately for any
 *            subsequent alignments.
 *            
 * Args:      r      : source of randomness
 *            gm     : configured profile to sample sequences from
 *            bg     : null model (for background residue frequencies)
 *            L      : mean length model for seq emission from profile
 *            N      : number of sequences to generate
 *            lambda : expected slope of the exponential tail (from p7_Lambda())
 *            tailp  : tail mass from which we will extrapolate mu
 *            ret_mu : RETURN: estimate for the Forward mu (base of exponential tail)
 *
 * Returns:   <eslOK> on success, and <*ret_fv> is the score difference
 *            in bits.
 *
 * Throws:    <eslEMEM> on allocation error, and <*ret_fv> is 0.
 */
int
p7_Tau(ESL_RANDOMNESS *r, P7_PROFILE *gm, P7_BG *bg, int L, int N, double lambda, double tailp, double *ret_tau)
{
  P7_GMX  *gx      = p7_gmx_Create(gm->M, L);	     /* DP matrix: L (rows) will grow as needed */
  ESL_DSQ *dsq     = NULL;
  double  *xv      = NULL;
  float    fsc, nullsc;		                  
  double   gmu, glam;
  int      status;
  int      i;

  ESL_ALLOC(xv,  sizeof(double)  * N);
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if (gx == NULL) { status = eslEMEM; goto ERROR; }

  p7_ReconfigLength(gm, L);
  p7_bg_SetLength(bg, L);

  for (i = 0; i < N; i++)
    {
      if ((status = esl_rsq_xfIID(r, bg->f, gm->abc->K, L, dsq)) != eslOK) goto ERROR;
      if ((status = p7_GForward(dsq, L, gm, gx, &fsc))           != eslOK) goto ERROR;
      if ((status = p7_bg_NullOne(bg, dsq, L, &nullsc))          != eslOK) goto ERROR;   
      xv[i] = (fsc - nullsc) / eslCONST_LOG2;
    }
  if ((status = esl_gumbel_FitComplete(xv, N, &gmu, &glam)) != eslOK) goto ERROR;

  /* Explanation of the eqn below: first find the x at which the Gumbel tail
   * mass is predicted to be equal to tailp. Then back up from that x
   * by log(tailp)/lambda to set the origin of the exponential tail to 1.0
   * instead of tailp.
   */
  *ret_tau =  esl_gumbel_invcdf(1.0-tailp, gmu, glam) + (log(tailp) / lambda);
  
  free(xv);
  free(dsq);
  p7_gmx_Destroy(gx);
  return eslOK;

 ERROR:
  *ret_tau = 0.;
  if (xv  != NULL) free(xv);
  if (dsq != NULL) free(dsq);
  if (gx  != NULL) p7_gmx_Destroy(gx);
  return status;
}


/*****************************************************************
 * 2. Experiment code.
 *****************************************************************/

#ifdef p7EVALUES_STATS
/* gcc -o stats -g -O2 -I. -L. -I../easel -L../easel -Dp7EVALUES_STATS evalues.c -lhmmer -leasel -lm
 * ./stats <hmmfile>
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_alphabet.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     NULL, NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed to <n>",        0 },
  { "-t",        eslARG_REAL,  "0.05", NULL, NULL,  NULL,  NULL, NULL, "set fitted tail probability to <x>",             0 },
  { "-L",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "set length to <n>",                              0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "set seq number to <n>",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "collect test statistics for E-value calculations";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  double          lambda;
  double          vmu;
  double          fmu;
  double          tailp = esl_opt_GetReal(go, "-t");
  int             L     = esl_opt_GetInteger(go, "-L");
  int             N     = esl_opt_GetInteger(go, "-N");
  int             status;

  if (esl_opt_IsDefault(go, "-s"))  r = esl_randomness_CreateTimeseeded();
  else                              r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslEOF) 
    {
      if (bg == NULL) bg = p7_bg_Create(abc);
      gm = p7_profile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);

      p7_Lambda(hmm, bg, &lambda);
      p7_Mu(r, gm, bg, L, N, lambda, &vmu);
      p7_Tau(r, gm, bg, L, N, tailp,  &fmu);
      
      printf("%s %.4f %.4f %.4f\n", hmm->name, lambda, vmu, fmu);

      p7_hmm_Destroy(hmm);      
      p7_profile_Destroy(gm);
    }

  p7_hmmfile_Close(hfp);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7EVALUES_STATS*/


#ifdef p7EXP_J1_135
/* Determining precision and accuracy of mu estimates. [xref J1/135]
 * gcc -o exp -g -O2 -I. -L. -I../easel -L../easel -Dp7EXP_J1_135 evalues.c -lhmmer -leasel -lm
 * ./exp <hmmfile>
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_alphabet.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-f",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "fit the Forward mu, not Viterbi mu",             0 },
  { "-s",        eslARG_INT,     NULL, NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed to <n>",        0 },
  { "-l",        eslARG_REAL,"0.6931", NULL, NULL,  NULL,  NULL, NULL, "set lambda param to <x>",                        0 },
  { "-t",        eslARG_REAL,  "0.05", NULL, NULL,  NULL,  NULL, NULL, "set fitted tail probability to <x>",             0 },
  { "-L",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "set length to <n>",                              0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "set seq number to <n>",                          0 },
  { "-Z",        eslARG_INT,   "1000", NULL, NULL,  NULL,  NULL, NULL, "set number of iterations to <n>",                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "collect test statistics for E-value calculations";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  double          lambda  = esl_opt_GetReal   (go, "-l");
  double          tailp   = esl_opt_GetReal   (go, "-t");
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  int             Z       = esl_opt_GetInteger(go, "-Z");
  double          mu;
  int             i;


  if (esl_opt_IsDefault(go, "-s"))  r = esl_randomness_CreateTimeseeded();
  else                              r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM from %s", hmmfile);
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);

  for (i = 0; i < Z; i++)
    {
      if (esl_opt_GetBoolean(go, "-f")) p7_Tau(r, gm, bg, L, N, tailp,  &mu);
      else                              p7_Mu(r, gm, bg, L, N, lambda, &mu);
      printf("%.4f\n", mu);
    }

  p7_hmm_Destroy(hmm);      
  p7_profile_Destroy(gm);
  p7_hmmfile_Close(hfp);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7EXP_J1_135*/


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
