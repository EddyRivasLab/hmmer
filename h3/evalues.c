/* Calculations and simulations relevant to E-value calculations.
 * 
 * SRE, Mon Aug  6 13:00:06 2007
 * SVN $Id$
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_gumbel.h"
#include "esl_random.h"

#include "hmmer.h"

/* Function:  p7_Lambda()
 * Synopsis:  Determines length-corrected local lambda parameter.
 * Incept:    SRE, Wed Aug  8 17:54:55 2007 [Janelia]
 *
 * Purpose:   Determine the effective scale parameter $\hat{\lambda}$ 
 *            for model <hmm>, applying both to Gumbel distributions and
 *            exponential tails.
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
 *            correction is based largely on \citep{Altschul01}.
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
 * Synopsis:  Determines the local Gumbel mu parameter for a model.
 * Incept:    SRE, Mon Aug  6 13:00:57 2007 [Janelia]
 *
 * Purpose:   Given model <gm> configured for local alignment (typically
 *            multihit, but may be unihit), determine the Gumbel location
 *            parameter $\mu$ by brief simulation. The simulation
 *            generates <N> random sequences of
 *            length <L> using background frequencies in the null
 *            model <bg> and the random number generator <r>; scores
 *            them with <gm> and <bg> with the Viterbi algorithm; and
 *            fitting the resulting distribution to a Gumbel of known
 *            <lambda>.
 *            
 *            ... TODO: insert typical N,L and accuracy here...
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
      if ((status = esl_rnd_xfIID(r, bg->f, gm->abc_r->K, L, dsq)) != eslOK) goto ERROR;
      if ((status = p7_GViterbi(dsq, L, gm, gx, &sc))              != eslOK) goto ERROR;
      if ((status = p7_bg_NullOne(bg, dsq, L, &nullsc))            != eslOK) goto ERROR;   
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


/* Function:  p7_FVOffset()
 * Synopsis:  Determine Forward-Viterbi score tail offset by brief simulation.
 * Incept:    SRE, Wed Aug  8 11:41:39 2007 [Janelia]
 *
 * Purpose:   Determine the offset of the high scoring tails of Viterbi
 *            and Forward score distributions by a brief
 *            simulation. Sequences are generated from the profile
 *            <gm>, which the caller should have already configured in
 *            unilocal or multilocal mode, accompanied by core model
 *            probabilities in <hmm> and background model
 *            probabilities in <bg>. The profile's length model 
 *            is set to mean length <L> for the purposes of generating
 *            the simulated sequences, and <N> sequences are sampled,
 *            using <r> as the source of randomness.
 *            
 * Args:      r      : source of randomness
 *            gm     : configured profile to sample sequences from
 *            hmm    : core model probabilities for profile
 *            bg     : null model (for background residue frequencies)
 *            L      : mean length model for seq emission from profile
 *            N      : number of sequences to generate
 *            ret_fv : RETURN: mean Forward-Viterbi score difference (in bits)
 *
 * Returns:   <eslOK> on success, and <*ret_fv> is the score difference
 *            in bits.
 *
 * Throws:    <eslEMEM> on allocation error, and <*ret_fv> is 0.
 */
#if 0
int
p7_FVOffset(ESL_RANDOMNESS *r, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, int L, int N, double *ret_fv)
{
  P7_GMX  *gx      = p7_gmx_Create(gm->M, L*2);	     /* DP matrix: L (rows) will grow as needed */
  ESL_SQ  *sq      = esl_sq_CreateDigital(hmm->abc); /* sampled digital sequence object: will grow as needed */
  float    fsc, vsc;		                     /* Viterbi, Forward scores         */
  double   sum = 0.;
  int      status;
  int      i;
  
  if (sq == NULL || gx == NULL) { status = eslEMEM; goto ERROR; }

  for (i = 0; i < N; i++)
    {
      if ((status = p7_ReconfigLength(gm, L))                  != eslOK) goto ERROR; /* configure for generating about L        */
      if ((status = p7_ProfileEmit(r, hmm, gm, bg, sq, NULL))  != eslOK) goto ERROR; /* NULL=unwanted trace. sq->dsq now set.   */
      if ((status = p7_gmx_GrowTo(gx, gm->M, sq->n))           != eslOK) goto ERROR;
      if ((status = p7_ReconfigLength(gm, sq->n))              != eslOK) goto ERROR; /* configure for scoring this seq length n */
      if ((status = p7_GViterbi(sq->dsq, sq->n, gm, gx, &vsc)) != eslOK) goto ERROR;
      if ((status = p7_GForward(sq->dsq, sq->n, gm, gx, &fsc)) != eslOK) goto ERROR;
      sum += fsc-vsc;
      /*printf("# %.4f %.4f  %.4f\n", fsc/eslCONST_LOG2, vsc/eslCONST_LOG2, (fsc-vsc)/eslCONST_LOG2);*/
    }

  *ret_fv = (sum / (double) N) / eslCONST_LOG2;
  esl_sq_Destroy(sq);
  p7_gmx_Destroy(gx);
  return eslOK;

 ERROR:
  *ret_fv = 0.;
  if (sq != NULL) esl_sq_Destroy(sq);
  if (gx != NULL) p7_gmx_Destroy(gx);
  return status;
}
#endif
int
p7_FVOffset(ESL_RANDOMNESS *r, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, int L, int N, double *ret_fv)
{
  P7_GMX  *gx      = p7_gmx_Create(gm->M, L);	  
  ESL_DSQ *dsq     = NULL;
  float    fsc, vsc, nullsc;		                  
  double   sum = 0.;
  int      status;
  int      i;

  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if (gx == NULL) { status = eslEMEM; goto ERROR; }

  p7_ReconfigLength(gm, L);
  p7_bg_SetLength(bg, L);

  for (i = 0; i < N; i++)
    {
      if ((status = esl_rnd_xfIID(r, bg->f, gm->abc_r->K, L, dsq)) != eslOK) goto ERROR;
      if ((status = p7_GViterbi(dsq, L, gm, gx, &vsc))             != eslOK) goto ERROR;
      if ((status = p7_GForward(dsq, L, gm, gx, &fsc))             != eslOK) goto ERROR;
      if ((status = p7_bg_NullOne(bg, dsq, L, &nullsc))            != eslOK) goto ERROR;   


      fsc  = (fsc - nullsc) / eslCONST_LOG2;
      vsc  = (vsc - nullsc) / eslCONST_LOG2;
      sum += fsc-vsc;
      printf("# %.4f %.4f %.4f\n", fsc, vsc, fsc-vsc);
    }

  *ret_fv = (sum / (double) N);
  free(dsq);
  p7_gmx_Destroy(gx);
  return eslOK;

 ERROR:
  *ret_fv = 0.;
  if (dsq != NULL) free(dsq);
  if (gx  != NULL) p7_gmx_Destroy(gx);
  return status;
}


/*****************************************************************
 * 
 *****************************************************************/

#ifdef p7EVALUES_STATS
/* mpicc -o stats -g -O2 -I. -L. -I../easel -L../easel -Dp7EVALUES_STATS evalues.c -lhmmer -leasel -lm
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
  double          mu;
  double          fv_offset;
  int             status;
  int             L = esl_opt_GetInteger(go, "-L");
  int             N = esl_opt_GetInteger(go, "-N");

  if (esl_opt_IsDefault(go, "-s"))  r = esl_randomness_CreateTimeseeded();
  else                              r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslEOF) 
    {
      if (bg == NULL) bg = p7_bg_Create(abc);
      gm = p7_profile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);

      p7_Lambda(hmm, bg, &lambda);
      p7_Mu(r, gm, bg, L, N, lambda, &mu);
      p7_FVOffset(r, gm, hmm, bg, L, N, &fv_offset);
      
      printf("%s %.4f %.4f %.4f\n", hmm->name, lambda, mu, fv_offset);

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
