/* Various summary statistics calculated for HMMs and profiles.
 * 
 * SRE, Fri May  4 11:43:20 2007 [Janelia]
 * SVN $Id$
 */

#include "p7_config.h"
#include <easel.h>
#include <esl_vectorops.h>
#include "hmmer.h"



/* Function:  p7_MeanMatchInfo()
 * Incept:    SRE, Fri May  4 11:43:56 2007 [Janelia]
 *
 * Purpose:   Calculate the mean information content per match state
 *            emission distribution, in bits:
 *            
 *            \[
 *              \frac{1}{M} \sum_{k=1}^{M} \left[ - \sum_x p_k(x) \log_2 p_k(x) + \sum_x f(x) \log_2 f(x) \right] 
 *            \]
 *            
 *            where $p_k(x)$ is emission probability for symbol $x$ from match state $k$,
 *            and $f(x)$ is the null model's background emission probability for $x$. 
 *            
 *            This statistic is used in "entropy weighting" to set the
 *            total sequence weight when model building.
 */
double
p7_MeanMatchInfo(const P7_HMM *hmm, const P7_BG *bg)
{
  return esl_vec_FEntropy(bg->f, hmm->abc->K) - p7_MeanMatchEntropy(hmm);
}

/* Function:  p7_MeanMatchEntropy()
 * Incept:    SRE, Fri May  4 13:37:15 2007 [Janelia]
 *
 * Purpose:   Calculate the mean entropy per match state emission
 *            distribution, in bits:
 *            
 *            \[
 *              \frac{1}{M} \sum_{k=1}^{M} -\sum_x p_k(x) \log_2 p_k(x)
 *            \]
 *       
 *            where $p_k(x)$ is emission probability for symbol $x$
 *            from match state $k$.
 */
double
p7_MeanMatchEntropy(const P7_HMM *hmm)
{
  int    k;
  double H = 0.;

  for (k = 1; k <= hmm->M; k++)
    H += esl_vec_FEntropy(hmm->mat[k], hmm->abc->K);
  H /= (double) hmm->M;
  return H;
}


/* Function:  p7_MeanMatchRelativeEntropy()
 * Incept:    SRE, Fri May 11 09:25:01 2007 [Janelia]
 *
 * Purpose:   Calculate the mean relative entropy per match state emission
 *            distribution, in bits:
 *            
 *            \[
 *              \frac{1}{M} \sum_{k=1}^{M} \sum_x p_k(x) \log_2 \frac{p_k(x)}{f(x)}
 *            \]
 *       
 *            where $p_k(x)$ is emission probability for symbol $x$
 *            from match state $k$, and $f(x)$ is the null model's 
 *            background emission probability for $x$. 
 */
double
p7_MeanMatchRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg)
{
  int    k;
  double KL = 0.;

#if 0
  p7_bg_Dump(stdout, hmm->bg);
  for (k = 1; k <= hmm->M; k++)
    printf("Match %d : %.2f %.2f\n", k, 
	   esl_vec_FRelEntropy(hmm->mat[k], hmm->bg->f, hmm->abc->K),
	   esl_vec_FEntropy(bg->f, hmm->abc->K) - esl_vec_FEntropy(hmm->mat[k], hmm->abc->K));
#endif

  for (k = 1; k <= hmm->M; k++)
    KL += esl_vec_FRelEntropy(hmm->mat[k], bg->f, hmm->abc->K);
  KL /= (double) hmm->M;
  return KL;
}



double
p7_MeanForwardScore(const P7_HMM *hmm, const P7_BG *bg)
{
  int             L   = 350;
  int             N   = 100;
  P7_PROFILE     *gm  = p7_profile_Create(hmm->M, hmm->abc);
  P7_GMX         *gx  = p7_gmx_Create(gm->M, L);
  ESL_SQ         *sq  = esl_sq_CreateDigital(hmm->abc);
  ESL_RANDOMNESS *r   = esl_randomness_CreateTimeseeded();
  int             fsc;
  int             nullsc;
  double          bitscore;
  double          sum = 0.;
  int             i;

  if (p7_ProfileConfig (hmm, bg, gm, p7_LOCAL) != eslOK) p7_Die("failed to configure profile");

  for (i = 0; i < N; i++)
    {
      if (p7_ReconfigLength(gm, L)                        != eslOK) p7_Die("failed to reconfig profile length");
      if (p7_ProfileEmit(r, gm, sq, NULL)                 != eslOK) p7_Die("failed to emit sequence");

      if (p7_ReconfigLength(gm, sq->n)                    != eslOK) p7_Die("failed to reconfig profile length");
      if (p7_gmx_GrowTo(gx, gm->M, sq->n)                 != eslOK) p7_Die("failed to grow the matrix");
      if (p7_GForward(sq->dsq, sq->n, gm, gx, &fsc)       != eslOK) p7_Die("failed to run Forward");
      if (p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc)      != eslOK) p7_Die("failed to run bg_NullOne()");
      bitscore = p7_SILO2Bitscore(fsc - nullsc);

      sum += bitscore;
    }
  
  esl_randomness_Destroy(r);
  esl_sq_Destroy(sq);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  return (sum / (double) N);
}


int
p7_MeanPositionRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg, double *ret_entropy)
{
  int     status;
  float  *mocc = NULL;
  int     k;
  double  mre, tre;
  double  xm, xi, xd;
  
  ESL_ALLOC(mocc, sizeof(float) * (hmm->M+1));
  if ((status = p7_hmm_CalculateOccupancy(hmm, mocc)) != eslOK) goto ERROR;
  
  /* The weighted relative entropy per match emission */
  for (mre = 0., k = 1; k <= hmm->M; k++)
    mre += mocc[k] * esl_vec_FRelEntropy(hmm->mat[k], bg->f, hmm->abc->K);
  mre /= esl_vec_FSum(mocc+1, hmm->M);

  /* The weighted relative entropy per match entry transition, 2..M 
   */
  for (tre = 0., k = 2; k <= hmm->M; k++)
    {
      xm = mocc[k-1]*hmm->t[k-1][p7_TMM] * log(hmm->t[k-1][p7_TMM] / bg->p1);
      xi = mocc[k-1]*hmm->t[k-1][p7_TMI] * (log(hmm->t[k-1][p7_TMM] / bg->p1) + log(hmm->t[k-1][p7_TIM] / bg->p1));
      xd = (1.-mocc[k-1])*hmm->t[k-1][p7_TDM] * log(hmm->t[k-1][p7_TDM] / bg->p1);
      tre += (xm+xi+xd) / eslCONST_LOG2;
    }
  tre /= esl_vec_FSum(mocc+2, hmm->M-1);

  free(mocc);
  *ret_entropy = mre+tre;
  return eslOK;

 ERROR:
  if (mocc != NULL) free(mocc);
  *ret_entropy = 0.;
  return status;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
