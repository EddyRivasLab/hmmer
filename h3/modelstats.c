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
p7_MeanMatchInfo(const P7_HMM *hmm)
{
  return esl_vec_FEntropy(hmm->bg->f, hmm->abc->K) - p7_MeanMatchEntropy(hmm);
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
p7_MeanMatchRelativeEntropy(const P7_HMM *hmm)
{
  int    k;
  double KL = 0.;

  for (k = 1; k <= hmm->M; k++)
    KL += esl_vec_FRelEntropy(hmm->mat[k], hmm->bg->f, hmm->abc->K);
  KL /= (double) hmm->M;
  return KL;
}





/*****************************************************************
 * @LICENSE@
 *****************************************************************/
