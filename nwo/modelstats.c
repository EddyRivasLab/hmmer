#include "h4_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "h4_profile.h"

/* Function:  h4_MeanMatchKL()
 * Synopsis:  Calculate and return mean relative entropy per match emission
 * Incept:    SRE, Wed 03 Apr 2019
 *
 * Purpose:   Calculate and return the mean expected score per match
 *            emission distribution, in bits:
 *
 *            \[
 *              \frac{1}{M} \sum_{k=1}^{M} \sum_x e_k(x) \log_2 \frac{e_k(x)}{f(x)}
 *            \]
 *
 *            This is also the mean relative entropy, or the mean Kullback-Leibler
 *            divergence, of the match emission distribution relative to background
 *            residue frequencies.
 *            
 * Returns:   mean score, in bits.
 *
 * Xref:      eweight.c : used as optimization target for entropy weighting
 */
float
h4_MeanMatchKL(H4_PROFILE *hmm)
{
  float  avgKL = 0.;
  int    k;

  for (k = 1; k <= hmm->M; k++)
    avgKL += esl_vec_FRelEntropy(hmm->e[k], hmm->f, hmm->abc->K); // in bits.
  return avgKL / (float) hmm->M;
}
  
  
  
