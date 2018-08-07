/* modelsample.c : sampling randomized profile HMMs
 * 
 * These routines are used primarily in model testing. Models are
 * often contrived and/or constrained in various ways that enable
 * particular tests.
 * 
 * Contents:
 *    1. Model sampling routines.
 *    
 */
#include "h4_config.h"

#include "esl_dirichlet.h"

#include "h4_profile.h"

#include "modelsample.h"

/* Function:  h4_modelsample()
 * Synopsis:  Sample a random profile HMM.
 * Incept:    SRE, Mon 06 Aug 2018
 */
int
h4_modelsample(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm)
{
  H4_PROFILE *hmm     = NULL;
  float       tII_max = 0.99;
  int         k;
  int         status;

  if (( hmm = h4_profile_Create(abc, M)) == NULL) { status = eslEMEM; goto ERROR; }

  /* transitions */
  hmm->t[0][h4_TMM] = esl_random(rng);
  hmm->t[0][h4_TMD] = 1. - hmm->t[0][h4_TMM];
  for (k = 1; k < M; k++)
    do {
      esl_dirichlet_FSampleUniform(rng, 3, hmm->t[k]);
      esl_dirichlet_FSampleUniform(rng, 3, hmm->t[k]+6);
      esl_dirichlet_FSampleUniform(rng, 3, hmm->t[k]+3);
    } while (hmm->t[k][h4_TII] > tII_max);                 // tII ~ 1.0 is too evil: watch out for infinite length seqs

  /* match emissions */
  for (k = 1; k <= M; k++)
    esl_dirichlet_FSampleUniform(rng, abc->K, hmm->e[k]);

  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  h4_profile_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}
