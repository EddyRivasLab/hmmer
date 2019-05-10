/* modelsample.c : sampling randomized profile HMMs
 * 
 * These routines are used primarily in model testing. Models are
 * often contrived and/or constrained in various ways that enable
 * particular tests.
 * 
 * Contents:
 *    1. Model sampling routines.
 *    2. Internal support functions.
 *    
 */
#include "h4_config.h"

#include "esl_dirichlet.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "h4_profile.h"

#include "modelsample.h"

static void zeropeppered_probvector(ESL_RANDOMNESS *rng, float *p, int n);

/*****************************************************************
 * 1. Model sampling
 *****************************************************************/

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
      esl_dirichlet_FSampleUniform(rng, 3, hmm->t[k]+3);
      esl_dirichlet_FSampleUniform(rng, 3, hmm->t[k]+6);
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

/* Function:  h4_modelsample_zeropeppered()
 * Synopsis:  Model, with occasional zeros. 
 * Incept:    SRE, Fri 10 May 2019 [Mike & The Mechanics, Silent Running]
 */
int
h4_modelsample_zeropeppered(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm)
{
  H4_PROFILE *hmm     = NULL;
  float       tII_max = 0.99;
  int         roll;
  int         k;
  int         status;

  if (( hmm = h4_profile_Create(abc, M)) == NULL) { status = eslEMEM; goto ERROR; }

  /* transitions */
  roll = esl_rnd_Roll(rng, 10);
  if      (roll == 0) hmm->t[0][h4_TMM] = 1.;               // 10%
  else if (roll == 1) hmm->t[0][h4_TMM] = 0.;               // 10%
  else                hmm->t[0][h4_TMM] = esl_random(rng);  // 80%
  hmm->t[0][h4_TMD] = 1. - hmm->t[0][h4_TMM];
  for (k = 1; k < M; k++)
    do {
      zeropeppered_probvector(rng, hmm->t[k],   3);
      zeropeppered_probvector(rng, hmm->t[k]+3, 3);
      zeropeppered_probvector(rng, hmm->t[k]+6, 3);
    } while (hmm->t[k][h4_TII] > tII_max);                 // tII ~ 1.0 is too evil: watch out for infinite length seqs

  /* match emissions */
  for (k = 1; k <= M; k++)
    zeropeppered_probvector(rng, hmm->e[k], abc->K);

  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  h4_profile_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}

/*****************************************************************
 * 2. internal support functions
 *****************************************************************/

static void
zeropeppered_probvector(ESL_RANDOMNESS *rng, float *p, int n)
{
  esl_dirichlet_FSampleUniform(rng, n, p);
  if (esl_rnd_Roll(rng, 2))	/* 50% of the time, we throw a zero into the sampled p's */
    {
      p[esl_rnd_Roll(rng, n)] = 0.0;
      esl_vec_FNorm(p, n);
    }
}
