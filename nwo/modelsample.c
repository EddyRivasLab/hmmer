/* modelsample.c : sampling randomized profile HMMs
 * 
 * These routines are used primarily in model testing. Models are
 * often contrived and/or constrained in various ways that enable
 * particular tests.
 * 
 * Contents:
 *    1. Model sampling routines.
 *    2. Internal support functions.
 *    3. Unit tests.
 *    4. Test driver.
 */
#include "h4_config.h"

#include "esl_dirichlet.h"
#include "esl_matrixops.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "h4_profile.h"

#include "standardize.h"
#include "vectorize.h"

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
  float       mtot    = 0.0;
  int         k;
  int         status;

  if (( hmm = h4_profile_Create(abc, M)) == NULL) { status = eslEMEM; goto ERROR; }

  /* transitions */
  while (mtot < 0.1)                                           // reject uber-pathological models that essentially never visit a match state
    {
      hmm->t[0][h4_TMM] = esl_random(rng);
      hmm->t[0][h4_TMD] = 1. - hmm->t[0][h4_TMM];
      for (k = 1; k < M; k++)
	do {
	  esl_dirichlet_FSampleUniform(rng, 3, hmm->t[k]);
	  esl_dirichlet_FSampleUniform(rng, 3, hmm->t[k]+3);
	  esl_dirichlet_FSampleUniform(rng, 3, hmm->t[k]+6);
	} while (hmm->t[k][h4_TII] > tII_max);                 // tII ~ 1.0 is uber-pathological too: gives infinite length seqs

      h4_profile_Occupancy(hmm, NULL, NULL, &mtot, NULL);      // only needs transition probs, so don't waste time setting emissions 'til we're good here
    }

  /* match emissions */
  for (k = 1; k <= M; k++)
    esl_dirichlet_FSampleUniform(rng, abc->K, hmm->e[k]);

  hmm->flags |= h4_HASPROBS;
  if (( status = h4_standardize(hmm)) != eslOK) goto ERROR;
  if (( status = h4_vectorize(hmm))   != eslOK) goto ERROR;
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
  float       mtot    = 0.0;
  int         roll;
  int         k;
  int         status;

  if (( hmm = h4_profile_Create(abc, M)) == NULL) { status = eslEMEM; goto ERROR; }

  /* transitions */
  while (mtot < 0.1)
    {
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

      h4_profile_Occupancy(hmm, NULL, NULL, &mtot, NULL);
    }

  /* match emissions */
  for (k = 1; k <= M; k++)
    zeropeppered_probvector(rng, hmm->e[k], abc->K);

  hmm->flags |= h4_HASPROBS;
  if (( status = h4_standardize(hmm)) != eslOK) goto ERROR;
  if (( status = h4_vectorize(hmm))   != eslOK) goto ERROR;
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  h4_profile_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}


/* Function:  h4_modelsample_Enumerable()
 * Synopsis:  Sample random HMM with no nonzero transitions to inserts.
 * Incept:    SRE, Sat 18 May 2019 [Nick Cave, People Ain't No Good]
 *
 * Purpose:   Sample a random HMM with random emission and transition
 *            probabilities with the exception that all transitions to
 *            insert are zero. This makes it possible to create a
 *            model with a finite, easily enumerable sequence space
 *            (all seqs of length 1..M).
 *            
 *            To achieve finite enumerability, the caller must also
 *            configure a unihit alignment mode with a target length
 *            of 0.
 *
 *            Useful for debugging and validating DP algorithms.
 *
 *            Compare <h4_modelsample_Enumerable2()>, which only makes
 *            tII transitions zero, and thus enumerates a larger space
 *            consisting of sequences of length 1..2M-1.
 *            
 * Returns:   <eslOK> on success. The newly allocated model is returned through
 *            <ret_hmm>. The caller is responsible for freeing it.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
h4_modelsample_Enumerable(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm)
{
  H4_PROFILE *hmm = NULL;
  float       mtot = 0.0;
  float       tmp[2];
  int         k;
  int         status;

  ESL_DASSERT1(( M > 0 ));

  if (( hmm = h4_profile_Create(abc, M)) == NULL) { status = eslEMEM; goto ERROR; }

  while (mtot < 0.1)
    {
      esl_dirichlet_FSampleUniform(rng, 2, tmp);
      hmm->t[0][h4_TMM] = tmp[0];
      hmm->t[0][h4_TMD] = tmp[1];

      for (k = 1; k < M; k++)
	{
	  esl_dirichlet_FSampleUniform(rng, 2, tmp); hmm->t[k][h4_TMM] = tmp[0]; hmm->t[k][h4_TMI] = 0.0f; hmm->t[k][h4_TMD] = tmp[1];
	  esl_dirichlet_FSampleUniform(rng, 2, tmp); hmm->t[k][h4_TIM] = tmp[0]; hmm->t[k][h4_TII] = 0.0f; hmm->t[k][h4_TID] = tmp[1]; // I transitions irrelevant, since I's are unreached.
	  esl_dirichlet_FSampleUniform(rng, 2, tmp); hmm->t[k][h4_TDM] = tmp[0]; hmm->t[k][h4_TDI] = 0.0f; hmm->t[k][h4_TDD] = tmp[1]; 
	}

      h4_profile_Occupancy(hmm, NULL, NULL, &mtot, NULL);
    }
      
  for (k = 1; k <= M; k++)
    esl_dirichlet_FSampleUniform(rng, abc->K, hmm->e[k]);    

  hmm->flags |= h4_HASPROBS;
  if (( status = h4_standardize(hmm)) != eslOK) goto ERROR;
  if (( status = h4_vectorize(hmm))   != eslOK) goto ERROR;
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  h4_profile_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}


/* Function:  h4_modelsample_Enumerable2()
 * Synopsis:  Sample an random HMM with no nonzero insert-insert transitions.
 * Incept:    SRE, Sat 18 May 2019
 *
 * Purpose:   Sample a random HMM with random emission and transition
 *            probabilities with the exception that all insert-insert
 *            transitions are zero. This makes it possible to create a
 *            model with a finite, easily enumerable sequence space
 *            (all seqs of length M=1..2M-1.)
 *            
 *            To achieve enumerability, the caller must also configure
 *            a unihit alignment mode with a target length of 0.
 *
 *            Useful for debugging and validating DP algorithms.
 *
 *            Compare <p7_modelsample_Enumerable()>, which makes all
 *            transitions to insert 0, and thus enumerates a smaller
 *            space.
 *            
 * Returns:   <eslOK> on success. The newly allocated model is returned through
 *            <ret_hmm>. The caller is responsible for freeing it.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
h4_modelsample_Enumerable2(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm)
{
  H4_PROFILE *hmm  = NULL;
  float       mtot = 0.0;
  float       tmp[2];
  int         k;
  int         status;

  ESL_DASSERT1(( M > 0 ));

  if (( hmm = h4_profile_Create(abc, M)) == NULL) { status = eslEMEM; goto ERROR; }

  while (mtot < 0.1) 
    {
      esl_dirichlet_FSampleUniform(rng, 2, tmp);
      hmm->t[0][h4_TMM] = tmp[0];
      hmm->t[0][h4_TMD] = tmp[1];

      for (k = 1; k < M; k++)
	{
	  esl_dirichlet_FSampleUniform(rng, 3, hmm->t[k]);
	  esl_dirichlet_FSampleUniform(rng, 2, tmp);
	  hmm->t[k][h4_TIM] = tmp[0];
	  hmm->t[k][h4_TII] = 0.0f;
	  hmm->t[k][h4_TID] = tmp[1];
	  esl_dirichlet_FSampleUniform(rng, 3, hmm->t[k]+6);
	}

      h4_profile_Occupancy(hmm, NULL, NULL, &mtot, NULL);
    }

  for (k = 1; k <= M; k++)
    esl_dirichlet_FSampleUniform(rng, abc->K, hmm->e[k]);

  hmm->flags |= h4_HASPROBS;
  if (( status = h4_standardize(hmm)) != eslOK) goto ERROR;
  if (( status = h4_vectorize(hmm))   != eslOK) goto ERROR;
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  h4_profile_Destroy(hmm);
  *ret_hmm = NULL;
  return status;

}


/* Function:  h4_modelsample_SinglePath()
 * Synopsis:  Sample a random profile HMM with only one P=1.0 path possible.
 * Incept:    SRE, Sun 19 May 2019 [Nick Cave and Warren Ellis]
 *
 * Purpose:   Sample a random profile HMM of length <M> in alphabet
 *            <abc> using random number generator <rng>, such that the
 *            profile has only a single possible path and sequence
 *            with probability 1; return it in <*ret_hmm>. All
 *            transition and emission distributions have one p=1.0
 *            component, with others p=0.0. (All II transitions have
 *            to be 0.0.) This creates a useful unambiguous case for
 *            debugging and testing.
 *
 *            The caller must also configure the alignment mode to
 *            unihit glocal and a target length of zero.
 */
int
h4_modelsample_SinglePath(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, H4_PROFILE **ret_hmm)
{
  H4_PROFILE *hmm  = NULL;
  float       mtot = 0.0;   // expected number of match states occupied
  int         k;
  int         status;

  if (( hmm = h4_profile_Create(abc, M)) == NULL) { status = eslEMEM; goto ERROR; }

  while (mtot < 0.1)
    {
      /* Initialize probability parameters to 0, then set edge conventions */
      esl_mat_FSet(hmm->t, M+1, h4_NT,   0.0f);
      esl_mat_FSet(hmm->e, M+1, abc->K,  0.0f);
      h4_profile_SetConventions(hmm);

      hmm->t[0][esl_rnd_Roll(rng, 2) ? h4_TMM : h4_TMD] = 1.0f;  // G->{M1|D1}

      for (k = 1; k < M; k++)
	{
	  hmm->t[k][esl_rnd_Roll(rng, 3)]                   = 1.0f;
	  hmm->t[k][esl_rnd_Roll(rng, 2) ? h4_TIM : h4_TID] = 1.0f;  // TII must be 0 else inf loop
	  hmm->t[k][esl_rnd_Roll(rng, 3) + 6]               = 1.0f;  // +6 = start of TDM|TDI|TDD
	}

      h4_profile_Occupancy(hmm, NULL, NULL, &mtot, NULL);
    }

  for (k = 1; k <= M; k++)
    hmm->e[k][esl_rnd_Roll(rng, abc->K)] = 1.0f;

  hmm->flags |= h4_HASPROBS;
  if (( status = h4_standardize(hmm)) != eslOK) goto ERROR;
  if (( status = h4_vectorize(hmm))   != eslOK) goto ERROR;
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  h4_profile_Destroy(hmm);
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


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef h4MODELSAMPLE_TESTDRIVE

#include "esl_sq.h"

#include "emit.h"
#include "h4_mode.h"
#include "h4_path.h"

/* utest_sanity()
 * 
 * Basic sanity checks: run the routines, make sure they don't crash, and that they
 * generate valid profile HMMs.
 */
static void
utest_sanity(ESL_RANDOMNESS *rng)
{
  char          msg[] = "modelsample sanity unit test failed";
  H4_PROFILE   *hmm   = NULL;
  ESL_ALPHABET *abc   = esl_alphabet_Create(eslCOINS);
  int           M     = 1 + esl_rnd_Roll(rng, 10);  // 1..10
  char          errbuf[eslERRBUFSIZE];
  
  if ( h4_modelsample(rng, abc, M, &hmm)              != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(hmm, errbuf)               != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  h4_profile_Destroy(hmm);

  if ( h4_modelsample_zeropeppered(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(hmm, errbuf)               != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  h4_profile_Destroy(hmm);

  if ( h4_modelsample_Enumerable(rng, abc, M, &hmm)   != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(hmm, errbuf)               != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  h4_profile_Destroy(hmm);

  if ( h4_modelsample_Enumerable2(rng, abc, M, &hmm)  != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(hmm, errbuf)               != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  h4_profile_Destroy(hmm);

  if ( h4_modelsample_SinglePath(rng, abc, M, &hmm)   != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(hmm, errbuf)               != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  h4_profile_Destroy(hmm);

  esl_alphabet_Destroy(abc);
}


/* utest_enumerated_length()
 * 
 * The "enumerable" profiles generate sequences only up to a fixed
 * length.  Test that that's true. <h4_modelsample_Enumerable()>
 * generates lengths 1..M. <h4_modelsample_Enumerable2()> generates
 * lengths 1..2M-1.
 */
static void
utest_enumerated_length(ESL_RANDOMNESS *rng)
{
  char          msg[] = "modelsample enumerated_length unit test failed";
  H4_PROFILE   *hmm   = NULL;
  H4_MODE      *mo    = h4_mode_Create();
  H4_PATH      *pi    = h4_path_Create();
  ESL_ALPHABET *abc   = esl_alphabet_Create(eslDICE);
  ESL_SQ       *sq    = esl_sq_CreateDigital(abc);
  int           M     = 1 + esl_rnd_Roll(rng, 10);
  int           nseq  = 100;
  int           i, roll;

  /* Enumerable profiles require unihit L=0 alignment mode */
  roll = esl_rnd_Roll(rng, 3);
  if      (roll == 0) h4_mode_SetUnihit(mo);
  else if (roll == 1) h4_mode_SetUnilocal(mo);
  else                h4_mode_SetUniglocal(mo);
  h4_mode_SetLength(mo, 0);

  if ( h4_modelsample_Enumerable(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  for (i = 0; i < nseq; i++)
    {
      if ( h4_emit(rng, hmm, mo, sq, pi) != eslOK) esl_fatal(msg);
      if (sq->n < 1 || sq->n > M)                  esl_fatal(msg);

      h4_path_Reuse(pi);
      esl_sq_Reuse(sq);
    }
  h4_profile_Destroy(hmm);

  if ( h4_modelsample_Enumerable2(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  for (i = 0; i < nseq; i++)
    {
      if ( h4_emit(rng, hmm, mo, sq, pi) != eslOK) esl_fatal(msg);
      if (sq->n < 1 || sq->n > 2*M-1)              esl_fatal(msg);

      h4_path_Reuse(pi);
      esl_sq_Reuse(sq);
    }
  h4_profile_Destroy(hmm);

  h4_path_Destroy(pi);
  h4_mode_Destroy(mo);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
}
#endif // h4MODELSAMPLE_TESTDRIVE

  

/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef h4MODELSAMPLE_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                          docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help summary",               0 },
  { "-s",         eslARG_INT,     "0", NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed",      0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv, "test driver for modelsample", "[-options]");
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng)); 

  utest_sanity(rng);
  utest_enumerated_length(rng);
  
  fprintf(stderr, "#  status   = ok\n");

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}  
#endif // h4MODELSAMPLE_TESTDRIVE
