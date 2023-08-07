/* modelsample.c : sampling randomized profile HMMs
 * 
 * These routines are used primarily in testing. They include various
 * ways to contrive constrained models for particular tests.
 * 
 * Contents:
 *    1. Profile sampling routines
 *    2. Profiles with a readily enumerable number of possible sequences
 *    3. Profiles with only a single possible path
 *    4. Profiles such that ASC paths = all possible paths
 *    5. Internal (static) functions
 *    6. Unit tests
 *    7. Test driver
 *    8. Example
 */
#include <h4_config.h>

#include "esl_dirichlet.h"
#include "esl_matrixops.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_vectorops.h"

#include "h4_anchorset.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"

#include "emit.h"
#include "standardize.h"
#include "vectorize.h"

#include "modelsample.h"

static void zeropeppered_probvector(ESL_RANDOMNESS *rng, float *p, int n);
static int  single_path_seq_engine (ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int do_multi,
                                    H4_PROFILE **ret_hmm, ESL_SQ **ret_sq,
                                    H4_MODE **opt_mo, H4_PATH **opt_pi, H4_ANCHORSET **opt_anch, float *opt_sc);
static int  anchored_ensemble_engine(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M,
                                     H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                                     H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_sc,
                                     int do_uni, int do_local, int do_multi);



/*****************************************************************
 * 1. Profile sampling routines
 *****************************************************************/

/* Function:  h4_modelsample()
 * Synopsis:  Sample a random profile HMM.
 * Incept:    SRE, Mon 06 Aug 2018
 *
 * Purpose:   Construct a profile HMM of <M> consensus positions, for alphabet <abc>,
 *            with randomly sampled parameters (using <rng>).
 *
 *            The profile is allocated here, and returned thru <ret_hmm>. Caller is
 *            responsible for free'ing, with <h4_profile_Destroy()>
 *
 * Args:      rng     - random number generator    [internal state changed]
 *            abc     - alphabet to use
 *            M       - model length in consensus positions
 *            ret_hmm - RETURN: new profile        [allocated here; caller frees]
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
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


/*****************************************************************
 * 2. Profiles with enumerable sequences
 *****************************************************************/


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


/*****************************************************************
 * 3. Profiles with only a single possible path
 *****************************************************************/

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
  *ret_hmm = NULL;
  return status;
}

/* Function:  h4_modelsample_SinglePathSeq()
 * Synopsis:  Create a profile/seq pair with only a single valid path.
 * Incept:    SRE, Mon 11 Jan 2021 [H9/133]
 *
 * Purpose:   Create a profile/sequence pair for which only a single
 *            path is possible. $P(\pi | x, \theta) = 1$, so $P(x, \pi
 *            | \theta) = P(x | \theta)$, so Fwd, Bck, Vit, ASC Fwd,
 *            ASC Bck, and ASC Vit scores are all the same.
 *
 *            Alignment mode must be uniglocal, but target length
 *            model is allowed to be L>0 as usual.
 *
 * Args:      input:
 *            rng      - random number generator
 *            abc      - digital alphabet 
 *            M        - length of model to sample (>= 1)
 *
 *            output:
 *            ret_hmm  - sampled profile \theta
 *            ret_sq   - sampled seq x
 *
 *            optional output:
 *            opt_mo   - comparison mode (uniglocal, with L=sq->n set)
 *            opt_pi   - sampled path \pi
 *            opt_anch - sampled anchorset A (for one domain D=1)
 *            opt_sc   - raw bitscore of the single path
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
h4_modelsample_SinglePathSeq(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, 
                             H4_PROFILE **ret_hmm, ESL_SQ **ret_sq,
                             H4_MODE **opt_mo, H4_PATH **opt_pi, H4_ANCHORSET **opt_anch, float *opt_sc)
{
  return single_path_seq_engine(rng, abc, M, /*do_multi=*/FALSE,
                                ret_hmm, ret_sq,
                                opt_mo, opt_pi, opt_anch, opt_sc);
}

/* Function:  h4_modelsample_SinglePathASC()
 * Synopsis:  Create a profile/seq/anchorset triplet with only a single valid path
 * Incept:    SRE, Mon 11 Jan 2021 [H9/133]
 *
 * Purpose:   Create a profile/sequence/anchorset triplet for which only
 *            a single path is possible. $P(\pi | x, A, \theta) = 1$,
 *            so $P(x, A, \pi | \theta) = P(x, A | \theta)$, so ASC
 *            Fwd, ASC Bck, and ASC Vit scores are all the same.
 *
 *            Alignment mode is multiglocal, and target length model
 *            is the usual L.
 *
 * Args:      input:
 *            rng      - random number generator
 *            abc      - digital alphabet 
 *            M        - length of model to sample (>= 1)
 *
 *            output:
 *            ret_hmm  - sampled profile \theta (new obj created here)
 *            ret_sq   - sampled seq x          (new obj created here)
 *            ret_anch - sampled anchorset A    (new obj created here)
 *
 *            optional output:
 *            opt_mo   - comparison mode (multiglocal, with L=sq->n set)
 *            opt_pi   - sampled path \pi
 *            opt_sc   - raw bitscore of the single path
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
h4_modelsample_SinglePathASC(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, 
                             H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                             H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_sc)
{
  return single_path_seq_engine(rng, abc, M, /*do_multi=*/TRUE,
                                ret_hmm, ret_sq, 
                                opt_mo, opt_pi, ret_anch, opt_sc);
}


/*****************************************************************
 * 4. Profiles such that ASC paths = all possible paths
 *****************************************************************/

/* Function:  h4_modelsample_AnchoredUni()
 * Synopsis:  Create profile/seq/anchorset triplet such that all paths use anchor
 * Incept:    SRE, Sun 31 Jan 2021
 *
 * Purpose:   One of three ways to sample a model, a sequence, and a
 *            single anchor, contrived such that all paths for the
 *            model/sequence comparison must use the anchor even when
 *            unconstrained by an ASC algorithm.
 *
 *            In this version (uni), we have a uniglocal profile, with
 *            an anchor Mk0 state that must emit a particular residue
 *            X with probability 1.0; all paths are forced to use Mk0 not
 *            Dk0; and the target sequence has only one X in it.
 *
 *            Unlike the other versions (local, multi), paths can
 *            contain inserts and N/C states, and the target length
 *            model is L>0.
 *
 * Args:      input:
 *            rng      - random number generator
 *            abc      - digital alphabet 
 *            M        - length of model to sample (>= 1)
 *
 *            output:
 *            ret_hmm  - sampled profile \theta
 *            ret_sq   - sampled seq x
 *            ret_anch - sampled anchorset A 
 *
 *            optional output:
 *            opt_mo   - comparison mode (uniglocal, with L=sq->n set)
 *            opt_pi   - path that <sq> was obtained from
 *            opt_tsc  - score of that path
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 * 
 * Xref:      See anchored_ensemble_engine() for more details.
 */
int
h4_modelsample_AnchoredUni(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, 
                           H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                           H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_tsc)
{
  return anchored_ensemble_engine(rng, abc, M, ret_hmm, ret_sq, ret_anch, opt_mo, opt_pi, opt_tsc, TRUE, FALSE, FALSE);
}


/* Function:  h4_modelsample_AnchoredLocal()
 * Synopsis:  Make profile/seq/anchorset such that all paths use anchor, allowing local paths.
 * Incept:    SRE, Mon 01 Feb 2021
 *
 * Purpose:   Another of three ways to sample a model, a sequence, and an
 *            anchor set such that the ensemble of all valid paths must
 *            use the anchor.
 *
 *            In this version (local), we make a unihit local/glocal
 *            profile that only emits from match states (L=0 and
 *            inserts prohibited), and only the anchor Mk0 state has
 *            e_k0(X)>0; the sequence has only one X in it, at the
 *            anchor.
 *
 *            This version differs from the others (uni, multi) in that
 *            tests local paths.
 *            
 * Args:      input:
 *            rng      - random number generator
 *            abc      - digital alphabet 
 *            M        - length of model to sample (>= 1)
 *
 *            output:
 *            ret_hmm  - sampled profile \theta
 *            ret_sq   - sampled seq x
 *            ret_anch - sampled anchorset A 
 *
 *            optional output:
 *            opt_mo   - comparison mode (unihit, L=0)
 *            opt_pi   - path that <sq> was obtained from
 *            opt_tsc  - score of that path
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      See anchored_ensemble_engine() for more details.
 */
int
h4_modelsample_AnchoredLocal(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M,
                             H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                             H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_tsc)
{
  return anchored_ensemble_engine(rng, abc, M, ret_hmm, ret_sq, ret_anch, opt_mo, opt_pi, opt_tsc, FALSE, TRUE, FALSE);
}

/* Function:  h4_modelsample_AnchoredMulti()
 * Synopsis:  Make profile/seq/anchorset such that all paths use anchorset
 * Incept:    SRE, Mon 01 Feb 2021
 *
 * Purpose:   Last of three ways to contrive a profile, sequence, and
 *            anchor set so that all valid paths must use the anchor
 *            set.
 *
 *            In this version (multi), we make a multiglocal profile that
 *            only emits from match states (L=0, inserts prohibited), 
 *            where the anchor Mk0 state must emit X and no other Mk can,
 *            and where paths are forced to use Mk0 not Dk0. The sequence
 *            contains D residues X, marking the anchors, one per domain.
 *
 *            (Forcing e_k0(X)=1 and transitions to use Mk0 not Dk0
 *            are both necessary, to prevent path from including
 *            additional domains that avoid an X.)
 *
 *            This version differs from the others by testing multiple
 *            domains per path.
 *
 * Args:      input:
 *            rng      - random number generator
 *            abc      - digital alphabet 
 *            M        - length of model to sample (>= 1)
 *
 *            output:
 *            ret_hmm  - sampled profile \theta
 *            ret_sq   - sampled seq x
 *            ret_anch - sampled anchorset A 
 *
 *            optional output:
 *            opt_mo   - comparison mode (multiglocal, L=0)
 *            opt_pi   - path that <sq> was obtained from
 *            opt_tsc  - score of that path
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
h4_modelsample_AnchoredMulti(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M,
                             H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                             H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_tsc)
{
  return anchored_ensemble_engine(rng, abc, M, ret_hmm, ret_sq, ret_anch, opt_mo, opt_pi, opt_tsc, FALSE, FALSE, TRUE);
}


/*****************************************************************
 * 5. Internal (static) functions
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

/* single_path_seq_engine()
 * SRE, Sun 10 Jan 2021
 *
 * Guts of h4_modelsample_SinglePathSeq(), h4_modelsample_SinglePathASC()
 * (controlled by <do_multi> FALSE or TRUE, respectively).
 *
 * With <do_multi=FALSE>:
 * Create a profile/sequence pair for which only one path is possible:
 * P(\pi | x, \theta) = 1, thus P(x, \pi | \theta) = P(x | \theta),
 * thus Fwd = Bck = Vit = ASC Fwd = ASC Bck = ASC Vit score.
 *
 * With <do_multi=TRUE>:
 * Create a profile/sequence/anchorset _triplet_ for which only one
 * path is possible:
 * P(\pi | x, A, \theta) = 1,
 * thus P(x, A, \pi | \theta) = P(x, A | \theta);
 * thus ASC Fwd = ASC Bck = ASC Vit score. 
 *
 * The idea of the first version is to force the model to visit a
 * specific set of match states per domain, by setting transition
 * probabilities such that we visit Mk with probability 1.0 or Dk with
 * probability 1.0, for all k=1..M; to set all the match emission
 * probabilities to have zero probability for a particular residue x;
 * and to set all NCJI-emitted positions in the sequence to x.
 *
 * Additionally, we need to avoid I->D->I; inserted x's must be
 * unambiguously associated with a particular Ik state. Thus, set all
 * tDI = 0.0. We also need to prohibit I's before the first match
 * state and after the last match state, to avoid ambiguity with N 
 * and C assignment; tDI=0 suffices for the first, but we additionally
 * set tMI=0 for the last occupied match state to achieve the second.
 *
 * So if a M=6 uniglocal profile has four visited match states M2 M3 M4 M6 that emit non-x 
 * residues a,b,c,d, and the sequence looks like:
 *     x x a b x x c x d x x x 
 * there is an unambiguous assignment to a single path:
 *     S N N N B G D1 M2 M3 I3 I3 M4 I4 D5 M6 E C C C C T
 *     . . x x . .  -  a  b  x  x  c  x  -  d . . x x x .
 *
 * The model must use unihit glocal-only mode for this to work, but it
 * can use any target length model L>0. With a multihit model, when
 * the sequence contains two (or more) domains, other paths become
 * possible that have fewer domains, using NCJI to account for the
 * remaining ones.
 *
 * The idea can be extended to multihit glocal if we use an anchor set
 * to constrain to the set of paths to include all domains, precluding
 * these paths that use NCJI to account for one or more of
 * them. That's what the <do_multi=TRUE> version does, by creating a
 * profile/sequence/anchorset triple, instead of just a
 * profile/sequence pair.
 *
 * input:
 *    rng  - random number generator
 *    abc  - digital alphabet
 *    M    - length of model to sample (>=1)
 * control options:
 *    do_multi - TRUE to do multihit profile/seq/anchorset version
 * output:
 *    ret_hmm - sampled profile \theta
 *    ret_sq  - sampled sequence x
 * optional output:
 *    opt_mo   - comparison mode (uniglocal or multiglocal, with target L set)
 *    opt_pi   - sampled path \pi
 *    opt_anch - sampled anchorset A (available in unihit version too, not just multihit)
 *    opt_sc   - raw bitscore of the single path
 */
static int
single_path_seq_engine(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int do_multi,
                       H4_PROFILE **ret_hmm, ESL_SQ **ret_sq,
                       H4_MODE **opt_mo, H4_PATH **opt_pi, H4_ANCHORSET **opt_anch, float *opt_sc)
{
  H4_PROFILE   *hmm  = NULL;
  ESL_SQ       *sq   = NULL;
  H4_PATH      *pi   = NULL;
  H4_MODE      *mo   = NULL;
  H4_ANCHORSET *anch = NULL;
  int          *mocc = NULL;    // (0) 1..M array of 0|1 flags for which nodes use Mk; at least one is up
  int           nm   = 0;       // counts how many match states we visit, so we can reject mute all-D paths
  int           lastm;          // k index of the last M node
  int           D;              // number of domains 
  int           i,k,z,r,s,d;
  ESL_DSQ       xx;             // prohibited residue for Mk emissions (0..K-1)
  float         sc;             // path score
  float         tII_max = 0.99; // don't let tII be ~1, or we get huge insert lengths
  int           status;

  ESL_DASSERT1(( M >= 1 ));

  if (( hmm = h4_profile_Create(abc, M)) == NULL) { status = eslEMEM; goto ERROR; }
  if (( sq  = esl_sq_CreateDigital(abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if (( pi  = h4_path_Create())          == NULL) { status = eslEMEM; goto ERROR; }
  if (( mo  = h4_mode_Create())          == NULL) { status = eslEMEM; goto ERROR; }

  /* choose the residue that match states can't emit */
  xx = esl_rnd_Roll(rng, abc->K);

  /* choose the pattern of M state occupancy */
  ESL_ALLOC(mocc, sizeof(int) * (M+1));
  mocc[0] = 0;
  while (nm == 0)
    {
      nm = 0;
      for (k = 1; k <= M; k++)
        if (esl_rnd_Roll(rng, 2)) { mocc[k] = 1; nm++; lastm = k; }
        else                      { mocc[k] = 0; }
    }
  
  /* sample the profile, with transitions and emissions constrained  
   */
  /* reset emission/transition probability params to zero, then set edge conventions */
  esl_mat_FSet(hmm->t, M+1, h4_NT,   0.0f);
  esl_mat_FSet(hmm->e, M+1, abc->K,  0.0f);
  h4_profile_SetConventions(hmm);

  /* Match emissions */
  for (k = 1; k <= M; k++)
    {
      esl_dirichlet_FSampleUniform(rng, abc->K-1, hmm->e[k]); // sample K-1 nonzero probs
      ESL_SWAP(hmm->e[k][abc->K-1], hmm->e[k][xx], float);    // swap 0.0 prob onto xx residue
    }

  /* Initialize G->M1|D1 transition so 50% of the time we start 100% in M1, else we start 100% in D1 */
  if (mocc[1]) hmm->t[0][h4_TMM] = 1.0; 
  else         hmm->t[0][h4_TMD] = 1.0;   

  for (k = 1; k < M; k++)
    {
      // code here ASSUMES TRANSITIONS IN THIS ORDER: MM | MI | MD | IM | II | ID | DM | DI | DD  (see h4_profile.h)
      if (mocc[k+1])                                                     // 50% of the time, make all paths go to Mk+1
        {                                                 
          esl_dirichlet_FSampleUniform(rng, 2, hmm->t[k]);               // tmd=0
          do {
            esl_dirichlet_FSampleUniform(rng, 2, hmm->t[k]+h4_TIM);      // tid = 0
          } while (hmm->t[k][h4_TII] > tII_max);                         //  ... but keep tii away from ~1.0
          hmm->t[k][h4_TDM] = 1.0;                                       //  ... and prohibit tdi to avoid alignment ambiguity
        }
      else                                                               // or, make all paths go to Dk+1
        {
          if (k < lastm) esl_dirichlet_FSampleUniform(rng, 2, hmm->t[k]+h4_TMI);  // tmm=0
          else hmm->t[k][h4_TMD] = 1.0;                                  // ... prohibit tmi for last match state to avoid ali ambiguity
          do {
            esl_dirichlet_FSampleUniform(rng, 2, hmm->t[k]+h4_TII);      // tim = 0
          } while (hmm->t[k][h4_TII] > tII_max);                         //  ... but keep tii away from ~1.0
          hmm->t[k][h4_TDD] = 1.0;                                       //  ... and prohibit tdi to avoid alignment ambiguity
        }
    }

  /* configure alignment mode */
  if (do_multi) { if (( status = h4_mode_SetGlocal(mo))    != eslOK) goto ERROR; }
  else          { if (( status = h4_mode_SetUniglocal(mo)) != eslOK) goto ERROR; }
  if (( status = h4_mode_SetLength(mo, 10)) != eslOK) goto ERROR;                  // set a small L=10 to avoid long NCJ stretches

  /* emit a sequence and path from profile */
  if (( status = h4_emit(rng, hmm, mo, sq, pi)) != eslOK) goto ERROR;

  /* give seq a name */
  if (( status = esl_sq_SetName(sq, "testseq")) != eslOK) goto ERROR;

  /* Change all N/C/J/I emitted residues to xx (the residue M states can't emit).
   * Also, while we're going through, count domains <D>, by counting G states. 
   * (No L states occur, because we're either uniglocal or multiglocal.)
   */
  i = 1;
  D = 0;
  for (z = 0; z < pi->Z; z++)
    if      (h4_path_IsX(pi->st[z])) for (r = 0; r < pi->rle[z]-1; r++) sq->dsq[i++] = xx;  // NJC: rle-1 because they emit on transition
    else if (h4_path_IsI(pi->st[z])) for (r = 0; r < pi->rle[z];   r++) sq->dsq[i++] = xx;
    else if (h4_path_IsM(pi->st[z])) i += pi->rle[z];
    else if (pi->st[z] == h4P_G)     D++;
  ESL_DASSERT1(( i == sq->n+1 ));
  
  /* Select an anchor set
   * We know that each domain uses <nm> match states, and any of them
   * is a suitable anchor; choose randomly for each domain.
   */
  if (( anch = h4_anchorset_Create(D, sq->n, hmm->M)) == NULL) { status = eslEMEM; goto ERROR; }
  i = 1;
  d = 0;
  for (z = 0; z < pi->Z; z++)
    {
      if      (pi->st[z] == h4P_G)     { d++; k = 1; s = esl_rnd_Roll(rng, nm); } // use s'th match state as anchor for domain d (s=0..nm-1)
      else if (h4_path_IsX(pi->st[z])) i += pi->rle[z]-1;
      else if (h4_path_IsI(pi->st[z])) i += pi->rle[z];
      else if (h4_path_IsD(pi->st[z])) k += pi->rle[z];
      else if (h4_path_IsM(pi->st[z]))
        for (r = 0; r < pi->rle[z]; r++)
          {
            if (s == 0) { anch->a[d].i0 = i; anch->a[d].k0 = k; }
            s--; k++; i++;
          }
    }

  /* We have a probability model; set standard and vectorized lod scores */
  hmm->flags |= h4_HASPROBS;
  if (( status = h4_standardize(hmm)) != eslOK) goto ERROR;
  if (( status = h4_vectorize(hmm))   != eslOK) goto ERROR;
 
  /* Configure length model before scoring the trace */
  if (( status = h4_mode_SetLength(mo, sq->n)) != eslOK) goto ERROR;
              
  /* Calculate the score of the path.
   * This has to come after length model config.
   */
  if (( status = h4_path_Score(pi, sq->dsq, hmm, mo, &sc)) != eslOK) goto ERROR;

  free(mocc);
  if (opt_mo)   *opt_mo   = mo;   else h4_mode_Destroy(mo);
  if (opt_pi)   *opt_pi   = pi;   else h4_path_Destroy(pi);
  if (opt_anch) *opt_anch = anch; else h4_anchorset_Destroy(anch);
  if (opt_sc)   *opt_sc   = sc;   
  *ret_hmm = hmm;
  *ret_sq  = sq;
  return eslOK;

 ERROR:
  if (opt_mo)   *opt_mo   = NULL;
  if (opt_pi)   *opt_pi   = NULL;
  if (opt_anch) *opt_anch = NULL;
  if (opt_sc)   *opt_sc   = -eslINFINITY;
  *ret_hmm = NULL;
  *ret_sq  = NULL;

  free(mocc);
  h4_mode_Destroy(mo);
  h4_path_Destroy(pi);
  h4_anchorset_Destroy(anch);
  h4_profile_Destroy(hmm);
  esl_sq_Destroy(sq);
  return status;
}


/* anchored_ensemble_engine()
 *
 * Used by h4_modelsample_AnchoredUni()   (with do_uni   = TRUE)
 * and     h4_modelsample_AnchoredLocal() (with do_local = TRUE)
 * and     h4_modelsample_AnchoredMulti() (with do_multi = TRUE)
 *
 * The idea of all of three versions is to choose a particular residue
 * in the alphabet, X, to put at the anchor(s), then force all
 * paths to align Mk0 to that X.  Now the ensemble of all valid paths
 * must pass through the anchor i0,k0, and:
 *    - standard F/B score = ASC F/B score
 *    - std and ASC Decoding matrices identical
 *    - std and ASC F/B matrices identical to a well-defined extent
 * 
 * In uni mode: a unihit glocal profile with anchor Mk0 state that
 * must emit X, transitions force use of Mk0 not Dk0, and a target
 * sequence that has only one X in it.
 *
 * In local mode: a unihit local/glocal profile that only emits from
 * match states (L=0, and inserts prohibited), only the anchor Mk0
 * state has e_k0(X) > 0; and the sequence only has one X in it.
 *
 * In multi mode: a multiglocal profile that only emits from match
 * states (L=0 and inserts prohibited), the anchor Mk0 state must emit
 * X and no other Mk state can, and paths are forced to use Mk0 not
 * Dk0. The sequence contains D residues X, marking the anchors for
 * D domains.
 *
 * Summary table [H10/22]:
 *                                                 uni         local       multi
 *                                              ---------   ----------   ------------
 * 1. Sample model of length M                  :::::::::::::: all ::::::::::::::::::
 * 2. Select k0, anchX                          :::::::::::::: all ::::::::::::::::::
 * 3. Prohibit insert use, tXI = 0                  no         YES          YES
 * 4. Force Mk0 visit, t_{k0-1}(XD) = 0            YES          no          YES
 * 5. All Mk k != k0 have e_k(X) = 0                -          YES          YES
 *                 ... or e_k(X) < 1               YES          -            -
 * 6. Mk0 has e_k0(X) = 1                          YES          -           YES
 *     ... or e_k0(X) > 0                           -          YES           -
 * 7. Alignment mode                             uni/glocal  uni/dual   multi/glocal
 * 8. Target length                                L>0          0            0
 * 9. Sampled path must include Mk0             (by const.)    YES       (by const.)
 * 10. Alter seq so Mk0 emits X                 (by const.)    YES       (by const.)
 * 11. Alter seq so no other Mk emits X            YES      (by const.)  (by const.)
 * 12. Make anchor set                         :::::::::::::: all ::::::::::::::::::
 */
static int
anchored_ensemble_engine(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M,
                         H4_PROFILE **ret_hmm, ESL_SQ **ret_sq, H4_ANCHORSET **ret_anch,
                         H4_MODE **opt_mo, H4_PATH **opt_pi, float *opt_tsc,
                         int do_uni, int do_local, int do_multi)
{
  H4_PROFILE   *hmm  = NULL;
  ESL_SQ       *sq   = NULL;
  H4_PATH      *pi   = NULL;
  H4_MODE      *mo   = NULL;
  H4_ANCHORSET *anch = NULL;
  int           k0;
  int           z,r,i,k,d;
  int           nanch, D;
  ESL_DSQ       anchX;
  int           status;

  ESL_DASSERT1(( M >= 1 ));
  ESL_DASSERT1(( abc->K  >= 2   ));          // avoid a pathological case
  ESL_DASSERT1(( do_uni   == TRUE || do_uni   == FALSE));
  ESL_DASSERT1(( do_local == TRUE || do_local == FALSE));
  ESL_DASSERT1(( do_multi == TRUE || do_multi == FALSE));
  ESL_DASSERT1(( do_uni + do_local + do_multi == 1 ));
    

  /* 1. Start with a randomly sampled model, that we'll tweak */
  if (( status = h4_modelsample(rng, abc, M, &hmm)) != eslOK) goto ERROR;

  if (( sq  = esl_sq_CreateDigital(abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if (( pi  = h4_path_Create())          == NULL) { status = eslEMEM; goto ERROR; }
  if (( mo  = h4_mode_Create())          == NULL) { status = eslEMEM; goto ERROR; }

  /* 2. Select anchor position k0 in model (1..M), and
   *    select residue that only the anchor Mk0 will be able to emit with e_k0>0.
   */
  k0    = 1 + esl_rnd_Roll(rng, M);
  anchX = esl_rnd_Roll(rng, abc->K);

  /* 3. Prohibit inserts, by making all tMI, tDI = 0  (local, multi; not uni)
   *    In local, multi versions, only match states generate residues.
   */
  if (do_local || do_multi)
    {
      for (k = 1; k < M; k++)
        {
          hmm->t[k][h4_TMM] = esl_random(rng);
          hmm->t[k][h4_TMI] = 0.0;
          hmm->t[k][h4_TMD] = 1.0 - hmm->t[k][h4_TMM];

          hmm->t[k][h4_TDM] = esl_random(rng);
          hmm->t[k][h4_TDI] = 0.0;
          hmm->t[k][h4_TDD] = 1.0 - hmm->t[k][h4_TDM];
        }
    }

  /* 4. Force Mk0 to be visited, not Dk0 (uni, multi; not local)
   */
  if (do_uni || do_multi)
    {
      hmm->t[k0-1][h4_TMM] += hmm->t[k0-1][h4_TMD];
      hmm->t[k0-1][h4_TMD]  = 0.0;
      hmm->t[k0-1][h4_TIM] += hmm->t[k0-1][h4_TID];
      hmm->t[k0-1][h4_TID]  = 0.0;
      hmm->t[k0-1][h4_TDM] += hmm->t[k0-1][h4_TDD];
      hmm->t[k0-1][h4_TDD]  = 0.0;
    }
  
  /* 5. All Mk!=k0 have e_k(X) = 0 (local, multi) 
   *                   or just < 1 (uni)
   */
  if (do_local || do_multi)
    {
      for (k = 1; k <= M; k++)
        if (k != k0)
          {
            hmm->e[k][abc->K-1] = 0.0;
            esl_dirichlet_FSampleUniform(rng, abc->K-1, hmm->e[k]);
            ESL_SWAP(hmm->e[k][abc->K-1], hmm->e[k][anchX], float);
          }
    }
  else 
    {
      for (k = 1; k <= M; k++)
        if (k != k0)
          while (hmm->e[k][anchX] > 0.5) // 0.5 is arbitrary; we need <1, and we also want to be able to sample y!=X efficiently below
            esl_dirichlet_FSampleUniform(rng, abc->K, hmm->e[k]);
    }

  /* 6. Mk has e_k0(X) = 1 (uni, multi)
   *           or just > 0 (local)
   */
  if (do_uni || do_multi)
    {
       esl_vec_FSet(hmm->e[k0], abc->K, 0.0);
       hmm->e[k0][anchX] = 1.0;
    }
  else 
    {
      while (hmm->e[k0][anchX] == 0.0f)
        esl_dirichlet_FSampleUniform(rng, abc->K, hmm->e[k0]);
    }
  
  /* 7. Configure alignment mode */
  if (do_uni)   { if ((status = h4_mode_SetUniglocal(mo)) != eslOK) goto ERROR; }
  if (do_local) { if ((status = h4_mode_SetUnihit(mo))    != eslOK) goto ERROR; }
  if (do_multi) { if ((status = h4_mode_SetGlocal(mo))    != eslOK) goto ERROR; }

  /* 8. Configure target length model */
  if (do_uni)   { if ((status = h4_mode_SetLength(mo, 10)) != eslOK) goto ERROR; }
  if (do_local) { if ((status = h4_mode_SetLength(mo,  0)) != eslOK) goto ERROR; }
  if (do_multi) { if ((status = h4_mode_SetLength(mo,  0)) != eslOK) goto ERROR; }

  /* 9. Emit a path and a sequence.
   *    In local mode, reject paths that don't use Mk0 at k0 anchor.
   *    (In uni and multi mode, all paths always use Mk0 anchor anyway.)
   */
  do {
    if (( status = h4_emit(rng, hmm, mo, sq, pi)) != eslOK) goto ERROR;

    nanch = D = 0;
    for (z = 0; z < pi->Z; z++)
      if      (pi->st[z] == h4P_G)     { D++; k = 1; }
      else if (pi->st[z] == h4P_L)     { D++; k = pi->rle[z]; }
      else if (h4_path_IsD(pi->st[z])) k += pi->rle[z];
      else if (h4_path_IsM(pi->st[z])) {
        if (k0 >= k && k0 < k+pi->rle[z]) nanch++;
        k += pi->rle[z];
      }
  } while (nanch != D);   // D is now the number of domains. We'll use it below.

  ESL_DASSERT1(( (do_uni && D==1) || (do_local && D==1) || (do_multi && D>0) ));

  if ((status = esl_sq_SetName(sq, "ensemble_test")) != eslOK) goto ERROR;

  /* 10. Alter sequence so Mk0 emits X (local; true by construction for uni, multi)
   *     and no other Mk emits X (uni; true by construction for local, multi)
   * 11. Make anchor set at the same time.    
   */
  if ((anch = h4_anchorset_Create(D, sq->n, hmm->M)) == NULL) goto ERROR;
  i = 1;
  d = 0;
  for (z = 0; z < pi->Z; z++) {  
    if      (pi->st[z] == h4P_G)     { d++; k =  1; }
    else if (pi->st[z] == h4P_L)     { d++; k =  pi->rle[z]; }
    else if (h4_path_IsD(pi->st[z])) {      k += pi->rle[z]; }
    else if (h4_path_IsX(pi->st[z]))
      for (r = 1; r < pi->rle[z]; r++) { // note starting at 1 for NJC's emitting on transition
        while (sq->dsq[i] == anchX) sq->dsq[i] = esl_rnd_FChoose(rng, hmm->f, abc->K);
        i++;
      }
    else if (h4_path_IsI(pi->st[z]))
      for (r = 0; r < pi->rle[z]; r++) {
        while (sq->dsq[i] == anchX) sq->dsq[i] = esl_rnd_FChoose(rng, hmm->f, abc->K);
        i++;
      }
    else if (h4_path_IsM(pi->st[z])) {
      for (r = 0; r < pi->rle[z]; r++) {
        if (k == k0) {
          anch->a[d].i0 = i;
          anch->a[d].k0 = k;
          sq->dsq[i]    = anchX;
        } else {
          while (sq->dsq[i] == anchX)  // must terminate because we already made sure e_k(y != x) has good prob
            sq->dsq[i] = esl_rnd_FChoose(rng, hmm->e[k], abc->K);
        }
        k++; i++;
      }
    }
  }
  ESL_DASSERT1(( d == D ));
  
  /* We've changed the probability model, so update lod scores */
  if (( status = h4_standardize(hmm)) != eslOK) goto ERROR;
  if (( status = h4_vectorize(hmm))   != eslOK) goto ERROR;

  /* For uni, which has a nonzero length model, configure length model
   * for actual seq length, before scoring path
   */
  if (do_uni) {
    if (( status = h4_mode_SetLength(mo, sq->n)) != eslOK) goto ERROR;
  }

  /* Obtain the path score, if caller wants it */
  if (opt_tsc && (status = h4_path_Score(pi, sq->dsq, hmm, mo, opt_tsc)) != eslOK) goto ERROR;

  if (opt_mo)   *opt_mo   = mo;   else h4_mode_Destroy(mo);
  if (opt_pi)   *opt_pi   = pi;   else h4_path_Destroy(pi);
  *ret_hmm  = hmm;
  *ret_sq   = sq;
  *ret_anch = anch;
  return eslOK;

 ERROR:
  if (opt_mo) *opt_mo = NULL; else h4_mode_Destroy(mo);
  if (opt_pi) *opt_pi = NULL; else h4_path_Destroy(pi);
  *ret_hmm  = NULL;                h4_profile_Destroy(hmm);
  *ret_sq   = NULL;                esl_sq_Destroy(sq);
  *ret_anch = NULL;                h4_anchorset_Destroy(anch);
  return status;
}


/*****************************************************************
 * 6. Unit tests
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
  ESL_SQ       *sq    = NULL;
  H4_MODE      *mo    = NULL;
  H4_PATH      *pi    = NULL;
  H4_ANCHORSET *anch  = NULL;
  float         sc;
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

  if ( h4_modelsample_SinglePathSeq(rng, abc, M, &hmm, &sq, &mo, &pi, &anch, &sc) != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(hmm, errbuf)       != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( esl_sq_Validate(sq, errbuf)            != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( h4_path_Validate(pi, M, sq->n, errbuf) != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( h4_anchorset_Validate(anch, errbuf)    != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( anch->D != 1 )                                   esl_fatal("h4_modelsample_SinglePathSeq is supposed to be unihit");
  h4_profile_Destroy(hmm);  esl_sq_Destroy(sq);  h4_mode_Destroy(mo);  h4_path_Destroy(pi);  h4_anchorset_Destroy(anch);

  if ( h4_modelsample_SinglePathSeq(rng, abc, M, &hmm, &sq, NULL, NULL, NULL, NULL) != eslOK) esl_fatal(msg); // optional args
  h4_profile_Destroy(hmm);  esl_sq_Destroy(sq);

  if ( h4_modelsample_SinglePathASC(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &sc) != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(hmm, errbuf)       != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( esl_sq_Validate(sq, errbuf)            != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( h4_path_Validate(pi, M, sq->n, errbuf) != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( h4_anchorset_Validate(anch, errbuf)    != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  h4_profile_Destroy(hmm);  esl_sq_Destroy(sq);  h4_anchorset_Destroy(anch);  h4_mode_Destroy(mo);  h4_path_Destroy(pi);  

  if ( h4_modelsample_SinglePathASC(rng, abc, M, &hmm, &sq, &anch, NULL, NULL, NULL) != eslOK) esl_fatal(msg);
  h4_profile_Destroy(hmm);  esl_sq_Destroy(sq);  h4_anchorset_Destroy(anch);

  if ( h4_modelsample_AnchoredUni(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &sc) != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(hmm, errbuf)       != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( esl_sq_Validate(sq, errbuf)            != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( h4_path_Validate(pi, M, sq->n, errbuf) != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( h4_anchorset_Validate(anch, errbuf)    != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  h4_profile_Destroy(hmm);  esl_sq_Destroy(sq);  h4_anchorset_Destroy(anch);  h4_mode_Destroy(mo);  h4_path_Destroy(pi);  

  if ( h4_modelsample_AnchoredLocal(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &sc) != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(hmm, errbuf)       != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( esl_sq_Validate(sq, errbuf)            != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( h4_path_Validate(pi, M, sq->n, errbuf) != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( h4_anchorset_Validate(anch, errbuf)    != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  h4_profile_Destroy(hmm);  esl_sq_Destroy(sq);  h4_anchorset_Destroy(anch);  h4_mode_Destroy(mo); h4_path_Destroy(pi); 

  if ( h4_modelsample_AnchoredMulti(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &sc) != eslOK) esl_fatal(msg);
  if ( h4_profile_Validate(hmm, errbuf)       != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( esl_sq_Validate(sq, errbuf)            != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( h4_path_Validate(pi, M, sq->n, errbuf) != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  if ( h4_anchorset_Validate(anch, errbuf)    != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
  h4_profile_Destroy(hmm);  esl_sq_Destroy(sq);  h4_anchorset_Destroy(anch);  h4_mode_Destroy(mo); h4_path_Destroy(pi); 

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
 * 7. Test driver
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



/*****************************************************************
 * 8. Example
 *****************************************************************/
#ifdef h4MODELSAMPLE_EXAMPLE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_anchorset.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"

#include "general.h"
#include "modelsample.h"


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                        docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help",             0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show HMMER version info",     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "example of using model sampling routines";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = h4_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(0);
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  int             M    = 10;
  H4_PROFILE     *hmm  = NULL;
  ESL_SQ         *sq   = NULL;
  H4_MODE        *mo   = NULL;
  H4_PATH        *pi   = NULL;
  H4_ANCHORSET   *anch = NULL;
  float           sc;

  // put whichever routine you want to play with here:
  //h4_modelsample_SinglePathSeq(rng, abc, M, &hmm, &sq, &mo, &pi, &anch, &sc);
  h4_modelsample_AnchoredUni(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &sc);
 
  h4_profile_Dump(stdout, hmm);
  esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, FALSE);
  h4_mode_Dump(stdout, mo);
  h4_path_Dump(stdout, pi);
  h4_anchorset_Dump(stdout, anch);
  printf("path score = %.3f\n", sc);

  h4_anchorset_Destroy(anch);
  h4_path_Destroy(pi);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}


#endif // h4MODELSAMPLE_EXAMPLE
