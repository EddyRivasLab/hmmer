/* Sampling random HMMs.
 * These routines are used primarily in unit testing.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_dirichlet.h"
#include "esl_random.h"

#include "base/p7_hmm.h"
#include "base/p7_prior.h"
#include "base/p7_hmm_sample.h"


/* Function:  p7_hmm_Sample()
 * Synopsis:  Sample an HMM at random.
 *
 * Purpose:   Creates a random HMM of length <M> nodes,
 *            for alphabet <abc>, obtaining randomness from
 *            <r>.
 * 
 *            Probably only useful for debugging.
 *            
 * Note:      Compare p7_hmm_Renormalize(), which has a similar
 *            structure, except it normalizes instead of
 *            sampling each probability vector.           
 *
 * Returns:   <eslOK> on success, and the new hmm is returned
 *            through <ret_hmm); caller is responsible for 
 *            freeing this object with <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_hmm_Sample(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[random HMM created by uniform sampling]";
  int     k;
  int     status;

  if ( (hmm = p7_hmm_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; }
  
  for (k = 0; k <= M; k++)
    {
      if (k > 0) esl_dirichlet_FSampleUniform(r, abc->K, hmm->mat[k]);
      esl_dirichlet_FSampleUniform(r, abc->K,    hmm->ins[k]);
      esl_dirichlet_FSampleUniform(r, p7H_NTMAT, P7H_TMAT(hmm, k));
      do {  /* put a cap on tII, because tII~1.0 leads to infinite-length sequence samples */
	esl_dirichlet_FSampleUniform(r, p7H_NTINS, P7H_TINS(hmm, k));
      } while (hmm->t[k][p7H_II] > p7H_II_SAMPLE_MAX);
      if (k > 0) esl_dirichlet_FSampleUniform(r, p7H_NTDEL, P7H_TDEL(hmm, k));
    }
  /* Node M is special: no transitions to D, transitions to M
   * are interpreted as transitions to E. Overwrite a little of
   * what we did in node M.
   */
  hmm->t[M][p7H_MD] = 0.0;
  esl_vec_FNorm(P7H_TMAT(hmm, M), p7H_NTMAT);
  
  esl_vec_FSet(P7H_TDEL(hmm, M), p7H_NTDEL, 0.0);
  hmm->t[M][p7H_DM] = 1.0;
  
  /* Add mandatory annotation, and some relevant optional annotation  */
  p7_hmm_SetName(hmm, "sampled-hmm");
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  p7_hmm_SetCtime(hmm);
  p7_hmm_SetConsensus(hmm, NULL);
  
  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}

/* Function:  p7_hmm_SamplePrior()
 * Synopsis:  Sample an HMM from a prior distribution.
 *
 * Purpose:   Creates a random HMM of length <M> nodes, for 
 *            alphabet <abc>, using random number generator
 *            <r>, by sampling from mixture Dirichlet priors
 *            given in <pri>.
 *            
 *            In general this should give more 'realistic' profile
 *            HMMs than <p7_hmm_Sample()> does.  In unit testing, we
 *            tend to use <p7_hmm_Sample()> for tests that should be
 *            able to deal with any profile, no matter how
 *            pathological. A few tests need to see typical,
 *            reasonable (high-posterior-probability) alignment paths,
 *            and we use <p7_hmm_SamplePrior()> in these cases.
 *            
 * Args:      r       : random number generator
 *            M       : length of profile HMM to sample
 *            abc     : alphabet
 *            pri     : mixture Dirichlet priors on emissions, transitions
 *            ret_hmm : RETURN: newly sampled profile HMM.
 *
 * Returns:   <eslOK> on success, and <*ret_hmm> contains the new profile
 *            HMM object. Caller is responsible for free'ing this 
 *            with <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_hmm_SamplePrior(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, const P7_PRIOR *pri, P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[random HMM created by sampling prior]";
  double  ep[p7_MAXABET];	/* tmp storage for sampled emission parameters as doubles */
  double  tp[p7H_NTMAX];	/* tmp storage, transitions */
  int     q,k;
  int     status;

  ESL_DASSERT1 ( ( M>0) );
  ESL_DASSERT1 ( ( pri->em->K == abc->K )   );
  ESL_DASSERT1 ( ( pri->ei->K == abc->K )   );
  ESL_DASSERT1 ( ( pri->tm->K == p7H_NTMAT) );
  ESL_DASSERT1 ( ( pri->ti->K == p7H_NTINS) );
  ESL_DASSERT1 ( ( pri->td->K == p7H_NTDEL) );
  
  if ( (hmm = p7_hmm_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; }

  for (k = 0; k <= M; k++)
    {
      if (k) {
	q = esl_rnd_DChoose(r, pri->em->pq, pri->em->N); 
	esl_dirichlet_DSample(r, pri->em->alpha[q], pri->em->K, ep); /* extra D2F step because Easel Dirichlet module is double-precision */
	esl_vec_D2F(ep, abc->K, hmm->mat[k]);
      }

      q = esl_rnd_DChoose(r, pri->ei->pq, pri->ei->N); 
      esl_dirichlet_DSample(r, pri->ei->alpha[q], pri->ei->K, ep);
      esl_vec_D2F(ep, abc->K, hmm->ins[k]);

      q = esl_rnd_DChoose(r, pri->tm->pq, pri->tm->N); 
      esl_dirichlet_DSample(r, pri->tm->alpha[q], pri->tm->K, tp);
      esl_vec_D2F(tp, p7H_NTMAT, P7H_TMAT(hmm,k));

      do {
	q = esl_rnd_DChoose(r, pri->ti->pq, pri->ti->N); 
	esl_dirichlet_DSample(r, pri->ti->alpha[q], pri->ti->K, tp);
	esl_vec_D2F(tp, p7H_NTINS, P7H_TINS(hmm,k));
      } while (hmm->t[k][p7H_II] > p7H_II_SAMPLE_MAX); /* put a cap on tII, because tII~1 gives us infinite-length sequence samples */

      if (k) {
	q = esl_rnd_DChoose(r, pri->td->pq, pri->td->N); 
	esl_dirichlet_DSample(r, pri->td->alpha[q], pri->td->K, tp);
	esl_vec_D2F(tp, p7H_NTDEL, P7H_TDEL(hmm,k));
      }
    }

  hmm->t[M][p7H_MD] = 0.0;
  esl_vec_FNorm(P7H_TMAT(hmm, M), p7H_NTMAT);
  
  esl_vec_FSet(P7H_TDEL(hmm, M), p7H_NTDEL, 0.0);
  hmm->t[M][p7H_DM] = 1.0;

  /* Add mandatory annotation, and some relevant optional annotation  */
  p7_hmm_SetName(hmm, "sampled-hmm");
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  p7_hmm_SetCtime(hmm);
  p7_hmm_SetConsensus(hmm, NULL);

  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  if (hmm) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}

/* Function:  p7_hmm_SampleUngapped()
 * Synopsis:  Sample a random HMM with no nonzero indel transitions.
 *
 * Purpose:   Same as <p7_hmm_Sample()>, except all 
 *            M $\rightarrow$ M transitions are 1.0:
 *            an ungapped model. Useful for testing 
 *            as a limit case.
 *            
 * Returns:   <eslOK> on success, and the new hmm is returned
 *            through <ret_hmm); caller is responsible for 
 *            freeing this object with <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      STL11/140
 */
int
p7_hmm_SampleUngapped(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  int     k;
  int     status;

  if ((status = p7_hmm_Sample(r, M, abc, &hmm)) != eslOK) goto ERROR;
  for (k = 0; k <= M; k++) {
    hmm->t[k][p7H_MM] = 1.0;
    hmm->t[k][p7H_MD] = 0.0;
    hmm->t[k][p7H_MI] = 0.0;
  }
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}

/* Function:  p7_hmm_SampleEnumerable()
 * Synopsis:  Sample an random HMM with no nonzero transitions to insert.
 *
 * Purpose:   Sample a random HMM with random emission and 
 *            transition probabilities with the exception that
 *            all transitions to insert are zero. This makes
 *            it possible to create a model with a finite,
 *            easily enumerable sequence space (all seqs of
 *            length $0..M$).
 *            
 *            To achieve finite enumerability in the profile as well
 *            as the core HMM, the caller must configure a unihit mode
 *            (<p7_profile_ConfigUnilocal(gm, hmm, bg, 0)> or
 *            <p7_profile_ConfigUniglocal()>, with a target length of
 *            zero.
 *            
 *            Useful for debugging and validating Forward/Viterbi
 *            algorithms.
 *
 *            Compare <p7_hmm_SampleEnumerable2()>, which only makes
 *            tII transitions zero, and thus enumerates a space consisting
 *            of sequences of length 0..2M+1 (and for an enumerable profile,
 *            1..2M-1).
 *            
 * Returns:   <eslOK> on success. The newly allocated hmm is returned through
 *            <ret_hmm>. The caller is responsible for freeing this object
 *            with <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_hmm_SampleEnumerable(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[random enumerable HMM created by sampling]";
  int     k;
  float   tmp[2];
  int     status;
  
  hmm = p7_hmm_Create(M, abc);
  if (hmm == NULL) { status = eslEMEM; goto ERROR; }

  for (k = 0; k <= M; k++)
    {
      if (k > 0) esl_dirichlet_FSampleUniform(r, abc->K, hmm->mat[k]); /* match emission probs  */
      esl_dirichlet_FSampleUniform(r, abc->K, hmm->ins[k]);            /* insert emission probs */
      esl_dirichlet_FSampleUniform(r, 2,      tmp);       
      hmm->t[k][p7H_MM] = tmp[0];
      hmm->t[k][p7H_MI] = 0.;
      hmm->t[k][p7H_MD] = tmp[1];
      hmm->t[k][p7H_IM] = 1.;                                          /* I transitions irrelevant since I's are unreached. */
      hmm->t[k][p7H_II] = 0.;
      if (k > 0) esl_dirichlet_FSampleUniform(r, 2,      hmm->t[k]+5); /* delete transitions to M,D */
    }

  /* Node M is special: no transitions to D, transitions to M
   * are interpreted as transitions to E. Overwrite a little of
   * what we did in node M.
   */
  hmm->t[M][p7H_MM] = 1.;
  hmm->t[M][p7H_MD] = 0.;	
  hmm->t[M][p7H_DM] = 1.;
  hmm->t[M][p7H_DD] = 0.;
  
  /* Add mandatory annotation
   */
  p7_hmm_SetName(hmm, "sampled-hmm");
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  p7_hmm_SetCtime(hmm);
  p7_hmm_SetConsensus(hmm, NULL);

#ifdef p7_DEBUGGING
  p7_hmm_Validate(hmm, NULL, 0.0001);
#endif

  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}



/* Function:  p7_hmm_SampleEnumerable2()
 * Synopsis:  Sample an random HMM with no nonzero insert-insert transitions.
 *
 * Purpose:   Sample a random HMM with random emission and 
 *            transition probabilities with the exception that
 *            all insert-insert transitions are zero. This makes
 *            it possible to create a model with a finite,
 *            easily enumerable sequence space (all seqs of
 *            length $\leq 2M+1).
 *            
 *            To achieve this in the profile as well as the core HMM,
 *            the caller must configure a unihit mode
 *            (<p7_profile_ConfigUnilocal(gm, hmm, bg, 0)> or
 *            <p7_profile_ConfigUniglocal()>, with a target length of
 *            zero.  The enumerable sequence space of the profile has
 *            $L=1..2M-1$, because I0 and IM are normalized away, and
 *            the B->D1..DM->E mute path (for a sequence of length 0)
 *            is also normalized away.
 *
 *            Useful for debugging and validating Forward/Viterbi
 *            algorithms. 
 *
 *            Compare <p7_hmm_SampleEnumerable()>, which makes all
 *            transitions to insert 0, and thus enumerates a smaller
 *            space.
 *            
 * Returns:   <eslOK> on success. The newly allocated hmm is returned through
 *            <ret_hmm>. The caller is responsible for freeing this object
 *            with <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_hmm_SampleEnumerable2(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm      = NULL;
  char   *logmsg   = "[random enumerable HMM, all II=0, created by sampling]";
  int     k;
  int     status;
  
  if (( hmm = p7_hmm_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; }

  for (k = 0; k < M; k++)
    {
      if (k > 0) esl_dirichlet_FSampleUniform(r, abc->K, hmm->mat[k]); /* match emission probs  */
      esl_dirichlet_FSampleUniform(r, abc->K, hmm->ins[k]);            /* insert emission probs */
      esl_dirichlet_FSampleUniform(r, 3,      hmm->t[k]);              /* match transitions     */
      hmm->t[k][p7H_IM] = 1.;	                                       /* tIM hardwired */
      hmm->t[k][p7H_II] = 0.;	                                       /* tII hardwired */
      if (k > 0) esl_dirichlet_FSampleUniform(r, 2,      hmm->t[k]+5); /* delete transitions */
    }

  /* Node M is special: no transitions to D, transitions to M
   * are interpreted as transitions to E. Overwrite a little of
   * what we just did in node M.
   */
  esl_dirichlet_FSampleUniform(r, abc->K, hmm->mat[M]); /* match emission probs  */
  esl_dirichlet_FSampleUniform(r, abc->K, hmm->ins[M]); /* insert emission probs */
  esl_dirichlet_FSampleUniform(r, 2, hmm->mat[M]);	/* DEPENDS ON ORDER OF TRANSITIONS: MM, MI, MD */
  hmm->t[M][p7H_MD] = 0.0;
  hmm->t[k][p7H_IM] = 1.;	                                       
  hmm->t[k][p7H_II] = 0.;	                                       
  hmm->t[M][p7H_DM] = 1.;
  hmm->t[M][p7H_DD] = 0.;
  
  /* Add mandatory annotation */
  p7_hmm_SetName(hmm, "sampled-hmm");
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  p7_hmm_SetCtime(hmm);
  p7_hmm_SetConsensus(hmm, NULL);

#ifdef p7_DEBUGGING
  p7_hmm_Validate(hmm, NULL, 0.0001);
#endif

  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}



/* Function:  p7_hmm_SampleUniform()
 * Synopsis:  Sample a model that uses fixed (given) transition probs.
 *
 * Purpose:   Sample a model that uses uniform transition probabilities,
 *            determined by <tmi>, <tii>, <tmd>, and <tdd>,
 *            the probabilistic equivalent of gap-open/gap-extend for
 *            inserts, deletes.
 *            
 *            Useful for testing expected behavior on single-sequence
 *            models, where transitions are position-independent.
 *
 * Returns:   <eslOK> on success, and the new hmm is returned
 *            through <ret_hmm); caller is responsible for 
 *            freeing this object with <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      J1/5.
 */
int
p7_hmm_SampleUniform(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, 
		     float tmi, float tii, float tmd, float tdd,
		     P7_HMM **ret_hmm)
{
  int     status;
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[HMM with uniform transitions, random emissions]";
  int     k;

  hmm = p7_hmm_Create(M, abc);
  if (hmm == NULL) { status = eslEMEM; goto ERROR; }
  
  for (k = 0; k <= M; k++)
    {
      if (k > 0) esl_dirichlet_FSampleUniform(r, abc->K, hmm->mat[k]);
      esl_dirichlet_FSampleUniform(r, abc->K, hmm->ins[k]);
      hmm->t[k][p7H_MM] = 1.0 - tmi - tmd;
      hmm->t[k][p7H_MI] = tmi;
      hmm->t[k][p7H_MD] = tmd;
      hmm->t[k][p7H_IM] = 1.0 - tii;
      hmm->t[k][p7H_II] = tii;
      hmm->t[k][p7H_DM] = 1.0 - tdd;
      hmm->t[k][p7H_DD] = tdd;
    }

  /* Deal w/ special stuff at node 0, M, overwriting some of what we
   * just did. 
   */
  hmm->t[M][p7H_MM] = 1.0 - tmi;
  hmm->t[M][p7H_MD] = 0.;
  hmm->t[M][p7H_DM] = 1.0;
  hmm->t[M][p7H_DD] = 0.;
  
  /* Add mandatory annotation
   */
  p7_hmm_SetName(hmm, "sampled-hmm");
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  p7_hmm_SetCtime(hmm);
  p7_hmm_SetConsensus(hmm, NULL);

  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}


/* Function:  p7_hmm_SampleSinglePathed()
 * Synopsis:  Sample a random HMM with only a single P=1.0 path possible.
 *
 * Purpose:   Sample a random HMM that only has 1.0 or 0.0 for its
 *            transition and emission probabilities. This makes it possible
 *            to create a profile that has only a single alignment path --
 *            a useful case for unambiguous testing and debugging.
 *            
 *            To do this, the caller also needs to configure unihit
 *            glocal mode and a target length of zero:
 *            <p7_profile_ConfigUniglocal(gm, hmm, bg, 0)>.
 *            
 * Returns:   <eslOK> on success, and the new HMM is returned through
 *            <ret_hmm>. Caller is responsible for freeing it with 
 *            <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_hmm_SampleSinglePathed(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  char   *logmsg   = "[random single-pathed HMM, all t/e probs 0 or 1, created by sampling]";
  P7_HMM *hmm      = NULL;
  int     nm       = 0;		/* make sure the HMM uses at least one M state */
  int     k;
  int     status;
  
  if ( (hmm = p7_hmm_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; }

  while (nm == 0) 
    {
      if ( (status = p7_hmm_Zero(hmm))  != eslOK) goto ERROR;

      /* At node k=0 boundary:
       *    t[0][MM,MI,MD] = B->{MID} distribution
       *    t[0][IM,II] = a normal I transition distribution (I0 exists)
       *    t[0][DM,DD] = 1.0,0.0 by convention (no D0 state)
       *    mat[0][...] = {1.0,0.0...} by convention (no M_0 state)
       *    ins[0][...] = a normal I emission distribution (I0 exists)
       * At node k=M boundary: 
       *    t[M][MM,MI,MD] = ME,MI, 0.0    
       *    t[M][IM,II]    = IE,II
       *    t[M][DM,DD]    = 1.0(DE), 0.0
       *    mat[M],ins[M] as usual
       * [xref hmmer.h::P7_HMM documentation]
       */
      for (k = 0; k <= M; k++)
	{
	  hmm->mat[k][(k==0 ? 0 : esl_rnd_Roll(r,abc->K))] = 1.0f;
	  hmm->ins[k][            esl_rnd_Roll(r,abc->K)]  = 1.0f;
      
	  /* code below relies on partial transition order, vectors 
	   * starting at 0, 3, 5: MM MI MD {IM II} DM DD 
	   */
	  hmm->t[k][esl_rnd_Roll(r, (k==M ? 2 : 3))]          = 1.0f; /* at k==M, MD=0.0 */
	  hmm->t[k][p7H_IM]                                   = 1.0f; /* can't set II to 1.0; infinite loop */
	  hmm->t[k][p7H_DM + (k==M ? 0 : esl_rnd_Roll(r, 2))] = 1.0f; /* at k==M, DM=1.0, DD=0.0 */
	}

      /* Make sure there's at least one M state in that path. 
       * Otherwise the HMM can't generate anything; only 
       * sequences of length zero are possible; and that just
       * makes life miserable in test code that uses single-pathed
       * HMMs. Just disallow that case.
       */
      if (hmm->t[0][p7H_MM]) nm++;
      for (k = 1; k < M; k++)
	if (hmm->t[k][p7H_DM]) nm++;
    }


  /* Add mandatory annotation */
  p7_hmm_SetName(hmm, "single-pathed-hmm");
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  p7_hmm_SetCtime(hmm);
  p7_hmm_SetConsensus(hmm, NULL);

#ifdef p7_DEBUGGING
  p7_hmm_Validate(hmm, NULL, 0.0001);
#endif

  *ret_hmm = hmm;
  return eslOK;

  return eslOK;

 ERROR:
  p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}

/************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 ************************************************************/
