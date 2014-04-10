/* Sampling profile HMMs.
 * These routines are used primarily in unit testing.
 * The models are often contrived and/or constrained to enable particular tests.
 * 
 * Contents:
 *   1. Model sampling routines.
 *   2. Unit tests.
 *   3. Test driver.
 *   4. Example driver.
 *   5. Copyright and license information.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_dirichlet.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_vectorops.h"

#include "base/p7_bg.h"
#include "base/p7_coords2.h"
#include "base/p7_hmm.h"
#include "base/p7_prior.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "misc/emit.h"

#include "search/modelconfig.h"

#include "build/modelsample.h"


/*****************************************************************
 * 1. Model sampling routines.
 *****************************************************************/

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
  p7_hmm_SetComposition(hmm);
  
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
  P7_HMM *hmm      = NULL;
  char   *logmsg   = "[random HMM created by sampling prior]";
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
  p7_hmm_SetComposition(hmm);

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
  p7_hmm_SetComposition(hmm);

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
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[random enumerable HMM, all II=0, created by sampling]";
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
   * are interpreted as transitions to E. 
   */
  esl_dirichlet_FSampleUniform(r, abc->K, hmm->mat[M]); /* match emission probs  */
  esl_dirichlet_FSampleUniform(r, abc->K, hmm->ins[M]); /* insert emission probs */
  esl_dirichlet_FSampleUniform(r, 2, hmm->t[M]+p7H_MM);	/* DEPENDS ON ORDER OF TRANSITIONS: MM, MI, MD */
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
  p7_hmm_SetComposition(hmm);

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
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[HMM with uniform transitions, random emissions]";
  int     k;
  int     status;

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
  p7_hmm_SetComposition(hmm);

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
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[random single-pathed HMM, all t/e probs 0 or 1, created by sampling]";
  int     nm     = 0;		/* make sure the HMM uses at least one M state */
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
  p7_hmm_SetComposition(hmm);

  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}

  
/* see comments on p7_hmm_SampleSinglePathedSeq() and p7_hmm_SampleSinglePathedASC(), below;
 * this is the engine they both share.
 */
static int
sample_single_pathed_seq_engine(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
				P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
				P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc,
				int do_asc_version)
{
  char       *hmmname   = (do_asc_version ? "single_pathed_asc" : "single_pathed_seq");
  char       *logmsg    = (do_asc_version ? "[Test model produced by p7_hmm_SampleSinglePathedASC()]" : "[Test model produced by p7_hmm_SampleSinglePathedSeq()]");
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_SQ     *sq        = NULL;
  ESL_DSQ    *dsq       = NULL;
  int         L         = 0;
  P7_TRACE   *tr        = NULL;
  P7_COORD2  *anch      = NULL;
  int         D         = 0;
  float       sc        = 0.0;
  int         nm        = 0;
  int         notX;
  int         k,s,z,d;
  int         status;

  if (( hmm = p7_hmm_Create    (M, bg->abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if ((  gm = p7_profile_Create(M, bg->abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if ((  sq = esl_sq_CreateDigital(bg->abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if ((  tr = p7_trace_Create())             == NULL) { status = eslEMEM; goto ERROR; }

  /* Choose the residue that match states can't emit. */
  notX = esl_rnd_Roll(rng, bg->abc->K);

  /* Sample our specially constructed HMM.
   * With a wrapper that avoids the pathological case of
   * using no match states.
   */
  while (nm == 0)
    {
      p7_hmm_Zero(hmm);

      /* At node k=0 boundary:
       *    t[0][MM,MI,MD] = B->{MID} distribution. Profile config makes G->M1 = tMM+tMI; G->D1 = tMD.
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
      
      /* Match emissions. */
      hmm->mat[0][0] = 1.0;
      for (k = 1; k <= M; k++) 
	{
	  esl_dirichlet_FSampleUniform(rng, bg->abc->K-1, hmm->mat[k]);  // this samples K-1 nonzero probabilities,
	  ESL_SWAP(hmm->mat[k][bg->abc->K-1], hmm->mat[k][notX], float); // and this swap puts the 0.0 prob into the notX residue.
	}

      /* Insert emissions. */
      for (k = 0; k <= M; k++)
	esl_dirichlet_FSampleUniform(rng, bg->abc->K, hmm->ins[k]);     // remember, profile ignores insert emission probs, assumes e_Ik(x) = f(x) */

      /* Match transitions */
      for (k = 0; k <= M; k++) 
	if (esl_rnd_Roll(rng, 2) || k==M) 
	  esl_dirichlet_FSampleUniform(rng, 2, hmm->t[k]+p7H_MM);  // 50% of the time (and at k=M) t_MD = 0; t_MM,t_MI can be any distribution
	else
	  hmm->t[k][p7H_MD] = 1.;                                  // 50% of the time, t_MD = 1; t_MM,t_MI = 0, so profile will have tG->M1 = 0, tG->D1 = 1.

      /* Insert transitions */
      for (k = 0; k <= M; k++)                                     // Insert transitions are unconstrained,
	do {                                                       // but put a cap on tII, because tII~1.0 leads to infinite-length sequence samples 
	  esl_dirichlet_FSampleUniform(rng, 2, hmm->t[k]+p7H_IM);
	} while (hmm->t[k][p7H_II] > p7H_II_SAMPLE_MAX);
      
      /* Delete transitions */
      for (k = 0; k <= M; k++)
	hmm->t[k][p7H_DM + (k==M || k == 0 ? 0 : esl_rnd_Roll(rng, 2))] = 1.0f;    // Delete transitions constrained so (t_DM = 0 | t_DD = 0), 	

      /* Count the number of match states visited */
      s = nm = 0;
      for (k = 0; k < M; k++)
	{
	  if (s == 0) { if (hmm->t[k][p7H_MM] > 0) nm++; else s = 1; }
	  else if (hmm->t[k][p7H_DM] > 0) { nm++; s = 0; }
	}
    }

  /* Add mandatory annotation */
  p7_hmm_SetName(hmm, hmmname);
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  p7_hmm_SetCtime(hmm);
  p7_hmm_SetConsensus(hmm, NULL);
  p7_hmm_SetComposition(hmm);

  /* Create a profile from it.
   * We don't yet know actual length of the emitted seq, so choose  
   * a smallish L suitable for sampling smallish path from profile.
   */
  if (do_asc_version) p7_profile_ConfigGlocal   (gm, hmm, bg, 10);
  else                p7_profile_ConfigUniglocal(gm, hmm, bg, 10);

  /* Emit a trace from the profile */
  p7_ProfileEmit(rng, hmm, gm, bg, sq, tr);
  p7_trace_Index(tr);

  /* Extract the <dsq> from <sq> container.
   * We only do this because p7_ProfileEmit() takes <sq>, not <dsq> as arg 
   */
  dsq = sq->dsq;
  L   = sq->n;
  sq->dsq = NULL;
  sq->n   = 0;
  esl_sq_Destroy(sq);

  /* Doctor the emissions so that any N/C/J/I residue is residue notX */
  for (z = 1; z < tr->N-1; z++)
    if (tr->i[z] && ! p7_trace_IsM(tr->st[z])) dsq[tr->i[z]] = notX;

  /* Select an anchor set 
   * We know that each domain uses <nm> match states, and
   * any of them are suitable anchors; choose randomly for each 
   * domain.
   */
  ESL_ALLOC(anch, sizeof(P7_COORD2) * tr->ndom);
  D = tr->ndom;
  z = 0;
  for (d = 0; d < D; d++)
    {
      s = esl_rnd_Roll(rng, nm);   // use the s'th match state as anchor
      while (s || ! p7_trace_IsM(tr->st[z])) {
	if (p7_trace_IsM(tr->st[z])) s--;
	z++;
      }
      anch[d].n1 = tr->i[z];
      anch[d].n2 = tr->k[z];
      z++;
    }
  
  /* Configure final length model of <gm> before returning */
  p7_profile_SetLength(gm, L);
  
  if (opt_hmm)  *opt_hmm  = hmm;  else p7_hmm_Destroy(hmm);
  if (opt_gm)   *opt_gm   = gm;   else p7_profile_Destroy(gm); 
  if (opt_dsq)  *opt_dsq  = dsq;  else free(dsq);
  if (opt_L)    *opt_L    = L;   
  if (opt_tr)   *opt_tr   = tr;   else p7_trace_Destroy(tr);
  if (opt_anch) *opt_anch = anch; else free(anch);
  if (opt_D)    *opt_D    = D;
  if (opt_sc)   *opt_sc   = sc;
  return eslOK;

 ERROR:
  if (hmm)  p7_hmm_Destroy(hmm);
  if (gm)   p7_profile_Destroy(gm);
  if (sq)   esl_sq_Destroy(sq);
  if (dsq)  free(dsq);
  if (tr)   p7_trace_Destroy(tr);
  if (anch) free(anch);
  if (opt_hmm)  *opt_hmm  = NULL;
  if (opt_gm)   *opt_gm   = NULL;
  if (opt_dsq)  *opt_dsq  = NULL;
  if (opt_L)    *opt_L    = 0;
  if (opt_tr)   *opt_tr   = NULL;
  if (opt_anch) *opt_anch = NULL;
  if (opt_D)    *opt_D    = 0;
  if (opt_sc)   *opt_sc   = 0.;
  return status;
}

/* Function:  p7_hmm_SampleSinglePathedSeq()
 * Synopsis:  Sample model/seq pair such that only one P=1.0 path exists for seq.
 *
 * Purpose:   Sample a model and a sequence that have been contrived
 *            such that there is only a single possible path that can
 *            generate that sequence. That is, $P(\pi | x, M) = 1$ for
 *            one $\pi$. Nonetheless, the sequence may contain
 *            nonhomologous residues and arbitrary
 *            length deletions and insertions.
 *            
 *            Because $P(\pi | x, M) = 1$, we have $P(\pi, x | M) = P(
 *            x | M)$. Therefore the Viterbi, Forward, and Backward
 *            scores of the profile/sequence comparison are identical.
 *            Moreover, for any anchor set constructed by choosing any
 *            single M state for the single domain in $\pi$, the
 *            anchor-set-constrained Viterbi, Forward, Backward scores
 *            are identical to the unconstrained scores, and all valid
 *            cells in the ASC DP matrices can be exactly compared to
 *            corresponding cells in standard DP matrices.  But at the
 *            same time, both $P(x | \pi, M)$ and $P(\pi | M)$ are
 *            essentially arbitrary probabilities, so the situation
 *            creates a nontrivial test case for our various DP
 *            routines. That is, the model can generate multiple
 *            different state paths, and each path can emit multiple
 *            different sequences, but for the given sequence, only
 *            one path is possible.
 *            
 *            As input, caller provides a random number generator
 *            <rng>, a model length <M> for the sampled model to have,
 *            and a background null model <bg> (which we need both for
 *            the digital alphabet information, and for configuring a
 *            returned profile).
 *            
 *            Other arguments are optional result output. You may
 *            <NULL> for any feature you don't need.
 *            
 *            <opt_hmm> is the model. If you configure a profile from
 *            it, that profile must be put in glocal-only mode.
 *            
 *            <opt_gm> is a glocal-only profile built from <*opt_hmm>,
 *            with its length model configured for length <L>.
 *            
 *            <opt_dsq> is the digital sequence, and <opt_L> is its
 *            length. You would generally retrieve both or neither.
 *            
 *            <opt_tr> is the unique path by which <gm> can generate
 *            <dsq>. This traceback path is indexed with
 *            <p7_trace_Index()>.
 *            
 *            <opt_anch> is an anchor set, and <opt_D> is the number
 *            of anchors (domains).
 *            
 *            <opt_sc> is the score of the sequence: Viterbi, Forward,
 *            Backward, ASC Forward, ASC Backward are all identical.
 *            
 *            <opt_hmm>, <opt_gm>, <opt_dsq>, <opt_tr>, and <opt_anch>
 *            are allocated here, and the caller becomes responsible
 *            for freeing whichever of them it retrieved.
 *            
 *            The results are interrelated, of course. You can
 *            generate the <gm> from the <hmm> and <L> by
 *            <p7_profile_ConfigUniglocal()>. You can get the <dsq> and
 *            <L> from the traceback <tr>. You can get <tr> by a
 *            Viterbi alignment of <gm> and <dsq>. You can get <sc> by
 *            either Viterbi, Forward, or Backward. You can get <anch>
 *            from <tr> by choosing any M state for each domain.
 *            
 *            Key idea is that we contrive an HMM that always
 *            generates a fixed # and sequence of match states
 *            (possibly involving D state usage), and match states
 *            cannot generate some residue X. Therefore for a sequence
 *            that has exactly one domain, that consists of X's for
 *            all N/C/I-emitted residues and nonX for all M-emitted
 *            residues, we can examine any sequence and uniquely
 *            assign match states to it, and once we've done that,
 *            we've uniquely defined the rest of the path.
 *            
 *            To guarantee a fixed # of match states, the profile must
 *            be configured in unihit glocal mode, and transitions are set
 *            such that (t_MD = 0 or t_MM = 0), and (t_DM = 0 or t_DD
 *            = 0). To guarantee that we can uniquely identify which
 *            residues must be generated by match states, match
 *            emission distributions have e_X = 0 for a randomly
 *            chosen residue X, and when we construct the target
 *            sequence, we make all N/C/I emitted residues X's. To
 *            avoid a mute cycle, we require that the path uses at
 *            least one match state.
 *            
 *
 * Args:      rng      ; random number generator
 *            M        : length of <hmm>/<gm> to sample
 *            bg       : background model
 *            opt_hmm  : optRETURN: model
 *            opt_gm   : optRETURN: profile (glocal, L)
 *            opt_dsq  : optRETURN: digital sequence
 *            opt_L    : optRETURN: length of <dsq>
 *            opt_tr   : optRETURN: traceback, state path and emitted seq
 *            opt_anch : optRETURN: array of anchors i0,k0 for each domain d
 *            opt_D    : optRETURN: number of domains d in <anch>, <tr>
 *            opt_sc   : optRETURN: raw lod score, in nats
 *            
 * Returns:   <eslOK> on success, and any results that were requested
 *            by passing non-<NULL> pointers for result arguments are
 *            now available.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            On any exception, any requested result pointer is returned <NULL>,
 *            and any requested return value is returned 0.0.
 */
int
p7_hmm_SampleSinglePathedSeq(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
			     P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
			     P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc)
{
  return sample_single_pathed_seq_engine(rng, M, bg, opt_hmm, opt_gm, opt_dsq, opt_L, opt_tr, opt_anch, opt_D, opt_sc, /*do_asc_version=*/FALSE);
}

/* Function:  p7_hmm_SampleSinglePathedASC()
 * Synopsis:  Sample profile/seq/anch triplet such that only one P=1.0 path exists for anchored seq comparison.
 *
 * Purpose:   Sample a profile, sequence, anchor set combination that have been
 *            contrived such that there is only a single possible path that
 *            can generate that sequence and satisfy the anchor set constraint.
 *            That is, $P(\pi | x, A, M) = 1$ for one and only one $\pi$.
 *            
 *            Essentially the same as <p7_hmm_SampleSinglePathedSeq(),
 *            so for more information, see documentation there.
 *            
 *            The difference is that this profile is in multiglocal
 *            mode, and the sequence may contain more than one domain.
 *            (Because the anchor set constraint defines the unique
 *            path for a multidomain sequence.) The standard Forward
 *            score will not equal the ASC Forward score for
 *            multidomain sequences.
 */
int
p7_hmm_SampleSinglePathedASC(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
			     P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
			     P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc)
{
  return sample_single_pathed_seq_engine(rng, M, bg, opt_hmm, opt_gm, opt_dsq, opt_L, opt_tr, opt_anch, opt_D, opt_sc, /*do_asc_version=*/TRUE);
}



/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef p7MODELSAMPLE_TESTDRIVE

static void
utest_sample(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char    failmsg[] = "modelsample::p7_hmm_Sample() unit test failed";
  P7_HMM *hmm       = NULL;
  int     ntrials   = 10;
  char    errmsg[eslERRBUFSIZE];
  int     i;
  int     status;
  
  for (i = 0; i < ntrials; i++)
    {
      if (( status = p7_hmm_Sample(rng, M, abc, &hmm))     != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001)) != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      p7_hmm_Destroy(hmm);
    }
}

static void
utest_sample_prior(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char      failmsg[] = "modelsample::p7_hmm_SamplePrior() unit test failed";
  P7_HMM   *hmm       = NULL;
  P7_PRIOR *pri       = NULL;
  int       ntrials   = 10;
  char      errmsg[eslERRBUFSIZE];
  int       i;
  int       status;

  if      (abc->type == eslAMINO) pri = p7_prior_CreateAmino();
  else if (abc->type == eslDNA)   pri = p7_prior_CreateNucleic();
  else if (abc->type == eslRNA)   pri = p7_prior_CreateNucleic();
  else                            pri = p7_prior_CreateLaplace(abc);
  if (pri == NULL) esl_fatal(failmsg);
  
  for (i = 0; i < ntrials; i++)
    {
      if (( status = p7_hmm_SamplePrior(rng, M, abc, pri, &hmm)) != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001))       != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      p7_hmm_Destroy(hmm);
    }

  p7_prior_Destroy(pri);
}

static void
utest_sample_ungapped(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char      failmsg[] = "modelsample::p7_hmm_SampleUngapped() unit test failed";
  P7_HMM   *hmm       = NULL;
  ESL_SQ   *sq        = esl_sq_CreateDigital(abc);
  P7_TRACE *tr        = p7_trace_Create();
  int       ntrials   = 10;
  char      errmsg[eslERRBUFSIZE];
  int       i,z;
  int       status;

  if (sq == NULL) esl_fatal(failmsg);
  if (tr == NULL) esl_fatal(failmsg);
  
  for (i = 0; i < ntrials; i++)
    {
      if (( status = p7_hmm_SampleUngapped(rng, M, abc, &hmm))    != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001))        != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      if (( status = p7_CoreEmit(rng, hmm, sq, tr))               != eslOK) esl_fatal(failmsg);
      if (( status = p7_trace_Validate(tr, abc, sq->dsq, errmsg)) != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      
      for (z = 0; z < tr->N; z++)
	if (p7_trace_IsI(tr->st[z]) || p7_trace_IsD(tr->st[z])) esl_fatal(failmsg);

      esl_sq_Reuse(sq);
      p7_trace_Reuse(tr);
      p7_hmm_Destroy(hmm);
    }

  p7_trace_Destroy(tr);
  esl_sq_Destroy(sq);
}


static void
utest_sample_enumerable(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char      failmsg[] = "modelsample::p7_hmm_SampleEnumerable() unit test failed";
  P7_HMM   *hmm       = NULL;
  ESL_SQ   *sq        = esl_sq_CreateDigital(abc);
  P7_TRACE *tr        = p7_trace_Create();
  int       ntrials   = 10;
  char      errmsg[eslERRBUFSIZE];
  int       i,z;
  int       status;

  if (sq == NULL) esl_fatal(failmsg);
  if (tr == NULL) esl_fatal(failmsg);
  
  for (i = 0; i < ntrials; i++)
    {
      if (( status = p7_hmm_SampleEnumerable(rng, M, abc, &hmm))  != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001))        != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      if (( status = p7_CoreEmit(rng, hmm, sq, tr))               != eslOK) esl_fatal(failmsg);
      if (( status = p7_trace_Validate(tr, abc, sq->dsq, errmsg)) != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      
      for (z = 0; z < tr->N; z++)
	if (p7_trace_IsI(tr->st[z])) esl_fatal(failmsg);

      esl_sq_Reuse(sq);
      p7_trace_Reuse(tr);
      p7_hmm_Destroy(hmm);
    }

  p7_trace_Destroy(tr);
  esl_sq_Destroy(sq);
}


static void
utest_sample_enumerable2(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char      failmsg[] = "modelsample::p7_hmm_Sample_Enumerable2() unit test failed";
  P7_HMM   *hmm       = NULL;
  ESL_SQ   *sq        = esl_sq_CreateDigital(abc);
  P7_TRACE *tr        = p7_trace_Create();
  int       ntrials   = 10;
  char      errmsg[eslERRBUFSIZE];
  int       i,z;
  int       status;

  if (sq == NULL) esl_fatal(failmsg);
  if (tr == NULL) esl_fatal(failmsg);
  
  for (i = 0; i < ntrials; i++)
    {
      if (( status = p7_hmm_SampleEnumerable2(rng, M, abc, &hmm)) != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001))        != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      if (( status = p7_CoreEmit(rng, hmm, sq, tr))               != eslOK) esl_fatal(failmsg);
      if (( status = p7_trace_Validate(tr, abc, sq->dsq, errmsg)) != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      
      for (z = 0; z < tr->N-1; z++)
	if (p7_trace_IsI(tr->st[z]) && p7_trace_IsI(tr->st[z+1])) esl_fatal(failmsg);

      esl_sq_Reuse(sq);
      p7_trace_Reuse(tr);
      p7_hmm_Destroy(hmm);
    }

  p7_trace_Destroy(tr);
  esl_sq_Destroy(sq);
}

static void
utest_sample_uniform(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char    failmsg[] = "modelsample::p7_hmm_SampleUniform() unit test failed";
  P7_HMM *hmm       = NULL;
  int     ntrials   = 10;
  char    errmsg[eslERRBUFSIZE];
  int     i;
  int     status;
  
  for (i = 0; i < ntrials; i++)
    {
      if (( status = p7_hmm_SampleUniform(rng, M, abc, 0.1, 0.4, 0.1, 0.4, &hmm)) != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001))                        != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      p7_hmm_Destroy(hmm);
    }
}

static void
utest_sample_singlepathed(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "modelsample::p7_hmm_SampleSinglePathed() unit test failed";
  P7_HMM     *hmm       = NULL;
  ESL_SQ     *sq1       = esl_sq_CreateDigital(abc);
  ESL_SQ     *sq2       = esl_sq_CreateDigital(abc);
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_PROFILE *gm        = p7_profile_Create(M, abc);
  P7_TRACE   *tr1       = p7_trace_Create();
  P7_TRACE   *tr2       = p7_trace_Create();
  int         nhmm      = 10;
  int         ntrace    = 10;
  int         h,t;
  char        errmsg[eslERRBUFSIZE];
  int         status;

  if (sq1  == NULL || sq2 == NULL || bg == NULL || gm == NULL || tr1 == NULL || tr2 == NULL) esl_fatal(failmsg);
  
  for (h = 0; h < nhmm; h++)
    {
      if (( status = p7_hmm_SampleSinglePathed(rng, M, abc, &hmm))      != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001))              != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      if (( status = p7_CoreEmit(rng, hmm, sq1, tr1))                   != eslOK) esl_fatal(failmsg);
      if (( status = p7_trace_Validate(tr1, abc, sq1->dsq, errmsg))     != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      for (t = 1; t < ntrace; t++)
	{
	  if (( status = p7_CoreEmit(rng, hmm, sq2, tr2))               != eslOK) esl_fatal(failmsg);
	  if (( status = p7_trace_Validate(tr2, abc, sq2->dsq, errmsg)) != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
	  if (( status = p7_trace_Compare(tr1, tr2, 0.0))               != eslOK) esl_fatal(failmsg);
      
	  esl_sq_Reuse(sq2);
	  p7_trace_Reuse(tr2);
	}
      esl_sq_Reuse(sq1);
      p7_trace_Reuse(tr1);

      if ((status = p7_profile_ConfigUniglocal(gm, hmm, bg, 0))         != eslOK) esl_fatal(failmsg);  // must be uniglocal L=0
      if ((status = p7_profile_Validate(gm, errmsg, 0.0001))            != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      if ((status = p7_ProfileEmit(rng, hmm, gm, bg, sq1, tr1))         != eslOK) esl_fatal(failmsg);
      if (( status = p7_trace_Validate(tr1, abc, sq1->dsq, errmsg))     != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      for (t = 1; t < ntrace; t++)
	{
	  if (( status = p7_ProfileEmit(rng, hmm, gm, bg, sq2, tr2))    != eslOK) esl_fatal(failmsg);
	  if (( status = p7_trace_Validate(tr2, abc, sq2->dsq, errmsg)) != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
	  if (( status = p7_trace_Compare(tr1, tr2, 0.0))               != eslOK) esl_fatal(failmsg);
      
	  esl_sq_Reuse(sq2);
	  p7_trace_Reuse(tr2);
	}
      esl_sq_Reuse(sq1);
      p7_trace_Reuse(tr1);

      p7_profile_Reuse(gm);
      p7_hmm_Destroy(hmm);
    }

  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_trace_Destroy(tr1);
  p7_trace_Destroy(tr2);
  esl_sq_Destroy(sq1);
  esl_sq_Destroy(sq2);
}



static void
utest_sample_singlepathed_seq(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int do_asc_version)
{
  char        failmsg[] = "modelsample::p7_hmm_SampleSinglePathedSeq() unit test failed";
  P7_HMM     *hmm       = NULL;
  ESL_DSQ    *dsq1      = NULL;
  ESL_SQ     *sq2       = esl_sq_CreateDigital(abc);
  int         L1;
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_PROFILE *gm        = NULL;
  P7_TRACE   *tr1       = NULL;
  P7_TRACE   *tr2       = p7_trace_Create();
  P7_COORD2  *anch      = NULL;
  int         D1;
  float       sc1;
  int         nhmm      = 10;
  int         ntrace    = 10;
  int        *observ    = malloc(sizeof(int) * (M+1));
  int        *expect    = malloc(sizeof(int) * (M+1));
  int         h,t,z,k;
  char        errmsg[eslERRBUFSIZE];
  int         status;

  if (sq2 == NULL || bg == NULL || tr2 == NULL || expect == NULL || observ == NULL) esl_fatal(failmsg);

  for (h = 0; h < nhmm; h++)
    {
      if (do_asc_version) status = p7_hmm_SampleSinglePathedASC(rng, M, bg, &hmm, &gm, &dsq1, &L1, &tr1, &anch, &D1, &sc1);
      else                status = p7_hmm_SampleSinglePathedSeq(rng, M, bg, &hmm, &gm, &dsq1, &L1, &tr1, &anch, &D1, &sc1);
      if (status != eslOK) esl_fatal(failmsg);
	
      if (( status = p7_hmm_Validate    (hmm, errmsg, 0.0001))            != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      if (( status = p7_profile_Validate(gm,  errmsg, 0.0001))            != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      if (( status = p7_trace_Validate  (tr1, abc, dsq1, errmsg))         != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      /* expect[k] must be -D or +D; -D means position k is used by delete, +D means match */
      if (! do_asc_version && D1 != 1) esl_fatal(failmsg); /* SinglePathedSeq requires a single domain and unihit mode */
      esl_vec_ISet(expect, M+1, 0.);
      for (z = 0; z < tr1->N; z++) 
	if      (p7_trace_IsM(tr1->st[z])) expect[tr1->k[z]]++;
	else if (p7_trace_IsD(tr1->st[z])) expect[tr1->k[z]]--;
      for (k = 1; k <= M; k++)
	if (abs(expect[k]) != D1) esl_fatal(failmsg);

      for (t = 1; t < ntrace; t++)
	{
	  if (( status = p7_ProfileEmit(rng, hmm, gm, bg, sq2, tr2))        != eslOK) esl_fatal(failmsg);
	  if (( status = p7_trace_Validate(tr2, abc, sq2->dsq, errmsg))     != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

	  p7_trace_Index(tr2);
	  esl_vec_ISet(observ, M+1, 0.);
	  for (z = 0; z < tr2->N; z++) 
	    if      (p7_trace_IsM(tr2->st[z])) observ[tr2->k[z]]++;
	    else if (p7_trace_IsD(tr2->st[z])) observ[tr2->k[z]]--;

	  for (k = 1; k <= M; k++) {
	    if ( abs(observ[k]) != tr2->ndom) esl_fatal(failmsg);
	    if ( (observ[k] < 0 && expect[k] > 0) || (observ[k] > 0 && expect[k] < 0)) esl_fatal(failmsg);
	  }

	  esl_sq_Reuse(sq2);
	  p7_trace_Reuse(tr2);
	}

      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
      free(dsq1);
      p7_trace_Destroy(tr1);
      free(anch);      
    }

  free(observ);
  free(expect);
  esl_sq_Destroy(sq2);
  p7_bg_Destroy(bg);
  p7_trace_Destroy(tr2);
}



#endif /*p7MODELSAMPLE_TESTDRIVE*/
/*---------------------- end of unit tests -----------------------*/


/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7MODELSAMPLE_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for modelsample.c, routines for sampling test models";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  int             M    = 10;

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_sample                 (rng, M, abc);
  utest_sample_prior           (rng, M, abc);
  utest_sample_ungapped        (rng, M, abc);
  utest_sample_enumerable      (rng, M, abc);
  utest_sample_enumerable2     (rng, M, abc);
  utest_sample_uniform         (rng, M, abc);
  utest_sample_singlepathed    (rng, M, abc);
  utest_sample_singlepathed_seq(rng, M, abc, FALSE); /* tests p7_hmm_SampleSinglePathedSeq(), which uses uniglocal model   */
  utest_sample_singlepathed_seq(rng, M, abc, TRUE);  /* tests p7_hmm_SampleSinglePathedASC(), which uses multiglocal model */

  fprintf(stderr, "#  status = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  exit(0); /* success */
}

#endif /*p7MODELSAMPLE_TESTDRIVE*/
/*-------------------- end of test driver ---------------------*/


/*****************************************************************
 * 4. Example driver
 *****************************************************************/
#ifdef p7MODELSAMPLE_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"


#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "example driver for sampling test models";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  int             M       = 10;
  ESL_ALPHABET   *abc     = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg      = p7_bg_Create(abc);
  P7_HMM         *hmm     = NULL;
  P7_PROFILE     *gm      = NULL;
  ESL_DSQ        *dsq     = NULL;
  int             L       = 0;
  P7_TRACE       *tr      = NULL;
  P7_COORD2      *anch    = NULL;
  int             D       = 0;
  float           sc      = 0.;
  char            errmsg[eslERRBUFSIZE];
  
  p7_hmm_SampleSinglePathedSeq(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc);

  if (p7_hmm_Validate(hmm, errmsg, 0.0001) != eslOK ) esl_fatal(errmsg);

  p7_hmm_Dump           (stdout, hmm);
  p7_profile_Dump       (stdout, gm);
  p7_trace_DumpAnnotated(stdout, tr, gm, dsq);

  if (anch) free(anch);
  if (tr)   p7_trace_Destroy(tr);
  if (dsq)  free(dsq);
  if (gm)   p7_profile_Destroy(gm);
  if (hmm)  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(rng);
  return 0;
}



#endif /*p7MODELSAMPLE_EXAMPLE*/


/************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 ************************************************************/
