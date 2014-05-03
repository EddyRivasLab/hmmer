/* Sampling profile HMMs.
 * 
 * These routines are used primarily in unit testing.  Models are
 * often contrived and/or constrained in various ways that enable
 * particular tests.
 * 
 * Contents:
 *   1. Model sampling routines.
 *   2. Models with finite enumerable path #
 *   3. Models with only a single valid path
 *   4. Models such that ASC paths == all paths 
 *   5. Internal (static) functions.
 *   6. Unit tests.
 *   7. Test driver.
 *   8. Example driver.
 *   9. Copyright and license information.
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

static int modelsample_engine(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm);
static int sample_single_pathed_seq_engine(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
					   P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
					   P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc,
					   int do_asc_version);
static int sample_anchored_engine(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
				  P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
				  P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc,
				  int do_local_version);

/*****************************************************************
 * 1. Model sampling routines.
 *****************************************************************/

/* Function:  p7_modelsample()
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
p7_modelsample(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[random HMM created by uniform sampling]";
  int     status;

  if (( status = modelsample_engine(r, M, abc, &hmm) ) != eslOK) goto ERROR;

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

/* Function:  p7_modelsample_Prior()
 * Synopsis:  Sample an HMM from a prior distribution.
 *
 * Purpose:   Creates a random HMM of length <M> nodes, for 
 *            alphabet <abc>, using random number generator
 *            <r>, by sampling from mixture Dirichlet priors
 *            given in <pri>.
 *            
 *            In general this should give more 'realistic' profile
 *            HMMs than <p7_modelsample()> does.  In unit testing, we
 *            tend to use <p7_modelsample()> for tests that should be
 *            able to deal with any profile, no matter how
 *            pathological. A few tests need to see typical,
 *            reasonable (high-posterior-probability) alignment paths,
 *            and we use <p7_modelsample_Prior()> in these cases.
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
p7_modelsample_Prior(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, const P7_PRIOR *pri, P7_HMM **ret_hmm)
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

/* Function:  p7_modelsample_Ungapped()
 * Synopsis:  Sample a random HMM with no nonzero indel transitions.
 *
 * Purpose:   Same as <p7_modelsample()>, except all 
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
p7_modelsample_Ungapped(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[ungapped HMM sample]";
  int     k;
  int     status;

  if ( (status = modelsample_engine(r, M, abc, &hmm)) != eslOK) goto ERROR;

  for (k = 0; k <= M; k++) {
    hmm->t[k][p7H_MM] = 1.0;
    hmm->t[k][p7H_MD] = 0.0;
    hmm->t[k][p7H_MI] = 0.0;
  }

  /* Add mandatory annotation, and some relevant optional annotation  */
  p7_hmm_SetName(hmm, "sampled-ungapped-hmm");
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


/* Function:  p7_modelsample_Uniform()
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
p7_modelsample_Uniform(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, 
		       float tmi, float tii, float tmd, float tdd,
		       P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[HMM with uniform transitions, random emissions]";
  int     k;
  int     status;

  if ((hmm = p7_hmm_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; }
  
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
/*-------------- end, basic model sampling ----------------------*/



/*****************************************************************
 * 2. Models with finite enumerable path #
 *****************************************************************/

/* Function:  p7_modelsample_Enumerable()
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
 *            Compare <p7_modelsample_Enumerable2()>, which only makes
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
p7_modelsample_Enumerable(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[random enumerable HMM created by sampling]";
  int     k;
  float   tmp[2];
  int     status;
  
  if ((hmm = p7_hmm_Create(M, abc)) == NULL)  { status = eslEMEM; goto ERROR; }

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



/* Function:  p7_modelsample_Enumerable2()
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
 *            Compare <p7_modelsample_Enumerable()>, which makes all
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
p7_modelsample_Enumerable2(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
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
/*--------- end, models with enumerable paths -------------------*/





/*****************************************************************
 * 3. Models with only a single valid path
 *****************************************************************/

/* Function:  p7_modelsample_SinglePathed()
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
p7_modelsample_SinglePathed(ESL_RANDOMNESS *r, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
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

/* Function:  p7_modelsample_SinglePathedSeq()
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
 *            such that (t_MD = 0 or t_MM,tMI = 0), and (t_DM = 0 or t_DD
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
p7_modelsample_SinglePathedSeq(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
			       P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
			       P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc)
{
  return sample_single_pathed_seq_engine(rng, M, bg, opt_hmm, opt_gm, opt_dsq, opt_L, opt_tr, opt_anch, opt_D, opt_sc, /*do_asc_version=*/FALSE);
}

/* Function:  p7_modelsample_SinglePathedASC()
 * Synopsis:  Sample profile/seq/anch triplet such that only one P=1.0 path exists for anchored seq comparison.
 *
 * Purpose:   Sample a profile, sequence, anchor set combination that have been
 *            contrived such that there is only a single possible path that
 *            can generate that sequence and satisfy the anchor set constraint.
 *            That is, $P(\pi | x, A, M) = 1$ for one and only one $\pi$.
 *            
 *            Essentially the same as <p7_modelsample_SinglePathedSeq(),
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
p7_modelsample_SinglePathedASC(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
			       P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
			       P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc)
{
  return sample_single_pathed_seq_engine(rng, M, bg, opt_hmm, opt_gm, opt_dsq, opt_L, opt_tr, opt_anch, opt_D, opt_sc, /*do_asc_version=*/TRUE);
}
/*---------- end, models w/ only one valid path -----------------*/


/*****************************************************************
 * 4. Models such that ASC paths == all paths
 *****************************************************************/

/* Function:  p7_modelsample_AnchoredUni()
 * Synopsis:  Model/seq/anchor triplet such that all paths use anchor.
 *
 * Purpose:   Sample a model, a sequence, and an anchor, contrived such
 *            that all paths for the model/sequence comparison must
 *            use the anchor set even when unconstrained by an ASC
 *            algorithm. Thus, the standard Forward, Backward, and
 *            Decoding algorithms give the same scores (and to a
 *            well-defined extent, the same DP matrix values) as ASC
 *            F/B/D. This enables unit tests of ASC algorithms,
 *            comparing to standard algorithms.
 *            
 *            <_AnchoredUni()> allows multiple paths to contribute to
 *            the Forward/Backward sums, and allows N/C/J and insert
 *            states to emit. However, the model is in uniglocal mode,
 *            so local transitions and multiple domains aren't used or
 *            tested.
 *            
 *            The key idea is to pick an anchor match state k0 that
 *            must be visited in all paths, and must generate residue
 *            X; and make a sequence with one X in it. Thus, path must
 *            assign X to Mk0, because Mk0 has to emit something in
 *            the path and that's the only residue it can.
 *            
 *            The HMM is a standard HMM sample (all
 *            transition/emission probabilities valid) except: Mk0
 *            must emit residue X with p=1, and transitions force all
 *            paths to include Mk0 by having t_{k0-1}(MD) = 0,
 *            t_{k0-1}(DM) = 1, t_{k0-1}(DD) = 0.
 *            
 *            The profile is in glocal mode, to force all paths to
 *            include the Mk0, and in unihit mode, to assure that
 *            residue X is uniquely aligned to Mk0, rather than to
 *            N/C/J/I states that can emit anything.
 *            
 *            The sequence is created such that it only contains 1 X
 *            residue, at the anchor position.
 *
 * Args:      rng      : random number generator
 *            M        : length of <hmm>/<gm> to sample
 *            bg       : background model (contains digital alphabet to use too)
 *            opt_hmm  : optRETURN: model
 *            opt_gm   : optRETURN: profile (uniglocal, L)
 *            opt_dsq  : optRETURN: digital sequence
 *            opt_L    : optRETURN: length of <dsq>
 *            opt_tr   : optRETURN: traceback, state path and emitted seq
 *            opt_anch : optRETURN: array of anchors i0,k0 for each domain d (always 1 of them)
 *            opt_D    : optRETURN: number of domains d in <anch>, <tr> (always 1)
 *            opt_sc   : optRETURN: raw lod score, in nats
 *
 * Returns:   <eslOK> on success,  and any results that were requested
 *            by passing non-<NULL> pointers for result arguments are
 *            now available.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            On any exception, any requested result pointer is returned <NULL>,
 *            and any requested return value is returned 0.0.
 *
 * Xref:      SRE:J13/51
 */
int
p7_modelsample_AnchoredUni(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
			   P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
			   P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc)
{
  char       *hmmname   = "anchored_uni";
  char       *logmsg    = "[Test model produced by p7_modelsample_AnchoredUni()]";
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  P7_TRACE   *tr        = NULL;
  ESL_SQ     *sq        = NULL;
  ESL_DSQ    *dsq       = NULL;
  P7_COORD2  *anch      = NULL;
  int         D         = 1;        // Because we're unihit, D must be 1.
  int         L;
  int         k,z;
  int         k0;                   // Anchor state Mk0; (1..M)
  ESL_DSQ     anchX;                // Residue emitted with p=1 by the anchor Mk0; (0..K-1 in bg->abc)
  float       sc;
  int         status;

  /* allocations that we know how to make already */
  if ((  gm = p7_profile_Create(M, bg->abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if ((  sq = esl_sq_CreateDigital(bg->abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if ((  tr = p7_trace_Create())             == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(anch, sizeof(P7_COORD2) * D);

  /* Start with a randomly sampled model, that we'll tweak slightly. */
  if (( status = modelsample_engine(rng, M, bg->abc, &hmm)) != eslOK) goto ERROR;

  /* Select anchor position k0 in model, 1..M */
  k0 = 1 + esl_rnd_Roll(rng, M);

  /* Tweak transitions to guarantee that Mk0 is occupied w/ P=1;
   * tweak emission so that Mk0 generates one residue w/ P=1.
   */
  hmm->t[k0-1][p7H_MM] += hmm->t[k0-1][p7H_MD];
  hmm->t[k0-1][p7H_MD]  = 0.0;
  hmm->t[k0-1][p7H_DM]  = 1.0;    // also works at k=0, where t0[DM]=1 is convention for nonexistent D0 state
  hmm->t[k0-1][p7H_DD]  = 0.0;

  anchX = esl_rnd_Roll(rng, bg->abc->K);
  esl_vec_FSet(hmm->mat[k0], bg->abc->K, 0.);
  hmm->mat[k0][anchX] = 1.0;
  for (k = 1; k <= M; k++) 
    if (k != k0)                            // avoid a pathological case where e_Mk(X) = 1 for a nonanchor position 
      while (hmm->mat[k][anchX] > 0.5)      // the 0.5 is arbitrary. It can't be 1, and we want to be able to resample a non-X residue efficiently.
	esl_dirichlet_FSampleUniform(rng, bg->abc->K, hmm->mat[k]);

  ESL_DASSERT1(( bg->abc->K  >= 2   ));
  ESL_DASSERT1(( bg->f[anchX] < 1.0 ));     // avoid a pathological case that would break a while loop below

  /* Add mandatory annotation */
  p7_hmm_SetName(hmm, hmmname);
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  p7_hmm_SetCtime(hmm);
  p7_hmm_SetConsensus(hmm, NULL);
  p7_hmm_SetComposition(hmm);

  /* Create a profile; emit a trace from it.
   * We don't yet know actual L of target seq, so here we're 
   * just choosing a smallish L suitable for sampling smallish
   * path from profile.
   */
  if (( status = p7_profile_ConfigUniglocal(gm, hmm, bg, 10)) != eslOK) goto ERROR;
  if (( status = p7_ProfileEmit(rng, hmm, gm, bg, sq, tr))    != eslOK) goto ERROR;
  if (( status = p7_trace_Index(tr))                          != eslOK) goto ERROR;

  /* Extract <dsq> from <sq> container. Doctor it so that the only X
   * residues in it were emitted by Mk0.  You can't just willy-nilly
   * change X's to any old Y at Mk != k0, in case e_k(Y) = 0; also,
   * we rely on having already made sure that there is no Mk != k0
   * that is forced to generate X.
   */
  dsq = sq->dsq;
  L   = sq->n;
  sq->dsq = NULL;
  sq->n   = 0;
  esl_sq_Destroy(sq);
  sq = NULL;

  for (z = 1; z < tr->N-1; z++)
    if (tr->i[z])
      {
	if (p7_trace_IsM(tr->st[z]) && tr->k[z] != k0)
	  while (dsq[tr->i[z]] == anchX)   // while loop must terminate, because we made sure e_Mk(y != X) has good probability
	    dsq[tr->i[z]] = esl_rnd_FChoose(rng, hmm->mat[tr->k[z]], bg->abc->K);
	else if (! p7_trace_IsM(tr->st[z]))
	  while (dsq[tr->i[z]] == anchX)   // while loop will terminate so long as bg->f distribution isn't pathological (p_X = 1), which we asserted above
	      dsq[tr->i[z]] = esl_rnd_FChoose(rng, bg->f, bg->abc->K);
      }

  /* Make the anchor "set", which we know has one anchor for 
   * one domain. We know k0, we just need to find i0.
   */
  ESL_DASSERT1(( tr->ndom == 1 ));
  anch[0].n2 = k0;
  for (z = 1; z < tr->N-1; z++)
    if (tr->k[z] == k0 && p7_trace_IsM(tr->st[z]))
      {
	anch[0].n1 = tr->i[z];
	break;
      }

  /* Finish up; set the length model, then find the trace score.
   */
  p7_profile_SetLength(gm, L);
  p7_trace_Score(tr, dsq, gm, &sc);

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


/* Function:  p7_modelsample_AnchoredLocal()
 * Synopsis:  Sample a model/seq/anchorset triplet, such that all paths use anchor.
 *
 * Purpose:   Sample a model, sequence, and anchor set, contrived such that
 *            all paths for the model/sequence comparison must use the anchor
 *            set even when unconstrained by the ASC algorithm. 
 *
 *            <_AnchoredLocal()> is unique in that it allows and tests
 *            local transitions; <_AnchoredUni()> and
 *            <_AnchoredMulti()> do not. <_AnchoredLocal()> prohibits
 *            N/C/J states from emitting residues, prohibits I states,
 *            and prohibits multiple domains. Paths essentially consist
 *            only of M/D states.
 *            
 *            The key idea is that only match states can generate
 *            residues, and only the anchor Mk0 can generate a special
 *            residue X with emission prob > 0. Thus, every X in the
 *            sequence must be assigned to the anchor. In
 *            <_AnchoredLocal()>, the target sequence contains exactly
 *            1 X residue, and therefore all paths (including local
 *            entry/exit) must assign that X to the anchor.
 *            
 *            The HMM has all t_k(MI) = 0; e_Mk0(X) > 0; and e_Mk(X) =
 *            0 for all k != k0.
 *            
 *            The profile is unihit L=0 dual-mode glocal/local. To
 *            allow local paths, it has to prohibit N/C emissions,
 *            because otherwise the X residue could be assigned to N
 *            or C, and a local path could simply avoid using the
 *            Mk0. (That is, with local paths, there is no way we can
 *            use transition probabilities to guarantee that Mk0 will
 *            be visited.) Similarly, the model must be unihit,
 *            because otherwise non-X residues can be assigned to
 *            various combinations of additional local hits.
 *            
 *            The sequence has only one X, and because only Mk0 can
 *            generate it, this is forced to be the anchor.
 *            
 * Args:      rng      : random number generator
 *            M        : length of <hmm>/<gm> to sample
 *            bg       : background model (contains digital alphabet to use too)
 *            opt_hmm  : optRETURN: model
 *            opt_gm   : optRETURN: profile (dual-mode local/glocal, L=0, unihit)
 *            opt_dsq  : optRETURN: digital sequence
 *            opt_L    : optRETURN: length of <dsq>
 *            opt_tr   : optRETURN: traceback, state path and emitted seq
 *            opt_anch : optRETURN: array of anchors i0,k0 for each domain d (always 1 of them)
 *            opt_D    : optRETURN: number of domains d in <anch>, <tr> (always 1)
 *            opt_sc   : optRETURN: raw lod score, in nats
 *
 * Returns:   <eslOK> on success,  and any results that were requested
 *            by passing non-<NULL> pointers for result arguments are
 *            now available.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            On any exception, any requested result pointer is returned <NULL>,
 *            and any requested return value is returned 0.0.
 */
int
p7_modelsample_AnchoredLocal(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
			     P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
			     P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc)
{
  return sample_anchored_engine(rng, M, bg, opt_hmm, opt_gm, opt_dsq, opt_L, opt_tr, opt_anch, opt_D, opt_sc, TRUE);
}

/* Function:  p7_modelsample_AnchoredMulti()
 * Synopsis:  Sample model/seq/anchorset triplet such that all paths use anchor
 *
 * Purpose:   Sample a model, sequence, and anchor set contrived such that all 
 *            paths for the model/sequence comparison must use the anchor set
 *            even when unconstrained by the ASC algorithm.
 *            
 *            Of the three <_Anchored*> sample strategies,
 *            <_AnchoredMulti()> is unique in that it allows and tests
 *            multiple hit mode. <_AnchoredUni()> and
 *            <_AnchoredLocal()> do not. <_AnchoredLocal()> and
 *            <_AnchoredMulti()> are similar to each other, and both
 *            prohibit N/C/J states from emitting residues and
 *            prohibit I states. Paths essentially consist only of M/D
 *            states.
 *            
 *            Like <_AnchoredLocal()>, the key idea is that only match
 *            states can generate residues, and only the anchor Mk0
 *            can generate a special residue X; here, with e_Mk0(X)
 *            strictly 1.  The transition probabilities are crafted
 *            such that Mk0 must be visited in every domain.  Thus,
 *            every X in the sequence must be assigned to the anchor
 *            Mk0, and every domain must contain an Mk0.  In
 *            <_AnchoredMulti()>, there can be more than one X, each
 *            of which corresponds to an anchor to Mk0.
 *            
 *            The HMM has all t_k(MI) = 0; e_Mk0(X) = 1; e_Mk(X) = 0
 *            for all k != k0; and Mk0 inclusion is forced with
 *            t_{k0-1}(MD) = 0, t_{k0-1}(DM) = 1, t_{k0-1}(DD) = 0.
 *            
 *            The profile is in multiglocal L=0 mode. It cannot allow
 *            local transitions while allowing multiple hits, because
 *            otherwise non-X residues could be aligned to various
 *            different combinations of local segments of the model;
 *            whereas in glocal mode, each X must be uniquely assigned
 *            to Mk0.
 *            
 *            The sequence has D X's, one for each of D anchor/domains.
 *
 * Args:      rng      : random number generator
 *            M        : length of <hmm>/<gm> to sample
 *            bg       : background model (contains digital alphabet to use too)
 *            opt_hmm  : optRETURN: model
 *            opt_gm   : optRETURN: profile (multiglocal, L=0)
 *            opt_dsq  : optRETURN: digital sequence
 *            opt_L    : optRETURN: length of <dsq>
 *            opt_tr   : optRETURN: traceback, state path and emitted seq (w/ index)
 *            opt_anch : optRETURN: array of anchors i0,k0 for each domain d 
 *            opt_D    : optRETURN: number of domains d in <anch>, <tr> 
 *            opt_sc   : optRETURN: raw lod score, in nats
 *
 * Returns:   <eslOK> on success,  and any results that were requested
 *            by passing non-<NULL> pointers for result arguments are
 *            now available.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            On any exception, any requested result pointer is returned <NULL>,
 *            and any requested return value is returned 0.0.
 */
int
p7_modelsample_AnchoredMulti(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
			     P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
			     P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc)
{
  return sample_anchored_engine(rng, M, bg, opt_hmm, opt_gm, opt_dsq, opt_L, opt_tr, opt_anch, opt_D, opt_sc, FALSE);
}





/*----------- end, models with ASC paths == all paths -----------*/



/*****************************************************************
 * 5. Internal (static) functions
 *****************************************************************/

static int
modelsample_engine(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm = NULL;
  int     k;
  int     status;

  if ( (hmm = p7_hmm_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; }

  /* _Create() has already initialized conventions of nonexistent states at k=0,
   * where M0 is treated as B.
   *   - mat[0] emission distribution is { 1.0, 0.0, ... }
   *   - t[0][DM] = 1.0 and DD is 0.0.
   * Now the rest of node 0.
   */
  esl_dirichlet_FSampleUniform(rng, abc->K,    hmm->ins[0]);
  esl_dirichlet_FSampleUniform(rng, p7H_NTMAT, P7H_TMAT(hmm, 0));
  hmm->t[0][p7H_II] = esl_random(rng) * p7H_II_SAMPLE_MAX;
  hmm->t[0][p7H_IM] = 1.0 - hmm->t[0][p7H_II];

  /* Main recursion over most nodes */
  for (k = 1; k < M; k++)
    {
      esl_dirichlet_FSampleUniform(rng, abc->K, hmm->mat[k]);
      esl_dirichlet_FSampleUniform(rng, abc->K, hmm->ins[k]);
      esl_dirichlet_FSampleUniform(rng, p7H_NTMAT, P7H_TMAT(hmm, k));
      hmm->t[k][p7H_II] = esl_random(rng) * p7H_II_SAMPLE_MAX;
      hmm->t[k][p7H_IM] = 1.0 - hmm->t[k][p7H_II];
      esl_dirichlet_FSampleUniform(rng, p7H_NTDEL, P7H_TDEL(hmm, k));
    }

  /* Node M is special: no transitions to D, transitions to M
   * are interpreted as transitions to E.
   */
  esl_dirichlet_FSampleUniform(rng, abc->K, hmm->mat[M]);
  esl_dirichlet_FSampleUniform(rng, abc->K, hmm->ins[M]);
  hmm->t[M][p7H_MM] = esl_random(rng); 
  hmm->t[M][p7H_MI] = 1.0 - hmm->t[M][p7H_MM];
  hmm->t[M][p7H_MD] = 0.0;
  hmm->t[M][p7H_II] = esl_random(rng) * p7H_II_SAMPLE_MAX;
  hmm->t[M][p7H_IM] = 1.0 - hmm->t[k][p7H_II];
  hmm->t[M][p7H_DM] = 1.0;
  hmm->t[M][p7H_DD] = 0.0;
  
  /* Do not add other annotation here. That's the caller's job. */
  
  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}


  
/* see comments on p7_modelsample_SinglePathedSeq() and p7_modelsample_SinglePathedASC(), below;
 * this is the engine they both share.
 */
static int
sample_single_pathed_seq_engine(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
				P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
				P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc,
				int do_asc_version)
{
  char       *hmmname   = (do_asc_version ? "single_pathed_asc" : "single_pathed_seq");
  char       *logmsg    = (do_asc_version ? "[Test model produced by p7_modelsample_SinglePathedASC()]" : "[Test model produced by p7_modelsample_SinglePathedSeq()]");
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
  sq = NULL;

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
      while (tr->st[z] != p7T_B) z++;
      s = esl_rnd_Roll(rng, nm);   // use the s'th match state as anchor
      while (s || ! p7_trace_IsM(tr->st[z])) {
	if (p7_trace_IsM(tr->st[z])) s--;
	z++;
      }
      anch[d].n1 = tr->i[z];
      anch[d].n2 = tr->k[z];
    }
  
  /* Configure final length model of <gm> before returning */
  p7_profile_SetLength(gm, L);
  
  /* Calculate the score of the trace.
   * This has to come after length model config.
   */
  p7_trace_Score(tr, dsq, gm, &sc);

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


/* Shared by p7_modelsample_AnchoredLocal(), p7_modelsample_AnchoredMulti().
 * See documentation on them, above.
 */
static int
sample_anchored_engine(ESL_RANDOMNESS *rng, int M, const P7_BG *bg,
		       P7_HMM **opt_hmm, P7_PROFILE **opt_gm, ESL_DSQ **opt_dsq, int *opt_L, 
		       P7_TRACE **opt_tr, P7_COORD2 **opt_anch, int *opt_D, float *opt_sc,
		       int do_local_version)
{
  char       *hmmname   = (do_local_version ? "anchored_local" : "anchored_multihit");
  char       *logmsg    = (do_local_version ? "[Test model produced by p7_modelsample_AnchoredLocal()]" : "[Test model produced by p7_modelsample_AnchoredMulti()]");
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  P7_TRACE   *tr        = NULL;
  ESL_SQ     *sq        = NULL;
  ESL_DSQ    *dsq       = NULL;
  P7_COORD2  *anch      = NULL;
  int         D,d;
  int         L;
  int         k,z;
  int         k0;                   // Anchor state Mk0; (1..M)
  ESL_DSQ     anchX;                // Residue emitted with p=1 by the anchor Mk0; (0..K-1 in bg->abc)
  int         nanch;
  float       sc;
  int         status;

  /* allocations that we know how to make already */
  if ((  gm = p7_profile_Create(M, bg->abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if ((  sq = esl_sq_CreateDigital(bg->abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if ((  tr = p7_trace_Create())             == NULL) { status = eslEMEM; goto ERROR; }

  /* Select anchor position k0 in model (1..M), and
   * select residue that only the anchor Mk0 will be able to emit with e_Mk0>0.
   */
  k0    = 1 + esl_rnd_Roll(rng, M);
  anchX = esl_rnd_Roll(rng, bg->abc->K);

  /* Sample an HMM in the usual way, which we'll then tweak as we need. 
   * It can have arbitrary transitions into node k0 -- even tXMk0 = 0 is ok,
   * because Mk0 can still be reached by a local start.
   */
  if (( status = modelsample_engine(rng, M, bg->abc, &hmm)) != eslOK) goto ERROR;

  /* Overwrite some transitions in sampled HMM now, so all tkMI = 0 */
  for (k = 0; k < M; k++)
    {
      hmm->t[k][p7H_MM] = esl_random(rng);
      hmm->t[k][p7H_MI] = 0.0;
      hmm->t[k][p7H_MD] = 1.0 - hmm->t[k][p7H_MM];
    }
  hmm->t[M][p7H_MM] = 1.0;
  hmm->t[M][p7H_MI] = 0.0;
  hmm->t[M][p7H_MD] = 0.0;

  /* In multihit version, we must force inclusion of Mk0 in every glocal domain */
  if (! do_local_version) 
    {
      hmm->t[k0-1][p7H_MM] += hmm->t[k0-1][p7H_MD];
      hmm->t[k0-1][p7H_MD]  = 0.0;
      hmm->t[k0-1][p7H_DM]  = 1.0;    // also works at k0=1, k=0, where t0[DM]=1 is convention for nonexistent D0 state
      hmm->t[k0-1][p7H_DD]  = 0.0;
    }

  /* Overwrite match emissions;
   * Only Mk0 can emit anch0 with e>0. No other Mk can.  
   * In multihit version, eMk0(X) = 1, not just >0.
   */
  if (do_local_version) { 
    while (hmm->mat[k0][anchX] == 0.0f)  esl_dirichlet_FSampleUniform(rng, bg->abc->K, hmm->mat[k0]);	  
  } else {
    esl_vec_FSet(hmm->mat[k0], bg->abc->K, 0.);
    hmm->mat[k0][anchX] = 1.0;
  }

  for (k = 1; k <= M; k++) 
    if (k != k0) 
      {
	hmm->mat[k][bg->abc->K-1] = 0.0;
	esl_dirichlet_FSampleUniform(rng, bg->abc->K-1, hmm->mat[k]);	  
	hmm->mat[k][bg->abc->K-1] = hmm->mat[k][anchX];
	hmm->mat[k][anchX]        = 0.0;
      }

  /* Add mandatory annotation */
  p7_hmm_SetName(hmm, hmmname);
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  p7_hmm_SetCtime(hmm);
  p7_hmm_SetConsensus(hmm, NULL);
  p7_hmm_SetComposition(hmm);

  /* Create a profile in L=0 mode.
   * Params are:  
   */
  if (do_local_version) { /* L=0: no N/C/J emission; nj=0.0: unihit;  pglocal=0.5: dual-mode glocal/local */
    if (( status = p7_profile_ConfigCustom(gm, hmm, bg, 0, 0.0, 0.5))  != eslOK) goto ERROR;  
  } else {                /* multihit glocal */
    if (( status = p7_profile_ConfigGlocal(gm, hmm, bg, 0))            != eslOK) goto ERROR;
  }

  /* Emit a trace from it.  When the profile is in dual mode, it can
   * generate local alignments that don't include the anchor Mk0.
   * Also, glocal or local paths can also pass thru Dk0 instead of
   * Mk0.  Reject such traces until we get one in which all domains
   * pass thru the anchor Mk0.
   */
  do {
    if (( status = p7_ProfileEmit(rng, hmm, gm, bg, sq, tr)) != eslOK) goto ERROR;
    if (( status = p7_trace_Index(tr))                       != eslOK) goto ERROR;
    nanch = 0;
    for (z = 0; z < tr->N; z++)
      if (tr->k[z] == k0 && p7_trace_IsM(tr->st[z])) nanch++;
  } while (nanch != tr->ndom);
  /* Note that this condition required that all the anchor k0's are M states */

  /* Extract <dsq> from <sq> container.
   */
  dsq = sq->dsq;
  L   = sq->n;
  sq->dsq = NULL;
  sq->n   = 0;
  esl_sq_Destroy(sq);
  sq = NULL;

  /* Make the anchor set, which we know has one anchor for 
   * one domain. We know k0, we just need to find i0's.
   * On every anchor, make the dsq[] residue == anchX.
   */
  D = tr->ndom;  // ok because we indexed the trace
  ESL_ALLOC(anch, sizeof(P7_COORD2) * D);

  d = 0;
  for (z = 1; z < tr->N-1; z++)
    if (tr->k[z] == k0 && p7_trace_IsM(tr->st[z]))
      {
	anch[d].n1    = tr->i[z];
	anch[d].n2    = k0;
	dsq[tr->i[z]] = anchX;
	d++;
      }

  /* Finish up by finding the trace score.
   */
  p7_trace_Score(tr, dsq, gm, &sc);

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



/*------------------- end, internal functions -------------------*/



/*****************************************************************
 * 6. Unit tests
 *****************************************************************/
#ifdef p7MODELSAMPLE_TESTDRIVE

static void
utest_modelsample(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char    failmsg[] = "modelsample::p7_modelsample() unit test failed";
  P7_HMM *hmm       = NULL;
  int     ntrials   = 10;
  char    errmsg[eslERRBUFSIZE];
  int     i;
  int     status;
  
  for (i = 0; i < ntrials; i++)
    {
      if (( status = p7_modelsample(rng, M, abc, &hmm))    != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001)) != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      p7_hmm_Destroy(hmm);
    }
}

static void
utest_prior(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char      failmsg[] = "modelsample::p7_modelsample_Prior() unit test failed";
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
      if (( status = p7_modelsample_Prior(rng, M, abc, pri, &hmm)) != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001))         != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      p7_hmm_Destroy(hmm);
    }

  p7_prior_Destroy(pri);
}

static void
utest_ungapped(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char      failmsg[] = "modelsample::p7_modelsample_Ungapped() unit test failed";
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
      if (( status = p7_modelsample_Ungapped(rng, M, abc, &hmm))  != eslOK) esl_fatal(failmsg);
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
utest_enumerable(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char      failmsg[] = "modelsample::p7_modelsample_Enumerable() unit test failed";
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
      if (( status = p7_modelsample_Enumerable(rng, M, abc, &hmm)) != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001))         != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      if (( status = p7_CoreEmit(rng, hmm, sq, tr))                != eslOK) esl_fatal(failmsg);
      if (( status = p7_trace_Validate(tr, abc, sq->dsq, errmsg))  != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      
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
utest_enumerable2(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char      failmsg[] = "modelsample::p7_modelsample__Enumerable2() unit test failed";
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
      if (( status = p7_modelsample_Enumerable2(rng, M, abc, &hmm)) != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001))          != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      if (( status = p7_CoreEmit(rng, hmm, sq, tr))                 != eslOK) esl_fatal(failmsg);
      if (( status = p7_trace_Validate(tr, abc, sq->dsq, errmsg))   != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      
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
utest_uniform(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char    failmsg[] = "modelsample::p7_modelsample_Uniform() unit test failed";
  P7_HMM *hmm       = NULL;
  int     ntrials   = 10;
  char    errmsg[eslERRBUFSIZE];
  int     i;
  int     status;
  
  for (i = 0; i < ntrials; i++)
    {
      if (( status = p7_modelsample_Uniform(rng, M, abc, 0.1, 0.4, 0.1, 0.4, &hmm)) != eslOK) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate(hmm, errmsg, 0.0001))                          != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      p7_hmm_Destroy(hmm);
    }
}

static void
utest_singlepathed(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "modelsample::p7_modelsample_SinglePathed() unit test failed";
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
      if (( status = p7_modelsample_SinglePathed(rng, M, abc, &hmm))    != eslOK) esl_fatal(failmsg);
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
utest_singlepathed_seq(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int do_asc_version)
{
  char        failmsg[] = "modelsample::p7_modelsample_SinglePathedSeq() unit test failed";
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
      if (do_asc_version) status = p7_modelsample_SinglePathedASC(rng, M, bg, &hmm, &gm, &dsq1, &L1, &tr1, &anch, &D1, &sc1);
      else                status = p7_modelsample_SinglePathedSeq(rng, M, bg, &hmm, &gm, &dsq1, &L1, &tr1, &anch, &D1, &sc1);
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

static void
utest_anchored_uni(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "modelsample::p7_modelsample_AnchoredUni() unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL;
  ESL_SQ     *sq        = esl_sq_CreateDigital(abc);
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_COORD2  *anch      = NULL;
  int         D;
  float       sc;
  int         nhmm      = 10;
  int         ntrace    = 10;
  int         h,t,z,i;
  int         anchX     = -1;
  int         nvisit;
  char        errmsg[eslERRBUFSIZE];
  int         status;

  if (bg == NULL || sq == NULL) esl_fatal(failmsg);

  for (h = 0; h < nhmm; h++)
    {
      if (( status = p7_modelsample_AnchoredUni(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

      if (D != 1) esl_fatal(failmsg);
      if (( status = p7_hmm_Validate    (hmm, errmsg, 0.0001))   != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      if (( status = p7_profile_Validate(gm,  errmsg, 0.0001))   != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      if (( status = p7_trace_Validate  (tr, abc, dsq, errmsg))  != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      /* We don't directly know what the magic anchor residue is, so 
       * figure that out first.
       */
      for (z = 1; z < tr->N-1; z++)
	if (tr->k[z] == anch[0].n2 && p7_trace_IsM(tr->st[z])) { anchX = dsq[tr->i[z]]; break; }

      if (hmm->mat[anch[0].n2][anchX] != 1.0)   esl_fatal(failmsg);
      for (i = 1; i <= L; i++)
	if (dsq[i] == anchX && i != anch[0].n1) esl_fatal(failmsg);  // Check that dsq has only one X, at the anchor position.

      /* Emit a few traces from the profile.
       * Verify that they all visit the anchor state Mk0, and that the anchor residue is X.
       */
      p7_trace_Reuse(tr);
      for (t = 0; t < ntrace; t++)
	{
	  if (( status = p7_ProfileEmit(rng, hmm, gm, bg, sq, tr))    != eslOK) esl_fatal(failmsg);
	  if (( status = p7_trace_Validate(tr, abc, sq->dsq, errmsg)) != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
	  
	  nvisit = 0;
	  for (z = 0; z < tr->N; z++)
	    if (tr->k[z] == anch[0].n2 && p7_trace_IsM(tr->st[z]))
	      {
		nvisit++;          
		if (sq->dsq[tr->i[z]] != anchX) esl_fatal(failmsg);  // Anchor i0 residue is an X. (In emitted sequences, there can be more than one X)
	      } 
	  if (nvisit != 1) esl_fatal(failmsg);                   // Check that all traces visit anchor i0,k0,M once.

	  esl_sq_Reuse(sq);
	  p7_trace_Reuse(tr);
	}

      free(dsq);
      free(anch);      
      p7_trace_Destroy(tr);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);
    }

  esl_sq_Destroy(sq);
  p7_bg_Destroy(bg);
}

static void
utest_anchored_local(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "modelsample::p7_modelsample_AnchoredLocal() unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL;
  ESL_SQ     *sq        = esl_sq_CreateDigital(abc);
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_COORD2  *anch      = NULL;
  int         D;
  float       sc;
  int         nhmm      = 10;
  int         h,z,i,d;
  int         anchX     = -1;
  int         k,k0;
  char        errmsg[eslERRBUFSIZE];
  int         status;

  if (bg == NULL || sq == NULL) esl_fatal(failmsg);

  for (h = 0; h < nhmm; h++)
    {
      if (( status = p7_modelsample_AnchoredLocal(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

      if (( status = p7_hmm_Validate    (hmm, errmsg, 0.0001))   != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      if (( status = p7_profile_Validate(gm,  errmsg, 0.0001))   != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      if (( status = p7_trace_Validate  (tr, abc, dsq, errmsg))  != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      if (D        != 1) esl_fatal(failmsg);
      if (tr->ndom != 1) esl_fatal(failmsg);

      /* We don't directly know what the magic anchor residue is, so figure that out first. */
      k0 = anch[0].n2;
      for (z = 1; z < tr->N-1; z++)
	if (tr->k[z] == k0 && p7_trace_IsM(tr->st[z])) { anchX = dsq[tr->i[z]]; break; }

      /* The emission prob of residue anchX at anchor state k0 should be >0,
       * and the only anchX's in the sequence come from Mk0
       */
      if (hmm->mat[k0][anchX] == 0.) esl_fatal(failmsg);
      for (k = 1; k <= hmm->M; k++)
	if (k != k0 &&  hmm->mat[k][anchX] != 0.)  esl_fatal(failmsg);

      /* All t_k(MI) are 0 */
      for (k = 0; k <= hmm->M; k++)
	if (hmm->t[k][p7H_MI] != 0.) esl_fatal(failmsg);

      /* The dsq has only one X, at the anchor position.
       */
      d = 0;
      for (i = 1; i <= L; i++)
	if (dsq[i] == anchX) {
	  if (i != anch[d++].n1) esl_fatal(failmsg); 
	}
      if (d != 1) esl_fatal(failmsg);

      /* The trace may not contain emitting N/C/J, nor any I */
      for (z = 1; z < tr->N; z++)
	{
	  if (p7_trace_IsI(tr->st[z])) esl_fatal(failmsg);
	  if (tr->st[z] == p7T_N && tr->st[z-1] == p7T_N) esl_fatal(failmsg);
	  if (tr->st[z] == p7T_J && tr->st[z-1] == p7T_J) esl_fatal(failmsg);
	  if (tr->st[z] == p7T_C && tr->st[z-1] == p7T_C) esl_fatal(failmsg); 
	}

      free(dsq);
      free(anch);      
      p7_trace_Destroy(tr);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);
    }

  p7_bg_Destroy(bg);
}

static void
utest_anchored_multi(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "modelsample::p7_modelsample_AnchoredMulti() unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL;
  ESL_SQ     *sq        = esl_sq_CreateDigital(abc);
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_COORD2  *anch      = NULL;
  int         D;
  float       sc;
  int         nhmm      = 10;
  int         h,z,i,d;
  int         anchX     = -1;
  int         k,k0;
  char        errmsg[eslERRBUFSIZE];
  int         status;

  if (bg == NULL || sq == NULL) esl_fatal(failmsg);

  for (h = 0; h < nhmm; h++)
    {
      if (( status = p7_modelsample_AnchoredMulti(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

      if (( status = p7_hmm_Validate    (hmm, errmsg, 0.0001))   != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      if (( status = p7_profile_Validate(gm,  errmsg, 0.0001))   != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);
      if (( status = p7_trace_Validate  (tr, abc, dsq, errmsg))  != eslOK) esl_fatal("%s\n  %s", failmsg, errmsg);

      if (D != tr->ndom) esl_fatal(failmsg);

      /* We don't directly know what the magic anchor residue is, so figure that out first. */
      k0 = anch[0].n2;
      for (z = 1; z < tr->N-1; z++)
	if (tr->k[z] == k0 && p7_trace_IsM(tr->st[z])) { anchX = dsq[tr->i[z]]; break; }

      /* The emission prob of residue anchX at anchor state k0 should be >0,
       * and the only anchX's in the sequence come from Mk0
       */
      if (hmm->mat[k0][anchX] == 0.) esl_fatal(failmsg);
      for (k = 1; k <= hmm->M; k++)
	if (k != k0 &&  hmm->mat[k][anchX] != 0.)  esl_fatal(failmsg);

      /* All t_k(MI) are 0 */
      for (k = 0; k <= hmm->M; k++)
	if (hmm->t[k][p7H_MI] != 0.) esl_fatal(failmsg);

      /* The dsq has D X's, one at each anchor position.
       */
      d = 0;
      for (i = 1; i <= L; i++)
	if (dsq[i] == anchX) {
	  if (i != anch[d++].n1) esl_fatal(failmsg); 
	}
      if (d != D) esl_fatal(failmsg);

      /* The trace may not contain emitting N/C/J, nor any I */
      for (z = 1; z < tr->N; z++)
	{
	  if (p7_trace_IsI(tr->st[z])) esl_fatal(failmsg);
	  if (tr->st[z] == p7T_N && tr->st[z-1] == p7T_N) esl_fatal(failmsg);
	  if (tr->st[z] == p7T_J && tr->st[z-1] == p7T_J) esl_fatal(failmsg);
	  if (tr->st[z] == p7T_C && tr->st[z-1] == p7T_C) esl_fatal(failmsg); 
	}

      free(dsq);
      free(anch);      
      p7_trace_Destroy(tr);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);
    }

  p7_bg_Destroy(bg);
}


#endif /*p7MODELSAMPLE_TESTDRIVE*/
/*---------------------- end of unit tests -----------------------*/


/*****************************************************************
 * 7. Test driver
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

  utest_modelsample     (rng, M, abc);
  utest_prior           (rng, M, abc);
  utest_ungapped        (rng, M, abc);
  utest_enumerable      (rng, M, abc);
  utest_enumerable2     (rng, M, abc);
  utest_uniform         (rng, M, abc);
  utest_singlepathed    (rng, M, abc);
  utest_singlepathed_seq(rng, M, abc, FALSE); /* tests p7_modelsample_SinglePathedSeq(), which uses uniglocal model   */
  utest_singlepathed_seq(rng, M, abc, TRUE);  /* tests p7_modelsample_SinglePathedASC(), which uses multiglocal model */
  utest_anchored_uni    (rng, M, abc);
  utest_anchored_local  (rng, M, abc);
  utest_anchored_multi  (rng, M, abc);

  fprintf(stderr, "#  status = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  exit(0); /* success */
}

#endif /*p7MODELSAMPLE_TESTDRIVE*/
/*-------------------- end of test driver ---------------------*/




/*****************************************************************
 * 8. Example driver
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
  
  //p7_modelsample_SinglePathedSeq(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc);
  //p7_modelsample_AnchoredUni(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc);
  p7_modelsample_AnchoredMulti(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc);

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
