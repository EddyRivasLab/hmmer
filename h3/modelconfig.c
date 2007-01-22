/* Model configuration: 
 * Converting a core model to a fully configured Plan7 search profile.
 * 
 * Contents:
 *     1. Routines in the exposed API.
 *     2. The four config_*() functions for specific algorithm modes.
 *     3. Wing retraction routines.
 * 
 * Revised May 2005: xref STL9/77-81.  (Uniform fragment length distribution)
 * Again, Sept 2005: xref STL10/24-26. (Inherent target length dependency)
 * Again, Jan 2007:  xref STL11/125.   (HMMER3)
 *
 * SRE, Mon May  2 10:55:16 2005 [St. Louis]
 * SRE, Fri Jan 12 08:06:33 2007 [Janelia] [Kate Bush, Aerial]
 * SVN $Id$
 */
#include "p7_config.h"

#include <math.h>
#include <float.h>

#include "easel.h"

#include "p7_hmm.h"
#include "p7_profile.h"

/*----------------------------------------------------------------------
 * Preamble.
 * 
 * There are four search modes:
 *                  single-hit              multi-hit
 *              --------------------  ------------------------
 *     local  |   sw (p7_UNILOCAL)          fs (p7_LOCAL)
 *    glocal  |    s (p7_UNIGLOCAL)         ls (p7_GLOCAL)
 *
 * Additionally, each search mode is configured for a particular
 * target length. Thus "LS/400" means a model configured for glocal,
 * multihit alignment of a target sequence of length 400.
 *
 *-----------------------------------------------------------------------
 * Exegesis. 
 * 
 * When you enter this module, you've got an HMM (P7_HMM) in "core"
 * probability form: t[], mat[], ins[] are all valid, normalized
 * probabilities. The routines here are used to create the "profile"
 * form (P7_PROFILE) of the model: tsc[], msc[], isc[], bsc[], esc[],
 * and xsc[] fields as integer log-odds scores.
 * 
 * Also in the process, xt[] are set to their algorithm-dependent
 * probabilities, though these probabilities are only for reference.
 * 
 * The configuration process breaks down into distinct conceptual steps:
 * 
 * 1. Algorithm configuration.
 *    An "algorithm mode" is chosen. This determines whether
 *    alignments will allow local entry/exit in the model, and sets
 *    the probabilities in xt[XTE], which determine
 *    multi-hit/single-hit behavior.  The "nj" value of the HMM is
 *    also set here (the expected # of times the J state will be used;
 *    0 for single-hit mode and 1 for the default parameterization of
 *    multihit modes).
 *    
 * 2. Wing retraction.
 *    In a profile, the D_1 and D_M states of the core model are
 *    removed. The probability of the paths B->D1...->Mk ("BMk") that
 *    enter D1 and use all D's before reaching M_k is treated instead
 *    as an additional dollop of B->Mk entry probability, and the
 *    probability of paths Mk->Dk+1...D_M->E ("MkE") is treated
 *    instead as an additional dollop of Mk->E exit probability.  The
 *    MkE path probability is subtracted from the Mk->Dk+1 transition.
 *    
 *    In local algorithm modes, these extra dollops are ignored, and
 *    the model is renormalized appropriately. That is, the algorithm
 *    overrides all B->DDDD->M and/or M->DDDD->E path probabilities
 *    with its own internal entry/exit probabilities.
 *    
 *    If the algorithm mode is "global" at either entry or exit, then
 *    the internal entries are set to BMk and internal exits are set
 *    to MkE, and the model is renormalized appropriately.  That is,
 *    the algorithm treats B->DDDD->M and/or M->DDDD->E path
 *    probabilities as internal entries/exits, instead of allowing
 *    dynamic programming algorithms to use the D_1 or D_M states.
 *    
 *    These two alternatives are represented differently in traces,
 *    where an X state is used to signal 'missing data' in a local
 *    alignment. Thus B->X->Mk indicates local entry, whereas B->Mk in
 *    a trace indicates a wing-retracted B->DDD->Mk entry with respect
 *    to the core HMM; similarly Mk->X->E indicates local exit, and
 *    Mk->E indicates a Mk->DDDD->E path in the core HMM.
 *    
 *    Wing retraction is a compulsive detail with two purposes. First,
 *    it removes a mute cycle from the model, B->D1 ...D_M->E, which
 *    cannot be correctly and efficiently dealt with by DP
 *    recursions. (A DP algorithm could just *ignore* that path
 *    though, and ignore the negligible amount of probability in it.)
 *    Second, wing retraction reconciles the algorithm-dependent
 *    entry/exit probabilities with the core model. For algorithms
 *    that impose local internal entry/exit, we don't want there to be
 *    any additional probability coming from "internal" B->DDD->M and
 *    M->DDD->E paths, so wing retraction takes it away.
 *  
 *  3. Local alignment D-path leveling.
 *    For fully local alignments, we want every fragment ij (starting
 *    at match i, ending from match j) to be equiprobable. There are
 *    M(M+1)/2 possible such fragments, so the probability of each
 *    one is 2/M(M+1). 
 *    
 *    Notionally, we imagine a "model" consisting of the M(M+1)/2
 *    possible fragments, with entry probability of 2/M(M+1) for each.
 *    
 *    Operationally, we achieve this by a trick inspired by a
 *    suggestion from Bill Bruno. Bill suggested that for a model with
 *    no delete states, if we set begin[k] = 1/(M-k+1) and end[k] =
 *    (M-k+1) / [M(M+1)/2], all fragments are equiprobable: the prob
 *    of any given fragment is
 *         b_i * e_j * \prod_{k=i}^{j-1} (1-e_k);
 *    that is, the fragment also includes (j-i) penalizing terms for
 *    *not* ending at i..j-1. Remarkably, this gives the result we
 *    want: this product is always 2/M(M+1), for any ij.
 *    
 *    However, D->D transitions throw a wrench into this trick,
 *    though. A local alignment that goes M_i->D...D->M_j, for
 *    example, only gets hit with one not-end penalty (for the
 *    M_i->D). This means that paths including deletions will be
 *    artifactually favored.
 *    
 *    A solution is to subtract log(1-e_k) from the deletion
 *    transition scores as well as the match transition scores.  Thus
 *    one log(1-e_k) penalty is always exacted upon transitioning from
 *    any node k->k+1. This is *not* part of the probabilistic model:
 *    it is a score accounting trick that forces the DP algorithms to
 *    associate a log(1-e_k) penalty for each node k->k+1 transition,
 *    which makes the DP calculations give the result desired for our
 *    *notional* probabilistic model with a single 2/M(M+1) transition
 *    for each possible fragment. (A similar accounting trick is the
 *    use of log-odds scoring, where we associate null model
 *    transitions and emissions with appropriate terms in the HMM, to
 *    assure that the final score of any path accounts for all the
 *    desired probability terms in an overall log-odds score). The
 *    overall score of any fragment can be rearranged such that there
 *    is one term consisting of a product of all these penalties * b_i
 *    * e_j = 2/M(M+1), and another term consisting of the actual
 *    model transition path score between i,j.
 *    
 * 4. Target length dependence. 
 *    Given a particular target sequence of length L, we want our HMM score
 *    to be as independent as possible of L. Otherwise, long sequences will
 *    give higher scores, even if they are nonhomologous. 
 *    
 *    The traditional solution to this is Karlin/Altschul statistics,
 *    which tells us that E(s=x) = KMNe^-{\lambda x}, so we expect to
 *    have to make a -1 bit score correction for every 2x increase in
 *    target sequence length (ignoring edge correction effects). K/A
 *    statistics have been proven for local Viterbi single-hit
 *    ungapped alignments. There is abundant literature showing they
 *    hold empirically for local Viterbi single-hit gapped
 *    alignments. In my hands the length dependence (though not the
 *    form of the distribution) holds for any single-hit alignment
 *    (local or glocal, Viterbi or forward) but it does not
 *    hold for multihit alignment modes.
 *    
 *    HMMER's solution is to build the length dependence right into
 *    the probabilistic model, so that we have a full probabilistic
 *    model of the target sequence. We match the expected lengths of
 *    the model M and the null model R by setting the p1, N, C, and J
 *    transitions appropriately. R has to emit the whole sequence, so
 *    it has a self-transition of L/(L+1). N, C, and J have to emit
 *    (L-(k+1)x) residues of the sequence, where x is the expected
 *    length of an alignment to the core model, and k is the expected
 *    number of times that we cycle through the J state. k=0 in sw
 *    mode, and k=1 in fs/ls mode w/ the standard [XTE][LOOP]
 *    probability of 0.5.
 *
 * 5. Conversion of probabilities to integer log-odds scores.
 *    This step incorporates the contribution of the null model,
 *    and converts floating-point probs to the scaled integer log-odds
 *    score values that are used by the DP alignment routines. 
 *
 * Step 1 is done by the main p7_ProfileConfig() function, which takes
 * a choice of algorithm mode as an argument.
 *
 * Step 2 is done by the *wing_retraction*() functions, which also
 *  go ahead and convert the affected transitions to log-odds scores;
 *  left wing retraction sets bsc[], right wing retraction sets
 *  esc[] and tsc[TM*].
 *  
 * Step 3 is carried out by one of two delete path accounting routines,
 *  which go ahead and set tsc[TD*].
 *  
 * Step 4 is carried out by the p7_ReconfigLength() routine.
 * 
 * Step 5 is carried out for all remaining scores by logoddsify_the_rest().   
 * 
 * Note that the profile never exists in a configured probability
 * form. The probability model for the search profile is implicit, not
 * explicit, because of the handling of local entry/exit transitions.
 * You can see this in more detail in emit.c:p7_ProfileEmit()
 * function, which samples sequences from the profile's probabilistic
 * model.
 *
 * So, overall, to find where the various scores and probs are set:
 *   bsc      :  wing retraction          (section 2)
 *   esc      :  wing retraction          (section 2)
 *   tsc[TM*] :  wing retraction          (section 2)
 *   tsc[TI*] :  logoddsify_the_rest()    (section 4)
 *   tsc[TD*] :  dpath leveling           (section 3)
 *   p1       :  target_ldependence()     (section 4)  
 *   xt[NCJ]  :  target_ldependence()     (section 4)  
 *   xsc (all):  logoddsify_the_rest()    (section 4)
 *   msc      :  logoddsify_the_rest()    (section 5)
 *   isc      :  logoddsify_the_rest()    (section 5)
 */



static int config_fs(P7_HMM *hmm, P7_PROFILE *gm);


/*****************************************************************
 * 1. Routines in the exposed API.
 *****************************************************************/
 
/* Function:  p7_ProfileConfig()
 * Incept:    SRE, Sun Sep 25 12:21:25 2005 [St. Louis]
 *
 * Purpose:   Given a model <hmm> with core probabilities set and null
 *            model probabilities set, and a desired <mode> (one of
 *            <p7_LOCAL>, <p7_GLOCAL>, <p7_UNILOCAL>, or <p7_UNIGLOCAL>);
 *            configure the profile <gm> into the appropriate search form for
 *            that algorithm mode. 
 *
 *            Often <gm> will be the one the <hmm> holds a reference
 *            to: that is, <p7_ProfileConfig(hmm, hmm->gm...)>.
 *            
 *            The model is configured for a default target length of
 *            350. It needs to be set to the actual length of each
 *            target sequence by a call to <p7_ReconfigLength()>.
 *            
 *            If necessary (for numerical reasons), the <p7_LCORRECT>
 *            flag will be raised on the model. This indicates that we
 *            lack sufficient numeric precision to represent transition scores
 *            for the unaligned residues, so their contribution (a total
 *            of ~1 bit for single-hit mode, ~2 bits for default multihit 
 *            mode) must instead be added post hoc to a sequence score.
 *            This score correction is calculated as needed by a call to
 *            <p7_ScoreCorrection()>.
 */
int
p7_ProfileConfig(P7_HMM *hmm, P7_PROFILE *gm, int mode)
{
  int status;

  switch (mode) {
  case p7_LOCAL:      status = config_fs(hmm, gm); break;
  case p7_GLOCAL:     status = config_ls(hmm, gm); break;
  case p7_UNILOCAL:   status = config_sw(hmm, gm); break;
  case p7_UNIGLOCAL:  status = config_s(hmm, gm);  break;
  default:            ESL_XEXCEPTION(eslESYNTAX, "no such mode")
  }
  return status;
 ERROR:
  return status;
}


/*****************************************************************
 * 2. The four config_*() functions for specific algorithm modes.
 *****************************************************************/

/*****************************************************************
 * Exegesis.
 *
 * The following functions are the Plan7 equivalent of choosing
 * different alignment styles (fully local, fully global,
 * global/local, multihit, etc.)
 * 
 * When you come into a configuration routine, the following
 * probabilities are valid in the model:
 *    1. t[1..M-1][0..6]: all the state transitions.
 *       (Node M is special: it has only a match and a delete state,
 *       no insert state, and M_M->E = 1.0 and D_M->E = 1.0 by def'n.)
 *    2. mat[1..M][]:  all the match emissions.
 *    3. ins[1..M-1][]: all the insert emissions. Note that there is
 *       no insert state in node M.
 *    4. tbd1: the B->D1 probability. The B->M1 probability is 1-tbd1.
 * These are the "data-dependent" probabilities in the model.
 * 
 * The configuration routine gets to set the "algorithm-dependent"
 * probabilities:
 *    1. xt[XTN][MOVE,LOOP] dist controls unaligned N-terminal seq.
 *       The higher xt[XTN][LOOP] is, the more unaligned seq we allow.
 *       Similarly, xt[XTC][MOVE,LOOP] dist controls unaligned C-terminal 
 *       seq, and xt[XTJ][MOVE,LOOP] dist controls length of unaligned sequence
 *       between multiple copies of a domain. Normally, if these are nonzero,
 *       they are all set to be equal to hmm->p1, the loop probability
 *       for the null hypothesis (see below).
 *    2. xt[XTE][MOVE,LOOP] distribution controls multihits. 
 *       Setting xt[XTE][LOOP] to 0.0 forces one hit per model.
 *    3. begin[1..M] controls entry probabilities. An algorithm 
 *       mode either imposes internal begin probabilities, or leaves begin[1] 
 *       as 1.0 and begin[k] = 0.0 for k>1.
 *    4. end[1..M] controls exit probabilities. An algorithm mode either
 *       imposes internal exit probabilities, or leaves end[M] = 1.0
 *       and end[k] = 0.0 for k<M.
 *    
 * The configuration routine then calls routines as appropriate to set
 * up all the model's scores, given these configured probabilities. When
 * the config routine returns, all scores are ready for alignment:
 * bsc, esc, tsc, msc, isc, and xsc.
 * 
 *****************************************************************
 *
 * SRE: REVISIT THE ISSUE BELOW. THE CONDITIONS ARE NO LONGER MET!
 *
 * There is (at least) one more issue worth noting.
 * If you want per-domain scores to sum up to per-sequence scores, which is
 * generally desirable if you don't want "bug" reports from vigilant users,
 * then one of the following two sets of conditions must be met:
 *   
 *   1) t(E->J) = 0    
 *      e.g. no multidomain hits
 *      
 *   2) t(N->N) = t(C->C) = t(J->J) = hmm->p1 
 *      e.g. unmatching sequence scores zero, and 
 *      N->B first-model score is equal to J->B another-model score.
 *      
 * These constraints are obeyed in the default Config() functions below,
 * but in the future (say, when HMM editing may be allowed) we'll have
 * to remember this. Non-equality of the summed domain scores and
 * the total sequence score is a really easy "red flag" for people to
 * notice and report as a bug, even if it may make probabilistic
 * sense not to meet either constraint for certain modeling problems.
 *****************************************************************
 */


/* config_fs()
 * Incept:   SRE, Fri Jan  2 15:34:40 1998 [StL]
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a Plan7 model to hmmfs (multihit Smith/Waterman) configuration.
 *           
 * Args:     hmm    - the Plan7 model w/ data-dep prob's valid
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 * 
 * Xref:     STL11/125: refactored for HMMER3.
 */
static int
config_fs(P7_HMM *hmm, P7_PROFILE *gm)
{
  int   k;			/* counter over states      */
  int   status;

  gm->xt[p7_XTE][p7_MOVE] = 0.5;  /* allow loops / multihits   */
  gm->xt[p7_XTE][p7_LOOP] = 0.5;

  /* Configure entry:   (M-k+1) / [M(M+1)/2]   (xref STL9/77)
   * (tbd1 is ignored)
   */
  for (k = 1; k <= hmm->M; k++)
    hmm->begin[k] = 2. * (float) (hmm->M-k+1) / (float) hmm->M / (float) (hmm->M+1);

  /* Configure exit:   1/(M-k+1)  (xref STL9/77)
   */
  for (k = 1; k <= hmm->M; k++)
    hmm->end[k] = 1. / (float) (hmm->M-k+1);

  if ((status = left_wing_retraction_imposed(hmm)) != eslOK) goto ERROR;
  if ((status = right_wing_retraction_imposed(hmm)) != eslOK) goto ERROR;
  if ((status = local_dpath_accounting(hmm))        != eslOK) goto ERROR;
  if ((status = p7_ReconfigLength(hmm, gm))         != eslOK) goto ERROR;
  if ((status = logoddsify_the_rest(hmm))           != eslOK) goto ERROR;

  hmm->mode   = p7_LOCAL; 	/* hmmfs mode: local, multihit */
  hmm->flags |= PLAN7_HASBITS;  /* we're configured */
  return eslOK;

 ERROR:
  return status;
}


/*****************************************************************
 * 3. Wing retraction routines.
 *****************************************************************/

/* _imposed() forms are for local alignment: the algorithm mode
 *      sets the entry or exit probabilities to predetermined values.
 *      
 * _added() forms are for global/glocal alignment: the alignment
 *      algorithm per se doesn't allow internal entry or exit, so
 *      the wings are folded into the entry/exit probs.
 *      
 * xref STL9/81.
 */



/* left_wing_retraction_imposed()
 * 
 * Wing retraction, when the B->M_k entry distribution is imposed by
 * the algorithm (sw, fs modes). No calculation is needed in this
 * case. The D_1 state is simply removed from the model.
 * 
 * bsc[1..M] scores are set. 
 * 
 * xref STL8/91.
 */
static int
left_wing_retraction_imposed(P7_PROFILE *gm)
{
  int k;

  for (k = 1; k <= gm->M; k++)
    gm->bsc[k] = p7_Prob2Score(gm->begin[k], gm->bg->p1);

  /* Virtual removal of D_1; assure transitions are impossible */
  if (gm->M > 1) {
    gm->tsc[p7_TDM][1] = p7_IMPOSSIBLE;
    gm->tsc[p7_TDD][1] = p7_IMPOSSIBLE;
  }
  return eslOK;
}


/* left_wing_retraction_added()
 * 
 * Wing retraction, where B->M_k entry distribution comes entirely
 * from retracted paths for k>1 (begin[1] = 1.0 from algorithm; ls
 * mode, for example).
 * 
 * Sets bsc[1..M] (the entry/begin scores), using the core
 * model and the algorithmic begin[] distribution.
 * 
 * xref STL9/81.
 */
static int
left_wing_retraction_added(P7_PROFILE *gm)
{
  int    k;
  float  cumul;
  float *bmk = NULL;		/* log B->D...D->M_k path probabilities, BMk */
  float  x;
  int    status;

  ESL_ALLOC(bmk, sizeof(float) * (gm->M+1));

  /* Calculate the log B->M_k path probabilities; xref STL9/81
   */
  bmk[1] = log(gm->hmm->t[0][p7_TMM]);
  cumul  = log(gm->hmm->t[0][p7_TMD]);
  for (k = 2; k <= gm->M; k++)
    {
      bmk[k]  = cumul + log(gm->hmm->t[k-1][p7_TDM]);
      cumul  += log(gm->hmm->t[k-1][p7_TDD]);
    }

  /* Renormalization (has little if any effect)
   * 
   * <cumul> now contains the log P of the B->D_1...D_M->E mute path
   * that we're removing. If (1-BE) is significantly different than
   * 1.0, renormalize the B distribution by dividing by (1-BE).  
   * Because log(1-x) = -x for x << 1, we know that subtracting
   * log(1-BE) from a log prob is only significant if logBE > log epsilon.
   */
  if (cumul > log(FLT_EPSILON)) { 
    x = log(1. - exp(cumul));
    for (k = 1; k <= gm->M; k++)
      bmk[k] -= x;
  }

  /* Conversion to scores. 
   * At this step, we're assuming that hmm->begin[k] = 0.0 for
   * k>1: the algorithm has no internal entries of its own, and
   * internal entry comes exclusively from paths through D states.
   */
  for (k = 1; k <= gm->M; k++)
     gm->bsc[k] = p7_LL2Score(bmk[k], gm->bg->p1);
  
  /* Virtual removal of D_1 state.
   */
  if (gm->M > 1) {
    gm->tsc[p7_TDM][1] = p7_IMPOSSIBLE;
    gm->tsc[p7_TDD][1] = p7_IMPOSSIBLE;
  }
  free(bmk);
  return eslOK;

 ERROR:
  if (bmk != NULL) free(bmk);
  return status;
}

/* right_wing_retraction_imposed()
 * 
 * Wing retraction for exits, for algorithms where M_k->E exit
 * probabilities are imposed by the algorithm (sw, fs modes). 
 * 
 * Sets esc[1..M] (the exit scores); also sets tsc[TM*][1..M-1], which
 * are affected by the presence of the new M_k->E probabilities.
 * 
 * xref STL9/81.
 */
static int
right_wing_retraction_imposed(P7_PROFILE *gm)
{
  int    k;
  float *mke = NULL;
  float  cumul;
  float  x;			/* temporary log prob */
  int    status;

  ESL_ALLOC(mke, sizeof(float) * (hmm->M+1));  

  /* The log prob of the wing-retracted M_k -> D...D -> E paths,
   * for k < M. (undefined for k == M).
   */
  cumul = 0.;
  for (k = gm->M-1; k >= 1; k--)
    {
      mke[k] = cumul + log(gm->hmm->t[k][p7_TMD]);
      cumul += log(gm->hmm->t[k][TDD]);
    }

  /* Set the esc[] and tsc[][TM*] scores.
   * 
   * The MkE path probability is subtracted from t[k][TMD] transition.
   * The match transitions are renormalized to account for the new
   * end[k] probability. The match transitions are also renormalized
   * to account for the newly missing MkE path probability from TMD.
   * (xref STL9/81 for details).
   */
  for (k = 1; k < gm->M; k++)
    {
      gm->esc[k] = p7_Prob2Score(gm->end[k], 1.0); /* end[k] is imposed. */

      x = log(gm->hmm->t[k][p7_TMM]);
      if (gm->end[k] > FLT_EPSILON)  x += log(1. - gm->end[k]);
      if (mke[k] > log(FLT_EPSILON)) x -= log(1. - exp(mke[k]));
      gm->tsc[p7_TMM][k] = p7_LL2Score(x, gm->bg->p1);

      x = log(gm->hmm->t[k][p7_TMI]);
      if (gm->end[k] > FLT_EPSILON)  x += log(1. - gm->end[k]);
      if (mke[k] > log(FLT_EPSILON)) x -= log(1. - exp(mke[k]));
      gm->tsc[p7_TMI][k] = p7_LL2Score(x, gm->bg->p1);

      x = log(gm->hmm->t[k][p7_TMD]);
      if (mke[k] - x > log(FLT_EPSILON)) x += log(1. - exp(mke[k] - x));
      if (gm->end[k] > FLT_EPSILON)      x += log(1. - gm->end[k]);         
      if (mke[k] > log(FLT_EPSILON))     x -= log(1. - exp(mke[k]));
      gm->tsc[TMD][k] = p7_LL2Score(x, 1.0);
    }
  gm->esc[hmm->M] = 0.0;	/* by definition */

  /* Note that node M isn't even explicitly represented in the
   * configured HMM scores -- tsc[][] only contains values for
   * 1..M-1. So there's no need to set any M_M or D_M transition
   * scores to 0 and -INFTY as we virtually remove D_M state;
   * the only other affected score is tsc[TDM][hmm->M-1], which
   * is going to be set when we deal with delete paths.
   * So, we're done.
   */
  free(mke);
  return eslOK;

 ERROR:
  if (mke != NULL) free(mke);
  return status;
}


/* right_wing_retraction_added()
 * 
 * Retract the right wing (remove the D_M state, and all paths through
 * it), for algorithms which have no M_k->E end[k] internal exit
 * probability. The Mk->Dk+1...DM->E path probabilities are therefore
 * subtracted from t[k][TMD] and added to end[k].
 *
 * Sets esc[1..M] (the exit scores); also sets tsc[TM*][1..M-1], which
 * are affected by the presence of the new M_k->E probabilities.
 * 
 * xref STL9/81.
 */
static void
right_wing_retraction_added(struct plan7_s *hmm)
{
  int    k;
  float *mke;
  float  cumul;
  float  x;			/* temporary log prob */

  mke = MallocOrDie(sizeof(float) * (hmm->M+1));  

  /* The log prob of the wing-retracted M_k -> D...D -> E paths,
   * for k < M. (undefined for k == M).
   */
  cumul = 0.;
  for (k = hmm->M-1; k >= 1; k--)
    {
      mke[k] = cumul + log(hmm->t[k][TMD]);
      cumul += log(hmm->t[k][TDD]);
    }

  /* Set the esc[] and tsc[TM*][] scores.
   * 
   * The end probability is assumed to be exclusively the MkE
   * path probability; algorithm has no internal exit prob of its own.
   * 
   * The MkE path probability is moved from the t[k][TMD] transition
   * to the end[k] probability. No renormalization is needed, because
   * prob is conserved: we assume that the algorithm added no
   * internal exit probability end[k] of its own.
   * (xref STL9/91 for details).
   */
  for (k = 1; k < hmm->M; k++)
    {
      hmm->esc[k] = LL2Score(mke[k], 1.0); /* M->E comes only thru terminal deletes */

      hmm->tsc[TMM][k] = Prob2Score(hmm->t[k][TMM], hmm->p1);
      hmm->tsc[TMI][k] = Prob2Score(hmm->t[k][TMI], hmm->p1);
      
      x = log(hmm->t[k][TMD]);
      if ((mke[k] - x) > log(FLT_EPSILON)) x += log(1. - exp(mke[k] - x));
      hmm->tsc[TMD][k] = LL2Score(x, 1.0);
    }
  hmm->esc[hmm->M] = 0.0;	/* by definition */

  /* Note that node M isn't even explicitly represented in the
   * configured HMM scores -- tsc[][] only contains values for
   * 1..M-1. So there's no need to set any M_M or D_M transition
   * scores to 0 and -INFTY as we virtually remove D_M state;
   * the only other affected score is tsc[TDM][hmm->M-1], which
   * is going to be set when we deal with delete paths.
   * So, we're done.
   */
  free(mke);
  return;
}
