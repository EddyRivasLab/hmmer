/* modelconfig.c
 * 
 * Model configuration: the process of converting a core model to a
 * fully configured Plan7 HMM with an "algorithm mode" built into it.
 * 
 * Revised and refactored May 2005: xref STL9/77-81.
 * 
 * SVN $Id$
 * SRE, Mon May  2 10:55:16 2005 [St. Louis]
 */

/******************************************************************
 * General notes:
 * 
 * When you enter this module, you've got an HMM in "core" probability
 * form: t[], mat[], ins[], and tbd1 are all valid, normalized
 * probabilities. The routines here are used to create the "configured"
 * score form of the model: tsc[], msc[], isc[], bsc[], esc[], and
 * xsc[] fields become valid integer log-odds scores. Also in the process,
 * xt[], begin[], and end[] are set to their algorithm-dependent
 * probabilities.
 * 
 * The configuration process breaks down into three conceptual steps:
 * 
 * 1. Algorithm configuration.
 *    An "algorithm mode" is chosen; this sets the probabilities in
 *    xt[], begin[], and end[], which specify local/global and
 *    multihit/single hit behavior in an HMM alignment to a target
 *    sequence.
 *    
 * 2. Wing retraction.
 *    In a configured model, the D_1 and D_M states are removed.
 *    The probability of the paths B->D1...->Mk ("BMk") that enter
 *    D1 and use all D's before reaching M_k is treated instead
 *    as an additional dollop of B->Mk entry probability, and the
 *    probability of paths Mk->Dk+1...D_M->E ("MkE") is treated
 *    instead as an additional dollop of Mk->E exit probability.
 *    The MkE path probability is subtracted from the Mk->Dk+1
 *    transition.
 *    
 *    If the algorithm mode is imposing its own local entry/exit
 *    probabilities, then these extra dollops are ignored, and the
 *    model is renormalized appropriately. That is, the algorithm
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
 *    DP traceback algorithms check whether the model is imposing
 *    internal entry/exit, and construct tracebacks accordingly;
 *    either leaving B->Mk/Mk->E transitions as is, or interpreting
 *    them as B->DDDD->Mk/Mk->DDDD->E paths.
 *    
 *    Wing retraction has two purposes. First, it removes a mute cycle
 *    from the model, B->D1 ...D_M->E, which cannot be correctly and
 *    efficiently dealt with by DP recursions. (A DP algorithm could
 *    simply *ignore* that path though, and throw away the negligible
 *    amount of probability in it. Wing retraction is an elaborate,
 *    overly stylistic alternative that keeps everything fully
 *    probabilistic. Second, it reconciles the algorithm-dependent
 *    entry/exit probabilities with the core model. For algorithms
 *    that impose particular internal entry/exit, you don't want there
 *    to be any additional probability coming from "internal"
 *    B->DDD->M and M->DDD->E paths, so wing retraction takes it away.
 *    
 * 3. Local alignment D-path leveling.
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
 *    D->D transitions throw a wrench into this trick, though. A local
 *    alignment that goes M_i->D...D->M_j, for example, only gets
 *    hit with one not-end penalty (for the M_i->D). This means that
 *    paths including deletions will be artifactually favored. 
 *    
 *    The solution is to subtract log(1-e_k) from the deletion
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
 * 4. Conversion of probabilities to integer log-odds scores.
 *    This step incorporates the contribution of the null model,
 *    and converts floating-point probs to the scaled integer log-odds
 *    score values that are used by the DP alignment routines.
 *    
 * Step 1 is done by the Plan7*Config() functions, which come in sw,
 * ls, fs, s, and "naked" flavors.
 *
 * Steps 2,3 are both done in P7Logoddsify(), which calls appropriate
 * *wing_retraction*() functions to accomplish step 2/3 for transitions
 * that are affected by wing retraction, then deals with everything
 * else itself.
 * 
 * The model never exists in a configured probability form. The core
 * model is in probabilities; the configured model is in scores. It's
 * possible to interpret the scores as probabilities by backcalculation
 * (Score2Prob()). 
 *****************************************************************
 */

#include "config.h"		/* must be included first */
#include "squidconf.h"

#include <stdio.h>
#include <string.h>

#include "squid.h"
#include "funcs.h"
#include "structs.h"
#include "plan7.h"		/* the model structure */



/*****************************************************************
 * Section 1. The API: plan7's algorithm mode configuration functions.
 *
 * The following few functions are the Plan7 equivalent of choosing
 * different alignment styles (fully local, fully global,
 * global/local, multihit, etc.)
 * 
 * When you come *into* a configuration routine, the following
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
 *    2. xt[XTE][MOVE,LOOP] dist controls multihits. 
 *       Setting xt[XTE][LOOP] to 0.0 forces one hit per model.
 *    3. begin[1..M] controls entry probabilities. An algorithm 
 *       mode either imposes internal begin probabilities, or leaves begin[1] as
 *       1.0 and begin[k] = 0.0 for k>1.
 *    4. end[1..M] controls exit probabilities. An algorithm mode either
 *       imposes internal exit probabilities, or leaves end[M] = 1.0
 *       and end[k] = 0.0 for k<M.
 *    
 * When you leave the configuration routine, these probs (xt, begin, end)
 * have been set. However, they have not been reconciled with the core model;
 * that happens later, in P7Logoddsify(), which reconciles wing retraction
 * with the begin/end probabilities.
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
 * but in the future (when HMM editing may be allowed) we'll have
 * to remember this. Non-equality of the summed domain scores and
 * the total sequence score is a really easy "red flag" for people to
 * notice and report as a bug, even if it may make probabilistic
 * sense not to meet either constraint for certain modeling problems.
 *****************************************************************
 */

/* Function: Plan7NakedConfig()
 * 
 * Purpose:  Set the alignment-independent, algorithm-dependent parameters
 *           of a Plan7 model so that no special states (N,C,J) emit anything:
 *           one simple, full global pass through the model.
 * 
 * Args:     hmm - the plan7 model
 *                 
 * Return:   (void)
 *           The HMM is modified; algorithm dependent parameters are set.
 *           Previous scores are invalidated if they existed.
 */
void
Plan7NakedConfig(struct plan7_s *hmm)                           
{
  hmm->xt[XTN][MOVE] = 1.;	      /* disallow N-terminal tail */
  hmm->xt[XTN][LOOP] = 0.;
  hmm->xt[XTE][MOVE] = 1.;	      /* only 1 domain/sequence ("global" alignment) */
  hmm->xt[XTE][LOOP] = 0.;
  hmm->xt[XTC][MOVE] = 1.;	      /* disallow C-terminal tail */
  hmm->xt[XTC][LOOP] = 0.;
  hmm->xt[XTJ][MOVE] = 0.;	      /* J state unused */
  hmm->xt[XTJ][LOOP] = 1.;

  FSet(hmm->begin+2, hmm->M-1, 0.);   /* disallow internal entries. */
  hmm->begin[1]    = 1.0;
  FSet(hmm->end+1,   hmm->M-1, 0.);   /* disallow internal exits. */
  hmm->end[hmm->M] = 1.0;

  hmm->flags       &= ~PLAN7_BIMPOSED; /* no internal entry -> entry comes from wing retraction  */
  hmm->flags       &= ~PLAN7_EIMPOSED; /* no internal exit  -> exit comes from wing retraction */
  hmm->flags       &= ~PLAN7_HASBITS;  /* reconfig invalidates log-odds scores */
}
   
/* Function: Plan7GlobalConfig()
 * 
 * Purpose:  Set the alignment-independent, algorithm-dependent parameters
 *           of a Plan7 model to global (Needleman/Wunsch) configuration:
 *           "s" mode (hmms).
 * 
 *           Like a non-looping hmmls, since we actually allow flanking
 *           N and C terminal sequence. 
 *           
 * Args:     hmm - the plan7 model
 *                 
 * Return:   (void)
 *           The HMM is modified; algorithm dependent parameters are set.
 *           Previous scores are invalidated if they existed.
 */
void
Plan7GlobalConfig(struct plan7_s *hmm)                           
{
  hmm->xt[XTN][MOVE] = 1. - hmm->p1;  /* allow N-terminal tail */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 1.;	      /* only 1 domain/sequence ("global" alignment) */
  hmm->xt[XTE][LOOP] = 0.;
  hmm->xt[XTC][MOVE] = 1. - hmm->p1;  /* allow C-terminal tail */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 0.;	      /* J state unused */
  hmm->xt[XTJ][LOOP] = 1.;
  FSet(hmm->begin+2, hmm->M-1, 0.);   /* disallow internal entries. */
  hmm->begin[1]    = 1.0;
  FSet(hmm->end+1,   hmm->M-1, 0.);   /* disallow internal exits. */
  hmm->end[hmm->M] = 1.0;

  hmm->flags       &= ~PLAN7_BIMPOSED; /* no internal entry -> entry comes from wing retraction  */
  hmm->flags       &= ~PLAN7_EIMPOSED; /* no internal exit  -> exit comes from wing retraction */
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}
   
/* Function: Plan7LSConfig()
 * 
 * Purpose:  Set the alignment independent parameters of a Plan7 model
 *           to hmmls (global in HMM, local in sequence) configuration.
 *           
 * Args:     hmm  - the plan7 model
 *                 
 * Return:   (void);
 *           the HMM probabilities are modified.
 */
void
Plan7LSConfig(struct plan7_s *hmm)
{
  hmm->xt[XTN][MOVE] = 1.-hmm->p1;    /* allow N-terminal tail */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 0.5;	      /* sets an expectation of 2 domains/seq, geometric dist */
  hmm->xt[XTE][LOOP] = 0.5;
  hmm->xt[XTC][MOVE] = 1.-hmm->p1;    /* allow C-terminal tail */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 1.-hmm->p1;    /* allow J junction state */
  hmm->xt[XTJ][LOOP] = hmm->p1;
  FSet(hmm->begin+2, hmm->M-1, 0.);   /* disallow internal entries */
  hmm->begin[1]    = 1.0;
  FSet(hmm->end+1,   hmm->M-1, 0.);   /* disallow internal exits */
  hmm->end[hmm->M] = 1.0;

  hmm->flags       &= ~PLAN7_BIMPOSED; /* no internal entry -> entry comes from wing retraction  */
  hmm->flags       &= ~PLAN7_EIMPOSED; /* no internal exit  -> exit comes from wing retraction */
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}  
                             

/* Function: Plan7SWConfig()
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a Plan7 model to hmmsw (Smith/Waterman) configuration.
 *           
 * Notes:    The desideratum for begin/end probs is that all fragments ij
 *           (starting at match i, ending at match j) are
 *           equiprobable -- there is no information in the choice of
 *           entry/exit. There are M(M+1)/2 possible choices of ij, so
 *           each must get a probability of 2/M(M+1). This prob is the
 *           product of a begin, an end, and all the not-end probs in
 *           the path between i,j. 
 *            
 *           Thus: entry/exit is asymmetric because of the left/right
 *           nature of the HMM/profile. Entry probability is distributed
 *           simply by assigning p_x = pentry / (M-1) to M-1 
 *           internal match states. However, the same approach doesn't
 *           lead to a flat distribution over exit points. Exit p's
 *           must be corrected for the probability of a previous exit
 *           from the model. Requiring a flat distribution over exit
 *           points leads to an easily solved piece of algebra, giving:
 *                      p_1 = pexit / (M-1)
 *                      p_x = p_1 / (1 - (x-1) p_1)
 *           
 * Args:     hmm    - the Plan7 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over M-1 match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over M-1 match states.
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
Plan7SWConfig(struct plan7_s *hmm, float pentry, float pexit)
{
  float basep;			/* p1 for exits: the base p */
  int   k;			/* counter over states      */

  /* Configure special states.
   */
  hmm->xt[XTN][MOVE] = 1-hmm->p1;    /* allow N-terminal tail */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 1.;	     /* disallow jump state   */
  hmm->xt[XTE][LOOP] = 0.;
  hmm->xt[XTC][MOVE] = 1-hmm->p1;    /* allow C-terminal tail */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 1.;           /* J is unused */
  hmm->xt[XTJ][LOOP] = 0.;

  /* Configure entry.
   * Note: hmm->tbd1 is ignored.
   */
  for (k = 1; k <= hmm->M; k++)
    hmm->begin[k] = 2. * (float) (hmm->M-k+1) / (float) hmm->M / (float) (hmm->M+1);


  hmm->begin[1] = (1. - pentry) * (1. - hmm->tbd1);
  FSet(hmm->begin+2, hmm->M-1, (pentry * (1.- hmm->tbd1)) / (float)(hmm->M-1));
  
  /* Configure exit.
   */
  hmm->end[hmm->M] = 1.0;
  basep = pexit / (float) (hmm->M-1);
  for (k = 1; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) (k-1));
  Plan7RenormalizeExits(hmm);
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}

/* Function: Plan7FSConfig()
 * Date:     SRE, Fri Jan  2 15:34:40 1998 [StL]
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a Plan7 model to hmmfs (multihit Smith/Waterman) configuration.
 *           
 *           See comments on Plan7SWConfig() for explanation of
 *           how pentry and pexit are used.
 *           
 * Args:     hmm    - the Plan7 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over M-1 match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over M-1 match states.
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
Plan7FSConfig(struct plan7_s *hmm, float pentry, float pexit)
{
  float basep;			/* p1 for exits: the base p */
  int   k;			/* counter over states      */

  /* Configure special states.
   */
  hmm->xt[XTN][MOVE] = 1-hmm->p1;    /* allow N-terminal tail     */
  hmm->xt[XTN][LOOP] = hmm->p1;
  hmm->xt[XTE][MOVE] = 0.5;	     /* allow loops / multihits   */
  hmm->xt[XTE][LOOP] = 0.5;
  hmm->xt[XTC][MOVE] = 1-hmm->p1;    /* allow C-terminal tail     */
  hmm->xt[XTC][LOOP] = hmm->p1;
  hmm->xt[XTJ][MOVE] = 1.-hmm->p1;   /* allow J junction between domains */
  hmm->xt[XTJ][LOOP] = hmm->p1;

  /* Configure entry.
   */
  hmm->begin[1] = (1. - pentry) * (1. - hmm->tbd1);
  FSet(hmm->begin+2, hmm->M-1, (pentry * (1.-hmm->tbd1)) / (float)(hmm->M-1));
  
  /* Configure exit.
   */
  hmm->end[hmm->M] = 1.0;
  basep = pexit / (float) (hmm->M-1);
  for (k = 1; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) (k-1));
  Plan7RenormalizeExits(hmm);
  hmm->flags       &= ~PLAN7_HASBITS; /* reconfig invalidates log-odds scores */
}







/************************************************************
 * @LICENSE@
 ************************************************************/
