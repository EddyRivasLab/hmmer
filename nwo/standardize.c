/* Set standard profile's scores, from probability parameters.
 * 
 * Contents:
 *    1. h4_standardize()
 *    2. Setting entry/exit scores
 *     
 * See also:
 *    parameterize.c : determine probability parameters from counts
 *    vectorize.c    : set striped vector scores from standard scores
 */
#include <h4_config.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "h4_profile.h"

static int set_local_entry (H4_PROFILE *hmm);
static int set_glocal_entry(H4_PROFILE *hmm);
static int set_glocal_exit (H4_PROFILE *hmm);

/*****************************************************************
 * 1. h4_standardize()
 *****************************************************************/

/* Function:  h4_standardize()
 * Synopsis:  Set the scores in a profile HMM, from the probabilities.
 * Incept:    SRE, Wed 08 May 2019
 */
int
h4_standardize(H4_PROFILE *hmm)
{
  float  sc[h4_MAXCODE];   // tmp space for calculating residue scores, including degeneracies
  int    k,x;

  ESL_DASSERT1(( hmm->flags & h4_HASPROBS ));

  set_local_entry(hmm);  //  [LM]    'uniform' local entry  -  sets tsc[(0)1..M][h4_LM]
  set_glocal_entry(hmm); //  [GM|GI] left wing retraction   -  sets tsc[0..M-1(M)][h4_GM] and tsc[(0)1..M-1(M)][h4_GI]
  set_glocal_exit(hmm);  //  [DGE]   right wing retraction  

  /* Correct boundary conditions on hmm->t[0], t[M] will result in
   * correct boundary conditions on tsc[0], tsc[M]... especially since
   * esl_log2f(0) = -eslINFINITY.
   */
  for (k = 0; k <= hmm->M; k++) {
    hmm->tsc[k][h4_MM] = esl_log2f(hmm->t[k][h4_TMM]);
    hmm->tsc[k][h4_IM] = esl_log2f(hmm->t[k][h4_TIM]);
    hmm->tsc[k][h4_DM] = esl_log2f(hmm->t[k][h4_TDM]);
    hmm->tsc[k][h4_MD] = esl_log2f(hmm->t[k][h4_TMD]);
    hmm->tsc[k][h4_ID] = esl_log2f(hmm->t[k][h4_TID]);
    hmm->tsc[k][h4_DD] = esl_log2f(hmm->t[k][h4_TDD]);
    hmm->tsc[k][h4_MI] = esl_log2f(hmm->t[k][h4_TMI]);
    hmm->tsc[k][h4_II] = esl_log2f(hmm->t[k][h4_TII]);
    hmm->tsc[k][h4_DI] = esl_log2f(hmm->t[k][h4_TDI]);
  }

  sc[hmm->abc->K]    = -eslINFINITY; // gap
  sc[hmm->abc->Kp-2] = -eslINFINITY; // nonresidue
  sc[hmm->abc->Kp-1] = -eslINFINITY; // missing data
  for (k = 0; k <= hmm->M; k++)                       // k=0 will use boundary conditions from e[0]
    {
      for (x = 0; x < hmm->abc->K; x++)
	sc[x] = esl_log2f(hmm->e[k][x] / hmm->f[x]);   // log-odds scores for canonical alphabet
      esl_abc_FExpectScVec(hmm->abc, sc, hmm->f);     //     ...  scores for degeneracies

      for (x = 0; x < hmm->abc->Kp; x++)              // stored transposed
	hmm->rsc[x][k] = sc[x];
    }

  hmm->flags |= h4_HASBITS;
  return eslOK;
}

/*****************************************************************
 * 2. Setting entry/exit scores
 *****************************************************************/

/* set_local_entry()
 * 
 * Local mode entry prob tLMk is approx. 2/(M(M+1)), with a
 * correction for match state occupancy probability [Eddy08]:
 *    L->Mk = occ[k] /( \sum_i occ[i] * (M-i+1))
 *    
 * We store these params off-by-one, with tLMk stored in
 * hmm->t[k-1][h4_LM], for DP efficiency reasons.
 *
 * Therefore hmm->tsc[0..M-1][h4_LM] have calculated scores for the
 * LMk=1..M entry distribution; hmm->tsc[M][h4_LM] is set to -inf.
 * 
 * We need space for an M+1 occ[k] array of match occupancy.  We save
 * a malloc by using hmm->rsc[0]'s space, which we know is >= M+1, and
 * which we guarantee (by order of calls in config calls) has not been
 * parameterized yet.
 * 
 * If the <hmm>'s <h4_SINGLE> flag is up, <occ[1..M]> is 1.0:
 * i.e. the match state occupancy term is only applied to profiles,
 * not to single-seq queries. 
 */
static int
set_local_entry(H4_PROFILE *hmm)
{
  float *occ = hmm->rsc[0];   // safe malloc-saving hack, so long as we set transition scores before emissions
  float  Z   = 0.;
  int    k;

  if (hmm->flags & h4_SINGLE)
    {
      for (k = 1; k <= hmm->M; k++)
	hmm->tsc[k-1][h4_LM] = esl_log2f( 2. / ((float) hmm->M * (float) (hmm->M+1)));
    }
  else
    {
      h4_profile_Occupancy(hmm, occ, NULL, NULL, NULL);
      for (k = 1; k <= hmm->M; k++)
	Z += occ[k] * (float) (hmm->M + 1 - k );
      for (k = 1; k <= hmm->M; k++)
	hmm->tsc[k-1][h4_LM] = esl_log2f(occ[k] / Z);  // watch out for occ[k] = 0. esl_log2f(0) = -inf
    }

  hmm->tsc[hmm->M][h4_LM] = -eslINFINITY;
  return eslOK;
}


/* set_glocal_entry()
 * 
 * glocal entry is "wing retracted" into G->(D1..Dk-1)->Mk and
 * G->D1..Dk->Ik entry distributions, which we use to enter on {M|I}k
 * states.
 *
 * Wing-retracted entry is used in sparse DP to avoid having to
 * calculate through supercells that aren't in the sparse mask.
 *
 * Wing retraction is also used to remove the B->G->D1..Dm->E->J->B
 * mute cycle, by removing the D1 state. Our DP algorithms neglect
 * (disallow) this mute cycle. As a result, a (usually negligible)
 * amount of probability is unaccounted for. If you need to know what
 * that neglected mass is, see <h4_profile_GetMutePathLogProb()>.
 *
 * Unretracted "base" transition parameters G->{M1,D1} are still in
 * <tsc[0][h4_MM | h4_MD]> if you need them, but DP algorithms must use
 * wing retraction or not: G->{MI}k or G->D1: not both.
 *
 * Wing-retracted entries G->{D1..Dk-1}->Mk are computed and stored
 * off-by-one in tsc[][GM] as follows:
 *      tsc[0][GM]   = tGM1 = log t(G->M1) 
 *      tsc[k-1][GM] = tGMk = log t(G->D1) + \sum_j={2..k-1} log t(Dj-1 -> Dj) + log t(Dk-1->Mk)
 *      tsc[M][GM]   = -inf
 *
 * Wing-retracted G->(D1..Dk)->Ik are computed as:
 *      tGI0 = -inf   (no I0 state)
 *      tGIk = log t(G->D1) + \sum_{j=2..k} log t(Dj-1 -> Dj) + log t(Dk->Ik)
 *      tGIm = -inf   (no Im state)
 * and stored normally at tsc[k]. 
 */
static int
set_glocal_entry(H4_PROFILE *hmm)
{
  float tmp;
  int   k;

  hmm->tsc[0][h4_GM] = esl_log2f(hmm->t[0][h4_TMM]); // i.e. G->M1 transition
  hmm->tsc[0][h4_GI] = -eslINFINITY;                 // (no I0 state)
  tmp                = esl_log2f(hmm->t[0][h4_TMD]); // i.e. G->D1
  for (k = 1; k < hmm->M; k++) 
    {
      hmm->tsc[k][h4_GM] = tmp + esl_log2f(hmm->t[k][h4_TDM]);  // because of the off-by-one storage, this is G->D1..Dk->Mk+1, stored at k
      hmm->tsc[k][h4_GI] = tmp + esl_log2f(hmm->t[k][h4_TDI]);  // ... whereas this is G->D1..Dk->Ik, stored at k
      tmp += esl_log2f(hmm->t[k][h4_TDD]);
    }
  hmm->tsc[hmm->M][h4_GM] = -eslINFINITY;  // no Mm+1 state (remember, off by one storage)
  hmm->tsc[hmm->M][h4_GI] = -eslINFINITY;  // no Im state      

  return eslOK;
}


/* set_glocal_exit()
 * 
 * Right wing retraction, DGE
 *   t[k][DGE] = t(Dk+1->...Dm->E) 
 *             = [\prod_j=k+1..m-1 t(Dj->Dj+1)] * Dm->E
 *             = \prod_j=k+1..m-1 t(Dj->Dj+1)              | because Dm->E=1.0
 *  valid for k=0..M, but k=0 is unused (mute path), so here k=1..M is stored
 *  boundaries: tsc[M][DGE]   = 0     (i.e. probability 1)
 *              tsc[M-1][DGE] = 0  
 *              tsc[0][DGE]   = -inf  (i.e. probability 0)
 * note off by one storage.  
 * to get the glocal exit path from a Mk: tsc[k][MD] + tsc[k][DGE]
 *                                ... Ik: tsc[k][ID] + tsc[k][DGE]
 *                                ... Dk: tsc[k][DD] + tsc[k][DGE]
 */
static int
set_glocal_exit(H4_PROFILE *hmm)
{
  float ppath = 0.0;		/* accumulated Dk+1..Dm part of the path */
  int   k;

  hmm->tsc[hmm->M][h4_DGE] = 0.;
  for (k = hmm->M-1; k >= 1; k--)
    {
      hmm->tsc[k][h4_DGE] = ppath;
      ppath += esl_log2f(hmm->t[k][h4_TDD]);
    }
  hmm->tsc[0][h4_DGE] = -eslINFINITY;
  return eslOK;
}

