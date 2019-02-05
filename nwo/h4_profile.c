/* H4_PROFILE: dual-mode local/glocal profile HMM
 * 
 * Contents:
 *    1. H4_PROFILE structure
 *    
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_matrixops.h"
#include "esl_vectorops.h"

#include "h4_profile.h"

#include "general.h"        // h4_AminoFrequencies()


/* Function:  h4_profile_Create()
 * Synopsis:  Allocate a new profile.
 * Incept:    SRE, Sun 08 Jul 2018 [JB1486 PIT-BOS]
 *
 * Purpose:   Allocates a new empty profile of length 
 *            <M> consensus positions, for alphabet <abc>.
 *
 * Args:      abc : digital alphabet
 *            M   : model length in nodes (consensus positions)
 *
 * Returns:   ptr to new <H4_PROFILE>
 *
 * Throws:    <NULL> on allocation failure.
 */
H4_PROFILE *
h4_profile_Create(const ESL_ALPHABET *abc, int M)
{
  H4_PROFILE *hmm = NULL;

  if ((hmm    = h4_profile_CreateShell()) == NULL)  return NULL;
  if ( h4_profile_CreateBody(hmm, abc, M) != eslOK) return NULL;
  return hmm;
}


/* Function:  h4_profile_CreateShell()
 * Synopsis:  First step of two-step allocation of a profile
 * Incept:    SRE, Sun 08 Jul 2018 [JB1486 PIT-BOS]
 *
 * Purpose:   Allocate a new profile except for things that depend on
 *            model size <M> or alphabet size <K>. 
 *            
 *            When we read models from files, we need be storing data
 *            before we've read <M> and <K>, so we provide for
 *            allocating in two steps.
 *            
 * Returns:   ptr to new <H4_PROFILE>
 *
 * Throws:    <NULL> on allocation failure
 */
H4_PROFILE *
h4_profile_CreateShell(void)
{
  H4_PROFILE *hmm = NULL;
  int         status;

  ESL_ALLOC(hmm, sizeof(H4_PROFILE));
  hmm->M   = 0;
  hmm->t   = NULL;
  hmm->e   = NULL;
  hmm->f   = NULL;
  hmm->tsc = NULL;
  hmm->rsc = NULL;
  hmm->abc = NULL;
  return hmm;

 ERROR:
  h4_profile_Destroy(hmm);
  return NULL;
}

/* Function:  h4_profile_CreateBody()
 * Synopsis:  Second step of two-step allocation of profile
 * Incept:    SRE, Sun 08 Jul 2018 [JB1486 PIT-BOS]
 *
 * Purpose:   Given a profile <hmm> shell allocated by 
 *            <h4_profile_CreateShell()>, now allocate stuff 
 *            that depends on knowing alphabet <abc>
 *            and profile length <M>. Initialize
 *            data-dependent probability parameters to 0,
 *            and set data-independent boundary conditions
 *            as described in {h4_profile.md}.
 *            
 *            <hmm->M> is set, and a copy of <abc> is kept.  Caller
 *            remains responsible for <abc>, but must not free it
 *            until this profile is freed.
 *
 * Returns:   <eslOK> on success, and allocations in
 *            <hmm> are done.
 *
 * Throws:    <eslEMEM> on allocation failure
 */
int
h4_profile_CreateBody(H4_PROFILE *hmm, const ESL_ALPHABET *abc, int M)
{
  int status;

  if ((hmm->t   = esl_mat_FCreate( M+1,     h4_NT))   == NULL) goto ERROR;
  if ((hmm->e   = esl_mat_FCreate( M+1,     abc->K))  == NULL) goto ERROR;
  if ((hmm->tsc = esl_mat_FCreate( M+1,     h4_NTSC)) == NULL) goto ERROR;
  if ((hmm->rsc = esl_mat_FCreate( abc->Kp, M+1))     == NULL) goto ERROR; 

  ESL_ALLOC(hmm->f, sizeof(float) * abc->K);

  /* Initialize probability parameters to 0. */
  esl_mat_FSet(hmm->t, M+1, h4_NT,   0.);
  esl_mat_FSet(hmm->e, M+1, abc->K,  0.);

  /* Boundary conditions on probability params. {see h4_profile.md} */
  hmm->e[0][0]      = 1.;
  hmm->t[0][h4_TIM] = 1.;
  hmm->t[0][h4_TDM] = 1.;
  hmm->t[M][h4_TMM] = 1.;
  hmm->t[M][h4_TIM] = 1.;
  hmm->t[M][h4_TDM] = 1.;

  if (abc->type == eslAMINO) h4_AminoFrequencies(hmm->f);
  else                       esl_vec_FSet(hmm->f, abc->K, 1. / (float) abc->K);

  hmm->M   = M;
  hmm->abc = abc;
  return eslOK;

 ERROR:
  h4_profile_Destroy(hmm);
  return eslEMEM;
}

/* Function:  h4_profile_Sizeof()
 * Synopsis:  Returns allocated size of a profile, in bytes.
 */
size_t
h4_profile_Sizeof(H4_PROFILE *hmm)
{
  size_t n = 0;
  n += sizeof(H4_PROFILE);
  n += esl_mat_FSizeof(hmm->M+1,       h4_NT);         // t
  n += esl_mat_FSizeof(hmm->M+1,       hmm->abc->K);   // e
  n += sizeof(float) * hmm->abc->K;                    // f
  n += esl_mat_FSizeof(hmm->M+1,       h4_NTSC);       // tsc
  n += esl_mat_FSizeof(hmm->abc->Kp+1, hmm->M+1);      // rsc
  return n;
}


/* Function:  h4_profile_Destroy()
 * Synopsis:  Frees a profile HMM.
 * Incept:    SRE, Sun 08 Jul 2018 [JB1486 PIT-BOS]
 */
void
h4_profile_Destroy(H4_PROFILE *hmm)
{
  if (hmm)
    {
      esl_mat_FDestroy(hmm->t);
      esl_mat_FDestroy(hmm->e);
      free(hmm->f);
      esl_mat_FDestroy(hmm->tsc);
      esl_mat_FDestroy(hmm->rsc);
      free(hmm);
    }
}



/* Function:  h4_profile_SetConventions()
 * Synopsis:  Set the fixed edge conditions in a H4 profile HMM.
 * Incept:    SRE, Mon 06 Aug 2018
 */
int
h4_profile_SetConventions(H4_PROFILE *hmm)
{
  esl_vec_FSet(hmm->e[0], hmm->abc->K, 0.0);         // e[0] is unused; we make it a valid probability vector anyway
  hmm->e[0][0] = 1.0;

  hmm->t[0][h4_TMI] = 0.;                            // at [0] TMM,TMD are G->{MD}1; no I0 state, so G->MI0 = 0

  hmm->t[0][h4_TIM] = hmm->t[hmm->M][h4_TIM] = 1.0;  // at [0] and [M], there is no insert state;
  hmm->t[0][h4_TII] = hmm->t[hmm->M][h4_TII] = 0.0;  //   we make the transitions valid prob vectors anyway
  hmm->t[0][h4_TID] = hmm->t[hmm->M][h4_TID] = 0.0;

  hmm->t[0][h4_TDM] = hmm->t[hmm->M][h4_TDM] = 1.0;  // at [0] there is no delete state; at [M], delete -> E.
  hmm->t[0][h4_TDI] = hmm->t[hmm->M][h4_TDI] = 0.0;  
  hmm->t[0][h4_TDD] = hmm->t[hmm->M][h4_TDD] = 0.0;  

  hmm->t[hmm->M][h4_TMM] = 1.0;                      // at [M], match state must go M->E.
  hmm->t[hmm->M][h4_TMI] = 0.0;
  hmm->t[hmm->M][h4_TMD] = 0.0;

  return eslOK;
}


int
h4_profile_Renormalize(H4_PROFILE *hmm)
{
  float sum;
  int   k;

  for (k = 1; k <= hmm->M; k++) esl_vec_FNorm(hmm->e[k], hmm->abc->K);
  for (k = 1; k <  hmm->M; k++) esl_vec_FNorm(hmm->t[k],   3);
  for (k = 1; k <  hmm->M; k++) esl_vec_FNorm(hmm->t[k]+3, 3);
  for (k = 1; k <  hmm->M; k++) esl_vec_FNorm(hmm->t[k]+6, 3);

  /* You have to be unusually careful with the t[0] match transitions,
   * which are the G->{MD}1 glocal entry transitions. t[0][TMI] must
   * remain 0.0, which it won't do if there are zero counts and you
   * call FNorm() on the vector.
   */
  sum = hmm->t[0][h4_TMM] + hmm->t[0][h4_TMD];
  if (sum > 0) esl_vec_FScale(hmm->t[0], 3, 1./sum);
  else         hmm->t[0][h4_TMM] = hmm->t[0][h4_TMD] = 0.5;

  h4_profile_SetConventions(hmm);

  return eslOK;
}



/* Function:  h4_profile_CalculateOccupancy()
 * Synopsis:  Calculate match occupancy and insert expected use count vectors.
 *
 * Purpose:   Calculate a vector <mocc[1..M]> containing probability
 *            that each match state is used in a sampled path through
 *            the model. Caller provides allocated space (<M+1> floats)
 *            for <mocc>.
 *            
 *            Caller may optionally provide an array <iocc[0..M]> as
 *            well, which (if provided) will be set to contain the
 *            expected number of times that a sampled path would contain
 *            each insert state.
 *            
 * Returns:   <eslOK> on success.
 * 
 * TODO:      unit testing. I'm worried about roundoff error accumulation
 *            especially on long models.
 */
int
h4_profile_CalculateOccupancy(const H4_PROFILE *hmm, float *mocc, float *iocc)
{
  int   k;
  float mshare;

  mocc[0] = 0.;			        // no M_0 state 
  mocc[1] = hmm->t[0][h4_TMM];          // initialize w/ G->M1 
  for (k = 2; k <= hmm->M; k++) {
    mshare   = hmm->t[k-1][h4_TIM] / (hmm->t[k-1][h4_TIM] + hmm->t[k-1][h4_TID]);
    mocc[k]  =     mocc[k-1]  * hmm->t[k-1][h4_TMM];           // M->M
    mocc[k] += (1.-mocc[k-1]) * hmm->t[k-1][h4_TDM];           // D->M
    mocc[k] +=     mocc[k-1]  * hmm->t[k-1][h4_TMI]  * mshare; // M->I->M
    mocc[k] += (1.-mocc[k-1]) * hmm->t[k-1][h4_TDI]  * mshare; // D->I->M
    mocc[k] = ESL_MIN(1.0, mocc[k]);                           // avoid floating pt roundoff error making mocc[k] 1+epsilon; mocc[k] is a probability
  }

  if (iocc) {
    iocc[0] = 0.;                     // no I0 state
    for (k = 1; k < hmm->M; k++) {    // ... also no Im state
      iocc[k]  =     mocc[k]  * hmm->t[k][h4_TMI];
      iocc[k] += (1.-mocc[k]) * hmm->t[k][h4_TDI];
      iocc[k] /= (1.-hmm->t[k][h4_TII]);
    }
  }

  return eslOK;
}


/*****************************************************************
 * x. Profile configuration: probabilities to scores
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
      h4_profile_CalculateOccupancy(hmm, occ, NULL);
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
 * glocal entry is "wing retracted" into a G->Mk entry distribution.
 *
 * G->Mk entry is used in sparse DP to avoid having to calculate
 * through supercells that aren't in the sparse mask.
 *
 * Wing retraction also removes the B->G->D1..Dm->E->J->B mute cycle,
 * by removing the D1 state. Our DP algorithms neglect (disallow) this
 * mute cycle. As a result, a (usually negligible) amount of
 * probability is unaccounted for. If you need to know what that
 * neglected mass is, see <p7_profile_GetMutePathLogProb()>.
 *
 * Unretracted "base" transition parameters G->{M1,D1} are still in
 * <tsc[0][h4_MM | h4_MD]> if you need them.
 *
 * We compute the wing-retracted entries G->{D1..Dk-1}->Mk as follows:
 *      tGM1 = log t(G->M1) 
 *      tGMk = log t(G->D1) + \sum_j={1..k-2} log t(Dj->Dj+1) + log t(Dk-1->Mk)
 * and for technical/efficiency reasons, these params are stored
 * off-by-one in the profile: tGMk is stored at <tsc[k-1][h4_GM]>.
 */
static int
set_glocal_entry(H4_PROFILE *hmm)
{
  float ppath;   
  int   k;

  hmm->tsc[0][h4_GM] = esl_log2f(hmm->t[0][h4_TMM]); // G->M1
  ppath              = esl_log2f(hmm->t[0][h4_TMD]); // G->D1
  for (k = 1; k < hmm->M; k++) 
    {
      hmm->tsc[k][h4_GM] = ppath + esl_log2f(hmm->t[k][h4_TDM]);
      ppath             += esl_log2f(hmm->t[k][h4_TDD]);
    }
  hmm->tsc[hmm->M][h4_GM] = -eslINFINITY;
  return eslOK;
}


/* set_glocal_exit()
 * 
 * Right wing retraction, DGE
 *   tsc[k][DGE] = t(Dk+1->...Dm->E) 
 *               = [\prod_j=k+1..m-1 t(Dj->Dj+1)] * Dm->E
 *               = \prod_j=k+1..m-1 t(Dj->Dj+1)              | because Dm->E=1.0
 *  valid for k=0..M, but k=0 is unused (mute path), so here k=1..M is stored
 *  boundaries: tsc[M][DGE]   = 0     (i.e. probability 1)
 *              tsc[M-1][DGE] = 0  
 *              tsc[0][DGE]   = -inf  (i.e. probability 0)
 * note off by one storage.  
 * to get the glocal exit path from a Mk: tsc[k][MD] + tsc[k][DGE]
 * to get the glocal exit path from a Dk: tsc[k][DD] + tsc[k][DGE]
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
  hmm->t[0][h4_DGE] = -eslINFINITY;
  return eslOK;
}



int
h4_profile_Config(H4_PROFILE *hmm)
{
  float  sc[h4_MAXCODE];   // tmp space for calculating residue scores, including degeneracies
  int    k,x;

  set_local_entry(hmm);  //  [LM] 'uniform' local alignment entry
  set_glocal_entry(hmm); //  [GM] left wing retraction
  set_glocal_exit(hmm);  // [DGE] right wing retraction

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

  return eslOK;
}





/*****************************************************************
 * x. Debugging and development tools
 *****************************************************************/


/* Function:  h4_profile_Dump()
 * Synopsis:  Dump contents of an H4_PROFILE for inspection
 * Incept:    SRE, Mon 16 Jul 2018 [Benasque]
 */
int
h4_profile_Dump(FILE *fp, H4_PROFILE *hmm)
{
  int k,a,z;

  fprintf(fp, "Emissions:\n");
  fprintf(fp, "     ");
  for (a = 0; a < hmm->abc->K; a++)
    fprintf(fp, "         %c%c", hmm->abc->sym[a], a == hmm->abc->K-1 ? '\n':' ');
  for (k = 1; k <= hmm->M; k++)
    {
      fprintf(fp, "%4d ", k);
      for (a = 0; a < hmm->abc->K; a++)
	fprintf(fp, "%10.2f%c", hmm->e[k][a], a == hmm->abc->K-1 ? '\n':' ');
    }

  fprintf(fp, "Transitions:\n");
  fprintf(fp, "     ");
  fprintf(fp, "%10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
	  "TMM", "TMI", "TMD", "TIM", "TII", "TID", "TDM", "TDI", "TDD");
  for (k = 0; k < hmm->M; k++)  // include t[0] because GM, GD (in MM, MD) are data-dependent
    {                           // exclude t[M] which has no data dependent transitions.
      fprintf(fp, "%4d ", k);
      for (z = 0; z < h4_NT; z++)
	fprintf(fp, "%10.2f%c", hmm->t[k][z], z == h4_NT-1 ? '\n':' ');
    }
  return eslOK;
}



/* Function:  h4_profile_Validate()
 * Synopsis:  Validate an <H4_PROFILE> data structure
 * Incept:    SRE, Mon 06 Aug 2018 [Nick Cave and Warren Ellis, Mountain Lion Mean]
 */
int
h4_profile_Validate(const H4_PROFILE *hmm, char *errbuf)
{
  int   k;
  float tol = 1e-4;

  if (hmm->M < 1) ESL_FAIL(eslFAIL, errbuf, "invalid model size M");
  if (! hmm->abc) ESL_FAIL(eslFAIL, errbuf, "no model alphabet");

  /* emissions and transitions */
  for (k = 0; k <= hmm->M; k++)
    if ( esl_vec_FValidate(hmm->e[k],   hmm->abc->K, tol, NULL) != eslOK ||
	 esl_vec_FValidate(hmm->t[k],   3,           tol, NULL) != eslOK ||
         esl_vec_FValidate(hmm->t[k]+3, 3,           tol, NULL) != eslOK ||
	 esl_vec_FValidate(hmm->t[k]+6, 3,           tol, NULL) != eslOK)
      ESL_FAIL(eslFAIL, errbuf, "something awry at state %d", k);

  /* edge conventions */
  if (hmm->e[0][0]           != 1. ||
      hmm->t[0][h4_TMI]      != 0. ||
      hmm->t[0][h4_TIM]      != 1. ||
      hmm->t[0][h4_TDM]      != 1. ||
      hmm->t[hmm->M][h4_TMM] != 1. ||
      hmm->t[hmm->M][h4_TIM] != 1. ||
      hmm->t[hmm->M][h4_TDM] != 1.)
    ESL_FAIL(eslFAIL, errbuf, "something awry in edge conventions");

  return eslOK;
}

/* Function:  h4_profile_Compare()
 * Synopsis:  Compare two <H4_PROFILE>'s for equality
 * Incept:    SRE, Mon 06 Aug 2018 [Colter Wall, Sleeping on the Blacktop]
 */
int
h4_profile_Compare(const H4_PROFILE *h1, const H4_PROFILE *h2)
{
  int   k;
  float tol = 1e-4;

  if (h1->abc->type != h2->abc->type) ESL_FAIL(eslFAIL, NULL, "different alphabets");
  if (h1->M         != h2->M)         ESL_FAIL(eslFAIL, NULL, "different M");

  for (k = 0; k <= h1->M; k++)
    if ( esl_vec_FCompare(h1->e[k], h2->e[k], h1->abc->K, tol) != eslOK)
      ESL_FAIL(eslFAIL, NULL, "difference in match emission vector %d", k);

  for (k = 0; k <= h1->M; k++) 
    if ( esl_vec_FCompare(h1->t[k], h2->t[k], h4_NT, tol) != eslOK)
      ESL_FAIL(eslFAIL, NULL, "difference in state transition vector %d", k);

  return eslOK;
}
