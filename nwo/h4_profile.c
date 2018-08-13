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
  hmm->e   = NULL;
  hmm->t   = NULL;
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
  if ((hmm->t = esl_mat_FCreate((M+1), h4_NTRANSITIONS)) == NULL) goto ERROR;
  if ((hmm->e = esl_mat_FCreate((M+1), abc->K))          == NULL) goto ERROR;

  /* Initialize to 0. */
  esl_mat_FSet(hmm->t, M+1, h4_NTRANSITIONS, 0.);
  esl_mat_FSet(hmm->e, M+1, abc->K,          0.);

  /* Boundary conditions. {see h4_profile.md} */
  hmm->e[0][0]      = 1.;
  hmm->t[0][h4_TIM] = 1.;
  hmm->t[0][h4_TDM] = 1.;
  hmm->t[M][h4_TMM] = 1.;
  hmm->t[M][h4_TIM] = 1.;
  hmm->t[M][h4_TDM] = 1.;

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
  n += esl_mat_FSizeof((hmm->M+1), h4_NTRANSITIONS);
  n += esl_mat_FSizeof((hmm->M+1), hmm->abc->K);
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
      for (z = 0; z < h4_NTRANSITIONS; z++)
	fprintf(fp, "%10.2f%c", hmm->t[k][z], z == h4_NTRANSITIONS-1 ? '\n':' ');
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
    if ( esl_vec_FCompare(h1->t[k], h2->t[k], h4_NTRANSITIONS, tol) != eslOK)
      ESL_FAIL(eslFAIL, NULL, "difference in state transition vector %d", k);

  return eslOK;
}
