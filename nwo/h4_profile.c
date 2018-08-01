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
h4_profile_Create(ESL_ALPHABET *abc, int M)
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



int
h4_profile_SetConventions(H4_PROFILE *hmm)
{
  esl_vec_FSet(hmm->e[0], hmm->abc->K, 0.0);
  hmm->e[0][0] = 1.0;

  hmm->t[0][h4_TMI] = 0.;        // at [0] TMM,TMD are G->{MD}1.

  hmm->t[hmm->M][h4_TMM] = 1.0; 
  hmm->t[hmm->M][h4_TMI] = 0.0;
  hmm->t[hmm->M][h4_TMD] = 0.0;

  hmm->t[0][h4_TIM] = hmm->t[hmm->M][h4_TIM] = 1.0;  // at [0] and [M], there is no insert state; 
  hmm->t[0][h4_TII] = hmm->t[hmm->M][h4_TII] = 0.0;  //   these settings are documented conventions.
  hmm->t[0][h4_TID] = hmm->t[hmm->M][h4_TID] = 0.0;

  hmm->t[0][h4_TDM] = hmm->t[hmm->M][h4_TDM] = 1.0;  // at [0] there is no delete state; 
  hmm->t[0][h4_TDI] = hmm->t[hmm->M][h4_TDI] = 0.0;  // at [M], delete -> E.
  hmm->t[0][h4_TDD] = hmm->t[hmm->M][h4_TDD] = 0.0;  // these settings are documented conventions.

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
