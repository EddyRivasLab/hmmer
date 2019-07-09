/* H4_PROFILE: dual-mode local/glocal profile HMM
 * 
 * Contents:
 *    1. H4_PROFILE structure
 *    2. Model estimation: counts to probabilities
 *    4. Debugging and development tools
 *    5. Unit tests
 *    6. Test driver
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_alloc.h"
#include "esl_alphabet.h"
#include "esl_matrixops.h"
#include "esl_vectorops.h"

#include "h4_profile.h"

#include "general.h"        // h4_AminoFrequencies()
#include "simdvec.h"

/*****************************************************************
 * 1. H4_PROFILE structure
 *****************************************************************/

/* Function:  h4_profile_Create()
 * Synopsis:  Allocate a new profile.
 * Incept:    SRE, Sun 08 Jul 2018 [JB1486 PIT-BOS]
 *
 * Purpose:   Allocates a new empty profile of length 
 *            <M> consensus positions, for alphabet <abc>.
 *            
 *            The <M> and <abc> fields in the returned profile are
 *            set. The <abc> field is only a copy of the <abc>
 *            ptr. Because of this use of a ptr copy, caller needs to
 *            keep <abc> (unchanged) while this profile is in play.
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
  hmm->M         = 0;
  hmm->V         = 0;
  hmm->Qb        = 0;
  hmm->Qw        = 0;
  hmm->Qf        = 0;
  hmm->f         = NULL;
  hmm->flags     = 0;
  hmm->abc       = NULL;

  hmm->t         = NULL;
  hmm->e         = NULL;

  hmm->tsc       = NULL;
  hmm->rsc       = NULL;

  hmm->rbv       = NULL;
  hmm->tauBM     = 0;

  hmm->rwv       = NULL;
  hmm->twv       = NULL;
  hmm->ddbound_w = 0;

  hmm->rfv       = NULL;
  hmm->tfv       = NULL;

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
  int x, q;
  int status;

  hmm->V  = h4_simdvec_width();
  hmm->Qb = H4_Q(M, hmm->V);
  hmm->Qw = H4_Q(M, hmm->V/2);
  hmm->Qf = H4_Q(M, hmm->V/4);

  ESL_ALLOC(hmm->f, sizeof(float) * abc->K);

  if ((hmm->t   = esl_mat_FCreate( M+1,     h4_NT))   == NULL) goto ERROR;
  if ((hmm->e   = esl_mat_FCreate( M+1,     abc->K))  == NULL) goto ERROR;
  if ((hmm->tsc = esl_mat_FCreate( M+1,     h4_NTSC)) == NULL) goto ERROR;
  if ((hmm->rsc = esl_mat_FCreate( abc->Kp, M+1))     == NULL) goto ERROR; 

  ESL_ALLOC(hmm->rbv, sizeof(int8_t *)  * abc->Kp);      hmm->rbv[0] = NULL;
  ESL_ALLOC(hmm->rwv, sizeof(int16_t *) * abc->Kp);      hmm->rwv[0] = NULL;
  ESL_ALLOC(hmm->twv, sizeof(int16_t *) * (hmm->Qw+1));  hmm->twv[0] = NULL;
  ESL_ALLOC(hmm->rfv, sizeof(float *)   * abc->Kp);      hmm->rfv[0] = NULL;
  ESL_ALLOC(hmm->tfv, sizeof(float *)   * (hmm->Qf+1));  hmm->tfv[0] = NULL;

  hmm->rbv[0] = esl_alloc_aligned(abc->Kp * (hmm->Qb + h4_EXTRA_SB) * hmm->V, hmm->V);      // allocation size in bytes; V is # of bytes per vector
  hmm->rwv[0] = esl_alloc_aligned(abc->Kp * hmm->Qw                 * hmm->V, hmm->V);
  hmm->rfv[0] = esl_alloc_aligned(abc->Kp * hmm->Qf                 * hmm->V, hmm->V);
  hmm->twv[0] = esl_alloc_aligned( h4_NVT * (hmm->Qw+1)             * hmm->V, hmm->V);
  hmm->tfv[0] = esl_alloc_aligned( h4_NVT * (hmm->Qf+1)             * hmm->V, hmm->V);

  for (x = 1; x < abc->Kp; x++) {
    hmm->rbv[x] = hmm->rbv[0] + (x * (hmm->V / sizeof(int8_t))  * (hmm->Qb + h4_EXTRA_SB)); // ptr arithmetic, though, is in # of *values* per row, not bytes: Vb/Vw/Vf, not V
    hmm->rwv[x] = hmm->rwv[0] + (x * (hmm->V / sizeof(int16_t)) * hmm->Qw);                 
    hmm->rfv[x] = hmm->rfv[0] + (x * (hmm->V / sizeof(float))   * hmm->Qf);                 
  }

  /* vectorized transition scores have a special layout; 
   * tv[q=0..Q-1] point to arrays of MM..DI, all but DD's;
   * tv[Q] points to array of Q DD's all by themselves 
   */
  for (q = 1; q <= hmm->Qw; q++)  hmm->twv[q] = hmm->twv[0] + (q * (hmm->V / sizeof(int16_t)) * (h4_NVT-1));
  for (q = 1; q <= hmm->Qf; q++)  hmm->tfv[q] = hmm->tfv[0] + (q * (hmm->V / sizeof(float))   * (h4_NVT-1));

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

  hmm->M     = M;
  hmm->abc   = abc;
  hmm->flags = 0;
  return eslOK;

 ERROR:
  h4_profile_Destroy(hmm);
  return eslEMEM;
}


/* Function:  h4_profile_Clone()
 * Synopsis:  Duplicate a profile into newly allocated space.
 * Incept:    SRE, Fri 05 Apr 2019
 *
 * Purpose:   Make a duplicate of <hmm> in freshly allocated
 *            space. Return a pointer to the new profile.
 *
 * Throws:    <NULL> on allocation failure.
 */
H4_PROFILE *
h4_profile_Clone(const H4_PROFILE *hmm)
{
  H4_PROFILE *hmm2 = NULL;

  if (( hmm2 = h4_profile_Create(hmm->abc, hmm->M)) == NULL) return NULL;
  h4_profile_Copy(hmm, hmm2);
  return hmm2;
}


/* Function:  h4_profile_Copy()
 * Synopsis:  Copy a profile into existing space.
 * Incept:    SRE, Fri 05 Apr 2019
 *
 * Purpose:   Copy profile <src> to <dst>, where space for
 *            <dst> is already allocated.
 *            
 *            Currently assumes that <dst->M> is exactly the same as
 *            <src->M>: the profiles are exactly the same allocated
 *            size. We might relax this in the future, allowing there
 *            to be a difference between allocated size and used size,
 *            and then we'd just need <src->M> <= <dst->allocM>.
 *            
 *            Also assumes that <dst->abc> is the same type as
 *            <src->abc>, so that they're identical and their sizes
 *            match. Usually, the ptrs themselves will be identical -
 *            pointing to the same one alphabet structure in the
 *            caller - but this isn't necessary, they could be two
 *            separately allocated but identical alphabets.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
h4_profile_Copy(const H4_PROFILE *src, H4_PROFILE *dst)
{
  ESL_DASSERT1(( src->M  == dst->M ));
  ESL_DASSERT1(( src->abc->type == dst->abc->type ));  // alphabet reference ptr is already in <dst>, matching <src>
  ESL_DASSERT1(( src->V  == dst->V ));
  ESL_DASSERT1(( src->Qb == dst->Qb ));
  ESL_DASSERT1(( src->Qw == dst->Qw ));
  ESL_DASSERT1(( src->Qf == dst->Qf ));

  /* null (background) residue freqs  */
  esl_vec_FCopy(src->f,                 src->abc->K, dst->f);

  /* probability model */
  esl_mat_FCopy(src->t,   src->M+1,     h4_NT,       dst->t);
  esl_mat_FCopy(src->e,   src->M+1,     src->abc->K, dst->e);

  /* standard scores */
  esl_mat_FCopy(src->tsc, src->M+1,     h4_NTSC,     dst->tsc);
  esl_mat_FCopy(src->rsc, src->abc->Kp, src->M+1,    dst->rsc);

  /* vectorized scores: SSV */
  esl_mat_BCopy(src->rbv, src->abc->Kp, src->Qb*src->V, dst->rbv);
  dst->tauBM = src->tauBM;

  /* vectorized scores: VF */
  esl_mat_WCopy(src->rwv,    src->abc->Kp, src->Qw*src->V/2, dst->rwv);
  esl_vec_WCopy(src->twv[0], src->Qw*h4_NVT*src->V/2,        dst->twv[0]);
  dst->ddbound_w = src->ddbound_w;
  
  /* vectorized scores: FB */
  esl_mat_FCopy(src->rfv, src->abc->Kp, src->Qf*src->V/4, dst->rfv);
  esl_vec_FCopy(src->tfv[0], src->Qf*h4_NVT*src->V/4,     dst->tfv[0]);

  dst->flags = src->flags;
  return eslOK;
}



/* Function:  h4_profile_Sizeof()
 * Synopsis:  Returns allocated size of a profile, in bytes.
 */
size_t
h4_profile_Sizeof(const H4_PROFILE *hmm)
{
  int    maxQb = H4_Q(hmm->M, h4_VMAX_SSV);   
  int    maxQw = H4_Q(hmm->M, h4_VMAX_VF);
  int    maxQf = H4_Q(hmm->M, h4_VMAX_FB);
  size_t n     = 0;

  n += sizeof(H4_PROFILE);
  n += esl_mat_FSizeof(hmm->M+1,       h4_NT);            // t
  n += esl_mat_FSizeof(hmm->M+1,       hmm->abc->K);      // e
  n += sizeof(float) * hmm->abc->K;                       // f
  n += esl_mat_FSizeof(hmm->M+1,       h4_NTSC);          // tsc
  n += esl_mat_FSizeof(hmm->abc->Kp+1, hmm->M+1);         // rsc 
  n += sizeof(int8_t *)  * hmm->abc->Kp;                  // rbv ptrs
  n += sizeof(int16_t *) * hmm->abc->Kp;                  // rwv ptrs
  n += sizeof(float *)   * hmm->abc->Kp;                  // rfv ptrs
  n += sizeof(int16_t *) * (maxQw+1);                     // twv ptrs
  n += sizeof(float *)   * (maxQf+1);                     // tfv ptrs
  n += hmm->abc->Kp * h4_VWIDTH * (maxQb + h4_EXTRA_SB);  // rbv[0]
  n += hmm->abc->Kp * h4_VWIDTH *  maxQw;                 // rwv[0]
  n += hmm->abc->Kp * h4_VWIDTH *  maxQf;                 // rfv[0]
  n += h4_NVT       * h4_VWIDTH * (maxQw + 1);            // twv[0]
  n += h4_NVT       * h4_VWIDTH * (maxQf + 1);            // tfv[0]
  
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



/*****************************************************************
 * 2. Model estimation: counts to probabilities.
 *****************************************************************/

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


/* Function:  h4_profile_Renormalize()
 * Synopsis:  Renormalize profile probability params.
 * Incept:    SRE, Fri 05 Apr 2019
 *
 * Purpose:   Renormalize the probability parameters in profile <hmm>,
 *            i.e. transitions <hmm->t> and emissions <hmm->e>.
 *
 *            Assures that the fixed boundary conditions are set properly.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
h4_profile_Renormalize(H4_PROFILE *hmm)
{
  float sum;
  int   k;

  for (k = 1; k <= hmm->M; k++) esl_vec_FNorm(hmm->e[k], hmm->abc->K); // M0 doesn't exist. Keep at conventions.
  for (k = 1; k <  hmm->M; k++) esl_vec_FNorm(hmm->t[k],   3); // M states. M0 doesn't exist, Mm forced to Mm->E. Keep them untouched, at conventions.
  for (k = 1; k <  hmm->M; k++) esl_vec_FNorm(hmm->t[k]+3, 3); // I states. I0, IM don't exist. Keep untouched, at conventions.
  for (k = 1; k <  hmm->M; k++) esl_vec_FNorm(hmm->t[k]+6, 3); // D states. D0 doesn't exist, Dm is forced to Dm->E. Keep untouched, at conventions.

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


/* Function:  h4_profile_Occupancy()
 * Synopsis:  Calculate match occupancy and insert expected use count vectors.
 *
 * Purpose: Calculate the expected number of times that each match and
 *            insert state is used in a single domain subpath
 *            (B->...->E) in profile <hmm>. (Since each match state
 *            can only be used once, this is also the probability of
 *            using the state. Inserts can be used more than once.)
 *            
 *            Caller can elect to get back either the total expected
 *            number of match or insert states used (in <*opt_mtot>
 *            and <*opt_itot>), and/or provide allocated space of
 *            at least (M+1) floats for <mocc> and <iocc> arrays
 *            to get expected counts per state.
 *
 *            <mocc[0]>, <iocc[0]>, and <iocc[M]> are 0 (there's no
 *            M0, I0, or Im state).
 *            
 *            Only depends on transition probabilities, not emissions.
 *            It's ok to abuse this function, calling it on a
 *            partially parameterized <hmm> that only has transition
 *            probs.  Some <h4_modelsample_*()> functions do this as a
 *            speed thing. Therefore, it deliberately doesn't check
 *            the model's <h4_HASPROBS> flag.
 * 
 *            <mtot>+<itot> is the expected length of one domain.
 *            
 * Returns:   <eslOK> on success.
 */
int
h4_profile_Occupancy(const H4_PROFILE *hmm, float *mocc, float *iocc, float *opt_mtot, float *opt_itot)
{
  float mprv, mcur, icur;
  float mtot = 0.0;
  float itot = 0.0;
  int   k;
  float mshare;

  if (mocc) mocc[0] = 0.;
  if (iocc) iocc[0] = 0.;

  for (k = 1; k <= hmm->M; k++)
    {
      if (k == 1)
	mcur = hmm->t[0][h4_TMM];
      else {
	mshare = hmm->t[k-1][h4_TIM] / (hmm->t[k-1][h4_TIM] + hmm->t[k-1][h4_TID]);  // factors out II
	mcur   =     mprv  * hmm->t[k-1][h4_TMM];           // M->M
	mcur  += (1.-mprv) * hmm->t[k-1][h4_TDM];           // D->M
	mcur  +=     mprv  * hmm->t[k-1][h4_TMI]  * mshare; // M->I->M
	mcur  += (1.-mprv) * hmm->t[k-1][h4_TDI]  * mshare; // D->I->M
	mcur   = ESL_MIN(1.0, mcur);                        // avoid floating pt roundoff error making mocc[k] 1+epsilon
      }
      if (mocc) mocc[k] = mcur;
      mtot += mcur;
      mprv  = mcur;

      if (k < hmm->M) {
	icur  =     mcur  * hmm->t[k][h4_TMI];  // M->I
	icur += (1.-mcur) * hmm->t[k][h4_TDI];  // D->I
	icur /= (1.-hmm->t[k][h4_TII]);         // account for expected I length 
      } else icur = 0.f;
      if (iocc) iocc[k] = icur;
      itot += icur;
    }
  
  if (opt_mtot) *opt_mtot = mtot;
  if (opt_itot) *opt_itot = itot;
  return eslOK;
}


/*****************************************************************
 * 3. Profile configuration: probabilities to scores
 *****************************************************************/




/*****************************************************************
 * 4. Debugging and development tools
 *****************************************************************/


/* Function:  h4_profile_Dump()
 * Synopsis:  Dump contents of an H4_PROFILE for inspection
 * Incept:    SRE, Mon 16 Jul 2018 [Benasque]
 */
int
h4_profile_Dump(FILE *fp, H4_PROFILE *hmm)
{
  int k,a,z;

  fprintf(fp, "Emission probabilities:\n");
  fprintf(fp, "     ");
  for (a = 0; a < hmm->abc->K; a++)
    fprintf(fp, "         %c%c", hmm->abc->sym[a], a == hmm->abc->K-1 ? '\n':' ');
  for (k = 1; k <= hmm->M; k++)
    {
      fprintf(fp, "%4d ", k);
      for (a = 0; a < hmm->abc->K; a++)
	fprintf(fp, "%10.4f%c", hmm->e[k][a], a == hmm->abc->K-1 ? '\n':' ');
    }

  fprintf(fp, "Transition probabilities:\n");
  fprintf(fp, "     ");
  fprintf(fp, "%10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
	  "TMM", "TMI", "TMD", "TIM", "TII", "TID", "TDM", "TDI", "TDD");
  for (k = 0; k < hmm->M; k++)  // include t[0] because GM, GD (in MM, MD) are data-dependent
    {                           // exclude t[M] which has no data dependent transitions.
      fprintf(fp, "%4d ", k);
      for (z = 0; z < h4_NT; z++)
	fprintf(fp, "%10.4f%c", hmm->t[k][z], z == h4_NT-1 ? '\n':' ');
    }

  fprintf(fp, "Emission scores:\n");
  fprintf(fp, "     ");
  for (a = 0; a < hmm->abc->K; a++)
    fprintf(fp, "         %c%c", hmm->abc->sym[a], a == hmm->abc->K-1 ? '\n':' ');
  for (k = 1; k <= hmm->M; k++)
    {
      fprintf(fp, "%4d ", k);
      for (a = 0; a < hmm->abc->K; a++)
	fprintf(fp, "%10.4f%c", hmm->rsc[a][k], a == hmm->abc->K-1 ? '\n':' ');
    }

  fprintf(fp, "Transition scores:\n");
  fprintf(fp, "     ");
  fprintf(fp, "%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
	  "MM", "IM", "DM", "LM", "GM", "MI", "II", "DI", "GI", "MD", "ID", "DD", "DGE");
  for (k = 0; k <= hmm->M; k++)  // include t[0] because GM, GD (in MM, MD) are data-dependent
    {             
      fprintf(fp, "%4d ", k);
      for (z = 0; z < h4_NTSC; z++)
	fprintf(fp, "%10.4f%c", hmm->tsc[k][z], z == h4_NTSC-1 ? '\n':' ');
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
  int   M   = hmm->M;
  float tol = 1e-4;

  if (M < 1)      ESL_FAIL(eslFAIL, errbuf, "invalid model size M");
  if (! hmm->abc) ESL_FAIL(eslFAIL, errbuf, "no model alphabet");

  /* Probability model */
  if (hmm->flags & h4_HASPROBS)
    {
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
	ESL_FAIL(eslFAIL, errbuf, "something awry in transition probability edge conventions");
    }


  /* Score model */
  if (hmm->flags & h4_HASBITS)
    {
      /* Check edge conventions on scores; see h4_profile.md */
      if (hmm->tsc[0][h4_IM]  != 0.           ||  // [0][MM] is log2(G->M1)
	  hmm->tsc[0][h4_DM]  != 0.           ||  
	  hmm->tsc[0][h4_MI]  != -eslINFINITY ||  // [0][LM] and [0][GM] are set
	  hmm->tsc[0][h4_II]  != -eslINFINITY ||
	  hmm->tsc[0][h4_DI]  != -eslINFINITY ||
	  hmm->tsc[0][h4_GI]  != -eslINFINITY ||
	  hmm->tsc[0][h4_ID]  != -eslINFINITY ||  // [0][MD] is log2(G->D1)
	  hmm->tsc[0][h4_DD]  != -eslINFINITY ||
	  hmm->tsc[0][h4_DGE] != -eslINFINITY)
	ESL_FAIL(eslFAIL, errbuf, "something awry in transition score conventions at k=0");

      if (M > 1 && hmm->tsc[M-1][h4_DGE] != 0.)  // if M=1, then tsc[M-1] = tsc[0], which we just checked above.
	ESL_FAIL(eslFAIL, errbuf, "something awry in transition score conventions at k=M-1");

      if (hmm->tsc[M][h4_MM]  != 0. ||
	  hmm->tsc[M][h4_IM]  != 0. ||
	  hmm->tsc[M][h4_DM]  != 0. ||
	  hmm->tsc[M][h4_LM]  != -eslINFINITY ||  // because LM, GM are stored off by one, they're -inf at k=M
	  hmm->tsc[M][h4_GM]  != -eslINFINITY ||
	  hmm->tsc[M][h4_MI]  != -eslINFINITY ||
	  hmm->tsc[M][h4_II]  != -eslINFINITY ||
	  hmm->tsc[M][h4_GI]  != -eslINFINITY ||
	  hmm->tsc[M][h4_MD]  != -eslINFINITY ||
	  hmm->tsc[M][h4_ID]  != -eslINFINITY ||
	  hmm->tsc[M][h4_DD]  != -eslINFINITY ||
	  hmm->tsc[M][h4_DGE] != 0.)
	ESL_FAIL(eslFAIL, errbuf, "something awry in transition score conventions at k=M");
    }

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


/* Function:  h4_profile_MutePathScore()
 * Synopsis:  Calculate log2 prob of G->D1..Dm->E path
 * Incept:    SRE, Sat 18 May 2019
 *
 * Purpose:   Calculate the log2 prob (i.e. in bits) of the mute glocal
 *            path G->D1..Dm->E, the probability mass that
 *            wing-retracted DP algorithms neglect. We sometimes need
 *            to add this term back, especially in some unit tests
 *            that check for perfect summation over paths or sequences.
 * 
 * Xref:      [SRE:J9/98,100]
 */
int
h4_profile_MutePathScore(const H4_PROFILE *hmm, float *ret_sc)
{
  int   k;
  float sc;
    
  sc = hmm->tsc[0][h4_MD];    // G->D1
  for (k = 1; k < hmm->M; k++)
    sc += hmm->tsc[k][h4_DD]; // D1..Dm
  *ret_sc = sc;
  return eslOK;
}


/* Function:  h4_profile_SameAsVF()
 * Synopsis:  Scale and round standard profile scores to match VF
 * Incept:    SRE, Sun 07 Jul 2019
 *
 * Purpose:   Make a copy of <hmm> with its standard scores scaled and rounded 
 *            to match Viterbi filter scores, by scaling by a
 *            factor of <h4_SCALE_W> and rounding to the nearest
 *            integer. Return the copy in <*ret_xhmm>.
 *            
 *            Some unit tests use this, for example, when testing that
 *            a local-only standard model alignment score matches a
 *            Viterbi filter score.
 * 
 *            To convert a standard raw Viterbi score <vsc> calculated
 *            with a "SameAsVF" profile to a raw bit score that should
 *            match Viterbi filter score (as long as the VF filter
 *            score doesn't overflow), do <(vsc / h4_SCALE_W) -
 *            h4_3NAT_APPROX>.
 */
int
h4_profile_SameAsVF(H4_PROFILE *hmm, H4_PROFILE **ret_xhmm)
{
  H4_PROFILE *xhmm = NULL;
  int k,s,x;
  int status;

  ESL_DASSERT1(( hmm->flags & h4_HASBITS ));
  
  if (( xhmm = h4_profile_Clone(hmm)) == NULL) { status = eslEMEM; goto ERROR; }

  for (x = 0; x < hmm->abc->Kp; x++)
    for (k = 0; k <= hmm->M; k++)
      xhmm->rsc[x][k] = (hmm->rsc[x][k] <= -eslINFINITY) ? -eslINFINITY : roundf(h4_SCALE_W * hmm->rsc[x][k]);

  for (k = 0; k <= hmm->M; k++)
    for (s = 0; s < h4_NTSC; s++)
      xhmm->tsc[k][s] = (hmm->tsc[k][s] <= -eslINFINITY) ? -eslINFINITY : roundf(h4_SCALE_W * hmm->tsc[k][s]);

  *ret_xhmm = xhmm;
  return eslOK;

 ERROR:
  h4_profile_Destroy(xhmm);
  *ret_xhmm = NULL;
  return status;
}

/*****************************************************************
 * 5. Unit tests
 *****************************************************************/
#ifdef h4PROFILE_TESTDRIVE

#include "esl_sq.h"

#include "h4_counts.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "emit.h"
#include "modelsample.h"

/* utest_occupancy()
 * 
 * The "occupancy" unit test checks that the analytical calculation
 * of match and insert state occupancy in h4_profile_Occupancy()
 * agrees with counts in a brute force simulation of <N> emitted
 * sequences.
 */
static void
utest_occupancy(FILE *diagfp, ESL_RANDOMNESS *rng, int alphatype, int M, int N)
{
  char           msg[]    = "h4_profile occupancy unit test failed";
  ESL_ALPHABET  *abc      = esl_alphabet_Create(alphatype);
  ESL_SQ        *sq       = esl_sq_CreateDigital(abc);
  H4_PROFILE    *hmm      = NULL;
  H4_COUNTS     *ctm      = h4_counts_Create(abc, M);
  H4_MODE       *mo       = h4_mode_Create();
  H4_PATH       *pi       = h4_path_Create();
  float         *mocc     = malloc(sizeof(float) * (M+1));
  float         *iocc     = malloc(sizeof(float) * (M+1));
  float         *sim_mocc = malloc(sizeof(float) * (M+1));
  float         *sim_iocc = malloc(sizeof(float) * (M+1));
  float          mtot, itot;
  double         exp_len  = 0.;
  double         sim_len  = 0.;
  int            k,idx;

  if (diagfp)
    esl_dataheader(stdout, 6, "k", 12, "mocc",   12, "sim_mocc", 12, "mocc_diff", 
		                   12, "iocc",   12, "sim_iocc", 12, "iocc_diff",
		                    9, "exp_len", 9, "sim_len",   9, "len_diff", 0); 
 

  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  if ( h4_mode_SetUniglocal(mo)          != eslOK) esl_fatal(msg);
  if ( h4_mode_SetLength(mo, 0)          != eslOK) esl_fatal(msg);
  //  h4_profile_Dump(stdout, hmm);
  
  for (idx = 0; idx < N; idx++)
    {
      if ( h4_emit(rng, hmm, mo, sq, pi)        != eslOK) esl_fatal(msg);
      if ( h4_path_Count(pi, sq->dsq, 1.0, ctm) != eslOK) esl_fatal(msg);
      sim_len += sq->n;
      h4_path_Reuse(pi);
      esl_sq_Reuse(sq);
    }
  sim_len /= (double) N;

  sim_mocc[0] = sim_iocc[0] = 0.;
  for (k = 1; k <= M; k++)
    {
      sim_mocc[k] = esl_vec_DSum(ctm->e[k], ctm->abc->K) / (double) N;
      if (k == M) sim_iocc[0] = 0.;
      else        sim_iocc[k] = (ctm->t[k][h4_TMI] + ctm->t[k][h4_TII] + ctm->t[k][h4_TDI]) / (double) N;
    }

  if ( h4_profile_Occupancy(hmm, mocc, iocc, &mtot, &itot) != eslOK) esl_fatal(msg);
  exp_len = mtot+itot;
  
  if (diagfp)
    {
      for (k = 1; k <= M; k++)
	fprintf(stdout, "%6d %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %9.4f %9.4f %9.4f\n",
		k, mocc[k], sim_mocc[k], (mocc[k]-sim_mocc[k]) / sim_mocc[k],
		iocc[k], sim_iocc[k], (iocc[k]-sim_iocc[k]) / sim_iocc[k],
		exp_len, sim_len, (exp_len-sim_len)/sim_len);
    }

  if (esl_DCompareNew(sim_len, exp_len, 0.01, 0.0) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++) if (esl_DCompareNew(sim_mocc[k], mocc[k], 0.05, 0.0) != eslOK) esl_fatal(msg);
  for (k = 1; k <  M; k++) if (esl_DCompareNew(sim_iocc[k], iocc[k], 0.10, 0.0) != eslOK) esl_fatal(msg);

  h4_path_Destroy(pi);
  h4_counts_Destroy(ctm);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
}


#endif // h4PROFILE_TESTDRIVE

/*****************************************************************
 * 6. Test driver
 *****************************************************************/
#ifdef h4PROFILE_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                          docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help summary",                  0 },
  { "-s",         eslARG_INT,    "42", NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed",         0 },
  { "-M",         eslARG_INT,    "20", NULL, NULL,  NULL,  NULL, NULL, "set test profile length",                  0 },
  { "-N",         eslARG_INT,"100000", NULL, NULL,  NULL,  NULL, NULL, "set # of sampled paths in occupancy test", 0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv, "test driver for h4_profile", "[-options]");
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  int             M   = esl_opt_GetInteger(go, "-M");
  int             N   = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng)); 

  utest_occupancy(stdout, rng, eslRNA, M, N);
  
  fprintf(stderr, "#  status   = ok\n");

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}  
#endif // h4PROFILE_TESTDRIVE
