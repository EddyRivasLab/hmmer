/* Model configuration: 
 * Converting a core model to a fully configured Plan7 search profile.
 * 
 * Contents:
 *   1. Profile configuration
 *   2. Internal functions
 *   3. Unit tests
 *   4. Test driver
 *   5. References
 *   6. Copyright and license information
 */
#include "p7_config.h"

#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "base/p7_bg.h"
#include "base/p7_hmm.h"
#include "base/p7_profile.h"

#include "search/modelconfig.h"


static int set_local_entry (P7_PROFILE *gm, const P7_HMM *hmm);
static int set_glocal_entry(P7_PROFILE *gm, const P7_HMM *hmm);
static int set_glocal_exit (P7_PROFILE *gm, const P7_HMM *hmm);


/*****************************************************************
 * 1. Profile configuration 
 *****************************************************************/
 
/* Function:  p7_profile_Config()
 * Synopsis:  Standard HMMER3 profile configuration.
 *
 * Purpose:   Parameterize search profile <gm>, from core
 *            probability model <hmm> and null model <bg>,
 *            using the standard HMMER3 configuration.
 *            
 *            The caller assures that <gm> has been allocated for a
 *            model of at least <hmm->M> nodes, and that the alphabet
 *            <gm->abc> is the same as <hmm->abc>.
 *
 *            The standard configuration sets <L=0> (no length model
 *            yet; caller will do <p7_profile_SetLength(gm, L)>),
 *            <nj=1.0> (multihit), and <pglocal=0.5> (dual-mode
 *            local/glocal).
 *
 * Returns:   <eslOK> on success; the profile <gm> now contains 
 *            scores. Before it is used, it still needs a length
 *            model; see <p7_profile_SetLength()>.
 *            
 * Throws:    <eslEMEM> on allocation error.           
 */
int
p7_profile_Config(P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg)
{
  return p7_profile_ConfigCustom(gm, hmm, bg, 0, 1.0, 0.5); 
  /* that's L=0         : caller will call p7_profile_SetLength() before first use;
   *        nj=1.0      : multihit;
   *        pglocal=0.5 : dual-mode local/glocal.
   */
}


int
p7_profile_ConfigLocal(P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L)
{
  return p7_profile_ConfigCustom(gm, hmm, bg, L, 1.0, 0.0); /* nj=1, pglocal=0: multihit local */
}

int
p7_profile_ConfigUnilocal(P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L)
{
  return p7_profile_ConfigCustom(gm, hmm, bg, L, 0.0, 0.0); /* nj=0, pglocal=0: unihit local */
}

int
p7_profile_ConfigGlocal(P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L)
{
  return p7_profile_ConfigCustom(gm, hmm, bg, L, 1.0, 1.0); /* nj=1, pglocal=1: multihit glocal */
}

int
p7_profile_ConfigUniglocal(P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L)
{
  return p7_profile_ConfigCustom(gm, hmm, bg, L, 0.0, 1.0); /* nj=0, pglocal=1: unihit glocal */
}

/* Function:  p7_profile_ConfigCustom()
 * Synopsis:  Configure a search profile from a core HMM.
 *
 * Purpose:   Parameterize search model <gm>, given a model <hmm> with
 *            core probabilities, the null1 model <bg>, and the free
 *            parameters <L>, <nj>, <pglocal> that determine search
 *            profile configuration (see below).
 *            
 *            The caller assures that <gm> has been allocated for a
 *            model of at least <hmm->M> nodes, and that the alphabet
 *            <gm->abc> is the same as <hmm->abc>.
 *            
 *            Parameter <L> controls the parameters of the target
 *            sequence length model. The length model must be (re)set
 *            for each new target sequence length. The function
 *            <p7_profile_SetLength()> is the fastest way to
 *            reconfigure the length model without reconfiguring the
 *            entire profile. It is common to initially config with
 *            <L=0>, knowing that <p7_profile_SetLength()> will also
 *            be called appropriately before the profile is first
 *            used.
 *            
 *            Parameter <nj> controls the unihit/multihit behavior of
 *            the profile. <nj> is the expected number of J segments.
 *            <nj=0> is unihit. <nj=1> is the standard multihit model,
 *            with <t_{EJ} = 0.5>.
 *            
 *            Parameter <pglocal> controls the local/glocal mode of
 *            the profile, by setting the <t_{BG}> transition
 *            probability. <pglocal=0> makes a local alignment
 *            profile.  <pglocal=1> makes a glocal alignment
 *            profile. <pglocal=0.5> makes the (now-standard)
 *            "dual-mode" profile.
 *            
 * Returns:   <eslOK> on success; the profile <gm> now contains 
 *            scores and is ready for searching target sequence(s)
 *            of length <L>.
 *            
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_profile_ConfigCustom(P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L, float nj, float pglocal)
{
  float *tp, *rp;
  float  sc[p7_MAXCODE];
  int    k, x, z;	/* counters over states, residues, annotation */
  int    status;

  /* Contract checks */
  if (gm->abc->type != hmm->abc->type) ESL_EXCEPTION(eslEINVAL, "HMM and profile alphabet don't match");
  if (hmm->M > gm->allocM)             ESL_EXCEPTION(eslEINVAL, "profile too small to hold HMM");
  if (! (hmm->flags & p7H_CONS))       ESL_EXCEPTION(eslEINVAL, "HMM must have a consensus to transfer to the profile");

  /* Copy annotation across from HMM  */
  if (gm->name) free(gm->name);   if ((status = esl_strdup(hmm->name,   -1, &(gm->name))) != eslOK) return status;
  if (gm->acc)  free(gm->acc);    if ((status = esl_strdup(hmm->acc,    -1, &(gm->acc)))  != eslOK) return status;
  if (gm->desc) free(gm->desc);   if ((status = esl_strdup(hmm->desc,   -1, &(gm->desc))) != eslOK) return status;

  if (hmm->flags & p7H_RF)    strcpy(gm->rf,        hmm->rf);
  if (hmm->flags & p7H_MMASK) strcpy(gm->mm,        hmm->mm);
  if (hmm->flags & p7H_CONS)  strcpy(gm->consensus, hmm->consensus); /* must be present, actually, so the flag test is just for symmetry w/ other optional HMM fields */
  if (hmm->flags & p7H_CS)    strcpy(gm->cs,        hmm->cs);

  for (z = 0; z < p7_NEVPARAM; z++) gm->evparam[z] = hmm->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) gm->cutoff[z]  = hmm->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) gm->compo[z]   = hmm->compo[z];

  gm->offs[p7_MOFFSET] = -1;
  gm->offs[p7_FOFFSET] = -1;
  gm->offs[p7_POFFSET] = -1;
  gm->roff             = -1;
  gm->eoff             = -1;

  gm->max_length       = hmm->max_length;
  gm->M                = hmm->M;

  gm->xsc[p7P_E][p7P_MOVE] = logf ( 1.0f / (1.0f + nj));       // E->C
  gm->xsc[p7P_E][p7P_LOOP] = logf ( nj / (1.0f + nj));         // E->J

  gm->xsc[p7P_B][0]        = logf(1.0 - pglocal);              // B->L
  gm->xsc[p7P_B][1]        = logf(pglocal);                    // B->G

  gm->xsc[p7P_G][0]        = logf(1.0 - hmm->t[0][p7H_MD]);    // not p7H_MM, because we're absorbing core model's I0 state and p7H_MI I
  gm->xsc[p7P_G][1]        = logf(hmm->t[0][p7H_MD]);

  /* Initialize tsc[k=M] boundary condition to -inf; 
   * this deletes Im state;
   * will overwrite as needed in set_ functions below
   */
  esl_vec_FSet( gm->tsc + gm->M*p7P_NTRANS, p7P_NTRANS, -eslINFINITY);

  /* Entry and exit distributions: L->MLk, G->MGk, MGk->E */
  if (( status = set_local_entry (gm, hmm)) != eslOK) return status;   // L->MLk params: uniform fragment length parameterization 
  if (( status = set_glocal_entry(gm, hmm)) != eslOK) return status;   // G->MGk params: left wing retraction : glocal entry
  if (( status = set_glocal_exit (gm, hmm)) != eslOK) return status;   // MGk->E params: right wing retraction : glocal exit

  /* Remaining transition scores (BLM, BGM, MGE entries/exits were
   * just set, above). 
   */
  for (k = 1; k < gm->M; k++) {
    tp = gm->tsc + k * p7P_NTRANS;
    tp[p7P_MM] = logf(hmm->t[k][p7H_MM]);
    tp[p7P_IM] = logf(hmm->t[k][p7H_IM]);
    tp[p7P_DM] = logf(hmm->t[k][p7H_DM]);
    tp[p7P_MD] = logf(hmm->t[k][p7H_MD]);
    tp[p7P_DD] = logf(hmm->t[k][p7H_DD]);
    tp[p7P_MI] = logf(hmm->t[k][p7H_MI]);
    tp[p7P_II] = logf(hmm->t[k][p7H_II]);
  }
  
  /* All transition scores at k=M are special cases for DP boundary conditions! */
  tp = gm->tsc + gm->M * (p7P_NTRANS);
  tp[p7P_MD] = 0.0;		/* glocal Forward and Backwards rely on this    */
  tp[p7P_DD] = 0.0;		/* ditto */

  /* match emission scores. */
  sc[hmm->abc->K]     = -eslINFINITY; /* gap character */
  sc[hmm->abc->Kp-2]  = -eslINFINITY; /* nonresidue character */
  sc[hmm->abc->Kp-1]  = -eslINFINITY; /* missing data character */
  for (k = 1; k <= hmm->M; k++) {
    for (x = 0; x < hmm->abc->K; x++) 
      sc[x] = logf(hmm->mat[k][x] / bg->f[x]);
    esl_abc_FExpectScVec(hmm->abc, sc, bg->f); 
    for (x = 0; x < hmm->abc->Kp; x++) {
      rp = gm->rsc[x] + k * p7P_NR;
      rp[p7P_M] = sc[x];
    }
  }
  
 /* Insert emission scores are currently hardwired to zero, an
  * implicit assumption that insert emissions are equal to
  * background. Benchmarking shows that setting inserts to informative
  * emission distributions causes more problems than it's worth: polar
  * biased composition hits driven by stretches of "insertion" occur,
  * and are difficult to correct for. This may change someday, if we
  * get a better null model for example.
  */
  for (x = 0; x < gm->abc->Kp; x++)
    {
      for (k = 1; k < hmm->M; k++) P7P_ISC(gm, k, x) = 0.0f;
      P7P_ISC(gm, hmm->M, x) = -eslINFINITY;   /* init I_M to impossible.   */
    }
  for (k = 1; k <= hmm->M; k++) P7P_ISC(gm, k, gm->abc->K)    = -eslINFINITY; /* gap symbol */
  for (k = 1; k <= hmm->M; k++) P7P_ISC(gm, k, gm->abc->Kp-2) = -eslINFINITY; /* nonresidue symbol */
  for (k = 1; k <= hmm->M; k++) P7P_ISC(gm, k, gm->abc->Kp-1) = -eslINFINITY; /* missing data symbol */

  /* Remaining specials, [NCJ][MOVE | LOOP] are set by length model */
  gm->nj      = nj;
  gm->pglocal = pglocal;
  return p7_profile_SetLength(gm, L);
}


/* Function:  p7_profile_SetLength()
 * Synopsis:  Set the target sequence length of a model.
 *
 * Purpose:   Given a model <gm> already configured for scoring, in some
 *            particular algorithm mode; reset the expected length
 *            distribution of the profile for a new mean of <L>.
 *
 *            This doesn't affect the length distribution of the null
 *            model. That must also be reset, using <p7_bg_SetLength()>.
 *            
 *            We want this routine to run as fast as possible, because
 *            the caller needs to dynamically reconfigure the model
 *            for the length of each target sequence in a database
 *            search.
 *
 * Returns:   <eslOK> on success; xsc[NCJ] scores are set here. These
 *            control the target length dependence of the model.
 */
int
p7_profile_SetLength(P7_PROFILE *gm, int L)
{
  float ploop, pmove;
  
  /* Configure N,J,C transitions so they bear L/(2+nj) of the total
   * unannotated sequence length L. 
   */
  pmove = (2.0f + gm->nj) / ((float) L + 2.0f + gm->nj); /* 2/(L+2) for unihit; 3/(L+3) for standard multihit */
  ploop = 1.0f - pmove;
  gm->xsc[p7P_N][p7P_LOOP] =  gm->xsc[p7P_C][p7P_LOOP] = gm->xsc[p7P_J][p7P_LOOP] = logf(ploop);
  gm->xsc[p7P_N][p7P_MOVE] =  gm->xsc[p7P_C][p7P_MOVE] = gm->xsc[p7P_J][p7P_MOVE] = logf(pmove);
  gm->L = L;
  return eslOK;
}
/*------------ end, initial profile config ----------------------*/



/*****************************************************************
 * 2. Internal functions
 *****************************************************************/

/* set_local_entry()
 * 
 * Local mode entry prob t_BLMk is approx. 2/(M(M+1)), with a
 * correction for match state occupancy probability [Eddy08]:
 *    L->Mk = occ[k] /( \sum_i occ[i] * (M-i+1))
 *    
 * and we store these params off-by-one, with tLMk stored in TSC(gm,
 * k-1, p7P_LM), for DP efficiency reasons.
 * 
 * We need space for an M+1 occ[k] array of match occupancy.  We save
 * a malloc by using gm->rsc[0]'s space, which we know is >= M+1, and
 * which we guarantee (by order of calls in config calls) has not been
 * parameterized yet.
 * 
 * If the <hmm>'s <p7H_SINGLE> flag is up, <occ[1..M]> is 1.0:
 * i.e. the match state occupancy term is only applied to profiles,
 * not to single-seq queries. 
 */
static int
set_local_entry(P7_PROFILE *gm, const P7_HMM *hmm)
{
  float *occ = gm->rsc[0];	/* a safe, malloc-saving hack; see note above  */
  float  Z   = 0.;
  int    k;
  int    status;

  if (hmm->flags & p7H_SINGLE)
    {
      for (k = 0; k < hmm->M; k++) 
	P7P_TSC(gm, k, p7P_LM) = logf( 2. / ((float) hmm->M * (float) (hmm->M+1))); 
    }
  else 
    {
      if (( status  = p7_hmm_CalculateOccupancy(hmm, occ, NULL)) != eslOK) return status; /* NULL=iocc[k], I state expected uses, which we don't need here */
      for (k = 1; k <= hmm->M; k++) 
	Z += occ[k] * (float) (hmm->M-k+1);
      for (k = 1; k <= hmm->M; k++) 
	P7P_TSC(gm, k-1, p7P_LM) = logf(occ[k] / Z); /* note off-by-one: entry at Mk stored as [k-1][LM] */
      /* leave tsc[M,LM] as -inf, as we'd already initialized it */
    }
  return eslOK;
}

/* set_glocal_entry()
 * 
 * glocal entry is "wing retracted" into a G->Mk entry distribution.
 *
 * Wing retraction removes the B->G->D1..Dm->E->J->B mute cycle, by
 * removing the D1 state. In a profile, all entries start at an Mk
 * state. As a result, a (usually negligible) amount of probability is
 * unaccounted for (see p7_profile_GetMutePathLogProb(), if you need
 * to know what that neglected mass is).
 *
 * Unretracted "base" transition parameters G->{M1,D1} are in
 * xsc[p7G][0=M,1=D].
 *
 * We compute the wing-retracted entries G->{D1..Dk-1}->Mk as follows:
 *      tGM1 = log t(G->M1) 
 *      tGMk = log t(G->D1) + \sum_j={1..k-2} log t(Dj->Dj+1) + log t(Dk-1->Mk)
 * and for technical/efficiency reasons, these params are stored
 * off-by-one in the profile: tGMk is stored at TSC(k-1, p7P_GM).
 */
static int
set_glocal_entry(P7_PROFILE *gm, const P7_HMM *hmm)
{
  float Z;   
  int   k;

  P7P_TSC(gm, 0, p7P_GM) = gm->xsc[p7P_G][0];  // i.e. G->M1, stored at TSC(0,GM)
  Z                      = gm->xsc[p7P_G][1];  // i.e. G->D1
  for (k = 1; k < hmm->M; k++) 
    {
      P7P_TSC(gm, k, p7P_GM) = Z + logf(hmm->t[k][p7H_DM]);
      Z += logf(hmm->t[k][p7H_DD]);
    }
  /* leave tsc[M,GM] as -inf, as we'd already initialized it */
  return eslOK;
}

/* set_glocal_exit()
 * 
 * Right wing retraction, DGE
 *   TSC(k,DGE) = t(Dk+1->...Dm->E) 
 *              = [\prod_j=k+1..m-1 t(Dj->Dj+1)] * Dm->E
 *              = \prod_j=k+1..m-1 t(Dj->Dj+1)              | because Dm->E=1.0
 *  valid for k=0..M
 *  boundaries: TSC(M,DGE)   = 0
 *              TSC(M-1,DGE) = 0
 * note off by one. 
 * to get the glocal exit path from an Mk: t(k,MD) + t(k,DGE)
 * to get the glocal exit path from a Dk:  t(k,DD) + t(k,DGE)
 */
static int
set_glocal_exit(P7_PROFILE *gm, const P7_HMM *hmm)
{
  float Z = 0.0;		/* accumulated Dk+1..Dm part of the path */
  int   k;

  P7P_TSC(gm, gm->M, p7P_DGE) = Z;
  for (k = hmm->M-1; k >= 1; k--)
    {
      P7P_TSC(gm, k, p7P_DGE) = Z;
      Z += logf(hmm->t[k][p7H_DD]);
    }
  return eslOK;
}
/*--------------- end, internal functions -----------------------*/




/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7MODELCONFIG_TESTDRIVE

/* The config test simply makes sure a random profile passes
 * a Validate() check.
 */
static void
utest_config(P7_HMM *hmm, P7_BG *bg)
{
  char       *msg = "modelconfig.c::p7_profile_Config() unit test failed";
  P7_PROFILE *gm  = NULL;

  if ((gm = p7_profile_Create(hmm->M, hmm->abc))    == NULL)   esl_fatal(msg);
  if (p7_profile_Config(gm, hmm, bg)                != eslOK)  esl_fatal(msg);
  if (p7_profile_Validate(gm, NULL, 0.0001)         != eslOK)  esl_fatal(msg);

  p7_profile_Destroy(gm);
  return;
}

/* Note that calculate_occupancy has moved to p7_hmm.c, but
 * unit tests over there aren't hooked up yet; so leave a copy of the unit test 
 * here for now.
 */
static void
utest_occupancy(P7_HMM *hmm)
{
  char  *msg = "modelconfig.c::calculate_occupancy() unit test failed";
  float *occ;
  float  x;

  occ = malloc(sizeof(float) * (hmm->M+1));
  p7_hmm_CalculateOccupancy(hmm, occ, NULL);
  x = esl_vec_FSum(occ+1, hmm->M) / (float) hmm->M;
  if (esl_FCompare(x, 0.6, 0.1) != eslOK)           esl_fatal(msg);
  free(occ);
  return;
}
#endif /*p7MODELCONFIG_TESTDRIVE*/



/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7MODELCONFIG_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
   /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                            0},
  {"-s",  eslARG_INT,       "0", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for modelconfig.c";

int
main(int argc, char **argv)
{  
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc    = NULL;
  P7_HMM         *hmm    = NULL;
  P7_BG          *bg     = NULL;
  int             M      = 10000;
  
  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create amino alphabet");
  if (p7_hmm_Sample(rng, M, abc, &hmm)      != eslOK) esl_fatal("failed to sample random HMM");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to created null model");

  utest_config(hmm, bg);
  utest_occupancy(hmm);

  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);

  fprintf(stderr, "#  status = ok\n");
  return eslOK;
}
#endif /*p7MODELCONFIG_TESTDRIVE*/


/*****************************************************************
 * 5. References.
 *****************************************************************/

/* Major revisions include:
 *   May 2005:  xref STL9/77-81         Uniform fragment distribution
 *   Sep 2005:  xref STL10/24-26        Target length dependency
 *   Jan 2007:  xref STL11/125,136-137  HMMER3
 *   Jul 2007:  xref J1/103             floating point ops
 *   Nov 2011:  xref J9/2               dual-mode local/glocal
 *   
 */


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
