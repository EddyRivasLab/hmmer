/* Model configuration: 
 * Converting a core model to a fully configured Plan7 search profile.
 * 
 * Contents:
 *   1. Profile configuration
 *   2. Internal functions
 *   3. Unit tests
 *   4. Test driver
 *   5. Statistics collection driver
 *   6. References
 *   7. Copyright and license information
 */
#include "p7_config.h"

#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

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
  return p7_profile_ConfigCustom(gm, hmm, bg, L, 1.0, 1.0); /* nj=1, pglocal=0: multihit glocal */
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
 */
static int
set_local_entry(P7_PROFILE *gm, const P7_HMM *hmm)
{
  float *occ = gm->rsc[0];	/* a safe, malloc-saving hack; see note above  */
  float  Z   = 0.;
  int    k;
  int    status;

  if (( status  = p7_hmm_CalculateOccupancy(hmm, occ, NULL)) != eslOK) return status; /* NULL=iocc[k], I state expected uses, which we don't need here */
  for (k = 1; k <= hmm->M; k++) 
    Z += occ[k] * (float) (hmm->M-k+1);
  for (k = 1; k <= hmm->M; k++) 
    P7P_TSC(gm, k-1, p7P_LM) = logf(occ[k] / Z); /* note off-by-one: entry at Mk stored as [k-1][LM] */
  /* leave tsc[M,LM] as -inf, as we'd already initialized it */
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

/* The Config test simply makes sure a random profile passes
 * a Validate() check.
 */
static void
utest_Config(P7_HMM *hmm, P7_BG *bg)
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

/* gcc -g -Wall -Dp7MODELCONFIG_TESTDRIVE -I. -I../easel -L. -L../easel -o modelconfig_utest modelconfig.c -lhmmer -leasel -lm
 * ./modelconfig_utest
 */
#include "easel.h"

#include "p7_config.h"
#include "hmmer.h"


int
main(int argc, char **argv)
{  
  ESL_ALPHABET   *abc    = NULL;
  ESL_RANDOMNESS *r      = NULL;
  P7_HMM         *hmm    = NULL;
  P7_BG          *bg     = NULL;
  int             M      = 10000;
  
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create amino alphabet");
  if ((r   = esl_randomness_CreateFast(0))  == NULL)  esl_fatal("failed to create randomness");
  if (p7_hmm_Sample(r, M, abc, &hmm)        != eslOK) esl_fatal("failed to sample random HMM");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to created null model");

  utest_Config(hmm, bg);
  utest_occupancy(hmm);

  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7MODELCONFIG_TESTDRIVE*/



/*****************************************************************
 * 5. Statistics collection driver.
 *****************************************************************/
#ifdef p7MODELCONFIG_STATS
/* gcc -g -Wall -Dp7MODELCONFIG_STATS -I. -I../easel -L. -L../easel -o statprog modelconfig.c -lhmmer -leasel -lm
 * ./statprog
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",       0 },
  { "-i",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "sample by two-step ideal rule, not from profile", 0},
  { "-m",        eslARG_INFILE,  NULL, NULL, NULL,      NULL,      NULL, "-u,-M", "input HMM from file <f> instead of sampling",0 },
  { "-n",        eslARG_INT, "100000", NULL, "n>0",     NULL,      NULL,    NULL, "number of seqs to sample",                   0 },
  { "-s",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-m,-u", "make sampled HMM uniform transitions, as S/W", 0},
  { "-u",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-m,-s", "make sampled HMM ungapped",                  0 },
  { "-L",        eslARG_INT,    "400", NULL,"n>=0",     NULL,      NULL,    NULL, "set expected length from profile to <n>",    0 },
  { "-M",        eslARG_INT,     "50", NULL, "n>0",     NULL,      NULL,    "-m", "set sampled model length to <n>",            0 },
  { "-2",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "emulate HMMER2 configuration",               0 },
  { "--ips",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,    NULL, "output PostScript mx of i endpoints to <f>", 0 },
  { "--kps",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,    NULL, "output PostScript mx of k endpoints to <f>", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[] = "./statprog [options]";

static int ideal_local_endpoints(ESL_RANDOMNESS *r, P7_HMM *hmm, ESL_SQ *sq, P7_TRACE *tr, int Lbins,
				 int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2);
static int profile_local_endpoints(ESL_RANDOMNESS *r, P7_HMM *core, P7_PROFILE *gm, ESL_SQ *sq, P7_TRACE *tr, int Lbins,
				   int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2);

int
main(int argc, char **argv)
{
  ESL_ALPHABET    *abc     = NULL;     /* sequence alphabet                       */
  ESL_GETOPTS     *go      = NULL;     /* command line processing                 */
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
  P7_HMM          *hmm     = NULL;     /* sampled HMM to emit from                */
  P7_HMM          *core    = NULL;     /* safe copy of the HMM, before config     */
  P7_BG           *bg      = NULL;     /* null model                              */
  ESL_SQ          *sq      = NULL;     /* sampled sequence                        */
  P7_TRACE        *tr      = NULL;     /* sampled trace                           */
  P7_PROFILE      *gm      = NULL;     /* profile                                 */
  int              i,j;
  int              i1,i2;
  int              k1,k2;
  int              iseq;
  FILE            *fp      = NULL;
  double           expected;

  int              do_ilocal;
  char            *hmmfile = NULL;
  int              nseq;
  int              do_swlike;
  int              do_ungapped;
  int              L;
  int              M;
  int              do_h2;
  char            *ipsfile = NULL;
  char            *kpsfile = NULL;
  ESL_DMATRIX     *imx     = NULL;
  ESL_DMATRIX     *kmx     = NULL;
  ESL_DMATRIX     *iref    = NULL; /* reference matrix: expected i distribution under ideality */
  int              Lbins;
  int              status;
  char             errbuf[eslERRBUFSIZE];
  
  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    puts(usage);
    puts("\n  where options are:\n");
    esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=all docgroups; 2 = indentation; 80=textwidth*/
    return eslOK;
  }
  do_ilocal   = esl_opt_GetBoolean(go, "-i");
  hmmfile     = esl_opt_GetString (go, "-m");
  nseq        = esl_opt_GetInteger(go, "-n");
  do_swlike   = esl_opt_GetBoolean(go, "-s");
  do_ungapped = esl_opt_GetBoolean(go, "-u");
  L           = esl_opt_GetInteger(go, "-L");
  M           = esl_opt_GetInteger(go, "-M");
  do_h2       = esl_opt_GetBoolean(go, "-2");
  ipsfile     = esl_opt_GetString (go, "--ips");
  kpsfile     = esl_opt_GetString (go, "--kps");

  if (esl_opt_ArgNumber(go) != 0) {
    puts("Incorrect number of command line arguments.");
    printf("Usage: %s [options]\n", argv[0]);
    return eslFAIL;
  }

  r = esl_randomness_CreateFast(0);

  if (hmmfile != NULL)
    {	/* Read the HMM (and get alphabet from it) */
      P7_HMMFILE      *hfp     = NULL;

      status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
      if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
      else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
      else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);  
    
      if ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslOK) {
	if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
	else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
	else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
	else                             esl_fatal("Unexpected error in reading HMMs");
      }
      M = hmm->M;
      p7_hmmfile_Close(hfp);
    }
  else
    {			/* Or sample the HMM (create alphabet first) */
      abc = esl_alphabet_Create(eslAMINO);    
      if      (do_ungapped) p7_hmm_SampleUngapped(r, M, abc, &hmm);
      else if (do_swlike)   p7_hmm_SampleUniform (r, M, abc, 0.05, 0.5, 0.05, 0.2, &hmm); /* tmi, tii, tmd, tdd */
      else                  p7_hmm_Sample        (r, M, abc, &hmm);
    }

  Lbins = M;
  imx  = esl_dmatrix_Create(Lbins, Lbins);
  iref = esl_dmatrix_Create(Lbins, Lbins);
  kmx  = esl_dmatrix_Create(M, M);
  esl_dmatrix_SetZero(imx);
  esl_dmatrix_SetZero(iref);
  esl_dmatrix_SetZero(kmx);
  tr    = p7_trace_Create();
  sq    = esl_sq_CreateDigital(abc);
  bg    = p7_bg_Create(abc);
  core  = p7_hmm_Clone(hmm);

  if (do_h2) {
    gm = p7_profile_Create(hmm->M, abc);
    p7_H2_ProfileConfig(hmm, bg, gm, p7_UNILOCAL);
  } else {
    gm = p7_profile_Create(hmm->M, abc);
    p7_profile_ConfigUnilocal(gm, hmm, bg, L);
    if (p7_hmm_Validate    (hmm, NULL, 0.0001) != eslOK) esl_fatal("whoops, HMM is bad!");
    if (p7_profile_Validate(gm,  NULL, 0.0001) != eslOK) esl_fatal("whoops, profile is bad!");
  }

  /* Sample endpoints.
   * Also sample an ideal reference distribution for i endpoints.  i
   * endpoints are prone to discretization artifacts, when emitted
   * sequences have varying lengths. Taking log odds w.r.t. an ideal
   * reference that is subject to the same discretization artifacts 
   * cancels out the effect.
   */
  for (iseq = 0; iseq < nseq; iseq++)
    {				
      if (do_ilocal) ideal_local_endpoints  (r, core,     sq, tr, Lbins, &i1, &i2, &k1, &k2);
      else           profile_local_endpoints(r, core, gm, sq, tr, Lbins, &i1, &i2, &k1, &k2);

      imx->mx[i1-1][i2-1] += 1.;
      kmx->mx[k1-1][k2-1] += 1.; 

      /* reference distribution for i */
      ideal_local_endpoints  (r, core, sq, tr, Lbins, &i1, &i2, &k1, &k2);
      iref->mx[i1-1][i2-1] += 1.;
    }


  /* Adjust both mx's to log_2(obs/exp) ratio */
  printf("Before normalization/log-odds:\n");
  printf("   i matrix values range from %f to %f\n", dmx_upper_min(imx), dmx_upper_max(imx));
  printf("   k matrix values range from %f to %f\n", dmx_upper_min(kmx), dmx_upper_max(kmx));
  printf("iref matrix values range from %f to %f\n", dmx_upper_min(iref), dmx_upper_max(iref));

  expected = (double) nseq * 2. / (double) (M*(M+1));
  for (i = 0; i < kmx->m; i++)
    for (j = i; j < kmx->n; j++)
      kmx->mx[i][j] = log(kmx->mx[i][j] / expected) / log(2.0);

  for (i = 0; i < imx->m; i++)
    for (j = i; j < imx->m; j++)
      if (iref->mx[i][j] == 0. && imx->mx[i][j] == 0.) 
	imx->mx[i][j] = 0.;
      else if (iref->mx[i][j] == 0.)
	imx->mx[i][j] = eslINFINITY;
      else if (imx->mx[i][j] == 0.)
	imx->mx[i][j] = -eslINFINITY;
      else
	imx->mx[i][j] = log(imx->mx[i][j] / iref->mx[i][j]) / log(2.0);
  
  /* Print ps files */
  if (kpsfile != NULL) {
    if ((fp = fopen(kpsfile, "w")) == NULL) esl_fatal("Failed to open output postscript file %s", kpsfile);
    dmx_Visualize(fp, kmx, -4., 5.);
    fclose(fp);
  }
  if (ipsfile != NULL) {
    if ((fp = fopen(ipsfile, "w")) == NULL) esl_fatal("Failed to open output postscript file %s", ipsfile);
    dmx_Visualize(fp, imx, -4., 5.); 
    /* dmx_Visualize(fp, imx, dmx_upper_min(imx), dmx_upper_max(imx)); */
    fclose(fp);
  }

  printf("After normalization/log-odds:\n");
  printf("i matrix values range from %f to %f\n", dmx_upper_min(imx), dmx_upper_max(imx));
  printf("k matrix values range from %f to %f\n", dmx_upper_min(kmx), dmx_upper_max(kmx));

  
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(core);
  p7_hmm_Destroy(hmm);
  p7_trace_Destroy(tr);
  esl_sq_Destroy(sq);
  esl_dmatrix_Destroy(imx);
  esl_dmatrix_Destroy(kmx);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}

/* ideal_local_endpoints()
 *
 * Purpose:  Implementation of the "two-step" fragment sampling
 *           algorithm, sampling a uniform local fragment w.r.t.
 *           sequence coords, by first sampling a complete
 *           sequence of length L from <hmm>; then choosing
 *           a random fragment <i1..i2> uniformly from all
 *           possible $\frac{L(L+1)/2}$ fragments;  then finding
 *           local alignment coordinates wrt model and sequence,
 *           using convention that local alignment starts/stops
 *           with match states. (Thus, if the initially selected
 *           i1 or i2 were generated by insert states, bounds
 *           are moved to reach first/last match state.)
 *           
 *           The caller also provides an allocated sequence <sq> and
 *           traceback <tr>, as storage to be provided to
 *           <p7_CoreEmit()>. They contain the generated global
 *           sequence and trace upon return (not a local trace, note).
 *           
 *           i endpoints are normalized/discretized to 1..<Lbins>, so
 *           we can collate i statistics from sampled sequences of
 *           varying L. Note this causes discretization artifacts,
 *           leading to underrepresentation of j=M and
 *           overrepresentation of i=1.
 *           
 *           This routine is only intended for collecting endpoint
 *           statistics (i1,i2,k1,k2); it does not generate a local
 *           alignment trace. (xref milestone 2, STL11/115).
 *           
 * Returns:  <eslOK> on success; returns normalized/binned sequence
 *           coords in <*ret_i1> and <*ret_i2> in range <1..Lbins> and
 *           the model entry/exit coords in <*ret_k1> and <*ret_k2> in
 *           range <1..M>. By internal def'n of local alignment endpoints,
 *           M_k1 emits residue x_i1, M_k2 emits residue x_i2.
 *           
 * Xref:     STL11/142-143 
 */
static int
ideal_local_endpoints(ESL_RANDOMNESS *r, P7_HMM *hmm, ESL_SQ *sq, P7_TRACE *tr, int Lbins,
		      int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2)
{
  int status;
  int tpos;
  int i1, i2, k1,k2, t1,t2;
  int all_insert;
  int failsafe = 0;		/* a failsafe timer for rejection sampling */

  do {
    if (failsafe++ == 1000) ESL_XEXCEPTION(eslENOHALT, "failed to obtain local alignment that wasn't all inserts");

    if ((status = p7_CoreEmit(r, hmm, sq, tr)) != eslOK) goto ERROR;

    /* a simple way to sample uniformly from upper triangle is by rejection 
     * this do/while cannot infinite loop, doesn't need failsafe 
     */
    do {
      i1 = 1 + esl_rnd_Roll(r, sq->n);
      i2 = 1 + esl_rnd_Roll(r, sq->n);
    } while (i1 > i2);

    /* Get initial k1,k2 coords: this step must work in a core model, 
     * i1/i2 were generated by an M or I. Also record t1,t2 endpoints
     * on core's trace.
     */
    for (tpos = 0; tpos < tr->N; tpos++)
      if (tr->i[tpos] == i1) { t1 = tpos; k1 = tr->k[tpos]; break; }
    for (tpos = tr->N-1; tpos >= 0; tpos--)
      if (tr->i[tpos] == i2) { t2 = tpos; k2 = tr->k[tpos]; break; }

    /* Enforce the definition of local alignment endpoints being
     * match-delimited - roll up any leading/trailing I states. 
     * Watch out for pathological case of a local fragment that
     * includes no M state at all.
     */
    all_insert = FALSE;
    for (; t1 <= t2; t1++) if (tr->st[t1] == p7T_M) break;
    for (; t2 >= t1; t2--) if (tr->st[t2] == p7T_M) break;
    if (t2 < t1) all_insert = TRUE; /* sufficient to check both. */
    i1 = tr->i[t1];  i2 = tr->i[t2];
    k1 = tr->k[t1];  k2 = tr->k[t2];
  } while (all_insert);

  /* Normalize sequence coords.
   * They're 1..L now; make them 1..Lbins
   */
  *ret_i1 = ((i1-1) * Lbins / sq->n) + 1;
  *ret_i2 = ((i2-1) * Lbins / sq->n) + 1;
  *ret_k1 = k1;
  *ret_k2 = k2;
  return eslOK;

 ERROR:
  *ret_i1 = 0.;
  *ret_i2 = 0.;
  *ret_k1 = 0;
  *ret_k2 = 0;
  return status;
}

/* profile_local_endpoints()
 *
 * Purpose:   Wrapper around <p7_ProfileEmit()>, sampling a local
 *            alignment fragment from the profile's probabilistic model
 *            (which may be the implicit model of HMMER3, or the
 *            Plan7 model of HMMER2), and reporting coordinates
 *            of the fragment w.r.t. both model and sequence.
 *            
 *            To simplify the implementation, the profile must be in
 *            <p7_UNILOCAL> mode, not <p7_LOCAL> mode, so we know we
 *            only have to deal with a single hit per sampled
 *            sequence. 
 *            
 *            We want <i1..i2> to be relative to the sequence coords
 *            of a complete (global) sampled sequence that we could
 *            have sampled this local alignment from; but the <i1..i2>
 *            we initially get are relative to our profile-sampled
 *            trace, so they are offset both by N-generated residues
 *            that occur in the profile and by residues that the
 *            profile's local entry skipped. To translate from
 *            profile/sequence coords to core model/sequence coords,
 *            we use rejection sampling: sample traces from the core
 *            model until we find one that uses the same statetypes
 *            at *initial* entry/exit points <k1>,<k2>, then use
 *            that sample's sequence to determine offsets and correct
 *            <i1..i2> reference frame.
 *            
 *            Local alignment endpoints are defined to be
 *            match-delimited. However, an H3 model allows exit on
 *            either a D or M state. Thus, the initially sampled end
 *            point k2 may need to be rolled back to last M state, to
 *            satisfy local alignment endpoint definition. Entries are
 *            not a problem; both H2 and H3 profiles can only enter on
 *            a M state. (This rollback has to occur after we've
 *            matched a core trace to the profile trace to determine
 *            i offsets.)
 *            
 *            Then, sampling from both the core model and the profile
 *            in the same routine introduces a complication:
 *            conceivably, profile configuration alters the transition
 *            probabilities in the core model (by adding <M->E>
 *            transitions and renormalizing the M transition
 *            distributions, for example; H2 configuration does this,
 *            though H3 does not). So you can't <CoreSample()> the
 *            <gm->hmm> safely. To avoid such things, the caller
 *            provides a clean copy of the core model in <core>.
 *            
 *           i endpoints are normalized/discretized to 1..<Lbins>, so
 *           we can collate i statistics from sampled sequences of
 *           varying L. Note this causes discretization artifacts,
 *           leading to underrepresentation of j=M and
 *           overrepresentation of i=1.
 *           
 * Returns:  <eslOK> on success; returns normalized sequence coords in
 *           <*ret_i1> and <*ret_i2>, and the model entry/exit coords
 *           in <*ret_k1> and <*ret_k2>. 
 *           
 * Xref:     STL11/142-143 
 */
static int
profile_local_endpoints(ESL_RANDOMNESS *r, P7_HMM *core, P7_PROFILE *gm, ESL_SQ *sq, P7_TRACE *tr, int Lbins,
			int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2)
{
  int status;
  int i1,i2;
  int k1,k2;
  int t1,t2;			/* entry/exit positions in local trace, tr */
  int tg1, tg2;			/* entry/exit positions in global trace, tr2 */
  int tpos;
  int nterm, cterm;		/* offsets at N, C terminus. */
  int L;			/* inferred length from 3-part patching */
  ESL_SQ *sq2   = NULL;
  P7_TRACE *tr2 = NULL;
  int failsafe  = 0;
  
  if (gm->mode != p7_UNILOCAL) ESL_XEXCEPTION(eslEINVAL, "profile must be unilocal");
  if ((sq2 = esl_sq_CreateDigital(gm->abc))  == NULL)   { status = eslEMEM; goto ERROR; }
  if ((tr  = p7_trace_Create())              == NULL)   { status = eslEMEM; goto ERROR; }

  /* sample local alignment from the implicit model */
  if (gm->h2_mode) {
    if ((status = p7_H2_ProfileEmit(r, gm, sq, tr)) != eslOK) goto ERROR;
  } else {
    if ((status = p7_ProfileEmit(r, gm, sq, tr)) != eslOK) goto ERROR;
  }
    
  /* Get initial trace coords */
  for (tpos = 0;       tpos < tr->N; tpos++)  if (tr->st[tpos] == p7T_B) { t1 = tpos+1; break; }
  for (tpos = tr->N-1; tpos >= 0;    tpos--)  if (tr->st[tpos] == p7T_E) { t2 = tpos-1; break; }
  
  /* Match a core trace to this local trace by rejection sampling;
   * this is to let us calculate sequence offsets; see comments above in preamble
   */
  do {
    if (failsafe++ == 100000) ESL_XEXCEPTION(eslENOHALT, "failed to match core,local traces in %d tries\n", failsafe);

    if ((status = p7_CoreEmit(r, core, sq2, tr2)) != eslOK) goto ERROR;
    for (tpos = 0; tpos < tr2->N; tpos++)
      if (tr2->k[tpos] == tr->k[t1]) { tg1 = tpos; break; }
    for (tpos = tr2->N-1; tpos >= 0; tpos--)
      if (tr2->k[tpos] == tr->k[t2]) { tg2 = tpos; break; }
  }  while (tr2->st[tg1] != tr->st[t1] && tr2->st[tg2] != tr->st[t2]);

  /* tg1..tg2 in core trace is now matched to t1..t2 in the profile trace.
   * Calculate # of residues preceding tg1 and following tg2 in the core trace.
   * A core trace can only generate residues from M or I states.
   */
  for (nterm = 0, tpos = 0; tpos < tg1; tpos++) 
    if (tr2->st[tpos] == p7T_M || tr2->st[tpos] == p7T_I) nterm++;
  for (cterm = 0, tpos = tr2->N-1; tpos > tg2; tpos--)
    if (tr2->st[tpos] == p7T_M || tr2->st[tpos] == p7T_I) cterm++;

  /* rectify the t2 endpoint, rolling back any trailing D path 
   */
  for (; t2 >= 0; t2--) if (tr->st[t2] == p7T_M) break;
  if (t2 < t1) ESL_XEXCEPTION(eslEINCONCEIVABLE, "this only happens on an all-D path through profile");  
  
  /* determine initial endpoint coords from t1 and t2 */
  i1 = tr->i[t1];  i2 = tr->i[t2];
  k1 = tr->k[t1];  k2 = tr->k[t2];

  /* offset the i coords. */
  L  = (i2-i1+1) + nterm + cterm;
  i2 = (i2-i1+1) + nterm;
  i1 = nterm+1;

  /* normalize the i coords into range 1..Lbins, instead of 1..L */
  i1 = ((i1-1) * Lbins / L) + 1;
  i2 = ((i2-1) * Lbins / L) + 1;

  *ret_i1 = i1;
  *ret_i2 = i2;
  *ret_k1 = k1;
  *ret_k2 = k2;
  p7_trace_Destroy(tr2);
  esl_sq_Destroy(sq2);
  return eslOK;

 ERROR:
  if (sq2 != NULL)  esl_sq_Destroy(sq2);
  if (tr2 != NULL)  p7_trace_Destroy(tr2);
  *ret_i1 = 0.;
  *ret_i2 = 0.;
  *ret_k1 = 0;
  *ret_k2 = 0;
  return status;
}
#endif /*p7MODELCONFIG_STATS*/

/*****************************************************************
 * 6. References.
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
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
