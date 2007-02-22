/* Model configuration: 
 * Converting a core model to a fully configured Plan7 search profile.
 * 
 * Contents:
 *     1. Routines in the exposed API.
 *     2. Private functions.
 *     3. Routines emulating HMMER2, for testing.
 *     4. Unit tests.
 *     5. Test driver.
 *     6. Statistics collection driver.
 *     7. Copyright and license
 * 
 * Revised May 2005: xref STL9/77-81.       (Uniform fragment distribution)
 * Again, Sept 2005: xref STL10/24-26.      (Inherent target length dependency)
 * Again, Jan 2007:  xref STL11/125,136-137 (HMMER3)
 *
 * SRE, Mon May  2 10:55:16 2005 [St. Louis]
 * SRE, Fri Jan 12 08:06:33 2007 [Janelia] [Kate Bush, Aerial]
 * SVN $Id$
 */
#include "p7_config.h"

#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

static int calculate_occupancy(P7_HMM *hmm, float *occ);
static int logoddsify(P7_HMM *hmm, P7_PROFILE *gm);





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
 *            
 * Returns:   <eslOK> on success; the profile <gm> has its scores filled
 *            in, its implicit probabilistic model probabilities filled in,
 *            and it contains copied reference pointers for the HMM's
 *            alphabet, null model, and the HMM itself.
 *            
 * Throws:    <eslECONTRACT> if the <hmm> doesn't have a null model assigned
 *            to it.
 */
int
p7_ProfileConfig(P7_HMM *hmm, P7_PROFILE *gm, int mode)
{
  int   k;			/* counter over states      */
  float *occ = NULL;
  float Z;
  int   status;

  /* Contract checks: HMM must have null model */
  if (hmm->bg == NULL)  
    ESL_XEXCEPTION(eslECONTRACT, "HMM needs to have a null model here");
  if (mode != p7_LOCAL && mode != p7_UNILOCAL) 
    ESL_XEXCEPTION(eslECONTRACT, "I'm not ready for any mode but local modes");

  /* Copy some pointer references and other info across from HMM
   */
  gm->M    = hmm->M;
  gm->abc  = hmm->abc;
  gm->hmm  = hmm;
  gm->bg   = hmm->bg;
  gm->mode = mode;

  /* E state loop/move probabilities: nonzero for MOVE allows loops/multihits
   * N,C,J transitions are set later by length config 
   */
  if (mode == p7_LOCAL) {
    gm->xt[p7_XTE][p7_MOVE] = 0.5;  
    gm->xt[p7_XTE][p7_LOOP] = 0.5;  
  } else {
    gm->xt[p7_XTE][p7_MOVE] = 1.0;  
    gm->xt[p7_XTE][p7_LOOP] = 0.0;  
  }
    
  /* Begin probabilities:   occ[k] /( \sum_i occ[i] * (M-k+1))
   * (Reduces to 2/(M(M+1)) for occupancies of 1.0)
   * These are only probabilistic w.r.t. the implicit model.
   */
  ESL_ALLOC(occ, sizeof(float) * (hmm->M+1));
  if ((status = calculate_occupancy(hmm, occ)) != eslOK) goto ERROR;
  for (Z = 0., k = 1; k <= hmm->M; k++) 
    Z += occ[k] * (float) (hmm->M-k+1);
  for (gm->begin[0] = 0., k=1; k<=hmm->M; k++)
    gm->begin[k] = occ[k] / Z;
  free(occ);

  /* Exit probabilities: 1.0
   * 
   * Only probabilistic w.r.t. the implicit model.  HMMER3 will use
   * implicit Dk->E and Mk->E exit probabilities of 1.0 but we use
   * end[k] for the Mk->E exits, at least temporarily, for regression
   * testing against HMMER2 configuration strategies.
   */
  gm->end[0] = 0.;
  esl_vec_FSet(gm->end+1, hmm->M, 1.0);

  /* Set probabilities and scores for N,C,J; then the rest of the model
   */
  if ((status = p7_ReconfigLength(gm, 350.))        != eslOK) goto ERROR;
  if ((status = logoddsify(hmm, gm))                != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  if (occ != NULL) free(occ);
  return status;
}


/* Function:  p7_ReconfigLength()
 * Incept:    SRE, Sun Sep 25 12:38:55 2005 [St. Louis]
 *
 * Purpose:   Given a model already configured for scoring, in some
 *            particular algorithm mode; reset the expected length
 *            distribution of both the HMM and the null model to a
 *            mean of <L>.
 *            
 *            Do this as quickly as possible, because the caller needs
 *            to dynamically reconfigure the model for the length of
 *            each target sequence in a database search.  
 *
 * Returns:   <eslOK> on success.
 *            p1, xt[NCJ] probabilities, and xsc[NCJ] scores are set 
 *            here. These control the target length dependence of the
 *            model. 
 *            
 * Throws:    <eslECONTRACT> if the <gm> does not contain a null model.           
 */
int
p7_ReconfigLength(P7_PROFILE *gm, int L)
{
  float ploop, pmove;
  float nj;
  int   status;

  /* Contract checks: profile must have null model */
  if (gm->bg == NULL) ESL_XEXCEPTION(eslECONTRACT, "profile needs to have a null model here");

  /* Configure p1 in the null model to an expected length of L */
  gm->bg->p1 = (float) L / (float) (L+1);
  
  /* Figure out the expected number of uses of the J state.
   */
  if (gm->xt[p7_XTE][p7_LOOP] == 0.) nj = 0.; /* make sure. */
  else                               nj = gm->xt[p7_XTE][p7_LOOP] / (1. - gm->xt[p7_XTE][p7_LOOP]);

  /* Configure N,J,C transitions so they bear L/(2+nj) of the total
   * unannotated sequence length L. 
   */
  pmove = (2. + nj) / ((float) L + 2. + nj); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1. - pmove;
  gm->xt[p7_XTN][p7_MOVE] = pmove;
  gm->xt[p7_XTN][p7_LOOP] = ploop;
  gm->xt[p7_XTC][p7_MOVE] = pmove;
  gm->xt[p7_XTC][p7_LOOP] = ploop;
  gm->xt[p7_XTJ][p7_MOVE] = pmove;	/* note, J unused if [XTE][LOOP] = 0. */
  gm->xt[p7_XTJ][p7_LOOP] = ploop;

  /* Set the N,J,C scores. (integer scaled lod scores) */
  gm->xsc[p7_XTN][p7_LOOP] = p7_Prob2Score(gm->xt[p7_XTN][p7_LOOP], gm->bg->p1);
  gm->xsc[p7_XTN][p7_MOVE] = p7_Prob2Score(gm->xt[p7_XTN][p7_MOVE], 1.0);
  gm->xsc[p7_XTC][p7_LOOP] = p7_Prob2Score(gm->xt[p7_XTC][p7_LOOP], gm->bg->p1);
  gm->xsc[p7_XTC][p7_MOVE] = p7_Prob2Score(gm->xt[p7_XTC][p7_MOVE], 1.0 - gm->bg->p1);
  gm->xsc[p7_XTJ][p7_LOOP] = p7_Prob2Score(gm->xt[p7_XTJ][p7_LOOP], gm->bg->p1);
  gm->xsc[p7_XTJ][p7_MOVE] = p7_Prob2Score(gm->xt[p7_XTJ][p7_MOVE], 1.0);

  /* Detect the "special" (actually common) case of the LOOP scores
   * having too small of a magnitude to keep track of properly.  We
   * will have trouble with the [NCJ][LOOP] scores for large L,
   * because these scores will be ~ 1/L, and we can only hold scores
   * down to 0.001 bits, if INTSCALE is at its default 1000. As a
   * workaround, we catch the case where the absolute value of
   * these scores would be <10. In this case, we set a PLAN7_LCORRECT
   * flag in the hmm, and store the difference between the expected
   * score per residue in floating pt versus integer lod scores as
   * hmm->lscore. We can then post-hoc correct an alignment score by
   * L' * hmm->lscore.  
   * 
   * This code assumes that all three scores (N,C,J)[LOOP] are equal,
   * which is how we set them above. gm->lscore becomes a correction
   * per unannotated residue that we add, when the do_lcorrect flag is
   * set.  Note that the correction is "soft": we try to use what
   * we've got as an integer score, and only add an expected
   * difference.
   * (xref STL10/26)
   */
  if (abs(gm->xsc[p7_XTN][p7_LOOP]) < 10) {
    gm->do_lcorrect = TRUE;
    /* the real cost per residue, as a float: */
    gm->lscore   = log(gm->xt[p7_XTN][p7_LOOP] / gm->bg->p1);
    /* minus what we're going to imprecisely calculate, in integer scores: */
    gm->lscore -= (float) gm->xsc[p7_XTN][p7_LOOP] / (float) p7_INTSCALE;
  } else {
    gm->do_lcorrect = FALSE;
    gm->lscore = 0.;
  }

  return eslOK;

 ERROR:
  return status;
}

/*****************************************************************
 * 2. Private functions
 *****************************************************************/

/* calculate_occupancy()
 * Incept:    SRE, Mon Jan 22 08:10:05 2007 [Janelia]
 *
 * Purpose:   Calculate a vector <occ[1..M]> containing probability
 *            that each match state is used in a sampled path through
 *            the model. Caller provides allocated space (<M+1> floats)
 *            for <occ>.
 *
 * Returns:   <eslOK> on success.
 */
static int
calculate_occupancy(P7_HMM *hmm, float *occ)
{
  int k;

  occ[0] = 0.;			/* no M_0 state */
  occ[1] = hmm->t[0][p7_TMM];	/* initialize w/ B->M_1 */
  for (k = 2; k <= hmm->M; k++)
    occ[k] = occ[k-1] * (hmm->t[k-1][p7_TMM] + hmm->t[k-1][p7_TMI]) +
      (1.0-occ[k-1]) * hmm->t[k-1][p7_TDM];
  return eslOK;
}

/* logoddsify()
 * Incept:    SRE, Mon Jan 22 08:45:57 2007 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
static int
logoddsify(P7_HMM *hmm, P7_PROFILE *gm)
{
  int k, x;
  int sc[p7_MAXCODE];

  /* p7_profile_Create() has already initialized unused memory */

  /* Match and insert transitions, 1..M-1.  */
  for (k = 1; k < gm->M; k++) {
    gm->tsc[p7_TMM][k] = p7_Prob2Score(hmm->t[k][p7_TMM], gm->bg->p1);
    gm->tsc[p7_TMI][k] = p7_Prob2Score(hmm->t[k][p7_TMI], gm->bg->p1);
    gm->tsc[p7_TMD][k] = p7_Prob2Score(hmm->t[k][p7_TMD], 1.0);
    gm->tsc[p7_TIM][k] = p7_Prob2Score(hmm->t[k][p7_TIM], gm->bg->p1);
    gm->tsc[p7_TII][k] = p7_Prob2Score(hmm->t[k][p7_TII], gm->bg->p1);
  }
  /* Delete transitions, 2..M-1. */
  for (k = 2; k < gm->M; k++) {
    gm->tsc[p7_TDM][k] = p7_Prob2Score(hmm->t[k][p7_TDM], gm->bg->p1);
    gm->tsc[p7_TDD][k] = p7_Prob2Score(hmm->t[k][p7_TDD], 1.0);
  }

  /* Match emissions (including degeneracies). 
   * Temp sc[x] vector because of the rearranged msc[x][k] array */
  sc[gm->abc->K] = p7_IMPOSSIBLE;
  for (k = 1; k < gm->M; k++)  {
    for (x = 0; x < gm->abc->K; x++)
      sc[x] = p7_Prob2Score(hmm->mat[k][x], gm->bg->f[x]); /* base */
    esl_abc_IExpectScVec(gm->abc, sc, gm->bg->f);             /* degens */
    for (x = 0; x < gm->abc->K; x++)
      gm->msc[x][k] = sc[x];	
  }
  
  /* Then the same for insert emissions */
  for (k = 1; k < gm->M; k++) {
    for (x = 0; x < gm->abc->K; x++)
      sc[x] = p7_Prob2Score(hmm->ins[k][x], gm->bg->f[x]); /* base */
    esl_abc_IExpectScVec(gm->abc, sc, gm->bg->f);             /* degens */
    for (x = 0; x < gm->abc->K; x++)
      gm->isc[x][k] = sc[x];	
  }
    
  /* B->Mk begin transitions 1..M */
  for (k = 1; k <= gm->M; k++)
    gm->bsc[k] = p7_Prob2Score(gm->begin[k], gm->bg->p1);
  
  /* Mk->E end transitions 1..M */
  for (k = 1; k <= gm->M; k++)
    gm->esc[k] = 0;		/* by construction */

  /* E state transitions (loop, move); N,C,J are done by p7_ReconfigLength() */
  gm->xsc[p7_XTE][p7_LOOP] = p7_Prob2Score(gm->xt[p7_XTE][p7_LOOP], 1.0);
  gm->xsc[p7_XTE][p7_MOVE] = p7_Prob2Score(gm->xt[p7_XTE][p7_MOVE], 1.0);

  return eslOK;
}

/*****************************************************************
 * 3. Routines emulating HMMER2, for testing.
 *****************************************************************/

/* Function:  p7_H2_ProfileConfig()
 * Incept:    SRE, Thu Jan 25 17:12:31 2007 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      HMMER2.3.2 plan7.c::Plan7FSConfig()
 */
int
p7_H2_ProfileConfig(P7_HMM *hmm, P7_PROFILE *gm, int mode)
{
  int   k;			/* counter over states      */
  float basep;
  float pentry = 0.5;
  float pexit  = 1.0;
  int   status;
  float d;

  /* Contract checks: HMM must have null model */
  if (hmm->bg == NULL)  ESL_XEXCEPTION(eslECONTRACT, "HMM needs to have a null model here");
  if (mode != p7_LOCAL && mode != p7_UNILOCAL) 
    ESL_XEXCEPTION(eslECONTRACT, "I'm not ready for any mode but local mode");

  /* Copy some pointer references and other info across from HMM */
  gm->M       = hmm->M;
  gm->abc     = hmm->abc;
  gm->hmm     = hmm;
  gm->bg      = hmm->bg;
  gm->mode    = mode;
  gm->h2_mode = TRUE;

  /* Configure special states */
  gm->xt[p7_XTN][p7_MOVE] = 1.0 - hmm->bg->p1;    /* allow N-terminal tail     */
  gm->xt[p7_XTN][p7_LOOP] = hmm->bg->p1;
  if (mode == p7_LOCAL) {
    gm->xt[p7_XTE][p7_MOVE] = 0.5; /* multihit */
    gm->xt[p7_XTE][p7_LOOP] = 0.5;  
  } else {
    gm->xt[p7_XTE][p7_MOVE] = 1.0; /* singlehit */
    gm->xt[p7_XTE][p7_LOOP] = 0.0;  
  }
  gm->xt[p7_XTC][p7_MOVE] = 1.0 - hmm->bg->p1;    /* allow C-terminal tail     */
  gm->xt[p7_XTC][p7_LOOP] = hmm->bg->p1;
  gm->xt[p7_XTJ][p7_MOVE] = 1.0 - hmm->bg->p1;    /* allow J junction between domains */
  gm->xt[p7_XTJ][p7_LOOP] = hmm->bg->p1;

  /* Configure entry. HMMER2 used (almost) uniform entry. */  
  /* gm->begin[1] = (1. - pentry) * (1. - hmm->t[0][p7_TMD]); */
  /* gm->begin[1] = (1. - pentry);
     esl_vec_FSet(gm->begin+2, hmm->M-1, (pentry * (1.-hmm->t[0][p7_TMD])) / (float) (hmm->M-1));
  */
  esl_vec_FSet(gm->begin+1, hmm->M, 1.0 / hmm->M);

  /* Configure exit, uniform exit given entry: which is 1/(M-k+1)  */
  gm->end[hmm->M] = 1.0;
  basep = pexit / (float) (hmm->M-1);
  for (k = 1; k < hmm->M; k++)
    gm->end[k] = basep / (1. - basep * (float) (k-1));
  /* renormalize exits (this is Plan7RenormalizeExits from H2, inlined */
  for  (k = 1; k < hmm->M; k++) {
    d = esl_vec_FSum(hmm->t[k], 3);
    esl_vec_FScale(hmm->t[k], 3, 1./(d + d*gm->end[k]));
  }

  /* logoddsify the stuff that ReconfigLength would've done  */
  gm->xsc[p7_XTN][p7_LOOP] = p7_Prob2Score(gm->xt[p7_XTN][p7_LOOP], gm->bg->p1);
  gm->xsc[p7_XTN][p7_MOVE] = p7_Prob2Score(gm->xt[p7_XTN][p7_MOVE], 1.0);
  gm->xsc[p7_XTC][p7_LOOP] = p7_Prob2Score(gm->xt[p7_XTC][p7_LOOP], gm->bg->p1);
  gm->xsc[p7_XTC][p7_MOVE] = p7_Prob2Score(gm->xt[p7_XTC][p7_MOVE], 1.0 - gm->bg->p1);
  gm->xsc[p7_XTJ][p7_LOOP] = p7_Prob2Score(gm->xt[p7_XTJ][p7_LOOP], gm->bg->p1);
  gm->xsc[p7_XTJ][p7_MOVE] = p7_Prob2Score(gm->xt[p7_XTJ][p7_MOVE], 1.0);

  /* logoddsify everything else (without HMMER2's wing retraction) */
  logoddsify(hmm, gm);
  return eslOK;

 ERROR:
  return status;
}

/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7MODELCONFIG_TESTDRIVE

/* The Config test simply makes sure a random profile passes
 * a Validate() check.
 */
static void
utest_Config(P7_HMM *hmm)
{
  char       *msg = "modelconfig.c::p7_ProfileConfig() unit test failed";
  P7_PROFILE *gm  = NULL;

  if ((gm = p7_profile_Create(hmm->M, hmm->abc)) == NULL)   esl_fatal(msg);
  if (p7_ProfileConfig(hmm, gm, p7_LOCAL)        != eslOK)  esl_fatal(msg);
  if (p7_profile_Validate(gm, 0.0001)            != eslOK)  esl_fatal(msg);
  return;
}

/* The occupancy test is based on the principle that
 * the stationary match occupancy probability in a random HMM 
 * converges to 0.6, for long enough M (STL11/138)
 */
static void
utest_occupancy(P7_HMM *hmm)
{
  char  *msg = "modelconfig.c::calculate_occupancy() unit test failed";
  float *occ;
  float  x;

  occ = malloc(sizeof(float) * (hmm->M+1));
  calculate_occupancy(hmm, occ);
  x = esl_vec_FSum(occ+1, hmm->M) / (float) hmm->M;
  if (esl_FCompare(x, 0.6, 0.1) != eslOK)           esl_fatal(msg);
  free(occ);
  return;
}



#endif /*p7MODELCONFIG_TESTDRIVE*/



/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7MODELCONFIG_TESTDRIVE

/* gcc -g -Wall -Dp7MODELCONFIG_TESTDRIVE -I. -I../easel -L. -L../easel -o testprog modelconfig.c -lhmmer -leasel -lm
 * ./testprog
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
  int             M      = 10000;
  
  if ((abc = esl_alphabet_Create(eslAMINO))     == NULL)  esl_fatal("failed to create amino alphabet");
  if ((r   = esl_randomness_CreateTimeseeded()) == NULL)  esl_fatal("failed to create randomness");
  if (p7_hmm_Sample(r, M, abc, &hmm)            != eslOK) esl_fatal("failed to sample random HMM");
  if ((hmm->bg = p7_bg_Create(abc))             == NULL)  esl_fatal("failed to created null model");

  utest_Config(hmm);
  utest_occupancy(hmm);


  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7MODELCONFIG_TESTDRIVE*/



/*****************************************************************
 * 6. Statistics collection driver.
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
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",       0 },
  { "-i",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "sample by two-step ideal rule, not from profile", 0},
  { "-m",        eslARG_INFILE,  NULL, NULL, NULL,      NULL,      NULL, "-u,-M", "input HMM from file <f> instead of sampling",0 },
  { "-n",        eslARG_INT, "100000", NULL, "n>0",     NULL,      NULL,    NULL, "number of seqs to sample",                   0 },
  { "-u",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    "-m", "make sampled HMM ungapped",                  0 },
  { "-L",        eslARG_INT,    "400", NULL,"n>=0",     NULL,      NULL,    NULL, "set expected length from profile to <n>",    0 },
  { "-M",        eslARG_INT,     "50", NULL, "n>0",     NULL,      NULL,    "-m", "set sampled model length to <n>",            0 },
  { "-2",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "emulate HMMER2 configuration",               0 },
  { "--ips",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,    NULL, "output PostScript mx of i endpoints to <f>", 0 },
  { "--kps",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,    NULL, "output PostScript mx of k endpoints to <f>", 0 },
  { "--ibins",   eslARG_INT,    "100", NULL, "n>0",     NULL,      NULL,    NULL, "set # bins in i endpoint plot to <n>",       0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[] = "./statprog [options]";

static int ideal_local_endpoints(ESL_RANDOMNESS *r, P7_HMM *hmm, ESL_SQ *sq, P7_TRACE *tr, int ibins,
			 int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2);
static int emitted_local_endpoints(ESL_RANDOMNESS *r, P7_HMM *core, P7_PROFILE *gm, ESL_SQ *sq, P7_TRACE *tr, int ibins,
			 int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2);

int
main(int argc, char **argv)
{
  int              status;
  ESL_ALPHABET    *abc     = NULL;     /* sequence alphabet                       */
  ESL_GETOPTS     *go      = NULL;     /* command line processing                 */
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
  P7_HMM          *hmm     = NULL;     /* sampled HMM to emit from                */
  P7_HMM          *core    = NULL;     /* safe copy of the HMM, before config     */
  ESL_SQ          *sq      = NULL;     /* sampled sequence                        */
  P7_TRACE        *tr      = NULL;     /* sampled trace                           */
  int              i,j;
  int              i1,i2;
  int              k1,k2;
  int              iseq;
  FILE            *fp      = NULL;

  int              do_ilocal;
  char            *hmmfile = NULL;
  int              nseq;
  int              do_ungapped;
  int              L;
  int              M;
  int              do_h2;
  char            *ipsfile = NULL;
  char            *kpsfile = NULL;
  ESL_DMATRIX     *imx     = NULL;
  ESL_DMATRIX     *kmx     = NULL;
  int              ibins;
  
  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  if (esl_opt_IsSet(go, "-h")) {
    puts(usage);
    puts("\n  where options are:\n");
    esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=all docgroups; 2 = indentation; 80=textwidth*/
    return eslOK;
  }
  esl_opt_GetBooleanOption(go, "-i",     &do_ilocal);
  esl_opt_GetStringOption (go, "-m",     &hmmfile);
  esl_opt_GetIntegerOption(go, "-n",     &nseq);
  esl_opt_GetBooleanOption(go, "-u",     &do_ungapped);
  esl_opt_GetIntegerOption(go, "-L",     &L);
  esl_opt_GetIntegerOption(go, "-M",     &M);
  esl_opt_GetBooleanOption(go, "-2",     &do_h2);
  esl_opt_GetStringOption (go, "--ips",  &ipsfile);
  esl_opt_GetStringOption (go, "--kps",  &kpsfile);
  esl_opt_GetIntegerOption(go, "--ibins", &ibins);

  if (esl_opt_ArgNumber(go) != 0) {
    puts("Incorrect number of command line arguments.");
    printf("Usage: %s [options]\n", argv[0]);
    return eslFAIL;
  }

  r   = esl_randomness_CreateTimeseeded();

  if (hmmfile != NULL)
    {	/* Read the HMM (and get alphabet from it) */
      P7_HMMFILE      *hfp     = NULL;

      status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
      if (status == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n", hmmfile);
      else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", hmmfile);
    
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
      if (do_ungapped) p7_hmm_SampleUngapped(r, M, abc, &hmm);
      else             p7_hmm_Sample        (r, M, abc, &hmm);
    }
  imx = esl_dmatrix_Create(ibins, ibins);
  kmx = esl_dmatrix_Create(M, M);
  esl_dmatrix_SetZero(imx);
  esl_dmatrix_SetZero(kmx);
  p7_trace_Create(256, &tr);
  sq = esl_sq_CreateDigital(abc);
  hmm->bg = p7_bg_Create(abc);
  core    = p7_hmm_Duplicate(hmm);

  if (do_h2) {
    hmm->gm = p7_profile_Create(hmm->M, abc);
    p7_H2_ProfileConfig(hmm, hmm->gm, p7_UNILOCAL);
  } else {
    hmm->gm = p7_profile_Create(hmm->M, abc);
    p7_ProfileConfig(hmm, hmm->gm, p7_UNILOCAL);
    p7_ReconfigLength(hmm->gm, L);
    if (p7_hmm_Validate    (hmm,     0.0001, NULL) != eslOK) esl_fatal("whoops, HMM is bad!");
    if (p7_profile_Validate(hmm->gm, 0.0001)       != eslOK) esl_fatal("whoops, profile is bad!");
  }

  for (iseq = 0; iseq < nseq; iseq++)
    {				
      if (do_ilocal) ideal_local_endpoints  (r, core,          sq, tr, ibins, &i1, &i2, &k1, &k2);
      else           emitted_local_endpoints(r, core, hmm->gm, sq, tr, ibins, &i1, &i2, &k1, &k2);

      imx->mx[i1-1][i2-1] += 1.;
      kmx->mx[k1-1][k2-1] += 1.; 
    }

  /* Adjust both mx's to log_2(obs/exp) ratio */
  printf("Before normalization/log-odds:\n");
  printf("i matrix values range from %f to %f\n", dmx_upper_min(imx), dmx_upper_max(imx));
  printf("k matrix values range from %f to %f\n", dmx_upper_min(kmx), dmx_upper_max(kmx));

  dmx_upper_norm(kmx);
  esl_dmx_Scale(kmx, (double) (M*(M+1)/2));
  for (i = 0; i < kmx->m; i++)
    for (j = i; j < kmx->n; j++)
      kmx->mx[i][j] = log(kmx->mx[i][j]) / log(2.0);

  dmx_upper_norm(imx);
  esl_dmx_Scale(imx, (double) (imx->m*(imx->m+1)/2));
  for (i = 0; i < imx->m; i++)
    for (j = i; j < imx->m; j++)
      imx->mx[i][j] = log(imx->mx[i][j]) / log(2.0);
  
  /* Print ps files */
  if (kpsfile != NULL) {
    if ((fp = fopen(kpsfile, "w")) == NULL) esl_fatal("Failed to open output postscript file %s", kpsfile);
    dmx_Visualize(fp, kmx, -2., 2.);
    fclose(fp);
  }
  if (ipsfile != NULL) {
    if ((fp = fopen(ipsfile, "w")) == NULL) esl_fatal("Failed to open output postscript file %s", ipsfile);
    dmx_Visualize(fp, imx, -2., 2.); 
    /* dmx_Visualize(fp, imx, dmx_upper_min(imx), dmx_upper_max(imx)); */
    fclose(fp);
  }

  printf("After normalization/log-odds:\n");
  printf("i matrix values range from %f to %f\n", dmx_upper_min(imx), dmx_upper_max(imx));
  printf("k matrix values range from %f to %f\n", dmx_upper_min(kmx), dmx_upper_max(kmx));

  
  p7_profile_Destroy(hmm->gm);
  p7_bg_Destroy(hmm->bg);
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
 * Incept:    SRE, Fri Jan 26 13:01:34 2007 [Janelia]
 *
 * Purpose:  Implementation of the "two-step" fragment sampling
 *           algorithm, sampling a uniform local fragment w.r.t.
 *           sequence coords, by first sampling a complete
 *           sequence of length L from <hmm>, then choosing
 *           a random fragment <i1..i2> uniformly from all
 *           possible $\frac{L(L+1)/2}$ fragments. 
 *           
 *           The caller also provides an allocated sequence <sq> and
 *           traceback <tr>, to be passed to <p7_CoreEmit()>. They
 *           contain the generated (global) sequence and trace upon
 *           return.
 *           
 *           i endpoints are normalized/discretized to 1..<ibins>, so
 *           we can collate i statistics from sampled sequences of
 *           varying L. However, this introduces a binning artifact
 *           leading to underrepresentation of j=M and
 *           overrepresentation of i=1.  
 *           
 *           Note that this routine is only intended for collecting
 *           endpoint statistics (i1,i2,k1,k2); it does not generate a
 *           local alignment trace. (xref milestone 2, STL11/115).
 *           
 * Returns:  <eslOK> on success; returns normalized/binned sequence coords in
 *           <*ret_i1> and <*ret_i2>, and the model entry/exit coords
 *           in <*ret_k1> and <*ret_k2>. Model coords are defined as
 *           the nodes that emitted residues <i1> and <i2>, so they
 *           must be I or M states, and cannot be D states.
 *           
 * Xref:     STL11/142-143 
 */
static int
ideal_local_endpoints(ESL_RANDOMNESS *r, P7_HMM *hmm, ESL_SQ *sq, P7_TRACE *tr, int ibins,
		      int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2)
{
  int status;
  int tpos;
  int i1, i2;

  if ((status = p7_CoreEmit(r, hmm, sq, tr)) != eslOK) goto ERROR;

  /* a simple way to sample uniformly from upper triangle is by rejection */
  do {
    i1 = 1 + esl_rnd_Choose(r, sq->n);
    i2 = 1 + esl_rnd_Choose(r, sq->n);
  } while (i1 > i2);

  /* Get model coords */
  for (tpos = 0; tpos < tr->N; tpos++)
    if (tr->i[tpos] == i1) { *ret_k1 = tr->k[tpos]; break; }
  for (tpos = tr->N-1; tpos >= 0; tpos--)
    if (tr->i[tpos] == i2) { *ret_k2 = tr->k[tpos]; break; }

  /* Normalize sequence coords.
   * They're 1..L now; make them 1..M
   */
  *ret_i1 = ((i1-1) * ibins / sq->n) + 1;
  *ret_i2 = ((i2-1) * ibins / sq->n) + 1;
  return eslOK;

 ERROR:
  *ret_i1 = 0.;
  *ret_i2 = 0.;
  *ret_k1 = 0;
  *ret_k2 = 0;
  return status;
}

/* emitted_local_endpoints()
 * Incept:    SRE, Fri Jan 26 13:16:09 2007 [Janelia]
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
 *            The model start/end coords <k1..k2> are well-defined, as
 *            the coord of the first/last node in the trace. In
 *            principle, this could be at any state type, but in
 *            practice, HMMER3 allows entry only on M states and exit
 *            from M or D, whereas HMMER2 allowed entry/exit only on M
 *            states.
 *            
 *            Normalized sequence coords <i1..i2> are not
 *            well-defined, because we want to report them relative to
 *            sequence coords of a complete (global) sampled sequence
 *            that we could have sampled this local alignment
 *            from. There are different ways we could try to calculate
 *            the missing sequence 'offsets' flanking our model entry
 *            points <k1> and <k2>. (Note that we do expect variation,
 *            and a dependency on both the statetypes at <k1> and
 *            <k2>.)  The method used here, rejection sampling, is
 *            inefficient: sample sequences from the core
 *            probabilistic model until we obtain one that uses that
 *            same statetypes at entry point <k1> and <k2>, and use
 *            that sample's sequence offsets. I believe this sampling
 *            procedure gives us P(offsets | k1,k2 used). 
 *            
 *            Then, sampling from both the core model and the profile
 *            in the same routine introduces a complication:
 *            conceivably, profile configuration alters the transition
 *            probabilities in the core model (by adding <M->E>
 *            transitions and renormalizing the M transition
 *            distributions, for example). So you can't <CoreSample()>
 *            the <gm->hmm> safely. Instead, the caller must provide a
 *            clean copy of the core model in <core>.
 *            
 * Returns:  <eslOK> on success; returns normalized sequence coords in
 *           <*ret_i1> and <*ret_i2>, and the model entry/exit coords
 *           in <*ret_k1> and <*ret_k2>. 
 *           
 * Xref:     STL11/142-143 
 */
static int
emitted_local_endpoints(ESL_RANDOMNESS *r, P7_HMM *core, P7_PROFILE *gm, ESL_SQ *sq, P7_TRACE *tr, int ibins,
			int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2)
{
  int i1,i2;
  int t1,t2;			/* entry/exit positions in local trace, tr */
  int tg1, tg2;			/* entry/exit positions in global trace, tr2 */
  int tpos;
  int nterm, cterm;		/* offsets at N, C terminus. */
  int L;			/* inferred length from 3-part patching */
  ESL_SQ *sq2   = NULL;
  P7_TRACE *tr2 = NULL;
  int status;
  
  if (gm->mode != p7_UNILOCAL) ESL_XEXCEPTION(eslECONTRACT, "profile must be unilocal");
  if ((sq2 = esl_sq_CreateDigital(gm->abc))  == NULL)   { status = eslEMEM; goto ERROR; }
  if ((status = p7_trace_Create(256, &tr2))  != eslOK)  goto ERROR;

  /* sample local alignment from the implicit model */
  if (gm->h2_mode) {
    if ((status = p7_H2_ProfileEmit(r, gm, sq, tr)) != eslOK) goto ERROR;
  } else {
    if ((status = p7_ProfileEmit(r, gm, sq, tr)) != eslOK) goto ERROR;
  }
    
  /* Get trace coords */
  for (tpos = 0; tpos < tr->N; tpos++)
    if (tr->st[tpos] == p7_STM || tr->st[tpos] == p7_STD) { t1 = tpos; break; }
  for (tpos = tr->N-1; tpos >= 0; tpos--)
    if (tr->st[tpos] == p7_STM || tr->st[tpos] == p7_STD) { t2 = tpos; break; }
  
  /* Convert to model coords (easy) */
  /*
  *ret_k1 = tr->k[t1];
  *ret_k2 = tr->k[t2];
  */

  /* Match a core trace to this local trace by rejection sampling;
   * this lets us calculate sequence offsets; see comments above in preamble
   */
  do {
    if ((status = p7_CoreEmit(r, core, sq2, tr2)) != eslOK) goto ERROR;
    for (tpos = 0; tpos < tr2->N; tpos++)
      if (tr2->k[tpos] == tr->k[t1]) { tg1 = tpos; break; }
    for (tpos = tr2->N-1; tpos >= 0; tpos--)
      if (tr2->k[tpos] == tr->k[t2]) { tg2 = tpos; break; }
  }  while (tr2->st[tg1] != tr->st[t1] && tr2->st[tg2] != tr->st[t2]);

  for (nterm = 0, tpos = 0; tpos < tg1; tpos++) 
    if (tr2->st[tpos] == p7_STM || tr2->st[tpos] == p7_STI) nterm++;
  for (cterm = 0, tpos = tr2->N-1; tpos > tg2; tpos--)
    if (tr2->st[tpos] == p7_STM || tr2->st[tpos] == p7_STI) cterm++;

  /* Now determine endpoint coords (t1 and/or t2 might be sitting on a D) */
  for (; t1 < tr->N; t1++) 
    if (tr->i[t1] > 0) { i1 = tr->i[t1]; *ret_k1 = tr->k[t1]; break; }
  for (; t2 >= 0; t2--) 
    if (tr->i[t2] > 0) { i2 = tr->i[t2]; *ret_k2 = tr->k[t2]; break; }
  if (t2 < t1) ESL_XEXCEPTION(eslEINCONCEIVABLE, "this only happens on an all-D path through profile");

  /* Now, offset and normalize the coords. */
  L  = (i2-i1+1) + nterm + cterm;
  i2 = (i2-i1+1) + nterm;
  i1 = nterm+1;
  /*      printf("reconstructed global length L=%d, model length %d\n", L, gm->M); */
    
  /*  printf("success: L=M, %d\n", L); */
  *ret_i1 = ((i1-1) * ibins / L) + 1;
  *ret_i2 = ((i2-1) * ibins / L) + 1;
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



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
