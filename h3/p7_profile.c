/* Routines for the P7_PROFILE structure - a Plan 7 search profile
 *                                         
 *    1. The P7_PROFILE object: allocation, initialization, destruction.
 *    2. Access methods.
 *    3. Debugging and development code.
 *    4. Unit tests.
 *    5. Test driver.
 *    
 * SRE, Thu Jan 11 15:16:47 2007 [Janelia] [Sufjan Stevens, Illinois]
 * SVN $Id$
 */

#include "p7_config.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"


/*****************************************************************
 * 1. The P7_PROFILE object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  p7_profile_Create()
 * Incept:    SRE, Thu Jan 11 15:53:28 2007 [Janelia]
 *
 * Purpose:   Creates a profile of <M> nodes, for digital alphabet <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    NULL on allocation error.
 *
 * Xref:      STL11/125.
 */
P7_PROFILE *
p7_profile_Create(int M, const ESL_ALPHABET *abc)
{
  P7_PROFILE *gm = NULL;
  int         x;
  int         status;

  /* level 0 */
  ESL_ALLOC(gm, sizeof(P7_PROFILE));
  gm->tsc   = gm->msc = gm->isc = NULL;
  gm->bsc   = gm->esc = NULL;
  gm->begin = gm->end = NULL;
  
  /* level 1 */
  ESL_ALLOC(gm->tsc, sizeof(int *) * 7);       gm->tsc[0] = NULL;
  ESL_ALLOC(gm->msc, sizeof(int *) * abc->Kp); gm->msc[0] = NULL;
  ESL_ALLOC(gm->isc, sizeof(int *) * abc->Kp); gm->isc[0] = NULL;
  ESL_ALLOC(gm->bsc, sizeof(int) * (M+1));      
  ESL_ALLOC(gm->esc, sizeof(int) * (M+1));      

  /* Begin, end may eventually disappear in production
   * code, but we need them in research code for now to
   * be able to emulate & test HMMER2 configurations.
   */
  ESL_ALLOC(gm->begin, sizeof(float) * (M+1));
  ESL_ALLOC(gm->end,   sizeof(float) * (M+1));
  
  /* level 2 */
  ESL_ALLOC(gm->tsc[0], sizeof(int) * 7*M);    
  ESL_ALLOC(gm->msc[0], sizeof(int) * abc->Kp * (M+1));
  ESL_ALLOC(gm->isc[0], sizeof(int) * abc->Kp * M);
  for (x = 1; x < abc->Kp; x++) {
    gm->msc[x] = gm->msc[0] + x * (M+1);
    gm->isc[x] = gm->isc[0] + x * M;
  }
  for (x = 0; x < 7; x++)
    gm->tsc[x] = gm->tsc[0] + x * M;

  /* Initialize some pieces of memory that are never used,
   * only there for indexing convenience.
   */
  gm->tsc[p7_TMM][0] = p7_IMPOSSIBLE; /* node 0 nonexistent, has no transitions  */
  gm->tsc[p7_TMI][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TMD][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TIM][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TII][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TDM][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TDD][0] = p7_IMPOSSIBLE;
  gm->tsc[p7_TDM][1] = p7_IMPOSSIBLE; /* delete state D_1 is wing-retracted */
  gm->tsc[p7_TDD][1] = p7_IMPOSSIBLE;
  for (x = 0; x < abc->Kp; x++) {     /* no emissions from nonexistent M_0, I_0 */
    gm->msc[x][0] = p7_IMPOSSIBLE;
    gm->isc[x][0] = p7_IMPOSSIBLE;
  }
  x = esl_abc_XGetGap(abc);	      /* no emission can emit/score gap characters */
  esl_vec_ISet(gm->msc[x], M+1, p7_IMPOSSIBLE);
  esl_vec_ISet(gm->isc[x], M,   p7_IMPOSSIBLE);
  x = esl_abc_XGetMissing(abc);	      /* no emission can emit/score missing data characters */
  esl_vec_ISet(gm->msc[x], M+1, p7_IMPOSSIBLE);
  esl_vec_ISet(gm->isc[x], M,   p7_IMPOSSIBLE);
  gm->bsc[0] = p7_IMPOSSIBLE;
  gm->esc[0] = p7_IMPOSSIBLE;

  /* Set remaining info
   */
  gm->mode        = p7_NO_MODE;
  gm->M           = M;
  gm->abc         = abc;
  gm->hmm         = NULL;
  gm->bg          = NULL;
  gm->do_lcorrect = FALSE;
  gm->lscore      = 0.;
  gm->h2_mode     = FALSE;
  return gm;

 ERROR:
  p7_profile_Destroy(gm);
  return NULL;
}

/* Function:  p7_profile_Clone()
 * Synopsis:  Duplicates a profile.
 * Incept:    SRE, Mon Jun 25 08:29:23 2007 [Janelia]
 *
 * Purpose:   Duplicate profile <gm>; return a pointer
 *            to the newly allocated copy.
 */
P7_PROFILE *
p7_profile_Clone(const P7_PROFILE *gm)
{
  P7_PROFILE *g2 = NULL;
  int x;

  if ((g2 = p7_profile_Create(gm->M, gm->abc)) == NULL) return NULL;
  g2->mode        = gm->mode;
  g2->bg          = gm->bg;
  g2->do_lcorrect = gm->do_lcorrect;
  g2->lscore      = gm->lscore;
  g2->h2_mode     = gm->h2_mode;

  for (x = 0; x < 7;           x++) esl_vec_ICopy(gm->tsc[x], gm->M,   g2->tsc[x]);
  for (x = 0; x < gm->abc->Kp; x++) esl_vec_ICopy(gm->msc[x], gm->M+1, g2->msc[x]);
  for (x = 0; x < gm->abc->Kp; x++) esl_vec_ICopy(gm->isc[x], gm->M,   g2->isc[x]);
  for (x = 0; x < 2;           x++) esl_vec_ICopy(gm->xsc[x], 2,       g2->xsc[x]);
  esl_vec_ICopy(gm->bsc, gm->M+1, g2->bsc);
  esl_vec_ICopy(gm->esc, gm->M+1, g2->esc);

  for (x = 0; x < 2;      x++) esl_vec_FCopy(gm->xt[x],  2,       g2->xt[x]);
  esl_vec_FCopy(gm->begin, gm->M+1, g2->begin);
  esl_vec_FCopy(gm->end,   gm->M+1, g2->end);

  return g2;
}



/* Function:  p7_profile_SetNullEmissions()
 * Synopsis:  Set all emission scores to zero (experimental).
 * Incept:    SRE, Mon Jun 25 08:12:06 2007 [Janelia]
 *
 * Purpose:   Set all emission scores in profile <gm> to zero.
 *            This makes the profile a null model, with all the same
 *            length distributions as the original model, but
 *            the emission probabilities of the background.
 *            
 *            Written to test the idea that score statistics will be
 *            even better behaved when using a null model with the
 *            same length distribution as the search model.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_profile_SetNullEmissions(P7_PROFILE *gm)
{
  int x;

  /* Canonicals */
  for (x = 0; x <= gm->abc->K; x++) {
    esl_vec_ISet(gm->msc[x], gm->M+1, 0);
    esl_vec_ISet(gm->isc[x], gm->M,   0);
  }
  /* Noncanonicals */
  for (x = gm->abc->K+1; x <= gm->abc->Kp-2; x++) {
    esl_vec_ISet(gm->msc[x], gm->M+1, 0);
    esl_vec_ISet(gm->isc[x], gm->M,   0);
  }
  return eslOK;
}



/* Function:  p7_profile_Destroy()
 * Incept:    SRE, Thu Jan 11 15:54:17 2007 [Janelia]
 *
 * Purpose:   Frees a profile <gm>.
 *
 * Returns:   (void).
 *
 * Xref:      STL11/125.
 */
void
p7_profile_Destroy(P7_PROFILE *gm)
{
  if (gm != NULL) {
    if (gm->tsc   != NULL && gm->tsc[0] != NULL) free(gm->tsc[0]);
    if (gm->msc   != NULL && gm->msc[0] != NULL) free(gm->msc[0]);
    if (gm->isc   != NULL && gm->isc[0] != NULL) free(gm->isc[0]);
    if (gm->tsc   != NULL) free(gm->tsc);
    if (gm->msc   != NULL) free(gm->msc);
    if (gm->isc   != NULL) free(gm->isc);
    if (gm->bsc   != NULL) free(gm->bsc);
    if (gm->esc   != NULL) free(gm->esc);
    if (gm->begin != NULL) free(gm->begin);
    if (gm->end   != NULL) free(gm->end);
  }
  free(gm);
  return;
}

/*****************************************************************
 * 2. Access methods.
 *****************************************************************/


/* Function:  p7_profile_GetTScore()
 * Incept:    SRE, Wed Apr 12 14:20:18 2006 [St. Louis]
 *
 * Purpose:   Convenience function that looks up a transition score in
 *            profile <gm> for a transition from state type <st1> in
 *            node <k1> to state type <st2> in node <k2>. For unique
 *            state types that aren't in nodes (<p7_STS>, for example), the
 *            <k> value is ignored, though it would be customarily passed as 0.
 *            Return the transition score in <ret_tsc>.
 *            
 * Returns:   <eslOK> on success, and <*ret_tsc> contains the requested
 *            transition score.            
 * 
 * Throws:    <eslEINVAL> if a nonexistent transition is requested. Now
 *            <*ret_tsc> is set to 0.
 */
int
p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1, char st2, int k2, int *ret_tsc)
{
  int status;
  int tsc    = 0;

  switch (st1) {
  case p7_STS:  break;
  case p7_STT:  break;

  case p7_STN:
    switch (st2) {
    case p7_STB: tsc =  gm->xsc[p7_XTN][p7_MOVE]; break;
    case p7_STN: tsc =  gm->xsc[p7_XTN][p7_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", 
				p7_hmm_DescribeStatetype(st1),
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STB:
    switch (st2) {
    case p7_STM: tsc = gm->bsc[k2]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", 
				p7_hmm_DescribeStatetype(st1),
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STM:
    switch (st2) {
    case p7_STM: tsc = gm->tsc[p7_TMM][k1]; break;
    case p7_STI: tsc = gm->tsc[p7_TMI][k1]; break;
    case p7_STD: tsc = gm->tsc[p7_TMD][k1]; break;
    case p7_STE: tsc = gm->esc[k1];         break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", 
				p7_hmm_DescribeStatetype(st1), k1,
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STD:
    switch (st2) {
    case p7_STM: tsc = gm->tsc[p7_TDM][k1]; break;
    case p7_STD: tsc = gm->tsc[p7_TDD][k1]; break;
    case p7_STE: tsc = gm->esc[k1];         break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", 
				p7_hmm_DescribeStatetype(st1), k1,
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STI:
    switch (st2) {
    case p7_STM: tsc = gm->tsc[p7_TIM][k1]; break;
    case p7_STI: tsc = gm->tsc[p7_TII][k1]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", 
				p7_hmm_DescribeStatetype(st1), k1,
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STE:
    switch (st2) {
    case p7_STC: tsc = gm->xsc[p7_XTE][p7_MOVE]; break;
    case p7_STJ: tsc = gm->xsc[p7_XTE][p7_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", 
				p7_hmm_DescribeStatetype(st1),
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STJ:
    switch (st2) {
    case p7_STB: tsc = gm->xsc[p7_XTJ][p7_MOVE]; break;
    case p7_STJ: tsc = gm->xsc[p7_XTJ][p7_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", 
				p7_hmm_DescribeStatetype(st1),
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7_STC:
    switch (st2) {
    case p7_STT:  tsc = gm->xsc[p7_XTC][p7_MOVE]; break;
    case p7_STC:  tsc = gm->xsc[p7_XTC][p7_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", 
				p7_hmm_DescribeStatetype(st1),
				p7_hmm_DescribeStatetype(st2));
    }
    break;

  default: ESL_XEXCEPTION(eslEINVAL, "bad state type %d in traceback", st1);
  }

  *ret_tsc = tsc;
  return eslOK;

 ERROR:
  *ret_tsc = 0;
  return status;
}


/*****************************************************************
 * 3. Debugging and development code.
 *****************************************************************/

/* Function:  p7_profile_Validate()
 * Incept:    SRE, Tue Jan 23 13:58:04 2007 [Janelia]
 *
 * Purpose:   Validates the internals of the generic profile structure
 *            <gm>. Probability vectors in the implicit profile
 *            probabilistic model are validated to sum to 1.0 +/- <tol>.
 *            
 *            TODO: currently this function only validates the implicit
 *            model's probabilities, nothing else.
 *            
 * Returns:   <eslOK> if <gm> internals look fine. Returns <eslFAIL>
 *            if something is wrong.
 */
int
p7_profile_Validate(const P7_PROFILE *gm, float tol)
{
  float sum;
  int k,i;

  /* begin[k] should sum to 1.0 over the M(M+1)/2 entries in
   * the implicit model
   */
  for (sum = 0., k = 1; k <= gm->M; k++)
    sum += gm->begin[k] * (gm->M - k + 1);
  if (esl_FCompare(sum, 1.0, tol) != eslOK) return eslFAIL;

  /* end[k] should all be 1.0 in the implicit model
   */
  for (k = 1; k <= gm->M; k++)
    if (gm->end[k] != 1.0) return eslFAIL;

  /* all four xt's should sum to 1.0
   */
  for (i = 0; i < 4; i++)
    if (esl_FCompare(gm->xt[i][p7_MOVE] + gm->xt[i][p7_LOOP], 1.0, tol) != eslOK) return eslFAIL;

  return eslOK;
}

/* Function:  p7_profile_Compare()
 * Synopsis:  Compare two profiles for equality.
 * Incept:    SRE, Thu Jun 21 17:57:56 2007 [Janelia]
 *
 * Purpose:   Compare two profiles <gm1> and <gm2> to each other.
 *            Return <eslOK> if they're identical, and <eslFAIL> if
 *            they differ. Floating-point probabilities are 
 *            compared for equality within a fractional tolerance
 *            <tol>. 
 */
int
p7_profile_Compare(P7_PROFILE *gm1, P7_PROFILE *gm2, float tol)
{
  int x;

  if (gm1->mode != gm2->mode) return eslFAIL;
  if (gm1->M    != gm2->M)    return eslFAIL;

  for (x = 0; x < 7; x++) 
    if (esl_vec_ICompare(gm1->tsc[x], gm2->tsc[x], gm1->M)   != eslOK) return eslFAIL;
  for (x = 0; x < gm1->abc->Kp; x++) {
    if (esl_vec_ICompare(gm1->msc[x], gm2->msc[x], gm1->M+1) != eslOK) return eslFAIL;
    if (esl_vec_ICompare(gm1->isc[x], gm2->isc[x], gm1->M)   != eslOK) return eslFAIL;
  }
  for (x = 0; x < 4; x++)
    if (esl_vec_ICompare(gm1->xsc[x], gm2->xsc[x], 2)        != eslOK) return eslFAIL;
  if (esl_vec_ICompare(gm1->bsc, gm2->bsc, gm1->M+1)         != eslOK) return eslFAIL;
  if (esl_vec_ICompare(gm1->esc, gm2->esc, gm1->M+1)         != eslOK) return eslFAIL;

  for (x = 0; x < 4; x++)
    if (esl_vec_FCompare(gm1->xt[x], gm2->xt[x], 2, tol)     != eslOK) return eslFAIL;
  if (esl_vec_FCompare(gm1->begin, gm2->begin, gm1->M+1, tol)!= eslOK) return eslFAIL;
  if (esl_vec_FCompare(gm1->end,   gm2->end,   gm1->M+1, tol)!= eslOK) return eslFAIL;
  return eslOK;
}


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7PROFILE_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_random.h"

static void
utest_Compare(void)
{
  ESL_RANDOMNESS *r    = esl_randomness_Create(42);
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm  = NULL;
  P7_BG          *bg   = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_PROFILE     *gm2  = NULL;
  int             M    = 200;
  int             L    = 400;

  p7_hmm_Sample(r, M, abc, &hmm); /* master and worker's sampled profiles are identical */
  bg  = p7_bg_Create(abc);
  gm  = p7_profile_Create(hmm->M, abc);
  gm2 = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm,  p7_LOCAL);
  p7_ProfileConfig(hmm, bg, gm2, p7_LOCAL);
  p7_ReconfigLength(gm,  L);
  p7_ReconfigLength(gm2, L);

  if (p7_profile_Compare(gm, gm2, 0.001) != eslOK) p7_Die("identical profile comparison failed");
  
  p7_profile_Destroy(gm);
  p7_profile_Destroy(gm2);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return;
}


#endif /*p7PROFILE_TESTDRIVE*/

/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7PROFILE_TESTDRIVE

/* gcc -o profile_utest -g -Wall -I../easel -L../easel -I. -L. -Dp7PROFILE_TESTDRIVE p7_profile.c -lhmmer -leasel -lm
 * ./profile_utest
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",              0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for p7_profile.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  utest_Compare();

  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7PROFILE_TESTDRIVE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
