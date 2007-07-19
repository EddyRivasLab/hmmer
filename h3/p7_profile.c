/* Routines for the P7_PROFILE structure - Plan 7's search profile
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
 * Synopsis:  Create a profile.
 * Incept:    SRE, Thu Jan 11 15:53:28 2007 [Janelia]
 *
 * Purpose:   Creates a profile of <M> nodes, for digital alphabet <abc>.
 *            
 *            Scores (and length model probabilities) are uninitialized;
 *            the <p7_ProfileConfig()> call is what sets these.
 *            The alignment mode is set to <p7_NO_MODE>. 
 *            The reference pointers <gm->hmm_r> and <gm->bg_r> are set to <NULL>.
 *            The reference pointer <gm->abc_r> is set to <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    <NULL> on allocation error.
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
  gm->tsc   = NULL;
  gm->rsc   = NULL;
  
  /* level 1 */
  ESL_ALLOC(gm->tsc,   sizeof(float)   * M * p7P_NTRANS); 
  ESL_ALLOC(gm->rsc,   sizeof(float *) * abc->Kp);
  gm->rsc[0] = NULL;
  
  /* level 2 */
  ESL_ALLOC(gm->rsc[0], sizeof(float) * abc->Kp * (M+1) * p7P_NR);
  for (x = 1; x < abc->Kp; x++) 
    gm->rsc[x] = gm->rsc[0] + x * (M+1) * p7P_NR;

  /* Initialize some edge pieces of memory that are never used,
   * and are only present for indexing convenience.
   */
  esl_vec_FSet(gm->tsc, p7P_NTRANS, -eslINFINITY);     /* node 0 nonexistent, has no transitions  */
  if (M > 1) {
    p7P_TSC(gm, 1, p7P_DM) = -eslINFINITY;             /* delete state D_1 is wing-retracted      */
    p7P_TSC(gm, 1, p7P_DD) = -eslINFINITY;
  }
  for (x = 0; x < abc->Kp; x++) {        
    p7P_MSC(gm, 0, x) = -eslINFINITY;                  /* no emissions from nonexistent M_0... */
    p7P_ISC(gm, 0, x) = -eslINFINITY;                  /* or I_0... */
    p7P_ISC(gm, M, x) = -eslINFINITY;                  /* or I_M.   */
  }
  x = esl_abc_XGetGap(abc);	                       /* no emission can emit/score gap characters */
  esl_vec_FSet(gm->rsc[x], (M+1)*p7P_NR, -eslINFINITY);
  x = esl_abc_XGetMissing(abc);	                      /* no emission can emit/score missing data characters */
  esl_vec_FSet(gm->rsc[x], (M+1)*p7P_NR, -eslINFINITY);

  /* Set remaining info  */
  gm->mode        = p7_NO_MODE;
  gm->M           = M;
  gm->abc_r       = abc;
  gm->hmm_r       = NULL;
  gm->bg_r        = NULL;
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

  if ((g2 = p7_profile_Create(gm->M, gm->abc_r)) == NULL) return NULL;
  g2->mode        = gm->mode;
  g2->hmm_r       = gm->hmm_r;
  g2->bg_r        = gm->bg_r;

  esl_vec_FCopy(gm->tsc, gm->M*p7P_NTRANS, g2->tsc);
  for (x = 0; x < gm->abc_r->Kp;  x++) esl_vec_FCopy(gm->rsc[x], (gm->M+1)*p7P_NR, g2->rsc[x]);
  for (x = 0; x < p7P_NXSTATES;   x++) esl_vec_FCopy(gm->xsc[x], p7P_NXTRANS,      g2->xsc[x]);
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
  for (x = 0; x <= gm->abc_r->K; x++)                  esl_vec_FSet(gm->rsc[x], (gm->M+1)*p7P_NR, 0.0);   /* canonicals    */
  for (x = gm->abc_r->K+1; x <= gm->abc_r->Kp-2; x++)  esl_vec_FSet(gm->rsc[x], (gm->M+1)*p7P_NR, 0.0);   /* noncanonicals */
  return eslOK;
}



/* Function:  p7_profile_Destroy()
 * Synopsis:  Frees a profile.
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
    if (gm->rsc   != NULL && gm->rsc[0] != NULL) free(gm->rsc[0]);
    if (gm->tsc   != NULL) free(gm->tsc);
    if (gm->rsc   != NULL) free(gm->rsc);
    free(gm);
  }
  return;
}


/*****************************************************************
 * 2. Access methods.
 *****************************************************************/

/* Function:  p7_profile_IsLocal()
 * Synopsis:  Return TRUE if profile is in a local alignment mode.
 * Incept:    SRE, Thu Jul 12 11:57:49 2007 [Janelia]
 *
 * Purpose:   Return <TRUE> if profile is in a local alignment mode.
 */
int
p7_profile_IsLocal(const P7_PROFILE *gm)
{
  if (gm->mode == p7_UNILOCAL || gm->mode == p7_LOCAL) return TRUE;
  return FALSE;
}

/* Function:  p7_profile_IsMultihit()
 * Synopsis:  Return TRUE if profile is in a multihit alignment mode.
 * Incept:    SRE, Thu Jul 12 11:58:58 2007 [Janelia]
 *
 * Purpose:   Return <TRUE> if profile is in a multihit alignment mode.
 */
int
p7_profile_IsMultihit(const P7_PROFILE *gm)
{
  if (gm->mode == p7_LOCAL || gm->mode == p7_GLOCAL) return TRUE;
  return FALSE;
}




/* Function:  p7_profile_GetTScore()
 * Incept:    SRE, Wed Apr 12 14:20:18 2006 [St. Louis]
 *
 * Purpose:   Convenience function that looks up a transition score in
 *            profile <gm> for a transition from state type <st1> in
 *            node <k1> to state type <st2> in node <k2>. For unique
 *            state types that aren't in nodes (<p7T_S>, for example), the
 *            <k> value is ignored, though it would be customarily passed as 0.
 *            Return the transition score in <ret_tsc>.
 *            
 * Returns:   <eslOK> on success, and <*ret_tsc> contains the requested
 *            transition score.            
 * 
 * Throws:    <eslEINVAL> if a nonexistent transition is requested. Now
 *            <*ret_tsc> is set to $-\infty$.
 */
int
p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1, char st2, int k2, float *ret_tsc)
{
  int   status;
  float tsc = 0.0f;

  switch (st1) {
  case p7T_S:  break;
  case p7T_T:  break;

  case p7T_N:
    switch (st2) {
    case p7T_B: tsc =  gm->xsc[p7P_N][p7P_MOVE]; break;
    case p7T_N: tsc =  gm->xsc[p7P_N][p7P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DescribeStatetype(st1), p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7T_B:
    switch (st2) {
    case p7T_M:  tsc = p7P_TSC(gm, k2-1, p7P_BM); break; /* remember, B->Mk is stored in [k-1][p7P_BM] */
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DescribeStatetype(st1), p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7T_M:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm, k1, p7P_MM); break;
    case p7T_I: tsc = p7P_TSC(gm, k1, p7P_MI); break;
    case p7T_D: tsc = p7P_TSC(gm, k1, p7P_MD); break;
    case p7T_E: 
      if (k1 != gm->M && ! p7_profile_IsLocal(gm)) ESL_EXCEPTION(eslEINVAL, "local end transition (M%d of %d) in non-local model", k1, gm->M);
      tsc = 0.0f;		/* by def'n in H3 local alignment */
      break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", p7_hmm_DescribeStatetype(st1), k1, p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7T_D:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm, k1, p7P_DM); break;
    case p7T_D: tsc = p7P_TSC(gm, k1, p7P_DD); break;
    case p7T_E: 
      if (k1 != gm->M && ! p7_profile_IsLocal(gm)) ESL_EXCEPTION(eslEINVAL, "local end transition (D%d of %d) in non-local model", k1, gm->M);
      tsc = 0.0f;		/* by def'n in H3 local alignment */
      break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", p7_hmm_DescribeStatetype(st1), k1, p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7T_I:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm, k1, p7P_IM); break;
    case p7T_I: tsc = p7P_TSC(gm, k1, p7P_II); break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", p7_hmm_DescribeStatetype(st1), k1, p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7T_E:
    switch (st2) {
    case p7T_C: tsc = gm->xsc[p7P_E][p7P_MOVE]; break;
    case p7T_J: tsc = gm->xsc[p7P_E][p7P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DescribeStatetype(st1), p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7T_J:
    switch (st2) {
    case p7T_B: tsc = gm->xsc[p7P_J][p7P_MOVE]; break;
    case p7T_J: tsc = gm->xsc[p7P_J][p7P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DescribeStatetype(st1), p7_hmm_DescribeStatetype(st2));
    }
    break;

  case p7T_C:
    switch (st2) {
    case p7T_T:  tsc = gm->xsc[p7P_C][p7P_MOVE]; break;
    case p7T_C:  tsc = gm->xsc[p7P_C][p7P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DescribeStatetype(st1), p7_hmm_DescribeStatetype(st2));
    }
    break;

  default: ESL_XEXCEPTION(eslEINVAL, "bad state type %d in traceback", st1);
  }

  *ret_tsc = tsc;
  return eslOK;

 ERROR:
  *ret_tsc = -eslINFINITY;
  return status;
}


/*****************************************************************
 * 3. Debugging and development code.
 *****************************************************************/

/* Function:  p7_profile_Validate()
 * Incept:    SRE, Tue Jan 23 13:58:04 2007 [Janelia]
 *
 * Purpose:   Validates the internals of the generic profile structure
 *            <gm>.
 *            
 *            TODO: currently this function is a no-op!
 *            
 * Returns:   <eslOK> if <gm> internals look fine. Returns <eslFAIL>
 *            if something is wrong.
 */
int
p7_profile_Validate(const P7_PROFILE *gm, float tol)
{
  int     status;
  int     k;
  double *pstart = NULL;

  ESL_ALLOC(pstart, sizeof(double) * (gm->M+1));
  pstart[0] = 0.0;

  /* Validate the entry distribution.
   * In a glocal model, this is an explicit probability distribution,
   * corresponding to left wing retraction.
   * In a local model, this is an implicit probability distribution,
   * corresponding to the implicit local alignment model, and we have
   * to calculate the M(M+1)/2 fragment probabilities accordingly.
   */
  if (p7_profile_IsLocal(gm))
    {				/* the code block below is also in emit.c:sample_endpoints */
      for (k = 1; k <= gm->M; k++)
	pstart[k] = exp(p7P_TSC(gm, k-1, p7P_BM)) * (gm->M - k + 1); /* multiply p_ij by the number of exits j */
    }
  else
    {
      for (k = 1; k <= gm->M; k++)
	pstart[k] = exp(p7P_TSC(gm, k-1, p7P_BM));
    }

  if (esl_vec_DValidate(pstart, gm->M+1, tol, NULL) != eslOK) { status = eslFAIL; goto ERROR; }
  free(pstart);
  return eslOK;

 ERROR:
  if (pstart != NULL) free(pstart);
  return eslFAIL;
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

  if (esl_vec_FCompare(gm1->tsc, gm2->tsc, gm1->M*p7P_NTRANS, tol)         != eslOK) return eslFAIL;
  for (x = 0; x < gm1->abc_r->Kp; x++) 
    if (esl_vec_FCompare(gm1->rsc[x], gm2->rsc[x], (gm1->M+1)*p7P_NR, tol) != eslOK) return eslFAIL;

  for (x = 0; x < p7P_NXSTATES; x++)
    if (esl_vec_FCompare(gm1->xsc[x], gm2->xsc[x], p7P_NXTRANS, tol)     != eslOK) return eslFAIL;

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
