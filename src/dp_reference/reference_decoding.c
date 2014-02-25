/* Reference implementation of posterior decoding, with the dual-mode
 * glocal/local model.
 * 
 * The reference implementation is for testing. It is not for use in
 * HMMER3's main executables. The (more complicated) sparse
 * implementation is the production version.
 * 
 * Contents:
 *   1. Posterior decoding.
 *   2. Unit tests.
 *   3. Test driver.
 *   4. Example.
 *   5. Copyright and license information.
 */

#include "p7_config.h"

#include <math.h>

#include "easel.h"

#include "base/p7_profile.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_decoding.h"

/*****************************************************************
 * 1. Posterior decoding
 *****************************************************************/

/* Function:  p7_ReferenceDecoding()
 * Synopsis:  Reference implementation of posterior decoding.
 *
 * Purpose:   Given previously calculated Forward and Backward matrices
 *            <fwd> and <bck>, for a comparison of query profile <gm>
 *            to digital sequence <dsq> of length <L>, perform
 *            posterior decoding for states responsible for residue
 *            emissions.  (That is, $P(s \mid x_i)$ for each state <s>
 *            that could account for each residue $x_i$.) The resulting
 *            posterior decoding matrix is left in <pp>, which has
 *            been allocated and provided by the caller.
 *            
 * Args:      dsq - digital sequence, 1..L
 *            L   - length of <dsq>
 *            gm  - query profile
 *            fwd - Forward matrix, provided by caller
 *            bck - Backward matrix, provided by caller
 *            pp  - RESULT: posterior decoding matrix
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if a reallocation fails.
 */
int
p7_ReferenceDecoding(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_REFMX *fwd, const P7_REFMX *bck, P7_REFMX *pp)
{
  const float *tsc = gm->tsc;	/* activates the TSC() convenience macro */
  const float *rsc;		/* ptr to current row's residue score x_i vector in <gm> */
  const int    M   = fwd->M;
  float     xJ, xC, xG;
  float    *fwdp;
  float    *bckp;
  float    *ppp;
  float    *ppp2;
  float     denom;
  float     sc;
  int       i;
  int       k;
  int       s;
  float     delta;		/* additions to DGk's, resulting from unfolding wing-retracted entry/exit paths */
  int       status;

  /* Contract checks / argument validation */
  ESL_DASSERT1( (fwd->type == p7R_FORWARD) );
  ESL_DASSERT1( (bck->type == p7R_BACKWARD) );
  ESL_DASSERT1( (fwd->L == bck->L && fwd->M == bck->M && fwd->M == gm->M) );

  /* Reallocation, if needed */
  if (bck != pp && ( status = p7_refmx_GrowTo(pp, gm->M, L)) != eslOK) return status;
  pp->M    = M;
  pp->L    = L;
  pp->type = p7R_DECODING;

  /* On row 0, all main states are 0; initialize them so. set
   * N is 1.0 by definition.
   * Some specials (BLG) can be reached, calculate w/ fwdp,bckp on specials.
   */
  ppp = pp->dp[0];
  for (s = 0; s < p7R_NSCELLS * (M+1); s++) *ppp++ = 0.0;
  fwdp = fwd->dp[0] + (p7R_NSCELLS * (M+1)); /* position on first special of row 0 */
  bckp = bck->dp[0] + (p7R_NSCELLS * (M+1)); 
  sc   = bckp[p7R_N]; /* bck[0][N] is the tot score. Use it to normalize */
  ppp[p7R_E]  = 0.0f;
  ppp[p7R_N]  = 1.0f;
  ppp[p7R_J]  = 0.0f;
  ppp[p7R_B]  = expf(fwdp[p7R_B] + bckp[p7R_B] - sc);
  ppp[p7R_L]  = expf(fwdp[p7R_L] + bckp[p7R_L] - sc);
  ppp[p7R_G]  = expf(fwdp[p7R_G] + bckp[p7R_G] - sc);
  ppp[p7R_C]  = 0.0f;
  ppp[p7R_JJ] = 0.0f;
  ppp[p7R_CC] = 0.0f;
  
  /* xJ/xC/xG hold the previous row i-1 values from forward matrix;
   * needed for decoding emit-on-transition, and G->D1..Dk-1->Mk entry
   */
  xJ = fwdp[p7R_J];     /* i.e. -inf */
  xC = fwdp[p7R_C];	/* i.e. -inf */
  xG = fwdp[p7R_G];

  /* main recursion */
  for (i = 1; i <= L; i++)
    {
      fwdp  = fwd->dp[i]  + M*p7R_NSCELLS;    /* <*fwdp> is on supercell M, we'll wing-retract back from there. */
      bckp  = bck->dp[i]  + M*p7R_NSCELLS;    /* <*bckp> ditto */
      ppp   = pp->dp[i]   + M*p7R_NSCELLS;    /* <*ppp> ditto. We haven't initialized ppp[0] supercell to 0.'s yet    */
      ppp2  = pp->dp[i-1] + M*p7R_NSCELLS;    /* <ppp2> is on PREVIOUS row (i-1) where we need to apply wing-retraction deltas */
      rsc   = gm->rsc[dsq[i]] + M * p7P_NR;   /* <rsc> now points at M,I scores for k=M, and we'll walk it back from there */
      denom = 0.0;

      /* Wing-retracted G->D1..Dk-1->MGk entry paths, right to left. The delta is applied to DGk's on the PREVIOUS row i-1 */
      /* Note we need the residue match emission score, *rsc */
      for (delta=0.0, k = M; k >= 1; k--, ppp2-= p7R_NSCELLS, ppp-= p7R_NSCELLS, fwdp-=p7R_NSCELLS, bckp -= p7R_NSCELLS, rsc -= p7P_NR)
	{
	  ppp2[p7R_DG] += delta;                                                  /* deltas are applied to all DGk on PREVIOUS row  */
	  ppp[p7R_DG]   = expf(fwdp[p7R_DG] + bckp[p7R_DG] - sc);                 /* because we allow decoding to overwrite bck, we may have just overwritten bckp[p7R_DG]. Do not access it again. */
	  delta        += expf(xG + TSC(p7P_GM, k-1) + *rsc + bckp[p7R_MG] - sc); /* the access of bckp[MG] forces this line before the ppp[MG] calculation below */
	}

      /* Those ptrs are now on supercell k=0. Initialize ppp[0] supercell, and move them all to supercell 1 */
      for (s = 0; s < p7R_NSCELLS; s++) ppp[s] = 0.0;
      ppp  += p7R_NSCELLS;
      fwdp += p7R_NSCELLS;
      bckp += p7R_NSCELLS;

      /* [ ML MG IL IG DL . ] for k=1..M */
      for (k = 1; k <= M; k++, ppp+= p7R_NSCELLS, fwdp+= p7R_NSCELLS, bckp += p7R_NSCELLS )
	{
	  ppp[p7R_ML] = expf(fwdp[p7R_ML] + bckp[p7R_ML] - sc); denom += ppp[p7R_ML]; 
	  ppp[p7R_MG] = expf(fwdp[p7R_MG] + bckp[p7R_MG] - sc); denom += ppp[p7R_MG]; 
	  ppp[p7R_IL] = expf(fwdp[p7R_IL] + bckp[p7R_IL] - sc); denom += ppp[p7R_IL]; // at k=M IL=0.0; made so because fwd/bck are -inf
	  ppp[p7R_IG] = expf(fwdp[p7R_IG] + bckp[p7R_IG] - sc); denom += ppp[p7R_IG];
	  ppp[p7R_DL] = expf(fwdp[p7R_DL] + bckp[p7R_DL] - sc);                  
	  /* DG were already done above, with unfolding the entry wing retraction */
	}
      /* now our pointers are on the specials */

      /* [ E N J B L G C JJ CC ] */
      /* JJ, CC must be done first; they access bckp[J,C], which ppp[J,C] calc will clobber */
      ppp[p7R_JJ] = expf(xJ + gm->xsc[p7P_J][p7P_LOOP] + bckp[p7R_J] - sc); xJ = fwdp[p7R_J]; denom += ppp[p7R_JJ];
      ppp[p7R_CC] = expf(xC + gm->xsc[p7P_C][p7P_LOOP] + bckp[p7R_C] - sc); xC = fwdp[p7R_C]; denom += ppp[p7R_CC];
      ppp[p7R_E]  = expf(fwdp[p7R_E] + bckp[p7R_E] - sc);                   
      ppp[p7R_N]  = expf(fwdp[p7R_N] + bckp[p7R_N] - sc); denom += ppp[p7R_N]; /* only NN is possible for i>=1, so N=NN */
      ppp[p7R_J]  = expf(fwdp[p7R_J] + bckp[p7R_J] - sc);
      ppp[p7R_B]  = expf(fwdp[p7R_B] + bckp[p7R_B] - sc);
      ppp[p7R_L]  = expf(fwdp[p7R_L] + bckp[p7R_L] - sc);
      ppp[p7R_G]  = expf(fwdp[p7R_G] + bckp[p7R_G] - sc);                  xG = fwdp[p7R_G];
      ppp[p7R_C]  = expf(fwdp[p7R_C] + bckp[p7R_C] - sc);


      /* Renormalization.
       * 
       * Roundoff error accumulation in F/B is significant. For large
       * target seqs, it isn't unusual to have a whole nat of
       * difference in overall fwd vs. bck raw score; for example, fn3
       * vs. TITIN_HUMAN. Since we only use the bck score (at i=0) for
       * normalization above, pp's would have very large systematic
       * error. 
       * 
       * To squash this error in production code, we renormalize.
       * Because renormalization can hide real errors, we don't
       * renormalize when we've compiled the code for unit testing.
       * Default unit tests don't run large enough M/L to create a
       * lot of error accumulation.
       */
#ifndef p7REFERENCE_DECODING_TESTDRIVE
      denom = 1.0 / denom;	/* multiplication faster than division... */
      ppp = pp->dp[i] + p7R_NSCELLS;
      for (s = 0; s < M*p7R_NSCELLS + p7R_NXCELLS; s++)
	*ppp++ *= denom;
#endif
    }
  
  return eslOK;
}
/*-------------- end, posterior decoding ------------------------*/


/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef p7REFERENCE_DECODING_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_random.h"

#include "base/p7_bg.h"

#include "build/modelsample.h"
#include "search/modelconfig.h"

#include "misc/emit.h"
#include "misc/logsum.h"

#include "dp_reference/reference_fwdback.h"
#include "dp_reference/reference_trace.h"

/* The "overwrite" utest verifies an important wrinkle in the API:
 * we're allowed to overwrite the input Backwards matrix with the new
 * posterior decoding matrix, thus saving a matrix allocation.  The
 * utest samples a random profile, compares it against a generated
 * sequence, decodes it both with and without overwriting, and
 * verifies that the resulting decoding matrices are valid and
 * identical.
 */
static void
utest_overwrite(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L)
{
  char          msg[] = "reference_decoding: overwrite unit test failed";
  ESL_SQ       *sq    = esl_sq_CreateDigital(abc);
  P7_HMM       *hmm   = NULL;
  P7_PROFILE   *gm    = p7_profile_Create(M, abc);
  P7_REFMX     *fwd   = NULL;
  P7_REFMX     *bck   = NULL; 
  P7_REFMX     *pp    = NULL;
  float         tol   = 0.0f;	/* exact match is expected! */
  char          errbuf[eslERRBUFSIZE];

  if ( p7_hmm_Sample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)   != eslOK) esl_fatal(msg);
  
  /* Emit sequence from model, using a length model of <L>;
   * restrict the emitted sequence length to 3x (L+M), arbitrarily, to 
   * keep it down to something reasonable.
   */
  if ( p7_profile_SetLength(gm, L)      != eslOK) esl_fatal(msg);
  do {
    esl_sq_Reuse(sq);
    if (p7_ProfileEmit(rng, hmm, gm, bg, sq, /*tr=*/NULL) != eslOK) esl_fatal(msg);
  } while (sq->n > (gm->M+gm->L) * 3); 

  /* Allocate matrices, and set length models */
  if (p7_profile_SetLength(gm, sq->n)        != eslOK) esl_fatal(msg);
  if ( (fwd = p7_refmx_Create(gm->M, sq->n))  == NULL) esl_fatal(msg);
  if ( (bck = p7_refmx_Create(gm->M, sq->n))  == NULL) esl_fatal(msg);
  if ( (pp  = p7_refmx_Create(gm->M, sq->n))  == NULL) esl_fatal(msg);
  
  /* F/B, then decode both ways */
  if (p7_ReferenceForward (sq->dsq, sq->n, gm, fwd, /*fsc=*/NULL)  != eslOK) esl_fatal(msg);
  if (p7_ReferenceBackward(sq->dsq, sq->n, gm, bck, /*bsc=*/NULL)  != eslOK) esl_fatal(msg);
  if (p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, bck, pp)       != eslOK) esl_fatal(msg); /* <pp> is decoded independently      */
  if (p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, bck, bck)      != eslOK) esl_fatal(msg); /* <bck> is overwritten with decoding */

  /* Tests. */
  if (p7_refmx_Validate(pp,  errbuf) != eslOK) esl_fatal("%s:\n%s\n", msg, errbuf);
  if (p7_refmx_Validate(bck, errbuf) != eslOK) esl_fatal("%s:\n%s\n", msg, errbuf);
  if (p7_refmx_Compare(pp, bck, tol) != eslOK) esl_fatal(msg);

  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_refmx_Destroy(pp);
  p7_refmx_Destroy(bck);
  p7_refmx_Destroy(fwd);
}

/* The "rowsum" utest verifies that the sum of posterior probs for all
 * emitters (M,I,NN,CC,JJ) on any row 1..L is ~1.0, for comparisons of
 * a randomly sampled profile to <N> randomly sampled homologs.
 */
static void
utest_rowsum(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char          msg[] = "reference_decoding: rowsum unit test failed";
  ESL_SQ       *sq    = esl_sq_CreateDigital(abc);
  P7_HMM       *hmm   = NULL;
  P7_PROFILE   *gm    = p7_profile_Create(M, abc);
  P7_REFMX     *fwd   = p7_refmx_Create(M, L);
  P7_REFMX     *pp    = p7_refmx_Create(M, L);
  float         tol   = ( p7_logsum_IsSlowExact() ? 0.005 : 0.03); /* tuned to the test's default <M>,<L>,<N> */
  int           idx;
  int           i,k;
  float         rowsum;
  float        *dpc;

  if ( p7_hmm_Sample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)   != eslOK) esl_fatal(msg);
  
  for (idx = 0; idx < N; idx++)
    {
      /* Emit sequence from model, using a length model of <L>;
       * restrict the emitted sequence length to 3x (L+M), arbitrarily, to 
       * keep it down to something reasonable.
       */
      if ( p7_profile_SetLength(gm, L)      != eslOK) esl_fatal(msg);
      do {
	esl_sq_Reuse(sq);
	if (p7_ProfileEmit(rng, hmm, gm, bg, sq, /*tr=*/NULL) != eslOK) esl_fatal(msg);
      } while (sq->n > (gm->M+gm->L) * 3); 


      /* F/B, then decode in place */
      if ( p7_profile_SetLength(gm, sq->n)        != eslOK) esl_fatal(msg);
      if (p7_ReferenceForward (sq->dsq, sq->n, gm, fwd, /*fsc=*/NULL) != eslOK) esl_fatal(msg);
      if (p7_ReferenceBackward(sq->dsq, sq->n, gm, pp, /*bsc=*/NULL)  != eslOK) esl_fatal(msg); /* <pp> is temporarily the Bck matrix */
      if (p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, pp, pp)       != eslOK) esl_fatal(msg); 

      for (i = 1; i <= sq->n; i++)
	{
	  for (rowsum=0., dpc=pp->dp[i], k=0; k <= gm->M; k++, dpc+=p7R_NSCELLS) /* k=0 included; why not, it's all 0's. */
	    rowsum += dpc[p7R_ML] + dpc[p7R_MG] + dpc[p7R_IL] + dpc[p7R_IG];
	  rowsum += dpc[p7R_N] + dpc[p7R_JJ] + dpc[p7R_CC];

	  if (esl_FCompareAbs(rowsum, 1.0, tol) != eslOK) esl_fatal(msg);
	}
    }

  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_refmx_Destroy(pp);
  p7_refmx_Destroy(fwd);
}


/* The "colsum" unit test verifies that for a uniglocal model,
 * which forces each consensus position k (1..M) to be visited
 * once and only once in Mk or Dk, the sum of posterior probs
 * over a column k is ~1.0; using a randomly sampled profile,
 * to <N> randomly sampled homologs.
 *
 * This (\sum_i = 1.0) is also true for E,B,G, which can only
 * be visited once in a uniglocal path. L, J and JJ must be 0.0,
 * because the model is uniglocal.
 */
static void
utest_colsum(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char          msg[] = "reference_decoding: rowsum unit test failed";
  ESL_SQ       *sq    = esl_sq_CreateDigital(abc);
  P7_HMM       *hmm   = NULL;
  P7_PROFILE   *gm    = p7_profile_Create(M, abc);
  P7_REFMX     *fwd   = p7_refmx_Create(M, L);
  P7_REFMX     *pp    = p7_refmx_Create(M, L);
  float         tol   = ( p7_logsum_IsSlowExact() ? 0.003 : 0.02); /* tuned to the test's default <M>,<L>,<N> */
  int           idx;
  int           i,k,s;
  float         xsum[p7R_NXCELLS];
  float        *colsum = malloc(sizeof(float) * (M+1));
  float        *dpc;

  if ( p7_hmm_Sample(rng, M, abc, &hmm)           != eslOK) esl_fatal(msg);
  if ( p7_profile_ConfigUniglocal(gm, hmm, bg, L) != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      /* Emit sequence from model, using a length model of <L>;
       * restrict the emitted sequence length to 3x (L+M), arbitrarily, to 
       * keep it down to something reasonable.
       */
      if ( p7_profile_SetLength(gm, L)      != eslOK) esl_fatal(msg);
      do {
	esl_sq_Reuse(sq);
	if (p7_ProfileEmit(rng, hmm, gm, bg, sq, /*tr=*/NULL) != eslOK) esl_fatal(msg);
      } while (sq->n > (gm->M+gm->L) * 3); 

      /* F/B, then decode in place */
      if ( p7_profile_SetLength(gm, sq->n) != eslOK) esl_fatal(msg);
      if (p7_ReferenceForward (sq->dsq, sq->n, gm, fwd, /*fsc=*/NULL) != eslOK) esl_fatal(msg);
      if (p7_ReferenceBackward(sq->dsq, sq->n, gm, pp, /*bsc=*/NULL)  != eslOK) esl_fatal(msg); /* <pp> is temporarily the Bck matrix */
      if (p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, pp, pp)       != eslOK) esl_fatal(msg); 

      /* Column sums */
      for (k = 0; k <= M;          k++) colsum[k] = 0.0f;
      for (s = 0; s < p7R_NXCELLS; s++) xsum[s] = 0.0f;
      for (i = 0; i <= sq->n; i++)
	{
	  for (dpc=pp->dp[i], k=0; k <= gm->M; k++, dpc+=p7R_NSCELLS) /* k=0 included; why not, it's all 0's. */
	    colsum[k] += dpc[p7R_ML] + dpc[p7R_MG] + dpc[p7R_DL] + dpc[p7R_DG];
	  for (s = 0; s < p7R_NXCELLS; s++)
	    xsum[s]   += dpc[s];
	}

      /* Tests */
      for (k = 1; k <= M; k++) 
	if (esl_FCompareAbs(colsum[k], 1.0, tol) != eslOK) esl_fatal(msg);
      if (esl_FCompareAbs(xsum[p7R_E], 1.0, tol) != eslOK) esl_fatal(msg);
      if (esl_FCompareAbs(xsum[p7R_B], 1.0, tol) != eslOK) esl_fatal(msg);
      if (esl_FCompareAbs(xsum[p7R_G], 1.0, tol) != eslOK) esl_fatal(msg);

      if (xsum[p7R_L]  != 0.0f)  esl_fatal(msg);
      if (xsum[p7R_J]  != 0.0f)  esl_fatal(msg);
      if (xsum[p7R_JJ] != 0.0f)  esl_fatal(msg);
    }

  free(colsum);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_refmx_Destroy(pp);
  p7_refmx_Destroy(fwd);
}


/* The 'approx-decoding' utest compares exact posterior decoding (via
 * p7_ReferenceDecoding()) to stochastic approximation (by a large
 * ensemble of stochastic tracebacks). It does this for a randomly
 * sampled profile HMM of length <M> compared against one homologous
 * (generated) sequence. (Only one sequence, not <N>, because the
 * stochastic tracebacks are computationally expensive.)
 * 
 * Tests:
 * 1. The two decoding approaches give identical matrices within 
 *    a given sampling error tolerance. (Additionally, cells that
 *    are exactly zero in exact posterior decoding must not be
 *    visited in any stochastic trace.) All this is checked by 
 *    <p7_refmx_CompareDecoding()>.
 * 2. The two decoding matrices both Validate().
 * 
 */
static void
utest_approx_decoding(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L)
{
  char           msg[]  = "reference_decoding: approx-decoding unit test failed";
  P7_HMM        *hmm    = NULL;
  P7_PROFILE    *gm     = p7_profile_Create(M, abc);
  ESL_SQ        *sq     = esl_sq_CreateDigital(abc);       /* space for generated (homologous) target seqs              */
  P7_TRACE      *tr     = NULL;
  P7_REFMX      *fwd    = NULL;
  P7_REFMX      *bck    = NULL;
  P7_REFMX      *ppe    = NULL;
  P7_REFMX      *ppa    = NULL;
  float         *wrk    = NULL;	/* reusable scratch workspace needed by stochastic trace */
  int            idx;
  int            ntr    = 100000;
  float          tol    = 0.01;	/* with utest's defaults, max diff will be ~+/- 0.004; exact vs. approx logsum seems to make no difference    */
  char           errbuf[eslERRBUFSIZE];

  /* Sample a profile. 
   * Config as usual: multihit dual-mode local/glocal, so all paths in it are valid.
   */
  if ( p7_hmm_Sample(rng, M, abc, &hmm)  != eslOK) esl_fatal(msg);  
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);
  if ( p7_profile_SetLength(gm, L)       != eslOK) esl_fatal(msg);  

  /* Generate (sample) a sequence from the profile */
  do {   
    esl_sq_Reuse(sq);
    p7_ProfileEmit(rng, hmm, gm, bg, sq, NULL);
  } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
  if ( p7_profile_SetLength(gm, sq->n)       != eslOK) esl_fatal(msg);

  /*  DP calculations, and exact posterior decoding */
  if ( (fwd = p7_refmx_Create(gm->M, sq->n)) == NULL)  esl_fatal(msg);
  if ( (bck = p7_refmx_Create(gm->M, sq->n)) == NULL)  esl_fatal(msg);
  if ( (ppe = p7_refmx_Create(gm->M, sq->n)) == NULL)  esl_fatal(msg);

  if ( p7_ReferenceForward (sq->dsq, sq->n, gm, fwd, NULL)     != eslOK) esl_fatal(msg);
  if ( p7_ReferenceBackward(sq->dsq, sq->n, gm, bck, NULL)     != eslOK) esl_fatal(msg);
  if ( p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, bck, ppe) != eslOK) esl_fatal(msg);

  /* Approximate decoding by stochastic traceback  */
  if ( (ppa = p7_refmx_Create(gm->M, sq->n)) == NULL) esl_fatal(msg);
  if ( (tr  = p7_trace_Create())             == NULL) esl_fatal(msg);
  if ( p7_refmx_Zero(ppa, gm->M, sq->n)     != eslOK) esl_fatal(msg);
  for (idx = 0; idx < ntr; idx++)
    {
      if ( p7_reference_trace_Stochastic(rng, &wrk, gm, fwd, tr) != eslOK) esl_fatal(msg);
      if ( p7_refmx_CountTrace(tr, ppa)                          != eslOK) esl_fatal(msg);
      if ( p7_trace_Reuse(tr)                                    != eslOK) esl_fatal(msg);
    }
  p7_refmx_Rescale(ppa, 1./(float)ntr);

  //p7_refmx_Dump(stdout, ppe);
  //p7_refmx_Dump(stdout, ppa);

  /* Tests */
  if ( p7_refmx_CompareDecoding(ppe, ppa, tol) != eslOK) esl_fatal(msg);
  if ( p7_refmx_Validate(ppe, errbuf)          != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
  if ( p7_refmx_Validate(ppa, errbuf)          != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
  
  if (wrk) free(wrk);
  p7_refmx_Destroy(ppa);
  p7_refmx_Destroy(ppe);
  p7_refmx_Destroy(bck);
  p7_refmx_Destroy(fwd);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  esl_sq_Destroy(sq);
}
#endif /*p7REFERENCE_DECODING_TESTDRIVE*/
/*---------------- end, unit tests ------------------------------*/


/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7REFERENCE_DECODING_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,     "20", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for sparse DP of dual-mode profile";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg   = p7_bg_Create(abc);
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  utest_approx_decoding(r, abc, bg, M, L);

  utest_overwrite      (r, abc, bg, M, L);
  utest_rowsum         (r, abc, bg, M, L, N);
  utest_colsum         (r, abc, bg, M, L, N);

  fprintf(stderr, "#  status = ok\n");

  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}

#endif /*p7REFERENCE_DECODING_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/


/*****************************************************************
 * 8. Example
 *****************************************************************/
#ifdef p7REFERENCE_DECODING_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",              0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                     0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Forward DP matrix for examination",            0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Backward DP matrix for examination",           0 },
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Decoding matrix for examination",              0 },
  { "-X",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump best decoding for each residue position",      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of posterior decoding, reference dual implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_HMM         *hmm     = NULL;
  ESL_SQ         *sq      = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_REFMX       *fwd     = NULL;
  P7_REFMX       *pp      = NULL;
  int           i,k;
  float         fsc, bsc;
  float         rowsum;
  float        *dpc;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Open sequence database */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Read in one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);
  esl_sqfile_Close(sqfp);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  //p7_profile_ConfigGlocal(gm, hmm, bg, sq->n);
  p7_bg_SetLength     (bg, sq->n);
  p7_profile_SetLength(gm, sq->n);

  /* Create matrices */
  fwd = p7_refmx_Create(gm->M, sq->n);
  pp  = p7_refmx_Create(gm->M, sq->n);

  /*  DP calculations */
  p7_ReferenceForward (sq->dsq, sq->n, gm, fwd, &fsc);     if (esl_opt_GetBoolean(go, "-F")) p7_refmx_Dump(stdout, fwd);
  p7_ReferenceBackward(sq->dsq, sq->n, gm, pp,  &bsc);     if (esl_opt_GetBoolean(go, "-B")) p7_refmx_Dump(stdout, pp);
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, pp, pp);   if (esl_opt_GetBoolean(go, "-D")) p7_refmx_Dump(stdout, pp);

  if (esl_opt_GetBoolean(go, "-X")) p7_refmx_DumpBestDecoding(stdout, sq->dsq, sq->n, gm, pp);

  printf("# Forward raw nat score = %9.3f\n", fsc);
  printf("# Backward raw nat score = %9.3f\n", bsc);
  
#if 0
  /* Row sums. */
  for (i = 1; i <= sq->n; i++)
    {
      for (rowsum=0., dpc=pp->dp[i], k=0; k <= gm->M; k++, dpc+=p7R_NSCELLS) /* k=0 included; why not, it's all 0's. */
	rowsum += dpc[p7R_ML] + dpc[p7R_MG] + dpc[p7R_IL] + dpc[p7R_IG];
      rowsum += dpc[p7R_N] + dpc[p7R_JJ] + dpc[p7R_CC];

      printf("row %3d : %f\n", i, rowsum-1.0f);
    }
#endif

  esl_sq_Destroy(sq);
  p7_refmx_Destroy(pp);
  p7_refmx_Destroy(fwd);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_DECODING_EXAMPLE*/
/*------------------ end, example driver ------------------------*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
