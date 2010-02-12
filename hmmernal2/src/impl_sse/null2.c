/* "null2" model, biased composition correction; SSE implementations.
 * 
 * Contents:
 *   1. Null2 estimation algorithms.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 *
 * SRE, Mon Aug 18 08:31:11 2008 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdlib.h>
#include <string.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse.h"

/*****************************************************************
 * 1. Null2 estimation algorithms.
 *****************************************************************/

/* Function:  p7_Null2_ByExpectation()
 * Synopsis:  Calculate null2 model from posterior probabilities.
 * Incept:    SRE, Mon Aug 18 08:32:55 2008 [Janelia]
 *
 * Purpose:   Identical to <p7_GNull2_ByExpectation()> except that
 *            <om>, <pp> are SSE optimized versions of the profile
 *            and the residue posterior probability matrix. See 
 *            <p7_GNull2_ByExpectation()>  documentation.
 */
int
p7_Null2_ByExpectation(const P7_OPROFILE *om, const P7_OMX *pp, float *null2)
{
  int      M    = om->M;
  int      Ld   = pp->L;
  int      Q    = p7O_NQF(M);
  float    norm;
  __m128  *rp;
  __m128   sv;
  float    xfactor;
  int      i,q,x;
  
  /* Calculate expected # of times that each emitting state was used
   * in generating the Ld residues in this domain.
   * The 0 row in <wrk> is used to hold these numbers.
   */
  memcpy(pp->dpf[0], pp->dpf[1], sizeof(__m128) * 3 * Q);
  pp->xmx[p7X_N][0] = pp->xmx[p7X_N][1];
  pp->xmx[p7X_C][0] = pp->xmx[p7X_C][1]; /* 0.0 */
  pp->xmx[p7X_J][0] = pp->xmx[p7X_J][1]; /* 0.0 */
  for (i = 2; i <= Ld; i++)
    {
      for (q = 0; q < Q; q++)
	{
	  pp->dpf[0][q*3 + p7X_M] = _mm_add_ps(pp->dpf[i][q*3 + p7X_M], pp->dpf[0][q*3 + p7X_M]);
	  pp->dpf[0][q*3 + p7X_I] = _mm_add_ps(pp->dpf[i][q*3 + p7X_I], pp->dpf[0][q*3 + p7X_I]);
	}
      pp->xmx[p7X_N][0] += pp->xmx[p7X_N][i];
      pp->xmx[p7X_C][0] += pp->xmx[p7X_C][i];
      pp->xmx[p7X_J][0] += pp->xmx[p7X_J][i];
    }

  /* Convert those expected #'s to frequencies, to use as posterior weights. */
  norm = 1.0 / (float) Ld;
  sv   = _mm_set1_ps(norm);
  for (q = 0; q < Q; q++)
    {
      pp->dpf[0][q*3 + p7X_M] = _mm_mul_ps(pp->dpf[0][q*3 + p7X_M], sv);
      pp->dpf[0][q*3 + p7X_I] = _mm_mul_ps(pp->dpf[0][q*3 + p7X_I], sv);
    }
  pp->xmx[p7X_N][0] *= norm;
  pp->xmx[p7X_C][0] *= norm;
  pp->xmx[p7X_J][0] *= norm;

  /* Calculate null2's emission odds, by taking posterior weighted sum
   * over all emission vectors used in paths explaining the domain.
   */
  xfactor =  pp->xmx[p7X_N][0] +   pp->xmx[p7X_C][0] +  pp->xmx[p7X_J][0];
  for (x = 0; x < om->abc->K; x++)
    {
      sv = _mm_setzero_ps();
      rp = om->rf[x];
      for (q = 0; q < Q; q++)
	{
	  sv = _mm_add_ps(sv, _mm_mul_ps(pp->dpf[0][q*3 + p7X_M], *rp)); rp++;
	  sv = _mm_add_ps(sv, _mm_mul_ps(pp->dpf[0][q*3 + p7X_I], *rp)); rp++;
	}
      esl_sse_hsum_ps(sv, &(null2[x]));
      null2[x] += xfactor;
    }
  /* now null2[x] = \frac{f_d(x)}{f_0(x)} for all x in alphabet,
   * 0..K-1, where f_d(x) are the ad hoc "null2" residue frequencies
   * for this envelope.
   */

  /* make valid scores for all degeneracies, by averaging the odds ratios. */
  esl_abc_FAvgScVec(om->abc, null2);

  /* ta-da */
  return eslOK;
}


/* Function:  p7_Null2_ByTrace()
 * Synopsis:  Assign null2 scores to an envelope by the sampling method.
 * Incept:    SRE, Mon Aug 18 10:22:49 2008 [Janelia]
 *
 * Purpose:   Identical to <p7_GNull2_ByTrace()> except that
 *            <om>, <wrk> are SSE optimized versions of the profile
 *            and the residue posterior probability matrix. See 
 *            <p7_GNull2_ByTrace()>  documentation.
 */
int
p7_Null2_ByTrace(const P7_OPROFILE *om, const P7_TRACE *tr, int zstart, int zend, P7_OMX *wrk, float *null2)
{
  union { __m128 v; float p[4]; } u;
  int    Q  = p7O_NQF(om->M);
  int    Ld = 0;
  float  norm;
  float  xfactor;
  __m128 sv;
  __m128 *rp;
  int    q, r, s;
  int    x;
  int    z;

  /* We'll use the i=0 row in wrk for working space: dp[0][] and xmx[][0]. */
  for (q = 0; q < Q; q++)
    {
      wrk->dpf[0][q*3 + p7X_M] = _mm_setzero_ps();
      wrk->dpf[0][q*3 + p7X_I] = _mm_setzero_ps();
    }
  wrk->xmx[p7X_N][0] = 0.0;
  wrk->xmx[p7X_C][0] = 0.0;
  wrk->xmx[p7X_J][0] = 0.0;

  /* Calculate emitting state usage in this particular trace segment */
  for (z = zstart; z <= zend; z++)
    {
      if (tr->i[z] == 0) continue; /* quick test for whether this trace elem emitted or not */
      Ld++;
      if (tr->k[z] > 0)	/* must be an M or I */
	{ /* surely there's an easier way? but our workspace is striped, interleaved quads... */
	  s = ( (tr->st[z] == p7T_M) ?  p7X_M : p7X_I);
	  q = p7X_NSCELLS * ( (tr->k[z] - 1) % Q) + p7X_M;
	  r = (tr->k[z] - 1) / Q;
	  u.v            = wrk->dpf[0][q];
	  u.p[r]        += 1.0;	/* all this to increment a count by one! */
	  wrk->dpf[0][q] = u.v;

	}
      else /* emitted an x_i with no k; must be an N,C,J */
	{
	  switch (tr->st[z]) {
	  case p7T_N: wrk->xmx[p7X_N][0] += 1.0; break;
	  case p7T_C: wrk->xmx[p7X_C][0] += 1.0; break;
	  case p7T_J: wrk->xmx[p7X_J][0] += 1.0; break;
	  }
	}
    }
  norm = 1.0 / (float) Ld;
  sv = _mm_set1_ps(norm);
  for (q = 0; q < Q; q++)
    {
      wrk->dpf[0][q*3 + p7X_M] = _mm_mul_ps(wrk->dpf[0][q*3 + p7X_M], sv);
      wrk->dpf[0][q*3 + p7X_I] = _mm_mul_ps(wrk->dpf[0][q*3 + p7X_I], sv);
    }
  wrk->xmx[p7X_N][0] *= norm;
  wrk->xmx[p7X_C][0] *= norm;
  wrk->xmx[p7X_J][0] *= norm;

  /* Calculate null2's emission odds, by taking posterior weighted sum
   * over all emission vectors used in paths explaining the domain.
   */
  xfactor =  wrk->xmx[p7X_N][0] +   wrk->xmx[p7X_C][0] +  wrk->xmx[p7X_J][0];
  for (x = 0; x < om->abc->K; x++)
    {
      sv = _mm_setzero_ps();
      rp = om->rf[x];
      for (q = 0; q < Q; q++)
	{
	  sv = _mm_add_ps(sv, _mm_mul_ps(wrk->dpf[0][q*3 + p7X_M], *rp)); rp++;
	  sv = _mm_add_ps(sv, _mm_mul_ps(wrk->dpf[0][q*3 + p7X_I], *rp)); rp++;
	}
      esl_sse_hsum_ps(sv, &(null2[x]));
      null2[x] += xfactor;
    }

  /* now null2[x] = \frac{f_d(x)}{f_0(x)} for all x in alphabet,
   * 0..K-1, where f_d(x) are the ad hoc "null2" residue frequencies
   * for this envelope.
   */

  /* make valid scores for all degeneracies, by averaging the odds ratios. */
  esl_abc_FAvgScVec(om->abc, null2);

  return eslOK;
}


/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/
#ifdef p7NULL2_BENCHMARK
/*
   icc  -O3 -static -o null2_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7NULL2_BENCHMARK null2.c -lhmmer -leasel -lm 
   ./null2_benchmark    <hmmfile>      Does the expectation version.
   ./null2_benchmark -t <hmmfile>      Does the stochastic-traceback-dependent version. 
                                       (This version isn't really dependent on M, so Mc/s may not be an appropriate measure.)

                       RRM_1 (M=72)       Caudal_act (M=136)     SMC_N (M=1151)
                     -----------------    ------------------     ---------------
        21 Aug 2008   3.00u (480 Mc/s)     5.45u (499 Mc/s)     77.56u (297 Mc/s)
    -t  21 Aug 2008  30.50u  (47 Mc/s)    44.96u  (61 Mc/s)  32.03u*10 ( 72 Mc/s)
             
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                  0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-t",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "benchmark the trace-dependent version of null2",   0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for null2 estimation, SSE version";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox1     = NULL;
  P7_OMX         *ox2     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  float           null2[p7_MAXCODE];
  int             i,j;
  int             nsamples = 200;
  float           fsc, bsc;
  double          Mcs;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);                 p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);    p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);    p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  ox1 = p7_omx_Create(gm->M, L, L);
  ox2 = p7_omx_Create(gm->M, L, L);

  esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  p7_Forward (dsq, L, om, ox1,      &fsc);

  if (esl_opt_GetBoolean(go, "-t"))
    {
      P7_TRACE *tr   = p7_trace_Create();
      float    *n2sc = malloc(sizeof(float) * (L+1));
      int       pos;

      esl_stopwatch_Start(w);
      for (i = 0; i < N; i++)
	{ /* This is approximately what p7_domaindef.c::region_trace_ensemble() is doing: */
	  for (j = 0; j < nsamples; j++)
	    {
	      p7_Null2_ByTrace(om, tr, 0, tr->N-1, ox2, null2);
	      for (pos = 1; pos <= L; pos++)
		n2sc[pos] += null2[dsq[pos]];
	    }
	  for (pos = 1; pos <= L; pos++)
	    n2sc[pos] = logf(n2sc[pos] / nsamples);
	}
      esl_stopwatch_Stop(w);

      free(n2sc);
      p7_trace_Destroy(tr);
    }
  else
    {
      p7_Backward(dsq, L, om, ox1, ox2, &bsc);
      p7_Decoding(om, ox1, ox2, ox2);              

      esl_stopwatch_Start(w);
      for (i = 0; i < N; i++)
	p7_Null2_ByExpectation(om, ox2, null2);
      esl_stopwatch_Stop(w);
    }


  Mcs = (double) N * (double) L * (double) gm->M * 1e-6 / (double) w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_omx_Destroy(ox1);
  p7_omx_Destroy(ox2);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7NULL2_BENCHMARK*/
/*------------------ end, benchmark driver ----------------------*/




/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7NULL2_TESTDRIVE

#endif /*p7NULL2_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7NULL2_TESTDRIVE

#endif /*p7NULL2_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/




/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7NULL2_EXAMPLE

#endif /*p7NULL2_EXAMPLE*/
/*------------------------ example ------------------------------*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/

