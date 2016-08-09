/* Viterbi filter implementation; SSE version.
 * 
 * This is a SIMD vectorized, striped, interleaved, one-row, reduced
 * precision (epi16) implementation of the Viterbi algorithm.
 * 
 * It calculates a close approximation of the Viterbi score, in
 * limited precision (signed words: 16 bits) and range. It may overflow on
 * high scoring sequences, but this indicates that the sequence is a
 * high-scoring hit worth examining more closely anyway.  It will not
 * underflow, in local alignment mode.
 * 
 * Contents:
 *   1. Viterbi filter implementation.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#if p7_CPU_ARCH == intel 
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#ifdef p7_build_AVX512
 #include <immintrin.h>
 #include "esl_avx_512.h"
#endif
#endif /* intel arch */

#if p7_CPU_ARCH == arm || p7_CPU_ARCH == arm64
#include <arm_neon.h>
#endif
#include "easel.h"
#include "esl_sse.h"

#ifdef p7_build_AVX2
 #include <immintrin.h>
 #include "esl_avx.h"
#endif 

#include "esl_gumbel.h"

#include "base/p7_hmmwindow.h"
#include "search/p7_pipeline.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"
#include "dp_vector/vitfilter.h"

/*****************************************************************
 * 1. Viterbi filter implementation.
 *****************************************************************/

/* Function:  p7_ViterbiFilter()
 * Synopsis:  Calculates Viterbi score, vewy vewy fast, in limited precision.
 *
 * Purpose:   Calculates an approximation of the Viterbi score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and a preallocated one-row DP matrix <ox>. Return the 
 *            estimated Viterbi score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow. 
 *            
 *            <ox> will be resized if needed. It's fine if it was just
 *            <_Reuse()'d> from a previous, smaller profile comparison.
 *            
 *            The model must be in a local alignment mode; other modes
 *            cannot provide the necessary guarantee of no underflow.
 *            
 *            This is a striped SIMD Viterbi implementation using Intel
 *            SSE/SSE2 integer intrinsics \citep{Farrar07}, in reduced
 *            precision (signed words, 16 bits).
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: Viterbi score (in nats)          
 *
 * Returns:   <eslOK> on success;
 *            <eslERANGE> if the score overflows; in this case
 *            <*ret_sc> is <eslINFINITY>, and the sequence can 
 *            be treated as a high-scoring hit.
 *            <ox> may be reallocated.
 *
 * Xref:      [Farrar07] for ideas behind striped SIMD DP.
 *            J2/46-47 for layout of HMMER's striped SIMD DP.
 *            J2/50 for single row DP.
 *            J2/60 for reduced precision (epu8)
 *            J2/65 for initial benchmarking
 *            J2/66 for precision maximization
 *            J4/138-140 for reimplementation in 16-bit precision
 *            J9/110-111 for reimplementation with P7_FILTERMX, memory share w/ checkpointed DP matrix
 *            J10/101 for separating P7_FILTERMX from P7_CHECKPTMX again: don't share these
 */
 /* This filter just dispatches to the appropriate SIMD version.  When possible, call that version directly to
  * save the dispatch overhead.
  */

int
p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
  switch(om->simd){
    case SSE:
      return p7_ViterbiFilter_sse(dsq, L, om, ox, ret_sc);
      break;
    case AVX:
      return p7_ViterbiFilter_avx(dsq, L, om, ox, ret_sc);
      break;
    case AVX512:
      return p7_ViterbiFilter_avx512(dsq, L, om, ox, ret_sc);
      break;
    case NEON:
      return p7_ViterbiFilter_neon(dsq, L, om, ox, ret_sc);
      break;
    case NEON64:
      return p7_ViterbiFilter_neon64(dsq, L, om, ox, ret_sc);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_ViterbiFilter");  
  }

}
/*---------------- end, p7_ViterbiFilter() ----------------------*/



/* Function:  p7_ViterbiFilter_longtarget()
 * Synopsis:  Finds windows within potentially long sequence blocks with Viterbi
 *            scores above threshold (vewy vewy fast, in limited precision)
 *
 * Purpose:   Calculates an approximation of the Viterbi score for regions
 *            of sequence <dsq>, using optimized profile <om>, and a pre-
 *            allocated one-row DP matrix <ox>, and captures the positions
 *            at which such regions exceed the score required to be
 *            significant in the eyes of the calling function (usually
 *            p=0.001).
 *
 *            The resulting landmarks are converted to subsequence
 *            windows by the calling function
 *
 *            The model must be in a local alignment mode; other modes
 *            cannot provide the necessary guarantee of no underflow.
 *
 *            This is a striped SIMD Viterbi implementation using Intel
 *            SSE/SSE2 integer intrinsics \citep{Farrar07}, in reduced
 *            precision (signed words, 16 bits).
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues
 *            om      - optimized profile
 *            ox      - DP matrix
 *            filtersc   - null or bias correction, required for translating a P-value threshold into a score threshold
 *            P          - p-value below which a region is captured as being above threshold
 *            windowlist - RETURN: preallocated array of hit windows (start and end of diagonal) for the above-threshold areas
 *
 * Returns:   <eslOK> on success;
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small, or if
 *            profile isn't in a local alignment mode. (Must be in local
 *            alignment mode because that's what helps us guarantee
 *            limited dynamic range.)
 *
 * Xref:      See p7_ViterbiFilter()
 */
int
p7_ViterbiFilter_longtarget(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox,
                            float filtersc, double P, P7_HMM_WINDOWLIST *windowlist)
{
  switch(om->simd){
    case SSE:
      return p7_ViterbiFilter_longtarget_sse(dsq, L, om, ox, filtersc, P, windowlist);
      break;
    case AVX:
      return p7_ViterbiFilter_longtarget_avx(dsq, L, om, ox, filtersc, P, windowlist);
      break;
    case AVX512:
      return p7_ViterbiFilter_longtarget_avx512(dsq, L, om, ox, filtersc, P, windowlist);
      break;
    case NEON:
      return p7_ViterbiFilter_longtarget_neon(dsq, L, om, ox, filtersc, P, windowlist);
      break;
    case NEON64:
      return p7_ViterbiFilter_longtarget_neon64(dsq, L, om, ox, filtersc, P, windowlist);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_ViterbiFilter_longtarget");  
  }

}


/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7VITFILTER_BENCHMARK
/* -c, -x are used for debugging, testing; see msvfilter.c for explanation 
   ./vitfilter_benchmark <hmmfile>          runs benchmark 
   ./vitfilter_benchmark -N100 -c <hmmfile> compare scores to generic impl
   ./vitfilter_benchmark -N100 -x <hmmfile> compare scores to exact emulation
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-x", "compare scores to generic implementation (debug)", 0 }, 
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-c", "equate scores to trusted implementation (debug)",  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for Viterbi filter";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_FILTERMX    *ox      = NULL;
  P7_REFMX       *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc1, sc2;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, L);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  if (esl_opt_GetBoolean(go, "-x")) p7_profile_SameAsVF(om, gm);

  ox = p7_filtermx_Create(om->M);
  gx = p7_refmx_Create(gm->M, L);

  /* Get a baseline time: how long it takes just to generate the sequences */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Run the benchmark */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_ViterbiFilter(dsq, L, om, ox, &sc1);   

      if (esl_opt_GetBoolean(go, "-c")) 
	{
	  p7_ReferenceViterbi(dsq, L, gm, gx, NULL, &sc2); 
	  printf("%.4f %.4f\n", sc1, sc2);  
	}

      if (esl_opt_GetBoolean(go, "-x"))
	{
	  p7_ReferenceViterbi(dsq, L, gm, gx, NULL, &sc2); 
	  sc2 /= om->scale_w;
	  if (om->mode == p7_UNILOCAL)   sc2 -= 2.0; /* that's ~ L \log \frac{L}{L+2}, for our NN,CC,JJ */
	  else if (om->mode == p7_LOCAL) sc2 -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
	  printf("%.4f %.4f\n", sc1, sc2);  
	}
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_filtermx_Destroy(ox);
  p7_refmx_Destroy(gx);
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
#endif /*p7VITFILTER_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/




/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7VITFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_viterbi.h"

/* utest_comparison()
 * 
 * Check against the reference Viterbi, after configuring  
 * a profile such that its scores will match the roundoffs in 
 * the ViterbiFilter -- p7_profile_SameAsVF().
 * 
 * Sample a random model of length <M>, and score <N> random
 * test sequences of length <L>.
 *
 * We assume that we don't accidentally generate a high-scoring random
 * sequence that overflows ViterbiFilter()'s limited range.
 * 
 */
static void
utest_comparison(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HARDWARE *hw;
  if ((hw = p7_hardware_Create ()) == NULL)  p7_Fail("Couldn't get HW information data structure"); 
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_FILTERMX *ox  = p7_filtermx_Create(M, hw->simd);
  P7_REFMX    *gx  = p7_refmx_Create(M, L);
  float sc1, sc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_profile_SameAsVF(om, gm);	/* round and scale the scores in <gm> the same as in <om> */

#if 0
  p7_oprofile_Dump(stdout, om);                   // dumps the optimized profile
  p7_filtermx_SetDumpMode(ox, stdout, TRUE);      // makes the fast DP algorithms dump their matrices
#endif

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_ViterbiFilter   (dsq, L, om, ox,       &sc1);
      p7_ReferenceViterbi(dsq, L, gm, gx, NULL, &sc2);

#if 0
      p7_refmx_Dump(stdout, gx);   // dumps a generic DP matrix
#endif
      
      sc2 /= om->scale_w;
      sc2 -= 3.0;

      if (fabs(sc1-sc2) > 0.001) esl_fatal("viterbi filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);
      
      p7_refmx_Reuse(gx);
      p7_filtermx_Reuse(ox);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_filtermx_Destroy(ox);
  p7_refmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7VITFILTER_TESTDRIVE*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7VITFILTER_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the SSE implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter() tests, DNA\n");
  utest_comparison(r, abc, bg, M, L, N);   
  utest_comparison(r, abc, bg, 1, L, 10);  
  utest_comparison(r, abc, bg, M, 1, 10);  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter() tests, protein\n");
  utest_comparison(r, abc, bg, M, L, N); 
  utest_comparison(r, abc, bg, 1, L, 10);
  utest_comparison(r, abc, bg, M, 1, 10);

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);

  fprintf(stderr, "#  status = ok\n");
  return eslOK;
}
#endif /*VITFILTER_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/



/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7VITFILTER_EXAMPLE
/* A minimal example.
   Also useful for debugging on small HMMs and sequences.
   ./vitfilter_example <hmmfile> <seqfile>
 */ 
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in one line awkable format",                0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in profmark format",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of Viterbi filter algorithm";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_FILTERMX    *ox      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           vfraw, nullsc, vfscore;
  double          P;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, sq->n);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  ox = p7_filtermx_Create(gm->M);

  /* Useful to place and compile in for debugging: 
     p7_oprofile_Dump(stdout, om);                      // dumps the optimized profile
     p7_filtermx_SetDumpMode(ox, stdout, TRUE);         // makes the fast DP algorithms dump their matrices
     p7_refmx_Dump(stdout, gx);                         // dumps a generic DP matrix
  */

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_bg_SetLength(bg,            sq->n);

      p7_ViterbiFilter  (sq->dsq, sq->n, om, ox, &vfraw);
      p7_bg_NullOne (bg, sq->dsq, sq->n, &nullsc);
      vfscore = (vfraw - nullsc) / eslCONST_LOG2;
      P        = esl_gumbel_surv(vfscore,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);

      if (esl_opt_GetBoolean(go, "-1"))
	{
	  printf("%-30s\t%-20s\t%9.2g\t%7.2f\n", sq->name, hmm->name, P, vfscore);
	}
      else if (esl_opt_GetBoolean(go, "-P"))
	{ /* output suitable for direct use in profmark benchmark postprocessors: */
	  printf("%g\t%.2f\t%s\t%s\n", P, vfscore, sq->name, hmm->name);
	}
      else
	{
	  printf("target sequence:      %s\n",        sq->name);
	  printf("vit filter raw score: %.2f nats\n", vfraw);
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", vfscore);
	  printf("P-value:              %g\n",        P);
	}
      
      p7_filtermx_Reuse(ox);
      esl_sq_Reuse(sq);
    }

  /* cleanup */
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_filtermx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7VITFILTER_EXAMPLE*/
/*-------------------- end, example -----------------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

