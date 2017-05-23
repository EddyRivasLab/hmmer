/* SSV (Single Segment Viterbi) filter
 * 
 * See ssvfilter.md for notes.
 *
 * This is only the CPU dispatch front end, which dispatches a
 * p7_SSVFilter() call to the appropriate vector code. See
 * ssvfilter_{sse,avx...} for the various available vector
 * implementations.

 * Contents:
 *   1. SSVFilter() API.
 *   2. CPU dispatching to vector implementations.
 *   3. Benchmark driver.
 *   4. Unit tests.
 *   5. Test driver.
 *   4. Example.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/ssvfilter.h"

static int ssvfilter_dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);



/*****************************************************************
 * 1. SSVFilter API call
 *****************************************************************/

/* Function:  p7_SSVFilter()
 * Synopsis:  The SSV filter, the first step in the acceleration pipeline
 * Incept:    SRE, Mon May  8 07:37:41 2017 [HHGTTG, The Campaign for Real Time]
 *
 * Purpose:   Calculates approximate SSV score for digital sequence <dsq>
 *            of length <L> residues, using vector profile <om>. Return
 *            the SSV score, in nats, in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring sequences,
 *            but will not underflow.
 *            
 *            The model <om> may be in any mode. Only its match
 *            emission scores will be used. The SSV filter inherently
 *            assumes a singlehit local mode, and uses its own special
 *            state transitions scores, not the scores in the profile.
 *
 * Args:      dsq    - digital target sequence 1..L
 *            L      - length of <dsq> in residues
 *            om     - optimized profile
 *            ret_sc - RETURN: SSV score in nats
 *
 * Returns:   <eslOK> on success, and <*ret_sc> is the SSV score.
 * 
 *            <eslERANGE> if score overflows limited range. In this
 *            case, this is a high-scoring hit that passes the filter,
 *            and <*ret_sc> is set to its maximum allowed value.
 *
 * Throws:    (no abnormal error conditions)
 */
int
(*p7_SSVFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc) = 
  ssvfilter_dispatcher;




/*****************************************************************
 * 2. CPU dispatching to vector implementations.
 *****************************************************************/

static int 
ssvfilter_dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc) 
{
#ifdef eslENABLE_AVX512  // Fastest first.
  if (esl_cpu_has_avx512())
    {
      p7_SSVFilter = p7_SSVFilter_sse;
      return p7_SSVFilter_sse(dsq, L, om, ret_sc);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      p7_SSVFilter = p7_SSVFilter_avx;
      return p7_SSVFilter_avx(dsq, L, om, ret_sc);
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse())
    {
      p7_SSVFilter = p7_SSVFilter_sse;
      return p7_SSVFilter_sse(dsq, L, om, ret_sc);
    }
#endif
  
#ifdef eslENABLE_NEON
  p7_SSVFilter = p7_SSVFilter_neon;
  return p7_SSVFilter_neon(dsq, L, om, ret_sc);
#endif

  //#ifdef eslENABLE_VMX
  //  p7_SSVFilter = p7_SSVFilter_vmx;
  //  return p7_SSVFilter_vmx(dsq, L, om, ret_sc);
  //#endif

  p7_Die("ssvfilter_dispatcher found no vector implementation - that shouldn't happen.");
}



/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/
#ifdef p7SSVFILTER_BENCHMARK

/* The benchmark driver an additional non-benchmarking options
 * to facilitate small-scale (by-eye) comparison of SSV scores against
 * other implementations, for debugging purposes.
 * The -x option compares against an emulation that should give
 * exactly the same scores. The emulation is achieved by jiggering the
 * fp scores in a generic profile to disallow gaps, have the same
 * rounding and precision as the int8_t's SSVFilter() is using, and
 * to make the same post-hoc corrections for the NN, CC, JJ
 * contributions to the final nat score; under these contrived
 * circumstances, p7_ReferenceViterbi() gives the same scores as
 * p7_SSVFilter().
 * 
 * For using -x, you probably also want to limit the number of
 * generated target sequences, using -N10 or -N100 for example.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_cpu.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "baseline version, not production version",         0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "equate scores to trusted implementation (debug)",  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT, "200000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for SSVFilter() implementation";

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
  P7_FILTERMX    *fx      = NULL;
  P7_REFMX       *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ       **dsq     = malloc(N * sizeof(ESL_DSQ *));
  int             i;
  float           sc1, sc2;
  double          Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigUnilocal(gm, hmm, bg, L);  // unilocal, to match SSV model.
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  if (esl_opt_GetBoolean(go, "-x")) p7_profile_SameAsSSV(gm, om->scale_b);

  fx = p7_filtermx_Create(gm->M);
  gx = p7_refmx_Create(gm->M, L);

  /* Generate the random seqs and hold them in memory. */
  for (i = 0; i < N; i++)
    {
      dsq[i] = malloc(sizeof(ESL_DSQ) * (L+2));
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq[i]);
    }

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      if (esl_opt_GetBoolean(go, "-b"))      p7_SSVFilter_base_sse(dsq[i], L, om, fx, &sc1);   
      else                                   p7_SSVFilter(dsq[i], L, om, &sc1);   

      /* -x option: compare generic to fast score in a way that should give exactly the same result */
      if (esl_opt_GetBoolean(go, "-x"))
	{
	  p7_ReferenceViterbi(dsq[i], L, gm, gx, NULL, &sc2); 
	  sc2 /= om->scale_b;
          sc2 -= 2.0;
	  printf("%.4f %.4f\n", sc1, sc2);  
	}
    }
  esl_stopwatch_Stop(w);
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) w->elapsed;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n", gm->M);
  printf("# %.1f Mc/s\n", Mcs);
  printf("# implementation: %s\n", esl_cpu_Get());

  for (i = 0; i < N; i++) free(dsq[i]);
  free(dsq);
  p7_filtermx_Destroy(fx);
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
#endif /*p7SSVFILTER_BENCHMARK*/
/*-------------- end, benchmark driver --------------------------*/

/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7SSVFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

#include "search/modelconfig.h"
#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_viterbi.h"

/* utest_comparison()
 * 
 * SSV is tested against reference Viterbi, after configuring a
 * generic profile such that its scores should match SSV scores,
 * using p7_profile_SameAsSSV().
 *
 * We sample a random model of length <M>, and score <N> random test
 * sequences of length <L>. 
 * 
 * Because SSV scores can overflow, we don't sample high-scoring
 * homologs for this test. If a random sequence happens to overflow,
 * that's ok, we skip the score comparison for those rarities.
 */
static void
utest_comparison(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_REFMX    *gx  = p7_refmx_Create(M, L);
  float        sc1, sc2;
  int          status;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_profile_ConfigUnilocal(gm, hmm, bg, L);            // unilocal, to match SSV model
  p7_profile_SameAsSSV(gm, om->scale_b);

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      status = p7_SSVFilter       (dsq, L, om,           &sc1);
      if      (status == eslERANGE) continue;
      else if (status != eslOK)     esl_fatal("SSVFilter unit test failed.");

      p7_ReferenceViterbi(dsq, L, gm, gx, NULL, &sc2);
      sc2 = sc2 / om->scale_b - 2.0f;
      if (fabs(sc1-sc2) > 0.001) esl_fatal("ssv filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);

      p7_refmx_Reuse(gx);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_refmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7SSVFILTER_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7SSVFILTER_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the SSE SSVFilter() implementation";

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

  if (esl_opt_GetBoolean(go, "-v")) printf("SSVFilter() tests, DNA\n");
  utest_comparison(r, abc, bg, M, L, N);   /* normal sized models */
  utest_comparison(r, abc, bg, 1, L, 10);  /* size 1 models       */
  utest_comparison(r, abc, bg, M, 1, 10);  /* size 1 sequences    */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("SSVFilter() tests, protein\n");
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
#endif /*p7SSVFILTER_TESTDRIVE*/




/*****************************************************************
 * 4. Example
 *****************************************************************/
#ifdef p7SSVFILTER_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_cpu.h"
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
static char banner[] = "example of SSV filter algorithm";

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
  P7_FILTERMX    *fx      = p7_filtermx_Create(100);
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           sfraw, nullsc, sfscore;
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

  //p7_oprofile_Dump(stdout, om);                      // dumps the optimized profile

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_bg_SetLength(bg,            sq->n);
      
      p7_SSVFilter  (sq->dsq, sq->n, om, &sfraw);
      //p7_SSVFilter_base_sse(sq->dsq, sq->n, om, fx, &sfraw);


      p7_bg_NullOne (bg, sq->dsq, sq->n, &nullsc);
      sfscore = (sfraw - nullsc) / eslCONST_LOG2;
      P       = esl_gumbel_surv(sfscore,  om->evparam[p7_SMU],  om->evparam[p7_SLAMBDA]);

      if (esl_opt_GetBoolean(go, "-1"))
	{
	  printf("%-30s\t%-20s\t%9.2g\t%7.2f\n", sq->name, hmm->name, P, sfscore);
	}
      else if (esl_opt_GetBoolean(go, "-P"))
	{ /* output suitable for direct use in profmark benchmark postprocessors: */
	  printf("%g\t%.2f\t%s\t%s\n", P, sfscore, sq->name, hmm->name);
	}
      else
	{
          printf("vector code used:     %s\n",        esl_cpu_Get());
	  printf("target sequence:      %s\n",        sq->name);
	  printf("SSV filter raw score: %.2f nats\n", sfraw);
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", sfscore);
	  printf("P-value:              %g\n",        P);
	}
      
      esl_sq_Reuse(sq);
    }

  /* cleanup */
  p7_filtermx_Destroy(fx);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SSVFILTER_EXAMPLE*/
/*-------------------- end, example ---------------------------*/
