/* The SSV (Single Segment Viterbi) filter.
 * 
 * SSV is the first vectorized DP filter in HMMER's acceleration
 * pipeline.  It computes the score of the optimal single-hit local
 * ungapped alignment, in limited (8-bit) numeric range and precision.
 * It works entirely in SIMD registers and does not require any memory
 * allocation.
 * 
 * The code here is only the runtime CPU dispatcher. The ISA-specific
 * SIMD implementations are in ssvfilter_{sse,avx,avx512...}. Documentation 
 * and testing code is here, rather than repeating it across each 
 * ISA-specific implementation.
 * 
 * Contents:
 *    1. h4_ssvfilter() stub
 *    2. Runtime dispatcher
 *    3. Benchmark driver
 *    4. Unit tests
 *    5. Test driver
 *    6. Example
 */
#include <h4_config.h>

#include "easel.h"
#include "esl_cpu.h"

#include "h4_profile.h"

#include "ssvfilter.h"

static int ssvfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc);


/*****************************************************************
 * 1. h4_ssvfilter() stub
 *****************************************************************/

/* Function:  h4_ssvfilter()
 * Synopsis:  The SSV filter, the first step in the acceleration pipeline
 *
 * Purpose:   Calculates approximate SSV score for digital sequence <dsq>
 *            of length <L> residues, using profile <hmm>. Return
 *            the SSV raw score, in bits, in <ret_sc>.
 * 
 *            The <hmm> has its striped vector parameters set (and its 
 *            <h4_HASVECS> flag).
 *            
 *            There is no alignment mode to set. SSV assumes
 *            single-hit local alignment by construction.  There isn't
 *            any DP matrix either; SSV works entirely in SIMD
 *            registers.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow. On overflow, returns
 *            <eslERANGE> and <*ret_sc> is set to the maximum
 *            representable value, which is a lower bound on the
 *            actual score.
 *            
 * Args:      dsq    - digital target sequence 1..L
 *            L      - length of <dsq> in residues
 *            hmm    - profile HMM, with striped vector params
 *            ret_sc - RETURN: SSV raw score in bits
 *
 * Returns:   <eslOK> on success, and <*ret_sc> is the SSV score.
 * 
 *            <eslERANGE> if score overflows limited range. In this
 *            case, this is a high-scoring hit that passes the filter,
 *            and <*ret_sc> is set to its maximum representable value.
 *
 * Throws:    (no abnormal error conditions)
 */
int
(*h4_ssvfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc) = 
  ssvfilter_dispatcher;




/*****************************************************************
 * 2. Runtime CPU dispatching to vector implementations.
 *****************************************************************/

static int 
ssvfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc) 
{
  h4_ssvfilter = h4_ssvfilter_sse;               // SRE FIXME: SSE only for now
  return h4_ssvfilter_sse(dsq, L, hmm, ret_sc);

#ifdef eslENABLE_AVX512  // Fastest first.
  if (esl_cpu_has_avx512())
    {
      h4_ssvfilter = h4_ssvfilter_avx512;
      return h4_ssvfilter_avx512(dsq, L, hmm, ret_sc);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      h4_ssvfilter = h4_ssvfilter_avx;
      return h4_ssvfilter_avx(dsq, L, hmm, ret_sc);
    }
#endif

#ifdef eslENABLE_SSE4
  if (esl_cpu_has_sse4())
    {
      h4_ssvfilter = h4_ssvfilter_sse;
      return h4_ssvfilter_sse(dsq, L, hmm, ret_sc);
    }
#endif
  
#ifdef eslENABLE_NEON
  h4_ssvfilter = h4_ssvfilter_neon;
  return h4_ssvfilter_neon(dsq, L, hmm, ret_sc);
#endif

#ifdef eslENABLE_VMX
  h4_ssvfilter = h4_ssvfilter_vmx;
  return h4_ssvfilter_vmx(dsq, L, hmm, ret_sc);
#endif

  esl_fatal("ssvfilter_dispatcher found no vector implementation - that shouldn't happen.");
}


/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/
#ifdef h4SSVFILTER_BENCHMARK
#include <h4_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"

#include "general.h"
#include "ssvfilter.h"

#define SIMDOPTS "--sse,--avx,--avx512"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range   toggles  reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,     NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,     NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0",    NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,"1000000", NULL, "n>0",    NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "--sse",     eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the SSE4 implementation",              0 },
  { "--avx",     eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the AVX2 implementation",              0 },
  { "--avx512",  eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the AVX512 implementation",            0 },
  { "--version", eslARG_NONE,   FALSE, NULL,  NULL,    NULL,  NULL, NULL, "show HMMER version info",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 1, argc, argv, "benchmark driver for SSV filter", "[-options] <hmmfile>");
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ       **dsq     = malloc(N * sizeof(ESL_DSQ *));
  int             i;
  float           sfsc;
  double          Mcs;

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);

  /* Overriding the CPU dispatcher */
  if      (esl_opt_GetBoolean(go, "--sse"))    h4_ssvfilter = h4_ssvfilter_sse;
  else if (esl_opt_GetBoolean(go, "--avx"))    h4_ssvfilter = h4_ssvfilter_avx;
  else if (esl_opt_GetBoolean(go, "--avx512")) h4_ssvfilter = h4_ssvfilter_avx512;

  /* generate and store the random test seqs */
  for (i = 0; i < N; i++)
    {
      dsq[i] = malloc(sizeof(ESL_DSQ) * (L+2));
      esl_rsq_xfIID(rng, hmm->f, abc->K, L, dsq[i]);
    }

  /* benchmark timing */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    h4_ssvfilter(dsq[i], L, hmm, &sfsc);
  esl_stopwatch_Stop(w);

  Mcs        = (double) N * (double) L * (double) hmm->M * 1e-6 / (double) w->elapsed;
  printf("# implementation: ");
  if      (esl_opt_GetBoolean(go, "--sse"))    printf("SSE\n");
  else if (esl_opt_GetBoolean(go, "--avx"))    printf("AVX\n");
  else if (esl_opt_GetBoolean(go, "--avx512")) printf("AVX512\n");
  else                                         printf("%s\n", esl_cpu_Get());
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n", hmm->M);
  printf("# %.1f Mc/s\n", Mcs);

  for (i = 0; i < N; i++) free(dsq[i]);
  free(dsq);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif // h4SSVFILTER_BENCHMARK




/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef h4SSVFILTER_TESTDRIVE

#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "h4_mode.h"
#include "h4_refmx.h"

#include "modelsample.h"
#include "reference_dp.h"
#include "simdvec.h"

static void
utest_compare_reference(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, int M, int L, int nseq)
{
  char        msg[] = "ssvfilter compare_reference unit test failed";
  ESL_DSQ    *dsq   = malloc(sizeof(ESL_DSQ) * (L+2));
  H4_PROFILE *hmm   = NULL;
  H4_PROFILE *xhmm  = NULL;
  H4_MODE    *mo    = h4_mode_Create();
  H4_MODE    *xmo   = NULL;
  H4_REFMX   *vit   = h4_refmx_Create(M, L);
  float       sfsc, vsc, minsc;
  int         status;

  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  if ( h4_mode_SetUnilocal(mo)           != eslOK) esl_fatal(msg);
  if ( h4_mode_SetLength(mo, L)          != eslOK) esl_fatal(msg);

  if ( h4_profile_SameAsSSV(hmm, &xhmm)  != eslOK) esl_fatal(msg);
  if ( h4_mode_SameAsSSV(mo, &xmo)       != eslOK) esl_fatal(msg);

  /* SSV filter assumes best local alignment diagonal has nonnegative
   * score. For an optimal ungapped local alignment of score D, the
   * SSV score reported is max(0,D) + tNB + tBM + tCT - 2 nats. It's
   * almost always true that there's at least one positive scoring
   * residue somewhere, but on edge cases (M ~ 1, L ~ 1) there may not
   * be. The unit test has to check for the SSV minimal score floor.
   */
  minsc = mo->xsc[h4_N][h4_MOVE] + hmm->tauBM + mo->xsc[h4_C][h4_MOVE] - h4_2NAT_APPROX;

  while (nseq--)
    {
      esl_rsq_xfIID(rng, hmm->f, hmm->abc->K, L, dsq);

      status = h4_ssvfilter(dsq, L, hmm,  &sfsc);
      if (status == eslOK) 
	{      
	  h4_reference_Viterbi(dsq, L, xhmm, xmo, vit, NULL, &vsc);
	  vsc = (vsc / h4_SCALE_B) - h4_2NAT_APPROX;
	  vsc = ESL_MAX(minsc, vsc);  // catch case where emulated ref Vit score is < minimum reportable SSV score.

	  //printf("M: %6d L: %6d sfsc: %8.2f vsc: %8.2f\n", M, L, sfsc, vsc);

	  if (fabs(sfsc - vsc) > 0.001) esl_fatal(msg);
	}
      else if (status != eslERANGE) esl_fatal(msg); // allow SSV filter to overflow, though this should be extremely rare on random profile/seq comparisons

      h4_refmx_Reuse(vit);
    }

  h4_refmx_Destroy(vit);
  h4_mode_Destroy(xmo);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(xhmm);
  h4_profile_Destroy(hmm);
  free(dsq);
}
#endif // h4SSVFILTER_TESTDRIVE


/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef h4SSVFILTER_TESTDRIVE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version info",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = h4_CreateDefaultApp(options, 0, argc, argv, "ssvfilter test driver", "[-options]");
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");

  utest_compare_reference(rng, abc, M, L, N);   
  utest_compare_reference(rng, abc, 1, L, 10);  
  utest_compare_reference(rng, abc, M, 1, 10);  

  esl_alphabet_Destroy(abc);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");

  utest_compare_reference(rng, abc, M, L, N); 
  utest_compare_reference(rng, abc, 1, L, 10);
  utest_compare_reference(rng, abc, M, 1, 10);

  esl_alphabet_Destroy(abc);

  fprintf(stderr, "#  status   = ok\n");

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(rng);
  return eslOK;
}
#endif // h4SSVFILTER_TESTDRIVE


/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef h4SSVFILTER_EXAMPLE
#include <h4_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_cpu.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"

#include "general.h"
#include "ssvfilter.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in one line awkable format",                0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version info",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv, "example of using the SSV filter", "[-options] <hmmfile> <seqfile>");
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();   // SSV doesn't use a mode. We only use this for mo->nullsc.
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  float           sfraw, sfscore;
  int             status;

  /* Read in one HMM */
  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);

  /* Open sequence file */
  if ( esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, /*env=*/NULL, &sqfp) != eslOK) esl_fatal("couldn't open sequence file %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  /* Loop over each sequence in <sqfp> */
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      h4_mode_SetLength(mo, sq->n); 

      h4_ssvfilter(sq->dsq, sq->n, hmm, &sfraw);
      sfscore = sfraw - mo->nullsc;

      if (esl_opt_GetBoolean(go, "-1"))
	{
	  printf("%-30s %10.2f\n", sq->name, sfscore);
	}
      else
	{
          printf("vector code used:     %s\n",        esl_cpu_Get());
	  printf("target sequence:      %s\n",        sq->name);
	  printf("SSV filter raw score: %.2f nats\n", sfraw);
	  printf("null score:           %.2f nats\n", mo->nullsc);
	  printf("per-seq score:        %.2f bits\n", sfscore);
	}      
 
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);

  /* cleanup */
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*h4SSVFILTER_EXAMPLE*/
/*-------------------- end, example ---------------------------*/
