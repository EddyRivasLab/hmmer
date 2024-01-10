/* The Viterbi filter (VF).
 * 
 * VF is the second vectorized DP filter in the HMMER's acceleration
 * pipeline, between SSV (single segment ungapped Viterbi filter) and
 * FB (the checkpointed Forward/Backward filter)
 * 
 * The Viterbi filter computes an approximate Viterbi optimal gapped
 * alignment score in limited (16-bit) numeric precision and range.
 * It uses striped SIMD vectorization in a memory-efficient single DP
 * row.
 * 
 * The code here is only a runtime dispatcher. Depending on the
 * processor the code is running on, the runtime dispatcher calls the
 * appropriate (fastest) vector ISA implementation. For the actual
 * VF implementations for each supported ISA, see
 * vitfilter_{sse,avx,avx512...}.c. Most of the documentation for VF
 * is here, rather than repeating it across each ISA implementation.
 * 
 * Contents:
 *    1. h4_vitfilter() stub
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
#include "h4_mode.h"
#include "h4_filtermx.h"

#include "vitfilter.h"

static int vitfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);


/*****************************************************************
 * 1. h4_vitfilter() stub
 *****************************************************************/

/* Function:  h4_vitfilter()
 * Synopsis:  Calculates Viterbi score, vewy vewy fast, in limited precision.
 * Incept:    SRE, Sun 30 Jun 2019 (Lawless, Dear God)
 *
 * Purpose:   Calculates an approximation of the Viterbi score for sequence
 *            <dsq> of length <L> residues, using profile <hmm>, alignment
 *            mode <mo>, and an allocated DP matrix <fx>. Return the 
 *            estimated Viterbi raw score (in bits) in <ret_sc>.
 *            
 *            The <hmm> has its striped vector parameters set (and its
 *            <h4_HASVECS) flag).
 *            
 *            The caller has set the alignment mode <mo>, including
 *            its length parameterization (see <h4_mode_SetLength()>).
 *            However, the Viterbi filter uses local alignment only.
 *            The local vs. glocal parameters of <mo> (i.e.  the B
 *            $\rightarrow$ L | G parameters) are ignored.  (Thus,
 *            <mo> can be a default dual-mode local/glocal multihit
 *            mode for length <L>, and the Viterbi filter will still do
 *            local multihit alignment with it.)
 *            
 *            Caller has allocated DP matrix <fx> to any size. It will
 *            be resized here as needed, using a
 *            <h4_filtermx_Reinit()> call. Caller does not have to
 *            call any reuse or reinit on the <fx> to reuse it for
 *            another DP calculation.
 *            
 *            The Viterbi filter uses a memory-efficient one-row
 *            dynamic programming algorithm, so <fx> doesn't really
 *            contain any useful information upon return. If you want
 *            to see the dynamic programming matrix (for debugging,
 *            say), use <h4_filtermx_SetDumpMode()> to provide an open
 *            stream, and the DP matrix will be written there row by
 *            row. Dumping is only available when the code is compiled
 *            in debugging mode (with a nonzero <eslDEBUGLEVEL>).
 *            
 *            The Viterbi filter has limited numeric range, in 16-bit
 *            integers. Scores will overflow on high-scoring
 *            sequences, which is fine for filtering purposes. On
 *            overflow, returns <eslERANGE> and <*ret_sc> is set to
 *            the maximum representable bitscore, which is a lower
 *            bound on the actual score.
 *
 *            The score is not supposed to underflow, but if it did,
 *            <*ret_sc> to -eslINFINITY and return <eslERANGE>.
 *
 *            Because <eslERANGE> is the return code for both overflow
 *            and a should-not-happen underflow, don't use it to
 *            detect one or the other without also checking <*ret_sc>.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            hmm     - profile HMM, with striped vector params
 *            mo      - alignment mode (algorithm-dependent parameters)
 *            fx      - viterbil filter DP matrix 
 *            ret_sc  - RETURN: Viterbi score (in bits)          
 *
 * Returns:   <eslOK> on success; <*ret_sc> is the bitscore; 
 *            <fx> is likely to have been reallocated.
 *
 *            <eslERANGE> if the score overflows. In this case,
 *            this is a high-scoring hit that passes the filter, 
 *            and <*ret_sc> is set to the maximum representable bitscore.
 *
 *            <eslERANGE> if the score underflows despite me thinking
 *            that I've guaranteed that it can't; <*ret_sc> is set to
 *            <-eslINFINITY>.
 *            
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
int
(*h4_vitfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc) =
  vitfilter_dispatcher;




/*****************************************************************
 * 2. CPU dispatching to vector implementations.
 *****************************************************************/

static int
vitfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc)
{
  /* SRE: for now */
  h4_vitfilter = h4_vitfilter_sse;
  return h4_vitfilter_sse(dsq, L, hmm, mo, fx, ret_sc);

#ifdef eslENABLE_AVX512  // Fastest first.
  if (esl_cpu_has_avx512())
    {
      h4_vitfilter = h4_vitfilter_avx512;
      return h4_vitfilter_avx512(dsq, L, hmm, mo, fx, ret_sc);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      h4_vitfilter = h4_vitfilter_avx;
      return h4_vitfilter_avx(dsq, L, hmm, mo, fx, ret_sc);
    }
#endif

#ifdef eslENABLE_SSE4
  if (esl_cpu_has_sse4())
    {
      h4_vitfilter = h4_vitfilter_sse;
      return h4_vitfilter_sse(dsq, L, hmm, mo, fx, ret_sc);
    }
#endif
  
#ifdef eslENABLE_NEON
  h4_vitfilter = h4_vitfilter_neon;
  return h4_vitfilter_neon(dsq, L, hmm, mo, fx, ret_sc);
#endif

#ifdef eslENABLE_VMX
  h4_vitfilter = h4_vitfilter_vmx;
  return h4_vitfilter_vmx(dsq, L, hmm, mo, fx, ret_sc);
#endif

  esl_fatal("vitfilter_dispatcher found no vector implementation - that shouldn't happen.");
}


/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/
#ifdef h4VITFILTER_BENCHMARK
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
#include "vitfilter.h"

#define SIMDOPTS "--sse,--avx,--avx512"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range   toggles  reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,     NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,     NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0",    NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT, "100000", NULL, "n>0",    NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "--sse",     eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the SSE4 implementation",              0 },
  { "--avx",     eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the AVX2 implementation",              0 },
  { "--avx512",  eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the AVX512 implementation",            0 },
  { "--version", eslARG_NONE,   FALSE, NULL,  NULL,    NULL,  NULL, NULL, "show HMMER version info",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 1, argc, argv, "benchmark driver for Viterbi filter", "[-options] <hmmfile>");
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  H4_FILTERMX    *fx      = h4_filtermx_Create(100);
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ       **dsq     = malloc(N * sizeof(ESL_DSQ *));
  int             i;
  float           vfsc;
  double          Mcs;

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);

  h4_mode_SetLength(mo, L);

  /* Overriding the CPU dispatcher */
  if      (esl_opt_GetBoolean(go, "--sse"))    h4_vitfilter = h4_vitfilter_sse;
  else if (esl_opt_GetBoolean(go, "--avx"))    h4_vitfilter = h4_vitfilter_avx;
  else if (esl_opt_GetBoolean(go, "--avx512")) h4_vitfilter = h4_vitfilter_avx512;

  /* generate and store the random test seqs */
  for (i = 0; i < N; i++)
    {
      dsq[i] = malloc(sizeof(ESL_DSQ) * (L+2));
      esl_rsq_xfIID(rng, hmm->f, abc->K, L, dsq[i]);
    }

  /* benchmark timing */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    h4_vitfilter(dsq[i], L, hmm, mo, fx, &vfsc);
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
  h4_filtermx_Destroy(fx);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif // h4VITFILTER_BENCHMARK


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef h4VITFILTER_TESTDRIVE

#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "h4_refmx.h"

#include "modelsample.h"
#include "reference_dp.h"
#include "simdvec.h"

/* utest_compare_reference()
 * 
 * The "compare_reference" unit test checks that the Viterbi filter
 * gives the same raw score as reference DP, when the model for
 * reference DP is configured to local-only mode and its parameters
 * are specially scaled/rounded to match the precision loss in VF 
 * (see h4_{mode,profile}_SameAsVF).
 * 
 * Randomly sample a profile HMM of length <M> for alphabet <abc>, and
 * compare it to <nseq> iid random sequences.
 * 
 * Assumes that we won't accidentally generate a high-scoring random
 * sequence that overflows VF's limited range. 
 */
static void
utest_compare_reference(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, int M, int L, int nseq)
{
  char         msg[] = "vitfilter compare_reference unit test failed";
  ESL_DSQ     *dsq   = malloc(sizeof(ESL_DSQ) * (L+2));
  H4_PROFILE  *hmm   = NULL;
  H4_PROFILE  *xhmm  = NULL;
  H4_MODE     *mo    = h4_mode_Create();
  H4_MODE     *xmo   = NULL;
  H4_PATH     *pi    = h4_path_Create();
  H4_REFMX    *vit   = h4_refmx_Create(M,L);
  H4_FILTERMX *fx    = h4_filtermx_Create(M);
  float       vsc, vfsc;

  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  if ( h4_mode_SetLocal(mo)              != eslOK) esl_fatal(msg);
  if ( h4_mode_SetLength(mo, L)          != eslOK) esl_fatal(msg);

  if ( h4_profile_SameAsVF(hmm, &xhmm)   != eslOK) esl_fatal(msg);
  if ( h4_mode_SameAsVF(mo, &xmo)        != eslOK) esl_fatal(msg);

  //h4_profile_Dump(stdout, hmm);
  //h4_mode_Dump(stdout, xmo);

  while (nseq--)
    {
      esl_rsq_xfIID(rng, hmm->f, hmm->abc->K, L, dsq);

      //h4_profile_Dump(stdout, hmm);
      //h4_filtermx_SetDumpMode(fx, stdout);

      h4_vitfilter        (dsq, L, hmm,  mo,  fx,      &vfsc);
      h4_reference_Viterbi(dsq, L, xhmm, xmo, vit, pi, &vsc);
      vsc = (vsc / h4_SCALE_W) - h4_3NAT_APPROX - xmo->nullsc;    // SRE note: when we make a planned change to have h4_reference_Viterbi return bitscore, not rawscore, delete the xmo->nullsc term here
      
      //h4_refmx_Dump(stdout, vit);
      //h4_refmx_DumpAsVF(stdout, vit);
      //h4_path_Dump(stdout, pi);
      //printf("M: %6d L: %6d vfsc: %8.2f vsc: %8.2f\n", M, L, vfsc, vsc);

      if (fabs(vfsc - vsc) > 0.001) esl_fatal(msg);

      h4_path_Reuse(pi);
      h4_filtermx_Reuse(fx);
      h4_refmx_Reuse(vit);
    }
  
  free(dsq);
  h4_path_Destroy(pi);
  h4_profile_Destroy(hmm);
  h4_profile_Destroy(xhmm);
  h4_mode_Destroy(mo);
  h4_mode_Destroy(xmo);
  h4_filtermx_Destroy(fx);
  h4_refmx_Destroy(vit);
}
#endif // h4VITFILTER_TESTDRIVE



/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef h4VITFILTER_TESTDRIVE

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
  ESL_GETOPTS    *go   = h4_CreateDefaultApp(options, 0, argc, argv, "vitfilter test driver", "[-options]");
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
#endif // h4VITFILTER_TESTDRIVE



/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef h4VITFILTER_EXAMPLE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_filtermx.h"
#include "h4_hmmfile.h"
#include "h4_profile.h"
#include "h4_mode.h"

#include "general.h"
#include "reference_dp.h"
#include "vitfilter.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version info",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv,
						"example of using the Viterbi filter",
						"[-options] <hmmfile> <seqfile>");
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_PROFILE     *xhmm    = NULL;             // clone of <hmm> with scaled/rounded VF-matching scores
  H4_MODE        *mo      = h4_mode_Create();
  H4_MODE        *xmo     = NULL;             // clone of <mo> with scaled/rounded VF-matching scores
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  H4_REFMX       *vit     = h4_refmx_Create(100,100);
  H4_FILTERMX    *fx      = h4_filtermx_Create(100);
  float           vsc, vfsc;
  int             status;

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);

  h4_mode_SetLocal(mo);

  //h4_profile_Dump(stdout, hmm);
  //h4_mode_Dump(stdout, mo);

  if ( esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, /*env=*/NULL, &sqfp) != eslOK) esl_fatal("couldn't open sequence file %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      h4_mode_SetLength(mo, sq->n);

      h4_profile_SameAsVF(hmm, &xhmm);
      h4_mode_SameAsVF(mo, &xmo);

      //h4_profile_Dump(stdout, xhmm);
      //h4_mode_Dump(stdout, xmo);

      //h4_filtermx_SetDumpMode(fx, stdout);

      h4_vitfilter(sq->dsq, sq->n, hmm, mo, fx, &vfsc);
      printf("%s vfsc %.6f\n", sq->name, vfsc);

      h4_reference_Viterbi(sq->dsq, sq->n, xhmm, xmo, vit, NULL, &vsc);
      printf("%s vsc %.6f\n", sq->name, (vsc / h4_SCALE_W) - h4_3NAT_APPROX);
      //h4_refmx_Dump(stdout, vit);

      h4_refmx_Reuse(vit);
      esl_sq_Reuse(sq);
      h4_mode_Destroy(xmo);
      h4_profile_Destroy(xhmm);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  h4_refmx_Destroy(vit);
  h4_filtermx_Destroy(fx);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
}  
#endif // h4VITFILTER_EXAMPLE

