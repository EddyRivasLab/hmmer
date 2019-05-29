/* Viterbi filter.
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
 *   1. ViterbiFilter() API.
 *   2. CPU dispatching to vector implementations.
 *   3. Benchmark driver.
 *   4. Unit tests.
 *   5. Test driver.
 *   6. Example.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "base/p7_hmmwindow.h"
#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"
#include "dp_vector/vitfilter.h"


static int vitfilter_dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);

/*****************************************************************
 * 1. ViterbiFilter() API
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
int
(*p7_ViterbiFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc) =
  vitfilter_dispatcher;







/*****************************************************************
 * 2. CPU dispatching to vector implementations.
 *****************************************************************/

static int
vitfilter_dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
#ifdef eslENABLE_AVX512  // Fastest first.
  if (esl_cpu_has_avx512())
    {
      p7_ViterbiFilter = p7_ViterbiFilter_avx512;
      return p7_ViterbiFilter_avx512(dsq, L, om, ox, ret_sc);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      p7_ViterbiFilter = p7_ViterbiFilter_avx;
      return p7_ViterbiFilter_avx(dsq, L, om, ox, ret_sc);
    }
#endif

#ifdef eslENABLE_SSE4
  if (esl_cpu_has_sse4())
    {
      p7_ViterbiFilter = p7_ViterbiFilter_sse;
      return p7_ViterbiFilter_sse(dsq, L, om, ox, ret_sc);
    }
#endif
  
#ifdef eslENABLE_NEON
  p7_ViterbiFilter = p7_ViterbiFilter_neon;
  return p7_ViterbiFilter_neon(dsq, L, om, ox, ret_sc);
#endif

#ifdef eslENABLE_VMX
  p7_ViterbiFilter = p7_ViterbiFilter_vmx;
  return p7_ViterbiFilter_vmx(dsq, L, om, ox, ret_sc);
#endif

  p7_Die((char *) "vitfilter_dispatcher found no vector implementation - that shouldn't happen.");
}




/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7VITFILTER_BENCHMARK
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

#define SIMDOPTS "--sse,--avx,--avx512"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range   toggles  reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,     NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,     NULL,  NULL, "-x", "compare scores to generic implementation (debug)", 0 }, 
  { "-s",        eslARG_INT,      "0", NULL, NULL,     NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,     NULL,  NULL, "-c", "equate scores to trusted implementation (debug)",  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0",    NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT, "100000", NULL, "n>0",    NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "--sse",     eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the SSE4 implementation",              0 },
  { "--avx",     eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the AVX2 implementation",              0 },
  { "--avx512",  eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the AVX512 implementation",            0 },

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
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, (char *) "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_FILTERMX    *ox      = NULL;
  P7_REFMX       *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, (char *) "-L");
  int             N       = esl_opt_GetInteger(go,(char *)  "-N");
  ESL_DSQ       **dsq     = (ESL_DSQ **) malloc(N * sizeof(ESL_DSQ *));
  int             i;
  float           sc1, sc2;
  double          Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail((char *) "Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail((char *) "Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, L);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  if (esl_opt_GetBoolean(go, (char *) "-x")) p7_profile_SameAsVF(gm, om->scale_w);

  /* Overriding the CPU dispatcher */
  if      (esl_opt_GetBoolean(go, "--sse"))    p7_ViterbiFilter = p7_ViterbiFilter_sse;
  else if (esl_opt_GetBoolean(go, "--avx"))    p7_ViterbiFilter = p7_ViterbiFilter_avx;
  else if (esl_opt_GetBoolean(go, "--avx512")) p7_ViterbiFilter = p7_ViterbiFilter_avx512;

  ox = p7_filtermx_Create(om->M);
  gx = p7_refmx_Create(gm->M, L);

  for (i = 0; i < N; i++)
    {
      dsq[i] = (ESL_DSQ *) malloc(sizeof(ESL_DSQ) * (L+2));
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq[i]);
    }


  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_ViterbiFilter(dsq[i], L, om, ox, &sc1);   

      if (esl_opt_GetBoolean(go, (char *) "-c")) 
	{
	  p7_ReferenceViterbi(dsq[i], L, gm, gx, NULL, &sc2); 
	  printf("%.4f %.4f\n", sc1, sc2);  
	}

      if (esl_opt_GetBoolean(go, (char *) "-x"))
	{
	  p7_ReferenceViterbi(dsq[i], L, gm, gx, NULL, &sc2); 
	  sc2 /= om->scale_w;
	  if (om->mode == p7_UNILOCAL)   sc2 -= 2.0; /* that's ~ L \log \frac{L}{L+2}, for our NN,CC,JJ */
	  else if (om->mode == p7_LOCAL) sc2 -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
	  printf("%.4f %.4f\n", sc1, sc2);  
	}
    }
  esl_stopwatch_Stop(w);
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) w->elapsed;

  printf("# implementation: ");
  if      (esl_opt_GetBoolean(go, "--sse"))    printf("SSE\n");
  else if (esl_opt_GetBoolean(go, "--avx"))    printf("AVX\n");
  else if (esl_opt_GetBoolean(go, "--avx512")) printf("AVX512\n");
  else                                         printf("%s\n", esl_cpu_Get());
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  printf("# M    = %d\n", gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  for (i = 0; i < N; i++) free(dsq[i]);
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
 * 4. Unit tests.
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
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = (ESL_DSQ *) malloc(sizeof(ESL_DSQ) * (L+2));
  P7_FILTERMX *ox  = p7_filtermx_Create(M);
  P7_REFMX    *gx  = p7_refmx_Create(M, L);
  float sc1, sc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_profile_SameAsVF(gm, om->scale_w);	/* round and scale the scores in <gm> the same as in <om> */

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

      if (fabs(sc1-sc2) > 0.001) esl_fatal("viterbi filter unit test failed: scores differ (%.3f, %.3f)", sc1, sc2);
      
      p7_refmx_Reuse(gx);
      //p7_filtermx_Reuse(ox);
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
 * 5. Test driver
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
  { (char *) "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,(char *)  "show brief help on version and usage",           0 },
  { (char *) "-s",        eslARG_INT,     (char *)  "0", NULL, NULL,  NULL,  NULL, NULL, (char *) "set random number seed to <n>",                  0 },
  { (char *) "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, (char *) "be verbose",                                     0 },
  { (char *) "-L",        eslARG_INT,    (char *) "200", NULL, NULL,  NULL,  NULL, NULL, (char *) "size of random sequences to sample",             0 },
  {(char *)  "-M",        eslARG_INT,    (char *) "145", NULL, NULL,  NULL,  NULL, NULL, (char *) "size of random models to sample",                0 },
  { (char *) "-N",        eslARG_INT,    (char *) "100", NULL, NULL,  NULL,  NULL, NULL, (char *) "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the SSE implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, (char *) "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  int             M    = esl_opt_GetInteger(go, (char *) "-M");
  int             L    = esl_opt_GetInteger(go, (char *) "-L");
  int             N    = esl_opt_GetInteger(go, (char *) "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, (char *) "-v")) printf("ViterbiFilter() tests, DNA\n");
  utest_comparison(r, abc, bg, M, L, N);   
  utest_comparison(r, abc, bg, 1, L, 10);  
  utest_comparison(r, abc, bg, M, 1, 10);  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go,(char *)  "-v")) printf("ViterbiFilter() tests, protein\n");
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
 * 6. Example
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
  { (char *) "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL,(char *)  "show brief help on version and usage",             0 },
  { (char *) "-1",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, (char *) "output in one line awkable format",                0 },
  { (char *) "-P",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, (char *) "output in profmark format",                        0 },
#if eslDEBUGLEVEL > 0
  { (char *) "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, (char *) "dump vector DP matrix for examination (verbose)", 0 },
#endif
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
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail((char *) "Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail((char *) "Failed to read HMM");

  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail((char *) "No such file.");
  else if (status == eslEFORMAT)   p7_Fail((char *) "Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail((char *) "Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail((char *) "Open failed, code %d.", status);

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, sq->n);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  ox = p7_filtermx_Create(gm->M);
#if eslDEBUGLEVEL > 0
  if (esl_opt_GetBoolean(go, (char *) "-D"))
      p7_filtermx_SetDumpMode(ox, stdout, TRUE);         // makes VF dump its DP matrix rows as it goes
#endif
  
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_bg_SetLength(bg,            sq->n);

      p7_ViterbiFilter  (sq->dsq, sq->n, om, ox, &vfraw);
      p7_bg_NullOne (bg, sq->dsq, sq->n, &nullsc);
      vfscore = (vfraw - nullsc) / eslCONST_LOG2;
      P        = esl_gumbel_surv(vfscore,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);

      if (esl_opt_GetBoolean(go, (char *) "-1"))
	{
	  printf("%-30s\t%-20s\t%9.2g\t%7.2f\n", sq->name, hmm->name, P, vfscore);
	}
      else if (esl_opt_GetBoolean(go, (char *) "-P"))
	{ /* output suitable for direct use in profmark benchmark postprocessors: */
	  printf("%g\t%.2f\t%s\t%s\n", P, vfscore, sq->name, hmm->name);
	}
      else
	{
          printf("vector code used:     %s\n",        esl_cpu_Get());
	  printf("target sequence:      %s\n",        sq->name);
	  printf("vit filter raw score: %.2f nats\n", vfraw);
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", vfscore);
	  printf("P-value:              %g\n",        P);
	}
      
      //p7_filtermx_Reuse(ox);
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

