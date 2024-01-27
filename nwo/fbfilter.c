/* The Forward/Backward filter (FB).
 *
 * FB is the third vectorized DP filter in the HMMER acceleration
 * pipeline, after SSV and VF.
 *
 * The Forward pass of the FB filter computes a local (only) ensemble
 * Forward score. If this score is deemed sufficient, then we do a
 * second pass with Backward and posterior decoding to identify i,k
 * lattice cells (sequence positions i, profile positions k) that have
 * more than a threshold amount of probability mass. These cells are
 * marked in a data structure (H4_SPARSEMASK) for subsequent sparse
 * dynamic programming with the full glocal/local dual-mode model.
 *
 * The FB filter is memory-efficient, using checkpointed dynamic
 * programming. It requires $O(M \sqrt L)$ memory for a profile of
 * length $M$ and a sequence of length $L$.
 *
 * Denormalized floating point computation must be turned off for best
 * performance on Intel/x86 platforms. See notes in simdvec.md.
 *
 * The code here is only the runtime dispatcher. The actual
 * implementations are in fbfilter_{sse,avx,avx512...}.c.
 *
 * See fwdfilter.md for notes.
 *
 * Contents:
 *    1. h4_fwdfilter(), h4_bckfilter()
 *    2. CPU dispatching to vector implementations.
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
#include "h4_checkptmx.h"
#include "h4_sparsemask.h"

#include "fbfilter.h"

static int fwdfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc);
static int bckfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh);


/*****************************************************************
 * 1. h4_fwdfilter(), h4_bckfilter() 
 *****************************************************************/

/* Function:  h4_fwdfilter()
 * Synopsis:  the Forward acceleration filter 
 *
 * Purpose:   Calculates the local Forward filter score for a comparison
 *            of profile <hmm> to digital sequence <dsq> of length
 *            <L>, using alignment mode <mo>, and checkpoint DP matrix
 *            <cpx>.  Returns the Forward filter score (in bits) in
 *            <ret_sc>.
 *
 *            The <hmm> has its striped vector parameters set (and its
 *            <h4_HASVECS) flag).
 *            
 *            The Forward filter uses multihit local alignment. The
 *            caller provides a multihit alignment mode <mo>, with its
 *            length parameterization set
 *            (<h4_mode_SetLength()>). Local alignment is hardcoded; B
 *            $\rightarrow$ L|G local/glocal parameterization is
 *            ignored, implicitly assuming B->L = 1.0.  Multihit mode
 *            is required to guarantee that scores will not underflow
 *            in sparse rescaling [Eddy11].
 *
 *            (Thus, <mo> can be our default dual-mode local/glocal
 *            multihit mode for length <L>, and the Forward filter
 *            will still do **local** multihit alignment with it.)
 *            
 *            Caller provides an allocated DP matrix <cpx> of any
 *            size. It is resized here as needed. Caller does not have
 *            to call any reuse or reinit on <cpx> to reuse it for
 *            another DP calculation.
 *    
 *            The Forward filter uses a memory-efficient checkpointed
 *            DP algorithm that requires $O(M \sqrt L)$ memory.
 *            
 *
 * Args:      dsq    : digital target sequence 1..L
 *            L      : length of <dsq>
 *            hmm    : profile HMM, with striped vector params
 *            mo     : alignment mode, multihit, with length L set
 *            cpx    : checkpointed DP matrix, any allocation is fine
 *            ret_sc : RETURN: Forward filter score (bits)
 *
 * Returns:   <eslOK> on success;
 *            <*ret_sc> is the Forward filter score in bits;
 *            <cpx> contains a checkpointed Forward DP matrix, ready for <h4_bckfilter()>.
 *
 * Throws:    <eslEMEM> if a DP matrix reallocation fails.
 *
 * Xref:      fbfilter.md for more notes
 */
int
(*h4_fwdfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc) =
  fwdfilter_dispatcher;

int
(*h4_bckfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh) =
  bckfilter_dispatcher;



/*****************************************************************
 * 2. CPU dispatching to vector implementations.
 *****************************************************************/

static int 
fwdfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc)
{
  /* SRE: for now */
  h4_fwdfilter = h4_fwdfilter_sse;
  return h4_fwdfilter_sse(dsq, L, hmm, mo, cpx, opt_sc);

  /* When a platform supports more than one vector implementation, 
   * put fastest one first to prefer enabling it over others.
   */
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512())
    {
      h4_fwdfilter = h4_fwdfilter_avx512;
      return h4_fwdfilter_avx512(dsq, L, hmm, mo, cpx, opt_sc);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      h4_fwdfilter = h4_fwdfilter_avx;
      return h4_fwdfilter_avx(dsq, L, hmm, mo, cpx, opt_sc);
    }
#endif

#ifdef eslENABLE_SSE4
  if (esl_cpu_has_sse4())
    {
      h4_fwdfilter = h4_fwdfilter_sse;
      return h4_fwdfilter_sse(dsq, L, hmm, mo, cpx, opt_sc);
    }
#endif
  
#ifdef eslENABLE_NEON
  h4_fwdfilter = h4_fwdfilter_neon;
  return h4_fwdfilter_neon(dsq, L, hmm, mo, cpx, opt_sc);
#endif

#ifdef eslENABLE_VMX
  h4_fwdfilter = h4_fwdfilter_vmx;
  return h4_fwdfilter_vmx(dsq, L, hmm, mo, cpx, opt_sc);
#endif

  esl_fatal("fwdfilter_dispatcher found no vector implementation - that shouldn't happen.");
}



static int
bckfilter_dispatcher(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh)
{
  /* SRE: for now */
  h4_bckfilter = h4_bckfilter_sse;
  return h4_bckfilter_sse(dsq, L, hmm, mo, cpx, sm, sm_thresh);

#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512())
    {
      h4_bckfilter = h4_bckfilter_avx512;
      return h4_bckfilter_avx512(dsq, L, hmm, mo, cpx, sm, sm_thresh);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      h4_bckfilter = h4_bckfilter_avx;
      return h4_bckfilter_avx(dsq, L, hmm, mo, cpx, sm, sm_thresh);
    }
#endif

#ifdef eslENABLE_SSE4
  if (esl_cpu_has_sse4())
    {
      h4_bckfilter = h4_bckfilter_sse;
      return h4_bckfilter_sse(dsq, L, hmm, mo, cpx, sm, sm_thresh);
    }
#endif
  
#ifdef eslENABLE_NEON
  h4_bckfilter = h4_bckfilter_neon;
  return h4_bckfilter_neon(dsq, L, hmm, mo, cpx, sm, sm_thresh);
#endif

#ifdef eslENABLE_VMX
  h4_bckfilter = h4_bckfilter_vmx;
  return h4_bckfilter_vmx(dsq, L, hmm, mo, cpx, sm, sm_thresh);
#endif

  esl_fatal("bckfilter_dispatcher found no vector implementation - that shouldn't happen.");
}




/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/
#ifdef h4FBFILTER_BENCHMARK
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
#include "fbfilter.h"

#define SIMDOPTS "--sse,--avx,--avx512"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range   toggles  reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,     NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,     NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,     NULL,  NULL, NULL, "only benchmark Forward",                           0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0",    NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "20000", NULL, "n>0",    NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "--sse",     eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the SSE4 implementation",              0 },
  { "--avx",     eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the AVX2 implementation",              0 },
  { "--avx512",  eslARG_NONE,    NULL, NULL,  NULL,SIMDOPTS,  NULL, NULL, "force using the AVX512 implementation",            0 },
  { "--version", eslARG_NONE,   FALSE, NULL,  NULL,    NULL,  NULL, NULL, "show HMMER version info",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go         = h4_CreateDefaultApp(options, 1, argc, argv, "benchmark driver for Viterbi filter", "[-options] <hmmfile>");
  char           *hmmfile    = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w          = esl_stopwatch_Create();
  ESL_RANDOMNESS *rng        = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc        = NULL;
  H4_HMMFILE     *hfp        = NULL;
  H4_PROFILE     *hmm        = NULL;
  H4_MODE        *mo         = h4_mode_Create();
  H4_CHECKPTMX   *cpx        = h4_checkptmx_Create(100,100,ESL_MBYTES(32));
  H4_SPARSEMASK  *sm         = h4_sparsemask_Create(100,100);
  int             L          = esl_opt_GetInteger(go, "-L");
  int             N          = esl_opt_GetInteger(go, "-N");
  int             do_fwdonly = esl_opt_GetBoolean(go, "-F");
  ESL_DSQ       **dsq        = malloc(N * sizeof(ESL_DSQ *));
  int             i;
  float           vfsc;
  double          Mcs;

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);

  h4_mode_SetLength(mo, L);

  /* Overriding the CPU dispatcher */
  if      (esl_opt_GetBoolean(go, "--sse"))    { h4_fwdfilter = h4_fwdfilter_sse;     /* h4_bckfilter = h4_bckfilter_sse;    */ }
  else if (esl_opt_GetBoolean(go, "--avx"))    { h4_fwdfilter = h4_fwdfilter_avx;     /* h4_bckfilter = h4_bckfilter_avx;    */ }
  else if (esl_opt_GetBoolean(go, "--avx512")) { h4_fwdfilter = h4_fwdfilter_avx512;  /* h4_bckfilter = h4_bckfilter_avx512; */ }

  /* generate and store the random test seqs */
  for (i = 0; i < N; i++)
    {
      dsq[i] = malloc(sizeof(ESL_DSQ) * (L+2));
      esl_rsq_xfIID(rng, hmm->f, abc->K, L, dsq[i]);
    }

  /* benchmark timing */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      h4_fwdfilter(dsq[i], L, hmm, mo, cpx, &vfsc);
      if (! do_fwdonly)
	{
	  h4_bckfilter(dsq[i], L, hmm, mo, cpx, sm, h4SPARSIFY_THRESH);
	  h4_sparsemask_Reuse(sm);
	}
    }
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
  h4_sparsemask_Destroy(sm);
  h4_checkptmx_Destroy(cpx);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif // h4FBFILTER_BENCHMARK



/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef h4FBFILTER_TESTDRIVE

#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sq.h"

#include "emit.h"
#include "logsum.h"
#include "modelsample.h"
#include "reference_dp.h"

static void
utest_compare_reference(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, int M, int L, int nseq, int64_t ramlimit)
{
  char           msg[]  = "fbfilter compare_reference unit test failed";
  ESL_DSQ       *dsq    = NULL;                             // will point either to sq->dsq or dsqmem for homologous vs. random seqs
  ESL_DSQ       *dsqmem = malloc(sizeof(ESL_DSQ) * (L+2));  // iid random seqs go here
  ESL_SQ        *sq     = esl_sq_CreateDigital(abc);        // homologous (emitted) seqs go here
  H4_PROFILE    *hmm    = NULL;
  H4_MODE       *mo     = h4_mode_Create();
  H4_REFMX      *fwd    = h4_refmx_Create(M,L);
  H4_REFMX      *bck    = h4_refmx_Create(M,L);
  H4_REFMX      *pp     = h4_refmx_Create(M,L);
  H4_CHECKPTMX  *cpx    = h4_checkptmx_Create(M,L,ramlimit);
  H4_SPARSEMASK *sm     = h4_sparsemask_Create(M,L);
  int            tL;
  float          fsc, fsc_ref, bsc_ref;
  float          ftol   = (h4_logsum_IsSlowExact() ?  0.001 : 0.1);
#if eslDEBUGLEVEL > 0
  float          pptol  = (h4_logsum_IsSlowExact() ? 0.0001 : 0.01);  /* posterior decoding values differ by no more than this */
#endif

  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg); 
  if ( h4_mode_SetLocal(mo)              != eslOK) esl_fatal(msg);  // FB filter is local-only. It ignores local/glocal setting of <mo>, but set <mo> for multilocal anyway.

  /* In debugging mode, set checkpoint mx to record full fwd/bck/pp matrices for
   * comparison to reference implementation.
   */
#if eslDEBUGLEVEL > 0
  cpx->fwd = h4_refmx_Create(M,L);
  cpx->bck = h4_refmx_Create(M,L);
  cpx->pp  = h4_refmx_Create(M,L);
#endif

  while (nseq--)
    {
      // test uses a 50:50 mix of emitted ("homologous") vs. iid random seqs
      if (esl_rnd_Roll(rng, 2))
	{
	  if ( esl_rsq_xfIID(rng, hmm->f, hmm->abc->K, L, dsqmem) != eslOK) esl_fatal(msg);
	  dsq = dsqmem;
	  tL  = L;
	}
      else
	{
	  if ( h4_mode_SetLength(mo, L)        != eslOK) esl_fatal(msg); // expected len, for emitting the target seq
	  do {
	    esl_sq_Reuse(sq);
	    if ( h4_emit(rng, hmm, mo, sq, NULL) != eslOK) esl_fatal(msg);
	  } while (sq->n > L*3);	// keep sequence length from getting ridiculous; also, long seqs do have higher expected absolute numerical error
	  dsq = sq->dsq;
	  tL  = sq->n;
	}

      if ( h4_mode_SetLength(mo, tL)                                  != eslOK) esl_fatal(msg);
      if ( h4_fwdfilter(dsq, tL, hmm, mo, cpx, &fsc)                  != eslOK) esl_fatal(msg);
      if ( h4_bckfilter(dsq, tL, hmm, mo, cpx, sm, h4SPARSIFY_THRESH) != eslOK) esl_fatal(msg);

      if ( h4_reference_Forward (dsq, tL, hmm, mo, fwd, &fsc_ref)     != eslOK) esl_fatal(msg);
      if ( h4_reference_Backward(dsq, tL, hmm, mo, bck, &bsc_ref)     != eslOK) esl_fatal(msg);
      if ( h4_reference_Decoding(dsq, tL, hmm, mo, fwd, bck, pp)      != eslOK) esl_fatal(msg);

      if ( fabs(fsc - fsc_ref) > ftol) esl_fatal(msg);

#if eslDEBUGLEVEL > 0
      if ( fabs(cpx->bcksc - bsc_ref) > ftol) esl_fatal(msg);

      if ( h4_refmx_Compare     (pp,  cpx->pp,  pptol) != eslOK) esl_fatal(msg);
      if ( h4_refmx_Compare     (fwd, cpx->fwd, ftol)  != eslOK) esl_fatal(msg);
      if ( h4_refmx_CompareLocal(bck, cpx->bck, ftol)  != eslOK) esl_fatal(msg);
#endif      
			
      h4_sparsemask_Reuse(sm);
      h4_refmx_Reuse(fwd);
      h4_refmx_Reuse(bck);
      h4_refmx_Reuse(pp);
    }
  
  free(dsqmem);
  esl_sq_Destroy(sq);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  h4_refmx_Destroy(fwd);
  h4_refmx_Destroy(bck);
  h4_refmx_Destroy(pp);
  h4_sparsemask_Destroy(sm);
  h4_checkptmx_Destroy(cpx);
}
#endif // h4FBFILTER_TESTDRIVE

/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef h4FBFILTER_TESTDRIVE

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
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random model to sample",                 0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version info",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = h4_CreateDefaultApp(options, 0, argc, argv, "fbfilter test driver", "[-options]");
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");

  utest_compare_reference(rng, abc, M, L, N,  ESL_MBYTES(32)); // normal use
  utest_compare_reference(rng, abc, M, L, N,  ESL_MBYTES(0));  // zero memory: force checkpointing
  utest_compare_reference(rng, abc, 1, L, 10, ESL_MBYTES(32)); // size 1 profiles
  utest_compare_reference(rng, abc, M, 1, 10, ESL_MBYTES(32)); // size 1 seqs

  esl_alphabet_Destroy(abc);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");

  utest_compare_reference(rng, abc, M, L, N,  ESL_MBYTES(32)); 
  utest_compare_reference(rng, abc, M, L, N,  ESL_MBYTES(0));  
  utest_compare_reference(rng, abc, 1, L, 10, ESL_MBYTES(32)); 
  utest_compare_reference(rng, abc, M, 1, 10, ESL_MBYTES(32)); 

  esl_alphabet_Destroy(abc);

  fprintf(stderr, "#  status   = ok\n");

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(rng);
  return eslOK;
}
#endif // h4FBFILTER_TESTDRIVE


/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef h4FBFILTER_EXAMPLE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_hmmfile.h"
#include "h4_profile.h"
#include "h4_mode.h"
#include "h4_refmx.h"
#include "h4_checkptmx.h"

#include "general.h"
#include "fbfilter.h"
#include "reference_dp.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version info",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv, "example of using the Fwd/Bck filter", "[-options] <hmmfile> <seqfile>");
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  H4_REFMX       *fwd     = NULL;
  H4_CHECKPTMX   *cpx     = NULL;
  H4_SPARSEMASK  *sm      = NULL;
  float           fsc, bsc, fsc_ref;
  int             status;

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);

  if ( esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, /*env=*/NULL, &sqfp) != eslOK)
    esl_fatal("couldn't open sequence file %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  cpx = h4_checkptmx_Create (hmm->M, 400, ESL_MBYTES(32));
  sm  = h4_sparsemask_Create(hmm->M, 400);
  fwd = h4_refmx_Create     (hmm->M, 400);

  h4_mode_SetLocal(mo);  // FB filter is implicitly local. We set local mode so reference Fwd scores match FB.
  
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      h4_mode_SetLength(mo, sq->n);

#if eslDEBUGLEVEL > 0
      h4_checkptmx_SetDumpMode(cpx, stdout);
#endif

      h4_fwdfilter(sq->dsq, sq->n, hmm, mo, cpx, &fsc);
      h4_bckfilter(sq->dsq, sq->n, hmm, mo, cpx, sm, h4SPARSIFY_THRESH);
      h4_reference_Forward(sq->dsq, sq->n, hmm, mo, fwd, &fsc_ref);

      bsc = 0.0;
#if eslDEBUGLEVEL > 0
      bsc = cpx->bcksc;   // backward score only available in debugging mode
#endif      

      //h4_refmx_Dump(stdout, fwd);

      printf("fwdfilter raw score (bits): %.2f\n", fsc);
      printf("bckfilter raw score (bits): %.2f\n", bsc);
      printf("reference Fwd score (bits): %.2f\n", fsc_ref);

      h4_refmx_Reuse(fwd);
      h4_sparsemask_Reuse(sm);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);


  h4_sparsemask_Destroy(sm);
  h4_checkptmx_Destroy(cpx);
  h4_refmx_Destroy(fwd);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif // h4FBFILTER_EXAMPLE
