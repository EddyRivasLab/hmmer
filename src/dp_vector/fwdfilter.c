/* Forwards/Backwards filters
 * 
 * See fwdfilter.md for notes.
 *
 * Contents:
 *    1. ForwardFilter() and BackwardFilter() API.
 *    2. CPU dispatching to vector implementations.
 *    3. Stats driver (memory requirements)
 *    4. Benchmark driver.
 *    5. Unit tests.
 *    6. Test driver.
 *    7. Example.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_sparse/p7_sparsemx.h"
#include "dp_vector/fwdfilter.h"


static int fwdfilter_dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc);
static int bckfilter_dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh);



/*****************************************************************
 * 1. ForwardFilter, BackwardFilter API calls
 *****************************************************************/

/* Function:  p7_ForwardFilter()
 * Synopsis:  Checkpointed striped vector Forward calculation, producing score.
 *
 * Purpose:   Calculate the Forward algorithm for target sequence <dsq>
 *            of <L> residues aligned to query profile <om>, using the
 *            checkpointed DP matrix <ox> provided by the caller. Upon
 *            successful return, <ox> contains the filled Forward
 *            matrix, and <*opt_sc> optionally contains the raw Forward
 *            score in nats.
 *            
 *            <ox> will be reallocated, if needed, for the MxL problem;
 *            caller does not need to call <p7_checkptmx_Reinit()> itself.
 *  
 * Args:      dsq    - digital target sequence, 1..L
 *            L      - length of dsq, residues
 *            om     - optimized profile (multihit local)
 *            ox     - checkpointed DP matrix
 *            opt_sc - optRETURN: raw Forward score (nats)
 *
 * Returns:   <eslOK> on success, <ox> contains the checkpointed
 *            Forward matrix calculation, and <*opt_sc> optionally
 *            has the raw Forward score in nats.
 *
 * Throws:    <eslEMEM> on reallocation error.
 * 
 * Xref:      For layout of checkpointed <ox> see exegesis in p7_checkptmx.h.
 */
int 
(*p7_ForwardFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc) =
  fwdfilter_dispatcher;


/* Function:  p7_BackwardFilter()
 * Synopsis:  Checkpointed striped vector Backward calculation, producing sparse mask.
 *
 * Purpose:   Given a target sequence <dsq> of length <L>, a query model
 *            <om>, and a DP matrix <ox> resulting from a successful
 *            call to <p7_ForwardFilter()>; calculate the Backward and
 *            posterior decoding algorithms. On each row <i=1..L>, use
 *            posterior decoding to determine which <i,k> cells pass a
 *            significance threshold for inclusion in the sparse DP
 *            mask. Store those sparse cells in <sm>, which the caller
 *            allocates (or reuses) and provides. <sm> will be resized
 *            here as needed; caller does not need to call 
 *            <p7_sparsemask_Reinit()> on it.
 *            
 *            
 * Args:      dsq       - digital target sequence, 1..L
 *            L         - length of dsq, residues
 *            om        - optimized profile (multihit local)
 *            ox        - checkpointed DP matrix, ForwardFilter already run
 *            sm        - allocated P7_SPARSEMASK structure to hold sparse DP mask
 *            sm_thresh - Threshold for determining 'significant' posterior
 *                        alignment probability, passed in to posterior_decode_row()
 *
 * Throws:    <eslEINVAL> if something's awry with a data structure's internals.
 *            <eslEMEM> on allocation failure.
 */
int
(*p7_BackwardFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh) =
  bckfilter_dispatcher;

/*----------- end forward/backward API calls --------------------*/



/*****************************************************************
 * 2. CPU dispatching to vector implementations.
 *****************************************************************/

static int 
fwdfilter_dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc)
{
  /* When a platform supports more than one vector implementation, 
   * put fastest one first to prefer enabling it over others.
   */
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512())
    {
      p7_ForwardFilter = p7_ForwardFilter_sse;
      return p7_ForwardFilter_sse(dsq, L, om, ox, opt_sc);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      p7_ForwardFilter = p7_ForwardFilter_avx;
      return p7_ForwardFilter_avx(dsq, L, om, ox, opt_sc);
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse())
    {
      p7_ForwardFilter = p7_ForwardFilter_sse;
      return p7_ForwardFilter_sse(dsq, L, om, ox, opt_sc);
    }
#endif
  
#ifdef eslENABLE_NEON
  p7_ForwardFilter = p7_ForwardFilter_neon;
  return p7_ForwardFilter_neon(dsq, L, om, ox, opt_sc);
#endif

  //#ifdef eslENABLE_VMX
  //  p7_ForwardFilter = p7_ForwardFilter_vmx;
  //  return p7_ForwardFilter_vmx(dsq, L, om, ox, opt_sc);
  //#endif

  p7_Die("fwdfilter_dispatcher found no vector implementation - that shouldn't happen.");
}




static int
bckfilter_dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512())
    {
      p7_BackwardFilter = p7_BackwardFilter_sse;
      return p7_BackwardFilter_sse(dsq, L, om, ox, sm, sm_thresh);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      p7_BackwardFilter = p7_BackwardFilter_avx;
      return p7_BackwardFilter_avx(dsq, L, om, ox, sm, sm_thresh);
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse())
    {
      p7_BackwardFilter = p7_BackwardFilter_sse;
      return p7_BackwardFilter_sse(dsq, L, om, ox, sm, sm_thresh);
    }
#endif
  
#ifdef eslENABLE_NEON
  p7_BackwardFilter = p7_BackwardFilter_neon;
  return p7_BackwardFilter_neon(dsq, L, om, ox, sm, sm_thresh);
#endif

  //#ifdef eslENABLE_VMX
  //  p7_BackwardFilter = p7_BackwardFilter_vmx;
  //  return p7_BackwardFilter_vmx(dsq, L, om, ox, sm, sm_thresh);
  //#endif

  p7_Die("bckfilter_dispatcher found no vector implementation - that shouldn't happen.");
}







/*****************************************************************
 * 3. Stats driver: memory requirement profiling 
 *****************************************************************/

/* In production pipeline, there are up to four sorts of DP matrices
 * in play:
 *    MSV   MSV filter
 *    VF    Viterbi filter
 *    CHK   checkpointed Forward/Backward
 *    SP    sparse DP
 * and the reference implementation has one, for our baseline:
 *    REF   reference DP implementation of everything
 * and, the sparse DP implementation also requires a sparse mask,
 * which is of on the same memory usage order as the DP matrix:
 *    SM    sparse DP mask
 * We count it separately here because we will typically have one
 * mask, for up to four (V/F/B/D) matrices in play.
 *    
 * MSV and VF are O(M).
 * CHK is O(M \sqrt L)
 * SP and SM are O(L) when a nonzero sparsity posterior threshold is used in the f/b filter.
 * REF is O(ML). 
 * 
 * The goal of this driver is to characterize the minimal memory
 * requirements for each, and in particular, for sparse DP.  MSV, VF,
 * CHK, and REF minimal requirements are completely determined by
 * dimensions M and L. The requirement of SP, though, is an empirical
 * question of how many sparse cells get included by the f/b filter.
 * This driver runs the f/b filter, using one query profile, against a
 * target seq db, and for each profile/seq comparison, it prints a
 * line summarizing all the minimal memory requirements.
 * 
 * Note that even though we're empirically characterizing the typical
 * minimum requirement of SP, it nonetheless does have a proven O(L)
 * memory complexity for a nonzero posterior prob threshold in f/b
 * filter.
 * 
 * Note also that this is characterizing the "minimal" requirement,
 * not the actual allocated size of the various matrices, which tend
 * to grow depending on previous requirements. We're not testing our
 * reallocation strategy here. We might use these statistics to help
 * refine that strategy.
 */
#ifdef p7FWDFILTER_STATS
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "stats driver, ForwardFilter()";

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
  P7_FILTERMX    *fx      = NULL;
  P7_CHECKPTMX   *ox      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fraw, nullsc, fsc, msvsc;
  float           msvmem, vfmem, fbmem, smmem, spmem, refmem;
  double          P;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  /* Open sequence file for reading */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);               
  gm = p7_profile_Create(hmm->M, abc); 
  p7_profile_Config(gm, hmm, bg);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  
  /* Initially allocate matrices for a default sequence size, 500; 
   * we resize as needed for each individual target seq 
   */
  fx  = p7_filtermx_Create  (gm->M);
  ox  = p7_checkptmx_Create (gm->M, 500, ESL_MBYTES(32));  
  sm  = p7_sparsemask_Create(gm->M, 500);

  printf("# %-28s %-30s %9s %9s %9s %9s %9s %9s %9s %9s\n",
	 "target seq", "query profile", "score", "P-value", "MSV (KB)", "VF (KB)", "CHK (MB)", "SM (MB)", "SP (MB)", "REF (MB)");
  printf("# %-28s %-30s %9s %9s %9s %9s %9s %9s %9s %9s\n",
	 "----------------------------", "------------------------------",  "---------", "---------", "---------", "---------", "---------", "---------", "---------", "---------");
  
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_bg_SetLength(bg,            sq->n);

      p7_checkptmx_GrowTo (ox, om->M, sq->n); 

      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);

      /* Filter insig hits, partially simulating the real pipeline */
      p7_MSVFilter(sq->dsq, sq->n, om, fx, &msvsc);
      msvsc = (msvsc - nullsc) / eslCONST_LOG2;
      P     =  esl_gumbel_surv(msvsc,  om->evparam[p7_SMU],  om->evparam[p7_MLAMBDA]);
      if (P > 0.02) goto NEXT_SEQ;

      p7_ForwardFilter (sq->dsq, sq->n, om, ox, &fraw);
      p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH);

      /* Calculate minimum memory requirements for each step */
      msvmem = (double) ( P7_Q(om->M, om->V) * sizeof(__m128i)) / 1024.;  
      vfmem  = (double) p7_filtermx_MinSizeof(om->M)            / 1024.;
      fbmem  = (double) p7_checkptmx_MinSizeof(om->M, sq->n)    / 1024. / 1024.;
      smmem  = (double) p7_sparsemask_MinSizeof(sm)             / 1024. / 1024.;
      spmem  = (double) p7_sparsemx_MinSizeof(sm)               / 1024. / 1024.;
      refmem = (double) p7_refmx_MinSizeof(om->M, sq->n)        / 1024. / 1024.;

      fsc  =  (fraw-nullsc) / eslCONST_LOG2;
      P    = esl_exp_surv(fsc,   om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);

      printf("%-30s %-30s %9.2f %9.2g %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	     sq->name,
	     hmm->name,
	     fsc, P,
	     msvmem,
	     vfmem,
	     fbmem, 
	     smmem,
	     spmem,
	     refmem);
	     
    NEXT_SEQ:
      esl_sq_Reuse(sq);
      p7_sparsemask_Reuse(sm);
      //p7_checkptmx_Reuse(ox);
      //p7_filtermx_Reuse(fx);
    }

  printf("# SPARSEMASK: kmem reallocs: %d\n", sm->n_krealloc);
  printf("#             seg reallocs:  %d\n", sm->n_srealloc);
  printf("#             row reallocs:  %d\n", sm->n_rrealloc);

  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(ox);
  p7_filtermx_Destroy(fx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}


#endif
/*--------------- end, stats driver -----------------------------*/


/*****************************************************************
 * 4. Benchmark
 *****************************************************************/

#ifdef p7FWDFILTER_BENCHMARK
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Forward",                         0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,   "2000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for checkpointed ForwardFilter()";

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
  P7_CHECKPTMX   *ox      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, L);

  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  ox  = p7_checkptmx_Create (om->M, L, ESL_MBYTES(32));
  sm  = p7_sparsemask_Create(om->M, L);

  /* Baseline time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_ForwardFilter(dsq, L, om, ox, &sc);
      if (! esl_opt_GetBoolean(go, "-F")) 
	{
	  p7_BackwardFilter(dsq, L, om, ox, sm, p7_SPARSIFY_THRESH);
	  esl_vec_IReverse(sm->kmem, sm->kmem, sm->ncells);
	}

      //p7_checkptmx_Reuse(ox);
      p7_sparsemask_Reuse(sm);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_checkptmx_Destroy(ox);
  p7_sparsemask_Destroy(sm);
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
#endif /*p7FWDFILTER_BENCHMARK*/
/*-------------------- end, benchmark ---------------------------*/


/*****************************************************************
 * 6. Unit tests
 *****************************************************************/
#ifdef p7FWDFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sqio.h"

#include "search/modelconfig.h"
#include "misc/logsum.h"
#include "misc/emit.h"

#include "dp_reference/reference_fwdback.h"
#include "dp_reference/reference_decoding.h"

/* Compare scores of Forward, Backward to those from the reference
 * implementation.
 * 
 */
static void
utest_scores(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int64_t ramlimit)
{
  char           msg[]  = "fwdfilter scores unit test failed";
  P7_HMM        *hmm    = NULL;
  P7_PROFILE    *gm     = NULL;
  P7_OPROFILE   *om     = NULL;
  ESL_DSQ       *dsqmem = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_DSQ       *dsq    = NULL;
  int            tL     = 0;
  ESL_SQ        *sq     = esl_sq_CreateDigital(abc);
  P7_CHECKPTMX  *ox     = p7_checkptmx_Create(M, L, ramlimit);
  P7_REFMX      *fwd    = p7_refmx_Create   (M, L);
  P7_REFMX      *bck    = p7_refmx_Create   (M, L);
  P7_REFMX      *pp     = p7_refmx_Create   (M, L);
  P7_SPARSEMASK *sm     = p7_sparsemask_Create(M, L);
  float          tol2   = ( p7_logsum_IsSlowExact() ? 0.001  : 0.1);   /* absolute agreement of reference (log-space) and vector (prob-space) depends on whether we're using LUT-based logsum() */
  float          fsc1, fsc2;
  float          bsc2;
#if eslDEBUGLEVEL > 0
  float          bsc1;
  float          tol1   = 0.0001;	                               /* forward and backward scores from same implementation type should agree with high tolerance */
  float          ptol   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01);  /* posterior decoding values differ by no more than this */
#endif

#if eslDEBUGLEVEL > 0
  /* We set the debugging tools to record full pp, fwd, bck matrices
   * for comparison to reference implementation: 
   */
  ox->pp  = p7_refmx_Create(M,L);	
  ox->fwd = p7_refmx_Create(M,L);	
  ox->bck = p7_refmx_Create(M,L);	
#endif

  if ( p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om) != eslOK) esl_fatal(msg);
  /* note the <gm> is config'ed local-only; both <om>,<gm> have length model set to <L> */

  while (N--)
    {
      p7_profile_SetLength(gm, L);       /* because it may have been reset to tL by last emitted seq */
      p7_oprofile_ReconfigLength(om, L); 

      /* A mix of generated (homologous) and random (nonhomologous) sequences */
      if (esl_rnd_Roll(r, 2)) 
	{
	  esl_rsq_xfIID(r, bg->f, abc->K, L, dsqmem);  
	  dsq = dsqmem;  
	  tL = L;     
	  /* fixed-length random emission: length config on <gm>,<om> don't need to change  */
	}
      else  
	{
	  do {
	    esl_sq_Reuse(sq);
	    p7_ProfileEmit(r, hmm, gm, bg, sq, NULL);
	  } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
	  dsq = sq->dsq; 
	  tL = sq->n; 
	  /* variable-length seq emission: length config on <gm>,<om> will change */
	}
	
      if ( p7_profile_SetLength(gm, tL)       != eslOK) esl_fatal(msg);
      if ( p7_oprofile_ReconfigLength(om, tL) != eslOK) esl_fatal(msg); 
      if ( p7_checkptmx_Reinit(ox,  M, tL)    != eslOK) esl_fatal(msg);
      
      p7_ForwardFilter (dsq, tL, om, ox, &fsc1);
      p7_BackwardFilter(dsq, tL, om, ox,  sm, p7_SPARSIFY_THRESH);

      p7_ReferenceForward (dsq, tL, gm, fwd,  &fsc2);
      p7_ReferenceBackward(dsq, tL, gm, bck,  &bsc2);
      p7_ReferenceDecoding(dsq, tL, gm, fwd, bck, pp);

#if eslDEBUGLEVEL > 0
      /* vector Forward and Backward scores should agree with high tolerance.
       * Backward score is only available in debugging mode 
       */
      bsc1 = ox->bcksc;
      if (fabs(fsc1-bsc1) > tol1)    esl_fatal(msg);
#endif

      /* reference scores should agree with tolerance depending on whether logsum compiled to use LUT or not */
      if (fabs(fsc2-bsc2) > tol2)    esl_fatal(msg);

      /* Reference and vector implementations should agree depending on logsum */
      if (fabs(fsc1-fsc2) > tol2) esl_fatal(msg);

#if eslDEBUGLEVEL > 0
      /* Compare all DP cell values to reference implementation,
       * in fwd, bck, and pp matrices. Note the need for CompareLocal()
       * for the Backward matrix comparison, because the zero B->G
       * transition isn't evaluated until the *end* of glocal paths;
       * thus Backward values along the glocal paths are finite for 
       * the reference implementation, whereas in the ForwardFilter
       * they're -inf by construction (ForwardFilter is local-only).
       */
      if (p7_refmx_Compare     (pp,  ox->pp,  ptol) != eslOK) esl_fatal(msg);
      if (p7_refmx_Compare     (fwd, ox->fwd, tol2) != eslOK) esl_fatal(msg);
      if (p7_refmx_CompareLocal(bck, ox->bck, tol2) != eslOK) esl_fatal(msg);
#endif
      esl_sq_Reuse(sq);
      p7_refmx_Reuse(fwd);
      p7_refmx_Reuse(bck);
      p7_refmx_Reuse(pp);
      //p7_checkptmx_Reuse(ox);
      p7_sparsemask_Reuse(sm);
    }

  free(dsqmem);
  esl_sq_Destroy(sq);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(ox);
  p7_refmx_Destroy(fwd);
  p7_refmx_Destroy(bck);
  p7_refmx_Destroy(pp);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7FWDFILTER_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/


/*****************************************************************
 * 7. Test driver
 *****************************************************************/
#ifdef p7FWDFILTER_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "hardware/hardware.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for checkpointed vector SSE Forward, Backward implementations";

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

  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(32));   /* normal sized models              */
  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(0));    /* zero memory: force checkpointing */
  utest_scores(r, abc, bg, 1, L, 10, ESL_MBYTES(32));   /* size 1 models                    */
  utest_scores(r, abc, bg, M, 1, 10, ESL_MBYTES(32));   /* size 1 sequences                 */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(32));  
  utest_scores(r, abc, bg, M, L, N,  ESL_MBYTES(0));   
  utest_scores(r, abc, bg, 1, L, 10, ESL_MBYTES(32));  
  utest_scores(r, abc, bg, M, 1, 10, ESL_MBYTES(32));  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);

  fprintf(stderr, "#  status = ok\n");
  return eslOK;
}
#endif /*p7FWDFILTER_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


/*****************************************************************
 * 8. Example
 *****************************************************************/
#ifdef p7FWDFILTER_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "output in one line awkable format",                0 },
#if eslDEBUGLEVEL > 0
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump vector DP matrices for examination (verbose)",0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump recorded forward matrix",                     0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump recorded backward matrix",                    0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump recorded posterior prob matrix",              0 },
  { "--diplot",  eslARG_OUTFILE, NULL, NULL, NULL,   NULL,  NULL, NULL, "save domain inference plot to <f>",                0 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example driver, ForwardFilter()";

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
  P7_REFMX       *gx      = NULL;
  P7_CHECKPTMX   *ox      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fraw, nullsc, fsc, bsc;
  float           gfraw, gbraw, gfsc;
  float           gmem, cmem, bmem;
  double          P, gP;
  int             status;
#if eslDEBUGLEVEL > 0
  int             store_pp   = FALSE;
  char           *diplotfile = esl_opt_GetString(go, "--diplot");
  FILE           *difp       = NULL;
#endif

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  /* Open sequence file for reading */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

#if eslDEBUGLEVEL > 0
  if (diplotfile && (difp = fopen(diplotfile, "w")) == NULL) p7_Fail("couldn't open %s for writing", diplotfile);
  if (difp || esl_opt_GetBoolean(go, "-P")) store_pp = TRUE;
#endif

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);               
  gm = p7_profile_Create(hmm->M, abc); 
  p7_profile_Config(gm, hmm, bg);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  /* p7_oprofile_Dump(stdout, om);  */

  /* Initially allocate matrices for a default sequence size, 500; 
   * we resize as needed for each individual target seq 
   */
  ox  = p7_checkptmx_Create (gm->M, 500, ESL_MBYTES(32));  
  gx  = p7_refmx_Create     (gm->M, 500);
  sm  = p7_sparsemask_Create(gm->M, 500);
#if eslDEBUGLEVEL > 0
  /* When the eslDEBUGLEVEL is nonzero, <ox> matrix has the ability to
   * record generic, complete <fwd>, <bck>, and <pp> matrices, for
   * comparison to reference implementation, even when checkpointing.
   */
  if (esl_opt_GetBoolean(go, "-D")) p7_checkptmx_SetDumpMode(ox, stdout, TRUE);
  if (esl_opt_GetBoolean(go, "-F")) ox->fwd = p7_refmx_Create(gm->M, 100);
  if (esl_opt_GetBoolean(go, "-B")) ox->bck = p7_refmx_Create(gm->M, 100);
  if (store_pp)                     ox->pp  = p7_refmx_Create(gm->M, 100);
#endif

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_bg_SetLength(bg,            sq->n);
      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);

      p7_checkptmx_GrowTo  (ox, om->M, sq->n); 
    
      p7_ForwardFilter (sq->dsq, sq->n, om, ox, &fraw);
      p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH);

      p7_ReferenceForward (sq->dsq, sq->n, gm, gx, &gfraw);
      p7_ReferenceBackward(sq->dsq, sq->n, gm, gx, &gbraw);

      bsc = 0.0;		/* Backward score only available in debugging mode */
#if eslDEBUGLEVEL > 0
      if (esl_opt_GetBoolean(go, "-F")) p7_refmx_Dump(stdout, ox->fwd);
      if (esl_opt_GetBoolean(go, "-B")) p7_refmx_Dump(stdout, ox->bck);
      if (esl_opt_GetBoolean(go, "-P")) p7_refmx_Dump(stdout, ox->pp);
      if (difp)                         p7_refmx_PlotDomainInference(difp, ox->pp, 1, sq->n, NULL);
      bsc  =  (ox->bcksc-nullsc) / eslCONST_LOG2;
#endif

      fsc  =  (fraw-nullsc) / eslCONST_LOG2;
      gfsc = (gfraw-nullsc) / eslCONST_LOG2;
      P  = esl_exp_surv(fsc,   om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
      gP = esl_exp_surv(gfsc,  gm->evparam[p7_FTAU],  gm->evparam[p7_FLAMBDA]);

      gmem = (float) p7_refmx_Sizeof(gx)     / 1024 / 1024;
      cmem = (float) p7_checkptmx_Sizeof(ox) / 1024 / 1024;
      bmem = (float) sm->ncells * 6. * sizeof(float) / 1024 / 1024;

      if (esl_opt_GetBoolean(go, "-1")) 
	printf("%-30s\t%-20s\t%9.2g\t%7.4f\t%7.4f\t%9.2g\t%6.1f\t%6.2fM\t%6.2fM\t%6.2fM\n", sq->name, hmm->name, P, fsc, bsc, gP, gfsc, gmem, cmem, bmem);
      else
	{
	  printf("query model:               %s\n",        hmm->name);
	  printf("target sequence:           %s\n",        sq->name);
	  printf("fwd filter raw score:      %.4f nats\n", fraw);
#if eslDEBUGLEVEL > 0
	  printf("bck filter raw score:      %.4f nats\n", ox->bcksc);
#endif
	  printf("null score:                %.2f nats\n", nullsc);
	  printf("per-seq score:             %.2f bits\n", fsc);
	  printf("P-value:                   %g\n",        P);
	  printf("Reference fwd raw score:   %.2f nats\n", gfraw);
	  printf("Reference bck raw score:   %.2f nats\n", gbraw);
	  printf("Reference Fwd bit score:   %.2f bits\n", gfsc);
	  printf("Reference Forward P-val:   %g\n",        gP);
	  printf("RAM, f/b filter:           %.2fM\n",    cmem);
	  printf("RAM, generic:              %.2fM\n",    gmem);
	  printf("RAM, sparse:               %.2fM\n",    bmem);
	}

      esl_sq_Reuse(sq);
      p7_refmx_Reuse(gx);
      p7_sparsemask_Reuse(sm);
      //p7_checkptmx_Reuse(ox);
    }

#if eslDEBUGLEVEL > 0
  if (difp) fclose(difp);
#endif
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(ox);
  p7_refmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7FWDFILTER_EXAMPLE*/
/*---------------------- end, example ---------------------------*/



                                          
