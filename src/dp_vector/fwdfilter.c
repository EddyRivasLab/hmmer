/* Forwards/Backwards filters
 *
 * This is code in the acceleration pipeline:
 * SSVFilter -> MSVFilter -> VitFilter -> ForwardFilter -> BackwardFilter
 *                                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 *                                                (you are here)
 *
 * Striped SIMD vector implementation [Farrar07], using checkpointing 
 * to guarantee O(M sqrt L) memory [Grice97,TarnasHughey98,Newberg08],
 * in probability space using sparse rescaling [Eddy11].
 *
 * Called in two pieces. p7_ForwardFilter() returns the Forward score
 * in nats. Caller chooses whether or not to proceed with Backward and
 * posterior decoding. If caller proceeds, p7_BackwardFilter() does
 * Backward and posterior decoding. Based on the posterior decoding
 * probabilities on each row i, it determines which cells are to be
 * added to the sparse DP mask for subsequent local/glocal
 * reprocessing.  The sparse DP mask is returned in a <P7_SPARSEMASK>
 * structure.
 * 
 * Any cell (i,k) with total posterior probability (i.e., summed over
 * M,D,I) >= sm_thresh is marked and included in the sparse mask.
 * Cells needed for glocal entry/exit delete paths (G->DDDD->Mk,
 * Mk->DDDD->E) are not marked, because the filter only uses local
 * alignment, not glocal. The default for sm_thresh is
 * <p7_SPARSIFY_THRESH>, in <p7_config.h.in>; currently set
 * to 0.01.
 *
 * ForwardFilter() and BackwardFilter() are a dependent pair, sharing
 * the same DP matrix object. They must be called sequentially,
 * because BackwardFilter() is doing posterior decoding on the fly,
 * and for this it needs both forward and backward scores; it uses
 * checkpointed Forward rows that have been stored by the preceding
 * ForwardFilter() call.
 * 
 * Prob-space implementations using sparse rescaling require multihit
 * local alignment mode for numeric range reasons. Unihit or glocal
 * will result in errors due to underflow.
 *
 * Contents:
 *    1. ForwardFilter() and BackwardFilter() API.
 *    2. Internal functions: inlined recursions.
 *    3. Debugging tools.
 *    4. Stats driver (memory requirements)
 *    5. Benchmark driver.
 *    6. Unit tests.
 *    7. Test driver.
 *    8. Example.
 *    9. Notes
 *       a. On debugging and testing methods.
 *       b. On running time, in theory and in practice.
 *       c. Copyright and license information.
 */

#include "p7_config.h"
#if p7_CPU_ARCH == intel 
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#endif
#include "easel.h"
#include "esl_vectorops.h"

#include "dp_reference/p7_refmx.h"
#include "dp_sparse/p7_sparsemx.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"


/* vectorized DP recursions are tediously lengthy, so for some
 * semblance of clarity, they're broken out into one-page-ish
 * chunks, using static inlined functions.
 */

#ifdef p7_DEBUGGING
static inline float backward_row_zero(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
extern inline float backward_row_zero_sse(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
extern inline float backward_row_zero_avx(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
extern inline float backward_row_zero_avx512(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
extern inline float backward_row_zero_neon(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
extern inline float backward_row_zero_neon64(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox);
static        void  save_debug_row_pp(P7_CHECKPTMX *ox,               debug_print *dpc, int i);
extern void  save_debug_row_pp_sse(P7_CHECKPTMX *ox,               debug_print *dpc, int i);
extern void  save_debug_row_pp_avx(P7_CHECKPTMX *ox,               debug_print *dpc, int i);
extern void  save_debug_row_pp_avx512(P7_CHECKPTMX *ox,               debug_print *dpc, int i);
extern void  save_debug_row_pp_neon(P7_CHECKPTMX *ox,               debug_print *dpc, int i);
extern void  save_debug_row_pp_neon64(P7_CHECKPTMX *ox,               debug_print *dpc, int i);
static void  save_debug_row_fb(P7_CHECKPTMX *ox, P7_REFMX *gx, debug_print *dpc, int i, float totscale);
extern void  save_debug_row_fb_sse(P7_CHECKPTMX *ox, P7_REFMX *gx, debug_print *dpc, int i, float totscale);
extern void  save_debug_row_fb_avx(P7_CHECKPTMX *ox, P7_REFMX *gx, debug_print *dpc, int i, float totscale);
extern void  save_debug_row_fb_avx512(P7_CHECKPTMX *ox, P7_REFMX *gx, debug_print *dpc, int i, float totscale);
extern void  save_debug_row_fb_neon(P7_CHECKPTMX *ox, P7_REFMX *gx, debug_print *dpc, int i, float totscale);
extern void  save_debug_row_fb_neon64(P7_CHECKPTMX *ox, P7_REFMX *gx, debug_print *dpc, int i, float totscale);

#endif

/*****************************************************************
 * 1. Forward and Backward API calls
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
 *            caller does not need to call <p7_checkptmx_GrowTo()> itself.
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

 /* Note: This is just a dispatch function that calls the correct SIMD version of the filter
  * Whenever possible, you should call the SIMD version directly to avoid the dispatch
     overhead.  */

int
p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, float *opt_sc)
{
switch(om->simd){
    case SSE:
      return p7_ForwardFilter_sse(dsq, L, om, ox, opt_sc);
      break;
    case AVX:
      return p7_ForwardFilter_avx(dsq, L, om, ox, opt_sc);
      break;
    case AVX512:
      return p7_ForwardFilter_avx512(dsq, L, om, ox, opt_sc);
      break;
    case NEON:
      return p7_ForwardFilter_neon(dsq, L, om, ox, opt_sc);
      break;
    case NEON64:
      return p7_ForwardFilter_neon64(dsq, L, om, ox, opt_sc);
      break;
   default:
      p7_Fail("Unrecognized SIMD type passed to p7_ForwardFilter");  
  }

}


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
p7_BackwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHECKPTMX *ox, P7_SPARSEMASK *sm, float sm_thresh)
{
switch(om->simd){
    case SSE:
      return p7_BackwardFilter_sse(dsq, L, om, ox, sm, sm_thresh);
      break;
    case AVX:
      return p7_BackwardFilter_avx(dsq, L, om, ox, sm, sm_thresh);
      break;
    case AVX512:
      return p7_BackwardFilter_avx512(dsq, L, om, ox, sm, sm_thresh);
      break;
    case NEON:
      return p7_BackwardFilter_neon(dsq, L, om, ox, sm, sm_thresh);
      break;
    case NEON64:
      return p7_BackwardFilter_neon64(dsq, L, om, ox, sm, sm_thresh);
      break;
   default:
      p7_Fail("Unrecognized SIMD type passed to p7_BackwardFilter");  
  }
   
}
/*----------- end forward/backward API calls --------------------*/



/*****************************************************************
 * 2. Internal functions: inlined recursions
 *****************************************************************/



/*****************************************************************
 * 3. Debugging tools
 *****************************************************************/
#ifdef p7_DEBUGGING

/* backward_row_zero()
 * 
 * Slightly peculiar but true: in production code we don't
 * need to calculate backward row 0, because we're only doing
 * posterior decoding on residues 1..L.  We only need row 0
 * if we need a complete Backward calculation -- which happens
 * when we're in debugging mode, and we're going to test the
 * Backwards score or compare the Backwards matrix to the
 * reference implementation.
 * 
 * x1  - residue dsq[1] on row 1
 * om  - query model
 * ox  - DP matrix; dpf[0] is the current Backward row
 * 
 * Upon return, ox->dpf[0] (backward row 0) has been calculated.
 * Returns the log of the value in the N cell at row 0. When this 
 * is added to the sum of the logs of all the scalefactors (which
 * the caller has been accumulating in <ox->bcksc>), then
 * the <ox->bcksc> value is finished and equal to the Backwards
 * raw score in nats.
 */
static inline float
backward_row_zero(ESL_DSQ x1, const P7_OPROFILE *om, P7_CHECKPTMX *ox)
{
 switch(om->simd){
    case SSE:
      return backward_row_zero_sse(x1, om, ox);
      break;
    case AVX:
      return backward_row_zero_avx(x1, om, ox);
      break;
    case AVX512:
      return backward_row_zero_avx512(x1, om, ox);
      break;
    case NEON:
      return backward_row_zero_neon(x1, om, ox);
      break;
    case NEON64:
      return backward_row_zero_neon64(x1, om, ox);
      break;
   default:
      p7_Fail("Unrecognized SIMD type passed to backward_row_zero");  
  }
}


/* save_debug_row_pp()
 * 
 * Debugging only. Transfer posterior decoding values from a
 * vectorized row to appropriate row of <ox->pp>, as probabilities.
 *
 * Zero wherever there's no emission:
 *   all values on row 0 
 *   all D states
 *   all E,B states
 * Zero in row positions that can't be reached in some state:
 *   all insert states on row 1
 *   J(1),C(1)
 *   N(L), J(L) on row L
 * Zero wherever there's no state:
 *   D1, IM states
 */
static void
save_debug_row_pp(P7_CHECKPTMX *ox, debug_print *dpc, int i)
{
  switch(ox->simd){
    case SSE:
      return save_debug_row_pp_sse(ox, dpc, i);
      break;
    case AVX:
      return save_debug_row_pp_avx(ox, dpc, i);
      break;
    case AVX512:
      return save_debug_row_pp_avx512(ox, dpc, i);
      break;
    case NEON:
      return save_debug_row_pp_neon(ox, dpc, i);
      break;
    case NEON64:
      return save_debug_row_pp_neon64(ox, dpc, i);
      break;
   default:
      p7_Fail("Unrecognized SIMD type passed to save_debug_row_pp");  
  }
}

/* save_debug_row_fb()
 * 
 * Debugging only. Transfer posterior decoding values (sparse scaled,
 * prob space) from a vectorized row, to appropriate row of <ox->fwd>
 * or <ox->bck> (log space, inclusive of partial sum of scalefactors);
 * <ox->fwd> and <ox->bck> should be identical (within numerical error
 * tolerance) to a reference implementation Forward/Backward in log
 * space.
 */
static void
save_debug_row_fb(P7_CHECKPTMX *ox, P7_REFMX *gx, debug_print *dpc, int i, float totscale)
{
  switch(ox->simd){
    case SSE:
      return save_debug_row_fb_sse(ox, gx, dpc, i, totscale);
      break;
    case AVX:
      return save_debug_row_fb_avx(ox, gx, dpc, i, totscale);
      break;
    case AVX512:
      return save_debug_row_fb_avx512(ox, gx, dpc, i, totscale);
      break;
    case NEON:
      return save_debug_row_fb_neon(ox, gx, dpc, i, totscale);
      break;
    case NEON64:
      return save_debug_row_fb_neon64(ox, gx, dpc, i, totscale);
      break;
   default:
      p7_Fail("Unrecognized SIMD type passed to save_debug_row_fb");  
  } 
}
#endif
/*---------------- end, debugging tools -------------------------*/


/*****************************************************************
 * 4. Stats driver: memory requirement profiling 
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
      p7_sparsemask_Reinit(sm, gm->M, sq->n);

      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);

      /* Filter insig hits, partially simulating the real pipeline */
      p7_MSVFilter(sq->dsq, sq->n, om, fx, &msvsc);
      msvsc = (msvsc - nullsc) / eslCONST_LOG2;
      P     =  esl_gumbel_surv(msvsc,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      if (P > 0.02) goto NEXT_SEQ;

      p7_ForwardFilter (sq->dsq, sq->n, om, ox, &fraw);
      p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH);

      /* Calculate minimum memory requirements for each step */
      msvmem = (double) ( P7_NVB(om->M) * sizeof(__m128i))    / 1024.;  
      vfmem  = (double) p7_filtermx_MinSizeof(om->M)          / 1024.;
      fbmem  = (double) p7_checkptmx_MinSizeof(om->M, sq->n)  / 1024. / 1024.;
      smmem  = (double) p7_sparsemask_MinSizeof(sm)           / 1024. / 1024.;
      spmem  = (double) p7_sparsemx_MinSizeof(sm)             / 1024. / 1024.;
      refmem = (double) p7_refmx_MinSizeof(om->M, sq->n)      / 1024. / 1024.;

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
      p7_checkptmx_Reuse(ox);
      p7_filtermx_Reuse(fx);
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
 * 5. Benchmark
 *****************************************************************/
/* Difference between gcc -g vs. icc -O3 is large! */

#ifdef p7FWDFILTER_BENCHMARK
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

      p7_sparsemask_Reinit(sm, om->M, L);

      p7_ForwardFilter(dsq, L, om, ox, &sc);
      if (! esl_opt_GetBoolean(go, "-F")) 
	{
	  p7_BackwardFilter(dsq, L, om, ox, sm, p7_SPARSIFY_THRESH);
	  esl_vec_IReverse(sm->kmem, sm->kmem, sm->ncells);
	}

      p7_checkptmx_Reuse(ox);
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
  char         msg[]  = "fbfilter scores unit test failed";
  P7_HARDWARE *hw;
  if ((hw = p7_hardware_Create ()) == NULL)  p7_Fail("Couldn't get HW information data structure"); 
  P7_HMM      *hmm    = NULL;
  P7_PROFILE  *gm     = NULL;
  P7_OPROFILE *om     = NULL;
  ESL_DSQ     *dsqmem = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_DSQ     *dsq    = NULL;
  int          tL     = 0;
  ESL_SQ      *sq     = esl_sq_CreateDigital(abc);
  P7_CHECKPTMX *ox    = p7_checkptmx_Create(M, L, ramlimit, hw->simd);
  P7_REFMX    *fwd    = p7_refmx_Create   (M, L);
  P7_REFMX    *bck    = p7_refmx_Create   (M, L);
  P7_REFMX    *pp     = p7_refmx_Create   (M, L);
  P7_SPARSEMASK *sm   = p7_sparsemask_Create(M, L, hw->simd);
  float        tol2   = ( p7_logsum_IsSlowExact() ? 0.001  : 0.1);   /* absolute agreement of reference (log-space) and vector (prob-space) depends on whether we're using LUT-based logsum() */
  float fsc1, fsc2;
  float bsc2;
#ifdef p7_DEBUGGING
  float        bsc1;
  float        tol1   = 0.0001;	                                     /* forward and backward scores from same implementation type should agree with high tolerance */
  float        ptol   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01);  /* posterior decoding values differ by no more than this */
#endif

#ifdef p7_DEBUGGING
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
      if ( p7_checkptmx_GrowTo(ox,  M, tL)    != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_Reinit(sm,M, tL)     != eslOK) esl_fatal(msg);

      p7_ForwardFilter (dsq, tL, om, ox, &fsc1);
      p7_BackwardFilter(dsq, tL, om, ox,  sm, p7_SPARSIFY_THRESH);

      p7_ReferenceForward (dsq, tL, gm, fwd,  &fsc2);
      p7_ReferenceBackward(dsq, tL, gm, bck,  &bsc2);
      p7_ReferenceDecoding(dsq, tL, gm, fwd, bck, pp);

#ifdef p7_DEBUGGING
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

#ifdef p7_DEBUGGING
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
      p7_checkptmx_Reuse(ox);
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
#ifdef p7_DEBUGGING
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
#ifdef p7_DEBUGGING
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

#ifdef p7_DEBUGGING
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
#ifdef p7_DEBUGGING
  /* When the p7_DEBUGGING flag is up, <ox> matrix has the ability to
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

      p7_checkptmx_GrowTo  (ox, om->M, sq->n); 
      p7_sparsemask_Reinit(sm, gm->M, sq->n);

      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);
    
      p7_ForwardFilter (sq->dsq, sq->n, om, ox, &fraw);
      p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH);

      p7_ReferenceForward (sq->dsq, sq->n, gm, gx, &gfraw);
      p7_ReferenceBackward(sq->dsq, sq->n, gm, gx, &gbraw);

      bsc = 0.0;		/* Backward score only available in debugging mode */
#ifdef p7_DEBUGGING
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
#ifdef p7_DEBUGGING
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
      p7_checkptmx_Reuse(ox);
    }

#ifdef p7_DEBUGGING
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

/*****************************************************************
 * 9. Notes
 *****************************************************************/

/* [a.]  On debugging and testing methods.
 * 
 *    When compiled with <p7_DEBUGGING> (specifically, when
 *    fwdfilter.c and p7_checkptmx.[ch] are thus compiled), the
 *    <P7_CHECKPTMX> structure is augmented with additional fields for
 *    debugging dumps and unit test comparisons.
 *    
 *   %% Dumping vector matrices for examination
 *      The values from the vectorized <P7_CHECKPTMX> can be dumped
 *      during DP calculations by calling <p7_checkptmx_SetDumpMode()>
 *      on the object. Dumping has to happen during DP, not after,
 *      because of the way checkpointing discards rows as it goes (and
 *      for posterior decoding, the implementation never stores a row
 *      at all). Dumped rows are prefixed by a tag "f1 O", "f1 X", "f2
 *      O", "f2 X", or "bck", indicating backward (bck), first pass
 *      Forward (f1), second pass Forward (f2), checkpointed rows (O)
 *      that get saved and recalled, and discarded rows (X) that get
 *      recalculated in the second Forward pass.
 *      
 *      This capability is most useful for examining small DP matrices
 *      by hand; see fwdfilter_example, -D option.
 *      
 *   %% Saving matrices for comparison to reference
 *      With the <p7_DEBUGGING> compile flag, a caller may
 *      additionally provide an allocated <P7_REFMX> to the
 *      <P7_CHECKPTMX>, to enable storage of all DP matrix values in a
 *      form suitable for a <p7_refmx_Compare()> call against <P7_REFMX>
 *      DP matrices calculated by the reference implementation.
 *      Caller does something like <ox->fwd = p7_refmx_Create(M,L)> to
 *      save a Forward matrix, and/or analogously for <ox->bck> and/or
 *      <ox->pp> for Backward and Decoding.
 *      
 *      This capability is most useful for unit tests and automated
 *      comparison to the reference implementation. See utest_scores().
 *      
 *   %% High-precision comparison to reference
 *      Normally the reference implementation uses a table-driven
 *      log-sum-exp approximation (see misc/logsum.c), in order to do
 *      stable numerical calculations in log space. This introduces
 *      nonnegligible numerical error into DP calculations, so
 *      comparisons between a probability space vector implementation
 *      and the reference implementation must allow a large amount of
 *      numeric slop. At a cost of about 20x in speed, if the
 *      p7_LOGSUM_SLOWEXACT flag is compiled in, the p7_FLogsum()
 *      function uses the (more) exact calculation, allowing DP cells
 *      values to be compared more stringently.
 *      
 *      This capability is reflected in unit tests that set tolerances
 *      for floating-point comparison, after checking the flag with a
 *      <p7_logsum_IsSlowExact()> call. See utest_scores() for
 *      example.
 */

/* [b.] Running time, in theory and in practice.
 *
 *    Checkpointing requires more time in theory, but in practice you
 *    probably won't notice. The checkpointing method requires
 *    recalculation of Forward rows that weren't checkpointed, meaning
 *    up to two Forwards passes (plus one Backwards and one posterior
 *    decoding pass) over the DP matrix, rather than one each. First,
 *    using the Forward score as a filter minimizes this penalty:
 *    relatively few sequences pass on to the BackwardFilter
 *    call. Second, memory management in the checkpointed P7_CHECKPTMX
 *    structure uses "partial checkpointing" to minimize the use of
 *    checkpointing; for most comparisons, of all but the longest
 *    query/target combinations, DP calculations will fit in the
 *    available memory, and checkpointing is not invoked. And finally,
 *    the implementation here has been further optimized, such that
 *    it's actually slightly faster (with checkpointing) than the
 *    original HMMER3.0 implementation of Forward and Backward without
 *    checkpointing.
 * 
 */

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/




                                          
