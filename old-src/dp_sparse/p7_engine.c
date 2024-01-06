#include <p7_config.h>

/* SIMD-vectorized acceleration filters, local only: */
#include "dp_vector/ssvfilter.h"          // SSV primary acceleration filter
#include "dp_vector/vitfilter.h"          // Viterbi secondary acceleration filter
#include "dp_vector/fwdfilter.h"          // Sparsification w/ checkpointed local Forward/Backward
#include "dp_vector/p7_filtermx.h"        // DP matrix for SSV, VF
#include "dp_vector/p7_checkptmx.h"       // DP matrix for FB

#include "search/modelconfig.h"

/* Sparse DP, dual-mode glocal/local:    */
#include "dp_sparse/sparse_fwdback.h"     // sparse Forward/Backward
#include "dp_sparse/sparse_viterbi.h"     // sparse Viterbi
#include "dp_sparse/sparse_decoding.h"    // sparse Decoding
#include "dp_sparse/sparse_anchors.h"     // most probable anchor set (MPAS) 

/* Sparse anchor-set-constrained (ASC): */
#include "dp_sparse/sparse_asc_fwdback.h"   // ASC Forward/Backward
#include "dp_sparse/sparse_envelopes.h"     // Envelope inference
#include "dp_sparse/sparse_null2.h"         // Null2 score correction
#include "dp_sparse/sparse_aec_align.h"     // anchor/envelope constrained alignment

#include "dp_sparse/p7_engine.h"  // FIXME: we'll move the engine somewhere else, I think

#include "server/p7_hitlist.h" //  This probably wants to move somewhere else once we figure out how command-line search will work
#include "easel.h"
#include "esl_random.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"

/*****************************************************************
 * 1. P7_ENGINE_PARAMS: config/control parameters for the Engine.
 *****************************************************************/

P7_ENGINE_PARAMS *
p7_engine_params_Create(P7_MPAS_PARAMS *mpas_params)
{
  P7_ENGINE_PARAMS *prm = NULL;
  int               status;

  ESL_ALLOC(prm, sizeof(P7_ENGINE_PARAMS));
  prm->rng_seed          = p7_ENGINE_FIXED_SEED;     // Defaults are set in p7_config.h.in
  prm->rng_reproducible  = p7_ENGINE_REPRODUCIBLE; 
  prm->sparsify_ramlimit = p7_SPARSIFY_RAMLIMIT;     
  prm->sparsify_thresh   = p7_SPARSIFY_THRESH;   
  prm->do_biasfilter     = p7_ENGINE_DO_BIASFILTER;
  prm->mpas_params       = (mpas_params ? mpas_params : NULL);
  return prm;

 ERROR:
  p7_engine_params_Destroy(prm);
  return NULL;
}

void
p7_engine_params_Destroy(P7_ENGINE_PARAMS *prm)
{
  if (prm) {
    if (prm->mpas_params) p7_mpas_params_Destroy(prm->mpas_params);
  }
  free(prm);
}


/*****************************************************************
 * 2. P7_ENGINE_STATS
 *****************************************************************/ 

P7_ENGINE_STATS *
p7_engine_stats_Create(void)
{
  P7_ENGINE_STATS *stats = NULL;
  int              status;

  ESL_ALLOC(stats, sizeof(P7_ENGINE_STATS));
  stats->n_past_ssv  = 0;
  stats->n_past_bias = 0;
  stats->n_ran_vit   = 0;
  stats->n_past_vit  = 0;
  stats->n_past_fwd  = 0;
  return stats;

 ERROR:
  p7_engine_stats_Destroy(stats);
  return NULL;
}

void
p7_engine_stats_Destroy(P7_ENGINE_STATS *stats)
{
  if (stats) free(stats);
}


/*****************************************************************
 * 3. P7_ENGINE
 *****************************************************************/

P7_ENGINE *
p7_engine_Create(const ESL_ALPHABET *abc, P7_ENGINE_PARAMS *prm, P7_ENGINE_STATS *stats, int M_hint, int L_hint)
{
  P7_ENGINE *eng               = NULL;
  uint32_t   rng_seed          = (prm ? prm->rng_seed          : p7_ENGINE_FIXED_SEED);
  int        sparsify_ramlimit = (prm ? prm->sparsify_ramlimit : p7_SPARSIFY_RAMLIMIT);
  int        status;

  /* level 0 */
  ESL_ALLOC(eng, sizeof(P7_ENGINE));
  eng->rng   = NULL;
  eng->fx    = NULL;
  eng->cx    = NULL;
  eng->sm    = NULL;

  eng->used_main = FALSE;
  eng->sxf   = NULL;
  eng->sxd   = NULL;
  eng->asf   = NULL;
  eng->asb   = NULL;
  eng->asd   = NULL;
  eng->anch  = NULL;
  eng->vanch = NULL;
  eng->ahash = NULL;
  eng->env   = NULL;
  eng->tr    = NULL;

  /* Use params if provided, else create defaults.
   * Add optional stats collection if provided.
   */
  eng->params = (prm   ? prm   : NULL);
  eng->stats  = (stats ? stats : NULL);
  eng->rng    = esl_randomness_CreateFast(rng_seed);

  /* All DP algorithms resize their input mx's as needed. 
   * The initial allocation can be anything.
   */
  eng->fx    = p7_filtermx_Create  (M_hint);
  eng->cx    = p7_checkptmx_Create (M_hint, L_hint, ESL_MBYTES(sparsify_ramlimit));
  eng->sm    = p7_sparsemask_Create(M_hint, L_hint);
  eng->sxf   = p7_sparsemx_Create(eng->sm);
  eng->sxd   = p7_sparsemx_Create(eng->sm);
  eng->asf   = p7_sparsemx_Create(eng->sm);
  eng->asb   = p7_sparsemx_Create(eng->sm);
  eng->asd   = p7_sparsemx_Create(eng->sm);
  eng->anch  = p7_anchors_Create();
  eng->vanch = p7_anchors_Create();
  eng->ahash = p7_anchorhash_Create();
  eng->env   = p7_envelopes_Create();
  eng->tr    = p7_trace_CreateWithPP();

  eng->wrkM  = NULL;  // wrkM is handled by bypass idiom; starts NULL, reallocated as needed.
  ESL_ALLOC(eng->wrkKp, abc->Kp * sizeof(float));

  eng->nullsc = 0.;
  eng->biassc = 0.;
  eng->sfsc   = 0.;
  eng->vfsc   = 0.;

  eng->ffsc   = 0.;
  eng->vsc    = 0.;
  eng->fsc    = 0.;
  eng->asc_f  = 0.;

  eng->F1     = 0.02;
  eng->F2     = 0.001;
  eng->F3     = 1e-5;

  return eng;

 ERROR:
  p7_engine_Destroy(eng);
  return NULL;
}


int
p7_engine_Reuse(P7_ENGINE *eng)
{
  int      rng_reproducible = (eng->params ? eng->params->rng_reproducible : p7_ENGINE_REPRODUCIBLE);
  uint32_t rng_seed         = (eng->params ? eng->params->rng_seed         : p7_ENGINE_FIXED_SEED);
  int status;

  if (rng_reproducible) 
    esl_randomness_Init(eng->rng, rng_seed);

  if ((status = p7_sparsemask_Reuse(eng->sm))    != eslOK) return status;

  /* Most Reuse()'s are cheap, but the p7_anchorhash_Reuse() is a little
   * expensive. That's why we avoid Reuse()'ing the structures that only
   * the main engine uses.
   */
  if (eng->used_main)
    {
      if ((status = p7_sparsemx_Reuse  (eng->sxf))   != eslOK) return status;
      if ((status = p7_sparsemx_Reuse  (eng->sxd))   != eslOK) return status;
      if ((status = p7_sparsemx_Reuse  (eng->asf))   != eslOK) return status;
      if ((status = p7_sparsemx_Reuse  (eng->asb))   != eslOK) return status;
      if ((status = p7_sparsemx_Reuse  (eng->asd))   != eslOK) return status;
      if ((status = p7_anchors_Reuse   (eng->anch))  != eslOK) return status;
      if ((status = p7_anchors_Reuse   (eng->vanch)) != eslOK) return status;
      if ((status = p7_anchorhash_Reuse(eng->ahash)) != eslOK) return status;  
      if ((status = p7_envelopes_Reuse (eng->env))   != eslOK) return status;  
      if ((status = p7_trace_Reuse     (eng->tr))    != eslOK) return status;
      /* wrkM and wrkKp are scratch workspaces, don't need to be reused/reinitialized */
    }
  eng->used_main = FALSE;

  eng->nullsc = 0.;
  eng->biassc = 0.;
  eng->sfsc   = 0.;
  eng->vfsc   = 0.;

  eng->ffsc   = 0.;
  eng->vsc    = 0.;
  eng->fsc    = 0.;
  eng->asc_f  = 0.;

  /* Don't clear out the current hit chunk -- that needs to stay around until we finish the entire search we're doing.
     Caveat is that the chunk needs to be cleaned manually between searches */

  /* F1, F2, F3 are constants, they don't need to be reset. */
  return eslOK;
}

void
p7_engine_Destroy(P7_ENGINE *eng)
{
  if (eng)
    {
      if (eng->rng)   esl_randomness_Destroy(eng->rng);
      if (eng->fx)    p7_filtermx_Destroy   (eng->fx);
      if (eng->cx)    p7_checkptmx_Destroy  (eng->cx);
      if (eng->sm)    p7_sparsemask_Destroy (eng->sm);
      if (eng->sxf)   p7_sparsemx_Destroy   (eng->sxf);
      if (eng->sxd)   p7_sparsemx_Destroy   (eng->sxd);
      if (eng->asf)   p7_sparsemx_Destroy   (eng->asf);
      if (eng->asb)   p7_sparsemx_Destroy   (eng->asb);
      if (eng->asd)   p7_sparsemx_Destroy   (eng->asd);
      if (eng->vanch) p7_anchors_Destroy    (eng->vanch);
      if (eng->anch)  p7_anchors_Destroy    (eng->anch);
      if (eng->ahash) p7_anchorhash_Destroy (eng->ahash);
      if (eng->env)   p7_envelopes_Destroy  (eng->env);
      if (eng->tr)    p7_trace_Destroy      (eng->tr);

      if (eng->wrkM)  free(eng->wrkM);
      if (eng->wrkKp) free(eng->wrkKp);

      if (eng->params) p7_engine_params_Destroy(eng->params);
      if (eng->stats)  p7_engine_stats_Destroy (eng->stats);
    }
  free(eng);
}


/*****************************************************************
 * 4. The engines themselves.
 *****************************************************************/

/* Function:  p7_engine_Overthruster()
 * Synopsis:  The acceleration filters.
 *
 * Purpose:   The acceleration filters: using
 *            engine <eng>, compare sequence <dsq> of length <L> to
 *            the vectorized profile <om>, calculating log-odds raw
 *            scores using null model <bg>.
 * 
 *            Return <eslOK> if the sequence passes the filters,
 *            <eslFAIL> if it doesn't.
 *            
 *            Caller must have configured length models in <om> and
 *            <bg> already.
 *            
 *            Upon return, the sparse mask <eng->sm> has been calculated,
 *            and the following score fields are set in the
 *            <eng>. Later steps may not have run, depending on earlier
 *            steps:
 *               <nullsc>  : null model raw score, nats;     = 0 if not reached.
 *               <biassc>  : ad hoc bias filter score, nats; = 0 if not reached.
 *               <sfsc>    : SSV filter raw score, nats;     = -eslINFINITY if not reached.
 *               <vfsc>    : Viterbi filter raw score, nats; = -eslINFINITY if not reached.
 *               <ffsc>    : Forward filter raw score, nats; = -eslINFINITY if not reached.
 * 
 *            If the bias filter is off, biassc = nullsc.
 *
 *            If the SSV filter score is so high that it satisfies
 *            both F1 and F2 thresholds, the Viterbi filter step is
 *            skipped. 
 *            
 *            Caller can recalculate output scores in bits
 *            by (raw_sc - biassc) / eslCONST_LOG2.
 *            This works even for calculations that aren't
 *            reached (in which case the score will come out as -inf).
 *            
 *            If the <eng> is collecting statistics in a non-NULL
 *            <eng->stats>, its <n_past_ssv>, <n_past_bias>,
 *            <n_past_vit>, <n_ran_vit>, <n_past_fwd> counters can advance here.
 *            
 *            The O(M) filter DP matrix <eng->fx> and the O(M sqrt L) checkpoint
 *            matrix <eng->cx> may be reallocated here.
 *            
 * Throws:    <eslEMEM> if a DP matrix reallocation fails.           
 *            
 */
int
p7_engine_Overthruster(P7_ENGINE *eng, ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_BG *bg)
{
  int   do_biasfilter    = (eng->params ? eng->params->do_biasfilter   : p7_ENGINE_DO_BIASFILTER);
  float sparsify_thresh  = (eng->params ? eng->params->sparsify_thresh : p7_SPARSIFY_THRESH);
  float seq_score;
  float P;
  int   status;

  if (L == 0) return eslFAIL;

  if ((status = p7_bg_NullOne(bg, dsq, L, &(eng->nullsc))) != eslOK) return status; 

  /* First level: SSV filter */
  status = p7_SSVFilter(dsq, L, om, &(eng->sfsc));
  if (status != eslOK && status != eslERANGE) return status;

  seq_score = (eng->sfsc - eng->nullsc) / eslCONST_LOG2;          
  P = esl_gumbel_surv(seq_score,  om->evparam[p7_SMU],  om->evparam[p7_SLAMBDA]);
  if (P > eng->F1) return eslFAIL;

  if (eng->stats) eng->stats->n_past_ssv++;

  /* Biased composition HMM, ad hoc, acts as a modified null */
  if  (do_biasfilter)
    {
      if ((status = p7_bg_FilterScore(bg, dsq, L, &(eng->biassc))) != eslOK) return status;
      seq_score = (eng->sfsc - eng->biassc) / eslCONST_LOG2;
      P = esl_gumbel_surv(seq_score,  om->evparam[p7_SMU],  om->evparam[p7_SLAMBDA]);
      if (P > eng->F1) return eslFAIL;
    }
  else eng->biassc = eng->nullsc;
  if (eng->stats) eng->stats->n_past_bias++;

  // TODO: in scan mode, you have to load the rest of the oprofile now,
  // configure its length model, and get GA/TC/NC thresholds.

  /* Second level: ViterbiFilter(), multihit with <om> */
  if (P > eng->F2)
    {
      if (eng->stats) eng->stats->n_ran_vit++;

      status = p7_ViterbiFilter(dsq, L, om, eng->fx, &(eng->vfsc));  
      if (status != eslOK && status != eslERANGE) return status;

      seq_score = (eng->vfsc - eng->biassc) / eslCONST_LOG2;
      P  = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
      if (P > eng->F2) return eslFAIL;
    }
  if (eng->stats) eng->stats->n_past_vit++;


  /* Checkpointed vectorized Forward, local-only.
   */
  status = p7_ForwardFilter(dsq, L, om, eng->cx, &(eng->ffsc));
  if (status != eslOK) return status;

  seq_score = (eng->ffsc - eng->biassc) / eslCONST_LOG2;
  P  = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  if (P > eng->F3) return eslFAIL;
  if (eng->stats) eng->stats->n_past_fwd++;

  /* Sequence has passed all acceleration filters.
   * Calculate the sparse mask, by checkpointed vectorized decoding.
   */
  p7_BackwardFilter(dsq, L, om, eng->cx, eng->sm, sparsify_thresh);

  return eslOK;
}


  // om is assumed to be complete, w/ GA/NC/TC thresholds set, and w/ length model set.
  // Use dsq, L -- not sq -- so subseqs can be processed (should help longtarget/nhmmer)
/* Function:  
 * Synopsis:  
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   <eslOK> on success, and the Engine <eng> contains results in eng->*:
 *               fsc  : sparse Forward raw score for the whole seq
 *               asc_f: sparse ASC Forward raw score, also for the whole sequence
 *               tr:  optimal alignment of the target sequence 
 *               env: envelope information:
 *                    D : number of domains
 *                    and for each domain 1..D in env->arr[d].*:
 *                       i0,k0    : anchor
 *                       ia,ib    : envelope on seq
 *                       ka,kb    : alignment start/end on model
 *                       oa,ob    : outer envelope 
 *                       env_sc   : envelope raw score
 *                       null2_sc : envelope null2 score correction
 *                  
 *
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
p7_engine_Main(P7_ENGINE *eng, ESL_DSQ *dsq, int L, P7_PROFILE *gm)
{
  P7_MPAS_PARAMS *mpas_params     = (eng->params && eng->params->mpas_params ? eng->params->mpas_params : NULL);
  eng->used_main = TRUE;  // This flag causes engine_Reuse() to reuse all of the engine, 
                          // not just the structures used by the Overthruster.

  /* First pass analysis 
   * Uses two sparse matrices: <sxf>, <sxd>,
   * and also gets the unconstrained Viterbi trace, <tr>
   */
  p7_SparseViterbi (dsq, L, gm, eng->sm,  eng->sxf, eng->tr, &(eng->vsc));
  p7_SparseForward (dsq, L, gm, eng->sm,  eng->sxf,          &(eng->fsc));
  p7_SparseBackward(dsq, L, gm, eng->sm,  eng->sxd,          /*bsc=*/NULL);
  p7_SparseDecoding(dsq, L, gm, eng->sxf, eng->sxd,          eng->sxd);

  /* MPAS algorithm for finding the anchor set */
  p7_sparse_anchors_SetFromTrace(eng->sxd, eng->tr, eng->vanch);
  p7_trace_Reuse(eng->tr);
  p7_sparse_Anchors(eng->rng, dsq, L, gm,
        eng->vsc, eng->fsc, eng->sxf, eng->sxd, eng->vanch,
        eng->tr, &(eng->wrkM), eng->ahash,  
        eng->asf, eng->anch, &(eng->asc_f), 
        mpas_params);

  /* Remaining ASC calculations. MPAS already did <asf> for us. */
  p7_sparse_asc_Backward(dsq, L, gm, eng->anch->a, eng->anch->D, eng->sm,    eng->asb, /*asc_b=*/NULL);
  p7_sparse_asc_Decoding(dsq, L, gm, eng->anch->a, eng->anch->D, eng->asc_f, eng->asf, eng->asb, eng->asd);

  /* Envelope determination */
  p7_sparse_Envelopes(dsq, L, gm, eng->anch->a, eng->anch->D, eng->asf, eng->asd, eng->env);

  /* null2 score corrections on each envelope.
   * Store them in <env>: env->arr[d].null2_sc.     ($r_d$, in our print documentation)
   */
  p7_sparse_Null2(dsq, L, gm, eng->asd, eng->env, &(eng->wrkM), eng->wrkKp);

  /* Optimal alignments for each envelope */
  p7_sparsemx_Reuse(eng->sxf);                                      // sxf overwritten with AEC DP matrix
  p7_trace_Reuse(eng->tr); 
  p7_sparse_aec_Align(gm, eng->asd, eng->env, eng->sxf, eng->tr);

  /* Pick up posterior probability annotation for the alignment */
  p7_sparsemx_TracePostprobs(eng->sxd, eng->tr);

  return eslOK;
}

/* Calls the engine to compare a sequence to an HMM.  Heavily cribbed from Seans 0226-px code*/
// returns 0 if no hit was found, 1 if a hit was found
int p7_engine_Compare_Sequence_HMM(P7_ENGINE *eng, ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg){

    int status;
    // reset the models for the length of this sequence
    p7_bg_SetLength(bg, L);           
    p7_oprofile_ReconfigLength(om, L);
    
    // First, the overthruster (filters)
    status = p7_engine_Overthruster(eng, dsq, L, om, bg);  
    if (status == eslFAIL) { // filters say no match
      p7_engine_Reuse(eng);
      return 0;
    }

    // if we get here, run the full comparison
    p7_profile_SetLength(gm, L);
    status = p7_engine_Main(eng, dsq, L, gm); 

    p7_engine_Reuse(eng);

    return(1); // for now, everything that reaches the main stage is a hit
}



/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7ENGINE_EXAMPLE

#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"


/* trace_dump_runlengths()
 * Used in testing zigar compression strategies.
 */
static int
trace_dump_runlengths(P7_TRACE *tr)
{
  int nrun      = 0;
  int z;

  for (z = 1; z < tr->N; z++)                                // starts at z=1 (S state); we use z-1.
    if ( p7_trace_IsMain(tr->st[z]) || tr->st[z] == p7T_E)   // we're in a domain; or ending one.
      {
	if (! p7_trace_IsMain(tr->st[z-1]) ) nrun = 1;       // start of a new domain? reset counter.
	else if (tr->st[z] == tr->st[z-1])   nrun++;         // extending a run of the same state? bump counter.
	else {                                               // finished previous run (including E)? print cigar element.
	  printf("%5d %-2s%c", nrun, p7_trace_DecodeStatetype(tr->st[z-1]), tr->st[z] == p7T_E ? '\n' : ' ');
	  nrun = 1;
	}
      }
  return eslOK;
}

static int
trace_dump_postprobs(P7_TRACE *tr)
{
  int z;

  for (z = 1; z < tr->N; z++)
    if      (p7_trace_IsM(tr->st[z]) || p7_trace_IsI(tr->st[z])) printf("%.6f ", tr->pp[z]);
    else if (tr->st[z] == p7T_E)                                 printf("\n");
  return eslOK;
}



static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",            0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                   0 },
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump trace structure(s)",                         0 },
  { "--Zcig",    eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump CIGAR strings",                              0 },  // xref projects/zigar
  { "--Zpp",     eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump posterior probability lines",                0 },  // ditto
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example driver for the engine";

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
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_ENGINE      *eng     = NULL;
  int             status;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);

  while ( p7_hmmfile_Read(hfp, &abc, &hmm) == eslOK) 
    {
      /* Configure a profile from the HMM */
      bg = p7_bg_Create(abc);
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_profile_Config   (gm, hmm, bg);
      p7_oprofile_Convert (gm, om);

      p7_bg_SetFilter(bg, om->M, om->compo);

      /* Open sequence file */
      sq     = esl_sq_CreateDigital(abc);
      status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
      if      (status == eslENOTFOUND) p7_Fail("No such file.");
      else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
      else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
      else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
      /* Create the comparison engine */
      eng   = p7_engine_Create(abc, /*params=*/NULL, /*stats=*/NULL, gm->M, /*L_hint=*/400);

      /* For each sequence in <sqfile>: */
      while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
	{ 
	  p7_bg_SetLength(bg, sq->n);
	  p7_oprofile_ReconfigLength(om, sq->n);

	  status = p7_engine_Overthruster(eng, sq->dsq, sq->n, om, bg);
	  if      (status == eslFAIL) { 
	    p7_engine_Reuse(eng);
	    esl_sq_Reuse(sq); 
	    continue; 
	  }
	  else if (status != eslOK)   p7_Fail("overthruster failed with code %d\n", status);

	  p7_profile_SetLength(gm, sq->n);
	  status = p7_engine_Main(eng, sq->dsq, sq->n, gm);

	  if (esl_opt_GetBoolean(go, "-T"))     p7_trace_DumpAnnotated(stdout, eng->tr, gm, sq->dsq);
	  if (esl_opt_GetBoolean(go, "--Zcig")) trace_dump_runlengths(eng->tr);
	  if (esl_opt_GetBoolean(go, "--Zpp"))  trace_dump_postprobs(eng->tr);

	  p7_engine_Reuse(eng);
	  esl_sq_Reuse(sq);
	}
      if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s)\n%s\n",	   sqfp->filename, sqfp->get_error(sqfp));     
      else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s", status, sqfp->filename);
  
      p7_engine_Destroy(eng);     eng  = NULL;
      esl_sqfile_Close(sqfp);     sqfp = NULL;
      esl_sq_Destroy(sq);         sq   = NULL;
      p7_oprofile_Destroy(om);    om   = NULL;
      p7_profile_Destroy(gm);     gm   = NULL;
      p7_bg_Destroy(bg);          bg   = NULL;
      p7_hmm_Destroy(hmm);        hmm  = NULL;
    }

  p7_hmmfile_Close(hfp);      
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}

#endif /*p7ENGINE_EXAMPLE*/



