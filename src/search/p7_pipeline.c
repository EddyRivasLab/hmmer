/* H3's accelerated seq/profile comparison pipeline
 *  
 * Contents:
 *   1. P7_PIPELINE: allocation, initiazation, destruction
 *   2. Setting pipeline for next target or model
 *   3. Testing bitscore/E-value against reporting/inclusion thresholds
 *   4. Statistics output from a completed pipeline.
 *   5. The main pipeline call: p7_Pipeline().
 *   6. The pipeline speciazed for nhmmer/longtarget.
 *   7. Example 1: search mode (in a sequence db)
 *   8. Example 2: scan mode (in an HMM db)
 *   9. Copyright and license information
 * 
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

#include "easel.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "base/p7_bg.h"
#include "base/p7_domain.h"
#include "base/p7_profile.h"   /* used by the hmmscan workaround */
#include "base/p7_scoredata.h"
#include "base/p7_tophits.h"
#include "base/p7_hmmwindow.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/io.h"
#include "dp_vector/msvfilter.h"
#include "dp_vector/vitfilter.h"
#include "dp_vector/fwdfilter.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_viterbi.h"
#include "dp_sparse/sparse_fwdback.h"
#include "dp_sparse/sparse_decoding.h"
#include "dp_sparse/sparse_masstrace.h"
#include "dp_sparse/sparse_envscore.h"
#include "dp_sparse/sparse_null2.h"


#include "misc/logsum.h"

#include "search/modelconfig.h" /* used by the hmmscan workaround */
#include "search/p7_pipeline.h"


/* SRE: FIXME 3.1 in progress */
static int workaround_get_profile(P7_PIPELINE *pli, const P7_OPROFILE *om, const P7_BG *bg, P7_PROFILE **ret_gm);

typedef struct {
  P7_BG            *bg;
  P7_PROFILE       *gm;
  P7_OPROFILE      *om;
  P7_TRACE         *trc;
  P7_SPARSEMASK    *sm;
  P7_SPARSEMX      *sxf;
  P7_SPARSEMX      *sxb;
  P7_SPARSEMX      *sxd;
  P7_SPARSEMX      *sxx;
  P7_MASSTRACE     *mt;
  float            *scores;

} P7_PIPELINE_TMP_OBJS;

/*****************************************************************
 * 1. The P7_PIPELINE object: allocation, initiazation, destruction.
 *****************************************************************/

/* Function:  p7_pipeline_Create()
 * Synopsis:  Create a new accelerated comparison pipeline.
 *
 * Purpose:   Given an application configuration structure <go>
 *            containing certain standardized options (described
 *            below), some initial guesses at the model size <M_hint>
 *            and sequence length <L_hint> that will be processed,
 *            an nhmmer flag <do_longtargets> (TRUE for nhmmer's
 *            long target sequence scanning; FALSE otherwise);
 *            and a <mode> that can be either <p7_SCAN_MODELS> or
 *            <p7_SEARCH_SEQS> depending on whether we're searching one sequence
 *            against a model database (hmmscan mode) or one model
 *            against a sequence database (hmmsearch mode); create new
 *            pipeline object.
 *
 *            In search mode, we would generally know the length of
 *            our query profile exactly, and would pass <om->M> as <M_hint>;
 *            in scan mode, we generally know the length of our query
 *            sequence exactly, and would pass <sq->n> as <L_hint>.
 *            Targets will come in various sizes as we read them,
 *            and the pipeline will resize any necessary objects as
 *            needed, so the other (unknown) length is only an
 *            initial allocation.
 *            
 *            The configuration <go> must include settings for the 
 *            following options:
 *            
 *            || option      ||            description                    || usually  ||
 *            | --noali      |  don't output alignments (smaller output)   |   FALSE   |
 *            | -E           |  report hits <= this E-value threshold      |    10.0   |
 *            | -T           |  report hits >= this bit score threshold    |    NULL   |
 *            | -Z           |  set initial hit search space size          |    NULL   |
 *            | --domZ       |  set domain search space size               |    NULL   |
 *            | --domE       |  report domains <= this E-value threshold   |    10.0   |
 *            | --domT       |  report domains <= this bit score threshold |    NULL   |
 *            | --incE       |  include hits <= this E-value threshold     |    0.01   |
 *            | --incT       |  include hits >= this bit score threshold   |    NULL   |
 *            | --incdomE    |  include domains <= this E-value threshold  |    0.01   |
 *            | --incdomT    |  include domains <= this score threshold    |    NULL   |
 *            | --cut_ga     |  model-specific thresholding using GA       |   FALSE   |
 *            | --cut_nc     |  model-specific thresholding using NC       |   FALSE   |
 *            | --cut_tc     |  model-specific thresholding using TC       |   FALSE   |
 *            | --max        |  turn all heuristic filters off             |   FALSE   |
 *            | --F1         |  Stage 1 (MSV) thresh: promote hits P <= F1 |    0.02   |
 *            | --F2         |  Stage 2 (Vit) thresh: promote hits P <= F2 |    1e-3   |
 *            | --F3         |  Stage 2 (Fwd) thresh: promote hits P <= F3 |    1e-5   |
 *            | --nobias     |  turn OFF composition bias filter HMM       |   FALSE   |
 *            | --nonull2    |  turn OFF biased comp score correction      |   FALSE   |
 *            | --seed       |  RNG seed (0=use arbitrary seed)            |      42   |
 *            | --acc        |  prefer accessions over names in output     |   FALSE   |
 *
 *            As a special case, if <go> is <NULL>, defaults are set as above.
 *            This shortcut is used in simplifying test programs and the like.
 *            
 * Returns:   ptr to new <P7_PIPELINE> object on success. Caller frees this
 *            with <p7_pipeline_Destroy()>.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_PIPELINE *
p7_pipeline_Create(ESL_GETOPTS *go, int M_hint, int L_hint, int do_longtargets, enum p7_pipemodes_e mode)
{
  P7_PIPELINE *pli  = NULL;
  int          status;

  /* Contract checks / arg validation */
  ESL_DASSERT1( (M_hint >= 1 && M_hint <= 100000) ); /* H3 design limit: M,L <= 100Kres */
  ESL_DASSERT1( (L_hint >= 1 && L_hint <= 100000) );

  ESL_ALLOC(pli, sizeof(P7_PIPELINE));
  pli->fx   = NULL;
  pli->cx   = NULL;
  pli->sm   = NULL;
  pli->sxf  = NULL;
  pli->sxb  = NULL;
  pli->sxd  = NULL;
  pli->sxx  = NULL;
  pli->tr   = NULL;
  pli->mt   = NULL;
  pli->n2sc = NULL;
  pli->wrk  = NULL;
  
  if ( (pli->fx  = p7_filtermx_Create  (M_hint))                                  == NULL) goto ERROR;
  if ( (pli->cx  = p7_checkptmx_Create (M_hint, L_hint, ESL_MBYTES(p7_RAMLIMIT))) == NULL) goto ERROR; /* p7_RAMLIMIT=256MB, in p7_config.h.in */
  if ( (pli->sm  = p7_sparsemask_Create(M_hint, L_hint))                          == NULL) goto ERROR;
  if ( (pli->sxf = p7_sparsemx_Create  (pli->sm))                                 == NULL) goto ERROR;
  if ( (pli->sxb = p7_sparsemx_Create  (pli->sm))                                 == NULL) goto ERROR;
  if ( (pli->sxd = p7_sparsemx_Create  (pli->sm))                                 == NULL) goto ERROR;
  if ( (pli->sxx = p7_sparsemx_Create  (pli->sm))                                 == NULL) goto ERROR;
  if ( (pli->tr  = p7_trace_CreateWithPP())                                       == NULL) goto ERROR;
  if ( (pli->mt  = p7_masstrace_Create (M_hint, L_hint))                          == NULL) goto ERROR;
  ESL_ALLOC(pli->n2sc, sizeof(float) * (L_hint+1));
  ESL_ALLOC(pli->wrk,  sizeof(float) * (M_hint+1));
  esl_vec_FSet(pli->n2sc, L_hint+1, 0.0f);

  /* Configure reporting thresholds */
  pli->by_E            = TRUE;
  pli->E               = (go ? esl_opt_GetReal(go, "-E") : 10.0);
  pli->T               = 0.0;
  pli->dom_by_E        = TRUE;
  pli->domE            = (go ? esl_opt_GetReal(go, "--domE") : 10.0);
  pli->domT            = 0.0;
  pli->use_bit_cutoffs = FALSE;
  if (go && esl_opt_IsOn(go, "-T"))     { pli->T    = esl_opt_GetReal(go, "-T");     pli->by_E     = FALSE; }
  if (go && esl_opt_IsOn(go, "--domT")) { pli->domT = esl_opt_GetReal(go, "--domT"); pli->dom_by_E = FALSE; }

  /* Configure inclusion thresholds */
  pli->inc_by_E           = TRUE;
  pli->incE               = (go ? esl_opt_GetReal(go, "--incE") : 0.01);
  pli->incT               = 0.0;
  pli->incdom_by_E        = TRUE;
  pli->incdomE            = (go ? esl_opt_GetReal(go, "--incdomE") : 0.01);
  pli->incdomT            = 0.0;
  if (go && esl_opt_IsOn(go, "--incT"))    { pli->incT    = esl_opt_GetReal(go, "--incT");    pli->inc_by_E    = FALSE; }
  if (go && esl_opt_IsOn(go, "--incdomT")) { pli->incdomT = esl_opt_GetReal(go, "--incdomT"); pli->incdom_by_E = FALSE; }

  /* Configure for one of the model-specific thresholding options */
  if (go && esl_opt_GetBoolean(go, "--cut_ga"))
    {
      pli->T        = pli->domT        = 0.0;
      pli->by_E     = pli->dom_by_E    = FALSE;
      pli->incT     = pli->incdomT     = 0.0;
      pli->inc_by_E = pli->incdom_by_E = FALSE;
      pli->use_bit_cutoffs = p7H_GA;
    }
  if (go && esl_opt_GetBoolean(go, "--cut_nc"))
    {
      pli->T        = pli->domT        = 0.0;
      pli->by_E     = pli->dom_by_E    = FALSE;
      pli->incT     = pli->incdomT     = 0.0;
      pli->inc_by_E = pli->incdom_by_E = FALSE;
      pli->use_bit_cutoffs = p7H_NC;
    }
  if (go && esl_opt_GetBoolean(go, "--cut_tc"))
    {
      pli->T        = pli->domT        = 0.0;
      pli->by_E     = pli->dom_by_E    = FALSE;
      pli->incT     = pli->incdomT     = 0.0;
      pli->inc_by_E = pli->incdom_by_E = FALSE;
      pli->use_bit_cutoffs = p7H_TC;
    }

  /* Configure search space sizes for E value calculations  */
  pli->Z       = pli->domZ       = 0.0;
  pli->Z_setby = pli->domZ_setby = p7_ZSETBY_NTARGETS;
  if (go && esl_opt_IsOn(go, "-Z")) 
    {
      pli->Z_setby = p7_ZSETBY_OPTION;
      pli->Z       = esl_opt_GetReal(go, "-Z");
    }
  if (go && esl_opt_IsOn(go, "--domZ")) 
    {
      pli->domZ_setby = p7_ZSETBY_OPTION;
      pli->domZ       = esl_opt_GetReal(go, "--domZ");
    }

  /* Configure acceleration pipeline thresholds */
  pli->do_max        = FALSE;
  pli->do_biasfilter = TRUE;
  pli->do_null2      = TRUE;
  pli->F1     = (go ? ESL_MIN(1.0, esl_opt_GetReal(go, "--F1")) : 0.02);
  pli->F2     = (go ? ESL_MIN(1.0, esl_opt_GetReal(go, "--F2")) : 1e-3);
  pli->F3     = (go ? ESL_MIN(1.0, esl_opt_GetReal(go, "--F3")) : 1e-5);
  if (do_longtargets) {
    pli->B1     = (go ? esl_opt_GetInteger(go, "--B1") : 100);
    pli->B2     = (go ? esl_opt_GetInteger(go, "--B2") : 240);
    pli->B3     = (go ? esl_opt_GetInteger(go, "--B3") : 1000);
  } else {
    pli->B1 = pli->B2 = pli->B3 = -1;
  }

  if (go && esl_opt_GetBoolean(go, "--max")) 
    {
      pli->do_max        = TRUE;
      pli->do_biasfilter = FALSE;

      pli->F2 = pli->F3 = 1.0;
      pli->F1 = (do_longtargets ? 0.3 : 1.0); // need to set some threshold for F1 even on long targets. Should this be tighter?
    }
  if (go && esl_opt_GetBoolean(go, "--nonull2")) pli->do_null2      = FALSE;
  if (go && esl_opt_GetBoolean(go, "--nobias"))  pli->do_biasfilter = FALSE;
  
  /* Zero the accounting counters */
  p7_pipeline_stats_Init(&(pli->stats));

  /* State config, flags. */
  pli->mode            = mode;
  pli->show_accessions = (go && esl_opt_GetBoolean(go, "--acc")   ? TRUE  : FALSE);
  pli->show_alignments = (go && esl_opt_GetBoolean(go, "--noali") ? FALSE : TRUE);
  pli->hfp             = NULL;

  /* Additional state config for nhmmer */
  pli->long_targets = do_longtargets;               /* TRUE | FALSE */
  pli->strand       = p7_STRAND_BOTH;               /* nhmmer's default */
  pli->W            = -1;                     /* p7_pipeline_NewModel() initiazes this, from profile */
  pli->block_length = p7_NHMMER_MAX_RESIDUE_COUNT;

  /* Diagnostic info */
  pli->errbuf[0]       = '\0';

  return pli;

 ERROR:
  p7_pipeline_Destroy(pli);
  return NULL;
}


/* Function:  p7_pipeline_Reuse()
 * Synopsis:  Reuse a pipeline for next target.
 *
 * Purpose:   Reuse <pli> for next target sequence (search mode)
 *            or model (scan mode). 
 *            
 *            Accounting statistics are not altered by this call;
 *            we accumulate statistics across many targets.
 *            
 *            May eventually need to distinguish from reusing pipeline
 *            for next query, but we're not really focused on multiquery
 *            use of hmmscan/hmmsearch/phmmer for the moment.
 */
int
p7_pipeline_Reuse(P7_PIPELINE *pli)
{
  p7_filtermx_Reuse  (pli->fx);
  p7_checkptmx_Reuse (pli->cx);
  p7_sparsemask_Reuse(pli->sm);
  p7_sparsemx_Reuse  (pli->sxf);
  p7_sparsemx_Reuse  (pli->sxb);
  p7_sparsemx_Reuse  (pli->sxd);
  p7_sparsemx_Reuse  (pli->sxx);
  p7_trace_Reuse     (pli->tr);
  p7_masstrace_Reuse (pli->mt);

  /* The rest of the state of the pipeline stays constant
   * as we move to the next target, including acct'ing stats.
   */
  return eslOK;
}

/* Function:  p7_pipeline_Destroy()
 * Synopsis:  Free a <P7_PIPELINE> object.
 *
 * Purpose:   Free a <P7_PIPELINE> object.
 */
void
p7_pipeline_Destroy(P7_PIPELINE *pli)
{
  if (pli) {
    p7_filtermx_Destroy  (pli->fx);
    p7_checkptmx_Destroy (pli->cx);
    p7_sparsemask_Destroy(pli->sm);
    p7_sparsemx_Destroy  (pli->sxf);
    p7_sparsemx_Destroy  (pli->sxb);
    p7_sparsemx_Destroy  (pli->sxd);
    p7_sparsemx_Destroy  (pli->sxx);
    p7_trace_Destroy     (pli->tr);
    p7_masstrace_Destroy (pli->mt);

    if (pli->n2sc) free(pli->n2sc);
    if (pli->wrk)  free(pli->wrk);
    free(pli);
  }
}
/*---------------- end, P7_PIPELINE object ----------------------*/


/*****************************************************************
 * 2. Setting pipeline for next target or model
 *****************************************************************/

/* Function:  p7_pipeline_NewModel()
 * Synopsis:  Prepare pipeline for a new model (target or query)
 *
 * Purpose:   Caller has a new model <om>. Prepare the pipeline <pli>
 *            to receive this model as either a query or a target.
 *
 *            If the "experimental" bias filter HMM is in use, this
 *            call resets it to use the new model's composition. This
 *            overwrites the bias filter HMM's expected length! You
 *            need to call <p7_bg_SetLength()> after a <NewModel()> call.
 *            (Failure to do this is bug #h85, 14 Dec 10.)
 *
 *            The pipeline may alter the null model <bg> in a model-specific
 *            way (if we're using a composition bias filter HMM in the
 *            pipeline).
 *
 * Returns:   <eslOK> on success.
 * 
 *            <eslEINVAL> if pipeline expects to be able to use a
 *            model's bit score thresholds, but this model does not
 *            have the appropriate ones set.
 *            
 *            <eslEMEM> if an allocation fails.
 */
int
p7_pipeline_NewModel(P7_PIPELINE *pli, const P7_OPROFILE *om, P7_BG *bg)
{
  int status = eslOK;

  pli->stats.nmodels++;
  pli->stats.nnodes += om->M;
  if (pli->Z_setby == p7_ZSETBY_NTARGETS && pli->mode == p7_SCAN_MODELS) pli->Z = pli->stats.nmodels;

  if (pli->do_biasfilter) p7_bg_SetFilter(bg, om->M, om->compo);

  if (pli->mode == p7_SEARCH_SEQS) 
    status = p7_pipeline_NewModelThresholds(pli, om);

  pli->W = om->max_length;
  return status;
}

/* Function:  p7_pipeline_NewModelThresholds()
 * Synopsis:  Set reporting and inclusion bit score thresholds on a new model.
 *
 * Purpose:   Set the bit score thresholds on a new model, if we're 
 *            using Pfam GA, TC, or NC cutoffs for reporting or
 *            inclusion.
 *            
 *            In a "search" pipeline, this only needs to be done once
 *            per query model, so <p7_pipeline_NewModelThresholds()> gets 
 *            called by <p7_pipeline_NewModel()>.
 *            
 *            In a "scan" pipeline, this needs to be called for each
 *            model, and it needs to be called after
 *            <p7_oprofile_ReadRest()>, because that's when the bit
 *            score thresholds get read.
 *
 * Returns:   <eslOK> on success. 
 *            
 *            <eslEINVAL> if pipeline expects to be able to use a
 *            model's bit score thresholds, but this model does not
 *            have the appropriate ones set.
 *
 * Xref:      Written to fix bug #h60.
 */
int
p7_pipeline_NewModelThresholds(P7_PIPELINE *pli, const P7_OPROFILE *om)
{

  if (pli->use_bit_cutoffs)
    {
      if (pli->use_bit_cutoffs == p7H_GA)
    {
      if (om->cutoff[p7_GA1] == p7_CUTOFF_UNSET) ESL_FAIL(eslEINVAL, pli->errbuf, "GA bit thresholds unavailable on model %s\n", om->name);
      pli->T    = pli->incT    = om->cutoff[p7_GA1];
      pli->domT = pli->incdomT = om->cutoff[p7_GA2];
    }
      else if  (pli->use_bit_cutoffs == p7H_TC)
    {
      if (om->cutoff[p7_TC1] == p7_CUTOFF_UNSET) ESL_FAIL(eslEINVAL, pli->errbuf, "TC bit thresholds unavailable on model %s\n", om->name);
      pli->T    = pli->incT    = om->cutoff[p7_TC1];
      pli->domT = pli->incdomT = om->cutoff[p7_TC2];
    }
      else if (pli->use_bit_cutoffs == p7H_NC)
    {
      if (om->cutoff[p7_NC1] == p7_CUTOFF_UNSET) ESL_FAIL(eslEINVAL, pli->errbuf, "NC bit thresholds unavailable on model %s\n", om->name);
      pli->T    = pli->incT    = om->cutoff[p7_NC1];
      pli->domT = pli->incdomT = om->cutoff[p7_NC2];
    }
    }

  return eslOK;
}


/* Function:  p7_pipeline_NewSeq()
 * Synopsis:  Prepare pipeline for a new sequence (target or query)
 *
 * Purpose:   Caller has a new sequence <sq>. Prepare the pipeline <pli>
 *            to receive this model as either a query or a target.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_pipeline_NewSeq(P7_PIPELINE *pli, const ESL_SQ *sq)
{
  if (!pli->long_targets) pli->stats.nseqs++; // if long_targets, sequence counting happens in the serial loop, which can track multiple windows for a single long sequence
  pli->stats.nres += sq->n;
  if (pli->Z_setby == p7_ZSETBY_NTARGETS && pli->mode == p7_SEARCH_SEQS) pli->Z = pli->stats.nseqs;
  return eslOK;
}
/*-------- end, setting pipeline for next target ----------------*/



/*****************************************************************
 * 3. Testing bitscore/E-value against reporting/inclusion threshold
 *****************************************************************/

/* Function:  p7_pipeline_TargetReportable
 * Synopsis:  Returns TRUE if target score meets reporting threshold.
 *
 * Purpose:   Returns <TRUE> if the bit score <score> and/or 
 *            log P-value <lnP> meet per-target reporting thresholds 
 *            for the processing pipeline.
 */
int
p7_pipeline_TargetReportable(P7_PIPELINE *pli, float score, double lnP)
{
  if (pli->by_E )
    {
      if ( !pli->long_targets  && exp(lnP) * pli->Z <= pli->E) return TRUE;
      if (  pli->long_targets  && exp(lnP) <= pli->E)          return TRUE; // database size is already built into the Pval if pli->targetlength == p7_TARGET_LONG
    }
  else if (! pli->by_E   && score >= pli->T) return TRUE;
  return FALSE;
}

/* Function:  p7_pipeline_DomainReportable
 * Synopsis:  Returns TRUE if domain score meets reporting threshold. 
 *
 * Purpose:   Returns <TRUE> if the bit score <score> and/or 
 *            log P-value <lnP> meet per-domain reporting thresholds 
 *            for the processing pipeline.
 */
int
p7_pipeline_DomainReportable(P7_PIPELINE *pli, float dom_score, double lnP)
{
  if ( pli->dom_by_E )
    {
      if ( !pli->long_targets  &&  exp(lnP) * pli->domZ <= pli->domE) return TRUE;
      if (  pli->long_targets  &&  exp(lnP) <= pli->domE) return TRUE;
    }
  else if (! pli->dom_by_E   && dom_score        >= pli->domT) return TRUE;
  return FALSE;
}

/* Function:  p7_pipeline_TargetIncludable()
 * Synopsis:  Returns TRUE if target score meets inclusion threshold.
 */
int
p7_pipeline_TargetIncludable(P7_PIPELINE *pli, float score, double lnP)
{
  if (pli->inc_by_E )
    {
      if ( !pli->long_targets && exp(lnP) * pli->Z <= pli->incE) return TRUE;
      if (  pli->long_targets && exp(lnP) <= pli->incE) return TRUE;
    }
  else if (! pli->inc_by_E   && score         >= pli->incT) return TRUE;
  return FALSE;
}

/* Function:  p7_pipeline_DomainIncludable()
 * Synopsis:  Returns TRUE if domain score meets inclusion threshold.
 */
int
p7_pipeline_DomainIncludable(P7_PIPELINE *pli, float dom_score, double lnP)
{
  if      (  pli->incdom_by_E   && exp(lnP) * pli->domZ <= pli->incdomE) return TRUE;
  else if (! pli->incdom_by_E   && dom_score        >= pli->incdomT) return TRUE;
  else return FALSE;
}
/*-------- end, testing score/E-value against thresholds  -------*/


/*****************************************************************
 * 4. Statistics output from a completed pipeline
 *****************************************************************/

/* Function:  p7_pipeline_stats_Init()
 * Synopsis:  Zero the accounting count accumulators.
 */
int
p7_pipeline_stats_Init(P7_PIPELINE_STATS *stats)
{
  stats->nmodels       = 0;
  stats->nseqs         = 0;
  stats->nres          = 0;
  stats->nnodes        = 0;

  stats->n_past_msv    = 0;
  stats->n_past_bias   = 0;
  stats->n_past_vit    = 0;
  stats->n_past_fwd    = 0;

  stats->n_output      = 0;
  stats->pos_past_msv  = 0;
  stats->pos_past_bias = 0;
  stats->pos_past_vit  = 0;
  stats->pos_past_fwd  = 0;
  stats->pos_output    = 0;
  return eslOK;
}

/* Function:  p7_pipeline_stats_Merge()
 * Synopsis:  Merge the pipeline statistics (among workers/threads)
 *
 * Purpose:   Merge the statistics of a subsidiary pipeline <stats> (from a worker thread
 *            or MPI process, say) into a master pipeline <p1>, prior to output
 *            of summary statistics.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_pipeline_stats_Merge(P7_PIPELINE *p1, const P7_PIPELINE_STATS *stats)
{
  /* if we are searching a sequence database, we need to keep track of the
   * number of sequences and residues processed.
   */
  if (p1->mode == p7_SEARCH_SEQS)
    {
      p1->stats.nseqs   += stats->nseqs;
      p1->stats.nres    += stats->nres;
    }
  else
    {
      p1->stats.nmodels += stats->nmodels;
      p1->stats.nnodes  += stats->nnodes;
    }

  p1->stats.n_past_msv  += stats->n_past_msv;
  p1->stats.n_past_bias += stats->n_past_bias;
  p1->stats.n_past_vit  += stats->n_past_vit;
  p1->stats.n_past_fwd  += stats->n_past_fwd;
  p1->stats.n_output    += stats->n_output;

  p1->stats.pos_past_msv  += stats->pos_past_msv;
  p1->stats.pos_past_bias += stats->pos_past_bias;
  p1->stats.pos_past_vit  += stats->pos_past_vit;
  p1->stats.pos_past_fwd  += stats->pos_past_fwd;
  p1->stats.pos_output    += stats->pos_output;

  if (p1->Z_setby == p7_ZSETBY_NTARGETS)
    {
      p1->Z += (p1->mode == p7_SCAN_MODELS) ? stats->nmodels : stats->nseqs;
    }
  return eslOK;
}

/* Function:  p7_pipeline_WriteStats()
 * Synopsis:  Final statistics output from a processing pipeline.
 *
 * Purpose:   Print a standardized report of the internal statistics of
 *            a finished processing pipeline <pli> to stream <ofp>.
 *            
 *            If stopped, non-<NULL> stopwatch <w> is provided for a
 *            stopwatch that was timing the pipeline, then the report
 *            also includes timing information.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_pipeline_WriteStats(FILE *ofp, P7_PIPELINE *pli, ESL_STOPWATCH *w)
{
  double ntargets; 

  fprintf(ofp, "Internal pipeline statistics summary:\n");
  fprintf(ofp, "-------------------------------------\n");
  if (pli->mode == p7_SEARCH_SEQS) {
    fprintf(ofp, "Query model(s):              %15" PRId64 "  (%" PRId64 " nodes)\n",              pli->stats.nmodels, pli->stats.nnodes);
    fprintf(ofp, "Target sequences:            %15" PRId64 "  (%" PRId64 " residues searched)\n",  pli->stats.nseqs,   pli->stats.nres);
    ntargets = pli->stats.nseqs;
  } else {
    fprintf(ofp, "Query sequence(s):           %15" PRId64 "  (%" PRId64 " residues searched)\n",  pli->stats.nseqs,   pli->stats.nres);
    fprintf(ofp, "Target model(s):             %15" PRId64 "  (%" PRId64 " nodes)\n",              pli->stats.nmodels, pli->stats.nnodes);
    ntargets = pli->stats.nmodels;
  }

  if (pli->long_targets) {
      fprintf(ofp, "Residues passing MSV filter:   %15" PRId64 "  (%.3g); expected (%.3g)\n",
        pli->stats.pos_past_msv,
        (double)pli->stats.pos_past_msv / (pli->stats.nres*pli->stats.nmodels) ,
        pli->F1);

      fprintf(ofp, "Residues passing bias filter:  %15" PRId64 "  (%.3g); expected (%.3g)\n",
        pli->stats.pos_past_bias,
        (double)pli->stats.pos_past_bias / (pli->stats.nres*pli->stats.nmodels) ,
        pli->F1);

      fprintf(ofp, "Residues passing Vit filter:   %15" PRId64 "  (%.3g); expected (%.3g)\n",
        pli->stats.pos_past_vit,
        (double)pli->stats.pos_past_vit / (pli->stats.nres*pli->stats.nmodels) ,
        pli->F2);

      fprintf(ofp, "Residues passing Fwd filter:   %15" PRId64 "  (%.3g); expected (%.3g)\n",
        pli->stats.pos_past_fwd,
        (double)pli->stats.pos_past_fwd / (pli->stats.nres*pli->stats.nmodels) ,
        pli->F3);

      fprintf(ofp, "Total number of hits:          %15d  (%.3g)\n",
          (int)pli->stats.n_output,
          (double)pli->stats.pos_output / (pli->stats.nres*pli->stats.nmodels) );

  } else { // typical case output

      fprintf(ofp, "Passed MSV filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",
          pli->stats.n_past_msv,
          (double) pli->stats.n_past_msv / ntargets,
          pli->F1 * ntargets,
          pli->F1);

      fprintf(ofp, "Passed bias filter:          %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",
          pli->stats.n_past_bias,
          (double) pli->stats.n_past_bias / ntargets,
          pli->F1 * ntargets,
          pli->F1);

      fprintf(ofp, "Passed Vit filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",
          pli->stats.n_past_vit,
          (double) pli->stats.n_past_vit / ntargets,
          pli->F2 * ntargets,
          pli->F2);

      fprintf(ofp, "Passed Fwd filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",
          pli->stats.n_past_fwd,
          (double) pli->stats.n_past_fwd / ntargets,
          pli->F3 * ntargets,
          pli->F3);

      fprintf(ofp, "Initial search space (Z):    %15.0f  %s\n", pli->Z,    pli->Z_setby    == p7_ZSETBY_OPTION ? "[as set by --Z on cmdline]"    : "[actual number of targets]");
      fprintf(ofp, "Domain search space  (domZ): %15.0f  %s\n", pli->domZ, pli->domZ_setby == p7_ZSETBY_OPTION ? "[as set by --domZ on cmdline]" : "[number of targets reported over threshold]");
  }

  if (w != NULL) {
    esl_stopwatch_Display(ofp, w, "# CPU time: ");
    fprintf(ofp, "# Mc/sec: %.2f\n", 
        (double) pli->stats.nres * (double) pli->stats.nnodes / (w->elapsed * 1.0e6));
  }

  return eslOK;
}
/*--------------- end, pipeline statistics ----------------------*/



/*****************************************************************
 * 5. The main pipeline call: p7_Pipeline()
 *****************************************************************/


/* Function:  p7_Pipeline()
 * Synopsis:  HMMER3's accelerated seq/profile comparison pipeline.
 *
 * Purpose:   Run H3's accelerated pipeline to compare profile <om>
 *            against sequence <sq>. If a significant hit is found,
 *            information about it is added to the <hitlist>. The pipeline 
 *            accumulates beancounting information about how many comparisons
 *            flow through the pipeline while it's active.
 *            
 * Returns:   <eslOK> on success. If a significant hit is obtained,
 *            its information is added to the growing <hitlist>. 
 *            
 *            <eslEINVAL> if (in a scan pipeline) we're supposed to
 *            set GA/TC/NC bit score thresholds but the model doesn't
 *            have any.
 *            
 *            <eslERANGE> on numerical overflow errors in the
 *            optimized vector implementations; particularly in
 *            posterior decoding. I don't believe this is possible for
 *            multihit local models, but I'm set up to catch it
 *            anyway. We may emit a warning to the user, but cleanly
 *            skip the problematic sequence and continue.
 *
 * Throws:    <eslEMEM> on allocation failure.
 * 
 * Future:    Right now it needs both <gm> and <om> because the filters
 *            use <om> and sparse DP uses <gm>. Not a problem in hmmsearch,
 *            which typically has both; more of a problem for hmmscan,
 *            which reads <om>. In the future, we need to optimize a bit
 *            more; perhaps we'll have P7_MODEL, holding annotation,
 *            hmm, profile params, and striped vector params for MSV,
 *            Vit, and Fwd/Back, with the ability to create vector params
 *            from profile and vice versa, on demand.
 */
int
p7_Pipeline(P7_PIPELINE *pli, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, P7_TOPHITS *hitlist)
{
  P7_DOMAIN       *dcl     = NULL;     /* array of domain data structures <0..ndom-1> */
  P7_HIT          *hit     = NULL;     /* ptr to the current hit output data          */
  float            usc, vitsc, fwdsc;  /* DP scores                                   */
  float            filtersc;           /* HMM null filter score                       */
  float            nullsc;             /* null model score                            */
  float            seqbias;  
  float            seq_score;          /* the corrected per-seq bit score */
  float            sum_score;           /* the corrected reconstruction score for the seq */
  float            pre_score, pre2_score; /* uncorrected bit scores for seq */
  double           P;                /* P-value of a hit */
  double           lnP;              /* log P-value of a hit */
  int              Ld;               /* # of residues in envelopes */
  int              d,z,i;
  float            null2[p7_MAXCODE];
  int              noverlaps;
  int              last_ibe;
  int              best_d;
  int              status;
  
  if (sq->n == 0) return eslOK;    /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */

  /* Base null model score (we could calculate this in NewSeq(), for a scan pipeline) */
  p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);

  /* First level filter: the SSV and MSV filters */
  p7_MSVFilter(sq->dsq, sq->n, om, pli->fx, &usc);
  seq_score = (usc - nullsc) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  if (P > pli->F1) return eslOK;
  pli->stats.n_past_msv++;

  /* biased composition HMM filtering */
  if (pli->do_biasfilter)
    {
      p7_bg_FilterScore(bg, sq->dsq, sq->n, &filtersc);
      seq_score = (usc - filtersc) / eslCONST_LOG2;
      P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      if (P > pli->F1) return eslOK;
    }
  else filtersc = nullsc;
  pli->stats.n_past_bias++;

  /* In scan mode, if it passes the MSV filter, read the rest of the profile */
  if (pli->mode == p7_SCAN_MODELS)
    {
      if (pli->hfp) p7_oprofile_ReadRest(pli->hfp, om);
      p7_oprofile_ReconfigRestLength(om, sq->n);
      if ((status = p7_pipeline_NewModelThresholds(pli, om)) != eslOK) return status; /* pli->errbuf has err msg set */
    }

  /* Second level filter: ViterbiFilter(), multihit with <om> */
  if (P > pli->F2)
    {
      p7_ViterbiFilter(sq->dsq, sq->n, om, pli->fx, &vitsc);  
      seq_score = (vitsc-filtersc) / eslCONST_LOG2;
      P  = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
      if (P > pli->F2) return eslOK;
    }
  pli->stats.n_past_vit++;

  /* Checkpointed Forward. Check score as a filter step before proceeding to backwards/decoding. */
  p7_ForwardFilter(sq->dsq, sq->n, om, pli->cx, &fwdsc);
  seq_score = (fwdsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  if (P > pli->F3) return eslOK;
  pli->stats.n_past_fwd++;

  /* ok, it's for real; passes vectorized local scoring filters.
   * Finish with Backwards and Decoding, defining a sparse mask in <sm>.
   */
  p7_BackwardFilter(sq->dsq, sq->n, om, pli->cx, pli->sm, p7_SPARSEMASK_THRESH_DEFAULT);

  /* FIXME 3.1
   * hmmscan needs <gm>; read from the hmm file (.h3m) on disk
   * expect that this is slow and inefficient; come back and refactor later.
   * (in hmmpgmd, hfp is NULL; all our profiles are already in a cache)
   */
  if (pli->mode == p7_SCAN_MODELS && pli->hfp)
    workaround_get_profile(pli, om, bg, &gm);

  /* Now we can hand it over to sparse DP, with the full glocal/local model */
  p7_SparseViterbi (sq->dsq, sq->n, gm, pli->sm,  pli->sxx, pli->tr, &vitsc);
  p7_SparseForward (sq->dsq, sq->n, gm, pli->sm,  pli->sxf,          &fwdsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, pli->sm,  pli->sxb,          NULL);
  p7_SparseDecoding(sq->dsq, sq->n, gm, pli->sxf, pli->sxb, pli->sxd);
  p7_sparsemx_TracePostprobs(pli->sxd, pli->tr); /* annotate the trace. */
  p7_trace_Index(pli->tr);             /* index domains in the trace */
  p7_sparsemx_Reuse(pli->sxx);                   /* reuse the Viterbi matrix for mass tracing, coming up next */
  

  dcl = p7_domain_Create(pli->tr->ndom);
  //ESL_ALLOC(dcl, sizeof(P7_DOMAIN) * pli->tr->ndom);
  ESL_REALLOC(pli->wrk,  sizeof(float) * (gm->M+1));
  ESL_REALLOC(pli->n2sc, sizeof(float) * (sq->n+1));
  esl_vec_FSet(pli->n2sc, sq->n+1, 0.0f);

  noverlaps  = 0;
  last_ibe   = -1;
  best_d     = 0;
  for (d = 0; d < pli->tr->ndom; d++)
    {

      /* Determine envelope coords by mass trace. */
      p7_SparseMasstrace(sq->dsq, sq->n, gm, pli->sxf, pli->sxb, pli->tr, pli->tr->anch[d], p7_MASSTRACE_THRESH_DEFAULT, pli->sxx, pli->mt, 
                         &(dcl[d].iae), &(dcl[d].ibe), &(dcl[d].kae), &(dcl[d].kbe));
      p7_sparsemx_Reuse (pli->sxx);
      p7_masstrace_Reuse(pli->mt);

      /* Keep track of overlaps */
      if (dcl[d].iae <= last_ibe) noverlaps++;
      last_ibe = dcl[d].ibe;

      /* Transfer alignment coords from trace. [Do we need to do this?] */
      dcl[d].ia = pli->tr->sqfrom[d];  
      dcl[d].ib = pli->tr->sqto[d];  
      dcl[d].ka = pli->tr->hmmfrom[d];  
      dcl[d].kb = pli->tr->hmmto[d];  

      /* Determine envelope score. [We have an approximation available, but we determined it rarely holds */
      p7_SparseEnvscore(sq->dsq, sq->n, gm, dcl[d].iae, dcl[d].ibe, dcl[d].kae, dcl[d].kbe, pli->sm, pli->sxx, &(dcl[d].envsc));
      p7_sparsemx_Reuse(pli->sxx);

      /* Determine null2 correction 
       * Upon return, null2[x] = \log f'(x)/f(x).
       * Then we record that in n2sc[iae..ibe]: seqbias correction mustn't overcount overlapped envelopes
       */
      if (pli->do_null2) 
      {
        p7_sparse_Null2ByExpectation(gm, pli->sxd, dcl[d].iae, dcl[d].ibe, dcl[d].kae, dcl[d].kbe, pli->wrk, null2);
        dcl[d].domcorrection = 0.;
        for (i = dcl[d].iae; i <= dcl[d].ibe; i++)
          {
            dcl[d].domcorrection += null2[ sq->dsq[i] ];
            pli->n2sc[i]          = null2[ sq->dsq[i] ];
          }
        dcl[d].dombias = p7_FLogsum(0., log(bg->omega) + dcl[d].domcorrection);
      }
      else 
      {
        dcl[d].domcorrection = 0.;
        dcl[d].dombias       = 0.;
      }

      /* alignment "accuracy score": expected number of correctly aligned positions */
      dcl[d].oasc = 0.;
      for (z = pli->tr->tfrom[d]; z <= pli->tr->tto[d]; z++)
        if (pli->tr->i[z]) dcl[d].oasc += pli->tr->pp[z];

      /* domain bit score = null2-corrected bit score of the sequence if this were the only envelope in it. */
      dcl[d].bitscore = (dcl[d].envsc - (nullsc + dcl[d].dombias)) / eslCONST_LOG2;
      dcl[d].lnP      = esl_exp_logsurv( dcl[d].bitscore, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

      if (dcl[d].bitscore > dcl[best_d].bitscore) best_d = d;

      /* Viterbi alignment of the domain */
      dcl[d].ad = p7_alidisplay_Create(pli->tr, d, om, sq);

      /* We're initiazing a P7_DOMAIN structure in dcl[d] by hand, without a Create().
       * We're responsible for initiazing all elements of this structure.
       */
      dcl[d].is_reported = FALSE; /* will get set later by p7_tophits_Threshold() */
      dcl[d].is_included = FALSE; /* ditto */
    }
  
  /* Calculate the null2-corrected per-seq score */
  seqbias = (pli->do_null2 ? p7_FLogsum(0.0, log(bg->omega) + esl_vec_FSum(pli->n2sc, sq->n+1)) : 0.0);
  pre_score =  (fwdsc - nullsc) / eslCONST_LOG2;              /* BITS */
  seq_score =  (fwdsc - (nullsc + seqbias)) / eslCONST_LOG2;  /* BITS */

  
  /* Calculate the "reconstruction score": estimated
   * per-sequence score as sum of individual domains,
   * discounting domains that don't contribute positive
   * score after being null2-corrected.
   */
  if (pli->do_null2) 
    {
      float baseline = sq->n * gm->xsc[p7P_C][p7P_LOOP] + gm->xsc[p7P_C][p7P_MOVE];
      sum_score = 0.0f;
      seqbias   = 0.0f;
      Ld        = 0;
      for (d = 0; d < pli->tr->ndom; d++)
      {
        if (dcl[d].envsc - dcl[d].dombias > baseline)
          {
            sum_score += dcl[d].envsc - (sq->n - (dcl[d].ibe - dcl[d].iae + 1)) * gm->xsc[p7P_C][p7P_LOOP] - gm->xsc[p7P_C][p7P_MOVE];
            Ld        += dcl[d].ibe - dcl[d].iae + 1;
            seqbias   += dcl[d].domcorrection; /* NATS */
          }
      }
      seqbias    = p7_FLogsum(0.0, log(bg->omega) + seqbias);      /* NATS */
      sum_score += (sq->n - Ld) * gm->xsc[p7P_C][p7P_LOOP] + gm->xsc[p7P_C][p7P_MOVE];
      pre2_score = (sum_score - nullsc) / eslCONST_LOG2;           /* BITS */
      sum_score  = (sum_score - (nullsc+seqbias)) / eslCONST_LOG2; /* BITS */
    }
  else 
    {
      pre2_score = pre_score;
      sum_score  = seq_score;
    }

  /* Let sum_score override the seq_score when it's better, and it includes at least 1 domain */
  if (Ld > 0 && sum_score > seq_score)
    {
      seq_score = sum_score;
      pre_score = pre2_score;
    }

  /* Apply thresholding and determine whether to put this
   * target into the hit list. E-value thresholding may
   * only be a lower bound for now, so this list may be longer
   * than eventually reported.
   */
  lnP =  esl_exp_logsurv (seq_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
  if (! p7_pipeline_TargetReportable(pli, seq_score, lnP))
    {
      p7_domain_Destroy(dcl, pli->tr->ndom);
    }
  else
    {
      p7_tophits_CreateNextHit(hitlist, &hit);
      if (pli->mode == p7_SEARCH_SEQS) {
        if (                       (status  = esl_strdup(sq->name, -1, &(hit->name)))  != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
        if (sq->acc[0]  != '\0' && (status  = esl_strdup(sq->acc,  -1, &(hit->acc)))   != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
        if (sq->desc[0] != '\0' && (status  = esl_strdup(sq->desc, -1, &(hit->desc)))  != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
      } else {
        if ((status  = esl_strdup(om->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure");
        if ((status  = esl_strdup(om->acc,  -1, &(hit->acc)))   != eslOK) esl_fatal("allocation failure");
        if ((status  = esl_strdup(om->desc, -1, &(hit->desc)))  != eslOK) esl_fatal("allocation failure");
      } 
      hit->window_length = 0;
      hit->sortkey       = pli->inc_by_E ? -lnP : seq_score; /* per-seq output sorts on bit score if inclusion is by score  */
      
      hit->score      = seq_score;
      hit->pre_score  = pre_score;
      hit->sum_score  = sum_score;

      hit->lnP        = lnP;
      hit->pre_lnP    = esl_exp_logsurv (hit->pre_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
      hit->sum_lnP    = esl_exp_logsurv (hit->sum_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

      p7_sparsemx_ExpectedDomains(pli->sxd, 1, sq->n, &(hit->nexpected));
      hit->noverlaps  = noverlaps;
      hit->ndom       = pli->tr->ndom;

      hit->flags       = 0;
      hit->nreported   = 0;
      hit->nincluded   = 0;
      hit->best_domain = best_d;

      hit->seqidx       = -1; /* nhmmer/longtarget only */
      hit->subseq_start = -1; /* nhmmer/longtarget only */
      
      hit->dcl        = dcl;
      dcl             = NULL;

      /* If we're using model-specific bit score thresholds (GA | TC |
       * NC) and we're in an hmmscan pipeline (mode = p7_SCAN_MODELS),
       * then we *must* apply those reporting or inclusion thresholds
       * now, because this model is about to go away; we won't have
       * its thresholds after all targets have been processed.
       * 
       * If we're using E-value thresholds and we don't know the
       * search space size (Z_setby or domZ_setby =
       * p7_ZSETBY_NTARGETS), we *cannot* apply those thresholds now,
       * and we *must* wait until all targets have been processed
       * (see p7_tophits_Threshold()).
       * 
       * For any other thresholding, it doesn't matter whether we do
       * it here (model-specifically) or at the end (in
       * p7_tophits_Threshold()). 
       * 
       * What we actually do, then, is to set the flags if we're using
       * model-specific score thresholds (regardless of whether we're
       * in a scan or a search pipeline); otherwise we leave it to 
       * p7_tophits_Threshold(). p7_tophits_Threshold() is always
       * responsible for *counting* the reported, included sequences.
       * 
       * [xref J5/92]
       */
      if (pli->use_bit_cutoffs)
      {
        if (p7_pipeline_TargetReportable(pli, hit->score, hit->lnP))
        {
          hit->flags |= p7_IS_REPORTED;
          if (p7_pipeline_TargetIncludable(pli, hit->score, hit->lnP))
            hit->flags |= p7_IS_INCLUDED;
        }

        for (d = 0; d < hit->ndom; d++)
        {
          if (p7_pipeline_DomainReportable(pli, hit->dcl[d].bitscore, hit->dcl[d].lnP))
          {
            hit->dcl[d].is_reported = TRUE;
            if (p7_pipeline_DomainIncludable(pli, hit->dcl[d].bitscore, hit->dcl[d].lnP))
              hit->dcl[d].is_included = TRUE;
          }
        }
      }
    }
  if (pli->mode == p7_SCAN_MODELS && pli->hfp) p7_profile_Destroy(gm); /* DON'T free it if we're in hmmpgmd, with all models cached */
  return eslOK;

 ERROR:
  if (dcl) free(dcl);
  return status;
}


/* Temporary workaround for a problem in the hmmscan version of the pipeline.
 * hmmscan reads vectorized profile in two pieces from .h3f, .h3p files.
 * In H3.0 we only needed vectorized profile in the pipeline.
 * In H3.1 we need both <om> and <gm>.
 * hmmscan doesn't have <gm>.
 * Options include:
 *   1. Convert <om> to <gm>.
 *      - We'd have to recalculate the BGMk and DGkE entry/exit wing retractions
 *        from the profile; normally modelconfig does this from the HMM.
 *      - This might be slow.
 *   2. Store the profile on disk too.
 *      - I like this option, but I'd rather do it together with some
 *        other reengineering. We have a lot of redundancy in HMM, PROFILE,
 *        and OPROFILE, particularly in annotation. Should consider
 *        having a P7_MODEL container around P7_MODELINFO (annotation),
 *        P7_HMM, P7_PROFILE, vector MSV, vector VF, vector FB subparts, with ways to
 *        convert amongst HMM, PROFILE, MSV, VF, and FB sections, or read
 *        them from disk. Support delayed read of annotation including
 *        name/desc, allow HMMs (and target seqs?) to be numbered for
 *        efficiency. Measure memory footprints and timing of read, MPI 
 *        transmit, conversion.
 *   3. As it happens, we do have the HMM on disk, in the .h3m file.
 *      We can read it, and convert it to a profile.
 *      This might be slow - especially since we need to alloc/dealloc
 *      the HMM and profile in the pipeline, rather than reusing them.
 *      
 * (3) is the fastest option to implement,
 * and right now the pressure is to get 3.1 compiling and running asap;
 * we can optimize/polish later from a baseline implementation.
 */
static int
workaround_get_profile(P7_PIPELINE *pli, const P7_OPROFILE *om, const P7_BG *bg, P7_PROFILE **ret_gm)
{
  
  P7_HMM     *hmm = NULL;
  P7_PROFILE *gm  = NULL;
  int         status;

  if ( (status = p7_hmmfile_Position(pli->hfp, om->offs[p7_MOFFSET])) != eslOK) goto ERROR; /* {eslESYS | eslEINVAL} */
  if ( (status = p7_hmmfile_Read(pli->hfp, (ESL_ALPHABET **) &(om->abc), &hmm)) != eslOK) goto ERROR; /* eslEOF | eslEINCOMPAT; {eslEMEM | eslESYS} */
  /* the ESL_ALPHABET ** cast was to get rid of a const; safe, but ugly */

  if ( (    gm = p7_profile_Create(hmm->M, om->abc))                  == NULL)  { status = eslEMEM; goto ERROR; }
  if ( (status = p7_profile_Config(gm, hmm, bg))                      != eslOK) goto ERROR;
  if ( (status = p7_profile_SetLength(gm, om->L))                     != eslOK) goto ERROR;
  *ret_gm = gm;

  p7_hmm_Destroy(hmm);
  return eslOK;
  
 ERROR:
  if (hmm) p7_hmm_Destroy(hmm);
  if (gm)  p7_profile_Destroy(gm);
  return status;
}



/*****************************************************************
 * 6. The pipeline, speciazed for nhmmer/longtarget
 *****************************************************************/
/* p7_pipeline_ExtendAndMergeWindows()
 * Synopsis:  Turns a list of ssv diagonals into windows, and merges
 *            overlapping windows.
 *
 * Purpose:   Accepts a <windowlist> of SSV diagonals, extends those
 *            to windows based on a combination of the max_length
 *            value from <om> and the prefix and suffix lengths stored
 *            in <data>, then merges (in place) windows that overlap
 *            by more than <pct_overlap> percent, ensuring that windows
 *            stay within the bounds of 1..<L>.
 *
 * Returns:   <eslOK>
 */
static int
p7_pipeline_ExtendAndMergeWindows (P7_OPROFILE *om, const P7_SCOREDATA *data, P7_HMM_WINDOWLIST *windowlist, int L, float pct_overlap)
{
  int i;
  P7_HMM_WINDOW   *prev_window = NULL;
  P7_HMM_WINDOW   *curr_window = NULL;
  int              window_start;
  int              window_end;
  int new_hit_cnt = 0;

  if (windowlist->count == 0) return eslOK;

  /* extend windows */
  for (i=0; i<windowlist->count; i++) {
    curr_window = windowlist->windows+i;

    // the 0.1 multiplier provides for a small buffer in excess of the predefined prefix/suffix lengths - one proportional to max_length
    window_start = ESL_MAX( 1,   curr_window->n - ( om->max_length * (0.1 + data->prefix_lengths[curr_window->k - curr_window->length + 1]  )) ) ;
    window_end   = ESL_MIN( L ,  curr_window->n + curr_window->length + (om->max_length * (0.1 + data->suffix_lengths[curr_window->k] ) ) )   ;

    curr_window->n = window_start;
    curr_window->length = window_end - window_start + 1;
  }


  /* merge overlapping windows, compressing list in place. */
  for (i=1; i<windowlist->count; i++) {
    prev_window = windowlist->windows+new_hit_cnt;
    curr_window = windowlist->windows+i;

    window_start = ESL_MAX(prev_window->n, curr_window->n);
    window_end   = ESL_MIN(prev_window->n+prev_window->length-1, curr_window->n+curr_window->length-1);

    if (  (float)(window_end-window_start+1)/ESL_MIN(prev_window->length, curr_window->length) > pct_overlap ) {
      //merge windows
      if (  curr_window->n + curr_window->length >  prev_window->n + prev_window->length )
        prev_window->length = curr_window->n + curr_window->length - prev_window->n;  //+1, -1 factored out

    } else {
      new_hit_cnt++;
      windowlist->windows[new_hit_cnt] = windowlist->windows[i];
    }
  }
  windowlist->count = new_hit_cnt+1;

  /* ensure that window start and end are within target sequence bounds */
  if ( windowlist->windows[0].n  <  1)
    windowlist->windows[0].n =  1;

  if ( windowlist->windows[windowlist->count-1].n + windowlist->windows[windowlist->count-1].length - 1  >  L)
    windowlist->windows[windowlist->count-1].length =  L - windowlist->windows[windowlist->count-1].n + 1;

  return eslOK;
}


/* Function:  p7_pipeline_computeAliScores()
 * Synopsis:  Compute per-position scores for the alignment for a domain
 *
 * Purpose:   Compute per-position (Viterbi) scores for the alignment for a domain,
 *            for the purpose of optionally printing these scores out in association
 *            with each alignment. Such scores can, for example, be used to detangle
 *            overlapping alignments (from different models)
 *
 * Args:      dom             - domain with the alignment for which we wish to compute scores
 *            seq             - sequence in which domain resides
 *            data            - contains model's emission and transition values in unstriped form
 *            K               - alphabet size
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
static int
p7_pipeline_computeAliScores(P7_DOMAIN *dom, ESL_DSQ *seq, const P7_SCOREDATA *data, int K)
{
  int status;
  int i, j, k;
  float sc;

  //Compute score contribution of each position in the alignment to the overall Viterbi score
  ESL_ALLOC( dom->scores_per_pos, sizeof(float) * dom->ad->N );
  for (i=0; i<dom->ad->N; i++)  dom->scores_per_pos[i] = 0.0;

  i = dom->ia - 1;        //sequence position
  j = dom->ad->hmmfrom - 1; //model position
  k = 0;
  while ( k<dom->ad->N) {
    if (dom->ad->model[k] != '.' && dom->ad->aseq[k] != '-') { //match
      i++;  j++;
      // Including the MM cost is a hack. The cost of getting to/from this match
      // state does matter, but an IM or DM transition would improperly deflate
      // the score of this column, so just give MM. That amount is offset out of
      // the score shown for preceding indels
      dom->scores_per_pos[k] = data->fwd_scores[K * j + seq[i]]
                             +  (j==1 ? 0 : log(data->fwd_transitions[p7O_MM][j]) );
      k++;
    } else if (dom->ad->model[k] == '.' ) { // insert
      //spin through the insert, accumulating cost;  only assign to final column in gap
      dom->scores_per_pos[k] = -eslINFINITY;

      sc = log(data->fwd_transitions[p7O_MI][j]);
      i++; k++;
      while (k<dom->ad->N && dom->ad->model[k] == '.') { //extend insert
        dom->scores_per_pos[k] = -eslINFINITY;
        sc += log(data->fwd_transitions[p7O_II][j]);
        i++; k++;
      }
      sc += log(data->fwd_transitions[p7O_IM][j+1]) - log(data->fwd_transitions[p7O_MM][j+1]);
      dom->scores_per_pos[k-1] = sc;

    } else if (dom->ad->aseq[k] == '-' ) { // delete
      dom->scores_per_pos[k] = -eslINFINITY;
      sc = log(data->fwd_transitions[p7O_MD][j]);
      j++; k++;
      while (k<dom->ad->N && dom->ad->aseq[k] == '-')  { //extend delete
        dom->scores_per_pos[k] = -eslINFINITY;
        sc += log(data->fwd_transitions[p7O_DD][j]);
        j++; k++;
      }
      sc += log(data->fwd_transitions[p7O_DM][j+1]) - log(data->fwd_transitions[p7O_MM][j+1]);
      dom->scores_per_pos[k-1] = sc;
    }
  }

  return eslOK;

ERROR:
  return eslEMEM;

}


/* parameterize_gm()
 *
 * Compute new (temporary) background priors based on a sequence
 * window, and, for a model that will be used narrowly within
 * the post-fwd portion of the longtarget pipeline, set match
 * state emission log-odds scores accordingly.
 *
 * Given the frequency distribution of a block of sequence
 * (<dsq>, of length <sq_len>), and a prior with which that
 * distribution can be mixed (<bg>), compute a new background
 * in the pre-allocated <bg_tmp> (used as a temporary holding
 * place). Use this to compute values for a pre-allocated
 * P&_PROFILE <gm_dest> based on a simple calculation involving
 * <gm>, <bg>, and <bg_tmp>.
 *
 */
static int
parameterize_gm (P7_BG *bg, P7_BG *bg_tmp, P7_PROFILE *gm_src, P7_PROFILE *gm_dest, const ESL_SQ *sq, int L, float *sc_tmp) {
  int     i, j;
  int     K   = gm_src->abc->K;
  int     Kp  = gm_src->abc->Kp;

  /* Fraction of new bg frequencies that comes from a prior determined by the sequence block.
   * This is 25% for long sequences, more for shorter sequences (e.g. 50% for sequences of length 50)
   */
  float   bg_smooth = 25.0 / (ESL_MIN(100,L));

  if (sq == NULL) return eslFAIL;


  /* Compute new background frequencies as a weighted avg of the
   * default bg and the observed frequencies dsq. Store in bg_tmp.
   */
  esl_vec_FSet (bg_tmp->f, gm_src->abc->K, 0);
  esl_sq_TallyCounts(sq, bg_tmp->f);
  esl_vec_FNorm(bg_tmp->f, gm_src->abc->K);

  esl_vec_FScale(bg_tmp->f, K, (1.0-bg_smooth));
  esl_vec_FAddScaled(bg_tmp->f, bg->f, bg_smooth, K);


  /* Compute new forward emission scores, based on default scores and
   * the new bg.
   *
   * The value at a given cell i,j in the rsc (fwd_score) matrix is
   *    sc[i][j]    = log(hmm->mat[i][j] / bg->f[j]);
   * We'd like to compute a new score, sc_tmp:
   *    sctmp[i][j] = log(hmm->mat[i][j] / bg_tmp->f[j]);
   *
   * Since log(hmm->mat[i][j] / bg->f[j])    == log(hmm->mat[i][j]) - log(bg->f[j])
   * and   log(hmm->mat[i][j] / bgtmp->f[j]) == log(hmm->mat[i][j]) - log(bgtmp->f[j])
   * we can set sctmp[i][j] = sc[i][j] + log(bg->f[j]) - log(bgtmp->f[j])
   *                        = sc[i][j] + log(bg->f[j]/bgtmp->f[j])
   */
  for (i = 1; i <=gm_src->M; i++) {

    for (j=0; j<K; j++) {
      if (gm_src->mm[0] != 0 && gm_src->mm[i] == 'm')
        sc_tmp[j] = 0;
      else
        sc_tmp[j] = gm_src->rsc[j][(i) * p7P_NR  + p7P_M] + log(bg->f[j]/bg_tmp->f[j]);
    }

    sc_tmp[K] = sc_tmp[Kp-2] = sc_tmp[Kp-1] = 0;
    esl_abc_FExpectScVec(gm_src->abc, sc_tmp, bg->f);

    for (j=0; j<Kp; j++)
      gm_dest->rsc[j][(i) * p7P_NR  + p7P_M] =  sc_tmp[j];

  }


  /* Set transitions */
  p7_profile_SetLength(gm_dest, L);

  return eslOK;
}


/* Function:  p7_pipeline_postViterbi_LongTarget()
 * Synopsis:  the part of the LongTarget P7 search Pipeline downstream
 *            of the Viterbi filter
 *
 * Purpose:   This is called by postMSV_LongTarget(), and runs the
 *            post-Viterbi part of H3's accelerated pipeline to
 *            compare profile <om> against sequence <sq>. If a
 *            significant hit is found, information about it is
 *            added to the <hitlist>.
 *            The pipeline accumulates beancounting information
 *            about how many comparisons (and residues) flow through
 *            the pipeline while it's active.
 *
 * Args:      pli             - the main pipeline object
 *            om              - optimized profile (query)
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin
 *            data         - for computing windows based on maximum prefix/suffix extensions
 *            seqidx          - the id # of the sequence from which the current window was extracted
 *            window_start    - the starting position of the extracted window (offset from the first
 *                              position of the block of a possibly longer sequence)
 *            window_len      - the length of the extracted window
 *            tmpseq          - a new or reused digital sequence object used for p7_alidisplay_Create() call
 *            subseq          - digital sequence of the extracted window
 *            seq_start       - first position of the sequence block passed in to the calling pipeline function
 *            seq_name        - name of the sequence the window comes from
 *            seq_source      - source of the sequence the window comes from
 *            seq_acc         - acc of the sequence the window comes from
 *            seq_desc        - desc of the sequence the window comes from
 *            nullsc          - score of the passed window vs the bg model
 *            usc             - msv score of the passed window
 *            complementarity - boolean; is the passed window sourced from a complementary sequence block
 *
 * Returns:   <eslOK> on success. If a significant hit is obtained,
 *            its information is added to the growing <hitlist>.
 *
 *            <eslERANGE> on numerical overflow errors in the
 *            optimized vector implementations; particularly in
 *            posterior decoding. I don't believe this is possible for
 *            multihit local models, but I'm set up to catch it
 *            anyway. We may emit a warning to the user, but cleanly
 *            skip the problematic sequence and continue.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 *
 */
static int
p7_pipeline_postViterbi_LongTarget(P7_PIPELINE *pli, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_TOPHITS *hitlist,
    const P7_SCOREDATA *data, int64_t seqidx, int window_start, int window_len, ESL_SQ *tmpseq, ESL_DSQ *subseq,
    int seq_start, char *seq_name, char *seq_source, char* seq_acc, char* seq_desc,
    int complementarity, int *overlap, P7_PIPELINE_TMP_OBJS *pli_tmp
)
{
  P7_DOMAIN        *dcl     = NULL;     /* array of domain data structures <0..ndom-1>, from the first trace */
  P7_DOMAIN        *dcl_tmp    = NULL;  /* After reparameterization, and re-tracing, store a working domain here, then either point the hit->dcl to it, or free it  */
  P7_HIT           *hit     = NULL;     /* ptr to the current hit output data      */
  float            fwdsc;   /* filter scores                           */
  float            nullsc;
  float            filtersc;           /* HMM null filter score                   */
  float            bias_filtersc;           /* HMM null filter score                   */
  float            seq_score;          /* the corrected per-seq bit score */
  double           P;               /* P-value of a hit */
  int              d, dd;
  int              status;
  int              nres;

//  int              Ld;               /* # of residues in envelopes */
  int              z,i;
  int              noverlaps;
  int              last_ibe;
  int              best_d;
  int              ii, zz;
  int              env_offset;


  int env_len, env_len2;
  int ali_len;
  float bitscore;
  float dom_score;
  double dom_lnP;
  ESL_DSQ *subseq_tmp;

  int max_env_extra = 20;
  int F3_L          = ESL_MIN( window_len,  pli->B3);

  p7_oprofile_ReconfigRestLength(om, window_len);
  p7_profile_SetLength(gm, window_len);
  p7_bg_SetLength(bg, window_len);

  p7_bg_NullOne  (bg, subseq, window_len, &nullsc);
  p7_bg_FilterScore(bg, subseq, window_len, &bias_filtersc);
  bias_filtersc -= nullsc;  //remove nullsc, so bias scaling can be done, then add it back on later

  /* ad hoc scaling down of bias, in case of long windows*/
  filtersc =  nullsc + (bias_filtersc * ( F3_L>window_len ? 1.0 : (float)F3_L/window_len) );


  /* Parse with Forward and obtain its real Forward score. */
  p7_checkptmx_Reuse (pli->cx);
  p7_ForwardFilter(subseq, window_len, om, pli->cx, &fwdsc);
  seq_score = (fwdsc - filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  if (P > pli->F3 ) return eslOK;

  pli->stats.pos_past_fwd += window_len - *overlap;

  *overlap = -1; // overload variable to tell calling function that this window passed fwd

  /* ok, it's for real; passes vectorized local scoring filters.
    * Finish with Backwards and Decoding, defining a sparse mask in <sm>.
    * In this case "domains" will end up being translated as independent "hits"
    */
   p7_BackwardFilter(subseq, window_len, om, pli->cx, pli->sm, /*p7_SPARSEMASK_THRESH_DEFAULT*/0.005);

   /* FIXME 3.1
    * hmmscan needs <gm>; read from the hmm file (.h3m) on disk
    * expect that this is slow and inefficient; come back and refactor later.
    * (in hmmpgmd, hfp is NULL; all our profiles are already in a cache)
    */
   if (pli->mode == p7_SCAN_MODELS && pli->hfp)
     workaround_get_profile(pli, om, bg, &gm);


   /* Now we can hand it over to sparse DP, with the full glocal/local model */
   p7_trace_Reuse   (pli->tr);
   p7_SparseViterbi (subseq, window_len, gm, pli->sm,  pli->sxx, pli->tr, &fwdsc); // don't care what the function returns, so used fwdsc just to make it happy
   p7_SparseForward (subseq, window_len, gm, pli->sm,  pli->sxf,          &fwdsc);
   p7_SparseBackward(subseq, window_len, gm, pli->sm,  pli->sxb,          NULL);
   p7_SparseDecoding(subseq, window_len, gm, pli->sxf, pli->sxb, pli->sxd);
   p7_sparsemx_TracePostprobs(pli->sxd, pli->tr); /* annotate the trace. */
   p7_trace_Index(pli->tr);                       /* index domains in the trace */
   p7_sparsemx_Reuse(pli->sxx);                   /* reuse the Viterbi matrix for mass tracing, coming up next */


   ESL_ALLOC(dcl, sizeof(P7_DOMAIN) * pli->tr->ndom);
   ESL_REALLOC(pli->wrk,  sizeof(float) * (gm->M+1));
   ESL_REALLOC(pli->n2sc, sizeof(float) * (window_len+1));
   esl_vec_FSet(pli->n2sc, window_len+1, 0.0f);

   noverlaps  = 0;
   last_ibe   = -1;
   best_d     = 0;


  /* For each domain range identified above,
   *  - determine an envelope,
   *  - reparameterize the model according to the (regularized)
   *    nucleotide frequencies in the envelope (this is done in
   *    a pre-allocated om and gm),
   *  - copy the appropriate part of the sparse matrix into a
   *    new smaller pre-allocated mask,
   *  - re-run the sparse trace and envelope pipeline
   *  - compute new score
   *  - stick the resulting domains onto the hit list (not
   *    all will be above threshold)
   *
   * The passing "domains" are stored as independent hits so
   * the remainder of the typical-case hit-merging process can
   * remain mostly intact.
   *
   * Note: we expect only a single reparameterized hit to come out
   * of a single envelope. But there may be cases of overlapping
   * envelopes from above stage, where the reparameterized
   * alignments overlap (one base is found in more than one hit).
   * That's not kosher. This can be avoided by merging overlapping
   * envelopes into a single longer envelope, then passing this
   * through the reparameterized pipeline. This has the potential
   * to result in multiple hits, so we must loop over all hits
   * in a sub-loop.
   *
   * Some of them may not pass eventual E-value thresholds. In
   * protein context, these would be reported as supplementary
   * data (domains contributing to a full-sequence score), but
   * in nhmmer context, they'll just get thrown away later, so
   * drop them now, if possible.
   *
   * So first ... merge overlapping envelopes :
   */
   dd=-1;
   for (d = 0; d < pli->tr->ndom; d++)
   {
     /* Determine envelope coords by mass trace. (with the default parameters) */
     p7_SparseMasstrace(subseq, window_len, gm, pli->sxf, pli->sxb, pli->tr, pli->tr->anch[d], p7_MASSTRACE_THRESH_DEFAULT, pli->sxx, pli->mt,
                        &(dcl[d].iae), &(dcl[d].ibe), &(dcl[d].kae), &(dcl[d].kbe));
     p7_sparsemx_Reuse (pli->sxx);
     p7_masstrace_Reuse(pli->mt);

     //if the envelope is very long, trim it
     dcl[d].iae = ESL_MAX(dcl[d].iae,  pli->tr->sqfrom[d] - max_env_extra );
     dcl[d].ibe = ESL_MIN(dcl[d].ibe,  pli->tr->sqto[d] + max_env_extra );


     if (d>0 && dcl[d].iae <= dcl[dd].ibe) {
       dcl[dd].ibe = ESL_MAX(dcl[d].ibe, dcl[dd].ibe);
     } else {
       dd++;
       dcl[dd].iae = dcl[d].iae;
       dcl[dd].ibe = dcl[d].ibe;
     }
   }


   pli->tr->ndom = dd+1;


   for (d = 0; d < pli->tr->ndom; d++)
   {

     /* I have an envelope based on the default null model; now update that
      * null model based on the contents of the current envelope (effectively
      * a partial null3 correction), then re-run everything required to get the
      * corrected score/trace/alignment for the domain
      */
     subseq_tmp = subseq + dcl[d].iae - 1;
     env_len    = dcl[d].ibe-dcl[d].iae+1;
     env_offset = dcl[d].iae;

     /* set up seq object required for parameterize_gm(). We only need
      * to have these three fields filled in for the parameterize_gm() call.
      * Will fill in the rest if necessary
      */
     tmpseq->n = env_len;
     tmpseq->seq = NULL;
     tmpseq->dsq = subseq_tmp; //using this instead of copying, I need to remember to set ->dsq to NULL before destroying tmpseq

     parameterize_gm (bg, pli_tmp->bg, gm, pli_tmp->gm, tmpseq, env_len, pli_tmp->scores);
     //pli_tmp->gm = p7_profile_Clone(gm); // this keeps the default bg-based scoring

     /* copy over the part of the sparse mask related to the current envelope*/
     p7_sparsemask_Reinit(pli_tmp->sm, gm->M, env_len);
     for (i=env_len; i >= 1; i--)
     {
       ii = i + env_offset - 1;
       if (pli->sm->n[ii] > 0) {
         if (  (status = p7_sparsemask_StartRow(pli_tmp->sm, i))                != eslOK) return status;
         for (z = pli->sm->n[ii]-1; z >= 0; z--) {
           zz = pli->sm->k[ii][z]-1;
           if ((status = p7_sparsemask_Add(pli_tmp->sm, zz%pli->sm->Q, zz/pli->sm->Q))   != eslOK) return status;
         }
         if (  (status = p7_sparsemask_FinishRow(pli_tmp->sm))                      != eslOK) return status;
       }
     }
     p7_sparsemask_Finish(pli_tmp->sm);

     p7_SparseViterbi (subseq_tmp, env_len, pli_tmp->gm, pli_tmp->sm,  pli_tmp->sxx, pli_tmp->trc, &fwdsc); // don't care what the function returns, so used fwdsc just to make it happy
     p7_SparseForward (subseq_tmp, env_len, pli_tmp->gm, pli_tmp->sm,  pli_tmp->sxf,               &fwdsc);
     p7_SparseBackward(subseq_tmp, env_len, pli_tmp->gm, pli_tmp->sm,  pli_tmp->sxb,               NULL);
     p7_SparseDecoding(subseq_tmp, env_len, pli_tmp->gm, pli_tmp->sxf, pli_tmp->sxb, pli_tmp->sxd);
     p7_sparsemx_TracePostprobs(pli_tmp->sxd, pli_tmp->trc);
     p7_trace_Index(pli_tmp->trc);
     p7_sparsemx_Reuse(pli_tmp->sxx);

     if (pli_tmp->trc->ndom > 0)
       p7_oprofile_Convert(pli_tmp->gm, pli_tmp->om);

     for (dd = 0; dd < pli_tmp->trc->ndom; dd++) {
       dcl_tmp = p7_domain_Create(1);

       p7_SparseMasstrace(subseq_tmp, env_len, pli_tmp->gm, pli_tmp->sxf, pli_tmp->sxb, pli_tmp->trc, pli_tmp->trc->anch[dd], p7_MASSTRACE_THRESH_DEFAULT, pli_tmp->sxx, pli_tmp->mt,
                             &(dcl_tmp->iae), &(dcl_tmp->ibe), &(dcl_tmp->kae), &(dcl_tmp->kbe));
       p7_sparsemx_Reuse (pli_tmp->sxx);
       p7_masstrace_Reuse(pli_tmp->mt);


       /* Transfer alignment coords from trace. [Do we need to do this?] */
       dcl_tmp->ia = pli_tmp->trc->sqfrom[dd] + env_offset - 1;
       dcl_tmp->ib = pli_tmp->trc->sqto[dd]   + env_offset - 1;
       dcl_tmp->ka = pli_tmp->trc->hmmfrom[dd];
       dcl_tmp->kb = pli_tmp->trc->hmmto[dd];
       p7_SparseEnvscore(subseq_tmp, env_len, pli_tmp->gm, dcl_tmp->iae, dcl_tmp->ibe, dcl_tmp->kae, dcl_tmp->kbe, pli_tmp->sm, pli_tmp->sxx, &(dcl_tmp->envsc));
       p7_sparsemx_Reuse(pli_tmp->sxx);

       /* this is an assessment of the null correction */
       p7_SparseEnvscore(subseq_tmp, env_len,          gm, dcl_tmp->iae, dcl_tmp->ibe, dcl_tmp->kae, dcl_tmp->kbe, pli_tmp->sm, pli_tmp->sxx, &(dcl_tmp->domcorrection));

       dcl_tmp->dombias = 0.0;
       p7_sparsemx_Reuse(pli_tmp->sxx);


       dcl_tmp->iae +=  env_offset - 1;
       dcl_tmp->ibe +=  env_offset - 1;

       /* Keep track of overlaps */
       if (dcl_tmp->iae <= last_ibe) noverlaps++;
       last_ibe   = dcl_tmp->ibe;


       /* alignment "accuracy score": expected number of correctly aligned positions */
       dcl_tmp->oasc = 0.;
       for (z = pli_tmp->trc->tfrom[dd]; z <= pli_tmp->trc->tto[dd]; z++)
         if (pli_tmp->trc->i[z]) dcl_tmp->oasc += pli_tmp->trc->pp[z];

       /* We're initiazing a P7_DOMAIN structure in dcl[d] by hand, without a Create().
        * We're responsible for initiazing all elements of this structure.
        */
       dcl_tmp->is_reported = FALSE; /* will get set later by p7_tophits_Threshold() */
       dcl_tmp->is_included = FALSE; /* ditto */



      /* note: the bitscore of a hit depends on the window_len of the
       * current window. Here, the score is modified (reduced) by treating
       * all passing windows as though they came from windows of length
       * om->max_length. For details, see
       * ~wheelert/archive/notebook/2012/0130_bits_v_evalues/00NOTES (Feb 1)
       */
       //adjust the score of a hit to account for the full length model - the characters outside the envelope but in the window
       env_len2 = dcl_tmp->ibe - dcl_tmp->iae + 1;
       ali_len  = dcl_tmp->ib - dcl_tmp->ia + 1;
       bitscore = dcl_tmp->envsc ;
       //For these modifications, see notes, ~/notebook/2010/0716_hmmer_score_v_eval_bug/, end of Thu Jul 22 13:36:49 EDT 2010
       bitscore -= 2 * log(2. / (env_len+2))          +   (env_len2-ali_len)            * log((float)env_len / (env_len+2));
       bitscore += 2 * log(2. / (om->max_length+2)) ;
       //the ESL_MAX test handles the extremely rare case that the env_len is actually larger than om->max_length
       bitscore +=  (ESL_MAX(om->max_length, env_len) - ali_len) * log((float)om->max_length / (float) (om->max_length+2));

      /*compute scores used to decide if we should keep this "domain" as a hit.
       */
       dom_score  = (bitscore - nullsc)  / eslCONST_LOG2;
       dom_lnP   = esl_exp_logsurv(dom_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

       //note: this test is conservative:
       // - if in scan mode, assume query sequence is at most the length of a single instance,
       // - otherwise, use just the nres from the current pipeline, while the final filters
       //   will add up nres from all threads, over all windows, which will increase stringency
       if (pli->mode == p7_SCAN_MODELS)  nres = ESL_MIN(env_len2,om->max_length);
       else                              nres = pli->stats.nres;

       if (!( p7_pipeline_TargetReportable(pli, dom_score, dom_lnP + log((float)nres / om->max_length) ) ) ){
         p7_domain_Destroy(dcl_tmp,1);
       } else {

          /*now that almost everything has been filtered away, set up seq object required
           * for the p7_alidisplay_Create() call, if it wasn't already done for an earlier domain
           * on the same sequence */
          if (tmpseq->name[0] !=  '\0') {
            if ((status = esl_sq_SetName     (tmpseq, seq_name))   != eslOK) goto ERROR;
            if ((status = esl_sq_SetSource   (tmpseq, seq_source)) != eslOK) goto ERROR;
            if ((status = esl_sq_SetAccession(tmpseq, seq_acc))    != eslOK) goto ERROR;
            if ((status = esl_sq_SetDesc     (tmpseq, seq_desc))   != eslOK) goto ERROR;
          }
          /* Viterbi alignment of the domain */
          dcl_tmp->ad = p7_alidisplay_Create(pli_tmp->trc, dd, pli_tmp->om, tmpseq); // it's ok to use the om instead of the gm here; it doesn't depend on updated fwd scores


          pli_tmp->trc->sqfrom[dd] += env_offset - 1;
          pli_tmp->trc->sqto[dd]   += env_offset - 1;
          dcl_tmp->ad->sqfrom      += env_offset - 1;
          dcl_tmp->ad->sqto        += env_offset - 1;

          //I think this is broken: the "scores" array should be with the reparameterized values
          p7_pipeline_computeAliScores(dcl_tmp, subseq, data, om->abc->Kp); // using subseq, because dcl[d] contains offsets into that sequence

          p7_tophits_CreateNextHit(hitlist, &hit);
          hit->dcl = dcl_tmp;

          hit->ndom        = 1;
          hit->best_domain = 0;

          hit->window_length = om->max_length;
          hit->seqidx = seqidx;
          hit->subseq_start = seq_start;

          if (complementarity == p7_NOCOMPLEMENT) {
            hit->dcl[0].iae += window_start - 1; // represents the real position within the sequence handed to the pipeline
            hit->dcl[0].ibe += window_start - 1;
            hit->dcl[0].ia += window_start - 1;
            hit->dcl[0].ib += window_start - 1;
            hit->dcl[0].ad->sqfrom += window_start - 1;
            hit->dcl[0].ad->sqto += window_start - 1;
          } else {

            hit->dcl[0].iae = window_start + window_len - 1 - hit->dcl[0].iae + 1; // represents the real position within the sequence handed to the pipeline
            hit->dcl[0].ibe = window_start + window_len - 1 - hit->dcl[0].ibe + 1;
            hit->dcl[0].ia = window_start + window_len - 1 - hit->dcl[0].ia + 1;
            hit->dcl[0].ib = window_start + window_len - 1 - hit->dcl[0].ib + 1;
            hit->dcl[0].ad->sqfrom = window_start + window_len - 1 - hit->dcl[0].ad->sqfrom + 1;
            hit->dcl[0].ad->sqto   = window_start + window_len - 1 - hit->dcl[0].ad->sqto + 1;
          }


          hit->pre_score = bitscore  / eslCONST_LOG2;
          hit->pre_lnP   = esl_exp_logsurv (hit->pre_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);


          //hit->dcl[0].dombias  = dcl_tmp->dombias;
          hit->sum_score  = hit->score  = hit->dcl[0].bitscore = dom_score;
          hit->sum_lnP    = hit->lnP    = hit->dcl[0].lnP  = dom_lnP;
          hit->sortkey       = pli->inc_by_E ? -dom_lnP : dom_score; /* per-seq output sorts on bit score if inclusion is by score  */


          if (pli->mode == p7_SEARCH_SEQS)
          {
            if (                       (status  = esl_strdup(seq_name, -1, &(hit->name)))  != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
            if (seq_acc[0]  != '\0' && (status  = esl_strdup(seq_acc,  -1, &(hit->acc)))   != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
            if (seq_desc[0] != '\0' && (status  = esl_strdup(seq_desc, -1, &(hit->desc)))  != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
          } else {
            if ((status  = esl_strdup(om->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure");
            if ((status  = esl_strdup(om->acc,  -1, &(hit->acc)))   != eslOK) esl_fatal("allocation failure");
            if ((status  = esl_strdup(om->desc, -1, &(hit->desc)))  != eslOK) esl_fatal("allocation failure");
          }


          /* If using model-specific thresholds, filter now.  See notes in front
           * of the analogous piece of code in p7_Pipeline() for further explanation
           * of timing.
           */
          hit->flags       = 0;
          if (pli->use_bit_cutoffs)
          {
            if (p7_pipeline_TargetReportable(pli, hit->score, hit->lnP))
            {
              hit->flags |= p7_IS_REPORTED;
              if (p7_pipeline_TargetIncludable(pli, hit->score, hit->lnP))
                hit->flags |= p7_IS_INCLUDED;
            }

            if (p7_pipeline_DomainReportable(pli, hit->dcl[0].bitscore, hit->dcl[0].lnP))
            {
              hit->dcl[0].is_reported = TRUE;
              if (p7_pipeline_DomainIncludable(pli, hit->dcl[0].bitscore, hit->dcl[0].lnP))
                hit->dcl[0].is_included = TRUE;
            }

          }


       }
     }
     tmpseq->dsq    = NULL;
     esl_sq_Reuse(tmpseq);

     p7_trace_Reuse   (pli_tmp->trc);
  }


  if (dcl != NULL) free (dcl);
  return eslOK;

ERROR:
  if (dcl != NULL) free (dcl);

  ESL_EXCEPTION(eslEMEM, "Error in LongTarget pipeline\n");
}


/* Function:  p7_pipeline_postMSV_LongTarget()
 * Synopsis:  the part of the LongTarget P7 search Pipeline downstream
 *            of the MSV filter
 *
 * Purpose:   This is called by either the standard (SIMD-MSV) long-target
 *            pipeline (p7_Pipeline_LongTarget) or the FM-index long-target
 *            pipeline (p7_Pipeline_FM), and runs the post-MSV part of H3's
 *            accelerated pipeline to compare profile <om> against
 *            sequence <sq>. If a significant hit is found,
 *            information about it is added to the <hitlist>.
 *            The pipeline accumulates beancounting information
 *            about how many comparisons (and residues) flow through
 *            the pipeline while it's active.
 *
 * Args:      pli             - the main pipeline object
 *            gm              - profile (query)
 *            om              - optimized profile (query)
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin
 *            data         - for computing windows based on maximum prefix/suffix extensions
 *            seqidx          - the id # of the sequence from which the current window was extracted
 *            window_start    - the starting position of the extracted window (offset from the first
 *                              position of the block of a possibly longer sequence)
 *            window_len      - the length of the extracted window
 *            tmpseq          - a new or reused digital sequence object used for "domain" definition
 *            subseq          - digital sequence of the extracted window
 *            seq_start       - first position of the sequence block passed in to the calling pipeline function
 *            seq_name        - name of the sequence the window comes from
 *            seq_source      - source of the sequence the window comes from
 *            seq_acc         - acc of the sequence the window comes from
 *            seq_desc        - desc of the sequence the window comes from
 *            nullsc          - score of the passed window vs the bg model
 *            usc             - msv score of the passed window
 *            complementarity - boolean; is the passed window sourced from a complementary sequence block
 *
 * Returns:   <eslOK> on success. If a significant hit is obtained,
 *            its information is added to the growing <hitlist>.
 *
 *            <eslERANGE> on numerical overflow errors in the
 *            optimized vector implementations; particularly in
 *            posterior decoding. I don't believe this is possible for
 *            multihit local models, but I'm set up to catch it
 *            anyway. We may emit a warning to the user, but cleanly
 *            skip the problematic sequence and continue.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 */
static int
p7_pipeline_postMSV_LongTarget(P7_PIPELINE *pli, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_TOPHITS *hitlist, const P7_SCOREDATA *data,
    int64_t seqidx, int window_start, int window_len, ESL_SQ *tmpseq,
    ESL_DSQ *subseq, int seq_start, char *seq_name, char *seq_source, char* seq_acc, char* seq_desc,
    float nullsc, float usc, int complementarity, P7_HMM_WINDOWLIST *vit_windowlist,
    P7_PIPELINE_TMP_OBJS *pli_tmp
)
{
  float            filtersc;           /* HMM null filter score                   */
  float            bias_filtersc;      /* HMM null filter score                   */
  float            seq_score;          /* the corrected per-seq bit score */
  double           P;               /* P-value of a hit */
  int i;
  int overlap;
  int new_n;
  int new_len;

  int   loc_window_len;  //used to re-parameterize to shorter target windows

  int max_window_len      = 80000;
  int smaller_window_len  = 40000;

  int F1_L = ESL_MIN( window_len,  pli->B1);
  int F2_L = ESL_MIN( window_len,  pli->B2);

  //initial bias filter, based on the input window_len
  if (pli->do_biasfilter) {
      p7_bg_SetLength(bg, window_len);
      p7_bg_FilterScore(bg, subseq, window_len, &bias_filtersc);
      bias_filtersc -= nullsc; // doing this because I'll be modifying the bias part of filtersc based on length, then adding nullsc back in.
      filtersc =  nullsc + (bias_filtersc * (float)(( F1_L>window_len ? 1.0 : (float)F1_L/window_len)));
      seq_score = (usc - filtersc) / eslCONST_LOG2;
      P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      if (P > pli->F1) return eslOK;
  } else {
    bias_filtersc = 0; // mullsc will be added in later
  }
  pli->stats.pos_past_bias += window_len;

  //establish a possibly shorter target window parameterization
  loc_window_len = ESL_MIN(window_len,om->max_length);

  //compute the new nullsc based on possibly shorter window
  p7_bg_SetLength(bg, loc_window_len);
  p7_bg_NullOne  (bg, subseq, loc_window_len, &nullsc);

  // bias_filtersc has already been reduced by nullsc based on window_len
  // We compute a --B2-scaled bias, then tack on the nullsc based on the new,
  // possibly shorter length model
  filtersc =  nullsc + (bias_filtersc * ( F2_L>window_len ? 1.0 : (float)F2_L/window_len) );


  //Then configure the model length based on the possibly shorter window length
  p7_oprofile_ReconfigRestLength(om, loc_window_len);

  /* Second level viterbi filter*/
  //use window_len instead of loc_window_len, because length parameterization is done, just need to loop over subseq
  p7_ViterbiFilter_longtarget(subseq, window_len, om, pli->fx, filtersc, pli->F2, vit_windowlist);


  p7_pipeline_ExtendAndMergeWindows (om, data, vit_windowlist, window_len, 0.5);


  // if a window is still too long (>80Kb), need to split it up to
  // ensure numeric stability in Fwd.
  for (i=0; i<vit_windowlist->count; i++) {

      if (vit_windowlist->windows[i].length > max_window_len) {

         //modify the current window to restrict length to 40K, then add
         //new windows with max length 40K, and MAXL overlap w/ preceding window
         new_n   = vit_windowlist->windows[i].n ;
         new_len = vit_windowlist->windows[i].length ;
         vit_windowlist->windows[i].length = smaller_window_len;

         do {
           new_n   +=  (smaller_window_len - om->max_length);
           new_len -=  (smaller_window_len - om->max_length);
           p7_hmmwindow_new(vit_windowlist, 0, new_n, 0, 0, ESL_MIN(smaller_window_len,new_len), 0.0, p7_NOCOMPLEMENT );
         } while (new_len > smaller_window_len);
      }
  }

  overlap = 0;
  for (i=0; i<vit_windowlist->count; i++) {
    pli->stats.pos_past_vit += vit_windowlist->windows[i].length;
    //remove overlap with preceding window
    if (i>0)
      pli->stats.pos_past_vit -= ESL_MAX(0,  vit_windowlist->windows[i-1].n + vit_windowlist->windows[i-1].length - vit_windowlist->windows[i].n );

    p7_pipeline_postViterbi_LongTarget(pli, gm, om, bg, hitlist, data, seqidx,
        window_start+vit_windowlist->windows[i].n-1, vit_windowlist->windows[i].length, tmpseq,
        subseq + vit_windowlist->windows[i].n - 1,
        seq_start, seq_name, seq_source, seq_acc, seq_desc, complementarity, &overlap,
        pli_tmp
    );

    if (overlap == -1 && i<vit_windowlist->count-1) {
      overlap = ESL_MAX(0,  vit_windowlist->windows[i].n + vit_windowlist->windows[i].length - vit_windowlist->windows[i+1].n );
    } else {
      //that window didn't pass Fwd
      overlap = 0;
    }

    pli->tr->ndom = 0;

  }


  return eslOK;

}

/* Function:  p7_Pipeline_LongTarget()
 * Synopsis:  HMMER3's accelerated seq/profile comparison pipeline, modified to use the scanning MSV/SSV filters.
 *
 * Purpose:   Run H3's accelerated pipeline to compare profile <om>
 *            against sequence <sq>. If a significant hit is found,
 *            information about it is added to the <hitlist>.
 *            This is a variant of p7_Pipeline that runs the
 *            versions of the MSV/SSV filters that scan a long
 *            sequence and find high-scoring regions (windows), then pass 
 *            those to the remainder of the pipeline. The pipeline
 *            accumulates beancounting information about how many comparisons
 *            flow through the pipeline while it's active.
 *
 * Returns:   <eslOK> on success. If a significant hit is obtained,
 *            its information is added to the growing <hitlist>.
 *
 *            <eslEINVAL> if (in a scan pipeline) we're supposed to
 *            set GA/TC/NC bit score thresholds but the model doesn't
 *            have any.
 *
 *            <eslERANGE> on numerical overflow errors in the
 *            optimized vector implementations; particularly in
 *            posterior decoding. I don't believe this is possible for
 *            multihit local models, but I'm set up to catch it
 *            anyway. We may emit a warning to the user, but cleanly
 *            skip the problematic sequence and continue.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      J4/25.
 */
int
p7_Pipeline_LongTarget(P7_PIPELINE *pli, P7_PROFILE *gm, P7_OPROFILE *om, P7_SCOREDATA *data, P7_BG *bg, const ESL_SQ *sq, P7_TOPHITS *hitlist, int64_t seqidx)
{
  int              i;
  int              status;
  float            nullsc;   /* null model score                        */
  float            usc;      /* msv score  */
  float            P;
  ESL_DSQ          *subseq;
  ESL_SQ           *tmpseq   = NULL;
  float            bias_filtersc;

  P7_HMM_WINDOWLIST msv_windowlist;
  P7_HMM_WINDOWLIST vit_windowlist;

  P7_PIPELINE_TMP_OBJS *pli_tmp = NULL;


  if (sq->n == 0) return eslOK;    /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */

  ESL_ALLOC(pli_tmp, sizeof(P7_PIPELINE_TMP_OBJS));
  pli_tmp->bg = p7_bg_Clone(bg);
  pli_tmp->gm = p7_profile_Clone(gm);
  pli_tmp->om = p7_oprofile_Create(gm->M, gm->abc);
  pli_tmp->sm = p7_sparsemask_Create(gm->M, 100);
  ESL_ALLOC(pli_tmp->scores, sizeof(float) * om->abc->Kp);
  if ( (pli_tmp->trc = p7_trace_CreateWithPP())            == NULL) goto ERROR;
  if ( (pli_tmp->sxf = p7_sparsemx_Create (pli_tmp->sm))   == NULL) goto ERROR;
  if ( (pli_tmp->sxb = p7_sparsemx_Create (pli_tmp->sm))   == NULL) goto ERROR;
  if ( (pli_tmp->sxd = p7_sparsemx_Create (pli_tmp->sm))   == NULL) goto ERROR;
  if ( (pli_tmp->sxx = p7_sparsemx_Create (pli_tmp->sm))   == NULL) goto ERROR;
  if ( (pli_tmp->mt  = p7_masstrace_Create (gm->M, 100))   == NULL) goto ERROR;

  msv_windowlist.windows = NULL;
  vit_windowlist.windows = NULL;
  p7_hmmwindow_init(&msv_windowlist);


  /* Set false target length. This is a conservative estimate of the length of window that'll
   * soon be passed on to later phases of the pipeline;  used to recover some bits of the score
   * that we would miss if we left length parameters set to the full target length */
  p7_oprofile_ReconfigMSVLength(om, om->max_length);

  /* First level filter: the SSV filter, with <om>.
   * This variant of SSV will scan a long sequence and find
   * short high-scoring regions.
   */
  p7_MSVFilter_longtarget(sq->dsq, sq->n, om, pli->fx, data, bg, pli->F1, &msv_windowlist);

  /* convert hits to windows, possibly filtering based on composition bias,
   * definitely merging neighboring windows, and
   * TODO: splitting overly-large windows
   */
  if ( msv_windowlist.count > 0 ) {

    /* In scan mode, if it passes the MSV filter, read the rest of the profile */
    if (pli->mode == p7_SCAN_MODELS)
    {
      if (om->base_w == 0 &&  om->scale_w == 0) { // we haven't already read this hmm (if we're on the second strand, we would've)
        if (pli->hfp) p7_oprofile_ReadRest(pli->hfp, om);
        if ((status = p7_pipeline_NewModelThresholds(pli, om)) != eslOK) goto ERROR;
      }
    }

    if (data->prefix_lengths == NULL) { //otherwise, already filled in
      p7_hmm_ScoreDataComputeRest(om, data);
    }

    p7_pipeline_ExtendAndMergeWindows (om, data, &msv_windowlist, sq->n, 0);


  /*
   * pass each remaining window on to the remaining pipeline
   */
    p7_hmmwindow_init(&vit_windowlist);
    tmpseq = esl_sq_CreateDigital(sq->abc);
    free (tmpseq->dsq);  //this ESL_SQ object is just a container that'll point to a series of other DSQs, so free the one we just created inside the larger SQ object

    for (i=0; i<msv_windowlist.count; i++){
      int window_len = msv_windowlist.windows[i].length;

      subseq = sq->dsq + msv_windowlist.windows[i].n - 1;
      p7_bg_SetLength(bg, window_len);
      p7_bg_NullOne  (bg, subseq, window_len, &nullsc);

      p7_bg_FilterScore(bg, subseq, window_len, &bias_filtersc);

      // Compute standard MSV to ensure that bias doesn't overcome SSV score when MSV
      // would have survived it
      p7_oprofile_ReconfigMSVLength(om, window_len);
      p7_MSVFilter(subseq, window_len, om, pli->fx, &usc);

      P = esl_gumbel_surv( (usc-nullsc)/eslCONST_LOG2,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      if (P > pli->F1 ) continue;

      pli->stats.pos_past_msv += window_len;

      status = p7_pipeline_postMSV_LongTarget(pli, gm, om, bg, hitlist, data, seqidx, msv_windowlist.windows[i].n, window_len, tmpseq,
                        subseq, sq->start, sq->name, sq->source, sq->acc, sq->desc, nullsc, usc, p7_NOCOMPLEMENT, &vit_windowlist,
                        pli_tmp
      );

      if (status != eslOK) goto ERROR;

    }

    tmpseq->dsq = NULL;  //It's just a pointer to another sequence's dsq entry, which has already been freed.
    esl_sq_Destroy(tmpseq);
    free (vit_windowlist.windows);
  }

  if (msv_windowlist.windows != NULL) free (msv_windowlist.windows);

  if (pli_tmp != NULL) {
    if (pli_tmp->bg != NULL)     p7_bg_Destroy(pli_tmp->bg);
    if (pli_tmp->gm != NULL)     p7_profile_Destroy(pli_tmp->gm);
    if (pli_tmp->om != NULL)     p7_oprofile_Destroy(pli_tmp->om);
    if (pli_tmp->sxf != NULL)    p7_sparsemx_Destroy(pli_tmp->sxf);
    if (pli_tmp->sxb != NULL)    p7_sparsemx_Destroy(pli_tmp->sxb);
    if (pli_tmp->sxd != NULL)    p7_sparsemx_Destroy(pli_tmp->sxd);
    if (pli_tmp->sxx != NULL)    p7_sparsemx_Destroy(pli_tmp->sxx);
    if (pli_tmp->sm != NULL)     p7_sparsemask_Destroy(pli_tmp->sm);
    if (pli_tmp->trc != NULL)    p7_trace_Destroy(pli_tmp->trc);
    if (pli_tmp->mt != NULL)     p7_masstrace_Destroy(pli_tmp->mt);
    if (pli_tmp->scores != NULL) free(pli_tmp->scores);
    free(pli_tmp);
  }

  return eslOK;





ERROR:
  if (tmpseq != NULL) esl_sq_Destroy(tmpseq);
  if (msv_windowlist.windows != NULL) free (msv_windowlist.windows);
  if (vit_windowlist.windows != NULL) free (vit_windowlist.windows);
  if (pli_tmp != NULL) {
    if (pli_tmp->bg != NULL)     p7_bg_Destroy(pli_tmp->bg);
    if (pli_tmp->gm != NULL)     p7_profile_Destroy(pli_tmp->gm);
    if (pli_tmp->om != NULL)     p7_oprofile_Destroy(pli_tmp->om);
    if (pli_tmp->sxf != NULL)    p7_sparsemx_Destroy(pli_tmp->sxf);
    if (pli_tmp->sxb != NULL)    p7_sparsemx_Destroy(pli_tmp->sxb);
    if (pli_tmp->sxd != NULL)    p7_sparsemx_Destroy(pli_tmp->sxd);
    if (pli_tmp->sxx != NULL)    p7_sparsemx_Destroy(pli_tmp->sxx);
    if (pli_tmp->sm != NULL)     p7_sparsemask_Destroy(pli_tmp->sm);
    if (pli_tmp->trc != NULL)    p7_trace_Destroy(pli_tmp->trc);
    if (pli_tmp->mt != NULL)     p7_masstrace_Destroy(pli_tmp->mt);
    if (pli_tmp->scores != NULL) free(pli_tmp->scores);
    free(pli_tmp);
  }

  return status;

}


/*------------- end, nhmmer/longtarget pipelines ----------------*/

/*****************************************************************
 * 7. Example 1: "search mode" in a sequence db
 *****************************************************************/

#ifdef p7PIPELINE_EXAMPLE
/* ./pipeline_example <hmmfile> <sqfile>
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles   reqs   incomp                             help                                                  docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,                          "show brief help on version and usage",                         0 },
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "E-value cutoff for reporting significant sequence hits",       0 },
  { "-T",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "bit score cutoff for reporting significant sequence hits",     0 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                          "set # of comparisons done, for E-value calculation",           0 },
  { "--domE",       eslARG_REAL,"1000.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "E-value cutoff for reporting individual domains",              0 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "bit score cutoff for reporting individual domains",            0 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                          "set # of significant seqs, for domain E-value calculation",    0 },
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--seqE,--seqT,--domE,--domT", "use GA gathering threshold bit score cutoffs in <hmmfile>",    0 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--seqE,--seqT,--domE,--domT", "use NC noise threshold bit score cutoffs in <hmmfile>",        0 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--seqE,--seqT,--domE,--domT", "use TC trusted threshold bit score cutoffs in <hmmfile>",      0 },
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--F1,--F2,--F3",               "Turn all heuristic filters off (less speed, more power)",      0 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,      NULL,  NULL, "--max",                        "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             0 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,      NULL,  NULL, "--max",                        "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             0 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,      NULL,  NULL, "--max",                        "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             0 },
  { "--nobias",     eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL, "--max",                        "turn off composition bias filter",                             0 },
  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL,  NULL,                          "turn off biased composition score corrections",                0 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",    NULL,  NULL,  NULL,                          "set RNG seed to <n> (if 0: one-time arbitrary seed)",          0 },
  { "--acc",        eslARG_NONE,  FALSE,  NULL, NULL,      NULL,  NULL,  NULL,                          "output target accessions instead of names if possible",        0 },
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqdb>";
static char banner[] = "example of using acceleration pipeline in search mode (seq targets)";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char         *hmmfile = esl_opt_GetArg(go, 1);
  char         *seqfile = esl_opt_GetArg(go, 2);
  int           format  = eslSQFILE_FASTA;
  P7_HMMFILE   *hfp     = NULL;
  ESL_ALPHABET *abc     = NULL;
  P7_BG        *bg      = NULL;
  P7_HMM       *hmm     = NULL;
  P7_PROFILE   *gm      = NULL;
  P7_OPROFILE  *om      = NULL;
  ESL_SQFILE   *sqfp    = NULL;
  ESL_SQ       *sq      = NULL;
  P7_PIPELINE  *pli     = NULL;
  P7_TOPHITS   *hitlist = NULL;
  int           h,d,namew;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Open a sequence file */
  if (esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp) != eslOK) p7_Fail("Failed to open sequence file %s\n", seqfile);
  sq = esl_sq_CreateDigital(abc);

  /* Create a pipeline and a top hits list */
  pli     = p7_pipeline_Create(go, hmm->M, 400, FALSE, p7_SEARCH_SEQS);
  hitlist = p7_tophits_Create(p7_TOPHITS_DEFAULT_INIT_ALLOC);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_oprofile_Convert(gm, om);   
  p7_pipeline_NewModel(pli, om, bg);

  /* Run each target sequence through the pipeline */
  while (esl_sqio_Read(sqfp, sq) == eslOK)
    { 
      p7_pipeline_NewSeq(pli, sq);
      p7_bg_SetLength(bg, sq->n);
      p7_profile_SetLength(gm, sq->n);
      p7_oprofile_ReconfigLength(om, sq->n);
  
      p7_Pipeline(pli, gm, om, bg, sq, hitlist);

      esl_sq_Reuse(sq);
      p7_pipeline_Reuse(pli);
    }

  /* Print the results. 
   * This example is a stripped version of hmmsearch's tabular output.
   */
  p7_tophits_SortBySortkey(hitlist);
  namew = ESL_MAX(8, p7_tophits_GetMaxNameLength(hitlist));
  for (h = 0; h < hitlist->N; h++)
    {
      d    = hitlist->hit[h]->best_domain;

      printf("%10.2g %7.1f %6.1f  %7.1f %6.1f %10.2g  %6.1f %5d  %-*s %s\n",
             exp(hitlist->hit[h]->lnP) * (double) pli->Z,
             hitlist->hit[h]->score,
             hitlist->hit[h]->pre_score - hitlist->hit[h]->score, /* bias correction */
             hitlist->hit[h]->dcl[d].bitscore,
             eslCONST_LOG2R * p7_FLogsum(0.0, log(bg->omega) + hitlist->hit[h]->dcl[d].domcorrection), /* print in units of bits */
             exp(hitlist->hit[h]->dcl[d].lnP) * (double) pli->Z,
             hitlist->hit[h]->nexpected,
             hitlist->hit[h]->nreported,
             namew,
             hitlist->hit[h]->name,
             hitlist->hit[h]->desc);
    }

  /* Done. */
  p7_tophits_Destroy(hitlist);
  p7_pipeline_Destroy(pli);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7PIPELINE_EXAMPLE*/
/*----------- end, search mode (seq db) example -----------------*/




/*****************************************************************
 * 8. Example 2: "scan mode" in an HMM db
 *****************************************************************/
#ifdef p7PIPELINE_EXAMPLE2
/* ./pipeline_example2 <hmmdb> <sqfile>
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles   reqs   incomp                             help                                                  docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,                          "show brief help on version and usage",                         0 },
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "E-value cutoff for reporting significant sequence hits",       0 },
  { "-T",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "bit score cutoff for reporting significant sequence hits",     0 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                          "set # of comparisons done, for E-value calculation",           0 },
  { "--domE",       eslARG_REAL,"1000.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "E-value cutoff for reporting individual domains",              0 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "bit score cutoff for reporting individual domains",            0 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                          "set # of significant seqs, for domain E-value calculation",    0 },
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--seqE,--seqT,--domE,--domT", "use GA gathering threshold bit score cutoffs in <hmmfile>",    0 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--seqE,--seqT,--domE,--domT", "use NC noise threshold bit score cutoffs in <hmmfile>",        0 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--seqE,--seqT,--domE,--domT", "use TC trusted threshold bit score cutoffs in <hmmfile>",      0 },
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--F1,--F2,--F3",               "Turn all heuristic filters off (less speed, more power)",      0 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,      NULL,  NULL, "--max",                        "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             0 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,      NULL,  NULL, "--max",                        "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             0 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,      NULL,  NULL, "--max",                        "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             0 },
  { "--nobias",     eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL, "--max",                        "turn off composition bias filter",                             0 },
  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL,  NULL,                          "turn off biased composition score corrections",                0 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",    NULL,  NULL,  NULL,                          "set RNG seed to <n> (if 0: one-time arbitrary seed)",          0 },
  { "--acc",        eslARG_NONE,  FALSE,  NULL, NULL,      NULL,  NULL,  NULL,                          "output target accessions instead of names if possible",        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of using acceleration pipeline in scan mode (HMM targets)";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char         *hmmfile = esl_opt_GetArg(go, 1);
  char         *seqfile = esl_opt_GetArg(go, 2);
  int           format  = eslSQFILE_FASTA;
  P7_HMMFILE   *hfp     = NULL;
  ESL_ALPHABET *abc     = NULL;
  P7_BG        *bg      = NULL;
  P7_OPROFILE  *om      = NULL;
  ESL_SQFILE   *sqfp    = NULL;
  ESL_SQ       *sq      = NULL;
  P7_PIPELINE  *pli     = NULL;
  P7_TOPHITS   *hitlist = p7_tophits_Create(p7_TOPHITS_DEFAULT_INIT_ALLOC);
  int           h,d,namew;

  /* Open a sequence file, read one seq from it.
   * Convert to digital later, after 1st HMM is input and abc becomes known 
   */
  sq = esl_sq_Create();
  if (esl_sqfile_Open(seqfile, format, NULL, &sqfp) != eslOK) p7_Fail("Failed to open sequence file %s\n", seqfile);
  if (esl_sqio_Read(sqfp, sq)                       != eslOK) p7_Fail("Failed to read sequence from %s\n", seqfile);
  esl_sqfile_Close(sqfp);

  /* Open the HMM db */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);

  /* Create a pipeline for the query sequence in scan mode */
  pli      = p7_pipeline_Create(go, 100, sq->n, FALSE, p7_SCAN_MODELS);
  p7_pipeline_NewSeq(pli, sq);
   
  /* Some additional config of the pipeline specific to scan mode */
  pli->hfp = hfp;
  if (! pli->Z_is_fixed && hfp->is_pressed) { pli->Z_is_fixed = TRUE; pli->Z = hfp->ssi->nprimary; }

  /* Read (partial) of each HMM in file */
  while (p7_oprofile_ReadMSV(hfp, &abc, &om) == eslOK) 
    {
      /* One time only initiazation after abc becomes known */
      if (bg == NULL) 
    {
      bg = p7_bg_Create(abc);
      if (esl_sq_Digitize(abc, sq) != eslOK) p7_Die("alphabet mismatch");
    }
      p7_pipeline_NewModel(pli, om, bg);
      p7_bg_SetLength(bg, sq->n); /* SetLength() call MUST follow NewModel() call, because NewModel() resets the filter HMM, including its default expected length; see bug #h85 */
      p7_oprofile_ReconfigLength(om, sq->n);

      p7_Pipeline(pli, om, bg, sq, hitlist);
      
      p7_oprofile_Destroy(om);
      p7_pipeline_Reuse(pli);
    } 

  /* Print the results. 
   * This example is a stripped version of hmmsearch's tabular output.
   */
  p7_tophits_SortBySortkey(hitlist);
  namew = ESL_MAX(8, p7_tophits_GetMaxNameLength(hitlist));
  for (h = 0; h < hitlist->N; h++)
    {
      d    = hitlist->hit[h]->best_domain;

      printf("%10.2g %7.1f %6.1f  %7.1f %6.1f %10.2g  %6.1f %5d  %-*s %s\n",
             exp(hitlist->hit[h]->lnP) * (double) pli->Z,
             hitlist->hit[h]->score,
             hitlist->hit[h]->pre_score - hitlist->hit[h]->score, /* bias correction */
             hitlist->hit[h]->dcl[d].bitscore,
             eslCONST_LOG2R * p7_FLogsum(0.0, log(bg->omega) + hitlist->hit[h]->dcl[d].domcorrection), /* print in units of BITS */
             exp(hitlist->hit[h]->dcl[d].lnP) * (double) pli->Z,
             hitlist->hit[h]->nexpected,
             hitlist->hit[h]->nreported,
             namew,
             hitlist->hit[h]->name,
             hitlist->hit[h]->desc);
    }

  /* Done. */
  p7_tophits_Destroy(hitlist);
  p7_pipeline_Destroy(pli);
  esl_sq_Destroy(sq);
  p7_hmmfile_Close(hfp);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7PIPELINE_EXAMPLE2*/
/*--------------- end, scan mode (HMM db) example ---------------*/



/*****************************************************************
 * @LICENSE@
 *
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
