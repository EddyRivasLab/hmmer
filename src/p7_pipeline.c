/* H3's accelerated seq/profile comparison pipeline
 *  
 * Contents:
 *   1. P7_PIPELINE: allocation, initialization, destruction
 *   2. Pipeline API
 *   3. Example 1: search mode (in a sequence db)
 *   4. Example 2: scan mode (in an HMM db)
 *   5. Copyright and license information
 * 
 * SRE, Fri Dec  5 10:09:39 2008 [Janelia] [BSG3, Bear McCreary]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

/*****************************************************************
 * 1. The P7_PIPELINE object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  p7_pipeline_Create()
 * Synopsis:  Create a new accelerated comparison pipeline.
 * Incept:    SRE, Fri Dec  5 10:11:31 2008 [Janelia]
 *
 * Purpose:   Given an application configuration structure <go>
 *            containing certain standardized options (described
 *            below), some initial guesses at the model size <M_hint>
 *            and sequence length <L_hint> that will be processed,
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
 * Returns:   ptr to new <P7_PIPELINE> object on success. Caller frees this
 *            with <p7_pipeline_Destroy()>.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_PIPELINE *
p7_pipeline_Create(ESL_GETOPTS *go, int M_hint, int L_hint, enum p7_pipemodes_e mode)
{
  P7_PIPELINE *pli  = NULL;
  int          seed = esl_opt_GetInteger(go, "--seed");
  int          status;

  ESL_ALLOC(pli, sizeof(P7_PIPELINE));

  if ((pli->fwd = p7_omx_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR;
  if ((pli->bck = p7_omx_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR;	
  if ((pli->oxf = p7_omx_Create(M_hint, 0,      L_hint)) == NULL) goto ERROR;
  if ((pli->oxb = p7_omx_Create(M_hint, 0,      L_hint)) == NULL) goto ERROR;     

  /* Normally, we reinitialize the RNG to the original seed every time we're
   * about to collect a stochastic trace ensemble. This eliminates run-to-run
   * variability. As a special case, if seed==0, we choose an arbitrary one-time 
   * seed: time() sets the seed, and we turn off the reinitialization.
   */
  pli->r                  =  esl_randomness_CreateFast(seed);
  pli->do_reseeding       = (seed == 0) ? FALSE : TRUE;
  pli->ddef               = p7_domaindef_Create(pli->r);
  pli->ddef->do_reseeding = pli->do_reseeding;

  /* Configure reporting thresholds */
  pli->by_E            = TRUE;
  pli->E               = esl_opt_GetReal(go, "-E");
  pli->T               = 0.0;
  pli->dom_by_E        = TRUE;
  pli->domE            = esl_opt_GetReal(go, "--domE");
  pli->domT            = 0.0;
  pli->use_bit_cutoffs = FALSE;
  if (esl_opt_IsOn(go, "-T")) 
    {
      pli->T    = esl_opt_GetReal(go, "-T"); 
      pli->by_E = FALSE;
    } 
  if (esl_opt_IsOn(go, "--domT")) 
    {
      pli->domT     = esl_opt_GetReal(go, "--domT"); 
      pli->dom_by_E = FALSE;
    }


  /* Configure inclusion thresholds */
  pli->inc_by_E           = TRUE;
  pli->incE               = esl_opt_GetReal(go, "--incE");
  pli->incT               = 0.0;
  pli->incdom_by_E        = TRUE;
  pli->incdomE            = esl_opt_GetReal(go, "--incdomE");
  pli->incdomT            = 0.0;
  if (esl_opt_IsOn(go, "--incT")) 
    {
      pli->incT     = esl_opt_GetReal(go, "--incT"); 
      pli->inc_by_E = FALSE;
    } 
  if (esl_opt_IsOn(go, "--incdomT")) 
    {
      pli->incdomT     = esl_opt_GetReal(go, "--incdomT"); 
      pli->incdom_by_E = FALSE;
    }


  /* Configure for one of the model-specific thresholding options */
  if (esl_opt_GetBoolean(go, "--cut_ga"))
    {
      pli->T        = pli->domT        = 0.0;
      pli->by_E     = pli->dom_by_E    = FALSE;
      pli->incT     = pli->incdomT     = 0.0;
      pli->inc_by_E = pli->incdom_by_E = FALSE;
      pli->use_bit_cutoffs = p7H_GA;
    }
  if (esl_opt_GetBoolean(go, "--cut_nc"))
    {
      pli->T        = pli->domT        = 0.0;
      pli->by_E     = pli->dom_by_E    = FALSE;
      pli->incT     = pli->incdomT     = 0.0;
      pli->inc_by_E = pli->incdom_by_E = FALSE;
      pli->use_bit_cutoffs = p7H_NC;
    }
  if (esl_opt_GetBoolean(go, "--cut_tc"))
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
  if (esl_opt_IsOn(go, "-Z")) 
    {
      pli->Z_setby = p7_ZSETBY_OPTION;
      pli->Z       = esl_opt_GetReal(go, "-Z");
    }
  if (esl_opt_IsOn(go, "--domZ")) 
    {
      pli->domZ_setby = p7_ZSETBY_OPTION;
      pli->domZ       = esl_opt_GetReal(go, "--domZ");
    }


  /* Configure acceleration pipeline thresholds */
  pli->do_max        = FALSE;
  pli->do_biasfilter = TRUE;
  pli->do_null2      = TRUE;
  pli->F1     = ESL_MIN(1.0, esl_opt_GetReal(go, "--F1"));
  pli->F2     = ESL_MIN(1.0, esl_opt_GetReal(go, "--F2"));
  pli->F3     = ESL_MIN(1.0, esl_opt_GetReal(go, "--F3"));
  if (esl_opt_GetBoolean(go, "--max")) 
    {
      pli->do_max        = TRUE;
      pli->do_biasfilter = FALSE;
      pli->F1 = pli->F2 = pli->F3 = 1.0; 
    }
  if (esl_opt_GetBoolean(go, "--nonull2")) pli->do_null2      = FALSE;
  if (esl_opt_GetBoolean(go, "--nobias"))  pli->do_biasfilter = FALSE;
  

  /* Accounting as we collect results */
  pli->nmodels     = 0;
  pli->nseqs       = 0;
  pli->nres        = 0;
  pli->nnodes      = 0;
  pli->n_past_msv  = 0;
  pli->n_past_bias = 0;
  pli->n_past_vit  = 0;
  pli->n_past_fwd  = 0;
  
  pli->mode            = mode;
  pli->show_accessions = (esl_opt_GetBoolean(go, "--acc")   ? TRUE  : FALSE);
  pli->show_alignments = (esl_opt_GetBoolean(go, "--noali") ? FALSE : TRUE);
  pli->hfp             = NULL;
  pli->errbuf[0]       = '\0';

  return pli;

 ERROR:
  p7_pipeline_Destroy(pli);
  return NULL;
}


/* Function:  p7_pipeline_Reuse()
 * Synopsis:  Reuse a pipeline for next target.
 * Incept:    SRE, Fri Dec  5 10:31:36 2008 [Janelia]
 *
 * Purpose:   Reuse <pli> for next target sequence (search mode)
 *            or model (scan mode). 
 *            
 *            May eventually need to distinguish from reusing pipeline
 *            for next query, but we're not really focused on multiquery
 *            use of hmmscan/hmmsearch/phmmer for the moment.
 */
int
p7_pipeline_Reuse(P7_PIPELINE *pli)
{
  p7_omx_Reuse(pli->oxf);
  p7_omx_Reuse(pli->oxb);
  p7_omx_Reuse(pli->fwd);
  p7_omx_Reuse(pli->bck);
  p7_domaindef_Reuse(pli->ddef);
  return eslOK;
}



/* Function:  p7_pipeline_Destroy()
 * Synopsis:  Free a <P7_PIPELINE> object.
 * Incept:    SRE, Fri Dec  5 10:30:23 2008 [Janelia]
 *
 * Purpose:   Free a <P7_PIPELINE> object.
 */
void
p7_pipeline_Destroy(P7_PIPELINE *pli)
{
  if (pli == NULL) return;
  
  p7_omx_Destroy(pli->oxf);
  p7_omx_Destroy(pli->oxb);
  p7_omx_Destroy(pli->fwd);
  p7_omx_Destroy(pli->bck);
  esl_randomness_Destroy(pli->r);
  p7_domaindef_Destroy(pli->ddef);
  free(pli);
}
/*---------------- end, P7_PIPELINE object ----------------------*/





/*****************************************************************
 * 2. The pipeline API.
 *****************************************************************/

/* Function:  p7_pli_TargetReportable
 * Synopsis:  Returns TRUE if target score meets reporting threshold.
 * Incept:    SRE, Tue Dec  9 08:57:26 2008 [Janelia]
 *
 * Purpose:   Returns <TRUE> if the bit score <score> and/or 
 *            P-value <Pval> meeds per-target reporting thresholds 
 *            for the processing pipeline.
 */
int
p7_pli_TargetReportable(P7_PIPELINE *pli, float score, float Pval)
{
  if      (  pli->by_E   && Pval * pli->Z <= pli->E) return TRUE;
  else if (! pli->by_E   && score         >= pli->T) return TRUE;
  else return FALSE;
}

/* Function:  p7_pli_DomainReportable
 * Synopsis:  Returns TRUE if domain score meets reporting threshold. 
 * Incept:    SRE, Tue Dec  9 09:01:01 2008 [Janelia]
 *
 * Purpose:   Returns <TRUE> if the bit score <score> and/or 
 *            P-value <Pval> meets per-domain reporting thresholds 
 *            for the processing pipeline.
 */
int
p7_pli_DomainReportable(P7_PIPELINE *pli, float dom_score, float Pval)
{
  if      (  pli->dom_by_E   && Pval * pli->domZ <= pli->domE) return TRUE;
  else if (! pli->dom_by_E   && dom_score        >= pli->domT) return TRUE;
  else return FALSE;
}

/* Function:  p7_pli_TargetIncludable()
 * Synopsis:  Returns TRUE if target score meets inclusion threshold.
 * Incept:    SRE, Fri Jan 16 11:18:08 2009 [Janelia]
 */
int
p7_pli_TargetIncludable(P7_PIPELINE *pli, float score, float Pval)
{
  if      (  pli->inc_by_E   && Pval * pli->Z <= pli->incE) return TRUE;
  else if (! pli->inc_by_E   && score         >= pli->incT) return TRUE;
  else return FALSE;
}

/* Function:  p7_pli_DomainIncludable()
 * Synopsis:  Returns TRUE if domain score meets inclusion threshold.
 * Incept:    SRE, Fri Jan 16 11:20:38 2009 [Janelia]
 */
int
p7_pli_DomainIncludable(P7_PIPELINE *pli, float dom_score, float Pval)
{
  if      (  pli->incdom_by_E   && Pval * pli->domZ <= pli->incdomE) return TRUE;
  else if (! pli->incdom_by_E   && dom_score        >= pli->incdomT) return TRUE;
  else return FALSE;
}




/* Function:  p7_pli_NewModel()
 * Synopsis:  Prepare pipeline for a new model (target or query)
 * Incept:    SRE, Fri Dec  5 10:35:37 2008 [Janelia]
 *
 * Purpose:   Caller has a new model <om>. Prepare the pipeline <pli>
 *            to receive this model as either a query or a target.
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
 */
int
p7_pli_NewModel(P7_PIPELINE *pli, const P7_OPROFILE *om, P7_BG *bg)
{
  int status = eslOK;

  pli->nmodels++;
  pli->nnodes += om->M;
  if (pli->Z_setby == p7_ZSETBY_NTARGETS && pli->mode == p7_SCAN_MODELS) pli->Z = pli->nmodels;

  if (pli->do_biasfilter) p7_bg_SetFilter(bg, om->M, om->compo);

  if (pli->mode == p7_SEARCH_SEQS) 
    status = p7_pli_NewModelThresholds(pli, om);

  return status;
}

/* Function:  p7_pli_NewModelThresholds()
 * Synopsis:  Set reporting and inclusion bit score thresholds on a new model.
 * Incept:    SRE, Sat Oct 17 12:07:43 2009 [Janelia]
 *
 * Purpose:   Set the bit score thresholds on a new model, if we're 
 *            using Pfam GA, TC, or NC cutoffs for reporting or
 *            inclusion.
 *            
 *            In a "search" pipeline, this only needs to be done once
 *            per query model, so <p7_pli_NewModelThresholds()> gets 
 *            called by <p7_pli_NewModel()>.
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
p7_pli_NewModelThresholds(P7_PIPELINE *pli, const P7_OPROFILE *om)
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


/* Function:  p7_pli_NewSeq()
 * Synopsis:  Prepare pipeline for a new sequence (target or query)
 * Incept:    SRE, Fri Dec  5 10:57:15 2008 [Janelia]
 *
 * Purpose:   Caller has a new sequence <sq>. Prepare the pipeline <pli>
 *            to receive this model as either a query or a target.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_pli_NewSeq(P7_PIPELINE *pli, const ESL_SQ *sq)
{
  pli->nseqs++;
  pli->nres += sq->n;
  if (pli->Z_setby == p7_ZSETBY_NTARGETS && pli->mode == p7_SEARCH_SEQS) pli->Z = pli->nseqs;
  return eslOK;
}

/* Function:  p7_pipeline_Merge()
 * Synopsis:  Merge the pipeline statistics
 * Incept:    
 *
 * Purpose:   Caller has a new model <om>. Prepare the pipeline <pli>
 *            to receive this model as either a query or a target.
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
 */
int
p7_pipeline_Merge(P7_PIPELINE *p1, P7_PIPELINE *p2)
{
  /* if we are searching a sequence database, we need to keep track of the
   * number of sequences and residues processed.
   */
  if (p1->mode == p7_SEARCH_SEQS)
    {
      p1->nseqs   += p2->nseqs;
      p1->nres    += p2->nres;
    }
  else
    {
      p1->nmodels += p2->nmodels;
      p1->nnodes  += p2->nnodes;
    }

  p1->n_past_msv  += p2->n_past_msv;
  p1->n_past_bias += p2->n_past_bias;
  p1->n_past_vit  += p2->n_past_vit;
  p1->n_past_fwd  += p2->n_past_fwd;

  if (p1->Z_setby == p7_ZSETBY_NTARGETS)
    {
      p1->Z += (p1->mode == p7_SCAN_MODELS) ? p2->nmodels : p2->nseqs;
    }
  else
    {
      p1->Z = p2->Z;
    }

  return eslOK;
}

/* Function:  p7_Pipeline()
 * Synopsis:  HMMER3's accelerated seq/profile comparison pipeline.
 * Incept:    SRE, Thu Dec  4 17:17:01 2008 [Janelia]
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
 * Xref:      J4/25.
 */
int
p7_Pipeline(P7_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, P7_TOPHITS *hitlist)
{
  P7_HIT          *hit     = NULL;     /* ptr to the current hit output data      */
  float            usc, vfsc, fwdsc;   /* filter scores                           */
  float            filtersc;           /* HMM null filter score                   */
  float            nullsc;             /* null model score                        */
  float            seqbias;  
  float            seq_score;          /* the corrected per-seq bit score */
  float            sum_score;	       /* the corrected reconstruction score for the seq */
  float            pre_score, pre2_score; /* uncorrected bit scores for seq */
  double           P;		       /* P-value of a hit */
  int              Ld;		       /* # of residues in envelopes */
  int              d;
  int              status;
  
  if (sq->n == 0) return eslOK;	/* silently skip length 0 seqs; they'd cause us all sorts of weird problems */

  p7_omx_GrowTo(pli->oxf, om->M, 0, sq->n);    /* expand the one-row omx if needed */

  /* Base null model score (we could calculate this in NewSeq(), for a scan pipeline) */
  p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);


  /* First level filter: the MSV filter, multihit with <om> */
  int* window_starts;
  int* window_ends;
  float invsurv = esl_gumbel_invsurv(pli->F1, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  float pthresh = nullsc  + (invsurv * eslCONST_LOG2);

  //int ret = p7_SSVFilter_longseq(sq->dsq, sq->n, om, pli->oxf, pthresh, &window_starts, &window_ends);
  int ret = p7_MSVFilter_longseq(sq->dsq, sq->n, om, pli->oxf, pthresh, &window_starts, &window_ends);
  if (ret == eslOK) return eslOK;



 // old msv filter
/*
  p7_MSVFilter(sq->dsq, sq->n, om, pli->oxf, &usc);
  seq_score = (usc - nullsc) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  if (P > pli->F1) return eslOK;
*/
  pli->n_past_msv++;




  /* biased composition HMM filtering */
  /*
  if (pli->do_biasfilter)
	{
	  p7_bg_FilterScore(bg, sq->dsq, sq->n, &filtersc);

	  //have to run msvfilter again.  This time, it should just be on the windows passed in from above
	  float pthresh2 = filtersc  + (invsurv * eslCONST_LOG2);
	  int ret = p7_MSVFilter_longseq(sq->dsq, sq->n, om, pli->oxf, pthresh2, &window_starts, &window_ends);
	  if (ret == eslOK) return eslOK;


//	  P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
//	  if (P > pli->F1) return eslOK;

	}
  else filtersc = nullsc;
  pli->n_past_bias++;
*/

  /* In scan mode, if it passes the MSV filter, read the rest of the profile */
  if (pli->hfp) {
	p7_oprofile_ReadRest(pli->hfp, om);
	p7_oprofile_ReconfigRestLength(om, sq->n);
	if ((status = p7_pli_NewModelThresholds(pli, om)) != eslOK) return status; /* pli->errbuf has err msg set */
  }


  /* Second level filter: ViterbiFilter(), multihit with <om> */
  //if (P > pli->F2 )
	{
	  p7_ViterbiFilter(sq->dsq, sq->n, om, pli->oxf, &vfsc);
	  seq_score = (vfsc-filtersc) / eslCONST_LOG2;
	  P  = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
	  if (P > pli->F2) return eslOK;
	}
  pli->n_past_vit++;


  /* Parse it with Forward and obtain its real Forward score. */
  p7_ForwardParser(sq->dsq, sq->n, om, pli->oxf, &fwdsc);
  seq_score = (fwdsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  if (P > pli->F3) return eslOK;
  pli->n_past_fwd++;

  /* ok, it's for real. Now a Backwards parser pass, and hand it to domain definition workflow */
  p7_omx_GrowTo(pli->oxb, om->M, 0, sq->n);
  p7_BackwardParser(sq->dsq, sq->n, om, pli->oxf, pli->oxb, NULL);

  status = p7_domaindef_ByPosteriorHeuristics(sq, om, pli->oxf, pli->oxb, pli->fwd, pli->bck, pli->ddef);
  if (status != eslOK) ESL_FAIL(status, pli->errbuf, "domain definition workflow failure"); /* eslERANGE can happen */
  if (pli->ddef->nregions   == 0) return eslOK; /* score passed threshold but there's no discrete domains here       */
  if (pli->ddef->nenvelopes == 0) return eslOK; /* rarer: region was found, stochastic clustered, no envelopes found */


  /* Calculate the null2-corrected per-seq score */
  if (pli->do_null2)
    {
      seqbias = esl_vec_FSum(pli->ddef->n2sc, sq->n+1);
      seqbias = p7_FLogsum(0.0, log(bg->omega) + seqbias);
    }
  else seqbias = 0.0;
  pre_score =  (fwdsc - nullsc) / eslCONST_LOG2; 
  seq_score =  (fwdsc - (nullsc + seqbias)) / eslCONST_LOG2;

  
  /* Calculate the "reconstruction score": estimated
   * per-sequence score as sum of individual domains,
   * discounting domains that aren't significant after they're
   * null-corrected.
   */
  sum_score = 0.0f;
  seqbias   = 0.0f;
  Ld        = 0;
  if (pli->do_null2) 
    {
      for (d = 0; d < pli->ddef->ndom; d++) 
	{
	  if (pli->ddef->dcl[d].envsc - pli->ddef->dcl[d].domcorrection > 0.0) 
	    {
	      sum_score += pli->ddef->dcl[d].envsc;
	      Ld        += pli->ddef->dcl[d].jenv  - pli->ddef->dcl[d].ienv + 1;
	      seqbias   += pli->ddef->dcl[d].domcorrection;
	    }
	}
      seqbias = p7_FLogsum(0.0, log(bg->omega) + seqbias);
    }
  else 
    {
      for (d = 0; d < pli->ddef->ndom; d++) 
	{
	  if (pli->ddef->dcl[d].envsc > 0.0) 
	    {
	      sum_score += pli->ddef->dcl[d].envsc;
	      Ld        += pli->ddef->dcl[d].jenv  - pli->ddef->dcl[d].ienv + 1;
	    }
	}
      seqbias = 0.0;
    }    
  sum_score += (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); 
  pre2_score = (sum_score - nullsc) / eslCONST_LOG2;
  sum_score  = (sum_score - (nullsc + seqbias)) / eslCONST_LOG2;

  /* A special case: let sum_score override the seq_score when it's better, and it includes at least 1 domain */
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
  P =  esl_exp_surv (seq_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
  if (p7_pli_TargetReportable(pli, seq_score, P))
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
      hit->ndom       = pli->ddef->ndom;
      hit->nexpected  = pli->ddef->nexpected;
      hit->nregions   = pli->ddef->nregions;
      hit->nclustered = pli->ddef->nclustered;
      hit->noverlaps  = pli->ddef->noverlaps;
      hit->nenvelopes = pli->ddef->nenvelopes;

      hit->pre_score  = pre_score;
      hit->pre_pvalue = esl_exp_surv (hit->pre_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

      hit->score      = seq_score;
      hit->pvalue     = P;
      hit->sortkey    = pli->inc_by_E ? -log(P) : seq_score; /* per-seq output sorts on bit score if inclusion is by score  */

      hit->sum_score  = sum_score;
      hit->sum_pvalue = esl_exp_surv (hit->sum_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

      /* Transfer all domain coordinates (unthresholded for
       * now) with their alignment displays to the hit list,
       * associated with the sequence. Domain reporting will
       * be thresholded after complete hit list is collected,
       * because we probably need to know # of significant
       * hits found to set domZ, and thence threshold and
       * count reported domains.
       */
      hit->dcl         = pli->ddef->dcl;
      pli->ddef->dcl   = NULL;
      hit->best_domain = 0;
      for (d = 0; d < hit->ndom; d++)
	{
	  Ld = hit->dcl[d].jenv - hit->dcl[d].ienv + 1;
	  hit->dcl[d].bitscore = hit->dcl[d].envsc + (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); 
	  hit->dcl[d].dombias  = (pli->do_null2 ? p7_FLogsum(0.0, log(bg->omega) + hit->dcl[d].domcorrection) : 0.0);
	  hit->dcl[d].bitscore = (hit->dcl[d].bitscore - (nullsc + hit->dcl[d].dombias)) / eslCONST_LOG2;
	  hit->dcl[d].pvalue   = esl_exp_surv (hit->dcl[d].bitscore,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
	  
	  if (hit->dcl[d].bitscore > hit->dcl[hit->best_domain].bitscore) hit->best_domain = d;
	}

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
	  if (p7_pli_TargetReportable(pli, hit->score, hit->pvalue))
	    {
	      hit->flags |= p7_IS_REPORTED;
	      if (p7_pli_TargetIncludable(pli, hit->score, hit->pvalue))
		hit->flags |= p7_IS_INCLUDED;
	    }

	  for (d = 0; d < hit->ndom; d++)
	    {
	      if (p7_pli_DomainReportable(pli, hit->dcl[d].bitscore, hit->dcl[d].pvalue))
		{
		  hit->dcl[d].is_reported = TRUE;
		  if (p7_pli_DomainIncludable(pli, hit->dcl[d].bitscore, hit->dcl[d].pvalue))
		    hit->dcl[d].is_included = TRUE;
		}
	    }	  
	}
    }

  return eslOK;
}


/* Function:  p7_pli_Statistics()
 * Synopsis:  Final statistics output from a processing pipeline.
 * Incept:    SRE, Tue Dec  9 10:19:45 2008 [Janelia]
 *
 * Purpose:   Print a standardized report of the internal statistics of
 *            a finished processing pipeline <pli> to stream <ofp>.
 *            
 *            If stopped, non-<NULL> stopwatch <w> is provided for a
 *            stopwatch that was timing the pipeline, then the report
 *            includes timing information.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_pli_Statistics(FILE *ofp, P7_PIPELINE *pli, ESL_STOPWATCH *w)
{
  double ntargets; 

  fprintf(ofp, "Internal pipeline statistics summary:\n");
  fprintf(ofp, "-------------------------------------\n");
  if (pli->mode == p7_SEARCH_SEQS) {
    fprintf(ofp, "Query model(s):              %15" PRId64 "  (%" PRId64 " nodes)\n",     pli->nmodels, pli->nnodes);
    fprintf(ofp, "Target sequences:            %15" PRId64 "  (%" PRId64 " residues)\n",  pli->nseqs,   pli->nres);
    ntargets = pli->nseqs;
  } else {
    fprintf(ofp, "Query sequence(s):           %15" PRId64 "  (%" PRId64 " residues)\n",  pli->nseqs,   pli->nres);
    fprintf(ofp, "Target model(s):             %15" PRId64 "  (%" PRId64 " nodes)\n",     pli->nmodels, pli->nnodes);
    ntargets = pli->nmodels;
  }

  fprintf(ofp, "Passed MSV filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n", 
	  pli->n_past_msv,
	  (double) pli->n_past_msv / ntargets,
	  pli->F1 * ntargets,
	  pli->F1);

  fprintf(ofp, "Passed bias filter:          %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n", 
	  pli->n_past_bias,
	  (double) pli->n_past_bias / ntargets,
	  pli->F1 * ntargets,
	  pli->F1);

  fprintf(ofp, "Passed Vit filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",   
	  pli->n_past_vit,
	  (double) pli->n_past_vit / ntargets,
	  pli->F2 * ntargets,
	  pli->F2);

  fprintf(ofp, "Passed Fwd filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",         
	  pli->n_past_fwd, 
	  (double) pli->n_past_fwd / ntargets,
	  pli->F3 * ntargets,
	  pli->F3);	  

  fprintf(ofp, "Initial search space (Z):    %15.0f  %s\n", pli->Z,    pli->Z_setby    == p7_ZSETBY_OPTION ? "[as set by --Z on cmdline]"    : "[actual number of targets]"); 
  fprintf(ofp, "Domain search space  (domZ): %15.0f  %s\n", pli->domZ, pli->domZ_setby == p7_ZSETBY_OPTION ? "[as set by --domZ on cmdline]" : "[number of targets reported over threshold]"); 

  if (w != NULL) {
    esl_stopwatch_Display(ofp, w, "# CPU time: ");
    fprintf(ofp, "# Mc/sec: %.2f\n", 
	    (double) pli->nres * (double) pli->nnodes / (w->elapsed * 1.0e6));
  }

  return eslOK;
}
/*------------------- end, pipeline API -------------------------*/




/*****************************************************************
 * 3. Example 1: "search mode" in a sequence db
 *****************************************************************/

#ifdef p7PIPELINE_EXAMPLE
/* gcc -o pipeline_example -g -Wall -I../easel -L../easel -I. -L. -Dp7PIPELINE_EXAMPLE p7_pipeline.c -lhmmer -leasel -lm
 * ./pipeline_example <hmmfile> <sqfile>
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
  ESL_GETOPTS  *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
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

  /* Don't forget this. Null2 corrections need FLogsum() */
  p7_FLogsumInit();

  /* Read in one HMM */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Open a sequence file */
  if (esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp) != eslOK) p7_Fail("Failed to open sequence file %s\n", seqfile);
  sq = esl_sq_CreateDigital(abc);

  /* Create a pipeline and a top hits list */
  pli     = p7_pipeline_Create(go, hmm->M, 400, p7_SEARCH_SEQS);
  hitlist = p7_tophits_Create();

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, 400, p7_LOCAL);
  p7_oprofile_Convert(gm, om);     /* <om> is now p7_LOCAL, multihit */
  p7_pli_NewModel(pli, om, bg);

  /* Run each target sequence through the pipeline */
  while (esl_sqio_Read(sqfp, sq) == eslOK)
    { 
      p7_pli_NewSeq(pli, sq);
      p7_bg_SetLength(bg, sq->n);
      p7_oprofile_ReconfigLength(om, sq->n);
  
      p7_Pipeline(pli, om, bg, sq, hitlist);

      esl_sq_Reuse(sq);
      p7_pipeline_Reuse(pli);
    }

  /* Print the results. 
   * This example is a stripped version of hmmsearch's tabular output.
   */
  p7_tophits_Sort(hitlist);
  namew = ESL_MAX(8, p7_tophits_GetMaxNameLength(hitlist));
  for (h = 0; h < hitlist->N; h++)
    {
      d    = hitlist->hit[h]->best_domain;

      printf("%10.2g %7.1f %6.1f  %7.1f %6.1f %10.2g  %6.1f %5d  %-*s %s\n",
	     hitlist->hit[h]->pvalue * (double) pli->Z,
	     hitlist->hit[h]->score,
	     hitlist->hit[h]->pre_score - hitlist->hit[h]->score, /* bias correction */
	     hitlist->hit[h]->dcl[d].bitscore,
	     p7_FLogsum(0.0, log(bg->omega) + hitlist->hit[h]->dcl[d].domcorrection),
	     hitlist->hit[h]->dcl[d].pvalue * (double) pli->Z,
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
 * 4. Example 2: "scan mode" in an HMM db
 *****************************************************************/
#ifdef p7PIPELINE_EXAMPLE2
/* gcc -o pipeline_example2 -g -Wall -I../easel -L../easel -I. -L. -Dp7PIPELINE_EXAMPLE2 p7_pipeline.c -lhmmer -leasel -lm
 * ./pipeline_example2 <hmmdb> <sqfile>
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
  ESL_GETOPTS  *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
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
  P7_TOPHITS   *hitlist = p7_tophits_Create();
  int           h,d,namew;

  /* Don't forget this. Null2 corrections need FLogsum() */
  p7_FLogsumInit();

  /* Open a sequence file, read one seq from it.
   * Convert to digital later, after 1st HMM is input and abc becomes known 
   */
  sq = esl_sq_Create();
  if (esl_sqfile_Open(seqfile, format, NULL, &sqfp) != eslOK) p7_Fail("Failed to open sequence file %s\n", seqfile);
  if (esl_sqio_Read(sqfp, sq)                       != eslOK) p7_Fail("Failed to read sequence from %s\n", seqfile);
  esl_sqfile_Close(sqfp);

  /* Open the HMM db */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);

  /* Create a pipeline for the query sequence in scan mode */
  pli      = p7_pipeline_Create(go, 100, sq->n, p7_SCAN_MODELS);
  p7_pli_NewSeq(pli, sq);
   
  /* Some additional config of the pipeline specific to scan mode */
  pli->hfp = hfp;
  if (! pli->Z_is_fixed && hfp->is_pressed) { pli->Z_is_fixed = TRUE; pli->Z = hfp->ssi->nprimary; }

  /* Read (partial) of each HMM in file */
  while (p7_oprofile_ReadMSV(hfp, &abc, &om) == eslOK) 
    {
      /* One time only initialization after abc becomes known */
      if (bg == NULL) 
	{
	  bg = p7_bg_Create(abc);
	  if (esl_sq_Digitize(abc, sq) != eslOK) p7_Die("alphabet mismatch");
	  p7_bg_SetLength(bg, sq->n);
	}
      p7_pli_NewModel(pli, om, bg);
      p7_oprofile_ReconfigLength(om, sq->n);

      p7_Pipeline(pli, om, bg, sq, hitlist);
      
      p7_oprofile_Destroy(om);
      p7_pipeline_Reuse(pli);
    } 

  /* Print the results. 
   * This example is a stripped version of hmmsearch's tabular output.
   */
  p7_tophits_Sort(hitlist);
  namew = ESL_MAX(8, p7_tophits_GetMaxNameLength(hitlist));
  for (h = 0; h < hitlist->N; h++)
    {
      d    = hitlist->hit[h]->best_domain;

      printf("%10.2g %7.1f %6.1f  %7.1f %6.1f %10.2g  %6.1f %5d  %-*s %s\n",
	     hitlist->hit[h]->pvalue * (double) pli->Z,
	     hitlist->hit[h]->score,
	     hitlist->hit[h]->pre_score - hitlist->hit[h]->score, /* bias correction */
	     hitlist->hit[h]->dcl[d].bitscore,
	     p7_FLogsum(0.0, log(bg->omega) + hitlist->hit[h]->dcl[d].domcorrection),
	     hitlist->hit[h]->dcl[d].pvalue * (double) pli->Z,
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
 *****************************************************************/
