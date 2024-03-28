/* H3's accelerated seq/profile comparison pipeline
 *  
 * Contents:
 *   1. P7_PIPELINE: allocation, initialization, destruction
 *   2. Pipeline API
 */
#include <p7_config.h> 
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

#include "easel.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "esl_sqio.h" //!!!!DEBUG

#include "e2_config.h"
#include "e1_rate.h"
#include "evohmmer.h"
#include "minimize.h"
#include "p7_evopipeline.h"
#include "ratematrix.h"

#define OPT_TIME_STEP 0.01
#define TMAX          5.0

/* Struct used to pass a collection of useful temporary objects around
 * within the LongTarget functions
 *  */
typedef struct {
  ESL_SQ           *tmpseq; // - a new or reused digital sequence object used for p7_alidisplay_Create() call
  P7_BG            *bg;
  P7_OPROFILE      *om;
  float            *scores;
  float            *fwd_emissions_arr;
} P7_PIPELINE_LONGTARGET_OBJS; 

static inline void   optimize_pack_paramvector        (double *p, struct optimize_data *data);
static inline void   optimize_unpack_paramvector      (double *p, struct optimize_data *data);

static inline double optimize_msvfilter_func          (double *p, int np, void *dptr);
static inline double optimize_viterbifilter_func      (double *p, int np, void *dptr);
static inline double optimize_forwardparser_func      (double *p, int np, void *dptr);
static inline double func_msvfilter    (ESL_RANDOMNESS *r, ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf,
					float time, int update_hmm, int calibrate);
static inline double func_viterbifilter(ESL_RANDOMNESS *r, ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf,
					float time, int update_hmm, int calibrate);
static inline double func_forwardparser(ESL_RANDOMNESS *r, ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf,
					float time, int update_hmm, int calibrate);
static inline void workaround_evolve_profile(ESL_RANDOMNESS *r, double time, int n, const P7_RATE *R, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om,
					     int calibrate);
static inline void workaround_recalibrate(ESL_RANDOMNESS *r, int len, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om);


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
 *            <eslETYPE> if <sq> is more than 100K long, which can
 *            happen when someone uses hmmsearch/hmmscan instead of
 *            nhmmer/nhmmscan on a genome DNA seq db.
 *
 * Xref:      J4/25.
 *
 * Note:      Error handling needs improvement. The <eslETYPE> exception
 *            was added as a late bugfix. It really should be an <eslEINVAL>
 *            normal error (because it's a user error). But then we need
 *            all our p7_Pipeline() calls to check their return status
 *            and handle normal errors appropriately, which we haven't 
 *            been careful enough about. [SRE H9/4]
 */
int
p7_EvoPipeline(P7_PIPELINE *pli, ESL_RANDOMNESS *r, P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm,
	       P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, const ESL_SQ *ntsq, P7_TOPHITS *hitlist,
	       float fixtime, int noevo, int recalibrate, int *ret_hmm_update)
{
  P7_HIT          *hit       = NULL;     /* ptr to the current hit output data      */
  float            time_star = 1.0;
  float            usc, vfsc, fwdsc;      /* filter scores                           */
  float            filtersc;              /* HMM null filter score                   */
  float            nullsc;                /* null model score                        */
  float            seqbias;  
  float            seq_score;             /* the corrected per-seq bit score */
  float            sum_score;             /* the corrected reconstruction score for the seq */
  float            pre_score, pre2_score; /* uncorrected bit scores for seq */
  double           P;                     /* P-value of a hit */
  double           lnP;                   /* log P-value of a hit */
  float            time;                  /* fwdfilter      time */
  float            spvtime;               /* sparse viterbi time */
  float            spftime;               /* sparse forward time */
  float            tol = 0.1;
  int              Ld;                    /* # of residues in envelopes */
  int              hmm_update = *ret_hmm_update; // do we need to update the hmm to a given time? 
  int              vit_opt    = FALSE;
  int              msv_opt    = FALSE;
  int              be_verbose = FALSE;
  int              d;
  int              vfsc_optimized;
  int              status;

  if (sq->n == 0) return eslOK;    /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */
  if (sq->n > 100000) ESL_EXCEPTION(eslETYPE, "Target sequence length > 100K, over comparison pipeline limit.\n(Did you mean to use nhmmer/nhmmscan?)");

  p7_omx_GrowTo(pli->oxf, om->M, 0, sq->n);    /* expand the one-row omx if needed */

  /* Base null model score (we could calculate this in NewSeq(), for a scan pipeline) */
  p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);

  /* First level filter: the MSV filter, multihit with <om> */
  // with noevo == TRUE
  time = time_star;
  if ((status = p7_OptimizeMSVFilter(r, sq->dsq, sq->n, &time, R, hmm, gm, om, bg, pli->oxf, &usc, fixtime, noevo, recalibrate,
				     msv_opt, hmm_update, tol)) != eslOK)      
    printf("\nsequence %s msvfilter did not optimize\n", sq->name);  
  //printf("^^MSV %s update? %d time %f usc %f\n", sq->name, hmm_update, time, usc);

  seq_score = (usc - nullsc) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  if (P > pli->F1) {
    *ret_hmm_update = (msv_opt)? TRUE : FALSE;
    return eslOK;
  }
  
  pli->n_past_msv++;

  /* biased composition HMM filtering */
  if (pli->do_biasfilter)
    {
      p7_bg_FilterScore(bg, sq->dsq, sq->n, &filtersc);
      seq_score = (usc - filtersc) / eslCONST_LOG2;
      P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      if (P > pli->F1) {
	*ret_hmm_update = (msv_opt)? TRUE : FALSE;
	return eslOK;
      }
    }
  else filtersc = nullsc;
  pli->n_past_bias++;

  /* In scan mode, if it passes the MSV filter, read the rest of the profile */
  if (pli->mode == p7_SCAN_MODELS)
    {
      if (pli->hfp) p7_oprofile_ReadRest(pli->hfp, om);
      p7_oprofile_ReconfigRestLength(om, sq->n);
      if ((status = p7_pli_NewModelThresholds(pli, om)) != eslOK) return status; /* pli->errbuf has err msg set */
    }

  /* Second level filter: ViterbiFilter(), multihit with <om> */
  vfsc_optimized = FALSE;
  if (P > pli->F2)
    {
      time       = time_star;
      hmm_update = (msv_opt)? TRUE:FALSE;
      if ((status = p7_OptimizeViterbiFilter(r, sq->dsq, sq->n, &time, R, hmm, gm, om, bg, pli->oxf, &vfsc, fixtime, noevo, recalibrate,
					     vit_opt, hmm_update, tol)) != eslOK) 
	printf("\nsequence %s vitfilter did not optimize\n", sq->name);
      if (vit_opt) {
	vfsc_optimized = TRUE;
	if (time == time_star) vfsc_optimized = FALSE;
      }

      //printf("^^VIT %s update? %d time %f vfsc %f\n", sq->name, hmm_update, time, vfsc);
      seq_score = (vfsc-filtersc) / eslCONST_LOG2;
      P  = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
      if (P > pli->F2) {
	*ret_hmm_update = (vit_opt)? TRUE : FALSE;
	return eslOK;
      }
    }
  pli->n_past_vit++;


  /* Parse it with Forward and obtain its real Forward score. */
  if (vfsc_optimized) {
    hmm_update = FALSE;
    fwdsc = func_forwardparser(NULL, sq->dsq, sq->n, hmm, R, gm, om, bg, pli->oxf, time, hmm_update, FALSE);
    //printf("^^FWD %s len %d time %f fwdsc %f\n", sq->name, sq->n, time, fwdsc);
  }
  else {
    time       = time_star;
    hmm_update = (vit_opt)? TRUE:FALSE;
    if ((status = p7_OptimizeForwardParser(r, sq->dsq, sq->n, &time, R, hmm, gm, om, bg, pli->oxf, &fwdsc, fixtime, noevo, recalibrate, hmm_update, tol)) != eslOK)      
      printf("\nsequence %s forwardparser did not optimize\n", sq->name);
    //printf("^^FWD OPT %s update? %d time %f fwdsc %f\n", sq->name, hmm_update, time, fwdsc);
  }

  seq_score = (fwdsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  if (P > pli->F3)  {
    *ret_hmm_update = TRUE;
    return eslOK;
  }
  pli->n_past_fwd++;
  

  /* ok, it's for real. Now a Backwards parser pass, and hand it to domain definition workflow */
  p7_omx_GrowTo(pli->oxb, om->M, 0, sq->n);
  p7_BackwardParser(sq->dsq, sq->n, om, pli->oxf, pli->oxb, NULL);

  status = p7_domaindef_ByPosteriorHeuristics(sq, ntsq, om, pli->oxf, pli->oxb, pli->fwd, pli->bck, pli->ddef, bg, FALSE, NULL, NULL, NULL);
  if (status != eslOK) ESL_FAIL(status, pli->errbuf, "domain definition workflow failure"); /* eslERANGE can happen  */
  if (pli->ddef->nregions   == 0) { *ret_hmm_update = TRUE; return eslOK; } /* score passed threshold but there's no discrete domains here       */
  if (pli->ddef->nenvelopes == 0) { *ret_hmm_update = TRUE; return eslOK; } /* rarer: region was found, stochastic clustered, no envelopes found */
  if (pli->ddef->ndom       == 0) { *ret_hmm_update = TRUE; return eslOK; } /* even rarer: envelope found, no domain identified {iss131}         */


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
	      sum_score += pli->ddef->dcl[d].envsc;         /* NATS */
	      Ld        += pli->ddef->dcl[d].jenv  - pli->ddef->dcl[d].ienv + 1;
	      seqbias   += pli->ddef->dcl[d].domcorrection; /* NATS */  
	    }
	}
      seqbias = p7_FLogsum(0.0, log(bg->omega) + seqbias);  /* NATS */
    }
  else 
    {
      for (d = 0; d < pli->ddef->ndom; d++) 
	{
	  if (pli->ddef->dcl[d].envsc > 0.0)
	    {
	      sum_score += pli->ddef->dcl[d].envsc;      /* NATS */
	      Ld        += pli->ddef->dcl[d].jenv  - pli->ddef->dcl[d].ienv + 1;
	    }
	}
      seqbias = 0.0;
    }    
  sum_score += (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); /* NATS */
  pre2_score = (sum_score - nullsc) / eslCONST_LOG2;                /* BITS */
  sum_score  = (sum_score - (nullsc + seqbias)) / eslCONST_LOG2;    /* BITS */

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
  lnP =  esl_exp_logsurv (seq_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
  if (p7_pli_TargetReportable(pli, seq_score, lnP))
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

      hit->pre_score  = pre_score; /* BITS */
      hit->pre_lnP    = esl_exp_logsurv (hit->pre_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

      hit->score      = seq_score; /* BITS */
      hit->lnP        = lnP;
      hit->sortkey    = pli->inc_by_E ? -lnP : seq_score; /* per-seq output sorts on bit score if inclusion is by score  */

      hit->sum_score  = sum_score; /* BITS */
      hit->sum_lnP    = esl_exp_logsurv (hit->sum_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

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
        hit->dcl[d].bitscore = hit->dcl[d].envsc + (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); /* NATS, for the moment... */
        hit->dcl[d].dombias  = (pli->do_null2 ? p7_FLogsum(0.0, log(bg->omega) + hit->dcl[d].domcorrection) : 0.0); /* NATS, and will stay so */
        hit->dcl[d].bitscore = (hit->dcl[d].bitscore - (nullsc + hit->dcl[d].dombias)) / eslCONST_LOG2; /* now BITS, as it should be */
        hit->dcl[d].lnP      = esl_exp_logsurv (hit->dcl[d].bitscore,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

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
        if (p7_pli_TargetReportable(pli, hit->score, hit->lnP))
        {
          hit->flags |= p7_IS_REPORTED;
          if (p7_pli_TargetIncludable(pli, hit->score, hit->lnP))
            hit->flags |= p7_IS_INCLUDED;
        }

        for (d = 0; d < hit->ndom; d++)
        {
          if (p7_pli_DomainReportable(pli, hit->dcl[d].bitscore, hit->dcl[d].lnP))
          {
            hit->dcl[d].is_reported = TRUE;
            if (p7_pli_DomainIncludable(pli, hit->dcl[d].bitscore, hit->dcl[d].lnP))
              hit->dcl[d].is_included = TRUE;
          }
        }
      }
	  
    }

  *ret_hmm_update = TRUE;
  return eslOK;
}



/* Function:  p7_pli_computeAliScores()
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
static inline int
p7_pli_computeAliScores(P7_DOMAIN *dom, ESL_DSQ *seq, const P7_SCOREDATA *data, int K)
{
  int status;
  int i, j, k;
  float sc;

  //Compute score contribution of each position in the alignment to the overall Viterbi score
  ESL_ALLOC( dom->scores_per_pos, sizeof(float) * dom->ad->N );
  for (i=0; i<dom->ad->N; i++)  dom->scores_per_pos[i] = 0.0;

  i = dom->iali - 1;        //sequence position
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







/*------------------- end, pipeline API -------------------------*/




int
p7_OptimizeMSVFilter(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
		     P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_usc,
		     float fixtime, int noevo, int recalibrate, int msv_opt, int hmm_update, float tol)
{
  struct optimize_data   data;
  ESL_MIN_CFG           *cfg   = esl_min_cfg_Create(1);
  ESL_MIN_DAT           *stats = esl_min_dat_Create(cfg);
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches       */
  double                 sc;
  double                 usc_init;
  float                  time_init;
  int                    isfixtime = (fixtime >= 0.0)? TRUE : FALSE;
  int                    np;
  int                    status;

  time_init = (isfixtime)? fixtime : *ret_time;
  usc_init = func_msvfilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, hmm_update, recalibrate);
  if (noevo || !msv_opt || isfixtime || usc_init == eslINFINITY) {
    *ret_usc  = usc_init;
    *ret_time = time_init;
    esl_min_cfg_Destroy(cfg);
    esl_min_dat_Destroy(stats);
    return eslOK;
  }
  
  // adjust tolerances
#if 1
  cfg->cg_rtol    = tol;
  cfg->cg_atol    = tol;
  cfg->brent_rtol = tol;
  cfg->brent_atol = tol;
  cfg->deriv_step = OPT_TIME_STEP;
#endif
  
  /* variables: time */
  np = 1;
  
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  
 /* Copy shared info into the "data" structure
   */
  data.time       = time_init;
  data.dsq        = (ESL_DSQ *)dsq;
  data.n          = n;
  data.R          = (P7_RATE *)R;
  data.hmm        = (P7_HMM *)hmm;
  data.gm         = (P7_PROFILE *)gm;
  data.om         = (P7_OPROFILE *)om;

  data.bg         = bg;
  data.oxf        = oxf;
  data.tol        = tol;
  data.errbuf     = NULL;
  data.be_verbose = FALSE;
 
  /* Create the parameter vector.
   */
  optimize_pack_paramvector(p, &data);
 
  /* pass problem to the optimizer
   */
  status = esl_min_ConjugateGradientDescent(cfg, p, np,
					    &optimize_msvfilter_func, NULL,
					    (void *) (&data), &sc, stats);
  if (0&&status == eslENOHALT) 
    printf("optimize_msvfilter(): bracket minimization did not converge. You may want to consider increasing the number of iterations\n");		
  if (status != eslENOHALT && status != eslOK) 
    esl_fatal("optimize_msvfilter(): bad bracket minimization");	
  
  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, &data);
  data.usc = -sc;
  //printf("END MSV OPTIMIZATION: time %f usc %f --> %f\n", data.time, usc_init, data.usc);

  if (usc_init > data.usc || data.usc == eslINFINITY) {
    *ret_usc  = usc_init;
    *ret_time = time_init;
    func_msvfilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, TRUE, recalibrate);
  }
  else {
    *ret_usc  = data.usc;
    *ret_time = data.time;
    if (recalibrate) p7_Calibrate(hmm, NULL, &r, &bg, NULL, NULL);
  }

  /* clean up */
  esl_min_cfg_Destroy(cfg);
  esl_min_dat_Destroy(stats);
  if (u != NULL) free(u);
  if (p != NULL) free(p);
  return eslOK;

 ERROR:
  if (cfg)   esl_min_cfg_Destroy(cfg);
  if (stats) esl_min_dat_Destroy(stats);
  if (p != NULL) free(p);
  if (u != NULL) free(u);
  return status;
}

int
p7_OptimizeViterbiFilter(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
			 P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_vfsc,
			 float fixtime, int noevo, int recalibrate, int vit_opt, int hmm_update, float tol)
{
  struct optimize_data   data;
  ESL_MIN_CFG           *cfg   = esl_min_cfg_Create(1);
  ESL_MIN_DAT           *stats = esl_min_dat_Create(cfg);
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                 sc;
  double                 vfsc_init;
  float                  time_init;
  int                    isfixtime = (fixtime >= 0.0)? TRUE : FALSE;
  int                    np;
  int                    status;

  time_init = (isfixtime)? fixtime : *ret_time;
  vfsc_init = func_viterbifilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, hmm_update, recalibrate);
  if (noevo || !vit_opt || isfixtime || vfsc_init == eslINFINITY) {
    *ret_vfsc = vfsc_init;
    *ret_time = time_init;
    esl_min_cfg_Destroy(cfg);
    esl_min_dat_Destroy(stats);
    return eslOK;
  }

  // adjust tolerances
#if 1
  cfg->cg_rtol    = tol;
  cfg->cg_atol    = tol;
  cfg->brent_rtol = tol;
  cfg->brent_atol = tol;
  cfg->deriv_step = OPT_TIME_STEP;
#endif
  
  /* variables: time */
  np = 1; 

  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));
  
 /* Copy shared info into the "data" structure
   */
  data.time       = time_init;
  data.dsq        = (ESL_DSQ *)dsq;
  data.n          = n;
  data.R          = (P7_RATE *)R;
  data.hmm        = (P7_HMM *)hmm;
  data.gm         = (P7_PROFILE *)gm;
  data.om         = (P7_OPROFILE *)om;
  data.bg         = bg;
  data.oxf        = oxf;
  data.tol        = tol;
  data.errbuf     = NULL;
  data.be_verbose = FALSE;
 
  /* Create the parameter vector.
   */
  optimize_pack_paramvector(p, &data);
 
  /* pass problem to the optimizer
   */
  status = esl_min_ConjugateGradientDescent(cfg, p, np,
					    &optimize_viterbifilter_func, NULL,
					    (void *) (&data), &sc, stats);
  if (0&&status == eslENOHALT) 
    printf("optimize_viterbiparser(): bracket minimization did not converge. You may want to consider increasing the number of iterations\n");		
  if (status != eslENOHALT && status != eslOK) 
    esl_fatal("optimize_viterbifilter(): bad bracket minimization");	
  
  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, &data);
  data.vfsc = -sc;
  printf("END VIT OPTIMIZATION: time %f vfsc %f --> %f\n", data.time, vfsc_init, data.vfsc);

  if (vfsc_init > data.vfsc || data.vfsc == eslINFINITY) {
    *ret_vfsc = vfsc_init;
    *ret_time = time_init;
    func_viterbifilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, TRUE, recalibrate);
  }
  else {
    *ret_vfsc = data.vfsc;
    *ret_time = data.time;
    if (recalibrate) p7_Calibrate(hmm, NULL, &r, &bg, NULL, NULL);
  }
  
  /* clean up */
  esl_min_cfg_Destroy(cfg);
  esl_min_dat_Destroy(stats);
  if (u != NULL) free(u);
  if (p != NULL) free(p);
  return eslOK;

 ERROR:
  if (cfg)   esl_min_cfg_Destroy(cfg);
  if (stats) esl_min_dat_Destroy(stats);
  if (p   != NULL) free(p);
  if (u   != NULL) free(u);
  return status;
}

int
p7_OptimizeForwardParser(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
			 P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_fwdsc,
			 float fixtime, int noevo, int recalibrate, int hmm_update, float tol)
{
  struct optimize_data   data;
  ESL_MIN_CFG           *cfg   = esl_min_cfg_Create(1);
  ESL_MIN_DAT           *stats = esl_min_dat_Create(cfg);
  double                *p = NULL;	       /* parameter vector                      */
  double                *u = NULL;             /* max initial step size vector          */
  double                 sc;
  double                 fwdsc_init;
  float                  time_init;
  int                    isfixtime = (fixtime >= 0.0)? TRUE : FALSE;
  int                    np;
  int                    status;

  time_init = (isfixtime)? fixtime : *ret_time;
  fwdsc_init = func_forwardparser(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, hmm_update, recalibrate);
  if (noevo || isfixtime || fwdsc_init == eslINFINITY) {
    *ret_fwdsc = fwdsc_init;
    *ret_time  = time_init;
    esl_min_cfg_Destroy(cfg);
    esl_min_dat_Destroy(stats);
    return eslOK;
  }
  
  // adjust tolerances
#if 1
  cfg->cg_rtol    = tol;
  cfg->cg_atol    = tol;
  cfg->brent_rtol = tol;
  cfg->brent_atol = tol;
  cfg->deriv_step = OPT_TIME_STEP;
#endif
  
  /* variables: time */
  np = 1;     
  
  /* allocate */
  ESL_ALLOC(p,   sizeof(double) * (np+1));
  ESL_ALLOC(u,   sizeof(double) * (np+1));

  /* Copy shared info into the "data" structure
   */
  data.time       = time_init;
  data.dsq        = (ESL_DSQ *)dsq;
  data.n          = n;
  data.R          = (P7_RATE *)R;
  data.hmm        = (P7_HMM *)hmm;
  data.gm         = (P7_PROFILE *)gm;
  data.om         = (P7_OPROFILE *)om;
  data.bg         = bg;
  data.oxf        = oxf;
  data.tol        = tol;
  data.errbuf     = NULL;
  data.be_verbose = FALSE;
 
  /* Create the parameter vector.
   */
  optimize_pack_paramvector(p, &data);

  /* pass problem to the optimizer
   */
  status = esl_min_ConjugateGradientDescent(cfg, p, np,
					    &optimize_forwardparser_func, NULL,
					    (void *) (&data), &sc, stats);
  
  if (0&&status == eslENOHALT) 
    printf("optimize_forwardparser(): bracket minimization did not converge. You may want to consider increasing the number of iterations\n");		
  if (status != eslENOHALT && status != eslOK) 
    esl_fatal("optimize_forwardparser(): bad bracket minimization. status %d tol %f", status, data.tol);		
  
  /* unpack the final parameter vector */
  optimize_unpack_paramvector(p, &data);
  data.fwdsc = -sc;
  //printf("END FWD OPTIMIZATION: time %f fwdsc %f --> %f\n", data.time, fwdsc_init, data.fwdsc);

  if (fwdsc_init > data.fwdsc || data.fwdsc == eslINFINITY) {
    *ret_fwdsc = fwdsc_init;
    *ret_time  = time_init;
    func_forwardparser(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, TRUE, recalibrate);
  }
  else {
    *ret_fwdsc = data.fwdsc;
    *ret_time  = data.time;
    if (recalibrate) p7_Calibrate(hmm, NULL, &r, &bg, NULL, NULL);
  }

  /* clean up */
  esl_min_cfg_Destroy(cfg);
  esl_min_dat_Destroy(stats);
  if (p != NULL) free(p);
  if (u != NULL) free(u);
  return eslOK;

 ERROR:
  if (cfg)   esl_min_cfg_Destroy(cfg);  
  if (stats) esl_min_dat_Destroy(stats);  
  if (p != NULL) free(p);
  if (u != NULL) free(u);
  return status;
}

static inline void
optimize_pack_paramvector(double *p, struct optimize_data *data)
{
  p[0] = (data->time < 1.0)? log(data->time) : data->time - 1.0;
}

static inline void
optimize_unpack_paramvector(double *p, struct optimize_data *data)
{
  data->time = (p[0] < 0.0)? exp(p[0]) : ((p[0] + 1.0 < TMAX)? p[0] + 1.0 : TMAX); 
}

static inline double
optimize_msvfilter_func(double *p, int np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  ESL_DSQ              *dsq = data->dsq;
  
  optimize_unpack_paramvector(p, data);
  
  data->usc = func_msvfilter(NULL, dsq, data->n, data->hmm, data->R, data->gm, data->om, data->bg, data->oxf, data->time, TRUE, FALSE);
  
  if (data->usc == eslINFINITY) data->usc = 1000.;
  return -(double)data->usc;
}

static inline double
optimize_viterbifilter_func(double *p, int np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  ESL_DSQ              *dsq = data->dsq;
  
  optimize_unpack_paramvector(p, data);
  
  data->vfsc = func_viterbifilter(NULL, dsq, data->n, data->hmm, data->R, data->gm, data->om, data->bg, data->oxf, data->time, TRUE, FALSE);
  
  if (data->vfsc == eslINFINITY) data->vfsc = 1000.;
  return -(double)data->vfsc;
}

static inline double
optimize_forwardparser_func(double *p, int np, void *dptr)
{
  struct optimize_data *data = (struct optimize_data *) dptr;
  ESL_DSQ              *dsq = data->dsq;
  
  optimize_unpack_paramvector(p, data);

  data->fwdsc = func_forwardparser(NULL, dsq, data->n, data->hmm, data->R, data->gm, data->om, data->bg, data->oxf, data->time, TRUE, FALSE);
  
  return -(double)data->fwdsc;
}



static inline double
func_msvfilter(ESL_RANDOMNESS *r, ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf,
	       float time, int hmm_update, int calibrate)
{
  float  usc;
  
  /* Construct the evolved profile */
  if (hmm_update) workaround_evolve_profile(r, (double)time, n, R, bg, hmm, gm, om, calibrate);
  
  p7_MSVFilter(dsq, n, om, oxf, &(usc));
  
#if 0
  printf("time %f usc %f\n", time, usc);
#endif

  return (double)usc;
 }

static inline double
func_viterbifilter(ESL_RANDOMNESS *r, ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf,
		   float time, int hmm_update, int calibrate)
{
  float  vfsc;
  
  /* Construct the evolved profile */
  if (hmm_update) workaround_evolve_profile(r, (double)time, n, R, bg, hmm, gm, om, calibrate);

  p7_ViterbiFilter(dsq, n, om, oxf, &(vfsc));
  
#if 0
  printf("time %f vfsc %f\n", time, vfsc);
#endif

  return (double)vfsc;
 }

static inline double
func_forwardparser(ESL_RANDOMNESS *r, ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf,
		   float time, int hmm_update, int calibrate)
{
  float   fwdsc;

  /* Construct the evolved profile */
  if (hmm_update) workaround_evolve_profile(r, (double)time, n, R, bg, hmm, gm, om, calibrate);
  
  p7_ForwardParser(dsq, n, om, oxf, &(fwdsc));

#if 0
  printf("time %f fwdsc %f\n", time, fwdsc);
#endif

  return (double)fwdsc;
 }

static inline void
workaround_evolve_profile(ESL_RANDOMNESS *r, double time, int len, const P7_RATE *R, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, int calibrate)
{  
  if (R == NULL) return;
  
  /* evolved HMM */
  p7_EvolveFromRate(NULL, hmm, R, bg, time, 0.0001, NULL, FALSE); 

  if (calibrate) p7_Calibrate(hmm, NULL, &r, &bg, NULL, NULL);
  
  /* evolved profiles gm and om */
  p7_ProfileConfig(hmm, bg, gm, len, p7_LOCAL);
  p7_oprofile_Convert(gm, om);
  
}

static inline void
workaround_recalibrate(ESL_RANDOMNESS *r, int len, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om)
{
  if (!r) return;
  
  /* calibrate HMM */
  p7_Calibrate(hmm, NULL, &r, &bg, NULL, NULL);

  /* evolved profiles gm and om */
  //p7_ProfileConfig(hmm, bg, gm, len, p7_LOCAL);
  //p7_oprofile_Convert(gm, om);
}


/*****************************************************************
 * @LICENSE@
 *
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/src/hmmer/trunk/src/search/p7_pipeline.c $
 * SVN $Id: p7_evopipeline.c 4617 2014-02-21 18:06:34Z eddys $
 *****************************************************************/
