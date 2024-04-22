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

#define OPT_TIME_STEP 0.5
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
					float time, int hmm_evolve, int calibrate);
static inline double func_viterbifilter(ESL_RANDOMNESS *r, ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf,
					float time, int hmm_evolve, int calibrate);
static inline double func_forwardparser(ESL_RANDOMNESS *r, ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf,
					float time, int hmm_evolve, int calibrate);
static inline void workaround_evolve_profile(ESL_RANDOMNESS *r, double time, int n, const P7_RATE *R, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om,
					     int calibrate);
static inline void workaround_calibrate_profile(ESL_RANDOMNESS *r, int len, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om);
static inline void workaround_restore_profile(float *evparam_star, int len, const P7_RATE *R, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om);



/* Function:  p7_Pipeline_Overthruster()
 * Synopsis:  The filters of HMMER3's accelerated seq/profile comparison pipeline.
 *
 * Purpose:   Run H3's comparison filters to determine if a sequence is sufficiently
 *            similar to an HMM to justify running the third stage.
 * 
 * Returns:   <eslOK> if the sequence and HMM are similar enough that the   
 *            main stage should be run.
 *            <eslFAIL> if the sequence and HMM are not similar enough that
 *            the main stage should be run.
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
 *            and handle normal errors appropriately, which we haven'ts
 *            been careful enough about. [SRE H9/4]
 */
extern int p7_EvoPipeline_Overthruster(P7_PIPELINE *pli, ESL_RANDOMNESS *r, float *evparam_star, EVOPIPE_OPT evopipe_opt,
				       P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq,
				       int *ret_hmm_restore, float *ret_fwdsc, float *ret_nullsc)
{
  float            time_star = 1.0;
  float            usc, vfsc, fwdsc;         /* filter scores                           */
  float            filtersc = -eslINFINITY;  /* HMM null filter score                   */
  float            nullsc;                   /* null model score                        */
  float            seq_score;                /* the corrected per-seq bit score */
  double           P;                        /* P-value of a hit */
  float            time;                     /* fwdfilter      time */
  float            spvtime;                  /* sparse viterbi time */
  float            spftime;                  /* sparse forward time */
  float            tol = 0.1;
  int              hmm_restore = *ret_hmm_restore; // do we need to restore the HMM to tstar? 
  int              be_verbose  = FALSE;
  int              hmm_evolve;
  int              vfsc_optimized;
  int              status;

  if (sq->n == 0)
    return eslFAIL; /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */
  if (sq->n > 100000) ESL_EXCEPTION(eslETYPE, "Target sequence length > 100K, over comparison pipeline limit.\n(Did you mean to use nhmmer/nhmmscan?)");

  p7_omx_GrowTo(pli->oxf, om->M, 0, sq->n);    /* expand the one-row omx if needed */

  /* Base null model score (we could calculate this in NewSeq(), for a scan pipeline) */
  p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

  // the bias composition score
  if (pli->do_biasfilter) p7_bg_FilterScore(bg, sq->dsq, sq->n, &filtersc);
 
  if (hmm_restore) 
    workaround_restore_profile(evparam_star, sq->n, R, bg, hmm, gm, om);
    
  /* First level filter: the MSV filter, multihit with <om> */
  time = time_star;
  hmm_evolve = FALSE;
  if ((status = p7_OptimizeMSVFilter(r, evopipe_opt, sq->dsq, sq->n, &time, R, hmm, gm, om, bg, pli->oxf, &usc, nullsc, filtersc, pli->F1, hmm_evolve, tol)) != eslOK)      
    printf("\nsequence %s msvfilter did not optimize\n", sq->name);  
  //printf("^^MSV %s evolve? %d time %f usc %f\n", sq->name, hmm_evolve, time, usc);

  seq_score = (usc - nullsc) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  if (P > pli->F1) {
    *ret_hmm_restore = (evopipe_opt.MSV_topt != TIMEOPT_NONE && time != time_star)? TRUE : FALSE;
    return eslFAIL;
  }
  pli->n_past_msv++;

  /* biased composition HMM filtering */
  if (pli->do_biasfilter)
    {
      seq_score = (usc - filtersc) / eslCONST_LOG2;
      P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      //printf("^^BIAS %s P %f F1 %f usc %f filtersc %f score %f\n", sq->name, P, pli->F1, usc, filtersc, seq_score);
      if (P > pli->F1) {
	*ret_hmm_restore = (evopipe_opt.MSV_topt != TIMEOPT_NONE && time != time_star)? TRUE : FALSE;
	return eslFAIL;
      }
    }
  else filtersc = nullsc;
  pli->n_past_bias++;

  /* In scan mode, if it passes the MSV filter, read the rest of the profile */
  if (pli->mode == p7_SCAN_MODELS)
    {
      if (pli->hfp)
	p7_oprofile_ReadRest(pli->hfp, om);
      p7_oprofile_ReconfigRestLength(om, sq->n);
      if ((status = p7_pli_NewModelThresholds(pli, om)) != eslOK)
	return status; /* pli->errbuf has err msg set */
    }
  
  /* Second level filter: ViterbiFilter(), multihit with <om> */
  vfsc_optimized = FALSE;
  if (P > pli->F2)
    {
      hmm_evolve = (evopipe_opt.MSV_topt != TIMEOPT_NONE)? TRUE:FALSE;
      if ((status = p7_OptimizeViterbiFilter(r, evopipe_opt, sq->dsq, sq->n, &time, R, hmm, gm, om, bg, pli->oxf, &vfsc, filtersc, pli->F2, hmm_evolve, tol)) != eslOK) 
	printf("\nsequence %s vitfilter did not optimize\n", sq->name);
      if (evopipe_opt.VIT_topt != TIMEOPT_NONE) {
	vfsc_optimized = TRUE;
	if (time == time_star) vfsc_optimized = FALSE;
      }

      //printf("^^VIT %s update? %d time %f vfsc %f\n", sq->name, hmm_evolve, time, vfsc);
      seq_score = (vfsc-filtersc) / eslCONST_LOG2;
      P  = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
      if (P > pli->F2) {
	*ret_hmm_restore = (evopipe_opt.VIT_topt != TIMEOPT_NONE && time != time_star)? TRUE : FALSE;
	return eslFAIL;
      }
    }
  pli->n_past_vit++;


  /* Parse it with Forward and obtain its real Forward score. */
  if (vfsc_optimized) {
    hmm_evolve = FALSE;
    fwdsc = func_forwardparser(NULL, sq->dsq, sq->n, hmm, R, gm, om, bg, pli->oxf, time, hmm_evolve, FALSE);
    //printf("^^FWD %s len %d time %f fwdsc %f\n", sq->name, sq->n, time, fwdsc);
  }
  else {
    hmm_evolve = (evopipe_opt.VIT_topt != TIMEOPT_NONE)? TRUE:FALSE;
    if ((status = p7_OptimizeForwardParser(r, evopipe_opt, sq->dsq, sq->n, &time, R, hmm, gm, om, bg, pli->oxf, &fwdsc, filtersc, pli->F3, hmm_evolve, tol)) != eslOK)      
      printf("\nsequence %s forwardparser did not optimize\n", sq->name);
    //printf("^^FWD OPT %s updated? %d time %f fwdsc %f filter %f score %f\n", sq->name, hmm_restore, time, fwdsc, filtersc, (fwdsc-filtersc) / eslCONST_LOG2);
  }

  seq_score = (fwdsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  //printf("^^FWD %s P %f time %f fwdsc %f filter %f score %f tau %f lambda %f\n", sq->name, P, time, fwdsc, filtersc, (fwdsc-filtersc) / eslCONST_LOG2, om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);

  if (P > pli->F3)  {
    *ret_hmm_restore = (time != time_star)? TRUE:FALSE;
    return eslFAIL;
  }
  pli->n_past_fwd++;
  
  *ret_hmm_restore = TRUE;
  *ret_fwdsc = fwdsc;
  *ret_nullsc = nullsc;
  return eslOK; // if we get this far, we passed all the filters and should proceed to the main stage
}

/* Function:  p7_Pipeline()
 * Synopsis:  HMMER3's accelerated seq/profile comparison pipeline.
 *
 * Purpose:   This function is now just a wrapper around the 
 *            p7_pipeline_Overthruster() and 
 *            p7_pipeline_Mainstage() calls: it first calls p7_pipeline_Overthruster(),
 *            to run HMMER's filter stages on the sequence/HMM comparison. If the 
 *            filters score highly enough that the main stage should be run, it
 *            calls p7_pipeline_Mainstage for final hit/miss determination and hitlist
 *            insertion.
 */
extern int p7_EvoPipeline(P7_PIPELINE *pli, ESL_RANDOMNESS *r, float *evparam_star, EVOPIPE_OPT evopipe_opt, P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm,
			  P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, const ESL_SQ *ntsq, P7_TOPHITS *hitlist, int *ret_hmm_restore)
{
  int status;
  float fwdsc;
  float nullsc;

  status = p7_EvoPipeline_Overthruster(pli, r, evparam_star, evopipe_opt, R, hmm, gm, om, bg, sq, ret_hmm_restore, &fwdsc, &nullsc);
  if (status == eslOK){ //run the main stage
    return (p7_Pipeline_Mainstage(pli, om, bg, sq, ntsq, hitlist, fwdsc, nullsc));
  }
  if (status == eslFAIL){  /* overthruster status of eslFAIL indicates that the thruster completed 
                              correctly, but the comparison did not score highly enough to proceed,
                              so return a correct completion of the pipeline */
    return eslOK;
  }
  else{
    return status;
  }
}

/*------------------- end, pipeline API -------------------------*/

extern int
p7_OptimizeMSVFilter(ESL_RANDOMNESS *r, EVOPIPE_OPT evopipe_opt, const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
		     P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_usc,
		     float nullsc, float filtersc, float F1, int hmm_evolve, float tol)
{
  struct optimize_data   data;
  ESL_MIN_CFG           *cfg   = esl_min_cfg_Create(1);
  ESL_MIN_DAT           *stats = esl_min_dat_Create(cfg);
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                *wrk = NULL;           /* 4 tmp vectors of length nbranches       */
  double                 sc;
  double                 usc_init;
  double                 usc;
  double                 seq_score;
  double                 P;
  float                  F1_brac = 2.0*F1;
  float                  F1_grad = 1.5*F1;
  float                  time_init;
  float                  time;
  enum timeopt_e         MSV_topt = evopipe_opt.MSV_topt;
  int                    isfixtime = (evopipe_opt.fixtime >= 0.0)? TRUE : FALSE;
  int                    np = 1;
  int                    status;

  time_init = (isfixtime)? evopipe_opt.fixtime : *ret_time;
  usc_init  = func_msvfilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, hmm_evolve, evopipe_opt.recalibrate);
  if (usc_init == eslINFINITY) MSV_topt = TIMEOPT_NONE;

  // this is a filter; if the score is already good enough, we don't need to optimize
  seq_score = (usc_init - ESL_MAX(nullsc,filtersc)) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
  if (P <= F1) MSV_topt = TIMEOPT_NONE;

  // adjust tolerances
#if 1
  cfg->cg_rtol    = tol;
  cfg->cg_atol    = tol;
  cfg->brent_rtol = tol;
  cfg->brent_atol = tol;
  cfg->deriv_step = OPT_TIME_STEP;
#endif

  switch(MSV_topt) {
  case TIMEOPT_NONE:
    *ret_usc  = usc_init;
    *ret_time = time_init;
    break;

  case TIMEOPT_BRAC:
    // pass the filter?
    if (P > F1_brac) {
      *ret_usc  = usc_init;
      *ret_time = time_init;
      break;
    }
    
    // check if the Rate need to be calculated
    p7_RateCalculate(hmm, bg, R, NULL, FALSE);

    // first time < tstar
    time = time_init - cfg->deriv_step;
    usc = func_msvfilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time, TRUE, evopipe_opt.recalibrate);
    if (usc > usc_init) {
      *ret_usc  = usc;
      *ret_time = time;
    }
    else {
      time = time_init + cfg->deriv_step;
      usc = func_msvfilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time, TRUE, evopipe_opt.recalibrate);
      if (usc > usc_init) {
	*ret_usc  = usc;
	*ret_time = time;
      }
      else {
	*ret_usc  = usc_init;
	*ret_time = time_init;
      }
    }
    
    // done with this brac optimization
    if (evopipe_opt.recalibrate) workaround_calibrate_profile(r, n, bg, hmm, gm, om);
    break;

  case TIMEOPT_GRAD:
    // pass the filter?
    if (P > F1_grad) {
      *ret_usc  = usc_init;
      *ret_time = time_init;
      break;
    }
    
    // check if the Rate need to be calculated
    p7_RateCalculate(hmm, bg, R, NULL, FALSE);
    
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
      func_msvfilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, TRUE, evopipe_opt.recalibrate);
    }
    else {
      *ret_usc  = data.usc;
      *ret_time = data.time;
      if (evopipe_opt.recalibrate) workaround_calibrate_profile(r, n, bg, hmm, gm, om);
    }

  case TIMEOPT_BOTH:
    break;
    
  default:
    esl_fatal("MSVFilter: could not find time optimization option");
    break;
  } // end of the switch over time optimization options
  
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


extern int
p7_OptimizeViterbiFilter(ESL_RANDOMNESS *r, EVOPIPE_OPT evopipe_opt, const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
			 P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_vfsc,
			 float filtersc, float F2, int hmm_evolve, float tol)
{
  struct optimize_data   data;
  ESL_MIN_CFG           *cfg   = esl_min_cfg_Create(1);
  ESL_MIN_DAT           *stats = esl_min_dat_Create(cfg);
  double                *p = NULL;	       /* parameter vector                        */
  double                *u = NULL;             /* max initial step size vector            */
  double                 sc;
  double                 vfsc_init;
  double                 vfsc;
  float                  time_init;
  float                  time;
  double                 seq_score;
  double                 P;
  float                  F2_brac = 5.0*F2;
  float                  F2_grad = 1.5*F2;
  enum timeopt_e         VIT_topt = evopipe_opt.VIT_topt;
  int                    isfixtime = (evopipe_opt.fixtime >= 0.0)? TRUE : FALSE;
  int                    np = 1;
  int                    status;

  time_init = (isfixtime)? evopipe_opt.fixtime : *ret_time;
  vfsc_init = func_viterbifilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, hmm_evolve, FALSE);
  if (vfsc_init == eslINFINITY) VIT_topt = TIMEOPT_NONE;

  // this is a filter; if the score is already good enough, we don't need to optimize
  seq_score = (vfsc_init - filtersc) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
  if (P < F2) VIT_topt = TIMEOPT_NONE;

  // adjust tolerances
#if 1
  cfg->cg_rtol    = tol;
  cfg->cg_atol    = tol;
  cfg->brent_rtol = tol;
  cfg->brent_atol = tol;
  cfg->deriv_step = OPT_TIME_STEP;
#endif

  switch(VIT_topt) {
  case TIMEOPT_NONE:
    *ret_vfsc = vfsc_init;
    *ret_time = time_init;
    break;

  case TIMEOPT_BRAC:
    if (P > F2_brac) {
      *ret_vfsc = vfsc_init;
      *ret_time = time_init;
      break;
     }

    // check if the Rate need to be calculated
    p7_RateCalculate(hmm, bg, R, NULL, FALSE);

    time = time_init - cfg->deriv_step;
    vfsc = func_viterbifilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time, TRUE, evopipe_opt.recalibrate);

    if (vfsc > vfsc_init) {
      *ret_vfsc = vfsc;
      *ret_time = time;
    }
    else {
      time = time_init + cfg->deriv_step;
      vfsc = func_viterbifilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time, TRUE, evopipe_opt.recalibrate);
      if (vfsc > vfsc_init) {
	*ret_vfsc  = vfsc;
	*ret_time = time;
      }
      else {
	*ret_vfsc = vfsc_init;
	*ret_time = time_init;
      }
    }
    
    // done with this brac optimization
    if (evopipe_opt.recalibrate) workaround_calibrate_profile(r, n, bg, hmm, gm, om);
    break;
    
  case TIMEOPT_GRAD:
    // pass the filter?
    if (P > F2_grad) {
      *ret_vfsc = vfsc_init;
      *ret_time = time_init;
      break;
    }

    // check if the Rate need to be calculated
    p7_RateCalculate(hmm, bg, R, NULL, FALSE);
    
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
    //printf("END VIT OPTIMIZATION: time %f vfsc %f --> %f\n", data.time, vfsc_init, data.vfsc);
    
    if (vfsc_init > data.vfsc || data.vfsc == eslINFINITY) {
      *ret_vfsc = vfsc_init;
      *ret_time = time_init;
      func_viterbifilter(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, TRUE, evopipe_opt.recalibrate);
    }
    else {
      *ret_vfsc = data.vfsc;
      *ret_time = data.time;
      if (evopipe_opt.recalibrate) workaround_calibrate_profile(r, n, bg, hmm, gm, om);
    }
    break;
    
  case TIMEOPT_BOTH:
    break;

  default:
    esl_fatal("ViterbiFilter: could not find time optimization option");
    break;
    
  } // end of the switch over time optimization options

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
p7_OptimizeForwardParser(ESL_RANDOMNESS *r, EVOPIPE_OPT evopipe_opt, const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
			 P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_fwdsc,
			 float filtersc, float F3, int hmm_evolve, float tol)
{
  struct optimize_data   data;
  ESL_MIN_CFG           *cfg   = esl_min_cfg_Create(1);
  ESL_MIN_DAT           *stats = esl_min_dat_Create(cfg);
  double                *p = NULL;	       /* parameter vector                      */
  double                *u = NULL;             /* max initial step size vector          */
  double                 sc;
  double                 seq_score;
  double                 P;
  double                 fwdsc_init;
  double                 fwdsc;
  float                  time_init;
  float                  time;
  float                  F3_brac = 2.0*F3;
  float                  F3_grad = 1.0;
  enum timeopt_e         FWD_topt = evopipe_opt.FWD_topt;
  int                    isfixtime = (evopipe_opt.fixtime >= 0.0)? TRUE : FALSE;
  int                    np = 1;
  int                    status;

  time_init = (isfixtime)? evopipe_opt.fixtime : *ret_time;
  fwdsc_init = func_forwardparser(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, hmm_evolve, FALSE);
  if (fwdsc_init == eslINFINITY) FWD_topt = TIMEOPT_NONE;

  // this is NOT a filter; if the score is already good enough, we don't need to optimize
  seq_score = (fwdsc_init - filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);

  // adjust tolerances
#if 1
  cfg->cg_rtol    = tol;
  cfg->cg_atol    = tol;
  cfg->brent_rtol = tol;
  cfg->brent_atol = tol;
  cfg->deriv_step = OPT_TIME_STEP;
#endif

  switch(FWD_topt) {
  case TIMEOPT_NONE:
    *ret_fwdsc = fwdsc_init;
    *ret_time = time_init;
    break;

  case TIMEOPT_BRAC:
    if (P > F3_brac) {
      *ret_fwdsc = fwdsc_init;
      *ret_time  = time_init;
      break;
    }

    // check if the Rate need to be calculated
    p7_RateCalculate(hmm, bg, R, NULL, FALSE);
	
    time = time_init - cfg->deriv_step;
    fwdsc = func_forwardparser(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time, TRUE, evopipe_opt.recalibrate);
    if (fwdsc > fwdsc_init) {
      *ret_fwdsc = fwdsc;
      *ret_time  = time;
    }
    else {
      time = time_init + cfg->deriv_step;
      fwdsc = func_forwardparser(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time, TRUE, evopipe_opt.recalibrate);
      if (fwdsc > fwdsc_init) {
	*ret_fwdsc = fwdsc;
	*ret_time  = time;
      }
      else {
	*ret_fwdsc = fwdsc_init;
	*ret_time  = time_init;
      }
    }
    
    // done with this bracket optimization
    if (evopipe_opt.recalibrate) workaround_calibrate_profile(r, n, bg, hmm, gm, om);
    break;
    
  case TIMEOPT_GRAD:
    if (P > F3_grad) {
      *ret_fwdsc = fwdsc_init;
      *ret_time  = time_init;
      break;
    }

    // check if the Rate need to be calculated
    p7_RateCalculate(hmm, bg, R, NULL, FALSE);

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
      func_forwardparser(r, (ESL_DSQ *)dsq, n, hmm, R, gm, om, bg, oxf, time_init, TRUE, evopipe_opt.recalibrate);
    }
    else {
      *ret_fwdsc = data.fwdsc;
      *ret_time  = data.time;
      if (evopipe_opt.recalibrate) workaround_calibrate_profile(r, n, bg, hmm, gm, om);
    }
    break;
    
  default:
    esl_fatal("ForwardParser: could not find time optimization option");
    break;
    
  } // end of the switch over time optimization options
 
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
	       float time, int hmm_evolve, int calibrate)
{
  float  usc;
  
  /* Construct the evolved profile */
  if (hmm_evolve) workaround_evolve_profile(r, (double)time, n, R, bg, hmm, gm, om, calibrate);
  
  p7_MSVFilter(dsq, n, om, oxf, &(usc));
  
#if 0
  printf("time %f usc %f\n", time, usc);
#endif

  return (double)usc;
 }

static inline double
func_viterbifilter(ESL_RANDOMNESS *r, ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf,
		   float time, int hmm_evolve, int calibrate)
{
  float  vfsc;
  
  /* Construct the evolved profile */
  if (hmm_evolve) workaround_evolve_profile(r, (double)time, n, R, bg, hmm, gm, om, calibrate);

  p7_ViterbiFilter(dsq, n, om, oxf, &(vfsc));
  
#if 0
  printf("time %f vfsc %f\n", time, vfsc);
#endif

  return (double)vfsc;
 }

static inline double
func_forwardparser(ESL_RANDOMNESS *r, ESL_DSQ *dsq, int n, P7_HMM *hmm, P7_RATE *R, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf,
		   float time, int hmm_evolve, int calibrate)
{
  float   fwdsc;

  /* Construct the evolved profile */
  if (hmm_evolve) workaround_evolve_profile(r, (double)time, n, R, bg, hmm, gm, om, calibrate);
  
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
  p7_EvolveFromRate(NULL, hmm, R, bg, time, NULL, FALSE); 

  if (calibrate) p7_Calibrate(hmm, NULL, &r, NULL, NULL, NULL);
  
  /* evolved profiles gm and om */
  p7_ProfileConfig(hmm, bg, gm, len, p7_LOCAL);
  p7_oprofile_Convert(gm, om);
  
}

static inline void
workaround_calibrate_profile(ESL_RANDOMNESS *r, int len, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om)
{  
  p7_Calibrate(hmm, NULL, &r, NULL, NULL, NULL);
  
  // assign the calibration parameters to the profiles
  om->evparam[p7_MLAMBDA] = hmm->evparam[p7_MLAMBDA];
  om->evparam[p7_VLAMBDA] = hmm->evparam[p7_VLAMBDA];
  om->evparam[p7_FLAMBDA] = hmm->evparam[p7_FLAMBDA];
  om->evparam[p7_MMU]     = hmm->evparam[p7_MMU];
  om->evparam[p7_VMU]     = hmm->evparam[p7_VMU];
  om->evparam[p7_FTAU]    = hmm->evparam[p7_FTAU];

  if (gm != NULL) {
    gm->evparam[p7_MLAMBDA] = hmm->evparam[p7_MLAMBDA];
    gm->evparam[p7_VLAMBDA] = hmm->evparam[p7_VLAMBDA];
    gm->evparam[p7_FLAMBDA] = hmm->evparam[p7_FLAMBDA];
    gm->evparam[p7_MMU]     = hmm->evparam[p7_MMU];
    gm->evparam[p7_VMU]     = hmm->evparam[p7_VMU];
    gm->evparam[p7_FTAU]    = hmm->evparam[p7_FTAU];
  }
}

static inline void
workaround_restore_profile(float *evparam_star, int len, const P7_RATE *R, P7_BG *bg, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om)
{  
  if (R == NULL) return;
  
  /* evolved HMM */
  p7_EvolveFromRate(NULL, hmm, R, bg, 1.0, NULL, FALSE); 

  /* evolved profiles gm and om */
  p7_ProfileConfig(hmm, bg, gm, len, p7_LOCAL);
  p7_oprofile_Convert(gm, om);

  // reassign the calibration parameters from the HMMstar
  hmm->evparam[p7_MLAMBDA] = om->evparam[p7_MLAMBDA] = evparam_star[p7_MLAMBDA];
  hmm->evparam[p7_VLAMBDA] = om->evparam[p7_VLAMBDA] = evparam_star[p7_VLAMBDA];
  hmm->evparam[p7_FLAMBDA] = om->evparam[p7_FLAMBDA] = evparam_star[p7_FLAMBDA];
  hmm->evparam[p7_MMU]     = om->evparam[p7_MMU]     = evparam_star[p7_MMU];
  hmm->evparam[p7_VMU]     = om->evparam[p7_VMU]     = evparam_star[p7_VMU];
  hmm->evparam[p7_FTAU]    = om->evparam[p7_FTAU]    = evparam_star[p7_FTAU];
  hmm->flags              |= p7H_STATS;

  if (gm != NULL) {
    gm->evparam[p7_MLAMBDA] = evparam_star[p7_MLAMBDA];
    gm->evparam[p7_VLAMBDA] = evparam_star[p7_MLAMBDA];
    gm->evparam[p7_FLAMBDA] = evparam_star[p7_MLAMBDA];
    gm->evparam[p7_MMU]     = evparam_star[p7_MMU];
    gm->evparam[p7_VMU]     = evparam_star[p7_VMU];
    gm->evparam[p7_FTAU]    = evparam_star[p7_FTAU];
  }
}



/*****************************************************************
 * @LICENSE@
 *
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/src/hmmer/trunk/src/search/p7_pipeline.c $
 * SVN $Id: p7_evopipeline.c 4617 2014-02-21 18:06:34Z eddys $
 *****************************************************************/
