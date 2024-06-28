/* P7_EVOPIPELINE is the standardized pipeline for one profile/sequence
 * comparison, from the fast filters down through domain postprocessing,
 * alignment, and scoring.
 */

#ifndef P7_EVOPIPELINE_INCLUDED
#define P7_EVOPIPELINE_INCLUDED

#include <p7_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"

#include "evohmmer.h"
#include "minimize.h"


enum timeopt_e {
  TIMEOPT_NONE = 0, // do not optimize the time 
  TIMEOPT_BRAC = 1, // do a simple bracketing by two points aroung tstar
  TIMEOPT_GRAD = 2, // gradient descent optimization
  TIMEOPT_BOTH = 3  // bracketing followed by gradient descent optimization
};

typedef struct {
  float          fixtime;
  int            noevo;
  int            recalibrate;
  enum timeopt_e MSV_topt;
  enum timeopt_e VIT_topt;
  enum timeopt_e FWD_topt;
} EVOPIPE_OPT;

struct optimize_data {
  float            time;
  int              unihit;
  ESL_DSQ         *dsq;
  int              n;
  P7_RATE         *R;
  P7_HMM          *hmm;
  P7_PROFILE      *gm;
  P7_OPROFILE     *om;
  P7_BG           *bg;
  P7_OMX          *oxf;
  float            usc;
  float            vfsc;
  float            fwdsc;
  double           tol;
  char            *errbuf;
  int              be_verbose;
};

extern int p7_OptimizeMSVFilter    (ESL_RANDOMNESS *r, ESL_MIN_CFG *cfg, ESL_MIN_DAT *stats, EVOPIPE_OPT evopipe_opt,
				    const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
				    P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_usc,
				    float nullsc, float filtersc, float F1, float tol);
extern int p7_OptimizeViterbiFilter(ESL_RANDOMNESS *r, ESL_MIN_CFG *cfg, ESL_MIN_DAT *stats, EVOPIPE_OPT evopipe_opt,
				    const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
				    P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_vfsc,
				    float filtersc, float F2, float tol);
extern int p7_OptimizeForwardParser(ESL_RANDOMNESS *r, ESL_MIN_CFG *cfg, ESL_MIN_DAT *stats, EVOPIPE_OPT evopipe_opt,
				    const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
				    P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_fwdsc,
				    float filtersc, float F3, float tol);
extern int p7_EvoPipeline_Overthruster(P7_PIPELINE *pli, ESL_RANDOMNESS *r, float *evparam_star, EVOPIPE_OPT evopipe_opt,
				       P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, const ESL_SQ *sq, const ESL_SQ *ntsq, 
				       int *ret_hmm_restore, float *ret_fwdsc, float *ret_nullsc, float *ret_time);
extern int p7_EvoPipeline_Mainstage(P7_PIPELINE * pli, float *evparam_star, P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om,
				    P7_BG * bg, const ESL_SQ *sq, const ESL_SQ *ntsq, P7_TOPHITS *hitlist,
				    float fwdsc, float nullsc, float time, int *ret_hmm_restore);
extern int p7_EvoPipeline(P7_PIPELINE *pli, ESL_RANDOMNESS *r, float *evparam_star, EVOPIPE_OPT evopipe_opt,
			  P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg,
			  const ESL_SQ *sq, const ESL_SQ *ntsq, P7_TOPHITS *hitlist, int *ret_hmm_restore);

#endif /*P7_EVOPIPELINE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
