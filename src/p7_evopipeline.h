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
  float           usc;
  float           vfsc;
  float           fwdsc;
  double          tol;
  char           *errbuf;
  int             be_verbose;
};

extern int p7_OptimizeMSVFilter    (ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
				    P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_usc,
				    float fixtime, int noevo, int recalibrate, int msv_opt, int hmm_evolve, float tol);
extern int p7_OptimizeViterbiFilter(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
				    P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_vfsc,
				    float fixtime, int noevo, int recalibrate, int vit_opt, int hmm_evolve, float tol);
extern int p7_OptimizeForwardParser(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int n, float *ret_time, P7_RATE *R,
				    P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg, P7_OMX *oxf, float *ret_fwdsc,
				    float fixtime, int noevo, int recalibrate,              int hmm_evolve, float tol);
extern int p7_EvoPipeline(P7_PIPELINE *pli, ESL_RANDOMNESS *r, float *evparam_star, P7_RATE *R, P7_HMM *hmm, P7_PROFILE *gm, P7_OPROFILE *om, P7_BG *bg,
			  const ESL_SQ *sq, const ESL_SQ *ntsq, P7_TOPHITS *hitlist, float fixtime, int noevo, int recalibrate, int *ret_hmm_restore);

#endif /*P7_EVOPIPELINE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
