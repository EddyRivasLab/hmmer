#ifndef p7MODELSTATS_INCLUDED
#define p7MODELSTATS_INCLUDED

#include <p7_config.h>

#include "base/p7_bg.h"
#include "base/p7_hmm.h"
#include "base/p7_hmmwindow.h"

extern double p7_MeanMatchInfo           (const P7_HMM *hmm, const P7_BG *bg);
extern double p7_MeanMatchEntropy        (const P7_HMM *hmm);
extern double p7_MeanMatchRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg);
extern int    p7_MeanPositionRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg, double *ret_entropy);
extern int    p7_hmm_CompositionKLDist(P7_HMM *hmm, P7_BG *bg, float *ret_KL, float **opt_avp);
extern int    p7_hmm_GetSimpleRepeats(P7_HMM *hmm, int maxK, int min_rep, int min_length, float relent_thresh, P7_HMM_WINDOWLIST *ranges);

#endif /*p7MODELSTATS_INCLUDED*/
