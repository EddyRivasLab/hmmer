#ifndef P7_EWEIGHT_INCLUDED
#define P7_EWEIGHT_INCLUDED

#include "p7_config.h"

#include "base/p7_hmm.h"
#include "base/p7_bg.h"
#include "base/p7_prior.h"

extern int p7_EntropyWeight(const P7_HMM *hmm, const P7_BG *bg, const P7_PRIOR *pri, double infotarget, double *ret_Neff);

#endif /*P7_EWEIGHT_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
