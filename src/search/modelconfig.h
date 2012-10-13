#ifndef P7_MODELCONFIG_INCLUDED
#define P7_MODELCONFIG_INCLUDED

#include "p7_config.h"

#include "base/p7_bg.h"
#include "base/p7_hmm.h"
#include "base/p7_profile.h"

extern int p7_profile_Config         (P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg);
extern int p7_profile_ConfigLocal    (P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L);
extern int p7_profile_ConfigUnilocal (P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L);
extern int p7_profile_ConfigGlocal   (P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L);
extern int p7_profile_ConfigUniglocal(P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L);
extern int p7_profile_ConfigCustom   (P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, int L, float nj, float pglocal);
extern int p7_profile_SetLength      (P7_PROFILE *gm, int L);

#endif /*P7_MODELCONFIG_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
