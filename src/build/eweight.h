#ifndef p7EWEIGHT_INCLUDED
#define p7EWEIGHT_INCLUDED

#include "p7_config.h"

#include "base/p7_hmm.h"
#include "base/p7_bg.h"
#include "base/p7_prior.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
extern int p7_EntropyWeight(const P7_HMM *hmm, const P7_BG *bg, const P7_PRIOR *pri, double infotarget, double *ret_Neff);
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7EWEIGHT_INCLUDED*/

