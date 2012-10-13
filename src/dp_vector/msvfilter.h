#ifndef P7_MSVFILTER_INCLUDED
#define P7_MSVFILTER_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_bg.h"
#include "dp_vector/p7_filtermx.h"
#include "dp_vector/p7_oprofile.h"
#include "p7_hmmwindow.h"
#include "p7_scoredata.h"

extern int p7_MSVFilter           (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_MSVFilter_longtarget(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_FILTERMX *ox, const P7_SCOREDATA *msvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);

#endif  /*P7_MSVFILTER_INCLUDED*/
