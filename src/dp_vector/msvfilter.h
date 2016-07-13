#ifndef p7MSVFILTER_INCLUDED
#define p7MSVFILTER_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_bg.h"
#include "base/p7_hmmwindow.h"
#include "base/p7_scoredata.h"
#include "dp_vector/p7_filtermx.h"
#include "dp_vector/p7_oprofile.h"

extern int p7_MSVFilter           (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_MSVFilter_sse          (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_MSVFilter_avx           (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_MSVFilter_avx512           (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);


extern int p7_SSVFilter_longtarget(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_FILTERMX *ox, const P7_SCOREDATA *msvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);

#endif  /*p7MSVFILTER_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
