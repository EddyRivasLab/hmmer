#ifndef P7_VITFILTER_INCLUDED
#define P7_VITFILTER_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"

#include "p7_hmmwindow.h"

extern int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_longtarget(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox,
                                        float filtersc, double P, P7_HMM_WINDOWLIST *windowlist);

#endif /*P7_VITFILTER_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
