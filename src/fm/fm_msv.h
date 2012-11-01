#ifndef P7_FM_MSV_INCLUDED
#define P7_FM_MSV_INCLUDED

#include "p7_config.h"

#include "base/p7_bg.h"
#include "base/p7_hmmwindow.h"
#include "base/p7_scoredata.h"

#include "dp_vector/p7_oprofile.h"

extern int p7_FM_MSV( P7_OPROFILE *om, void *gx, float nu, P7_BG *bg, double F1,
         const FM_DATA *fmf, const FM_DATA *fmb, FM_CFG *fm_cfg, const P7_SCOREDATA *scoredata,
         P7_HMM_WINDOWLIST *windowlist);

#endif /*P7_FM_MSV_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
