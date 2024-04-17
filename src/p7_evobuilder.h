/* P7_EVOPIPELINE is the standardized pipeline for one profile/sequence
 * comparison, from the fast filters down through domain postprocessing,
 * alignment, and scoring.
 */

#ifndef P7_EVOBUILDER_INCLUDED
#define P7_EVOBUILDER_INCLUDED

#include <p7_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"

#include "evohmmer.h"

/* p7_evobuilder.c */
extern int p7_EvoBuilder      (P7_BUILDER *bld, ESL_MSA *msa, P7_BG *bg, HMMRATE *hmmrate, 
			       P7_HMM **opt_hmm, P7_TRACE ***opt_trarr, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om, ESL_MSA **opt_postmsa,
			       int noevo_msv, int noevo_vit, int noevo_fwd);
extern int p7_EvoSingleBuilder(P7_BUILDER *bld, ESL_SQ *sq,   P7_BG *bg, HMMRATE *hmmrate,
			       P7_HMM **opt_hmm, P7_TRACE  **opt_tr,    P7_PROFILE **opt_gm, P7_OPROFILE **opt_om,
			       int noevo_msv, int noevo_vit, int noevo_fwd); 


#endif /*P7_EVOBUILDER_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
