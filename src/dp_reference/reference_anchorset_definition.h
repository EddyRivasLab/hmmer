#ifndef p7REFERENCE_ANCHORSET_DEFINITION_INCLUDED
#define p7REFERENCE_ANCHORSET_DEFINITION_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "base/p7_coords2.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"

#include "search/p7_mpas.h"

extern int p7_reference_Anchors(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
				const P7_REFMX *rxf, const P7_REFMX *rxd,
				P7_TRACE *tr,  float **byp_wrk,  P7_COORDS2_HASH *hashtbl,
				P7_REFMX *afu, P7_REFMX *afd, P7_COORDS2 *anch,  float *ret_asc,
				P7_MPAS_PARAMS *prm, P7_MPAS_STATS *stats);

#endif /*p7REFERENCE_ANCHORSET_DEFINITION_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
