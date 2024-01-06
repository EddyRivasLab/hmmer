#ifndef p7REFERENCE_ANCHORS_INCLUDED
#define p7REFERENCE_ANCHORS_INCLUDED

#include <p7_config.h>

#include "easel.h"
#include "esl_random.h"

#include "base/p7_anchors.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"

#include "search/p7_mpas.h"

extern int p7_reference_Anchors(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
				const P7_REFMX *rxf, const P7_REFMX *rxd,
				P7_TRACE *tr,  float **byp_wrk,  P7_ANCHORHASH *ah,
				P7_REFMX *afu, P7_REFMX *afd, P7_ANCHORS *anch,  float *ret_asc,
				P7_MPAS_PARAMS *prm, P7_MPAS_STATS *stats);

extern int p7_reference_anchors_SetFromTrace(const P7_REFMX *pp, const P7_TRACE *tr, P7_ANCHORS *anch);

#endif /*p7REFERENCE_ANCHORS_INCLUDED*/

