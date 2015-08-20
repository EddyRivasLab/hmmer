#ifndef p7SPARSE_ANCHORS_INCLUDED
#define p7SPARSE_ANCHORS_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "base/p7_anchors.h"
#include "base/p7_anchorhash.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_sparse/p7_sparsemx.h"

#include "search/p7_mpas.h"

extern int p7_sparse_Anchors(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,
			     float vsc, float fsc, const P7_SPARSEMX *sxf, const P7_SPARSEMX *sxd, const P7_ANCHORS *vanch,
			     P7_TRACE *tr, float **byp_wrk, P7_ANCHORHASH *ah, 
			     P7_SPARSEMX *asf, P7_ANCHORS *anch, float *ret_asc,
			     P7_MPAS_PARAMS *prm);

extern int p7_sparse_AnchorsGlobal(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,
				   float vsc, float fsc, const P7_SPARSEMX *sxf, const P7_SPARSEMX *sxd, const P7_ANCHORS *vanch,
				   P7_TRACE *tr, float **byp_wrk, P7_ANCHORHASH *ah,
				   P7_SPARSEMX *asf, P7_ANCHORS *anch, float *ret_asc,
				   P7_MPAS_PARAMS *prm);


extern int p7_sparse_anchors_SetFromTrace(const P7_SPARSEMX *sxd, const P7_TRACE *tr, P7_ANCHORS *anch);


#endif /*p7SPARSE_ANCHORS_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
