#ifndef p7SPARSE_ENVELOPES_INCLUDED
#define p7SPARSE_ENVELOPES_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_anchors.h"
#include "base/p7_envelopes.h"
#include "dp_sparse/p7_sparsemx.h"


extern int p7_sparse_Envelopes(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,
			       const P7_ANCHOR *anch, int D, 
			       P7_SPARSEMX *asf, const P7_SPARSEMX *asd,
			       P7_ENVELOPES *env);



#endif /*p7SPARSE_AEC_ALIGN_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
