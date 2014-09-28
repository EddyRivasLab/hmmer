#ifndef P7_SPARSE_AEC_ALIGN_INCLUDED
#define P7_SPARSE_AEC_ALIGN_INCLUDED

#include "p7_config.h"

#include "base/p7_profile.h"
#include "base/p7_envelopes.h"
#include "base/p7_trace.h"

#include "dp_sparse/p7_sparsemx.h"

extern int p7_sparse_aec_Align(const P7_PROFILE *gm, const P7_SPARSEMX *asd, 
			       P7_ENVELOPES *env, P7_SPARSEMX *aec, P7_TRACE *tr);




#endif /*P7_SPARSE_AEC_ALIGN_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
