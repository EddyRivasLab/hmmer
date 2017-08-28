#ifndef p7SPARSE_AEC_ALIGN_INCLUDED
#define p7SPARSE_AEC_ALIGN_INCLUDED

#include "p7_config.h"

#include "base/p7_profile.h"
#include "base/p7_envelopes.h"
#include "base/p7_trace.h"

#include "dp_sparse/p7_sparsemx.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
	
extern int p7_sparse_aec_Align(const P7_PROFILE *gm, const P7_SPARSEMX *asd, 
			       P7_ENVELOPES *env, P7_SPARSEMX *aec, P7_TRACE *tr);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7SPARSE_AEC_ALIGN_INCLUDED*/

