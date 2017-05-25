#ifndef p7SPARSE_MASSTRACE_INCLUDED
#define p7SPARSE_MASSTRACE_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_trace.h"
#include "base/p7_masstrace.h"

#include "dp_sparse/p7_sparsemx.h"

extern int p7_SparseMasstrace(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *fwd, const P7_SPARSEMX *bck, const P7_TRACE *tr, int z, float massthresh, 
			      P7_SPARSEMX *mass, P7_MASSTRACE *mt,
			      int *ret_iae, int *ret_ibe, int *ret_kae, int *ret_kbe);

#endif /*p7SPARSE_MASSTRACE_INCLUDED*/
