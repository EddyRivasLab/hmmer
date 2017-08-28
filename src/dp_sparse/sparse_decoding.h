#ifndef p7SPARSE_DECODING_INCLUDED
#define p7SPARSE_DECODING_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "dp_sparse/p7_sparsemx.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif


extern int p7_SparseDecoding(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *sxf, P7_SPARSEMX *sxb, P7_SPARSEMX *sxd);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif

#endif /*p7SPARSE_DECODING_INCLUDED*/

