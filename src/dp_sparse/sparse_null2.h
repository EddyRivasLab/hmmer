#include "p7_config.h"

#include "hmmer.h"
#include "p7_sparsemx.h"

extern int p7_sparse_Null2ByExpectation(const P7_PROFILE *gm, const P7_SPARSEMX *sxd, 
					int iae, int ibe, int kae, int kbe,
					float *wrk, float *null2);
