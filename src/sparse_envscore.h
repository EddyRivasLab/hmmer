
#include "hmmer.h"
#include "p7_sparsemx.h"

extern int p7_SparseEnvScore(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
			     int iae, int ibe, int kae, int kbe,
			     P7_SPARSEMX *sx, float *ret_envsc);

