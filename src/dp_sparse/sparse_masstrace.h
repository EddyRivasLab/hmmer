#include "easel.h"
#include "hmmer.h"
#include "p7_masstrace.h"
#include "p7_sparsemx.h"

extern int
p7_SparseMasstrace(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *fwd, const P7_SPARSEMX *bck, const P7_TRACE *tr, int z, float massthresh, 
		   P7_SPARSEMX *mass, P7_MASSTRACE *mt,
		   int *ret_iae, int *ret_ibe, int *ret_kae, int *ret_kbe);


