#include "easel.h"
#include "p7_sparsemx.h"

extern int p7_sparse_masstrace_Down(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *fwd, P7_SPARSEMX *mass, P7_TRACE *tr, int z, float massthresh, int *ret_ibe, int *ret_kbe);
extern int p7_sparse_masstrace_Up  (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *fwd, P7_SPARSEMX *mass, P7_TRACE *tr, int z, float massthresh, int *ret_iae, int *ret_kae);
