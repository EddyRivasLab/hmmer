/* Sparse dynamic programming, Forward/Backward.
 */
#ifndef SPARSE_FWDBACK_INCLUDED
#define SPARSE_FWDBACK_INCLUDED

#include "easel.h"
#include "hmmer.h"
#include "p7_sparsemx.h"

/* Sparse DP routines {sparse_fwdback.c} */
extern int p7_SparseForward (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_SPARSEMX *sxf, float *opt_sc);
extern int p7_SparseBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_SPARSEMX *sxb, float *opt_sc);

#endif
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
