/* Sparse dynamic programming, Viterbi/Forward/Backward/Decoding.
 */
#ifndef SPARSE_FWDBACK_INCLUDED
#define SPARSE_FWDBACK_INCLUDED

#include "easel.h"
#include "hmmer.h"
#include "p7_sparsemx.h"

/* Sparse DP routines {sparse_fwdback.c} */
extern int p7_SparseViterbi (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_SPARSEMX *sx, P7_TRACE *opt_tr, float *opt_sc);
extern int p7_SparseForward (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_SPARSEMX *sx,                   float *opt_sc);
extern int p7_SparseBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_SPARSEMX *sx,                   float *opt_sc);
extern int p7_SparseDecoding(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *sxf, P7_SPARSEMX *sxb, P7_SPARSEMX *sxd);

#endif
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
