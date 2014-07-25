#ifndef P7_SPARSE_ASC_FWDBACK_INCLUDED
#define P7_SPARSE_ASC_FWDBACK_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_anchors.h"
#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spascmx.h"


extern int p7_sparse_asc_Forward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D, 
				 const P7_SPARSEMASK *sm, P7_SPARSEMX *asf, float *opt_sc);

extern int p7_sparse_asc_ForwardSeg(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
				    const P7_ANCHOR     *anch, int D, int d,
				    const P7_SPARSEMASK *sm,   int g,
				    P7_SPARSEMX *asf, float xN, float xJ, float xC, float *dpc, float *xc,
				    float *ret_xN, float *ret_xJ, float *ret_xC, int *ret_d, float **ret_dpc, float **ret_xc, float *opt_asc);

#endif /*P7_SPARSE_ASC_FWDBACK_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
