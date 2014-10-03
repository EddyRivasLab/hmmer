#ifndef p7SPARSE_ASC_FWDBACK_INCLUDED
#define p7SPARSE_ASC_FWDBACK_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_anchors.h"
#include "dp_sparse/p7_sparsemx.h"


extern int p7_sparse_asc_Forward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D, 
				 const P7_SPARSEMASK *sm, P7_SPARSEMX *asf, float *opt_sc);

extern int p7_sparse_asc_ForwardSeg(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
				    const P7_ANCHOR     *anch, int D, int d,
				    const P7_SPARSEMASK *sm,   int ngap, int ia, int ib,
				    P7_SPARSEMX *asf, float xN, float xJ, float xC, float *dpc, float *xc,
				    float *opt_xN, float *opt_xJ, float *opt_xC, int *opt_d, float **opt_dpc, float **opt_xc, float *opt_asc);

extern int p7_sparse_asc_Backward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
				  const P7_ANCHOR *anch, int D, const P7_SPARSEMASK *sm,
				  P7_SPARSEMX *asb, float *opt_sc);

extern int p7_sparse_asc_Decoding(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
				  const P7_ANCHOR *anch, int D,
				  float totsc, const P7_SPARSEMX *asf, P7_SPARSEMX *asb, P7_SPARSEMX *asd);



#endif /*p7SPARSE_ASC_FWDBACK_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
