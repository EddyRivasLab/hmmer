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


#endif /*P7_SPARSE_ASC_FWDBACK_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
