#ifndef P7_SPARSE_NULL2_INCLUDED
#define P7_SPARSE_NULL2_INCLUDED

#include "p7_config.h"

#include "base/p7_profile.h"
#include "dp_sparse/p7_sparsemx.h"

extern int p7_sparse_Null2ByExpectation(const P7_PROFILE *gm, const P7_SPARSEMX *sxd, 
					int iae, int ibe, int kae, int kbe,
					float *wrk, float *null2);

#endif /*P7_SPARSE_NULL2_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
