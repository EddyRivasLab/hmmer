#ifndef p7SPARSE_NULL2_INCLUDED
#define p7SPARSE_NULL2_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "base/p7_envelopes.h"
#include "base/p7_profile.h"
#include "dp_sparse/p7_sparsemx.h"

extern int p7_sparse_Null2(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *apd, P7_ENVELOPES *env, float *wrk, float *null2);

/* legacy: */
extern int p7_sparse_Null2ByExpectation(const P7_PROFILE *gm, const P7_SPARSEMX *sxd, 
					int iae, int ibe, int kae, int kbe,
					float *wrk, float *null2);

#endif /*p7SPARSE_NULL2_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
