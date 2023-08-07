#ifndef p7SPARSE_ENVSCORE_INCLUDED
#define p7SPARSE_ENVSCORE_INCLUDED

#include <p7_config.h>

#include <stdio.h>

#include "easel.h"

#include "base/p7_profile.h"
#include "dp_sparse/p7_sparsemx.h"


extern int p7_SparseEnvscore(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, int iae, int ibe, int kae, int kbe,
			     const P7_SPARSEMASK *sm, P7_SPARSEMX *sx, float *ret_envsc);

extern int p7_SparseEnvscoreApprox(P7_PROFILE *gm, P7_SPARSEMX *sxf, int iae, int ibe, float *ret_envsc);

/* Debugging tools */
extern int p7_sparse_envscore_Dump(FILE *ofp, P7_SPARSEMX *sxe, int iae, int ibe, int kae, int kbe);
extern int p7_sparse_envscore_IntersectedMask(P7_SPARSEMASK *osm, int iae, int ibe, int kae, int kbe, P7_SPARSEMASK **ret_sm);

#endif /*p7SPARSE_ENVSCORE_INCLUDED*/


