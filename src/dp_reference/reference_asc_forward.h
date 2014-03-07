#ifndef p7REFERENCE_ASC_FORWARD_INCLUDED
#define p7REFERENCE_ASC_FORWARD_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_coords2.h"

#include "dp_reference/p7_refmx.h"

extern int p7_ReferenceASCForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_COORD2 *anch, int D,
				  P7_REFMX *mxu, P7_REFMX *mxd, float *opt_sc);

#endif /*p7REFERENCE_ASC_FORWARD_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
