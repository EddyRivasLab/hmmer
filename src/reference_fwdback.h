#ifndef p7REFERENCE_FWDBACK_INCLUDED
#define p7REFERENCE_FWDBACK_INCLUDED

#include "p7_config.h"

#include "hmmer.h"
#include "p7_refmx.h"

extern int p7_ReferenceForward (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, float *opt_sc);
extern int p7_ReferenceBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, float *opt_sc);
extern int p7_ReferenceAlign   (const P7_PROFILE *gm, float gamma, const P7_REFMX *pp,  P7_REFMX *rmx, P7_TRACE *tr, float *opt_gain);

#endif /* p7REFERENCE_FWDBACK_INCLUDED */
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
