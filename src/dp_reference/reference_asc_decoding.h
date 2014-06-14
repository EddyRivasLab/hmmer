#ifndef p7REFERENCE_ASC_DECODING_INCLUDED
#define p7REFERENCE_ASC_DECODING_INCLUDED

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_anchors.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_asc_fwdback.h"

extern int p7_ReferenceASCDecoding(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D, 
				   const P7_REFMX *afu, const P7_REFMX *afd,
				         P7_REFMX *abu,       P7_REFMX *abd, 
				         P7_REFMX *apu,       P7_REFMX *apd);

#endif /*p7REFERENCE_ASC_DECODING_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
