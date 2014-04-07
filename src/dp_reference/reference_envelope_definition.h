
#ifndef p7REFERENCE_ENVELOPE_DEFINITION_INCLUDED
#define p7REFERENCE_ENVELOPE_DEFINITION_INCLUDED
#include "p7_config.h"

#include "easel.h"

#include "base/p7_envelopes.h" 
#include "base/p7_coords2.h"
#include "base/p7_profile.h"

#include "dp_reference/p7_refmx.h"

extern int p7_reference_Envelopes(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_COORD2 *anch, int D,
				  const P7_REFMX *apu, const P7_REFMX *apd,
				  P7_REFMX *afu, P7_REFMX *afd, P7_ENVELOPE *env);


#endif /*p7REFERENCE_ENVELOPE_DEFINITION_INCLUDED*/
