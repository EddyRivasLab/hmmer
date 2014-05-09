#ifndef p7REFERENCE_AEC_TRACE_INCLUDED
#define p7REFERENCE_AEC_TRACE_INCLUDED
#include "p7_config.h"

#include "base/p7_envelopes.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"

extern int p7_reference_aec_trace_MEG(const P7_PROFILE *gm, P7_ENVELOPES *env, const P7_REFMX *apd, const P7_REFMX *mx, P7_TRACE *tr);

#endif  /*p7REFERENCE_AEC_TRACE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
