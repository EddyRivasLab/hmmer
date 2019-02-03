#ifndef h4REFERENCE_TRACE_INCLUDED
#define h4REFERENCE_TRACE_INCLUDED

#include "h4_config.h"

#include "easel.h"

#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_path.h"
#include "h4_refmx.h"

extern int h4_reference_trace_Viterbi(const H4_PROFILE *hmm, const H4_MODE *mo, const H4_REFMX *rmx, H4_PATH *pi);

#endif /*h4REFERENCE_TRACE_INCLUDED*/
