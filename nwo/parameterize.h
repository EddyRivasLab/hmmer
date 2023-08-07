#ifndef h4PARAMETERIZE_INCLUDED
#define h4PARAMETERIZE_INCLUDED
#include <h4_config.h>

#include "h4_counts.h"
#include "h4_prior.h"
#include "h4_profile.h"

extern int h4_parameterize(const H4_COUNTS *ctm, const H4_PRIOR *pri, H4_PROFILE *hmm);

#endif /* h4PARAMETERIZE_INCLUDED */
