#ifndef h4EWEIGHT_INCLUDED
#define h4EWEIGHT_INCLUDED
#include "h4_config.h"

#include "h4_counts.h"
#include "h4_prior.h"

extern int h4_EntropyWeight(const H4_COUNTS *ctm, const H4_PRIOR *pri, int nseq, float etarget, float *ret_Neff);

#endif // h4EWEIGHT_INCLUDED
