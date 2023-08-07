#ifndef h4FWDFILTER_INCLUDED
#define h4FWDFILTER_INCLUDED
#include <h4_config.h>

#include "easel.h"

#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_checkptmx.h"
#include "h4_sparsemask.h"

extern int (*h4_fwdfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc);
extern int (*h4_bckfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh);

extern int h4_fwdfilter_sse   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc);
extern int h4_fwdfilter_avx   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc);
extern int h4_fwdfilter_avx512(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc);
extern int h4_fwdfilter_neon  (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc);
extern int h4_fwdfilter_vmx   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, float *opt_sc);

extern int h4_bckfilter_sse   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh);
extern int h4_bckfilter_avx   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh);
extern int h4_bckfilter_avx512(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh);
extern int h4_bckfilter_neon  (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh);
extern int h4_bckfilter_vmx   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_CHECKPTMX *cpx, H4_SPARSEMASK *sm, float sm_thresh);

#endif /*h4FWDFILTER_INCLUDED*/

