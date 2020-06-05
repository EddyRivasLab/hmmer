#ifndef h4VITFILTER_INCLUDED
#define h4VITFILTER_INCLUDED
#include "h4_config.h"

#include "easel.h"

#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_filtermx.h"

extern int (*h4_vitfilter)(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);

/* These are deliberately not protected by eslENABLE_SSE4, etc. We
 * always provide an implementation; but when an ISA isn't enabled,
 * the implementation issues a fatal error. See additonal comments in
 * vitfilter_sse.c, etc.
 */
extern int h4_vitfilter_sse   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);
extern int h4_vitfilter_avx   (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);
extern int h4_vitfilter_avx512(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc);

#endif /*p7VITFILTER_INCLUDED*/
