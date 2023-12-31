#ifndef h4SPASCMX_INCLUDED
#define h4SPASCMX_INCLUDED
#include <h4_config.h>

#include "h4_anchorset.h"
#include "h4_sparsemask.h"
#include "h4_sparsemx.h"

extern int    h4_spascmx_Reinit(H4_SPARSEMX *asx, const H4_ANCHORSET *anch, const H4_SPARSEMASK *sm);
extern size_t h4_spascmx_MinSizeof(const H4_ANCHORSET *anch, const H4_SPARSEMASK *sm, int64_t *opt_dalloc, int *opt_xalloc);


#endif //h4SPASCMX_INCLUDED
