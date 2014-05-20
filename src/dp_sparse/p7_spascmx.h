#ifndef P7_SPASCMX_INCLUDED
#define P7_SPASCMX_INCLUDED

#include "p7_config.h"

#include "base/p7_coords2.h"
#include "dp_sparse/p7_sparsemx.h"

extern int    p7_spascmx_Reinit(P7_SPARSEMX *asf, const P7_SPARSEMASK *sm, const P7_COORD2 *anch, int D);
extern size_t p7_spascmx_MinSizeof(const P7_SPARSEMASK *sm, const P7_COORD2 *anch, int D, int64_t *opt_dalloc, int *opt_xalloc);
extern int    p7_spascmx_Dump(FILE *fp, const P7_SPARSEMX *asx, const P7_COORD2 *anch, int D);

#endif /*P7_SPASCMX_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
