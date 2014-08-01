/* There is no P7_SPASCMX structure. Rather, the p7_spascmx code
 * provides for using P7_SPARSEMX and P7_SPARSEMASK in anchor set
 * constrained DP calculations.
 */
#ifndef P7_SPASCMX_INCLUDED
#define P7_SPASCMX_INCLUDED

#include "p7_config.h"

#include "base/p7_anchors.h"
#include "dp_sparse/p7_sparsemx.h"


/* To create a sparse ASC matrix: see p7_sparsemx_Create()
 * To reuse it:                   see p7_sparsemx_Reuse()
 * To destroy it:                 see p7_sparsemx_Destroy()
 *
 * Only _Reinit() needs to be specialized, because that's where we
 * need to know the exact # of cells in the sparse ASC matrix.
 */

extern int    p7_spascmx_Resize(P7_SPARSEMX *asx, const P7_SPARSEMASK *sm, const P7_ANCHOR *anch, int D);
extern size_t p7_spascmx_MinSizeof(const P7_SPARSEMASK *sm, const P7_ANCHOR *anch, int D, int64_t *opt_dalloc, int *opt_xalloc);
extern int    p7_spascmx_Dump(FILE *fp, const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D);
extern int    p7_spascmx_CompareReference(const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D, const P7_REFMX *rxu, const P7_REFMX *rxd, float tol);

#endif /*P7_SPASCMX_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
