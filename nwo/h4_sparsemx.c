/* H4_SPARSEMX: DP matrix for production implementation
 *
 * Contents:
 *   1. H4_SPARSEMX object
 */
#include "h4_config.h"

#include "easel.h"

#include "h4_sparsemask.h"
#include "h4_sparsemx.h"

/*****************************************************************
 * 1. H4_SPARSEMX object
 *****************************************************************/

/* Function:  h4_sparsemx_Create()
 * Synopsis:  Create a new sparse DP matrix.
 *
 * Purpose:   Create a new sparse DP matrix, defined by the sparse DP
 *            mask <sm>.
 *            
 *            If <sm> is <NULL>, allocate for a reasonable-sized
 *            matrix. Matrix won't be usable until caller calls
 *            <h4_sparsemx_Reinit()>. This style works well in loops
 *            over sequences, where you can create the object once
 *            before any sequences are known, then
 *            <_Reinit()/_Reuse()> it in the loop.
 */
H4_SPARSEMX *
h4_sparsemx_Create(const H4_SPARSEMASK *sm)
{
  H4_SPARSEMX *sx            = NULL;
  int64_t      default_ncell = 4096; /* if <sm> is NULL, what <dalloc> should be (stored sparsified cells) */
  int          default_nx    = 256;  /* if <sm> is NULL, what <xalloc> should be (stored rows of specials) */
  int          status;

  ESL_ALLOC(sx, sizeof(H4_SPARSEMX));
  sx->dp   = NULL;
  sx->xmx  = NULL;
  sx->sm   = sm;
  sx->type = h4S_UNSET;

  /* We must avoid zero-sized mallocs. If there are no rows or cells, alloc the default sizes */
  sx->dalloc = ( (sm && sm->ncells)         ? sm->ncells       : default_ncell);
  sx->xalloc = ( (sm && (sm->nrow + sm->S)) ? sm->nrow + sm->S : default_nx);

  ESL_ALLOC(sx->dp,  sizeof(float) * h4S_NSCELLS * sx->dalloc);
  ESL_ALLOC(sx->xmx, sizeof(float) * h4S_NXCELLS * sx->xalloc);
  return sx;

 ERROR:
  h4_sparsemx_Destroy(sx);
  return NULL;
}

/* Function:  h4_sparsemx_Reinit()
 * Synopsis:  Reinitialize, reallocate sparse DP matrix for new calculation.
 *
 * Purpose:   Reinitialize an existing sparse matrix <sx> to use the
 *            sparse mask <sm>. Equivalent to a call to
 *            <h4_sparsemx_Create()> except we reuse as much memory as
 *            possible in the preexisting <sx>, to minimize
 *            realloc/free calls.
 *
 * Returns:   <eslOK> on success, and <sx> is now ready for sparse DP
 *            calculations defined by the mask <sm>.
 *
 * Throws:    <eslEMEM> on (re-)allocation failure.
 */
int
h4_sparsemx_Reinit(H4_SPARSEMX *sx, const H4_SPARSEMASK *sm)
{
  int64_t dalloc_req = sm->ncells;
  int     xalloc_req = sm->nrow + sm->S;
  int     status;

  if (sx->dalloc < dalloc_req) {
    ESL_REALLOC(sx->dp, sizeof(float) * h4S_NSCELLS * dalloc_req);
    sx->dalloc = dalloc_req;
  }
  if (sx->xalloc < xalloc_req) {
    ESL_REALLOC(sx->xmx, sizeof(float) * h4S_NXCELLS * xalloc_req);
    sx->xalloc = xalloc_req;
  }
  sx->sm   = sm;
  sx->type = h4S_UNSET;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  h4_sparsemx_Sizeof()
 * Synopsis:  Returns current allocated size of a H4_SPARSEMX, in bytes.
 *
 * Purpose:   Returns allocated size of <sx>, in bytes.  Does not
 *            include the size of its <H4_SPARSEMASK>.
 */
size_t
h4_sparsemx_Sizeof(const H4_SPARSEMX *sx)
{
  size_t n = sizeof(H4_SPARSEMX);
  n += sizeof(float) * h4S_NSCELLS * sx->dalloc;
  n += sizeof(float) * h4S_NXCELLS * sx->xalloc;
  return n;
}

/* Function:  h4_sparsemx_MinSizeof()
 * Synopsis:  Returns minimal required allocation size of a H4_SPARSEMX, in bytes.
 *
 * Purpose:   Calculate and return the minimum required size, in 
 *            bytes, of a sparse DP matrix, for a profile/sequence
 *            comparison using the sparse mask in <sm>.
 *            
 *            Taking a <H4_SPARSEMASK> as the arg, not <H4_SPARSEMX>,
 *            is not a typo. Does not require having an actual DP
 *            matrix allocated.  We use this function when
 *            planning/profiling memory allocation strategies.
 */
size_t
h4_sparsemx_MinSizeof(const H4_SPARSEMASK *sm)
{
  size_t n = sizeof(H4_SPARSEMX);
  n += sizeof(float) * h4S_NSCELLS * sm->ncells;          // dp[]
  n += sizeof(float) * h4S_NXCELLS * (sm->nrow + sm->S);  // xmx[]; for each seg ia..ib, ia-1..ib has special cells
  return n;
}

/* Function:  h4_sparsemx_Destroy()
 * Synopsis:  Free a sparse DP matrix.
 */
void
h4_sparsemx_Destroy(H4_SPARSEMX *sx)
{
  if (sx) {
    free(sx->dp);
    free(sx->xmx);
    /* sx->sm is a reference copy. caller remains responsible for it. */
    free(sx);
  }
}
/*----------- end, H4_SPARSEMX implementation -------------------*/

