/* The h4_spascmx module does not provide an independent data
 * structure, but rather a set of routines for using H4_SPARSEMX and
 * H4_SPARSEMASK for sparse anchor set contrained DP calculations.
 * 
 * To be included in a sparse ASC matrix, a main supercell i,k must
 * satisfy all of the following three conditions:
 * 
 * 1. It's in the sparse mask. 
 * 2. It's in an ASC UP or DOWN sector.
 * 3. It is "connected to" its anchor.
 * 
 * In more detail:
 *    
 * 1. Row i is in a sparse mask segment g if seg[g].ia <= i <= seg[g].ib;
 *    Cell (i,k) is in the mask if for some z,  sm->k[i][z] == k
 *    
 * 2. Suppose d is the index of the next anchor on or after row i:
 *       d = argmin_d i0[d] >= i
 *    Then (i,k) is in UP(d) if k < k0[d]
 *    and (i,k) is in DOWN(d-1) if k >= k0[d-1]
 *    (It is possible to be in both UP and DOWN.)
 *    
 * 3. In principle, we only need to calculate (i,k,s) cell if there is
 *    a valid DP path that connects it to "its anchor" anch[d]: i.e.
 *    the anchor down and right from i,k for an UP cell, or up/left
 *    for a DOWN cell. However, calculating whether an individual
 *    i,k,s cell is connected to its anchor appears to require a DP
 *    calculation all of its own, so we don't do that.
 *    
 *    Instead, we use a less stringent criterion for "connected to its
 *    anchor", based only on whether rows i,i0 are connected, rather
 *    than a detailed path through cells. For cell (i,k),
 *    row i and its anchor i0 must be in the same segment g:
 *          seg[g].ia <= (i, i0) <= seg[g].ib
 *          
 * Another way to think about why we have the 3rd criterion for
 * connectedness, instead of just taking a pure intersection of the
 * sparse mask and the ASC UP/DOWN sectors: Consider a segment that
 * contains 0 anchors. Here we don't have to do anything (not even
 * storing specials), because no path through any cells in this
 * segment can pass thru an anchor. Consider a segment that contains 1
 * anchor. Now we'll just have an UP and a DOWN sector, with no rows
 * that have both UP and DOWN cells. The only time we can have a row
 * with both UP and DOWN cells in it is when there's 2 or more anchors
 * in the same sparse segment.
 * 
 * For a simple example of traversing a sparse ASC matrix, see
 * h4_spascmx_MinSizeof().
 * 
 * Contents:
 *    1. Using H4_SPARSEMX for ASC calculations.
 *    2. Debugging and development tools.
 *    3. Debugging and dev: validation routine.
 */
#include <h4_config.h>

#include "easel.h"

#include "h4_anchorset.h"
#include "h4_sparsemask.h"
#include "h4_spascmx.h"

/*****************************************************************
 * 1. Using H4_SPARSEMX for ASC calculations
 *****************************************************************/


/* Function:  h4_spascmx_Reinit()
 * Synopsis:  Reinitialize, reallocate sparse ASC matrix for new DP problem.
 *
 * Purpose:   Reinitialize and, if necessary, reallocate an existing sparse
 *            matrix <asx> for a new ASC DP calculation that will be
 *            dually constrained by anchor set <anch> and sparse mask <sm>.
 *            
 *            In an MPAS algorithm, where we're optimizing the anchor
 *            set, we don't know what it is yet (more precisely: we're
 *            going to try different anchor sets) so we want to do an
 *            allocation that will work for any anchor set. In this
 *            case, pass <NULL> for <anch>, and <_Reinit()> will
 *            allocate space sufficient for any anchor set.
 *            
 *            <asx> keeps an internal pointer to <sm>, so the caller
 *            must not modify <sm> while <asx> remains in use.
 *
 * Args:      asx   : sparse ASC matrix to reinitialize
 *            anch  : anchor set
 *            sm    : sparse mask
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error. Now <asx> is in an undefined state,
 *            and can only be <_Destroy>'ed.
 *            
 * Notes:     When <anch> is <NULL>, we allocate the upper bound on space:
 *            2x sparse matrices, where typically we end up using less
 *            than 1x.  If memory use becomes an issue, we can explore
 *            more space-efficient memory reallocation strategies. In
 *            the current design, we can reallocate at most once in
 *            MPAS per profile/seq comparison, a design that
 *            prioritizes minimum reallocation (i.e. speed) over memory
 *            use.
 */
int
h4_spascmx_Reinit(H4_SPARSEMX *asx, const H4_ANCHORSET *anch, const H4_SPARSEMASK *sm)
{
  int64_t dalloc_req;    // denominated in i,k supercells, each w/ h4_NSCELLS; number of i,k main supercells stored by a sparse ASC DP calculation
  int     xalloc_req;    // denominated in i supercells, each w/ h4_NXCELLS:   number of rows that have specials stored
  int     status;

  if (anch) 
    h4_spascmx_MinSizeof(sm, anch, &dalloc_req, &xalloc_req);
  else
    {
      dalloc_req = sm->ncells * 2;     // unknown anchor set: upper bound on space is 2x sparse main matrices (UP/DOWN)
      xalloc_req = (sm->S + sm->nrow); //  ... but still only need one set of specials.
    }

  /* <sm> could be completely empty! 
   * Avoid zero length mallocs by assigning arbitrary minimum allocation sizes. 
   */
  if (dalloc_req == 0) dalloc_req = 16;
  if (xalloc_req == 0) xalloc_req = 16;
      
  if (dalloc_req > asx->dalloc) {
    ESL_REALLOC(asx->dp, sizeof(float) * h4S_NSCELLS * dalloc_req);
    asx->dalloc = dalloc_req;
  }
  if (xalloc_req > asx->xalloc) {
    ESL_REALLOC(asx->xmx, sizeof(float) * h4S_NXCELLS * xalloc_req);
    asx->xalloc = xalloc_req;
  }
  asx->sm   = sm;
  return eslOK;
  
 ERROR: 
  return status;
}

/* Function:  h4_spascmx_MinSizeof()
 * Synopsis:  Calculate minimum allocation size needed for sparse ASC DP calculation.
 *
 * Purpose:   For a sparse ASC DP calculation dually constrained by
 *            sparse mask <sm> and anchorset <anch>, calculate minimum
 *            required allocation size. Return the total size in
 *            bytes.
 *            
 *            Optionally, return in <opt_dalloc> the number of (i,k)
 *            main supercells that need to be stored, and in
 *            <opt_xalloc> the number of i rows for which specials
 *            need to be stored. <h4_spascmx_Reinit()> uses these
 *            numbers when it reuses and reallocates an existing
 *            matrix for a new DP problem.
 *            
 *            This routine also makes a good example of how to
 *            traverse a sparse ASC matrix.
 *
 * Args:      sm         : sparse mask that constrains the DP calc
 *            anch       : anchor set 
 *            opt_dalloc : optRETURN: number of main i,k supercells that need to be stored
 *            opt_xalloc : optRETURN: number of rows for which specials need to be stored
 *
 * Returns:   Minimum allocation size required for the complete <H4_SPARSEMX>
 *            structure, in bytes.
 */
size_t
h4_spascmx_MinSizeof(const H4_ANCHORSET *anch, const H4_SPARSEMASK *sm, int64_t *opt_dalloc, int *opt_xalloc)
{
  size_t  n       = sizeof(H7_SPARSEMX);
  int     g       = 1;		// index of next or current segment. When we enter it, and while we're in it, in_seg = TRUE.
  int     in_seg  = FALSE;      //   ... this bumps to TRUE when we see ia(g), first row of segment; to FALSE on ib(g), last row.
  int     d       = 1;          // index of next anchor we will reach; thus a current cell may be in UP(d) or DOWN(d-1) sector.
  int     ndown   = 0;          //   ... this bumps to 1 when i reaches an anchor, then counts rows in the DOWN sector, then goes to 0 when we leave seg g. Using counter allows knowing when we're on top DOWN row.
  int64_t dalloc  = 0;          // number of supercells in main DP matrix
  int     xalloc  = 0;          // number of special supercells 
  int     i,z;

  for (i = 0; i <= sm->L; i++)
    {
      if      (i == anch->a[d].i0)  { ndown = 1;  d++;     }    // when i reaches next anchor; bump d to next domain index, and DOWN sector is active...
      else if (sm->n[i] == 0)       { ndown = 0;           }    //  ... until when we reach end of segment, when DOWN becomes inactive again.
      else if (ndown)               { ndown++;             }    // counting ndown lets us easily test if we're on the special top row. Not used in this routine; here for pedagogy, since this is a traversal example

      if      (i >  sm->seg[g].ib)  { in_seg = FALSE; g++; }    // g bumps, to start expecting to see start of segment <g> next.
      else if (i == sm->seg[g].ia)  { in_seg = TRUE;       }    // g might be S+1, but this is safe because of sentinel seg[S+1].ia=ib=L+2      

      if (ndown)                                                // if i is in a DOWN sector:
	{
	  for (z  = 0; z < sm->n[i]; z++)
	    if (sm->k[i][z] >= anch->a[d-1].k0) break;             
	  dalloc += (sm->n[i] - z);                             // z is now on first cell in DOWN row; remainder of line is valid
	}

      if ( (i >= sm->seg[g].ia-1 && anch->a[d].i0 <= sm->seg[g].ib) ||  // if we need to store specials for row i because of UP xB...
	   (ndown))                                                     //   .. or because of DOWN xE ...
	xalloc++;
      
      if (in_seg && anch[d].i0 <= sm->seg[g].ib)               // if i is in an UP sector:
	{
	  for (z = 0; z < sm->n[i]; z++)
	    if (sm->k[i][z] >= anch->a[d].k0) break;   	       
	  dalloc += z;                                        // z is now +1 past the last sparse cell on the UP row
	}
    }

  n += dalloc * sizeof(float) * h4S_NSCELLS;
  n += xalloc * sizeof(float) * h4S_NXCELLS;
  if (opt_dalloc) *opt_dalloc  = dalloc;
  if (opt_xalloc) *opt_xalloc  = xalloc;
  return n;
}
