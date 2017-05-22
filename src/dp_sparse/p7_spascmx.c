/* The p7_spascmx module does not provide an independent data
 * structure, but rather a set of routines for using P7_SPARSEMX and
 * P7_SPARSEMASK for sparse anchor set contrained DP calculations.
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
 * 1. Row i is in the sparse mask segment g if  seg[g].ia <= i <= seg[g].ib;
 *    Cell (i,k) is in the mask if for some z,  sm->k[i][z] == k
 *    
 * 2. Suppose d is the index of the anchor on or after i:
 *       d = argmin_d anch[d].i0 >= i
 *    Then (i,k) is in UP(d) if k < anch[d].k0;
 *    and (i,k) is in DOWN(d-1) if k >= anch[d-1].k0
 *    (It is possible to be in both UP and DOWN.)
 *    
 * 3. In principle, we only need to calculate (i,k,s) cell if
 *    there is a valid DP path that connects it to "its anchor":
 *    where "its anchor" means anch[d].i0,k0 for a cell in UP(d)
 *    (the anchor down and right from i,k), or anch[d-1].i0,k0 
 *    for a cell in DOWN(d-1) (anchor up/left from i,k). However,
 *    calculating whether an individual i,k,s cell is connected
 *    to its anchor appears to require a DP calculation all of its
 *    own, so we don't do that.
 *    
 *    Instead, we use a less stringent criterion for "connected to its
 *    anchor", based only on whether rows i,i0 are connected, rather
 *    than a detailed path through cells. For cell (i,k),
 *    row i and its anchor i0 must be in the same segment g:
 *          seg[g].ia <= (i, i0) <= seg[g].ib
 *          
 * Another way to think about why we have the 3rd criterion for
 * connectedness, instead of just taking the intersection of the
 * sparse mask and the ASC UP/DOWN sectors: Consider a segment that
 * contains 0 anchors. Here we don't have to do anything (not even
 * storing specials), because no path through any cells in this
 * segment can pass thru an anchor. Consider a segment that contains 1
 * anchor. Now we'll just have an UP and a DOWN sector, with no rows
 * that have both UP and DOWN cells. The only time we can have a row
 * with both UP and DOWN cells in it is when there's 2 or more anchors
 * in the same sparse segment.
 * 
 *          
 * For a simple example of traversing a sparse ASC matrix, see
 * p7_spascmx_MinSizeof().
 * 
 * 
 * Contents:
 *    1. Using P7_SPARSEMX for ASC calculations.
 *    2. Debugging and development tools.
 *    3. Debugging and dev: Validation routine.
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"

#include "base/p7_anchors.h"
#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spascmx.h"


/*****************************************************************
 * 1. Using P7_SPARSEMX for ASC calculations
 *****************************************************************/


/* Function:  p7_spascmx_Resize()
 * Synopsis:  Reinitialize, reallocate sparse ASC matrix for new DP problem.
 *
 * Purpose:   Reinitialize and, if necessary, reallocate an existing sparse
 *            matrix <asx> for a new ASC DP calculation that will be
 *            constrained both by the sparse mask <sm> and the 
 *            anchor array <anch>.
 *            
 *            In an MPAS algorithm, where we're optimizing the anchor
 *            set, we don't know what it is yet (more precisely: we're
 *            going to try different anchor sets) so we want to do an
 *            allocation that will work for any anchor set. In this
 *            case, pass <NULL> for <anch> (and anything for <D>; 0 is
 *            fine).  Then <_Resize()> will allocate space sufficient
 *            for any anchor set.
 *            
 *            <asx> keeps an internal pointer to <sm>, so the caller
 *            must not modify <sm> while <asx> remains in use.
 *
 * Args:      asx   : sparse ASC matrix to reinitialize
 *            sm    : sparse mask that will constrain the new ASC DP calculation
 *            anch  : anchor array (1..D, with sentinels) that constrains new ASC DP calc
 *            D     : number of domains defined by anchors in <anch>
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
p7_spascmx_Resize(P7_SPARSEMX *asx, const P7_SPARSEMASK *sm, const P7_ANCHOR *anch, int D)
{
  int64_t dalloc_req;    // denominated in i,k supercells, each w/ p7_NSCELLS; number of i,k main supercells stored by a sparse ASC DP calculation
  int     xalloc_req;    // denominated in i supercells, each w/ p7_NXCELLS:   number of rows that have specials stored
  int     status;

  if (anch) 
    p7_spascmx_MinSizeof(sm, anch, D, &dalloc_req, &xalloc_req);
  else
    {
      dalloc_req = sm->ncells * 2;   // unknown anchor set: upper bound on space is 2x sparse matrices
      xalloc_req = (sm->S + sm->nrow);
    }

  /* <sm> could be completely empty. 
   * Avoid zero length mallocs by assigning arbitrary minimum allocation sizes. 
   */
  if (dalloc_req == 0) dalloc_req = 16;
  if (xalloc_req == 0) xalloc_req = 16;
      
  if (dalloc_req > asx->dalloc) {
    ESL_REALLOC(asx->dp, sizeof(float) * p7S_NSCELLS * dalloc_req);
    asx->dalloc = dalloc_req;
  }
  if (xalloc_req > asx->xalloc) {
    ESL_REALLOC(asx->xmx, sizeof(float) * p7S_NXCELLS * xalloc_req);
    asx->xalloc = xalloc_req;
  }
  asx->sm   = sm;
  return eslOK;
  
 ERROR: 
  return status;
}


/* Function:  p7_spascmx_MinSizeof()
 * Synopsis:  Calculate minimum allocation size needed for sparse ASC DP calculation.
 *
 * Purpose:   For a sparse ASC DP calculation constrained by sparse mask
 *            <sm> and anchor array <anch> for <1..D> domains,
 *            calculate minimum required allocation size. Return the
 *            total size in bytes. You might use this for collecting
 *            statistics on how much memory is required by sparse ASC
 *            DP calculations.
 *            
 *            Optionally, return in <opt_dalloc> the number of (i,k)
 *            main supercells that need to be stored, and in
 *            <opt_xalloc> the number of i rows for which specials
 *            need to be stored. <p7_spascmx_Resize()> uses these
 *            numbers when it reallocates a matrix for a new DP problem.
 *            
 *            This routine also makes a good example of how to
 *            traverse a sparse ASC matrix.
 *
 * Args:      sm         : sparse mask that constrains the DP calc
 *            anch       : anchor set array 1..D (with sentinels 0,D+1) 
 *            D          : number of domains/anchors in <anch>
 *            opt_dalloc : optRETURN: number of main i,k supercells that need to be stored
 *            opt_xalloc : optRETURN: number of rows for which specials need to be stored
 *
 * Returns:   Minimum allocation size required for the complete <P7_SPARSEMX>
 *            structure, in bytes.
 *
 * Throws:    (no abnormal error conditions)
 */
size_t
p7_spascmx_MinSizeof(const P7_SPARSEMASK *sm, const P7_ANCHOR *anch, int D, int64_t *opt_dalloc, int *opt_xalloc)
{
  size_t  n       = sizeof(P7_SPARSEMX);
  int     g       = 1;		// index of next or current segment. When we enter it, and while we're in it, in_seg = TRUE.
  int     in_seg  = FALSE;      //   ... this bumps to TRUE when we see ia(g), first row of segment; to FALSE on ib(g), last row.
  int     d       = 1;          // index of next anchor we will reach; thus a current cell may be in UP(d) or DOWN(d-1) sector.
  int     ndown   = 0;          //   ... this bumps to 1 when i reaches an anchor, then counts rows in the DOWN sector, then goes to 0 when we leave seg g. Using counter allows knowing when we're on top DOWN row.
  int64_t dalloc  = 0;          // number of supercells in main DP matrix
  int     xalloc  = 0;          // number of special supercells 
  int     i,z;

  for (i = 0; i <= sm->L; i++)
    {
      if      (i == anch[d].i0)     { ndown = 1;  d++;     }    // when i reaches next anchor; bump d to next domain index, and DOWN sector is active...
      else if (sm->n[i] == 0)       { ndown = 0;           }    //  ... until when we reach end of segment, when DOWN becomes inactive again.
      else if (ndown)               { ndown++;             }    // counting ndown lets us easily test if we're on the special top row. Not used in this routine; here for pedagogy, since this is a traversal example

      if      (i >  sm->seg[g].ib)  { in_seg = FALSE; g++; }    // g bumps, to start expecting to see start of segment <g> next.
      else if (i == sm->seg[g].ia)  { in_seg = TRUE;       }    // g might be S+1, but this is safe because of sentinel seg[S+1].ia=ib=L+2      

      if (ndown)                                                // if i is in a DOWN sector:
	{
	  for (z  = 0; z < sm->n[i]; z++)
	    if (sm->k[i][z] >= anch[d-1].k0) break;             
	  dalloc += (sm->n[i] - z);                             // z is now on first cell in DOWN row; remainder of line is valid
	}

      if ( (i >= sm->seg[g].ia-1 && anch[d].i0 <= sm->seg[g].ib) ||  // if we need to store specials for row i because of UP xB...
	   (ndown))                                                  //   .. or because of DOWN xE ...
	xalloc++;
      
      if (in_seg && anch[d].i0 <= sm->seg[g].ib)               // if i is in an UP sector:
	{
	  for (z = 0; z < sm->n[i]; z++)
	    if (sm->k[i][z] >= anch[d].k0) break;   	       
	  dalloc += z;                                        // z is now +1 past the last sparse cell on the UP row
	}
    }

  n += dalloc * sizeof(float) * p7S_NSCELLS;
  n += xalloc * sizeof(float) * p7S_NXCELLS;
  if (opt_dalloc) *opt_dalloc  = dalloc;
  if (opt_xalloc) *opt_xalloc  = xalloc;
  return n;
}


int
p7_spascmx_MinSizeofSeg(const P7_SPARSEMASK *sm, const P7_ANCHOR *anch, int D, int d, int g, 
			int64_t *ret_dalloc, int *ret_xalloc)
{
  int     in_down = FALSE;
  int64_t dalloc  = 0;
  int     xalloc  = 0;
  int i,z;
  
  if (anch[d].i0 <= sm->seg[g].ib)  // If there's no anchor in this segment: no storage needed.
    {
      xalloc = sm->seg[g].ib - sm->seg[g].ia + 2;

      for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)
	{
	  if (i == anch[d].i0) { in_down = TRUE; d++; }
	  
	  if (in_down)
	    {
	      for (z = 0; z < sm->n[i]; z++)
		if (sm->k[i][z] >= anch[d-1].k0) break;
	      dalloc += (sm->n[i] - z);
	    }
      
	  if (anch[d].i0 <= sm->seg[g].ib)
	    {
	      for (z = 0; z < sm->n[i]; z++)
		if (sm->k[i][z] >= anch[d].k0) break;
	      dalloc += z;;
	    }
	}
    }

  *ret_dalloc = dalloc;
  *ret_xalloc = xalloc;
  return eslOK;
}


int
p7_spascmx_Renormalize(P7_SPARSEMX *asd, const P7_ANCHOR *anch, int D)
{
  const P7_SPARSEMASK *sm  = asd->sm;
  float *dpc = asd->dp;
  float *xc  = asd->xmx;
  float *dpc0;
  int d = 1;
  float  norm;
  int i,g,x,z;

  for (g = 1; g <= sm->S; g++)
    {
      if (anch[d].i0 > sm->seg[g].ib) continue;

      norm = xc[p7S_N] + xc[p7S_JJ] + xc[p7S_CC];
      for (x = 0; x < p7S_NXCELLS; x++)
	xc[x] /= norm;
      xc += p7S_NXCELLS;

      for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)
	{
	  if (i == anch[d].i0) d++;

	  dpc0 = dpc;
	  norm = xc[p7S_N] + xc[p7S_JJ] + xc[p7S_CC];

	  if (anch[d-1].i0 >= sm->seg[g].ia) // in DOWN
	    {
	      z = 0; while (z < sm->n[i] && sm->k[i][z] < anch[d-1].k0) z++;
	      for (; z < sm->n[i]; z++, dpc += p7S_NSCELLS) 
		{
		  norm += dpc[p7S_ML] + dpc[p7S_MG];
		  norm += dpc[p7S_IL] + dpc[p7S_IG];
		}
	    }

	  if (anch[d].i0 <= sm->seg[g].ib) // in UP
	    {
	      for (z = 0; z < sm->n[i] && sm->k[i][z] < anch[d].k0; z++, dpc+=p7S_NSCELLS)
		{
		  norm += dpc[p7S_ML] + dpc[p7S_MG];
		  norm += dpc[p7S_IL] + dpc[p7S_IG];
		}
	    }

	  /* Renormalize all of row i */
	  for (; dpc0 < dpc; dpc0++)        *dpc0 /= norm;
	  for (x = 0; x < p7S_NXCELLS; x++) xc[x] /= norm;

	  xc += p7S_NXCELLS;
	}
    }
  return eslOK;
}


/*****************************************************************
 * 2. Debugging and development tools
 *****************************************************************/


static void
dump_up_header(FILE *fp, int k1, int k2)
{
  int  width     = 9;
  int k;

  fprintf(fp, "\n# UP component(s) of sparse ASC matrix\n");
  fprintf(fp, "       ");
  for (k = k1; k <= k2;         k++) fprintf(fp, "%*d ", width, k);
  fprintf(fp, "\n");

  fprintf(fp, "       ");
  for (k = k1; k <= k2;        k++) fprintf(fp, "%*.*s ", width, width, "----------");
  fprintf(fp, "\n");
}

static void
dump_down_header(FILE *fp, int k1, int k2)
{
  int  width     = 9;
  int k,x;

  fprintf(fp, "\n# DOWN component(s) of sparse ASC matrix\n");
  fprintf(fp, "       ");
  for (k = k1; k <= k2;         k++) fprintf(fp, "%*d ", width, k);
  for (x = 0;  x < p7S_NXCELLS; x++) fprintf(fp, "%*s ", width, p7_sparsemx_DecodeSpecial(x));
  fprintf(fp, "\n");

  fprintf(fp, "       ");
  for (k = k1; k <= k2;        k++) fprintf(fp, "%*.*s ", width, width, "----------");
  for (x = 0; x < p7S_NXCELLS; x++) fprintf(fp, "%*.*s ", width, width, "----------");
  fprintf(fp, "\n");
}


static void 
dump_up_row(FILE *fp, int i, const P7_SPARSEMASK *sm, const float *dpc, int k0, int k1, int k2, int s)
{
  int  width     = 9;
  int  precision = 4;
  int  k,z;

  /* as a useful side effect, if you pass -1 for k0,
   * the first loop will do nothing, and the second loop 
   * will dump a full blank line from k1..k2: hence, for
   * rows i that aren't connected to an anchor, we pass k0=-1.
   */
  fprintf(fp, "%3d %2s ", i, p7_sparsemx_DecodeState(s));
  for (z = 0, k = k1; k <= k2 && k < k0; k++) {               // from k1..k0-1, you're in potentially valid UP cells
    while (z < sm->n[i] && sm->k[i][z]  < k) z++;
    if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*.*f ", width, precision, *(dpc + z*p7S_NSCELLS + s));  
    else                                     fprintf(fp, "%*s ",   width, "......");
  }
  for ( ; k <= k2; k++) {                                    // from k0..k2, you're outside the UP sector
    while (z < sm->n[i] && sm->k[i][z]  < k) z++;
    if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*s ", width, "...xx.");
    else                                     fprintf(fp, "%*s ", width, "......");
  }
  fprintf(fp, "\n");
}

static void
dump_down_row(FILE *fp, int i, const P7_SPARSEMASK *sm, const float *dpc, int k0, int k1, int k2, int s)
{
  int width     = 9;
  int precision = 4;
  int k,z0,z;

  /* Similar to the side effect mentioned in dump_up_row,
   * if you pass k0 == M+1 and k2 = M+1, you'll dump a completely
   * blank row; hence, that's what we do for rows i that are
   * unconnected to an anchor.
   */
  fprintf(fp, "%3d %2s ", i, p7_sparsemx_DecodeState(s));
  for (z = 0, k = k1; k <= k2 && k < k0; k++) {
    while (z < sm->n[i] && sm->k[i][z]  < k) z++;
    if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*s ", width, "...xx.");
    else                                     fprintf(fp, "%*s ", width, "......");
  }

  z0 = z;
  while (z0 < sm->n[i] && sm->k[i][z0] < k0) z0++;

  for ( ; k <= k2; k++) {
    while (z < sm->n[i] && sm->k[i][z]  < k) z++;
    if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(fp, "%*.*f ", width, precision, *(dpc + (z-z0)*p7R_NSCELLS + s));
    else                                     fprintf(fp, "%*s ",   width, "......");
  }
}


/* Function:  p7_spascmx_Dump()
 * Synopsis:  Dump a sparse ASC matrix for examination.
 */
int
p7_spascmx_Dump(FILE *fp, const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D)
{
  const P7_SPARSEMASK *sm  = asx->sm;
  const float         *dpc;
  const float         *xc;
  int   width     = 9;
  int   precision = 4;
  int   i1        = 0;
  int   i2        = sm->L;
  int   k1        = 0;
  int   k2        = sm->M;
  int   in_seg, ndown, in_x;
  int   i,g,d,z,s;
  int   k0;
  int   in_up;

  /* First pass: print the UP sectors, skip DOWN; do nothing to xc yet  */
  dump_up_header(fp, k1, k2);
  g        = 1;
  in_seg   = FALSE;
  d        = 1;
  ndown    = FALSE;
  dpc      = asx->dp;
  for (i = 0; i <= sm->L; i++)
    {
      if      (i == anch[d].i0)     { ndown = TRUE;  d++;  }  
      else if (sm->n[i] == 0)       { ndown = FALSE;       }  

      if      (i >  sm->seg[g].ib)  { in_seg = FALSE; g++; }
      else if (i == sm->seg[g].ia)  { in_seg = TRUE;       }      

      if (in_seg && anch[d].i0 <= sm->seg[g].ib)  { in_up = TRUE;  k0 = anch[d].k0; }
      else                                        { in_up = FALSE; k0 = -1;         }

      if (ndown)
	{
	  for (z  = 0; z < sm->n[i]; z++)
	    if (sm->k[i][z] >= anch[d-1].k0) break;  // d-1 safe because you can't be in DOWN sector with d=0.
	  dpc += (sm->n[i] - z) * p7S_NSCELLS;       // skip past DOWN rows in the matrix in this pass
	}

      if (i >= i1 && i <= i2) { 
	for (s = 0; s < p7S_NSCELLS; s++) dump_up_row(fp, i, sm, dpc, k0, k1, k2, s);
	fprintf(fp, "\n");
      }
	  
      if (in_up) { // if in UP... skip past those rows
	for (z = 0; z < sm->n[i]; z++)
	  if (sm->k[i][z] >= anch[d].k0) break;    
	dpc += z * p7S_NSCELLS;
      }
    }
  
  /* Second pass: print DOWN sectors and xc rows */
  dump_down_header(fp, k1, k2);
  g        = 1;
  in_seg   = FALSE;
  d        = 1;
  ndown    = FALSE;
  dpc      = asx->dp;
  xc       = asx->xmx;
  for (i = 0; i <= sm->L; i++)
    {
      if      (i == anch[d].i0) { ndown = TRUE;  d++;  }  
      else if (sm->n[i] == 0)   { ndown = FALSE;       }  

      if      (i >  sm->seg[g].ib)  { in_seg = FALSE; g++; }
      else if (i == sm->seg[g].ia)  { in_seg = TRUE;       }      

      k0 = (ndown ? anch[d-1].k0 : sm->M+1);

      if ( (i >= sm->seg[g].ia-1 && anch[d].i0 <= sm->seg[g].ib) ||  ndown) in_x = TRUE;
      else in_x = FALSE;

      if (i >= i1 && i <= i2) {
	dump_down_row(fp, i, sm, dpc, k0, k1, k2, p7S_ML);
	if (in_x) for (s = 0; s < p7S_NXCELLS; s++) fprintf(fp, "%*.*f ", width, precision, xc[s]);
	else      for (s = 0; s < p7S_NXCELLS; s++) fprintf(fp, "%*s ", width, ".....");
	fprintf(fp, "\n");
	for (s = 1; s < p7S_NSCELLS; s++) {
	  dump_down_row(fp, i, sm, dpc, k0, k1, k2, s);
	  fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
      }

      if (ndown) {
	for (z  = 0; z < sm->n[i]; z++)
	  if (sm->k[i][z] >= anch[d-1].k0) break;  
	dpc += (sm->n[i] - z) * p7S_NSCELLS;       // skip past DOWN rows in the matrix in this pass
      }

      if (in_x) 
	xc += p7S_NXCELLS;

      if (in_seg && anch[d].i0 <= sm->seg[g].ib) {
	for (z = 0; z < sm->n[i]; z++)
	  if (sm->k[i][z] >= anch[d].k0) break;    
	dpc += z * p7S_NSCELLS;
      }
    }

  return eslOK;
}


int
p7_spascmx_DumpSpecials(FILE *fp, const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D)
{
  const P7_SPARSEMASK *sm        = asx->sm;
  const float         *xc        = asx->xmx;
  int                  width     = 9;
  int                  precision = 4;
  int d,g,x,i;

  fprintf(fp, "     ");
  for (x = 0; x < p7S_NXCELLS; x++)
    fprintf(fp, "%*s ", width, p7_sparsemx_DecodeSpecial(x));
  fprintf(fp, "\n");

  fprintf(fp, "     ");
  for (x = 0; x < p7S_NXCELLS; x++) 
    fprintf(fp, "%*.*s ", width, width, "----------");
  fprintf(fp, "\n");

  d = 1;
  for (g = 1; g <= sm->S; g++)
    {
      if (anch[d].i0 > sm->seg[g].ib) continue;

      fprintf(fp, "\n");
      fprintf(fp, "%4d ", sm->seg[g].ia-1);
      for (x = 0; x < p7S_NXCELLS; x++)
	fprintf(fp, "%*.*f ", width, precision, xc[x]);
      fprintf(fp, "\n");
      xc += p7S_NXCELLS;
      
      for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)
	{
	  fprintf(fp, "%4d ", i);
	  for (x = 0; x < p7S_NXCELLS; x++)
	    fprintf(fp, "%*.*f ", width, precision, xc[x]);
	  fprintf(fp, "\n");
	  xc += p7S_NXCELLS;
	}
    }
  return eslOK;
}


/* Function:  p7_spascmx_CompareReference()
 * Synopsis:  Compare sparse and reference ASC DP calculations, cell by cell.
 *
 * Purpose:   Compare each cell in the sparse ASC matrix <asx> to the
 *            corresponding value in the reference ASC UP and DOWN
 *            matrices <rxu> and <rxd>, allowing an absolute
 *            difference of up to <tol>. If any value differs by more
 *            than that, return <eslFAIL>.
 *            
 *            You only expect a sparse ASC matrix to compare
 *            identically to a reference ASC calculation if the sparse
 *            mask was set full, marking all cells, as in some unit
 *            tests.
 *            
 *            In some rare cases, it makes sense to compare a sparse
 *            ASC matrix to a normal (non-ASC) reference matrix
 *            <rx>. To do that, you can pass the same matrix ptr <rx>
 *            for <rxu> and <rxd>. An example is in the singlesingle
 *            unit test of sparse_asc_fwdback.
 */
int
p7_spascmx_CompareReference(const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D, const P7_REFMX *rxu, const P7_REFMX *rxd, float tol)
{
  char                 msg[] = "failed comparison of sparse ASC and reference ASC matrices";
  const P7_SPARSEMASK *sm    = asx->sm;
  const float         *dpc   = asx->dp;
  const float         *xc    = asx->xmx;
  int i;                              // index over rows in DP matrices, 0.1..L
  int g      = 1;                     // idx of next or current segment, 1..S, with sentinels. When we enter it, & while we're in it, in_seg = TRUE
  int in_seg = FALSE;                 //  ... => TRUE when starting ia(g), 1st row of segment; => FALSE when we pass ib(g).
  int d      = 1;                     // idx of next domain anchor we will see. Row i may be in sector UP(d), DOWN(d-1). 
  int ndown  = 0;                     // row # of DOWN sector, 1..; becomes 1 when i reaches anchor.
  int z;                              // index over sparse cell list for a row
  int k,s;


  ESL_DASSERT1(( p7R_ML      == p7S_ML      ));  // We assume that the main states come in the same order in both sparse, reference matrices. Spot check that.
  ESL_DASSERT1(( p7R_E       == p7S_E       ));
  ESL_DASSERT1(( p7R_NSCELLS == p7S_NSCELLS ));
  ESL_DASSERT1(( p7R_NXCELLS == p7S_NXCELLS ));

  for (i = 0; i <= sm->L; i++)
    {
      if      (i == anch[d].i0)    { ndown = 1; d++;      }
      else if (sm->n[i] == 0)      { ndown = 0;           }
      else if (ndown)              { ndown++;             }

      if      (i > sm->seg[g].ib)  { in_seg = FALSE; g++; }
      else if (i == sm->seg[g].ia) { in_seg = TRUE;       }

      if (ndown)
	{
	  for (z = 0; z < sm->n[i]; z++)              // Skip ahead in sparse cell list to first z in DOWN sector (anch[d-1].k0).
	    if (sm->k[i][z] >= anch[d-1].k0) break;
	  
	  for (; z < sm->n[i]; z++)
	    {
	      k = sm->k[i][z];
	      for (s = 0; s < p7R_NSCELLS; s++, dpc++)
		if (esl_FCompareAbs(*dpc, P7R_MX(rxd,i,k,s), tol) == eslFAIL) 
		  ESL_FAIL(eslFAIL, NULL, msg);
	    }
	}

      if ( (i >= sm->seg[g].ia-1 && anch[d].i0 <= sm->seg[g].ib) || ndown)
	{
	  for (s = 0; s < p7R_NXCELLS; s++, xc++)
	    if (esl_FCompareAbs(*xc, P7R_XMX(rxd,i,s), tol) == eslFAIL)
	      ESL_FAIL(eslFAIL, NULL, msg);
	}
	   
      /* Cells in UP sector on row i, if any  */
      if (in_seg && anch[d].i0 <= sm->seg[g].ib)
	{
	  for (z = 0; z < sm->n[i]; z++)
	    {
	      k = sm->k[i][z];
	      if (k >= anch[d].k0) break;

	      for (s = 0; s < p7R_NSCELLS; s++, dpc++)
		if (esl_FCompareAbs(*dpc, P7R_MX(rxu,i,k,s), tol) == eslFAIL) 
		  ESL_FAIL(eslFAIL, NULL, msg);
	    }
	}
    }
  return eslOK;
}


/*****************************************************************
 * 3. Debugging and dev: Validation routine
 *****************************************************************/ 

/* Validation of a SPARSEMX structure used for ASC calculations:
 *   we check for the correct pattern of unreached cells 
 *   (see notes in p7_sparsemx.h);
 *   for no NaN's;
 *   and in posterior decoding, that values are probabilities, 0 <= p <= 1.+tol
 *   
 * Validation routines, by Easel spec, must return normal errors (eslFAIL), 
 *   with an error message.
 */

static int
supercell_sc(const float *dpc, int M_used, int I_used, int D_used)
{
  int s;
  for (s = 0; s < p7S_NSCELLS; s++)
    if (isnan(dpc[s])) ESL_FAIL(eslFAIL, NULL, NULL);

  if (! M_used && (dpc[p7S_ML] != -eslINFINITY || dpc[p7S_MG] != -eslINFINITY)) ESL_FAIL(eslFAIL, NULL, NULL);
  if (! I_used && (dpc[p7S_IL] != -eslINFINITY || dpc[p7S_IG] != -eslINFINITY)) ESL_FAIL(eslFAIL, NULL, NULL);
  if (! D_used && (dpc[p7S_DL] != -eslINFINITY || dpc[p7S_DG] != -eslINFINITY)) ESL_FAIL(eslFAIL, NULL, NULL);
  
  return eslOK;
}


static int
supercell_pp(const float *dpc, int M_used, int I_used, int D_used, float tol)
{
  int s;
  for (s = 0; s < p7S_NSCELLS; s++) 
    {
      if (! isfinite(dpc[s]))             ESL_FAIL(eslFAIL, NULL, NULL);
      if (dpc[s] < 0. || dpc[s] > 1.+tol) ESL_FAIL(eslFAIL, NULL, NULL);
    }

  if (! M_used && (dpc[p7S_ML] != 0. || dpc[p7S_MG] != 0.)) ESL_FAIL(eslFAIL, NULL, NULL);
  if (! I_used && (dpc[p7S_IL] != 0. || dpc[p7S_IG] != 0.)) ESL_FAIL(eslFAIL, NULL, NULL);
  if (! D_used && (dpc[p7S_DL] != 0. || dpc[p7S_DG] != 0.)) ESL_FAIL(eslFAIL, NULL, NULL);

  return eslOK;
}

static int
xcell_sc(float xv, int X_used)
{
  if (isnan(xv))                      ESL_FAIL(eslFAIL, NULL, NULL);
  if (! X_used && xv != -eslINFINITY) ESL_FAIL(eslFAIL, NULL, NULL);
  return eslOK;
}

static int
xcell_pp(float xv, int X_used, float tol)
{
  if (! isfinite(xv))          ESL_FAIL(eslFAIL, NULL, NULL);
  if ( xv < 0. || xv > 1.+tol) ESL_FAIL(eslFAIL, NULL, NULL);
  if (! X_used && xv != 0.)    ESL_FAIL(eslFAIL, NULL, NULL);
  return eslOK;
}

static int
decoding_rowsums(const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D)
{
  const P7_SPARSEMASK *sm  = asx->sm;
  const float         *dpc = asx->dp;
  const float         *xc  = asx->xmx;
  float                rowsum;
  float                tol = 0.01;
  int                  d,g,i,z;
  
  d = 1;
  for (g = 1; g <= sm->S; g++)
    {
      if (anch[d].i0 > sm->seg[g].ib) continue;

      /* Row ia-1, just before a segment. */
      rowsum = xc[p7S_N] + xc[p7S_JJ] + xc[p7S_CC];
      if (esl_FCompare(rowsum, 1.0, tol) != eslOK) ESL_FAIL(eslFAIL, NULL, NULL);
      xc += p7S_NXCELLS;

      for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)
	{
	  if (i == anch[d].i0) d++; 

	  rowsum = xc[p7S_N] + xc[p7S_JJ] + xc[p7S_CC];
	  xc += p7S_NXCELLS;

	  if (anch[d-1].i0 >= sm->seg[g].ia)
	    {
	      z = 0; while (z < sm->n[i] && sm->k[i][z] < anch[d-1].k0) z++;
	      for (; z < sm->n[i]; z++, dpc += p7S_NSCELLS)
		{
		  rowsum += dpc[p7S_ML] + dpc[p7S_MG];
		  rowsum += dpc[p7S_IL] + dpc[p7S_IG];
		}
	    }
	  
	  if (anch[d].i0 <= sm->seg[g].ib)
	    {
	      for (z = 0; z < sm->n[i] && sm->k[i][z] < anch[d].k0; z++, dpc += p7S_NSCELLS)
		{
		  rowsum += dpc[p7S_ML] + dpc[p7S_MG];
		  rowsum += dpc[p7S_IL] + dpc[p7S_IG];
		}
	    }

	  if (esl_FCompare(rowsum, 1.0, tol) != eslOK) ESL_FAIL(eslFAIL, NULL, NULL);
	}
    }
  return eslOK;
}
	  
/* decoding_colsums()
 * 
 * Iff profile is in glocal-only mode, then we know that it visits M/D
 * for every position k=1..M, for each of the D domains. Therefore in
 * a full matrix the sum of posterior probability in each column k
 * must be D.
 * 
 * However (alas) in sparse matrices, including sparse ASC, decoding
 * only decodes the DG mass for G->DG1...DGk-1->Mk entry for cells
 * that happen to be in the sparse mask. So in sparse DP, because of
 * this lossage in DG states, colsum tests are less powerful than when
 * we have a complete matrix (as in the reference
 * implementation). Here we only test the bound, colsum <= D, rather
 * than equality, colsum == D.
 * 
 * For any profile, regardless of mode, \sum_i = D for E, B, and L+G.
 */
static int
decoding_colsums(const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D)
{
  const P7_SPARSEMASK *sm       = asx->sm;
  const float         *dpc      = asx->dp;
  const float         *xc       = asx->xmx;
  int                  M        = sm->M;
  float               *colsum   = malloc(sizeof(float) * (M+1));
  float                tol      = 0.01;
  float                xsum[p7S_NXCELLS];
  int   k,s,d,i,z,g;
  int   status;

  for (k = 0; k <= M;          k++) colsum[k] = 0.0f;
  for (s = 0; s < p7S_NXCELLS; s++) xsum[s]   = 0.0f;

  d = 1;
  for (g = 1; g <= sm->S; g++)
    {
      if (anch[d].i0 > sm->seg[g].ib) continue;

      for (s = 0; s < p7S_NXCELLS; s++) xsum[s] += xc[s];
      xc += p7S_NXCELLS;

      for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)
	{
	  if (i == anch[d].i0) d++;

	  for (s = 0; s < p7S_NXCELLS; s++) xsum[s] += xc[s];
	  xc += p7S_NXCELLS;

	  if (anch[d-1].i0 >= sm->seg[g].ia)
	    {
	      z = 0; while (z < sm->n[i] && sm->k[i][z] < anch[d-1].k0) z++;
	      for (; z < sm->n[i]; z++, dpc += p7S_NSCELLS)
		{
		  colsum[sm->k[i][z]] += dpc[p7S_ML] + dpc[p7S_MG];
		  colsum[sm->k[i][z]] += dpc[p7S_DL] + dpc[p7S_DG];
		}
	    }

	  if (anch[d].i0 <= sm->seg[g].ib)
	    {
	      for (z = 0; z < sm->n[i] && sm->k[i][z] < anch[d].k0; z++, dpc += p7S_NSCELLS)
		{
		  colsum[sm->k[i][z]] += dpc[p7S_ML] + dpc[p7S_MG];
		  colsum[sm->k[i][z]] += dpc[p7S_DL] + dpc[p7S_DG];
		}
	    }
	}
    }

  if (colsum[0]   != 0.) ESL_XFAIL(eslFAIL, NULL, NULL);
  for (k = 1; k <= M; k++)
    if (colsum[k] < 0.0 || colsum[k] > (float) D + tol)                    ESL_XFAIL(eslFAIL, NULL, NULL);
  if (esl_FCompareAbs(xsum[p7S_E],               (float) D, tol) != eslOK) ESL_XFAIL(eslFAIL, NULL, NULL);
  if (esl_FCompareAbs(xsum[p7S_B],               (float) D, tol) != eslOK) ESL_XFAIL(eslFAIL, NULL, NULL);
  if (esl_FCompareAbs(xsum[p7S_G] + xsum[p7S_L], (float) D, tol) != eslOK) ESL_XFAIL(eslFAIL, NULL, NULL);

  free(colsum);
  return eslOK;

 ERROR:
  free(colsum);
  return status;
}



static int
validate_decoding(const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D)
{
  const P7_SPARSEMASK *sm  = asx->sm;
  const float         *dpc = asx->dp;
  const float         *xc  = asx->xmx;
  int                  M   = sm->M;
  float                tol = 0.01;
  int in_up, in_down;
  int u1, u2;
  int i0, d2;
  int k0;
  int d;
  int g;
  int k;
  int z;
  int i;

  if (decoding_rowsums(asx, anch, D) != eslOK) return eslFAIL;
  if (decoding_colsums(asx, anch, D) != eslOK) return eslFAIL;

  d = 1;
  for (g = 1; g <= sm->S; g++)
    {
      if (anch[d].i0 > sm->seg[g].ib) continue;

      in_up   = TRUE;
      in_down = FALSE;
      u1      = sm->seg[g].ia;
      u2      = anch[d].i0-1;
      i0      = -1;
      d2      = -1;

      /* SPECIALS on row ia-1:
       *   E only reachable when i is in DOWN sector, so E(ia-1) must be 0.
       *   N only reachable above first anchor
       *   J only reachable between 1st, last anchor
       *   C not reachable above an anchor
       *   B, L, G always reachable in UP sector and its initialization
       *   If J is used, all its mass is in JJ.
       *   CC not reachable because C isn't.
       */
      if (xcell_pp(xc[p7S_E],  0,               tol) != eslOK) return eslFAIL;
      if (xcell_pp(xc[p7S_N],  (d==1),          tol) != eslOK) return eslFAIL;
      if (xcell_pp(xc[p7S_J],  (d>1 && d <= D), tol) != eslOK) return eslFAIL;
      if (xcell_pp(xc[p7S_C],  0,               tol) != eslOK) return eslFAIL;
      if (xcell_pp(xc[p7S_B],  1,               tol) != eslOK) return eslFAIL;
      if (xcell_pp(xc[p7S_L],  1,               tol) != eslOK) return eslFAIL;
      if (xcell_pp(xc[p7S_G],  1,               tol) != eslOK) return eslFAIL;
      if (xcell_pp(xc[p7S_JJ], (d>1 && d <= D), tol) != eslOK) return eslFAIL;
      if (xcell_pp(xc[p7S_CC], 0,               tol) != eslOK) return eslFAIL;
      xc += p7S_NXCELLS;

      for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)
	{
	  if (i == anch[d].i0) 
	    {
	      in_down = TRUE;
	      i0      = anch[d].i0;
	      d2      = ESL_MIN(sm->seg[g].ib, anch[d+1].i0 - 1);

	      if (anch[d+1].i0 > sm->seg[g].ib)
		{
		  in_up = FALSE;
		  u1 = u2 = -1;
		}
	      else 
		{
		  u1 = anch[d].i0+1;
		  u2 = anch[d+1].i0-1;
		}
	      d++;
	    }

	  /* SPECIALS */
	  if (xcell_pp(xc[p7S_E],  in_down,                      tol) != eslOK) return eslFAIL;
	  if (xcell_pp(xc[p7S_N],  (d==1),                       tol) != eslOK) return eslFAIL;  
	  if (xcell_pp(xc[p7S_J],  (d>1 && d <= D),              tol) != eslOK) return eslFAIL;  
	  if (xcell_pp(xc[p7S_C],  (d == D+1),                   tol) != eslOK) return eslFAIL;  
	  if (xcell_pp(xc[p7S_B],  (i >= u1-1 && i <= u2),       tol) != eslOK) return eslFAIL;
	  if (xcell_pp(xc[p7S_L],  (i >= u1-1 && i <= u2),       tol) != eslOK) return eslFAIL;
	  if (xcell_pp(xc[p7S_G],  (i >= u1-1 && i <= u2),       tol) != eslOK) return eslFAIL;
	  if (xcell_pp(xc[p7S_JJ], (i != i0 && (d>1 && d <= D)), tol) != eslOK) return eslFAIL;
	  if (xcell_pp(xc[p7S_CC], (d == D+1),                   tol) != eslOK) return eslFAIL;
	  xc += p7S_NXCELLS;

	  /* DOWN sector */
	  if (in_down)
	    {
	      k0 = anch[d-1].k0;
	      z  = 0; while (z < sm->n[i] && sm->k[i][z] < k0) z++;        
	      for (; z < sm->n[i]; z++, dpc += p7S_NSCELLS) 
		{
		  k  = sm->k[i][z];
		  if        (i == i0)
		    {                           
		      if        (k == k0)   { if (supercell_pp(dpc, 1, 0, 0, tol) != eslOK) return eslFAIL; }
		      else                  { if (supercell_pp(dpc, 0, 0, 1, tol) != eslOK) return eslFAIL; }
		    } 
		  else if (i == d2)    
		    {
		      if      (k == k0)     { if (supercell_pp(dpc, 0, 0, 0, tol) != eslOK) return eslFAIL; }
		      else if (k == k0+1)   { if (supercell_pp(dpc, 1, 0, 0, tol) != eslOK) return eslFAIL; }
		      else                  { if (supercell_pp(dpc, 1, 0, 1, tol) != eslOK) return eslFAIL; }
		    } 
		  else if (i == i0+1) 
		    { 
		      if      (k == k0)     { if (supercell_pp(dpc, 0, 1, 0, tol) != eslOK) return eslFAIL; }
		      else if (k == k0+1)   { if (supercell_pp(dpc, 1, 0, 0, tol) != eslOK) return eslFAIL; }
		      else                  { if (supercell_pp(dpc, 1, 0, 1, tol) != eslOK) return eslFAIL; }
		    }
		  else
		    {
		      if      (k == k0)     { if (supercell_pp(dpc, 0, 1, 0, tol) != eslOK) return eslFAIL; }
		      else if (k == k0+1)   { if (supercell_pp(dpc, 1, 1, 0, tol) != eslOK) return eslFAIL; }
		      else if (k == M)      { if (supercell_pp(dpc, 1, 0, 1, tol) != eslOK) return eslFAIL; }
		      else                  { if (supercell_pp(dpc, 1, 1, 1, tol) != eslOK) return eslFAIL; }
		    }
		}
	    }
		    
	  /* UP sector */
	  if (in_up)
	    {
	      k0 = anch[d].k0;
	      for (z = 0; z < sm->n[i] && sm->k[i][z] < k0; z++, dpc+=p7S_NSCELLS)
		{
		  k = sm->k[i][z];
		  if      (i == u1-1)    
		    {
		      if      (k == k0-1) { if (supercell_pp(dpc, 0, 0, (i == u2), tol) != eslOK) return eslFAIL; }  // i==u2 is tricksy. Possible for UP to have just the u1-1 row, 
		      else                { if (supercell_pp(dpc, 0, 0, 1,         tol) != eslOK) return eslFAIL; }  //   when anchors occur on adjacent rows. In that case, G->DDD->Mk can 
		    }                                                                                                //   pass thru k0-1 cell of u1-1 row to the anchor cell.
		  else if (i == u2)
		    { 
		      if      (k == k0-1) { if (supercell_pp(dpc, 1, 1, 1, tol) != eslOK) return eslFAIL; }
		      else                { if (supercell_pp(dpc, 1, 0, 1, tol) != eslOK) return eslFAIL; }
		    }
		  else if (i == u1)
		    {
		      if      (k == k0-1) { if (supercell_pp(dpc, 1, 0, 0, tol) != eslOK) return eslFAIL; }
		      else                { if (supercell_pp(dpc, 1, 0, 1, tol) != eslOK) return eslFAIL; }
		    }
		  else 
		    {
		      if      (k == k0-1) { if (supercell_pp(dpc, 1, 1, 0, tol) != eslOK) return eslFAIL; }
		      else                { if (supercell_pp(dpc, 1, 1, 1, tol) != eslOK) return eslFAIL; }
		    }
		}
	    }
	} // end loop over i=ia..ib in one segment g
    } // end loop over segments g
  return eslOK;
}

static int
validate_forward(const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D)
{
  const P7_SPARSEMASK *sm  = asx->sm;
  const float         *dpc = asx->dp;
  const float         *xc  = asx->xmx;
  int                  M   = sm->M;
  int in_up, in_down;
  int u1, u2;
  int i0;
  int d,g,k,z,i,k0;

  d = 1;
  for (g = 1; g <= sm->S; g++)
    {
      if (anch[d].i0 > sm->seg[g].ib) continue;

      in_up   = TRUE;
      in_down = FALSE;
      u1      = sm->seg[g].ia;
      u2      = anch[d].i0-1;
      i0      = -1;

      if (xcell_sc(xc[p7S_E],  0)      != eslOK) return eslFAIL;
      if (xcell_sc(xc[p7S_N],  (d==1)) != eslOK) return eslFAIL;
      if (xcell_sc(xc[p7S_J],  (d>1))  != eslOK) return eslFAIL; // In Fwd direction, J may be finite even after last anchor
      if (xcell_sc(xc[p7S_C],  (d>1))  != eslOK) return eslFAIL; // and similarly, C may be finite even with anchors to come
      if (xcell_sc(xc[p7S_B],  1)      != eslOK) return eslFAIL;
      if (xcell_sc(xc[p7S_L],  1)      != eslOK) return eslFAIL;
      if (xcell_sc(xc[p7S_G],  1)      != eslOK) return eslFAIL;
      if (xcell_sc(xc[p7S_JJ], 0)      != eslOK) return eslFAIL; // JJ, CC only used in Decoding; -inf in F/B
      if (xcell_sc(xc[p7S_CC], 0)      != eslOK) return eslFAIL;
      xc += p7S_NXCELLS;

      for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)
	{
	  if (i == anch[d].i0) 
	    {
	      in_down = TRUE;
	      i0      = anch[d].i0;

	      if (anch[d+1].i0 > sm->seg[g].ib)
		{
		  in_up = FALSE;
		  u1 = u2 = -1;
		}
	      else 
		{
		  u1 = anch[d].i0+1;
		  u2 = anch[d+1].i0-1;
		}
	      d++;
	    }

	  if (xcell_sc(xc[p7S_E],  in_down)                != eslOK) return eslFAIL;
	  if (xcell_sc(xc[p7S_N],  (d==1))                 != eslOK) return eslFAIL;  
	  if (xcell_sc(xc[p7S_J],  (d>1))                  != eslOK) return eslFAIL;  
	  if (xcell_sc(xc[p7S_C],  (d>1))                  != eslOK) return eslFAIL;  
	  if (xcell_sc(xc[p7S_B],  1)                      != eslOK) return eslFAIL; // in Fwd direction, can always enter {NJ}->B...
	  if (xcell_sc(xc[p7S_L],  1)                      != eslOK) return eslFAIL; //   ... and thence B->{LG}, even though Bck will be -inf on DOWN-only rows
	  if (xcell_sc(xc[p7S_G],  1)                      != eslOK) return eslFAIL;
	  if (xcell_sc(xc[p7S_JJ], 0)                      != eslOK) return eslFAIL;
	  if (xcell_sc(xc[p7S_CC], 0)                      != eslOK) return eslFAIL;
	  xc += p7S_NXCELLS;

	  if (in_down)
	    {
	      k0 = anch[d-1].k0;
	      z  = 0; while (z < sm->n[i] && sm->k[i][z] < k0) z++;
	      for (; z < sm->n[i]; z++, dpc += p7S_NSCELLS) 
		{
		  k  = sm->k[i][z];
		  if        (i == i0)
		    {                           
		      if        (k == k0)   { if (supercell_sc(dpc, 1, 0, 0) != eslOK) return eslFAIL; }
		      else                  { if (supercell_sc(dpc, 0, 0, 1) != eslOK) return eslFAIL; }
		    } 
		  else if (i == i0+1) 
		    { 
		      if      (k == k0)     { if (supercell_sc(dpc, 0, 1, 0) != eslOK) return eslFAIL; }   // if k == M too, it'd be 0 0 0
		      else if (k == k0+1)   { if (supercell_sc(dpc, 1, 0, 0) != eslOK) return eslFAIL; }
		      else                  { if (supercell_sc(dpc, 1, 0, 1) != eslOK) return eslFAIL; }
		    }
		  else 
		    {
		      if      (k == k0)     { if (supercell_sc(dpc, 0, 1, 0) != eslOK) return eslFAIL; }   // if k == k0   && k == M, it'd be 0 0 0; so we're underchecking that case
		      else if (k == k0+1)   { if (supercell_sc(dpc, 1, 1, 0) != eslOK) return eslFAIL; }   // if k == k0+1 && k == M, it'd be 1 0 0 
		      else if (k == M)      { if (supercell_sc(dpc, 1, 0, 1) != eslOK) return eslFAIL; }
		      else                  { if (supercell_sc(dpc, 1, 1, 1) != eslOK) return eslFAIL; }
		    } 
		}
	    }
		    
	  /* UP sector */
	  if (in_up)
	    {
	      k0 = anch[d].k0;
	      for (z = 0; z < sm->n[i] && sm->k[i][z] < k0; z++, dpc+=p7S_NSCELLS)
		{
		  k = sm->k[i][z];
		  if      (i == u1-1)     { if (supercell_sc(dpc, 0, 0, 0) != eslOK) return eslFAIL; }
		  else if (i == u2)
		    { 
		      if      (k == 1)    { if (supercell_sc(dpc, 1, 1, 0) != eslOK) return eslFAIL; }
		      else                { if (supercell_sc(dpc, 1, 1, 1) != eslOK) return eslFAIL; }
		    }
		  else if (i == u1)
		    {
		      if      (k == 1)    { if (supercell_sc(dpc, 1, 0, 0) != eslOK) return eslFAIL; }
		      else                { if (supercell_sc(dpc, 1, 0, 1) != eslOK) return eslFAIL; }
		    }
		  else 
		    {
		      if      (k == 1)    { if (supercell_sc(dpc, 1, 1, 0) != eslOK) return eslFAIL; }
		      else                { if (supercell_sc(dpc, 1, 1, 1) != eslOK) return eslFAIL; }
		    }
		}
	    }
	} // end loop over i=ia..ib in one segment g
    } // end loop over segments g
  return eslOK;
}


static int
validate_backward(const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D)
{
  const P7_SPARSEMASK *sm  = asx->sm;
  const float         *dpc = asx->dp;
  const float         *xc  = asx->xmx;
  int                  M   = sm->M;
  int in_up, in_down;
  int u1, u2;
  int d2;
  int d,g,k,z,i,k0;

  d = 1;
  for (g = 1; g <= sm->S; g++)
    {
      if (anch[d].i0 > sm->seg[g].ib) continue;

      in_up   = TRUE;
      in_down = FALSE;
      u1      = sm->seg[g].ia;
      u2      = anch[d].i0-1;
      d2      = -1;

      if (xcell_sc(xc[p7S_E],  1)        != eslOK) return eslFAIL; // E->{JC} means E always reachable in Backwards direction, like B in the Fwd direction
      if (xcell_sc(xc[p7S_N],  (d<=D))   != eslOK) return eslFAIL; // N is reachable after we pass 1st anchor on the way back...
      if (xcell_sc(xc[p7S_J],  (d<=D))   != eslOK) return eslFAIL; //   ... as is J.
      if (xcell_sc(xc[p7S_C],  (d==D+1)) != eslOK) return eslFAIL; // C only reachable below last anchor
      if (xcell_sc(xc[p7S_B],  1)        != eslOK) return eslFAIL;
      if (xcell_sc(xc[p7S_L],  1)        != eslOK) return eslFAIL;
      if (xcell_sc(xc[p7S_G],  1)        != eslOK) return eslFAIL;
      if (xcell_sc(xc[p7S_JJ], 0)        != eslOK) return eslFAIL; // JJ, CC only used in Decoding; -inf in F/B
      if (xcell_sc(xc[p7S_CC], 0)        != eslOK) return eslFAIL;
      xc += p7S_NXCELLS;

      for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)
	{
	  if (i == anch[d].i0) 
	    {
	      in_down = TRUE;
	      d2      = ESL_MIN(sm->seg[g].ib, anch[d+1].i0 - 1);

	      if (anch[d+1].i0 > sm->seg[g].ib)
		{
		  in_up = FALSE;
		  u1 = u2 = -1;
		}
	      else 
		{
		  u1 = anch[d].i0+1;
		  u2 = anch[d+1].i0-1;
		}
	      d++;
	    }

	  if (xcell_sc(xc[p7S_E],  1)                      != eslOK) return eslFAIL;
	  if (xcell_sc(xc[p7S_N],  (d<=D))                 != eslOK) return eslFAIL;  
	  if (xcell_sc(xc[p7S_J],  (d<=D))                 != eslOK) return eslFAIL;  
	  if (xcell_sc(xc[p7S_C],  (d==D+1))               != eslOK) return eslFAIL;  
	  if (xcell_sc(xc[p7S_B],  (i >= u1-1 && i <= u2)) != eslOK) return eslFAIL;
	  if (xcell_sc(xc[p7S_L],  (i >= u1-1 && i <= u2)) != eslOK) return eslFAIL;
	  if (xcell_sc(xc[p7S_G],  (i >= u1-1 && i <= u2)) != eslOK) return eslFAIL;
	  if (xcell_sc(xc[p7S_JJ], 0)                      != eslOK) return eslFAIL;
	  if (xcell_sc(xc[p7S_CC], 0)                      != eslOK) return eslFAIL;
	  xc += p7S_NXCELLS;

	  if (in_down)
	    {
	      k0 = anch[d-1].k0;
	      z  = 0; while (z < sm->n[i] && sm->k[i][z] < k0) z++;
	      for (; z < sm->n[i]; z++, dpc += p7S_NSCELLS) 
		{
		  k  = sm->k[i][z];
		  if (i == d2)              { if (supercell_sc(dpc, 1, 0, 1) != eslOK) return eslFAIL; }
		  else                      
		    {
		      if      (k == M)      { if (supercell_sc(dpc, 1, 0, 1) != eslOK) return eslFAIL; } 
		      else                  { if (supercell_sc(dpc, 1, 1, 1) != eslOK) return eslFAIL; }
		    } 
		}
	    }
		    
	  /* UP sector */
	  if (in_up)
	    {
	      k0 = anch[d].k0;
	      for (z = 0; z < sm->n[i] && sm->k[i][z] < k0; z++, dpc+=p7S_NSCELLS)
		{
		  k = sm->k[i][z];
		  if      (i == u1-1)     { if (supercell_sc(dpc, 0, 0, 0) != eslOK) return eslFAIL; }
		  else if (i == u2)
		    { 
		      if      (k == k0-1) { if (supercell_sc(dpc, 1, 1, 1) != eslOK) return eslFAIL; }
		      else                { if (supercell_sc(dpc, 1, 0, 1) != eslOK) return eslFAIL; }
		    }
		  else
		    {
		      if      (k == k0-1) { if (supercell_sc(dpc, 1, 1, 0) != eslOK) return eslFAIL; }
		      else                { if (supercell_sc(dpc, 1, 1, 1) != eslOK) return eslFAIL; }
		    }
		}
	    }
	} // end loop over i=ia..ib in one segment g
    } // end loop over segments g
  return eslOK;
}

/* Function:  p7_spascmx_Validate()
 * Synopsis:  Validate a sparse ASC DP matrix.
 *
 * Purpose:   
 * 
 *            Caller may optionally provide <gm>, a pointer to the
 *            profile that was used in calculating the <asx> matrix.
 *            If <gm> is provided (non-<NULL>), validation tests that
 *            require the profile to be in a particular mode can be
 *            performed. Currently there's only one such validation
 *            test: the "colsum" test for a decoding matrix.
 *
 * Args:      gm  - profile that <asx> was computed with. Optional; can pass <NULL>.
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_spascmx_Validate(const P7_SPARSEMX *asx, const P7_ANCHOR *anch, int D, char *errbuf)
{
  if (asx->type == p7S_ASC_DECODE) 
    {
      if (validate_decoding(asx, anch, D) != eslOK) 
	ESL_FAIL(eslFAIL, errbuf, "Sparse ASC posterior decoding matrix failed validation");
    } 
  else if (asx->type == p7S_ASC_FWD)
    {
      if (validate_forward(asx, anch, D) != eslOK) 
	ESL_FAIL(eslFAIL, errbuf, "Sparse ASC forward matrix failed validation");
    }
  else if (asx->type == p7S_ASC_BCK)
    {
      if (validate_backward(asx, anch, D) != eslOK) 
	ESL_FAIL(eslFAIL, errbuf, "Sparse ASC forward matrix failed validation");
    }
  else 
    ESL_FAIL(eslFAIL, errbuf, "validation of other sparse ASC matrix types unimplemented");
  
  return eslOK;
}


/*****************************************************************
 * x. Statistics collection driver.
 *****************************************************************/
#ifdef p7SPASCMX_STATS
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type           default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,   FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",                   0 },
  { "-s",          eslARG_INT,      "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "statistics collection on sparse ASC DP matrices";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_ANCHORS     *anch    = p7_anchors_Create();
  P7_REFMX       *rxf     = NULL;
  P7_REFMX       *rxd     = NULL;
  P7_REFMX       *afu     = NULL;
  P7_REFMX       *afd     = NULL;
  P7_TRACE       *vtr     = NULL;
  P7_FILTERMX    *fx      = p7_filtermx_Create(100);
  P7_CHECKPTMX   *cx      = p7_checkptmx_Create(100, 100, ESL_MBYTES(p7_SPARSIFY_RAMLIMIT));
  P7_SPARSEMASK  *sm      = p7_sparsemask_Create(100, 100);
  float          *wrk     = NULL;
  P7_ANCHORHASH  *ah      = p7_anchorhash_Create();
  float           fsc, vsc, asc;
  int64_t         dalloc;
  int             xalloc, spascmxsize;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_oprofile_Convert(gm, om);

  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Allocate DP matrices and tracebacks */
  vtr = p7_trace_Create();
  rxf = p7_refmx_Create(gm->M, 100);
  rxd = p7_refmx_Create(gm->M, 100);
  afu = p7_refmx_Create(gm->M, 100);
  afd = p7_refmx_Create(gm->M, 100);

  /* For each target sequence... */
  while (( status = esl_sqio_Read(sqfp, sq)) == eslOK) 
    {
      /* Set the profile and null model's target length models */
      p7_bg_SetLength           (bg, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_oprofile_ReconfigLength(om, sq->n);

      if (( status = p7_pipeline_AccelerationFilter(sq->dsq, sq->n, om, bg, fx, cx, sm)) == eslOK)
	{
	  /* First pass analysis */
	  p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf, vtr, &vsc);
	  p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc);
	  p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd, NULL);   
	  p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd);   

	  /* Find most probable anchor set */
	  p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, vtr, &wrk, ah,
			       afu, afd, anch, &asc, NULL, NULL);

      
	  printf("%-15s  %-30s  ", gm->name, sq->name);

	  /* Reference matrix size: cells, special rows, total in bytes */
	  printf("%10d %10d %10d ",
		 (int) (sq->n+1) * (gm->M+1),
		 (int) (sq->n+1),
		 (int) p7_refmx_MinSizeof(gm->M, sq->n));
		 
	  /* Sparse matrix size: cells, special rows, total in bytes */
	  printf("%10d %10d %10d ",
		 (int) sm->ncells,
		 (int) (sm->nrow + sm->S),
		 (int) p7_sparsemx_MinSizeof(sm));

	  /* Sparse ASC matrix size: UP cells, DOWN cells, total cells, special rows, total in bytes */
	  spascmxsize = p7_spascmx_MinSizeof(sm, anch->a, anch->D, &dalloc, &xalloc);
	  printf("%10d %10d %10d\n",
		 (int) dalloc,
		 (int) xalloc,
		 (int) spascmxsize);

	  p7_trace_Reuse(vtr);
	  p7_refmx_Reuse(rxf);   p7_refmx_Reuse(rxd);
	  p7_refmx_Reuse(afu);   p7_refmx_Reuse(afd);
	  p7_anchorhash_Reuse(ah);
	  p7_anchors_Reuse(anch);
	}

      esl_sq_Reuse(sq);
      p7_sparsemask_Reuse(sm);
    }
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslEOF)     p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  esl_sqfile_Close(sqfp);
  p7_filtermx_Destroy(fx);
  p7_checkptmx_Destroy(cx);
  p7_sparsemask_Destroy(sm);
  p7_anchorhash_Destroy(ah);
  if (wrk) free(wrk);
  p7_trace_Destroy(vtr);
  p7_refmx_Destroy(afd);  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(rxd);  p7_refmx_Destroy(rxf);
  p7_anchors_Destroy(anch);
  esl_sq_Destroy(sq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPASCMX_STATS*/

