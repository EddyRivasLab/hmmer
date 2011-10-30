/* P7_GMXCHK implementation
 * Checkpointed forward/backward dynamic programming matrix.
 * Derived from P7_GMX structure; see p7_gmx.[ch]
 * 
 * Contents:
 *   1. The <P7_GMXCHK> object.
 *   2. Internal (static) routines.
 *   x. Copyright and license information.
 *   
 * References:
 *   SRE J8/109-112, Oct 2011: implementation plan.
 */
#include "p7_config.h"

#include <math.h>

#include "hmmer.h"
#include "p7_gmxchk.h"


static int    set_row_layout  (P7_GMXCHK *gxc, int allocL, int maxR);
static void   set_full        (P7_GMXCHK *gxc, int L);
static void   set_checkpointed(P7_GMXCHK *gxc, int L, int R);
static void   set_redlined    (P7_GMXCHK *gxc, int L, double minR);
static double minimum_rows     (int L);
static double checkpointed_rows(int L, int R);


/*****************************************************************
 *= 1. The <P7_GMXCHK> object
 *****************************************************************/

/* Function:  p7_gmxchk_Create()
 * Synopsis:  Allocate a new <P7_GMXCHK> matrix.
 *
 * Purpose:   Allocate a new <P7_GMXCHK> matrix sufficient to
 *            hold the checkpointed forward and two-row backwards
 *            DP matrices for a comparison of a query model of
 *            length <M> and a target sequence of length <L>.
 *            
 *            Try to keep allocation within <ramlimit> MB in memory.
 *            <ramlimit=128>, for example, means a recommended RAM
 *            limit of 128MB. The allocation can exceed this, if the
 *            <MxL> comparison requires it (if so, the matrix is fully
 *            checkpointed, minimizing the allocation); but any
 *            subsequent <p7_gmxchk_GrowTo()> call attempting to reuse
 *            the matrix will try to reallocate it back downwards to
 *            the <ramlimit>.
 *            
 *            Choice for <ramlimit> should take into account how many
 *            parallel threads there are, each with its own
 *            <P7_GMXCHK> matrix allocation.
 *            
 *            By design spec, <M> and <L> are $\leq$ 100000.  Because
 *            some size calculations are in int32, <ramlimit> must be
 *            $\leq$ 2000 (2 GB).
 *
 * Args:      M         - query profile size in consensus positions (<=100000)
 *            L         - target sequence length in residues (<=100000)
 *            ramlimit  - recommended memory limit in MB (<=2000)
 *
 * Returns:   ptr to new <P7_GMXCHK> object on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_GMXCHK *
p7_gmxchk_Create(int M, int L, int ramlimit)
{
  P7_GMXCHK *gxc = NULL;
  int        maxR;
  int        r;
  int        status;

  /* Validity of integer variable ranges may depend on design spec:                   */
  ESL_DASSERT1( (M        <= 100000) ); /* design spec says, model length M <= 100000 */
  ESL_DASSERT1( (L        <= 100000) ); /*           ... and,  seq length L <= 100000 */
  ESL_DASSERT1( (ramlimit <= 2000)   ); /* this, because calculations are in int32    */

  /* Level 1 allocation: the structure itself */
  ESL_ALLOC(gxc, sizeof(P7_GMXCHK));
  gxc->dp      = NULL;
  gxc->xmx     = NULL;
  gxc->dp_mem  = NULL;

  /* Set allocR, R{abc}, L{abc} fields: row layout, full vs. checkpointed */
  gxc->R0          = 3;	                  /* fwd[0]; bck[prv,cur] */
  gxc->allocW      = M + 1;               /* allocM+1 because we're [0..allocM] in DP columns in generic mx. */
  gxc->ncell_limit = (ramlimit * 1048576L) / (p7G_NSCELLS *  sizeof(float));
  maxR             = gxc->ncell_limit / gxc->allocW;  /* this int32 calculation is why ramlimit must be <= 2000 (2GB) */
  set_row_layout(gxc, L, maxR);
  gxc->allocR      = gxc->R0 + gxc->Ra + gxc->Rb + gxc->Rc;
  gxc->validR      = gxc->allocR;

  /* Level 2 allocations: row pointers and dp cell memory */
  gxc->ncells = gxc->allocR * gxc->allocW;
  ESL_ALLOC( gxc->dp_mem, sizeof(float)   * gxc->ncells * p7G_NSCELLS);
  ESL_ALLOC( gxc->dp,     sizeof(float *) * gxc->allocR);
  for (r = 0; r < gxc->allocR; r++) 
    gxc->dp[r] = gxc->dp_mem + (r * gxc->allocW * p7G_NSCELLS);

  ESL_ALLOC( gxc->xmx, sizeof(float) * (L+1) * p7G_NXCELLS);
  gxc->allocXR = L+1;

  gxc->M      = 0;
  gxc->L      = 0;
  gxc->R      = 0;
  return gxc;

 ERROR:
  if (gxc) p7_gmxchk_Destroy(gxc);
  return NULL;
}

/* Function:  p7_gmxchk_GrowTo()
 * Synopsis:  Resize a checkpointed matrix structure for new comparison.
 *
 * Purpose:   Given an existing checkpointed matrix structure <gxc>,
 *            and the dimensions <M> and <L> of a new comparison,
 *            reallocate and reinitialize <gxc>.
 *
 *            Essentially the same as free'ing the previous matrix
 *            and creating a new one -- but trying to minimize 
 *            expensive memory allocation/reallocation calls.
 *            
 *            If <gxc> is redlined (over its recommended allocation)
 *            and the new problem size <M,L> can fit in the recommended
 *            allocation, <gxc> is reallocated to the smaller size.
 *            
 * Args:      gxc   - existing checkpointed matrix
 *            M     - new query profile length
 *            L     - new target sequence length         
 * 
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if an allocation fails. The state of <gxc> is
 *            now undefined, and the caller should not use it. 
 */
int
p7_gmxchk_GrowTo(P7_GMXCHK *gxc, int M, int L)
{
  int minR_chk      = (int) ceil(minimum_rows(L)) + gxc->R0; /* minimum number of DP rows needed  */
  int reset_dp_ptrs = FALSE;
  int maxR;
  int r;
  int status;

  /* Are the current allocations sufficient for the new problem? */
  if (M + 1 <= gxc->allocW &&  L + 1 <= gxc->allocXR &&  gxc->ncells <= gxc->ncell_limit)
    {
      if      (L + gxc->R0 <= gxc->validR) { set_full        (gxc, L);              return eslOK; }
      else if (minR_chk   <=  gxc->validR) { set_checkpointed(gxc, L, gxc->validR); return eslOK; }
    }

  /* Do individual matrix rows need to expand? */
  if (M+1 > gxc->allocW) 
    {
      gxc->allocW   = M+1;
      gxc->validR   = gxc->ncells / gxc->allocW; /* validR must be <= allocR */
      reset_dp_ptrs = TRUE;
    }

  /* Does matrix need reallocation, either up or down? */
  maxR  = gxc->ncell_limit / gxc->allocW;                       /* max rows if we use up to the recommended allocation size.      */
  if ( (gxc->ncells > gxc->ncell_limit && minR_chk <= maxR) ||  /* we were redlined, and recommended alloc will work: so downsize */
       minR_chk > gxc->validR)				        /* not enough memory for needed rows: so upsize                   */
    {
      set_row_layout(gxc, L, maxR); 
      gxc->validR = gxc->R0 + gxc->Ra + gxc->Rb + gxc->Rc; /* this may be > allocR now; we'll reallocate dp[] next, if so     */
      gxc->ncells = gxc->validR * gxc->allocW;
      ESL_REALLOC(gxc->dp_mem, sizeof(float) * gxc->ncells * p7G_NSCELLS);
      reset_dp_ptrs = TRUE;
    }
  
  /* Does the array of row ptrs need reallocation? */
  if (gxc->validR > gxc->allocR)
    {
      ESL_REALLOC(gxc->dp, sizeof(float *) * gxc->validR);
      gxc->allocR = gxc->validR;
      reset_dp_ptrs = TRUE;
    }

  /* Do the row ptrs need to be reset? */
  if (reset_dp_ptrs)
    for (r = 0; r < gxc->validR; r++)
      gxc->dp[r] = gxc->dp_mem + (r * gxc->allocW + p7G_NSCELLS);

  /* Does the xmx special array need reallocation? */
  if (L+1 > gxc->allocXR)
    {
      ESL_REALLOC( gxc->xmx, sizeof(float) * (L+1) * p7G_NXCELLS);
      gxc->allocXR = L+1;
    }

  gxc->M = 0;
  gxc->L = 0;
  gxc->R = 0;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_gmxchk_Destroy()
 * Synopsis:  Frees a checkpointed generic DP matrix.
 *
 * Purpose:   Free the checkpointed generic DP matrix <gxc>.
 *            <gxc> may be <NULL> or incompletely allocated.
 */
void
p7_gmxchk_Destroy(P7_GMXCHK *gxc)
{
  if (gxc)
    {
      if (gxc->dp_mem) free(gxc->dp_mem);
      if (gxc->dp)     free(gxc->dp);
      if (gxc->xmx)    free(gxc->xmx);
      free(gxc);
    }
}
/*-------------- end, P7_GMXCHK object --------------------------*/



/*****************************************************************
 * 2. Internal routines
 *****************************************************************/

/* set_row_layout()
 * Determines the size of the a,b,c regions ("all", "between",
 * "checkpointed") of rows in the DP matrix. Sets the R{0abc},
 * L{abc}, <allocW> and <allocR> fields in the <gxc> structure.
 *
 * Caller has already set gx->allocW, the row width, to
 * something >= allocM+1.
 * Caller has also set gxc->R0, the number of extra rows
 * it needs (two bck[cur,prv], one boundary fwd[0]).
 * 
 * Design spec says <allocM> <= 100000, <allocL> <= 100000.
 * 
 * <ramlimit> is the # of MB of RAM we're "allowed" to use.
 * We can exceed this for one comparison if necessary,
 * but the next <_Reuse()> call will bring the allocation 
 * back down. It should not exceed 2000 (2 GB).
 * 
 * Three possibilities:
 *  1. A full matrix fits into our recommended max memory use.
 *  2. A checkpointed matrix fits into our recommended memory.
 *     Allocate right up to our redline.
 *  3. We can't satisfy the recommended max memory. Allocate a
 *     checkpointed matrix anyway, temporarily.  The next _Reuse()
 *     pulls the allocation back down to redline.
 */
static int
set_row_layout(P7_GMXCHK *gxc, int allocL, int maxR)
{
  double Rbc      = minimum_rows(allocL);               
  int    minR_chk = gxc->R0 + (int) ceil(Rbc);	          /* min # rows we need for checkpointing          */
  int    minR_all = gxc->R0 + allocL;		          /* min # rows we need for full matrix            */

  if      (minR_all <= maxR) set_full        (gxc, allocL);
  else if (minR_chk <= maxR) set_checkpointed(gxc, allocL, maxR);
  else                       set_redlined    (gxc, allocL, Rbc);

  return eslOK;
}

static void
set_full(P7_GMXCHK *gxc, int L)
{
  gxc->Ra     = L;
  gxc->Rb     = 0;
  gxc->Rc     = 0;
  gxc->La     = L;
  gxc->Lb     = 0;
  gxc->Lc     = 0;
}

static void
set_checkpointed(P7_GMXCHK *gxc, int L, int R)
{
  double Rbc = checkpointed_rows(L, R-gxc->R0);
  double Rc  = floor(Rbc);

  gxc->Rc     = (int) Rc;
  gxc->Rb     = (Rbc > Rc ? 1 : 0);
  gxc->Ra     = R - gxc->Rb - gxc->Rc - gxc->R0;
  gxc->Lc     = ((gxc->Rc + 2) * (gxc->Rc + 1)) / 2 - 1;
  gxc->La     = gxc->Ra;
  gxc->Lb     = L - gxc->La - gxc->Lc;
}   

static void
set_redlined(P7_GMXCHK *gxc, int L, double minR)
{
  double Rc = floor(minR);

  gxc->Rc     = (int) Rc;
  gxc->Rb     = (minR > Rc ? 1 : 0);
  gxc->Ra     = 0;
  gxc->Lc     = ((gxc->Rc + 2) * (gxc->Rc + 1)) / 2 - 1;
  gxc->La     = 0;
  gxc->Lb     = L - gxc->La - gxc->Lc;
}

/* minimum_rows()
 * 
 * Calculate the minimum number of rows needed to checkpoint a
 * forward matrix for a sequence of length L, exclusive of 
 * other constant row overhead (R0: two backwards rows, fwd[0]
 * initial row).
 * 
 * This value is a float; if it has a fraction, a partial checkpoint
 * segment is needed, as in this typical code:
 *    Rbc  = minimum_rows(L);
 *    Rc   = floor(Rbc);
 *    Rb   = (Rbc > Rc ? 1 : 0);
 *    minR = (int) ceil(Rbc);    // or, Rc+Rb
 */
static double 
minimum_rows(int L)
{
  return (sqrt(9. + 8. * (double) L) - 3.) / 2.;
}

static double
checkpointed_rows(int L, int R)
{
  return (sqrt(1. + 8. * (double) (L - R)) - 1.) / 2.;
}

/*----------------- end, internals ------------------------------*/

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

