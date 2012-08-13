/* 
 * Data structures used by sparse dynamic programming.
 *
 * P7_SPARSEMASK and P7_SPARSEMX work together. P7_SPARSEMASK defines
 * which cells are included in the sparse DP matrix, and P7_SPARSEMX
 * is the matrix. The caller constructs a P7_SPARSEMASK first, then
 * uses it to create or reinitialize a P7_SPARSEMX, before running any
 * sparse DP routines (see sparse_fwdback.c).
 * 
 * For comparison see the reference DP implementation (p7_refmx.[ch],
 * reference_fwdback.c). Both the reference and the sparse DP
 * implementations work with a dual-mode profile <P7_PROFILE>, either
 * in the production code's back end (downstream of our fast vector
 * filters) or in testing.
 * 
 * Contents:
 *   1. P7_SPARSEMASK; defines cells to be included in sparse DP matrix
 *   2. API for constructing a sparse DP mask
 *   3. Debugging tools for P7_SPARSEMASK
 *   4. P7_SPARSEMX; sparse DP matrix, as specified by a given mask
 *   5. API for extracting information from a sparse DP matrix
 *   6. Debugging tools for P7_SPARSEMX
 *   7. Validation of a P7_SPARSEMX
 *   8. Copyright and license information  
 */

#include "p7_config.h"

#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_sparsemx.h"

/*****************************************************************
 * 1. P7_SPARSEMASK; defines cells to be included in sparse DP matrix
 *****************************************************************/

/* Function:  p7_sparsemask_Create()
 * Synopsis:  Creates a new P7_SPARSEMASK object.
 *
 * Purpose:   Create a new <P7_SPARSEMASK> for a comparison of a profile
 *            of length <M> to a sequence of length <L>. Return a ptr to the
 *            new object.
 *
 *            The allocation will generally be for (much) less than <ML> cells;
 *            the API for creating the sparse mask will grow the structure
 *            appropriately. The structure does require at least $O(M)$ cells
 *            of temporary storage, in four "slots" used to sort input
 *            from striped vector code.  
 *
 * Args:      M       - model length
 *            L       - sequence length
 *
 * Returns:   a pointer to a new <P7_SPARSEMASK>
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_SPARSEMASK *
p7_sparsemask_Create(int M, int L)
{
  P7_SPARSEMASK *sm             = NULL;
  int            default_ialloc = 4;
  int64_t        default_kalloc = 4096;
  int            i,r;
  int            status;

  ESL_ALLOC(sm, sizeof(P7_SPARSEMASK));
  sm->L      = L;
  sm->M      = M;
  sm->Q      = p7O_NQF(M);  // approx M/4, for striped vectors of four floats
		
  sm->i      = NULL;
  sm->k      = NULL;
  sm->n      = NULL;
  sm->kmem   = NULL;

  sm->nseg    = 0;
  sm->nrow    = 0;
  sm->ncells  = 0;
  sm->last_i  = L+1;     	 // sentinel to assure StartRow() is called in reverse L..1 order 
  for (r = 0; r < p7_VNF; r++) 
    sm->last_k[r]  = -1;        // sentinels to assure StartRow() is called before Add() 
  /* sn[] are initialized for each sparse row by _StartRow() */

  /* if Ws is really large, we might already know we need a
   * bigger-than-default allocation, just to enable the slots.
   * Rather than allocating the default and letting StartRow()
   * reallocate for the slots, go ahead and figure this out now.
   */
  sm->kalloc = default_kalloc;
  while (sm->kalloc < p7_VNF*sm->Q) sm->kalloc *= 2;

  sm->ralloc   = L+1;		
  sm->ialloc   = default_ialloc;

  sm->n_krealloc = 0;
  sm->n_irealloc = 0;
  sm->n_rrealloc = 0;

  ESL_ALLOC(sm->i,    sizeof(int)   * 2 * sm->ialloc); // *2 because ia,ib pairs 
  ESL_ALLOC(sm->k,    sizeof(int *) * sm->ralloc);
  ESL_ALLOC(sm->n,    sizeof(int)   * sm->ralloc);
  ESL_ALLOC(sm->kmem, sizeof(int)   * sm->kalloc);

  sm->k[0]   = NULL;		// always. 
  for (i = 0; i <= L; i++)	// n[0] will always be 0; n[i=1..L] initialized to 0, then count as cells are added 
    sm->n[i] = 0;
  return sm;

 ERROR:
  p7_sparsemask_Destroy(sm);
  return NULL;
}

/* Function:  p7_sparsemask_Reinit()
 * Synopsis:  Reinitialize an existing P7_SPARSEMASK for a new comparison.
 *
 * Purpose:   Same as a <_Create()>, but reusing an existing 
 *            <P7_SPARSEMASK> to minimize reallocation calls.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_sparsemask_Reinit(P7_SPARSEMASK *sm, int M, int L)
{
  int i,r;
  int status;

  sm->L  = L;
  sm->M  = M;
  sm->Q  = p7O_NQF(M);
		
  /* i[], kmem stay at their previous ialloc, kalloc
   * but do we need to reallocate rows for k[] and n[]? 
   */
  if (sm->ralloc < L+1) {
    ESL_REALLOC(sm->k, sizeof(int *) * (L+1));
    ESL_REALLOC(sm->n, sizeof(int)   * (L+1));
    sm->ralloc = L+1;
    sm->n_rrealloc++;
  }

  sm->nseg    = 0;
  sm->nrow    = 0;
  sm->ncells  = 0;
  sm->last_i  = sm->L+1;
  for (r = 0; r < p7_VNF; r++) 
    sm->last_k[r]  = -1; 
  /* sn[] are initialized for each sparse row by _StartRow() */

  /* The realloc counters are NOT reset. They keep accumulating during
   * the life of the object. 
   */
  for (i = 1; i <= L; i++)	/* n[0] will always be 0, but reinit n[1..L] */
    sm->n[i] = 0;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_sparsemask_Sizeof()
 * Synopsis:  Returns current allocated size of a <P7_SPARSEMASK>, in bytes.
 */
size_t
p7_sparsemask_Sizeof(const P7_SPARSEMASK *sm)
{
  size_t n = sizeof(P7_SPARSEMASK);
  n += sm->ialloc * sizeof(int)   * 2; // <i>; *2 = ia,ib pairs 
  n += sm->ralloc * sizeof(int *);     // <k>                   
  n += sm->ralloc * sizeof(int);       // <n>                   
  n += sm->kalloc * sizeof(int);       // <kmem>                
  return n;
}

/* Function:  p7_sparsemask_MinSizeof()
 * Synopsis:  Returns minimum required size of a <P7_SPARSEMASK>, in bytes.
 */
size_t
p7_sparsemask_MinSizeof(const P7_SPARSEMASK *sm)
{
  size_t n = sizeof(P7_SPARSEMASK);
  n += sm->nseg   * sizeof(int) * 2;  // <i>
  n += (sm->L+1)  * sizeof(int *);    // <k>
  n += (sm->L+1)  * sizeof(int);      // <n>
  n += sm->ncells * sizeof(int);      // <kmem>
  return n;
}


/* Function:  p7_sparsemask_Reuse()
 * Synopsis:  Reuse a <P7_SPARSEMASK>.
 *
 * Purpose:   Essentially the equivalent of <_Destroy()>, but
 *            we will reuse the existing allocations. 
 *
 *            Currently a no-op. <_Reinit()> does all the work of
 *            reusing a <P7_SPARSEMASK>. Someday, though, we may want
 *            to have <_Reuse()> reallocate downwards, if a problem
 *            has temporarily required an unusually large amount of
 *            memory.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_sparsemask_Reuse(P7_SPARSEMASK *sm)
{
  return eslOK;
}

/* Function:  p7_sparsemask_Destroy()
 * Synopsis:  Destroy a <P7_SPARSEMASK>.
 */
void
p7_sparsemask_Destroy(P7_SPARSEMASK *sm)
{
  if (sm) {
    if (sm->i)    free(sm->i);
    if (sm->k)    free(sm->k);
    if (sm->n)    free(sm->n);
    if (sm->kmem) free(sm->kmem);
    free(sm);
  }
}
/*-------------- end, P7_SPARSEMASK -----------------------------*/




/*****************************************************************
 * 2. API for constructing a sparse DP matrix
 *****************************************************************/

/* Function:  p7_sparsemask_AddAll()
 * Synopsis:  Add all cells to the mask.
 *
 * Purpose:   For sparse mask <sm> that has been created or
 *            reinitialized for an <M> by <L> comparison with
 *            posterior probability <0.0>; include every cell.
 *            
 *            We could do this more efficiently since we know we're
 *            going to add exactly ML cells, but since we use AddAll()
 *            only in test code, it is advantageous to implement it
 *            using the standard API.  
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation error. Now <sm> may only
 *            be destroyed; otherwise its contents are undefined.
 * 
 *            <eslEINVAL> on coding errors, failures of contract checks.
 */
int
p7_sparsemask_AddAll(P7_SPARSEMASK *sm)
{
  int  i,k;
  int  status;

  for (i = sm->L; i >= 1; i--)
    {
      if (  (status = p7_sparsemask_StartRow(sm, i))                   != eslOK) return status;
      for (k = sm->M; k >= 1; k--)
	if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) return status;
      if (  (status = p7_sparsemask_FinishRow(sm))                     != eslOK) return status;
    }
  return p7_sparsemask_Finish(sm);
}



/* Function:  p7_sparsemask_StartRow()
 * Synopsis:  Prepare to store sparse cells on a new row i.
 *
 * Purpose:   Here we set up the vector "slots" (generally four of them,
 *            for vectors of four floats; the actual number is
 *            configured in <p7_VNF>) for temporarily storing sorted
 *            lists of k indices for each striped vector segment.
 *            Each slot must allow up to Q entries; so kmem must be
 *            allocated for at least ncells+p7_VNF*Q]. Reallocate
 *            <kmem> (by doubling) if needed. 
 * 
 *            Why do we need these shenanigans? The f/b filter uses
 *            striped vectors, which are not in the M..1 order that
 *            sparse DP needs. Temporary storage in four sorted slots,
 *            followed by their concatenation, is one way to rearrange
 *            efficiently. See note [3] in p7_sparsemx.h.
 *
 *            Remember, the whole <kmem> array is in reverse order during
 *            collection, so slot n[3] is first, n[0] is last; see note [3] in
 *            p7_sparsemx.h.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure. Now the contents of <sm> 
 *            are undefined, except that <_Destroy()> may be safely called 
 *            on it.
 *
 *            <eslEINVAL> on coding errors, failures of contract checks.
 */
int
p7_sparsemask_StartRow(P7_SPARSEMASK *sm, int i)
{
  int r;
  int status;
  
#ifdef p7_DEBUGGING
  if (i < 1 || i > sm->L) ESL_EXCEPTION(eslEINVAL, "i is 1..L: sequence position");
  if (sm->last_i <= i)    ESL_EXCEPTION(eslEINVAL, "rows need to be added in reverse order L..1");
#endif

  /* Make sure kmem has enough memory; if not, double it.
   * Because we know the original allocation was enough to hold
   * the slots, we know that doubling (even if ncells has filled
   * the current kalloc) is sufficient.
   */
  if (sm->ncells + p7_VNF*sm->Q > sm->kalloc)
    {
      int64_t kalloc_req = sm->kalloc * 2;
      ESL_REALLOC(sm->kmem, sizeof(int) * kalloc_req);
      sm->kalloc = kalloc_req;
      sm->n_krealloc++;
    }
  
  for (r = 0; r < p7_VNF; r++)
    {
      sm->s[p7_VNF-r-1] = sm->kmem + sm->ncells + r*sm->Q;
      sm->sn[r]         = 0;
    }
  sm->last_i = i;
  for (r = 0; r < p7_VNF; r++) 
    sm->last_k[r] = sm->M+1;		/* sentinel to be sure that Add() is called in reverse order M..1 */
  return eslOK;
  
 ERROR:
  return status;
}



/* Function:  p7_sparsemask_Add()
 * Synopsis:  Add a vector cell q,r to k list on current row.
 *
 * Purpose:   For a row in progress (we've already called <_StartRow()>),
 *            add a striped vector cell <q,r>. We must call <q>'s
 *            in descending order <Q-1..0>. The <r>'s in a given vector
 *            may be called in any order. 
 *            
 *            This is obviously designed for the production-code case,
 *            where we are constructing the <P7_SPARSEMASK> in the
 *            vectorized f/b filter's linear-memory backwards pass, we
 *            have striped vector coords q,r, and where we have to get
 *            the k indices back into order from their striped
 *            indexing. For testing and debugging code, where we want
 *            to store normal k indices (as opposed to q,r vector
 *            coordinates), the magic is to convert k to simulated
 *            vector coords: that's <q = (k-1)%sm->Q>, <r =
 *            (k-1)/sm->Q>. Though that uglifies calls to this
 *            function from our test code, it's better to do such
 *            contortions so our test code goes through the same API
 *            as production uses.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> on coding errors, failures of contract checks. 
 */
int
p7_sparsemask_Add(P7_SPARSEMASK *sm, int q, int r)
{
  int     k = r*sm->Q+q+1;

  //printf("Adding i=%d q=%d r=%d k=%d M=%d\n", sm->last_i, q, r, k, sm->M);

#ifdef p7_DEBUGGING 
  if (q < 0 || q >= sm->Q)  ESL_EXCEPTION(eslEINVAL, "q is 0..Q-1; striped vector index");
  if (r < 0 || r >= p7_VNF) ESL_EXCEPTION(eslEINVAL, "r is 0..%d-1; segment index in vectors", p7_VNF);
  if (sm->last_k[r] == -1)  ESL_EXCEPTION(eslEINVAL, "need to StartRow() before calling Add()");
  if (sm->last_k[r] <= k)   ESL_EXCEPTION(eslEINVAL, "cells must be added in reverse order M..1");
#endif

  sm->s[r][sm->sn[r]] = k;
  sm->sn[r]++;
  sm->last_k[r] = k;
  return eslOK;
}

/* Function:  p7_sparsemask_FinishRow()
 * Synopsis:  Done adding sparse cells on a current row.
 *
 * Purpose:   This collapses (concatenates) the slots in <kmem[]>,
 *            increments <ncells>, and sets n[i] for the current
 *            row. No more cells can be added until <_StartRow()>
 *            initializes slots for a new row. The kmem[] for
 *            this row is now contiguous, albeit in reverse order.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_sparsemask_FinishRow(P7_SPARSEMASK *sm)
{
  int *p;
  int  r;

  /* s[3] is already where it belongs; so we start at p = kmem + ncells + sn[3] */
  p = sm->kmem + sm->ncells + sm->sn[p7_VNF-1];
  sm->n[sm->last_i] = sm->sn[p7_VNF-1];
  for (r = p7_VNF-2; r >= 0; r--)
    {
      memmove(p, sm->s[r], sizeof(int) * sm->sn[r]);
      p += sm->sn[r];
      sm->n[sm->last_i] += sm->sn[r];
    }
  /* now the slots are invalid; the next StartRow() will reset them */
  sm->ncells += sm->n[sm->last_i];
  for (r = 0; r < p7_VNF; r++) 
    sm->last_k[r]  = -1;	/* that'll suffice to prevent Add() from being called after FinishRow(). */
  return eslOK;
}

/* Function:  p7_sparsemask_Finish()
 * Synopsis:  Done adding cells to the mask.
 *
 * Purpose:   We have finished collecting an entire <P7_SPARSEMASK> for
 *            a given profile/sequence comparison we want to do using
 *            sparse DP routines. <kmem>, <n[i]>, and <ncells> fields
 *            are valid, though <kmem> is reversed.
 *            
 *            Here we reverse <kmem>, then set the rest: k[i] pointers
 *            into <kmem>, the <i[]> coord pairs ia,ib for each
 *            segment, the <nseg> counter for segments, and the <nrow>
 *            counter for the nonzero <n[i]>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. Now <sm> may only be
 *            safely destroyed; its contents are otherwise undefined.
 */
int
p7_sparsemask_Finish(P7_SPARSEMASK *sm)
{
  int i,r;
  int *p;
  int status;

  /* Reverse kmem. */
  esl_vec_IReverse(sm->kmem, sm->kmem, sm->ncells);

  /* Set the k[] pointers; count <nseg> and <nrow> */
  p = sm->kmem;
  sm->nseg = sm->nrow = 0;
  for (i = 1; i <= sm->L; i++)
    if (sm->n[i]) 
      {
	sm->nrow++;
	sm->k[i] = p;
	p       += sm->n[i];
	if (sm->n[i-1] == 0) sm->nseg++;
      } 
    else 
      {
	sm->k[i] = NULL;
      }

  /* Reallocate i[] if needed. */
  if (sm->nseg > sm->ialloc) 
    {
      ESL_REALLOC(sm->i, sizeof(int) * 2 * sm->nseg); /* *2, for ia,ib pairs */
      sm->ialloc = sm->nseg;
      sm->n_irealloc++;
    }
      
  /* Set i[] coord pairs. */
  p = sm->i;
  for (i = 1; i <= sm->L; i++) 
    {
      if (sm->n[i-1] > 0 && sm->n[i] == 0) *p++ = i-1; /* end of a previous segment, ib coord */
      if (sm->n[i] > 0 && sm->n[i-1] == 0) *p++ = i;   /* start a new segment, ia coord */
    }
  if (sm->n[sm->L]) *p = sm->L;	/* last segment ended exactly at L. */

  sm->last_i = -1;
  for (r = 0; r < p7_VNF; r++) 
    sm->last_k[r] = -1;
  return eslOK;

 ERROR:
  return eslEMEM;
}
/*----------------- end, P7_SPARSEMASK API ----------------------*/




/*****************************************************************
 * 3. Debugging tools for P7_SPARSEMASK
 *****************************************************************/

/* Function:  p7_sparsemask_Dump()
 * Synopsis:  Dump contents of a P7_SPARSEMASK
 *
 * Purpose:   Dump the contents of the sparse mask <sm> to 
 *            open stream <ofp> for debugging, diagnostics.
 */
int
p7_sparsemask_Dump(FILE *ofp, P7_SPARSEMASK *sm)
{
  int i,k,z;

  fprintf(ofp, "# sparse mask: M=%d L=%d Q=%d\n", sm->M, sm->L, sm->Q);
  fputs("     ", ofp);  for (k = 1; k <= sm->M; k++) fprintf(ofp, "%3d ", k);  fputs(" n \n", ofp);
  fputs("     ", ofp);  for (k = 1; k <= sm->M; k++) fputs("--- ", ofp);       fputs("---\n", ofp);

  for (i = 1; i <= sm->L; i++)
    {
      fprintf(ofp, "%3d: ", i);
      for (z = 0, k = 1; k <= sm->M; k++)
	{
	  while (z < sm->n[i] && sm->k[i][z] < k)  z++;
	  if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "  X ");
	  else                                     fprintf(ofp, "  . ");
	}
      fprintf(ofp, "%3d\n", sm->n[i]);
    }
  return eslOK;
}
/*-------- end, P7_SPARSEMASK debugging tools -------------------*/



/*****************************************************************
 * 4. P7_SPARSEMX; sparse DP matrix, as specified by a given mask
 *****************************************************************/

/* Function:  p7_sparsemx_Create()
 * Synopsis:  Create a new sparse DP matrix.
 *
 * Purpose:   Create a new sparse DP matrix, defined by the sparse DP
 *            mask <sm>.
 *            
 *            If <sm> is <NULL>, allocate for a reasonable-sized
 *            matrix. Matrix won't be usable until caller calls
 *            <p7_sparsemx_Reinit()>. This style works well in loops
 *            over sequences, where you can create the object once
 *            before any sequences are known, then
 *            <_Reinit()/_Reuse()> it in the loop.
 */
P7_SPARSEMX *
p7_sparsemx_Create(P7_SPARSEMASK *sm)
{
  P7_SPARSEMX *sx            = NULL;
  int64_t      default_ncell = 4096; /* if <sm> is NULL, what <dalloc> should be (stored sparsified cells) */
  int          default_nx    = 256;  /* if <sm> is NULL, what <xalloc> should be (stored rows of specials) */
  int          status;

  ESL_ALLOC(sx, sizeof(P7_SPARSEMX));
  sx->dp   = NULL;
  sx->xmx  = NULL;
  sx->sm   = sm;
  sx->type = p7S_UNSET;

  sx->dalloc = (sm ? sm->ncells          : default_ncell);
  sx->xalloc = (sm ? sm->nrow + sm->nseg : default_nx);

  ESL_ALLOC(sx->dp,  sizeof(float) * p7S_NSCELLS * sx->dalloc);
  ESL_ALLOC(sx->xmx, sizeof(float) * p7S_NXCELLS * sx->xalloc);
  return sx;

 ERROR:
  p7_sparsemx_Destroy(sx);
  return NULL;
}


/* Function:  p7_sparsemx_Reinit()
 * Synopsis:  Reinitialize a sparse DP matrix.
 *
 * Purpose:   Reinitialize an existing sparse matrix <sx> to use the
 *            sparse mask <sm>. Equivalent to a call to
 *            <p7_sparsemx_Create()> except we reuse as much memory as
 *            possible in the preexisting <sx>, to minimize
 *            realloc/free calls.
 *
 *            Caller should have also called <p7_sparsemx_Reuse()> on
 *            <sx>, but at least at present, this is not strictly
 *            necessary.
 *
 * Returns:   <eslOK> on success, and <sx> is now ready for sparse DP
 *            calculations defined by the mask <sm>.
 *
 * Throws:    <eslEMEM> on (re-)allocation failure.
 */
int
p7_sparsemx_Reinit(P7_SPARSEMX *sx, P7_SPARSEMASK *sm)
{
  int64_t dalloc_req = sm->ncells;
  int     xalloc_req = sm->nrow + sm->nseg;
  int     status;

  if (sx->dalloc < dalloc_req) {
    ESL_REALLOC(sx->dp, sizeof(float) * p7S_NSCELLS * dalloc_req);
    sx->dalloc = dalloc_req;
  }
  if (sx->xalloc < xalloc_req) {
    ESL_REALLOC(sx->xmx, sizeof(float) * p7S_NXCELLS * xalloc_req);
    sx->xalloc = xalloc_req;
  }
  sx->sm   = sm;
  sx->type = p7S_UNSET;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_sparsemx_Sizeof()
 * Synopsis:  Returns current allocated size of a P7_SPARSEMX, in bytes.
 *
 * Purpose:   Returns allocated size of <sx>, in bytes.  Does not
 *            include the size of its <P7_SPARSEMASK>.
 */
size_t
p7_sparsemx_Sizeof(const P7_SPARSEMX *sx)
{
  size_t n = sizeof(P7_SPARSEMX);
  n += sizeof(float) * p7S_NSCELLS * sx->dalloc;
  n += sizeof(float) * p7S_NXCELLS * sx->xalloc;
  return n;
}

/* Function:  p7_sparsemx_MinSizeof()
 * Synopsis:  Returns minimal required allocation size of a P7_SPARSEMX, in bytes.
 *
 * Purpose:   Calculate and return the minimum required size, in 
 *            bytes, of a sparse DP matrix, for a profile/sequence
 *            comparison using the sparse mask in <sm>.
 *            
 *            Does not require having an actual DP matrix allocated.
 *            We use this function when planning/profiling memory
 *            allocation strategies.
 */
size_t
p7_sparsemx_MinSizeof(const P7_SPARSEMASK *sm)
{
  size_t n = sizeof(P7_SPARSEMX);
  n += sizeof(float) * p7S_NSCELLS * sm->ncells;             // dp[]
  n += sizeof(float) * p7S_NXCELLS * (sm->nrow + sm->nseg);  // xmx[]; for each seg ia..ib, ia-1..ib has special cells
  return n;
}


/* Function:  p7_sparsemx_GetSpecial()
 * Synopsis:  Look up and return cell value for a special state
 *
 * Purpose:   Returns the value for state <s> on row (sequence position) <i> 
 *            in the sparse matrix <sx>. 
 *            
 *            The caller must know that specials for this row <i> is
 *            stored; <i> must lie in <ia-1..ib> for a stored segment
 *            <g>. If caller requests an invalid (unstored) <i> this
 *            functions will fail immediately with an <esl_fatal()>
 *            coding error.
 *
 *            The <P7_SPARSEMX> is designed for fast sweeps Forward or
 *            Backward, not for random access of its values. Random
 *            access of one value with a <p7_sparsemx_Get*()> requires
 *            a partial sweep, and so is relatively expensive. If you
 *            need to do many random accesses, it is probably better
 *            to implement your own ordered sweep across the values you
 *            need, rather than calling <p7_sparsemx_Get*()> many times.
 *
 * Args:      sx - sparse matrix (of any type except p7S_ENVSCORE)
 *            i  - which row, 1..L, must be ia-1..ib in some stored segment
 *            s  - which special state; example: p7S_C
 *
 * Returns:   the looked-up value.
 */
float
p7_sparsemx_GetSpecial(const P7_SPARSEMX *sx, int i, int s)
{
  float *xc = sx->xmx;
  int g, ia, ib;

  for (g = 0; g < sx->sm->nseg; g++)
    {
      ia = sx->sm->i[2*g];
      ib = sx->sm->i[2*g+1];
      
      if      (i > ib)    { xc += p7S_NXCELLS * (ib-ia+2); continue;     }
      else if (i >= ia-1) { xc += p7S_NXCELLS * (i-ia+1);  return xc[s]; } // normal return
      else    esl_fatal("i=%d was not stored - why are you asking for it?");
    }
  esl_fatal("i=%d was not stored - why are you asking for it?");
  /*NOTREACHED*/
  return -eslINFINITY;
}


/* Function:  p7_sparsemx_Reuse()
 * Synopsis:  Reuse (rather than free) a sparse DP matrix.
 *
 * Purpose:   Reuse the memory allocation in <sx>. This is an
 *            alternative to free'ing the structure entirely 
 *            (using <p7_sparsemx_Destroy()>.) The caller must
 *            reinitialize the matrix with a new mask using
 *            <p7_sparsemx_Reinit()> before using it again.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_sparsemx_Reuse(P7_SPARSEMX *sx)
{
  sx->sm   = NULL;
  sx->type = p7S_UNSET;
  return eslOK;
}

/* Function:  p7_sparsemx_Destroy()
 * Synopsis:  Free a sparse DP matrix.
 */
void
p7_sparsemx_Destroy(P7_SPARSEMX *sx)
{
  if (sx) {
    if (sx->dp)  free(sx->dp);
    if (sx->xmx) free(sx->xmx);
    /* sx->sm is a reference copy. caller remains responsible for it. */
    free(sx);
  }
}
/*----------- end, P7_SPARSEMX implementation -------------------*/



/*****************************************************************
 * 5. API for extracting information from a sparse DP matrix
 *****************************************************************/

int
p7_sparsemx_TracePostprobs(P7_SPARSEMX *sxd, P7_TRACE *tr)
{
  const P7_SPARSEMASK *sm  = sxd->sm;
  const float         *dpc = sxd->dp;   /* ptr that steps through stored main supercells */
  const float         *xc  = sxd->xmx;  /* ptr that steps through stored special rows, including ia-1 segment edges */
  int i  = 0;       /* important to start at i=0, <xc> can start there */
  int k;	    /* index in main MI states, 0,1..M */
  int v;	    /* index that steps through sparse cell list on a row */
  int z;	    /* index that steps through trace */
  int last_ib;	    /* last stored row index - all i>last_ib are emitted by CC with postprob 1.0 */

#ifdef p7_DEBUGGING
  if (sxd->type != p7S_DECODING) ESL_EXCEPTION(eslEINVAL, "<sxd> must be a decoding matrix");
  if (! tr->pp)                  ESL_EXCEPTION(eslEINVAL, "<tr> must be allocated for including posterior probs");
#endif

  /* <last_ib> addresses a boundary case; i>last_ib 
   * must be generated by CC with postprob = 1.0.
   */
  for (last_ib = sm->L; last_ib >=1; last_ib--)
    if (sm->n[last_ib]) break;

  for (z = 0; z < tr->N; z++)
    {
      if (! tr->i[z]) { tr->pp[z] = 0.0; continue; }  // only storing pp's for emissions. Only M,I, or NN/CC/JJ have i=1..L.

      k = tr->k[z];         // will be 0 for s={NJC}; 1..M for s={MID}
      while (i < tr->i[z])  // this happens to also accomodate the case that a trace somehow starts with i>1
	{                   // i is < L here, so n[i+1] is ok below
	  dpc += p7S_NSCELLS * sm->n[i]; 
	  if (sm->n[i] || sm->n[i+1]) xc += p7S_NXCELLS; // special rows also store ia-1, hence the move on sm->n[i+1] too
	  i++;
	}
      /*  ok, now: 
       *   dpc is either on first supercell of row i (if i is stored, n[i]>0); 
       *       or it's waiting on next stored row;
       *       or it's just off area (and can't be deferenced) if there are no more stored rows.
       *       (If it's off area, we'll get v=0 == n[i] in logic below, and avoid dereference.)
       *   xc  is on special rows of row i (if i specials stored, including ia-1 start-segment edge cases);
       *       or it's waiting on the next stored special row;
       *       or it's just off area and i>last_ib test will avoid deferencing it.
       *   
       *   in unstored intersegment regions ib+1..ia-2, the postprob of the downstream {CJ}(ia-1)
       *   is propagated up. The logic below exploits the fact that xc is waiting on the next
       *   stored row, the downstream ia-1, when i is in unstored intersegment region.
       */

      /* Look for sparse cell i,k in the list. If not found: v==sm->n[i], and tests below will catch this */
      if (k) { v = 0; while (v < sm->n[i] && sm->k[i][v] < k) v++; }
      
      /* and finally, in logic of neutron-star-like-density: */
      switch (tr->st[z]) {
      case p7T_ML: case p7T_MG: tr->pp[z] = ((v < sm->n[i] && sm->k[i][v] == k) ? *(dpc + p7S_NSCELLS * v + p7S_ML) + *(dpc + p7S_NSCELLS * v + p7S_MG) : 0.0); break;
      case p7T_IL: case p7T_IG: tr->pp[z] = ((v < sm->n[i] && sm->k[i][v] == k) ? *(dpc + p7S_NSCELLS * v + p7S_IL) + *(dpc + p7S_NSCELLS * v + p7S_IG) : 0.0); break;
      case p7T_N:   tr->pp[z] = ((i <= last_ib) ? xc[p7S_N]  : 0.0);  break;
      case p7T_J:   tr->pp[z] = ((i <= last_ib) ? xc[p7S_JJ] : 0.0);  break;
      case p7T_C:   tr->pp[z] = ((i <= last_ib) ? xc[p7S_CC] : 1.0);  break;
      default:      ESL_EXCEPTION(eslEINCONCEIVABLE, "sorry. after all that, that's rather embarrassing.");
      }
    }
  return eslOK;
}


int
p7_sparsemx_ExpectedDomains(P7_SPARSEMX *sxd, int iae, int ibe, float *ret_ndom_expected)
{
  const P7_SPARSEMASK *sm = sxd->sm;
  float  expB = 0.;
  float  expE = 0.;
  int    i;
  float *xc   = sxd->xmx;

  for (i = 0; i <= ibe; i++)
    {
      if (! (sm->n[i] || (i < sm->L && sm->n[i+1]))) continue; // skip unstored special rows
      if (i >= iae-1 && i < ibe) expB += xc[p7S_B];            // B(i-1) is pp of starting on i 
      if (i >= iae)              expE += xc[p7S_E];
      xc += p7S_NXCELLS;
    }

  /* Over the whole sequence, expB=expE, within roundoff error
   * accumulation.  But on an envelope iae..ibe of a full-sequence
   * Decoding matrix, expE is not necessarily equal to expB; one end
   * may be better determined and entirely within the envelope, for
   * example. So which do we use as the expected # of domains in the
   * envelope? Since the main use of this function is to decide
   * whether an envelope is simply a single domain or not -- whether
   * we can directly take a per-domain score from the whole-seq
   * Forward matrix -- we report the max(B,E), to be conservative
   * about flagging when it looks like a second (or more) domain
   * is impinging on this envelope.
   */
  *ret_ndom_expected = ESL_MAX(expE, expB);
  return eslOK;
}


/* Function:  p7_sparsemx_ApproxEnvScore()
 * Synopsis:  Envelope score by fast Forward endpoint approximation method.
 *
 * Purpose:   Implements envelope scoring by the fast Forward endpoint
 *            approximation method. Returns the envelope score of
 *            envelope <iae..ibe> in the target sequence for which
 *            we've calculated a Forward matrix <sxf>, comparing it to
 *            profile <gm>. The "envelope score" is defined as the
 *            Forward score for the entire sequence <1..L>, subject to
 *            the condition that there is only a single domain that lies
 *            somewhere within <iae..ibe> on the target sequence.
 *
 *            <*ret_envsc> is a raw Forward score in nats.  The caller
 *            still needs to include the null model(s) contribution to
 *            this score and convert it to bits, before calculating
 *            P-values or reporting it in output.
 *            
 *            Caller has determined that the envelope <iae..ibe>
 *            contains only a single domain (negligible probability
 *            mass in any more than that); otherwise, this technique
 *            does not work. 
 *            
 *            The model is assumed to have <tCC=tJJ=tNN> equal
 *            transition scores, or this method does not work.
 *            
 *            For a single-domain envelope <iae..ibe>, with equal
 *            CC,JJ,TT transition scores, the score <delta>
 *            contributed by the <iae..ibe> region is
 *            
 *            $\delta = C_j 
 *                      + \log \[ 1 - e^{C_{i-1} + L_e \tau_{CC} - C_j} \]
 *                      - \log \[ e^{N_{i-1}} + e^{C_{i-1}} \]$
 *            
 *            where $L_e$ is the envelope length <ibe-iae+1>. We then
 *            add the <S->N...> contribution for <iae-1> flanking
 *            residues and <...C->T> contribution for <L-ibe> flanking
 *            residues to convert <delta> to a raw Forward score for
 *            the whole sequence, constrained to a single domain in
 *            the envelope.
 *
 * Args:      gm        - profile used to create the Forward matrix
 *            sxf       - Forward matrix calculated for gm x seq comparison
 *            iae       - left endpoint of envelope, 1..L on target seq
 *            ibe       - right endpoint of envelope, 1..L on target seq
 *            ret_envsc - RETURN: raw envelope Forward score, in nats
 *
 * Returns:   <eslOK> on success, and <*ret_envsc> contains the raw 
 *            envelope Forward score in nats.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      SRE:J10/43-45.
 */
int
p7_sparsemx_ApproxEnvScore(P7_PROFILE *gm, P7_SPARSEMX *sxf, int iae, int ibe, float *ret_envsc)
{
  const float *xc = sxf->xmx;
  float Ci1, Ni1, Cj;
  float deltaC;
  float Le;
  int   g, ia, ib, last_ib;
  float envsc;

  /* Accessing the {CN}(iae-1) and C(ibe) entries is complicated because of sparse storage; forgive. */
  /* Find first stored special row i>=iae-1; get C(i-1), N(i-1) from there; reset iae=i */
  for (g = 0; g < sxf->sm->nseg; g++)
    {
      ia = sxf->sm->i[g*2];
      ib = sxf->sm->i[g*2+1];
      
      if      (iae < ia-1) {                             iae = ia;  Ci1 = xc[p7S_C]; Ni1 = xc[p7S_N]; break; }
      else if (iae <= ib)  { xc += (iae-ia)*p7S_NXCELLS; ia  = iae; Ci1 = xc[p7S_C]; Ni1 = xc[p7S_N]; break; }
      else    xc += (ib-ia+2) * p7S_NXCELLS;
    }
  /* now xc is on row ia-1, in segment g, which ends at ib. check if ib is in that segment. */
  
  /* Find last stored row j<=ibe, get C(j) from there; reset ibe=j */
  if (ib < iae)  { *ret_envsc = -eslINFINITY; return eslOK; }
  if (ibe <= ib) {
    xc += (ibe - ia + 1) * p7S_NXCELLS;
    Cj = xc[p7S_C];
  } else {
    xc     += (ib-ia+1) * p7S_NXCELLS; // now xc is on last supercell (ibe) in previous segment
    last_ib = ib;
    for (g = g+1; g < sxf->sm->nseg; g++)
      {
	ia = sxf->sm->i[g*2];
	ib = sxf->sm->i[g*2+1];
	
	if      (ibe < ia-1) { ibe = last_ib;                 Cj = xc[p7S_C]; break; }
	else if (ibe <= ib)  { xc += (ibe-ia+1)*p7S_NXCELLS;  Cj = xc[p7S_C]; break; }
	else  {
	  xc += (ib-ia+2)*p7S_NXCELLS;
	  last_ib = ib;
	}
      }
    if (g == sxf->sm->nseg) { ibe = last_ib; Cj = xc[p7S_C]; }
  }
  /* now iae,ibe may have been moved, to correspond to outermost stored rows in the envelope */


  /* First part of envsc is Cj + log(1-exp(-deltaC)), and the log term needs
   * special numerical handling; using lim x->0 1-e^-x = x for small deltaC,
   * lim x->0 log (1-x) = -x for large deltaC
   */
  envsc  = Cj;			/* first term. */
  Le     = ibe - iae + 1;
  deltaC = Cj - Ci1 - Le * gm->xsc[p7P_C][p7P_LOOP];

  if      (deltaC  < 0) ESL_EXCEPTION(eslEINCONCEIVABLE, "no, something's wrong, this intermediate term is >= 0 by construction");
  if      (deltaC == 0)                  envsc = -eslINFINITY;
  else if (deltaC < eslSMALLX1)          envsc += logf(deltaC);        // logf(deltaC) -> -eslINFINITY, for small deltaC ->0
  else if (exp(-1.*deltaC) < eslSMALLX1) envsc -= expf(-1.*deltaC);    // expf(-deltaC) -> 0, for large deltaC -> inf
  else                                   envsc += logf(1.-exp(-1.*deltaC));

  /* third term is a standard log-sum-exp-2, may as well use our fast approx */
  envsc -= p7_FLogsum(Ci1, Ni1);
  
  /* left flank: envsc already includes ...N->B->...; add S->N..N flank. S->N is 1.0.  */
  envsc += gm->xsc[p7P_N][p7P_LOOP] * (iae-1);
  
  /* right flank: envsc includes E->C; add C..C->T flank */
  envsc += gm->xsc[p7P_C][p7P_LOOP] * (sxf->sm->L - ibe);
  envsc += gm->xsc[p7P_C][p7P_MOVE];

  *ret_envsc = envsc;
  return eslOK;
}



/*****************************************************************
 * 6. Debugging tools for P7_SPARSEMX
 *****************************************************************/ 

char *
p7_sparsemx_DecodeState(int type)
{
  switch (type) {
  case p7S_ML: return "ML";
  case p7S_MG: return "MG";
  case p7S_IL: return "IL";
  case p7S_IG: return "IG";
  case p7S_DL: return "DL";
  case p7S_DG: return "DG";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_REFMX main state code %d\n", type);
  return NULL;
}

char *
p7_sparsemx_DecodeSpecial(int type)
{
  switch (type) {
  case p7S_E:  return "E";
  case p7S_N:  return "N";
  case p7S_J:  return "J";
  case p7S_B:  return "B";
  case p7S_L:  return "L";
  case p7S_G:  return "G";
  case p7S_C:  return "C";
  case p7S_JJ: return "JJ";
  case p7S_CC: return "CC";
  default:     break;
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_SPARSEMX special state type code %d\n", type);
  return NULL;
}

int
p7_sparsemx_Dump(FILE *ofp, P7_SPARSEMX *sx)
{
  return p7_sparsemx_DumpWindow(ofp, sx, 0, sx->sm->L, 0, sx->sm->M);
}

int
p7_sparsemx_DumpWindow(FILE *ofp, P7_SPARSEMX *sx, int i1, int i2, int k1, int k2)
{
  P7_SPARSEMASK *sm  = sx->sm;
  float         *dpc = sx->dp;
  float         *xc  = sx->xmx;
  int width          = 9;
  int precision      = 4;
  int i,k,x,z;

   /* Header */
  fprintf(ofp, "       ");
  for (k = k1; k <= k2;         k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0;  x < p7S_NXCELLS; x++) fprintf(ofp, "%*s ", width, p7_sparsemx_DecodeSpecial(x));
  fprintf(ofp, "\n");

  fprintf(ofp, "       ");
  for (k = k1; k <= k2;        k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < p7S_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");

  /* Skipping ahead in matrix, over rows we're not dumping: */
  for (i = 1; i < i1; i++) 
    if (sm->n[i]) {
      if (sm->n[i-1] == 0) xc += p7R_NXCELLS; /* skip an extra chunk of specials (ia-1) before each segment start on ia */
      dpc += sm->n[i] * p7R_NSCELLS;          /* skip over rows we're not dumping */
      xc  += p7R_NXCELLS;		      /* skip specials on sparsified rows */
    }

  for (i = i1; i <= i2; i++)
    {
      /* If current row has no cells...  */
      if (sm->n[i] == 0) {
	if (i > 1     && sm->n[i-1] > 0) fputs("...\n\n", ofp);  /* ... if last row did, then we ended a segment. Print an ellipsis. */
	if (i < sm->L && sm->n[i+1] > 0)                          /* if next row does, then we're about to start a segment; print specials on an ia-1 row   */
	  {
	    fprintf(ofp, "%3d -- ", i);
	    for (k = k1; k <= k2;         k++) fprintf(ofp, "%*s ", width, ".....");
	    for (x = 0;  x < p7S_NXCELLS; x++) fprintf(ofp, "%*.*f ", width, precision, xc[x]);
	    fputs("\n\n", ofp);
	    xc += p7R_NXCELLS;	  
	  }
	continue;			
      }

      fprintf(ofp, "%3d ML ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
	if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_ML));
	else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      for (x = 0; x < p7S_NXCELLS; x++) fprintf(ofp, "%*.*f ", width, precision, xc[x]);
      fputc('\n', ofp);

      fprintf(ofp, "%3d IL ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k)  z++;
	if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_IL));
	else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d DL ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
	if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_DL));
	else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d MG ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
	if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_MG));
	else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d IG ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
	if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_IG));
	else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d DG ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
	if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_DG));
	else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputs("\n\n", ofp);

      dpc += sm->n[i] * p7R_NSCELLS;
      xc  += p7R_NXCELLS;
    }
  return eslOK;
}

int
p7_sparsemx_Copy2Reference(P7_SPARSEMX *sx, P7_REFMX *rx)
{
  P7_SPARSEMASK *sm = sx->sm;
  int             M = sm->M;
  int             L = sm->L;
  int             W = (M+1)*p7R_NSCELLS + p7R_NXCELLS; /* total width of a reference DP row, in cells */
  float           vimp;
  int            *imem;
  float          *xc  = sx->xmx;
  float          *dpc = sx->dp;
  int g,i,ia,ib,k,x,z;

  /* clearest/most robust thing to do is first -inf the entire <rx>;
   * we're not worried too much about time here, no need to be fancy.
   */
  vimp = (sx->type == p7S_DECODING ? 0.0f : -eslINFINITY);
  for (i = 0; i <= L; i++)
    for (x = 0; x < W; x++) 
      rx->dp[i][x] = vimp;

  /* These must be copied first: <rx> needs to know L,M dimensions for P7R_XMX,P7R_MX macros to work */
  rx->type = sx->type;
  rx->L    = sm->L;
  rx->M    = sm->M;

  /* now traverse the sparse <sx> */
  imem = sm->i;
  for (g = 0; g < sm->nseg; g++)
    {
      ia = *imem++;
      ib = *imem++;

      /* specials on row ia-1 before each segment are stored. */
      /* don't rely on the order being the same in ref, sparse */
      P7R_XMX(rx, ia-1, p7R_E)  = *xc++;
      P7R_XMX(rx, ia-1, p7R_N)  = *xc++;
      P7R_XMX(rx, ia-1, p7R_J)  = *xc++;
      P7R_XMX(rx, ia-1, p7R_B)  = *xc++;
      P7R_XMX(rx, ia-1, p7R_L)  = *xc++;
      P7R_XMX(rx, ia-1, p7R_G)  = *xc++;
      P7R_XMX(rx, ia-1, p7R_C)  = *xc++;
      P7R_XMX(rx, ia-1, p7R_JJ) = *xc++;
      P7R_XMX(rx, ia-1, p7R_CC) = *xc++;

      for (i = ia; i <= ib; i++)
	{
	  for (z = 0; z < sm->n[i]; z++)
	    {
	      k = sm->k[i][z];
	      P7R_MX(rx, i, k, p7R_ML) = *dpc++;
	      P7R_MX(rx, i, k, p7R_MG) = *dpc++;
	      P7R_MX(rx, i, k, p7R_IL) = *dpc++;
	      P7R_MX(rx, i, k, p7R_IG) = *dpc++;
	      P7R_MX(rx, i, k, p7R_DL) = *dpc++;
	      P7R_MX(rx, i, k, p7R_DG) = *dpc++;
	    }
	  P7R_XMX(rx, i, p7R_E)  = *xc++;
	  P7R_XMX(rx, i, p7R_N)  = *xc++;
	  P7R_XMX(rx, i, p7R_J)  = *xc++;
	  P7R_XMX(rx, i, p7R_B)  = *xc++;
	  P7R_XMX(rx, i, p7R_L)  = *xc++;
	  P7R_XMX(rx, i, p7R_G)  = *xc++;
	  P7R_XMX(rx, i, p7R_C)  = *xc++;
	  P7R_XMX(rx, i, p7R_JJ) = *xc++;
	  P7R_XMX(rx, i, p7R_CC) = *xc++;
	}
    }

  return eslOK;
}


/* Function:  p7_sparsemx_CompareReference()
 * Synopsis:  Test sparse DP matrix for tolerable equality to reference
 *
 * Purpose:   Compare all the values in sparse DP matrix <sx> to the
 *            corresponding values in reference DP matrix <rx> for
 *            equality within the absolute epsilon <tol>, using
 *            <esl_FCompareAbs()> calls. Return <eslOK> if comparison
 *            succeeds; return <eslFAIL> otherwise.
 *            
 *            In generaly, this is only going to succeed if <sx> was
 *            calculated as a 'full' sparse matrix, with all cells
 *            marked, as in <p7_sparsemask_AddAll()>. For truly sparse
 *            calculations, you will more likely want to use
 *            <p7_sparsemx_CompareReferenceAsBound()>, which treats
 *            the reference calculation's matrix values as an upper
 *            bound on the corresponding sparse cell values.
 *            
 *            The comparison tests two different traversal mechanisms
 *            for sparse main cells.  Should be completely redundant,
 *            unless we've corrupted the structure.  The routine is
 *            intended for debugging, not speed.
 *            
 *            As also noted in <p7_refmx_Compare()>, absolute
 *            difference comparison is preferred over relative
 *            differences. Numerical error accumulation in DP scales
 *            more with the number of terms than their magnitude. DP
 *            cells with values close to zero (and hence small
 *            absolute differences) may reasonably have large relative
 *            differences.
 *            
 *            If <p7_DEBUGGING> compile-time flag is up, indicating
 *            that we're working in development code (as opposed to
 *            production code or a distro), we <abort()> immediately
 *            on any element comparison failure, to facilitate 
 *            finding errant elements during debugging.
 *            
 *            We assume that <p7S_NSCELLS> and <p7S_NXCELLS> main and
 *            special elements are in the same order and number in
 *            <p7R_NSCELLS> and <p7R_NXCELLS>.
 *
 * Args:      sx  - sparse DP matrix
 *            rx  - reference DP matrix
 *            tol - absolute floating point comparison tolerance
 *                 
 * Returns:   <eslOK> on success.
 *            <eslFAIL> if comparison fails.
 */
int
p7_sparsemx_CompareReference(const P7_SPARSEMX *sx, const P7_REFMX *rx, float tol)
{
  const P7_SPARSEMASK *sm  = sx->sm;
  int            killmenow = FALSE;
  const float   *dpc, *dpc2;
  const float   *xc,  *xc2;
  int            g,i,s,z;
  int            ia,ib;
#ifdef p7_DEBUGGING
  killmenow = TRUE;
#endif
  
  if (sx->type != rx->type) { if (killmenow) abort(); return eslFAIL; }
  if (sm->M    != rx->M)    { if (killmenow) abort(); return eslFAIL; }
  if (sm->L    != rx->L)    { if (killmenow) abort(); return eslFAIL; }

  /* First traversal way: sm->dp[] is just one big array */
  dpc = sx->dp;
  xc  = sx->xmx;
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1])         /* ia-1 specials at a segment start */
	{
	  for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	    if (esl_FCompareAbs(*xc, P7R_XMX(rx,i-1,s), tol) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	}

      for (z = 0; z < sm->n[i]; z++)       /* sparse cells */
	{
	  for (s = 0; s < p7S_NSCELLS; dpc++, s++) 
	    if (esl_FCompareAbs(*dpc, P7R_MX(rx,i,sm->k[i][z],s), tol) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	}
  
      if (sm->n[i])       /* specials */
	{
	  for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	    if (esl_FCompareAbs(*xc, P7R_XMX(rx,i,s), tol) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	}
    }


  /* Second way: "segments" */
  dpc2 = sx->dp;
  xc2  = sx->xmx;
  for (g = 0; g < sm->nseg; g++)
    {
      ia = sm->i[g*2];
      ib = sm->i[g*2+1];

      for (s = 0; s < p7S_NXCELLS; xc2++, s++)        /* ia-1 specials at segment start */
	if (esl_FCompareAbs(*xc2, P7R_XMX(rx,ia-1,s), tol) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }

      for (i = ia; i <= ib; i++) 
	{
	  for (z = 0; z < sm->n[i]; z++)      	  /* sparse main cells */
	    {
	      for (s = 0; s < p7S_NSCELLS; dpc2++, s++) 
		if (esl_FCompareAbs(*dpc2, P7R_MX(rx,i,sm->k[i][z],s), tol) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	    }
	  for (s = 0; s < p7S_NXCELLS; xc2++, s++)  	  /* specials */
	    if (esl_FCompareAbs(*xc2, P7R_XMX(rx,i,s), tol) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	}
    }
  
  /* Both ways must reach the same end */
  if (dpc != dpc2)  { if (killmenow) abort(); return eslFAIL; }
  if (xc  != xc2)   { if (killmenow) abort(); return eslFAIL; }
  
  return eslOK;
}


/* Function:  p7_sparsemx_CompareReferenceAsBound()
 * Synopsis:  Test sparse DP matrix cell values are upper-bounded by reference
 *
 * Purpose:   Check that all values in sparse matrix <sx> are bounded
 *            (less than or equal to) the corresponding values in reference
 *            matrix <rx>, allowing for a small amount <tol> of floating-point
 *            roundoff error accumulation. That is, <v_s <= v_r+tol> for all
 *            matrix cell values <v>. Return <eslOK> if comparison succeeds;
 *            return <eslFAIL> otherwise.
 *            
 *            This is a variant of <p7_sparsemx_CompareReference()>;
 *            see additional documentation there, with the exception
 *            that this routine only uses one traversal mechanism in
 *            the sparse matrix.
 *
 * Args:      sx  - sparse DP matrix
 *            rx  - reference DP matrix
 *            tol - absolute floating point comparison tolerance
 *                 
 * Returns:   <eslOK> on success.
 *            <eslFAIL> if comparison fails.
 */
int
p7_sparsemx_CompareReferenceAsBound(const P7_SPARSEMX *sx, const P7_REFMX *rx, float tol)
{
  const P7_SPARSEMASK *sm  = sx->sm;
  int            killmenow = FALSE;
  const float   *dpc       = sx->dp;
  const float   *xc        = sx->xmx;
  int            i,s,z;
#ifdef p7_DEBUGGING
  killmenow = TRUE;
#endif

  if (sx->type != rx->type) { if (killmenow) abort(); return eslFAIL; }
  if (sm->M    != rx->M)    { if (killmenow) abort(); return eslFAIL; }
  if (sm->L    != rx->L)    { if (killmenow) abort(); return eslFAIL; }
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {    /* ia-1 specials at a segment start */
	for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	  if (*xc > P7R_XMX(rx,i-1,s)+tol)           { if (killmenow) abort(); return eslFAIL; }
      }
      for (z = 0; z < sm->n[i]; z++)    /* sparse cells */
	for (s = 0; s < p7S_NSCELLS; dpc++, s++) 
	  if (*dpc > P7R_MX(rx,i,sm->k[i][z],s)+tol) { if (killmenow) abort(); return eslFAIL; }
  
      if (sm->n[i]) {    /* specials */
	for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	  if (*xc > P7R_XMX(rx,i,s)+tol)             { if (killmenow) abort(); return eslFAIL;  }
      }
    }
  return eslOK;
}
/*------------ end, P7_SPARSEMX debugging tools -----------------*/


/*****************************************************************
 * 7. Validation of a P7_SPARSEMX
 *****************************************************************/
/* also a debugging tool, but in its own section because it's
 * fiddly and complicated
 */

static int
validate_dimensions(P7_SPARSEMX *sx, char *errbuf)
{
  P7_SPARSEMASK *sm     = sx->sm;
  int            g      = 0;
  int            r      = 0;
  int            ncells = 0;
  int            i;

  if ( sm->M <= 0)               ESL_FAIL(eslFAIL, errbuf, "nonpositive M");
  if ( sm->L <= 0)               ESL_FAIL(eslFAIL, errbuf, "nonpositive L");
  if ( sm->Q <  p7O_NQF(sm->M))  ESL_FAIL(eslFAIL, errbuf, "insufficient Q");

  for (r=0, g=0, i = 1; i <= sm->L; i++) {
    if (sm->n[i] && !sm->n[i-1]) g++; /* segment count */
    if (sm->n[i])                r++; /* sparse row count */
    ncells += sm->n[i];
  }
  if (g      != sm->nseg)       ESL_FAIL(eslFAIL, errbuf, "nseg is wrong");
  if (r      != sm->nrow)       ESL_FAIL(eslFAIL, errbuf, "nrow is wrong");
  if (ncells != sm->ncells)     ESL_FAIL(eslFAIL, errbuf, "ncells is wrong");

  if (sm->L+1    > sm->ralloc)  ESL_FAIL(eslFAIL, errbuf, "k[] row allocation too small");
  if (sm->ncells > sm->kalloc)  ESL_FAIL(eslFAIL, errbuf, "kmem[] cell allocation too small");
  if (sm->nseg   > sm->ialloc)  ESL_FAIL(eslFAIL, errbuf, "i[] segment allocation too small");
  return eslOK;
}


static int
validate_no_nan(P7_SPARSEMX *sx, char *errbuf)
{
  P7_SPARSEMASK *sm  = sx->sm;
  float         *dpc = sx->dp;
  float         *xc  = sx->xmx;
  int            i,k,z,s;

  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1])          /* ia-1 specials at a segment start */
	{
	  for (s = 0; s < p7S_NXCELLS; s++) {
	    if (isnan(*xc)) ESL_FAIL(eslFAIL, errbuf, "nan at i=%d, %s", i, p7_sparsemx_DecodeSpecial(s));
	    xc++;
	  }
	}
      for (z = 0; z < sm->n[i]; z++)       /* sparse main cells */
	{
	  k = sm->k[i][z];
	  for (s = 0; s < p7S_NSCELLS; s++) {
	    if (isnan(*dpc)) ESL_FAIL(eslFAIL, errbuf, "nan at i=%d, k=%d, %s", i, k, p7_sparsemx_DecodeState(s));
	    dpc++;
	  }
	}

      if (sm->n[i])                       /* specials on sparse row */
	{
	  for (s = 0; s < p7S_NXCELLS; s++) {
	    if (isnan(*xc)) ESL_FAIL(eslFAIL, errbuf, "nan at i=%d, %s", i, p7_sparsemx_DecodeSpecial(s));
	    xc++;
	  }
	}
    }
  return eslOK;
}

static int
validate_fwdvit(P7_SPARSEMX *sx, char *errbuf)
{
  P7_SPARSEMASK *sm  = sx->sm;
  float         *dpc = sx->dp;
  float         *xc  = sx->xmx;
  int            i,z;

  /* Check special cases prohibited in the first ia-1 presegment specials: */
  if ( xc[p7S_J] !=  -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "first J not -inf");
  if ( xc[p7S_C] !=  -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "first C not -inf");
  
  /* Sweep, checking for (the most easily spotchecked) prohibited values (must be -inf) */
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {       /* ia-1 specials at a segment start */
	if ( xc[p7S_E]  != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "E seg start for ia=%d not -inf", i);
	if ( xc[p7S_JJ] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "JJ seg start for ia=%d not -inf", i);
	if ( xc[p7S_CC] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "CC seg start for ia=%d not -inf", i);
	xc += p7S_NXCELLS;
      }
      for (z = 0; z < sm->n[i]; z++)       /* sparse main cells */
	{
	  /* if k-1 supercell doesn't exist, can't reach D's */
	  if ((z == 0 || sm->k[i][z] != sm->k[i][z-1]+1) && dpc[p7S_DL] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "first DL on i=%d not -inf", i);
	  if ((z == 0 || sm->k[i][z] != sm->k[i][z-1]+1) && dpc[p7S_DG] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "first DG on i=%d not -inf", i);
	  if (   sm->k[i][z] == sm->M                    && dpc[p7S_IL] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "IL on i=%d,k=M not -inf", i);
	  if (   sm->k[i][z] == sm->M                    && dpc[p7S_IG] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "IG on i=%d,k=M not -inf", i);
	  dpc += p7S_NSCELLS;
	  /* there are other conditions where I(i,k) values must be zero but this is more tedious to check */
	}
      if (sm->n[i]) {                     
	if ( xc[p7S_JJ] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "JJ at i=%d not -inf", i);
	if ( xc[p7S_CC] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "CC at i=%d not -inf", i);
	xc += p7S_NXCELLS;
      }
    }
  return eslOK;
}

static int
validate_backward(P7_SPARSEMX *sx, char *errbuf)
{
  P7_SPARSEMASK *sm     = sx->sm;
  float         *dpc    = sx->dp  + (sm->ncells-1)*p7S_NSCELLS;		// last supercell in dp  
  float         *xc     = sx->xmx + (sm->nrow + sm->nseg - 1)*p7S_NXCELLS; // last supercell in xmx 
  int            last_n = 0;
  int            i,z;

  /* Backward sweep; many of our prohibits are on ib segment-end rows */
  /* first: some special cases on absolute final stored row ib */
  if (xc[p7S_N] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "N on last row not 0");
  if (xc[p7S_J] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "N on last row not 0");
  /* sweep: */
  for (i = sm->L; i >= 1; i--)
    {
      if (sm->n[i]) {                   /* specials on stored row i */
	if (               xc[p7S_JJ] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "JJ on row i=%d not -inf", i);
	if (               xc[p7S_CC] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "CC on row i=%d not -inf", i);
	if (last_n == 0 && xc[p7S_B]  != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "B on end-seg row ib=%d not -inf", i);
	if (last_n == 0 && xc[p7S_L]  != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "L on end-seg row ib=%d not -inf", i);
	if (last_n == 0 && xc[p7S_G]  != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "G on end-seg row ib=%d not -inf", i);
	xc -= p7S_NXCELLS;
      }
      for (z = sm->n[i]-1; z >= 0; z--) /* sparse main cells */
	{
	  if (sm->k[i][z] == sm->M && dpc[p7S_IL] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "IL on i=%d,k=M not -inf", i);
	  if (sm->k[i][z] == sm->M && dpc[p7S_IL] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "IL on i=%d,k=M not -inf", i);
	  if (     last_n == 0     && dpc[p7S_IL] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "IL on end-segment row ib=%d,k=%d not -inf", i, sm->k[i][z]);
	  if (     last_n == 0     && dpc[p7S_IG] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "IG on end-segment row ib=%d,k=%d not -inf", i, sm->k[i][z]);
	  dpc -= p7S_NSCELLS;
	}
      if (sm->n[i] && sm->n[i-1] == 0) xc -= p7S_NXCELLS; /* specials on ia-1 row before a start-segment */
      last_n = sm->n[i];
    }
  return eslOK;
}

static int
is_prob(float val, float tol)
{
  if (val < 0.0-tol || val > 1.0+tol) return FALSE; 
  return TRUE;
}

static int
validate_decoding(P7_SPARSEMX *sx, char *errbuf)
{
  P7_SPARSEMASK *sm  = sx->sm;
  float         *dpc = sx->dp;
  float         *xc  = sx->xmx;
  int            i,z,s;
  int            last_n;
  float          tol = 0.001;

  /* Check special cases prohibited in the first ia-1 presegment specials: */
  if ( esl_FCompareAbs(xc[p7S_N],  1.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first N not 1");
  if ( esl_FCompareAbs(xc[p7S_J],  0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first J not 0");
  if ( esl_FCompareAbs(xc[p7S_C],  0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first C not 0");
  if ( esl_FCompareAbs(xc[p7S_JJ], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first JJ not 0");
  if ( esl_FCompareAbs(xc[p7S_CC], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first CC not 0");

  /* Sweep, checking for (the most easily spotchecked) prohibited values (must be 0.0's) */
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {       /* ia-1 specials at a segment start */
	if ( esl_FCompareAbs(xc[p7S_E], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "E seg start for ia=%d not 0", i);
	xc += p7S_NXCELLS;
      }
      for (z = 0; z < sm->n[i]; z++)       /* sparse main cells */
	{
	  /* if k-1 supercell doesn't exist, can't reach D's */
	  if ((z == 0 || sm->k[i][z] != sm->k[i][z-1]+1) && esl_FCompareAbs(dpc[p7S_DL], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first DL on i=%d not 0", i);
	  if ((z == 0 || sm->k[i][z] != sm->k[i][z-1]+1) && esl_FCompareAbs(dpc[p7S_DG], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first DG on i=%d not 0", i);
	  if (   sm->k[i][z] == sm->M                    && esl_FCompareAbs(dpc[p7S_IL], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "IL on i=%d,M not 0", i);
	  if (   sm->k[i][z] == sm->M                    && esl_FCompareAbs(dpc[p7S_IG], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "IG on i=%d,M not 0", i);
	  dpc += p7S_NSCELLS;
	  /* there are other conditions where I(i,k) values must be zero but this is more tedious to check */
	}
      if (sm->n[i]) xc += p7S_NXCELLS;
    }

  /* Backwards sweep, looking only at ib end rows. */
  dpc    = sx->dp  + (sm->ncells-1)*p7S_NSCELLS;		 // last supercell in dp  
  xc     = sx->xmx + (sm->nrow + sm->nseg - 1)*p7S_NXCELLS; // last supercell in xmx 
  last_n = 0;
  /* special cases on absolute final stored row ib: */
  if (esl_FCompareAbs(xc[p7S_N], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "N on last row not 0");
  if (esl_FCompareAbs(xc[p7S_J], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "J on last row not 0");
  if (esl_FCompareAbs(xc[p7S_B], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "B on last row not 0");
  if (esl_FCompareAbs(xc[p7S_L], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "L on last row not 0");
  if (esl_FCompareAbs(xc[p7S_G], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "G on last row not 0");
  /* sweep: */
  for (i = sm->L; i >= 1; i--)
    {
      if (sm->n[i]) xc -= p7S_NXCELLS; /* specials on stored row i */

      for (z = sm->n[i]-1; z >= 0; z--)
	{ // last_n == 0 checks if we're on an end-segment row ib
	  if (last_n == 0 && esl_FCompareAbs(dpc[p7S_IL], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "IL on end-seg row ib=%d not 0", i);
	  if (last_n == 0 && esl_FCompareAbs(dpc[p7S_IG], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "IG on end-seg row ib=%d not 0", i);
	  dpc -= p7S_NSCELLS;
	}

      if (sm->n[i] && sm->n[i-1] == 0) xc -= p7S_NXCELLS; /* specials on ia-1 row before a start-segment */
      last_n = sm->n[i];
    }

  /* Sweep again; check all values are probabilities, 0<=x<=1, allowing a bit of numerical slop. */
  dpc = sx->dp;
  xc  = sx->xmx;
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {       /* ia-1 specials at a segment start */
	for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	  if (! is_prob(*xc, tol)) ESL_FAIL(eslFAIL, errbuf, "bad decode prob %f for %s, seg start, ia=%d\n", *xc, p7_sparsemx_DecodeSpecial(s), i); 
      }
      for (z = 0; z < sm->n[i]; z++)       /* sparse main cells */
	for (s = 0; s < p7S_NSCELLS; dpc++, s++) 
	  if (! is_prob(*dpc, tol)) ESL_FAIL(eslFAIL, errbuf, "bad decode prob %f at i=%d,k=%d,%s", *xc, i,sm->k[i][z], p7_sparsemx_DecodeState(s));
      if (sm->n[i]) {                      /* specials on sparse row */
	for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	  if (! is_prob(*xc, tol)) { abort(); ESL_FAIL(eslFAIL, errbuf, "bad decode prob %f at i=%d,%s", *xc, i, p7_sparsemx_DecodeSpecial(s)); }
      }
    }
  return eslOK;
}

/* Function:  p7_sparsemx_Validate()
 * Synopsis:  Validate a sparse DP matrix.
 *
 * Purpose:   Validate the contents of sparse DP matrix <sx>.
 *            Return <eslOK> if it passes. Return <eslFAIL> if
 *            it fails, and set <errbuf> to contain an 
 *            explanation, if caller provides a non-<NULL>
 *            <errbuf>.
 *            
 *            Currently validation is only implemented for
 *            Forward, Backward, Viterbi, and Decoding matrix
 *            types; not for Masstrace or Envscore.
 *
 * Args:      sx      - sparse DP matrix to validate
 *            errbuf  - char[eslERRBUFSIZE] space for error msg, or NULL.
 *
 * Returns:   <eslOK> on success.
 *            <eslFAIL> on failure, with an error message in <errbuf>
 *            if <errbuf> was provided.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_sparsemx_Validate(P7_SPARSEMX *sx, char *errbuf)
{
  int status;

  if ( (status = validate_dimensions(sx, errbuf)) != eslOK) return status;
  if ( (status = validate_no_nan    (sx, errbuf)) != eslOK) return status;

  switch (sx->type) {
  case p7S_UNSET:      ESL_FAIL(eslFAIL, errbuf, "validating an unset sparse DP matrix? probably not what you meant");
  case p7S_FORWARD:    if ( (status = validate_fwdvit  (sx, errbuf)) != eslOK) return status; break;
  case p7S_BACKWARD:   if ( (status = validate_backward(sx, errbuf)) != eslOK) return status; break;
  case p7S_DECODING:   if ( (status = validate_decoding(sx, errbuf)) != eslOK) return status; break;
  case p7S_ALIGNMENT:  ESL_FAIL(eslFAIL, errbuf, "unimplemented");
  case p7S_VITERBI:    if ( (status = validate_fwdvit  (sx, errbuf)) != eslOK) return status; break;
  case p7S_MASSTRACE:  ESL_FAIL(eslFAIL, errbuf, "unimplemented");
  case p7S_ENVSCORE:   ESL_FAIL(eslFAIL, errbuf, "unimplemented");
  default:             ESL_FAIL(eslFAIL, errbuf, "no such spare DP matrix type %d", sx->type);
  }
  return eslOK;
}
/*----------------- end, P7_SPARSEMX validation -----------------*/



/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
