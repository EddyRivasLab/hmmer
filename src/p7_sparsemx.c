/* 
 * Data structures used by sparse dynamic programming.
 *
 * P7_SPARSEMASK and P7_SPARSEMX work together. P7_SPARSEMASK defines
 * which cells are included in the sparse DP matrix. P7_SPARSEMX is
 * the sparse DP matrix. The caller constructs a P7_SPARSEMASK first,
 * then uses it to create or reinitialize a P7_SPARSEMX, before
 * running any sparse DP routines (see sparse_fwdback.c). 
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
 *   5. Debugging tools for P7_SPARSEMX
 *   6. Copyright and license information  
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

int
p7_sparsemask_Dump(FILE *ofp, P7_SPARSEMASK *sm)
{
  int i,k,z;

  fprintf(ofp, "# sparse mask: M=%d L=%d\n", sm->M, sm->L);
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


int
p7_sparsemx_Reuse(P7_SPARSEMX *sx)
{
  sx->sm   = NULL;
  sx->type = p7S_UNSET;
  return eslOK;
}


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
 * 5. Debugging tools for P7_SPARSEMX
 *****************************************************************/ 

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
/*------------ end, P7_SPARSEMX debugging tools -----------------*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
