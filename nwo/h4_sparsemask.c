/* H4_SPARSEMASK: marks cells included in sparse dynamic programming
 *
 * See h4_sparsemask.md for notes.
 *
 * Contents:
 *    1. H4_SPARSEMASK object
 *    2. API for constructing a sparse DP mask (backwards)
 *    3. Debugging and development tools
 *    4. Unit tests
 *    5. Test driver
 *    6. Example
 */
#include "h4_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  // memmove()

#include "easel.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "simdvec.h"

#include "h4_path.h"
#include "h4_sparsemask.h"

/*****************************************************************
 * 1. H4_SPARSEMASK object
 *****************************************************************/

/* Function:  h4_sparsemask_Create()
 * Synopsis:  Creates a new H4_SPARSEMASK object.
 *
 * Purpose:   Create a new <H4_SPARSEMASK> for a comparison of a profile
 *            of length <M> to a sequence of length <L>. Return a ptr to
 *            the new object.
 *            
 *            The allocation will generally be for (much) less than <ML> cells;
 *            the API for creating the sparse mask will grow the structure
 *            appropriately. The structure does require at least $O(M)$ cells
 *            of temporary storage, in <V> "slots" used to sort input
 *            from striped vector code.  
 *
 * Args:      M       - model length
 *            L       - sequence length
 *
 * Returns:   a pointer to a new <H4_SPARSEMASK>
 *
 * Throws:    <NULL> on allocation failure.
 */
H4_SPARSEMASK *
h4_sparsemask_Create(int M, int L)
{
  H4_SPARSEMASK *sm             = NULL;
  int            default_salloc = 8;
  int64_t        default_kalloc = 4096;
  int            i,r;
  int            status;

  ESL_ALLOC(sm, sizeof(H4_SPARSEMASK));
  sm->L      = L;
  sm->M      = M;

  sm->V      = H4_V_FB;            
  sm->Q      = H4_Q(M, sm->V);   // V*Q ~ M; row of len M is striped into Q vectors of width V floats.
    
  sm->seg    = NULL;
  sm->k      = NULL;
  sm->n      = NULL;
  sm->kmem   = NULL;

  sm->S       = 0;
  sm->nrow    = 0;
  sm->ncells  = 0;
  sm->last_i  = L+1;       // sentinel to assure StartRow() is called in reverse L..1 order 

  for (r = 0; r < sm->V; r++) 
    sm->last_k[r]  = -1;        // sentinels to assure StartRow() is called before Add() 
  /* sn[] are initialized for each sparse row by _StartRow() */

  /* if Ws is really large, we might already know we need a
   * bigger-than-default allocation, just to enable the slots.
   * Rather than allocating the default and letting StartRow()
   * reallocate for the slots, go ahead and figure this out now.
   */
  sm->kalloc = default_kalloc;
  while (sm->kalloc < sm->V*sm->Q) sm->kalloc *= 2;

  sm->ralloc   = L+1;   
  sm->salloc   = default_salloc;

  ESL_ALLOC(sm->seg,  sm->salloc * sizeof(struct h4_sparsemask_seg_s)); // salloc is the actual allocation, inclusive of +2 for sentinels
  ESL_ALLOC(sm->k,    sm->ralloc * sizeof(int *));
  ESL_ALLOC(sm->n,    sm->ralloc * sizeof(int));
  ESL_ALLOC(sm->kmem, sm->kalloc * sizeof(int));

  sm->k[0]   = NULL;        // always. 
  for (i = 0; i <= L; i++)  // n[0] will always be 0; n[i=1..L] initialized to 0, then count as cells are added 
    sm->n[i] = 0;

  sm->n_krealloc = 0;
  sm->n_rrealloc = 0;
  sm->n_srealloc = 0;
  return sm;

 ERROR:
  h4_sparsemask_Destroy(sm);
  return NULL;
}


/* Function:  h4_sparsemask_Reinit()
 * Synopsis:  Reinitialize an existing H4_SPARSEMASK for a new comparison.
 *
 * Purpose:   Same as a <_Create()>, but reusing an existing 
 *            <H4_SPARSEMASK> to minimize reallocation calls.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
h4_sparsemask_Reinit(H4_SPARSEMASK *sm, int M, int L)
{
  int i,r;
  int status;

  sm->L  = L;
  sm->M  = M; 
  sm->V  = H4_V_FB;
  sm->Q  = H4_Q(M, sm->V);
    
  /* seg[], kmem stay at their previous salloc, kalloc
   * but check if we need to reallocate rows for k[] and n[] 
   */
  if (sm->ralloc < L+1) {
    ESL_REALLOC(sm->k, sizeof(int *) * (L+1));
    ESL_REALLOC(sm->n, sizeof(int)   * (L+1));
    sm->ralloc = L+1;
    sm->n_rrealloc++;
  }

  sm->S       = 0;
  sm->nrow    = 0;
  sm->ncells  = 0;
  sm->last_i  = sm->L+1;
  for (r = 0; r < sm->V; r++) 
    sm->last_k[r]  = -1; 
  /* sn[] are initialized for each sparse row by _StartRow() */

  /* The realloc counters are NOT reset. They keep accumulating during
   * the life of the object. 
   */
  for (i = 1; i <= L; i++)  /* n[0] will always be 0, but reinit n[1..L] */
    sm->n[i] = 0;

  return eslOK;

 ERROR:
  return status;
}

/* Function:  h4_sparsemask_Sizeof()
 * Synopsis:  Returns current allocated size of a <H4_SPARSEMASK>, in bytes.
 */
size_t
h4_sparsemask_Sizeof(const H4_SPARSEMASK *sm)
{
  size_t n = sizeof(H4_SPARSEMASK);   
  n += sm->salloc * sizeof(struct h4_sparsemask_seg_s); // <seg>
  n += sm->ralloc * sizeof(int *);                      // <k>                   
  n += sm->ralloc * sizeof(int);                        // <n>                   
  n += sm->kalloc * sizeof(int);                        // <kmem>                
  return n;
}


/* Function:  h4_sparsemask_MinSizeof()
 * Synopsis:  Returns minimum required size of a <H4_SPARSEMASK>, in bytes.
 * 
 * Purpose:   As opposed to <h4_sparsemask_Sizeof()>, which calculates
 *            the actual current allocated size of the structure
 *            (including overallocations), <_MinSizeof()> returns the
 *            amount of memory that's actually used and needed -- the
 *            minimum allocation that would hold the same information.
 *            
 *            This gets used when we're characterizing empirical
 *            memory performance, and we don't want to consider the
 *            overhead created by the reallocation-minimization
 *            strategies.
 */
size_t
h4_sparsemask_MinSizeof(const H4_SPARSEMASK *sm)
{
  size_t n = sizeof(H4_SPARSEMASK);
  n += (sm->S+2)  * sizeof(struct h4_sparsemask_seg_s);  // <seg>; includes sentinels at 0,S+1
  n += (sm->L+1)  * sizeof(int *);                       // <k>
  n += (sm->L+1)  * sizeof(int);                         // <n>
  n += sm->ncells * sizeof(int);                         // <kmem>
  return n;
}


/* Function:  h4_sparsemask_Reuse()
 * Synopsis:  Reuse a <H4_SPARSEMASK>.
 *
 * Purpose:   Essentially the equivalent of <_Destroy()>, but
 *            we will reuse the existing allocations. 
 *
 *            Currently a no-op. <_Reinit()> does all the work of
 *            reusing a <H4_SPARSEMASK>. Someday, though, we may want
 *            to have <_Reuse()> reallocate downwards, if a problem
 *            has temporarily required an unusually large amount of
 *            memory.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
h4_sparsemask_Reuse(H4_SPARSEMASK *sm)
{
  return eslOK;
}



/* Function:  h4_sparsemask_Destroy()
 * Synopsis:  Destroy a <H4_SPARSEMASK>.
 */
void
h4_sparsemask_Destroy(H4_SPARSEMASK *sm)
{
  if (sm) {
    if (sm->seg)  free(sm->seg);
    if (sm->k)    free(sm->k);
    if (sm->n)    free(sm->n);
    if (sm->kmem) free(sm->kmem);
    free(sm);
  }
}
/*-------------- end, H4_SPARSEMASK -----------------------------*/





/*****************************************************************
 * 2. API for constructing a sparse DP matrix
 *****************************************************************/


/* Function:  h4_sparsemask_StartRow()
 * Synopsis:  Prepare to store sparse cells on a new row i.
 *
 * Purpose:   Here we set up the vector "slots" (<sm->V> of them,
 *            for vectors V floats wide) for temporarily storing sorted
 *            lists of k indices for each striped vector segment.
 *            Each slot must allow up to Q entries; so kmem must be
 *            allocated for at least ncells+V*Q]. Reallocate
 *            <kmem> (by doubling) if needed. 
 * 
 *            Why do we need these shenanigans? The F/B filter uses
 *            striped vectors, which are not in the M..1 order that
 *            sparse DP needs. Temporary storage in V sorted slots,
 *            followed by their concatenation, is one way to rearrange
 *            efficiently. See note [3] in h4_sparsemx.md.
 *
 *            Remember, the whole <kmem> array is in reverse order during
 *            collection, so slot n[V-1] is first, n[0] is last; see note [1] in
 *            h4_sparsemx.md.
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
h4_sparsemask_StartRow(H4_SPARSEMASK *sm, int i)
{
  int r;
  int status;
  
#if eslDEBUGLEVEL > 0
  if (i < 1 || i > sm->L) ESL_EXCEPTION(eslEINVAL, "i must be 1..L: sequence position");
  if (sm->last_i <= i)    ESL_EXCEPTION(eslEINVAL, "rows need to be added in reverse order L..1");
#endif

  /* Make sure kmem has enough memory; if not, double it.
   * Because we know the original allocation was enough to hold
   * the slots, we know that doubling (even if ncells has filled
   * the current kalloc) is sufficient.
   */
  if (sm->ncells + sm->V*sm->Q > sm->kalloc)
    {
      int64_t kalloc_req = sm->kalloc * 2;
      ESL_REALLOC(sm->kmem, sizeof(int) * kalloc_req);
      sm->kalloc = kalloc_req;
      sm->n_krealloc++;
    }
  
  for (r = 0; r < sm->V; r++)
    {
      sm->s[sm->V-r-1] = sm->kmem + sm->ncells + r*sm->Q;
      sm->sn[r]        = 0;
    }
  sm->last_i = i;

  for (r = 0; r < sm->V; r++) 
    sm->last_k[r] = sm->M+1;            /* sentinel to be sure that Add() is called in reverse order M..1 */
  return eslOK;
  
 ERROR:
  return status;
}


/* Function:  h4_sparsemask_Add()
 * Synopsis:  Add a vector cell q,r to k list on current row.
 *
 * Purpose:   For a row in progress (we've already called <_StartRow()>),
 *            add a striped vector cell <q,r>. We must call <q>'s
 *            in descending order <Q-1..0>. The <r>'s in a given vector
 *            may be called in any order. 
 *            
 *            This is obviously designed for the production-code case,
 *            where we are constructing the <H4_SPARSEMASK> in the
 *            vectorized F/B filter's linear-memory backwards pass, we
 *            have striped vector coords q,r, and where we have to get
 *            the k indices back into order from their striped
 *            indexing. For testing/debugging code, where we want
 *            to store normal k indices (as opposed to q,r vector
 *            coordinates), the magic is to convert k to simulated
 *            vector coords: that's <q = (k-1)%sm->Q>, <r =
 *            (k-1)/sm->Q>. Though that uglifies calls to this
 *            function from our test code, it's better to do such
 *            contortions so that our test code goes through the same API
 *            as production uses.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
h4_sparsemask_Add(H4_SPARSEMASK *sm, int q, int r)
{
  int     k = r*sm->Q+q+1;

  ESL_DASSERT1(( q >= 0 && q < sm->Q )); // q is 0..Q-1, striped vector index
  ESL_DASSERT1(( r >= 0 && r < sm->V )); // r is 0..V-1, segment index in vectors
  ESL_DASSERT1(( sm->last_k[r] != -1 )); // need to StartRow() before calling Add()
  ESL_DASSERT1(( sm->last_k[r] >  k  )); // cells are added in reverse order M..1

  sm->s[r][sm->sn[r]] = k;
  sm->sn[r]++;
  sm->last_k[r] = k;
  return eslOK;
} 


/* Function:  h4_sparsemask_FinishRow()
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
h4_sparsemask_FinishRow(H4_SPARSEMASK *sm)
{
  int *p;
  int  r;

  /* s[V-1] is already where it belongs; so we start at p = kmem + ncells + sn[V-1] */
  p = sm->kmem + sm->ncells + sm->sn[sm->V-1];
  sm->n[sm->last_i] = sm->sn[sm->V-1];
  for (r = sm->V-2; r >= 0; r--)
    {
      memmove(p, sm->s[r], sizeof(int) * sm->sn[r]);
      p += sm->sn[r];
      sm->n[sm->last_i] += sm->sn[r];
    }

  /* now the slots are invalid; the next StartRow() will reset them */
  sm->ncells += sm->n[sm->last_i];
  for (r = 0; r < sm->V; r++) 
    sm->last_k[r]  = -1;        /* that'll suffice to prevent Add() from being called after FinishRow(). */
  return eslOK;
}


/* Function:  h4_sparsemask_Finish()
 * Synopsis:  Done adding cells to the mask.
 *
 * Purpose:   We have finished collecting an entire <H4_SPARSEMASK> for
 *            a given profile/sequence comparison we want to do using
 *            sparse DP routines. <kmem>, <n[i]>, and <ncells> fields
 *            are valid, though <kmem> is reversed.
 *            
 *            Here we reverse <kmem>, then set the rest: k[i] pointers
 *            into <kmem>, the <seg[]> coord pairs ia,ib for each
 *            segment, the <S> counter for segment number, and the <nrow>
 *            counter for the nonzero <n[i]>.
 *            
 *            Sentinels for <seg>: <seg[0].ia> and <seg[0].ib> are set
 *            to -1. For some boundary conditions, we need seg[0].ib < 
 *            seg[1].ia to always succeed; for another, we need seg[0].ib +1
 *            to be 0 (so seg[s].ib+1 always starts us at a valid i in n[i]).
 *            <seg[S+1].ia> and <seg[S+1].ib> are set to <L+2>; for some boundary
 *            conditions, we need a test of <i < seg[S+1].ia-1> to
 *            fail for all i up to L.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. Now <sm> may only be
 *            safely destroyed; its contents are otherwise undefined.
 */
int
h4_sparsemask_Finish(H4_SPARSEMASK *sm)
{
  int i,r;
  int s;
  int *p;
  int status;

  esl_vec_IReverse(sm->kmem, sm->kmem, sm->ncells);

  /* Set the k[] pointers; count <S> and <nrow> */
  p = sm->kmem;
  sm->S = sm->nrow = 0;

  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i]) 
        {
          sm->nrow++;
          sm->k[i] = p;
          p       += sm->n[i];
          if (sm->n[i-1] == 0) sm->S++;
        } 
      else 
        sm->k[i] = NULL;
    }

  /* Reallocate seg[] if needed. */
  if ( (sm->S+2) > sm->salloc) 
    {
      ESL_REALLOC(sm->seg, (sm->S+2) * sizeof(struct h4_sparsemask_seg_s)); // +2, for sentinels 
      sm->salloc = sm->S + 2;                                               // also inclusive of sentinels
      sm->n_srealloc++;
    }
      
  /* Set seg[] coord pairs. */
  sm->seg[0].ia = sm->seg[0].ib = -1;
  for (s = 1, i = 1; i <= sm->L; i++)
    {
      if (sm->n[i]   && sm->n[i-1] == 0)                 sm->seg[s].ia   = i; 
      if (sm->n[i]   && (i == sm->L || sm->n[i+1] == 0)) sm->seg[s++].ib = i; 
    }
  ESL_DASSERT1(( s == sm->S+1 ));
  sm->seg[s].ia = sm->seg[s].ib = sm->L+2;

   sm->last_i = -1;
   for (r = 0; r < sm->V; r++) 
    sm->last_k[r] = -1;

  return eslOK;

 ERROR:
  return eslEMEM;
}

/* Function:  h4_sparsemask_AddAll()
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
h4_sparsemask_AddAll(H4_SPARSEMASK *sm)
{
  int  i,k;
  int  status;

  for (i = sm->L; i >= 1; i--)
    {
      if (  (status = h4_sparsemask_StartRow(sm, i))                   != eslOK) return status;
      for (k = sm->M; k >= 1; k--)
        if ((status = h4_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) return status;
      if (  (status = h4_sparsemask_FinishRow(sm))                     != eslOK) return status;
    }
  return h4_sparsemask_Finish(sm);
}
/*----------------- end, H4_SPARSEMASK API ----------------------*/




/*****************************************************************
 * 3. Debugging tools for H4_SPARSEMASK
 *****************************************************************/

/* Function:  h4_sparsemask_Dump()
 * Synopsis:  Dump contents of a H4_SPARSEMASK
 *
 * Purpose:   Dump the contents of the sparse mask <sm> to 
 *            open stream <ofp> for debugging, diagnostics.
 */
int
h4_sparsemask_Dump(FILE *ofp, H4_SPARSEMASK *sm)
{
  int i,k,z;

  fprintf(ofp, "# sparse mask: M=%d L=%d Q=%d V=%d\n", sm->M, sm->L, sm->Q, sm->V);
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

/* Function:  h4_sparsemask_Compare()
 * Synopsis:  Compare two sparse masks for equality.
 *
 * Purpose:   Compare <sm1> and <sm2>; return <eslOK> if they
 *            are equal, <eslFAIL> if they are not.
 *            
 *            This does not look at the "slots", which are only used
 *            during sparse mask construction; only at the completed
 *            sparse mask. Therefore it is independent of
 *            vectorization (specifically, of vector width). You can
 *            compare masks created by different vector
 *            implementations.
 */
int
h4_sparsemask_Compare(const H4_SPARSEMASK *sm1, const H4_SPARSEMASK *sm2)
{
  char msg[] = "H4_SPARSEMASK comparison failed";
  int  i;
  int  s;

  if ( (sm1->L      != sm2->L)      ||
       (sm1->M      != sm2->M)      ||
       (sm1->S      != sm2->S)      ||
       (sm1->nrow   != sm2->nrow)   ||
       (sm1->ncells != sm2->ncells)) 
    ESL_FAIL(eslFAIL, NULL, msg);

  for (s = 0; s <= sm1->S+1; s++)
    {
      if (sm1->seg[s].ia != sm2->seg[s].ia)   ESL_FAIL(eslFAIL, NULL, msg);
      if (sm1->seg[s].ib != sm2->seg[s].ib)   ESL_FAIL(eslFAIL, NULL, msg);
    }
  if ( esl_vec_ICompare(sm1->n, sm2->n, sm1->L+1)    != eslOK)  ESL_FAIL(eslFAIL, NULL, msg);
  for (i = 0; i <= sm1->L; i++)
    if ( esl_vec_ICompare(sm1->k[i], sm2->k[i], sm1->n[i]) != eslOK) ESL_FAIL(eslFAIL, NULL, msg);

  return eslOK;
}



/* Function:  h4_sparsemask_Validate()
 * Synopsis:  Validate a H4_SPARSEMASK sparse DP mask.
 *
 * Purpose:   Validate the contents of sparse mask <sm>. 
 *            Return <eslOK> if it passes. Return <eslFAIL>
 *            if it fails, and set <errbuf> to contain an
 *            explanation, if caller provides a non-<NULL>
 *            <errbuf>.
 *
 * Args:      sm      - sparse DP mask to validate
 *            errbuf  - [eslERRBUFSIZE] space for an error msg; or NULL      
 *
 * Returns:   <eslOK> on success; <errbuf>, if provided, is set
 *            to an empty string "\0".
 *            
 *            <eslFAIL> on failure; <errbuf>, if provided, contains an
 *            informative error message.
 *            
 * Note:      We don't check for all possible invalidity; the goal of a
 *            Validate() is primarily to catch any future problems
 *            similar to past problems that we've already run across
 *            in debugging/testing.
 */
int
h4_sparsemask_Validate(const H4_SPARSEMASK *sm, char *errbuf)
{
  int g, i;

  if (errbuf) errbuf[0] = '\0';

  if ( sm->L < 1) ESL_FAIL(eslFAIL, errbuf, "L must be >=1");
  if ( sm->M < 1) ESL_FAIL(eslFAIL, errbuf, "M must be >=1");
  if ( sm->S < 0) ESL_FAIL(eslFAIL, errbuf, "S must be >=0");

  for (g = 1; g <= sm->S; g++)
    {
      if (sm->seg[g-1].ib >= sm->seg[g].ia)           ESL_FAIL(eslFAIL, errbuf, "seg %d overlaps with previous one", g);  // Note boundary condition, seg[0].ib=-1
      if (sm->seg[g].ia   >  sm->seg[g].ib)           ESL_FAIL(eslFAIL, errbuf, "ia..ib are not in order for seg %d", g);
      if (sm->seg[g].ia < 1 || sm->seg[g].ia > sm->L) ESL_FAIL(eslFAIL, errbuf, "ia[%d] is invalid", g);
      if (sm->seg[g].ib < 1 || sm->seg[g].ib > sm->L) ESL_FAIL(eslFAIL, errbuf, "ib[%d] is invalid", g);

      for (i = sm->seg[g-1].ib+1; i < sm->seg[g].ia; i++)   // Note boundary condition. Sentinel seg[0].ib == -1, so (i = seg[0]+1) means 0
        if (sm->n[i] != 0) ESL_FAIL(eslFAIL, errbuf, "n[i] != 0 for i unmarked, not in sparse segment");
      for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)
        if (sm->n[i] == 0) ESL_FAIL(eslFAIL, errbuf, "n[i] == 0 for i supposedly marked in sparse seg");
    }
  for (i = sm->seg[sm->S].ib+1; i <= sm->L; i++)
    if (sm->n[i] != 0) ESL_FAIL(eslFAIL, errbuf, "n[i] != 0 for i unmarked, not in sparse segment");

  return eslOK;
}

/* Function:  h4_sparsemask_SetFromTrace()
 * Synopsis:  Set a sparse mask to contain cells in given trace, plus random scatter of others.
 *
 * Purpose:   Add every supercell <i,k> in trace <tr> to the sparse mask <sm>.
 *            
 *            Dirty option: if <rng> is provided (i.e. non-<NULL>), on
 *            rows <i> with at least one such cell, and on 20% of
 *            empty rows, also mark random sparse supercells with 50%
 *            probability each. This creates a sparse mask in which
 *            the path defined by <tr> is marked and can be scored by
 *            sparse DP routines; plus additional random cells, to try
 *            to exercise possible failure modes.
 */
int
h4_sparsemask_SetFromTrace(H4_SPARSEMASK *sm, ESL_RANDOMNESS *rng, const H4_PATH *pi)
{
  int  *kend     = NULL;
  float cellprob = 0.3;
  float rowprob  = 0.2;
  int   i,k,kx,z,r;
  int   status;

  /* The API for building a sparse mask works backwards from i=L..1,
   * but path information is stored asymmetrically and is most easily
   * traversed i=1..L.  Precompute kend[0..Z-1]: k index of the last
   * M|I|D state in a run in run-length=encoded path. kend[z] = -1 for
   * other states z.
   */
  ESL_ALLOC(kend, sizeof(int) * pi->Z);
  for (z = 0; z < pi->Z; z++)
    {
      kend[z] = -1;
      if (pi->st[z] == h4P_L  || pi->st[z] == h4P_G)
	k  = pi->rle[z];   // enter at k. G->{MD}1 for G, simple enough; but we need L->Mk entry information to know k for following states
      else if (h4_path_IsI(pi->st[z]))
	kend[z] = k-1;     // k is index of next M|D state, so for an I state, use k-1.
      else if (h4_path_IsM(pi->st[z]) || h4_path_IsD(pi->st[z]))
	{
	  k += pi->rle[z]; // advance k by the runlength to be k index of the next M|D, if there is one
	  kend[z]  = k-1;  // ... then kend[z] for *this* M|D is k-1.
	}
    }

  /* Now we can build the sparse mask */
  i = sm->L;
  for (z = pi->Z-1; z >= 0; z--)  // sparsemask api requires building backwards
    {
      if (pi->st[z] == h4P_C || pi->st[z] == h4P_J || pi->st[z] == h4P_N)  // if i was emitted by N|C|J, we only add random cells
	{
	  for (r = 1; r < pi->rle[z]; r++)  // start at 1 not 0 because N/C/J account for rle-1 residues
	    {
	      if ( (status = h4_sparsemask_StartRow(sm, i)) != eslOK) return status;
	      if (rng && esl_random(rng) < rowprob)
		for (k = sm->M; k >= 1; k--)
		  if (rng && esl_random(rng) < cellprob)
		    if ((status = h4_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) return status; 
	      if ((status = h4_sparsemask_FinishRow(sm)) != eslOK) return status;
	      i--;
	    }
	}
      else if (h4_path_IsM(pi->st[z]) || h4_path_IsI(pi->st[z]))    // if i was emitted by M|I, we add that (i,k) cell plus random cells.
	{
	  for (r = 0; r < pi->rle[z]; r++, i--)
	    {
	      if ( (status = h4_sparsemask_StartRow(sm, i)) != eslOK) return status;

	      // initialize k at start of run and handle any following D's, as follows:
	      if (r == 0) {                         // we also need to add (i,k) cells for D, when D's follow M|I. (G->D is handled below)
		if (h4_path_IsD(pi->st[z+1])) z++;  // only the last M|I in a run can be followed by D's. bump z to the D run.
		k = kend[z];                        // k initializes here at start of a run - either the M|I run, or the D run if there is one.
	      }

	      // random cells at end of row i, down to k+1
	      for (kx = sm->M; kx > k; kx--)
		if (rng && esl_random(rng) < cellprob)
		  if ((status = h4_sparsemask_Add(sm, (kx-1)%sm->Q, (kx-1)/sm->Q)) != eslOK) return status;

	      // run of k's for a DDD path (r=0 only, where we bumped z up to the D's)
	      if (h4_path_IsD(pi->st[z]))
		{
		  for (r = 0; r < pi->rle[z]; r++, k--)
		    if ((status = h4_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) return status;
		  z--; // go back to the M|I state. k is now k for that state.
		}

	      // k is now k for M|I state in profile 1..M coords; z is now back on that M|I state even if we detoured over to D's.
	      if ((status = h4_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) return status;
	  
	      // remainder of row at start. On  G->DDD->{MI}, add *all* Dk; else, add random cells as we did on the end of the row.
	      if (pi->st[z-1] == h4P_DG && pi->st[z-2] == h4P_G)  // if z is M|I, we know z-1 and z-2 must always exist in path. (worst case: N-{GL}-{MI})
		{
		  for (kx = k-1; kx >= 1; kx--)
		    if ((status = h4_sparsemask_Add(sm, (kx-1)%sm->Q, (kx-1)/sm->Q)) != eslOK) return status;
		}
	      else
		{
		  for (kx = k-1; kx >= 1; kx--)
		    if (rng && esl_random(rng) < cellprob)
		      if ((status = h4_sparsemask_Add(sm, (kx-1)%sm->Q, (kx-1)/sm->Q)) != eslOK) return status;
		}

	      if ((status = h4_sparsemask_FinishRow(sm)) != eslOK) return status;
	      if (h4_path_IsM(pi->st[z])) k--;
	    }
	}
    }

  if ( (status = h4_sparsemask_Finish(sm)) != eslOK) return status;
  free(kend);
  return eslOK;

 ERROR:
  return status;
}
/*-------- end, H4_SPARSEMASK debugging tools -------------------*/


/***************************************************************** 
 * 4. Unit tests
 *****************************************************************/
#ifdef h4SPARSEMASK_TESTDRIVE

#include "esl_alphabet.h"
#include "esl_sq.h"

#include "h4_profile.h"
#include "h4_mode.h"

#include "emit.h"
#include "modelsample.h"

/* utest_sanity() 
 *
 * Basic sanity checks that stuff doesn't obviously crash.
 */
static void
utest_sanity(ESL_RANDOMNESS *rng)
{
  char           msg[] = "h4_sparsemask sanity unit test failed";
  ESL_ALPHABET  *abc   = esl_alphabet_Create(eslCOINS);  // unusual alphabet, just to add stress
  H4_PROFILE    *hmm   = NULL;
  H4_MODE       *mo    = h4_mode_Create();               // default multihit dual-mode
  H4_PATH       *pi    = h4_path_Create();
  H4_SPARSEMASK *sm    = h4_sparsemask_Create(10, 10);   // initial sizes irrelevant; we'll Reinit()
  ESL_SQ       *sq     = esl_sq_CreateDigital(abc);
  int            M     = 1 + esl_rnd_Roll(rng, 20);      // profile length 1..20
  int            nseq  = 10;
  char           errbuf[eslERRBUFSIZE];

  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  if ( h4_mode_SetLength(mo, 10)         != eslOK) esl_fatal(msg);

  while (nseq--)
    {
      if ( h4_emit(rng, hmm, mo, sq, pi)           != eslOK) esl_fatal(msg);
      if ( h4_sparsemask_Reinit(sm, hmm->M, sq->n) != eslOK) esl_fatal(msg);
      if ( h4_sparsemask_SetFromTrace(sm, rng, pi) != eslOK) esl_fatal(msg);

      if ( h4_sparsemask_Validate(sm, errbuf)      != eslOK) esl_fatal("%s\n  %s", msg, errbuf);

      h4_sparsemask_Reuse(sm);
      h4_path_Reuse(pi);
      esl_sq_Reuse(sq);
    }

  esl_sq_Destroy(sq);
  h4_sparsemask_Destroy(sm);
  h4_path_Destroy(pi);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
}

#endif // h4SPARSEMASK_TESTDRIVE


/***************************************************************** 
 * 5. Test driver
 *****************************************************************/
#ifdef h4SPARSEMASK_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                          docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help summary",                 0 },
  { "-s",         eslARG_INT,     "0", NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed to <n>", 0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv, "test driver for h4_sparsemask", "[-options]");
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng)); 

  utest_sanity(rng);
  
  fprintf(stderr, "#  status   = ok\n");

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif // h4SPARSEMASK_TESTDRIVE


/***************************************************************** 
 * 6. Example
 *****************************************************************/
#ifdef h4SPARSEMASK_EXAMPLE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "h4_profile.h"
#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "h4_sparsemask.h"

#include "emit.h"
#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                              docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",   0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL, NULL, NULL, NULL, "set random number generator seed",       0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show HMMER version info",                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 1, argc, argv, "example of using sparse mask", "[-options] <hmmfile>");
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  ESL_SQ         *sq      = NULL;
  H4_MODE        *mo      = h4_mode_Create();  // leave in default multihit dual-mode
  H4_PATH        *pi      = h4_path_Create();
  H4_SPARSEMASK  *sm      = NULL;
  char            errbuf[eslERRBUFSIZE];

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);

  h4_mode_SetLength(mo, 10);
  sq = esl_sq_CreateDigital(abc);

  h4_emit(rng, hmm, mo, sq, pi);

  sm = h4_sparsemask_Create(hmm->M, sq->n);
  h4_sparsemask_SetFromTrace(sm, rng, pi);

  if ( h4_sparsemask_Validate(sm, errbuf) != eslOK)
    esl_fatal("sparsemask validation failed:\n%s\n", errbuf);

  h4_sparsemask_Dump(stdout, sm);
  h4_path_Dump(stdout, pi);

  h4_hmmfile_Close(hfp);
  h4_sparsemask_Destroy(sm);
  h4_path_Destroy(pi);
  h4_mode_Destroy(mo);
  esl_sq_Destroy(sq);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  return 0;
}
#endif // h4SPARSEMASK_EXAMPLE
