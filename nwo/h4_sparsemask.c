
#include "h4_config.h"

#include "easel.h"

#include "h4_sparsemask.h"

/*****************************************************************
 * 1. H4_SPARSEMASK; defines cells to be included in sparse DP matrix
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
