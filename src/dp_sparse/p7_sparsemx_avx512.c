///AVX-512 routines for manipulating P7_SPARSEMX data structures

/* 
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "base/p7_trace.h"

#include "dp_vector/simdvec.h"	    /* #define's of SIMD vector sizes */

#include "dp_reference/p7_refmx.h"
#include "dp_sparse/p7_sparsemx.h"


/*****************************************************************
 * 2. API for constructing a sparse DP matrix
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
p7_sparsemask_Create_avx512(int M, int L)
{
#ifdef HAVE_AVX512  
  P7_SPARSEMASK *sm             = NULL;
  int            default_salloc = 8;
  int64_t        default_kalloc = 4096;
  int            i,r;
  int            status;

  ESL_ALLOC(sm, sizeof(P7_SPARSEMASK));
  sm->L      = L;
  sm->M      = M;
  sm->simd = AVX512; 

  sm->Q_AVX_512      = P7_NVF_AVX_512(M);  // approx M/4, for striped vectors of four floats
    
  sm->seg_AVX_512    = NULL;
  sm->k_AVX_512      = NULL;
  sm->n_AVX_512      = NULL;
  sm->kmem_AVX_512   = NULL;

  sm->S_AVX_512       = 0;
  sm->nrow_AVX_512    = 0;
  sm->ncells_AVX_512  = 0;
  sm->last_i_AVX_512  = L+1;       // sentinel to assure StartRow() is called in reverse L..1 order 
  for (r = 0; r < p7_VNF_AVX_512; r++) 
    sm->last_k_AVX_512[r]  = -1;        // sentinels to assure StartRow() is called before Add() 
  /* sn[] are initialized for each sparse row by _StartRow() */

  /* if Ws is really large, we might already know we need a
   * bigger-than-default allocation, just to enable the slots.
   * Rather than allocating the default and letting StartRow()
   * reallocate for the slots, go ahead and figure this out now.
   */
  sm->kalloc_AVX_512 = default_kalloc;
  while (sm->kalloc_AVX_512 < p7_VNF_AVX_512*sm->Q_AVX_512) sm->kalloc_AVX_512 *= 2;

  sm->ralloc_AVX_512   = L+1;   
  sm->salloc_AVX_512   = default_salloc;

 

  ESL_ALLOC(sm->seg_AVX_512,  sm->salloc_AVX_512 * sizeof(p7_sparsemask_seg_s)); // salloc is the actual allocation, inclusive of +2 for sentinels
  ESL_ALLOC(sm->k_AVX_512,    sm->ralloc_AVX_512 * sizeof(int *));
  ESL_ALLOC(sm->n_AVX_512,    sm->ralloc_AVX_512 * sizeof(int));
  ESL_ALLOC(sm->kmem_AVX_512, sm->kalloc_AVX_512 * sizeof(int));

  sm->k_AVX_512[0]   = NULL;    // always. 
  for (i = 0; i <= L; i++)  // n[0] will always be 0; n[i=1..L] initialized to 0, then count as cells are added 
    sm->n_AVX_512[i] = 0;

  sm->n_krealloc = 0;
  sm->n_rrealloc = 0;
  sm->n_srealloc = 0;
  return sm;

 ERROR:
  p7_sparsemask_Destroy(sm);
  return NULL;
#endif //HAVE_AVX512
#ifndef  HAVE_AVX512
  return NULL
#endif  
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
p7_sparsemask_Reinit_avx512(P7_SPARSEMASK *sm, int M, int L)
{
 #ifdef HAVE_AVX512 
  int i,r;
  int status;

  sm->L  = L;
  sm->M  = M;
  sm->Q_AVX_512  = P7_NVF_AVX_512(M);
    
  /* seg[], kmem stay at their previous salloc, kalloc
   * but do we need to reallocate rows for k[] and n[]? 
   */
  if (sm->ralloc_AVX_512 < L+1) {
    ESL_REALLOC(sm->k_AVX_512, sizeof(int *) * (L+1));
    ESL_REALLOC(sm->n_AVX_512, sizeof(int)   * (L+1));
    sm->ralloc_AVX_512 = L+1;
    sm->n_rrealloc++;
  }

  sm->S_AVX_512       = 0;
  sm->nrow_AVX_512    = 0;
  sm->ncells_AVX_512  = 0;
  sm->last_i_AVX_512  = sm->L+1;
  for (r = 0; r < p7_VNF_AVX_512; r++) 
    sm->last_k_AVX_512[r]  = -1; 
  /* sn[] are initialized for each sparse row by _StartRow() */

  /* The realloc counters are NOT reset. They keep accumulating during
   * the life of the object. 
   */
  for (i = 1; i <= L; i++)  /* n[0] will always be 0, but reinit n[1..L] */
    sm->n_AVX_512[i] = 0;

  return eslOK;

 ERROR:
  return status;
#endif //HAVE_AVX512
#ifndef HAVE_AVX512
return eslENORESULT;
#endif  
}

/* Function:  p7_sparsemask_Sizeof()
 * Synopsis:  Returns current allocated size of a <P7_SPARSEMASK>, in bytes.
 */
size_t
p7_sparsemask_Sizeof_avx512(const P7_SPARSEMASK *sm)
{
  size_t n = sizeof(P7_SPARSEMASK);
#ifdef HAVE_AVX512 
  n += sm->salloc_AVX_512 * sizeof(p7_sparsemask_seg_s); // <seg>
  n += sm->ralloc_AVX_512 * sizeof(int *);                      // <k>                   
  n += sm->ralloc_AVX_512 * sizeof(int);                        // <n>                   
  n += sm->kalloc_AVX_512 * sizeof(int);                        // <kmem>                
#endif  
  return n;
}

/* Function:  p7_sparsemask_MinSizeof()
 * Synopsis:  Returns minimum required size of a <P7_SPARSEMASK>, in bytes.
 */
size_t
p7_sparsemask_MinSizeof_avx512(const P7_SPARSEMASK *sm)
{
  size_t n = sizeof(P7_SPARSEMASK);
#ifdef HAVE_AVX512  
  n += (sm->S_AVX_512+2)  * sizeof(struct p7_sparsemask_seg_s);  // <seg>; includes sentinels at 0,S+1
  n += (sm->L+1)  * sizeof(int *);                       // <k>
  n += (sm->L+1)  * sizeof(int);                         // <n>
  n += sm->ncells_AVX_512 * sizeof(int);                         // <kmem>
#endif 
  return n;
}

/* Function:  p7_sparsemask_Destroy()
 * Synopsis:  Destroy a <P7_SPARSEMASK>.
 */
void
p7_sparsemask_Destroy_avx512(P7_SPARSEMASK *sm)
{
  if (sm) {
#ifdef HAVE_AVX512   
    if (sm->seg_AVX_512)  free(sm->seg_AVX_512);
    if (sm->k_AVX_512)    free(sm->k_AVX_512);
    if (sm->n_AVX_512)    free(sm->n_AVX_512);
    if (sm->kmem_AVX_512) free(sm->kmem_AVX_512);
#endif      
    free(sm);
  }
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
p7_sparsemask_StartRow_avx512(P7_SPARSEMASK *sm, int i)
{
 #ifdef HAVE_AVX512 
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
  if (sm->ncells_AVX_512 + p7_VNF_AVX_512*sm->Q_AVX_512 > sm->kalloc_AVX_512)
    {
      int64_t kalloc_req = sm->kalloc_AVX_512 * 2;
      ESL_REALLOC(sm->kmem_AVX_512, sizeof(int) * kalloc_req);
      sm->kalloc_AVX_512 = kalloc_req;
      sm->n_krealloc++;
    }
  
  for (r = 0; r < p7_VNF_AVX_512; r++)
    {
      sm->s_AVX_512[p7_VNF_AVX_512-r-1] = sm->kmem_AVX_512 + sm->ncells_AVX_512 + r*sm->Q_AVX_512;
      sm->sn_AVX_512[r]         = 0;
    }
  sm->last_i_AVX_512 = i;
  for (r = 0; r < p7_VNF_AVX_512; r++) 
    sm->last_k_AVX_512[r] = sm->M+1;    /* sentinel to be sure that Add() is called in reverse order M..1 */
  return eslOK;
  
 ERROR:
  return status;
  #endif //HAVE_AVX512
  #ifndef HAVE_AVX512
  return eslENORESULT;
  #endif
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
p7_sparsemask_Add_avx512(P7_SPARSEMASK *sm, int q, int r)
{
  #ifdef HAVE_AVX512
  int     k = r*sm->Q_AVX_512+q+1;

  //printf("Sparsemask_AVX_512 Adding i=%d q=%d r=%d k=%d M=%d\n", sm->last_i, q, r, k, sm->M);

#ifdef p7_DEBUGGING 
  if (q < 0 || q >= sm->Q_AVX_512)  ESL_EXCEPTION(eslEINVAL, "q is 0..Q_AVX-1; striped vector index");
  if (r < 0 || r >= p7_VNF_AVX_512) ESL_EXCEPTION(eslEINVAL, "r is 0..%d-1; segment index in vectors", p7_VNF_AVX);
  if (sm->last_k_AVX_512[r] == -1)  ESL_EXCEPTION(eslEINVAL, "need to StartRow_AVX() before calling Add_AVX()");
  if (sm->last_k_AVX_512[r] <= k)   ESL_EXCEPTION(eslEINVAL, "cells must be added in reverse order M..1");
#endif

  sm->s_AVX_512[r][sm->sn_AVX_512[r]] = k;
  sm->sn_AVX_512[r]++;
  sm->last_k_AVX_512[r] = k;
  return eslOK;
#endif //HAVE_AVX512
#ifndef HAVE_AVX512
  return eslENORESULT;
#endif
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
p7_sparsemask_FinishRow_avx512(P7_SPARSEMASK *sm)
{
 #ifdef HAVE_AVX512 
  int *p;
  int  r;
//printf("sm->ncells_AVX = %li\n", sm->ncells_AVX);
  /* s[7] is already where it belongs; so we start at p = kmem + ncells + sn[7] */
  p = sm->kmem_AVX_512 + sm->ncells_AVX_512 + sm->sn_AVX_512[p7_VNF_AVX_512-1];
  sm->n_AVX_512[sm->last_i_AVX_512] = sm->sn_AVX_512[p7_VNF_AVX_512-1];
  for (r = p7_VNF_AVX_512-2; r >= 0; r--)
    {
      memmove(p, sm->s_AVX_512[r], sizeof(int) * sm->sn_AVX_512[r]);
      p += sm->sn_AVX_512[r];
      sm->n_AVX_512[sm->last_i_AVX_512] += sm->sn_AVX_512[r];
    }
  /* now the slots are invalid; the next StartRow() will reset them */
  sm->ncells_AVX_512 += sm->n_AVX_512[sm->last_i_AVX_512];
  for (r = 0; r < p7_VNF_AVX_512; r++) 
    sm->last_k_AVX_512[r]  = -1;  /* that'll suffice to prevent Add() from being called after FinishRow(). */
  return eslOK;

#endif // HAVE_AVX512
#ifndef HAVE_AVX512
return eslENORESULT;
#endif
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
 *            into <kmem>, the <seg[]> coord pairs ia,ib for each
 *            segment, the <S> counter for segment number, and the <nrow>
 *            counter for the nonzero <n[i]>.
 *            
 *            Sentinels for <seg>: <seg[0].ia> and <seg[0].ib> are set
 *            to -1. For some boundary conditions, we need seg[0].ib < 
 *            seg[1].ia to always succeed; for another, we need seg[0].ib +1
 *            to be 0 (so seg[s].ib+1 always starts us at a valid i in n[i]).
 *            <seg[S+1].ia> and
 *            <seg[S+1].ib> are set to <L+2>; for some boundary
 *            conditions, we need a test of <i < seg[S+1].ia-1> to
 *            fail for all i up to L.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. Now <sm> may only be
 *            safely destroyed; its contents are otherwise undefined.
 */

int
p7_sparsemask_Finish_avx512(P7_SPARSEMASK *sm)
{
#ifdef HAVE_AVX512  
  int i,r;
  int s;
  int status;
//printf("calling p7_sparsemask_Finish, sm->ncells = %li\n", sm->ncells);

  /* Reverse kmem. */

  int *p_AVX_512;
  esl_vec_IReverse(sm->kmem_AVX_512, sm->kmem_AVX_512, sm->ncells_AVX_512);

  /* Set the k[] pointers; count <S> and <nrow> */
  p_AVX_512 = sm->kmem_AVX_512;
  sm->S_AVX_512 = sm->nrow_AVX_512 = 0;

  for (i = 1; i <= sm->L; i++){
     if (sm->n_AVX_512[i]) 
      {
  sm->nrow_AVX_512++;
  sm->k_AVX_512[i] = p_AVX_512;
  p_AVX_512       += sm->n_AVX_512[i];
  if (sm->n_AVX_512[i-1] == 0) sm->S_AVX_512++;
      } 
    else 
      sm->k_AVX_512[i] = NULL;
  }

 
  /* Reallocate seg[] if needed. */

  if ( (sm->S_AVX_512+2) > sm->salloc_AVX_512) 
    {
      ESL_REALLOC(sm->seg_AVX_512, (sm->S_AVX_512+2) * sizeof(p7_sparsemask_seg_s)); /* +2, for sentinels */
      sm->salloc_AVX_512 = sm->S_AVX_512 + 2; // inclusive of sentinels
      sm->n_srealloc++;
    }
      
  /* Set seg[] coord pairs. */
  sm->seg_AVX_512[0].ia = sm->seg_AVX_512[0].ib = -1;
  for (s = 1, i = 1; i <= sm->L; i++)
    {
      if (sm->n_AVX_512[i]   && sm->n_AVX_512[i-1] == 0)                 sm->seg_AVX_512[s].ia   = i; 
      if (sm->n_AVX_512[i]   && (i == sm->L || sm->n_AVX_512[i+1] == 0)) sm->seg_AVX_512[s++].ib = i; 
    }
  ESL_DASSERT1(( s == sm->S_AVX_512+1 ));
  sm->seg_AVX_512[s].ia = sm->seg_AVX_512[s].ib = sm->L+2;

   sm->last_i_AVX_512 = -1;
  for (r = 0; r < p7_VNF_AVX_512; r++) 
    sm->last_k_AVX_512[r] = -1;
 
 // if we're running AVX-512 code and not SSE, need to copy some values into the SSE data structure
  // so the downstream code will see them
  sm->seg = sm->seg_AVX_512;
  sm->k = sm->k_AVX_512;
  sm->n = sm->n_AVX_512;
  sm->kmem = sm->kmem_AVX_512;
  sm->S = sm->S_AVX_512;
  sm->nrow = sm->nrow_AVX_512;
  sm->ncells = sm->ncells_AVX_512; 

  return eslOK;

 ERROR:
  return eslEMEM;
  #endif
#ifndef HAVE_AVX512
return eslENORESULT;
#endif
}
/*----------------- end, P7_SPARSEMASK API ----------------------*/
/* Function:  p7_sparsemask_Dump()
 * Synopsis:  Dump contents of a P7_SPARSEMASK
 *
 * Purpose:   Dump the contents of the sparse mask <sm> to 
 *            open stream <ofp> for debugging, diagnostics.
 */
int
p7_sparsemask_Dump_avx512(FILE *ofp, P7_SPARSEMASK *sm)
{
  int i,k,z;
#ifdef HAVE_AVX512
  fprintf(ofp, "# sparse mask: M=%d L=%d Q=%d\n", sm->M, sm->L, sm->Q);
  fputs("     ", ofp);  for (k = 1; k <= sm->M; k++) fprintf(ofp, "%3d ", k);  fputs(" n \n", ofp);
  fputs("     ", ofp);  for (k = 1; k <= sm->M; k++) fputs("--- ", ofp);       fputs("---\n", ofp);

  for (i = 1; i <= sm->L; i++)
    {
      fprintf(ofp, "%3d: ", i);
      for (z = 0, k = 1; k <= sm->M; k++)
  {
    while (z < sm->n_AVX_512[i] && sm->k_AVX_512[i][z] < k)  z++;
    if    (z < sm->n_AVX_512[i] && sm->k_AVX_512[i][z] == k) fprintf(ofp, "  X ");
    else                                     fprintf(ofp, "  . ");
  }
      fprintf(ofp, "%3d\n", sm->n_AVX_512[i]);
    }
  return eslOK;
#endif
#ifndef HAVE_AVX512
return eslENORESULT;
#endif
}

/* Function:  p7_sparsemask_Compare()
 * Synopsis:  Compare two sparse masks for equality.
 *
 * Purpose:   Compare <sm1> and <sm2>; return <eslOK> if they
 *            are equal, <eslFAIL> if they are not.
 */
int
p7_sparsemask_Compare_avx512(const P7_SPARSEMASK *sm1, const P7_SPARSEMASK *sm2)
{
#ifdef HAVE_AVX512  
  char msg[] = "P7_SPARSEMASK comparison failed";
  int  i;
  int  s;
 if(sm2->simd != SSE){
    ESL_FAIL(eslFAIL, NULL, "Can't compare sparsemasks generated for different SIMD instruction sets");
  }
  if ( (sm1->L      != sm2->L)      ||
       (sm1->M      != sm2->M)      ||
       (sm1->S_AVX_512     != sm2->S_AVX_512)      ||
       (sm1->nrow_AVX_512   != sm2->nrow_AVX_512)   ||
       (sm1->ncells_AVX_512 != sm2->ncells_AVX_512)) 
    ESL_FAIL(eslFAIL, NULL, msg);

  for (s = 0; s <= sm1->S_AVX_512+1; s++)
    {
      if (sm1->seg_AVX_512[s].ia != sm2->seg_AVX_512[s].ia)   ESL_FAIL(eslFAIL, NULL, msg);
      if (sm1->seg_AVX_512[s].ib != sm2->seg_AVX_512[s].ib)   ESL_FAIL(eslFAIL, NULL, msg);
    }
  if ( esl_vec_ICompare(sm1->n_AVX_512, sm2->n_AVX_512, sm1->L+1)    != eslOK)  ESL_FAIL(eslFAIL, NULL, msg);
  for (i = 0; i <= sm1->L; i++)
    if ( esl_vec_ICompare(sm1->k_AVX_512[i], sm2->k_AVX_512[i], sm1->n_AVX_512[i]) != eslOK) ESL_FAIL(eslFAIL, NULL, msg);
  return eslOK;
  #endif
#ifndef HAVE_AVX512
return eslENORESULT;
#endif
}


/* Function:  p7_sparsemask_Validate()
 * Synopsis:  Validate a P7_SPARSEMASK sparse DP mask.
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
p7_sparsemask_Validate_avx512(const P7_SPARSEMASK *sm, char *errbuf)
{
#ifdef HAVE_AVX512
  int g, i;

  if (errbuf) errbuf[0] = '\0';

  if ( sm->L < 1) ESL_FAIL(eslFAIL, errbuf, "L must be >=1");
  if ( sm->M < 1) ESL_FAIL(eslFAIL, errbuf, "M must be >=1");
  if ( sm->S_AVX_512 < 0) ESL_FAIL(eslFAIL, errbuf, "S must be >=0");

  for (g = 1; g <= sm->S_AVX_512; g++)
    {
      if (sm->seg_AVX_512[g-1].ib >= sm->seg_AVX_512[g].ia)           ESL_FAIL(eslFAIL, errbuf, "seg %d overlaps with previous one", g);  // Note boundary condition, seg[0].ib=-1
      if (sm->seg_AVX_512[g].ia   >  sm->seg_AVX_512[g].ib)           ESL_FAIL(eslFAIL, errbuf, "ia..ib are not in order for seg %d", g);
      if (sm->seg_AVX_512[g].ia < 1 || sm->seg_AVX_512[g].ia > sm->L) ESL_FAIL(eslFAIL, errbuf, "ia[%d] is invalid", g);
      if (sm->seg_AVX_512[g].ib < 1 || sm->seg_AVX_512[g].ib > sm->L) ESL_FAIL(eslFAIL, errbuf, "ib[%d] is invalid", g);

      for (i = sm->seg_AVX_512[g-1].ib+1; i < sm->seg_AVX_512[g].ia; i++)   // Note boundary condition. Sentinel seg[0].ib == -1, so (i = seg[0]+1) means 0
  if (sm->n_AVX_512[i] != 0) ESL_FAIL(eslFAIL, errbuf, "n[i] != 0 for i unmarked, not in sparse segment");
      for (i = sm->seg_AVX_512[g].ia; i <= sm->seg_AVX_512[g].ib; i++)
  if (sm->n_AVX_512[i] == 0) ESL_FAIL(eslFAIL, errbuf, "n[i] == 0 for i supposedly marked in sparse seg");
    }
  for (i = sm->seg_AVX_512[sm->S_AVX_512].ib+1; i <= sm->L; i++)
    if (sm->n_AVX_512[i] != 0) ESL_FAIL(eslFAIL, errbuf, "n[i] != 0 for i unmarked, not in sparse segment");

  return eslOK;
  #endif
#ifndef HAVE_AVX512
return eslENORESULT;
#endif
}

/* Function:  p7_sparsemask_SetFromTrace()
 * Synopsis:  Set a sparse mask to contain cells in given trace, plus random scatter of others.
 *
 * Purpose:   Add every supercell <i,k> in trace <tr> to the sparse mask <sm>.
 *            
 *            If <rng> is provided (i.e. non-<NULL>), on rows <i> with
 *            at least one such cell, and on 20% of empty rows, also
 *            mark random sparse supercells with 50% probability
 *            each. This creates a sparse mask in which the path
 *            defined by <tr> is marked and can be scored by sparse DP
 *            routines; plus additional random cells, to try to
 *            exercise possible failure modes.
 *            
 */
int
p7_sparsemask_SetFromTrace_avx512(P7_SPARSEMASK *sm, ESL_RANDOMNESS *rng, const P7_TRACE *tr)
{
#ifdef HAVE_AVX512  
  float cellprob = 0.5;
  float rowprob  = 0.2;
  int   i,k,z;
  int   status;
  z = tr->N-1;
  for (i = sm->L; i >= 1; i--) /* sparsemask api requires building it backwards */
    {
      while (tr->i[z] != i) z--; /* find trace position that generated this residue. */
    
      if ( (status = p7_sparsemask_StartRow(sm, i)) != eslOK) return status;

      /* If this residue was emitted by the model, at least that cell
       * must be present; thus the row must be present. 
       * Tricky: in addition to the actual emitting cell i,k, we may
       * also need to add one or more delete cells i,k-1... 
       */
      if (p7_trace_IsM(tr->st[z]) || p7_trace_IsI(tr->st[z])) 
  {
    while (p7_trace_IsD(tr->st[z+1])) z++;
    
    for (k = sm->M; k > tr->k[z]; k--) 
      if (rng && esl_random(rng) < cellprob)  
        if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q_AVX_512, (k-1)/sm->Q_AVX_512)) != eslOK) return status;

    while (p7_trace_IsD(tr->st[z])) {
      k = tr->k[z]; 
      if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q_AVX_512, (k-1)/sm->Q_AVX_512)) != eslOK) return status;
      z--;
    }

    k = tr->k[z];
    if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q_AVX_512, (k-1)/sm->Q_AVX_512)) != eslOK) return status;
    
    for (k = k-1; k >= 1; k--)
      if (rng && esl_random(rng) < cellprob)
        if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q_AVX_512, (k-1)/sm->Q_AVX_512)) != eslOK) return status;
  }
      else
  {
    if (rng && esl_random(rng) < rowprob)
      for (k = sm->M; k >= 1; k--)
        if (rng && esl_random(rng) < cellprob)
    if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q_AVX_512, (k-1)/sm->Q_AVX_512)) != eslOK) return status; /* append to k[i] list, increment n[i] count, reallocating as needed; doesn't deal w/ segments (nrow,nseg,i[]) */
  }

      if ((status = p7_sparsemask_FinishRow(sm)) != eslOK) return status;
    }
  if ( (status = p7_sparsemask_Finish(sm)) != eslOK) return status;

  return eslOK;
  #endif
#ifndef HAVE_AVX512
return eslENORESULT;
#endif
}
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
