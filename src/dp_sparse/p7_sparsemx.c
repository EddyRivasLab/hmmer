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
#include "hardware/hardware.h"


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
p7_sparsemask_Create(int M, int L, SIMD_TYPE simd)
{
switch(simd){
    case SSE:
      return p7_sparsemask_Create_sse(M, L);
      break;
    case AVX:
      return p7_sparsemask_Create_avx(M, L);
      break;
    case AVX512:
      return p7_sparsemask_Create_avx512(M, L);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_Create");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_Create");  
  }
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
 switch(sm->simd){
    case SSE:
      return p7_sparsemask_Reinit_sse(sm, M, L);
      break;
    case AVX:
      return p7_sparsemask_Reinit_avx(sm, M, L);
      break;
    case AVX512:
      return p7_sparsemask_Reinit_avx512(sm, M, L);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_Reinit");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_Reinit");  
  }
}

/* Function:  p7_sparsemask_Sizeof()
 * Synopsis:  Returns current allocated size of a <P7_SPARSEMASK>, in bytes.
 */
size_t
p7_sparsemask_Sizeof(const P7_SPARSEMASK *sm)
{
 switch(sm->simd){
    case SSE:
      return p7_sparsemask_Sizeof_sse(sm);
      break;
    case AVX:
      return p7_sparsemask_Sizeof_avx(sm);
      break;
    case AVX512:
      return p7_sparsemask_Sizeof_avx512(sm);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_Sizeof");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_Sizeof");  
  }
}

/* Function:  p7_sparsemask_MinSizeof()
 * Synopsis:  Returns minimum required size of a <P7_SPARSEMASK>, in bytes.
 */
size_t
p7_sparsemask_MinSizeof(const P7_SPARSEMASK *sm)
{
 switch(sm->simd){
    case SSE:
      return p7_sparsemask_MinSizeof_sse(sm);
      break;
    case AVX:
      return p7_sparsemask_MinSizeof_avx(sm);
      break;
    case AVX512:
      return p7_sparsemask_MinSizeof_avx512(sm);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_MinSizeof");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_MinSizeof");  
  }
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
  switch(sm->simd){
    case SSE:
      p7_sparsemask_Destroy_sse(sm);
      break;
    case AVX:
      p7_sparsemask_Destroy_avx(sm);
      break;
    case AVX512:
      p7_sparsemask_Destroy_avx512(sm);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_Destroy");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_Destroy");  
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

// Not sure what needs to be done to make this work for AVX, need to think about it more
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

// Separate versions of these functions for each ISA because they're called from 
// ISA-specific functions

int
p7_sparsemask_StartRow(P7_SPARSEMASK *sm, int i)
{
switch(sm->simd){
    case SSE:
      return p7_sparsemask_StartRow_sse(sm, i);
      break;
    case AVX:
      return p7_sparsemask_StartRow_avx(sm, i);
      break;
    case AVX512:
      return p7_sparsemask_StartRow_avx512(sm, i);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_StartRow");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_StartRow");  
  }
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
switch(sm->simd){
    case SSE:
      return p7_sparsemask_Add_sse(sm, q, r);
      break;
    case AVX:
      return p7_sparsemask_Add_avx(sm, q, r);
      break;
    case AVX512:
      return p7_sparsemask_Add_avx512(sm, q, r);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_Add");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_Add");  
  }
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
 switch(sm->simd){
    case SSE:
      return p7_sparsemask_FinishRow_sse(sm);
      break;
    case AVX:
      return p7_sparsemask_FinishRow_avx(sm);
      break;
    case AVX512:
      return p7_sparsemask_FinishRow_avx512(sm);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_FinishRow");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_FinishRow");  
  }
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
p7_sparsemask_Finish(P7_SPARSEMASK *sm)
{
 switch(sm->simd){
    case SSE:
      return p7_sparsemask_Finish_sse(sm);
      break;
    case AVX:
      return p7_sparsemask_Finish_avx(sm);
      break;
    case AVX512:
      return p7_sparsemask_Finish_avx512(sm);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_Finish");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_Finish");  
  }
}

//Old version
#ifdef NEVER_DEFINE_THIS 
int
p7_sparsemask_Finish(P7_SPARSEMASK *sm)
{
  int i,r;
  int s;
  int status;
//printf("calling p7_sparsemask_Finish, sm->ncells = %li\n", sm->ncells);

  /* Reverse kmem. */
#ifdef p7_build_SSE 
  int *p;
  esl_vec_IReverse(sm->kmem, sm->kmem, sm->ncells);

  /* Set the k[] pointers; count <S> and <nrow> */
  p = sm->kmem;
  sm->S = sm->nrow = 0;
#endif
#ifdef p7_build_AVX2 
  int *p_AVX;
  esl_vec_IReverse(sm->kmem_AVX, sm->kmem_AVX, sm->ncells_AVX);

  /* Set the k[] pointers; count <S> and <nrow> */
  p_AVX = sm->kmem_AVX;
  sm->S_AVX = sm->nrow_AVX = 0;
#endif
#ifdef p7_build_AVX512 
  int *p_AVX_512;
  esl_vec_IReverse(sm->kmem_AVX_512, sm->kmem_AVX_512, sm->ncells_AVX_512);

  /* Set the k[] pointers; count <S> and <nrow> */
  p_AVX_512 = sm->kmem_AVX_512;
  sm->S_AVX_512 = sm->nrow_AVX_512 = 0;
#endif

#ifdef p7_build_check_AVX2
if (sm->ncells != sm->ncells_AVX){
  printf("p7_sparsemask_Finish found different numbers of cells: %li (SSE) vs %li (AVX)\n", sm->ncells, sm->ncells_AVX);
}
int check_index;
for(check_index = 0; check_index < sm->ncells; check_index++){
  if(sm->kmem[check_index] != sm->kmem_AVX[check_index]){
    printf("Sparsemask index miss-match at position %d, %i (SSE) vs. %i (AVX)\n", check_index, sm->kmem[check_index], sm->kmem_AVX[check_index]);
  }
}
#endif  
#ifdef p7_build_check_AVX512
//printf("Checking AVX-512 sparsemask\n");
if (sm->ncells != sm->ncells_AVX_512){
  printf("p7_sparsemask_Finish found different numbers of cells: %li (SSE) vs %li (AVX-512)\n", sm->ncells, sm->ncells_AVX_512);
}
int check_index;
for(check_index = 0; check_index < sm->ncells; check_index++){
  if(sm->kmem[check_index] != sm->kmem_AVX_512[check_index]){
    printf("Sparsemask index miss-match at position %d, %i (SSE) vs. %i (AVX-512)\n", check_index, sm->kmem[check_index], sm->kmem_AVX_512[check_index]);
  }
}
#endif  
  for (i = 1; i <= sm->L; i++){
#ifdef p7_build_SSE  
    if (sm->n[i]) 
      {
  sm->nrow++;
  sm->k[i] = p;
  p       += sm->n[i];
  if (sm->n[i-1] == 0) sm->S++;
      } 
    else 
      sm->k[i] = NULL;
#endif
#ifdef p7_build_AVX2
     if (sm->n_AVX[i]) 
      {
  sm->nrow_AVX++;
  sm->k_AVX[i] = p_AVX;
  p_AVX       += sm->n_AVX[i];
  if (sm->n_AVX[i-1] == 0) sm->S_AVX++;
      } 
    else 
      sm->k_AVX[i] = NULL;
#endif
#ifdef p7_build_AVX512
     if (sm->n_AVX_512[i]) 
      {
  sm->nrow_AVX_512++;
  sm->k_AVX_512[i] = p_AVX_512;
  p_AVX_512       += sm->n_AVX_512[i];
  if (sm->n_AVX_512[i-1] == 0) sm->S_AVX_512++;
      } 
    else 
      sm->k_AVX_512[i] = NULL;
#endif
  }

 
  /* Reallocate seg[] if needed. */
#ifdef p7_build_SSE
  if ( (sm->S+2) > sm->salloc) 
    {
      ESL_REALLOC(sm->seg, (sm->S+2) * sizeof(p7_sparsemask_seg_s)); /* +2, for sentinels */
      sm->salloc = sm->S + 2; // inclusive of sentinels
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
#endif
#ifdef p7_build_AVX2
  if ( (sm->S_AVX+2) > sm->salloc_AVX) 
    {
      ESL_REALLOC(sm->seg_AVX, (sm->S_AVX+2) * sizeof(p7_sparsemask_seg_s)); /* +2, for sentinels */
      sm->salloc_AVX = sm->S_AVX + 2; // inclusive of sentinels
      sm->n_srealloc++;
    }
      
  /* Set seg[] coord pairs. */
  sm->seg_AVX[0].ia = sm->seg_AVX[0].ib = -1;
  for (s = 1, i = 1; i <= sm->L; i++)
    {
      if (sm->n_AVX[i]   && sm->n_AVX[i-1] == 0)                 sm->seg_AVX[s].ia   = i; 
      if (sm->n_AVX[i]   && (i == sm->L || sm->n_AVX[i+1] == 0)) sm->seg_AVX[s++].ib = i; 
    }
  ESL_DASSERT1(( s == sm->S_AVX+1 ));
  sm->seg_AVX[s].ia = sm->seg_AVX[s].ib = sm->L+2;
#endif  
#ifdef p7_build_AVX512
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
#endif  
#ifdef p7_build_check_AVX2
  if(sm->S != sm->S_AVX){
    printf("SSE reports %d segments, AVX reports %d\n", sm->S, sm->S_AVX);
  }
  for(i = 1; i <= sm->S; i++){
    if (sm->seg[i].ia != sm->seg_AVX[i].ia){
      printf("ia miss-match at seg[%d]: %d (SSE) vs %d (AVX)\n", i, sm->seg[i].ia, sm->seg_AVX[i].ia);
    }
    if (sm->seg[i].ib != sm->seg_AVX[i].ib){
      printf("ib miss-match at seg[%d]: %d (SSE) vs %d (AVX)\n", i, sm->seg[i].ib, sm->seg_AVX[i].ib);
    }
  }  
#endif   
 #ifdef p7_build_check_AVX512
 // printf("Checking AVX-512 sparsemask segments\n");
  if(sm->S != sm->S_AVX_512){
    printf("SSE reports %d segments, AVX-512 reports %d\n", sm->S, sm->S_AVX_512);
  }
  for(i = 1; i <= sm->S; i++){
    if (sm->seg[i].ia != sm->seg_AVX_512[i].ia){
      printf("ia miss-match at seg[%d]: %d (SSE) vs %d (AVX-512)\n", i, sm->seg[i].ia, sm->seg_AVX_512[i].ia);
    }
    if (sm->seg[i].ib != sm->seg_AVX_512[i].ib){
      printf("ib miss-match at seg[%d]: %d (SSE) vs %d (AVX-512)\n", i, sm->seg[i].ib, sm->seg_AVX_512[i].ib);
    }
  }  
#endif   
 #ifdef p7_build_SSE 
   sm->last_i = -1;
  for (r = 0; r < p7_VNF; r++) 
    sm->last_k[r] = -1;
 #endif  
#ifdef p7_build_AVX2
   sm->last_i_AVX = -1;
  for (r = 0; r < p7_VNF_AVX; r++) 
    sm->last_k_AVX[r] = -1;
 #endif  
#ifdef p7_build_AVX512
   sm->last_i_AVX_512 = -1;
  for (r = 0; r < p7_VNF_AVX_512; r++) 
    sm->last_k_AVX_512[r] = -1;
 #endif  
#ifndef p7_build_SSE
#ifdef p7_build_AVX2
  // if we're running AVX code and not SSE, need to copy some values into the SSE data structure
  // so the downstream code will see them
  sm->seg = sm->seg_AVX;
  sm->k = sm->k_AVX;
  sm->n = sm->n_AVX;
  sm->kmem = sm->kmem_AVX;
  sm->S = sm->S_AVX;
  sm->nrow = sm->nrow_AVX;
  sm->ncells = sm->ncells_AVX; 
#endif
#ifndef p7_build_AVX2
#ifdef p7_build_AVX512
 // if we're running AVX-512 code and not SSE, need to copy some values into the SSE data structure
  // so the downstream code will see them
  sm->seg = sm->seg_AVX_512;
  sm->k = sm->k_AVX_512;
  sm->n = sm->n_AVX_512;
  sm->kmem = sm->kmem_AVX_512;
  sm->S = sm->S_AVX_512;
  sm->nrow = sm->nrow_AVX_512;
  sm->ncells = sm->ncells_AVX_512; 
#endif
#endif  
#endif

  return eslOK;

 ERROR:
  return eslEMEM;
}
#endif // NEVER_DEFINE_THIS
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
 switch(sm->simd){
    case SSE:
      return p7_sparsemask_Dump_sse(ofp, sm);
      break;
    case AVX:
      return p7_sparsemask_Dump_avx(ofp, sm);
      break;
    case AVX512:
      return p7_sparsemask_Dump_avx512(ofp, sm);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_Dump");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_Dump");  
  }
}

/* Function:  p7_sparsemask_Compare()
 * Synopsis:  Compare two sparse masks for equality.
 *
 * Purpose:   Compare <sm1> and <sm2>; return <eslOK> if they
 *            are equal, <eslFAIL> if they are not.
 */
int
p7_sparsemask_Compare(const P7_SPARSEMASK *sm1, const P7_SPARSEMASK *sm2)
{
switch(sm1->simd){
    case SSE:
      return p7_sparsemask_Compare_sse(sm1, sm2);
      break;
    case AVX:
      return p7_sparsemask_Compare_avx(sm1, sm2);
      break;
    case AVX512:
      return p7_sparsemask_Compare_avx512(sm1, sm2);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_Compare");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_Compare");  
  }
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
p7_sparsemask_Validate(const P7_SPARSEMASK *sm, char *errbuf)
{
  switch(sm->simd){
    case SSE:
      return p7_sparsemask_Validate_sse(sm, errbuf);
      break;
    case AVX:
      return p7_sparsemask_Validate_avx(sm, errbuf);
      break;
    case AVX512:
      return p7_sparsemask_Validate_avx512(sm, errbuf);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_Validate");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_Validate");  
  }
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
p7_sparsemask_SetFromTrace(P7_SPARSEMASK *sm, ESL_RANDOMNESS *rng, const P7_TRACE *tr)
{
  switch(sm->simd){
    case SSE:
      return p7_sparsemask_SetFromTrace_sse(sm, rng, tr);
      break;
    case AVX:
      return p7_sparsemask_SetFromTrace_avx(sm, rng, tr);
      break;
    case AVX512:
      return p7_sparsemask_SetFromTrace_avx512(sm, rng, tr);
      break;
    case NEON:
      p7_Fail("Neon support not yet integrated into p7_sparsemask_SetFromTrace");
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_sparsemask_SetFromTrace");  
  }
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

  /* We must avoid zero-sized mallocs. If there are no rows or cells, alloc the default sizes */
  sx->dalloc = ( (sm && sm->ncells)         ? sm->ncells       : default_ncell);
  sx->xalloc = ( (sm && (sm->nrow + sm->S)) ? sm->nrow + sm->S : default_nx);

  ESL_ALLOC(sx->dp,  sizeof(float) * p7S_NSCELLS * sx->dalloc);
  ESL_ALLOC(sx->xmx, sizeof(float) * p7S_NXCELLS * sx->xalloc);
  return sx;

 ERROR:
  p7_sparsemx_Destroy(sx);
  return NULL;
}


/* Function:  p7_sparsemx_Reinit()
 * Synopsis:  Reinitialize, reallocate sparse DP matrix for new calculation.
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
p7_sparsemx_Reinit(P7_SPARSEMX *sx, const P7_SPARSEMASK *sm)
{
  int64_t dalloc_req = sm->ncells;
  int     xalloc_req = sm->nrow + sm->S;
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


int
p7_sparsemx_Zero(P7_SPARSEMX *sx)
{
  esl_vec_FSet(sx->dp,   sx->sm->ncells*p7S_NSCELLS,          0.0f);
  esl_vec_FSet(sx->xmx, (sx->sm->nrow+sx->sm->S)*p7S_NXCELLS, 0.0f);
  return eslOK;
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
  n += sizeof(float) * p7S_NSCELLS * sm->ncells;          // dp[]
  n += sizeof(float) * p7S_NXCELLS * (sm->nrow + sm->S);  // xmx[]; for each seg ia..ib, ia-1..ib has special cells
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

  for (g = 1; g <= sx->sm->S; g++)
    {
      ia = sx->sm->seg[g].ia;
      ib = sx->sm->seg[g].ib;
      
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

/* Function:  p7_sparsemx_TracePostprobs()
 * Synopsis:  Annotate trace with posterior probs from sparse decoding matrix
 *
 * Purpose:   Annotate trace <tr> with posterior probabilities taken from 
 *            the sparse decoding matrix in <sxd>.
 *            
 *            The trace <tr> must have an allocated <tr->pp> field;
 *            for example it should have been allocated with
 *            <p7_trace_CreateWithPP()>.
 *            
 *            The posterior probability is marginalized over
 *            local/glocal decision. x_i in a domain can be emitted by
 *            Mk or Ik.
 *
 * Returns:   <eslOK> on success, and now <tr> contains posterior probabilities
 *            for each emitted residue in <tr->pp>.
 */
int
p7_sparsemx_TracePostprobs(const P7_SPARSEMX *sxd, P7_TRACE *tr)
{
  const P7_SPARSEMASK *sm  = sxd->sm;
  const float         *dpc = sxd->dp;   // ptr that steps through stored main supercells 
  const float         *xc  = sxd->xmx;  // ptr that steps through stored special rows, including ia-1 segment edges 
  int i  = 0;                           // important to start at i=0, <xc> can start there 
  int k;	                        // index in main MI states, 0,1..M 
  int v  = 0;	                        // index that steps through sparse cell list on a row. (Initialization is to silence static code checkers.) 
  int z;	                        // index that steps through trace 
  int last_ib;	                        // last stored row index - all i>last_ib are emitted by CC with pp 1.0 

  /* Contract checks */
  ESL_DASSERT1(( sxd->type == p7S_DECODING )); 
  ESL_DASSERT1(( tr->pp    != NULL ));

  /* <last_ib> addresses a boundary case; i>last_ib 
   * must be generated by CC with postprob = 1.0.
   */
  for (last_ib = sm->L; last_ib >=1; last_ib--)
    if (sm->n[last_ib]) break;

  for (z = 0; z < tr->N; z++)
    {
      if (! tr->i[z]) { tr->pp[z] = 0.0; continue; }  // only storing pp's for emissions. Only M,I, or NN/CC/JJ have i=1..L.

      k = tr->k[z];         // will be 0 for s={NJC}; 1..M for s={MID}
      while (i < tr->i[z])  // this happens to also accommodate the case that a trace somehow starts with i>1
	{                   // i is < L here, so n[i+1] is ok below
	  dpc += p7S_NSCELLS * sm->n[i]; 
	  if (sm->n[i] || sm->n[i+1]) xc += p7S_NXCELLS; // special rows also store ia-1, hence the move on sm->n[i+1] too
	  i++;
	}
      /*  ok, now: 
       *   dpc is either on first supercell of row i (if i is stored, n[i]>0); 
       *       or it's waiting on next stored row;
       *       or it's just off area (and can't be dereferenced) if there are no more stored rows.
       *       (If it's off area, we'll get v=0 == n[i] in logic below, and avoid dereference.)
       *   xc  is on special rows of row i (if i specials stored, including ia-1 start-segment edge cases);
       *       or it's waiting on the next stored special row;
       *       or it's just off area and i>last_ib test will avoid dereferencing it.
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
p7_sparsemx_CountTrace(const P7_TRACE *tr, P7_SPARSEMX *sxd)
{
  static int sttype[p7T_NSTATETYPES] = { -1, p7S_ML, p7S_MG, p7S_IL, p7S_IG, p7S_DL, p7S_DG, -1, p7S_N, p7S_B, p7S_L, p7S_G, p7S_E, p7S_C, p7S_J, -1 }; /* sttype[] translates trace idx to DP matrix idx*/
  const P7_SPARSEMASK *sm  = sxd->sm;
  float               *dpc = sxd->dp;   /* ptr that steps thru stored main supercells i,k */
  float               *xc  = sxd->xmx;	/* ptr that steps thru stored special rows, including ia-1 seg edges */
  int   i   = 0;		/* row index (seq position) */
  int   y   = 0;		/* index in sparse cell list on a row */
  int   z;			/* index in trace positions */

  if (sxd->type == p7S_UNSET) sxd->type = p7S_DECODING; /* first time */
  
  for (z = 1; z < tr->N-1; z++)	/* z=0 is S; z=N-1 is T; neither is represented in DP matrix, so avoid them */
    {
      if (tr->i[z])		/* we'll emit i, so we need to update our matrix ptrs */
	{
	  while (y < sm->n[i]) { y++; dpc += p7S_NSCELLS; }  /* skip remainder of prev sparse row */
	  if (sm->n[i] || sm->n[i+1])  xc += p7S_NXCELLS;    /* increment specials. note, i+1 is safe here; i<L at this point, because tr->i[z] <= L and i < tr->i[z] */
	  i = tr->i[z];
	  y = 0;
	}

      if (tr->k[z])		/* if it's a main state, {MID}{LG} */
	{
	  while (y < sm->n[i] && sm->k[i][y] <  tr->k[z]) { y++; dpc += p7S_NSCELLS; } /* try to find sparse cell for i,k. note, if row isn't stored at all (n[i]==0) nothing happens here, because y==n[i] */
	  if    (y < sm->n[i] && sm->k[i][y] == tr->k[z]) dpc[sttype[(int) tr->st[z]]] += 1.0;  /* did we find it? then increment +1 */
	  else if ( tr->st[z] != p7T_DG )       	  ESL_EXCEPTION(eslEINCONCEIVABLE, "failed to find i,k supercell, and this cannot be a DG wing retraction");
	    /* Because of wing-retraction, it's possible to have DGk states in a trace that can't be stored in sparse cells,
	     * because they were implied by a wing-retracted entry/exit.
	     * for y == sm->n[i], then it must be a DG on a MGk->Dk+1...E wing-retracted exit. 
	     * for k[i][z] > tr->k[z], then it must be a DG on a G->...Dk-1->Mk wing-retracted entry.
	     * but if it's not a DG, fail.
	     */
	}
      else 			/* if it's a special state, {ENJBLGC} */
	{
	  if (sm->n[i] || (i < sm->L && sm->n[i+1])) {
	    xc[sttype[(int) tr->st[z]]] += 1.0;
	    if (tr->st[z] == p7T_C && tr->i[z]) xc[p7S_CC] += 1.0;
	    if (tr->st[z] == p7T_J && tr->i[z]) xc[p7S_JJ] += 1.0;
	  } else {   /* else, make sure it's something that's allowed to be implicit */
	    if ( ! (i == 0     && tr->st[z] == p7T_N) &&    /* S->N... on row 0 can be implicit */
		 ! (tr->i[z]))                        	    /* NN/CC/JJ emissions are implied on unstored rows; if (tr->i[z]) is used to check for NN/CC/JJ  */
	      ESL_EXCEPTION(eslEINCONCEIVABLE, "failed to find stored xmx[i] specials for this trace position");
	  }
	}
    }
  return eslOK;
}




int
p7_sparsemx_ExpectedDomains(const P7_SPARSEMX *sxd, int iae, int ibe, float *ret_ndom_expected)
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
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_SPARSEMX main state code %d\n", type);
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
p7_sparsemx_DumpWindow(FILE *ofp, const P7_SPARSEMX *sx, int i1, int i2, int k1, int k2)
{
  const P7_SPARSEMASK *sm  = sx->sm;
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

      fprintf(ofp, "%3d MG ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
	if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_MG));
	else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d IL ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k)  z++;
	if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_IL));
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

      fprintf(ofp, "%3d DL ", i);
      for (z = 0, k = k1; k <= k2; k++) { 
	while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
	if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*p7R_NSCELLS + p7R_DL));
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
p7_sparsemx_Copy2Reference(const P7_SPARSEMX *sx, P7_REFMX *rx)
{
  const P7_SPARSEMASK *sm = sx->sm;
  int             M = sm->M;
  int             L = sm->L;
  int             W = (M+1)*p7R_NSCELLS + p7R_NXCELLS; /* total width of a reference DP row, in cells */
  float           vimp;
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
  for (g = 1; g <= sm->S; g++)
    {
      ia = sm->seg[g].ia;
      ib = sm->seg[g].ib;

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

int
p7_sparsemx_Compare(const P7_SPARSEMX *sx1, const P7_SPARSEMX *sx2, float tol)
{
  char msg[] = "P7_SPARSEMX comparison failed";

  if ( sx1->type != sx2->type)                            ESL_FAIL(eslFAIL, NULL, msg);                          
  if ( p7_sparsemask_Compare(sx1->sm, sx2->sm) != eslOK)  ESL_FAIL(eslFAIL, NULL, msg);
  if ( esl_vec_FCompare(sx1->dp,  sx2->dp,   p7S_NSCELLS*sx1->sm->ncells,            tol) != eslOK)  ESL_FAIL(eslFAIL, NULL, msg);
  if ( esl_vec_FCompare(sx1->xmx, sx2->xmx,  p7S_NXCELLS*(sx1->sm->nrow+sx1->sm->S), tol) != eslOK)  ESL_FAIL(eslFAIL, NULL, msg);
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
 *            In general, this is only going to succeed if <sx> was
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
  char                 msg[] = "comparison of P7_SPARSEMX to P7_REFMX failed";
  const P7_SPARSEMASK *sm    = sx->sm;
  const float   *dpc, *dpc2;
  const float   *xc,  *xc2;
  int            g,i,s,z;
  int            ia,ib;
  
  if (sx->type != rx->type) ESL_FAIL(eslFAIL, NULL, msg);
  if (sm->M    != rx->M)    ESL_FAIL(eslFAIL, NULL, msg);
  if (sm->L    != rx->L)    ESL_FAIL(eslFAIL, NULL, msg);

  /* First traversal way: sm->dp[] is just one big array */
  dpc = sx->dp;
  xc  = sx->xmx;
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1])         /* ia-1 specials at a segment start */
	{
	  for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	    if (esl_FCompareAbs(*xc, P7R_XMX(rx,i-1,s), tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
	}

      for (z = 0; z < sm->n[i]; z++)       /* sparse cells */
	{
	  for (s = 0; s < p7S_NSCELLS; dpc++, s++) 
	    if (esl_FCompareAbs(*dpc, P7R_MX(rx,i,sm->k[i][z],s), tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
	}
  
      if (sm->n[i])       /* specials */
	{
	  for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	    if (esl_FCompareAbs(*xc, P7R_XMX(rx,i,s), tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
	}
    }


  /* Second way: "segments" */
  dpc2 = sx->dp;
  xc2  = sx->xmx;
  for (g = 1; g <= sm->S; g++)
    {
      ia = sm->seg[g].ia;
      ib = sm->seg[g].ib;

      for (s = 0; s < p7S_NXCELLS; xc2++, s++)        /* ia-1 specials at segment start */
	if (esl_FCompareAbs(*xc2, P7R_XMX(rx,ia-1,s), tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);

      for (i = ia; i <= ib; i++) 
	{
	  for (z = 0; z < sm->n[i]; z++)      	  /* sparse main cells */
	    {
	      for (s = 0; s < p7S_NSCELLS; dpc2++, s++) 
		if (esl_FCompareAbs(*dpc2, P7R_MX(rx,i,sm->k[i][z],s), tol) == eslFAIL)  ESL_FAIL(eslFAIL, NULL, msg);
	    }
	  for (s = 0; s < p7S_NXCELLS; xc2++, s++)  	  /* specials */
	    if (esl_FCompareAbs(*xc2, P7R_XMX(rx,i,s), tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
	}
    }
  
  /* Both ways must reach the same end */
  if (dpc != dpc2) ESL_FAIL(eslFAIL, NULL, msg);
  if (xc  != xc2)  ESL_FAIL(eslFAIL, NULL, msg);
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
  char                 msg[] = "failed comparison of P7_SPARSEMX to upper bound in P7_REFMX"; 
  const P7_SPARSEMASK *sm    = sx->sm;
  const float         *dpc   = sx->dp;
  const float         *xc    = sx->xmx;
  int   i,s,z;

  if (sx->type != rx->type) ESL_FAIL(eslFAIL, NULL, msg);
  if (sm->M    != rx->M)    ESL_FAIL(eslFAIL, NULL, msg);
  if (sm->L    != rx->L)    ESL_FAIL(eslFAIL, NULL, msg);
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {    /* ia-1 specials at a segment start */
	for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	  if (*xc > P7R_XMX(rx,i-1,s)+tol) 
	    ESL_FAIL(eslFAIL, NULL, msg);
      }
      for (z = 0; z < sm->n[i]; z++)    /* sparse cells */
	for (s = 0; s < p7S_NSCELLS; dpc++, s++) 
	  if (*dpc > P7R_MX(rx,i,sm->k[i][z],s)+tol) 
	    ESL_FAIL(eslFAIL, NULL, msg);
  
      if (sm->n[i]) {    /* specials */
	for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	  if (*xc > P7R_XMX(rx,i,s)+tol) 
	    ESL_FAIL(eslFAIL, NULL, msg);
      }
    }
  return eslOK;
}

/* Function:  p7_sparsemx_CompareDecoding()
 * Synopsis:  Compare exact, approximate posterior decoding matrices.
 *
 * Purpose:   Compare exact sparse decoding matrix <sxe> (calculated by 
 *            <p7_SparseDecoding()> to an approximate one <sxa> 
 *            (calculated by sampling lots of stochastic traces).
 *            Make sure that no nonzero value in <sxe> differs by
 *            more than absolute difference <tol> in <sxa>, and make
 *            sure that an exact zero value in <sxe> is also zero
 *            in <sxa>. Return <eslOK> if these tests succeed; <eslFAIL>
 *            if they do not.
 *            
 *            This comparison is used in the main unit test of
 *            posterior decoding. See
 *            <sparse_fwdback.c::utest_approx_decoding()>.
 *            
 * Args:      sxe  - exact posterior decoding matrix, from p7_SparseDecoding
 *            sxa  - approximate decoding mx, from stochastic trace ensemble
 *            tol  - absolute difference to tolerate per cell value
 *
 * Returns:   <eslOK> if all comparisons pass. <eslFAIL> if not.
 */
int
p7_sparsemx_CompareDecoding(const P7_SPARSEMX *sxe, const P7_SPARSEMX *sxa, float tol)
{
  char                 msg[] = "failed comparison of exact and sampled sparse decoding matrices";
  const P7_SPARSEMASK *sm  = sxe->sm;
  const float         *dpe = sxe->dp;
  const float         *dpa = sxa->dp;
  const float         *xce = sxe->xmx;
  const float         *xca = sxa->xmx;
  int   i,s,z;
  float diff;
  float max = 0.;

  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {    /* ia-1 specials at a segment start */
	for (s = 0; s < p7S_NXCELLS; xce++, xca++, s++) 
	  {
	    diff = *xce-*xca;
	    max  = ESL_MAX(fabs(diff), max);
	    //printf(" diff %9.4f  i=%d %s\n", diff, i-1, p7_sparsemx_DecodeSpecial(s));
	    if ( (*xce == 0. && *xca != 0.) || fabs(*xce-*xca) > tol) 
	      ESL_FAIL(eslFAIL, NULL, msg);
	  }	    
      }
      for (z = 0; z < sm->n[i]; z++)    /* sparse cells */
	for (s = 0; s < p7S_NSCELLS; dpe++, dpa++, s++)
	  {
	    diff = *dpe-*dpa;
	    max  = ESL_MAX(fabs(diff),max);
	    //printf(" diff %9.4f  i=%d k=%d %s  %9.4f %9.4f\n", diff, i, sm->k[i][z], p7_sparsemx_DecodeState(s), *dpe, *dpa);
	    if ( (*dpe == 0. && *dpa != 0.) || fabs(*dpe-*dpa) > tol) 
	      ESL_FAIL(eslFAIL, NULL, msg);
	  }
      if (sm->n[i]) {		        /* specials */
	for (s = 0; s < p7S_NSCELLS; xce++, xca++, s++)
	  {
	    diff = *xce-*xca;
	    max  = ESL_MAX(fabs(diff),max);
	    //printf(" diff %9.4f  i=%d %s\n", diff, i, p7_sparsemx_DecodeSpecial(s));
	    if ( (*xce == 0. && *xca != 0.) || fabs(*xce-*xca) > tol) 
	      ESL_FAIL(eslFAIL, NULL, msg);
	  }
      }
    }
  //printf("max diff %f\n", max);
  return eslOK;
}

/* Function:  p7_sparsemx_PlotDomainInference()
 * Synopsis:  Plot graph of domain inference data, in xmgrace XY format
 *
 * Purpose:   Given posterior decoding matrix <sxd>, plot a 2D graph of
 *            various quantities that could be useful for inferring 
 *            sequence coordinates of a domain location. Output is
 *            in XMGRACE xy format, and is sent to open writable
 *            stream <ofp>.
 *            
 *            Restrict the plot to sequence coords <ia..ib>, which can
 *            range from <1..sxd->L>. To get the full sequence, pass
 *            <ia=1>, <ib=sxd->L>.
 *            
 *            At least four datasets are plotted. Set 0 is
 *            P(homology): the posterior probability that this residue
 *            is 'in' the model (as opposed to being emitted by N,C,
 *            or J states).  Set 1 is P(B), the posterior probability
 *            of the BEGIN state. Set 2 is P(E), the posterior
 *            probability of the END state. Set 3 plots the maximum
 *            posterior probability of any emitting model state (M or
 *            I, glocal or local).
 *            
 *            Optionally, caller may also provide an indexed trace
 *            <tr>, showing where domains have been defined (as
 *            horizontal lines above the plot); each domain is an
 *            xmgrace dataset (so, sets 4..4+ndom-1 are these
 *            horizontal lines).
 *
 * Xref:      p7_refmx_PlotDomainInference() is the same thing,
 *            for P7_REFMX, with sets in the same order.
 */
int
p7_sparsemx_PlotDomainInference(FILE *ofp, const P7_SPARSEMX *sxd, int ia, int ib, const P7_TRACE *tr)
{
  const P7_SPARSEMASK *sm = sxd->sm;
  const float  tr_height  = 1.2;
  const float *xc;
  const float *dpc;
  float  val;
  int    i,z,d;


  /* Set #0: P(homology), via 1 - (NN+CC+JJ) */
  xc = sxd->xmx;
  for (i = 0; i <= ib; i++)
    {
      if ( sm->n[i] || (i < sm->L && sm->n[i+1]) ) /* either i is a stored row in segment, or an ia-1 special row before segment start */
	{
	  val = 1.0 - (xc[p7S_N] + xc[p7S_JJ] + xc[p7S_CC]); /* this might get calculated for i=0, and be wrong, because N(0) is not what we want. but we won't use <val> on i=0; see below */
	  xc += p7S_NXCELLS;
	}
      else val = 0.0;

      if (i >= ia) fprintf(ofp, "%-6d %.5f\n", i, val);
    }
  fprintf(ofp, "&\n");
  
  /* Set #1: P(B), begin */
  xc = sxd->xmx;
  for (i = 0; i <= ib; i++)
    {
      if ( sm->n[i] || (i < sm->L && sm->n[i+1]) ) /* either i is a stored row in segment, or an ia-1 special row before segment start */
	{
	  val = xc[p7S_B];
	  xc += p7S_NXCELLS;
	}
      else val = 0.0;
      if (i >= ia) fprintf(ofp, "%-6d %.5f\n", i, val);
    }
  fprintf(ofp, "&\n");
  
  /* Set #2: P(E), end */
  xc = sxd->xmx;
  for (i = 0; i <= ib; i++)
    {
      if ( sm->n[i] || (i < sm->L && sm->n[i+1]) ) /* either i is a stored row in segment, or an ia-1 special row before segment start */
	{
	  val = xc[p7S_E];
	  xc += p7S_NXCELLS;
	}
      else val = 0.0;
      if (i >= ia) fprintf(ofp, "%-6d %.5f\n", i, val);
    }
  fprintf(ofp, "&\n");

  /* Set #3: max P(i,k,s) for s=M,I emitting states (both G and L) */
  for (dpc = sxd->dp, i = 1; i < ia; i++) 
    dpc += sm->n[i] * p7S_NSCELLS;
  for (; i <= ib; i++)
    {
      for (val = 0., z = 0; z < sm->n[i]; z++)
	{
	  val = ESL_MAX(val, 
			ESL_MAX( ESL_MAX(dpc[p7S_ML], dpc[p7S_MG]),
				 ESL_MAX(dpc[p7S_IL], dpc[p7S_IG])));
	  dpc += p7S_NSCELLS;
	}
      fprintf(ofp, "%-6d %.5f\n", i, val);
    }
  fprintf(ofp, "&\n");

  if (tr && tr->ndom)  /* Remaining sets are horiz lines, representing individual domains appearing in the optional <tr> */
    {
      for (d = 0; d < tr->ndom; d++)
	{
	  fprintf(ofp, "%-6d %.5f\n", tr->sqfrom[d], tr_height);
	  fprintf(ofp, "%-6d %.5f\n", tr->sqto[d],   tr_height);
	  fprintf(ofp, "&\n");
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
validate_dimensions(const P7_SPARSEMX *sx, char *errbuf)
{
  const P7_SPARSEMASK *sm  = sx->sm;
  int   g      = 0;
  int   r      = 0;
  int   ncells = 0;
  int   i;

  if ( sm->M <= 0)              ESL_FAIL(eslFAIL, errbuf, "nonpositive M");
  if ( sm->L <= 0)              ESL_FAIL(eslFAIL, errbuf, "nonpositive L");
  if ( sm->Q <  P7_NVF(sm->M))  ESL_FAIL(eslFAIL, errbuf, "insufficient Q");

  for (r=0, g=0, i = 1; i <= sm->L; i++) {
    if (sm->n[i] && !sm->n[i-1]) g++; /* segment count */
    if (sm->n[i])                r++; /* sparse row count */
    ncells += sm->n[i];
  }
  if (g      != sm->S)          ESL_FAIL(eslFAIL, errbuf, "S (nseg) is wrong");
  if (r      != sm->nrow)       ESL_FAIL(eslFAIL, errbuf, "nrow is wrong");
  if (ncells != sm->ncells)     ESL_FAIL(eslFAIL, errbuf, "ncells is wrong");

  if (sm->L+1    > sm->ralloc)  ESL_FAIL(eslFAIL, errbuf, "k[] row allocation too small");
  if (sm->ncells > sm->kalloc)  ESL_FAIL(eslFAIL, errbuf, "kmem[] cell allocation too small");
  if (sm->S+2    > sm->salloc)  ESL_FAIL(eslFAIL, errbuf, "seg[] segment allocation too small");
  return eslOK;
}

static int
validate_no_nan(const P7_SPARSEMX *sx, char *errbuf)
{
  const P7_SPARSEMASK *sm  = sx->sm;
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
validate_fwdvit(const P7_SPARSEMX *sx, char *errbuf)
{
  const P7_SPARSEMASK *sm  = sx->sm;
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
validate_backward(const P7_SPARSEMX *sx, char *errbuf)
{
  const P7_SPARSEMASK *sm     = sx->sm;
  float         *dpc    = sx->dp  + (sm->ncells-1)*p7S_NSCELLS;		// last supercell in dp  
  float         *xc     = sx->xmx + (sm->nrow + sm->S - 1)*p7S_NXCELLS; // last supercell in xmx 
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
  if (val < 0. || val > 1.+tol) return FALSE; 
  return TRUE;
}

static int
validate_decoding(const P7_SPARSEMX *sx, char *errbuf)
{
  const P7_SPARSEMASK *sm  = sx->sm;
  float         *dpc = sx->dp;
  float         *xc  = sx->xmx;
  int            i,z,s;
  int            last_n;
  float          tol = 0.001;

  /* Check special cases prohibited in the first ia-1 presegment specials: */
  /* Note that the N on the first ia-1 row is *not* necessarily 1, because ia-1 is not necessarily row 0 */
  if ( esl_FCompareAbs(xc[p7S_J],  0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first J not 0");
  if ( esl_FCompareAbs(xc[p7S_C],  0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first C not 0");
  if ( esl_FCompareAbs(xc[p7S_JJ], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first JJ not 0");
  if ( esl_FCompareAbs(xc[p7S_CC], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first CC not 0");

  /* Sweep, checking for (the most easily spotchecked) prohibited values (must be 0.0's) */
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {       /* ia-1 specials at a segment start */
	if ( esl_FCompareAbs(xc[p7S_E], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "E seg start for ia=%d not 0", i);
	if ( xc[p7S_J]+tol < xc[p7S_JJ])                    ESL_FAIL(eslFAIL, errbuf, "JJ>J at seg start for ia=%d ", i);
	if ( xc[p7S_C]+tol < xc[p7S_CC])                    ESL_FAIL(eslFAIL, errbuf, "CC>C at seg start for ia=%d ", i);
	xc += p7S_NXCELLS;
      }
      for (z = 0; z < sm->n[i]; z++)       /* sparse main cells */
	{
	  /* if k-1 supercell doesn't exist, can't reach DL. But all DGk are reachable, because of wing-retracted entry/exit */
	  if ((z == 0 || sm->k[i][z] != sm->k[i][z-1]+1) && esl_FCompareAbs(dpc[p7S_DL], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first DL on i=%d not 0", i);
	  if (   sm->k[i][z] == sm->M                    && esl_FCompareAbs(dpc[p7S_IL], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "IL on i=%d,M not 0", i);
	  if (   sm->k[i][z] == sm->M                    && esl_FCompareAbs(dpc[p7S_IG], 0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "IG on i=%d,M not 0", i);
	  dpc += p7S_NSCELLS;
	  /* there are other conditions where I(i,k) values must be zero but this is more tedious to check */
	}
      if (sm->n[i]) {
	if ( xc[p7S_J]+tol < xc[p7S_JJ])                    ESL_FAIL(eslFAIL, errbuf, "JJ>J at i=%d ", i);
	if ( xc[p7S_C]+tol < xc[p7S_CC])                    ESL_FAIL(eslFAIL, errbuf, "CC>C at i=%d ", i);
	xc += p7S_NXCELLS;
      }
    }

  /* Backwards sweep, looking only at ib end rows. */
  dpc    = sx->dp  + (sm->ncells-1)*p7S_NSCELLS;		 // last supercell in dp  
  xc     = sx->xmx + (sm->nrow + sm->S - 1)*p7S_NXCELLS; // last supercell in xmx 
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
	  if (! is_prob(*dpc, tol)) ESL_FAIL(eslFAIL, errbuf, "bad decode prob %f at i=%d,k=%d,%s", *dpc, i,sm->k[i][z], p7_sparsemx_DecodeState(s));
      if (sm->n[i]) {                      /* specials on sparse row */
	for (s = 0; s < p7S_NXCELLS; xc++, s++) 
	  if (! is_prob(*xc, tol))  ESL_FAIL(eslFAIL, errbuf, "bad decode prob %f at i=%d,%s", *xc, i, p7_sparsemx_DecodeSpecial(s)); 
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
p7_sparsemx_Validate(const P7_SPARSEMX *sx, char *errbuf)
{
  int status;

  if (errbuf) errbuf[0] = '\0';

  if ( (status = validate_dimensions(sx, errbuf)) != eslOK) return status;
  if ( (status = validate_no_nan    (sx, errbuf)) != eslOK) return status;

  switch (sx->type) {
  case p7S_UNSET:      ESL_FAIL(eslFAIL, errbuf, "validating an unset sparse DP matrix? probably not what you meant");
  case p7S_FORWARD:    if ( (status = validate_fwdvit  (sx, errbuf)) != eslOK) return status; break;
  case p7S_BACKWARD:   if ( (status = validate_backward(sx, errbuf)) != eslOK) return status; break;
  case p7S_DECODING:   if ( (status = validate_decoding(sx, errbuf)) != eslOK) return status; break;
  case p7S_VITERBI:    if ( (status = validate_fwdvit  (sx, errbuf)) != eslOK) return status; break;
  default:             ESL_FAIL(eslFAIL, errbuf, "no such sparse DP matrix type %d", sx->type);
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
