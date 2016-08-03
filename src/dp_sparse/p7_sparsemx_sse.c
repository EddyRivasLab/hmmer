// SSE versions of functions to manipulate P7_SPARSEMASK data structures

/* 
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
p7_sparsemask_Create_sse(int M, int L)
{
#ifdef HAVE_SSE2  
  P7_SPARSEMASK *sm             = NULL;
  int            default_salloc = 8;
  int64_t        default_kalloc = 4096;
  int            i,r;
  int            status;

  ESL_ALLOC(sm, sizeof(P7_SPARSEMASK));
  sm->L      = L;
  sm->M      = M;
  sm->simd = SSE; 
  sm->Q      = P7_NVF(M);  // approx M/4, for striped vectors of four floats
    
  sm->seg    = NULL;
  sm->k      = NULL;
  sm->n      = NULL;
  sm->kmem   = NULL;

  sm->S       = 0;
  sm->nrow    = 0;
  sm->ncells  = 0;
  sm->last_i  = L+1;       // sentinel to assure StartRow() is called in reverse L..1 order 
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
  sm->salloc   = default_salloc;

 

  ESL_ALLOC(sm->seg,  sm->salloc * sizeof(struct p7_sparsemask_seg_s)); // salloc is the actual allocation, inclusive of +2 for sentinels
  ESL_ALLOC(sm->k,    sm->ralloc * sizeof(int *));
  ESL_ALLOC(sm->n,    sm->ralloc * sizeof(int));
  ESL_ALLOC(sm->kmem, sm->kalloc * sizeof(int));

  sm->k[0]   = NULL;    // always. 
  for (i = 0; i <= L; i++)  // n[0] will always be 0; n[i=1..L] initialized to 0, then count as cells are added 
    sm->n[i] = 0;
  sm->n_krealloc = 0;
  sm->n_rrealloc = 0;
  sm->n_srealloc = 0;
  return sm;

 ERROR:
  p7_sparsemask_Destroy(sm);
#endif //HAVE_SSE2
#ifndef HAVE_SSE2
return NULL;
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
p7_sparsemask_Reinit_sse(P7_SPARSEMASK *sm, int M, int L)
{
 // printf("Calling p7_sparsemask_Reinit with M= %d, L = %d\n", M, L);
  /* Implement one version of this function that handles all possible vector ISAs because we only 
  call it from one place */
#ifdef HAVE_SSE2
  int i,r;
  int status;

  sm->L  = L;
  sm->M  = M; 
  sm->Q  = P7_NVF(M);
    
  /* seg[], kmem stay at their previous salloc, kalloc
   * but do we need to reallocate rows for k[] and n[]? 
   */
  if (sm->ralloc < L+1) {
   // printf("reallocating base\n");
    ESL_REALLOC(sm->k, sizeof(int *) * (L+1));
    ESL_REALLOC(sm->n, sizeof(int)   * (L+1));
    sm->ralloc = L+1;
    sm->n_rrealloc++;
  }

  sm->S       = 0;
  sm->nrow    = 0;
  sm->ncells  = 0;
  sm->last_i  = sm->L+1;
  for (r = 0; r < p7_VNF; r++) 
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
#endif //HAVE_SSE2
#ifndef HAVE_SSE2
return eslENORESULT;
#endif  
}

/* Function:  p7_sparsemask_Sizeof()
 * Synopsis:  Returns current allocated size of a <P7_SPARSEMASK>, in bytes.
 */
size_t
p7_sparsemask_Sizeof_sse(const P7_SPARSEMASK *sm)
{
  size_t n = sizeof(P7_SPARSEMASK);   
#ifdef HAVE_SSE2  
  n += sm->salloc * sizeof(p7_sparsemask_seg_s); // <seg>
  n += sm->ralloc * sizeof(int *);                      // <k>                   
  n += sm->ralloc * sizeof(int);                        // <n>                   
  n += sm->kalloc * sizeof(int);                        // <kmem>                
  return n;
 #endif
}

/* Function:  p7_sparsemask_MinSizeof()
 * Synopsis:  Returns minimum required size of a <P7_SPARSEMASK>, in bytes.
 */
size_t
p7_sparsemask_MinSizeof_sse(const P7_SPARSEMASK *sm)
{
  size_t n = sizeof(P7_SPARSEMASK);
#ifdef HAVE_SSE2 
  n += (sm->S+2)  * sizeof(struct p7_sparsemask_seg_s);  // <seg>; includes sentinels at 0,S+1
  n += (sm->L+1)  * sizeof(int *);                       // <k>
  n += (sm->L+1)  * sizeof(int);                         // <n>
  n += sm->ncells * sizeof(int);                         // <kmem>
#endif 
  return n;
}

/* Function:  p7_sparsemask_Destroy()
 * Synopsis:  Destroy a <P7_SPARSEMASK>.
 */
void
p7_sparsemask_Destroy_sse(P7_SPARSEMASK *sm)
{
  if (sm) {
#ifdef HAVE_SSE2    
    if (sm->seg)  free(sm->seg);
    if (sm->k)    free(sm->k);
    if (sm->n)    free(sm->n);
    if (sm->kmem) free(sm->kmem);
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

// Separate versions of these functions for each ISA because they're called from 
// ISA-specific functions

int
p7_sparsemask_StartRow_sse(P7_SPARSEMASK *sm, int i)
{
#ifdef HAVE_SSE2  
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
  #endif // HAVE_SSE2
#ifndef HAVE_SSE2
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
p7_sparsemask_Add_sse(P7_SPARSEMASK *sm, int q, int r)
{
#ifdef HAVE_SSE2  
  int     k = r*sm->Q+q+1;

  //printf("Sparsemask Adding i=%d q=%d r=%d k=%d M=%d\n", sm->last_i, q, r, k, sm->M);

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
  #endif // HAVE_SSE2
#ifndef HAVE_SSE2
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
p7_sparsemask_FinishRow_sse(P7_SPARSEMASK *sm)
{
#ifdef HAVE_SSE2  
  int *p;
  int  r;
//  printf("sm->ncells = %li\n", sm->ncells);
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
#endif // HAVE_SSE2
#ifndef HAVE_SSE2
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
 #ifdef HAVE_SSE2
int
p7_sparsemask_Finish_sse(P7_SPARSEMASK *sm)
{
  int i,r;
  int s;
  int status;
//printf("calling p7_sparsemask_Finish, sm->ncells = %li\n", sm->ncells);

  int *p;
  esl_vec_IReverse(sm->kmem, sm->kmem, sm->ncells);

  /* Set the k[] pointers; count <S> and <nrow> */
  p = sm->kmem;
  sm->S = sm->nrow = 0;

  for (i = 1; i <= sm->L; i++){
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

   sm->last_i = -1;
  for (r = 0; r < p7_VNF; r++) 
    sm->last_k[r] = -1;

  return eslOK;

 ERROR:
  return eslEMEM;
  #endif // HAVE_SSE2
#ifndef HAVE_SSE2
  return eslENORESULT;
  #endif
}
/*----------------- end, P7_SPARSEMASK API ----------------------*/

int
p7_sparsemask_Dump_sse(FILE *ofp, P7_SPARSEMASK *sm)
{
#ifdef  HAVE_SSE2
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
#endif // HAVE_SSE2
#ifndef HAVE_SSE2
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
p7_sparsemask_Compare_sse(const P7_SPARSEMASK *sm1, const P7_SPARSEMASK *sm2)
{
  #ifdef HAVE_SSE2
  char msg[] = "P7_SPARSEMASK comparison failed";
  int  i;
  int  s;
  if(sm2->simd != SSE){
    ESL_FAIL(eslFAIL, NULL, "Can't compare sparsemasks generated for different SIMD instruction sets");
  }

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
  #endif //HAVE_SSE2
#ifndef HAVE_SSE2
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
p7_sparsemask_Validate_sse(const P7_SPARSEMASK *sm, char *errbuf)
{
#ifdef HAVE_SSE2  
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
  #endif //HAVE_SSE2
#ifndef HAVE_SSE2
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
p7_sparsemask_SetFromTrace_sse(P7_SPARSEMASK *sm, ESL_RANDOMNESS *rng, const P7_TRACE *tr)
{
#ifdef HAVE_SSE2
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
        if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) return status;

    while (p7_trace_IsD(tr->st[z])) {
      k = tr->k[z]; 
      if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) return status;
      z--;
    }

    k = tr->k[z];
    if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) return status;
    
    for (k = k-1; k >= 1; k--)
      if (rng && esl_random(rng) < cellprob)
        if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) return status;
  }
      else
  {
    if (rng && esl_random(rng) < rowprob)
      for (k = sm->M; k >= 1; k--)
        if (rng && esl_random(rng) < cellprob)
    if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) return status; /* append to k[i] list, increment n[i] count, reallocating as needed; doesn't deal w/ segments (nrow,nseg,i[]) */
  }

      if ((status = p7_sparsemask_FinishRow(sm)) != eslOK) return status;
    }
  if ( (status = p7_sparsemask_Finish(sm)) != eslOK) return status;

  return eslOK;
   #endif //HAVE_SSE2
#ifndef HAVE_SSE2
  return eslENORESULT;
  #endif
}

/*****************************************************************
 * 7. Validation of a P7_SPARSEMX
 *****************************************************************/
/* also a debugging tool, but in its own section because it's
 * fiddly and complicated
 */
static int
validate_dimensions_sse(const P7_SPARSEMX *sx, char *errbuf)
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
validate_no_nan_sse(const P7_SPARSEMX *sx, char *errbuf)
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
validate_fwdvit_sse(const P7_SPARSEMX *sx, char *errbuf)
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
validate_backward_sse(const P7_SPARSEMX *sx, char *errbuf)
{
  const P7_SPARSEMASK *sm     = sx->sm;
  float         *dpc    = sx->dp  + (sm->ncells-1)*p7S_NSCELLS;   // last supercell in dp  
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
validate_decoding_sse(const P7_SPARSEMX *sx, char *errbuf)
{
  const P7_SPARSEMASK *sm  = sx->sm;
  float         *dpc = sx->dp;
  float         *xc  = sx->xmx;
  int            i,z,s;
  int            last_n;
  float          tol = 0.001;

  /* Check special cases prohibited in the first ia-1 presegment specials: */
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
  dpc    = sx->dp  + (sm->ncells-1)*p7S_NSCELLS;     // last supercell in dp  
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
p7_sparsemx_Validate_sse(const P7_SPARSEMX *sx, char *errbuf)
{
  int status;

  if (errbuf) errbuf[0] = '\0';

  if ( (status = validate_dimensions_sse(sx, errbuf)) != eslOK) return status;
  if ( (status = validate_no_nan_sse    (sx, errbuf)) != eslOK) return status;

  switch (sx->type) {
  case p7S_UNSET:      ESL_FAIL(eslFAIL, errbuf, "validating an unset sparse DP matrix? probably not what you meant");
  case p7S_FORWARD:    if ( (status = validate_fwdvit_sse  (sx, errbuf)) != eslOK) return status; break;
  case p7S_BACKWARD:   if ( (status = validate_backward_sse(sx, errbuf)) != eslOK) return status; break;
  case p7S_DECODING:   if ( (status = validate_decoding_sse(sx, errbuf)) != eslOK) return status; break;
  case p7S_VITERBI:    if ( (status = validate_fwdvit_sse  (sx, errbuf)) != eslOK) return status; break;
  default:             ESL_FAIL(eslFAIL, errbuf, "no such sparse DP matrix type %d", sx->type);
  }
  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
