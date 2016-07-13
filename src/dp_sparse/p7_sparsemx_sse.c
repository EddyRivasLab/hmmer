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
#ifdef HAVE_SSE2
int
p7_sparsemask_StartRow_sse(P7_SPARSEMASK *sm, int i)
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
#endif


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
#ifdef HAVE_SSE2
 int
p7_sparsemask_Add_sse(P7_SPARSEMASK *sm, int q, int r)
{
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
}
#endif

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
#ifdef HAVE_SSE2
int
p7_sparsemask_FinishRow_sse(P7_SPARSEMASK *sm)
{
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
}
#endif

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
}
#endif
/*----------------- end, P7_SPARSEMASK API ----------------------*/



/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
