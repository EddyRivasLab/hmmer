/* H4_ANCHORSET: an array of (i0,k0) anchor points that define domain locations
 *
 * Contents:
 *   1. H4_ANCHORSET object
 *   2. Debugging, development tools
 */
#include <h4_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "h4_anchorset.h"
#include "h4_path.h"


/*****************************************************************
 * 1. H4_ANCHORSET object.
 *****************************************************************/

/* Function:  h4_anchorset_Create()
 * Synopsis:  Create a new <H4_ANCHORSET>
 *
 * Purpose:   Returns a pointer to a newly created <H4_ANCHORSET>
 *            object, for <D> anchors (thus defining <D> domains) for
 *            a comparison of a sequence of length <L> to a profile of
 *            length <M>.
 *
 *            <D> can be 0, to create an initially empty anchor set
 *            that will be grown later.
 *
 *            If <D> is 0, <L> and <M> can also be passed as 0, if
 *            they aren't known yet. In this case, the caller will
 *            need to call <h4_anchorset_SetSentinels()>. 
 *            
 * Returns:   pointer to the new <H4_ANCHORSET>.
 *
 * Throws:    <NULL> on allocation failure.
 */
H4_ANCHORSET *
h4_anchorset_Create(int D, int L, int M)
{
  H4_ANCHORSET *anch             = NULL;
  int32_t       default_nalloc   = 8;     // i.e. for up to D=6;  +2 sentinels
  int32_t       default_nredline = 64;    // i.e. for up to D=62; +2 sentinels
  int           status;

  ESL_DASSERT1(( ( D >= 0 && L > 0 && M > 0 ) ||
                 ( D == 0 && L == 0 && M == 0 ) ));

  ESL_ALLOC(anch, sizeof(H4_ANCHORSET));
  anch->a        = NULL;
  anch->D        = D;
  anch->nalloc   = (D+2 > default_nalloc ? D+2 : default_nalloc);
  anch->nredline = default_nredline;

  ESL_ALLOC(anch->a, sizeof(H4_ANCHOR) * (anch->nalloc));
  anch->a[0].i0   = 0;
  anch->a[0].k0   = (M == 0 ? 0 : M+1);
  anch->a[D+1].i0 = (L == 0 ? 0 : L+1);
  anch->a[D+1].k0 = 0;
  return anch;

 ERROR:
  h4_anchorset_Destroy(anch);
  return NULL;
}


/* Function:  h4_anchorset_Add()
 * Synopsis:  Add an anchor to an anchor set, left to right.
 * Incept:    SRE, Tue 05 Jan 2021
 *
 * Purpose:   Add anchor <i0>,<k0> to an anchor set <anch>.  Anchors
 *            must be added in order, left to right across the
 *            sequence, from smallest to largest <i>.
 *
 *            <anch> is reallocated if it needs to grow.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 *            <eslEINVAL> if <i0> isn't after the previous anchor.
 *            <eslERANGE> if <anch->D> gets ridiculously large
 */
int
h4_anchorset_Add(H4_ANCHORSET *anch, int i0, int k0)
{
  int status;

  if (( status = h4_anchorset_GrowFor(anch, anch->D+1)) != eslOK) return status;

  if ( i0 <= anch->a[anch->D].i0 )
    ESL_EXCEPTION(eslEINVAL, "anchors must be added in order, left to right across sequence");

  anch->D++;
  anch->a[anch->D+1].i0 = anch->a[anch->D].i0;  // move the right sentinel over one
  anch->a[anch->D+1].k0 = anch->a[anch->D].k0;
  anch->a[anch->D].i0   = i0;
  anch->a[anch->D].k0   = k0;
  return eslOK;
}


/* Function:  h4_anchorset_GetSentinels()
 * Synopsis:  Get L,M from existing sentinels
 *
 * Purpose:   Use the existing sentinels at <0,D+1> to retrieve
 *            sequence length <L> and profile length <M>. This gets
 *            used when we need to change the anchor set length
 *            and set new sentinels.
 */
int
h4_anchorset_GetSentinels(const H4_ANCHORSET *anch, int *ret_L, int *ret_M)
{
  *ret_L = anch->a[anch->D+1].i0 - 1;
  *ret_M = anch->a[0].k0 - 1;
  return eslOK;
}

/* Function:  h4_anchorset_SetSentinels()
 * Synopsis:  Initialize sentinel values in an H4_ANCHORSET
 *
 * Purpose:   Initialize sentinels <0> and <D+1> in a <H4_ANCHORSET>
 *            array <anch> defining <1..D> domains, for a (sub)sequence
 *            of length <L> and a profile of length <M>.
 *
 *            <anch> must have <anch->D> already set.
 *            
 *            <anch->a[0]> is set to <(i0,k0) = (0, M+1)>.
 *            <anch->a[D+1] is set to <(i0,k0) = (L+1, 0)>.
 *
 * Returns:   <eslOK> on success.
 */
int
h4_anchorset_SetSentinels(H4_ANCHORSET *anch, int L, int M)
{
  anch->a[0].i0   = 0;         // UP(d) sector runs i0(d-1)..i0(d)-1: so for UP(1), i0(0) must evaluate to exactly 0
  anch->a[anch->D+1].i0 = L+1; // DOWN(d) sector runs i0(d)..i0(d+1)-1; so for DOWN(D), i0(D+1) must evaluate to exactly L+1

  anch->a[0].k0   = M+1;       // DOWN(d) sector runs k0(d)..M. For access patterns that do DOWN(d-1),UP(d) on row i, must have k0(0) > M so DOWN(0) isn't evaluated;    M+1 suffices
  anch->a[anch->D+1].k0 = 0;   // UP(d) secor runs 1..k0(d)-1.  For access patterns that do DOWN(d-1),UP(d) on row i, must have k0(D+1) <= 1 so UP(D+1) isn't evaluated; 0 suffices
  return eslOK;
}


/* Function:  h4_anchorset_GrowFor()
 * Synopsis:  Reallocate H4_ANCHORSET object, if necessary
 *
 * Purpose:   Make sure that <anch> can hold an array of
 *            at least <D> anchors.
 *            
 *            Does not alter any data that are already stored
 *            in <anch>, so it's safe to resize an anchor
 *            array that we're growing incrementally (as in
 *            segmental divide and conquer MPAS algorithm).
 *            
 *            D=0 is a valid argument and may occur in normal use; it
 *            results in a no-op, because the structure is always big
 *            enough to hold zero anchors.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslERANGE> if <D> is so ridiculously large that doubling
 *            allocation size could overflow an integer.
 *            
 * Xref:      First example of a new pattern for how we 
 *            can handle reallocation/reuse strategy, 
 *            replacing _Reinit() and _Grow() interfaces.
 *            [SRE:J14/1]
 */
int
h4_anchorset_GrowFor(H4_ANCHORSET *anch, int D)
{
  int nalloc;
  int status;
  
  /* Contract checks, argument validation */
  ESL_DASSERT1(( anch->nalloc   >= 2 ));               
  ESL_DASSERT1(( anch->nredline >= 2 ));               
  if (D+2 > INT32_MAX / 2) ESL_XEXCEPTION(eslERANGE, "domain number too large, can't reallocate anchorset"); // safeguard, because we realloc by doubling

  if      (D+2 <= anch->nalloc) return eslOK;       // If we're big enough already, do nothing;
  else if (D+2 <  anch->nredline || anch->D > 0)    // If we're under the redline max, or if it looks like
    {                                               //   we're building the anchor array incrementally (because anch->D is already >0)
      nalloc = anch->nalloc;                        //   we reallocate by doubling, trying to minimize 
      while (nalloc < D+2) nalloc *= 2;             //   the need for more reallocations soon.
    }                                               // If we're over redline AND it looks like we're
  else nalloc = D+2;                                //   starting an empty object, allocate exactly.
                                                    //   Now nalloc will probably not be a multiple of two -- 
                                                    //   but the next _Reuse() call will pull it back
                                                    //   to the redline, which is.
  ESL_REALLOC(anch->a, sizeof(H4_ANCHOR) * nalloc); // Do not alter the data, or the sentinels. 
  anch->nalloc = nalloc;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  h4_anchorset_Copy()
 * Synopsis:  Copy one H4_ANCHORS object to another.
 *
 * Purpose:   Copy <src> to <dst>. <dst> is reallocated
 *            if necessary. 
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
h4_anchorset_Copy(const H4_ANCHORSET *src, H4_ANCHORSET *dst)
{
  int32_t d;
  int     status;

  if (( status = h4_anchorset_GrowFor(dst, src->D) ) != eslOK) goto ERROR;
  
  for (d = 0; d <= src->D+1; d++)       // inclusive of sentinels!
    {
      dst->a[d].i0 = src->a[d].i0;
      dst->a[d].k0 = src->a[d].k0;
    }
  dst->D    = src->D;
  return eslOK;

 ERROR:
  return status;
}




/* Function:  h4_anchorset_Reuse()
 * Synopsis:  Reuse an H4_ANCHORSET
 * Incept:    SRE, Sat 02 Jan 2021
 *
 * Purpose:   Reuse the anchor set <anch>: reinitialize the object
 *            as if it's been newly created, reusing the memory that's
 *            already allocated.
 *
 *            Sentinels are preserved. <anch> is returned for D=0,
 *            with i0[1] = L+1, k0[1] = M+1.
 *
 *            If the allocation size exceeds the redline, reallocate
 *            it down to redline.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on allocation failure. (Shouldn't happen; it can only
 *            downsize the allocation.)
 */
int
h4_anchorset_Reuse(H4_ANCHORSET *anch)
{
  int     status;

  ESL_DASSERT1(( anch->nredline >= 2 ));  // we assume we at least have an allocation for 2 sentinels

  /* Empty the object and move the D+1 sentinel back */
  anch->a[1].i0 = anch->a[anch->D+1].i0;
  anch->a[1].k0 = anch->a[anch->D+1].k0;
  anch->D       = 0;

  /* Downsize the allocation if it's over redline */
  if (anch->nalloc > anch->nredline) {
    ESL_REALLOC(anch->a, sizeof(H4_ANCHOR) * anch->nredline);
    anch->nalloc = anch->nredline;
  }
  return eslOK;

 ERROR:
  return status;
}




/* Function:  h4_anchorset_Destroy()
 * Synopsis:  Free an H4_ANCHORSET
 * Incept:    SRE, Sat 02 Jan 2021
 */
void
h4_anchorset_Destroy(H4_ANCHORSET *anch)
{
  if (anch) {
    free(anch->a);
    free(anch);
  }
  return;
}


/*****************************************************************
 * 2. Debugging, development tools
 *****************************************************************/


/* Function:  h4_anchorset_Dump()
 * Synopsis:  Dump an H4_ANCHORSET
 * Incept:    SRE, Mon 11 Jan 2021 [H9/133]
 */
int
h4_anchorset_Dump(FILE *fp, const H4_ANCHORSET *anch)
{
  int d;

  fprintf(fp, "d    i0    k0\n");
  fprintf(fp, "---- ----- -----\n");

  for (d = 1; d <= anch->D; d++)
    fprintf(fp, "%4d %5d %5d\n", d, anch->a[d].i0, anch->a[d].k0);

  fprintf(fp, "\n# D        = %d\n", anch->D);
  fprintf(fp,   "# nalloc   = %d\n", anch->nalloc);
  fprintf(fp,   "# nredline = %d\n", anch->nredline);
  return eslOK;
}

int
h4_anchorset_DumpOneLine(FILE *ofp, const H4_ANCHORSET *anch)
{
  int d;
  
  fprintf(ofp, "%2d ", anch->D);
  for (d = 1; d <= anch->D; d++)
    fprintf(ofp, "| %4d %4d ", anch->a[d].i0, anch->a[d].k0);
  fprintf(ofp, "\n");
  return eslOK;
}



/* Function:  h4_anchorset_Validate()
 * Synopsis:  Validate an H4_ANCHORSET
 * Incept:    SRE, Tue 12 Jan 2021 [H9/136]
 */
int
h4_anchorset_Validate(const H4_ANCHORSET *anch, char *errmsg)
{
  int M = anch->a[0].k0 - 1;          // the 0 sentinel is (0, M+1);  or it can be 0,0 (not set yet) if D=0, so M=-1 can happen here
  int L = anch->a[anch->D+1].i0 - 1;  // the D+1 sentinel is (L+1,0); or it can be 0,0 (not set yet) if D=0, so L=-1 can happen here
  int D = anch->D;
  int d;

  if (anch->a[0].i0   != 0) ESL_FAIL(eslFAIL, errmsg, "i0(0) sentinel not 0");
  if (anch->a[D+1].k0 != 0) ESL_FAIL(eslFAIL, errmsg, "k0(D+1) sentinel not 0");

  if (D == 0)  // the anchorset is still empty. Now sentinels can be (0,0),(0,0), or they can be (0,M+1),(L+1,0).
    {
      if (L < -1 || L == 0) ESL_FAIL(eslFAIL, errmsg, "sentinels improperly set in empty anchorset"); 
      if (M < -1 || M == 0) ESL_FAIL(eslFAIL, errmsg, "sentinels improperly set in empty anchorset"); 
    }
  else
    {
      if (anch->D < 1)     ESL_FAIL(eslFAIL, errmsg, "number of domains D must be >= 1");
      if (L < 1 || M < 1)  ESL_FAIL(eslFAIL, errmsg, "sentinels improperly set in anchorset");
      for (d = 1; d <= anch->D; d++)
        {
          if (anch->a[d].i0 < 1 || anch->a[d].i0 > L) ESL_FAIL(eslFAIL, errmsg, "bad i0. anchors are i0,k0 for i0=1..L, k0=1..M");
          if (anch->a[d].k0 < 1 || anch->a[d].k0 > M) ESL_FAIL(eslFAIL, errmsg, "bad k0. anchors are i0,k0 for i0=1..L, k0=1..M");
          if (anch->a[d].i0 <= anch->a[d-1].i0)       ESL_FAIL(eslFAIL, errmsg, "anchors need to be ordered in i; i0(d) > i0(d-1)");
        }
    }
  if (anch->nalloc   < anch->D+2) ESL_FAIL(eslFAIL, errmsg, "anchorset allocation is too small (D+2 sentinels)");
  if (anch->nredline <= 2)        ESL_FAIL(eslFAIL, errmsg, "anchorset redline needs to be at least 2, for sentinels");
  return eslOK;
}

/* Function:  h4_anchorset_Compare()
 * Synopsis:  Compare two H4_ANCHORSET's
 * Incept:    SRE, Fri 19 Feb 2021
 */
int
h4_anchorset_Compare(const H4_ANCHORSET *anch1, const H4_ANCHORSET *anch2)
{
  char msg[] = "H4_ANCHORSET comparison failed";
  int D = anch1->D;
  int d;

  if (anch1->D != anch2->D) ESL_FAIL(eslFAIL, NULL, msg);
  for (d = 0; d <= D+1; d++)  // inclusive of sentinels
    {
      if (anch1->a[d].i0 != anch2->a[d].i0) ESL_FAIL(eslFAIL, NULL, msg);
      if (anch1->a[d].k0 != anch2->a[d].k0) ESL_FAIL(eslFAIL, NULL, msg);
    }
  return eslOK;
}



/* Function:  h4_anchorset_Sample()
 * Synopsis:  Sample a randomized anchor set, for testing.
 *
 * Purpose:   Randomly generate an anchor set for a profile of 
 *            length <M> compared to a sequence of length <L>, 
 *            with a random number of up to <maxD> anchors.
 *            Caller provides an allocated but perhaps empty
 *            <anch>, which will be reallocated and reinitialized
 *            here as needed.
 */
int
h4_anchorset_Sample(ESL_RANDOMNESS *rng, int L, int M, int maxD, H4_ANCHORSET *anch)
{
  int      D   = 1 + esl_rnd_Roll(rng, maxD);
  int32_t *tmp = NULL;
  int      i,d,r;
  int      status;

  if ((status = h4_anchorset_GrowFor(anch, D)) != eslOK) goto ERROR;

  /* A reservoir sampling algorithm samples <D> i0 anchors w/o replacement; then sorts them */
  ESL_ALLOC(tmp, sizeof(int32_t) * D);
  for (i = 0; i < L; i++)
    {
      if (i < D) tmp[i] = i+1;   
      else {
	r = esl_rnd_Roll(rng, L);
	if (r < D) tmp[r] = i+1;
      }
    }
  esl_vec_ISortIncreasing(tmp, D);
  
  // then we sample random k0 anchors 1..M
  for (d = 1; d <= D; d++) {
    anch->a[d].i0 = tmp[d-1];                   // the <D> i0's are sorted
    anch->a[d].k0 = 1 + esl_rnd_Roll(rng, M);   // k0's are independent, uniform on 1..M
  } 

  anch->D = D;
  h4_anchorset_SetSentinels(anch, L, M);

  free(tmp);
  return eslOK;

 ERROR:
  free(tmp);
  return status;
}




/* Function:  h4_anchorset_SampleFromPath()
 * Synopsis:  Make a random anchorset that includes a given path.
 * Incept:    SRE, Wed 20 Jan 2021
 *
 * Purpose:   Using <rng>, sample a random anchorset that includes path
 *            <pi>.  For each domain in <pi>, a random match state is
 *            sampled uniformly as the anchor for that domain.
 *
 *            You provide a fresh <anch> structure to use, which
 *            you've initialized to <D=0>, with <L> and <M>, so the
 *            anchorset's sentinels are already set. For example,
 *            <anch = h4_anchorset_Create(0, L, M)>.
 *
 *            If one or more domains in <pi> have no match states,
 *            no anchorset that includes this path is possible.
 *            Return <eslEINVAL> in this case. 
 *
 * Returns:   <eslOK> on success; the new anchorset is in <anch>.
 *            
 *            <eslEINVAL> if one or more domains in <pi> have no match
 *            states; in this case <anch> is reset (with _Reuse())
 *            before returning.
 */
int
h4_anchorset_SampleFromPath(ESL_RANDOMNESS *rng, const H4_PATH *pi, H4_ANCHORSET *anch)
{
  int           i,k;
  int           i0,k0;
  int           nm;
  int           z,r;

  i = 1;
  for (z = 0; z < pi->Z; z++)
    {
      if      (pi->st[z] == h4P_N) i += pi->rle[z]-1;
      else if (pi->st[z] == h4P_G) { nm = 0; k = 1;          }
      else if (pi->st[z] == h4P_L) { nm = 0; k = pi->rle[z]; }
      else if (h4_path_IsM(pi->st[z]))
        {
          for (r = 0; r < pi->rle[z]; r++)
            {
              nm++;
              if (esl_rnd_Roll(rng, nm) == 0) { i0 = i; k0 = k; } // this is a mini reservoir sampling algorithm, w/ reservoir of 1
              i++; k++;
            }
        }
      else if (h4_path_IsI(pi->st[z])) { i += pi->rle[z]; }
      else if (h4_path_IsD(pi->st[z])) { k += pi->rle[z]; }
      else if (pi->st[z] == h4P_J || pi->st[z] == h4P_C)
        {
          if (nm == 0) { h4_anchorset_Reuse(anch); return eslEINVAL; }
          h4_anchorset_Add(anch, i0, k0);
          i += pi->rle[z]-1;
        }
    }
  return eslOK;
}
