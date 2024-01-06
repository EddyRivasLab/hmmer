/* P7_ANCHORS maintains a resizeable array of i0,k0 anchor points that
 * define domain locations.
 * 
 * Contents:
 *   1. P7_ANCHOR array, either naked or in P7_ANCHORS
 *   2. P7_ANCHORS object
 *   3. Debugging and development tools
 *   4. Unit tests
 *   5. Test driver
 */
#include <p7_config.h>

#include <stdlib.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "base/p7_anchors.h"

/*****************************************************************
 * 1. P7_ANCHOR arrays, either naked or as part of P7_ANCHORS
 *****************************************************************/

/* Function:  p7_anchor_GetSentinels()
 * Synopsis:  Get L,M from existing sentinels
 *
 * Purpose:   Use the existing sentinels at <0,D+1> to retrieve
 *            sequence length <L> and profile length <M>. This gets
 *            used when we need to change the anchor set length
 *            and set new sentinels.
 */
int
p7_anchor_GetSentinels(P7_ANCHOR *arr, int D, int *ret_L, int *ret_M)
{
  *ret_L = arr[D+1].i0 - 1;
  *ret_M = arr[0].k0   - 1;
  return eslOK;
}


/* Function:  p7_anchor_SetSentinels()
 * Synopsis:  Initialize sentinel values for a P7_ANCHOR array.
 *
 * Purpose:   Initialize sentinels <0> and <D+1> in a <P7_ANCHOR>
 *            array <anch> defining <1..D> domains, for a (sub)sequence
 *            of length <L> and a profile of length <M>.
 *            
 *            <anch[0]> is set to <(i0,k0) = (0, M+1)>.
 *            <anch[D+1] is set to <(i0,k0) = (L+1, 0)>.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_anchor_SetSentinels(P7_ANCHOR *anch, int D, int L, int M)
{
  anch[0].i0   = 0;      // UP(d) sector runs i0(d-1)..i0(d)-1: so for UP(1), i0(0) must evaluate to exactly 0
  anch[D+1].i0 = L+1;    // DOWN(d) sector runs i0(d)..i0(d+1)-1; so for DOWN(D), i0(D+1) must evaluate to exactly L+1

  anch[0].k0   = M+1;    // DOWN(d) sector runs k0(d)..M. For access patterns that do DOWN(d-1),UP(d) on row i, must have k0(0) > M so DOWN(0) isn't evaluated;    M+1 suffices
  anch[D+1].k0 = 0;      // UP(d) secor runs 1..k0(d)-1.  For access patterns that do DOWN(d-1),UP(d) on row i, must have k0(D+1) <= 1 so UP(D+1) isn't evaluated; 0 suffices
  return eslOK;
}


/*****************************************************************
 * 2. P7_ANCHORS object.
 *****************************************************************/

/* Function:  p7_anchors_Create()
 * Synopsis:  Create a new <P7_ANCHORS>
 *
 * Purpose:   Returns a pointer to a newly created <P7_ANCHORS> 
 *            object. 
 *            
 *            The object is allocated with a default initial size,
 *            which happens to be <D=8> right now. The caller will
 *            resize it as needed using <p7_anchors_Resize()>.
 *            
 *            Sentinels are not set. Caller needs to call 
 *            <p7_anchor_SetSentinels(anch->a, D, L, M)> when 
 *            it knows <D>, <L>, and <M>.
 *
 * Args:      (none)
 *
 * Returns:   pointer to the new <P7_ANCHORS>.
 *
 * Throws:    <NULL> on allocation failure.
 *
 * Notes:     If needed, we can add a <p7_anchors_CreateCustom()> call later,
 *            and allow caller to pass D and/or D_redline.
 */
P7_ANCHORS *
p7_anchors_Create(void)
{
  P7_ANCHORS *anch             = NULL;
  int32_t     default_nalloc   = 8;     // i.e. for up to D=6;  +2 sentinels
  int32_t     default_nredline = 64;    // i.e. for up to D=62; +2 sentinels
  int         status;

  ESL_ALLOC(anch, sizeof(P7_ANCHORS));
  anch->a        = NULL;
  anch->D        = 0;
  anch->nalloc   = default_nalloc;
  anch->nredline = default_nredline;

  ESL_ALLOC(anch->a, sizeof(P7_ANCHOR) * (anch->nalloc));
  return anch;

 ERROR:
  p7_anchors_Destroy(anch);
  return NULL;
}


/* Function:  p7_anchors_Resize()
 * Synopsis:  Reallocate a P7_ANCHORS object, if necessary
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
 * Xref:      First example of a new pattern for how we 
 *            can handle reallocation/reuse strategy, 
 *            replacing _Reinit() and _Grow() interfaces.
 *            [SRE:J14/1]
 */
int
p7_anchors_Resize(P7_ANCHORS *anch, int D)
{
  int nalloc;
  int status;
  
  /* Contract checks, argument validation */
  ESL_DASSERT1(( anch->nalloc > 0 ));               

  if      (D+2 <= anch->nalloc) return eslOK;       // If we're big enough already, do nothing;
  else if (D+2 <  anch->nredline || anch->D > 0)    // If we're under the redline max, or if it looks like
    {                                               //   we're building the anchor array incrementally,
      nalloc = anch->nalloc;                        //   we reallocate by doubling, trying to minimize 
      while (nalloc < D+2) nalloc *= 2;             //   the need for more reallocations soon.
    }                                               // If we're over redline AND it looks like we're
  else nalloc = D+2;                                //   starting an empty object, allocate exactly.
                                                    //   Now nalloc will probably not be a multiple of two -- 
                                                    //   but the next _Reuse() call will pull it back
                                                    //   to the redline, which is.
  ESL_REALLOC(anch->a, sizeof(P7_ANCHOR) * nalloc);
  anch->nalloc = nalloc;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_anchors_Copy()
 * Synopsis:  Copy one P7_ANCHORS object to another.
 *
 * Purpose:   Copy <src> to <dst>. <dst> is reallocated
 *            if necessary. 
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_anchors_Copy(const P7_ANCHORS *src, P7_ANCHORS *dst)
{
  int32_t d;
  int     status;

  if (( status = p7_anchors_Resize(dst, src->D) ) != eslOK) goto ERROR;
  
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


/* Function:  p7_anchors_Catenate()
 * Synopsis:  Append anchors to an anchor set, possibly with some overwriting.
 *
 * Purpose:   Append the <Dg> anchors in <arr[1..Dg]> to the anchor set
 *            in <anch->a[1..D0]> to form a new anchor set in <anch> of
 *            <D=D0+Dg> anchors; reset the sentinels at <0,D0+Dg+1>. 
 *
 *            This is used in segmental divide and conquer MPAS
 *            algorithm, where the optimal anchor set is determined
 *            segment by segment.  We have a prefix of <D0> anchors
 *            that are already optimal, and we are testing different
 *            suffixes, overwriting suffixes as we try them.
 *            
 *            <anch> is reallocated, if necessary, to hold the new
 *            anchor set.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on reallocation failure.           
 */
int
p7_anchors_Catenate(P7_ANCHORS *anch, int D0, P7_ANCHOR *arr, int Dg)
{
  int L;
  int M;
  int d;
  int status;

  if (( status = p7_anchors_Resize(anch, D0+Dg)) != eslOK) return status;

  p7_anchor_GetSentinels(anch->a, anch->D, &L, &M);
  for (d = 1; d <= Dg; d++)
    {
      anch->a[D0+d].i0 = arr[d].i0;
      anch->a[D0+d].k0 = arr[d].k0;
    }
  anch->D = D0 + Dg;
  p7_anchor_SetSentinels(anch->a, anch->D, L, M);
  return eslOK;
}




int
p7_anchors_Reuse(P7_ANCHORS *anch)
{
  int status;

  if (anch->nalloc > anch->nredline) {
    ESL_REALLOC(anch->a, sizeof(P7_ANCHOR) * anch->nredline);
    anch->nalloc = anch->nredline;
  }
  anch->D = 0;
  return eslOK;

 ERROR:
  return status;
}

void
p7_anchors_Destroy(P7_ANCHORS *anch)
{
  if (anch) {
    if (anch->a) free(anch->a);
    free(anch);
  }
  return;
}


/*****************************************************************
 * 3. Debugging and development tools
 *****************************************************************/

int
p7_anchors_Dump(FILE *fp, const P7_ANCHORS *anch)
{
  int d;

  fprintf(fp, "# anchors: for D=%d\n", anch->D);
  for (d = 0; d <= anch->D+1; d++)
    fprintf(fp, "%-4d %6d %6d\n", d, anch->a[d].i0, anch->a[d].k0);
  return eslOK;
}

int
p7_anchors_DumpOneLine(FILE *ofp, const P7_ANCHORS *anch)
{
  int d;
  
  fprintf(ofp, "%2d ", anch->D);
  for (d = 1; d <= anch->D; d++)
    fprintf(ofp, "%4d %4d ", anch->a[d].i0, anch->a[d].k0);
  fprintf(ofp, "\n");
  return eslOK;
}


int
p7_anchors_Compare(P7_ANCHORS *anch1, P7_ANCHORS *anch2)
{
  int D = anch1->D;
  int d;

  if (anch1->D != anch2->D) return eslFAIL;
  for (d = 0; d <= D+1; d++)  // inclusive of sentinels
    {
      if (anch1->a[d].i0 != anch2->a[d].i0) return eslFAIL;
      if (anch1->a[d].k0 != anch2->a[d].k0) return eslFAIL;
    }
  return eslOK;
}

/* Function:  p7_anchors_Validate()
 * Synopsis:  Validate an anchor set object. 
 *
 * Purpose:   Validates an anchor set object.
 * 
 *            If <M>,<L> dimensions are provided, then the sentinels
 *            at <0> and <D+1> are validated too. If <M> or <L> are
 *            unknown they can be passed as 0, and the sentinels in
 *            <anch> will be used to determine them -- which of course
 *            depends on the sentinels being valid, so is less strong.
 *
 * Args:      anch   - anchors to validate
 *            L      - sequence length if known; else 0
 *            M      - profile length if known; else 0
 *            errbuf - optional error message, allocated for eslERRBUFSIZE; or NULL
 *
 * Returns:   <eslOK> on success.
 *            <eslFAIL> on failure, and if <errbuf> was provided, it contains
 *            an informative error message.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_anchors_Validate(P7_ANCHORS *anch, int L, int M, char *errbuf)
{
  int D = anch->D;
  int d;

  /* If M or L aren't provided, set them from the sentinels */
  if (!L) L = anch->a[D+1].i0 - 1;
  if (!M) M = anch->a[0].k0 - 1;

  for (d = 0; d <= D; d++)
    if (! (anch->a[d].i0 < anch->a[d+1].i0)) ESL_FAIL(eslFAIL, errbuf, "i0 anchors not sorted");

  for (d = 1; d <= D; d++) {
    if (! (anch->a[d].i0 >= 1 && anch->a[d].i0 <= L)) ESL_FAIL(eslFAIL, errbuf, "i0 %d not in range 1..L", d);
    if (! (anch->a[d].k0 >= 1 && anch->a[d].k0 <= M)) ESL_FAIL(eslFAIL, errbuf, "k0 %d not in range 1..M", d);
  }

  if (anch->a[0].i0   != 0   || anch->a[0].k0   != M+1) ESL_FAIL(eslFAIL, errbuf, "sentinel 0 invalid");
  if (anch->a[D+1].i0 != L+1 || anch->a[D+1].k0 != 0)   ESL_FAIL(eslFAIL, errbuf, "sentinel D+1 invalid");

  return eslOK;
}




/* Function:  p7_anchors_Sample()
 * Synopsis:  Sample a randomized anchor set, for testing.
 *
 * Purpose:   Randomly generate an anchor set for a profile of 
 *            length <M> compared to a sequence of length <L>, 
 *            with a random number of up to <maxD> anchors.
 */
int
p7_anchors_Sample(ESL_RANDOMNESS *rng, int L, int M, int maxD, P7_ANCHORS *anch)
{
  int      D   = 1 + esl_rnd_Roll(rng, maxD);
  int32_t *tmp = NULL;
  int      i,d,r;
  int      status;

  if ((status = p7_anchors_Resize(anch, D)) != eslOK) goto ERROR;

  /* A reservoir sort like algorithm samples a combination of <D> i0 anchors, w/o replacement */
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
  
  for (d = 1; d <= D; d++) {
    anch->a[d].i0 = tmp[d-1];                   // the <D> i0's are sorted
    anch->a[d].k0 = 1 + esl_rnd_Roll(rng, M);   // k0's are independent, uniform on 1..M
  } 

  p7_anchor_SetSentinels(anch->a, D, L, M);
  anch->D = D;

  free(tmp);
  return eslOK;

 ERROR:
  if (tmp) free(tmp);
  return status;
}


/* Function:  p7_anchors_SampleFromTrace()
 * Synopsis:  Make a reasonable anchor set from a trace.
 *
 * Purpose:   Make a reasonable anchor set from trace <tr>, by
 *            randomly sampling a match state in each domain.
 *            Return the anchor set in <anch>, which will be
 *            reallocated if needed.
 *            
 *            <tr> must be indexed by the caller with <p7_trace_Index()>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_anchors_SampleFromTrace(P7_ANCHORS *anch, ESL_RANDOMNESS *rng, const P7_TRACE *tr)
{
  int D = tr->ndom;
  int d,z,w;
  int nM;
  int status;

  if ((status = p7_anchors_Resize(anch, D)) != eslOK) goto ERROR;
  
  for (d = 1; d <= D; d++)
    {
      for (nM = 0, z = tr->tfrom[d-1]; z <= tr->tto[d-1]; z++)   // P7_TRACE numbers domains 0..D-1, off by one from P7_ANCHORS
	if (p7_trace_IsM(tr->st[z])) nM++;
      ESL_DASSERT1(( nM ));

      w = 1+esl_rnd_Roll(rng, nM);               // w = 1..nM : choice of which M state to make the anchor
      
      for ( z = tr->tfrom[d-1]; w; z++)          // when w reaches 0, tr->st[z] is the M state we want to make the anchor, and we break out; there's a final z++, so the state we want ends up being z-1
	if (p7_trace_IsM(tr->st[z])) w--;   
      ESL_DASSERT1(( p7_trace_IsM(tr->st[z-1]) )); // since the logic above is overly elegant... better doublecheck.
      
      anch->a[d].i0 = tr->i[z-1];
      anch->a[d].k0 = tr->k[z-1];
    }
  
  p7_anchor_SetSentinels(anch->a, D, tr->L, tr->M);
  anch->D = D;
  return eslOK;

 ERROR:
  return status;
}


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7ANCHORS_TESTDRIVE

#include "easel.h"
#include "esl_random.h"

#include "hmmer.h"

static void
utest_sampling(ESL_RANDOMNESS *rng)
{
  char        msg[]    = "p7_anchors.c sampling unit test failed";
  P7_ANCHORS *anch     = p7_anchors_Create();
  P7_ANCHORS *anch2    = p7_anchors_Create();
  int         nsamples = 100;
  int         L        = 400;
  int         M        = 400;
  int         maxD     = 100;
  int         s;
  char        errmsg[eslERRBUFSIZE];

  for (s = 0; s < nsamples; s++)
    {
      if ( p7_anchors_Sample(rng, L, M, maxD, anch) != eslOK) esl_fatal(msg);
      if ( p7_anchors_Copy(anch, anch2)             != eslOK) esl_fatal(msg);
      if ( p7_anchors_Compare(anch, anch2)          != eslOK) esl_fatal(msg);
      if ( p7_anchors_Validate(anch, L, M, errmsg)  != eslOK) esl_fatal("%s:\n  %s", msg, errmsg);
      if ( p7_anchors_Reuse(anch)                   != eslOK) esl_fatal(msg);
    }
  p7_anchors_Destroy(anch);
  p7_anchors_Destroy(anch2);
}

#endif /*p7ANCHORS_TESTDRIVE*/




/*****************************************************************
 * 5. Test driver
 *****************************************************************/

#ifdef p7ANCHORS_TESTDRIVE
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for p7_anchors.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_sampling(rng);

  fprintf(stderr, "#  status = ok\n");

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  exit(0); /* success */
}

#endif /*p7ANCHORS_TESTDRIVE*/
/*-------------------- end of test driver ---------------------*/


