/* P7_ANCHORS maintains a resizeable array of i0,k0 anchor points that
 * define domain locations.
 * 
 * Contents:
 *   1. P7_ANCHOR array, either naked or in P7_ANCHORS
 *   2. P7_ANCHORS object
 *   3. Debugging and development tools
 *   4. Unit tests
 *   5. Test driver
 *   6. Copyright and license information
 */
#include "p7_config.h"

#include <stdlib.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "base/p7_anchors.h"

/*****************************************************************
 * 1. P7_ANCHOR arrays, either naked or as part of P7_ANCHORS
 *****************************************************************/

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
 *            resize it as needed, either to a known fixed size <D>
 *            using <p7_anchors_Reinit()>, or by incrementally
 *            expanding the allocation using <p7_anchors_Grow()>.
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
  int32_t     default_nalloc   = 8;     // i.e. for up to D=6
  int32_t     default_nredline = 64;    // i.e. for up to D=62
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

/* Function:  p7_anchors_Reinit()
 * Synopsis:  Reinitialize (and reallocate), for a known number of anchors.
 *
 * Purpose:   Reinitialize <anch> to hold up to <D> anchors. Reallocate
 *            if needed. 
 *
 *            State of the sentinels is undefined. Just like after <_Create()>,
 *            caller needs to call <p7_anchor_SetSentinels(anch->a, D, L, M)> when 
 *            it knows <D>, <L>, and <M>.
 *            
 *            This only reinitializes an object big enough to contain
 *            <D> anchors. It does not set <anch->D>, because the
 *            object does not contain valid anchor data yet. Rather,
 *            as in a <Create()> call, <anch->D> is initialized to
 *            zero. Caller will set D, probably when it sets the
 *            sentinels.
 *
 * Args:      anch  : <P7_ANCHORS> object
 *            D     : number of anchors to prepare for
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_anchors_Reinit(P7_ANCHORS *anch, int D)
{
  int status;

  if (D+2 > anch->nalloc) { // grow?
    ESL_REALLOC(anch->a, sizeof(P7_ANCHOR) * (D+2));  
    anch->nalloc = D+2;
  } 
  else if (anch->nalloc > anch->nredline && D+2 <= anch->nredline) {  // shrink?
    ESL_REALLOC(anch->a, sizeof(P7_ANCHOR) * anch->nredline);
    anch->nalloc = anch->nredline;
  }

  anch->D       = 0;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_anchors_Grow()
 * Synopsis:  Increase allocation for anchors, if needed.
 *
 * Purpose:   Check if there's enough space in <anch> to hold
 *            one new anchor. If not, increase the allocation
 *            in <anch> by doubling it.
 *
 * Args:      anch  : the anchors object
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_anchors_Grow(P7_ANCHORS *anch)
{
  int status;

  if (anch->D+2 < anch->nalloc) return eslOK;                  // +2 because of the sentinels

  ESL_REALLOC(anch->a, sizeof(P7_ANCHOR) * anch->nalloc * 2);  // reallocation may exceed redline - until we Reuse() or Reinit().
  anch->nalloc = anch->nalloc * 2;                            
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

  if (src->D+2 > dst->nalloc &&
      (status = p7_anchors_Reinit(dst, src->D)) != eslOK) goto ERROR;
  
  for (d = 0; d <= src->D+1; d++)       // inclusive of sentinels
    {
      dst->a[d].i0 = src->a[d].i0;
      dst->a[d].k0 = src->a[d].k0;
    }
  dst->D    = src->D;
  return eslOK;

 ERROR:
  return status;
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
p7_anchors_Dump(FILE *fp, P7_ANCHORS *anch)
{
  int d;

  fprintf(fp, "# anchors: for D=%d\n", anch->D);
  for (d = 0; d <= anch->D+1; d++)
    fprintf(fp, "%-4d %6d %6d\n", d, anch->a[d].i0, anch->a[d].k0);
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

int
p7_anchors_Validate(P7_ANCHORS *anch, char *errbuf)
{
  int D = anch->D;
  int L = anch->a[D+1].i0 - 1;
  int M = anch->a[0].k0 - 1;
  int d;

  for (d = 0; d <= D; d++)
    if (! (anch->a[d].i0 < anch->a[d+1].i0)) ESL_FAIL(eslFAIL, errbuf, "i0 anchors not sorted");

  for (d = 1; d <= D; d++) {
    if (! (anch->a[d].i0 >= 1 && anch->a[d].i0 <= L)) ESL_FAIL(eslFAIL, errbuf, "i0 %d not in range 1..L", d);
    if (! (anch->a[d].k0 >= 1 && anch->a[d].k0 <= M)) ESL_FAIL(eslFAIL, errbuf, "k0 %d not in range 1..M", d);
  }
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

  if ((status = p7_anchors_Reinit(anch, D)) != eslOK) goto ERROR;

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
p7_anchors_SampleFromTrace(ESL_RANDOMNESS *rng, const P7_TRACE *tr, P7_ANCHORS *anch)
{
  int D = tr->ndom;
  int d,z,w;
  int nM;
  int status;

  if ((status = p7_anchors_Reinit(anch, D)) != eslOK) goto ERROR;
  
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
      if ( p7_anchors_Validate(anch, errmsg)        != eslOK) esl_fatal("%s:\n  %s", msg, errmsg);
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
#include "p7_config.h"

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


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

