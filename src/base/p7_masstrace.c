/* The "mass trace" algorithm determines the bounds of the probability
 * envelope around each domain. The algorithm needs some auxiliary
 * data storage, which is implemented in the P7_MASSTRACE object.
 * 
 * Contents:
 *    1. P7_MASSTRACE object
 *    2. Debugging tools
 */

#include "p7_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "base/p7_trace.h"
#include "base/p7_masstrace.h"

/*****************************************************************
 * 1. P7_MASSTRACE object
 *****************************************************************/

static P7_MASSTRACE *
masstrace_create_engine(int M_hint, int L_hint, int do_slim)
{
  P7_MASSTRACE *mt = NULL;
  int           status;

  ESL_ALLOC(mt, sizeof(P7_MASSTRACE));
  mt->kmass  = NULL;
  mt->imass  = NULL;
  mt->kalloc = 0;
  mt->ialloc = 0;

  /* i0,k0,L,M are set by someone who fills in imass[0.1..L.L+1],
   * kmass[0.1..M.M+1]. Until then, they are 0; if the object
   * is discarded (by a Reuse()), they are reset to 0.
   */
  mt->i0    = 0;   
  mt->k0    = 0;
  mt->st0   = 0; /* p7T_BOGUS */
  mt->L     = 0;   
  mt->M     = 0;   

  ESL_ALLOC(mt->kmass, sizeof(float) * (M_hint+2));
  mt->kalloc = M_hint+2;

  if (!do_slim) {
    ESL_ALLOC(mt->imass, sizeof(float) * (L_hint+2));
    mt->ialloc = L_hint+2;
  }

  return mt;
  
 ERROR:
  p7_masstrace_Destroy(mt);
  return NULL;
}  

P7_MASSTRACE *
p7_masstrace_Create(int M_hint, int L_hint)
{ 
  return masstrace_create_engine(M_hint, L_hint, FALSE); 
}


P7_MASSTRACE *
p7_masstrace_CreateSlim(int M_hint, int L_hint)
{
  return masstrace_create_engine(M_hint, L_hint, TRUE);
}

int
p7_masstrace_GrowTo(P7_MASSTRACE *mt, int M, int L)
{
  int status;
  
  if (mt->imass && mt->ialloc < L+2) {
    ESL_REALLOC(mt->imass, sizeof(float) * (L+2));
    mt->ialloc = L+2;
  }
  if (mt->kalloc < M+2) {
    ESL_REALLOC(mt->kmass, sizeof(float) * (M+2));
    mt->kalloc = M+2;
  }
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_masstrace_Zero()
 * Synopsis:  Initialize cumulative endpoint distributions to zeros.
 *
 * Purpose:   Zero the cumulative distributions in <mt>, preparing to
 *            collect masstrace endpoint data for a sequence of length
 *            <L> and a profile of length <M>.
 *
 * Args:      mt - mass trace object to collect endpoint data in
 *            M  - profile length
 *            L  - sequence length
 *
 * Returns:   <eslOK> on success.
 */
int
p7_masstrace_Zero(P7_MASSTRACE *mt, int M, int L)
{
  /* contract checks / argument validation */
  ESL_DASSERT1( (mt->imass == NULL || L+2 <= mt->ialloc ) );
  ESL_DASSERT1( (M+2 <= mt->kalloc) );

  if (mt->imass) esl_vec_FSet(mt->imass, L+2, 0.0f);
  esl_vec_FSet(mt->kmass, M+2, 0.0f);
  mt->L = L;
  mt->M = M;
  return eslOK;
}

int
p7_masstrace_Reuse(P7_MASSTRACE *mt)
{
  mt->i0  = 0;
  mt->k0  = 0;
  mt->st0 = 0;
  mt->L   = 0;
  mt->M   = 0;
  return eslOK;
}

void
p7_masstrace_Destroy(P7_MASSTRACE *mt)
{
  if (mt) {
    if (mt->imass) free(mt->imass);
    if (mt->kmass) free(mt->kmass);
    free(mt);
  }
  return;
}


/*****************************************************************
 * 2. Debugging tools 
 *****************************************************************/

/* Function:  p7_masstrace_CountTrace()
 * Synopsis:  Count domain endpoints into endpoint distributions.
 *
 * Purpose:   Given a traceback <tr>, determine if it contains a domain
 *            specified by the anchor point <i0>, <k0>, <st0>. If it
 *            doesn't, return without doing anything. If it does,
 *            count that domain's start/end points into the <mt>
 *            structure, and update the <*ntr> count by one.
 *            
 *            This function is useful in unit tests that approximate
 *            the mass trace calculation using a large ensemble of
 *            stochastic tracebacks.
 *            
 *            Before the first <_CountTrace()> call on an <mt>, you
 *            call <p7_masstrace_Zero()> on it to initialize it.
 *            
 *            The counts in <mt> are collected as a histogram; after
 *            an entire ensemble has been collected, <mt> needs to be
 *            converted to a cumulative distribution.
 *            
 *            A special case arises exactly at the 'midpoints' in the
 *            <mt> vectors <kmass> and <imass>. <kmass[k]> will be the
 *            start point cumulative distribution P(ka <= k) for
 *            k<=k0, and the end point cumulative distribution P(kb >=
 *            k) for k>=k0. In a cumulative distributoin, it's ok that
 *            kmass[k0] is defined as both the start and end value,
 *            because it's 1.0 in both cases. But in a histogram, we
 *            would have to distinguish whether kmass[k0] has seen a
 *            start ka versus an end kb. Instead of doing something
 *            special to handle this, instead we don't count kmass[k0]
 *            (or imass[i0]) at all; and when we convert to a
 *            cumulative distribution, we'll set these to 1.0.
 *        
 *            Because this is counting in single-precision floating
 *            point arithmetic, it can't accurately count an ensemble
 *            of more than about $10^7$ traces.
 *
 * Args:      tr  - trace structure 
 *            i0  - i sequence position coord of domain anchor
 *            k0  - k model position coord of domain anchor
 *            st0 - a main model {MID}{LG} state type of domain anchor
 *            mt  - mass trace object to count endpoints in
 *            ntr - updated number of traces that contain the anchor
 *
 * Returns:   <eslOK> on success. "Success" includes ignoring a trace
 *            that does not contain the anchor <i0>,<k0>,<st0>. If the
 *            trace does contain the anchor, start/endpoint counts in
 *            <mt> are incremented by one, and <*ntr> is incremented
 *            by one.
 */
int
p7_masstrace_CountTrace(const P7_TRACE *tr, int i0, int k0, int st0, P7_MASSTRACE *mt, int *ntr)
{
  int i,z0,z;
  int ia,ib, ka,kb;
  int foundit = FALSE;

  /* Contract checks on arguments */
  ESL_DASSERT1( ( i0>=1 && i0 <= tr->L) );
  ESL_DASSERT1( ( k0>=1 && k0 <= tr->M) );
  ESL_DASSERT1( ( p7_trace_IsMain(st0)) );
  ESL_DASSERT1( ( mt->i0  == 0 || mt->i0  == i0) );
  ESL_DASSERT1( ( mt->k0  == 0 || mt->k0  == k0) );
  ESL_DASSERT1( ( mt->st0 == 0 || mt->st0 == st0) );
  ESL_DASSERT1( ( mt->L  == tr->L) );
  ESL_DASSERT1( ( mt->M  == tr->M) );

  /* Find the anchor, if it's there. */
  for (i=0, z0 = 0; z0 < tr->N; z0++)
    {
      if (tr->i[z0]) i = tr->i[z0]; /* update i. only emitting states have tr->i[z] set */
      if (i > i0 )  break;	  /* failed to find anchor. */
      if (i == i0 && tr->st[z0] == st0 && tr->k[z0] == k0) { foundit = TRUE; break; }
    }
  if (! foundit) return eslOK;	/* If no anchor: successful return, ignoring this trace. */

  /* Find leftmost bounds of domain */
  for (ia = i0, ka = k0, z = z0; z >= 0 && tr->st[z] != p7T_B; z--) 
    {
      if (tr->i[z]) ia = tr->i[z];
      if (tr->k[z]) ka = tr->k[z];
    }
  ESL_DASSERT1( ( tr->st[z] == p7T_B) );

  /* Find rightmost bounds of domain */
  for (ib = i0, kb = k0, z = z0; z < tr->N && tr->st[z] != p7T_E; z++)
    {
      if (tr->i[z]) ib = tr->i[z];
      if (tr->k[z]) kb = tr->k[z];
    }
  ESL_DASSERT1( ( tr->st[z] == p7T_E) );

  /* Increment counters */
  if (ka < k0)              mt->kmass[ka] += 1.; /* note the guards against incrementing the overlapped start/end at k0,i0 */
  if (kb > k0)              mt->kmass[kb] += 1.;
  if (mt->imass && ia < i0) mt->imass[ia] += 1.; /* also, guard for the optional <imass> data in <mt> */
  if (mt->imass && ib > i0) mt->imass[ib] += 1.;
  *ntr += 1;

  /* Make sure i0,k0,st0 are set. */
  mt->i0  = i0;
  mt->k0  = k0;
  mt->st0 = st0;
  return eslOK;
}

/* Function:  p7_masstrace_FinishCount()
 * Synopsis:  Convert counted histograms to cumulative endpoint prob distributions.
 *
 * Purpose:   We've finished collecting endpoints from traces with
 *            <_CountTrace()> in <mt>, <ntr> of which had the
 *            specified domain anchor; now convert the counts to
 *            <mt>'s cumulative probability distributions.  
 *
 * Args:      mt  - mass trace object we've collected endpoint counts in
 *            ntr - number of traces we counted into <mt> that contained the domain anchor
 *
 * Returns:   <eslOK> on success; <mt> is now a valid <P7_MASSTRACE> object
 *            containing envelope endpoint cumulative probability distributions.
 */
int
p7_masstrace_FinishCount(P7_MASSTRACE *mt, int ntr)
{
  int i,k;

  ESL_DASSERT1( (ntr > 0) );
  ESL_DASSERT1( (mt->i0)  );
  ESL_DASSERT1( (mt->k0)  );
  ESL_DASSERT1( (p7_trace_IsMain(mt->st0)) );

  if (mt->imass)
    {
      for (i = 1;     i < mt->i0; i++) mt->imass[i] += mt->imass[i-1];
      for (i = mt->L; i > mt->i0; i--) mt->imass[i] += mt->imass[i+1];
      esl_vec_FScale(mt->imass+1, mt->L, 1./(float) ntr);
      mt->imass[mt->i0] = 1.;
    }
  for (k = 1;     k < mt->k0; k++) mt->kmass[k] += mt->kmass[k-1];
  for (k = mt->M; k > mt->k0; k--) mt->kmass[k] += mt->kmass[k+1];
  esl_vec_FScale(mt->kmass+1, mt->M, 1./(float) ntr);
  mt->kmass[mt->k0] = 1.;
  return eslOK;
}


int
p7_masstrace_Dump(FILE *ofp, P7_MASSTRACE *mt)
{
  int i,k;

  if (mt->imass) {
    fprintf(ofp, "# i (sequence) cumulative endpoint distributions\n");
    for (i = 0; i <= mt->L+1; i++)
	fprintf(ofp, "%-6d %.5f%c\n", i, mt->imass[i], (i==0 || i==mt->L+1 || i==mt->i0) ? '*' : ' ');
  }

  fprintf(ofp, "# k (model) cumulative endpoint distributions\n");
  for (k = 0; k <= mt->M+1; k++)
    fprintf(ofp, "%-6d %.5f%c\n", k, mt->kmass[k], (k==0 || k==mt->M+1 || k==mt->k0) ? '*' : ' ');

  fprintf(ofp, "anchor = %s\n", p7_trace_DecodeStatetype(mt->st0));
  return eslOK;
}


int
p7_masstrace_PlotImass(FILE *ofp, P7_MASSTRACE *mt)
{
  int i;
  for (i = 1; i <= mt->L; i++)
    fprintf(ofp, "%f\n", mt->imass[i]);
  fprintf(ofp, "&\n");
  return eslOK;
}

int
p7_masstrace_PlotKmass(FILE *ofp, P7_MASSTRACE *mt)
{
  int k;
  for (k = 1; k <= mt->M; k++)
    fprintf(ofp, "%f\n", mt->kmass[k]);
  fprintf(ofp, "&\n");
  return eslOK;
}

int
p7_masstrace_Compare(const P7_MASSTRACE *mte, const P7_MASSTRACE *mta, float tol)
{
  char msg[] = "masstrace object comparison failed";
  int i,k;

  if (mte->L   != mta->L)   ESL_FAIL(eslFAIL, NULL, msg);
  if (mte->M   != mta->M)   ESL_FAIL(eslFAIL, NULL, msg);
  if (mte->i0  != mta->i0)  ESL_FAIL(eslFAIL, NULL, msg);
  if (mte->k0  != mta->k0)  ESL_FAIL(eslFAIL, NULL, msg);
  if (mte->st0 != mta->st0) ESL_FAIL(eslFAIL, NULL, msg);

  if (mte->imass && mta->imass)
    {
      for (i = 1; i <= mte->L; i++)
	{
	  if (mte->imass[i] == 0.0 && mta->imass[i] > 0.0)                  ESL_FAIL(eslFAIL, NULL, msg);
	  if (esl_FCompareAbs(mte->imass[i], mta->imass[i], tol) != eslOK)  ESL_FAIL(eslFAIL, NULL, msg);
	}
    }
  for (k = 1; k <= mte->M; k++)
    {
      if (mte->kmass[k] == 0.0 && mta->kmass[k] > 0.0)                  ESL_FAIL(eslFAIL, NULL, msg);
      if (esl_FCompareAbs(mte->kmass[k], mta->kmass[k], tol) != eslOK)  ESL_FAIL(eslFAIL, NULL, msg);
    }
  return eslOK;
}
  
float
p7_masstrace_GetMaxAbsDiff(const P7_MASSTRACE *mte, const P7_MASSTRACE *mta)
{
  int i,k;
  float diff;
  float max = 0.; 

  if (mte->imass && mta->imass)
    {
      for (i = 1; i <= mte->L; i++)
	{
	  diff = fabs(mte->imass[i] - mta->imass[i]); 
	  max  = ESL_MAX(diff, max);
	}
    }
  for (k = 1; k <= mte->M; k++)
    {
      diff = fabs(mte->kmass[k] - mta->kmass[k]); 
      max  = ESL_MAX(diff, max);
    }
  return max;
}


int
p7_masstrace_Validate(const P7_MASSTRACE *mt, char *errbuf)
{
  float tol = 1e-3;
  int i,k;

  if (mt->L  <= 0)                           ESL_FAIL(eslFAIL, errbuf, "L=0");
  if (mt->M  <= 0)                           ESL_FAIL(eslFAIL, errbuf, "L=0");
  if (mt->i0 < 1 || mt->i0 > mt->L)          ESL_FAIL(eslFAIL, errbuf, "i0 range");
  if (mt->k0 < 1 || mt->k0 > mt->M)          ESL_FAIL(eslFAIL, errbuf, "k0 range");
  if ( ! p7_trace_IsMain(mt->st0))           ESL_FAIL(eslFAIL, errbuf, "st0 not {MID}{LG}");
  if (mt->imass && mt->imass[0]       != 0.) ESL_FAIL(eslFAIL, errbuf, "imass[0] not 0");
  if (mt->imass && mt->imass[mt->i0]  != 1.) ESL_FAIL(eslFAIL, errbuf, "imass[i0] not 1");
  if (mt->imass && mt->imass[mt->L+1] != 0.) ESL_FAIL(eslFAIL, errbuf, "imass[L+1] not 0");
  if (mt->kmass[0]       != 0.)              ESL_FAIL(eslFAIL, errbuf, "kmass[0] not 0");
  if (mt->kmass[mt->k0]  != 1.)              ESL_FAIL(eslFAIL, errbuf, "kmass[k0] not 1");
  if (mt->kmass[mt->M+1] != 0.)              ESL_FAIL(eslFAIL, errbuf, "kmass[M+1] not 0");
  if (mt->imass) {
    for (i = 0; i <= mt->L+1; i++) 
      if (!isfinite(mt->imass[i]) || mt->imass[i] < 0.0 || mt->imass[i] > 1+tol) 
	ESL_FAIL(eslFAIL, errbuf, "imass[%d] isn't a probability: %f\n", i, mt->imass[i]);
  }
  for (k = 0; k <= mt->M+1; k++) 
    if (!isfinite(mt->kmass[k]) || mt->kmass[k] < 0.0 || mt->kmass[k] > 1+tol) 
      ESL_FAIL(eslFAIL, errbuf, "kmass[%d] isn't a probability: %f\n", k, mt->kmass[k]);
  return eslOK;
}
