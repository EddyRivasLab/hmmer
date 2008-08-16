/* P7_TRACE, the traceback structure.
 *
 * Contents:
 *   1. The P7_TRACE structure and API
 *   2. Copyright and license information.
 * 
 * Stylistic note: paths are indexed by z.
 * 
 * SRE, Tue Jan  2 2007 [Casa de Gatos] 
 * SVN $Id$
 */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "hmmer.h"

/* Function:  p7_trace_Create()
 * Synopsis:  Allocates a (growable, reusable) traceback.
 *
 * Purpose:   Allocates a traceback. 
 *  
 *            Tracebacks are growable. A reasonable initial internal
 *            allocation is made here, and routines that generate
 *            tracebacks will dynamically grow the trace as needed.
 *            
 *            Tracebacks are reusable. Usually a routine only
 *            allocates one, and reuses its memory over and over as
 *            new target sequences are aligned.
 *
 * Returns:   a pointer to the new <P7_TRACE> structure on success.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_TRACE *
p7_trace_Create(void)
{
  P7_TRACE *tr                = NULL;
  int       initial_nalloc    = 256;
  int       initial_ndomalloc = 16;
  int       status;

  ESL_ALLOC(tr, sizeof(P7_TRACE));
  tr->st = NULL;
  tr->k  = NULL;
  tr->i  = NULL;
  tr->tfrom   = tr->tto   = NULL;
  tr->sqfrom  = tr->sqto  = NULL;
  tr->hmmfrom = tr->hmmto = NULL;

  /* The trace data itself */
  ESL_ALLOC(tr->st, sizeof(char) * initial_nalloc);
  ESL_ALLOC(tr->k,  sizeof(int)  * initial_nalloc);
  ESL_ALLOC(tr->i,  sizeof(int)  * initial_nalloc);
  tr->N      = 0;
  tr->nalloc = initial_nalloc;

  /* The trace's index: table of domain start/stop coords */
  ESL_ALLOC(tr->tfrom,   sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->tto,     sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->sqfrom,  sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->sqto,    sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->hmmfrom, sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->hmmto,   sizeof(int) * initial_ndomalloc);
  tr->ndom      = 0;
  tr->ndomalloc = initial_ndomalloc;
  return tr;

 ERROR:
  if (tr != NULL) p7_trace_Destroy(tr);
  return NULL;
}


/* Function:  p7_trace_Reuse()
 * Synopsis:  Prepare a trace for reuse.
 * Incept:    SRE, Tue Jan  9 13:02:34 2007 [Janelia]
 *
 * Purpose:   Reinitializes an existing trace object, reusing its
 *            memory.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      STL11/124
 */
int
p7_trace_Reuse(P7_TRACE *tr)
{
  tr->N    = 0;
  tr->ndom = 0;
  return eslOK;
}

/* Function:  p7_trace_Grow()
 * Synopsis:  Grow the allocation for trace data.
 *
 * Purpose:   If <tr> can't fit another state, double its allocation for
 *            traceback data.
 *            
 *            This doesn't reallocate the domain index; see
 *            <p7_trace_GrowIndex()> or <p7_trace_GrowIndexTo()> for
 *            that.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, the data in
 *            <tr> are unaffected.
 */
int
p7_trace_Grow(P7_TRACE *tr)
{
  void *tmp;
  int   status;
  
  if (tr->N < tr->nalloc) return eslOK;

  ESL_RALLOC(tr->st, tmp, sizeof(char) *2*tr->nalloc);
  ESL_RALLOC(tr->k,  tmp, sizeof(int)  *2*tr->nalloc);
  ESL_RALLOC(tr->i,  tmp, sizeof(int)  *2*tr->nalloc);
  tr->nalloc *= 2;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_trace_GrowIndex()
 * Synopsis:  Grows the allocation of the trace's domain index.
 * Incept:    SRE, Fri Jan  4 10:40:02 2008 [Janelia]
 *
 * Purpose:   If <tr> can't fit another domain in its index,
 *            double the allocation of the index in <tr>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, the 
 *            data in <tr> are unaffected.
 */
int
p7_trace_GrowIndex(P7_TRACE *tr)
{
  void *p;
  int   status;

  if (tr->ndom < tr->ndomalloc) return eslOK;

  ESL_RALLOC(tr->tfrom,   p, sizeof(int)*2*tr->ndomalloc);
  ESL_RALLOC(tr->tto,     p, sizeof(int)*2*tr->ndomalloc);
  ESL_RALLOC(tr->sqfrom,  p, sizeof(int)*2*tr->ndomalloc);
  ESL_RALLOC(tr->sqto,    p, sizeof(int)*2*tr->ndomalloc);
  ESL_RALLOC(tr->hmmfrom, p, sizeof(int)*2*tr->ndomalloc);
  ESL_RALLOC(tr->hmmto,   p, sizeof(int)*2*tr->ndomalloc);
  tr->ndomalloc *= 2;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_trace_GrowTo()
 * Synopsis:  Reallocates trace to a given minimum size.
 *
 * Purpose:   Reallocates a trace structure <tr> to hold a trace
 *            of at least length <N> states.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; the data in <tr>
 *            are unaffected by failure.
 */
int
p7_trace_GrowTo(P7_TRACE *tr, int N)
{
  int status;
  void *tmp;

  if (N < tr->nalloc) return eslOK; /* no-op */
  
  ESL_RALLOC(tr->st, tmp, sizeof(char) *N);
  ESL_RALLOC(tr->k,  tmp, sizeof(int)  *N);
  ESL_RALLOC(tr->i,  tmp, sizeof(int)  *N);
  tr->nalloc = N;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_trace_GrowIndexTo()
 * Synopsis:  Reallocates domain index for a given minimum number.
 * Incept:    SRE, Fri Jan  4 10:47:43 2008 [Janelia]
 *
 * Purpose:   Reallocates the domain index in <tr> to index
 *            at least <ndom> domains.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, in which case
 *            the data in <tr> are unaffected.
 */
int
p7_trace_GrowIndexTo(P7_TRACE *tr, int ndom)
{
  void *p;
  int   status;

  if (ndom < tr->ndomalloc) return eslOK;

  ESL_RALLOC(tr->tfrom,   p, sizeof(int)*ndom);
  ESL_RALLOC(tr->tto,     p, sizeof(int)*ndom);
  ESL_RALLOC(tr->sqfrom,  p, sizeof(int)*ndom);
  ESL_RALLOC(tr->sqto,    p, sizeof(int)*ndom);
  ESL_RALLOC(tr->hmmfrom, p, sizeof(int)*ndom);
  ESL_RALLOC(tr->hmmto,   p, sizeof(int)*ndom);
  tr->ndomalloc = ndom;
  return eslOK;
  
 ERROR:
  return status;
}


/* Function:  p7_trace_Destroy()
 * Synopsis:  Frees a trace.
 *
 * Purpose:   Frees a trace structure <tr>.
 *
 * Returns:   (void)
 */
void 
p7_trace_Destroy(P7_TRACE *tr)
{
  if (tr == NULL) return;
  if (tr->st      != NULL) free(tr->st);
  if (tr->k       != NULL) free(tr->k);
  if (tr->i       != NULL) free(tr->i);
  if (tr->tfrom   != NULL) free(tr->tfrom);
  if (tr->tto     != NULL) free(tr->tto);
  if (tr->sqfrom  != NULL) free(tr->sqfrom);
  if (tr->sqto    != NULL) free(tr->sqto);
  if (tr->hmmfrom != NULL) free(tr->hmmfrom);
  if (tr->hmmto   != NULL) free(tr->hmmto);
  free(tr);
  return;
}


/* Function:  p7_trace_DestroyArray()
 *
 * Purpose:   Frees an array of <N> trace structures, <tr[0..N-1]>.
 *
 * Returns:   (void)
 */
void 
p7_trace_DestroyArray(P7_TRACE **tr, int N)
{
  int idx;

  if (tr == NULL) return;
  for (idx = 0; idx < N; idx++)
    {
      if (tr[idx] == NULL) continue;
      p7_trace_Destroy(tr[idx]);
    }
  free(tr);
  return;
}

/* Function:  p7_trace_Validate()
 * Incept:    SRE, Fri Jan  5 09:17:24 2007 [Janelia]
 *
 * Purpose:   Validate the internal data in a trace structure <tr>
 *            representing an alignment of an HMM to a 
 *            digital sequence <sq>. The digital sequence may be either
 *            unaligned (usually) or aligned (in the case of "fake"
 *            tracebacks generated from an MSA during a
 *            model construction process). 
 *            
 *            We don't pass the HMM that the trace is associated with,
 *            because we might have constructed the trace during
 *            HMM construction when we don't have an HMM yet; but 
 *            we always have a digital sequence.
 *
 *            Intended for debugging, development, and testing
 *            purposes.
 *            
 * Args:      tr     - trace to validate
 *            abc    - alphabet corresponding to sequence <sq>
 *            sq     - digital sequence that <tr> is explaining
 *            errbuf - NULL, or an error message buffer allocated
 *                     for at least eslERRBUFSIZE chars.           
 *
 * Returns:   <eslOK> if trace appears fine.
 *            Returns <eslFAIL> if a problem is detected; if <errbuf> is
 *            provided (non-<NULL>), an informative message is formatted
 *            there to indicate the reason for the failure.
 */
int
p7_trace_Validate(P7_TRACE *tr, ESL_ALPHABET *abc, ESL_DSQ *dsq, char *errbuf)
{
  int  status;
  int  z;			/* position in trace    */
  int  i;			/* position in sequence */
  int  k;			/* position in model */
  char prv;			/* type of the previous state */
  int  is_core;			/* TRUE if trace is a core trace, not profile */

  /* minimum trace length is a core's B->Mk->E. If we don't have at least that,
   * we're definitely in trouble
   */
  if (tr->N < 3)          ESL_XFAIL(eslFAIL, errbuf, "trace is too short");
  if (tr->N > tr->nalloc) ESL_XFAIL(eslFAIL, errbuf, "N of %d isn't sensible", tr->N);

  /* Determine if this is a core trace or a profile trace, so we can
   * construct validation tests appropriately.
   */
  if      (tr->st[0] == p7T_B) is_core = TRUE;
  else if (tr->st[0] == p7T_S) is_core = FALSE;
  else    ESL_XFAIL(eslFAIL, errbuf, "first state neither S nor B");

  /* Verify "sentinels", the final states of the trace
   * (before we start looking backwards and forwards from each state in 
   * our main validation loop)
   */
  if (is_core  && tr->st[tr->N-1] != p7T_E) ESL_XFAIL(eslFAIL, errbuf, "last state not E");
  if (!is_core && tr->st[tr->N-1] != p7T_T) ESL_XFAIL(eslFAIL, errbuf, "last state not T");
  if (tr->k[0]        != 0)                  ESL_XFAIL(eslFAIL, errbuf, "first state shouldn't have k set");
  if (tr->i[0]        != 0)                  ESL_XFAIL(eslFAIL, errbuf, "first state shouldn't have i set");
  if (tr->k[tr->N-1]  != 0)                  ESL_XFAIL(eslFAIL, errbuf, "last state shouldn't have k set");
  if (tr->i[tr->N-1]  != 0)                  ESL_XFAIL(eslFAIL, errbuf, "last state shouldn't have i set");

  /* Main validation loop.
   */
  k = 0; 
  i = 1;
  for (z = 1; z < tr->N-1; z++)
    {
      for (; dsq[i] != eslDSQ_SENTINEL; i++) /* find next non-gap residue in dsq */
	if (esl_abc_XIsResidue(abc, dsq[i])) break;

      /* Handle missing data states: can only be one.
       * prv state might have to skip over one (but not more) missing data states
       */
      if (tr->st[z]   == p7T_X) continue; /* the main tests below will catch X->X cases. */
      prv = (tr->st[z-1] == p7T_X)? tr->st[z-2] : tr->st[z-1];

      switch (tr->st[z]) {
      case p7T_S:
	ESL_XFAIL(eslFAIL, errbuf, "S must be first state");
	break;
	
      case p7T_N:
	if (is_core)          ESL_XFAIL(eslFAIL, errbuf, "core trace can't contain N");
	if (tr->k[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "no N should have k set");
	if (prv == p7T_S) { /* 1st N doesn't emit */
	  if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "first N shouldn't have i set");
	} else if (prv == p7T_N) { /* subsequent N's do */
	  if (tr->i[z] != i) ESL_XFAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	  i++;
	} else ESL_XFAIL(eslFAIL, errbuf, "bad transition to N; expected {S,N}->N");
	break;

      case p7T_B:
	if (tr->k[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "B shouldn't have k set");
	if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "B shouldn't have i set");
	if (prv != p7T_N && prv != p7T_J) 
	  ESL_XFAIL(eslFAIL, errbuf, "bad transition to B; expected {N,J}->B");
	break;

      case p7T_M:
	if (! is_core && prv == p7T_B) k = tr->k[z]; else k++; /* on a B->Mk entry, trust k; else verify */

	if (tr->k[z] != k) ESL_XFAIL(eslFAIL, errbuf, "expected k doesn't match trace's k");
	if (tr->i[z] != i) ESL_XFAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	if (prv != p7T_B && prv != p7T_M && prv != p7T_D && prv != p7T_I)
	  ESL_XFAIL(eslFAIL, errbuf, "bad transition to M; expected {B,M,D,I}->M");
	i++;
	break;

      case p7T_D:
	k++;
	if (tr->k[z] != k) ESL_XFAIL(eslFAIL, errbuf, "expected k doesn't match trace's k");
	if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "D shouldn't have i set");
	if (is_core) {
	  if (prv != p7T_M && prv != p7T_D && prv != p7T_B)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to D; expected {B,M,D}->D");
	} else {
	  if (prv != p7T_M && prv != p7T_D)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to D; expected {M,D}->D");
	}
	break;
	
      case p7T_I:
	if (tr->k[z] != k) ESL_XFAIL(eslFAIL, errbuf, "expected k doesn't match trace's k");
	if (tr->i[z] != i) ESL_XFAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	if (is_core) {
	  if (prv != p7T_B && prv != p7T_M && prv != p7T_I)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to I; expected {B,M,I}->I");
	} else {
	  if (prv != p7T_M && prv != p7T_I)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to I; expected {M,I}->I");
	}
	i++;
	break;

      case p7T_E:
	if (tr->k[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "E shouldn't have k set");
	if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "E shouldn't have i set");
	if (is_core) {
	  if (prv != p7T_M && prv != p7T_D && prv != p7T_I)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to E; expected {M,D,I}->E");
	} else {
	  if (prv != p7T_M && prv != p7T_D)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to E; expected {M,D}->E");
	}
	break;
	
      case p7T_J:
	if (tr->k[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "no J should have k set");
	if (prv == p7T_E) { /* 1st J doesn't emit */
	  if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "first J shouldn't have i set");
	} else if (prv == p7T_J) { /* subsequent J's do */
	  if (tr->i[z] != i) ESL_XFAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	  i++;
	} else ESL_XFAIL(eslFAIL, errbuf, "bad transition to J; expected {E,J}->J");
	break;

      case p7T_C:
	if (is_core)          ESL_XFAIL(eslFAIL, errbuf, "core trace can't contain C");
	if (tr->k[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "no C should have k set");
	if (prv == p7T_E) { /* 1st C doesn't emit */
	  if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "first C shouldn't have i set");
	} else if (prv == p7T_C) { /* subsequent C's do */
	  if (tr->i[z] != i) ESL_XFAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	  i++;
	} else ESL_XFAIL(eslFAIL, errbuf, "bad transition to C; expected {E,C}->C");
	break;
	
      case p7T_T:
	ESL_XFAIL(eslFAIL, errbuf, "T must be last state");
	break;	
      }
    }

  /* Trace should have accounted for all residues in the dsq
   */
  for (; dsq[i] != eslDSQ_SENTINEL; i++) 
    if (esl_abc_XIsResidue(abc, dsq[i])) 
      ESL_XFAIL(eslFAIL, errbuf, "trace didn't account for all residues in the sq");

  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_trace_Dump()
 * Incept:    SRE, Fri Jan  5 09:26:04 2007 [Janelia]
 *
 * Purpose:   Dumps internals of a traceback structure <tr> to <fp>.
 *            If <gm> is non-NULL, also prints transition/emission scores.
 *            If <dsq> is non-NULL, also prints residues (using alphabet
 *            in the <gm>).
 *            
 * Args:      fp   - stream to dump to (often stdout)
 *            tr   - trace to dump
 *            gm   - NULL, or score profile corresponding to trace
 *            dsq  - NULL, or digitized seq corresponding to trace        
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEINVAL> if trace contains something corrupt or invalid;
 *            in this case, dump will be aborted, possibly after partial
 *            output.
 */
int
p7_trace_Dump(FILE *fp, P7_TRACE *tr, P7_PROFILE *gm, ESL_DSQ *dsq) /* replace void w/ P7_PROFILE */
{
  int j;		/* counter for trace position */

  if (tr == NULL) { fprintf(fp, " [ trace is NULL ]\n"); return eslOK; }

  if (gm == NULL) 
    {
      fprintf(fp, "st   k      i   - traceback len %d\n", tr->N);
      fprintf(fp, "--  ----   ----\n");
      for (j = 0; j < tr->N; j++) {
	fprintf(fp, "%1s  %4d %6d\n", 
		p7_hmm_DescribeStatetype(tr->st[j]),
		tr->k[j],
		tr->i[j]);
      } 
    } 
  else 
    {
      int   status;
      float sc = 0.0f;
      float tsc;
      int   xi;


      fprintf(fp, "st   k     i      transit emission - traceback len %d\n", tr->N);
      fprintf(fp, "--  ---- ------  -------- --------\n");
      for (j = 0; j < tr->N; j++) 
	{
	  if (j < tr->N-1) 
	    {
	      status = p7_profile_GetT(gm, tr->st[j], tr->k[j], tr->st[j+1], tr->k[j+1], &tsc);
	      if (status != eslOK) return status;
	    }
	  else tsc = 0.0f;

	  fprintf(fp, "%1s  %4d %6d  %8.4f", p7_hmm_DescribeStatetype(tr->st[j]),  tr->k[j], tr->i[j], tsc);
	  sc += tsc;
	  
	  if (dsq != NULL) {
	    xi = dsq[tr->i[j]];

	    if (tr->st[j] == p7T_M) {
	      fprintf(fp, " %8.4f %c", p7P_MSC(gm, tr->k[j], xi), gm->abc->sym[xi]);
	      sc += p7P_MSC(gm, tr->k[j], xi);
	    } 
	    else if (tr->st[j] == p7T_I) {
	      fprintf(fp, " %8.4f %c", p7P_ISC(gm, tr->k[j], xi), (char) tolower((int) gm->abc->sym[xi]));
	      sc += p7P_ISC(gm, tr->k[j], xi);
	    }
	    else if ((tr->st[j] == p7T_N && tr->st[j-1] == p7T_N) ||
		     (tr->st[j] == p7T_C && tr->st[j-1] == p7T_C) ||
		     (tr->st[j] == p7T_J && tr->st[j-1] == p7T_J))  {
	      fprintf(fp, " %8d %c", 0, (char) tolower((int) gm->abc->sym[xi]));
	    }
	  } 
	  else fprintf(fp, " %8s %c", "-", '-');
	  fputs("\n", fp);
	}
      fprintf(fp, "                 -------- --------\n");
      fprintf(fp, "         total: %8.4f\n\n", sc);
    }

  return eslOK;
}


/* Function:  p7_trace_Append()
 * Synopsis:  Add an element (state/residue) to a growing trace.
 *
 * Purpose:   Adds an element to a trace <tr> that is growing
 *            left-to-right. The element is defined by a state type
 *            <st> (such as <p7T_M>); a node index <k> (1..M for
 *            M,D,I main states; else 0); and a dsq position <i> (1..L
 *            for emitters, else 0).
 *            
 *            For CNJ states, which emit on transition, by convention
 *            we associate the emission with the downstream state; therefore
 *            the first state in any run of CNJ states has i=0. 
 *            
 *            Reallocates the trace (by doubling) if necessary.
 *            
 *            Caller can grow a trace right-to-left too, if it
 *            plans to call <p7_trace_Reverse()>. 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure. The element is successfully
 *            added, but no more elements can be added before this trace is
 *            destroyed.
 *            
 *            <eslEINVAL> if you try to add an element to a trace whose
 *            reallocation has already failed.
 */
int
p7_trace_Append(P7_TRACE *tr, char st, int k, int i)
{
  int status;

  if ((status = p7_trace_Grow(tr)) != eslOK) return status;

  switch (st) {
    /* Emit-on-transition states: */
  case p7T_N: 
  case p7T_C: 
  case p7T_J: if (tr->st[tr->N-1] == st) tr->i[tr->N-1] = i; tr->i[tr->N] = 0; tr->k[tr->N] = 0; break;
    /* Nonemitting states, outside main model: */
  case p7T_X:
  case p7T_S:
  case p7T_B:
  case p7T_E:
  case p7T_T: tr->i[tr->N] = 0; tr->k[tr->N] = 0; break;
    /* Nonemitting, but in main model (k valid) */
  case p7T_D: tr->i[tr->N] = 0; tr->k[tr->N] = k; break;
    /* Emitting states, with valid k position in model: */
  case p7T_M: 
  case p7T_I: tr->i[tr->N] = i; tr->k[tr->N] = k; break;
  default:    ESL_EXCEPTION(eslEINVAL, "no such state; can't append");
  }

  tr->st[tr->N] = st;
  tr->N++;
  return eslOK;
}


/* Function: p7_trace_Reverse()
 * Synopsis: Reverse the arrays in a traceback structure.
 * 
 * Purpose:  Reverse the arrays in a traceback structure.  Tracebacks
 *           from DP algorithms are collected backwards, and they call this
 *           function when they're done.
 *           
 *           At least for now, this invalidates any domain index
 *           table, if it exists. The expectd order of invocation is
 *           to create the traceback backwards, <Reverse()> it, then
 *           <IndexDomains()> it.
 *           
 * Args:     tr - the traceback to reverse. tr->N must be set.
 *                
 * Return:   <eslOK> on success; <tr> is modified.
 */                
int
p7_trace_Reverse(P7_TRACE *tr)
{
  int    z;
  int    tmp;

  /* In runs of CCCCC, NNNNN, JJJJJ, the first one doesn't emit (i=0),
   * and the rest do. Reassign if needed, so after reversal they'll be
   * as they should be.
   */
  for (z = 0; z < tr->N; z++)
    {
      if ( (tr->st[z] == p7T_N && tr->st[z+1] == p7T_N) ||
	   (tr->st[z] == p7T_C && tr->st[z+1] == p7T_C) ||
	   (tr->st[z] == p7T_J && tr->st[z+1] == p7T_J))
	{
	  if (tr->i[z] == 0 && tr->i[z+1] > 0) 
	    { 
	      tr->i[z] = tr->i[z+1]; tr->i[z+1] = 0; 
	    }
	}
    }

  /* Reverse the trace in place. */
  for (z = 0; z < tr->N/2; z++)
    {
      tmp = tr->st[tr->N-z-1];  tr->st[tr->N-z-1] = tr->st[z];   tr->st[z] = tmp;
      tmp = tr->k[tr->N-z-1];   tr->k[tr->N-z-1]  = tr->k[z];    tr->k[z]  = tmp;
      tmp = tr->i[tr->N-z-1];   tr->i[tr->N-z-1]  = tr->i[z];    tr->i[z]  = tmp;
    }
  /* don't worry about the middle residue in odd-length N, since we're in-place  */
  return eslOK;
}


/* Function:  p7_trace_Index()
 * Synopsis:  Internally index the domains in a trace.
 * Incept:    SRE, Fri Jan  4 11:12:24 2008 [Janelia]
 *
 * Purpose:   Create an internal index of the domains in <tr>.
 *            This makes calls to <GetDomainCount()> and 
 *            <GetDomainCoords()> more efficient, and it is
 *            a necessary prerequisite for creating alignments
 *            of any individual domains in a multidomain trace with
 *            <p7_alidisplay_Create()>.
 *
 * Returns:   <eslOK> on success. 
 *
 * Throws:    <eslEMEM> on allocation failure, in which case the
 *            data in the trace is still fine, but the domain index
 *            table isn't constructed.
 */
int
p7_trace_Index(P7_TRACE *tr)
{
  int z;
  int status;

  tr->ndom = 0;
  for (z = 0; z < tr->N; z++)
    {
      switch (tr->st[z]) {
      case p7T_B:
	if ((status = p7_trace_GrowIndex(tr)) != eslOK) goto ERROR;
	tr->tfrom[tr->ndom]   = z;
	tr->sqfrom[tr->ndom]  = 0;
	tr->hmmfrom[tr->ndom] = 0;
	break;

      case p7T_M:
	if (tr->sqfrom[tr->ndom]  == 0) tr->sqfrom[tr->ndom]  = tr->i[z];
	if (tr->hmmfrom[tr->ndom] == 0) tr->hmmfrom[tr->ndom] = tr->k[z];
	tr->sqto[tr->ndom]  = tr->i[z];
	tr->hmmto[tr->ndom] = tr->k[z];
	break;

      case p7T_E:
	tr->tto[tr->ndom]   = z;
	tr->ndom++;
	break;
      }
    }
  return eslOK;
  
 ERROR:
  return status;
}



/* Function: p7_trace_Count()
 * 
 * Purpose:  Count a traceback into a count-based core HMM structure.
 *           (Usually as part of a model parameter re-estimation.)
 *           
 *           The traceback may either be a core traceback (as in model
 *           construction) or a profile traceback (as in model
 *           reestimation).
 *           
 *           If it is a profile traceback, we have to be careful how
 *           we translate an internal entry path from a score profile
 *           back to the core model. Sometimes a B->M_k transition is
 *           an internal entry from local alignment, and sometimes it
 *           is a wing-folded B->D_1..DDM_k alignment to the core
 *           model.
 *           
 *           This is one of the purposes of the special p7T_X
 *           'missing data' state in tracebacks. Local alignment entry
 *           is indicated by a B->X->M_k 'missing data' path, and
 *           direct B->M_k or M_k->E transitions in a traceback are
 *           interpreted as wing retraction in a glocal model.
 * 
 *           The STX state is also used in core traces in model
 *           construction literally to mean missing data, in the
 *           treatment of sequence fragments.
 *
 * Args:     hmm   - counts-based HMM to count <tr> into
 *           tr    - alignment of seq to HMM
 *           dsq   - digitized sequence that traceback aligns to the HMM (1..L)
 *                   (or can be an ax, aligned digital seq)
 *           wt    - weight on this sequence
 *           
 * Return:   <eslOK> on success.
 *           Weighted count events are accumulated in hmm's mat[][], ins[][],
 *           t[][] fields: the core probability model.
 *           
 * Throws:   <eslEINVAL> if something's corrupt in the trace; effect on hmm
 *           counts is undefined, because it may abort at any point in the trace.
 */
int
p7_trace_Count(P7_HMM *hmm, ESL_DSQ *dsq, float wt, P7_TRACE *tr)
{
  int z;                     /* position in tr */
  int i;			/* symbol position in seq */
  int st,st2;     		/* state type (cur, nxt) */
  int k,k2,ktmp;		/* node index (cur, nxt) */
  
  for (z = 0; z < tr->N-1; z++) 
    {
      if (tr->st[z] == p7T_X) continue; /* skip missing data */

      /* pull some info into tmp vars for notational clarity later. */
      st  = tr->st[z]; st2 = tr->st[z+1];
      k   = tr->k[z];  k2  = tr->k[z+1];
      i   = tr->i[z];

      /* Emission counts. */
      if      (st == p7T_M) esl_abc_FCount(hmm->abc, hmm->mat[k], dsq[i], wt);
      else if (st == p7T_I) esl_abc_FCount(hmm->abc, hmm->ins[k], dsq[i], wt);

      /* Transition counts */
      if (st2 == p7T_X) continue; /* ignore transition to missing data */

      if (st == p7T_B) {
	if (st2 == p7T_M && k2 > 1)   /* wing-retracted B->DD->Mk path */
	  {
	    hmm->t[0][p7H_MD] += wt;                
	    for (ktmp = 1; ktmp < k2-1; ktmp++) 
	      hmm->t[ktmp][p7H_DD] += wt;
	    hmm->t[ktmp][p7H_DM] += wt;
	  }
	else  {
	  switch (st2) {
	  case p7T_M: hmm->t[0][p7H_MM] += wt; break;
	  case p7T_I: hmm->t[0][p7H_MI] += wt; break;
	  case p7T_D: hmm->t[0][p7H_MD] += wt; break;
	  default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	  }
	}
      }
      else if (st == p7T_M) {
     	switch (st2) {
	case p7T_M: hmm->t[k][p7H_MM] += wt; break;
	case p7T_I: hmm->t[k][p7H_MI] += wt; break;
	case p7T_D: hmm->t[k][p7H_MD] += wt; break;
	case p7T_E: hmm->t[k][p7H_MM] += wt; break; /* k==M. A local alignment would've been Mk->X->E. */
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
      else if (st == p7T_I) {
	switch (st2) {
	case p7T_M: hmm->t[k][p7H_IM] += wt; break;
	case p7T_I: hmm->t[k][p7H_II] += wt; break;
	case p7T_E: hmm->t[k][p7H_IM] += wt; break; /* k==M. */
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
      else if (st == p7T_D) {
	switch (st2) {
	case p7T_M: hmm->t[k][p7H_DM] += wt; break;
	case p7T_D: hmm->t[k][p7H_DD] += wt; break;
	case p7T_E: hmm->t[k][p7H_DM] += wt; break; /* k==M. A local alignment would've been Dk->X->E. */
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
    } /* end loop over trace position */
  return eslOK;
}

/* Function:  p7_trace_Score()
 * Incept:    SRE, Tue Mar  6 14:40:34 2007 [Janelia]
 *
 * Purpose:   Score path <tr> for digital target sequence <dsq> 
 *            using profile <gm>. Return the lod score in
 *            <ret_sc>.
 *
 * Args:      tr     - traceback path to score
 *            dsq    - digitized sequence
 *            gm     - score profile
 *            ret_sc - RETURN: lod score of trace <tr>
 *
 * Returns:   <eslOK> on success, and <*ret_sc> contains the
 *            lod score for the trace.
 *
 * Throws:    <eslEINVAL> if something's wrong with the trace.
 *            Now <*ret_sc> is returned as $-\infty$.
 */
int 
p7_trace_Score(P7_TRACE *tr, ESL_DSQ *dsq, P7_PROFILE *gm, float *ret_sc)
{
  float  sc;		/* total lod score   */
  float tsc;		/* a transition score */
  int    z;             /* position in tr */
  int    xi;		/* digitized symbol in dsq */
  int  status;

  sc = 0.0f;
  for (z = 0; z < tr->N-1; z++) {
    xi = dsq[tr->i[z]];

    if      (tr->st[z] == p7T_M) sc += p7P_MSC(gm, tr->k[z], xi);
    else if (tr->st[z] == p7T_I) sc += p7P_ISC(gm, tr->k[z], xi);

    if ((status = p7_profile_GetT(gm, tr->st[z], tr->k[z], 
				  tr->st[z+1], tr->k[z+1], &tsc)) != eslOK) goto ERROR;
    sc += tsc;
  }

  *ret_sc = sc;
  return eslOK;

 ERROR:
  *ret_sc = -eslINFINITY;
  return status;
}


/* Function:  p7_trace_GetDomainCount()
 * Incept:    SRE, Tue Feb 27 13:11:43 2007 [Janelia]
 *
 * Purpose:   Determine the number of hits in the trace <tr> -- that is,
 *            the number of times the trace traverses the model from
 *            B...E.  Return that number in <ret_ndom>.
 *            
 *            Done simply by counting the number of B states used in
 *            the trace.
 *            
 *            Only sensible on profile traces. Core traces have 1
 *            domain by definition.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_trace_GetDomainCount(P7_TRACE *tr, int *ret_ndom)
{
  int z;
  int ndom = 0;

  if (tr->ndom > 0) 
    ndom = tr->ndom; /* if we already indexed the domains, we know the answer */
  else {
    for (z = 0; z < tr->N; z++)
      if (tr->st[z] == p7T_B) ndom++;
  }
  *ret_ndom = ndom;
  return eslOK;
}

/* Function:  p7_trace_StateUseCounts()
 * Incept:    SRE, Sun May 27 10:30:13 2007 [Janelia]
 *
 * Purpose:   Accumulate counts of each different state type in trace <tr>. 
 *
 *            <counts[]> is allocated for at least <p7T_NSTATETYPES>
 *            integers, indexed by statetype. Upon return,
 *            <counts[p7T_M]> contains the number of match states
 *            in the trace, for example.
 */
int
p7_trace_StateUseCounts(const P7_TRACE *tr, int *counts)
{
  int x,z;

  for (x = 0; x < p7T_NSTATETYPES; x++) counts[x] = 0;

  for (z = 0; z < tr->N; z++) {
    x = tr->st[z];
    if (x < 0 || x >= p7T_NSTATETYPES) ESL_EXCEPTION(eslEINVAL, "bad state type");
    counts[x]++;
  }
  return eslOK;
}

/* Function:  p7_trace_GetDomainCoords()
 * Incept:    SRE, Tue Feb 27 13:08:32 2007 [Janelia]
 *
 * Purpose:   Retrieve the bounds of domain alignment number <which> in
 *            traceback <tr>. <which> starts from 0. The total number
 *            of domains in a trace can be obtained from
 *            <p7_trace_GetDomainCount()>, or caller can just loop
 *            an increasing <which> count until <eslEOD> is returned.
 *            
 *            Start/end in the sequence are returned in <ret_i1>,
 *            <ret_i2>. Start/end in the model are returned in <ret_k1>, 
 *            <ret_k2>.
 *
 *            It only makes sense to call this function on profile
 *            traces.
 *            
 *            By local alignment bounds convention, the domain
 *            alignment is defined as bounded by match states, so <k1>
 *            and <k2> are the coords of the first and last match
 *            state (in range 1..M), and <i1> and <i2> are the coords
 *            of the residues aligned to those match states. Profiles
 *            do allow a Mk->DDD->E trailer; nonetheless, if such a
 *            trailer occurs, the k2 coord still refers to the last
 *            match state's coordinate. Note that such trailers would
 *            only occur in generated or sampled paths, not Viterbi
 *            paths; in Viterbi alignments with exit probabilities of
 *            1.0, the direct Mk->E path will always have higher
 *            probability than a Mk->DDD->E path.
 *
 * Returns:   <eslOK> on success, and the coords are returned.
 *            <eslEOD> if the trace doesn't contain a <which>'th
 *            domain, and the coords are all returned as 0.
 *            
 * Throws:    <eslEINVAL> if you stupidly pass a <which> less than 0;
 *            <eslECORRUPT> if something is grievously wrong with <tr>.           
 */
int
p7_trace_GetDomainCoords(P7_TRACE *tr, int which,
			 int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2)
{
  int status;
  int z;

  if (which < 0) ESL_XEXCEPTION(eslEINVAL, "bad which < 0");

  if (tr->ndom) 		/* do we have an index? then this'll be easy */
    {
      if (which >= tr->ndom) { status = eslEOD; goto ERROR; }
      *ret_i1 = tr->sqfrom[which];
      *ret_i2 = tr->sqto[which];
      *ret_k1 = tr->hmmfrom[which];
      *ret_k2 = tr->hmmto[which];
      return eslOK;
    }

  /* else, the hard way.
   * skip z to one state past the which'th B state. 
   */
  for (z = 0; which >= 0 && z < tr->N; z++)
    if (tr->st[z] == p7T_B) which--;
  if (z == tr->N) { status = eslEOD; goto ERROR; }
  
  /* skip to the first M state and record i1,k1: 
   * in a profile trace, this must be the next state.
   */
  if (tr->st[z] != p7T_M) ESL_XEXCEPTION(eslECORRUPT, "not a profile trace?");
  *ret_i1 = tr->i[z];
  *ret_k1 = tr->k[z];

  /* skip to the end E, then look back at the last M, record i2,k2.
   */
  for (; z < tr->N; z++)
    if (tr->st[z] == p7T_E) break;
  if (z == tr->N)         ESL_EXCEPTION(eslECORRUPT, "invalid trace: no E for a B");
  do { z--; } while (tr->st[z] == p7T_D); /* roll back over any D trailer */
  if (tr->st[z] != p7T_M) ESL_EXCEPTION(eslECORRUPT, "invalid trace: no M");
  *ret_i2 = tr->i[z];
  *ret_k2 = tr->k[z];
  return eslOK;

 ERROR:
  *ret_i1 = 0;
  *ret_i2 = 0;
  *ret_k1 = 0;
  *ret_k2 = 0;
  return status;
}
			 

/************************************************************
 * @LICENSE@
 ************************************************************/

