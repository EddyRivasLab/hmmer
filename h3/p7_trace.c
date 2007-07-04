/* P7_TRACE, the traceback structure.
 *
 * Contents:
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
 *
 * Purpose:   Allocate a traceback of length <N> states. 
 *
 *            <N> only needs to be an initial guess. Routines that
 *            make tracebacks will dynamically grow the trace as
 *            needed.
 *
 * Returns:   a pointer to the new <P7_TRACE> structure on success.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_TRACE *
p7_trace_Create(int N)
{
  P7_TRACE *tr = NULL;
  int       status;

  ESL_ALLOC(tr, sizeof(P7_TRACE));
  tr->st = NULL;
  tr->k  = NULL;
  tr->i  = NULL;
  ESL_ALLOC(tr->st, sizeof(char) * N);
  ESL_ALLOC(tr->k,  sizeof(int)  * N);
  ESL_ALLOC(tr->i,  sizeof(int)  * N);
  tr->N      = 0;
  tr->nalloc = N;
  return tr;

 ERROR:
  if (tr != NULL) p7_trace_Destroy(tr);
  return NULL;
}


/* Function:  p7_trace_Reuse()
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
  /* At least right now, reusing a trace is as simple as: */
  tr->N = 0;
  return eslOK;
}



/* Function:  p7_trace_Grow()
 *
 * Purpose:   Doubles the allocation in a trace structure <tr>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; the data in
 *            <tr> are unaffected by failure.
 */
int
p7_trace_Grow(P7_TRACE *tr)
{
  void *tmp;
  int   status;
  
  ESL_RALLOC(tr->st, tmp, sizeof(char) *2*tr->nalloc);
  ESL_RALLOC(tr->k,  tmp, sizeof(int)  *2*tr->nalloc);
  ESL_RALLOC(tr->i,  tmp, sizeof(int)  *2*tr->nalloc);
  tr->nalloc *= 2;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_trace_GrowTo()
 *
 * Purpose:   Reallocates a trace structure <tr> to hold a trace
 *            of length <N> states.
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


/* Function:  p7_trace_Destroy()
 *
 * Purpose:   Frees a trace structure <tr>.
 *
 * Returns:   (void)
 */
void 
p7_trace_Destroy(P7_TRACE *tr)
{
  if (tr == NULL) return;
  if (tr->st != NULL) free(tr->st);
  if (tr->k  != NULL) free(tr->k);
  if (tr->i  != NULL) free(tr->i);
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
      if (tr[idx]->st != NULL) free(tr[idx]->st);
      if (tr[idx]->k  != NULL) free(tr[idx]->k);
      if (tr[idx]->i  != NULL) free(tr[idx]->i);
      free(tr[idx]);
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
  if      (tr->st[0] == p7_STB) is_core = TRUE;
  else if (tr->st[0] == p7_STS) is_core = FALSE;
  else    ESL_XFAIL(eslFAIL, errbuf, "first state neither S nor B");

  /* Verify "sentinels", the final states of the trace
   * (before we start looking backwards and forwards from each state in 
   * our main validation loop)
   */
  if (is_core  && tr->st[tr->N-1] != p7_STE) ESL_XFAIL(eslFAIL, errbuf, "last state not E");
  if (!is_core && tr->st[tr->N-1] != p7_STT) ESL_XFAIL(eslFAIL, errbuf, "last state not T");
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
      if (tr->st[z]   == p7_STX) continue; /* the main tests below will catch X->X cases. */
      prv = (tr->st[z-1] == p7_STX)? tr->st[z-2] : tr->st[z-1];

      switch (tr->st[z]) {
      case p7_STS:
	ESL_XFAIL(eslFAIL, errbuf, "S must be first state");
	break;
	
      case p7_STN:
	if (is_core)          ESL_XFAIL(eslFAIL, errbuf, "core trace can't contain N");
	if (tr->k[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "no N should have k set");
	if (prv == p7_STS) { /* 1st N doesn't emit */
	  if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "first N shouldn't have i set");
	} else if (prv == p7_STN) { /* subsequent N's do */
	  if (tr->i[z] != i) ESL_XFAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	  i++;
	} else ESL_XFAIL(eslFAIL, errbuf, "bad transition to N; expected {S,N}->N");
	break;

      case p7_STB:
	if (tr->k[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "B shouldn't have k set");
	if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "B shouldn't have i set");
	if (prv != p7_STN && prv != p7_STJ) 
	  ESL_XFAIL(eslFAIL, errbuf, "bad transition to B; expected {N,J}->B");
	break;

      case p7_STM:
	if (! is_core && prv == p7_STB) k = tr->k[z]; else k++; /* on a B->Mk entry, trust k; else verify */

	if (tr->k[z] != k) ESL_XFAIL(eslFAIL, errbuf, "expected k doesn't match trace's k");
	if (tr->i[z] != i) ESL_XFAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	if (prv != p7_STB && prv != p7_STM && prv != p7_STD && prv != p7_STI)
	  ESL_XFAIL(eslFAIL, errbuf, "bad transition to M; expected {B,M,D,I}->M");
	i++;
	break;

      case p7_STD:
	k++;
	if (tr->k[z] != k) ESL_XFAIL(eslFAIL, errbuf, "expected k doesn't match trace's k");
	if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "D shouldn't have i set");
	if (is_core) {
	  if (prv != p7_STM && prv != p7_STD && prv != p7_STB)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to D; expected {B,M,D}->D");
	} else {
	  if (prv != p7_STM && prv != p7_STD)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to D; expected {M,D}->D");
	}
	break;
	
      case p7_STI:
	if (tr->k[z] != k) ESL_XFAIL(eslFAIL, errbuf, "expected k doesn't match trace's k");
	if (tr->i[z] != i) ESL_XFAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	if (is_core) {
	  if (prv != p7_STB && prv != p7_STM && prv != p7_STI)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to I; expected {B,M,I}->I");
	} else {
	  if (prv != p7_STM && prv != p7_STI)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to I; expected {M,I}->I");
	}
	i++;
	break;

      case p7_STE:
	if (tr->k[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "E shouldn't have k set");
	if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "E shouldn't have i set");
	if (is_core) {
	  if (prv != p7_STM && prv != p7_STD && prv != p7_STI)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to E; expected {M,D,I}->E");
	} else {
	  if (prv != p7_STM && prv != p7_STD)
	    ESL_XFAIL(eslFAIL, errbuf, "bad transition to E; expected {M,D}->E");
	}
	break;
	
      case p7_STJ:
	if (tr->k[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "no J should have k set");
	if (prv == p7_STE) { /* 1st J doesn't emit */
	  if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "first J shouldn't have i set");
	} else if (prv == p7_STJ) { /* subsequent J's do */
	  if (tr->i[z] != i) ESL_XFAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	  i++;
	} else ESL_XFAIL(eslFAIL, errbuf, "bad transition to J; expected {E,J}->J");
	break;

      case p7_STC:
	if (is_core)          ESL_XFAIL(eslFAIL, errbuf, "core trace can't contain C");
	if (tr->k[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "no C should have k set");
	if (prv == p7_STE) { /* 1st C doesn't emit */
	  if (tr->i[z] != 0) ESL_XFAIL(eslFAIL, errbuf, "first C shouldn't have i set");
	} else if (prv == p7_STC) { /* subsequent C's do */
	  if (tr->i[z] != i) ESL_XFAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	  i++;
	} else ESL_XFAIL(eslFAIL, errbuf, "bad transition to C; expected {E,C}->C");
	break;
	
      case p7_STT:
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
      int status;
      int tsc;
      int xi;
      int sc = 0;


      fprintf(fp, "st   k     i     transit emission - traceback len %d\n", tr->N);
      fprintf(fp, "--  ---- ------  ------- --------\n");
      for (j = 0; j < tr->N; j++) 
	{
	  if (j < tr->N-1) 
	    {
	      status = p7_profile_GetT(gm, tr->st[j], tr->k[j],
				       tr->st[j+1], tr->k[j+1], &tsc);
	      if (status != eslOK) return status;
	    }
	  else tsc = 0;

	  fprintf(fp, "%1s  %4d %6d  %7d", p7_hmm_DescribeStatetype(tr->st[j]),
		  tr->k[j], tr->i[j], tsc);
	  sc += tsc;
	  
	  if (dsq != NULL) {
	    xi = dsq[tr->i[j]];

	    if (tr->st[j] == p7_STM) {
	      fprintf(fp, " %8d %c", gm->msc[xi][tr->k[j]], gm->abc->sym[xi]);
	      sc += gm->msc[xi][tr->k[j]];
	    } 
	    else if (tr->st[j] == p7_STI) {
	      fprintf(fp, " %8d %c", gm->isc[xi][tr->k[j]], 
		      (char) tolower((int) gm->abc->sym[xi]));
	      sc += gm->isc[xi][tr->k[j]];
	    }
	    else if ((tr->st[j] == p7_STN && tr->st[j-1] == p7_STN) ||
		     (tr->st[j] == p7_STC && tr->st[j-1] == p7_STC) ||
		     (tr->st[j] == p7_STJ && tr->st[j-1] == p7_STJ))  {
	      fprintf(fp, " %8d %c", 0, (char) tolower((int) gm->abc->sym[xi]));
	    }
	  } 
	  else fprintf(fp, " %8s %c", "-", '-');
	  fputs("\n", fp);
	}
      fprintf(fp, "                 ------- --------\n");
      fprintf(fp, "           total: %6d\n\n", sc);
    }

  return eslOK;
}


/* Function:  p7_trace_Append()
 *
 * Purpose:   Adds an element to a trace <tr> that is growing
 *            left-to-right. The element is defined by a state type
 *            <st> (such as <p7_STM>); a node index <k> (1..M for
 *            M,D,I main states; else 0); and a dsq position <i> (1..L
 *            for emitters, else 0).
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
  if (tr->N >= tr->nalloc) ESL_EXCEPTION(eslEINVAL, "no space in trace");

  tr->st[tr->N] = st;
  tr->k[tr->N]  = k;
  tr->i[tr->N]  = i;
  tr->N++;
  if (tr->N == tr->nalloc) return p7_trace_Grow(tr);
  return eslOK;
}


/* Function: p7_trace_Reverse()
 * 
 * Purpose:  Reverse the arrays in a traceback structure.  Tracebacks
 *           from DP algorithms are collected backwards, and call this
 *           function when they're done.
 *           
 *           It's possible to reverse the arrays in place more
 *           efficiently; but the realloc/copy strategy has the
 *           advantage of reallocating the trace into the right size of
 *           memory. (Tracebacks overallocate.)
 *           
 * Args:     tr - the traceback to reverse. tr->N must be set.
 *                
 * Return:   <eslOK> on success; <tr> is modified.
 * 
 * Throws:   <eslEMEM> on allocation failure; in which case, <tr>
 *           is left unmodified, but you're probably already in trouble.
 */                
int
p7_trace_Reverse(P7_TRACE *tr)
{
  int    status;
  char  *st = NULL;
  int   *k  = NULL;
  int   *i  = NULL;
  int    op, np;

  /* Allocate
   */
  ESL_ALLOC(st, sizeof(char)* tr->N);
  ESL_ALLOC(k,  sizeof(int) * tr->N);
  ESL_ALLOC(i,  sizeof(int) * tr->N);
  
  /* Reverse the trace.
   */
  for (op = tr->N-1, np = 0; np < tr->N; np++, op--)
    {
      st[np] = tr->st[op];
      k[np]  = tr->k[op];
      i[np]  = tr->i[op];
    }

  /* Swap old, new arrays.
   */
  free(tr->st); tr->st = st;
  free(tr->k);  tr->k  = k;
  free(tr->i);  tr->i  = i;
  tr->nalloc = tr->N;
  return eslOK;

 ERROR:
  if (st != NULL) free(st);
  if (k  != NULL) free(k);
  if (i  != NULL) free(i);
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
 *           This is one of the purposes of the special p7_STX
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
      if (tr->st[z] == p7_STX) continue; /* skip missing data */

      /* pull some info into tmp vars for notational clarity later. */
      st  = tr->st[z]; st2 = tr->st[z+1];
      k   = tr->k[z];  k2  = tr->k[z+1];
      i   = tr->i[z];

      /* Emission counts. */
      if      (st == p7_STM) esl_abc_FCount(hmm->abc, hmm->mat[k], dsq[i], wt);
      else if (st == p7_STI) esl_abc_FCount(hmm->abc, hmm->ins[k], dsq[i], wt);

      /* Transition counts */
      if (st2 == p7_STX) continue; /* ignore transition to missing data */

      if (st == p7_STB) {
	if (st2 == p7_STM && k2 > 1)   /* wing-retracted B->DD->Mk path */
	  {
	    hmm->t[0][p7_TMD] += wt;                
	    for (ktmp = 1; ktmp < k2-1; ktmp++) 
	      hmm->t[ktmp][p7_TDD] += wt;
	    hmm->t[ktmp][p7_TDM] += wt;
	  }
	else  {
	  switch (st2) {
	  case p7_STM: hmm->t[0][p7_TMM] += wt; break;
	  case p7_STI: hmm->t[0][p7_TMI] += wt; break;
	  case p7_STD: hmm->t[0][p7_TMD] += wt; break;
	  default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	  }
	}
      }
      else if (st == p7_STM) {
     	switch (st2) {
	case p7_STM: hmm->t[k][p7_TMM] += wt; break;
	case p7_STI: hmm->t[k][p7_TMI] += wt; break;
	case p7_STD: hmm->t[k][p7_TMD] += wt; break;
	case p7_STE: hmm->t[k][p7_TMM] += wt; break; /* k==M. A local alignment would've been Mk->X->E. */
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
      else if (st == p7_STI) {
	switch (st2) {
	case p7_STM: hmm->t[k][p7_TIM] += wt; break;
	case p7_STI: hmm->t[k][p7_TII] += wt; break;
	case p7_STE: hmm->t[k][p7_TIM] += wt; break; /* k==M. */
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
      else if (st == p7_STD) {
	switch (st2) {
	case p7_STM: hmm->t[k][p7_TDM] += wt; break;
	case p7_STD: hmm->t[k][p7_TDD] += wt; break;
	case p7_STE: hmm->t[k][p7_TDM] += wt; break; /* k==M. A local alignment would've been Dk->X->E. */
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
 *            using profile <gm>. Return the SILO score in
 *            <ret_sc>.
 *
 * Args:      tr     - traceback path to score
 *            dsq    - digitized sequence
 *            gm     - score profile
 *            ret_sc - RETURN: SILO score of trace <tr>
 *
 * Returns:   <eslOK> on success, and <*ret_sc> contains the
 *            SILO score for the trace.
 *
 * Throws:    <eslEINVAL> if something's wrong with the trace.
 *            Now <*ret_sc> is returned as <p7_IMPOSSIBLE>.
 */
int 
p7_trace_Score(P7_TRACE *tr, ESL_DSQ *dsq, P7_PROFILE *gm, int *ret_sc)
{
  int  sc;		/* total SILO score */
  int  z;               /* position in tr */
  int  xi;		/* digitized symbol in dsq */
  int  tsc;
  int  status;

  sc = 0;
  for (z = 0; z < tr->N-1; z++) {
    xi = dsq[tr->i[z]];

    if      (tr->st[z] == p7_STM) sc += gm->msc[xi][tr->k[z]];
    else if (tr->st[z] == p7_STI) sc += gm->isc[xi][tr->k[z]];

    if ((status = p7_profile_GetT(gm, tr->st[z], tr->k[z], 
				  tr->st[z+1], tr->k[z+1], &tsc)) != eslOK) goto ERROR;
    sc += tsc;
  }

  *ret_sc = sc;
  return eslOK;

 ERROR:
  *ret_sc = p7_IMPOSSIBLE;
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

  for (z = 0; z < tr->N; z++)
    if (tr->st[z] == p7_STB) ndom++;
  *ret_ndom = ndom;
  return eslOK;
}

/* Function:  p7_trace_StateUseCounts()
 * Incept:    SRE, Sun May 27 10:30:13 2007 [Janelia]
 *
 * Purpose:   Accumulate counts of each different state type in trace <tr>. 
 *
 *            <counts[]> is allocated for at least <p7_NSTATETYPES>
 *            integers, indexed by statetype. Upon return,
 *            <counts[p7_STM]> contains the number of match states
 *            in the trace, for example.
 */
int
p7_trace_StateUseCounts(const P7_TRACE *tr, int *counts)
{
  int x,z;

  for (x = 0; x < p7_NSTATETYPES; x++) counts[x] = 0;

  for (z = 0; z < tr->N; z++) {
    x = tr->st[z];
    if (x < 0 || x >= p7_NSTATETYPES) ESL_EXCEPTION(eslEINVAL, "bad state type");
    counts[x]++;
  }
  return eslOK;
}

/* Function:  p7_trace_GetDomainCoords()
 * Incept:    SRE, Tue Feb 27 13:08:32 2007 [Janelia]
 *
 * Purpose:   Retrieve the bounds of domain alignment number <which> in
 *            traceback <tr>. <which> starts from 1. The total number
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
 *            allow a Mk->DDD->E trailer; if such a trailer occurs,
 *            the k2 coord refers to the match state's
 *            coordinate. Note that such trailers would only occur in
 *            generated or sampled paths, not Viterbi paths; in
 *            Viterbi alignments with exit probabilities of 1.0, the
 *            direct Mk->E path will always have higher probability
 *            than a Mk->DDD->E path.
 *
 * Returns:   <eslOK> on success, and the coords are returned.
 *            <eslEOD> if the trace doesn't contain a <which>'th
 *            domain, and the coords are all returned as 0.
 */
int
p7_trace_GetDomainCoords(P7_TRACE *tr, int which, int *ret_i1, int *ret_i2,
			 int *ret_k1, int *ret_k2)
{
  int status;
  int z;

  if (which < 1) ESL_XEXCEPTION(eslEINVAL, "bad which < 1");

  /* skip z to one state past the which'th B state. */
  for (z = 0; which > 0 && z < tr->N; z++)
    if (tr->st[z] == p7_STB) which--;
  if (z == tr->N) { status = eslEOD; goto ERROR; }
  
  /* skip to the first M state and record i1,k1: 
   * in a profile trace, this must be the next state.
   */
  if (tr->st[z] != p7_STM) ESL_XEXCEPTION(eslECORRUPT, "not a profile trace?");
  *ret_i1 = tr->i[z];
  *ret_k1 = tr->k[z];

  /* skip to the end E, then look back at the last M, record i2,k2.
   */
  for (; z < tr->N; z++)
    if (tr->st[z] == p7_STE) break;
  if (z == tr->N)          ESL_EXCEPTION(eslECORRUPT, "invalid trace: no E for a B");
  do { z--; } while (tr->st[z] == p7_STD); /* roll back over any D trailer */
  if (tr->st[z] != p7_STM) ESL_EXCEPTION(eslECORRUPT, "invalid trace: no M");
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

