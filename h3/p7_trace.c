/* Support for P7_TRACE, the traceback structure.
 * 
 * SVN $Id$
 */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"

#include "p7_trace.h"


/* Function:  p7_trace_Create()
 *
 * Purpose:   Allocate a traceback of length <N> states (inclusive
 *            of S, T states); return it via <ret_tr>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <ret_tr> is returned NULL.
 */
int
p7_trace_Create(int N, P7_TRACE **ret_tr)
{
  int       status;
  P7_TRACE *tr = NULL;

  *ret_tr = NULL;
  
  ESL_ALLOC(tr,     sizeof(P7_TRACE));
  tr->st = NULL;
  tr->k  = NULL;
  tr->i  = NULL;
  ESL_ALLOC(tr->st, sizeof(char) * N);
  ESL_ALLOC(tr->k,  sizeof(int)  * N);
  ESL_ALLOC(tr->i,  sizeof(int)  * N);
  tr->N      = 0;
  tr->nalloc = N;

  *ret_tr = tr;
  return eslOK;

 ERROR:
  if (tr != NULL) p7_trace_Destroy(tr);
  return status;
}


/* Function:  p7_trace_Expand()
 *
 * Purpose:   Doubles the allocation in a trace structure <tr>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; the data in
 *            <tr> are unaffected by failure.
 */
int
p7_trace_Expand(P7_TRACE *tr)
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

/* Function:  p7_trace_ExpandTo()
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
p7_trace_ExpandTo(P7_TRACE *tr, int N)
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
  if (tr->N == tr->nalloc) return p7_trace_Expand(tr);
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
 * Note:     A traceback is relative to a profile (it does not contain
 *           B->D1 and Dm->E transitions), so we have to be careful
 *           how we translate internal/exit paths from a score profile
 *           back to the core model. Sometimes a B->M_k or M_k->E 
 *           transition is an internal entry/exit from local alignment,
 *           and sometimes it is a wing-folded B->D_1..DDM_k or
 *           M_k->DDD_M->E global alignment to the core model. 
 *           
 *           This is the main purpose of the special p7_STX 'missing
 *           data' state in tracebacks. Local alignment entry/exits
 *           are indicated by B->X->M_k and M_k->X->E 'missing data'
 *           paths, and direct B->M_k or M_k->E transitions in a
 *           traceback are interpreted as wing foldings.
 * 
 *           In hmmbuild, sequence fragments are treated as local
 *           alignments and dealt with by converting flanking gaps to
 *           'missing data' symbols, and fake traceback construction
 *           then converts a transition from missing data to an M_k to
 *           a B->X->M_k path (or M_k to missing data to M_k -> X ->
 *           E). As a side effect, other missing data would also be
 *           converted to an X.
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
  int tpos;                     /* position in tr */
  int i;			/* symbol position in seq */
  int st,st2,sttmp;		/* state type (cur, nxt) */
  int k,k2,ktmp;		/* node index (cur, nxt) */
  
  for (tpos = 0; tpos < tr->N-1; tpos++) 
    {
      if (tr->st[tpos] == p7_STX) continue; /* skip missing data */

      /* pull some info into tmp vars for notational clarity later. */
      st  = tr->st[tpos]; st2 = tr->st[tpos+1];
      k   = tr->k[tpos];  k2  = tr->k[tpos+1];
      i   = tr->i[tpos];

      /* Emission counts. */
      if      (st == p7_STM) esl_abc_FCount(hmm->abc, hmm->mat[k], dsq[i], wt);
      else if (st == p7_STI) esl_abc_FCount(hmm->abc, hmm->ins[k], dsq[i], wt);

      /* Transition counts */
      if (st2 == p7_STX) continue; /* ignore transition to missing data */

      if (st == p7_STB) {
	if (st2 != p7_STM) ESL_EXCEPTION(eslEINVAL, "can't happen: bad trace");

	if (k2 == 1) hmm->t[0][p7_TMM] += wt;  /* B->M1 */
	else {                                 /* wing-retracted B->DD->Mk path */
	  hmm->t[0][p7_TMD] += wt;                
	  for (ktmp = 1; ktmp < k2-1; k++) 
	    hmm->t[ktmp][p7_TDD] += wt;
	  hmm->t[ktmp][p7_TDM] += wt;
	}
      }
      else if (st == p7_STM) {
     	switch (st2) {
	case p7_STM: hmm->t[k][p7_TMM] += wt; break;
	case p7_STI: hmm->t[k][p7_TMI] += wt; break;
	case p7_STD: hmm->t[k][p7_TMD] += wt; break;
	case p7_STE:            /* M_M->E is a no-op, because hmm->t[M] is implicit, 1.0's to E */
	  if (k < hmm->M) {	/* M_k->E interp'ed as wing-retracted M_k->DD->E path */
	    hmm->t[k][p7_TMD] += wt;
	    for (ktmp = k+1; ktmp < hmm->M; ktmp++)
	      hmm->t[ktmp][p7_TDD] += wt;
	    /* hmm->t[M] implicit 1.0's to E */
	  } 
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
      else if (st == p7_STI) {
	switch (st2) {
	case p7_STM: hmm->t[k][p7_TIM] += wt; break;
	case p7_STI: hmm->t[k][p7_TII] += wt; break;
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
      else if (st == p7_STD) {
	switch (st2) {
	case p7_STM: hmm->t[k][p7_TDM] += wt; break;
	case p7_STD: hmm->t[k][p7_TDD] += wt; break;
	case p7_STE:                          break; /* ignore; must be D_M, and p(D_M->E) is implicitly 1.0 */
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
    } /* end loop over trace position */
  return eslOK;
}

