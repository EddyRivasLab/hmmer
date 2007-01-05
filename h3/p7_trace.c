/* Support for P7_TRACE, the traceback structure.
 * 
 * SRE, Tue Jan  2 2007 [Casa de Gatos] 
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

/* Function:  p7_trace_Validate()
 * Incept:    SRE, Fri Jan  5 09:17:24 2007 [Janelia]
 *
 * Purpose:   Validate the internal data in a trace structure <tr>.
 *            Intended for debugging, development, and testing
 *            purposes.
 *
 * Returns:   <eslOK> if trace appears fine.
 *
 * Throws:    <eslFAIL> if a problem is detected.
 */
int
p7_trace_Validate(P7_TRACE *tr)
{
  int  tpos;			/* position in trace    */
  int  i;			/* position in sequence */
  int  k;			/* position in model */
  char prv;			/* type of the previous state */
  int  status;

  /* Verify that N is sensible (before we start using it in loops!)
   * Minimum is S->N->B->M_k->E->C->T, for a seq of length 1: 7 states.
   */
  if (tr->N < 7 || tr->N > tr->nalloc) ESL_XEXCEPTION(eslFAIL, "N of %d isn't sensible", tr->N);

  /* Verify "sentinels", the S and T at each end of the trace
   * (before we start looking backwards and forwards from each state in 
   * our main validation loop)
   */
  if (tr->st[0]       != p7_STS)   ESL_XEXCEPTION(eslFAIL, "first state not S");
  if (tr->k[0]        != 0)        ESL_XEXCEPTION(eslFAIL, "first state shouldn't have k set");
  if (tr->i[0]        != 0)        ESL_XEXCEPTION(eslFAIL, "first state shouldn't have i set");
  if (tr->st[tr->N-1] != p7_STT)   ESL_XEXCEPTION(eslFAIL, "last state not N");
  if (tr->k[tr->N-1]  != 0)        ESL_XEXCEPTION(eslFAIL, "last state shouldn't have k set");
  if (tr->i[tr->N-1]  != 0)        ESL_XEXCEPTION(eslFAIL, "last state shouldn't have i set");

  /* Main validation loop.
   */
  i = 1;
  k = 0; 
  for (tpos = 1; tpos < tr->N-1; tpos++)
    {
      /* Handle missing data states: can only be one.
       * prv state might have to skip over one (but not more) missing data states
       */
      if (tr->st[tpos]   == p7_STX) continue; /* the main tests below will catch X->X cases. */
      prv = (tr->st[tpos-1] == p7_STX)? tr->st[tpos-2] : tr->st[tpos-1];

      switch (tr->st[tpos]) {
      case p7_STS:
	ESL_XEXCEPTION(eslFAIL, "S must be first state");
	break;
	
      case p7_STN:
	if (tr->k[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "no N should have k set");
	if (prv == p7_STS) { /* 1st N doesn't emit */
	  if (tr->i[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "first N shouldn't have i set");
	} else if (prv == p7_STN) { /* subsequent N's do */
	  if (tr->i[tpos] != i) ESL_XEXCEPTION(eslFAIL, "expected i doesn't match trace's i");
	  i++;
	} else ESL_XEXCEPTION(eslFAIL, "bad transition to N; expected {S,N}->N");
	break;

      case p7_STB:
	if (tr->k[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "B shouldn't have k set");
	if (tr->i[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "B shouldn't have i set");
	if (prv != p7_STN && prv != p7_STJ) 
	  ESL_XEXCEPTION(eslFAIL, "bad transition to B; expected {N,J}->B");
	break;

      case p7_STM:
	if (prv == p7_STB) k = tr->k[tpos]; else k++; /* on a B->Mk entry, trust k; else verify */

	if (tr->k[tpos] != k) ESL_XEXCEPTION(eslFAIL, "expected k doesn't match trace's k");
	if (tr->i[tpos] != i) ESL_XEXCEPTION(eslFAIL, "expected i doesn't match trace's i");
	if (prv != p7_STB && prv != p7_STM && prv != p7_STD && prv != p7_STI)
	  ESL_XEXCEPTION(eslFAIL, "bad transition to M; expected {B,M,D,I}->M");
	i++;
	break;

      case p7_STD:
	k++;
	if (tr->k[tpos] != k) ESL_XEXCEPTION(eslFAIL, "expected k doesn't match trace's k");
	if (tr->i[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "D shouldn't have i set");
	if (prv != p7_STM && prv != p7_STD)
	  ESL_XEXCEPTION(eslFAIL, "bad transition to D; expected {M,D}->D");
	break;
	
      case p7_STI:
	if (tr->k[tpos] != k) ESL_XEXCEPTION(eslFAIL, "expected k doesn't match trace's k");
	if (tr->i[tpos] != i) ESL_XEXCEPTION(eslFAIL, "expected i doesn't match trace's i");
	if (prv != p7_STM && prv != p7_STI)
	  ESL_XEXCEPTION(eslFAIL, "bad transition to I; expected {M,I}->I");
	i++;
	break;

      case p7_STE:
	if (tr->k[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "E shouldn't have k set");
	if (tr->i[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "E shouldn't have i set");
	if (prv != p7_STM)
	  ESL_XEXCEPTION(eslFAIL, "bad transition to E; expected M->E");
	break;
	
      case p7_STJ:
	if (tr->k[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "no J should have k set");
	if (prv == p7_STE) { /* 1st J doesn't emit */
	  if (tr->i[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "first J shouldn't have i set");
	} else if (prv == p7_STJ) { /* subsequent J's do */
	  if (tr->i[tpos] != i) ESL_XEXCEPTION(eslFAIL, "expected i doesn't match trace's i");
	  i++;
	} else ESL_XEXCEPTION(eslFAIL, "bad transition to J; expected {E,J}->J");
	break;

      case p7_STC:
	if (tr->k[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "no C should have k set");
	if (prv == p7_STE) { /* 1st C doesn't emit */
	  if (tr->i[tpos] != 0) ESL_XEXCEPTION(eslFAIL, "first C shouldn't have i set");
	} else if (prv == p7_STC) { /* subsequent C's do */
	  if (tr->i[tpos] != i) ESL_XEXCEPTION(eslFAIL, "expected i doesn't match trace's i");
	  i++;
	} else ESL_XEXCEPTION(eslFAIL, "bad transition to C; expected {E,C}->C");
	break;
	
      case p7_STT:
	ESL_XEXCEPTION(eslFAIL, "T must be last state");
	break;	
      }
    }
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
p7_trace_Dump(FILE *fp, P7_TRACE *tr, void *gm, ESL_DSQ *dsq) /* replace void w/ P7_PROFILE */
{
  int j;		/* counter for trace position */

  if (tr == NULL) { fprintf(fp, " [ trace is NULL ]\n"); return eslOK; }

  if (gm == NULL) {
    fprintf(fp, "st   k      i   - traceback len %d\n", tr->N);
    fprintf(fp, "--  ----   ----\n");
    for (j = 0; j < tr->N; j++) {
      fprintf(fp, "%1s  %4d %6d\n", 
	      p7_hmm_DescribeStatetype(tr->st[j]),
	      tr->k[j],
	      tr->i[j]);
    } 
  } else {
#if 0
    sc = 0;
    fprintf(fp, "st  node   rpos  transit emission - traceback len %d\n", tr->N);
    fprintf(fp, "--  ---- ------  ------- --------\n");
    for (j = 0; j < tr->N; j++) {
      if (j < tr->N-1) {
	status = p7_profile_GetTScore(gm, tr->st[j], tr->k[j],
				      tr->st[j+1], tr->k[j+1], &tsc);
	if (status != eslOK) return status;
      }
      else tsc = 0;

      fprintf(fp, "%1s  %4d %6d  %7d", p7_hmm_Statetype(tr->st[j]),
	      tr->k[j], tr->i[j], tsc);
      sc += tsc;

      if (dsq != NULL) {
	xi = dsq[tr->i[j]];

	if (tr->st[j] == p7_STM) {
	  fprintf(fp, " %8d %c", gm->msc[xi][tr->k[j]], 
		  gm->abc->sym[xi]);
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
#endif
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
  int st,st2;     		/* state type (cur, nxt) */
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
	  for (ktmp = 1; ktmp < k2-1; ktmp++) 
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
	  break; 	
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


/************************************************************
 * @LICENSE@
 ************************************************************/

