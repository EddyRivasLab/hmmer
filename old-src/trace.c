/* trace.c
 * The Plan 7 traceback data structure, P7_TRACE.
 * 
 * SVN $Id: trace.c 1387 2005-05-13 20:50:29Z eddy $
 */
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <easel.h>

#include "hmmer.h"

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
  *ret_tr      = NULL;
  
  ESL_ALLOC(tr,     sizeof(P7_TRACE));
  tr->st = NULL;
  tr->k  = NULL;
  tr->i  = NULL;
  ESL_ALLOC(tr->st, sizeof(char) * N);
  ESL_ALLOC(tr->k,  sizeof(int)  * N);
  ESL_ALLOC(tr->i,  sizeof(int)  * N);
  tr->N      = 0;
  tr->nalloc = N;
  status = eslOK;
 CLEANEXIT:
  if (status != eslOK) p7_trace_Destroy(tr);
  else                 *ret_tr = tr;
  return eslOK;
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
  int status;
  
  ESL_RALLOC(tr->st, tmp, sizeof(char) *2*tr->nalloc);
  ESL_RALLOC(tr->k,  tmp, sizeof(int)  *2*tr->nalloc);
  ESL_RALLOC(tr->i,  tmp, sizeof(int)  *2*tr->nalloc);
  tr->nalloc *= 2;
  status = eslOK;
 CLEANEXIT:
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
  status = eslOK;
 CLEANEXIT:
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


/* Function:  p7_trace_Dump()
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
p7_trace_Dump(FILE *fp, P7_TRACE *tr, P7_PROFILE *gm, char *dsq)
{
  int j;		/* counter for trace position */
  int xi;		/* digital residue            */
  int sc;		/* total score of the trace   */
  int tsc;		/* one transition score       */
  int status;


  if (tr == NULL) { fprintf(fp, " [ trace is NULL ]\n"); return; }

  if (gm == NULL) {
    fprintf(fp, "st  node   rpos  - traceback len %d\n", tr->N);
    fprintf(fp, "--  ---- ------\n");
    for (j = 0; j < tr->N; j++) {
      fprintf(fp, "%1s  %4d %6d\n", 
	      p7_hmm_Statetype(tr->st[j]),
	      tr->k[j],
	      tr->i[j]);
    } 
  } else {
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
  }
}

/* Function:  p7_trace_Append()
 * Incept:    SRE, Tue Apr 11 17:11:04 2006 [St. Louis]
 *
 * Purpose:   Adds an element to a trace <tr> that is growing
 *            left-to-right. The element is defined by a state type
 *            <st> (such as <p7_STM>); a node index <k> (1..M for
 *            M,D,I main states; else 0); and a seq position <i> (1..L
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
  
  if (tr->N >= tr->nalloc) ESL_ERROR(eslEINVAL, "no space in trace");

  tr->st[tr->N] = st;
  tr->k[tr->N]  = k;
  tr->i[tr->N]  = i;
  tr->N++;
  if (tr->N == tr->nalloc) return p7_trace_Expand(tr);
  return eslOK;
}



/* Function: p7_trace_Reverse()
 * Date:     SRE, Mon Aug 25 12:57:29 1997; Denver CO. 
 * 
 * Purpose:  Reverse the arrays in a traceback structure.
 *           Tracebacks from Forward() and Viterbi() are
 *           collected backwards, and call this function
 *           when they're done.
 *           
 *           It's possible to reverse the arrays in place
 *           more efficiently; but the realloc/copy strategy
 *           has the advantage of reallocating the trace
 *           into the right size of memory. (Tracebacks
 *           overallocate.)
 *           
 * Args:     tr - the traceback to reverse. tr->N must be set.
 *                
 * Return:   <eslOK> on success; <tr> is modified.
 * 
 * Throws:   <eslEMEM> on allocation failure; in which case, <tr>
 *           is left unmodified, but you're probably in trouble.
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

  status = eslOK;
 CLEANEXIT:
  if (status != eslOK) 
    {
      if (st != NULL) free(st);
      if (k  != NULL) free(k);
      if (i  != NULL) free(i);
    }
  return status;
}



/* Function: p7_trace_Count()
 * 
 * Purpose:  Count a traceback into a count-based core HMM structure.
 *           (Usually as part of a model parameter re-estimation.)
 *           
 * Note:     A traceback is relative to a score profile (it does not
 *           contain B->D1 and Dm->E transitions) so we have to define
 *           how we translate internal/exit paths from a score profile
 *           back to a profile HMM.
 *           
 *           That is, where are the terminal probabilities in the core
 *           model going to come from, if tracebacks never contain
 *           them?
 * 
 *           In hmmbuild, we want to deal with two different situations.
 *           Some sequences are considered to be fragments, and some
 *           sequences are considered to be full length. For fragments,
 *           we'll ignore internal entry/exit events. For full length
 *           sequences, we'll count internal entry/exits as paths through
 *           D1/Dm.
 *
 *           In DP alignments, including
 *           hmmsearch/hmmpfam/hmmcalibrate, local mode and glocal
 *           mode have disjoint treatments of entry/exit; in glocal
 *           mode, all internal entry/exit is actually a path through
 *           D1/Dm, whereas local mode's uniform fragment length
 *           distribution is imposed on the model in a way that
 *           ignored the paths through D1/Dm.
 *           
 *           Thus we can deal with both situations consistently, by
 *           passing the <mode> argument. For traces in a local mode
 *           (FS,SW) we ignore internal entry/exit; for traces in
 *           glocal mode (LS,S) we count them as BDDMk and MkDDE
 *           paths. That is, B and M->E transitions are only counted
 *           if alignment is glocal.
 *           
 * Args:     hmm   - counts-based HMM to count <tr> into
 *           tr    - alignment of seq to HMM
 *           dsq   - digitized sequence that traceback aligns to the HMM (1..L)
 *           wt    - weight on this sequence
 *           mode  - alignment mode, such as p7_LOCAL
 *           
 * Return:   <eslOK> on success.
 *           Weighted count events are accumulated in hmm's mat[][], ins[][],
 *           t[][] fields: the core probability model.
 *           
 * Throws:   <eslEINVAL> if something's corrupt in the trace; effect on hmm
 *           counts is undefined, because it may abort at any point in the trace.
 */
int
p7_trace_Count(P7_HMM *hmm, char *dsq, float wt, P7_TRACE *tr, int mode)
{
  int tpos;                     /* position in tr */
  int i;			/* symbol position in seq */
  int st,st2,sttmp;		/* state type (cur, nxt) */
  int k,k2,ktmp;		/* node index (cur, nxt) */
  
  for (tpos = 0; tpos < tr->N-1; tpos++) 
    {
      /* pull some info into tmp vars for notational clarity later. */
      st  = tr->st[tpos]; st2 = tr->st[tpos+1];
      k   = tr->k[tpos];  k2  = tr->k[tpos+1];
      i   = tr->i[tpos];

      /* Emission counts. */
      if      (st == p7_STM) esl_abc_FCount(hmm->abc, hmm->mat[k], dsq[i], wt);
      else if (st == p7_STI) esl_abc_FCount(hmm->abc, hmm->ins[k], dsq[i], wt);

      /* Transition counts */
      if (st == p7_STB && (mode == p7_GLOCAL || mode == p7_UNIGLOCAL)) {
	if (k2 == 1) hmm->t[0][p7_TMM] += wt;     /* B->M1 */
	else {                                 /* retracted B->DD->Mk path */
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
	case p7_STE: /* local:ignore; glocal Mm->E: implicit 1.0 */ break;
	default:     ESL_ERROR(eslEINVAL, "bad transition in trace");
	}
      }
      else if (st == p7_STI) {
	switch (st2) {
	case p7_STM: hmm->t[k][p7_TIM] += wt; break;
	case p7_STI: hmm->t[k][p7_TII] += wt; break;
	default:     ESL_ERROR(eslEINVAL, "bad transition in trace");
	}
      }
      else if (st == p7_STM) {
	switch (st2) {
	case p7_STM: hmm->t[k][p7_TDM] += wt; break;
	case p7_STD: hmm->t[k][p7_TDD] += wt; break;
	case p7_STE: /* ignore; must be D_M, and p(D_M->E) = 1.0 */ break;
	default:      ESL_ERROR(eslEINVAL, "bad transition in trace");
	}
      }

    } /* end loop over trace position */
  return eslOK;
}


/* Function:  p7_trace_Score()
 * Incept:    SRE, Thu Apr 13 08:57:55 2006 [St. Louis]
 *
 * Purpose:   Score traceback <tr> using profile <gm>, for digital target
 *            sequence <dsq>; return the SILO score in <ret_sc>.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEINVAL> if there's something corrupt or invalid in the trace;
 *            in this case, <ret_sc> comes back as <p7_IMPOSSIBLE>.           
 */
int
p7_trace_Score(P7_PROFILE *gm, char *dsq, P7_TRACE *tr, int *ret_sc)
{
  int  sc;			/* total SILO score */
  int  tpos;                    /* position in tr */
  char xi;			/* digitized symbol in dsq */
  int  i;			/* an extra tr position index, for lookahead */
  int  tsc;
  int  status = eslOK;

  sc = 0;
  for (tpos = 0; tpos < tr->N-1; tpos++)
    {
      xi = dsq[tr->i[tpos]];

      if      (tr->st[tpos] == p7_STM) 	sc += gm->msc[xi][tr->k[tpos]];
      else if (tr->st[tpos] == p7_STI) 	sc += gm->isc[xi][tr->k[tpos]];

      status = p7_profile_GetTScore(gm, tr->st[tpos], tr->k[tpos], 
				    tr->st[tpos+1], tr->k[tpos+1], &tsc);
      if (status != eslOK) break;
      sc += tsc;
    }

  if (status == eslOK) *ret_sc = sc;
  else                 *ret_sc = p7_IMPOSSIBLE;
  return status;
}



/* Function: p7_trace_DomainCount()
 * 
 * Purpose:  Returns the number of domains in a single Plan7 trace --
 *           i.e. the number of times we traverse the core model --
 *           equivalent to counting the number of begin states in it.
 */
int
p7_trace_DomainCount(P7_TRACE *tr)
{
  int i;
  int ndom = 0;
  
  for (i = 0; i < tr->N; i++)
    if (tr->st[i] == p7_STB) ndom++;
  return ndom;
}




/* Function: p7_trace_Decompose()
 * Date:     Sat Aug 30 11:18:40 1997 (Denver CO)
 * 
 * Purpose:  Decompose a long multi-hit trace into domains - one or more
 *           traces without N,C,J transitions - for consistent
 *           scoring and statistical evaluation of single domain
 *           hits.
 *           
 * Args:     otr    - original trace structure
 *           ret_tr - RETURN: array of simpler traces        
 *           ret_ntr- RETURN: number of traces.
 *           
 * Returns:  <eslOK> on success.
 *           <ret_tr> alloc'ed here; free individuals with <p7_trace_Destroy()>.
 *
 * Throws:   <eslEMEM> on allocation failure; <otr> is untouched, 
 *           <ret_tr> is returned NULL, and <ret_ntr> is returned 0.
 * 
 *           <eslEINVAL> if trace is corrupted.
 */            
int
p7_trace_Decompose(P7_TRACE *otr, P7_TRACE ***ret_tr, int *ret_ntr)
{
  P7_TRACE **tr;        /* array of new traces          */
  int ntr;              /* number of traces             */
  int i,j;		/* position counters in traces  */
  int idx;		/* index over ntr subtraces     */
  int status;

  ntr = p7_trace_DomainCount(otr);
  if (ntr < 1) ESL_ERROR(eslEINVAL, "can't happen");

  ESL_ALLOC(tr, sizeof(P7_TRACE *) * ntr);
  for (idx = 0; idx < ntr; idx++) tr[idx] = NULL;

  for (idx = 0, i = 0; i < otr->N; i++) /* i = position in old trace */
    if (otr->st[i] == p7_STB)	/* start of a domain in old trace */
      {
	for (j = i+1; j < otr->N; j++) /* j = tmp; get length of subtrace */
	  if (otr->st[j] == p7_STE) break;
			/* trace = S-N-(B..E)-C-T : len + 4 : j-i+1 + 4*/
	status = p7_trace_Create(j-i+5, &(tr[idx]));
	if (status != eslOK) goto CLEANEXIT;

	esl_trace_Append(tr[idx], p7_STS, 0, 0);
	esl_trace_Append(tr[idx], p7_STN, 0, 0);
	for (; i < otr->N && otr->st[i] != p7_STE; i++)
	  esl_trace_Append(tr[idx], otr->st[i], otr->k[i], otr->i[i]);
	if (otr->st[i] != p7_STE) ESL_DIE(eslEINVAL, "invalid trace");
	esl_trace_Append(tr[idx], p7_STE, 0, 0);
	esl_trace_Append(tr[idx], p7_STC, 0, 0);
	esl_trace_Append(tr[idx], p7_STJ, 0, 0);
	idx++;
      }

 CLEANEXIT:
  if (status == eslOK) {
    *ret_tr  = tr;
    *ret_ntr = ntr;
  }
  else {
    if (tr != NULL) {
      for (idx = 0; idx < ntr; idx++) 
	if (tr[idx] != NULL) free(tr[idx]);
      free(tr);
    }
    *ret_tr  = NULL;
    *ret_ntr = 0;
  }
  return status;
}


/* Function:  p7_trace_GetDomainCoords()
 * Incept:    SRE, Fri May 13 08:54:43 2005 [St. Louis]
 *
 * Purpose:   Retrieve the bounds of domain alignment number <which> in
 *            a traceback <tr>. Usually, <which>==1 (the first and
 *            only domain), but for a multihit trace, there may be
 *            more than one domain in the trace.
 *            
 *            The total number of domains in a trace can be calculated
 *            by a call to <p7_trace_DomainCount()>. To retrieve coords
 *            for all N domains, you can loop <which> from 1..N.
 *            
 * Args:      tr     - traceback structure, containing N>=1 domains
 *            which  - which domain to get coords for, 1..N
 *            ret_i1 - optRETURN: start point on sequence [1..L] (or NULL)
 *            ret_i2 - optRETURN: end point on sequence   [1..L] (or NULL)
 *            ret_k1 - optRETURN: start point on model    [1..M] (or NULL)
 *            ret_k2 - optRETURN: end point on model      [1..M] (or NULL)
 *            ret_avlen - optRETURN: the average length of the alignment,
 *                        ((k2-k1+1)+(i2-i1+1))/2.
 *                        edge correction of E-value statistics.
 *            
 * Returns:   <eslOK> on success.
 *            <eslEOD> if there isn't a domain <which>.
 *            
 * Throws:    <eslEINVAL> if trace is corrupt.
 */
int
p7_trace_GetDomainCoords(P7_TRACE *tr, int which,
			 int *ret_i1, int *ret_i2, int *ret_k1, int *ret_k2,
			 int *ret_avlen)
{
  int status;
  int i1, i2, k1, k2;
  int tpos;

  /* skip tpos to the which'th B state.
   */
  if (which < 1) { status == eslEOD; goto CLEANEXIT; }
  for (tpos = 0; which > 0 && tpos < tr->N; tpos++)
    if (tr->st[tpos] == p7_STB) which--;
  if (tpos == tr->N) { status = eslEOD; goto CLEANEXIT; }

  /* skip to the first M state and record i1,k1: this must be the next state.
   */
  tpos++;
  if (tr->st[tpos] != p7_STM) ESL_DIE(eslEINVAL, "bad trace");
  k1 = tr->k[tpos];
  i1 = tr->i[tpos];

  /* skip to the end E, then look back at the last M, record i2,k2.
   */
  for (; tpos < tr->N; tpos++)
    if (tr->st[tpos] == p7_STE) break;
  if (tpos == tr->N) ESL_DIE(eslEINVAL, "invalid trace");
  if (tr->st[tpos-1] != p7_STM) ESL_DIE(eslEINVAL, "invalid trace");
  k2 = tr->k[tpos-1];
  i2 = tr->i[tpos-1];
  status = eslOK;

 CLEANEXIT:
  if (status == eslOK) {
    if (ret_k1 != NULL)    *ret_k1    = k1;
    if (ret_i1 != NULL)    *ret_i1    = i1;
    if (ret_k2 != NULL)    *ret_k2    = k2;
    if (ret_i2 != NULL)    *ret_i2    = i2;
    if (ret_avlen != NULL) *ret_avlen = ((k2-k1+1) + (i2-i1+1)) / 2;
  } else {
    if (ret_k1 != NULL)    *ret_k1    = 0;
    if (ret_i1 != NULL)    *ret_i1    = 0;
    if (ret_k2 != NULL)    *ret_k2    = 0;
    if (ret_i2 != NULL)    *ret_i2    = 0;
    if (ret_avlen != NULL) *ret_avlen = 0;
  }
  return status;
}

/************************************************************
 * @LICENSE@
 ************************************************************/

