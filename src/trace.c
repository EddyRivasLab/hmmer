/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* trace.c
 * SRE, Sat Nov 16 12:34:57 1996
 * RCS $Id$
 * 
 * Support for Plan 7 traceback data structure, p7trace_s.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "structs.h"
#include "config.h"
#include "squid.h"
#include "funcs.h"
#include "version.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static void rightjustify(char *s, int n);

/* Function: P7AllocTrace(), P7ReallocTrace(), P7FreeTrace()
 * 
 * Purpose:  allocation and freeing of traceback structures
 */
void
P7AllocTrace(int tlen, struct p7trace_s **ret_tr)
{
  struct p7trace_s *tr;
  
  tr = (struct p7trace_s *) MallocOrDie (sizeof(struct p7trace_s));
  tr->statetype = (enum p7stype *)  MallocOrDie (sizeof(enum p7stype) * tlen);
  tr->nodeidx   = (int *)   MallocOrDie (sizeof(int)  * tlen);
  tr->pos       = (int *)   MallocOrDie (sizeof(int)  * tlen);
  *ret_tr = tr;
}
void
P7ReallocTrace(struct p7trace_s *tr, int tlen)
{
  tr->statetype = (enum p7stype *) ReallocOrDie (tr->statetype, tlen * sizeof(enum p7stype));
  tr->nodeidx   = (int *)  ReallocOrDie (tr->nodeidx,   tlen * sizeof(int));
  tr->pos       = (int *)  ReallocOrDie (tr->pos,       tlen * sizeof(int));
}
void 
P7FreeTrace(struct p7trace_s *tr)
{
  free(tr->pos);
  free(tr->nodeidx);
  free(tr->statetype);
  free(tr);
}

/* Function: TraceSet()
 * Date:     SRE, Sun Mar  8 12:39:00 1998 [St. Louis]
 *
 * Purpose:  Convenience function; set values at position tpos
 *           in a trace. 
 *           
 *
 * Args:     tr   - trace object to write to
 *           tpos - ptr to position in trace to set
 *           type - statetype e.g. STS, etc.
 *           idx  - nodeidx 1..M or 0
 *           pos  - seq position 1..L or 0
 *
 * Returns:  void
 */
void
TraceSet(struct p7trace_s *tr, int tpos, enum p7stype type, int idx, int pos)
{
  tr->statetype[tpos] = type;
  tr->nodeidx[tpos]   = idx;
  tr->pos[tpos]       = pos;
}


/* Function: MergeTraceArrays()
 * Date:     SRE, Sun Jul  5 15:09:10 1998 [St. Louis]
 *
 * Purpose:  Combine two arrays of traces into a single array.
 *           Used in hmmalign to merge traces from a fixed alignment
 *           with traces from individual unaligned seqs.
 * 
 *           t1 traces always precede t2 traces in the resulting array.
 *
 * Args:     t1 - first set of traces
 *           n1 - number of traces in t1
 *           t2 - second set of traces
 *           n2 - number of traces in t2
 *
 * Returns:  pointer to new array of traces.
 *           Both t1 and t2 are free'd here! Do not reuse.
 */
struct p7trace_s **
MergeTraceArrays(struct p7trace_s **t1, int n1, struct p7trace_s **t2, int n2)
{
  struct p7trace_s **tr;
  int i;			/* index in traces */

  tr = MallocOrDie(sizeof(struct p7trace_s *) * (n1+n2));
  for (i = 0; i < n1; i++) tr[i]    = t1[i];
  for (i = 0; i < n2; i++) tr[n1+i] = t2[i];
  free(t1);
  free(t2);
  return tr;
}



/* Function: P7ReverseTrace()
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
 * Args:     tr - the traceback to reverse. tr->tlen must be set.
 *                
 * Return:   (void)
 *           tr is modified.
 */                
void
P7ReverseTrace(struct p7trace_s *tr)
{
  enum p7stype *statetype;
  int          *nodeidx;
  int          *pos;
  int           opos, npos;

  /* Allocate
   */
  statetype = (enum p7stype *) MallocOrDie (sizeof(enum p7stype) * tr->tlen);
  nodeidx   = (int *) MallocOrDie (sizeof(int) * tr->tlen);
  pos       = (int *) MallocOrDie (sizeof(int) * tr->tlen);
  
  /* Reverse the trace.
   */
  for (opos = tr->tlen-1, npos = 0; npos < tr->tlen; npos++, opos--)
    {
      statetype[npos] = tr->statetype[opos];
      nodeidx[npos]   = tr->nodeidx[opos];
      pos[npos]       = tr->pos[opos];
    }

  /* Swap old, new arrays.
   */
  free(tr->statetype);
  free(tr->nodeidx);
  free(tr->pos);
  tr->statetype = statetype;
  tr->nodeidx   = nodeidx;
  tr->pos       = pos;
}



/* Function: P7TraceCount()
 * 
 * Purpose:  Count a traceback into a count-based HMM structure.
 *           (Usually as part of a model parameter re-estimation.)
 *           
 * Args:     hmm   - counts-based HMM
 *           dsq   - digitized sequence that traceback aligns to the HMM (1..L)
 *           wt    - weight on the sequence
 *           tr    - alignment of seq to HMM
 *           
 * Return:   (void)
 */
void
P7TraceCount(struct plan7_s *hmm, char *dsq, float wt, struct p7trace_s *tr)
{
  int tpos;                     /* position in tr */
  int i;			/* symbol position in seq */
  
  for (tpos = 0; tpos < tr->tlen; tpos++)
    {
      i = tr->pos[tpos];

      /* Emission counts. 
       * Don't bother counting null states N,J,C.
       */
      if (tr->statetype[tpos] == STM) 
	P7CountSymbol(hmm->mat[tr->nodeidx[tpos]], dsq[i], wt);
      else if (tr->statetype[tpos] == STI) 
	P7CountSymbol(hmm->ins[tr->nodeidx[tpos]], dsq[i], wt);

      /* State transition counts
       */
      switch (tr->statetype[tpos]) {
      case STS:
	break;			/* don't bother; P=1 */
      case STN:
	switch (tr->statetype[tpos+1]) {
	case STB: hmm->xt[XTN][MOVE] += wt; break;
	case STN: hmm->xt[XTN][LOOP] += wt; break;
	default:
	  Die("illegal state transition %s->%s in traceback", 
	      Statetype(tr->statetype[tpos]),
	      Statetype(tr->statetype[tpos+1]));
	}
	break;
      case STB:
	switch (tr->statetype[tpos+1]) {
	case STM: hmm->begin[tr->nodeidx[tpos+1]] += wt; break;
	case STD: hmm->tbd1 += wt;                       break;
	default:      
	  Die("illegal state transition %s->%s in traceback", 
	      Statetype(tr->statetype[tpos]),
	      Statetype(tr->statetype[tpos+1]));
	}
	break;
      case STM:
	switch (tr->statetype[tpos+1]) {
	case STM: hmm->t[tr->nodeidx[tpos]][TMM] += wt; break;
	case STI: hmm->t[tr->nodeidx[tpos]][TMI] += wt; break;
	case STD: hmm->t[tr->nodeidx[tpos]][TMD] += wt; break;
	case STE: hmm->end[tr->nodeidx[tpos]]    += wt; break;
	default:    
	  Die("illegal state transition %s->%s in traceback", 
	      Statetype(tr->statetype[tpos]),
	      Statetype(tr->statetype[tpos+1]));
	}
	break;
      case STI:
	switch (tr->statetype[tpos+1]) {
	case STM: hmm->t[tr->nodeidx[tpos]][TIM] += wt; break;
	case STI: hmm->t[tr->nodeidx[tpos]][TII] += wt; break;
	default:    
	  Die("illegal state transition %s->%s in traceback", 
	      Statetype(tr->statetype[tpos]),
	      Statetype(tr->statetype[tpos+1]));
	}
	break;
      case STD:
	switch (tr->statetype[tpos+1]) {
	case STM: hmm->t[tr->nodeidx[tpos]][TDM] += wt; break;
	case STD: hmm->t[tr->nodeidx[tpos]][TDD] += wt; break;
	case STE: /* ignore; p(D->E) = 1.0 */           break;
	default: 
	  Die("illegal state transition %s->%s in traceback", 
	      Statetype(tr->statetype[tpos]),
	      Statetype(tr->statetype[tpos+1]));
	}
	break;
      case STE:
	switch (tr->statetype[tpos+1]) {
	case STC: hmm->xt[XTE][MOVE] += wt; break;
	case STJ: hmm->xt[XTE][LOOP] += wt; break;
	default:     
	  Die("illegal state transition %s->%s in traceback", 
	      Statetype(tr->statetype[tpos]),
	      Statetype(tr->statetype[tpos+1]));
	}
	break;
     case STJ:
	switch (tr->statetype[tpos+1]) {
	case STB: hmm->xt[XTJ][MOVE] += wt; break;
	case STJ: hmm->xt[XTJ][LOOP] += wt; break;
	default:     
	  Die("illegal state transition %s->%s in traceback", 
	      Statetype(tr->statetype[tpos]),
	      Statetype(tr->statetype[tpos+1]));
	}
	break;
      case STC:
	switch (tr->statetype[tpos+1]) {
	case STT: hmm->xt[XTC][MOVE] += wt; break;
	case STC: hmm->xt[XTC][LOOP] += wt; break;
	default:     
	 Die("illegal state transition %s->%s in traceback", 
	     Statetype(tr->statetype[tpos]),
	     Statetype(tr->statetype[tpos+1])); 
	}
	break;
      case STT:
	break;			/* T is the last. It makes no transitions. */
      default:
	Die("illegal state %s in traceback", 
	    Statetype(tr->statetype[tpos]));
      }
    }
}


/* Function: P7TraceScore()
 * 
 * Purpose:  Score a traceback and return the score in scaled bits.
 *           
 * Args:     hmm   - HMM with valid log odds scores.
 *           dsq   - digitized sequence that traceback aligns to the HMM (1..L)
 *           tr    - alignment of seq to HMM
 *           
 * Return:   (void)
 */
float
P7TraceScore(struct plan7_s *hmm, char *dsq, struct p7trace_s *tr)
{
  int score;			/* total score as a scaled integer */
  int tpos;                     /* position in tr */
  int sym;			/* digitized symbol in dsq */
  
  /*  P7PrintTrace(stdout, tr, hmm, dsq); */
  score = 0;
  for (tpos = 0; tpos < tr->tlen-1; tpos++)
    {
      sym = (int) dsq[tr->pos[tpos]];

      /* Emissions.
       * Don't bother counting null states N,J,C.
       */
      if (tr->statetype[tpos] == STM) 
	score += hmm->msc[sym][tr->nodeidx[tpos]];
      else if (tr->statetype[tpos] == STI) 
	score += hmm->isc[sym][tr->nodeidx[tpos]];

      /* State transitions.
       */
      score += TransitionScoreLookup(hmm, 
			     tr->statetype[tpos], tr->nodeidx[tpos],
			     tr->statetype[tpos+1], tr->nodeidx[tpos+1]);
    }
  return Scorify(score);
}



/* Function: P7Traces2Alignment()
 * 
 * Purpose:  Convert an array of traceback structures for a set
 *           of sequences into a new multiple alignment. 
 *           
 *           Insertions are put into lower case and 
 *           are not aligned; instead, Nterm is right-justified,
 *           Cterm is left-justified, and internal insertions
 *           are split in half and the halves are justified in
 *           each direction (the objective being to increase
 *           the chances of getting insertions aligned well enough
 *           for them to become a match). SAM gap char conventions
 *           are used: - in match columns, . in insert columns
 * 
 * NOTE:     Does not recognize J state.
 *
 * Args:     dsq        - digitized unaligned sequences 
 *           sqinfo     - array of info about the sequences
 *           wgt        - weights on seqs
 *           nseq       - number of sequences
 *           mlen       - length of model (number of match states)
 *           tr         - array of tracebacks
 *           matchonly  - TRUE if we don't print insert-generated symbols at all
 *           ret_aseqs  - RETURN: multiple sequence alignment           
 *           ainfo      - RETURN: optional info about alignment 
 *           
 * Return:   1 on success, 0 on failure.
 *           ret_aseqs is alloc'ed here and ainfo is filled in; 
 *           FreeAlignment(aseqs, &ainfo).
 */          
void
P7Traces2Alignment(char **dsq, SQINFO *sqinfo, float *wgt, int nseq, int mlen, 
		   struct p7trace_s **tr, int matchonly, 
		   char ***ret_aseqs, AINFO *ainfo)
{
  char **aseqs;                 /* RETURN: aligned sequence set */
  int    idx;                   /* counter for sequences */
  int    alen;                  /* width of alignment */
  int   *inserts;               /* array of max gaps between aligned columns */
  int   *matmap;                /* matmap[k] = apos of match k [1..M] */
  int    nins;                  /* counter for inserts */
  int    apos;                  /* position in aligned sequence (0..alen-1)*/
  int    rpos;                  /* position in raw digital sequence (1..L)*/
  int    tpos;                  /* position counter in traceback */
  int    statetype;		/* type of current state, e.g. STM */
  int    k;                     /* counter over states in model */

  /* Here's the problem. We want to align the match states in columns,
   * but some sequences have inserted symbols in them; we need some
   * sort of overall knowledge of where the inserts are and how long
   * they are in order to create the alignment.
   * 
   * Here's our trick. inserts[] is a 0..hmm->M array; inserts[i] stores
   * the maximum number of times insert substate i was used. This
   * is the maximum number of gaps to insert between canonical 
   * column i and i+1.  inserts[0] is the N-term tail; inserts[M] is
   * the C-term tail.
   * 
   * Remember that N and C emit on transition, hence the check for an
   * N->N or C->C transition before bumping nins. 
   */
  inserts = (int *) MallocOrDie (sizeof(int) * (mlen+1));
  for (k = 0; k <= mlen; k++)
    inserts[k] = 0;
  for (idx = 0; idx < nseq; idx++) {
    nins = 0;
    for (tpos = 0; tpos < tr[idx]->tlen; tpos++) {
      switch (tr[idx]->statetype[tpos]) { 
      case STI: nins++; break;
      case STN: if (tr[idx]->statetype[tpos-1] == STN) nins++; break;
      case STC: if (tr[idx]->statetype[tpos-1] == STC) nins++; break;
      case STM:
      case STD:		/* M,D: record max. reset ctr. */
	if (nins > inserts[tr[idx]->nodeidx[tpos]-1])
	  inserts[tr[idx]->nodeidx[tpos]-1] = nins;
	nins = 0;
	break;
      case STB:		/* B; record N-tail max, reset ctr */
	if (nins > inserts[0])
	  inserts[0] = nins;
	nins = 0;
	break;
      case STT:		/* T: record C-tail max */
	if (nins > inserts[mlen])
	  inserts[mlen] = nins;
	break;
      case STS: case STE: break; /* ignore other states */
      case STJ:
	Die("yo! you don't support J in Traces2Alignment(), remember?");
      default:
	Die("Traces2Alignment reports unrecognized statetype %c", 
	    Statetype(tr[idx]->statetype[tpos]));
      }
    }
  }

				/* Insert compression option. */
  if (matchonly) 
    for (k = 0; k <= mlen; k++)
      if (inserts[k] > 1) 
	inserts[k] = 1;

  /***********************************************
   * Construct the alignment
   ***********************************************/
				/* calculate alignment length and matmap */
  matmap= (int *)   MallocOrDie (sizeof(int) * (mlen+1));
  matmap[0] = -1;
  alen = inserts[0];
  for (k = 1; k <= mlen ; k++) {
    matmap[k] = alen;
    alen += inserts[k] + 1;
  }
                                /* allocations for new alignment */
  AllocAlignment(nseq, alen, &aseqs, ainfo);

  for (idx = 0; idx < nseq; idx++) {
				/* blank an aseq */
    for (apos = 0; apos < alen; apos++)
      aseqs[idx][apos] = '.';
    for (k = 1; k <= mlen; k++)
      aseqs[idx][matmap[k]] = '-';
    aseqs[idx][alen] = '\0';
				/* align the sequence */
    apos = 0;
    for (tpos = 0; tpos < tr[idx]->tlen; tpos++) {
      statetype = tr[idx]->statetype[tpos]; /* just for clarity */
      rpos      = tr[idx]->pos[tpos]; 
      k         = tr[idx]->nodeidx[tpos];

      if (statetype == STM) {
	apos = matmap[k];
	aseqs[idx][apos] = Alphabet[(int) dsq[idx][rpos]];
	apos++;
      }
      else if (statetype == STI) {
	if (matchonly) 
	  aseqs[idx][apos] = '*'; /* insert compression option */
	else {
	  aseqs[idx][apos] = tolower(Alphabet[(int) dsq[idx][rpos]]);
	  apos++;
	}
      }
      else if ((statetype == STN || statetype == STC) && rpos > 0) {
	if (matchonly)
	  aseqs[idx][apos] = '*'; /* insert compression option */
	else {
	  aseqs[idx][apos] = tolower(Alphabet[(int) dsq[idx][rpos]]);
	  apos++;
	}
      }
      else if (statetype == STE)
	apos = matmap[mlen]+1;	/* set position for C-term tail */
    }

  /* N-terminal extension is right-justified.
   * Internal inserts are split in half, and C-term is right-justified.
   * C-terminal extension remains left-justified.
   */
    if (! matchonly) {
      rightjustify(aseqs[idx], inserts[0]);

      for (k = 1; k < mlen; k++) 
	if (inserts[k] > 1) {
	  for (nins = 0, apos = matmap[k]+1; islower(aseqs[idx][apos]); apos++)
	    nins++;
	  nins /= 2;		/* split the insertion in half */
	  rightjustify(aseqs[idx]+matmap[k]+1+nins, inserts[k]-nins);
	}
    }

  }
    
  /***********************************************
   * Build ainfo
   ***********************************************/
        
  sprintf(ainfo->au, "HMMER %s", RELEASE);
  ainfo->flags |= AINFO_AUTH;
				/* copy sqinfo array and weights */
  for (idx = 0; idx < nseq; idx++)
    {
      SeqinfoCopy(&(ainfo->sqinfo[idx]), &(sqinfo[idx]));
      ainfo->wgt[idx] = wgt[idx];
    }

  /* #=RF annotation: x for match column, . for insert column
   */
  ainfo->flags |= AINFO_RF;
  ainfo->rf = (char *) MallocOrDie (sizeof(char) * (alen+1));
  for (apos = 0; apos < alen; apos++)
    ainfo->rf[apos] = '.';
  for (k = 1; k <= mlen; k++)
    ainfo->rf[matmap[k]] = 'x';
  ainfo->rf[alen] = '\0';

    /* Currently, we produce no consensus structure. 
     * #=CS, generated from HMM structural annotation, would go here.
     */

  free(inserts);
  free(matmap);
  *ret_aseqs = aseqs;
  return;
}

/* Function: TransitionScoreLookup()
 * 
 * Purpose:  Convenience function used in PrintTrace() and TraceScore();
 *           given state types and node indices for a transition,
 *           return the integer score for that transition. 
 */
int
TransitionScoreLookup(struct plan7_s *hmm, enum p7stype st1, int k1, 
		      enum p7stype st2, int k2)
{
  switch (st1) {
  case STS: return 0;	/* S never pays */
  case STN:
    switch (st2) {
    case STB: return hmm->xsc[XTN][MOVE]; 
    case STN: return hmm->xsc[XTN][LOOP]; 
    default:      Die("illegal %s->%s transition", Statetype(st1), Statetype(st2));
    }
    break;
  case STB:
    switch (st2) {
    case STM: return hmm->bsc[k2]; 
    case STD: return Prob2Score(hmm->tbd1, 1.);
    default:      Die("illegal %s->%s transition", Statetype(st1), Statetype(st2));
    }
    break;
  case STM:
    switch (st2) {
    case STM: return hmm->tsc[k1][TMM];
    case STI: return hmm->tsc[k1][TMI];
    case STD: return hmm->tsc[k1][TMD];
    case STE: return hmm->esc[k1];
    default:      Die("illegal %s->%s transition", Statetype(st1), Statetype(st2));
    }
    break;
  case STI:
    switch (st2) {
    case STM: return hmm->tsc[k1][TIM];
    case STI: return hmm->tsc[k1][TII];
    default:      Die("illegal %s->%s transition", Statetype(st1), Statetype(st2));
    }
    break;
  case STD:
    switch (st2) {
    case STM: return hmm->tsc[k1][TDM]; 
    case STD: return hmm->tsc[k1][TDD];
    case STE: return 0;	/* D_m->E has probability 1.0 by definition in Plan7 */
    default:      Die("illegal %s->%s transition", Statetype(st1), Statetype(st2));
    }
    break;
  case STE:
    switch (st2) {
    case STC: return hmm->xsc[XTE][MOVE]; 
    case STJ: return hmm->xsc[XTE][LOOP]; 
    default:      Die("illegal %s->%s transition", Statetype(st1), Statetype(st2));
    }
    break;
  case STJ:
    switch (st2) {
    case STB: return hmm->xsc[XTJ][MOVE]; 
    case STJ: return hmm->xsc[XTJ][LOOP]; 
    default:      Die("illegal %s->%s transition", Statetype(st1), Statetype(st2));
    }
    break;
  case STC:
    switch (st2) {
    case STT: return hmm->xsc[XTC][MOVE]; 
    case STC: return hmm->xsc[XTC][LOOP]; 
    default:      Die("illegal %s->%s transition", Statetype(st1), Statetype(st2));
    }
    break;
  case STT:   return 0;	/* T makes no transitions */
  default:        Die("illegal state %s in traceback", Statetype(st1));
  }
  /*NOTREACHED*/
  return 0;
}


/* Function: CreateFancyAli()
 * Date:     SRE, Mon Oct 27 06:49:44 1997 [Sanger Centre UK]
 * 
 * Purpose:  Output of an HMM/sequence alignment, using a
 *           traceback structure. Deliberately similar to 
 *           the output of BLAST, to make it easier for
 *           people to adapt their Perl parsers (or what have
 *           you) from BLAST to HMMER.
 *           
 * Args:     tr  - traceback structure that gives the alignment      
 *           hmm - the model 
 *           dsq - the sequence (digitized form)       
 *           name- name of the sequence  
 *                 
 * Return:   allocated, filled fancy alignment structure.
 */
struct fancyali_s *
CreateFancyAli(struct p7trace_s *tr, struct plan7_s *hmm,
	       char *dsq, char *name)
{
  struct fancyali_s *ali;       /* alignment to create                */
  int   tpos;			/* position in trace and alignment    */
  int   bestsym;		/* index of best symbol at this pos   */
  float mthresh;		/* above this P(x), display uppercase */

  /* Allocate and initialize the five lines of display
   */
  ali         = AllocFancyAli();
  ali->rfline = NULL;
  ali->csline = NULL;
  ali->model  = (char *) MallocOrDie (sizeof(char) * (tr->tlen+1));
  ali->mline  = (char *) MallocOrDie (sizeof(char) * (tr->tlen+1));
  ali->aseq   = (char *) MallocOrDie (sizeof(char) * (tr->tlen+1));

  memset(ali->model,  ' ', tr->tlen);
  memset(ali->mline,  ' ', tr->tlen);
  memset(ali->aseq,   ' ', tr->tlen);

  if (hmm->flags & PLAN7_RF)
    {
      ali->rfline = (char *) MallocOrDie (sizeof(char) * (tr->tlen+1));
      memset(ali->rfline, ' ', tr->tlen);
    }
  if (hmm->flags & PLAN7_CS)
    {
      ali->csline = (char *) MallocOrDie (sizeof(char) * (tr->tlen+1));
      memset(ali->csline, ' ', tr->tlen);  
    }

  ali->query  = Strdup(hmm->name);
  ali->target = Strdup(name);

  if (Alphabet_type == hmmAMINO) mthresh = 0.5;
  else                           mthresh = 0.9;

  /* Find first, last seq position
   * HMM start/end positions currently not recorded, because there
   * might be multiple HMM hits per sequence.
   */
  for (tpos = 0; tpos < tr->tlen; tpos++) 
    if (tr->pos[tpos] > 0) {
	ali->sqfrom = tr->pos[tpos];
	break;
    }
  for (tpos = tr->tlen-1; tpos >= 0; tpos--)
    if (tr->pos[tpos] > 0) {
      ali->sqto = tr->pos[tpos];
      break;
    }

  /* Fill in the five lines of display
   */
  for (tpos = 0; tpos < tr->tlen; tpos++) {
    switch (tr->statetype[tpos]) {
    case STS: 
    case STT:
      ali->model[tpos] = '*';
      break;

    case STN:
    case STJ:
    case STC:
      ali->model[tpos] = '-';
      if (tr->pos[tpos] > 0) { 
	ali->aseq[tpos] = tolower(Alphabet[(int) dsq[tr->pos[tpos]]]);
      }
      break;

    case STB: 
      ali->model[tpos] = '>';
      break;

    case STE:
      ali->model[tpos] = '<';
      break;

    case STM:
      if (hmm->flags & PLAN7_RF) ali->rfline[tpos] = hmm->rf[tr->nodeidx[tpos]];
      if (hmm->flags & PLAN7_CS) ali->csline[tpos] = hmm->cs[tr->nodeidx[tpos]];
      bestsym = FMax(hmm->mat[tr->nodeidx[tpos]], Alphabet_size);
      ali->model[tpos] = Alphabet[bestsym];
      if (hmm->mat[tr->nodeidx[tpos]][bestsym] < mthresh)
	ali->model[tpos] = tolower(ali->model[tpos]);
      if (dsq[tr->pos[tpos]] == bestsym)
	{
	  ali->mline[tpos] = Alphabet[(int) dsq[tr->pos[tpos]]];
	  if (hmm->mat[tr->nodeidx[tpos]][bestsym] < mthresh)
	    ali->mline[tpos] = tolower(ali->mline[tpos]);
	}
      else if (hmm->msc[(int) dsq[tr->pos[tpos]]] [tr->nodeidx[tpos]] > 0)
	ali->mline[tpos] = '+';
      ali->aseq[tpos]  = Alphabet[(int) dsq[tr->pos[tpos]]];
      break;

    case STD:
      if (hmm->flags & PLAN7_RF) ali->rfline[tpos] = hmm->rf[tr->nodeidx[tpos]];
      if (hmm->flags & PLAN7_CS) ali->csline[tpos] = hmm->cs[tr->nodeidx[tpos]];
      bestsym = FMax(hmm->mat[tr->nodeidx[tpos]], Alphabet_size);
      ali->model[tpos] = Alphabet[bestsym];
      if (hmm->mat[tr->nodeidx[tpos]][bestsym] < mthresh)
	ali->model[tpos] = tolower(ali->model[tpos]);
      ali->aseq[tpos]  = '-';
      break;

    case STI:
      ali->model[tpos] = '.';
      if (hmm->isc[(int) dsq[tr->pos[tpos]]] [tr->nodeidx[tpos]] > 0)
	ali->mline[tpos] = '+';
      ali->aseq[tpos]  = tolower(Alphabet[(int) dsq[tr->pos[tpos]]]);
      break;

    default:
      Die("bogus statetype");
    } /* end switch over statetypes */
  }  /* end loop over tpos */

  ali->len          = tpos;
  if (hmm->flags & PLAN7_RF) ali->rfline[tpos] = '\0';
  if (hmm->flags & PLAN7_CS) ali->csline[tpos] = '\0';
  ali->model[tpos]  = '\0';
  ali->mline[tpos]  = '\0';
  ali->aseq[tpos]   = '\0';
  return ali;
} 


/* Function: PrintFancyAli()
 * Date:     SRE, Mon Oct 27 06:56:42 1997 [Sanger Centre UK]
 * 
 * Purpose:  Print an HMM/sequence alignment from a fancyali_s 
 *           structure. Line length controlled by ALILENGTH in
 *           config.h (set to 50).
 *           
 * Args:     fp  - where to print it (stdout or open FILE)
 *           ali - alignment to print 
 *                 
 * Return:   (void)                
 */
void
PrintFancyAli(FILE *fp, struct fancyali_s *ali)
{
  char  buffer[ALILENGTH+1];	/* output line buffer                 */
  int   starti, endi;
  int   pos;
  int   i;

  buffer[ALILENGTH] = '\0';
  endi = ali->sqfrom - 1;
  for (pos = 0; pos < ali->len; pos += ALILENGTH)
    {
				/* coords of target seq for this line */
      starti = endi + 1;
      for (i = pos; ali->aseq[i] != '\0' && i < pos + ALILENGTH; i++)
	if (!isgap(ali->aseq[i])) endi++;
	    
      if (ali->csline != NULL) {
	strncpy(buffer, ali->csline+pos, ALILENGTH);
	fprintf(fp, "  %16s %s\n", "CS", buffer);
      }
      if (ali->rfline != NULL) {
	strncpy(buffer, ali->rfline+pos, ALILENGTH);
	fprintf(fp, "  %16s %s\n", "RF", buffer);
      }
      if (ali->model  != NULL) {
	strncpy(buffer, ali->model+pos, ALILENGTH);
	fprintf(fp, "  %16s %s\n", " ", buffer);
      }
      if (ali->mline  != NULL) {
	strncpy(buffer, ali->mline+pos, ALILENGTH);
	fprintf(fp, "  %16s %s\n", " ", buffer);
      }
      if (ali->aseq   != NULL) { 
	strncpy(buffer, ali->aseq+pos, ALILENGTH);
	if (endi >= starti)
	  fprintf(fp, "  %10.10s %5d %s %-5d\n\n", ali->target, starti, buffer, endi);
	else
	  fprintf(fp, "  %10.10s %5s %s %-5s\n\n", ali->target, "-", buffer, "-"); 
      }
    }

  /* Cleanup and return
   */
  fflush(fp);
  return;
} 



/* Function: TraceDecompose()
 * Date:     Sat Aug 30 11:18:40 1997 (Denver CO)
 * 
 * Purpose:  Decompose a long multi-hit trace into zero or more
 *           traces without N,C,J transitions: for consistent
 *           scoring and statistical evaluation of single domain
 *           hits.
 *           
 * Args:     otr    - original trace structure
 *           ret_tr - RETURN: array of simpler traces        
 *           ret_ntr- RETURN: number of traces.
 *           
 * Return:   (void)
 *           ret_tr alloc'ed here; free individuals with FreeTrace().
 */            
void
TraceDecompose(struct p7trace_s *otr, struct p7trace_s ***ret_tr, int *ret_ntr)
{
  struct p7trace_s **tr;        /* array of new traces          */
  int ntr;			/* number of traces             */
  int i,j;			/* position counters in traces  */
  int idx;			/* index over ntr subtraces     */

  /* First pass: count begin states to get ntr.
   */
  for (ntr = 0, i = 0; i < otr->tlen; i++)
    if (otr->statetype[i] == STB) ntr++;

  /* Allocations.
   */
  if (ntr == 0) {
    *ret_ntr = 0;
    *ret_tr  = NULL;
    return;
  }
  tr = (struct p7trace_s **) MallocOrDie (sizeof(struct p7trace_s *) * ntr);

  for (idx = 0, i = 0; i < otr->tlen; i++) /* i = position in old trace */
    if (otr->statetype[i] == STB)
      {
	for (j = i+1; j < otr->tlen; j++) /* j = tmp; get length of subtrace */
	  if (otr->statetype[j] == STE) break;
			/* trace = S-N-(B..E)-C-T : len + 4 : j-i+1 + 4*/
	P7AllocTrace(j-i+5, &(tr[idx]));
	tr[idx]->tlen = j-i+5;

	tr[idx]->statetype[0] = STS;
	tr[idx]->nodeidx[0]   = 0;
	tr[idx]->pos[0]       = 0;
	tr[idx]->statetype[1] = STN;
	tr[idx]->nodeidx[1]   = 0;
	tr[idx]->pos[1]       = 0;
	j = 2;			/* now j = position in new subtrace */
	while (1)		/* copy subtrace */
	  {
	    tr[idx]->statetype[j] = otr->statetype[i];
	    tr[idx]->nodeidx[j]   = otr->nodeidx[i];
	    tr[idx]->pos[j]       = otr->pos[i];
	    if (otr->statetype[i] == STE) break;
	    i++; j++;
	  }
	j++;
	tr[idx]->statetype[j] = STC;
	tr[idx]->nodeidx[j]   = 0;
	tr[idx]->pos[j]       = 0;
	j++;
	tr[idx]->statetype[j] = STT;
	tr[idx]->nodeidx[j]   = 0;
	tr[idx]->pos[j]       = 0;
	idx++;
      }

  *ret_tr  = tr;
  *ret_ntr = ntr;
  return;
}


/* Function: TraceDomainNumber()
 * 
 * Purpose:  Count how many times we traverse the
 *           model in a single Plan7 trace -- equivalent
 *           to counting the number of domains.
 *           
 *           (A weakness is that we might discard some of
 *           those domains because they have low scores
 *           below E or T threshold.)
 */
int
TraceDomainNumber(struct p7trace_s *tr)
{
  int i;
  int ndom = 0;
  
  for (i = 0; i < tr->tlen; i++)
    if (tr->statetype[i] == STB) ndom++;
  return ndom;
}


/* Function: TraceSimpleBounds()
 * 
 * Purpose:  For a trace that contains only a single
 *           traverse of the model (i.e. something that's
 *           come from TraceDecompose(), or a global
 *           alignment), determine the bounds of
 *           the match on both the sequence [1..L] and the
 *           model [1..M].
 *           
 * Args:     tr   - trace to look at
 *           i1   - RETURN: start point in sequence [1..L]
 *           i2   - RETURN: end point in sequence [1..L]
 *           k1   - RETURN: start point in model [1..M]
 *           k2   - RETURN: end point in model [1..M]
 */
void
TraceSimpleBounds(struct p7trace_s *tr, int *ret_i1, int *ret_i2, 
		  int *ret_k1,  int *ret_k2)
{
  int i1, i2, k1, k2, tpos;

  i1 = k1 = i2 = k2 = -1;

				/* Look forwards to find start of match */
  for (tpos = 0; tpos < tr->tlen; tpos++)
    {
      if (k1 == -1 && (tr->statetype[tpos] == STM || tr->statetype[tpos] == STD))
	k1 = tr->nodeidx[tpos];
      if (tr->statetype[tpos] == STM)
	{
	  i1 = tr->pos[tpos];
	  break;
	}
    }
  if (tpos == tr->tlen || i1 == -1 || k1 == -1)
    Die("sanity check failed: didn't find a match state in trace");

				/* Look backwards to find end of match */
  for (tpos = tr->tlen-1; tpos >= 0; tpos--)
    {
      if (k2 == -1 && (tr->statetype[tpos] == STM || tr->statetype[tpos] == STD))
	k2 = tr->nodeidx[tpos];
      if (tr->statetype[tpos] == STM)
	{
	  i2 = tr->pos[tpos];
	  break;
	}
    }
  if (tpos == tr->tlen || i2 == -1 || k2 == -1)
    Die("sanity check failed: didn't find a match state in trace");

  *ret_k1 = k1;
  *ret_i1 = i1;
  *ret_k2 = k2;
  *ret_i2 = i2;
}


/* Function: MasterTraceFromMap()
 * Date:     SRE, Tue Jul  7 18:51:11 1998 [St. Louis]
 *
 * Purpose:  Convert an alignment map (e.g. hmm->map) to
 *           a master trace. Used for mapping an alignment
 *           onto an HMM. Generally precedes a call to
 *           ImposeMasterTrace(). Compare P7ViterbiAlignAlignment(),
 *           which aligns an alignment to the model using a
 *           Viterbi algorithm to get a master trace. 
 *           MasterTraceFromMap() only works if the alignment
 *           is exactly the one used to train the model.
 *
 * Args:     map  - the map (usually hmm->map is passed) 1..M
 *           M    - length of map (model; usually hmm->M passed)
 *           alen - length of alignment that map refers to
 *
 * Returns:  ptr to master trace
 *           Caller must free: P7FreeTrace().
 */
struct p7trace_s *
MasterTraceFromMap(int *map, int M, int alen)
{
  struct p7trace_s *tr;         /* RETURN: master trace */ 
  int tpos;			/* position in trace */
  int apos;			/* position in alignment, 1..alen */
  int k;			/* position in model */

  /* Allocate for the trace.
   * S-N-B- ... - E-C-T  : 6 states + alen is maximum trace,
   * because each of alen columns is an N*, M*, I*, or C* metastate.
   * No D* metastates possible.
   */
  P7AllocTrace(alen+6, &tr);

  /* Initialize the trace
   */
  tpos = 0;
  TraceSet(tr, tpos, STS, 0, 0); tpos++;
  TraceSet(tr, tpos, STN, 0, 0); tpos++;

  /* Leading N's
   */
  for (apos = 1; apos < map[1]; apos++) {
    TraceSet(tr, tpos, STN, 0, apos); tpos++;
  } /* now apos == map[1] */
  TraceSet(tr, tpos, STB, 0, 0); tpos++;

  for (k = 1; k < M; k++)
    {
      TraceSet(tr, tpos, STM, k, apos); tpos++;
      apos++;

      for (; apos < map[k+1]; apos++) {
	TraceSet(tr, tpos, STI, k, apos); tpos++;
      }
    } /* now apos == map[M] and k == M*/
      
  TraceSet(tr, tpos, STM, M, apos); tpos++;
  apos++;

  /* Trailing C's
   */
  TraceSet(tr, tpos, STE, 0, 0); tpos++;
  TraceSet(tr, tpos, STC, 0, 0); tpos++;
  for (; apos <= alen; apos++) {
    TraceSet(tr, tpos, STC, 0, apos); tpos++;
  }

  /* Terminate and return
   */
  TraceSet(tr, tpos, STT, 0, 0); tpos++;
  tr->tlen = tpos;
  return tr;
}



/* Function: ImposeMasterTrace()
 * Date:     SRE, Sun Jul  5 14:27:16 1998 [St. Louis]
 *
 * Purpose:  Goes with P7ViterbiAlignAlignment(), which gives us
 *           a "master trace" for a whole alignment. Now, given
 *           the alignment and the master trace, construct individual
 *           tracebacks for each sequence. Later we'll hand these
 *           (and presumably other traces) to P7Traces2Alignment().
 *           
 *           It is possible to generate individual traces that
 *           are not consistent with Plan7 (e.g. D->I and I->D 
 *           transitions may be present). P7Traces2Alignment()
 *           can handle such traces; other functions may not.
 *           See modelmaker.c:trace_doctor() if this is a problem.
 * 
 *           Akin to modelmaker.c:fake_tracebacks().
 *
 * Args:     aseq  - aligned seqs
 *           nseq  - number of aligned seqs 
 *           mtr   - master traceback
 *           ret_tr- RETURN: array of individual tracebacks, one for each aseq
 *
 * Returns:  (void)
 */
void
ImposeMasterTrace(char **aseq, int nseq, struct p7trace_s *mtr, struct p7trace_s ***ret_tr)
{
  struct p7trace_s **tr;
  int  idx;			/* counter over sequences */
  int  i;                       /* position in raw sequence (1..L) */
  int  tpos;			/* position in traceback           */
  int  mpos;			/* position in master trace        */

  tr = (struct p7trace_s **) MallocOrDie (sizeof(struct p7trace_s *) * nseq);
  
  for (idx = 0; idx < nseq; idx++)
    {
      P7AllocTrace(mtr->tlen, &tr[idx]); /* we're guaranteed that individuals len < master len */
      
      tpos = 0;
      i    = 1;
      for (mpos = 0; mpos < mtr->tlen; mpos++)
	{
	  switch (mtr->statetype[mpos]) 
	    {
	    case STS:		/* straight copies w/ no emission: S, B, D, E, T*/
	    case STB:
	    case STD:
	    case STE:
	    case STT:
	      TraceSet(tr[idx], tpos, mtr->statetype[mpos], mtr->nodeidx[mpos], 0);
	      tpos++;
	      break;

	    case STM:		/* M* implies M or D */
	      if (isgap(aseq[idx][mtr->pos[mpos]-1])) 
		TraceSet(tr[idx], tpos, STD, mtr->nodeidx[mpos], 0);
	      else {
		TraceSet(tr[idx], tpos, STM, mtr->nodeidx[mpos], i);
		i++;
	      }
	      tpos++;
	      break;

	    case STI:		/* I* implies I or nothing */
	      if (!isgap(aseq[idx][mtr->pos[mpos]-1])) {
		TraceSet(tr[idx], tpos, STI, mtr->nodeidx[mpos], i);
		i++;
		tpos++;
	      }
	      break;

	    case STJ:		/* N,J,C: first N* -> N. After that, N* -> N or nothing. */
	    case STN:
	    case STC:
	      if (mtr->pos[mpos] == 0) { 
		TraceSet(tr[idx], tpos, mtr->statetype[mpos], 0, 0);
		tpos++;
	      } else if (!isgap(aseq[idx][mtr->pos[mpos]-1])) {
		TraceSet(tr[idx], tpos, mtr->statetype[mpos], 0, i);
		i++;
		tpos++; 
	      }
	      break;

	    case STBOGUS:
	      Die("never happens. Trust me.");
	    }
	}
      tr[idx]->tlen = tpos;
    }	  
  *ret_tr = tr;
}


/* Function: rightjustify()
 * 
 * Purpose:  Given a gap-containing string of length n,
 *           pull all the non-gap characters as far as
 *           possible to the right, leaving gaps on the
 *           left side. Used to rearrange the positions
 *           of insertions in HMMER alignments.
 */
static void
rightjustify(char *s, int n)
{
  int npos;
  int opos;

  npos = n-1;
  opos = n-1;
  while (opos >= 0) {
    if (isgap(s[opos])) opos--;
    else                s[npos--]=s[opos--];  
  }
  while (npos >= 0) 
    s[npos--] = '.';
}


