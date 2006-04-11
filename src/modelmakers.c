/* modelmakers.c
 * Construction of models from multiple alignments. Three versions:
 *    p7_Handmodelmaker() -- use #=RF annotation to indicate match columns
 *    p7_Fastmodelmaker() -- Krogh/Haussler heuristic
 *    p7_Maxmodelmaker()  -- MAP model construction algorithm (Eddy, 
 *                           unpublished)
 *                          
 * The meat of the model construction code is in matassign2hmm().
 * The three model construction strategies simply label which columns
 * are supposed to be match states, and then hand this info to
 * matassign2hmm().
 * 
 * Two wrinkles to watch for:
 * 1) The alignment is assumed to contain sequence fragments. Look in
 *    fake_tracebacks() for how internal entry/exit points are handled.
 * 2) Plan7 disallows DI and ID transitions, but an alignment may
 *    imply these. Look in trace_doctor() for how DI and ID transitions
 *    are removed.
 *
 * SVN $Id: modelmakers.c 1384 2005-05-06 23:04:04Z eddy $
 * SRE, Fri Nov 15 10:00:04 1996
 */

#include "p7config.h"

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#include <easel.h>
#include <esl_msa.h>

#include "plan7.h"
#include "p7_trace.h"

/* flags used for matassign[] arrays -- 
 *   assignment of aligned columns to match/insert states
 */
#define ASSIGN_MATCH      (1<<0) 
#define FIRST_MATCH       (1<<1)
#define LAST_MATCH        (1<<2)
#define ASSIGN_INSERT     (1<<3)
#define EXTERNAL_INSERT_N (1<<4)
#define EXTERNAL_INSERT_C (1<<5) 

static void matassign2hmm(ESL_MSA *msa, char **dsq, char *isfrag,
			  int *matassign, P7_HMM **ret_hmm,
			  P7_TRACE ***ret_tr);
static void fake_tracebacks(char **aseq, int nseq, int alen, char *isfrag,
			    int *matassign, P7_TRACE ***ret_tr);
static void trace_doctor(P7_TRACE *tr, int M, int *ret_ndi, int *ret_nid);
static void annotate_model(P7_HMM *hmm, int *matassign, ESL_MSA *msa);
#if 0
static void print_matassign(int *matassign, int alen);
#endif



/* Function: p7_Handmodelmaker()
 * 
 * Purpose:  Manual model construction:
 *           Construct an HMM from an alignment, where the #=RF line
 *           of a HMMER alignment file is given to indicate
 *           the columns assigned to matches vs. inserts.
 *           
 *           NOTE: Handmodelmaker() will slightly revise the alignment
 *           if necessary, if the assignment of columns implies
 *           DI and ID transitions.
 *           
 *           Returns both the HMM in counts form (ready for applying
 *           Dirichlet priors as the next step), and fake tracebacks
 *           for each aligned sequence.
 *           
 * Args:     msa  - multiple sequence alignment          
 *           dsq  - digitized unaligned aseq's
 *           isfrag  - [0..nseq-1] flags for candidate seq frags
 *           ret_hmm - RETURN: counts-form HMM
 *           ret_tr  - RETURN: array of tracebacks for aseq's
 *           
 * Return:   (void)
 *           ret_hmm and ret_tr alloc'ed here. 
 */            
void
p7_Handmodelmaker(ESL_MSA *msa, char **dsq, char *isfrag,
		  P7_HMM **ret_hmm, P7_TRACE ***ret_tr)
{
  int     *matassign;           /* MAT state assignments if 1; 1..alen */
  int      apos;                /* counter for aligned columns         */

  /* Make sure we have all the info about the alignment that we need */
  if (msa->rf == NULL)
    Die("Alignment must have RF annotation to hand-build an HMM");

				/* Allocation */
  matassign = (int *) MallocOrDie (sizeof(int) * (msa->alen+1));
  
  /* Determine match assignment from optional annotation
   */
  matassign[0] = 0;
  for (apos = 0; apos < msa->alen; apos++)
    {
      matassign[apos+1] = 0;
      if (!isgap(msa->rf[apos])) 
	matassign[apos+1] |= ASSIGN_MATCH;
      else 
	matassign[apos+1] |= ASSIGN_INSERT;
    }

  /* Hand matassign off for remainder of model construction
   */
  /*   print_matassign(matassign, msa->alen); */
  matassign2hmm(msa, dsq, isfrag, matassign, ret_hmm, ret_tr);

  free(matassign);
  return;
}


/* Function: P7Fastmodelmaker()
 * 
 * Purpose:  Heuristic model construction:
 *           Construct an HMM from an alignment using a heuristic,
 *           based on the fractional occupancy of each columns w/
 *           residues vs gaps. Roughly, any column w/ a fractional
 *           occupancy of >= <thresh> is assigned as a MATCH column;
 *           for instance, if thresh = 0.5, columns w/ >= 50% 
 *           residues are assigned to match... roughly speaking.
 *           
 *           "Roughly speaking" because sequences are weighted; we
 *           guard against fragments counting against us by not
 *           counting frags at all unless we have to; and we guard
 *           against highly gappy alignments by reweighting the threshold
 *           according to the most occupied column in the alignment.
 *           
 *           The detailed calculation is:
 *           
 *           Count the weighted fraction r_i of symbols in column i
 *           (weighted by relative sequence weights): 0<=r_i<=1, with
 *           r_i = 1.0 in fully occupied columns. Only sequences that
 *           are not flagged as candidate fragments are counted toward
 *           the r_i calculation. (If *all* sequences are candidate
 *           fragments, then revert to counting them, rather than
 *           giving up w/ undefined r_i's.)
 *           
 *           Determine the fraction of symbols in the most occupied
 *           column: R = max_i r_i. Normally R=1.0, but if we're given
 *           a gappy alignment, R may be less than that.
 *           
 *           Then given threshold t:
 *           
 *           if r_i >= tR, column is assigned as a MATCH column;
 *           else it's an INSERT column.
 *
 *           NOTE: Fastmodelmaker() will slightly revise the
 *           alignment if the assignment of columns implies
 *           DI and ID transitions.
 *           
 *           Returns the HMM in counts form (ready for applying Dirichlet
 *           priors as the next step). Also returns fake traceback
 *           for each training sequence.
 *           
 * Args:     msa       - multiple sequence alignment
 *           dsq       - digitized unaligned aseq's
 *           isfrag    - [0..nseq-1] T/F flags for candidate seq frags
 *           symfrac   - threshold for residue occupancy
 *           ret_hmm   - RETURN: counts-form HMM
 *           ret_tr    - RETURN: array of tracebacks for aseq's
 *           
 * Return:   (void)
 *           ret_hmm and ret_tr alloc'ed here; FreeTrace(tr[i]), free(tr),
 *           FreeHMM(hmm).       
 */
void
P7Fastmodelmaker(MSA *msa, unsigned char **dsq, char *isfrag, float symfrac,
		 struct plan7_s **ret_hmm, struct p7trace_s ***ret_tr)
{
  int     *matassign;        /* MAT state assignments if 1; 1..alen */
  int      idx;              /* counter over sequences              */
  int      apos;             /* counter for aligned columns         */
  float   *r;		     /* weighted frac of gaps in column     */
  float    totwgt;	     /* total non-fragment seq weight       */
  float    maxR;     	     /* maximum r_i                         */
  int      incfrags;		/* TRUE to include candidate frags  */

  /* Allocations: matassign is 1..alen array of bit flags;
   *              gapfrac is 1..alen array of fractions 0<=gapfrac[i]<=1
   */
  matassign = MallocOrDie (sizeof(int)   * (msa->alen+1));
  r         = MallocOrDie (sizeof(float) * (msa->alen+1));
  
  /* Determine total non-frag weight, just once.
   */
  incfrags = FALSE;
  totwgt   = 0.;
  for (idx = 0; idx < msa->nseq; idx++) 
    if (! isfrag[idx]) totwgt += msa->wgt[idx];

  /* Fallback position, if we don't have any non-candidate frags:
   * count all seqs.
   */
  if (totwgt == 0.) /* yes, this fp compare is safe */
    {
      totwgt = FSum(msa->wgt, msa->nseq);
      incfrags = TRUE;
    }

  /* Determine weighted sym freq in each column; keep track of max.
   * Mind the off-by-one (r is 1..alen, msa->aseq is 0..alen-1)
   */
  for (apos = 0; apos < msa->alen; apos++) 
    {  
      r[apos+1] = 0.;
      for (idx = 0; idx < msa->nseq; idx++) 
	if ((incfrags || ! isfrag[idx]) 
	    && ! isgap(msa->aseq[idx][apos])) 
	  r[apos+1] += msa->wgt[idx];
      r[apos+1] /= totwgt;
    }
  maxR = FMax(r+1, msa->alen);
  
  /* Determine match assignment. (Both matassign and r are 1..alen)
   */
  matassign[0] = 0;
  for (apos = 1; apos <= msa->alen; apos++) 
    {
      matassign[apos] = 0;
      if (r[apos] >= symfrac * maxR)
	matassign[apos] |= ASSIGN_MATCH;
      else
	matassign[apos] |= ASSIGN_INSERT;
    }

  /* Once we have matassign calculated, modelmakers behave
   * the same; matassign2hmm() does this stuff (traceback construction,
   * trace counting) and sets up ret_hmm and ret_tr.
   */
  matassign2hmm(msa, dsq, isfrag, matassign, ret_hmm, ret_tr);

  free(r);
  free(matassign);
  return;
}
  


/* Function: matassign2hmm()
 * 
 * Purpose:  Given an assignment of alignment columns to match vs.
 *           insert, finish the final part of the model construction 
 *           calculation that is constant between model construction
 *           algorithms.
 *           
 * Args:     msa       - multiple sequence alignment
 *           dsq       - digitized unaligned aseq's
 *           isfrag    - 0..nseq-1 T/F flags on candidate fragments
 *           matassign - 1..alen bit flags for column assignments
 *           ret_hmm   - RETURN: counts-form HMM
 *           ret_tr    - RETURN: array of tracebacks for aseq's
 *                         
 * Return:   (void)
 *           ret_hmm and ret_tr alloc'ed here for the calling
 *           modelmaker function.
 */
static void
matassign2hmm(MSA *msa, unsigned char **dsq, char *isfrag, int *matassign, 
	      struct plan7_s **ret_hmm, struct p7trace_s ***ret_tr)
{
  struct plan7_s    *hmm;       /* RETURN: new hmm                     */
  struct p7trace_s **tr;        /* fake tracebacks for each seq        */
  int      M;                   /* length of new model in match states */
  int      idx;                 /* counter over sequences              */
  int      apos;                /* counter for aligned columns         */

				/* how many match states in the HMM? */
  M = 0;
  for (apos = 1; apos <= msa->alen; apos++) {
    if (matassign[apos] & ASSIGN_MATCH) 
      M++;
  }

  if (M == 0) 
    Die("No conserved consensus columns found; aborting construction!\n\
This is an unusual situation. If you're using default --fast heuristic\n\
construction, reexamine your alignment; it is probably unusually full of\n\
gaps, or lots of sequence fragments. You may be able to force HMMER to\n\
model it by using the --symfrac option. If you're using --hand construction,\n\
then you didn't mark any match columns in your reference line annotation;\n\
see the documentation.");

				/* delimit N-terminal tail */
  for (apos=1; matassign[apos] & ASSIGN_INSERT && apos <= msa->alen; apos++)
    matassign[apos] |= EXTERNAL_INSERT_N;
  if (apos <= msa->alen) matassign[apos] |= FIRST_MATCH;

				/* delimit C-terminal tail */
  for (apos=msa->alen; matassign[apos] & ASSIGN_INSERT && apos > 0; apos--)
    matassign[apos] |= EXTERNAL_INSERT_C;
  if (apos > 0) matassign[apos] |= LAST_MATCH;

  /* print_matassign(matassign, msa->alen);  */

				/* make fake tracebacks for each seq */
  fake_tracebacks(msa->aseq, msa->nseq, msa->alen, isfrag, matassign, &tr);

				/* build count model from tracebacks */
  hmm = AllocPlan7(M);
  ZeroPlan7(hmm);
  for (idx = 0; idx < msa->nseq; idx++) {
    /* P7PrintTrace(stdout, tr[idx], NULL, NULL);   */
    P7TraceCount(hmm, dsq[idx], msa->wgt[idx], tr[idx]);
  }
				/* annotate new model */
  annotate_model(hmm, matassign, msa);

  /* Set #=RF line of alignment to reflect our assignment
   * of match, delete. matassign is valid from 1..alen and is off
   * by one from msa->rf.
   */
  if (msa->rf != NULL) free(msa->rf);
  msa->rf = (char *) MallocOrDie (sizeof(char) * (msa->alen + 1));
  for (apos = 0; apos < msa->alen; apos++)
    msa->rf[apos] = matassign[apos+1] & ASSIGN_MATCH ? 'x' : '.';
  msa->rf[msa->alen] = '\0';

				/* Cleanup and return. */
  if (ret_tr != NULL) *ret_tr = tr;
  else   { for (idx = 0; idx < msa->nseq; idx++) P7FreeTrace(tr[idx]); free(tr); }
  if (ret_hmm != NULL) *ret_hmm = hmm; else FreePlan7(hmm);
  return;
}
  


/* Function: fake_tracebacks()
 * 
 * Purpose:  From a consensus assignment of columns to MAT/INS, construct fake
 *           tracebacks for each individual sequence.
 *           
 * Note:     Fragment tolerant by default. Internal entries are 
 *           B->M_x, instead of B->D1->D2->...->M_x; analogously
 *           for internal exits. 
 *           
 * Args:     aseqs     - alignment [0..nseq-1][0..alen-1]
 *           nseq      - number of seqs in alignment
 *           alen      - length of alignment in columns
 *           isfrag    - T/F flags for candidate fragments
 *           matassign - assignment of column; [1..alen] (off one from aseqs)
 *           ret_tr    - RETURN: array of tracebacks
 *           
 * Return:   (void)
 *           ret_tr is alloc'ed here. Caller must free.
 */          
static void
fake_tracebacks(char **aseq, int nseq, int alen, char *isfrag, int *matassign,
		struct p7trace_s ***ret_tr)
{
  struct p7trace_s **tr;
  int  idx;                     /* counter over sequences          */
  int  i;                       /* position in raw sequence (1..L) */
  int  k;                       /* position in HMM                 */
  int  apos;                    /* position in alignment columns   */
  int  tpos;			/* position in traceback           */

  tr = (struct p7trace_s **) MallocOrDie (sizeof(struct p7trace_s *) * nseq);
  
  for (idx = 0; idx < nseq; idx++)
    {
      P7AllocTrace(alen+6, &tr[idx]); /* allow room for S,N,B,E,C,T */
      
				/* all traces start with S state... */
      tr[idx]->statetype[0] = STS;
      tr[idx]->nodeidx[0]   = 0;
      tr[idx]->pos[0]       = 0;
				/* ...and transit to N state; N-term tail
				   is emitted on N->N transitions */
      tr[idx]->statetype[1] = STN;
      tr[idx]->nodeidx[1]   = 0;
      tr[idx]->pos[1]       = 0;
      
      i = 1;
      k = 0;
      tpos = 2;
      for (apos = 0; apos < alen; apos++)
        {
	  tr[idx]->statetype[tpos] = STBOGUS; /* bogus, deliberately, to debug */

	  if (matassign[apos+1] & FIRST_MATCH)
	    {			/* BEGIN */
	      tr[idx]->statetype[tpos] = STB;
	      tr[idx]->nodeidx[tpos]   = 0;
	      tr[idx]->pos[tpos]       = 0;
	      tpos++;
	    }

	  if (matassign[apos+1] & ASSIGN_MATCH && ! isgap(aseq[idx][apos]))
	    {			/* MATCH */
	      k++;		/* move to next model pos */
	      tr[idx]->statetype[tpos] = STM;
	      tr[idx]->nodeidx[tpos]   = k;
	      tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
	    }	      
          else if (matassign[apos+1] & ASSIGN_MATCH)
            {                   /* DELETE */
		/* being careful about S/W transitions.
                 * Count B->D1 transitions only if we're not
                 * a candidate fragment seq.
		 */
	      k++;		/* *always* move on model when ASSIGN_MATCH */
	      if (tr[idx]->statetype[tpos-1] != STB || ! isfrag[idx])
		{
		  tr[idx]->statetype[tpos] = STD;
		  tr[idx]->nodeidx[tpos]   = k;
		  tr[idx]->pos[tpos]       = 0;
		  tpos++;
		}
            }
          else if (matassign[apos+1] & EXTERNAL_INSERT_N &&
		   ! isgap(aseq[idx][apos]))
            {                   /* N-TERMINAL TAIL */
              tr[idx]->statetype[tpos] = STN;
              tr[idx]->nodeidx[tpos]   = 0;
              tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
            }
	  else if (matassign[apos+1] & EXTERNAL_INSERT_C &&
		   ! isgap(aseq[idx][apos]))
	    {			/* C-TERMINAL TAIL */
	      tr[idx]->statetype[tpos] = STC;
              tr[idx]->nodeidx[tpos]   = 0;
              tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
	    }
	  else if (! isgap(aseq[idx][apos]))
	    {			/* INSERT */
	      tr[idx]->statetype[tpos] = STI;
              tr[idx]->nodeidx[tpos]   = k;
              tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
	    }

	  if (matassign[apos+1] & LAST_MATCH)
	    {			/* END */
	      /* be careful about S/W transitions; if we're 
               * a candidate sequence fragment, we need to roll
	       * back over some D's because there's no D->E transition
               * for fragments. k==M right now; don't change it.
	       */
	      if (isfrag[idx])
		while (tr[idx]->statetype[tpos-1] == STD) 
		  tpos--;
	      tr[idx]->statetype[tpos] = STE;
	      tr[idx]->nodeidx[tpos]   = 0;
	      tr[idx]->pos[tpos]       = 0;
	      tpos++;
				/* and then transit E->C;
				   alignments that use J are undefined;
				   C-term tail is emitted on C->C transitions */
	      tr[idx]->statetype[tpos] = STC;
	      tr[idx]->nodeidx[tpos]   = 0;
	      tr[idx]->pos[tpos]       = 0;
	      tpos++;
	    }
        }
                                /* all traces end with T state */
      tr[idx]->statetype[tpos] = STT;
      tr[idx]->nodeidx[tpos]   = 0;
      tr[idx]->pos[tpos]       = 0;
      tr[idx]->tlen = ++tpos;
				/* deal with DI, ID transitions */
				/* k == M here */
      trace_doctor(tr[idx], k, NULL, NULL);

    }    /* end for sequence # idx */

  *ret_tr = tr;
  return;
}

/* Function: trace_doctor()
 * 
 * Purpose:  Plan 7 disallows D->I and I->D "chatter" transitions.
 *           However, these transitions may be implied by many
 *           alignments for hand- or heuristic- built HMMs.
 *           trace_doctor() collapses I->D or D->I into a
 *           single M position in the trace. 
 *           Similarly, B->I and I->E transitions may be implied
 *           by an alignment.
 *           
 *           trace_doctor does not examine any scores when it does
 *           this. In ambiguous situations (D->I->D) the symbol
 *           will be pulled arbitrarily to the left, regardless
 *           of whether that's the best column to put it in or not.
 *           
 * Args:     tr      - trace to doctor
 *           M       - length of model that traces are for 
 *           ret_ndi - number of DI transitions doctored
 *           ret_nid - number of ID transitions doctored
 * 
 * Return:   (void)
 *           tr is modified
 */               
static void
trace_doctor(struct p7trace_s *tr, int mlen, int *ret_ndi, int *ret_nid)
{
  int opos;			/* position in old trace                 */
  int npos;			/* position in new trace (<= opos)       */
  int ndi, nid;			/* number of DI, ID transitions doctored */

				/* overwrite the trace from left to right */
  ndi  = nid  = 0;
  opos = npos = 0;
  while (opos < tr->tlen) {
      /* fix implied D->I transitions; D transforms to M, I pulled in */
    if (tr->statetype[opos] == STD && tr->statetype[opos+1] == STI) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos]; /* D transforms to M      */
      tr->pos[npos]       = tr->pos[opos+1];   /* insert char moves back */
      opos += 2;
      npos += 1;
      ndi++;
    } /* fix implied I->D transitions; D transforms to M, I is pushed in */
    else if (tr->statetype[opos]== STI && tr->statetype[opos+1]== STD) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos+1];/* D transforms to M    */
      tr->pos[npos]       = tr->pos[opos];      /* insert char moves up */
      opos += 2;
      npos += 1;
      nid++; 
    } /* fix implied B->I transitions; pull I back to its M */
    else if (tr->statetype[opos]== STI && tr->statetype[opos-1]== STB) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos]; /* offending I transforms to M */
      tr->pos[npos]       = tr->pos[opos];
      opos++;
      npos++;
    } /* fix implied I->E transitions; push I to next M */
    else if (tr->statetype[opos]== STI && tr->statetype[opos+1]== STE) {
      tr->statetype[npos] = STM;
      tr->nodeidx[npos]   = tr->nodeidx[opos]+1;/* offending I transforms to M */
      tr->pos[npos]       = tr->pos[opos];
      opos++;
      npos++;
    } /* rare: N-N-B-E becomes N-B-M_1-E (swap B,N) */
    else if (tr->statetype[opos]==STB && tr->statetype[opos+1]==STE
	     && tr->statetype[opos-1]==STN && tr->pos[opos-1] > 0) {   
      tr->statetype[npos]   = STM;
      tr->nodeidx[npos]     = 1;
      tr->pos[npos]         = tr->pos[opos-1];
      tr->statetype[npos-1] = STB;
      tr->nodeidx[npos-1]   = 0;
      tr->pos[npos-1]       = 0;
      opos++;
      npos++;
    } /* rare: B-E-C-C-x becomes B-M_M-E-C-x (swap E,C) */
    else if (tr->statetype[opos]==STE && tr->statetype[opos-1]==STB
	     && tr->statetype[opos+1]==STC 
	     && tr->statetype[opos+2]==STC) { 
      tr->statetype[npos]   = STM;
      tr->nodeidx[npos]     = mlen;
      tr->pos[npos]         = tr->pos[opos+2];
      tr->statetype[npos+1] = STE;
      tr->nodeidx[npos+1]   = 0;
      tr->pos[npos+1]       = 0;
      tr->statetype[npos+2] = STC; /* first C must be a nonemitter  */
      tr->nodeidx[npos+2]   = 0;
      tr->pos[npos+2]       = 0;
      opos+=3;
      npos+=3;
    } /* everything else is just copied */
    else {
      tr->statetype[npos] = tr->statetype[opos];
      tr->nodeidx[npos]   = tr->nodeidx[opos];
      tr->pos[npos]       = tr->pos[opos];
      opos++;
      npos++;
    }
  }
  tr->tlen = npos;

  if (ret_ndi != NULL) *ret_ndi = ndi;
  if (ret_nid != NULL) *ret_nid = nid;
  return;
}


/* Function: annotate_model()
 * 
 * Purpose:  Add rf, cs optional annotation to a new model.
 * 
 * Args:     hmm       - new model
 *           matassign - which alignment columns are MAT; [1..alen]
 *           msa       - alignment, including annotation to transfer
 *           
 * Return:   (void)
 */
static void
annotate_model(struct plan7_s *hmm, int *matassign, MSA *msa)
{                      
  int   apos;			/* position in matassign, 1.alen  */
  int   k;			/* position in model, 1.M         */
  char *pri;			/* X-PRM, X-PRI, X-PRT annotation */

  /* Transfer reference coord annotation from alignment,
   * if available
   */
  if (msa->rf != NULL) {
    hmm->rf[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH) /* ainfo is off by one from HMM */
	hmm->rf[k++] = (msa->rf[apos-1] == ' ') ? '.' : msa->rf[apos-1];
    hmm->rf[k] = '\0';
    hmm->flags |= PLAN7_RF;
  }

  /* Transfer consensus structure annotation from alignment, 
   * if available
   */
  if (msa->ss_cons != NULL) {
    hmm->cs[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH)
	hmm->cs[k++] = (msa->ss_cons[apos-1] == ' ') ? '.' : msa->ss_cons[apos-1];
    hmm->cs[k] = '\0';
    hmm->flags |= PLAN7_CS;
  }

  /* Transfer surface accessibility annotation from alignment,
   * if available
   */
  if (msa->sa_cons != NULL) {
    hmm->ca[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH)
	hmm->ca[k++] = (msa->sa_cons[apos-1] == ' ') ? '.' : msa->sa_cons[apos-1];
    hmm->ca[k] = '\0';
    hmm->flags |= PLAN7_CA;
  }

  /* Store the alignment map
   */
  for (apos = k = 1; apos <= msa->alen; apos++)
    if (matassign[apos] & ASSIGN_MATCH)
      hmm->map[k++] = apos;
  hmm->flags |= PLAN7_MAP;

  /* Translate and transfer X-PRM annotation. 
   * 0-9,[a-zA-Z] are legal; translate as 0-9,10-35 into hmm->mpri.
   * Any other char is translated as -1, and this will be interpreted
   * as a flag that means "unknown", e.g. use the normal mixture Dirichlet
   * procedure for this column.
   */
  if ((pri = MSAGetGC(msa, "X-PRM")) != NULL)
    {
      hmm->mpri = MallocOrDie(sizeof(int) * (hmm->M+1));
      for (apos = k = 1; apos <= msa->alen; apos++)
	if (matassign[apos] & ASSIGN_MATCH)
	  {
	    if      (isdigit((int) pri[apos-1])) hmm->mpri[k] = pri[apos-1] - '0';
	    else if (islower((int) pri[apos-1])) hmm->mpri[k] = pri[apos-1] - 'a' + 10;
	    else if (isupper((int) pri[apos-1])) hmm->mpri[k] = pri[apos-1] - 'A' + 10;
	    else hmm->mpri[k] = -1;
	    k++;
	  }
    }
  /* And again for X-PRI annotation on insert priors:
   */
  if ((pri = MSAGetGC(msa, "X-PRI")) != NULL)
    {
      hmm->ipri = MallocOrDie(sizeof(int) * (hmm->M+1));
      for (apos = k = 1; apos <= msa->alen; apos++)
	if (matassign[apos] & ASSIGN_MATCH)
	  {
	    if      (isdigit((int) pri[apos-1])) hmm->ipri[k] = pri[apos-1] - '0';
	    else if (islower((int) pri[apos-1])) hmm->ipri[k] = pri[apos-1] - 'a' + 10;
	    else if (isupper((int) pri[apos-1])) hmm->ipri[k] = pri[apos-1] - 'A' + 10;
	    else hmm->ipri[k] = -1;
	    k++;
	  }
    }
  /* And one last time for X-PRT annotation on transition priors:
   */
  if ((pri = MSAGetGC(msa, "X-PRT")) != NULL)
    {
      hmm->tpri = MallocOrDie(sizeof(int) * (hmm->M+1));
      for (apos = k = 1; apos <= msa->alen; apos++)
	if (matassign[apos] & ASSIGN_MATCH)
	  {
	    if      (isdigit((int) pri[apos-1])) hmm->tpri[k] = pri[apos-1] - '0';
	    else if (islower((int) pri[apos-1])) hmm->tpri[k] = pri[apos-1] - 'a' + 10;
	    else if (isupper((int) pri[apos-1])) hmm->tpri[k] = pri[apos-1] - 'A' + 10;
	    else hmm->tpri[k] = -1;
	    k++;
	  }
    }

}

#if 0
static void
print_matassign(int *matassign, int alen)
{
  int apos;

  for (apos = 0; apos <= alen; apos++) {
    printf("%3d %c %c %c\n", 
	   apos,
	   (matassign[apos] & ASSIGN_MATCH) ? 'x':' ',
	   (matassign[apos] & FIRST_MATCH || matassign[apos] & LAST_MATCH) ? '<' : ' ',
	   (matassign[apos] & EXTERNAL_INSERT_N ||
	    matassign[apos] & EXTERNAL_INSERT_C) ? '|':' ');
  }
}
#endif

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

