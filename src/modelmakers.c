/* modelmakers.c
 * Construction of models from multiple alignments. Two versions:
 *    p7_Handmodelmaker() -- use #=RF annotation to indicate match columns
 *    p7_Fastmodelmaker() -- Krogh/Haussler heuristic
 *                          
 * The meat of the model construction code is in matassign2hmm().
 * The two model construction strategies simply label which columns
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

#include "hmmer.h"
#include "plan7.h"
#include "p7_trace.h"

/* flags used for matassign[] arrays -- 
 *   assignment of aligned columns to match/insert states
 */
#define p7_ASSIGN_MATCH      (1<<0) 
#define p7_FIRST_MATCH       (1<<1)
#define p7_LAST_MATCH        (1<<2)
#define p7_ASSIGN_INSERT     (1<<3)
#define p7_EXTERNAL_INSERT_N (1<<4)
#define p7_EXTERNAL_INSERT_C (1<<5) 

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
 *           Construct an HMM from an alignment, where the <#=RF> line
 *           of the alignment file is used to indicate
 *           the columns assigned to matches vs. inserts.
 *           
 *           NOTE: <p7_Handmodelmaker()> will slightly revise the
 *           alignment if necessary, if the assignment of columns
 *           implies DI and ID transitions.
 *           
 *           Returns both the HMM in counts form (ready for applying
 *           Dirichlet priors as the next step), and fake tracebacks
 *           for each aligned sequence.
 *           
 * Args:     msa     - multiple sequence alignment          
 *           abc     - symbol alphabet to use
 *           dsq     - digitized unaligned aseq's
 *           ret_hmm - RETURN: counts-form HMM
 *           ret_tr  - RETURN: array of tracebacks for aseq's
 *           
 * Return:   <eslOK> on success.
 *
 *           <p7ERR_RF> if alignment doesn't have RF column annotation.
 *           ret_hmm and ret_tr alloc'ed here. 
 *           
 * Throws:   <eslEMEM> on allocation failure.          
 */            
int
p7_Handmodelmaker(ESL_MSA *msa, ESL_ALPHABET *abc, char **dsq, 
		  P7_HMM **ret_hmm, P7_TRACE ***ret_tr)
{
  int  status;
  int *matassign = NULL;    /* MAT state assignments if 1; 1..alen */
  int  apos;                /* counter for aligned columns         */
  *ret_hmm = NULL;
  *ret_tr  = NULL;

  if (msa->rf == NULL) { status = p7_ERR_NO_RF; goto CLEANEXIT; }
  ESL_ALLOC(matassign, sizeof(int) * (msa->alen+1));
  esl_vec_ISet(matassign, msa->alen+1, 0);
  
  for (apos = 0; apos < msa->alen; apos++)
    if (!esl_abc_CIsGap(abc, msa->rf[apos]))
      matassign[apos+1] |= p7_ASSIGN_MATCH;
    else
      matassign[apos+1] |= p7_ASSIGN_INSERT;

  status = matassign2hmm(msa, abc, dsq, matassign, ret_hmm, ret_tr);
  /* matassign2hmm leaves ret_hmm, ret_tr in their proper state */
 CLEANEXIT:
  if (matassign != NULL) free(matassign);
  return status;
}


/* Function: p7_Fastmodelmaker()
 * 
 * Purpose:  Heuristic model construction:
 *           Construct an HMM from an alignment using a heuristic,
 *           based on the fractional occupancy of each columns w/
 *           residues vs gaps. Roughly, any column w/ a fractional
 *           occupancy of $\geq$ <thresh> is assigned as a MATCH column;
 *           for instance, if thresh = 0.5, columns w/ $\geq$ 50\% 
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
 *           NOTE: p7_Fastmodelmaker() will slightly revise the
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
 * Return:   <eslOK> on success.
 *           ret_hmm and ret_tr alloc'ed here; FreeTrace(tr[i]), free(tr),
 *           FreeHMM(hmm).       
 *           
 * Throws:   <eslEMEM> on allocation failure.          
 */
int
p7_Fastmodelmaker(ESL_MSA *msa, char **dsq, char *isfrag, float symfrac,
		  P7_HMM **ret_hmm, P7_TRACE ***ret_tr)
{
  int      status;	     /* return status flag                  */
  int     *matassign = NULL; /* MAT state assignments if 1; 1..alen */
  int      idx;              /* counter over sequences              */
  int      apos;             /* counter for aligned columns         */
  float   *r = NULL;	     /* weighted frac of gaps in column     */
  float    totwgt;	     /* total non-fragment seq weight       */
  float    maxR;     	     /* maximum r_i                         */
  int      incfrags;	     /* TRUE to include candidate frags     */
  *ret_hmm = NULL;
  *ret_tr  = NULL;

  /* Allocations: matassign is 1..alen array of bit flags;
   *              gapfrac is 1..alen array of fractions 0<=gapfrac[i]<=1
   */
  ESL_ALLOC(matassign, sizeof(int)     * (msa->alen+1));
  ESL_ALLOC(r,         sizeof(float)   * (msa->alen+1));

  esl_vec_ISet(matassign, msa->alen+1, 0);  
  esl_vec_FSet(r,         msa->alen+1, 0.);

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
      totwgt = esl_vec_FSum(msa->wgt, msa->nseq);
      incfrags = TRUE;
    }

  /* Determine weighted sym freq in each column; keep track of max.
   * Mind the off-by-one (r is 1..alen, msa->aseq is 0..alen-1)
   */
  for (apos = 0; apos < msa->alen; apos++) 
    {  
      for (idx = 0; idx < msa->nseq; idx++) 
	if ((incfrags || ! isfrag[idx]) 
	    && ! isgap(msa->aseq[idx][apos])) 
	  r[apos+1] += msa->wgt[idx];
      r[apos+1] /= totwgt;
    }
  maxR = esl_vec_FMax(r+1, msa->alen);
  
  /* Determine match assignment. (Both matassign and r are 1..alen)
   */
  for (apos = 1; apos <= msa->alen; apos++) 
    if (r[apos] >= symfrac * maxR) matassign[apos] |= p7_ASSIGN_MATCH;
    else	                   matassign[apos] |= p7_ASSIGN_INSERT;

  /* Once we have matassign calculated, modelmakers behave
   * the same; matassign2hmm() does this stuff (traceback construction,
   * trace counting) and sets up ret_hmm and ret_tr.
   */
  status = matassign2hmm(msa, dsq, isfrag, matassign, ret_hmm, ret_tr);
 CLEANEXIT:
  if (r != NULL)         free(r);
  if (matassign != NULL) free(matassign);
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
 *           matassign - 1..alen bit flags for column assignments
 *           ret_hmm   - RETURN: counts-form HMM
 *           ret_tr    - RETURN: array of tracebacks for aseq's
 *                         
 * Return:   <eslOK> on success.
 *           <p7_ERR_NO_CONSENSUS> if no consensus columns are identified.
 *
 *           ret_hmm and ret_tr alloc'ed here.
 */
static int
matassign2hmm(ESL_MSA *msa, char **dsq, int *matassign, 
	      P7_HMM **ret_hmm, P7_TRACE ***ret_tr)
{
  int        status;		/* return status                       */
  P7_HMM    *hmm = NULL;        /* RETURN: new hmm                     */
  P7_TRACE **tr  = NULL;        /* RETURN: 0..nseq-1 fake traces       */
  int      M;                   /* length of new model in match states */
  int      idx;                 /* counter over sequences              */
  int      apos;                /* counter for aligned columns         */
  *ret_hmm = NULL;
  *ret_tr  = NULL;

  /* How many match states in the HMM? */
  for (M=0,apos = 1; apos <= msa->alen; apos++) 
    if (matassign[apos] & p7_ASSIGN_MATCH) M++;
  if (M == 0) { status = p7ERR_NO_CONSENSUS; goto CLEANEXIT; }

  /* Delimit N-terminal tail */
  for (apos=1; matassign[apos] & p7_ASSIGN_INSERT && apos <= msa->alen; apos++)
    matassign[apos] |= p7_EXTERNAL_INSERT_N;
  if (apos <= msa->alen) matassign[apos] |= p7_FIRST_MATCH;

  /* Delimit C-terminal tail */
  for (apos=msa->alen; matassign[apos] & p7_ASSIGN_INSERT && apos > 0; apos--)
    matassign[apos] |= p7_EXTERNAL_INSERT_C;
  if (apos > 0) matassign[apos] |= p7_LAST_MATCH;

  /* Make fake tracebacks for each seq */
  ESL_TRY( fake_tracebacks(msa->aseq, msa->nseq, msa->alen, matassign, &tr) );

  /* Build count model from tracebacks */
  hmm = p7_hmm_Create(M, global_abc);
  if (hmm == NULL) { status = eslEMEM; goto CLEANEXIT; }
  p7_hmm_ZeroCounts(hmm);
  for (idx = 0; idx < msa->nseq; idx++) {
    /* p7_trace_Dump(stdout, tr[idx], NULL, NULL);   */
    p7_trace_Count(hmm, dsq[idx], msa->wgt[idx], tr[idx]);
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
 CLEANEXIT:
  if (ret_tr != NULL) *ret_tr = tr;
  else   { for (idx = 0; idx < msa->nseq; idx++) P7FreeTrace(tr[idx]); free(tr); }
  if (ret_hmm != NULL) *ret_hmm = hmm; else FreePlan7(hmm);
  return status;
}
  


/* Function: fake_tracebacks()
 * 
 * Purpose:  From a consensus assignment of columns to MAT/INS, construct fake
 *           tracebacks for each individual sequence.
 *           
 *           Tracebacks are always relative to the search model (with
 *           internal entry/exit, B->Mk and Mk->E).
 *           
 * Args:     aseqs     - alignment [0..nseq-1][0..alen-1]
 *           nseq      - number of seqs in alignment
 *           alen      - length of alignment in columns
 *           matassign - assignment of column; [1..alen] (off one from aseqs)
 *           ret_tr    - RETURN: array of tracebacks
 *           
 * Return:   (void)
 *           ret_tr is alloc'ed here. Caller must free.
 */          
static int
fake_tracebacks(char **aseq, int nseq, int alen, int *matassign,
		P7_TRACE ***ret_tr)
{
  int  status;
  P7_TRACE **tr = NULL;
  int  idx;                     /* counter over sequences          */
  int  i;                       /* position in raw sequence (1..L) */
  int  k;                       /* position in HMM                 */
  int  apos;                    /* position in alignment columns   */

  ESL_ALLOC(tr, sizeof(P7_TRACE *) * nseq);
  for (idx = 0; idx < nseq; idx++) tr[idx] = NULL;
  
  for (idx = 0; idx < nseq; idx++)
    {
      ESL_TRY( p7_trace_Create(alen+6, &tr[idx]) ); /* +6 = S,N,B,E,C,T */
      
      ESL_TRY( p7_trace_Append(tr, p7_STS, 0, 0) ); /* traces start with S... */
      ESL_TRY( p7_trace_Append(tr, p7_STN, 0, 0) ); /*...and transit to N.    */

      i = 1;
      k = 0;
      for (apos = 0; apos < alen; apos++)
        {
	  if (matassign[apos+1] & p7_FIRST_MATCH) 	/* BEGIN */
	    ESL_TRY( p7_trace_Append(tr, p7_STB, 0, 0) );

	  if (matassign[apos+1] & p7_ASSIGN_MATCH && ! esl_abc_CIsGap(aseq[idx][apos]))
	    {			/* MATCH */
	      k++;		/* move to next model pos */
	      ESL_TRY( p7_trace_Append(tr, p7_STM, k, i) );
	      i++;
	    }	      
          else if (matassign[apos+1] & p7_ASSIGN_MATCH)
            {                   /* DELETE */ /* No B->D transitions */
	      k++;		/* *always* move on model when ASSIGN_MATCH */
	      if (tr[idx]->st[tpos-1] != STB)
		ESL_TRY( p7_trace_Append(tr, p7_STD, k, 0) );
            }
          else if (matassign[apos+1] & p7_EXTERNAL_INSERT_N &&
		   ! esl_abc_CIsGap(aseq[idx][apos]))
            {                   /* N-TERMINAL TAIL */
	      ESL_TRY( p7_trace_Append(tr, p7_STN, 0, i) );
	      i++;
            }
	  else if (matassign[apos+1] & p7_EXTERNAL_INSERT_C &&
		   ! esl_abc_CIsGap(aseq[idx][apos]))
	    {			/* C-TERMINAL TAIL */
	      ESL_TRY( p7_trace_Append(tr, p7_STC, 0, i) );
	      i++;
	    }
	  else if (! esl_abc_CIsGap(aseq[idx][apos]))
	    {			/* INSERT */
	      ESL_TRY( p7_trace_Append(tr, p7_STI, k, i) );
	      i++;
	    }

	  if (matassign[apos+1] & p7_LAST_MATCH)
	    {			/* END */
	      /* be careful about S/W transitions; if we're 
               * a candidate sequence fragment, we need to roll
	       * back over some D's because there's no D->E transition
               * for fragments. k==M right now; don't change it.
	       */
	      while (tr[idx]->st[tr->N-1] == p7_STD) 
		tr->N--;	/* gratuitous chumminess with trace */
	      ESL_TRY( p7_trace_Append(tr, p7_STE, 0, 0) );
	      ESL_TRY( p7_trace_Append(tr, p7_STC, 0, 0) );
	    }
        }
      ESL_TRY( p7_trace_Append(tr, p7_STT, 0, 0) );

      /* deal with DI, ID transitions and other plan7 impossibilities */
      /* (k == M at this point) */
      ESL_TRY ( trace_doctor(tr[idx], k, NULL, NULL) );
    } 

  status = eslOK;
 CLEANEXIT:
  if (status != eslOK) { esl_Free2D((void **) tr, nseq); tr = NULL; }
  *ret_tr = tr;
  return status;
}

/* Function: trace_doctor()
 * 
 * Purpose:  Plan 7 disallows D->I and I->D "chatter" transitions.
 *           However, these transitions may be implied by many
 *           alignments. trace_doctor() collapses I->D or D->I into a
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
 * Return:   <eslOK>
 *           tr is modified.
 */               
static int
trace_doctor(P7_TRACE *tr, int mlen, int *ret_ndi, int *ret_nid)
{
  int opos;			/* position in old trace                 */
  int npos;			/* position in new trace (<= opos)       */
  int ndi, nid;			/* number of DI, ID transitions doctored */

				/* overwrite the trace from left to right */
  ndi  = nid  = 0;
  opos = npos = 0;
  while (opos < tr->N) {
      /* fix implied D->I transitions; D transforms to M, I pulled in */
    if (tr->st[opos] == p7_STD && tr->st[opos+1] == p7_STI) {
      tr->st[npos] = p7_STM;
      tr->k[npos]  = tr->k[opos];     /* D transforms to M      */
      tr->i[npos]  = tr->i[opos+1];   /* insert char moves back */
      opos += 2;
      npos += 1;
      ndi++;
    } /* fix implied I->D transitions; D transforms to M, I is pushed in */
    else if (tr->st[opos]== p7_STI && tr->st[opos+1]== p7_STD) {
      tr->st[npos] = p7_STM;
      tr->k[npos]  = tr->k[opos+1];    /* D transforms to M    */
      tr->i[npos]  = tr->i[opos];      /* insert char moves up */
      opos += 2;
      npos += 1;
      nid++; 
    } /* fix implied B->I transitions; pull I back to its M */
    else if (tr->st[opos]== p7_STI && tr->st[opos-1]== p7_STB) {
      tr->st[npos] = p7_STM;
      tr->k[npos]  = tr->k[opos];    /* offending I transforms to M */
      tr->i[npos]  = tr->i[opos];
      opos++;
      npos++;
    } /* fix implied I->E transitions; push I to next M */
    else if (tr->st[opos]== p7_STI && tr->st[opos+1]== p7_STE) {
      tr->st[npos] = p7_STM;
      tr->k[npos]  = tr->k[opos]+1;  /* offending I transforms to M */
      tr->i[npos]  = tr->i[opos];
      opos++;
      npos++;
    } /* rare: N-N-B-E becomes N-B-M_1-E (swap B,N) */
    else if (tr->st[opos]==p7_STB && tr->st[opos+1]==p7_STE
	     && tr->st[opos-1]==p7_STN && tr->i[opos-1] > 0) {   
      tr->st[npos]   = p7_STM;
      tr->k[npos]    = 1;
      tr->i[npos]    = tr->i[opos-1];
      tr->st[npos-1] = p7_STB;
      tr->k[npos-1]  = 0;
      tr->i[npos-1]  = 0;
      opos++;
      npos++;
    } /* rare: B-E-C-C-x becomes B-M_M-E-C-x (swap E,C) */
    else if (tr->st[opos]==p7_STE && tr->st[opos-1]==p7_STB
	     && tr->st[opos+1]==p7_STC 
	     && tr->st[opos+2]==p7_STC) { 
      tr->st[npos]   = p7_STM;
      tr->k[npos]    = mlen;
      tr->i[npos]    = tr->pos[opos+2];
      tr->st[npos+1] = p7_STE;
      tr->k[npos+1]  = 0;
      tr->i[npos+1]  = 0;
      tr->st[npos+2] = p7_STC; /* first C must be a nonemitter  */
      tr->k[npos+2]  = 0;
      tr->i[npos+2]  = 0;
      opos+=3;
      npos+=3;
    } /* everything else is just copied */
    else {
      tr->st[npos] = tr->st[opos];
      tr->k[npos]  = tr->k[opos];
      tr->i[npos]  = tr->i[opos];
      opos++;
      npos++;
    }
  }
  tr->N = npos;

  if (ret_ndi != NULL) *ret_ndi = ndi;
  if (ret_nid != NULL) *ret_nid = nid;
  return eslOK;
}


/* Function: annotate_model()
 * 
 * Purpose:  Transfer rf, cs, and other optional annotation from the alignment
 *           to the new model.
 * 
 * Args:     hmm       - new model
 *           matassign - which alignment columns are MAT; [1..alen]
 *           msa       - alignment, including annotation to transfer
 *           
 * Return:   <eslOK> on success.
 */
static int
annotate_model(P7_HMM *hmm, int *matassign, ESL_MSA *msa)
{                      
  int   apos;			/* position in matassign, 1.alen  */
  int   k;			/* position in model, 1.M         */
  char *pri;			/* X-PRM, X-PRI, X-PRT annotation */

  /* Reference coord annotation  */
  if (msa->rf != NULL) {
    hmm->rf[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH) /* apos is off by one from HMM */
	hmm->rf[k++] = msa->rf[apos-1];
    hmm->rf[k] = '\0';
    hmm->flags |= p7_RF;
  }

  /* Consensus structure annotation */
  if (msa->ss_cons != NULL) {
    hmm->cs[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH)
	hmm->cs[k++] = msa->ss_cons[apos-1];
    hmm->cs[k] = '\0';
    hmm->flags |= p7_CS;
  }

  /* Surface accessibility annotation */
  if (msa->sa_cons != NULL) {
    hmm->ca[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos] & ASSIGN_MATCH)
	hmm->ca[k++] = msa->sa_cons[apos-1];
    hmm->ca[k] = '\0';
    hmm->flags |= p7_CA;
  }

  /* The alignment map (1..M in model, 1..alen in alignment) */
  for (apos = k = 1; apos <= msa->alen; apos++)
    if (matassign[apos] & ASSIGN_MATCH)
      hmm->map[k++] = apos;
  hmm->flags |= p7_MAP;

  return eslOK;
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

