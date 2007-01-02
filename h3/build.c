/* Construction of new HMMs from multiple alignments.
 * 
 * Two versions: 
 *    p7_Handmodelmaker() -- use #=RF annotation to indicate match columns
 *    p7_Fastmodelmaker() -- Krogh/Haussler heuristic
 * 
 * The maximum likelihood model construction algorithm that was in previous
 * HMMER versions has been deprecated, at least for the moment.
 * 
 * The meat of the model construction code is in matassign2hmm().
 * The two model construction strategies simply label which columns
 * are supposed to be match states, and then hand this info to
 * matassign2hmm().
 * 
 * 
 * Contents:
 *    1. Exported API: model construction routines.
 *    2. Private functions used in constructing models.
 *    3.
 * 
 * SVN $Id$
 */

#include "p7_config.h"

#include "p7_hmm.h"
#include "p7_trace.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

/* flags used for matassign[] arrays -- 
 *   assignment of aligned columns to match/insert states
 */
#define p7_ASSIGN_MATCH      (1<<0) 
#define p7_FIRST_MATCH       (1<<1)
#define p7_LAST_MATCH        (1<<2)
#define p7_ASSIGN_INSERT     (1<<3)
#define p7_EXTERNAL_INSERT_N (1<<4)
#define p7_EXTERNAL_INSERT_C (1<<5) 

static int matassign2hmm(ESL_MSA *msa, ESL_ALPHABET *abc, int *matassign, 
			 P7_HMM **ret_hmm, P7_TRACE ***ret_tr);

/*****************************************************************
 * 1. Exported API: model construction routines.
 *****************************************************************/

/* Function: p7_Handmodelmaker()
 * 
 * Purpose:  Manual model construction.
 *           Construct an HMM from a digital alignment, where the
 *           <#=RF> line of the alignment file is used to indicate the
 *           columns assigned to matches vs. inserts.
 *           
 *           The <msa> must be in digital mode, and it must have
 *           a reference annotation line.
 *           
 *           NOTE: <p7_Handmodelmaker()> will slightly revise the
 *           alignment if necessary, if the assignment of columns
 *           implies DI and ID transitions.
 *           
 *           Returns both the HMM in counts form (ready for applying
 *           Dirichlet priors as the next step), and fake tracebacks
 *           for each aligned sequence. 
 *           
 *           Models must have at least one node, so if the <msa> defined 
 *           no consensus columns, a <eslENORESULT> error is returned.
 *           
 * Args:     msa     - multiple sequence alignment          
 *           ret_hmm - RETURN: counts-form HMM
 *           ret_tr  - RETURN: array of tracebacks for aseq's
 *           
 * Return:   <eslOK> on success. <ret_hmm> and <ret_tr> are allocated 
 *           here, and must be free'd by caller.
 *
 *           Returns <eslENORESULT> if no consensus columns were annotated;
 *           in this case, <ret_hmm> and <ret_tr> are returned NULL.
 *           
 * Throws:   <eslEMEM> on allocation failure. <eslECONTRACT> if the <msa>
 *           does not have a reference annotation line, or if it isn't
 *           in digital mode.
 */            
int
p7_Handmodelmaker(ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***ret_tr)
{
  int  status;
  int *matassign = NULL;    /* MAT state assignments if 1; 1..alen */
  int  apos;                /* counter for aligned columns         */

  *ret_hmm = NULL;
  *ret_tr  = NULL;

  if (! (msa->flags & eslMSA_DIGITAL)) ESL_XEXCEPTION(eslECONTRACT, "need a digital msa");
  if (msa->rf == NULL)                 ESL_XEXCEPTION(eslECONTRACT, "need an RF line");

  ESL_ALLOC(matassign, sizeof(int) * (msa->alen+1));
  esl_vec_ISet(matassign, msa->alen+1, 0);
 
  /* Watch for off-by-one. rf is [0..alen-1]; matassign is [1..alen] */
  for (apos = 0; apos < msa->alen; apos++)
    if (!esl_abc_CIsGap(msa->abc, msa->rf[apos]))
      matassign[apos+1] |= p7_ASSIGN_MATCH;
    else
      matassign[apos+1] |= p7_ASSIGN_INSERT;

  /* matassign2hmm leaves ret_hmm, ret_tr in their proper state: */
  if ((status = matassign2hmm(msa, matassign, ret_hmm, ret_tr)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  if (matassign != NULL) free(matassign);
  return status;
}

/* Function: p7_Fastmodelmaker()
 * 
 * Purpose:  Heuristic model construction.
 *           Construct an HMM from an alignment by a simple rule,
 *           based on the fractional occupancy of each columns w/
 *           residues vs gaps. Any column w/ a fractional
 *           occupancy of $\geq$ <symfrac> is assigned as a MATCH column;
 *           for instance, if thresh = 0.5, columns w/ $\geq$ 50\% 
 *           residues are assigned to match... roughly speaking.
 *           
 *           "Roughly speaking" because sequences may be weighted
 *           in the input <msa>, and because missing data symbols are
 *           ignored, in order to deal with sequence fragments.
 *
 *           The <msa> must be in digital mode. 
 *
 *           If the caller wants to designate any sequences as
 *           fragments, it does so by converting all N-terminal and
 *           C-terminal flanking gap symbols to missing data symbols.
 *
 *           NOTE: p7_Fastmodelmaker() will slightly revise the
 *           alignment if the assignment of columns implies
 *           DI and ID transitions.
 *           
 *           Returns the HMM in counts form (ready for applying Dirichlet
 *           priors as the next step). Also returns fake traceback
 *           for each training sequence.
 *           
 *           Models must have at least one node, so if the <msa> defined 
 *           no consensus columns, a <eslENORESULT> error is returned.
 *           
 * Args:     msa       - multiple sequence alignment
 *           symfrac   - threshold for residue occupancy; >= assigns MATCH
 *           ret_hmm   - RETURN: counts-form HMM
 *           ret_tr    - RETURN: array of tracebacks for aseq's
 *           
 * Return:   <eslOK> on success. ret_hmm and ret_tr allocated here,
 *           and must be free'd by the caller (FreeTrace(tr[i]), free(tr),
 *           FreeHMM(hmm)).       
 *
 *           Returns <eslENORESULT> if no consensus columns were annotated;
 *           in this case, <ret_hmm> and <ret_tr> are returned NULL.
 *           
 * Throws:   <eslEMEM> on allocation failure; <eslECONTRACT> if the 
 *           <msa> isn't in digital mode.
 */
int
p7_Fastmodelmaker(ESL_MSA *msa, float symfrac, P7_HMM **ret_hmm, P7_TRACE ***ret_tr)
{
  int      status;	     /* return status flag                  */
  int     *matassign = NULL; /* MAT state assignments if 1; 1..alen */
  int      idx;              /* counter over sequences              */
  int      apos;             /* counter for aligned columns         */
  float    r;		     /* weighted residue count              */
  float    totwgt;	     /* weighted residue+gap count          */

  *ret_hmm = NULL;
  *ret_tr  = NULL;
  if (! (msa->flags & eslMSA_DIGITAL)) ESL_XEXCEPTION(eslECONTRACT, "need digital MSA");

  /* Allocations: matassign is 1..alen array of bit flags.
   */
  ESL_ALLOC(matassign, sizeof(int)     * (msa->alen+1));
  esl_vec_ISet(matassign, msa->alen+1, 0);  

  /* Determine weighted sym freq in each column, set matassign[] accordingly.
   */
  for (apos = 1; apos <= msa->alen; apos++) 
    {  
      r = totwgt = 0.;
      for (idx = 0; idx < msa->nseq; idx++) 
	{
	  if       (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) { r += msa->wgt[idx]; totwgt += msa->wgt[idx]; }
	  else if  (esl_abc_XIsGap(msa->abc,     msa->ax[idx][apos])) {                     totwgt += msa->wgt[idx]; }
	  else if  (esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos])) continue;
	}
      if (r > 0. && r / totwgt >= symfrac) matassign[apos] |= p7_ASSIGN_MATCH;
      else                                 matassign[apos] |= p7_ASSIGN_INSERT;
    }

  /* Once we have matassign calculated, modelmakers behave
   * the same; matassign2hmm() does this stuff (traceback construction,
   * trace counting) and sets up ret_hmm and ret_tr.
   */
  if ((status = matassign2hmm(msa, matassign, ret_hmm, ret_tr)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  if (matassign != NULL) free(matassign);
  return status;
}

/*-------------------- end, exported API -------------------------*/




/*****************************************************************
 * 2. Private functions used in constructing models.
 *****************************************************************/ 

/* Function: matassign2hmm()
 * 
 * Purpose:  Given an assignment of alignment columns to match vs.
 *           insert, finish the final part of the model construction 
 *           calculation that is constant between model construction
 *           algorithms.
 *           
 * Args:     msa       - multiple sequence alignment
 *           matassign - 1..alen bit flags for column assignments
 *           ret_hmm   - RETURN: counts-form HMM
 *           ret_tr    - RETURN: array of tracebacks for aseq's
 *                         
 * Return:   <eslOK> on success.
 *           <eslENORESULT> if no consensus columns are identified.
 *
 *           ret_hmm and ret_tr alloc'ed here.
 */
static int
matassign2hmm(ESL_MSA *msa, int *matassign, P7_HMM **ret_hmm, P7_TRACE ***ret_tr)
{
  int        status;		/* return status                       */
  P7_HMM    *hmm = NULL;        /* RETURN: new hmm                     */
  P7_TRACE **tr  = NULL;        /* RETURN: 0..nseq-1 fake traces       */
  int      M;                   /* length of new model in match states */
  int      idx;                 /* counter over sequences              */
  int      apos;                /* counter for aligned columns         */
  int      mode;

  *ret_hmm = NULL;
  *ret_tr  = NULL;

  /* How many match states in the HMM? */
  for (M=0,apos = 1; apos <= msa->alen; apos++) 
    if (matassign[apos] & p7_ASSIGN_MATCH) M++;
  if (M == 0) { status = eslENORESULT; goto ERROR; }

  /* Delimit N-terminal tail */
  for (apos=1; matassign[apos] & p7_ASSIGN_INSERT && apos <= msa->alen; apos++)
    matassign[apos] |= p7_EXTERNAL_INSERT_N;
  if (apos <= msa->alen) matassign[apos] |= p7_FIRST_MATCH;

  /* Delimit C-terminal tail */
  for (apos=msa->alen; matassign[apos] & p7_ASSIGN_INSERT && apos > 0; apos--)
    matassign[apos] |= p7_EXTERNAL_INSERT_C;
  if (apos > 0) matassign[apos] |= p7_LAST_MATCH;

  /* Make fake tracebacks for each seq */
  if ((status = fake_tracebacks(msa, matassign, &tr)) != eslOK) goto ERROR;

  /* Build count model from tracebacks */
  if ((hmm    = p7_hmm_Create(M, msa->abc)) == NULL)  { status = eslEMEM; goto ERROR; }
  if ((status = p7_hmm_Zero(hmm))           != eslOK) goto ERROR;
  for (idx = 0; idx < msa->nseq; idx++) {
    if (tr[idx] == NULL) continue; /* skip rare examples of empty sequences */
    if ((status = p7_trace_Count(hmm, msa->ax[idx], msa->wgt[idx], tr[idx])) != eslOK) goto ERROR;
  }

  /* Transfer annotation from the MSA to the new model
   */
  if ((status = annotate_model(hmm, matassign, msa)) != eslOK) goto ERROR;

  /* Reset #=RF line of alignment to reflect our assignment
   * of match, delete. matassign is valid from 1..alen and is off
   * by one from msa->rf.
   */
  if (msa->rf == NULL)  ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen + 1));
  for (apos = 0; apos < msa->alen; apos++)
    msa->rf[apos] = matassign[apos+1] & p7_ASSIGN_MATCH ? 'x' : '.';
  msa->rf[msa->alen] = '\0';

  if (ret_tr  != NULL) *ret_tr  = tr;
  if (ret_hmm != NULL) *ret_hmm = hmm;
  return eslOK;

 ERROR:
  if (tr  != NULL) p7_trace_DestroyArray(tr, msa->nseq);
  if (hmm != NULL) p7_hmm_Destroy(hmm);
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
 *           Rarely, an individual traceback may be left as NULL. This
 *           happens when a sequence has no residues at all. Plan7 models
 *           cannot emit L=0 sequences. Caller must watch out for NULL
 *           traces in the trace array.
 *           
 * Args:     msa       - digital alignment
 *           matassign - assignment of column; [1..alen] 
 *           ret_tr    - RETURN: array of tracebacks
 *           
 * Return:   <eslOK> on success.
 *           ret_tr array is alloc'ed here. Caller must free.
 */          
static int
fake_tracebacks(ESL_MSA *msa, int *matassign, P7_TRACE ***ret_tr)
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
      
      ESL_TRY( p7_trace_Append(tr[idx], p7_STS, 0, 0) ); /* traces start with S... */
      ESL_TRY( p7_trace_Append(tr[idx], p7_STN, 0, 0) ); /*...and transit to N.    */

      i = 1;
      k = 0;
      for (apos = 1; apos <= alen; apos++)
        {
	  if (matassign[apos] & p7_FIRST_MATCH) 	/* BEGIN */
	    ESL_TRY( p7_trace_Append(tr[idx], p7_STB, 0, 0) );

	  if (matassign[apos+1] & p7_ASSIGN_MATCH && ! esl_abc_XIsGap(abc, msa->ax[idx][apos]))
	    {			/* MATCH */
	      k++;		/* move to next model pos */
	      ESL_TRY( p7_trace_Append(tr[idx], p7_STM, k, i) );
	      i++;
	    }	      
          else if (matassign[apos+1] & p7_ASSIGN_MATCH)
            {                   /* DELETE */ /* No B->D transitions */
	      k++;		/* *always* move on model when ASSIGN_MATCH */
	      if (tr[idx]->st[tr[idx]->N-1] != p7_STB)
		ESL_TRY( p7_trace_Append(tr[idx], p7_STD, k, 0) );
            }
          else if (matassign[apos+1] & p7_EXTERNAL_INSERT_N &&
		   ! esl_abc_CIsGap(abc, aseq[idx][apos]))
            {                   /* N-TERMINAL TAIL */
	      ESL_TRY( p7_trace_Append(tr[idx], p7_STN, 0, i) );
	      i++;
            }
	  else if (matassign[apos+1] & p7_EXTERNAL_INSERT_C &&
		   ! esl_abc_CIsGap(abc, aseq[idx][apos]))
	    {			/* C-TERMINAL TAIL */
	      ESL_TRY( p7_trace_Append(tr[idx], p7_STC, 0, i) );
	      i++;
	    }
	  else if (! esl_abc_CIsGap(abc, aseq[idx][apos]))
	    {			/* INSERT */
	      ESL_TRY( p7_trace_Append(tr[idx], p7_STI, k, i) );
	      i++;
	    }

	  if (matassign[apos+1] & p7_LAST_MATCH)
	    {			/* END */
	      /* be careful about S/W transitions; if we're 
               * a candidate sequence fragment, we need to roll
	       * back over some D's because there's no D->E transition
               * for fragments. k==M right now; don't change it.
	       */
	      while (tr[idx]->st[tr[idx]->N-1] == p7_STD) 
		tr[idx]->N--;	/* gratuitous chumminess with trace */
	      ESL_TRY( p7_trace_Append(tr[idx], p7_STE, 0, 0) );
	      ESL_TRY( p7_trace_Append(tr[idx], p7_STC, 0, 0) );
	    }
        }
      ESL_TRY( p7_trace_Append(tr[idx], p7_STT, 0, 0) );

      /* deal with DI, ID transitions and other plan7 impossibilities */
      /* (k == M at this point) */
      ESL_TRY ( trace_doctor(tr[idx], k, NULL, NULL) );

      /* What happens if we had a length 0 sequence? 
       * At this point, we'd have an invalid SNBECT trace.
       * Check for this; if so, free the trace and leave it NULL.
       */
      if (tr[idx]->st[2] == p7_STB && tr[idx]->st[3] == p7_STE)
	{
	  p7_trace_Destroy(tr[idx]);
	  tr[idx] = NULL;
	}
    } 

  status = eslOK;
 CLEANEXIT:
  if (status != eslOK) { esl_Free2D((void **) tr, nseq); tr = NULL; }
  *ret_tr = tr;
  return status;
}
