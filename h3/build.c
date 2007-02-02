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
 *    3. Unit tests.
 *    4. Test driver.
 *    5. Copyright and license.
 * 
 * SRE, Tue Jan 2 2007 [Casa de Gatos]
 * SVN $Id$
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

#include "hmmer.h"


static int matassign2hmm(ESL_MSA *msa, int *matassign, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);
static int fake_tracebacks(ESL_MSA *msa, int *matassign, P7_TRACE ***ret_tr);
static int trace_doctor(P7_TRACE *tr, int mlen, int *ret_ndi, int *ret_nid);
static int annotate_model(P7_HMM *hmm, int *matassign, ESL_MSA *msa);

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
 *           ret_tr  - optRETURN: array of tracebacks for aseq's
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
  int        status;
  int       *matassign = NULL;    /* MAT state assignments if 1; 1..alen */
  int        apos;                /* counter for aligned columns         */

  if (! (msa->flags & eslMSA_DIGITAL)) ESL_XEXCEPTION(eslECONTRACT, "need a digital msa");
  if (msa->rf == NULL)                 ESL_XEXCEPTION(eslECONTRACT, "need an RF line");

  ESL_ALLOC(matassign, sizeof(int) * (msa->alen+1));
 
  /* Watch for off-by-one. rf is [0..alen-1]; matassign is [1..alen] */
  for (apos = 1; apos <= msa->alen; apos++)
    matassign[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos-1])? FALSE : TRUE);

  /* matassign2hmm leaves ret_hmm, ret_tr in their proper state: */
  if ((status = matassign2hmm(msa, matassign, ret_hmm, ret_tr)) != eslOK) goto ERROR;

  free(matassign);
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
 *           ret_tr    - optRETURN: array of tracebacks for aseq's
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

  if (! (msa->flags & eslMSA_DIGITAL)) ESL_XEXCEPTION(eslECONTRACT, "need digital MSA");

  /* Allocations: matassign is 1..alen array of bit flags.
   */
  ESL_ALLOC(matassign, sizeof(int)     * (msa->alen+1));

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
      if (r > 0. && r / totwgt >= symfrac) matassign[apos] = TRUE;
      else                                 matassign[apos] = FALSE;
    }

  /* Once we have matassign calculated, modelmakers behave
   * the same; matassign2hmm() does this stuff (traceback construction,
   * trace counting) and sets up ret_hmm and ret_tr.
   */
  if ((status = matassign2hmm(msa, matassign, ret_hmm, ret_tr)) != eslOK) goto ERROR;

  free(matassign);
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
 *           ret_tr    - optRETURN: array of tracebacks for aseq's
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

  /* How many match states in the HMM? */
  for (M = 0, apos = 1; apos <= msa->alen; apos++) 
    if (matassign[apos]) M++;
  if (M == 0) { status = eslENORESULT; goto ERROR; }

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
  for (apos = 1; apos <= msa->alen; apos++)
    msa->rf[apos-1] = matassign[apos] ? 'x' : '.';
  msa->rf[msa->alen] = '\0';

  if (ret_tr  != NULL) *ret_tr  = tr; 
  else                  p7_trace_DestroyArray(tr, msa->nseq);
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  if (tr     != NULL) p7_trace_DestroyArray(tr, msa->nseq);
  if (hmm    != NULL) p7_hmm_Destroy(hmm);
  if (ret_tr != NULL) *ret_tr = NULL;
  *ret_hmm = NULL;
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
 *           These faked tracebacks use sequence indices <i> relative
 *           to the aligned sequences in msa, not unaligned seqs as normal.
 *           
 *           Rarely, an individual traceback may be left as NULL. This
 *           happens when a sequence has no residues assigned to the
 *           model (or, more rarely, an entirely empty sequence of
 *           length 0, which is possible in many alignment formats).
 *           Plan7 models cannot emit L=0 sequences. Caller must watch
 *           out for NULL traces in the trace array.
 *           
 * Args:     msa       - digital alignment
 *           matassign - assignment of columns; [1..alen] 
 *           ret_tr    - RETURN: array of tracebacks
 *           
 * Return:   <eslOK> on success.
 *           ret_tr array is alloc'ed here. Caller must free.
 *           
 * Throws:   <eslEMEM> on allocation error.
 */          
static int
fake_tracebacks(ESL_MSA *msa, int *matassign, P7_TRACE ***ret_tr)
{
  int  status;
  P7_TRACE **tr = NULL;
  int  idx;                     /* counter over sequences          */
  int  k;                       /* position in HMM                 */
  int  apos;                    /* position in alignment columns   */
  int  k1,k2;			/* first, last match state used    */
  int  c1,c2;			/* entered/exited column           */
  int  cfirst,clast;		/* loc's of p7_{FIRST,LAST}_MATCH  */
  int  M;

  /* Allocate for nseq traces.
   */
  ESL_ALLOC(tr, sizeof(P7_TRACE *) * msa->nseq);
  for (idx = 0; idx < msa->nseq; idx++) tr[idx] = NULL;
  
  /* Precalculate where the first and last match columns are in matassign.
   * We know M>0; matassign2hmm already tested for M=0 pathology.
   */
  cfirst = clast = M = 0;
  for (apos = 1; apos <= msa->alen; apos++) {
    if (matassign[apos]) { 
      M++;
      clast = apos;
      if (cfirst == 0) cfirst = apos;
    }
  }

  /* Construct traces for each sequence. A trace may end up being NULL,
   * if it would have implied a B->E empty sequence that Plan7 can't do.
   */
  for (idx = 0; idx < msa->nseq; idx++)
    {
      /* Preprocess the sequence to identify position of its first
       * match state (node k1, column c1) and last match state
       * (node k2, column c2).
       */
      for (k1 = 0, c1 = cfirst; c1 <= clast; c1++) {
	if (matassign[c1]) {
	  k1++;
	  if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][c1])) break;
	}
      }
      for (k2 = M+1, c2 = clast; c2 >= cfirst; c2--) {
	if (matassign[c2]) {
	  k2--;
	  if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][c2])) break;
	}
      }
      if (c1 > clast || c2 < cfirst) continue; /* rare: empty sequence; leave trace NULL. */

      /* Now build the trace from left to right in three sections 1..c1-1, c1..c2, c2+1..alen */
      /* First section: left flank, residues assigned to N's.
       */
      if ((status = p7_trace_Create(msa->alen+6, &tr[idx])) != eslOK) goto ERROR; /* +6 = S,N,B,E,C,T */
      if ((status = p7_trace_Append(tr[idx], p7_STS, 0, 0)) != eslOK) goto ERROR; /* traces start with S... */
      if ((status = p7_trace_Append(tr[idx], p7_STN, 0, 0)) != eslOK) goto ERROR; /*...and transit to N.    */

      for (apos = 1; apos < c1; apos++)
	if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) 
	  {
	    if ((status = p7_trace_Append(tr[idx], p7_STN, 0, apos)) != eslOK) goto ERROR;
	  }
      if ((status = p7_trace_Append(tr[idx], p7_STB, 0, 0)) != eslOK) goto ERROR;
      if (esl_abc_XIsMissing(msa->abc, msa->ax[idx][cfirst]))
	if ((status = p7_trace_Append(tr[idx], p7_STX, 0, 0)) != eslOK) goto ERROR;
      
      /* Second section: the model itself. */
      k = k1; /* k = index of next node to deal with */
      for (apos = c1; apos <= c2; apos++)
	{
	  if (matassign[apos]) {
	    if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos]))
	      { /* on c1 and c2 endpoints themselves, this must be the case, by def'n of c1/c2 */
		if ((status = p7_trace_Append(tr[idx], p7_STM, k, apos)) != eslOK) goto ERROR;
		k++;
	      }
	    else if (esl_abc_XIsGap(msa->abc, msa->ax[idx][apos]))
	      {
		if ((status = p7_trace_Append(tr[idx], p7_STD, k, 0)) != eslOK) goto ERROR;
		k++;
	      }
	    else if (esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]))
	      {
		if (tr[idx]->st[tr[idx]->N-1] != p7_STX)
		  if ((status == p7_trace_Append(tr[idx], p7_STX, 0, 0)) != eslOK) goto ERROR;
	      }
	    else
	      ESL_XEXCEPTION(eslEINCONCEIVABLE, "can't happen");

	  } else { /* else, an insert-assigned column: */
	    if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos]))
	      {
		if ((status = p7_trace_Append(tr[idx], p7_STI, k-1, apos)) != eslOK) goto ERROR;
	      }
	    else if (esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]))
	      {
		if (tr[idx]->st[tr[idx]->N-1] != p7_STX)
		  if ((status = p7_trace_Append(tr[idx], p7_STX, 0, 0)) != eslOK) goto ERROR;
	      }
	    else if (! esl_abc_XIsGap(msa->abc, msa->ax[idx][apos]))
	      ESL_XEXCEPTION(eslEINCONCEIVABLE, "can't happen");
	  }
	}

      /* Third section: the C-terminal flank
       */
      if (esl_abc_XIsMissing(msa->abc, msa->ax[idx][clast]))
	if ((status = p7_trace_Append(tr[idx], p7_STX, 0, 0)) != eslOK) goto ERROR;
      if ((status = p7_trace_Append(tr[idx], p7_STE, 0, 0)) != eslOK) goto ERROR;
      if ((status = p7_trace_Append(tr[idx], p7_STC, 0, 0)) != eslOK) goto ERROR;

      for (apos = c2+1; apos <= msa->alen; apos++)
	if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) 
	  {
	    if ((status = p7_trace_Append(tr[idx], p7_STC, 0, apos)) != eslOK) goto ERROR;
	  }	
      if ((status = p7_trace_Append(tr[idx], p7_STT, 0, 0)) != eslOK) goto ERROR;

      /* finally, clean up: deal with DI, ID transitions and other plan7 impossibilities */
      if ((status = trace_doctor(tr[idx], M, NULL, NULL)) != eslOK) goto ERROR;

      /* Validate it.
       */
      /*
      p7_trace_Dump(stdout, tr[idx], NULL, NULL); 
      */
      if (p7_trace_Validate(tr[idx], msa->abc, msa->ax[idx], NULL) != eslOK) 
	ESL_XEXCEPTION(eslFAIL, "validation failed");
    } 

  *ret_tr = tr;
  return eslOK;

 ERROR:
  p7_trace_DestroyArray(tr, msa->nseq);
  *ret_tr = NULL;
  return status;
}



/* Function: trace_doctor()
 * 
 * Purpose:  Plan 7 disallows D->I and I->D "chatter" transitions.
 *           However, these transitions will be implied by many
 *           alignments. trace_doctor() arbitrarily collapses I->D or
 *           D->I into a single M position in the trace.
 *           
 *           trace_doctor does not examine any scores when it does
 *           this. In ambiguous situations (D->I->D) the symbol
 *           will be pulled arbitrarily to the left, regardless
 *           of whether that's the best column to put it in or not.
 *           
 * Args:     tr      - [M] trace to doctor
 *           M       - length of model that traces are for 
 *           ret_ndi - optRETURN: number of DI transitions doctored
 *           ret_nid - optRETURN: number of ID transitions doctored
 * 
 * Return:   <eslOK> on success, and the trace <tr> is modified.
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
 * Args:     hmm       - [M] new model to annotate
 *           matassign - which alignment columns are MAT; [1..alen]
 *           msa       - alignment, including annotation to transfer
 *           
 * Return:   <eslOK> on success.
 *
 * Throws:   <eslEMEM> on allocation error.
 */
static int
annotate_model(P7_HMM *hmm, int *matassign, ESL_MSA *msa)
{                      
  int   apos;			/* position in matassign, 1.alen  */
  int   k;			/* position in model, 1.M         */
  int   status;

  /* Reference coord annotation  */
  if (msa->rf != NULL) {
    ESL_ALLOC(hmm->rf, sizeof(char) * (hmm->M+2));
    hmm->rf[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++) 
      if (matassign[apos]) hmm->rf[k++] = msa->rf[apos-1]; /* watch off-by-one in msa's rf */
    hmm->rf[k] = '\0';
    hmm->flags |= p7_RF;
  }

  /* Consensus structure annotation */
  if (msa->ss_cons != NULL) {
    ESL_ALLOC(hmm->cs, sizeof(char) * (hmm->M+2));
    hmm->cs[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos]) hmm->cs[k++] = msa->ss_cons[apos-1];
    hmm->cs[k] = '\0';
    hmm->flags |= p7_CS;
  }

  /* Surface accessibility annotation */
  if (msa->sa_cons != NULL) {
    ESL_ALLOC(hmm->ca, sizeof(char) * (hmm->M+2));
    hmm->ca[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos]) hmm->ca[k++] = msa->sa_cons[apos-1];
    hmm->ca[k] = '\0';
    hmm->flags |= p7_CA;
  }

  /* The alignment map (1..M in model, 1..alen in alignment) */
  ESL_ALLOC(hmm->map, sizeof(int) * (hmm->M+1));
  hmm->map[0] = 0;
  for (apos = k = 1; apos <= msa->alen; apos++)
    if (matassign[apos]) hmm->map[k++] = apos;
  hmm->flags |= p7_MAP;

  return eslOK;

 ERROR:
  return status;
}

/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7BUILD_TESTDRIVE

#if 0
static void
utest_foo(void)
{

  return;
}
#endif

#endif /*p7BUILD_TESTDRIVE*/
/*---------------------- end of unit tests -----------------------*/




/*****************************************************************
 * 5. Test driver.
 *****************************************************************/

#ifdef p7BUILD_TESTDRIVE
/* gcc -g -Wall -Dp7BUILD_TESTDRIVE -I. -I../easel -L. -L../easel -o build_test build.c -lhmmer -leasel -lm
 */
#include "easel.h"

#include "p7_config.h"
#include "hmmer.h"

static int write_test_msa(FILE *ofp);

int
main(int argc, char **argv)
{  
  char          msafile[10] = "eslXXXXXX"; /* tmpfile template */
  FILE         *fp;        
  ESL_ALPHABET *abc;
  ESL_MSAFILE  *afp;
  ESL_MSA      *msa;
  P7_HMM       *hmm;
  int           nali;
  float         symfrac = 0.5;
  int           status;


  if ((status = esl_tmpfile_named(msafile, &fp)) != eslOK) esl_fatal("tmpfile creation failed");
  if ((status = write_test_msa(fp))              != eslOK) esl_fatal("test msa write failed");
  fclose(fp);

  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("alphabet creation failed");

  status = esl_msafile_OpenDigital(abc, msafile, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) esl_fatal("msa open failed");

  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;
      
      status = p7_Fastmodelmaker(msa, symfrac, &hmm, NULL);
      if (status != eslOK) esl_fatal("Model construction failed.");

      p7_hmm_Destroy(hmm);
      esl_msa_Destroy(msa);
    }
  if (status != eslEOF) esl_fatal("Alignment read failed");

  esl_msafile_Close(afp);
  esl_alphabet_Destroy(abc);
  remove(msafile);
  return eslOK;
}


/* 
 * An MSA to ex{e,o}rcise past demons.
 *   1. seq2 gives an I->end transition, that construction must
 *      recognize as a C tail instead.
 *   2. seq1 contains degenerate Z,X, exercising symbol counting
 *      of degenerate residues.
 */
static int
write_test_msa(FILE *ofp)
{
  fprintf(ofp, "# STOCKHOLM 1.0\n");
  fprintf(ofp, "#=GC RF --xxxxxxxxxxxxxxxx-xxx-x--\n");
  fprintf(ofp, "seq1    --ACDEFGHIKLMNPZXS-TVW-Yyy\n");
  fprintf(ofp, "seq2    aaACDEFGHIKLMNPQRS-TVWw---\n");
  fprintf(ofp, "seq3    aaAC-EFGHIKLMNPQRS-TVW-Y--\n");
  fprintf(ofp, "seq4    aaAC-EFGHIKLMNPQRS-TVW-Y--\n");
  fprintf(ofp, "//\n");
  return eslOK;
}


#endif /*p7BUILD_TESTDRIVE*/
/*-------------------- end of test driver ---------------------*/


/************************************************************
 * @LICENSE@
 ************************************************************/


