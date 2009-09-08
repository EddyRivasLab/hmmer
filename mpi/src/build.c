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


static int matassign2hmm(ESL_MSA *msa, int *matassign, P7_HMM **ret_hmm, P7_TRACE ***opt_tr);
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
 *           opt_tr  - optRETURN: array of tracebacks for aseq's
 *           
 * Return:   <eslOK> on success. <ret_hmm> and <opt_tr> are allocated 
 *           here, and must be free'd by caller.
 *
 *           Returns <eslENORESULT> if no consensus columns were annotated;
 *           in this case, <ret_hmm> and <opt_tr> are returned NULL. 
 *           
 *           Returns <eslEFORMAT> if the <msa> doesn't have a reference
 *           annotation line.
 *           
 * Throws:   <eslEMEM> on allocation failure. Throws <eslEINVAL> if the <msa>
 *           isn't in digital mode.
 */            
int
p7_Handmodelmaker(ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***opt_tr)
{
  int        status;
  int       *matassign = NULL;    /* MAT state assignments if 1; 1..alen */
  int        apos;                /* counter for aligned columns         */

  if (! (msa->flags & eslMSA_DIGITAL)) ESL_XEXCEPTION(eslEINVAL, "need a digital msa");
  if (msa->rf == NULL)                 return eslEFORMAT;

  ESL_ALLOC(matassign, sizeof(int) * (msa->alen+1));
 
  /* Watch for off-by-one. rf is [0..alen-1]; matassign is [1..alen] */
  for (apos = 1; apos <= msa->alen; apos++)
    matassign[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos-1])? FALSE : TRUE);

  /* matassign2hmm leaves ret_hmm, opt_tr in their proper state: */
  if ((status = matassign2hmm(msa, matassign, ret_hmm, opt_tr)) != eslOK) goto ERROR;

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
 *           opt_tr    - optRETURN: array of tracebacks for aseq's
 *           
 * Return:   <eslOK> on success. ret_hmm and opt_tr allocated here,
 *           and must be free'd by the caller (FreeTrace(tr[i]), free(tr),
 *           FreeHMM(hmm)).       
 *
 *           Returns <eslENORESULT> if no consensus columns were annotated;
 *           in this case, <ret_hmm> and <opt_tr> are returned NULL.
 *           
 * Throws:   <eslEMEM> on allocation failure; <eslEINVAL> if the 
 *           <msa> isn't in digital mode.
 */
int
p7_Fastmodelmaker(ESL_MSA *msa, float symfrac, P7_HMM **ret_hmm, P7_TRACE ***opt_tr)
{
  int      status;	     /* return status flag                  */
  int     *matassign = NULL; /* MAT state assignments if 1; 1..alen */
  int      idx;              /* counter over sequences              */
  int      apos;             /* counter for aligned columns         */
  float    r;		     /* weighted residue count              */
  float    totwgt;	     /* weighted residue+gap count          */

  if (! (msa->flags & eslMSA_DIGITAL)) ESL_XEXCEPTION(eslEINVAL, "need digital MSA");

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
   * trace counting) and sets up ret_hmm and opt_tr.
   */
  if ((status = matassign2hmm(msa, matassign, ret_hmm, opt_tr)) != eslOK) goto ERROR;

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
 *           opt_tr    - optRETURN: array of tracebacks for aseq's
 *                         
 * Return:   <eslOK> on success.
 *           <eslENORESULT> if no consensus columns are identified.
 *
 *           ret_hmm and opt_tr alloc'ed here.
 */
static int
matassign2hmm(ESL_MSA *msa, int *matassign, P7_HMM **ret_hmm, P7_TRACE ***opt_tr)
{
  int        status;		/* return status                       */
  P7_HMM    *hmm = NULL;        /* RETURN: new hmm                     */
  P7_TRACE **tr  = NULL;        /* RETURN: 0..nseq-1 fake traces       */
  int      M;                   /* length of new model in match states */
  int      idx;                 /* counter over sequences              */
  int      apos;                /* counter for aligned columns         */
  char errbuf[eslERRBUFSIZE];

  /* How many match states in the HMM? */
  for (M = 0, apos = 1; apos <= msa->alen; apos++) 
    if (matassign[apos]) M++;
  if (M == 0) { status = eslENORESULT; goto ERROR; }

  /* Make fake tracebacks for each seq */
  ESL_ALLOC(tr, sizeof(P7_TRACE *) * msa->nseq);
  if ((status = p7_trace_FauxFromMSA(msa, matassign, p7_MSA_COORDS, tr))        != eslOK) goto ERROR;
  for (idx = 0; idx < msa->nseq; idx++)
    {
      if ((status = p7_trace_Doctor(tr[idx], NULL, NULL))                       != eslOK) goto ERROR;
      if ((status = p7_trace_Validate(tr[idx], msa->abc, msa->ax[idx], errbuf)) != eslOK) 
	ESL_XEXCEPTION(eslFAIL, "validation failed: %s", errbuf);
    }

  /* Build count model from tracebacks */
  if ((hmm    = p7_hmm_Create(M, msa->abc)) == NULL)  { status = eslEMEM; goto ERROR; }
  if ((status = p7_hmm_Zero(hmm))           != eslOK) goto ERROR;
  for (idx = 0; idx < msa->nseq; idx++) {
    if (tr[idx] == NULL) continue; /* skip rare examples of empty sequences */
    if ((status = p7_trace_Count(hmm, msa->ax[idx], msa->wgt[idx], tr[idx])) != eslOK) goto ERROR;
  }
  hmm->nseq     = msa->nseq;
  hmm->eff_nseq = msa->nseq;

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

  if (opt_tr  != NULL) *opt_tr  = tr; 
  else                  p7_trace_DestroyArray(tr, msa->nseq);
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  if (tr     != NULL) p7_trace_DestroyArray(tr, msa->nseq);
  if (hmm    != NULL) p7_hmm_Destroy(hmm);
  if (opt_tr != NULL) *opt_tr = NULL;
  *ret_hmm = NULL;
  return status;
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
    hmm->flags |= p7H_RF;
  }

  /* Consensus structure annotation */
  if (msa->ss_cons != NULL) {
    ESL_ALLOC(hmm->cs, sizeof(char) * (hmm->M+2));
    hmm->cs[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos]) hmm->cs[k++] = msa->ss_cons[apos-1];
    hmm->cs[k] = '\0';
    hmm->flags |= p7H_CS;
  }

  /* Surface accessibility annotation */
  if (msa->sa_cons != NULL) {
    ESL_ALLOC(hmm->ca, sizeof(char) * (hmm->M+2));
    hmm->ca[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos]) hmm->ca[k++] = msa->sa_cons[apos-1];
    hmm->ca[k] = '\0';
    hmm->flags |= p7H_CA;
  }

  /* The alignment map (1..M in model, 1..alen in alignment) */
  ESL_ALLOC(hmm->map, sizeof(int) * (hmm->M+1));
  hmm->map[0] = 0;
  for (apos = k = 1; apos <= msa->alen; apos++)
    if (matassign[apos]) hmm->map[k++] = apos;
  hmm->flags |= p7H_MAP;

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
  char          msafile[16] = "p7tmpXXXXXX"; /* tmpfile template */
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
 *   1. seq2 gives an I->end transition.
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


