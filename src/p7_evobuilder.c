/* Standardized pipeline for construction of new HMMs.
 * 
 * Contents:
 *    1. P7_BUILDER: allocation, initialization, destruction
 *    2. Standardized model construction API.
 *    3. Internal functions.
 */   
#include <p7_config.h>

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msacluster.h"
#include "esl_msaweight.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "evo_evalues.h"


/*****************************************************************
 * 2. Standardized model construction API.
 *****************************************************************/

static int    validate_msa         (P7_BUILDER *bld, ESL_MSA *msa);
static int    relative_weights     (P7_BUILDER *bld, ESL_MSA *msa);
static int    build_model          (P7_BUILDER *bld, ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***opt_tr);
static int    effective_seqnumber  (P7_BUILDER *bld, const ESL_MSA *msa, P7_HMM *hmm, const P7_BG *bg);
static int    parameterize         (P7_BUILDER *bld, P7_HMM *hmm);
static int    annotate             (P7_BUILDER *bld, const ESL_MSA *msa, P7_HMM *hmm);
static int    evo_calibrate        (P7_BUILDER *bld, P7_RATE *R, P7_HMM *hmm, P7_BG *bg, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om,
				    int noevo_msv, int noevo_vit, int noevo_fwd,float fixtime, float tol);
static int    make_post_msa        (P7_BUILDER *bld, const ESL_MSA *premsa, const P7_HMM *hmm, P7_TRACE **tr, ESL_MSA **opt_postmsa);

/* Function:  p7_EvoBuilder()
 * Synopsis:  Build a new HMM from an MSA.
 *
 * Purpose:   Take the multiple sequence alignment <msa> and a build configuration <bld>,
 *            and build a new HMM. 
 * 
 *            Effective sequence number determination and calibration steps require
 *            additionally providing a null model <bg>.
 *
 * Args:      bld         - build configuration
 *            msa         - multiple sequence alignment
 *            bg          - null model
 *            opt_hmm     - optRETURN: new HMM
 *            opt_trarr   - optRETURN: array of faux tracebacks, <0..nseq-1>
 *            opt_gm      - optRETURN: profile corresponding to <hmm>
 *            opt_om      - optRETURN: optimized profile corresponding to <gm>
 *            opt_postmsa - optRETURN: RF-annotated, possibly modified MSA 
 *
 * Returns:   <eslOK> on success. The new HMM is optionally returned in
 *            <*opt_hmm>, along with optional returns of an array of faux tracebacks
 *            for each sequence in <*opt_trarr>, the annotated MSA used to construct
 *            the model in <*opt_postmsa>, a configured search profile in 
 *            <*opt_gm>, and an optimized search profile in <*opt_om>. These are
 *            all optional returns because the caller may, for example, be interested
 *            only in an optimized profile, or may only be interested in the HMM.
 *            
 *            Returns <eslENORESULT> if no consensus columns were annotated.
 *            Returns <eslEFORMAT> on MSA format problems, such as a missing RF annotation
 *            line in hand architecture construction. On any returned error,
 *            <bld->errbuf> contains an informative error message.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if relative weights couldn't be calculated from <msa>.
 *
 * Xref:      J4/30.
 */
int
p7_EvoBuilder(P7_BUILDER *bld, ESL_MSA *msa, P7_BG *bg, HMMRATE *hmmrate,
	      P7_HMM **opt_hmm, P7_TRACE ***opt_trarr, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om,
	      ESL_MSA **opt_postmsa, int noevo_msv, int noevo_vit, int noevo_fwd)
{
  int i,j;
  uint32_t    checksum = 0;	/* checksum calculated for the input MSA. hmmalign --mapali verifies against this. */
  P7_HMM     *hmm      = NULL;
  P7_RATE    *R        = NULL;
  P7_TRACE  **tr       = NULL;
  P7_TRACE ***tr_ptr   = (opt_trarr != NULL || opt_postmsa != NULL) ? &tr : NULL;
  int         status;
  if ((status =  validate_msa         (bld, msa))                       != eslOK) goto ERROR;
  if ((status =  esl_msa_Checksum     (msa, &checksum))                 != eslOK) ESL_XFAIL(status, bld->errbuf, "Failed to calculate checksum"); 
  if ((status =  relative_weights     (bld, msa))                       != eslOK) goto ERROR;
  if ((status =  esl_msa_MarkFragments_old(msa, bld->fragthresh))       != eslOK) goto ERROR;
  if ((status =  build_model          (bld, msa, &hmm, tr_ptr))         != eslOK) goto ERROR;
  
  //Ensures that the weighted-average I->I count <=  bld->max_insert_len
  //(MI currently contains the number of observed insert-starts)
  if (bld->max_insert_len>0)
    for (i=1; i<hmm->M; i++ )
      hmm->t[i][p7H_II] = ESL_MIN(hmm->t[i][p7H_II], bld->max_insert_len*hmm->t[i][p7H_MI]);

  if ((status =  effective_seqnumber  (bld, msa, hmm, bg))                                    != eslOK) goto ERROR;
  if ((status =  parameterize         (bld, hmm))                                             != eslOK) goto ERROR;
  if ((status =  annotate             (bld, msa, hmm))                                        != eslOK) goto ERROR;

  if (!noevo_fwd) 
    if ((status = p7_RateConstruct(hmm, bg, hmmrate, &R, NULL, FALSE)) != eslOK) goto ERROR;

  if ((status =  evo_calibrate(bld, R, hmm, bg, opt_gm, opt_om, noevo_msv, noevo_vit, noevo_fwd, hmmrate->fixtime, hmmrate->tol)) != eslOK) goto ERROR;
  if ((status =  make_post_msa(bld, msa, hmm, tr, opt_postmsa))                                                                   != eslOK) goto ERROR;

  //force masked positions to background  (it'll be close already, so no relevant impact on weighting)
  if (hmm->mm != NULL)
    for (i=1; i<hmm->M; i++ )
      if (hmm->mm[i] == 'm')
        for (j=0; j<hmm->abc->K; j++)
          hmm->mat[i][j] = bg->f[j];

  if ( bld->abc->type == eslDNA ||  bld->abc->type == eslRNA ) {
	  if (bld->w_len > 0)           hmm->max_length = bld->w_len;
	  else if (bld->w_beta == 0.0)  hmm->max_length = hmm->M *4;
	  else if ( (status =  p7_Builder_MaxLength(hmm, bld->w_beta)) != eslOK) goto ERROR;
  }

  hmm->checksum = checksum;
  hmm->flags   |= p7H_CHKSUM;

  if (opt_hmm   != NULL) *opt_hmm   = hmm; else p7_hmm_Destroy(hmm);
  if (opt_trarr != NULL) *opt_trarr = tr;  else p7_trace_DestroyArray(tr, msa->nseq);
  p7_RateDestroy(R);
  return eslOK;

 ERROR:
  p7_hmm_Destroy(hmm);
  p7_RateDestroy(R);
  p7_trace_DestroyArray(tr, msa->nseq);
  if (opt_gm    != NULL) p7_profile_Destroy(*opt_gm);
  if (opt_om    != NULL) p7_oprofile_Destroy(*opt_om);
  return status;
}


/* Function:  p7_SingleBuilder()
 * Synopsis:  Build a new HMM from a single sequence.
 *
 * Purpose:   Take the sequence <sq> and a build configuration <bld>, and
 *            build a new HMM.
 *            
 *            The single sequence scoring system in the <bld>
 *            configuration must have been previously initialized by
 *            <p7_builder_SetScoreSystem()>.
 *            
 * Args:      bld       - build configuration
 *            sq        - query sequence
 *            bg        - null model (needed to parameterize insert emission probs)
 *            opt_hmm   - optRETURN: new HMM
 *            opt_tr    - optRETURN: faux trace of <sq> against core <hmm>
 *            opt_gm    - optRETURN: profile corresponding to <hmm>
 *            opt_om    - optRETURN: optimized profile corresponding to <gm>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if <bld> isn't properly configured somehow.
 */
int
p7_EvoSingleBuilder(P7_BUILDER *bld, ESL_SQ *sq, P7_BG *bg, HMMRATE *hmmrate,
		    P7_HMM **opt_hmm, P7_TRACE **opt_tr, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om,
		    int noevo_msv, int noevo_vit, int noevo_fwd)
{
  P7_HMM   *hmm = NULL;
  P7_RATE  *R   = NULL;
  P7_TRACE *tr  = NULL;
  int       k;
  int       status;
  
  bld->errbuf[0] = '\0';
  if (! bld->Q) ESL_XEXCEPTION(eslEINVAL, "score system not initialized");

  if ((status = p7_Seqmodel(bld->abc, sq->dsq, sq->n, sq->name, bld->Q, bg->f, bld->popen, bld->pextend, &hmm))   != eslOK) goto ERROR;
  if ((status = p7_hmm_SetComposition(hmm))                                                                       != eslOK) goto ERROR;
  if ((status = p7_hmm_SetConsensus(hmm, sq))                                                                     != eslOK) goto ERROR;
  if (!noevo_fwd) 
    if ((status = p7_RateConstruct(hmm, bg, hmmrate, &R, NULL, FALSE)) != eslOK) goto ERROR;
  if ((status = evo_calibrate(bld, R, hmm, bg, opt_gm, opt_om, noevo_msv, noevo_vit, noevo_fwd, hmmrate->fixtime, hmmrate->tol)) != eslOK) goto ERROR;

  if ( bld->abc->type == eslDNA ||  bld->abc->type == eslRNA ) {
    if (bld->w_len > 0)           hmm->max_length = bld->w_len;
    else if (bld->w_beta == 0.0)  hmm->max_length = hmm->M *4;
    else if ( (status =  p7_Builder_MaxLength(hmm, bld->w_beta)) != eslOK) goto ERROR;
  }


  /* build a faux trace: relative to core model (B->M_1..M_L->E) */
  if (opt_tr != NULL) 
    {
      if ((tr = p7_trace_Create())                      == NULL)  goto ERROR;
      if ((status = p7_trace_Append(tr, p7T_B, 0, 0))   != eslOK) goto ERROR; 
      for (k = 1; k <= sq->n; k++)
        if ((status = p7_trace_Append(tr, p7T_M, k, k)) != eslOK) goto ERROR;
      if ((status = p7_trace_Append(tr, p7T_E, 0, 0))   != eslOK) goto ERROR; 
      tr->M = sq->n;
      tr->L = sq->n;
    }

  /* note that <opt_gm> and <opt_om> were already set by calibrate() call above. */
  if (opt_hmm   != NULL) *opt_hmm = hmm; else p7_hmm_Destroy(hmm);
  if (opt_tr    != NULL) *opt_tr  = tr;

  p7_RateDestroy(R);
  return eslOK;

 ERROR:
  p7_RateDestroy(R);
  p7_hmm_Destroy(hmm);
  if (tr        != NULL) p7_trace_Destroy(tr);
  if (opt_gm    != NULL) p7_profile_Destroy(*opt_gm);
  if (opt_om    != NULL) p7_oprofile_Destroy(*opt_om);
  return status;
}



/*------------- end, model construction API ---------------------*/




/*****************************************************************
 * 3. Internal functions
 *****************************************************************/


/* validate_msa:
 * SRE, Thu Dec  3 16:10:31 2009 [J5/119; bug #h70 fix]
 * 
 * HMMER uses a convention for missing data characters: they
 * indicate that a sequence is a fragment.  (See
 * esl_msa_MarkFragments_old()).
 *
 * Because of the way these fragments will be handled in tracebacks,
 * we reject any alignment that uses missing data characters in any
 * other way.
 * 
 * This validation step costs negligible time.
 */
static int
validate_msa(P7_BUILDER *bld, ESL_MSA *msa)
{
  int     idx;
  int64_t apos;

  for (idx = 0; idx < msa->nseq; idx++)
    {
      apos = 1;
      while (  esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]) && apos <= msa->alen) apos++;
      while (! esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]) && apos <= msa->alen) apos++;
      while (  esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]) && apos <= msa->alen) apos++;
      if (apos != msa->alen+1) ESL_FAIL(eslEINVAL, bld->errbuf, "msa %s; sequence %s\nhas missing data chars (~) other than at fragment edges", msa->name, msa->sqname[idx]);
    }
  
  return eslOK;
}


/* set_relative_weights():
 * Set msa->wgt vector, using user's choice of relative weighting algorithm.
 */
static int
relative_weights(P7_BUILDER *bld, ESL_MSA *msa)
{
  ESL_MSAWEIGHT_CFG *cfg    = esl_msaweight_cfg_Create();
  int                status = eslOK;

  cfg->ignore_rf = (bld->arch_strategy == p7_ARCH_HAND ? FALSE : TRUE);  // PB weights use RF consensus columns only if --hand is set. [iss #180]

  if      (bld->wgt_strategy == p7_WGT_NONE)                    { esl_vec_DSet(msa->wgt, msa->nseq, 1.); }
  else if (bld->wgt_strategy == p7_WGT_GIVEN)                   ;
  else if (bld->wgt_strategy == p7_WGT_PB)                      status = esl_msaweight_PB_adv(cfg, msa, /*ESL_MSAWEIGHT_DAT=*/ NULL); 
  else if (bld->wgt_strategy == p7_WGT_GSC)                     status = esl_msaweight_GSC(msa); 
  else if (bld->wgt_strategy == p7_WGT_BLOSUM)                  status = esl_msaweight_BLOSUM(msa, bld->wid); 
  else ESL_EXCEPTION(eslEINCONCEIVABLE, "no such weighting strategy");

  if (status != eslOK) ESL_FAIL(status, bld->errbuf, "failed to set relative weights in alignment");
  esl_msaweight_cfg_Destroy(cfg);
  return eslOK;
}


/* build_model():
 * Given <msa>, choose HMM architecture, collect counts;
 * upon return, <*ret_hmm> is newly allocated and contains
 * relative-weighted observed counts.
 * Optionally, caller can request an array of inferred traces for
 * the <msa> too.
 */
static int
build_model(P7_BUILDER *bld, ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***opt_tr)
{
  int status;

  if      (bld->arch_strategy == p7_ARCH_FAST)
    {
      status = p7_Fastmodelmaker( msa, bld->symfrac, bld, ret_hmm, opt_tr);
      if      (status == eslENORESULT) ESL_XFAIL(status, bld->errbuf, "Alignment %s has no consensus columns w/ > %d%% residues - can't build a model.\n", msa->name != NULL ? msa->name : "", (int) (100 * bld->symfrac));
      else if (status == eslEMEM)      ESL_XFAIL(status, bld->errbuf, "Memory allocation failure in model construction.\n");
      else if (status != eslOK)        ESL_XFAIL(status, bld->errbuf, "internal error in model construction.\n");      
    }
  else if (bld->arch_strategy == p7_ARCH_HAND)
    {
      status = p7_Handmodelmaker( msa, bld, ret_hmm, opt_tr);
      if      (status == eslENORESULT) ESL_XFAIL(status, bld->errbuf, "Alignment %s has no annotated consensus columns - can't build a model.\n", msa->name != NULL ? msa->name : "");
      else if (status == eslEFORMAT)   ESL_XFAIL(status, bld->errbuf, "Alignment %s has no reference annotation line\n", msa->name != NULL ? msa->name : "");            
      else if (status == eslEMEM)      ESL_XFAIL(status, bld->errbuf, "Memory allocation failure in model construction.\n");
      else if (status != eslOK)        ESL_XFAIL(status, bld->errbuf, "internal error in model construction.\n");
    }
  return eslOK;

 ERROR:
  return status;
}



/* effective_seqnumber()
 *
 * <hmm> comes in with weighted observed counts. It goes out with
 * those observed counts rescaled to sum to the "effective sequence
 * number". 
 *
 * <msa> is needed because we may need to see the sequences in order 
 * to determine effective seq #. (for --eclust)
 *
 * <prior> is needed because we may need to parameterize test models
 * looking for the right relative entropy. (for --eent, the default)
 */
static int
effective_seqnumber(P7_BUILDER *bld, const ESL_MSA *msa, P7_HMM *hmm, const P7_BG *bg)
{
  int    status;
  int    i;

  if (bld->effn_strategy == p7_EFFN_ENTROPY_EXP) {
      double etarget; 
      double eff_nseq = 0.0;
      double exp;
      etarget = (bld->esigma - eslCONST_LOG2R * log( 2.0 / ((double) hmm->M * (double) (hmm->M+1)))) / (double) hmm->M; /* xref J5/36. */
      etarget = ESL_MAX(bld->re_target, etarget);

      status = p7_EntropyWeight_exp(hmm, bg, bld->prior, etarget, &exp);
      if      (status == eslEMEM) ESL_XFAIL(status, bld->errbuf, "memory allocation failed");
      else if (status != eslOK)   ESL_XFAIL(status, bld->errbuf, "internal failure in entropy weighting algorithm");

      p7_hmm_ScaleExponential(hmm, exp);

      for (i = 1; i <= hmm->M; i++)
        eff_nseq +=  esl_vec_FSum(hmm->mat[i], hmm->abc->K);


      eff_nseq /= hmm->M;
      hmm->eff_nseq = eff_nseq;

  } else {

    if      (bld->effn_strategy == p7_EFFN_NONE)    hmm->eff_nseq = msa->nseq;
    else if (bld->effn_strategy == p7_EFFN_SET)     hmm->eff_nseq = bld->eset;
    else if (bld->effn_strategy == p7_EFFN_CLUST)
    {
        int nclust;

        status = esl_msacluster_SingleLinkage(msa, bld->eid, NULL, NULL, &nclust);
        if      (status == eslEMEM) ESL_XFAIL(status, bld->errbuf, "memory allocation failed");
        else if (status != eslOK)   ESL_XFAIL(status, bld->errbuf, "single linkage clustering algorithm (at %d%% id) failed", (int)(100 * bld->eid));

        hmm->eff_nseq = (double) nclust;
    }
    else if (bld->effn_strategy == p7_EFFN_ENTROPY)
    {
        double etarget;
        double eff_nseq;
        etarget = (bld->esigma - eslCONST_LOG2R * log( 2.0 / ((double) hmm->M * (double) (hmm->M+1)))) / (double) hmm->M; /* xref J5/36. */
        etarget = ESL_MAX(bld->re_target, etarget);

        status = p7_EntropyWeight(hmm, bg, bld->prior, etarget, &eff_nseq);
        if      (status == eslEMEM) ESL_XFAIL(status, bld->errbuf, "memory allocation failed");
        else if (status != eslOK)   ESL_XFAIL(status, bld->errbuf, "internal failure in entropy weighting algorithm");
        hmm->eff_nseq = eff_nseq;
    }

    p7_hmm_Scale(hmm, hmm->eff_nseq / (double) hmm->nseq);

  }

  return eslOK;

 ERROR:
  return status;
}


/* parameterize()
 * Converts counts to probability parameters.
 */
static int
parameterize(P7_BUILDER *bld, P7_HMM *hmm)
{
  int status;

  if ((status = p7_ParameterEstimation(hmm, bld->prior)) != eslOK) ESL_XFAIL(status, bld->errbuf, "parameter estimation failed");

  return eslOK;

 ERROR:
  return status;
}



/* annotate()
 * Transfer annotation information from MSA to new HMM.
 * Also sets model-specific residue composition (hmm->compo).
 */
static int
annotate(P7_BUILDER *bld, const ESL_MSA *msa, P7_HMM *hmm)
{
  int status;

  /* Name. */
  if (msa->name) p7_hmm_SetName(hmm, msa->name);  
  else ESL_XFAIL(eslEINVAL, bld->errbuf, "Unable to name the HMM.");

  if ((status = p7_hmm_SetAccession  (hmm, msa->acc))           != eslOK) ESL_XFAIL(status, bld->errbuf, "Failed to record MSA accession");
  if ((status = p7_hmm_SetDescription(hmm, msa->desc))          != eslOK) ESL_XFAIL(status, bld->errbuf, "Failed to record MSA description");
  //  if ((status = p7_hmm_AppendComlog(hmm, go->argc, go->argv))   != eslOK) ESL_XFAIL(status, errbuf, "Failed to record command log");
  if ((status = p7_hmm_SetCtime(hmm))                           != eslOK) ESL_XFAIL(status, bld->errbuf, "Failed to record timestamp");
  if ((status = p7_hmm_SetComposition(hmm))                     != eslOK) ESL_XFAIL(status, bld->errbuf, "Failed to determine model composition");
  if ((status = p7_hmm_SetConsensus(hmm, NULL))                 != eslOK) ESL_XFAIL(status, bld->errbuf, "Failed to set consensus line");

  if (msa->cutset[eslMSA_GA1]) {
    hmm->cutoff[p7_GA1] = msa->cutoff[eslMSA_GA1];
    hmm->flags |= p7H_GA;
    if (msa->cutset[eslMSA_GA2])
      hmm->cutoff[p7_GA2] = msa->cutoff[eslMSA_GA2];
  }
  if (msa->cutset[eslMSA_TC1]) {
    hmm->cutoff[p7_TC1] = msa->cutoff[eslMSA_TC1];
    hmm->flags |= p7H_TC;
    if (msa->cutset[eslMSA_TC2])
      hmm->cutoff[p7_TC2] = msa->cutoff[eslMSA_TC2];
  }
  if (msa->cutset[eslMSA_NC1]) {
    hmm->cutoff[p7_NC1] = msa->cutoff[eslMSA_NC1];
    hmm->flags |= p7H_NC;
    if (msa->cutset[eslMSA_NC2])
      hmm->cutoff[p7_NC2] = msa->cutoff[eslMSA_NC2];
  }

  return eslOK;

 ERROR:
  return status;
}

/* calibrate()
 * 
 * Sets the E value parameters of the model with two short simulations.
 * A profile and an oprofile are created here. If caller wants to keep either
 * of them, it can pass non-<NULL> <opt_gm>, <opt_om> pointers.
 */
static int
evo_calibrate(P7_BUILDER *bld, P7_RATE *R, P7_HMM *hmm, P7_BG *bg, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om,
	      int noevo_msv, int noevo_vit, int noevo_fwd, float fixtime, float tol)
{
  int status;

  if (opt_gm != NULL) *opt_gm = NULL;
  if (opt_om != NULL) *opt_om = NULL;

  if ((status = p7_EvoCalibrate(R, hmm, bld, &(bld->r), &bg, opt_gm, opt_om, noevo_msv, noevo_vit, noevo_fwd, fixtime, tol)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  return status;
}


/* make_post_msa()
 * 
 * Optionally, we can return the alignment we actually built the model
 * from (including RF annotation on assigned consensus columns, and any
 * trace doctoring to enforce Plan7 consistency). 
 */
static int
make_post_msa(P7_BUILDER *bld, const ESL_MSA *premsa, const P7_HMM *hmm, P7_TRACE **tr, ESL_MSA **opt_postmsa)
{
  ESL_MSA  *postmsa  = NULL;
  int       optflags = p7_DEFAULT;
  int       status;

  if (opt_postmsa == NULL) return eslOK;

  /* someday we might want to transfer more info from HMM to postmsa */
  if ((status = p7_tracealign_MSA(premsa, tr, hmm->M, optflags, &postmsa)) != eslOK) goto ERROR;
  
  *opt_postmsa = postmsa;
  return eslOK;
  
 ERROR:
  if (postmsa != NULL) esl_msa_Destroy(postmsa);
  return status;
}
/*---------------- end, internal functions ----------------------*/

