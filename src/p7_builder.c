/* Standardized pipeline for construction of new HMMs.
 * 
 * Contents:
 *    1. P7_BUILDER: allocation, initialization, destruction
 *    2. Standardized model construction API.
 *    3. Internal functions.
 *    4. Copyright and license information
 *    
 * SRE, Thu Dec 11 08:44:58 2008 [Janelia] [Requiem for a Dream]
 * SVN $Id$
 */   
#include "p7_config.h"

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

/*****************************************************************
 * 1. P7_BUILDER: allocation, initialization, destruction
 *****************************************************************/

/* Function:  p7_builder_Create()
 * Synopsis:  Create a default HMM construction configuration.
 * Incept:    SRE, Thu Dec 11 13:14:21 2008 [Janelia]
 *
 * Purpose:   Create a construction configuration for building
 *            HMMs in alphabet <abc>, and return a pointer to it.
 *            
 *            An application configuration <go> may optionally be
 *            provided. If <go> is <NULL>, default parameters are
 *            used. If <go> is non-<NULL>, it must include appropriate
 *            settings for all 24 ``standard build options'':
 *            
 *            Model construction:   --fast --hand --symfrac --fragthresh
 *            Relative weighting:   --wgsc --wblosum --wpb --wgiven --wid
 *            Effective seq #:      --eent --eclust --enone --eset --ere --esigma --eid
 *            E-val calibration:    --EmL --EmN --EvL --EvN --EfL --EfN --Eft
 *            run-to-run variation: --seed
 *            
 *            See <hmmbuild.c> or other big users of the build
 *            pipeline for an example of appropriate <ESL_GETOPTS>
 *            initializations of these 24 options.
 */
P7_BUILDER *
p7_builder_Create(const ESL_GETOPTS *go, const ESL_ALPHABET *abc)
{
  P7_BUILDER *bld = NULL;
  int         seed;
  int         status;


  ESL_ALLOC(bld, sizeof(P7_BUILDER));
  bld->prior        = NULL;
  bld->r            = NULL;
  bld->S            = NULL;
  bld->Q            = NULL;
  bld->eset         = -1.0;	/* -1.0 = unset; must be set if effn_strategy is p7_EFFN_SET */
  bld->re_target    = -1.0;

  if (go == NULL) 
    {
      bld->arch_strategy = p7_ARCH_FAST;
      bld->wgt_strategy  = p7_WGT_PB;
      bld->effn_strategy = p7_EFFN_ENTROPY;
      seed               = 0;
    }
  else 
    {
      if      (esl_opt_GetBoolean(go, "--fast"))    bld->arch_strategy = p7_ARCH_FAST;
      else if (esl_opt_GetBoolean(go, "--hand"))    bld->arch_strategy = p7_ARCH_HAND;

      

      if      (esl_opt_GetBoolean(go, "--wpb"))     bld->wgt_strategy = p7_WGT_PB;
      else if (esl_opt_GetBoolean(go, "--wgsc"))    bld->wgt_strategy = p7_WGT_GSC;
      else if (esl_opt_GetBoolean(go, "--wblosum")) bld->wgt_strategy = p7_WGT_BLOSUM;
      else if (esl_opt_GetBoolean(go, "--wnone"))   bld->wgt_strategy = p7_WGT_NONE;
      else if (esl_opt_GetBoolean(go, "--wgiven"))  bld->wgt_strategy = p7_WGT_GIVEN;

      if      (esl_opt_GetBoolean(go, "--eent"))    bld->effn_strategy = p7_EFFN_ENTROPY;
      else if (esl_opt_GetBoolean(go, "--eclust"))  bld->effn_strategy = p7_EFFN_CLUST;
      else if (esl_opt_GetBoolean(go, "--enone"))   bld->effn_strategy = p7_EFFN_NONE;
      else if (esl_opt_IsOn      (go, "--eset"))  { bld->effn_strategy = p7_EFFN_SET;      bld->eset = esl_opt_GetReal(go, "--eset"); }

      seed = esl_opt_GetInteger(go, "--seed");
    }

  /* The default RE target is alphabet dependent. */
  if (go != NULL &&  esl_opt_IsOn (go, "--ere")) 
    bld->re_target = esl_opt_GetReal(go, "--ere");
  else {
    switch (abc->type) {
    case eslAMINO:  bld->re_target = p7_ETARGET_AMINO; break;
    case eslDNA:    bld->re_target = p7_ETARGET_DNA;   break;
    case eslRNA:    bld->re_target = p7_ETARGET_DNA;   break;
    default:        bld->re_target = p7_ETARGET_OTHER; break;
    }
  }

  bld->symfrac    = (go != NULL) ?  esl_opt_GetReal   (go, "--symfrac")    : 0.5; 
  bld->fragthresh = (go != NULL) ?  esl_opt_GetReal   (go, "--fragthresh") : 0.5; 
  bld->wid        = (go != NULL) ?  esl_opt_GetReal   (go, "--wid")        : 0.62;
  bld->esigma     = (go != NULL) ?  esl_opt_GetReal   (go, "--esigma")     : 45.0;
  bld->eid        = (go != NULL) ?  esl_opt_GetReal   (go, "--eid")        : 0.62;
  bld->EmL        = (go != NULL) ?  esl_opt_GetInteger(go, "--EmL")        : 200;
  bld->EmN        = (go != NULL) ?  esl_opt_GetInteger(go, "--EmN")        : 200;
  bld->EvL        = (go != NULL) ?  esl_opt_GetInteger(go, "--EvL")        : 200;
  bld->EvN        = (go != NULL) ?  esl_opt_GetInteger(go, "--EvN")        : 200;
  bld->EfL        = (go != NULL) ?  esl_opt_GetInteger(go, "--EfL")        : 100;
  bld->EfN        = (go != NULL) ?  esl_opt_GetInteger(go, "--EfN")        : 200;
  bld->Eft        = (go != NULL) ?  esl_opt_GetReal   (go, "--Eft")        : 0.04;
    
  /* Normally we reinitialize the RNG to original seed before calibrating each model.
   * This eliminates run-to-run variation.
   * As a special case, seed==0 means choose an arbitrary seed and shut off the
   * reinitialization; this allows run-to-run variation.
   */
  bld->r            = esl_randomness_CreateFast(seed);
  bld->do_reseeding = (seed == 0) ? FALSE : TRUE;

  switch (abc->type) {
  case eslAMINO: bld->prior = p7_prior_CreateAmino();      break;
  case eslDNA:   bld->prior = p7_prior_CreateNucleic();    break;
  case eslRNA:   bld->prior = p7_prior_CreateNucleic();    break;
  default:       bld->prior = p7_prior_CreateLaplace(abc); break;
  }
  if (bld->prior == NULL) goto ERROR;

  bld->abc       = abc;
  bld->errbuf[0] = '\0';
  return bld;
  
 ERROR:
  p7_builder_Destroy(bld);
  return NULL;
}




/* Function:  p7_builder_SetScoreSystem()
 * Synopsis:  Initialize score system for single sequence queries.
 * Incept:    SRE, Fri Dec 12 10:04:36 2008 [Janelia]
 *
 * Purpose:   Initialize the builder <bld> to be able to parameterize
 *            single sequence queries.
 *            
 *            Read a standard substitution score matrix from file
 *            <mxfile>. If <mxfile> is <NULL>, default to BLOSUM62
 *            scores. If <mxfile> is "-", read score matrix from
 *            <stdin> stream. If <env> is non-<NULL> and <mxfile> is
 *            not found in the current working directory, look for
 *            <mxfile> in colon-delimited directory list contained in
 *            environment variable <env>.
 *            
 *            Set the gap-open and gap-extend probabilities to
 *            <popen>, <pextend>, respectively.
 *
 *
 * Args:      bld      - <P7_BUILDER> to initialize
 *            mxfile   - score matrix file to use, or NULL for BLOSUM62 default
 *            env      - env variable containing directory list where <mxfile> may reside
 *            popen    - gap open probability
 *            pextend  - gap extend probability
 *
 * Returns:   <eslOK> on success.
 *            
 *            <eslENOTFOUND> if <mxfile> can't be found or opened, even
 *            in any of the directories specified by the <env> variable.   
 *            
 *            <eslEINVAL> if the score matrix can't be converted into
 *            conditional probabilities by the Yu and Altschul method,
 *            either because it isn't a symmetric matrix or because
 *            the Yu/Altschul numerical method fails to converge. 
 * 
 *            On either error, <bld->errbuf> contains a useful error message
 *            for the user.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_builder_SetScoreSystem(P7_BUILDER *bld, const char *mxfile, const char *env, double popen, double pextend)
{
  ESL_FILEPARSER  *efp      = NULL;
  double          *fa       = NULL;
  double          *fb       = NULL;
  double           slambda;
  int              a,b;
  int              status;


  bld->errbuf[0] = '\0';

  /* If a score system is already set, delete it. */
  if (bld->S != NULL) esl_scorematrix_Destroy(bld->S);
  if (bld->Q != NULL) esl_dmatrix_Destroy(bld->Q);

  /* Get the scoring matrix */
  if ((bld->S  = esl_scorematrix_Create(bld->abc)) == NULL) { status = eslEMEM; goto ERROR; }
  if (mxfile == NULL) 
    {
      if ((status = esl_scorematrix_SetBLOSUM62(bld->S)) != eslOK) goto ERROR;
    } 
  else 
    {
      if ((status = esl_fileparser_Open(mxfile, env, &efp)) != eslOK) ESL_XFAIL(status, bld->errbuf, "Failed to find or open matrix file %s", mxfile);
      if ((status = esl_sco_Read(efp, bld->abc, &(bld->S))) != eslOK) ESL_XFAIL(status, bld->errbuf, "Failed to read matrix from %s:\n%s", mxfile, efp->errbuf);
      esl_fileparser_Close(efp); efp = NULL;
    }
  if (! esl_scorematrix_IsSymmetric(bld->S)) 
    ESL_XFAIL(eslEINVAL, bld->errbuf, "Matrix isn't symmetric");
  if ((status = esl_sco_Probify(bld->S, &(bld->Q), &fa, &fb, &slambda)) != eslOK) 
    ESL_XFAIL(eslEINVAL, bld->errbuf, "Yu/Altschul method failed to backcalculate probabilistic basis of score matrix");

  for (a = 0; a < bld->abc->K; a++)
    for (b = 0; b < bld->abc->K; b++)
      bld->Q->mx[a][b] /= fa[a];	/* Q->mx[a][b] is now P(b | a) */

  bld->popen   = popen;
  bld->pextend = pextend;

  free(fa);
  free(fb);
  return eslOK;

 ERROR:
  if (efp != NULL) esl_fileparser_Close(efp);
  if (fa  != NULL) free(fa);
  if (fb  != NULL) free(fb);
  return status;
}


/* Function:  p7_builder_Destroy()
 * Synopsis:  Free a <P7_BUILDER>
 * Incept:    SRE, Thu Dec 11 13:15:45 2008 [Janelia]
 *
 * Purpose:   Frees a <P7_BUILDER> object.
 */
void
p7_builder_Destroy(P7_BUILDER *bld)
{
  if (bld == NULL) return;

  if (bld->prior   != NULL) p7_prior_Destroy(bld->prior);
  if (bld->r       != NULL) esl_randomness_Destroy(bld->r);
  if (bld->Q       != NULL) esl_dmatrix_Destroy(bld->Q);
  if (bld->S       != NULL) esl_scorematrix_Destroy(bld->S);

  free(bld);
  return;
}
/*------------------- end, P7_BUILDER ---------------------------*/




/*****************************************************************
 * 2. Standardized model construction API.
 *****************************************************************/

static int    validate_msa         (P7_BUILDER *bld, ESL_MSA *msa);
static int    relative_weights     (P7_BUILDER *bld, ESL_MSA *msa);
static int    build_model          (P7_BUILDER *bld, ESL_MSA *msa, P7_HMM **ret_hmm, P7_TRACE ***opt_tr);
static int    effective_seqnumber  (P7_BUILDER *bld, const ESL_MSA *msa, P7_HMM *hmm, const P7_BG *bg);
static int    parameterize         (P7_BUILDER *bld, P7_HMM *hmm);
static int    annotate             (P7_BUILDER *bld, const ESL_MSA *msa, P7_HMM *hmm);
static int    calibrate            (P7_BUILDER *bld, P7_HMM *hmm, P7_BG *bg, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om);
static int    make_post_msa        (P7_BUILDER *bld, const ESL_MSA *premsa, const P7_HMM *hmm, P7_TRACE **tr, ESL_MSA **opt_postmsa);

/* Function:  p7_Builder()
 * Synopsis:  Build a new HMM from an MSA.
 * Incept:    SRE, Thu Dec 11 13:24:38 2008 [Janelia]
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
 *            opt_postmsa - optRETURN: RF-annotated, possibly modified MSA 
 *            opt_gm      - optRETURN: profile corresponding to <hmm>
 *            opt_om      - optRETURN: optimized profile corresponding to <gm>
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
p7_Builder(P7_BUILDER *bld, ESL_MSA *msa, P7_BG *bg,
	   P7_HMM **opt_hmm, P7_TRACE ***opt_trarr, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om,
	   ESL_MSA **opt_postmsa)
{
  P7_HMM     *hmm    = NULL;
  P7_TRACE  **tr     = NULL;
  P7_TRACE ***tr_ptr = (opt_trarr != NULL || opt_postmsa != NULL) ? &tr : NULL;
  int         status;

  if ((status =  validate_msa       (bld, msa))                       != eslOK) goto ERROR;
  if ((status =  relative_weights   (bld, msa))                       != eslOK) goto ERROR;
  if ((status =  esl_msa_MarkFragments(msa, bld->fragthresh))         != eslOK) goto ERROR;
  if ((status =  build_model        (bld, msa, &hmm, tr_ptr))         != eslOK) goto ERROR;
  if ((status =  effective_seqnumber(bld, msa, hmm, bg))              != eslOK) goto ERROR;
  if ((status =  parameterize       (bld, hmm))                       != eslOK) goto ERROR;
  if ((status =  annotate           (bld, msa, hmm))                  != eslOK) goto ERROR;
  if ((status =  calibrate          (bld, hmm, bg, opt_gm, opt_om))   != eslOK) goto ERROR;
  if ((status =  make_post_msa      (bld, msa, hmm, tr, opt_postmsa)) != eslOK) goto ERROR;

  if (opt_hmm   != NULL) *opt_hmm   = hmm; else p7_hmm_Destroy(hmm);
  if (opt_trarr != NULL) *opt_trarr = tr;  else p7_trace_DestroyArray(tr, msa->nseq);
  return eslOK;

 ERROR:
  p7_hmm_Destroy(hmm);
  p7_trace_DestroyArray(tr, msa->nseq);
  if (opt_gm    != NULL) p7_profile_Destroy(*opt_gm);
  if (opt_om    != NULL) p7_oprofile_Destroy(*opt_om);
  return status;
}


/* Function:  p7_SingleBuilder()
 * Synopsis:  Build a new HMM from a single sequence.
 * Incept:    SRE, Fri Dec 12 10:52:45 2008 [Janelia]
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
 *            bg        - null model (needed to paramaterize insert emission probs)
 *            opt_hmm   - optRETURN: new HMM
 *            opt_gm    - optRETURN: profile corresponding to <hmm>
 *            opt_om    - optRETURN: optimized profile corresponding to <gm>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if <bld> isn't properly configured somehow.
 */
int
p7_SingleBuilder(P7_BUILDER *bld, ESL_SQ *sq, P7_BG *bg, P7_HMM **opt_hmm,
		 P7_TRACE **opt_tr, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om)
{
  P7_HMM   *hmm = NULL;
  P7_TRACE *tr  = NULL;
  int       k;
  int       status;
  
  bld->errbuf[0] = '\0';
  if (! bld->Q) ESL_XEXCEPTION(eslEINVAL, "score system not initialized");

  if ((status = p7_Seqmodel(bld->abc, sq->dsq, sq->n, sq->name, bld->Q, bg->f, bld->popen, bld->pextend, &hmm)) != eslOK) goto ERROR;
  if ((status = calibrate(bld, hmm, bg, opt_gm, opt_om))                                                        != eslOK) goto ERROR;

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

  if (opt_hmm   != NULL) *opt_hmm = hmm; else p7_hmm_Destroy(hmm);
  if (opt_tr    != NULL) *opt_tr  = tr;
  return eslOK;

 ERROR:
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
 * esl_msa_MarkFragments()).
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
  int status;

  if      (bld->wgt_strategy == p7_WGT_NONE)                    { esl_vec_DSet(msa->wgt, msa->nseq, 1.); status = eslOK; }
  else if (bld->wgt_strategy == p7_WGT_GIVEN)                   status = eslOK;
  else if (bld->wgt_strategy == p7_WGT_PB)                      status = esl_msaweight_PB(msa); 
  else if (bld->wgt_strategy == p7_WGT_GSC)                     status = esl_msaweight_GSC(msa); 
  else if (bld->wgt_strategy == p7_WGT_BLOSUM)                  status = esl_msaweight_BLOSUM(msa, bld->wid); 

  if (status != eslOK) ESL_FAIL(status, bld->errbuf, "failed to set relative weights in alignment");
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
      status = p7_Fastmodelmaker(msa, bld->symfrac, ret_hmm, opt_tr);
      if      (status == eslENORESULT) ESL_XFAIL(status, bld->errbuf, "Alignment %s has no consensus columns w/ > %d%% residues - can't build a model.\n", msa->name != NULL ? msa->name : "", (int) (100 * bld->symfrac));
      else if (status == eslEMEM)      ESL_XFAIL(status, bld->errbuf, "Memory allocation failure in model construction.\n");
      else if (status != eslOK)        ESL_XFAIL(status, bld->errbuf, "internal error in model construction.\n");      
    }
  else if (bld->arch_strategy == p7_ARCH_HAND)
    {
      status = p7_Handmodelmaker(msa, ret_hmm, opt_tr);
      if      (status == eslENORESULT) ESL_XFAIL(status, bld->errbuf, "Alignment %s has no annotated consensus columns - can't build a model.\n", msa->name != NULL ? msa->name : "");
      else if (status == eslEFORMAT)   ESL_XFAIL(status, bld->errbuf, "Alignment %s has no reference annotation line\n", msa->name != NULL ? msa->name : "");            
      else if (status == eslEMEM)      ESL_XFAIL(status, bld->errbuf, "Memory allocation failure in model construction.\n");
      else if (status != eslOK)        ESL_XFAIL(status, bld->errbuf, "internal error in model construction.\n");
    }
  return eslOK;

 ERROR:
  return status;
}


/* set_effective_seqnumber()
 * Incept:    SRE, Fri May 11 08:14:57 2007 [Janelia]
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
  if ((status = esl_msa_Checksum(msa, &(hmm->checksum)))        != eslOK) ESL_XFAIL(status, bld->errbuf, "Failed to record checksum"); 
  if ((status = p7_hmm_SetComposition(hmm))                     != eslOK) ESL_XFAIL(status, bld->errbuf, "Failed to determine model composition");
  hmm->flags |= p7H_CHKSUM;
  hmm->flags |= p7H_COMPO;

  if (msa->cutset[eslMSA_GA1] && msa->cutset[eslMSA_GA2]) { hmm->cutoff[p7_GA1] = msa->cutoff[eslMSA_GA1]; hmm->cutoff[p7_GA2] = msa->cutoff[eslMSA_GA2]; hmm->flags |= p7H_GA; }
  if (msa->cutset[eslMSA_TC1] && msa->cutset[eslMSA_TC2]) { hmm->cutoff[p7_TC1] = msa->cutoff[eslMSA_TC1]; hmm->cutoff[p7_TC2] = msa->cutoff[eslMSA_TC2]; hmm->flags |= p7H_TC; }
  if (msa->cutset[eslMSA_NC1] && msa->cutset[eslMSA_NC2]) { hmm->cutoff[p7_NC1] = msa->cutoff[eslMSA_NC1]; hmm->cutoff[p7_NC2] = msa->cutoff[eslMSA_NC2]; hmm->flags |= p7H_NC; }

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
calibrate(P7_BUILDER *bld, P7_HMM *hmm, P7_BG *bg, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om)
{
  int status;

  if (opt_gm != NULL) *opt_gm = NULL;
  if (opt_om != NULL) *opt_om = NULL;

  if ((status = p7_Calibrate(hmm, bld, &(bld->r), &bg, opt_gm, opt_om)) != eslOK) goto ERROR;
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





/*****************************************************************
 * @LICENSE@
 *****************************************************************/
