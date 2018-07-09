/* Building profile HMMs from alignments.
 *
 * Contents: 
 *    x. h4_Build():      build new profile HMM from alignment
 *    x. H4_BUILD_CONFIG: customization of h4_Build()
 *    x. _experiment:     save counts files for training priors
 *    x. _experiment2:    compare old vs. new fragment marking
 * 
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_msa.h"
#include "esl_msaweight.h"
#include "esl_vectorops.h"


#include "h4_path.h"
#include "h4_profile.h"
#include "h4_build.h"


static int consensus_by_symfrac(const ESL_MSA *msa, float symfrac, const int8_t *fragassign, int8_t *matassign);
static int consensus_by_hand   (const ESL_MSA *msa, int8_t *matassign, char *errbuf);
static int stabilized_alen     (const ESL_MSA *msa, float occthresh, float *ret_alen);
static int mark_fragments      (const ESL_MSA *msa, float occthresh, float fragthresh, int8_t *fragassign);
static int collect_counts      (const ESL_MSA *msa, const int8_t *fragassign, const int8_t *matassign, H4_PROFILE *hmm);


int
h4_Build(H4_BUILD_CONFIG *cfg, ESL_MSA *msa, H4_PROFILE **ret_hmm, char *errbuf)
{
  H4_PROFILE *hmm        = NULL;
  int8_t     *fragassign = NULL;
  int8_t     *matassign  = NULL;
  int         M          = 0;
  int         apos;
  int         status;

  ESL_DASSERT1(( msa->flags && eslMSA_DIGITAL ));

  ESL_ALLOC(fragassign, sizeof(int8_t) * msa->nseq);
  ESL_ALLOC(matassign,  sizeof(int8_t) * (msa->alen + 1));
  matassign[0] = 0; // unused; convention.

  // validate_msa()
  // esl_msa_Checksum()

  /* x. Set relative weights, in msa->wgt.
   */
  status = eslOK;
  if      (! cfg || cfg->wgt_strategy == h4_WGT_PB)     status = esl_msaweight_PB(msa);
  else if (         cfg->wgt_strategy == h4_WGT_GIVEN)  status = eslOK;
  else if (         cfg->wgt_strategy == h4_WGT_GSC)    status = esl_msaweight_GSC(msa);
  else if (         cfg->wgt_strategy == h4_WGT_BLOSUM) status = esl_msaweight_BLOSUM(msa, cfg->wid);
  else if (         cfg->wgt_strategy == h4_WGT_NONE)   esl_vec_DSet(msa->wgt, msa->nseq, 1.);
  else    ESL_EXCEPTION(eslEINCONCEIVABLE, "no such weighting strategy");
  if (status != eslOK) goto ERROR;

  /* x. Define which sequences are considered to be fragments (local alignments).
   */
  status = mark_fragments(msa,
			  cfg ? cfg->alen_occthresh : h4_DEFAULT_ALEN_OCCTHRESH,
			  cfg ? cfg->fragthresh     : h4_DEFAULT_FRAGTHRESH,
			  fragassign);
  if (status != eslOK) goto ERROR;

  /* x. Define which columns are considered to be consensus.
   */
  if      (! cfg || cfg->arch_strategy == h4_ARCH_SYMFRAC)  status = consensus_by_symfrac(msa, cfg ? cfg->symfrac : h4_DEFAULT_SYMFRAC, fragassign, matassign);
  else if (         cfg->arch_strategy == h4_ARCH_GIVEN)    status = consensus_by_hand(msa, matassign, errbuf);
  else    ESL_EXCEPTION(eslEINCONCEIVABLE, "no such profile architecture strategy");
  if (status != eslOK) goto ERROR;

  /* Allocate the new profile.
   */
  for (apos = 1; apos <= msa->alen; apos++) if (matassign[apos]) M++;
  hmm = h4_profile_Create(msa->abc, M);

  /* x. Collect observed (relative-weighted) counts from alignment in hmm->t[] and ->e[]
   */
  if ((status = collect_counts(msa, fragassign, matassign, hmm)) != eslOK) goto ERROR;

  // mask_columns();

  // effective_seqnumber();
  // parameterize();
  // annotate();
  // calibrate();
  // make_post_msa()

  // *ret_hmm = hmm;
  free(fragassign);
  free(matassign);
  return eslOK;

 ERROR:
  free(fragassign);
  free(matassign);
  return status;
}




/* consensus_by_symfrac()
 * Sets <matassign[1..alen]> to 1/0 flags, defining consensus columns.
 * 
 * Args:    msa        : multiple sequence alignment
 *          symfrac    : define col as consensus if weighted residue fraction >= symfrac
 *          fragassign : [0..nseq-1] 1/0 flags marking fragments (local alignments)
 *          matassign  : RETURN: [1..alen] 1/0 flags marking consensus columns
 * 
 * Returns: <eslOK> on success;
 *          and matassign[(0)1..alen] has 1/0 for consensus vs. non.
 *          
 * Throws:  <eslEMEM> on allocation error.
 *          Now state of <matassign> is undefined.         
 */
static int
consensus_by_symfrac(const ESL_MSA *msa, float symfrac, const int8_t *fragassign, int8_t *matassign)
{
  float *r      = NULL;  // weighted residue count for each column 1..alen
  float *totwgt = NULL;  // weighted residue+gap count for each column 1..alen (not constant across cols, because of fragments)
  int    apos, idx;
  int    lpos, rpos;
  int    status;

  ESL_ALLOC(r,      sizeof(float) * (msa->alen+1));
  ESL_ALLOC(totwgt, sizeof(float) * (msa->alen+1));
  esl_vec_FSet(r,      msa->alen+1, 0.);
  esl_vec_FSet(totwgt, msa->alen+1, 0.);

  for (idx = 0; idx < msa->nseq; idx++)
    {
      lpos = 1;
      rpos = msa->alen;
      if (fragassign[idx])
	{
	  for (lpos = 1;         lpos <= msa->alen; lpos++) if (! esl_abc_XIsGap(msa->abc, msa->ax[idx][lpos])) break;  // find first residue of this seq
	  for (rpos = msa->alen; rpos >= 1;         rpos--) if (! esl_abc_XIsGap(msa->abc, msa->ax[idx][rpos])) break;  //  ... and last.
	  // if sequence is empty: lpos = msa->alen+1, rpos = 0, and the loop body below doesn't execute.
	}
      for (apos = lpos; apos <= rpos; apos++)
	{
	  if      (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) { r[apos] += msa->wgt[idx]; totwgt[apos] += msa->wgt[idx]; }
	  else if (esl_abc_XIsGap(msa->abc,     msa->ax[idx][apos])) {                           totwgt[apos] += msa->wgt[idx];  }
	  // missing ~, nonresidue * don't count either way toward the calculation.
	}
    }

  for (apos = 1; apos <= msa->alen; apos++)
    matassign[apos] = (totwgt[apos] > 0. && r[apos] / totwgt[apos] >= symfrac) ? 1 : 0;

  free(r);
  free(totwgt);
  return eslOK;

 ERROR:
  free(r);
  free(totwgt);
  return status;
}


/* consensus_by_hand()
 * Define consensus columns using provided alignment annotation (#=GC RF or seq_cons)
 */
static int
consensus_by_hand(const ESL_MSA *msa, int8_t *matassign, char *errbuf)
{
  int apos;
  
  if (! msa->rf)
    ESL_FAIL(eslEFORMAT, errbuf, "no consensus column (#=GC RF, #=GC seq_cons) annotation on MSA");
  for (apos = 1; apos <= msa->alen; apos++)
    matassign[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos-1])? 0 : 1);  // watch off-by-one. rf is 0..alen-1, matassign is 1..alen
  return eslOK;
}



/* stabilized_alen()
 * Calculates an nseq-stabilized alignment length for mark_fragments()
 * SRE, Fri 06 Jul 2018 [JB1685 BOS-PIT]
 * 
 * Calculates an alignment length for <mark_fragments()> to use.
 * 
 * The issue here is that the length of a multiple sequence alignment
 * (in columns) scales with the number of sequences in it, as
 * insertions in individual sequences accumulate. A <mark_fragments()>
 * rule based on rlen/alen can end up calling all sequences fragments,
 * for deep sequence alignments.
 * 
 * To mitigate this, we calculate the (weighted) fraction f_c of
 * residues in each column c. If f_c >= occthresh, the column is
 * counted toward the "stabilized" alignment length.  A reasonable
 * value of <occthresh> might be around 0.01, meaning the alignment
 * length if there were ~100 sequences.
 */
static int
stabilized_alen(const ESL_MSA *msa, float occthresh, float *ret_alen)
{
  float *colwgt = NULL;
  float  alen   = 0;
  int    idx,col;
  int    status;

  ESL_ALLOC(colwgt, sizeof(float) * (msa->alen+1));
  esl_vec_FSet(colwgt, msa->alen+1, 0.);
  for (idx = 0; idx < msa->nseq; idx++)
    for (col = 1; col <= msa->alen; col++)
      if (! esl_abc_XIsGap(msa->abc, msa->ax[idx][col])) // * and ~ count towards occupancy here
	colwgt[col] += msa->wgt[idx];
  esl_vec_FScale(colwgt+1, msa->alen, 1./(float) msa->nseq);

  for (col = 1; col <= msa->alen; col++)
    if (colwgt[col] >= occthresh) alen += 1.0;

  *ret_alen = alen;
  free(colwgt);
  return eslOK;

 ERROR:
  free (colwgt);
  return status;
}


/* mark_fragments()
 * Set fragassign[i] TRUE | FALSE to mark local alignment frags
 * SRE, Fri 06 Jul 2018 [JB1685 BOS-PIT]
 * 
 * Heuristically define sequence fragments (as opposed to "full
 * length" sequences) in <msa>. Set <fragassign[i]> to TRUE if seq <i>
 * is a fragment, else FALSE.
 * 
 * Args:    msa        : msa with msa->nseq seqs
 *          occthresh  : parameter passed to stabilized_alen()
 *          fragthresh : if rlen/alen <= this, seq is a fragment
 *          fragassign : RESULT: fragassign[i]=1|0 if seq i is|isn't a fragment.
 *                       Caller provides allocation for <msa->nseq> flags.
 * 
 * Let r = raw sequence length (in residues)
 *     a = stabilized alignment length (in columns)
 *     t = <fragthresh> parameter
 * If r/a <= t, define seq as a fragment.    
 * 
 * See <stabilized_alen()> for explanation of what "stabilization"
 * means.
 * 
 * Why not a rule based on raw unaligned length alone?
 * An example of an edge case this needs to deal with (see
 * the i007-hmmbuild-fragments test):
 *   seq1 ACDEFGHIKL------------------------------
 *   seq2 ----------MNPQRSTVWY--------------------
 *   seq3 --------------------ACDEFGHIKL----------
 *   seq4 ------------------------------MNPQRSTVWY
 * All four sequences are fragments of equal length 10, in an ungapped
 * sequence alignment of length 40. This sort of thing arises on
 * metagenomic sequence sample alignments.
 *
 * Adapted from (and supersedes) esl_msa_MarkFragments(), which does
 * not use stabilization, and which marked fragments by introducing
 * leading/trailing ~ characters into the alignment itself.
 */
static int
mark_fragments(const ESL_MSA *msa, float occthresh, float fragthresh, int8_t *fragassign)
{
  int   idx;
  float alen;

  stabilized_alen(msa, occthresh, &alen);
  for (idx = 0; idx < msa->nseq; idx++)
    fragassign[idx] = 
      ((float) esl_abc_dsqrlen(msa->abc, msa->ax[idx]) / alen <= fragthresh) ?  
      TRUE : FALSE;
  return eslOK;
}






static int
collect_counts(const ESL_MSA *msa, const int8_t *fragassign, const int8_t *matassign, H4_PROFILE *hmm)
{
  H4_PATH *pi = h4_path_Create();
  int      idx;
  int      status;
  
  for (idx = 0; idx < msa->nseq; idx++)
    {
      if (fragassign[idx]) { if ((status = h4_path_InferLocal (msa->abc, msa->ax[idx], msa->alen, matassign, pi)) != eslOK) goto ERROR; }
      else                 { if ((status = h4_path_InferGlocal(msa->abc, msa->ax[idx], msa->alen, matassign, pi)) != eslOK) goto ERROR; }

      if ((status = h4_path_Count(pi, msa->ax[idx], msa->wgt[idx], hmm)) != eslOK) goto ERROR;

      h4_path_Reuse(pi);
    }

  h4_path_Destroy(pi);
  return eslOK;

 ERROR:
  h4_path_Destroy(pi);
  return status;
}


#if 0
static int
annotate(H4_BUILD_CONFIG *cfg, ESL_MSA *msa, H4_PROFILE *hmm, char *errbuf, int8_t *matassign)
{
  /* transfer from MSA to new HMM: */
  /* msa->rf   */
  /* msa->mm   */
  /* msa->ss_cons */
  /* msa->sa_cons */
  /* map */
  

  return eslOK;
}
#endif


/*****************************************************************
 * x. H4_BUILD_CONFIG
 *****************************************************************/

H4_BUILD_CONFIG *
h4_build_config_Create(void)
{
  H4_BUILD_CONFIG *cfg = NULL;
  int              status;

  ESL_ALLOC(cfg, sizeof(H4_BUILD_CONFIG));

  cfg->arch_strategy  = h4_ARCH_SYMFRAC;
  cfg->wgt_strategy   = h4_WGT_PB;
  cfg->symfrac        = h4_DEFAULT_SYMFRAC;
  cfg->alen_occthresh = h4_DEFAULT_ALEN_OCCTHRESH;
  cfg->fragthresh     = h4_DEFAULT_FRAGTHRESH;
  cfg->wid            = h4_DEFAULT_WID;
  return cfg;

 ERROR:
  h4_build_config_Destroy(cfg);
  return NULL;
}

void
h4_build_config_Destroy(H4_BUILD_CONFIG *cfg)
{
  free(cfg);
}
/*----------------- end, H4_BUILD_CONFIG ------------------------*/




/*****************************************************************
 * x. _experiment: save counts files for training priors
 *****************************************************************/
#ifdef h4BUILD_EXPERIMENT

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "utility for saving counts files for training priors";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);




  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /* h4BUILD_EXPERIMENT */



/*****************************************************************
 * x. _experiment2: compare old vs. new fragment marking 
 *****************************************************************/

#ifdef h4BUILD_EXPERIMENT2
#include "h4_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"




static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",        0 },
  { "--dna",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use DNA alphabet",                            0 },
  { "--rna",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use RNA alphabet",                            0 },
  { "--amino",   eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use protein alphabet",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msafile>";
static char banner[] = "testing old v. new fragment-marking strategy";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char        *msafile = esl_opt_GetArg(go, 1);
  int          infmt   = eslMSAFILE_UNKNOWN;
  ESL_ALPHABET *abc    = NULL;
  ESL_MSAFILE  *afp    = NULL;
  ESL_MSA      *msa    = NULL;
  int           nali   = 0;
  float         occthresh = 0.01;
  float         fragthresh = 0.5;
  int           nold, nnew;
  float         s_alen;
  int           idx;
  int           status;
  
  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO); 

  status = esl_msafile_Open(&abc, msafile, NULL, infmt, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  esl_dataheader(stdout,
		 20, "name",  10,  "nseq",  10, "alen",     10, "s_alen", 
		 10, "n_old", 10, "n_new",  10, "frac_old", 10, "frac_new", 0);

  while ((status = esl_msafile_Read(afp, &msa)) == eslOK)
    {
      int8_t *old_fragassign = NULL;
      int8_t *new_fragassign = NULL;

      ESL_ALLOC(old_fragassign, sizeof(int8_t) * msa->nseq);
      ESL_ALLOC(new_fragassign, sizeof(int8_t) * msa->nseq);
      
      /* New calculation with stabilized alen */
      if ((status = mark_fragments(msa, occthresh, fragthresh, new_fragassign)) != eslOK) goto ERROR;
      if ((status = stabilized_alen(msa, occthresh, &s_alen))                   != eslOK) goto ERROR;

      /* Reproduce the H3 calculation, while setting flag instead of marking ~ in the msa */
      for (idx = 0; idx < msa->nseq; idx++)
	old_fragassign[idx] = 
	  ((float) esl_abc_dsqrlen(msa->abc, msa->ax[idx]) / (float) msa->alen <= fragthresh) ?  
	  TRUE : FALSE;

      for (idx = 0, nold = 0, nnew = 0; idx < msa->nseq; idx++)
	{
	  if (new_fragassign[idx]) nnew++;
	  if (old_fragassign[idx]) nold++;
	}
      
      printf("%20s %10d %10d %10.1f %10d %10d %10.4f %10.4f\n",
	     msa->name, msa->nseq, (int) msa->alen, s_alen, nold, nnew,
	     (float) nold / (float) msa->nseq,
	     (float) nnew / (float) msa->nseq);
      
      nali++;
      esl_msa_Destroy(msa);  msa            = NULL;
      free(old_fragassign);  old_fragassign = NULL;
      free(new_fragassign);  new_fragassign = NULL;
    }
  
  if (nali == 0 || status != eslEOF)
    esl_msafile_ReadFailure(afp, status);


  esl_msafile_Close(afp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;

 ERROR:
  return status;
}


#endif /* h4BUILD_EXPERIMENT2 */

/*--------------- end, experiment driver ------------------------*/
