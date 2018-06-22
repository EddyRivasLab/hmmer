

int
h4_Build(H4_BUILD_CONFIG *cfg, ESL_MSA *msa, H4_PROFILE **ret_hmm, char *errbuf)
{

  ESL_DASSERT1(( msa->flags && eslMSA_DIGITAL ));

  relative_weights();
  define_consensus();
  mask_columns();
  collect_counts();
  effective_seqnumber();
  parameterize();
  annotate();
  calibrate();

  *ret_hmm = hmm;
  return eslOK;

  

}


static int
define_consensus(H4_BUILD_CONFIG *cfg, ESL_MSA *msa, char *errbuf, int8_t **ret_matassign)
{
  int8_t *matassign = NULL;
  int     apos, idx;

  ESL_ALLOC(matassign, sizeof(int8_t) * (msa->alen + 1));

  if      (cfg->arch_strategy == h4_ARCH_SYMFRAC)  
    {
      float r      = 0.;  // weighted residue count for a column
      float totwgt = 0.;  // weighted residue+gap count for a column

      for (apos = 1; apos <= msa->alen; apos++)
	{
	  for (idx = 0; idx < msa->nseq; idx++) 
	    {
	      if       (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) { r += msa->wgt[idx]; totwgt += msa->wgt[idx]; }
	      else if  (esl_abc_XIsGap(msa->abc,     msa->ax[idx][apos])) {                     totwgt += msa->wgt[idx]; }
	      else if  (esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos])) continue;
	    }
	  if (r > 0. && r / totwgt >= symfrac) matassign[apos] = TRUE;
	  else                                 matassign[apos] = FALSE;
	}
    }

  else if (cfg->arch_strategy == h4_ARCH_GIVEN)
    {
      if (! msa->rf) ESL_XFAIL(eslEFORMAT, errbuf, "message")

      for (apos = 1; apos <= msa->alen; apos++)
	matassign[apos] = (esl_abc_CIsGap(msa->abc, msa->rf[apos-1])? FALSE : TRUE);
    }

  else
    {

    }

  *ret_matassign = matassign;
  return eslOK;

 ERROR:
  free(matassign);
  return status;
}



static int
collect_counts(H4_BUILD_CONFIG *cfg, ESL_MSA *msa, H4_PROFILE *hmm, char *errbuf, int8_t *matassign)
{
  H4_TRACE *tr = h4_trace_Create();
  int idx;

  for (idx = 0; idx < msa->nseq; idx++)
    {
      h4_trace_Faux (msa, idx, matassign, tr);
      h4_trace_Count(tr, hmm, msa->ax[idx], msa->wgt[idx], tr);

      h4_trace_Reuse();
    }

  h4_trace_Destroy();
  return eslOK;
}


static int
annotate(H4_BUILD_CONFIG *cfg, ESL_MSA *msa, H4_PROFILE *hmm, char *errbuf, int8_t *matassign)
{
  /* transfer from MSA to new HMM: */
  /* msa->rf   */
  /* msa->mm   */
  /* msa->ss_cons */
  /* msa->sa_cons */
  /* map */
  


}



/*****************************************************************
 * x. H4_BUILD_CONFIG
 *****************************************************************/

H4_BUILD_CONFIG *
h4_build_config_Create(void)
{

  cfg->arch_strategy = h4_ARCH_SYMFRAC;
  cfg->symfrac       = 0.5;
  cfg->fragthresh    = 0.5;
  return cfg;
}

void
h4_build_config_Destroy(H4_BUILD_CONFIG *cfg)
{

}
