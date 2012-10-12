/* sim-alignment: 
 * simulation testbed for evaluating alignment accuracy.
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_refmx.h"
#include "p7_sparsemx.h"
#include "p7_trace_metrics.h"

#include "reference_viterbi.h"
#include "sparse_fwdback.h"
#include "sparse_viterbi.h"

#include "impl_sse/impl_sse.h"
#include "impl_sse/p7_filtermx.h"


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",                          0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                                 0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of sequences sampled per input model",                   0 },
  { "--gpglocal",eslARG_REAL,   "1.0", NULL, NULL,  NULL,  NULL, NULL, "local/glocal prob param, generation; 0.5=dual-mode; 0=local",   0 },
  { "--gnj",     eslARG_REAL,   "1.0", NULL, NULL,  NULL,  NULL, NULL, "uni/multihit expected # param, generation; 1=multi; 0=uni",     0 },
  { "--gL",      eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "length config param, generation",                               0 },
  { "--apglocal",eslARG_REAL,   "0.5", NULL, NULL,  NULL,  NULL, NULL, "local/glocal prob param, alignment; 0.5=dual-mode; 0=local",    0 },
  { "--anj",     eslARG_REAL,   "1.0", NULL, NULL,  NULL,  NULL, NULL, "uni/multihit expected # param, alignment; 1=multi; 0=uni",      0 },
  { "--vit",     eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "test Viterbi alignment, not banded MEG alignment",              0 },
  { "--full",    eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "unrestricted sparse DP",                                        0 },
  { "--dumpseqs",eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump the generated sequences",                                  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <gen-hmmfile> <ali-hmmfile>";
static char banner[] = "simulation testbed for evaluating alignment accuracy";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS      *go       = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS   *rng      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET     *abc      = NULL;
  char             *ghmmfile = esl_opt_GetArg(go, 1); /* HMMs parameterized for sequence generation */
  char             *ahmmfile = esl_opt_GetArg(go, 2); /* HMMs parameterized for alignment */
  int               N        = esl_opt_GetInteger(go, "-N");
  P7_HMMFILE       *ghfp     = NULL;
  P7_HMMFILE       *ahfp     = NULL;
  P7_HMM           *ghmm     = NULL;
  P7_HMM           *ahmm     = NULL;
  P7_PROFILE       *ggm      = NULL;
  P7_PROFILE       *agm      = NULL;
  P7_OPROFILE      *aom      = NULL;
  P7_BG            *bg       = NULL;
  ESL_SQ           *sq       = NULL;
  P7_TRACE         *reftr    = p7_trace_Create();
  P7_TRACE         *testtr   = p7_trace_Create();
  P7_TRACE_METRICS *tmetrics = p7_trace_metrics_Create();
  P7_REFMX         *rmx      = p7_refmx_Create(100,100);
  //  P7_FILTERMX      *ox       = NULL;
  P7_SPARSEMASK    *sm       = p7_sparsemask_Create(100, 100);
  P7_SPARSEMX      *sxv      = p7_sparsemx_Create(NULL);
  int               idx;
  char              errbuf[eslERRBUFSIZE];
  int               status;
  
  /* Initialize log-sum calculator */
  impl_Init();
  p7_FLogsumInit();

  /* open HMM file containing models parameterized for generation (sampling) of seqs */
  status = p7_hmmfile_OpenE(ghmmfile, NULL, &ghfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", ghmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                ghmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",                       status, ghmmfile, errbuf);  

  /* open HMM file containing models parameterized for alignment (may be the same as ghmmfile) */
  status = p7_hmmfile_OpenE(ahmmfile, NULL, &ahfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", ahmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                ahmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",                       status, ahmmfile, errbuf);  

  while ( (status = p7_hmmfile_Read(ghfp, &abc, &ghmm)) == eslOK) /* <abc> gets set on first read  */
    {
      /* read the counterpart HMM from <ahfp> */
      status = p7_hmmfile_Read(ahfp, &abc, &ahmm);
      if      (status == eslEFORMAT)   p7_Fail("Bad file format in HMM file %s:\n%s\n",          ahfp->fname, ahfp->errbuf);
      else if (status == eslEINCOMPAT) p7_Fail("HMM in %s is not in the expected %s alphabet\n", ahfp->fname, esl_abc_DecodeType(abc->type));
      else if (status == eslEOF)       p7_Fail("Empty HMM file %s? No HMM data found.\n",        ahfp->fname);
      else if (status != eslOK)        p7_Fail("Unexpected error in reading HMMs from %s\n",     ahfp->fname);

      /* try to validate that they're the "same" */
      if (ahmm->M != ghmm->M || strcmp(ahmm->name, ghmm->name) != 0) p7_Fail("<gen-hmmfile>, <ali-hmmfile> contain different set or order of models");

      /* deferred one-time creation of structures that need to know the alphabet */
      if (!bg) bg = p7_bg_Create(abc);
      if (!sq) sq = esl_sq_CreateDigital(abc);

      ggm = p7_profile_Create(ghmm->M,  abc);
      agm = p7_profile_Create(ahmm->M,  abc);
      aom = p7_oprofile_Create(ahmm->M, abc);

      p7_profile_ConfigCustom(ggm, ghmm, bg, esl_opt_GetInteger(go, "--gL"), esl_opt_GetReal(go, "--gnj"), esl_opt_GetReal(go, "--gpglocal"));
      p7_profile_ConfigCustom(agm, ahmm, bg, 100,                            esl_opt_GetReal(go, "--anj"), esl_opt_GetReal(go, "--apglocal"));
      p7_oprofile_Convert(agm, aom);

      for (idx = 1; idx <= N; idx++)
	{
	  p7_ProfileEmit(rng, ghmm, ggm, bg, sq, reftr);

	  if (esl_opt_GetBoolean(go, "--dumpseqs")) {
	    esl_sq_FormatName(sq, "seq%d", idx);
	    esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, FALSE);
	  }

	  p7_bg_SetLength(bg, sq->n);
	  p7_profile_SetLength(agm, sq->n);
	  p7_sparsemask_Reinit(sm, agm->M, sq->n);
	  p7_sparsemask_AddAll(sm);

	  if (esl_opt_GetBoolean(go, "--vit"))  p7_ReferenceViterbi(sq->dsq, sq->n, agm,     rmx, testtr, /*opt_vsc=*/NULL);
	  else                         	        p7_SparseViterbi   (sq->dsq, sq->n, agm, sm, sxv, testtr, /*opt_vsc=*/NULL);

	  p7_trace_metrics(reftr, testtr, tmetrics);

	  p7_sparsemask_Reuse(sm);
	  p7_sparsemx_Reuse(sxv);
	  //p7_filtermx_Reuse(ox);
	  p7_refmx_Reuse(rmx);
	  esl_sq_Reuse(sq);
	  p7_trace_Reuse(reftr);
	  p7_trace_Reuse(testtr);
	}

      p7_oprofile_Destroy(aom);
      p7_profile_Destroy(ggm);
      p7_profile_Destroy(agm);
      p7_hmm_Destroy(ghmm);
      p7_hmm_Destroy(ahmm);
    }
  /* we leave the loop with <status> set by a p7_hmmfile_Read() on ghfp; if all is well, status=eslEOF */
  if      (status == eslEFORMAT)   p7_Fail("Bad file format in HMM file %s:\n%s\n",          ghfp->fname, ghfp->errbuf);
  else if (status == eslEINCOMPAT) p7_Fail("HMM in %s is not in the expected %s alphabet\n", ghfp->fname, esl_abc_DecodeType(abc->type));
  else if (status != eslEOF)       p7_Fail("Unexpected error in reading HMMs from %s\n",     ghfp->fname);
  
  p7_trace_metrics_Dump(stdout, tmetrics);

  p7_hmmfile_Close(ghfp);  
  p7_hmmfile_Close(ahfp);
  //  p7_filtermx_Destroy(ox);
  p7_sparsemask_Destroy(sm);
  p7_sparsemx_Destroy(sxv);
  p7_refmx_Destroy(rmx);
  p7_trace_metrics_Destroy(tmetrics);
  p7_trace_Destroy(testtr);
  p7_trace_Destroy(reftr);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
}

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
