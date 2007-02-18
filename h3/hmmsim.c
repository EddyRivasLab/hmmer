/* main() for scoring a profile HMM against simulated random sequences
 * 
 * SRE, Thu Feb  8 14:45:23 2007 [UA 1130 from Vancouver]
 * SVN $Id$
 */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0",     NULL,      NULL,    NULL, "number of random seqs",                  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0",     NULL,      NULL,    NULL, "length of random seqs",                  0 },
  { "--h2",      eslARG_NONE,    NULL, NULL, NULL,      NULL,      "-p",    NULL, "configure profile in old HMMER2 style",  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "hmmsim [-options] <hmmfile input>";

int
main(int argc, char **argv)
{
  int status;			       /* status of a function call               */
  ESL_ALPHABET    *abc     = NULL;     /* sequence alphabet                       */
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                 */
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
  char            *hmmfile = NULL;     /* file to read HMM(s) from                */
  P7_HMMFILE      *hfp     = NULL;     /* open hmmfile                            */
  P7_HMM          *hmm     = NULL;     /* HMM to emit from                        */
  P7_TRACE        *tr      = NULL;     /* sampled trace                           */
  P7_GMX          *mx      = NULL;     /* DP matrix                               */
  ESL_DSQ         *dsq     = NULL;     /* sampled digital sequence                */
  int              N;		       /* number of seqs to generate              */
  int              L;		       /* length of generated seqs                */
  int              do_oldconfig;       /* TRUE to use H2 exit/entry configuration */
  float            sc;		       /* a Viterbi score                         */
  int              nseq;	       /* counter over sequences                  */

  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  if (esl_opt_IsSet(go, "-h")) {
    puts(usage);
    puts("\n  where options are:\n");
    esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=all docgroups; 2 = indentation; 80=textwidth*/
    return eslOK;
  }

  esl_opt_GetIntegerOption(go, "-N",        &N);
  esl_opt_GetIntegerOption(go, "-L",        &L);
  esl_opt_GetBooleanOption(go, "--h2",      &do_oldconfig);

  if (esl_opt_ArgNumber(go) != 1) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return eslFAIL;
  }
  hmmfile = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL); /* NULL=no range checking */


  /*****************************************************************
   * Initializations, including opening and reading the HMM file 
   *****************************************************************/

  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if (status == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n", hmmfile);
  else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", hmmfile);

  if ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslOK) {
    if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
    else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
    else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
    else                             esl_fatal("Unexpected error in reading HMMs");
  }

  if ((r       = esl_randomness_CreateTimeseeded()) == NULL)  esl_fatal("Failed to create rng");
  if ((status  = p7_trace_Create(256, &tr))         != eslOK) esl_fatal("Failed to allocate trace");
  if ((mx      = p7_gmx_Create(hmm->M, L))          == NULL)  esl_fatal("failed to create dp matrix");
  if ((hmm->bg = p7_bg_Create(abc))                 == NULL)  esl_fatal("failed to create null model");
  if ((hmm->gm = p7_profile_Create(hmm->M, abc))    == NULL)  esl_fatal("failed to create profile");
  if ((dsq     = malloc(sizeof(ESL_DSQ) * (L+2)))   == NULL)  esl_fatal("failed to create dsq");

  if (do_oldconfig) {
    if (p7_H2_ProfileConfig(hmm, hmm->gm, p7_LOCAL) != eslOK) esl_fatal("failed to configure profile");
  } else {
    if (p7_ProfileConfig(hmm, hmm->gm, p7_LOCAL)    != eslOK) esl_fatal("failed to configure profile");
    if (p7_ReconfigLength(hmm->gm,  L)              != eslOK) esl_fatal("failed to reconfig profile L");
    if (p7_hmm_Validate    (hmm,     0.0001)        != eslOK) esl_fatal("whoops, HMM is bad!");
    if (p7_profile_Validate(hmm->gm, 0.0001)        != eslOK) esl_fatal("whoops, profile is bad!");
  }

  for (nseq = 1; nseq <= N; nseq++) 
    {
      if (esl_rnd_xfIID(r, hmm->bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");
      if (p7_Viterbi(dsq, L, hmm->gm, mx, tr, &sc)     != eslOK) esl_fatal("viterbi failed");
      printf("score = %6.2f bits\n", sc);
    }

  p7_profile_Destroy(hmm->gm);
  p7_bg_Destroy(hmm->bg);
  p7_gmx_Destroy(mx);
  free(dsq);
  p7_trace_Destroy(tr);
  esl_randomness_Destroy(r);
  p7_hmmfile_Close(hfp);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;


}
