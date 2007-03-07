/* main() for emitting sequences from a profile HMM.
 * 
 * SRE, Tue Jan  9 13:22:53 2007 [Janelia] [Verdi, Requiem]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",      0 },
  { "-n",        eslARG_INT,      "1", NULL, "n>0",     NULL,      NULL,    NULL, "number of seqs to sample",                  0 },
  { "-p",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "sample from profile, not core model",       0 },
  { "-L",        eslARG_INT,    "400", NULL, NULL,      NULL,      "-p",    NULL, "set expected length from profile to <n>",   0 },
  { "--h2",      eslARG_NONE,    NULL, NULL, NULL,      NULL,      "-p",    NULL, "configure profile in old HMMER2 style",     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "hmmemit [-options] <hmmfile (single)>";

struct emitcfg_s {
  int   nseq;			/* how many sequences to emit */
  int   do_profile;		/* TRUE to emit from implicit profile model, not core model */
  int   do_oldconfig;		/* TRUE to do old-style HMMER2 entry/exit configuration */
  int   L;			/* expected length from a profile */
};


int
main(int argc, char **argv)
{
  struct emitcfg_s cfg;		       /* command-line configurations             */
  int              status;	       /* status of a function call               */
  ESL_ALPHABET    *abc     = NULL;     /* sequence alphabet                       */
  ESL_GETOPTS     *go      = NULL;     /* command line processing                 */
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
  char            *hmmfile = NULL;     /* file to read HMM(s) from                */
  P7_HMMFILE      *hfp     = NULL;     /* open hmmfile                            */
  P7_HMM          *hmm     = NULL;     /* HMM to emit from                        */
  P7_TRACE        *tr      = NULL;     /* sampled trace                           */
  ESL_SQ          *sq      = NULL;     /* sampled digital sequence                */
  char             sqname[64];
  int              nseq;
  int              fmt;

  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  fmt = eslSQFILE_FASTA;

  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  if (esl_opt_IsSet(go, "-h")) {
    puts(usage);
    puts("\n  where options are:\n");
    esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=all docgroups; 2 = indentation; 80=textwidth*/
    return eslOK;
  }

  esl_opt_GetIntegerOption(go, "-n",        &(cfg.nseq));
  esl_opt_GetBooleanOption(go, "-p",        &(cfg.do_profile));
  esl_opt_GetIntegerOption(go, "-L",        &(cfg.L));
  esl_opt_GetBooleanOption(go, "--h2",      &(cfg.do_oldconfig));

  if (esl_opt_ArgNumber(go) != 1) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return eslFAIL;
  }
  hmmfile = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL); /* NULL=no range checking */

  /*****************************************************************
   * Initializations, including opening the HMM file
   *****************************************************************/

  if ((r = esl_randomness_CreateTimeseeded()) == NULL)
    esl_fatal("Failed to create random number generator: probably out of memory");

  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if (status == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n", hmmfile);
  else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", hmmfile);
    
  status = p7_trace_Create(256, &tr);
  if (status != eslOK) esl_fatal("Failed to allocate trace\n");

  if ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslOK) {
    if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
    else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
    else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
    else                             esl_fatal("Unexpected error in reading HMMs");
  }
   
  /* init sq (need to know abc to do this */
  if (sq == NULL) sq = esl_sq_CreateDigital(abc);
  if (sq == NULL) esl_fatal("Failed to allocated sequence");

  if ((hmm->bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to created null model");
  if ((hmm->gm = p7_profile_Create(hmm->M, abc)) == NULL)  esl_fatal("failed to create profile");

  if (cfg.do_oldconfig) {
    if (p7_H2_ProfileConfig(hmm, hmm->gm, p7_LOCAL) != eslOK) esl_fatal("failed to configure profile");
  } else {
    if (p7_ProfileConfig(hmm, hmm->gm, p7_LOCAL)    != eslOK) esl_fatal("failed to configure profile");
    if (p7_ReconfigLength(hmm->gm, cfg.L)           != eslOK) esl_fatal("failed to reconfig profile L");
    if (p7_hmm_Validate    (hmm,     0.0001, NULL)  != eslOK) esl_fatal("whoops, HMM is bad!");
    if (p7_profile_Validate(hmm->gm, 0.0001)       != eslOK) esl_fatal("whoops, profile is bad!");
  }

  for (nseq = 1; nseq <= cfg.nseq; nseq++) 
    {
      if (cfg.do_profile && cfg.do_oldconfig) {
	status = p7_H2_ProfileEmit(r, hmm->gm, sq, tr);
	if (status != eslOK) esl_fatal("Failed to emit sequence from hmm\n");
      } else if (cfg.do_profile) {
	status = p7_ProfileEmit(r, hmm->gm, sq, tr);
	if (status != eslOK) esl_fatal("Failed to emit sequence from hmm\n");
      } else {
	status = p7_CoreEmit(r, hmm, sq, tr);
	if (status != eslOK) esl_fatal("Failed to emit sequence from hmm\n");
      }
      
      sprintf(sqname, "%s-sample%d", hmm->name, nseq);
      status = esl_sq_SetName(sq, sqname);
      if (status != eslOK) esl_fatal("Failed to set sequence name\n");

      status = esl_sqio_Write(stdout, sq, eslSQFILE_FASTA);
      if (status != eslOK) esl_fatal("Failed to write sequence\n");
    }

  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  p7_hmmfile_Close(hfp);
  p7_hmm_Destroy(hmm);
  return eslOK;
}


