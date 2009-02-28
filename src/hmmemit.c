/* hmmemit: sample sequence(s) from a profile HMM.
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
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",                0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    "-p", "emit simple consensus sequence",                      0 },
  { "-o",        eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,    NULL, "send sequence output to file <f>, not stdout",        0 },
  { "-p",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    "-c", "sample from profile, not core model",                 0 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0",     NULL,      NULL,    "-c", "number of seqs to sample",                            0 },
  { "-L",        eslARG_INT,    "400", NULL, NULL,      NULL,      "-p",    "-c", "set expected length from profile to <n>",             0 },
  { "--seed",    eslARG_INT,      "0", NULL, "n>=0",    NULL,      NULL,    NULL, "set RNG seed to <n>",                                 0 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile (single)>";
static char banner[] = "sample sequence(s) from a profile HMM";

struct emitcfg_s {
  FILE *ofp;			/* output stream */
  int   nseq;			/* how many sequences to emit */
  int   do_profile;		/* TRUE to emit from implicit profile model, not core model */
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
  P7_PROFILE      *gm      = NULL;     /* profile HMM (scores)                    */
  P7_BG           *bg      = NULL;     /* null model                              */
  P7_TRACE        *tr      = NULL;     /* sampled trace                           */
  ESL_SQ          *sq      = NULL;     /* sampled digital sequence                */
  int              nseq;
  int              fmt;

  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  fmt = eslSQFILE_FASTA;

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h")) {
    p7_banner(stdout, argv[0], banner);
    esl_usage(stdout, argv[0], usage);
    puts("\nwhere options are:\n");
    esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=all docgroups; 2 = indentation; 80=textwidth*/
    return eslOK;
  }

  cfg.do_profile   = esl_opt_GetBoolean(go, "-p");
  cfg.nseq         = esl_opt_GetInteger(go, "-N");
  cfg.L            = esl_opt_GetInteger(go, "-L");

  if ( esl_opt_IsOn(go, "-o"))
    {
      cfg.ofp = fopen(esl_opt_GetString(go, "-o"), "w");
      if (cfg.ofp == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
    }
  else cfg.ofp = stdout;

  if (esl_opt_ArgNumber(go) != 1) {
    puts("Incorrect number of command line arguments.");
    esl_usage(stdout, argv[0], usage);
    return eslFAIL;
  }
  hmmfile = esl_opt_GetArg(go, 1);
  if (hmmfile == NULL) esl_fatal("failed to get <hmmfile> on cmdline: %s\n", go->errbuf);

  /*****************************************************************
   * Initializations, including opening the HMM file
   *****************************************************************/

  r = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));

  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open HMM file %s for reading.\n",                   hmmfile);
  else if (status == eslEFORMAT)   p7_Fail("File %s does not appear to be in a recognized HMM format.\n", hmmfile);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n",       status, hmmfile);  

  else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", hmmfile);
    
  if ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslOK) {
    if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
    else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
    else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
    else                             esl_fatal("Unexpected error in reading HMMs");
  }
   
  /* init sq (need to know abc to do this */
  if (sq == NULL) sq = esl_sq_CreateDigital(abc);
  if (sq == NULL) esl_fatal("Failed to allocated sequence");

  if (esl_opt_GetBoolean(go, "-c")) 
    {
      if (p7_emit_SimpleConsensus(hmm, sq)              != eslOK) esl_fatal("failed to create simple consensus seq");
      if (esl_sq_SetName(sq, "%s-consensus", hmm->name) != eslOK) esl_fatal("Failed to set sequence name");
      if (esl_sqio_Write(cfg.ofp, sq, eslSQFILE_FASTA)   != eslOK) esl_fatal("Failed to write sequence");
    }
  else
    {
      if ((tr = p7_trace_Create())              == NULL)  esl_fatal("Failed to allocate trace");
      if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");
      if ((gm = p7_profile_Create(hmm->M, abc)) == NULL)  esl_fatal("failed to create profile");

      if (p7_ProfileConfig(hmm, bg, gm, cfg.L, p7_LOCAL) != eslOK) esl_fatal("failed to configure profile");
      if (p7_bg_SetLength(bg, cfg.L)                     != eslOK) esl_fatal("failed to reconfig null model length");
      if (p7_hmm_Validate    (hmm, NULL, 0.0001)         != eslOK) esl_fatal("whoops, HMM is bad!");
      if (p7_profile_Validate(gm,  NULL, 0.0001)         != eslOK) esl_fatal("whoops, profile is bad!");

      for (nseq = 1; nseq <= cfg.nseq; nseq++) 
	{
	  if (cfg.do_profile) {
	    status = p7_ProfileEmit(r, hmm, gm, bg, sq, tr);
	    if (status != eslOK) esl_fatal("Failed to emit sequence from hmm\n");
	  } else {
	    status = p7_CoreEmit(r, hmm, sq, tr);
	    if (status != eslOK) esl_fatal("Failed to emit sequence from hmm\n");
	  }
      
	  status = esl_sq_SetName(sq, "%s-sample%d", hmm->name, nseq);
	  if (status != eslOK) esl_fatal("Failed to set sequence name\n");

	  status = esl_sqio_Write(cfg.ofp, sq, eslSQFILE_FASTA);
	  if (status != eslOK) esl_fatal("Failed to write sequence\n");
	}

      p7_profile_Destroy(gm);
      p7_bg_Destroy(bg);
      p7_trace_Destroy(tr);
    }

  if (esl_opt_IsOn(go, "-o")) { fclose(cfg.ofp); }
  esl_sq_Destroy(sq);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  p7_hmmfile_Close(hfp);
  p7_hmm_Destroy(hmm);
  return eslOK;
}


