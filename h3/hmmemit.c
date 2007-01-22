/* main() for emitting sequences from a profile HMM.
 * 
 * SRE, Tue Jan  9 13:22:53 2007 [Janelia] [Verdi, Requiem]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",     0 },
  { "-n",        eslARG_INT,      "1", NULL, "n>0",     NULL,      NULL,    NULL, "number of seqs to sample",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "hmmemit [-options] <hmmfile>";

static struct emitcfg_s {
  int do_nseq;			/* how many sequences to emit */
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

  esl_opt_GetIntegerOption(go, "-n", &(cfg.do_nseq));

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

  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) == eslOK)  
    {
      /* First time through, init sq (need to know abc to do this */
      if (sq == NULL) sq = esl_sq_CreateDigital(abc);
      if (sq == NULL) esl_fatal("Failed to allocated sequence");

      if (p7_hmm_Validate(hmm, 0.0001) != eslOK) esl_fatal("whoops, HMM is bad!");

      for (nseq = 1; nseq <= cfg.do_nseq; nseq++) 
	{
	  status = p7_CoreEmit(r, hmm, sq, tr);
	  if (status != eslOK) esl_fatal("Failed to emit sequence from hmm\n");
      
	  sprintf(sqname, "%s-sample%d", hmm->name, nseq);
	  status = esl_sq_SetName(sq, sqname);
	  if (status != eslOK) esl_fatal("Failed to set sequence name\n");

	  status = esl_sqio_Write(stdout, sq, eslSQFILE_FASTA);
	  if (status != eslOK) esl_fatal("Failed to write sequence\n");
	}
    }
  if (status != eslEOF) {
    if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
    else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
    else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
    else                             esl_fatal("Unexpected error in reading HMMs");
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
