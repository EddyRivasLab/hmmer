/* SSI indexing of an HMM file.
 * 
 * SRE, Sun Apr 22 09:07:30 2007 [Janelia]
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
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "hmmindex [-options] <hmmfile>";

int
main(int argc, char **argv)
{
  int status;
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                 */
  char            *hmmfile = NULL;     /* file to read HMM(s) from                */
  P7_HMMFILE      *hfp     = NULL;    
  ESL_ALPHABET    *abc     = NULL;
  P7_HMM          *hmm     = NULL;

  /* Process command line options.
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    puts(usage);
    puts("\ngeneral options are:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 2 = indentation; 80=textwidth*/
    return 0;
  }
  hmmfile = esl_opt_GetArg(go, 1);
  if (hmmfile == NULL)        esl_fatal("Failed to get <hmmfile> argument on command line.");

  /* Initializations, including opening the HMM file.
   */
  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if (status == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n", hmmfile);
  else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", hmmfile);

  /* Main loop over the HMMs
   */
  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) == eslOK)
    {
      printf("Read %s: %d nodes\n", hmm->name, hmm->M);
      p7_hmm_Destroy(hmm);
    }
  if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
  else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
  else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
  else if (status != eslEOF)       esl_fatal("Unexpected error in reading HMMs from %s", hmmfile);

  
  esl_alphabet_Destroy(abc);
  p7_hmmfile_Close(hfp);
  return 0;
}
