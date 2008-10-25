/* hmmconvert: converting profile HMM files to HMMER3 HMM format.
 * 
 * SRE, Thu Oct 16 08:57:43 2008 [janelia] [Enigma MCMXC a.D.]
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
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",          0 },
  { "-a",        eslARG_NONE,"default",NULL, NULL,   "-a,-b",      NULL,    NULL, "ascii:  output models in HMMER3 ASCII format",  0 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,   "-a,-b",      NULL,    NULL, "binary: output models in HMMER3 binary format", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "convert a profile HMM file to a HMMER3 format";


int 
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_ALPHABET  *abc     = NULL;
  char          *hmmfile = esl_opt_GetArg(go, 1);
  P7_HMMFILE    *hfp     = NULL;
  P7_HMM        *hmm     = NULL;
  FILE          *ofp     = stdout;
  int            status;

  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open HMM file %s for reading.\n",     hmmfile);
  else if (status == eslEFORMAT)   p7_Fail("File %s does not appear to be in a recognized HMM format.\n", hmmfile);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n", status, hmmfile);  

  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) == eslOK)
    {
      if      (esl_opt_GetBoolean(go, "-a") == TRUE) p7_hmmfile_WriteASCII (ofp, hmm);
      else if (esl_opt_GetBoolean(go, "-b") == TRUE) p7_hmmfile_WriteBinary(ofp, hmm);

      p7_hmm_Destroy(hmm);
    }
  if      (status == eslEFORMAT)   p7_Fail("bad file format in HMM file %s",             hmmfile);
  else if (status == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets",   hmmfile);
  else if (status != eslEOF)       p7_Fail("Unexpected error in reading HMMs from %s",   hmmfile);

  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
