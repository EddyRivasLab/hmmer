/* hmmconvert: converting profile HMM files to HMMER3 HMM format.
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { (char *) "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,       NULL,    NULL,(char *)  "show brief help on version and usage",                             0 },
  { (char *) "-a",        eslARG_NONE,(char *) "default",NULL, NULL,(char *)  "-a,-b,-2",      NULL,    NULL, (char *) "ascii:  output models in HMMER3 ASCII format",                     0 },
  { (char *) "-b",        eslARG_NONE,   FALSE, NULL, NULL,(char *)  "-a,-b,-2",      NULL,    NULL,(char *)  "binary: output models in HMMER3 binary format",                    0 },
  { (char *) "-2",        eslARG_NONE,   FALSE, NULL, NULL,(char *)  "-a,-b,-2",      NULL,    NULL, (char *) "HMMER2: output backward compatible HMMER2 ASCII format (ls mode)", 0 },
  {(char *)  "--outfmt",  eslARG_STRING, NULL,  NULL, NULL,      NULL,       NULL,  (char *)   "-2", (char *) "choose output legacy 3.x file formats by name, such as '3/a'",     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "convert profile file to a HMMER format";


int 
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_ALPHABET  *abc     = NULL;
  char          *hmmfile = esl_opt_GetArg(go, 1);
  P7_HMMFILE    *hfp     = NULL;
  P7_HMM        *hmm     = NULL;
  FILE          *ofp     = stdout;
  char          *outfmt  = esl_opt_GetString(go, (char *) "--outfmt");
  int            fmtcode = -1;	/* -1 = write the current default format */
  int            status;
  char           errbuf[eslERRBUFSIZE];

  if (outfmt != NULL) {
    if      (strcmp(outfmt, "3/a") == 0) fmtcode = p7_HMMFILE_3a;
    else if (strcmp(outfmt, "3/b") == 0) fmtcode = p7_HMMFILE_3b;
    else if (strcmp(outfmt, "3/c") == 0) fmtcode = p7_HMMFILE_3c;
    else    p7_Fail((char *)  "No such 3.x output format code %s.\n", outfmt);
  }

  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail((char *) "File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail((char *) "File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail((char *) "Unexpected error %d in opening HMM file %s.\n%s\n",                       status, hmmfile, errbuf);  

  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) == eslOK)
    {
      if      (esl_opt_GetBoolean(go,(char *)  "-a") == TRUE) p7_hmmfile_WriteASCII (ofp, fmtcode, hmm);
      else if (esl_opt_GetBoolean(go,(char *)  "-b") == TRUE) p7_hmmfile_WriteBinary(ofp, fmtcode, hmm);
      else if (esl_opt_GetBoolean(go,(char *)  "-2") == TRUE) p7_h2io_WriteASCII    (ofp, hmm);

      p7_hmm_Destroy(hmm);
    }
  if      (status == eslEFORMAT)   p7_Fail((char *) "bad file format in HMM file %s",             hmmfile);
  else if (status == eslEINCOMPAT) p7_Fail((char *) "HMM file %s contains different alphabets",   hmmfile);
  else if (status != eslEOF)       p7_Fail((char *) "Unexpected error in reading HMMs from %s",   hmmfile);

  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
