#include "h4_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_subcmd.h"

#include "h4_profile.h"
#include "modelio.h"

#include "h4_build.h"


static ESL_OPTIONS build_options[] = {
  /* name             type          default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,       FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",        0 },
  { "--informat",  eslARG_STRING,      NULL,  NULL, NULL,  NULL,  NULL, NULL, "specify the input MSA file is in format <s>", 0 }, 
  { "--dna",       eslARG_NONE,       FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use DNA alphabet",                            0 },
  { "--rna",       eslARG_NONE,       FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use RNA alphabet",                            0 },
  { "--amino",     eslARG_NONE,       FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use protein alphabet",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
h4_cmd_build(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS  *go      = esl_subcmd_CreateDefaultApp(topcmd, sub, build_options, argc, argv);
  char         *msafile = esl_opt_GetArg(go, 1);
  //  char         *hmmfile = esl_opt_GetArg(go, 2);
  int           infmt   = eslMSAFILE_UNKNOWN;
  ESL_ALPHABET *abc     = NULL;
  ESL_MSAFILE  *afp     = NULL;
  ESL_MSA      *msa     = NULL;
  H4_PROFILE   *hmm     = NULL;
  int           nali    = 0;
  int           status;

  if      (esl_opt_GetBoolean(go, "--rna"))   abc = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna"))   abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--amino")) abc = esl_alphabet_Create(eslAMINO); 

  if (esl_opt_IsOn(go, "--informat") &&
      (infmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"))) == eslMSAFILE_UNKNOWN)
    esl_fatal("%s is not a valid MSA file format for --informat", esl_opt_GetString(go, "--informat"));
  // SRE: check the esl_fatal().

  status = esl_msafile_Open(&abc, msafile, NULL, infmt, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  while ((status = esl_msafile_Read(afp, &msa)) == eslOK)
    {
      nali++;

      h4_Build(NULL, msa, &hmm, NULL);
  
      h4_modelio_Write(stdout, hmm);

      esl_msa_Destroy(msa);
      h4_profile_Destroy(hmm);
    }
  if (nali == 0 || status != eslEOF) esl_msafile_ReadFailure(afp, status); /* a convenience, like esl_msafile_OpenFailure() */
  esl_msafile_Close(afp);
  
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}


