#include <h4_config.h>

#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_cpu.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_subcmd.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"

#include "ssvfilter.h"

static ESL_OPTIONS kiteline_options[] = {
   /* name             type      default                         env  range      toggles reqs incomp  help                                                docgroup*/
  { "-h",           eslARG_NONE,  FALSE,                        NULL, NULL,      NULL, NULL, NULL,     "show brief help on version and usage",                     1 },
};


int
h4_cmd_kiteline(const char *topcmd, const ESL_SUBCMD *sub, int argc, char **argv)
{
  ESL_GETOPTS  *go      = esl_subcmd_CreateDefaultApp(topcmd, sub, kiteline_options, argc, argv);
  char         *hmmfile = esl_opt_GetArg(go, 1);
  char         *seqfile = esl_opt_GetArg(go, 2);
  int           infmt   = eslSQFILE_UNKNOWN;
  H4_HMMFILE   *hfp     = NULL;
  H4_PROFILE   *hmm     = NULL;
  H4_MODE        *mo    = h4_mode_Create();   // SSV doesn't use a mode. We only use this for mo->nullsc.
  ESL_ALPHABET *abc     = NULL;
  ESL_SQFILE   *sqfp    = NULL;
  ESL_SQ       *sq      = NULL;
  //H4_FILTERMX  *fx      = h4_filtermx_Create(100);
  float         sfraw, sfscore;
  //float         vfraw, vfscore;
  int           status;

  status = h4_hmmfile_Open(hmmfile, NULL, &hfp);
  if (status != eslOK) esl_fatal("Error: failed to open %s for reading profile HMM(s)\n%s\n", strcmp(hmmfile, "-") == 0? "<stdin>" : hmmfile, hfp->errmsg);

  status = h4_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT)   esl_fatal("Parse failed, bad profile HMM file format in %s:\n   %s",  strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status == eslEOF)       esl_fatal("Empty input? No profile HMM found in %s\n",                strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, hfp->errmsg);
  else if (status != eslOK)        esl_fatal("Unexpected error reading profile HMM from %s (code %d)\n", strcmp(hmmfile, "-") == 0 ? "<stdin>" : hmmfile, status);

  status = esl_sqfile_OpenDigital(abc, seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open %s for reading sequences\n",             strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile);
  else if (status == eslEFORMAT)   esl_fatal("Format of sequence file %s unrecognized\n",             strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile);
  else if (status != eslOK)        esl_fatal("Unexpected error opening sequence file %s (code %d)\n", strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, status);

  sq = esl_sq_CreateDigital(abc);

  while ((status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) esl_fatal("Error: failed to parse sequence format of %s\n%s\n",    strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, esl_sqfile_GetErrorBuf(sqfp));
      else if (status != eslOK)      esl_fatal("Unexpected error reading sequence from %s (code %d)\n", strcmp(seqfile, "-") == 0 ? "<stdin>" : seqfile, status);

      h4_mode_SetLength(mo, sq->n);

      h4_ssvfilter(sq->dsq, sq->n, hmm, &sfraw);
      sfscore = sfraw - mo->nullsc;

      //      h4_vitfilter(sq->dsq, sq->n, hmm, mo, fx, &vfraw);
      //vfscore = vfraw = mo->nullsc;

      printf("%-30s %10.2f\n", sq->name, sfscore);

      esl_sq_Reuse(sq);
    }
         

  esl_sqfile_Close(sqfp);
  h4_hmmfile_Close(hfp);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}

