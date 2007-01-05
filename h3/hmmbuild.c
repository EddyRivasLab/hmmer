/* main() for profile HMM construction from a multiple sequence alignment
 * 
 * SRE, Wed Jan  3 11:03:47 2007 [Janelia] [The Chemical Brothers, Push the Button]
 * SVN $Id$
 */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "hmmbuild [-options] <hmmfile output> <alignment file input>";


int
main(int argc, char **argv)
{
  int              status;	/* status of a function call               */
  ESL_ALPHABET    *abc;		/* sequence alphabet                       */
  ESL_GETOPTS     *go;		/* command line processing                 */
  char            *alifile;     /* seqfile to read alignment from          */
  int              fmt;	        /* format of the alifile                   */
  ESL_MSAFILE     *afp;         /* open alignment file                     */
  ESL_MSA         *msa;         /* a multiple sequence alignment           */
  int              nali;	/* count number of alignments/HMMs         */
  char            *hmmfile;     /* file to write HMM to                    */
  FILE            *hmmfp;       /* HMM output file handle                  */
  P7_HMM          *hmm;         /* constructed HMM; written to hmmfile     */
  float            symfrac;	/* controls fast modelmaking               */

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
    return 0;
  }
  if (esl_opt_ArgNumber(go) != 2) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return 1;
  }
  hmmfile = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL); /* NULL=no range checking */
  alifile = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL);
  fmt     = eslMSAFILE_UNKNOWN;  /* autodetect alignment format by default. */

  symfrac = 0.5;

  /*****************************************************************
   * Set up the alphabet
   *****************************************************************/
  
  abc = esl_alphabet_Create(eslAMINO);

  /*****************************************************************
   * Open the alignment file (it might have >1 alignment)
   *****************************************************************/

  status = esl_msafile_OpenDigital(abc, alifile, fmt, NULL, &afp); /* NULL= no database dir from the environment */
  if      (status == eslENOTFOUND) esl_fatal("Alignment file %s doesn't exist or isn't readable.\n",     alifile);
  else if (status == eslEFORMAT)   esl_fatal("Couldn't determine format of alignment file %s.\n",        alifile);
  else if (status != eslOK)        esl_fatal("Alignment file open unexpectedly failed with error %d.\n", status);


  /*****************************************************************
   * Open the HMM output file.
   *****************************************************************/

  hmmfp = fopen(hmmfile, "w");
  if (hmmfp == NULL) esl_fatal("Failed to open HMM file %s for writing", hmmfile);

  /*****************************************************************
   * Read alignments one at a time, build HMMs, and save them.
   *****************************************************************/

  nali = 0;
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;

      printf("Working on %s...\n", msa->name);

      status = p7_Fastmodelmaker(msa, symfrac, &hmm, NULL);
      if (status != eslOK) esl_fatal("Model construction failed.");

      status = p7_hmmfile_Write(hmmfp, hmm);
      if (status != eslOK) esl_fatal("Failed to write model to disk.");

      p7_hmm_Destroy(hmm);
      esl_msa_Destroy(msa);
    }
  if (status == eslEFORMAT) 
    esl_fatal("\
Alignment file parse error, line %d of file %s:\n\
%s\n\
Offending line is:\n\
%s\n", afp->linenumber, afp->fname, afp->errbuf, afp->buf);
  else if (status != eslEOF)
    esl_fatal("Alignment file read unexpectedly failed with code %d\n", status);
      
  esl_msafile_Close(afp);
  esl_getopts_Destroy(go);
  esl_alphabet_Destroy(abc);
  return 0;
}
