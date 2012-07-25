/* hmmstat: display summary statistics for an HMM database.
 * 
 * Example:
 *  ./hmmstat Pfam
 *  
 * SRE, Thu May 24 11:18:20 2007
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_exponential.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type       default   env  range    toggles    reqs       incomp  help   docgroup*/
  { "-h",        eslARG_NONE,    FALSE,  NULL, NULL,    NULL,  NULL,           NULL, "show brief help on version and usage",            0 },

  { "--eval2score",  eslARG_NONE, FALSE, NULL, NULL,    NULL,  NULL,           NULL,            "compute score required to get E-value (E) for database of (Z) sequences",     0 },
  { "-Z",            eslARG_INT,    "1", NULL, "n>0",   NULL,  "--eval2score", NULL,            "database size, by default in # sequences , for --eval2score (default 1)",     0 },
  { "--rescntZ",    eslARG_NONE,   FALSE, NULL, NULL,   NULL,  "--eval2score", NULL,            "for --eval2score, -Z is in millions of residues (DNA models only)",          0 },
  { "-E",           eslARG_REAL,  "0.01", NULL, NULL,   NULL,  "--eval2score", NULL,            "E-value threshold, for --eval2score",                                         0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "display summary statistics for a profile file";


static int
output_header(FILE *ofp, const ESL_GETOPTS *go)
{
  p7_banner(ofp, go->argv[0], banner);

  if (esl_opt_IsUsed(go, "--eval2score"))  {
     if (  fprintf(ofp, "# show scores required to reach E-value:    %.2g\n",        esl_opt_GetReal(go, "-E"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
     if (esl_opt_IsUsed(go, "--rescntZ") ) {
       if (  fprintf(ofp, "# assume database residue count:            %d Kb\n",     esl_opt_GetInteger(go, "-Z")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
     } else {
       if (  fprintf(ofp, "# assume database sequence count:           %d\n",        esl_opt_GetInteger(go, "-Z")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
     }
  }

  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}



int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go	   = NULL;      /* command line processing                   */
  ESL_ALPHABET    *abc     = NULL;
  char            *hmmfile = NULL;
  P7_HMMFILE      *hfp     = NULL;
  P7_HMM          *hmm     = NULL;
  P7_BG           *bg      = NULL;
  int              nhmm;	
  double           x;
  float            KL;
  int              status;
  char             errbuf[eslERRBUFSIZE];

  int              do_eval2score = 0;
  int              z_val;
  float            e_val;

  /* Process the command line options.
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || 
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nOptions:");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 90); /* 0=docgroup, 2 = indentation; 80=textwidth*/


      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if ((hmmfile = esl_opt_GetArg(go, 1)) == NULL) 
    {
      puts("Failed to read <hmmfile> argument from command line.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  output_header(stdout, go);

  if ( esl_opt_IsOn(go, "--eval2score") ) {
    do_eval2score = TRUE;
    z_val         =  esl_opt_GetInteger(go, "-Z");
    e_val         =  esl_opt_GetReal(go, "-E");
  }

  /* Initializations: open the HMM file
   */
  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);  

  /* Main body: read HMMs one at a time, print one line of stats
   */
  printf("#\n");
  printf("# %-4s %-20s %-12s %8s %8s %6s %6s %6s %6s %6s", "idx",  "name",                 "accession",    "nseq",     "eff_nseq", "M",      "relent", "info",   "p relE", "compKL");
  if (do_eval2score)
    printf (" %6s %6.2g", "sc for", e_val);
  printf("\n");
  printf("# %-4s %-20s %-12s %8s %8s %6s %6s %6s %6s %6s", "----", "--------------------", "------------", "--------", "--------", "------", "------", "------", "------", "------");
  if (do_eval2score)
    printf (" %13s", "-------------");
  printf("\n");


  nhmm = 0;
  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslEOF) 
    {
      if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
      else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s",             hmmfile);
      else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets",   hmmfile);
      else if (status != eslOK)        esl_fatal("Unexpected error in reading HMMs from %s",   hmmfile);
      nhmm++;

      if ( esl_opt_IsOn(go, "--eval2score") ) {
        if (esl_opt_IsUsed(go, "--rescntZ") ) {
          if ( hmm->abc->type != eslRNA   && hmm->abc->type != eslDNA) {
            puts("The flag --rescntZ can't be used with non-nucleotide models.");
            esl_usage(stdout, argv[0], usage);
            printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
            exit(1);
          }
        } else if ( hmm->abc->type != eslAMINO  && hmm->abc->type != eslRNA && hmm->abc->type != eslDNA) {
          puts("The flag --eval2score can't be used with non-sequence models.");
          esl_usage(stdout, argv[0], usage);
          printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
          exit(1);
        }
      }


      if (bg == NULL) bg = p7_bg_Create(abc);

      p7_MeanPositionRelativeEntropy(hmm, bg, &x); 
      p7_hmm_CompositionKLDist(hmm, bg, &KL, NULL);

      printf("%-6d %-20s %-12s %8d %8.2f %6d %6.2f %6.2f %6.2f %6.2f",
	     nhmm,
	     hmm->name,
	     hmm->acc == NULL ? "-" : hmm->acc,
	     hmm->nseq,
	     hmm->eff_nseq,
	     hmm->M,
	     p7_MeanMatchRelativeEntropy(hmm, bg),
	     p7_MeanMatchInfo(hmm, bg),
	     x,
	     KL);


      if ( esl_opt_IsOn(go, "--eval2score") ) {
        float nseq;
        float sc;
        if (esl_opt_IsUsed(go, "--rescntZ") )
          nseq = (float)((long)z_val*1000000) / (float)(hmm->max_length);
        else
          nseq = (float)z_val;

        sc = esl_exp_invsurv( e_val / nseq ,  hmm->evparam[p7_FTAU],  hmm->evparam[p7_FLAMBDA]);

        printf("%8.1f", sc);

      }


      printf("\n");

	     /*	     p7_MeanForwardScore(hmm, bg)); */

      p7_hmm_Destroy(hmm);
    }

  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  p7_hmmfile_Close(hfp);
  esl_getopts_Destroy(go);
  exit(0);
}
