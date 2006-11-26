/* 
   gcc -O2 -I ~/src/hmmer/src   -L ~/src/hmmer/src\
           -I ~/src/hmmer/squid -L ~/src/hmmer/squid\
           -I ~/src/hmmer/easel -L ~/src/hmmer/easel -o iobench\
    iobench.c -lhmmer -lsquid -leasel -lm
*/

#include "config.h"		/* compile-time HMMER configuration constants */

#include <stdio.h>

/* Easel
 */
#include <easel.h>
#include <esl_getopts.h>
#include <esl_stopwatch.h>
#include <esl_random.h>

/* HMMER
 */
#include "plan7.h"		/* plan 7 profile HMM structure         */
#include "structs.h"
#include "funcs.h"
#include "globals.h"		/* alphabet global variables            */

static char banner[] = "iobench :: various HMMER i/o speed benchmarks";

static char usage[] = "Usage: iobench [-options] <HMM file>";

static ESL_OPTIONS options[] = {
   /* name          type        default   env   range  togs           reqs  incompat    help                          docgrouptag*/
  { "-h",         eslARG_NONE,   FALSE,   NULL, NULL,  NULL,          NULL, NULL,  "show help and usage",                      1},
  { "-n",         eslARG_INT,       "1",  NULL, "n>0", NULL,          NULL, NULL,  "set number of repetitions to <n>",         1},  
  { 0,0,0,0,0,0,0,0,0,0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go;		/* application configuration */
  ESL_STOPWATCH  *w;
  char           *hmmfile;	/* HMM file to read */
  HMMFILE        *hmmfp;
  struct plan7_s *hmm;
  int             i;
  int             Mtotal;
  int             nhmm;

  int             show_help;	/* TRUE to show usage and exit    */
  int             ntrials;	/* number of trials to aggregate  */

  /*****************************************************************
   * Parse the command line
   *****************************************************************/

  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);

  esl_opt_GetBooleanOption(go, "-h",         &show_help);
  esl_opt_GetIntegerOption(go, "-n",         &ntrials);

  if (show_help)
    {
      puts(usage);
      puts("\n  where available options are:");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 80);
      exit(EXIT_SUCCESS);
    }

  if (esl_opt_ArgNumber(go) != 1) 
    esl_fatal("Incorrect number of command line arguments.\n%s\n", usage); 

  hmmfile = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL);

  w = esl_stopwatch_Create();

  /*****************************************************************
   * Get HMMs, possibly many cycles, while timing.
   *****************************************************************/

  esl_stopwatch_Start(w);

  Mtotal = 0;
  nhmm = 0;
  for (i = 0; i < ntrials; i++)
    {
      if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
	Die("Failed to open HMM file %s\n%s", hmmfile, usage);

      while (HMMFileRead(hmmfp, &hmm)) 
	{
	  if (hmm == NULL) 
	    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
	  Mtotal += hmm->M;
	  nhmm++;

	  FreePlan7(hmm);
	}
      HMMFileClose(hmmfp);
    }
  
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "CPU Time: ");

  printf("total HMMs: %d\n", nhmm);
  printf("total M:    %d\n", Mtotal);


  /*****************************************************************
   * Fini.
   *****************************************************************/
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}
  
