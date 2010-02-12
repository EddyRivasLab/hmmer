/* 
   gcc -O2 -I ~/src/hmmer/src   -L ~/src/hmmer/src\
           -I ~/src/hmmer/squid -L ~/src/hmmer/squid\
           -I ~/src/hmmer/easel -L ~/src/hmmer/easel -o speedbench\
    speedbench.c -lhmmer -lsquid -leasel -lm
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

static char banner[] = "\
speedbench :: various HMMER speed benchmarks";

static char usage[] = "Usage: speedbench [-options] <HMM>";

static ESL_OPTIONS options[] = {
   /* name          type        default   env   range  togs           reqs  incompat    help                          docgrouptag*/
  { "-h",         eslARG_NONE,   FALSE,   NULL, NULL,  NULL,          NULL, NULL,  "show help and usage",                      1},
  { "-n",         eslARG_INT,       "1",  NULL, "n>0", NULL,          NULL, NULL,  "set number of repetitions to <n>",         1},  
  { "-l",         eslARG_INT,    "350",   NULL, "n>0", NULL,          NULL, NULL,  "set simulated sequence length to <n>",     1},
  { "-F",         eslARG_NONE,   FALSE,   NULL, NULL,  "-C,-L,-S,-V", NULL, NULL,  "the Forward() algorithm",                  2},
  { "-S",         eslARG_NONE,   FALSE,   NULL, NULL,  "-C,-F,-L,-V", NULL, NULL,  "the SmallViterbi() algorithm",             2},
  { "-V",         eslARG_NONE,   FALSE,   NULL, NULL,  "-C,-F,-L,-S", NULL, NULL,  "the Viterbi() algorithm",                  2},
  { "-C",         eslARG_NONE,   FALSE,   NULL, NULL,  "-F,-L,-S,-V", NULL, NULL,  "model configuration, P7Config()",          2},
  { "-L",         eslARG_NONE,   FALSE,   NULL, NULL,  "-C,-F,-S,-V", NULL, NULL,  "target length config, P7ReconfigLength()", 2},
  { 0,0,0,0,0,0,0,0,0,0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go;		/* application configuration                  */
  ESL_STOPWATCH  *w;
  ESL_RANDOMNESS *r;		
  struct dpmatrix_s *mx;
  char           *hmmfile;	/* HMM file to read */
  HMMFILE        *hmmfp;
  struct plan7_s *hmm;
  int             i;
  double          ranp[MAXABET];
  char           *seq;
  unsigned char  *dsq;
  double          dpcells;

  int             show_help;	/* TRUE to show usage and exit    */
  int             ntrials;	/* number of trials to aggregate  */
  int             do_forward;	
  int             do_small;
  int             do_viterbi;
  int             do_config;
  int             do_length;
  int             L;

  /*****************************************************************
   * Parse the command line
   *****************************************************************/

  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);

  esl_opt_GetBooleanOption(go, "-h",         &show_help);
  esl_opt_GetIntegerOption(go, "-n",         &ntrials);
  esl_opt_GetIntegerOption(go, "-l",         &L);
  esl_opt_GetBooleanOption(go, "-F",         &do_forward);
  esl_opt_GetBooleanOption(go, "-S",         &do_small);
  esl_opt_GetBooleanOption(go, "-V",         &do_viterbi);
  esl_opt_GetBooleanOption(go, "-C",         &do_config);
  esl_opt_GetBooleanOption(go, "-L",         &do_length);

  if (show_help)
    {
      puts(usage);
      puts("\n  one of the following options must be set:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      puts("\n  and other available options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      exit(EXIT_SUCCESS);
    }

  if (esl_opt_ArgNumber(go) != 1) 
    esl_fatal("Incorrect number of command line arguments.\n%s\n", usage); 

  hmmfile = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL);

  /*****************************************************************
   * Get one HMM; synthesize one sequence.
   *****************************************************************/

  r = esl_randomness_CreateTimeseeded();

  if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  if (hmm == NULL) 
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  HMMFileClose(hmmfp);
  P7Config(hmm, P7_SW_MODE);  
  P7ReconfigLength(hmm, L);

  for (i = 0; i < Alphabet_size; i++) ranp[i] = hmm->null[i];
  seq = esl_rnd_IID(r, Alphabet, ranp, Alphabet_size, L);
  dsq = DigitizeSequence(seq, L);

  mx = CreatePlan7Matrix(L, hmm->M, 25, 0);

  w = esl_stopwatch_Create();

  /*****************************************************************
   * Run an alignment algorithm, however many times, while timing.
   *****************************************************************/

  esl_stopwatch_Start(w);
  for (i = 0; i < ntrials; i++)
    {
      if      (do_forward) P7Forward(dsq, L, hmm, NULL);
      else if (do_small)   P7SmallViterbi(dsq, L, hmm, mx, NULL);
      else if (do_viterbi) P7Viterbi(dsq, L, hmm, mx, NULL);
      else if (do_config)  P7Config(hmm, P7_SW_MODE); 
      else if (do_length)  P7ReconfigLength(hmm, 400);
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "CPU Time: ");

  dpcells = hmm->M * L * ntrials;
  
  if (do_forward || do_viterbi || do_small)
    printf("DP:      %.2f Mcells/sec\n", dpcells/(1e6*(w->user+w->sys)));
  else
    printf("config:  %.2g msec/call\n",  (w->user+w->sys)*1000./ntrials);

  /*****************************************************************
   * Fini.
   *****************************************************************/
  
  esl_randomness_Destroy(r);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);

  FreePlan7Matrix(mx);
  free(seq);
  free(dsq);
  FreePlan7(hmm);
  return 0;
}
  
