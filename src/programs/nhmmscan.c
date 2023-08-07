/* nhmmscan: search sequence(s) against a profile HMM database, using nhmmer pipeline
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#ifdef HAVE_PTHREAD
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HAVE_PTHREAD*/

#include "hmmer.h"

//#ifdef P7_IMPL_DUMMY_INCLUDED
//#include "esl_vectorops.h"
//#define NHMMER_MAX_RESIDUE_COUNT (1024 * 100)
//#else
//#define NHMMER_MAX_RESIDUE_COUNT (10000000)  /* 10 million bases */
//#endif



typedef struct {
#ifdef HAVE_PTHREAD
  ESL_WORK_QUEUE   *queue;
#endif /*HAVE_PTHREAD*/
  ESL_SQ           *qsq;
  P7_BG            *bg;	         /* null model                              */
  P7_PIPELINE      *pli;         /* work pipeline                           */
  P7_TOPHITS       *th;          /* top hit results                         */
//  P7_OPROFILE      *om;          /* optimized query profile                 */
//  P7_SCOREDATA     *scoredata;   /* hmm-specific data used by nhmmer */
//  P7_HMM           *hmm;
} WORKER_INFO;

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#define CPUOPTS     NULL
#define MPIOPTS     NULL

//#define DAEMONOPTS  "-o,--tblout,--domtblout"

static ESL_OPTIONS options[] = {
  /* name           type          default  env  range toggles  reqs   incomp                         help                                           docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                          1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                         2 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-sequence hits to file <s>",         2 },
  { "--dfamtblout", eslARG_OUTFILE,      NULL, NULL, NULL,    NULL,  NULL,  NULL,              "save table of hits to file, in Dfam format <s>",               2 },
  { "--longinsertout", eslARG_OUTFILE,   NULL, NULL, NULL,    NULL,  NULL,  NULL,              "save table of long inserts to file <s>",                       2 },
  { "--insertlength", eslARG_INT,        "40", NULL, "n>=5",  NULL,  "--longinsertout",  NULL, "min length of insert printed to --longinsertout",              2 },
  { "--aliscoresout", eslARG_OUTFILE,   NULL, NULL, NULL,    NULL,  NULL,  NULL,               "save of scores for each position in each alignment to <s>",    2 },
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                        2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                 2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                          2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                      2 },
  /* Control of scoring system */
  { "--popen",      eslARG_REAL,       "0.02", NULL, "0<=x<0.5",NULL,  NULL,  NULL,          "gap open probability",                                         3 },
  { "--pextend",    eslARG_REAL,        "0.4", NULL, "0<=x<1",  NULL,  NULL,  NULL,          "gap extend probability",                                       3 },

  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report models <= this E-value threshold in output",             4 },
  { "-T",           eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report models >= this score threshold in output",               4 },
  /* Control of inclusion (significance) thresholds: */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider models <= this E-value threshold as significant",      5 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider models >= this score threshold as significant",        5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's GA gathering cutoffs to set all thresholding",    6 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's NC noise cutoffs to set all thresholding",        6 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's TC trusted cutoffs to set all thresholding",      6 },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",       7 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,    NULL,  NULL, "--max",          "MSV threshold: promote hits w/ P <= F1",                        7 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,    NULL,  NULL, "--max",          "Vit threshold: promote hits w/ P <= F2",                        7 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,    NULL,  NULL, "--max",          "Fwd threshold: promote hits w/ P <= F3",                        7 },
  { "--nobias",     eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--max",          "turn off composition bias filter",                              7 },
  { "--B1",         eslARG_INT,         "110", NULL, NULL,    NULL,  NULL, "--max,--nobias", "window length for biased-composition modifier (MSV)",          7 },
  { "--B2",         eslARG_INT,         "240", NULL, NULL,    NULL,  NULL, "--max,--nobias", "window length for biased-composition modifier (Vit)",          7 },
  { "--B3",         eslARG_INT,        "1000", NULL, NULL,    NULL,  NULL, "--max,--nobias", "window length for biased-composition modifier (Fwd)",          7 },

  /* Other options */
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,             "assert input <seqfile> is in format <s>: no autodetection",     12 },
  { "--nonull2",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  NULL,             "turn off biased composition score corrections",                12 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,             "set # of comparisons done, for E-value calculation",           12 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",  NULL,  NULL,  NULL,             "set RNG seed to <n> (if 0: one-time arbitrary seed)",          12 },
  { "--w_beta",     eslARG_REAL,    NULL, NULL, NULL,    NULL,  NULL,           NULL,    "tail mass at which window length is determined",                12 },
  { "--w_length",   eslARG_INT,     NULL, NULL, NULL,    NULL,  NULL,           NULL,    "window length - essentially max expected hit length ",                                                12 },
  { "--toponly",     eslARG_NONE,   NULL, NULL, NULL,    NULL,  NULL,   "--bottomonly",  "only search the top strand",                                    12 },
  { "--bottomonly",  eslARG_NONE,   NULL, NULL, NULL,    NULL,  NULL,      "--toponly",  "only search the bottom strand",                                    12 },



  /* Not used, but retained because esl option-handling code errors if it isn't kept here.  Placed in group 99 so it doesn't print to help*/
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  DOMREPOPTS,      "report domains <= this E-value threshold in output",            99 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  DOMREPOPTS,      "report domains >= this score cutoff in output",                 99 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCDOMOPTS,      "consider domains <= this E-value threshold as significant",     99 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCDOMOPTS,      "consider domains >= this score threshold as significant",       99 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,           "set # of significant seqs, for domain E-value calculation",      99 },

//  { "--daemon",     eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  DAEMONOPTS,      "run program as a daemon",                                      12 },


#ifdef HAVE_PTHREAD
  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>=0",NULL,  NULL,  CPUOPTS,         "number of parallel CPU workers to use for multithreads",       12 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char            *seqfile;           /* query sequence file                             */
  char            *hmmfile;           /* database HMM file                               */

  int              do_mpi;            /* TRUE if we're doing MPI parallelization         */
  int              nproc;             /* how many MPI processes, total                   */
  int              my_rank;           /* who am I, in 0..nproc-1                         */
};

static char usage[]  = "[-options] <hmmdb> <seqfile>";
static char banner[] = "search DNA sequence(s) against a DNA profile database";

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, P7_HMMFILE *hfp);
#ifdef HAVE_PTHREAD
#define BLOCK_SIZE 1
static int  thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, P7_HMMFILE *hfp);
static void pipeline_thread(void *arg);
#endif /*HAVE_PTHREAD*/


/* process_commandline()
 * 
 * Processes the commandline, filling in fields in <cfg> and creating and returning
 * an <ESL_GETOPTS> options structure. The help page (hmmsearch -h) is formatted
 * here.
 */
static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_hmmfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      if (puts("\nBasic options:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

      if (puts("\nOptions controlling output:")                              < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      if (puts("\nOptions controlling reporting thresholds:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      if (puts("\nOptions controlling inclusion (significance) thresholds:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      if (puts("\nOptions for model-specific thresholding:")                 < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 

      if (puts("\nOptions controlling acceleration heuristics:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      if (puts("\nOther expert options:")                                    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)      { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { if (puts("Failed to get <hmmdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_hmmfile, "-") == 0) 
    { if (puts("nhmmscan cannot read <hmm database> from stdin stream, because it must have hmmpress'ed auxfiles") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");   goto FAILURE;  }

  *ret_go = go;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  if (puts("\nwhere most common options are:")                                 < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}


static int
output_header(FILE *ofp, ESL_GETOPTS *go, char *hmmfile, char *seqfile)
{
  p7_banner(ofp, go->argv[0], banner);
  
  if (fprintf(ofp, "# query sequence file:             %s\n", seqfile)                                                                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# target HMM database:             %s\n", hmmfile)                                                                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-o")          && fprintf(ofp, "# output directed to file:         %s\n",            esl_opt_GetString(go, "-o"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--tblout")    && fprintf(ofp, "# per-seq hits tabular output:     %s\n",            esl_opt_GetString(go, "--tblout"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--dfamtblout")    && fprintf(ofp, "# hits output in Dfam format:      %s\n",            esl_opt_GetString(go, "--dfamtblout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--longinsertout") && fprintf(ofp, "# long inserts output:             %s\n",            esl_opt_GetString(go, "--longinsertout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--insertlength")  && fprintf(ofp, "# long inserts minimum length:     %d\n",            esl_opt_GetInteger(go, "--insertlength")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--aliscoresout")  && fprintf(ofp, "# alignment scores output:         %s\n",            esl_opt_GetString(go, "--longinsertout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--acc")       && fprintf(ofp, "# prefer accessions over names:    yes\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--noali")     && fprintf(ofp, "# show alignments in output:       no\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notextw")   && fprintf(ofp, "# max ASCII text line length:      unlimited\n")                                           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--textw")     && fprintf(ofp, "# max ASCII text line length:      %d\n",            esl_opt_GetInteger(go, "--textw"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--popen")      && fprintf(ofp, "# gap open probability:            %f\n",             esl_opt_GetReal   (go, "--popen"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--pextend")    && fprintf(ofp, "# gap extend probability:          %f\n",             esl_opt_GetReal   (go, "--pextend"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "-E")          && fprintf(ofp, "# profile reporting threshold:     E-value <= %g\n", esl_opt_GetReal(go, "-E"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-T")          && fprintf(ofp, "# profile reporting threshold:     score >= %g\n",   esl_opt_GetReal(go, "-T"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incE")      && fprintf(ofp, "# profile inclusion threshold:     E-value <= %g\n", esl_opt_GetReal(go, "--incE"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incT")      && fprintf(ofp, "# profile inclusion threshold:     score >= %g\n",   esl_opt_GetReal(go, "--incT"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_ga")    && fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_nc")    && fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_tc")    && fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--max")       && fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n")                      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F1")        && fprintf(ofp, "# MSV filter P threshold:       <= %g\n",            esl_opt_GetReal(go, "--F1"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F2")        && fprintf(ofp, "# Vit filter P threshold:       <= %g\n",            esl_opt_GetReal(go, "--F2"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F3")        && fprintf(ofp, "# Fwd filter P threshold:       <= %g\n",            esl_opt_GetReal(go, "--F3"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nobias")    && fprintf(ofp, "# biased composition HMM filter:   off\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--B1")         && fprintf(ofp, "# biased comp MSV window len:      %d\n",             esl_opt_GetInteger(go, "--B1"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--B2")         && fprintf(ofp, "# biased comp Viterbi window len:  %d\n",             esl_opt_GetInteger(go, "--B2"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--B3")         && fprintf(ofp, "# biased comp Forward window len:  %d\n",             esl_opt_GetInteger(go, "--B3"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");


  if (esl_opt_IsUsed(go, "--nonull2")   && fprintf(ofp, "# null2 bias corrections:          off\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--toponly")    && fprintf(ofp, "# search only top strand:          on\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--bottomonly") && fprintf(ofp, "# search only bottom strand:       on\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "-Z")          && fprintf(ofp, "# sequence search space set to:    %.0f\n",          esl_opt_GetReal(go, "-Z"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed")==0 && fprintf(ofp, "# random number seed:              one-time arbitrary\n")                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if (                                  fprintf(ofp, "# random number seed set to:       %d\n",        esl_opt_GetInteger(go, "--seed"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--qformat")   && fprintf(ofp, "# input seqfile format asserted:   %s\n",            esl_opt_GetString(go, "--qformat"))   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--w_beta")     && fprintf(ofp, "# window length beta value:        %g\n",             esl_opt_GetReal(go, "--w_beta"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--w_length")   && fprintf(ofp, "# window length :                  %d\n",             esl_opt_GetInteger(go, "--w_length")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
//  if (esl_opt_IsUsed(go, "--daemon")    && fprintf(ofp, "run as a daemon process\n")                                                                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#ifdef HAVE_PTHREAD
  if (esl_opt_IsUsed(go, "--cpu")       && fprintf(ofp, "# number of worker threads:        %d\n",            esl_opt_GetInteger(go, "--cpu"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");  
#endif
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go  = NULL;	
  struct cfg_s     cfg;         
  int              status   = eslOK;

  p7_Init();			

  /* Initialize what we can in the config structure (without knowing the alphabet yet) */
  cfg.hmmfile    = NULL;
  cfg.seqfile    = NULL;
  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */

  process_commandline(argc, argv, &go, &cfg.hmmfile, &cfg.seqfile);    

  status = serial_master(go, &cfg);

  esl_getopts_Destroy(go);
  return status;
}


/* serial_master()
 * The serial version of hmmsearch.
 * For each query HMM in <hmmdb> search the database for hits.
 * 
 * A master can only return if it's successful. All errors are handled
 * immediately and fatally with p7_Fail().
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;	         /* output file for results (default stdout)        */
  FILE            *tblfp    = NULL;		 /* output stream for tabular per-seq (--tblout)    */
  FILE            *dfamtblfp    = NULL;            /* output stream for tabular Dfam format (--dfamtblout)    */
  FILE            *inserttblfp  = NULL;            /* output stream for table of long inserts (--longinsertout)    */
  FILE            *aliscoresfp  = NULL;            /* output stream for alignment scores (--aliscoresout)    */

//  P7_HMM          *hmm        = NULL;              /* one HMM query                                   */
//  P7_SCOREDATA    *scoredata  = NULL;

  int              seqfmt   = eslSQFILE_UNKNOWN; /* format of seqfile                               */
  ESL_SQFILE      *sqfp     = NULL;              /* open seqfile                                    */
  P7_HMMFILE      *hfp      = NULL;		 /* open HMM database file                          */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                               */
  P7_OPROFILE     *om       = NULL;		 /* target profile                                  */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                          */
  ESL_SQ          *qsq      = NULL;		 /* query sequence                                  */
  int              nquery   = 0;
  int              textw;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;

  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#ifdef HAVE_PTHREAD
  P7_OM_BLOCK     *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char             errbuf[eslERRBUFSIZE];

  double window_beta = -1.0 ;
  int window_length  = -1;
  if (esl_opt_IsUsed(go, "--w_beta")) { if (  ( window_beta   = esl_opt_GetReal(go, "--w_beta") )  < 0 || window_beta > 1  ) esl_fatal("Invalid window-length beta value\n"); }
  if (esl_opt_IsUsed(go, "--w_length")) { if (( window_length = esl_opt_GetInteger(go, "--w_length")) < 4  ) esl_fatal("Invalid window length value\n"); }


  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* If caller declared an input format, decode it */
  if (esl_opt_IsOn(go, "--qformat")) {
    seqfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (seqfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }

  /* validate options if running as a daemon */
//  if (esl_opt_IsOn(go, "--daemon")) {
    /* running as a daemon, the input format must be type daemon */
//    if (seqfmt != eslSQFILE_UNKNOWN && seqfmt != eslSQFILE_DAEMON)
//      esl_fatal("Input format %s not supported.  Must be daemon\n", esl_opt_GetString(go, "--qformat"));
//    seqfmt = eslSQFILE_DAEMON;

//    if (strcmp(cfg->seqfile, "-") != 0) esl_fatal("Query sequence file must be '-'\n");
//  }

  /* Open the target profile database to get the sequence alphabet */
  status = p7_hmmfile_OpenE(cfg->hmmfile, p7_HMMDBENV, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem, trying to open HMM file %s.\n%s\n",                  cfg->hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  
  if (! hfp->is_pressed)           p7_Fail("Failed to open binary auxfiles for %s: use hmmpress first\n",             hfp->fname);

  hstatus = p7_oprofile_ReadSSV(hfp, &abc, &om);
  if      (hstatus == eslEFORMAT)   p7_Fail("bad format, binary auxfiles, %s:\n%s",     cfg->hmmfile, hfp->errbuf);
  else if (hstatus == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets", cfg->hmmfile);
  else if (hstatus != eslOK)        p7_Fail("Unexpected error in reading HMMs from %s", cfg->hmmfile); 

  p7_oprofile_Destroy(om);
  p7_hmmfile_Close(hfp);

  /* Open the query sequence database */
  status = esl_sqfile_OpenDigital(abc, cfg->seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);
  if (sqfp->format > 100) // breaking the law!  That range is reserved for msa, for aligned formats
    p7_Fail("%s contains a multiple sequence alignment; expect unaligned sequences, like FASTA\n",   cfg->seqfile);
  qsq = esl_sq_CreateDigital(abc);


  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"),          "w")) == NULL)  esl_fatal("Failed to open output file %s for writing\n",                 esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp")); }
  if (esl_opt_IsOn(go, "--dfamtblout"))    { if ((dfamtblfp    = fopen(esl_opt_GetString(go, "--dfamtblout"),"w"))    == NULL)  esl_fatal("Failed to open tabular dfam output file %s for writing\n", esl_opt_GetString(go, "--dfamtblout")); }
  if (esl_opt_IsOn(go, "--longinsertout")) { if ((inserttblfp  = fopen(esl_opt_GetString(go, "--longinsertout"),"w")) == NULL)  esl_fatal("Failed to open long-insert output file %s for writing\n", esl_opt_GetString(go, "--longinsertout")); }
  if (esl_opt_IsOn(go, "--aliscoresout"))  { if ((aliscoresfp  = fopen(esl_opt_GetString(go, "--aliscoresout"),"w")) == NULL)  esl_fatal("Failed to open alignment scores output file %s for writing\n", esl_opt_GetString(go, "--aliscoresout")); }
 
  output_header(ofp, go, cfg->hmmfile, cfg->seqfile);

#ifdef HAVE_PTHREAD
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                           esl_threads_CPUCount(&ncpus);

  if (ncpus > 0)
    {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
    }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  for (i = 0; i < infocnt; ++i)
    {
      info[i].bg    = p7_bg_Create(abc);
#ifdef HAVE_PTHREAD
      info[i].queue = queue;
#endif
    }

#ifdef HAVE_PTHREAD
  for (i = 0; i < ncpus * 2; ++i)
    {
      block = p7_oprofile_CreateBlock(BLOCK_SIZE);
      if (block == NULL)    esl_fatal("Failed to allocate sequence block");

      status = esl_workqueue_Init(queue, block);
      if (status != eslOK)  esl_fatal("Failed to add block to work queue");
    }
#endif

  /* Outside loop: over each query sequence in <seqfile>. */
  while ((sstatus = esl_sqio_Read(sqfp, qsq)) == eslOK)
  {
      if (sstatus == eslEMEM)                 p7_Fail("Memory allocation error reading sequence file\n", status);
      if (sstatus == eslEINCONCEIVABLE)       p7_Fail("Unexpected error %d reading sequence file\n", status);
     // if (qsq->L > NHMMER_MAX_RESIDUE_COUNT)  p7_Fail("Input sequence %s in file %s exceeds maximum length of %d bases.\n",  qsq->name, cfg->seqfile, NHMMER_MAX_RESIDUE_COUNT);

      nquery++;
      esl_stopwatch_Start(w);	                          

      /* Open the target profile database */
      status = p7_hmmfile_OpenE(cfg->hmmfile, p7_HMMDBENV, &hfp, NULL);
      if (status != eslOK)        p7_Fail("Unexpected error %d in opening hmm file %s.\n",           status, cfg->hmmfile);  
  
#ifdef HAVE_PTHREAD
      /* if we are threaded, create a lock to prevent multiple readers */
      if (ncpus > 0)
      {
        status = p7_hmmfile_CreateLock(hfp);
        if (status != eslOK) p7_Fail("Unexpected error %d creating lock\n", status);
      }
#endif

      if (fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (qsq->acc[0]  != 0 && fprintf(ofp, "Accession:   %s\n", qsq->acc)     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (qsq->desc[0] != 0 && fprintf(ofp, "Description: %s\n", qsq->desc)    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      for (i = 0; i < infocnt; ++i)
      {
        /* Create processing pipeline and hit list */
        info[i].th  = p7_tophits_Create(p7_TOPHITS_DEFAULT_INIT_ALLOC);
        info[i].pli = p7_pipeline_Create(go, 100, 100, TRUE, p7_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
        info[i].pli->hfp = hfp;  /* for two-stage input, pipeline needs <hfp> */


        p7_pipeline_NewSeq(info[i].pli, qsq);
        info[i].qsq = qsq;

        if (  esl_opt_IsUsed(go, "--toponly") )
          info[i].pli->strand = p7_STRAND_TOPONLY;
        else if (  esl_opt_IsUsed(go, "--bottomonly") )
          info[i].pli->strand = p7_STRAND_BOTTOMONLY;
        else
          info[i].pli->strand = p7_STRAND_BOTH;


#ifdef HAVE_PTHREAD
          if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
      }

#ifdef HAVE_PTHREAD
      if (ncpus > 0)  hstatus = thread_loop(threadObj, queue, hfp);
      else	      hstatus = serial_loop(info, hfp);
#else
      hstatus = serial_loop(info, hfp);
#endif
      switch(hstatus)
      {
        case eslEFORMAT:   p7_Fail("bad file format in HMM file %s",             cfg->hmmfile);	  break;
        case eslEINCOMPAT: p7_Fail("HMM file %s contains different alphabets",   cfg->hmmfile);	  break;
        case eslEOF:
        case eslOK:   /* do nothing */
          break;
        default: 	   p7_Fail("Unexpected error in reading HMMs from %s",   cfg->hmmfile);
      }



      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i)
      {
        p7_tophits_Merge(info[0].th, info[i].th);
        p7_pipeline_stats_Merge(info[0].pli, &(info[i].pli.stats));

        p7_pipeline_Destroy(info[i].pli);
        p7_tophits_Destroy(info[i].th);
      }


      /* modify e-value to account for number of models */
      for (i = 0; i < info->th->N ; i++)
      {
        info->th->unsrt[i].lnP         += log((float)info->pli->stats.nmodels);
        info->th->unsrt[i].dcl[0].lnP   = info->th->unsrt[i].lnP;
        info->th->unsrt[i].sortkey      = -1.0 * info->th->unsrt[i].lnP;
      }


      /* it's possible to have duplicates based on how viterbi ranges can overlap */
      p7_tophits_SortByModelnameAndAlipos(info->th);
      p7_tophits_RemoveDuplicates(info->th, info->pli->use_bit_cutoffs);

      /* Print results */
      p7_tophits_SortBySortkey(info->th);
      p7_tophits_Threshold(info->th, info->pli);

      //tally up total number of hits and target coverage
      info->pli->stats.n_output = info->pli->stats.pos_output = 0;
      for (i = 0; i < info->th->N; i++) {
          if ( (info->th->hit[i]->flags & p7_IS_REPORTED) || info->th->hit[i]->flags & p7_IS_INCLUDED) {
              info->pli->stats.n_output++;
              info->pli->stats.pos_output += abs(info->th->hit[i]->dcl[0].jali - info->th->hit[i]->dcl[0].iali) + 1;
          }
      }




      p7_tophits_Targets(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      p7_tophits_Domains(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      if (tblfp)     p7_tophits_TabularTargets(tblfp,    qsq->name, qsq->acc, info->th, info->pli, (nquery == 1));
      if (dfamtblfp) p7_tophits_TabularXfam(dfamtblfp,   qsq->name, NULL, info->th, info->pli);
      if (inserttblfp) p7_tophits_LongInserts(inserttblfp, qsq->name, qsq->acc, info->th, info->pli, esl_opt_GetInteger(go, "--insertlength") );
      if (aliscoresfp) p7_tophits_AliScores(aliscoresfp, qsq->name, info->th );


      esl_stopwatch_Stop(w);
      info->pli->stats.nseqs = 1;
      p7_pipeline_WriteStats(ofp, info->pli, w);
      if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      fflush(ofp);

      p7_hmmfile_Close(hfp);
      p7_pipeline_Destroy(info->pli);
      p7_tophits_Destroy(info->th);
      esl_sq_Reuse(qsq);
  }



  if      (sstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					    sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (sstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    sstatus, sqfp->filename);

  /* Terminate outputs - any last words?
   */
  if (tblfp)    p7_tophits_TabularTail(tblfp,    "hmmscan", p7_SCAN_MODELS, cfg->seqfile, cfg->hmmfile, go);
  if (ofp)      { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for successful exit
   */
  for (i = 0; i < infocnt; ++i)
    p7_bg_Destroy(info[i].bg);

#ifdef HAVE_PTHREAD
  if (ncpus > 0)
    {
      esl_workqueue_Reset(queue);
      while (esl_workqueue_Remove(queue, (void **) &block) == eslOK)
        p7_oprofile_DestroyBlock(block);
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
    }
#endif

  free(info);

  esl_sq_Destroy(qsq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);

  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);
  if (dfamtblfp)     fclose(dfamtblfp);
  if (inserttblfp)   fclose(inserttblfp);
  if (aliscoresfp)   fclose(aliscoresfp);

  return eslOK;

 ERROR:
  return status;
}


static int
serial_loop(WORKER_INFO *info, P7_HMMFILE *hfp)
{
  int            status;
  int i;
  int seq_len = 0;
  int prev_hit_cnt = 0;
  P7_OPROFILE   *om        = NULL;
  P7_SCOREDATA  *scoredata = NULL;   /* hmm-specific data used by nhmmer */
  ESL_ALPHABET  *abc = NULL;
  P7_DOMAIN *dcl;

#ifdef eslAUGMENT_ALPHABET
  ESL_SQ        *sq_revcmp = NULL;
  if (info->pli->strand != p7_STRAND_TOPONLY && info->qsq->abc->complement != NULL ) {
    sq_revcmp =  esl_sq_CreateDigital(info->qsq->abc);
    esl_sq_Copy(info->qsq,sq_revcmp);
    esl_sq_ReverseComplement(sq_revcmp);

    info->pli->stats.nres += info->qsq->n;
  }
#endif /*eslAUGMENT_ALPHABET*/


  /* Main loop: */
  while ((status = p7_oprofile_ReadSSV(hfp, &abc, &om)) == eslOK)
  {
      seq_len = 0;

      p7_pipeline_NewModel(info->pli, om, info->bg);
      p7_bg_SetLength(info->bg, info->qsq->n);
      p7_oprofile_ReconfigLength(om, info->qsq->n);

      scoredata = p7_hmm_ScoreDataCreate(om, FALSE);

#ifdef eslAUGMENT_ALPHABET
      //reverse complement
      if (info->pli->strand != p7_STRAND_TOPONLY && info->qsq->abc->complement != NULL )
      {

        p7_Pipeline_LongTarget(info->pli, om, scoredata, info->bg, sq_revcmp, info->th, 0);
        p7_pipeline_Reuse(info->pli); // prepare for next search
        seq_len = info->qsq->n;
        for (i = prev_hit_cnt; i < info->th->N ; i++)
        {
          dcl = info->th->unsrt[i].dcl;
          // modify hit positions to account for the position of the window in the full sequence
          dcl->ienv = seq_len - dcl->ienv + 1;
          dcl->jenv = seq_len - dcl->jenv + 1;
          dcl->iali = seq_len - dcl->iali + 1;
          dcl->jali = seq_len - dcl->jali + 1;
          dcl->ad->sqfrom = seq_len - dcl->ad->sqfrom + 1;
          dcl->ad->sqto = seq_len - dcl->ad->sqto + 1;
        }

      }
#endif


      if (info->pli->strand != p7_STRAND_BOTTOMONLY) {
        p7_Pipeline_LongTarget(info->pli, om, scoredata, info->bg, info->qsq, info->th, 0);
        p7_pipeline_Reuse(info->pli);
        seq_len += info->qsq->n;
      }

      for (i = prev_hit_cnt; i < info->th->N ; i++)
      {
        info->th->unsrt[i].lnP         += log((float)seq_len / (float)om->max_length);
        info->th->unsrt[i].dcl[0].lnP   = info->th->unsrt[i].lnP;
        info->th->unsrt[i].sortkey      = -1.0 * info->th->unsrt[i].lnP;
        info->th->unsrt[i].dcl[0].ad->L =  om->M;
      }

      prev_hit_cnt = info->th->N;

      p7_oprofile_Destroy(om);
      p7_hmm_ScoreDataDestroy(scoredata);

  }

  esl_alphabet_Destroy(abc);
#ifdef eslAUGMENT_ALPHABET
  esl_sq_Destroy(sq_revcmp);
#endif
  return status;
}

#ifdef HAVE_PTHREAD
static int
thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, P7_HMMFILE *hfp)
{
  int  status   = eslOK;
  int  sstatus  = eslOK;
  int  eofCount = 0;
  P7_OM_BLOCK   *block;
  ESL_ALPHABET  *abc = NULL;
  void          *newBlock;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
      
  /* Main loop: */
  while (sstatus == eslOK)
  {
      block = (P7_OM_BLOCK *) newBlock;
      sstatus = p7_oprofile_ReadBlockSSV(hfp, &abc, block);

      if (sstatus == eslEOF)
      {
        if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
        ++eofCount;
      }
	  
      if (sstatus == eslOK)
      {
        status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
        if (status != eslOK) esl_fatal("Work queue reader failed");
      }
  }

  status = esl_workqueue_ReaderUpdate(queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF)
  {
      /* wait for all the threads to complete */
      esl_threads_WaitForFinish(obj);
      esl_workqueue_Complete(queue);  
  }
  
  esl_alphabet_Destroy(abc);
  return sstatus;
}

static void 
pipeline_thread(void *arg)
{
  int i, j;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  P7_OM_BLOCK   *block;
  void          *newBlock;
  P7_OPROFILE   *om        = NULL;
  P7_SCOREDATA  *scoredata = NULL;   /* hmm-specific data used by nhmmer */

  P7_DOMAIN *dcl;
  int seq_len = 0;
  int prev_hit_cnt = 0;

#ifdef eslAUGMENT_ALPHABET
  ESL_SQ        *sq_revcmp = NULL;
#endif /*eslAUGMENT_ALPHABET*/
  
  p7_Init();

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");

#ifdef eslAUGMENT_ALPHABET
  //reverse complement
  if (info->pli->strand != p7_STRAND_TOPONLY && info->qsq->abc->complement != NULL ) {
    sq_revcmp =  esl_sq_CreateDigital(info->qsq->abc);
    esl_sq_Copy(info->qsq,sq_revcmp);
    esl_sq_ReverseComplement(sq_revcmp);
    info->pli->stats.nres += info->qsq->n;
  }
#endif /*eslAUGMENT_ALPHABET*/


  /* loop until all blocks have been processed */
  block = (P7_OM_BLOCK *) newBlock;
  while (block->count > 0)
  {
      /* Main loop: */
      for (i = 0; i < block->count; ++i)
      {
        om = block->list[i];
        seq_len = 0;

        p7_pipeline_NewModel(info->pli, om, info->bg);
        p7_bg_SetLength(info->bg, info->qsq->n);
        p7_oprofile_ReconfigLength(om, info->qsq->n);

        scoredata = p7_hmm_ScoreDataCreate(om, FALSE);


#ifdef eslAUGMENT_ALPHABET
        //reverse complement
        if (info->pli->strand != p7_STRAND_TOPONLY && info->qsq->abc->complement != NULL )
        {
          p7_Pipeline_LongTarget(info->pli, om, scoredata, info->bg, sq_revcmp, info->th, 0);
          p7_pipeline_Reuse(info->pli); // prepare for next search

          seq_len = info->qsq->n;
          for (j = prev_hit_cnt; j < info->th->N ; j++)
          {
            dcl = info->th->unsrt[j].dcl;
            // modify hit positions to account for the position of the window in the full sequence
            dcl->ienv = seq_len - dcl->ienv + 1;
            dcl->jenv = seq_len - dcl->jenv + 1;
            dcl->iali = seq_len - dcl->iali + 1;
            dcl->jali = seq_len - dcl->jali + 1;
            dcl->ad->sqfrom = seq_len - dcl->ad->sqfrom + 1;
            dcl->ad->sqto = seq_len - dcl->ad->sqto + 1;
          }

        }
#endif
        if (info->pli->strand != p7_STRAND_BOTTOMONLY) {
          p7_Pipeline_LongTarget(info->pli, om, scoredata, info->bg, info->qsq, info->th, 0);
          p7_pipeline_Reuse(info->pli);
          seq_len += info->qsq->n;
        }

        for (j = prev_hit_cnt; j < info->th->N ; j++)
        {
          info->th->unsrt[j].lnP         += log((float)seq_len / (float)om->max_length);
          info->th->unsrt[j].dcl[0].lnP   = info->th->unsrt[j].lnP;
          info->th->unsrt[j].sortkey      = -1.0 * info->th->unsrt[j].lnP;
          info->th->unsrt[j].dcl[0].ad->L = om->M;
        }

        prev_hit_cnt = info->th->N;
        p7_hmm_ScoreDataDestroy(scoredata);
        p7_oprofile_Destroy(om);
        block->list[i] = NULL;
      }


      status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      block = (P7_OM_BLOCK *) newBlock;
  }

#ifdef eslAUGMENT_ALPHABET
  esl_sq_Destroy(sq_revcmp);
#endif


  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
  return;
}
#endif   /* HAVE_PTHREAD */


