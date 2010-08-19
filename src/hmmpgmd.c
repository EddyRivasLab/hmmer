/* hmmpgmd: hmmer deamon searchs against a sequence database.
 * 
 * MSF, Thu Aug 12, 2010 [Janelia]
 * SVN $Id: hmmsearch.c 3324 2010-07-07 19:30:12Z wheelert $
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "esl_mpi.h"
#endif /*HAVE_MPI*/

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"

#define MAX_BUFFER (4*1024)

typedef struct {
  char             *pgm;         /* program name               */
  
  P7_HMM           *hmm;         /* query HMM                  */
  ESL_SQ           *seq;         /* query sequence             */
  ESL_ALPHABET     *abc;         /* digital alphabet           */
  ESL_GETOPTS      *opts;        /* search specific options    */

  int               buf_size;
  char             *buffer;
  char             *opts_ptr;
  char             *query_ptr;
} QUERY_INFO;

static int read_QueryInfo(QUERY_INFO *queryInfo, FILE *qfp);
static int process_QueryInfo(QUERY_INFO *info);
static int destroy_QueryInfo(QUERY_INFO *info);
static int reuse_QueryInfo(QUERY_INFO *info);

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  P7_HMM           *hmm;         /* query HMM                               */
  ESL_SQ           *seq;         /* query sequence                          */
  ESL_ALPHABET     *abc;         /* digital alphabet                        */
  ESL_GETOPTS      *opts;
  /* Structure created and populated by the individual threads.
   * The main thread is responsible for freeing up the memory.
   */
  P7_PIPELINE      *pli;         /* work pipeline                           */
  P7_TOPHITS       *th;          /* top hit results                         */
} WORKER_INFO;

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#ifdef HAVE_MPI
#define DAEMONOPTS  "--mpi,--stall"
#else
#define DAEMONOPTS  NULL
#endif

#if defined (HMMER_THREADS) && defined (HAVE_MPI)
#define CPUOPTS     "--mpi"
#define MPIOPTS     "--cpu"
#else
#define CPUOPTS     NULL
#define MPIOPTS     NULL
#endif

static ESL_OPTIONS cmdlineOpts[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                         1 },
  /* Control of output */
  //m/{ "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                        2 },
  //m/{ "-A",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save multiple alignment of all hits to file <s>",              2 },
  //m/{ "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-sequence hits to file <s>",        2 },
  //m/{ "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-domain hits to file <s>",          2 },
  /* Other options */
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert target <seqfile> is in format <s>: no autodetection",  12 },
  { "--daemon",     eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  DAEMONOPTS,      "run program as a daemon",                                     12 },

#ifdef HMMER_THREADS 
  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>=0",NULL,  NULL,  CPUOPTS,         "number of parallel CPU workers to use for multithreads",      12 },
#endif
#ifdef HAVE_MPI
  { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--mpi", NULL,            "arrest after start: for debugging MPI under gdb",             12 },  
  { "--mpi",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  MPIOPTS,         "run as an MPI parallel program",                              12 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static ESL_OPTIONS searchOpts[] = {
  /* Control of output */
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                       2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                         2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                     2 },
  /* Control of scoring system */
  { "--popen",      eslARG_REAL,  "0.02", NULL, "0<=x<0.5",NULL,  NULL,  NULL,          "gap open probability",                                         3 },
  { "--pextend",    eslARG_REAL,   "0.4", NULL, "0<=x<1",  NULL,  NULL,  NULL,          "gap extend probability",                                       3 },
  { "--mxfile",     eslARG_INFILE,  NULL, NULL, NULL,      NULL,  NULL,  NULL,          "substitution score matrix [default: BLOSUM62]",                3 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report sequences <= this E-value threshold in output",         4 },
  { "-T",           eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report sequences >= this score threshold in output",           4 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  DOMREPOPTS,      "report domains <= this E-value threshold in output",           4 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  DOMREPOPTS,      "report domains >= this score cutoff in output",                4 },
  /* Control of inclusion (significance) thresholds */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider sequences >= this score threshold as significant",    5 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCDOMOPTS,      "consider domains <= this E-value threshold as significant",    5 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCDOMOPTS,      "consider domains >= this score threshold as significant",      5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's GA gathering cutoffs to set all thresholding",   6 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's NC noise cutoffs to set all thresholding",       6 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's TC trusted cutoffs to set all thresholding",     6 },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",      7 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,    NULL,  NULL, "--max",          "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,    NULL,  NULL, "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,    NULL,  NULL, "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--nobias",     eslARG_NONE,   NULL,  NULL, NULL,    NULL,  NULL, "--max",          "turn off composition bias filter",                             7 },
  /* Control of E-value calibration */
  { "--EmL",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,          "length of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EmN",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,          "number of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EvL",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,          "length of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EvN",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,          "number of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EfL",        eslARG_INT,    "100", NULL,"n>0",      NULL,  NULL,  NULL,          "length of sequences for Forward exp tail tau fit",            11 },   
  { "--EfN",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,          "number of sequences for Forward exp tail tau fit",            11 },   
  { "--Eft",        eslARG_REAL,  "0.04", NULL,"0<x<1",    NULL,  NULL,  NULL,          "tail mass for Forward exponential tail tau fit",              11 },   
  /* Other options */
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,    NULL,  NULL,  NULL,            "turn off biased composition score corrections",               12 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of comparisons done, for E-value calculation",          12 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of significant seqs, for domain E-value calculation",   12 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[options] <query file> <target seqfile>";
static char banner[] = "search against a sequence database";

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char            *dbfile;            /* target sequence database file                   */
  char            *hmmfile;           /* query HMM file                                  */

  int              do_mpi;            /* TRUE if we're doing MPI parallelization         */
  int              nproc;             /* how many MPI processes, total                   */
  int              my_rank;           /* who am I, in 0..nproc-1                         */
};

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ESL_SQFILE *dbfp);
#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/

#ifdef HAVE_MPI
static int  mpi_master   (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  mpi_worker   (ESL_GETOPTS *go, struct cfg_s *cfg);
#endif /*HAVE_MPI*/


static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_seqfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(cmdlineOpts)) == NULL)     p7_Die("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n",  go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n",  go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);

      puts("\nBasic options:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

      puts("\nOther expert options:");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                  != 1)     { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }
  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere most common options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}

static void
process_searchline(char *pgm, char *cmdstr, ESL_GETOPTS *go)
{
  if (esl_getopts_Reuse(go)            != eslOK)    p7_Die("Internal failure reusing options object");
  if (esl_opt_ProcessSpoof(go, cmdstr) != eslOK)  { printf("Failed to parse options string: %s\n",  go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)         != eslOK)  { printf("Failed to parse options string: %s\n",  go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    p7_banner(stdout, pgm, banner);
    esl_usage(stdout, pgm, usage);
    puts("\nBasic options:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

    puts("\nOptions controlling scoring system:");
    esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 

    puts("\nOptions controlling reporting thresholds:");
    esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

    puts("\nOptions controlling inclusion (significance) thresholds:");
    esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

    puts("\nOptions controlling model-specific thresholding:");
    esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 

    puts("\nOptions controlling acceleration heuristics:");
    esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

    puts("\nOptions controlling E value calibration:");
    esl_opt_DisplayHelp(stdout, go, 11, 2, 80); 

    puts("\nOther expert options:");
    esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
    exit(0);
  }

  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, pgm, usage);
  puts("\nwhere most common options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", pgm);
  exit(1);  
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *copt, const ESL_GETOPTS *sopt, char *seqfile)
{
  p7_banner(ofp, copt->argv[0], banner);
  
  fprintf(ofp, "# target sequence database:        %s\n", seqfile);
  //m/if (esl_opt_IsUsed(copt, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",      esl_opt_GetString(copt, "-o"));
  //m/if (esl_opt_IsUsed(copt, "-A"))          fprintf(ofp, "# MSA of all hits saved to file:   %s\n",      esl_opt_GetString(copt, "-A"));
  //m/if (esl_opt_IsUsed(copt, "--tblout"))    fprintf(ofp, "# per-seq hits tabular output:     %s\n",      esl_opt_GetString(copt, "--tblout"));
  //m/if (esl_opt_IsUsed(copt, "--domtblout")) fprintf(ofp, "# per-dom hits tabular output:     %s\n",      esl_opt_GetString(copt, "--domtblout"));
  if (esl_opt_IsUsed(sopt, "--acc"))       fprintf(ofp, "# prefer accessions over names:    yes\n");
  if (esl_opt_IsUsed(sopt, "--noali"))     fprintf(ofp, "# show alignments in output:       no\n");
  if (esl_opt_IsUsed(sopt, "--notextw"))   fprintf(ofp, "# max ASCII text line length:      unlimited\n");
  if (esl_opt_IsUsed(sopt, "--textw"))     fprintf(ofp, "# max ASCII text line length:      %d\n",             esl_opt_GetInteger(sopt, "--textw"));  
  if (esl_opt_IsUsed(sopt, "--popen"))     fprintf(ofp, "# gap open probability:            %f\n",             esl_opt_GetReal  (sopt, "--popen"));
  if (esl_opt_IsUsed(sopt, "--pextend"))   fprintf(ofp, "# gap extend probability:          %f\n",             esl_opt_GetReal  (sopt, "--pextend"));
  if (esl_opt_IsUsed(sopt, "--mxfile"))    fprintf(ofp, "# subst score matrix:              %s\n",             esl_opt_GetString(sopt, "--mxfile"));
  if (esl_opt_IsUsed(sopt, "-E"))          fprintf(ofp, "# sequence reporting threshold:    E-value <= %g\n",  esl_opt_GetReal(sopt, "-E"));
  if (esl_opt_IsUsed(sopt, "-T"))          fprintf(ofp, "# sequence reporting threshold:    score >= %g\n",    esl_opt_GetReal(sopt, "-T"));
  if (esl_opt_IsUsed(sopt, "--domE"))      fprintf(ofp, "# domain reporting threshold:      E-value <= %g\n",  esl_opt_GetReal(sopt, "--domE"));
  if (esl_opt_IsUsed(sopt, "--domT"))      fprintf(ofp, "# domain reporting threshold:      score >= %g\n",    esl_opt_GetReal(sopt, "--domT"));
  if (esl_opt_IsUsed(sopt, "--incE"))      fprintf(ofp, "# sequence inclusion threshold:    E-value <= %g\n",  esl_opt_GetReal(sopt, "--incE"));
  if (esl_opt_IsUsed(sopt, "--incT"))      fprintf(ofp, "# sequence inclusion threshold:    score >= %g\n",    esl_opt_GetReal(sopt, "--incT"));
  if (esl_opt_IsUsed(sopt, "--incdomE"))   fprintf(ofp, "# domain inclusion threshold:      E-value <= %g\n",  esl_opt_GetReal(sopt, "--incdomE"));
  if (esl_opt_IsUsed(sopt, "--incdomT"))   fprintf(ofp, "# domain inclusion threshold:      score >= %g\n",    esl_opt_GetReal(sopt, "--incdomT"));
  if (esl_opt_IsUsed(sopt, "--cut_ga"))    fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n"); 
  if (esl_opt_IsUsed(sopt, "--cut_nc"))    fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n"); 
  if (esl_opt_IsUsed(sopt, "--cut_tc"))    fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n"); 
  if (esl_opt_IsUsed(sopt, "--max"))       fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n");
  if (esl_opt_IsUsed(sopt, "--F1"))        fprintf(ofp, "# MSV filter P threshold:       <= %g\n", esl_opt_GetReal(sopt, "--F1"));
  if (esl_opt_IsUsed(sopt, "--F2"))        fprintf(ofp, "# Vit filter P threshold:       <= %g\n", esl_opt_GetReal(sopt, "--F2"));
  if (esl_opt_IsUsed(sopt, "--F3"))        fprintf(ofp, "# Fwd filter P threshold:       <= %g\n", esl_opt_GetReal(sopt, "--F3"));
  if (esl_opt_IsUsed(sopt, "--nobias"))    fprintf(ofp, "# biased composition HMM filter:   off\n");
  if (esl_opt_IsUsed(sopt, "--nonull2"))   fprintf(ofp, "# null2 bias corrections:          off\n");
  if (esl_opt_IsUsed(sopt, "--EmL") )      fprintf(ofp, "# seq length, MSV Gumbel mu fit:   %d\n",     esl_opt_GetInteger(sopt, "--EmL"));
  if (esl_opt_IsUsed(sopt, "--EmN") )      fprintf(ofp, "# seq number, MSV Gumbel mu fit:   %d\n",     esl_opt_GetInteger(sopt, "--EmN"));
  if (esl_opt_IsUsed(sopt, "--EvL") )      fprintf(ofp, "# seq length, Vit Gumbel mu fit:   %d\n",     esl_opt_GetInteger(sopt, "--EvL"));
  if (esl_opt_IsUsed(sopt, "--EvN") )      fprintf(ofp, "# seq number, Vit Gumbel mu fit:   %d\n",     esl_opt_GetInteger(sopt, "--EvN"));
  if (esl_opt_IsUsed(sopt, "--EfL") )      fprintf(ofp, "# seq length, Fwd exp tau fit:     %d\n",     esl_opt_GetInteger(sopt, "--EfL"));
  if (esl_opt_IsUsed(sopt, "--EfN") )      fprintf(ofp, "# seq number, Fwd exp tau fit:     %d\n",     esl_opt_GetInteger(sopt, "--EfN"));
  if (esl_opt_IsUsed(sopt, "--Eft") )      fprintf(ofp, "# tail mass for Fwd exp tau fit:   %f\n",     esl_opt_GetReal   (sopt, "--Eft"));
  if (esl_opt_IsUsed(sopt, "-Z"))          fprintf(ofp, "# sequence search space set to:    %.0f\n",    esl_opt_GetReal(sopt, "-Z"));
  if (esl_opt_IsUsed(sopt, "--domZ"))      fprintf(ofp, "# domain search space set to:      %.0f\n",    esl_opt_GetReal(sopt, "--domZ"));
  if (esl_opt_IsUsed(sopt, "--seed"))  {
    if (esl_opt_GetInteger(sopt, "--seed") == 0) fprintf(ofp, "# random number seed:              one-time arbitrary\n");
    else                                       fprintf(ofp, "# random number seed set to:       %d\n", esl_opt_GetInteger(sopt, "--seed"));
  }
  //m/if (esl_opt_IsUsed(copt, "--tformat"))   fprintf(ofp, "# targ <seqfile> format asserted:  %s\n", esl_opt_GetString(copt, "--tformat"));
  //m/if (esl_opt_IsUsed(copt, "--daemon"))    fprintf(ofp, "run as a daemon process\n");
#ifdef HMMER_THREADS
  //m/if (esl_opt_IsUsed(copt, "--cpu"))       fprintf(ofp, "# number of worker threads:        %d\n", esl_opt_GetInteger(copt, "--cpu"));  
#endif
#ifdef HAVE_MPI
  //m/if (esl_opt_IsUsed(copt, "--mpi"))       fprintf(ofp, "# MPI:                             on\n");
#endif
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}

int
main(int argc, char **argv)
{
  int              status   = eslOK;

  ESL_GETOPTS     *go  = NULL;	/* command line processing                 */
  struct cfg_s     cfg;         /* configuration data                      */

  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet) 
   */
  cfg.hmmfile    = NULL;
  cfg.dbfile     = NULL;

  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */

  /* Initializations */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */
  process_commandline(argc, argv, &go, &cfg.dbfile);    

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI
  /* pause the execution of the programs execution until the user has a
   * chance to attach with a debugger and send a signal to resume execution
   * i.e. (gdb) signal SIGCONT
   */
  if (esl_opt_GetBoolean(go, "--stall")) pause();

  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if (cfg.my_rank > 0)  status = mpi_worker(go, &cfg);
      else 		    status = mpi_master(go, &cfg);

      MPI_Finalize();
    }
  else
#endif /*HAVE_MPI*/
    {
      status = serial_master(go, &cfg);
    }

  esl_getopts_Destroy(go);

  return status;
}


/* serial_master()
 * The serial version of hmmsearch.
 * For each query HMM in <hmmfile> search the database for hits.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with p7_Fail().
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;              /* output stream for tabular per-seq (--domtblout) */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w;
  int              textw    = 0;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;

  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif

  QUERY_INFO       queryInfo;
 
  w = esl_stopwatch_Create();

  if (esl_opt_IsOn(go, "--tformat")) {
    dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  

  //m//* Open the query profile HMM file */
  //m/status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp);
  //m/if      (status == eslENOTFOUND) p7_Fail("Failed to open hmm file %s for reading.\n",                      cfg->hmmfile);
  //m/else if (status == eslEFORMAT)   p7_Fail("Unrecognized format, trying to open hmm file %s for reading.\n", cfg->hmmfile);
  //m/else if (status != eslOK)        p7_Fail("Unexpected error %d in opening hmm file %s.\n", status,          cfg->hmmfile);  

  //m//* Open the results output files */
  //m/if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  //m/if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
  //m/if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp")); }
  //m/if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblfp")); }

#ifdef HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                   esl_threads_CPUCount(&ncpus);

  if (ncpus > 0) {
    threadObj = esl_threads_Create(&pipeline_thread);
    queue = esl_workqueue_Create(ncpus * 2);
  }
#endif

  /* seqfile may need to be rewound (multiquery mode) */
  if (!esl_sqfile_IsRewindable(dbfp)) {
    esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);
  }

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  /* initialize the query information struction */
  queryInfo.buf_size = 0;
  queryInfo.buffer   = NULL;
  queryInfo.pgm      = go->argv[0];
  queryInfo.opts     = esl_getopts_Create(searchOpts);

  queryInfo.abc      = esl_alphabet_Create(eslAMINO);
  queryInfo.seq      = NULL;
  queryInfo.hmm      = NULL;

  dbfp->abc          = queryInfo.abc;      //ReadBlock requires knowledge of the alphabet to decide how best to read blocks

#ifdef HMMER_THREADS
  for (i = 0; i < infocnt; ++i)	{
    info[i].queue = queue;
  }

  for (i = 0; i < ncpus * 2; ++i) {
    block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, queryInfo.abc);
    if (block == NULL) 	      esl_fatal("Failed to allocate sequence block");

    status = esl_workqueue_Init(queue, block);
    if (status != eslOK)      esl_fatal("Failed to add block to work queue");
  }
#endif

  /* read query hmm/sequence */
  while (read_QueryInfo(&queryInfo, stdin) == eslOK) {

    /* process any run-time options */
    process_QueryInfo(&queryInfo);

    textw = (esl_opt_GetBoolean(queryInfo.opts, "--notextw")) ? 0 : esl_opt_GetInteger(queryInfo.opts, "--textw");

    output_header(ofp, go, queryInfo.opts, cfg->dbfile);

    esl_stopwatch_Start(w);

    if (queryInfo.hmm == NULL) {
      fprintf(ofp, "Query:       %s  [L=%ld]\n", queryInfo.seq->name, (long) queryInfo.seq->n);
      if (queryInfo.seq->acc[0]  != '\0') fprintf(ofp, "Accession:   %s\n", queryInfo.seq->acc);
      if (queryInfo.seq->desc[0] != '\0') fprintf(ofp, "Description: %s\n", queryInfo.seq->desc);  
    } else {
      fprintf(ofp, "Query:       %s  [M=%d]\n", queryInfo.hmm->name, queryInfo.hmm->M);
      if (queryInfo.hmm->acc)  fprintf(ofp, "Accession:   %s\n", queryInfo.hmm->acc);
      if (queryInfo.hmm->desc) fprintf(ofp, "Description: %s\n", queryInfo.hmm->desc);
    }

    /* Create processing pipeline and hit list */
    for (i = 0; i < infocnt; ++i) {
      info[i].abc   = queryInfo.abc;
      info[i].hmm   = queryInfo.hmm;
      info[i].seq   = queryInfo.seq;
      info[i].opts  = queryInfo.opts;

      info[i].th    = NULL;
      info[i].pli   = NULL;

#ifdef HMMER_THREADS
      if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
    }

#ifdef HMMER_THREADS
    if (ncpus > 0)  
      sstatus = thread_loop(threadObj, queue, dbfp);
    else            
#endif
      sstatus = serial_loop(info, dbfp);
    switch(sstatus){
    case eslEFORMAT: 
      esl_fatal("Parse failed (sequence file %s):\n%s\n", dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
      break;
    case eslEOF:
      /* do nothing */
      break;
    default:
      esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
    }

    /* merge the results of the search results */
    for (i = 1; i < infocnt; ++i) {
      p7_tophits_Merge(info[0].th, info[i].th);
      p7_pipeline_Merge(info[0].pli, info[i].pli);

      p7_pipeline_Destroy(info[i].pli);
      p7_tophits_Destroy(info[i].th);
    }

    /* Print the results.  */
    p7_tophits_Sort(info->th);
    p7_tophits_Threshold(info->th, info->pli);
    p7_tophits_Targets(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
    p7_tophits_Domains(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");

    if (tblfp)    p7_tophits_TabularTargets(tblfp,    queryInfo.hmm->name, queryInfo.hmm->acc, info->th, info->pli, TRUE);
    if (domtblfp) p7_tophits_TabularDomains(domtblfp, queryInfo.hmm->name, queryInfo.hmm->acc, info->th, info->pli, TRUE);
  
    esl_stopwatch_Stop(w);
    p7_pli_Statistics(ofp, info->pli, w);
    fprintf(ofp, "//\n");

    /* Output the results in an MSA (-A option) */
    if (afp) {
      ESL_MSA *msa = NULL;

      if (p7_tophits_Alignment(info->th, queryInfo.abc, NULL, NULL, 0, p7_ALL_CONSENSUS_COLS, &msa) == eslOK) {
        if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
        else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);
	  
        fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
      } 
      else fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n");
	  
      esl_msa_Destroy(msa);
    }

    p7_pipeline_Destroy(info->pli);
    p7_tophits_Destroy(info->th);

    esl_sqfile_Position(dbfp, 0);

    reuse_QueryInfo(&queryInfo);
  } /* end outer loop over query HMMs */

  if (hstatus != eslOK) p7_Fail("Unexpected error (%d) parsing input buffer", hstatus);

#ifdef HMMER_THREADS
  if (ncpus > 0) {
    esl_workqueue_Reset(queue);
    while (esl_workqueue_Remove(queue, (void **) &block) == eslOK) {
      esl_sq_DestroyBlock(block);
    }
    esl_workqueue_Destroy(queue);
    esl_threads_Destroy(threadObj);
  }
#endif

  destroy_QueryInfo(&queryInfo);

  free(info);

  esl_sqfile_Close(dbfp);
  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);

  return eslOK;

 ERROR:
  return eslFAIL;
}

static int
read_QueryInfo(QUERY_INFO *info, FILE *qfp)
{
  int    eoq;
  int    len;
  int    status;

  int    remaining;
  
  char   buffer[MAX_BUFFER];
  char  *tmp;

  if (info->buffer == NULL) {
    info->buf_size = 2 * MAX_BUFFER;
    ESL_ALLOC(info->buffer, info->buf_size);
  }
  remaining = info->buf_size;
  info->buffer[0] = 0;

  info->opts_ptr  = NULL;
  info->query_ptr = NULL;

  eoq = 0;
  while (!eoq) {
    if (fgets(buffer, MAX_BUFFER, qfp) == NULL) return eslEFORMAT;
    len = strlen(buffer);

    /* save the line into the query buffer */
    if (len >= remaining) {
      remaining += info->buf_size;
      info->buf_size *= 2;
      ESL_RALLOC(info->buffer, tmp, info->buf_size);
    }
    strcat(info->buffer, buffer);
    remaining -= len;

    /* check for the end of the query */
    eoq = (buffer[0] == '/' && buffer[1] == '/');
  }

  return eslOK;
 ERROR:
  return eslEMEM;
}

static int
process_QueryInfo(QUERY_INFO *info)
{
  char         *ptr;

  int           status     = eslOK;

  P7_HMMFILE      *hfp     = NULL;              /* open input HMM file                             */
  ESL_SQ          *seq     = NULL;              /* one target sequence (digital)                   */

  /* skip all leading white spaces */
  ptr = info->buffer;
  while (*ptr && isspace(*ptr)) ++ptr;

  /* process search specific options */
  if (*ptr == '@') {
    info->opts_ptr = ++ptr;

    /* skip to the end of the line */
    while (*ptr && (*ptr != '\n' && *ptr == '\r')) ++ptr;

    /* skip remaining white spaces */
    if (*ptr) {
      *ptr++ = 0;
      while (*ptr && isspace(*ptr)) ++ptr;
    }

    process_searchline(info->pgm, info->opts_ptr, info->opts);
  }

  if (*ptr) {
    info->query_ptr = ptr;

    /* try to parse the input buffer as a sequence */
    seq = esl_sq_CreateDigital(info->abc);
    status = esl_sqio_Parse(ptr, strlen(ptr), info->seq, eslSQFILE_DAEMON);
    if (status == eslOK) {
      info->seq = seq;
    } else {
      esl_sq_Destroy(seq);
    }

    /* now try to parse the buffer as an hmm */
    if (status != eslOK) {
      if ((status = p7_hmmfile_OpenBuffer(ptr, strlen(ptr), &hfp)) != eslOK) return status;
      if ((status = p7_hmmfile_Read(hfp, &info->abc,  &info->hmm)) != eslOK) return status;
      p7_hmmfile_Close(hfp);
   }
  }

  return status;
}

static int
reuse_QueryInfo(QUERY_INFO *info)
{
  if (info->hmm)  p7_hmm_Destroy(info->hmm);
  if (info->seq)  esl_sq_Destroy(info->seq);

  info->seq       = NULL;
  info->hmm       = NULL;

  if (info->buffer) info->buffer[0] = 0;
  info->opts_ptr  = NULL;
  info->query_ptr = NULL;

  return eslOK;
}

static int
destroy_QueryInfo(QUERY_INFO *info)
{
  if (info->buffer) {
    free(info->buffer);
    info->buffer = NULL;
  }

  esl_getopts_Destroy(info->opts);

  if (info->abc != NULL) esl_alphabet_Destroy(info->abc);
  if (info->hmm != NULL) p7_hmm_Destroy(info->hmm);
  if (info->seq != NULL) esl_sq_Destroy(info->seq);

  info->opts      = NULL;
  info->abc       = NULL;
  info->seq       = NULL;
  info->hmm       = NULL;

  info->opts_ptr  = NULL;
  info->query_ptr = NULL;

  free(info);

  return eslOK;
}

#ifdef HAVE_MPI

/* Define common tags used by the MPI master/slave processes */
#define HMMER_ERROR_TAG          1
#define HMMER_HMM_TAG            2
#define HMMER_SEQUENCE_TAG       3
#define HMMER_BLOCK_TAG          4
#define HMMER_PIPELINE_TAG       5
#define HMMER_TOPHITS_TAG        6
#define HMMER_HIT_TAG            7
#define HMMER_TERMINATING_TAG    8
#define HMMER_READY_TAG          9

/* mpi_failure()
 * Generate an error message.  If the clients rank is not 0, a
 * message is created with the error message and sent to the
 * master process for handling.
 */
static void
mpi_failure(char *format, ...)
{
  va_list  argp;
  int      status = eslFAIL;
  int      len;
  int      rank;
  char     str[512];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* format the error mesg */
  va_start(argp, format);
  len = vsnprintf(str, sizeof(str), format, argp);
  va_end(argp);

  /* make sure the error string is terminated */
  str[sizeof(str)-1] = '\0';

  /* if the caller is the master, print the results and abort */
  if (rank == 0)
    {
      fprintf(stderr, "\nError: ");
      fprintf(stderr, "%s", str);
      fprintf(stderr, "\n");
      fflush(stderr);

      MPI_Abort(MPI_COMM_WORLD, status);
      exit(1);
    }
  else
    {
      MPI_Send(str, len, MPI_CHAR, 0, HMMER_ERROR_TAG, MPI_COMM_WORLD);
      pause();
    }
}

#define MAX_BLOCK_SIZE (512*1024)

typedef struct {
  uint64_t  offset;
  uint64_t  length;
  uint64_t  count;
} SEQ_BLOCK;

typedef struct {
  int        complete;
  int        size;
  int        current;
  int        last;
  SEQ_BLOCK *blocks;
} BLOCK_LIST;

/* this routine parses the database keeping track of the blocks
 * offset within the file, number of sequences and the length
 * of the block.  These blocks are passed as work units to the
 * MPI workers.  If multiple hmm's are in the query file, the
 * blocks are reused without parsing the database a second time.
 */
int next_block(ESL_SQFILE *sqfp, ESL_SQ *sq, BLOCK_LIST *list, SEQ_BLOCK *block)
{
  int      status   = eslOK;

  /* if the list has been calculated, use it instead of parsing the database */
  if (list->complete)
    {
      if (list->current == list->last)
	{
	  block->offset = 0;
	  block->length = 0;
	  block->count  = 0;

	  status = eslEOF;
	}
      else
	{
	  int inx = list->current++;

	  block->offset = list->blocks[inx].offset;
	  block->length = list->blocks[inx].length;
	  block->count  = list->blocks[inx].count;

	  status = eslOK;
	}

      return status;
    }

  block->offset = 0;
  block->length = 0;
  block->count = 0;

  esl_sq_Reuse(sq);
  while (block->length < MAX_BLOCK_SIZE && (status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK)
    {
      if (block->count == 0) block->offset = sq->roff;
      block->length = sq->eoff - block->offset + 1;
      block->count++;
      esl_sq_Reuse(sq);
    }

  if (status == eslEOF && block->count > 0) status = eslOK;
  if (status == eslEOF) list->complete = 1;

  /* add the block to the list of known blocks */
  if (status == eslOK)
    {
      int inx;

      if (list->last >= list->size)
	{
	  void *tmp;
	  list->size += 500;
	  ESL_RALLOC(list->blocks, tmp, sizeof(SEQ_BLOCK) * list->size);
	}

      inx = list->last++;
      list->blocks[inx].offset = block->offset;
      list->blocks[inx].length = block->length;
      list->blocks[inx].count  = block->count;
    }

  return status;

 ERROR:
  return eslEMEM;
}

/* mpi_master()
 * The MPI version of hmmbuild.
 * Follows standard pattern for a master/worker load-balanced MPI program (J1/78-79).
 * 
 * A master can only return if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * 
 * Recoverable errors include all worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers.
 * 
 * Unrecoverable errors are master-side errors that may affect MPI
 * communication, meaning we cannot count on being able to reach the
 * workers and shut them down. Unrecoverable errors result in immediate
 * p7_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static int
mpi_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;              /* output stream for tabular per-seq (--domtblout) */
  P7_BG           *bg       = NULL;	         /* null model                                      */
  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  ESL_SQ          *dbsq     = NULL;              /* one target sequence (digital)                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w;
  int              textw    = 0;
  int              nquery   = 0;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              dest;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures */
  int              mpi_size = 0;                 /* size of the allocated buffer */
  BLOCK_LIST      *list     = NULL;
  SEQ_BLOCK        block;

  int              i;
  int              size;
  MPI_Status       mpistatus;

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  if (esl_opt_IsOn(go, "--tformat")) {
    dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbfmt == eslSQFILE_UNKNOWN) mpi_failure("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  

  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open hmm file %s for reading.\n",                      cfg->hmmfile);
  else if (status == eslEFORMAT)   mpi_failure("Unrecognized format, trying to open hmm file %s for reading.\n", cfg->hmmfile);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening hmm file %s.\n", status,          cfg->hmmfile);  

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o") && (ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL)
    mpi_failure("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o"));

  if (esl_opt_IsOn(go, "-A") && (afp = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) 
    mpi_failure("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A"));

  if (esl_opt_IsOn(go, "--tblout") && (tblfp = fopen(esl_opt_GetString(go, "--tblout"), "w")) == NULL)
    mpi_failure("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp"));

  if (esl_opt_IsOn(go, "--domtblout") && (domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)
    mpi_failure("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblfp"));

  ESL_ALLOC(list, sizeof(SEQ_BLOCK));
  list->complete = 0;
  list->size     = 0;
  list->current  = 0;
  list->last     = 0;
  list->blocks   = NULL;

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, cfg->hmmfile, cfg->dbfile);
      dbsq = esl_sq_CreateDigital(abc);
      bg = p7_bg_Create(abc);
    }
  
  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hstatus == eslOK) 
    {
      P7_PROFILE      *gm    = NULL;
      P7_OPROFILE     *om    = NULL;       /* optimized query profile                  */
      P7_PIPELINE     *pli   = NULL;
      P7_TOPHITS      *th    = NULL;

      nquery++;
      esl_stopwatch_Start(w);

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1) list->current = 0;

      fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M);
      if (hmm->acc)  fprintf(ofp, "Accession:   %s\n", hmm->acc);
      if (hmm->desc) fprintf(ofp, "Description: %s\n", hmm->desc);

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL);
      p7_oprofile_Convert(gm, om);

      /* Create processing pipeline and hit list */
      th  = p7_tophits_Create(); 
      pli = p7_pipeline_Create(go, hmm->M, 100, FALSE, p7_SEARCH_SEQS);
      p7_pli_NewModel(pli, om, bg);

      /* Main loop: */
      while ((sstatus = next_block(dbfp, dbsq, list, &block)) == eslOK)
	{
	  if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	    mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

	  MPI_Get_count(&mpistatus, MPI_PACKED, &size);
	  if (mpi_buf == NULL || size > mpi_size) {
	    void *tmp;
	    ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	    mpi_size = size; 
	  }

	  dest = mpistatus.MPI_SOURCE;
	  MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

	  if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	    mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
	  if (mpistatus.MPI_TAG != HMMER_READY_TAG)
	    mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
      
	  MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);
	}
      switch(sstatus)
	{
	case eslEFORMAT: 
	  mpi_failure("Parse failed (sequence file %s):\n%s\n", dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
	  break;
	case eslEOF:
	  break;
	default:
	  mpi_failure("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
	}

      block.offset = 0;
      block.length = 0;
      block.count  = 0;

      /* wait for all workers to finish up their work blocks */
      for (i = 1; i < cfg->nproc; ++i)
	{
	  if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	    mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

	  MPI_Get_count(&mpistatus, MPI_PACKED, &size);
	  if (mpi_buf == NULL || size > mpi_size) {
	    void *tmp;
	    ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	    mpi_size = size; 
	  }

	  dest = mpistatus.MPI_SOURCE;
	  MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

	  if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	    mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
	  if (mpistatus.MPI_TAG != HMMER_READY_TAG)
	    mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
	}

      /* merge the results of the search results */
      for (dest = 1; dest < cfg->nproc; ++dest)
	{
	  P7_PIPELINE     *mpi_pli   = NULL;
	  P7_TOPHITS      *mpi_th    = NULL;

	  /* send an empty block to signal the worker they are done */
	  MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);

	  /* wait for the results */
	  if ((status = p7_tophits_MPIRecv(dest, HMMER_TOPHITS_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, &mpi_th)) != eslOK)
	    mpi_failure("Unexpected error %d receiving tophits from %d", status, dest);

	  if ((status = p7_pipeline_MPIRecv(dest, HMMER_PIPELINE_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, go, &mpi_pli)) != eslOK)
	    mpi_failure("Unexpected error %d receiving pipeline from %d", status, dest);

	  p7_tophits_Merge(th, mpi_th);
	  p7_pipeline_Merge(pli, mpi_pli);

	  p7_pipeline_Destroy(mpi_pli);
	  p7_tophits_Destroy(mpi_th);
	}

      /* Print the results.  */
      p7_tophits_Sort(th);
      p7_tophits_Threshold(th, pli);
      p7_tophits_Targets(ofp, th, pli, textw); fprintf(ofp, "\n\n");
      p7_tophits_Domains(ofp, th, pli, textw); fprintf(ofp, "\n\n");

      if (tblfp)    p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, th, pli, (nquery == 1));
      if (domtblfp) p7_tophits_TabularDomains(domtblfp, hmm->name, hmm->acc, th, pli, (nquery == 1));
  
      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, pli, w);
      fprintf(ofp, "//\n");

      /* Output the results in an MSA (-A option) */
      if (afp) {
	ESL_MSA *msa = NULL;

	if (p7_tophits_Alignment(th, abc, NULL, NULL, 0, p7_ALL_CONSENSUS_COLS, &msa) == eslOK)
	  {
	    if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
	    else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);
	  
	    fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
	  } 
	else fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n");
	  
	esl_msa_Destroy(msa);
      }

      p7_pipeline_Destroy(pli);
      p7_tophits_Destroy(th);
      p7_hmm_Destroy(hmm);

      hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
    } /* end outer loop over query HMMs */

  switch(hstatus)
    {
    case eslEOD:
      mpi_failure("read failed, HMM file %s may be truncated?", cfg->hmmfile);
      break;
    case eslEFORMAT:
      mpi_failure("bad file format in HMM file %s", cfg->hmmfile);
      break;
    case eslEINCOMPAT:
      mpi_failure("HMM file %s contains different alphabets", cfg->hmmfile);
      break;
    case eslEOF:
      break;
    default:
      mpi_failure("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->hmmfile);
    }

  /* monitor all the workers to make sure they have ended */
  for (i = 1; i < cfg->nproc; ++i)
    {
      if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

      MPI_Get_count(&mpistatus, MPI_PACKED, &size);
      if (mpi_buf == NULL || size > mpi_size) {
	void *tmp;
	ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	mpi_size = size; 
      }

      dest = mpistatus.MPI_SOURCE;
      MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

      if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
      if (mpistatus.MPI_TAG != HMMER_TERMINATING_TAG)
	mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
    }

  free(list);
  if (mpi_buf != NULL) free(mpi_buf);

  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);

  p7_bg_Destroy(bg);
  esl_sq_Destroy(dbsq);
  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);

  return eslOK;

 ERROR:
  return eslEMEM;
}


static int
mpi_worker(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  ESL_SQ          *dbsq     = NULL;              /* one target sequence (digital)                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  P7_BG           *bg       = NULL;	         /* null model                                      */
  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures           */
  int              mpi_size = 0;                 /* size of the allocated buffer                    */

  MPI_Status       mpistatus;

  w = esl_stopwatch_Create();

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  

  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open hmm file %s for reading.\n",                      cfg->hmmfile);
  else if (status == eslEFORMAT)   mpi_failure("Unrecognized format, trying to open hmm file %s for reading.\n", cfg->hmmfile);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening hmm file %s.\n", status,          cfg->hmmfile);  

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      dbsq = esl_sq_CreateDigital(abc);
      bg = p7_bg_Create(abc);
    }
  
  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hstatus == eslOK) 
    {
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */
      P7_PIPELINE     *pli     = NULL;
      P7_TOPHITS      *th      = NULL;

      SEQ_BLOCK        block;

      esl_stopwatch_Start(w);

      status = 0;
      MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL);
      p7_oprofile_Convert(gm, om);

      th  = p7_tophits_Create(); 
      pli = p7_pipeline_Create(go, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
      p7_pli_NewModel(pli, om, bg);

      /* receive a sequence block from the master */
      MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
      while (block.count > 0)
	{
	  uint64_t length = 0;
	  uint64_t count  = block.count;

	  status = esl_sqfile_Position(dbfp, block.offset);
	  if (status != eslOK) mpi_failure("Cannot position sequence database to %ld\n", block.offset);

	  while (count > 0 && (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
	    {
	      length = dbsq->eoff - block.offset + 1;

	      p7_pli_NewSeq(pli, dbsq);
	      p7_bg_SetLength(bg, dbsq->n);
	      p7_oprofile_ReconfigLength(om, dbsq->n);
      
	      p7_Pipeline(pli, om, bg, dbsq, th);

	      esl_sq_Reuse(dbsq);
	      p7_pipeline_Reuse(pli);

	      --count;
	    }

	  /* lets do a little bit of sanity checking here to make sure the blocks are the same */
	  if (count > 0)              mpi_failure("Block count mismatch - expected %ld found %ld at offset %ld\n",  block.count,  block.count - count, block.offset);
	  if (block.length != length) mpi_failure("Block length mismatch - expected %ld found %ld at offset %ld\n", block.length, length,              block.offset);

	  /* inform the master we need another block of sequences */
	  status = 0;
	  MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);

	  /* wait for the next block of sequences */
	  MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
	}

      esl_stopwatch_Stop(w);

      /* Send the top hits back to the master. */
      p7_tophits_MPISend(th, 0, HMMER_TOPHITS_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);
      p7_pipeline_MPISend(pli, 0, HMMER_PIPELINE_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);

      p7_pipeline_Destroy(pli);
      p7_tophits_Destroy(th);
      p7_oprofile_Destroy(om);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);

      hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
    } /* end outer loop over query HMMs */

  switch(hstatus)
    {
    case eslEOF:
      /* do nothing */
      break;
    case eslEFORMAT:
      mpi_failure("bad file format in HMM file %s", cfg->hmmfile);
      break;
    case eslEINCOMPAT:
      mpi_failure("HMM file %s contains different alphabets", cfg->hmmfile);
      break;
    default:
      mpi_failure("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->hmmfile);
    }

  status = 0;
  MPI_Send(&status, 1, MPI_INT, 0, HMMER_TERMINATING_TAG, MPI_COMM_WORLD);

  if (mpi_buf != NULL) free(mpi_buf);

  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);

  p7_bg_Destroy(bg);
  esl_sq_Destroy(dbsq);
  esl_stopwatch_Destroy(w);

  return eslOK;
}
#endif /*HAVE_MPI*/

static int
serial_loop(WORKER_INFO *info, ESL_SQFILE *dbfp)
{
  int               seed;
  int               status;
  ESL_SQ           *dbsq     = NULL;         /* one target sequence (digital)  */

  P7_BG            *bg       = NULL;	     /* null model                     */
  P7_PIPELINE      *pli      = NULL;         /* work pipeline                  */
  P7_TOPHITS       *th       = NULL;         /* top hit results                */
  P7_PROFILE       *gm       = NULL;         /* generic model                  */
  P7_OPROFILE      *om       = NULL;         /* optimized query profile        */

  P7_BUILDER       *bld      = NULL;         /* HMM construction configuration */

  dbsq = esl_sq_CreateDigital(info->abc);

  /* Convert to an optimized model */
  bg = p7_bg_Create(info->abc);

  /* process a query sequence or hmm */
  if (info->seq != NULL) {
    bld = p7_builder_Create(NULL, info->abc);
    if ((seed = esl_opt_GetInteger(info->opts, "--seed")) > 0) {
      esl_randomness_Init(bld->r, seed);
      bld->do_reseeding = TRUE;
    }
    bld->EmL = esl_opt_GetInteger(info->opts, "--EmL");
    bld->EmN = esl_opt_GetInteger(info->opts, "--EmN");
    bld->EvL = esl_opt_GetInteger(info->opts, "--EvL");
    bld->EvN = esl_opt_GetInteger(info->opts, "--EvN");
    bld->EfL = esl_opt_GetInteger(info->opts, "--EfL");
    bld->EfN = esl_opt_GetInteger(info->opts, "--EfN");
    bld->Eft = esl_opt_GetReal   (info->opts, "--Eft");
    status = p7_builder_SetScoreSystem(bld, esl_opt_GetString(info->opts, "--mxfile"), NULL, esl_opt_GetReal(info->opts, "--popen"), esl_opt_GetReal(info->opts, "--pextend"));
    if (status != eslOK) esl_fatal("Failed to set single query seq score system:\n%s\n", bld->errbuf);

    p7_SingleBuilder(bld, info->seq, bg, NULL, NULL, NULL, &om); /* bypass HMM - only need model */
  } else {
    gm = p7_profile_Create (info->hmm->M, info->abc);
    om = p7_oprofile_Create(info->hmm->M, info->abc);
    p7_ProfileConfig(info->hmm, bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
    p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */
  }

  /* Create processing pipeline and hit list */
  th  = p7_tophits_Create(); 
  pli = p7_pipeline_Create(info->opts, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
  p7_pli_NewModel(pli, om, bg);

  /* Main loop: */
  while ((status = esl_sqio_Read(dbfp, dbsq)) == eslOK) {
    p7_pli_NewSeq(pli, dbsq);
    p7_bg_SetLength(bg, dbsq->n);
    p7_oprofile_ReconfigLength(om, dbsq->n);
      
    p7_Pipeline(pli, om, bg, dbsq, th);
	  
    esl_sq_Reuse(dbsq);
    p7_pipeline_Reuse(pli);
  }

  /* make available the pipeline objects to the main thread */
  info->th = th;
  info->pli = pli;

  /* clean up */
  p7_bg_Destroy(bg);
  p7_oprofile_Destroy(om);

  if (gm  != NULL) p7_profile_Destroy(gm);
  if (bld != NULL) p7_builder_Destroy(bld);

  esl_sq_Destroy(dbsq);

  return status;
}

#ifdef HMMER_THREADS
static int
thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp)
{
  int  status  = eslOK;
  int  sstatus = eslOK;
  int  eofCount = 0;
  ESL_SQ_BLOCK *block;
  void         *newBlock;


#ifdef HAVE_FLUSH_ZERO_MODE
  /* In order to avoid the performance penalty dealing with sub-normal
   * values in the floating point calculations, set the processor flag
   * so sub-normals are "flushed" immediately to zero.
   * On OS X, need to reset this flag for each thread
   * (see TW notes 05/08/10 for details)
   */
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
      
  /* Main loop: */
  while (sstatus == eslOK)
    {
      block = (ESL_SQ_BLOCK *) newBlock;
      sstatus = esl_sqio_ReadBlock(dbfp, block, -1);
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

  return sstatus;
}

static void 
pipeline_thread(void *arg)
{
  int               i;
  int               seed;
  int               status;
  int               workeridx;
  WORKER_INFO      *info;
  ESL_THREADS      *obj;

  ESL_SQ_BLOCK     *block = NULL;
  void             *newBlock;
  
  P7_BG            *bg       = NULL;	     /* null model                     */
  P7_PIPELINE      *pli      = NULL;         /* work pipeline                  */
  P7_TOPHITS       *th       = NULL;         /* top hit results                */
  P7_PROFILE       *gm       = NULL;         /* generic model                  */
  P7_OPROFILE      *om       = NULL;         /* optimized query profile        */

  P7_BUILDER       *bld      = NULL;         /* HMM construction configuration */

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  /* Convert to an optimized model */
  bg = p7_bg_Create(info->abc);

  /* process a query sequence or hmm */
  if (info->seq != NULL) {
    bld = p7_builder_Create(NULL, info->abc);
    if ((seed = esl_opt_GetInteger(info->opts, "--seed")) > 0) {
      esl_randomness_Init(bld->r, seed);
      bld->do_reseeding = TRUE;
    }
    bld->EmL = esl_opt_GetInteger(info->opts, "--EmL");
    bld->EmN = esl_opt_GetInteger(info->opts, "--EmN");
    bld->EvL = esl_opt_GetInteger(info->opts, "--EvL");
    bld->EvN = esl_opt_GetInteger(info->opts, "--EvN");
    bld->EfL = esl_opt_GetInteger(info->opts, "--EfL");
    bld->EfN = esl_opt_GetInteger(info->opts, "--EfN");
    bld->Eft = esl_opt_GetReal   (info->opts, "--Eft");
    status = p7_builder_SetScoreSystem(bld, esl_opt_GetString(info->opts, "--mxfile"), NULL, esl_opt_GetReal(info->opts, "--popen"), esl_opt_GetReal(info->opts, "--pextend"));
    if (status != eslOK) esl_fatal("Failed to set single query seq score system:\n%s\n", bld->errbuf);

    p7_SingleBuilder(bld, info->seq, bg, NULL, NULL, NULL, &om); /* bypass HMM - only need model */
  } else {
    gm = p7_profile_Create (info->hmm->M, info->abc);
    om = p7_oprofile_Create(info->hmm->M, info->abc);
    p7_ProfileConfig(info->hmm, bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
    p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */
  }

  /* Create processing pipeline and hit list */
  th  = p7_tophits_Create(); 
  pli = p7_pipeline_Create(info->opts, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
  p7_pli_NewModel(pli, om, bg);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  /* loop until all blocks have been processed */
  block = (ESL_SQ_BLOCK *) newBlock;
  while (block->count > 0)
    {
      /* Main loop: */
      for (i = 0; i < block->count; ++i)
	{
	  ESL_SQ *dbsq = block->list + i;

	  p7_pli_NewSeq(pli, dbsq);
	  p7_bg_SetLength(bg, dbsq->n);
	  p7_oprofile_ReconfigLength(om, dbsq->n);
	  
	  p7_Pipeline(pli, om, bg, dbsq, th);
	  
	  esl_sq_Reuse(dbsq);
	  p7_pipeline_Reuse(pli);
	}

      status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      block = (ESL_SQ_BLOCK *) newBlock;
    }

  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  /* make available the pipeline objects to the main thread */
  info->th = th;
  info->pli = pli;

  /* clean up */
  p7_bg_Destroy(bg);
  p7_oprofile_Destroy(om);

  if (gm != NULL)  p7_profile_Destroy(gm);
  if (bld != NULL) p7_builder_Destroy(bld);

  esl_threads_Finished(obj, workeridx);

  return;
}
#endif   /* HMMER_THREADS */
 

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

