/* ssv-hits: given a profile HMM and a nucleotide sequence, find the shortest diagonal that meets SSV filter threshold
 * 
 * TJW, Mon Oct  5 06:53:36 EDT 2010 [Janelia]
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_scorematrix.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"


#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"


/* set the max residue count to 1 meg when reading a block */
#ifdef P7_IMPL_DUMMY_INCLUDED
#define NHMMER_MAX_RESIDUE_COUNT (1024 * 100)
#else
#define NHMMER_MAX_RESIDUE_COUNT MAX_RESIDUE_COUNT   /*from esl_sqio_(ascii|ncbi).c*/
#endif

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  P7_BG            *bg;             /* null model                              */
  P7_PIPELINE      *pli;         /* work pipeline                           */
  P7_TOPHITS       *th;          /* top hit results                         */
  P7_OPROFILE      *om;          /* optimized query profile                 */
} WORKER_INFO;

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#define CPUOPTS     NULL
#define MPIOPTS     NULL


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                        2 },
  { "-A",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save multiple alignment of all hits to file <s>",              2 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-sequence hits to file <s>",        2 },
//  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-domain hits to file <s>",          2 },
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                       2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                         2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                     2 },
  /* Control of scoring system */
  { "--popen",      eslARG_REAL,  "0.02", NULL, "0<=x<0.5",NULL,  NULL,  NULL,              "gap open probability",                                         3 },
  { "--pextend",    eslARG_REAL,   "0.4", NULL, "0<=x<1",  NULL,  NULL,  NULL,              "gap extend probability",                                       3 },
  { "--mxfile",     eslARG_INFILE,  NULL, NULL, NULL,      NULL,  NULL,  NULL,              "substitution score matrix [default: BLOSUM62]",                3 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report sequences <= this E-value threshold in output",         4 },
  { "-T",           eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report sequences >= this score threshold in output",           4 },
  /* Control of inclusion (significance) thresholds */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider sequences >= this score threshold as significant",    5 },
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
/* Other options */
  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,    NULL,  NULL,  NULL,            "turn off biased composition score corrections",               12 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set database size in *Mb* to <x> for E-value calculations",   12 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--w_beta",     eslARG_REAL,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "tail mass at which window length is determined",               12 },
  { "--w_length",   eslARG_INT, NULL,  NULL, NULL,   NULL,  NULL,  NULL,            "window length ",                                              12 },
  { "--single",     eslARG_NONE, NULL,  NULL, NULL,   NULL,  NULL,  NULL,            "don't search reverse complement of database sequences ",         12 },

/* Not used, but retained because esl option-handling code errors if it isn't kept here.  Placed in group 99 so it doesn't print to help*/
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "Not used",   99 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  DOMREPOPTS,      "Not used",   99 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  DOMREPOPTS,      "Not used",   99 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCDOMOPTS,      "Not used",    99 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCDOMOPTS,      "Not used",      99 },
/* will eventually bring these back, but store in group 99 for now, so they don't print to help*/
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert query <seqfile> is in format <s>: no autodetection",   99 },
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert target <seqfile> is in format <s>>: no autodetection", 99 },


#ifdef HMMER_THREADS 
  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>=0",NULL,  NULL,  CPUOPTS,         "number of parallel CPU workers to use for multithreads",      12 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


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

//static char usage[]  = "[options] <query hmmfile|alignfile> <target seqfile>";
//static char banner[] = "search a DNA model or alignment against a DNA database";
static char usage[]  = "[options] <query hmmfile> <target seqfile>";
static char banner[] = "given a profile HMM and a nucleotide sequence, find the shortest diagonal that meets SSV filter threshold";


static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ESL_SQFILE *dbfp);
#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(WORKER_INFO *info, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/


static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_hmmfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     p7_Die("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nBasic options:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/

      puts("\nOptions directing output:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

//      puts("\nOptions controlling scoring system:");
//     esl_opt_DisplayHelp(stdout, go, 3, 2, 80);

      puts("\nOptions controlling reporting thresholds:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      puts("\nOptions controlling inclusion (significance) thresholds:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      puts("\nOptions controlling acceleration heuristics:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);

      puts("\nOther expert options:");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
      exit(0);

  }

  if (esl_opt_ArgNumber(go)                  != 2)     { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <hmmfile> argument on command line"); goto ERROR; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }
  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *hmmfile, char *seqfile)
{
  p7_banner(ofp, go->argv[0], banner);
  
  fprintf(ofp, "# query HMM file:                  %s\n", hmmfile);
  fprintf(ofp, "# target sequence database:        %s\n", seqfile);
  if (esl_opt_IsUsed(go, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",      esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "-A"))          fprintf(ofp, "# MSA of all hits saved to file:   %s\n",      esl_opt_GetString(go, "-A"));
  if (esl_opt_IsUsed(go, "--tblout"))    fprintf(ofp, "# per-seq hits tabular output:     %s\n",      esl_opt_GetString(go, "--tblout"));
  if (esl_opt_IsUsed(go, "--acc"))       fprintf(ofp, "# prefer accessions over names:    yes\n");
  if (esl_opt_IsUsed(go, "--noali"))     fprintf(ofp, "# show alignments in output:       no\n");
  if (esl_opt_IsUsed(go, "--notextw"))   fprintf(ofp, "# max ASCII text line length:      unlimited\n");
  if (esl_opt_IsUsed(go, "--textw"))     fprintf(ofp, "# max ASCII text line length:      %d\n",             esl_opt_GetInteger(go, "--textw"));
  if (esl_opt_IsUsed(go, "--popen"))     fprintf(ofp, "# gap open probability:            %f\n",             esl_opt_GetReal  (go, "--popen"));
  if (esl_opt_IsUsed(go, "--pextend"))   fprintf(ofp, "# gap extend probability:          %f\n",             esl_opt_GetReal  (go, "--pextend"));
  if (esl_opt_IsUsed(go, "--mxfile"))    fprintf(ofp, "# subst score matrix:              %s\n",             esl_opt_GetString(go, "--mxfile"));
  if (esl_opt_IsUsed(go, "-E"))          fprintf(ofp, "# sequence reporting threshold:    E-value <= %g\n",  esl_opt_GetReal(go, "-E"));
  if (esl_opt_IsUsed(go, "-T"))          fprintf(ofp, "# sequence reporting threshold:    score >= %g\n",    esl_opt_GetReal(go, "-T"));
  if (esl_opt_IsUsed(go, "--domE"))      fprintf(ofp, "# domain reporting threshold:      E-value <= %g\n",  esl_opt_GetReal(go, "--domE"));
  if (esl_opt_IsUsed(go, "--domT"))      fprintf(ofp, "# domain reporting threshold:      score >= %g\n",    esl_opt_GetReal(go, "--domT"));
  if (esl_opt_IsUsed(go, "--incE"))      fprintf(ofp, "# sequence inclusion threshold:    E-value <= %g\n",  esl_opt_GetReal(go, "--incE"));
  if (esl_opt_IsUsed(go, "--incT"))      fprintf(ofp, "# sequence inclusion threshold:    score >= %g\n",    esl_opt_GetReal(go, "--incT"));
  if (esl_opt_IsUsed(go, "--incdomE"))   fprintf(ofp, "# domain inclusion threshold:      E-value <= %g\n",  esl_opt_GetReal(go, "--incdomE"));
  if (esl_opt_IsUsed(go, "--incdomT"))   fprintf(ofp, "# domain inclusion threshold:      score >= %g\n",    esl_opt_GetReal(go, "--incdomT"));
//  if (esl_opt_IsUsed(go, "--cut_ga"))    fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n");
//  if (esl_opt_IsUsed(go, "--cut_nc"))    fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n");
//  if (esl_opt_IsUsed(go, "--cut_tc"))    fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n");
  if (esl_opt_IsUsed(go, "--max"))       fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n");
  if (esl_opt_IsUsed(go, "--F1"))        fprintf(ofp, "# MSV filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F1"));
  if (esl_opt_IsUsed(go, "--F2"))        fprintf(ofp, "# Vit filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F2"));
  if (esl_opt_IsUsed(go, "--F3"))        fprintf(ofp, "# Fwd filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F3"));
  if (esl_opt_IsUsed(go, "--nobias"))    fprintf(ofp, "# biased composition HMM filter:   off\n");
  if (esl_opt_IsUsed(go, "--nonull2"))   fprintf(ofp, "# null2 bias corrections:          off\n");
  if (esl_opt_IsUsed(go, "--single"))    fprintf(ofp, "# search reverse complement:       off\n");

  if (esl_opt_IsUsed(go, "-Z"))          fprintf(ofp, "# database size is set to:         %.1f Mb\n",    esl_opt_GetReal(go, "-Z"));
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0) fprintf(ofp, "# random number seed:              one-time arbitrary\n");
    else                                       fprintf(ofp, "# random number seed set to:       %d\n", esl_opt_GetInteger(go, "--seed"));
  }
//  if (esl_opt_IsUsed(go, "--qformat"))   fprintf(ofp, "# query <seqfile> format asserted: %s\n",     esl_opt_GetString(go, "--qformat"));
  if (esl_opt_IsUsed(go, "--tformat"))   fprintf(ofp, "# targ <seqfile> format asserted:  %s\n", esl_opt_GetString(go, "--tformat"));
  if (esl_opt_IsUsed(go, "--w_beta"))
                                           fprintf(ofp, "# window length beta value:        %g\n", esl_opt_GetReal(go, "--w_beta"));
  if (esl_opt_IsUsed(go, "--w_length") )
                                         fprintf(ofp, "# window length :                  %d\n", esl_opt_GetInteger(go, "--w_length"));
  #ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu"))       fprintf(ofp, "# number of worker threads:        %d\n", esl_opt_GetInteger(go, "--cpu"));  
#endif
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}

int
main(int argc, char **argv)
{
  int              status   = eslOK;

  ESL_GETOPTS     *go  = NULL;    /* command line processing                 */
  struct cfg_s     cfg;         /* configuration data                      */

  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet)
   */
  cfg.hmmfile    = NULL;
  cfg.dbfile     = NULL;

  cfg.do_mpi     = FALSE;               /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;                   /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;                   /* this gets reset below, if we init MPI */

  /* Initializations */
  p7_FLogsumInit();        /* we're going to use table-driven Logsum() approximations at times */
  process_commandline(argc, argv, &go, &cfg.hmmfile, &cfg.dbfile);    

  status = serial_master(go, &cfg);

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
  //FILE            *domtblfp = NULL;              /* output stream for tabular per-seq (--domtblout) */

  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  /*either the two above or ones below will be used*/
  //int              qhformat  = eslSQFILE_UNKNOWN;  /* format of qfile                                  */
  //ESL_SQFILE      *qfp      = NULL;          /* open qfile                                       */
  //ESL_SQ          *qsq      = NULL;               /* query sequence                                   */

  int              dbformat = eslSQFILE_UNKNOWN;  /* format of dbfile                                 */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */

  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  //int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w;
  //int              seed;
  int              textw    = 0;
  int              nquery   = 0;
  int              status   = eslOK;
  int              qhstatus  = eslOK;
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

  double window_beta = -1.0 ;
  int window_length  = -1;
  if (esl_opt_IsUsed(go, "--w_beta")) { if (  ( window_beta   = esl_opt_GetReal(go, "--w_beta") )  < 0 || window_beta > 1  ) esl_fatal("Invalid window-length beta value\n"); }
  if (esl_opt_IsUsed(go, "--w_length")) { if (( window_length = esl_opt_GetInteger(go, "--w_length")) < 4  ) esl_fatal("Invalid window length value\n"); }

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");


  /* If caller declared input formats, decode them */
/*  if (esl_opt_IsOn(go, "--qformat")) {
    qhformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (qhformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }
  */
  if (esl_opt_IsOn(go, "--tformat")) {
    dbformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbformat, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open target sequence database %s for reading\n",      cfg->dbfile);
  else if (status == eslEFORMAT)   esl_fatal("Target sequence database file %s is empty or misformatted\n",   cfg->dbfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        esl_fatal("Unexpected error %d opening target sequence database file %s\n", status, cfg->dbfile);


  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n",                      cfg->hmmfile);
  else if (status == eslEFORMAT)   esl_fatal("Unrecognized format, trying to open hmm file %s for reading.\n", cfg->hmmfile);
  else if (status != eslOK)        esl_fatal("Unexpected error %d in opening hmm file %s.\n", status,          cfg->hmmfile);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp")); }

#ifdef HMMER_THREADS
  /* initialize thread data */
//  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
//  else                                   esl_threads_CPUCount(&ncpus);


  if (ncpus > 0) {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
  }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  /* <abc> is not known 'til first HMM is read.   Could be DNA or RNA*/
  qhstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (qhstatus == eslOK) {
      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, cfg->hmmfile, cfg->dbfile);
      dbfp->abc = abc;

      for (i = 0; i < infocnt; ++i)    {
          info[i].pli   = NULL;
          info[i].th    = NULL;
          info[i].om    = NULL;
          info[i].bg    = p7_bg_Create(abc);
#ifdef HMMER_THREADS
          info[i].queue = queue;
#endif
      }

#ifdef HMMER_THREADS
      for (i = 0; i < ncpus * 2; ++i) {
          block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abc);
          if (block == NULL)           esl_fatal("Failed to allocate sequence block");

          status = esl_workqueue_Init(queue, block);
          if (status != eslOK)          esl_fatal("Failed to add block to work queue");
      }
#endif
  }

  /* Outer loop: over each query HMM in <hmmfile>. */
  while (qhstatus == eslOK) {

      if (window_length > 0)
          hmm->max_length = window_length;
      else if (window_beta > 0)
          p7_Builder_MaxLength(hmm, window_beta);
      else if (hmm->max_length == 0 ) {
          p7_Builder_MaxLength(hmm, p7_DEFAULT_WINDOW_BETA);
      }

      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */

      nquery++;
      esl_stopwatch_Start(w);

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1) {
          if (! esl_sqfile_IsRewindable(dbfp))
            esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);
          esl_sqfile_Position(dbfp, 0);
      }

      fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M);
      if (hmm->acc)  fprintf(ofp, "Accession:   %s\n", hmm->acc);
      if (hmm->desc) fprintf(ofp, "Description: %s\n", hmm->desc);

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, info->bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
      p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */



      for (i = 0; i < infocnt; ++i) {
          /* Create processing pipeline and hit list */
          info[i].th  = p7_tophits_Create();
          info[i].om  = p7_oprofile_Clone(om);
          info[i].pli = p7_pipeline_Create(go, om->M, 100, TRUE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
          p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);
          info[i].pli->single_strand = esl_opt_IsUsed(go, "--single");

#ifdef HMMER_THREADS
          if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
      }

#ifdef HMMER_THREADS
      if (ncpus > 0)  sstatus = thread_loop(info, threadObj, queue, dbfp);
      else            sstatus = serial_loop(info, dbfp);
#else
      sstatus = serial_loop(info, dbfp);
#endif
      switch(sstatus) {
        case eslEFORMAT:
          esl_fatal("Parse failed (sequence file %s):\n%s\n",
                    dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
          break;
        case eslEOF:
          /* do nothing */
          break;
        default:
          esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
      }

      //need to re-compute e-values before merging (when list will be sorted)
      double resCnt = 0.0;
      if (esl_opt_IsUsed(go, "-Z"))
    	  resCnt = 1000000*esl_opt_GetReal(go, "-Z");
      else
    	  for (i = 0; i < infocnt; ++i)
    		  resCnt += info[i].pli->nres;

      for (i = 0; i < infocnt; ++i)
          p7_tophits_ComputeNhmmerEvalues(info[i].th, resCnt);


      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i) {
          p7_tophits_Merge(info[0].th, info[i].th);
          p7_pipeline_Merge(info[0].pli, info[i].pli);

          p7_pipeline_Destroy(info[i].pli);
          p7_tophits_Destroy(info[i].th);
          p7_oprofile_Destroy(info[i].om);
      }


      /* Print the results.  */
      p7_tophits_Sort(info->th);
      p7_tophits_RemoveDuplicates(info->th);

      p7_tophits_Threshold(info->th, info->pli);

      //tally up total number of hits and target coverage
      info->pli->n_output = info->pli->pos_output = 0;
      for (i = 0; i < info->th->N; i++) {
          if (info->th->hit[i]->dcl[0].is_reported || info->th->hit[i]->dcl[0].is_included) {
              info->pli->n_output++;
              info->pli->pos_output += abs(info->th->hit[i]->dcl[0].jali - info->th->hit[i]->dcl[0].iali) + 1;
          }
      }


      p7_tophits_Targets(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
      p7_tophits_Domains(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");

      if (tblfp)    p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, info->th, info->pli, (nquery == 1));
  
      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, info->pli, w);
      fprintf(ofp, "//\n");

      /* Output the results in an MSA (-A option) */
      if (afp) {
          ESL_MSA *msa = NULL;

          if (p7_tophits_Alignment(info->th, abc, NULL, NULL, 0, p7_DEFAULT, &msa) == eslOK) {
            if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
            else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);

            fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
          }  else {
              fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n");
          }

          esl_msa_Destroy(msa);
      }

      p7_pipeline_Destroy(info->pli);
      p7_tophits_Destroy(info->th);
      p7_oprofile_Destroy(info->om);
      p7_oprofile_Destroy(om);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);

      qhstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  } /* end outer loop over query HMMs */

  switch(qhstatus) {
    case eslEOD:
      p7_Fail("read failed, HMM file %s may be truncated?", cfg->hmmfile);
      break;
    case eslEFORMAT:
      p7_Fail("bad file format in HMM file %s", cfg->hmmfile);
      break;
    case eslEINCOMPAT:
      p7_Fail("HMM file %s contains different alphabets", cfg->hmmfile);
      break;
    case eslEOF:
      /* do nothing */
      break;
    default:
      p7_Fail("Unexpected error (%d) in reading HMMs from %s", qhstatus, cfg->hmmfile);
  }

  for (i = 0; i < infocnt; ++i) {
      p7_bg_Destroy(info[i].bg);
  }

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

  free(info);

  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  //if (domtblfp)      fclose(domtblfp);

  return eslOK;

 ERROR:
  return eslFAIL;
}

//TODO: MPI code needs to be added here


int
get_shortest_ssv_hit(const ESL_DSQ *dsq, int L, P7_PROFILE *gm, float nu, P7_BG *bg, double P) {
	  int 	  	   status;
	  float        tloop = logf((float) gm->max_length / (float) (gm->max_length+3));
	  float        tmove = logf(     3.0f / (float) (gm->max_length+3));
	  float        tbmk  = logf(     2.0f / ((float) gm->M * (float) (gm->M+1)));
	  float        tec   = logf(1.0f / nu);
	  int          i,j,k, x;

	  float 	   tloop_total = tloop * gm->max_length;
	  float nullsc;
	  float invP = esl_gumbel_invsurv(P, gm->evparam[p7_MMU],  gm->evparam[p7_MLAMBDA]);
	  p7_bg_SetLength(bg, gm->max_length);
	  p7_ReconfigLength(gm, gm->max_length);
	  p7_bg_NullOne  (bg, dsq, gm->max_length, &nullsc);
	  float sc_thresh =   nullsc  + (invP * eslCONST_LOG2) - tmove - tloop_total - tmove - tbmk - tec;


	  int test = 1;
	  /*gather values from gm->rsc into a succinct 2D array*/
	  if (test == 1) {
		  printf ("         ");
		  for (k = 1; k <= gm->M; k++) {
			  printf (" %3d   ", k);
		  }
		  printf("\n");
	  }
	  float **scores;
	  ESL_ALLOC(scores, (L+1) * sizeof(float*));
	  for (i=1; i<=L; i++) {
		  ESL_ALLOC(scores[i], (gm->M + 1) * sizeof(float));
		  if (test == 1)  printf ("%3d (%c): ", i, gm->abc->sym[dsq[i]]);
		  for (k = 1; k <= gm->M; k++) {
			  scores[i][k] = gm->rsc[dsq[i]][(k) * p7P_NR     + p7P_MSC];
			  if (test == 1) printf ("%5.2f  ", scores[i][k]);
		  }
		  if (test == 1) printf ("\n");
	  }

	  //For each position in the model, start the query there, and scan along for high-scoring runs.
	  //Keep two pointers (head, tail): keep both together until a positive score is seen, then
	  //extend head while score remains positive. Bring tail up to head if score goes negative. If score
	  //exceeds sc_thresh, note the length (keeping the shortest),
	  float sc;
	  float run_sc;
	  float max_run_sc;
	  int  max_run_i;
	  int  max_run_j;
	  int  max_run_k;
	  int shortest_len = 100000;
	  int shortest_len_i ;
	  int shortest_len_j ;
	  int shortest_len_k ;
	  for (k = 1; k <= gm->M; k++) {
		  i = 1;
		  run_sc = 0;
		  for (j = 1; j<=L && j+k-1<=gm->M; j++)
		  {

			  sc = scores[j][j+k-1];
			  run_sc += sc;
			  if (run_sc > max_run_sc) {
				  max_run_sc = run_sc;
				  max_run_i = i;
				  max_run_j = j;
				  max_run_k = k;
			  }
			  if (run_sc <= 0 ){
				  i=j+1; // j will catch up in a moment
				  run_sc = 0;
			  } else {
				  if (run_sc >= sc_thresh) {
					  //now scan between i and j, checking if the score of any suffix of that subalignment is still above threshold
					  float tmp_run_sc = run_sc;

					  for (x=i; x<j; x++) {

						  sc = scores[x][x+k-1];
						  tmp_run_sc -= sc;
						  if (tmp_run_sc >= sc_thresh) {
							  i = x+1;
							  run_sc = tmp_run_sc;
						  }
					  }
					  if (j-i+1 < shortest_len) {
						  shortest_len = j-i+1;
						  shortest_len_i = i;
						  shortest_len_j = j;
						  shortest_len_k = k;
					  }
				  }
			  }
		  }
	  }

	  printf ("Max: %.2f (%.2f) (i=%d, j=%d, k=%d)\n", max_run_sc, sc_thresh, max_run_i, max_run_j, max_run_k);
	  if (shortest_len < 100000) {
		  printf ("Shortest: %d\n", shortest_len);
		  printf ("i=%d, j=%d, k=%d\n", shortest_len_i, shortest_len_j, shortest_len_k );
	  }
	  exit(0);

	ERROR:
	  return eslFAIL;

}


static int
serial_loop(WORKER_INFO *info, ESL_SQFILE *dbfp)
{

  int      wstatus;
  int i;
  int prev_hit_cnt;
  ESL_SQ   *dbsq   =  esl_sq_CreateDigital(info->om->abc);
#ifdef eslAUGMENT_ALPHABET
  ESL_SQ   *dbsq_revcmp;
  if (dbsq->abc->complement != NULL)
      dbsq_revcmp =  esl_sq_CreateDigital(info->om->abc);
#endif /*eslAUGMENT_ALPHABET*/


  wstatus = esl_sqio_ReadWindow(dbfp, 0, NHMMER_MAX_RESIDUE_COUNT, dbsq);

  while (wstatus == eslOK ) {

      p7_pli_NewSeq(info->pli, dbsq);
      info->pli->nres -= dbsq->C; // to account for overlapping region of windows
      prev_hit_cnt = info->th->N;

      //p7_Pipeline_LongTarget(info->pli, info->om, info->bg, dbsq, info->th);
      //Here - compare dbsq to om (in dummy mode, thats a gm object)
      get_shortest_ssv_hit(dbsq->dsq, dbsq->n, info->om, 2.0, info->bg, info->pli->F1);



      p7_pipeline_Reuse(info->pli); // prepare for next search



      P7_DOMAIN *dcl;
      // modify hit positions to account for the position of the window in the full sequence
      for (i=prev_hit_cnt; i < info->th->N ; i++) {
          dcl = info->th->unsrt[i].dcl;
          dcl->ienv += dbsq->start - 1;
          dcl->jenv += dbsq->start - 1;
          dcl->iali += dbsq->start - 1;
          dcl->jali += dbsq->start - 1;
          dcl->ad->sqfrom += dbsq->start - 1;
          dcl->ad->sqto += dbsq->start - 1;
      }

#ifdef eslAUGMENT_ALPHABET
      //reverse complement
      if (!info->pli->single_strand && dbsq->abc->complement != NULL )
      {
          prev_hit_cnt = info->th->N;
          esl_sq_Copy(dbsq,dbsq_revcmp);
          esl_sq_ReverseComplement(dbsq_revcmp);
          p7_Pipeline_LongTarget(info->pli, info->om, info->bg, dbsq_revcmp, info->th);
          p7_pipeline_Reuse(info->pli); // prepare for next search

          for (i=prev_hit_cnt; i < info->th->N ; i++) {
              dcl = info->th->unsrt[i].dcl;
              // modify hit positions to account for the position of the window in the full sequence
              dcl->ienv = dbsq_revcmp->start - dcl->ienv + 1;
              dcl->jenv = dbsq_revcmp->start - dcl->jenv + 1;
              dcl->iali = dbsq_revcmp->start - dcl->iali + 1;
              dcl->jali = dbsq_revcmp->start - dcl->jali + 1;
              dcl->ad->sqfrom = dbsq_revcmp->start - dcl->ad->sqfrom + 1;
              dcl->ad->sqto = dbsq_revcmp->start - dcl->ad->sqto + 1;

          }

          info->pli->nres += dbsq_revcmp->W;

      }
#endif /*eslAUGMENT_ALPHABET*/

      wstatus = esl_sqio_ReadWindow(dbfp, info->om->max_length, NHMMER_MAX_RESIDUE_COUNT, dbsq);
      if (wstatus == eslEOD) { // no more left of this sequence ... move along to the next sequence.
          info->pli->nseqs++;
          esl_sq_Reuse(dbsq);
          wstatus = esl_sqio_ReadWindow(dbfp, 0, NHMMER_MAX_RESIDUE_COUNT, dbsq);
      }
    }

  esl_sq_Destroy(dbsq);

  return wstatus;

}

#ifdef HMMER_THREADS
static int
thread_loop(WORKER_INFO *info, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp)
{

  //int      wstatus, wstatus_next;
  int  status  = eslOK;
  int  sstatus = eslOK;
  int  eofCount = 0;
  ESL_SQ_BLOCK *block;
  void         *newBlock;
  int i;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
  ((ESL_SQ_BLOCK *)newBlock)->complete = TRUE;

  /* Main loop: */
  while (sstatus == eslOK ) {
      block = (ESL_SQ_BLOCK *) newBlock;

      //reset block as an empty vessel, possibly keeping the first sq intact for reading in the next window
      if (block->count > 0 && block->complete)
          esl_sq_Reuse(block->list);
      for (i=1; i<block->count; i++)
          esl_sq_Reuse(block->list + i);

      sstatus = esl_sqio_ReadBlock(dbfp, block, NHMMER_MAX_RESIDUE_COUNT);

      info->pli->nseqs += block->count - (block->complete ? 0 : 1);// if there's an incomplete sequence read into the block wait to count it until it's complete.

      if (sstatus == eslEOF) {
          if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
          ++eofCount;
      }

      if (sstatus == eslOK) {
          status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
          if (status != eslOK) esl_fatal("Work queue reader failed");
      }

      //newBlock needs all this information so the next ReadBlock call will know what to do
      ((ESL_SQ_BLOCK *)newBlock)->complete = block->complete;
      if (!block->complete) {
          // the final sequence on the block was a probably-incomplete window of the active sequence,
          // so prep the next block to read in the next window
          esl_sq_Copy(block->list + (block->count - 1) , ((ESL_SQ_BLOCK *)newBlock)->list);
          ((ESL_SQ_BLOCK *)newBlock)->list->C = info->om->max_length;
      }

  }

  status = esl_workqueue_ReaderUpdate(queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF) {
      /* wait for all the threads to complete */
      esl_threads_WaitForFinish(obj);
      esl_workqueue_Complete(queue);  
    }

  return sstatus;
}

static void 
pipeline_thread(void *arg)
{
  int prev_hit_cnt;
  int i, j;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;

  ESL_SQ_BLOCK  *block = NULL;
  void          *newBlock;
  
#ifdef HAVE_FLUSH_ZERO_MODE
  /* In order to avoid the performance penalty dealing with sub-normal
   * values in the floating point calculations, set the processor flag
   * so sub-normals are "flushed" immediately to zero.
   * On OS X, need to reset this flag for each thread
   * (see TW notes 05/08/10 for details)
   */
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

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

      p7_pli_NewSeq(info->pli, dbsq);
      info->pli->nres -= dbsq->C; // to account for overlapping region of windows

      prev_hit_cnt = info->th->N;
      p7_Pipeline_LongTarget(info->pli, info->om, info->bg, dbsq, info->th);
      p7_pipeline_Reuse(info->pli); // prepare for next search


      P7_DOMAIN *dcl;
      // modify hit positions to account for the position of the window in the full sequence
      for (j=prev_hit_cnt; j < info->th->N ; ++j) {
          dcl = info->th->unsrt[j].dcl;
          dcl->ienv += dbsq->start - 1;
          dcl->jenv += dbsq->start - 1;
          dcl->iali += dbsq->start - 1;
          dcl->jali += dbsq->start - 1;
          dcl->ad->sqfrom += dbsq->start - 1;
          dcl->ad->sqto += dbsq->start - 1;
      }


#ifdef eslAUGMENT_ALPHABET
      //reverse complement
      if (!info->pli->single_strand && dbsq->abc->complement != NULL)
      {
          prev_hit_cnt = info->th->N;
          esl_sq_ReverseComplement(dbsq);
          p7_Pipeline_LongTarget(info->pli, info->om, info->bg, dbsq, info->th);
          p7_pipeline_Reuse(info->pli); // prepare for next search

          for (j=prev_hit_cnt; j < info->th->N ; ++j) {
              dcl = info->th->unsrt[j].dcl;
              // modify hit positions to account for the position of the window in the full sequence
              dcl->ienv = dbsq->start - dcl->ienv + 1;
              dcl->jenv = dbsq->start - dcl->jenv + 1;
              dcl->iali = dbsq->start - dcl->iali + 1;
              dcl->jali = dbsq->start - dcl->jali + 1;
              dcl->ad->sqfrom = dbsq->start - dcl->ad->sqfrom + 1;
              dcl->ad->sqto = dbsq->start - dcl->ad->sqto + 1;

          }

          info->pli->nres += dbsq->W;
      }

#endif /*eslAUGMENT_ALPHABET*/

    }

      status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      block = (ESL_SQ_BLOCK *) newBlock;
    }

  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
  return;
}
#endif   /* HMMER_THREADS */
 

/*****************************************************************
 * @LICENSE@
 *****************************************************************/



