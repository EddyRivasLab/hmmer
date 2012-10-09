/* nhmmer: search profile HMM(s) against a nucleotide sequence database.
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_scorematrix.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"


#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"


/* set the max residue count to 1/4 meg when reading a block */
#ifdef P7_IMPL_DUMMY_INCLUDED
#include "esl_vectorops.h"
#define NHMMER_MAX_RESIDUE_COUNT (1024 * 100)
#else
#define NHMMER_MAX_RESIDUE_COUNT (1024 * 256)  /* 1/4 Mb */
#endif

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  P7_BG            *bg;             /* null model                              */
  P7_PIPELINE      *pli;         /* work pipeline                           */
  P7_TOPHITS       *th;          /* top hit results                         */
  P7_OPROFILE      *om;          /* optimized query profile                 */
  FM_CFG           *fm_cfg;      /* global data for FM-index for fast MSV */
  P7_SCOREDATA     *scoredata;   /* hmm-specific data used by nhmmer */
  P7_HMM           *hmm;
} WORKER_INFO;

typedef struct {
  int    id;         /* internal sequence ID  */
  int    length;     /* length of sequence */
} ID_LENGTH;

typedef struct {
  ID_LENGTH  *id_lengths;
  int        count;
  int        size;
} ID_LENGTH_LIST;


static ID_LENGTH_LIST* init_id_length( int size );
static void            destroy_id_length( ID_LENGTH_LIST *list );
static int             add_id_length(ID_LENGTH_LIST *list, int id, int L);
static int             assign_Lengths(P7_TOPHITS *th, ID_LENGTH_LIST *id_length_list);

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#define CPUOPTS     NULL
#define MPIOPTS     NULL


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE,      NULL, NULL, NULL,    NULL,  NULL,  NULL,              "direct output to file <f>, not stdout",                        2 },
  { "-A",           eslARG_OUTFILE,      NULL, NULL, NULL,    NULL,  NULL,  NULL,              "save multiple alignment of all hits to file <s>",              2 },
  { "--tblout",     eslARG_OUTFILE,      NULL, NULL, NULL,    NULL,  NULL,  NULL,              "save parseable table of hits to file <s>",                     2 },
  { "--dfamtblout", eslARG_OUTFILE,      NULL, NULL, NULL,    NULL,  NULL,  NULL,              "save table of hits to file, in Dfam format <s>",               2 },
  { "--longinsertout", eslARG_OUTFILE,   NULL, NULL, NULL,    NULL,  NULL,  NULL,              "save table of long inserts to file <s>",                       2 },
  { "--insertlength", eslARG_INT,        "40", NULL, "n>=5",  NULL,  "--longinsertout",  NULL, "min length of insert printed to --longinsertout",              2 },
  { "--aliscoresout", eslARG_OUTFILE,   NULL, NULL, NULL,    NULL,  NULL,  NULL,               "save of scores for each position in each alignment to <s>",    2 },

  { "--acc",        eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                       2 },
  { "--noali",      eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
  { "--notextw",    eslARG_NONE,         NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                         2 },
  { "--textw",      eslARG_INT,         "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                     2 },
  /* Control of scoring system */
  { "--popen",      eslARG_REAL,       "0.02", NULL, "0<=x<0.5",NULL,  NULL,  NULL,          "gap open probability",                                         3 },
  { "--pextend",    eslARG_REAL,        "0.4", NULL, "0<=x<1",  NULL,  NULL,  NULL,          "gap extend probability",                                       3 },
/*  { "--mx",         eslARG_STRING, "BLOSUM62", NULL, NULL,      NULL,  NULL,  "--mxfile",    "substitution score matrix choice (of some built-in matrices)", 3 },*/
/*  { "--mxfile",     eslARG_INFILE,       NULL, NULL, NULL,      NULL,  NULL,  "--mx",        "read substitution score matrix from file <f>",                 3 },*/
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,       "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report sequences <= this E-value threshold in output",         4 },
  { "-T",           eslARG_REAL,        FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report sequences >= this score threshold in output",           4 },
  /* Control of inclusion (significance) thresholds */
  { "--incE",       eslARG_REAL,       "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",       eslARG_REAL,        FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider sequences >= this score threshold as significant",    5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's GA gathering cutoffs to set all thresholding",   6 },
  { "--cut_nc",     eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's NC noise cutoffs to set all thresholding",       6 },
  { "--cut_tc",     eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's TC trusted cutoffs to set all thresholding",     6 },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,        FALSE, NULL, NULL,    NULL,  NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",      7 },
  { "--F1",         eslARG_REAL,       "0.02", NULL, NULL,    NULL,  NULL, "--max",          "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F2",         eslARG_REAL,       "1e-3", NULL, NULL,    NULL,  NULL, "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F3",         eslARG_REAL,       "1e-5", NULL, NULL,    NULL,  NULL, "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--nobias",     eslARG_NONE,         NULL, NULL, NULL,    NULL,  NULL, "--max",          "turn off composition bias filter",                             7 },
  { "--B1",         eslARG_INT,         "110", NULL, NULL,    NULL,  NULL, "--max,--nobias", "window length for biased-composition modifier (MSV)",          7 },
  { "--B2",         eslARG_INT,         "240", NULL, NULL,    NULL,  NULL, "--max,--nobias", "window length for biased-composition modifier (Vit)",          7 },
  { "--B3",         eslARG_INT,        "1000", NULL, NULL,    NULL,  NULL, "--max,--nobias", "window length for biased-composition modifier (Fwd)",          7 },

  /* Control of FM pruning/extension */
  { "--fm_msv_length",   eslARG_INT,          "70", NULL, NULL,    NULL,  NULL, NULL,          "max length used when extending seed for MSV",                8 },
  { "--fm_max_depth",    eslARG_INT,          "16", NULL, NULL,    NULL,  NULL, NULL,          "seed length at which bit threshold must be met",             8 },
  { "--fm_max_neg_len",  eslARG_INT,           "4", NULL, NULL,    NULL,  NULL, NULL,          "maximum number consecutive negative scores in seed",         8 },
  { "--fm_req_pos",      eslARG_INT,           "5", NULL, NULL,    NULL,  NULL, NULL,          "minimum number consecutive positive scores in seed" ,        8 },
  { "--fm_sc_ratio",     eslARG_REAL,       "0.45", NULL, NULL,    NULL,  NULL, NULL,          "seed must maintain this bit ratio from one of two ends",     8 },
  { "--fm_max_scthresh", eslARG_REAL,       "10.5", NULL, NULL,    NULL,  NULL, NULL,          "max total bits required in seed of length fm_max_depth",     8 },

/* Other options */
  { "--tformat",    eslARG_STRING,       NULL, NULL, NULL,    NULL,  NULL,           NULL,     "assert target <seqdb> is in format <s>: no autodetection",      12 },
  { "--nonull2",    eslARG_NONE,         NULL, NULL, NULL,    NULL,  NULL,           NULL,     "turn off biased composition score corrections",                 12 },
  { "--usenull3",   eslARG_NONE,         NULL, NULL, NULL,    NULL,  NULL,        "--nonull2", "also use null3 for biased composition score corrections",       12 },
  { "-Z",           eslARG_REAL,        FALSE, NULL, "x>0",   NULL,  NULL,           NULL,     "set database size (Megabases) to <x> for E-value calculations", 12 },
  { "--seed",       eslARG_INT,          "42", NULL, "n>=0",  NULL,  NULL,           NULL,     "set RNG seed to <n> (if 0: one-time arbitrary seed)",           12 },
  { "--w_beta",     eslARG_REAL,         NULL, NULL, NULL,    NULL,  NULL,           NULL,     "tail mass at which window length is determined",                12 },
  { "--w_length",   eslARG_INT,          NULL, NULL, NULL,    NULL,  NULL,           NULL,     "window length - essentially max expected hit length ",                                                12 },
  { "--block_length", eslARG_INT,        NULL, NULL, "n>=50000", NULL, NULL,         NULL,     "length of blocks read from target database (threaded) ",        12 },
  { "--toponly",     eslARG_NONE,         NULL, NULL, NULL,    NULL,  NULL,   "--bottomonly",  "only search the top strand",                                    12 },
  { "--bottomonly",  eslARG_NONE,         NULL, NULL, NULL,    NULL,  NULL,      "--toponly",  "only search the bottom strand",                                 12 },


/* Not used, but retained because esl option-handling code errors if it isn't kept here.  Placed in group 99 so it doesn't print to help*/
  { "--domZ",       eslARG_REAL,        FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "Not used",   99 },
  { "--domE",       eslARG_REAL,       "10.0", NULL, "x>0",   NULL,  NULL,  DOMREPOPTS,      "Not used",   99 },
  { "--domT",       eslARG_REAL,        FALSE, NULL, NULL,    NULL,  NULL,  DOMREPOPTS,      "Not used",   99 },
  { "--incdomE",    eslARG_REAL,       "0.01", NULL, "x>0",   NULL,  NULL,  INCDOMOPTS,      "Not used",    99 },
  { "--incdomT",    eslARG_REAL,        FALSE, NULL, NULL,    NULL,  NULL,  INCDOMOPTS,      "Not used",      99 },
/* will eventually bring these back, but store in group 99 for now, so they don't print to help*/
  { "--qformat",    eslARG_STRING,       NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert query <seqfile> is in format <s>: no autodetection",   99 },



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
static char usage[]  = "[options] <hmmfile> <seqdb>";
static char banner[] = "search a DNA model against a DNA database";


static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_SQFILE *dbfp);
static int  serial_loop_FM(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_SQFILE *dbfp);

#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/


static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_hmmfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n", go->errbuf)  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n", go->errbuf)  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      if (puts("\nBasic options:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 100); /* 1= group; 2 = indentation; 120=textwidth*/

      if (puts("\nOptions directing output:")                                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 100);

//    if (puts("\nOptions controlling scoring system:")                      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
//     esl_opt_DisplayHelp(stdout, go, 3, 2, 80);

      if (puts("\nOptions controlling reporting thresholds:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 100);

      if (puts("\nOptions controlling inclusion (significance) thresholds:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 100);

      if (puts("\nOptions controlling model-specific thresholding:")         < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);

      if (puts("\nOptions controlling acceleration heuristics:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 100);

//      if (puts("\nOptions controlling trimming thresholds:")         < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
//      esl_opt_DisplayHelp(stdout, go, 9, 2, 100);

//      if (puts("\nControl of FM pruning and extension:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
//      esl_opt_DisplayHelp(stdout, go, 8, 2, 100);

      if (puts("\nOther expert options:")                                    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 100);
      exit(0);

  }

  if (esl_opt_ArgNumber(go)                  != 2)     { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { if (puts("Failed to get <hmmfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { if (puts("Failed to get <seqdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_hmmfile, "-") == 0 && strcmp(*ret_seqfile, "-") == 0) 
    { if (puts("Either <hmmfile> or <seqdb> may be '-' (to read from stdin), but not both.") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  *ret_go = go;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  if (puts("\nwhere basic options are:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *hmmfile, char *seqfile)
{
  p7_banner(ofp, go->argv[0], banner);
  
  if (fprintf(ofp, "# query HMM file:                  %s\n", hmmfile)                                                                                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# target sequence database:        %s\n", seqfile)                                                                                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-o")              && fprintf(ofp, "# output directed to file:         %s\n",            esl_opt_GetString(go, "-o"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-A")              && fprintf(ofp, "# MSA of all hits saved to file:   %s\n",            esl_opt_GetString(go, "-A"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--tblout")        && fprintf(ofp, "# hits tabular output:             %s\n",            esl_opt_GetString(go, "--tblout"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--dfamtblout")    && fprintf(ofp, "# hits output in Dfam format:      %s\n",            esl_opt_GetString(go, "--dfamtblout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--longinsertout") && fprintf(ofp, "# long inserts output:             %s\n",            esl_opt_GetString(go, "--longinsertout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--insertlength")  && fprintf(ofp, "# long inserts minimum length:     %d\n",            esl_opt_GetInteger(go, "--insertlength")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--aliscoresout")  && fprintf(ofp, "# alignment scores output:         %s\n",            esl_opt_GetString(go, "--longinsertout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--acc")        && fprintf(ofp, "# prefer accessions over names:    yes\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--noali")      && fprintf(ofp, "# show alignments in output:       no\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notextw")    && fprintf(ofp, "# max ASCII text line length:      unlimited\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--textw")      && fprintf(ofp, "# max ASCII text line length:      %d\n",             esl_opt_GetInteger(go, "--textw"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--popen")      && fprintf(ofp, "# gap open probability:            %f\n",             esl_opt_GetReal   (go, "--popen"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--pextend")    && fprintf(ofp, "# gap extend probability:          %f\n",             esl_opt_GetReal   (go, "--pextend"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
//if (esl_opt_IsUsed(go, "--mx")         && fprintf(ofp, "# subst score matrix (built-in):   %s\n",             esl_opt_GetString (go, "--mx"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
//if (esl_opt_IsUsed(go, "--mxfile")     && fprintf(ofp, "# subst score matrix (file):       %s\n",             esl_opt_GetString (go, "--mxfile"))   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-E")           && fprintf(ofp, "# sequence reporting threshold:    E-value <= %g\n",  esl_opt_GetReal   (go, "-E"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-T")           && fprintf(ofp, "# sequence reporting threshold:    score >= %g\n",    esl_opt_GetReal   (go, "-T"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incE")       && fprintf(ofp, "# sequence inclusion threshold:    E-value <= %g\n",  esl_opt_GetReal   (go, "--incE"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incT")       && fprintf(ofp, "# sequence inclusion threshold:    score >= %g\n",    esl_opt_GetReal   (go, "--incT"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_ga")     && fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_nc")     && fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_tc")     && fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--max")        && fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n")                       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F1")         && fprintf(ofp, "# MSV filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F1"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F2")         && fprintf(ofp, "# Vit filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F2"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F3")         && fprintf(ofp, "# Fwd filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F3"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nobias")     && fprintf(ofp, "# biased composition HMM filter:   off\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--B1")         && fprintf(ofp, "# biased comp MSV window len:      %d\n",             esl_opt_GetInteger(go, "--B1"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--B2")         && fprintf(ofp, "# biased comp Viterbi window len:  %d\n",             esl_opt_GetInteger(go, "--B2"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--B3")         && fprintf(ofp, "# biased comp Forward window len:  %d\n",             esl_opt_GetInteger(go, "--B3"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");


  if (esl_opt_IsUsed(go, "--fm_msv_length")   && fprintf(ofp, "# seed extension max length:       %d\n",             esl_opt_GetInteger(go, "--fm_msv_length"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fm_max_depth")    && fprintf(ofp, "# seed length:                     %d\n",             esl_opt_GetInteger(go, "--fm_max_depth"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fm_max_neg_len")  && fprintf(ofp, "# max consec neg score in seed:    %d\n",             esl_opt_GetInteger(go, "--fm_max_neg_len"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fm_req_pos")      && fprintf(ofp, "# min consec pos score in seed:    %d\n",             esl_opt_GetInteger(go, "--fm_req_pos"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fm_sc_ratio")     && fprintf(ofp, "# min seed bit ratio:              %.2f\n",           esl_opt_GetReal(go, "--fm_sc_ratio"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fm_max_scthresh") && fprintf(ofp, "# max bits in seed:                %.2f\n",           esl_opt_GetReal(go, "--fm_max_scthresh"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--nonull2")    && fprintf(ofp, "# null2 bias corrections:          off\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--usenull3")   && fprintf(ofp, "# null3 bias corrections:          on\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--toponly")    && fprintf(ofp, "# search only top strand:          on\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--bottomonly") && fprintf(ofp, "# search only bottom strand:       on\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-Z")           && fprintf(ofp, "# database size is set to:         %.1f Mb\n",        esl_opt_GetReal(go, "-Z"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0 && fprintf(ofp, "# random number seed:              one-time arbitrary\n")                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if                              (  fprintf(ofp, "# random number seed set to:       %d\n",             esl_opt_GetInteger(go, "--seed"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
//if (esl_opt_IsUsed(go, "--qformat")    && fprintf(ofp, "# query <seqfile> format asserted: %s\n",             esl_opt_GetString(go, "--qformat"))   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--tformat")    && fprintf(ofp, "# targ <seqfile> format asserted:  %s\n",             esl_opt_GetString(go, "--tformat"))   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--w_beta")     && fprintf(ofp, "# window length beta value:        %g\n",             esl_opt_GetReal(go, "--w_beta"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--w_length")   && fprintf(ofp, "# window length :                  %d\n",             esl_opt_GetInteger(go, "--w_length")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--block_length")&&fprintf(ofp, "# block length :                   %d\n",             esl_opt_GetInteger(go, "--block_length")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  #ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu")        && fprintf(ofp, "# number of worker threads:        %d\n",             esl_opt_GetInteger(go, "--cpu"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#endif
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}



int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go       = NULL;  
  struct cfg_s     cfg;         
  int              status   = eslOK;

  impl_Init();             /* processor specific initialization */
  p7_FLogsumInit();        /* we're going to use table-driven Logsum() approximations at times */

  /* Initialize what we can in the config structure (without knowing the alphabet yet)
   */
  cfg.hmmfile    = NULL;
  cfg.dbfile     = NULL;
  cfg.do_mpi     = FALSE;               /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;                   /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;                   /* this gets reset below, if we init MPI */

  process_commandline(argc, argv, &go, &cfg.hmmfile, &cfg.dbfile);    

  status = serial_master(go, &cfg);

  esl_getopts_Destroy(go);
  return status;
}

/* serial_master()
 * The serial version of hmmsearch.
 * For each query HMM in <hmmfile> search the database for hits.
 * 
 * A master can only return if it's successful. All errors are handled
 * immediately and fatally with p7_Fail().
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp          = stdout;          /* results output file (-o)                        */
  FILE            *afp          = NULL;            /* alignment output file (-A)                      */
  FILE            *tblfp        = NULL;            /* output stream for tabular  (--tblout)    */
  FILE            *dfamtblfp    = NULL;            /* output stream for tabular Dfam format (--dfamtblout)    */
  FILE            *inserttblfp  = NULL;            /* output stream for table of long inserts (--longinsertout)    */
  FILE            *aliscoresfp  = NULL;            /* output stream for alignment scores (--aliscoresout)    */


  P7_HMMFILE      *hfp        = NULL;              /* open input HMM file                             */
  P7_HMM          *hmm        = NULL;              /* one HMM query                                   */
  P7_SCOREDATA    *scoredata  = NULL;
  /*either the two above or ones below will be used*/
  //int              qhformat  = eslSQFILE_UNKNOWN;  /* format of qfile                                  */
  //ESL_SQFILE      *qfp      = NULL;          /* open qfile                                       */
  //ESL_SQ          *qsq      = NULL;               /* query sequence                                   */

  //int              dbformat = eslSQFILE_FMINDEX; // eslSQFILE_UNKNOWN;  /* format of dbfile                                 */
  int              dbformat =  eslSQFILE_UNKNOWN;  /* format of dbfile                                 */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */

  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  ESL_STOPWATCH   *w;
  //int              seed;
  int              textw    = 0;
  int              nquery   = 0;
  int              status   = eslOK;
  int              qhstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;
  double           resCnt = 0;
  /* used to keep track of the lengths of the sequences that are processed */
  ID_LENGTH_LIST  *id_length_list = NULL;


  /* these variables are only used if db type is FM-index*/
  void        *fm_cfg_mem      = NULL; //used to ensure cfg is 16-byte aligned, which matters since, for sse/vmx implementations, elements within cfg need to be aligned thusly
  FM_CFG      *fm_cfg       = NULL;
  FM_METADATA *fm_meta      = NULL;
  fpos_t       fm_basepos;
  /* end FM-index-specific variables */

  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char   errbuf[eslERRBUFSIZE];
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

  if (dbformat == eslSQFILE_FMINDEX) {
    //For now, this is a separate path from the typical esl_sqfile_Open() function call
    //TODO: create esl_sqio_fmindex.c, analogous to esl_sqio_ascii.c,

    fm_configAlloc(&fm_cfg_mem, &fm_cfg);
    fm_meta = fm_cfg->meta;

    if((fm_meta->fp = fopen(cfg->dbfile, "rb")) == NULL)
      esl_fatal("Cannot open file `%s': ", cfg->dbfile);

    fm_readFMmeta(fm_meta);
    fm_initConfig(fm_cfg, go);
    fm_createAlphabet(fm_meta, NULL); // don't override charBits
    fgetpos( fm_meta->fp, &fm_basepos);

  } else {
    /* Open the target sequence database */
    status = esl_sqfile_Open(cfg->dbfile, dbformat, p7_SEQDBENV, &dbfp);
    if      (status == eslENOTFOUND) p7_Fail("Failed to open target sequence database %s for reading\n",      cfg->dbfile);
    else if (status == eslEFORMAT)   p7_Fail("Target sequence database file %s is empty or misformatted\n",   cfg->dbfile);
    else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
    else if (status != eslOK)        p7_Fail("Unexpected error %d opening target sequence database file %s\n", status, cfg->dbfile);

    if (dbfp->format > 100) // breaking the law!  That range is reserved for msa, for aligned formats
      p7_Fail("%s contains a multiple sequence alignment; expect unaligned sequences, like FASTA\n",   cfg->dbfile);
  }

  /* Open the query profile HMM file */
  status = p7_hmmfile_OpenE(cfg->hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))              { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "-A"))              { if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
  if (esl_opt_IsOn(go, "--tblout"))        { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }
  if (esl_opt_IsOn(go, "--dfamtblout"))    { if ((dfamtblfp    = fopen(esl_opt_GetString(go, "--dfamtblout"),"w"))    == NULL)  esl_fatal("Failed to open tabular dfam output file %s for writing\n", esl_opt_GetString(go, "--dfamtblout")); }
  if (esl_opt_IsOn(go, "--longinsertout")) { if ((inserttblfp  = fopen(esl_opt_GetString(go, "--longinsertout"),"w")) == NULL)  esl_fatal("Failed to open long-insert output file %s for writing\n", esl_opt_GetString(go, "--longinsertout")); }
  if (esl_opt_IsOn(go, "--aliscoresout"))  { if ((aliscoresfp  = fopen(esl_opt_GetString(go, "--aliscoresout"),"w")) == NULL)  esl_fatal("Failed to open alignment scores output file %s for writing\n", esl_opt_GetString(go, "--aliscoresout")); }


#ifdef HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                   esl_threads_CPUCount(&ncpus);


  if (ncpus > 0) {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
  }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  /* <abc> is not known 'til first HMM is read.   Could be DNA or RNA*/
  qhstatus = p7_hmmfile_Read(hfp, &abc, &hmm);

  if (! (abc->type == eslRNA || abc->type == eslDNA))
    p7_Fail("Invalid alphabet type in hmm for nhmmer. Expect DNA or RNA\n");
  /* if using FM-index, verify that the alphabets are in agreement */
  if (dbformat == eslSQFILE_FMINDEX) {
    if ( ! (fm_meta->alph_type == fm_DNA ||
            fm_meta->alph_type == fm_DNA_full  ||
            fm_meta->alph_type == fm_RNA ||
            fm_meta->alph_type == fm_RNA_full  ) )
      p7_Fail("Alphabet type of hmm and database disagree.\n");
  }

  if (qhstatus == eslOK) {
      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, cfg->hmmfile, cfg->dbfile);

      if (dbformat != eslSQFILE_FMINDEX)
        dbfp->abc = abc;

      for (i = 0; i < infocnt; ++i)    {
          info[i].pli    = NULL;
          info[i].th     = NULL;
          info[i].om     = NULL;
          info[i].bg     = p7_bg_Create(abc);

          info[i].bg->use_null3  = esl_opt_IsUsed(go, "--usenull3");
          info[i].fm_cfg = NULL;
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
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */

      if      (window_length > 0)     hmm->max_length = window_length;
      else if (window_beta   > 0)     p7_Builder_MaxLength(hmm, window_beta);
      else if (hmm->max_length == -1 ) p7_Builder_MaxLength(hmm, p7_DEFAULT_WINDOW_BETA);

      nquery++;
      resCnt = 0;
      esl_stopwatch_Start(w);

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1) {
        if (dbformat == eslSQFILE_FMINDEX) {
          //rewind
          if (fsetpos(fm_meta->fp, &fm_basepos) != 0)  ESL_EXCEPTION(eslESYS, "rewind via fsetpos() failed");

        } else {
          if (! esl_sqfile_IsRewindable(dbfp))
            esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);
          esl_sqfile_Position(dbfp, 0);
        }
      }

      if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (hmm->acc  && fprintf(ofp, "Accession:   %s\n", hmm->acc)     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (hmm->desc && fprintf(ofp, "Description: %s\n", hmm->desc)    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");


      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_profile_ConfigLocal(gm, hmm, info->bg, 100); /* 100 is a dummy length for now; and MSVFilter requires local mode */
      p7_oprofile_Convert(gm, om);                    /* <om> is now p7_LOCAL, multihit */

      scoredata = p7_hmm_ScoreDataCreate(om, FALSE);

      for (i = 0; i < infocnt; ++i) {
          /* Create processing pipeline and hit list */
          info[i].th  = p7_tophits_Create();
          info[i].om = p7_oprofile_Copy(om);
          info[i].hmm = p7_hmm_Clone(hmm);
          info[i].pli = p7_pipeline_Create(go, om->M, 100, TRUE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
          p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);

          if (  esl_opt_IsUsed(go, "--toponly") )
            info[i].pli->strand = p7_STRAND_TOPONLY;
          else if (  esl_opt_IsUsed(go, "--bottomonly") )
            info[i].pli->strand = p7_STRAND_BOTTOMONLY;
          else
            info[i].pli->strand = p7_STRAND_BOTH;


          if (  esl_opt_IsUsed(go, "--block_length") )
            info[i].pli->block_length = esl_opt_GetInteger(go, "--block_length");
          else
            info[i].pli->block_length = NHMMER_MAX_RESIDUE_COUNT;

          info[i].fm_cfg = fm_cfg;
          info[i].scoredata = p7_hmm_ScoreDataClone(scoredata, om->abc->Kp);

#ifdef HMMER_THREADS
          if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
      }

      /* establish the id_lengths data structutre */
      id_length_list = init_id_length(1000);


#ifdef HMMER_THREADS
      if (dbformat == eslSQFILE_FMINDEX) {
        if (ncpus > 0)  sstatus = thread_loop(info, id_length_list, threadObj, queue, dbfp);
        else            sstatus = serial_loop_FM(info, id_length_list, dbfp);
      } else {
        if (ncpus > 0)  sstatus = thread_loop(info, id_length_list, threadObj, queue, dbfp);
        else            sstatus = serial_loop(info, id_length_list, dbfp);
      }
#else
      if (dbformat == eslSQFILE_FMINDEX)
        sstatus = serial_loop_FM(info, id_length_list, dbfp);
      else
        sstatus = serial_loop(info, id_length_list, dbfp);
#endif


      switch(sstatus) {
        case eslEFORMAT:
          esl_fatal("Parse failed (sequence file %s):\n%s\n",
                    dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
          break;
        case eslEOF:
        case eslOK:
          /* do nothing */
          break;
        default:
          esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
      }

      //need to re-compute e-values before merging (when list will be sorted)
      if (esl_opt_IsUsed(go, "-Z")) {
    	  resCnt = 1000000*esl_opt_GetReal(go, "-Z");

    	  if ( info[0].pli->strand == p7_STRAND_BOTH)
    	    resCnt *= 2;

      } else {
    	  for (i = 0; i < infocnt; ++i)
    		  resCnt += info[i].pli->nres;
      }

      for (i = 0; i < infocnt; ++i)
          p7_tophits_ComputeNhmmerEvalues(info[i].th, resCnt, info[i].om->max_length);


      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i) {
          p7_tophits_Merge(info[0].th, info[i].th);
          p7_pipeline_Merge(info[0].pli, info[i].pli);

          p7_pipeline_Destroy(info[i].pli);
          p7_tophits_Destroy(info[i].th);
          p7_oprofile_Destroy(info[i].om);
      }



      /* Print the results.  */
      p7_tophits_SortBySeqidxAndAlipos(info->th);
      assign_Lengths(info->th, id_length_list);
      p7_tophits_RemoveDuplicates(info->th, info->pli->use_bit_cutoffs);


      p7_tophits_SortBySortkey(info->th);
      p7_tophits_Threshold(info->th, info->pli);



      //tally up total number of hits and target coverage
      info->pli->n_output = info->pli->pos_output = 0;
      for (i = 0; i < info->th->N; i++) {
          if ( (info->th->hit[i]->flags & p7_IS_REPORTED) || info->th->hit[i]->flags & p7_IS_INCLUDED) {
              info->pli->n_output++;
              info->pli->pos_output += abs(info->th->hit[i]->dcl[0].jali - info->th->hit[i]->dcl[0].iali) + 1;
          }
      }


      p7_tophits_Targets(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      p7_tophits_Domains(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      if (tblfp)     p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, info->th, info->pli, (nquery == 1));
      if (dfamtblfp) p7_tophits_TabularXfam(dfamtblfp,   hmm->name, hmm->acc, info->th, info->pli);
      if (inserttblfp) p7_tophits_LongInserts(inserttblfp, hmm->name, hmm->acc, info->th, info->pli, esl_opt_GetInteger(go, "--insertlength") );
      if (aliscoresfp) p7_tophits_AliScores(aliscoresfp, hmm->name, info->th );


      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, info->pli, w);
      if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      /* Output the results in an MSA (-A option) */
      if (afp) {
          ESL_MSA *msa = NULL;

          if (p7_tophits_Alignment(info->th, abc, NULL, NULL, 0, p7_DEFAULT, &msa) == eslOK) {
            if (textw > 0) eslx_msafile_Write(afp, msa, eslMSAFILE_STOCKHOLM);
            else           eslx_msafile_Write(afp, msa, eslMSAFILE_PFAM);

            if (fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
          }  else {
	    if (fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
          }
          esl_msa_Destroy(msa);
      }

      for (i = 0; i < infocnt; ++i) {
        p7_hmm_ScoreDataDestroy(info[i].scoredata);
        //TODO: these cause errors, not sure why
        //p7_oprofile_Destroy(info[i].om);
        //p7_pipeline_Destroy(info[i].pli);
        //p7_tophits_Destroy(info[i].th);
      }

      p7_hmm_ScoreDataDestroy(scoredata);
      p7_pipeline_Destroy(info->pli);
      p7_tophits_Destroy(info->th);
      p7_oprofile_Destroy(om);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);
      destroy_id_length(id_length_list);

      qhstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  } /* end outer loop over query HMMs */

  switch(qhstatus) {
    case eslEOD:        p7_Fail("read failed, HMM file %s may be truncated?", cfg->hmmfile);      break;
    case eslEFORMAT:    p7_Fail("bad file format in HMM file %s",             cfg->hmmfile);      break;
    case eslEINCOMPAT:  p7_Fail("HMM file %s contains different alphabets",   cfg->hmmfile);      break;
    case eslEOF:        /* do nothing; EOF is what we expect here */                              break;
    default:            p7_Fail("Unexpected error (%d) in reading HMMs from %s", qhstatus, cfg->hmmfile);
  }

 /* Terminate outputs - any last words?
   */
  if (tblfp)    p7_tophits_TabularTail(tblfp,    "nhmmer", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (ofp)      { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for successful exit
   */
  for (i = 0; i < infocnt; ++i)
    p7_bg_Destroy(info[i].bg);

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

  if (dbformat == eslSQFILE_FMINDEX) {
    fm_destroyConfig(fm_cfg);
    free (fm_cfg->meta);
    free(fm_cfg_mem); //16-byte aligned memory in which cfg is found
  }

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  if (dfamtblfp)     fclose(dfamtblfp);
  if (inserttblfp)   fclose(inserttblfp);
  if (aliscoresfp)   fclose(aliscoresfp);

  return eslOK;

 ERROR:

   if (ofp != stdout) fclose(ofp);
   if (afp)           fclose(afp);
   if (tblfp)         fclose(tblfp);
   if (dfamtblfp)     fclose(dfamtblfp);
   if (inserttblfp)   fclose(inserttblfp);
   if (aliscoresfp)   fclose(aliscoresfp);

   return eslFAIL;
}

//TODO: MPI code needs to be added here
static int
serial_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_SQFILE *dbfp)
{

  int      wstatus;
  int i;
  int prev_hit_cnt;
  int seq_id = 0;
  P7_DOMAIN *dcl;
  ESL_SQ   *dbsq   =  esl_sq_CreateDigital(info->om->abc);
#ifdef eslAUGMENT_ALPHABET
  ESL_SQ   *dbsq_revcmp;
  if (dbsq->abc->complement != NULL)
    dbsq_revcmp =  esl_sq_CreateDigital(info->om->abc);
#endif /*eslAUGMENT_ALPHABET*/


  wstatus = esl_sqio_ReadWindow(dbfp, 0, info->pli->block_length, dbsq);


  while (wstatus == eslOK ) {
      dbsq->idx = seq_id;

      p7_pli_NewSeq(info->pli, dbsq);

      if (info->pli->strand != p7_STRAND_BOTTOMONLY) {

        info->pli->nres -= dbsq->C; // to account for overlapping region of windows
        prev_hit_cnt = info->th->N;

        p7_Pipeline_LongTarget(info->pli, info->om, info->scoredata, info->bg, dbsq, info->th, info->pli->nseqs);

        p7_pipeline_Reuse(info->pli); // prepare for next search

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
      } else {
        info->pli->nres -= dbsq->n;
      }
#ifdef eslAUGMENT_ALPHABET
      //reverse complement
      if (info->pli->strand != p7_STRAND_TOPONLY && dbsq->abc->complement != NULL )
      {
          prev_hit_cnt = info->th->N;
          esl_sq_Copy(dbsq,dbsq_revcmp);
          esl_sq_ReverseComplement(dbsq_revcmp);
          p7_Pipeline_LongTarget(info->pli, info->om, info->scoredata, info->bg, dbsq_revcmp, info->th, info->pli->nseqs);
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

      wstatus = esl_sqio_ReadWindow(dbfp, info->om->max_length, info->pli->block_length, dbsq);
      if (wstatus == eslEOD) { // no more left of this sequence ... move along to the next sequence.
          add_id_length(id_length_list, dbsq->idx, dbsq->L);

          info->pli->nseqs++;
          esl_sq_Reuse(dbsq);
          wstatus = esl_sqio_ReadWindow(dbfp, 0, info->pli->block_length, dbsq);

          seq_id++;

      }
    }

  if (dbsq) esl_sq_Destroy(dbsq);
  if (dbsq_revcmp) esl_sq_Destroy(dbsq_revcmp);

  return wstatus;

}

static int
serial_loop_FM(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_SQFILE *dbfp)
{

  int      wstatus = eslOK;
  int i;

  FM_DATA  fmf;
  FM_DATA  fmb;

  FM_METADATA *meta = info->fm_cfg->meta;

  for ( i=0; i<info->fm_cfg->meta->block_count; i++ ) {

    wstatus = fm_readFM( &fmf, meta, 1 );
    if (wstatus != eslOK) return wstatus;
    wstatus = fm_readFM( &fmb, meta, 0 );
    if (wstatus != eslOK) return wstatus;

    fmb.SA = fmf.SA;
    fmb.T  = fmf.T;

    wstatus = p7_Pipeline_FM(info->pli, info->om, info->scoredata, info->bg,
        info->th, info->pli->nseqs,  &fmf, &fmb, info->fm_cfg );
    if (wstatus != eslOK) return wstatus;

  }

  info->pli->nres = 2 * meta->char_count;

  return wstatus;

}

#ifdef HMMER_THREADS
static int
thread_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp)
{

  int i;
  int          status  = eslOK;
  int          sstatus = eslOK;
  int          eofCount = 0;
  ESL_SQ_BLOCK *block;
  void         *newBlock;
  int          seqid = -1;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
  ((ESL_SQ_BLOCK *)newBlock)->complete = TRUE;


  /* Main loop: */
  while (sstatus == eslOK ) {
      block = (ESL_SQ_BLOCK *) newBlock;

      sstatus = esl_sqio_ReadBlock(dbfp, block, info->pli->block_length, TRUE);


      block->first_seqidx = info->pli->nseqs;
      info->pli->nseqs += block->count  - (block->complete ? 0 : 1);// if there's an incomplete sequence read into the block wait to count it until it's complete.


      seqid = block->first_seqidx;
      for (i=0; i<block->count; i++) {
        block->list[i].idx = seqid;
        add_id_length(id_length_list, seqid, block->list[i].L);
        seqid++;
      }


      if (sstatus == eslEOF) {
          if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
          ++eofCount;
      }

      if (sstatus == eslOK) {
          status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
          if (status != eslOK) esl_fatal("Work queue reader failed");

          //newBlock needs all this information so the next ReadBlock call will know what to do
          ((ESL_SQ_BLOCK *)newBlock)->complete = block->complete;
          if (!block->complete) {
              // the final sequence on the block was an incomplete window of the active sequence,
              // so prep the next block to read in the next window
              esl_sq_Copy(block->list + (block->count - 1) , ((ESL_SQ_BLOCK *)newBlock)->list);

              if (  ((ESL_SQ_BLOCK *)newBlock)->list->n < info->om->max_length ) {
                //no reason to search the final partial sequence on the block, as the next block will search this whole chunk
                ((ESL_SQ_BLOCK *)newBlock)->list->C = ((ESL_SQ_BLOCK *)newBlock)->list->n;
                (((ESL_SQ_BLOCK *)newBlock)->count)--;
              } else {
                ((ESL_SQ_BLOCK *)newBlock)->list->C = info->om->max_length;
              }

          }
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
  P7_DOMAIN *dcl;
  ESL_SQ_BLOCK  *block = NULL;
  void          *newBlock;
  
  impl_ThreadInit();

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

      if (info->pli->strand != p7_STRAND_BOTTOMONLY) {
        info->pli->nres -= dbsq->C; // to account for overlapping region of windows

        prev_hit_cnt = info->th->N;
        p7_Pipeline_LongTarget(info->pli, info->om, info->scoredata, info->bg, dbsq, info->th, block->first_seqidx + i);
        p7_pipeline_Reuse(info->pli); // prepare for next search


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
      } else {
        info->pli->nres -= dbsq->n;
      }

#ifdef eslAUGMENT_ALPHABET
      //reverse complement
      if (info->pli->strand != p7_STRAND_TOPONLY && dbsq->abc->complement != NULL)
      {
          prev_hit_cnt = info->th->N;
          esl_sq_ReverseComplement(dbsq);
          p7_Pipeline_LongTarget(info->pli, info->om, info->scoredata, info->bg, dbsq, info->th, block->first_seqidx + i);
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



/* helper functions for tracking id_lengths */

static ID_LENGTH_LIST *
init_id_length( int size )
{
  int status;
  ID_LENGTH_LIST *list;

  ESL_ALLOC (list, sizeof(ID_LENGTH_LIST));
  list->count = 0;
  list->size  = size;
  list->id_lengths = NULL;

  ESL_ALLOC (list->id_lengths, size * sizeof(ID_LENGTH));

  return list;

ERROR:
  return NULL;
}

static void
destroy_id_length( ID_LENGTH_LIST *list )
{

  if (list != NULL) {
    if (list->id_lengths != NULL) free (list->id_lengths);
    free (list);
  }

}


static int
add_id_length(ID_LENGTH_LIST *list, int id, int L)
{
   int status;

   if (list->count > 0 && list->id_lengths[list->count-1].id == id) {
     // the last time this gets updated, it'll have the sequence's actual length
     list->id_lengths[list->count-1].length = L;
   } else {

     if (list->count == list->size) {
       list->size *= 10;
       ESL_REALLOC(list->id_lengths, list->size * sizeof(ID_LENGTH));
     }

     list->id_lengths[list->count].id     = id;
     list->id_lengths[list->count].length = L;

     list->count++;
   }
   return eslOK;

ERROR:
   return status;
}
 
static int
assign_Lengths(P7_TOPHITS *th, ID_LENGTH_LIST *id_length_list) {

  int i;
  int j = 0;

  for (i=0; i<th->N; i++) {
    while (th->hit[i]->seqidx != id_length_list->id_lengths[j].id) { j++;   }
    th->hit[i]->dcl[0].ad->L = id_length_list->id_lengths[j].length;
  }

  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/



