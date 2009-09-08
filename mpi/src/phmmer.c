/* phmmer: search a protein sequence against a protein database
 * 
 * SRE, Sun Dec  7 11:19:05 2008 [Janelia] [Requiem for a Dream]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_scorematrix.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif

#include "hmmer.h"

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif
  P7_BG            *bg;
  P7_PIPELINE      *pli;
  P7_TOPHITS       *th;
  P7_OPROFILE      *om;
} WORKER_INFO;

#ifdef HMMER_THREADS
#define BLOCK_SIZE 2500

static int  threadedLoop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp);
static void pipelineThread(void *arg);

#else
static int serialLoop(WORKER_INFO *info, ESL_SQFILE *dbfp);
#endif

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles   reqs   incomp                             help                                                  docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,                          "show brief help on version and usage",                         1 },
/* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL,                          "direct output to file <f>, not stdout",                        7 },
  { "-A",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL,                          "save multiple alignment of hits to file <s>",                  7 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL,                          "save parseable table of per-sequence hits to file <s>",        7 },
  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL,                          "save parseable table of per-domain hits to file <s>",          7 },
/* Control of scoring system */
  { "--popen",      eslARG_REAL,  "0.02", NULL, "0<=x<0.5",NULL,  NULL,  NULL,                          "gap open probability",                                         2 },
  { "--pextend",    eslARG_REAL,   "0.4", NULL, "0<=x<1",  NULL,  NULL,  NULL,                          "gap extend probability",                                       2 },
  { "--mxfile",     eslARG_INFILE,  NULL, NULL, NULL,      NULL,  NULL,  NULL,                          "substitution score matrix [default: BLOSUM62]",                2 },
/* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "report sequences <= this E-value threshold in output",         3 },
  { "-T",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "report sequences >= this score threshold in output",           3 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                          "set # of comparisons done, for E-value calculation",           3 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "report domains <= this E-value threshold in output",           3 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",  "report domains >= this score cutoff in output",                3 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                          "set # of significant seqs, for domain E-value calculation",    3 },
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "-E,-T,--domE,--domT",         "use profile's GA gathering cutoffs to set -T, --domT",         99 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "-E,-T,--domE,--domT",         "use profile's NC noise cutoffs to set -T, --domT",             99 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "-E,-T,--domE,--domT",         "use profile's TC trusted cutoffs to set -T, --domT",           99 },
/* Control of inclusion thresholds */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include sequences <= this E-value threshold in output ali",    4 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include sequences >= this score threshold in output ali",      4 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include domains <= this E-value threshold in output ali",      4 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include domains >= this score threshold in output ali",        4 },
  { "--inc_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  "-A",  "--incE,--incT,--incdomE,--incdomT", "use profile's GA gathering cutoffs to set --incT, --incdomT",  99 },
  { "--inc_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  "-A",  "--incE,--incT,--incdomE,--incdomT", "use profile's NC noise cutoffs to set --incT, --incdomT",      99 },
  { "--inc_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  "-A",  "--incE,--incT,--incdomE,--incdomT", "use profile's TC trusted cutoffs to set --incT, --incdomT",    99 },
/* Control of filter pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--F1,--F2,--F3",               "Turn all heuristic filters off (less speed, more power)",      5 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,      NULL,  NULL, "--max",                        "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             5 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,      NULL,  NULL, "--max",                        "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             5 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,      NULL,  NULL, "--max",                        "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             5 },
  { "--nobias",     eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL, "--max",                        "turn off composition bias filter",                             5 },
  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL,  NULL,                          "turn off biased composition score corrections",                5 },
/* Control of E-value calibration */
  { "--EmL",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,                          "length of sequences for MSV Gumbel mu fit",                    6 },   
  { "--EmN",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,                          "number of sequences for MSV Gumbel mu fit",                    6 },   
  { "--EvL",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,                          "length of sequences for Viterbi Gumbel mu fit",                6 },   
  { "--EvN",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,                          "number of sequences for Viterbi Gumbel mu fit",                6 },   
  { "--EfL",        eslARG_INT,    "100", NULL,"n>0",      NULL,  NULL,  NULL,                          "length of sequences for Forward exp tail tau fit",             6 },   
  { "--EfN",        eslARG_INT,    "200", NULL,"n>0",      NULL,  NULL,  NULL,                          "number of sequences for Forward exp tail tau fit",             6 },   
  { "--Eft",        eslARG_REAL,  "0.04", NULL,"0<x<1",    NULL,  NULL,  NULL,                          "tail mass for Forward exponential tail tau fit",               6 },   
/* other options */
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",    NULL,  NULL,    NULL,                        "set RNG seed to <n> (if 0: one-time arbitrary seed)",       8 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",  NULL,  NULL,  "--notextw",                   "set max width of ASCII text output lines",                  8 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,      NULL,  NULL,  "--textw",                     "unlimit ASCII text output line width",                      8 },
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,      NULL,  NULL,   NULL,                         "assert query <seqfile> is in format <s>: no autodetection", 8 },
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,      NULL,  NULL,   NULL,                         "assert target <seqdb> is in format <s>>: no autodetection", 8 },
#ifdef HMMER_THREADS
  { "--cpu",        eslARG_INT,     NULL,"HMMER_NCPU", "n>0",      NULL,    NULL,   NULL,                       "number of worker threads",                                  8 },
#endif
#ifdef HAVE_MPI
  // { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "arrest after start: for debugging MPI under gdb",          4 },  
  // { "--mpi",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "run as an MPI parallel program",                           4 },
#endif 
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <query seqfile> <target seqdb>";
static char banner[] = "search a protein sequence against a protein database";



/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static void 
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_qfile, char **ret_dbfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     esl_fatal("Internal failure creating options object");
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/

      puts("\nOptions directing output:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      puts("\nOptions controlling scoring system:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      puts("\nOptions controlling significance thresholds for reporting:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 

      puts("\nOptions controlling significance thresholds for inclusion in output alignment (-A):");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      puts("\nOptions controlling acceleration heuristics:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      puts("\nOptions controlling E value calibration:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 

      puts("\nOther options:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)    { puts("Incorrect number of command line arguments.");    goto ERROR; }
  if ((*ret_qfile  = esl_opt_GetArg(go, 1)) == NULL) { puts("Failed to get <qfile> argument on command line"); goto ERROR; }
  if ((*ret_dbfile = esl_opt_GetArg(go, 2)) == NULL) { puts("Failed to get <seqdb> argument on command line"); goto ERROR; }

  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}


static int
output_header(FILE *ofp, ESL_GETOPTS *go, char *qfile, char *dbfile)
{
  p7_banner(ofp, go->argv[0], banner);
  
  fprintf(ofp, "# query sequence file:             %s\n", qfile);
  fprintf(ofp, "# target sequence database:        %s\n", dbfile);
  if (esl_opt_IsUsed(go, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",      esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "-A"))          fprintf(ofp, "# MSA of hits saved to file:       %s\n",      esl_opt_GetString(go, "-A"));
  if (esl_opt_IsUsed(go, "--tblout"))    fprintf(ofp, "# per-seq hits tabular output:     %s\n",      esl_opt_GetString(go, "--tblout"));
  if (esl_opt_IsUsed(go, "--domtblout")) fprintf(ofp, "# per-dom hits tabular output:     %s\n",      esl_opt_GetString(go, "--domtblout"));
  if (esl_opt_IsUsed(go, "--popen"))     fprintf(ofp, "# gap open probability:            %f\n",      esl_opt_GetReal  (go, "--popen"));
  if (esl_opt_IsUsed(go, "--pextend"))   fprintf(ofp, "# gap extend probability:          %f\n",      esl_opt_GetReal  (go, "--pextend"));
  if (esl_opt_IsUsed(go, "--mxfile"))    fprintf(ofp, "# subst score matrix:              %s\n",      esl_opt_GetString(go, "--mxfile"));
  if (esl_opt_IsUsed(go, "-E"))          fprintf(ofp, "# sequence E-value threshold:   <= %g\n",      esl_opt_GetReal(go, "-E"));
  if (esl_opt_IsUsed(go, "-T"))          fprintf(ofp, "# sequence bit score threshold: >= %g\n",      esl_opt_GetReal(go, "-T"));
  if (esl_opt_IsUsed(go, "-Z"))          fprintf(ofp, "# sequence search space set to:    %.0f\n",    esl_opt_GetReal(go, "-Z"));
  if (esl_opt_IsUsed(go, "--domE"))      fprintf(ofp, "# domain E-value threshold:     <= %g\n",      esl_opt_GetReal(go, "--domE"));
  if (esl_opt_IsUsed(go, "--domT"))      fprintf(ofp, "# domain bit score threshold:   >= %g\n",      esl_opt_GetReal(go, "--domT"));
  if (esl_opt_IsUsed(go, "--domZ"))      fprintf(ofp, "# domain search space set to:      %.0f\n",    esl_opt_GetReal(go, "--domZ"));
//if (esl_opt_IsUsed(go, "--cut_ga"))    fprintf(ofp, "# set reporting thresholds to:     GA cutoffs\n"); 
//if (esl_opt_IsUsed(go, "--cut_nc"))    fprintf(ofp, "# set reporting thresholds to:     NC cutoffs\n"); 
//if (esl_opt_IsUsed(go, "--cut_tc"))    fprintf(ofp, "# set reporting thresholds to:     TC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--incE"))      fprintf(ofp, "# seq inclusion E-val thresh:   <= %g\n",      esl_opt_GetReal(go, "--incE"));
  if (esl_opt_IsUsed(go, "--incT"))      fprintf(ofp, "# seq inclusion score thresh:   >= %g\n",      esl_opt_GetReal(go, "--incT"));
  if (esl_opt_IsUsed(go, "--incdomE"))   fprintf(ofp, "# dom inclusion E-val thresh:   <= %g\n",      esl_opt_GetReal(go, "--incdomE"));
  if (esl_opt_IsUsed(go, "--incdomT"))   fprintf(ofp, "# dom inclusion score thresh:   >= %g\n",      esl_opt_GetReal(go, "--incdomT"));
//if (esl_opt_IsUsed(go, "--inc_ga"))    fprintf(ofp, "# set inclusion thresholds to:     GA cutoffs\n"); 
//if (esl_opt_IsUsed(go, "--inc_nc"))    fprintf(ofp, "# set inclusion thresholds to:     NC cutoffs\n"); 
//if (esl_opt_IsUsed(go, "--inc_tc"))    fprintf(ofp, "# set inclusion thresholds to:     TC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--max"))       fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n");
  if (esl_opt_IsUsed(go, "--F1"))        fprintf(ofp, "# MSV filter P threshold:       <= %g\n",      esl_opt_GetReal(go, "--F1"));
  if (esl_opt_IsUsed(go, "--F2"))        fprintf(ofp, "# Vit filter P threshold:       <= %g\n",      esl_opt_GetReal(go, "--F2"));
  if (esl_opt_IsUsed(go, "--F3"))        fprintf(ofp, "# Fwd filter P threshold:       <= %g\n",      esl_opt_GetReal(go, "--F3"));
  if (esl_opt_IsUsed(go, "--nobias"))    fprintf(ofp, "# biased composition HMM filter:   off\n");
  if (esl_opt_IsUsed(go, "--nonull2"))   fprintf(ofp, "# null2 bias corrections:          off\n");
  if (esl_opt_IsUsed(go, "--EmL") )      fprintf(ofp, "# seq length, MSV Gumbel mu fit:   %d\n",     esl_opt_GetInteger(go, "--EmL"));
  if (esl_opt_IsUsed(go, "--EmN") )      fprintf(ofp, "# seq number, MSV Gumbel mu fit:   %d\n",     esl_opt_GetInteger(go, "--EmN"));
  if (esl_opt_IsUsed(go, "--EvL") )      fprintf(ofp, "# seq length, Vit Gumbel mu fit:   %d\n",     esl_opt_GetInteger(go, "--EvL"));
  if (esl_opt_IsUsed(go, "--EvN") )      fprintf(ofp, "# seq number, Vit Gumbel mu fit:   %d\n",     esl_opt_GetInteger(go, "--EvN"));
  if (esl_opt_IsUsed(go, "--EfL") )      fprintf(ofp, "# seq length, Fwd exp tau fit:     %d\n",     esl_opt_GetInteger(go, "--EfL"));
  if (esl_opt_IsUsed(go, "--EfN") )      fprintf(ofp, "# seq number, Fwd exp tau fit:     %d\n",     esl_opt_GetInteger(go, "--EfN"));
  if (esl_opt_IsUsed(go, "--Eft") )      fprintf(ofp, "# tail mass for Fwd exp tau fit:   %f\n",     esl_opt_GetReal   (go, "--Eft"));
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0) fprintf(ofp, "# random number seed:              one-time arbitrary\n");
    else                                       fprintf(ofp, "# random number seed set to:       %d\n", esl_opt_GetInteger(go, "--seed"));
  }
  if (esl_opt_IsUsed(go, "--textw"))     fprintf(ofp, "# max ASCII text line length:      %d\n",     esl_opt_GetInteger(go, "--textw"));
  if (esl_opt_IsUsed(go, "--notextw"))   fprintf(ofp, "# max ASCII text line length:      unlimited\n");
  if (esl_opt_IsUsed(go, "--qformat"))   fprintf(ofp, "# query <seqfile> format asserted: %s\n",     esl_opt_GetString(go, "--qformat"));
  if (esl_opt_IsUsed(go, "--tformat"))   fprintf(ofp, "# target <seqdb> format asserted:  %s\n",     esl_opt_GetString(go, "--tformat"));
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go       = NULL;               /* application configuration options                */
  FILE            *ofp      = stdout;             /* output file for results (default stdout)         */
  FILE            *afp      = NULL;               /* alignment output file (-A option)                */
  FILE            *tblfp    = NULL;		  /* output stream for tabular per-seq (--tblout)     */
  FILE            *domtblfp = NULL;		  /* output stream for tabular per-seq (--domtblout)  */
  char            *qfile    = NULL;               /* file to read query sequence from                 */
  int              qformat  = eslSQFILE_UNKNOWN;  /* format of qfile                                  */
  ESL_SQFILE      *qfp      = NULL;		  /* open qfile                                       */
  ESL_SQ          *qsq      = NULL;               /* query sequence                                   */
  char            *dbfile   = NULL;               /* file to read sequence(s) from                    */
  int              dbformat = eslSQFILE_UNKNOWN;  /* format of dbfile                                 */
  ESL_SQFILE      *dbfp     = NULL;               /* open dbfile                                      */
  ESL_SQ          *dbsq     = NULL;               /* target sequence                                  */
  ESL_ALPHABET    *abc      = NULL;               /* sequence alphabet                                */
  P7_BUILDER      *bld      = NULL;               /* HMM construction configuration                   */
  ESL_STOPWATCH   *w        = NULL;               /* for timing                                       */
  int              nquery   = 0;
  int              seed;
  int              textw;
  int              status   = eslOK;
  int              qstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;

  int              ncpus    = 1;

  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif

  /* Initializations */
  process_commandline(argc, argv, &go, &qfile, &dbfile);    
  p7_FLogsumInit();
  abc     = esl_alphabet_Create(eslAMINO);
  w       = esl_stopwatch_Create();
  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");
  esl_stopwatch_Start(w);

  /* If caller declared input formats, decode them */
  if (esl_opt_IsOn(go, "--qformat")) {
    qformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (qformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }
  if (esl_opt_IsOn(go, "--tformat")) {
    dbformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Initialize a default builder configuration,
   * then set only the options we need for single sequence search
   */
  bld = p7_builder_Create(NULL, abc);
  if ((seed = esl_opt_GetInteger(go, "--seed")) > 0)
    {				/* a little wasteful - we're blowing a couple of usec by reinitializing */
      esl_randomness_Init(bld->r, seed);
      bld->do_reseeding = TRUE;
    }
  bld->EmL = esl_opt_GetInteger(go, "--EmL");
  bld->EmN = esl_opt_GetInteger(go, "--EmN");
  bld->EvL = esl_opt_GetInteger(go, "--EvL");
  bld->EvN = esl_opt_GetInteger(go, "--EvN");
  bld->EfL = esl_opt_GetInteger(go, "--EfL");
  bld->EfN = esl_opt_GetInteger(go, "--EfN");
  bld->Eft = esl_opt_GetReal   (go, "--Eft");
  status = p7_builder_SetScoreSystem(bld, esl_opt_GetString(go, "--mxfile"), NULL, esl_opt_GetReal(go, "--popen"), esl_opt_GetReal(go, "--pextend"));
  if (status != eslOK) esl_fatal("Failed to set single query seq score system:\n%s\n", bld->errbuf);


  /* Open results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"),          "w")) == NULL)  esl_fatal("Failed to open output file %s for writing\n",                 esl_opt_GetString(go, "-o")); } 
  if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"),          "w")) == NULL)  esl_fatal("Failed to open alignment output file %s for writing\n",       esl_opt_GetString(go, "-A")); } 
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp")); }
  if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblfp")); }
    

  /* Open the target sequence database for sequential access. */
  status =  esl_sqfile_OpenDigital(abc, dbfile, dbformat, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open target sequence database %s for reading\n",      dbfile);
  else if (status == eslEFORMAT)   esl_fatal("Target sequence database file %s is empty or misformatted\n",   dbfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        esl_fatal("Unexpected error %d opening target sequence database file %s\n", status, dbfile);
  dbsq = esl_sq_CreateDigital(abc);

  /* Open the query sequence file  */
  status = esl_sqfile_OpenDigital(abc, qfile, qformat, NULL, &qfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open sequence file %s for reading\n",      qfile);
  else if (status == eslEFORMAT)   esl_fatal("Sequence file %s is empty or misformatted\n",        qfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        esl_fatal ("Unexpected error %d opening sequence file %s\n", status, qfile);
  qsq  = esl_sq_CreateDigital(abc);

#ifdef HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                           esl_threads_CPUCount(&ncpus);

  threadObj = esl_threads_Create(&pipelineThread);
  queue = esl_workqueue_Create(ncpus * 2);
#else
  ncpus = 1;
#endif

  ESL_ALLOC(info, sizeof(*info) * ncpus);

  /* Show header output */
  output_header(ofp, go, qfile, dbfile);

  for (i = 0; i < ncpus; ++i)
    {
      info[i].pli   = NULL;
      info[i].th    = NULL;
      info[i].om    = NULL;
      info[i].bg    = p7_bg_Create(abc);
#ifdef HMMER_THREADS
      info[i].queue = queue;
#endif
    }

#ifdef HMMER_THREADS
  for (i = 0; i < ncpus * 2; ++i)
    {
      block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abc);
      if (block == NULL) 
	{
	  esl_fatal("Failed to allocate sequence block");
	}

      status = esl_workqueue_Init(queue, block);
      if (status != eslOK) 
	{
	  esl_fatal("Failed to add block to work queue");
	}
    }
#endif

  /* Outer loop over sequence queries */
  while ((qstatus = esl_sqio_Read(qfp, qsq)) == eslOK)
    {
      P7_OPROFILE     *om       = NULL;           /* optimized query profile                  */

      nquery++;
      if (qsq->n == 0) continue; /* skip zero length seqs as if they aren't even present */

      esl_stopwatch_Start(w);

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1)
	{
	  if (! esl_sqfile_IsRewindable(dbfp)) esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", dbfile);
	  esl_sqfile_Position(dbfp, 0);
	}

      fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n);
      if (qsq->acc[0]  != '\0') fprintf(ofp, "Accession:   %s\n", qsq->acc);
      if (qsq->desc[0] != '\0') fprintf(ofp, "Description: %s\n", qsq->desc);  


      /* Build the model */
      p7_SingleBuilder(bld, qsq, info->bg, NULL, NULL, NULL, &om); /* bypass HMM - only need model */

      for (i = 0; i < ncpus; ++i)
	{
	  /* Create processing pipeline and hit list */
	  info[i].th  = p7_tophits_Create(); 
	  info[i].om  = p7_oprofile_Clone(om);
	  info[i].pli = p7_pipeline_Create(go, om->M, 100, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
	  p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);

#ifdef HMMER_THREADS
	  esl_threads_AddThread(threadObj, &info[i]);
#endif
	}

#if 0
      /* Run each target sequence through the pipeline */
      while ((sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
	{ 
	  p7_pli_NewSeq(pli, dbsq);
	  p7_bg_SetLength(bg, dbsq->n);
	  p7_oprofile_ReconfigLength(om, dbsq->n);
  
	  p7_Pipeline(pli, om, bg, dbsq, th);

	  esl_sq_Reuse(dbsq);
	  p7_pipeline_Reuse(pli);
	}
      if      (sstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s line %" PRId64 "):\n%s\n",
						dbfp->filename, dbfp->linenumber, dbfp->errbuf);     
      else if (sstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
						    sstatus, dbfp->filename);
#endif
#ifdef HMMER_THREADS
      sstatus = threadedLoop(threadObj, queue, dbfp);
#else
      sstatus = serialLoop(info, dbfp);
#endif
      switch(sstatus)
	{
	case eslEFORMAT: 
	  esl_fatal("Parse failed (sequence file %s line %" PRId64 "):\n%s\n",
		    dbfp->filename, dbfp->linenumber, dbfp->errbuf);
	  break;
	case eslEOF:
	  /* do nothing */
	  break;
	default:
	  esl_fatal("Unexpected error %d reading sequence file %s",
		    sstatus, dbfp->filename);
	}


      /* merge the results of the search results */
      for (i = 1; i < ncpus; ++i)
	{
	  p7_tophits_Merge(info[0].th, info[i].th);
	  p7_pipeline_Merge(info[0].pli, info[i].pli);

	  p7_pipeline_Destroy(info[i].pli);
	  p7_tophits_Destroy(info[i].th);
	  p7_oprofile_Destroy(info[i].om);
	}

      /* Print the results.  */
      p7_tophits_Sort(info->th);
      p7_tophits_Threshold(info->th, info->pli);
      p7_tophits_Targets(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
      p7_tophits_Domains(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
  
      if (tblfp)    p7_tophits_TabularTargets(tblfp,    qsq->name, info->th, info->pli, (nquery == 1));
      if (domtblfp) p7_tophits_TabularDomains(domtblfp, qsq->name, info->th, info->pli, (nquery == 1));

      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, info->pli, w);
      fprintf(ofp, "//\n");

      /* Output the results in an MSA (-A option) */
      if (afp) {
	ESL_MSA *msa = NULL;

	if ( p7_tophits_Alignment(info->th, abc, NULL, NULL, 0, p7_DEFAULT, &msa) == eslOK) 
	  {
	    if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
	    else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);

	    fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
	  }
	else fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n");
	  
	esl_msa_Destroy(msa);
      }

      p7_tophits_Destroy(info->th);
      p7_pipeline_Destroy(info->pli);
      p7_oprofile_Destroy(info->om);
      p7_oprofile_Destroy(om);
      esl_sq_Reuse(qsq);
    } /* end outer loop over query sequences */
  if      (qstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s line %" PRId64 "):\n%s\n",
					    qfp->filename, qfp->linenumber, qfp->errbuf);     
  else if (qstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    qstatus, qfp->filename);

  for (i = 0; i < ncpus; ++i)
    {
      p7_bg_Destroy(info[i].bg);
    }

#ifdef HMMER_THREADS
  esl_workqueue_Reset(queue);
  while ((block = esl_workqueue_Remove(queue)) != NULL)
    {
      esl_sq_DestroyBlock(block);
    }
  esl_workqueue_Destroy(queue);
  esl_threads_Destroy(threadObj);
#endif

  free(info);

  esl_sqfile_Close(dbfp);
  esl_sqfile_Close(qfp);
  esl_stopwatch_Destroy(w);
  esl_sq_Destroy(dbsq);
  esl_sq_Destroy(qsq);
  p7_builder_Destroy(bld);
  esl_alphabet_Destroy(abc);

  if (ofp      != stdout) fclose(ofp);
  if (afp      != NULL)   fclose(afp);
  if (tblfp    != NULL)   fclose(tblfp);
  if (domtblfp != NULL)   fclose(domtblfp);
  esl_getopts_Destroy(go);
  return eslOK;

 ERROR:
  return eslFAIL;
}

#ifdef HMMER_THREADS
static int
threadedLoop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp)
{
  int  status  = eslOK;
  int  sstatus = eslOK;
  int  eofCount = 0;
  ESL_SQ_BLOCK *block;
  void         *newBlock;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
      
  /* Main loop: */
  while (sstatus == eslOK)
    {
      block = (ESL_SQ_BLOCK *) newBlock;
      sstatus = esl_sqio_ReadBlock(dbfp, block);
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
pipelineThread(void *arg)
{
  int i;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  ESL_SQ_BLOCK  *block = NULL;
  void          *newBlock;
  
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
	  p7_bg_SetLength(info->bg, dbsq->n);
	  p7_oprofile_ReconfigLength(info->om, dbsq->n);
	  
	  p7_Pipeline(info->pli, info->om, info->bg, dbsq, info->th);
	  
	  esl_sq_Reuse(dbsq);
	  p7_pipeline_Reuse(info->pli);
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
#else
static int
serialLoop(WORKER_INFO *info, ESL_SQFILE *dbfp)
{
  int      sstatus;
  ESL_SQ   *dbsq     = NULL;   /* one target sequence (digital)  */

  dbsq = esl_sq_CreateDigital(info->om->abc);

  /* Main loop: */
  while ((sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
    {
      p7_pli_NewSeq(info->pli, dbsq);
      p7_bg_SetLength(info->bg, dbsq->n);
      p7_oprofile_ReconfigLength(info->om, dbsq->n);
      
      p7_Pipeline(info->pli, info->om, info->bg, dbsq, info->th);
	  
      esl_sq_Reuse(dbsq);
      p7_pipeline_Reuse(info->pli);
    }

  esl_sq_Destroy(dbsq);

  return sstatus;
}
#endif


/*****************************************************************
 * @LICENSE@
 *****************************************************************/
