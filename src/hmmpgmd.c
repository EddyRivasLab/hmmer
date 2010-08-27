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

#ifdef HMMER_THREADS
#include <unistd.h>
#include <pthread.h>
#include "esl_threads.h"
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

static int read_QueryInfo(QUERY_INFO *info, FILE *qfp);
static int process_QueryInfo(QUERY_INFO *info);
static int destroy_QueryInfo(QUERY_INFO *info);
static int reuse_QueryInfo(QUERY_INFO *info);

typedef struct {
#ifdef HMMER_THREADS
  ESL_SQ           *sq_list;     /* list of sequences to process            */
  int               sq_cnt;      /* number of sequences                     */

  pthread_mutex_t   mutex;       /* protect data                            */
  int              *sq_inx;      /* next sequence to process                */
#endif /*HMMER_THREADS*/

  P7_HMM           *hmm;         /* query HMM                               */
  ESL_SQ           *seq;         /* query sequence                          */
  ESL_ALPHABET     *abc;         /* digital alphabet                        */
  ESL_GETOPTS      *opts;

  double            elapsed;

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
  { "--daemon",     eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  NULL,            "run program as a daemon",                                     12 },

#ifdef HMMER_THREADS 
  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>=0",NULL,  NULL,  NULL,            "number of parallel CPU workers to use for multithreads",      12 },
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

static char usage[]  = "[options] <target seqfile>";

static char banner[] = "search query against a sequence database";
static char banner_seq[] = "search a protein sequence against a sequence database";
static char banner_hmm[] = "search profile(s) against a sequence database";

static int  serial_loop  (WORKER_INFO *info, ESL_SQCACHE *cache);
#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(ESL_THREADS *obj);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/


static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go)
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

  if (esl_opt_ArgNumber(go) < 1)     { puts("Incorrect number of command line arguments.");      goto ERROR; }

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

  /* the options string can handle an optional database */
  if (esl_opt_ArgNumber(go) > 1)                  { puts("Incorrect number of command line arguments.");         goto ERROR; }

  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, pgm, usage);
  puts("\nwhere most common options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", pgm);
  exit(1);  
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *copt, const ESL_GETOPTS *sopt, char *seqfile, char *banner)
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

  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}

/* sort routines */
static int
sort_seq_size(const void *p1, const void *p2)
{
  int cmp;

  cmp  = (((ESL_SQ *)p1)->n < ((ESL_SQ *)p2)->n);
  cmp -= (((ESL_SQ *)p1)->n > ((ESL_SQ *)p2)->n);

  return cmp;
}

static int
sort_seq_inx(const void *p1, const void *p2)
{
  int cmp;

  cmp  = (((ESL_SQ *)p1)->idx > ((ESL_SQ *)p2)->idx);
  cmp -= (((ESL_SQ *)p1)->idx < ((ESL_SQ *)p2)->idx);

  return cmp;
}

int
main(int argc, char **argv)
{
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;              /* output stream for tabular per-seq (--domtblout) */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_ALPHABET    *abc;                          /* digital alphabet                                */
  ESL_STOPWATCH   *w;                            /* timer used for profiling statistics             */
  ESL_GETOPTS     *go       = NULL;              /* command line processing                         */
  int              textw    = 0;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i, j;

  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  ESL_THREADS     *threadObj= NULL;
  pthread_mutex_t  mutex;
  int              current_index;
#endif

  ESL_SQCACHE    **cache;
  QUERY_INFO      *query;

  /* Set processor specific flags */
  impl_Init();

  /* Initializations */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */

  process_commandline(argc, argv, &go);    

  w = esl_stopwatch_Create();

  abc = esl_alphabet_Create(eslAMINO);

  if (esl_opt_IsOn(go, "--tformat")) {
    dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* cache all the databases into memory */
  ESL_ALLOC(cache, sizeof(ESL_SQCACHE *) * esl_opt_ArgNumber(go));
  for (i = 0; i < esl_opt_ArgNumber(go); ++i) {
    ESL_RANDOMNESS *rnd  = NULL;
    ESL_SQ         *sq   = NULL;

    char *dbfile = esl_opt_GetArg(go, i+1);

    status = esl_sqfile_Cache(abc, dbfile, dbfmt, p7_SEQDBENV, &cache[i]);
    if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          dbfile);
    else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            dbfile);
    else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
    else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, dbfile);

    sq  = cache[i]->sq_list;

    /* sort the cache on sequence size */
    qsort(sq, cache[i]->seq_count, sizeof(ESL_SQ), sort_seq_size);

    /* jumble up the top 2/3 of the database.  This will leave the largest sequences at
     * the beginning of the cache.  then jumble up the bottom 1/3 of the database mixing
     * up the smaller sequences.  the reason is for load balancing the threads.  as we
     * process the database, smaller and smaller blocks of sequences will be processed
     * to try eleminate the case where one thread dominates the execution time.
     */
    rnd = esl_randomness_CreateFast(cache[i]->seq_count);
    for (j = 0 ; j < cache[i]->seq_count; ++j) {
      rnd->x = rnd->x * 69069 + 1;
      sq[j].idx = rnd->x;
    }
    esl_randomness_Destroy(rnd);

    j = cache[i]->seq_count / 3 * 2;
    qsort(sq, j, sizeof(ESL_SQ), sort_seq_inx);
    qsort(sq + j, cache[i]->seq_count - j, sizeof(ESL_SQ), sort_seq_inx);
    for (j = 0 ; j < cache[i]->seq_count; ++j) sq[j].idx = j;
  }

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
    if (pthread_mutex_init(&mutex, NULL) != 0) p7_Fail("mutex init failed");
  }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  /* initialize the query information struction */
  ESL_ALLOC(query, sizeof(QUERY_INFO));
  query->buf_size = 0;
  query->buffer   = NULL;
  query->pgm      = go->argv[0];
  query->opts     = esl_getopts_Create(searchOpts);

  query->abc      = abc;
  query->seq      = NULL;
  query->hmm      = NULL;

  /* read query hmm/sequence */
  while (read_QueryInfo(query, stdin) == eslOK) {
    int dbx;

    /* process any run-time options */
    if (process_QueryInfo(query) == eslEOF) break;

    /* figure out which cached database to use */
    dbx = 0;
    if (esl_opt_ArgNumber(query->opts) == 1) {
      char *db  = esl_opt_GetArg(query->opts, 1);
      int   len = strlen(db);
      for (i = 0; i < esl_opt_ArgNumber(go); ++i) {
        int n = strlen(cache[i]->filename);
        if (n >= len) {
          n = n - len;
          if (strcmp(cache[i]->filename + n, db) == 0) {
            dbx = i;
            break;
          }
        }
      }
      if (i >= esl_opt_ArgNumber(go)) {
        /* TODO report back an error that the db cannot be found */
      }
    }

    textw = (esl_opt_GetBoolean(query->opts, "--notextw")) ? 0 : esl_opt_GetInteger(query->opts, "--textw");

    esl_stopwatch_Start(w);

    if (query->hmm == NULL) {
      output_header(ofp, go, query->opts, cache[dbx]->filename, banner_seq);
      fprintf(ofp, "Query:       %s  [L=%ld]\n", query->seq->name, (long) query->seq->n);
      if (query->seq->acc[0]  != '\0') fprintf(ofp, "Accession:   %s\n", query->seq->acc);
      if (query->seq->desc[0] != '\0') fprintf(ofp, "Description: %s\n", query->seq->desc);  
    } else {
      output_header(ofp, go, query->opts, cache[dbx]->filename, banner_hmm);
      fprintf(ofp, "Query:       %s  [M=%d]\n", query->hmm->name, query->hmm->M);
      if (query->hmm->acc)  fprintf(ofp, "Accession:   %s\n", query->hmm->acc);
      if (query->hmm->desc) fprintf(ofp, "Description: %s\n", query->hmm->desc);
    }

    /* Create processing pipeline and hit list */
    for (i = 0; i < infocnt; ++i) {
      info[i].abc   = query->abc;
      info[i].hmm   = query->hmm;
      info[i].seq   = query->seq;
      info[i].opts  = query->opts;

      info[i].th    = NULL;
      info[i].pli   = NULL;

#ifdef HMMER_THREADS
      if (ncpus > 0) {
        info[i].sq_list = cache[dbx]->sq_list;
        info[i].sq_cnt  = cache[dbx]->seq_count;
        info[i].mutex   = mutex;
        info[i].sq_inx  = &current_index;

        esl_threads_AddThread(threadObj, &info[i]);
      }
#endif
    }

#ifdef HMMER_THREADS
    if (ncpus > 0) {
      current_index = 0;
      sstatus = thread_loop(threadObj);
    } else            
#endif
      sstatus = serial_loop(info, cache[dbx]);
    if (status != eslOK) {
      esl_fatal("Unexpected error %d reading sequence file %s", sstatus, cache[dbx]->filename);
    }

    /* merge the results of the search results */
    for (i = 1; i < infocnt; ++i) {
      p7_tophits_Merge(info[0].th, info[i].th);
      p7_pipeline_Merge(info[0].pli, info[i].pli);
    }

    /* Print the results.  */
    p7_tophits_Sort(info->th);
    p7_tophits_Threshold(info->th, info->pli);
    p7_tophits_Targets(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
    p7_tophits_Domains(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");

    if (tblfp)    p7_tophits_TabularTargets(tblfp,    query->hmm->name, query->hmm->acc, info->th, info->pli, TRUE);
    if (domtblfp) p7_tophits_TabularDomains(domtblfp, query->hmm->name, query->hmm->acc, info->th, info->pli, TRUE);
  
    esl_stopwatch_Stop(w);
    p7_pli_Statistics(ofp, info->pli, w);

#if 0
    fprintf (ofp, "   Sequences  Residues                        Elapsed\n");
    for (i = 1; i < infocnt; ++i) {
      char buf1[16];
      int h, m, s, hs;
      P7_PIPELINE *pli = info[i].pli;
      double elapsed;

      elapsed = info[i].elapsed;
      h  = (int) (elapsed / 3600.);
      m  = (int) (elapsed / 60.) - h * 60;
      s  = (int) (elapsed) - h * 3600 - m * 60;
      hs = (int) (elapsed * 100.) - h * 360000 - m * 6000 - s * 100;
      sprintf(buf1, "%02d:%02d.%02d", m,s,hs);

      fprintf (ofp, "%2d %9" PRId64 " %9" PRId64 " %7" PRId64 " %7" PRId64 " %6" PRId64 " %5" PRId64 " %s\n",
               i, pli->nseqs, pli->nres, pli->n_past_msv, pli->n_past_bias, pli->n_past_vit, pli->n_past_fwd, buf1);
    }
#endif

    for (i = 1; i < infocnt; ++i) {
      p7_pipeline_Destroy(info[i].pli);
      p7_tophits_Destroy(info[i].th);
    }
    fprintf(ofp, "//\n"); fflush(ofp);

    /* Output the results in an MSA (-A option) */
    if (afp) {
      ESL_MSA *msa = NULL;

      if (p7_tophits_Alignment(info->th, query->abc, NULL, NULL, 0, p7_ALL_CONSENSUS_COLS, &msa) == eslOK) {
        if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
        else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);
	  
        fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
      } 
      else fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n");
	  
      esl_msa_Destroy(msa);
    }

    p7_pipeline_Destroy(info->pli);
    p7_tophits_Destroy(info->th);

    reuse_QueryInfo(query);
  } /* end outer loop over query HMMs */

  if (hstatus != eslOK) p7_Fail("Unexpected error (%d) parsing input buffer", hstatus);

#ifdef HMMER_THREADS
  if (ncpus > 0) {
    pthread_mutex_destroy(&mutex);
    esl_threads_Destroy(threadObj);
  }
#endif

  for (i = 0; i < esl_opt_ArgNumber(go); ++i) {
    esl_sqfile_Free(cache[i]);
    cache[i] = NULL;
  }
  free(cache);

  destroy_QueryInfo(query);

  free(info);

  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);

  esl_getopts_Destroy(go);

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

  int           status  = eslOK;

  P7_HMMFILE   *hfp     = NULL;              /* open input HMM file                             */
  ESL_SQ       *seq     = NULL;              /* one target sequence (digital)                   */

  /* skip all leading white spaces */
  ptr = info->buffer;
  while (*ptr && isspace(*ptr)) ++ptr;

  /* process search specific options */
  if (*ptr == '@') {
    info->opts_ptr = ++ptr;

    /* skip to the end of the line */
    while (*ptr && (*ptr != '\n' && *ptr != '\r')) ++ptr;

    /* skip remaining white spaces */
    if (*ptr) {
      *ptr++ = 0;
      while (*ptr && isspace(*ptr)) ++ptr;
    }

    process_searchline(info->pgm, info->opts_ptr, info->opts);
  }

  if (*ptr) {
    info->query_ptr = ptr;

    if (strncmp(ptr, "//", 2) == 0) return eslEOF;

    /* try to parse the input buffer as a sequence */
    seq = esl_sq_CreateDigital(info->abc);
    status = esl_sqio_Parse(ptr, strlen(ptr), seq, eslSQFILE_DAEMON);
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

static int
serial_loop(WORKER_INFO *info, ESL_SQCACHE *cache)
{
  int               i;
  int               seed;
  int               status;
  ESL_SQ           *dbsq     = NULL;         /* one target sequence (digital)  */

  P7_BG            *bg       = NULL;	     /* null model                     */
  P7_PIPELINE      *pli      = NULL;         /* work pipeline                  */
  P7_TOPHITS       *th       = NULL;         /* top hit results                */
  P7_PROFILE       *gm       = NULL;         /* generic model                  */
  P7_OPROFILE      *om       = NULL;         /* optimized query profile        */

  P7_BUILDER       *bld      = NULL;         /* HMM construction configuration */

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
  for (i = 0; i < cache->seq_count; ++i) {
    dbsq = cache->sq_list + i;

    p7_pli_NewSeq(pli, dbsq);
    p7_bg_SetLength(bg, dbsq->n);
    p7_oprofile_ReconfigLength(om, dbsq->n);
      
    p7_Pipeline(pli, om, bg, dbsq, th);
	  
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

  return eslOK;
}

#ifdef HMMER_THREADS
static int
thread_loop(ESL_THREADS *obj)
{
  impl_Init();

  esl_threads_WaitForStart(obj);
  esl_threads_WaitForFinish(obj);

  return eslOK;
}

static void 
pipeline_thread(void *arg)
{
  int               i;
  int               seed;
  int               count;
  int               status;
  int               workeridx;
  WORKER_INFO      *info;
  ESL_THREADS      *obj;

  ESL_STOPWATCH    *w;

  P7_BG            *bg       = NULL;	     /* null model                     */
  P7_PIPELINE      *pli      = NULL;         /* work pipeline                  */
  P7_TOPHITS       *th       = NULL;         /* top hit results                */
  P7_PROFILE       *gm       = NULL;         /* generic model                  */
  P7_OPROFILE      *om       = NULL;         /* optimized query profile        */

  P7_BUILDER       *bld      = NULL;         /* HMM construction configuration */

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  w = esl_stopwatch_Create();
  esl_stopwatch_Start(w);

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

  /* loop until all sequences have been processed */
  count = 1;
  while (count > 0) {
    int     inx;
    ESL_SQ *dbsq;

#if 1
    /* grab the next block of sequences */
    if (pthread_mutex_lock(&info->mutex) != 0) p7_Fail("mutex lock failed");
    inx = *info->sq_inx;
    *info->sq_inx += BLOCK_SIZE;
    if (pthread_mutex_unlock(&info->mutex) != 0) p7_Fail("mutex unlock failed");

    dbsq = info->sq_list + inx;

    count = info->sq_cnt - inx;
    if (count > BLOCK_SIZE) count = BLOCK_SIZE;
    //printf("THREAD %08x: %d %d\n", workeridx, inx, count);
#endif

#if 0
    /* grab the next block of sequences */
    if (pthread_mutex_lock(&info->mutex) != 0) p7_Fail("mutex lock failed");
    inx = *info->sq_inx;
    count = info->sq_cnt - inx;
    count = count >> 7;
    if (count > 2500) count = 2500;
    if (count < 1000) count = 1000;
    *info->sq_inx += count;
    if (pthread_mutex_unlock(&info->mutex) != 0) p7_Fail("mutex unlock failed");

    dbsq = info->sq_list + inx;

    if (info->sq_cnt - inx < count) count = info->sq_cnt - inx;
    //printf("THREAD %08x: %d %d\n", workeridx, inx, count);
#endif

    /* Main loop: */
    for (i = 0; i < count; ++i, ++dbsq) {
      p7_pli_NewSeq(pli, dbsq);
      p7_bg_SetLength(bg, dbsq->n);
      p7_oprofile_ReconfigLength(om, dbsq->n);
	  
      p7_Pipeline(pli, om, bg, dbsq, th);
	  
      p7_pipeline_Reuse(pli);
    }
  }

  /* make available the pipeline objects to the main thread */
  info->th = th;
  info->pli = pli;

  /* clean up */
  p7_bg_Destroy(bg);
  p7_oprofile_Destroy(om);

  if (gm != NULL)  p7_profile_Destroy(gm);
  if (bld != NULL) p7_builder_Destroy(bld);

  esl_stopwatch_Stop(w);
  info->elapsed = w->elapsed;

  esl_stopwatch_Destroy(w);

  esl_threads_Finished(obj, workeridx);

  return;
}
#endif   /* HMMER_THREADS */
 

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

