/* hmmpgmd: hmmer deamon searchs against a sequence database.
 * 
 * MSF, Thu Aug 12, 2010 [Janelia]
 * SVN $Id: hmmsearch.c 3324 2010-07-07 19:30:12Z wheelert $
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>
#include <pthread.h>
#include <setjmp.h>
#include <sys/socket.h>
#include <arpa/inet.h>

#ifndef HMMER_THREADS
#error "Program requires pthreads be enabled."
#endif /*HMMER_THREADS*/

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"
#include "esl_threads.h"

#include "hmmer.h"
#include "hmmpgmd.h"

#define MAX_BUFFER 4096

typedef struct queue_data_s {
  P7_HMM             *hmm;         /* query HMM                      */
  ESL_SQ             *seq;         /* query sequence                 */
  ESL_ALPHABET       *abc;         /* digital alphabet               */
  ESL_GETOPTS        *opts;        /* search specific options        */
  HMMD_SEARCH_CMD    *cmd;         /* workers search command         */

  int                 sock;        /* socket descriptor of client    */

  int                 retry;
  int                 retry_cnt;

  int                 dbx;         /* database index to search       */
  int                 sq_cnt;
  int                 sq_inx;

  struct queue_data_s *next;
  struct queue_data_s *prev;
} QUEUE_DATA;

typedef struct {
  pthread_mutex_t  queue_mutex;
  pthread_cond_t   queue_cond;
  QUEUE_DATA      *head;
  QUEUE_DATA      *tail;
} SEARCH_QUEUE;

typedef struct {
  char           *pgm;
  int             sock_fd;
  SEARCH_QUEUE   *queue;
  ESL_SQCACHE   **cache;
  int             cache_size;
} CLIENTSIDE_ARGS;

typedef struct {
  char            *pgm;
  int              sock_fd;

  pthread_mutex_t  work_mutex;
  pthread_cond_t   start_cond;
  pthread_cond_t   worker_cond;
  pthread_cond_t   complete_cond;

  int              started;
  struct worker_s *head;
  struct worker_s *tail;

  int              completed;
  int              errors;
} WORKERSIDE_ARGS;

typedef struct worker_s {
  int                   sock_fd;
  
  HMMD_SEARCH_STATS     stats;
  HMMD_SEARCH_STATUS    status;
  char                 *err_buf;
  P7_HIT              **hit;
  void                **hit_data;
  int                   total;

  WORKERSIDE_ARGS      *parent;

  struct worker_s      *next;
  struct worker_s      *prev;

  int                   started;
  int                   completed;
  int                   terminated;
  HMMD_SEARCH_CMD      *cmd;
} WORKER_DATA;

typedef struct {
  ESL_SQ           *sq_list;     /* list of sequences to process     */
  int               sq_cnt;      /* number of sequences              */

  pthread_mutex_t   inx_mutex;   /* protect data                     */
  int              *sq_inx;      /* next sequence to process         */

  P7_HMM           *hmm;         /* query HMM                        */
  ESL_SQ           *seq;         /* query sequence                   */
  ESL_ALPHABET     *abc;         /* digital alphabet                 */
  ESL_GETOPTS      *opts;        /* search specific options          */
  P7_BUILDER       *bld;         /* HMM construction configuration   */

  double            elapsed;

  /* Structure created and populated by the individual threads.
   * The main thread is responsible for freeing up the memory.
   */
  P7_PIPELINE      *pli;         /* work pipeline                           */
  P7_TOPHITS       *th;          /* top hit results                         */
} WORKER_INFO;

static void free_QueueData(QUEUE_DATA *data);
static void pop_Queue(SEARCH_QUEUE *queue);
static void push_Queue(QUEUE_DATA *data, SEARCH_QUEUE *queue);
static QUEUE_DATA *read_Queue(SEARCH_QUEUE *queue);

static QUEUE_DATA *read_QueryCmd(int fd);

static int  setup_masterside_comm(ESL_GETOPTS *opts);
static void setup_clientside_comm(ESL_GETOPTS *opts, CLIENTSIDE_ARGS  *args);
static void setup_workerside_comm(ESL_GETOPTS *opts, WORKERSIDE_ARGS  *args);

static void clear_results(WORKERSIDE_ARGS *comm);
static void forward_results(QUEUE_DATA *query, ESL_STOPWATCH *w, WORKERSIDE_ARGS *comm);
static void send_results(int fd, ESL_STOPWATCH *w, WORKER_INFO *info);

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

static ESL_OPTIONS cmdlineOpts[] = {
  /* name           type         default  env   range  toggles  reqs   incomp           help                                                     docgroup */
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                         1 },
  /* Other options */
  { "--master",     eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  "--worker",      "run program as the master server",                            12 },
  { "--worker",     eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  "--master",      "run program as a worker with server at <s>",                  12 },
  { "--cport",      eslARG_INT,  "41139", NULL, "n>1024",NULL,  NULL,  "--worker",      "port to use for client/server communication",                 12 },
  { "--wport",      eslARG_INT,  "41023", NULL, "n>1024",NULL,  NULL,  NULL,            "port to use for server/worker communication",                 12 },
  { "--ccncts",     eslARG_INT,     "16", NULL, "n>0",   NULL,  NULL,  "--worker",      "maximum number of client side connections to accept",         12 },
  { "--wcncts",     eslARG_INT,     "32", NULL, "n>0",   NULL,  NULL,  "--worker",      "maximum number of worker side connections to accept",         12 },

  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>0", NULL,  NULL,  "--master",      "number of parallel CPU workers to use for multithreads",      12 },

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

#define BLOCK_SIZE 1000
static void pipeline_thread(void *arg);


typedef void sig_func(int);

sig_func *
signal(int signo, sig_func *fn)
{
  struct sigaction act;
  struct sigaction oact;

  act.sa_handler = fn;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  if (signo == SIGALRM) {
#ifdef SA_INTERRUMP
    act.sa_flags |= SA_INTERRUPT;  /* SunOS 4.x */
#endif
  } else {
#ifdef SA_RESTART
    act.sa_flags |= SA_RESTART;  /* SVR4, 4.4BSD */
#endif
  }
  if (sigaction(signo, &act, &oact) < 0) {
    return SIG_ERR;
  }

  return oact.sa_handler;
}

static void
print_client_error(int fd, int status, char *format, va_list ap)
{
  int   n;
  char  ebuf[512];

  HMMD_SEARCH_STATUS s;

  s.status  = status;
  s.err_len = vsnprintf(ebuf, sizeof(ebuf), format, ap);

  /* send back an unsuccessful status message */
  n = sizeof(s);
  if (writen(fd, &s, n) != n) {
    fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
    exit(1);
  }

  writen(fd, ebuf, s.err_len);
  if (writen(fd, ebuf, s.err_len) != s.err_len) {
    fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", s.err_len, errno, strerror(errno));
    exit(1);
  }
}

static void
client_error(int fd, int status, char *format, ...)
{
  va_list ap;

  va_start(ap, format);
  print_client_error(fd, status, format, ap);
  va_end(ap);
}

static void
client_error_longjmp(int fd, int status, jmp_buf *env, char *format, ...)
{
  va_list ap;

  va_start(ap, format);
  print_client_error(fd, status, format, ap);
  va_end(ap);

  longjmp(*env, 1);
}

static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go)
{
  int n;
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(cmdlineOpts)) == NULL)    p7_Die("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK) { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) { printf("Failed to parse command line: %s\n",  go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK) { printf("Failed to parse command line: %s\n",  go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    p7_banner(stdout, argv[0], banner);
    esl_usage(stdout, argv[0], usage);

    puts("\nBasic options:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

    puts("\nOther expert options:");
    esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
    exit(0);
  }

  n = esl_opt_ArgNumber(go);
  if (n < 0) { puts("Incorrect number of command line arguments."); goto ERROR; }

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
process_searchline(int fd, char *cmdstr, ESL_GETOPTS *go)
{
  int status;

  if ((status = esl_opt_ProcessSpoof(go, cmdstr)) != eslOK) client_error(fd, status, "Failed to parse options string: %s", go->errbuf);
  if ((status = esl_opt_VerifyConfig(go))         != eslOK) client_error(fd, status, "Failed to parse options string: %s", go->errbuf);

  /* the options string can handle an optional database */
  if (esl_opt_ArgNumber(go) > 1)                 client_error(fd, eslFAIL, "Incorrect number of command line arguments.");
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

static ESL_SQCACHE **
cache_dbs(ESL_GETOPTS *go, ESL_ALPHABET *abc)
{
  int i, j;
  int status;

  ESL_SQCACHE **cache;

  /* cache all the databases into memory */
  ESL_ALLOC(cache, sizeof(ESL_SQCACHE *) * esl_opt_ArgNumber(go));
  for (i = 0; i < esl_opt_ArgNumber(go); ++i) {
    ESL_RANDOMNESS *rnd  = NULL;
    ESL_SQ         *sq   = NULL;

    char *dbfile = esl_opt_GetArg(go, i+1);

    status = esl_sqfile_Cache(abc, dbfile, eslSQFILE_FASTA, p7_SEQDBENV, &cache[i]);
    if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          dbfile);
    else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            dbfile);
    else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
    else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, dbfile);

    sq  = cache[i]->sq_list;
#if 0
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
#endif

    rnd = esl_randomness_CreateFast(cache[i]->seq_count);
    for (j = 0 ; j < cache[i]->seq_count; ++j) {
      rnd->x = rnd->x * 69069 + 1;
      sq[j].idx = rnd->x;
    }
    esl_randomness_Destroy(rnd);
    qsort(sq, cache[i]->seq_count, sizeof(ESL_SQ), sort_seq_inx);
    for (j = 0 ; j < cache[i]->seq_count; ++j) sq[j].idx = j;

#if 1
    fprintf (stdout, "%2d: %8d %5" PRId64 "M %5" PRId64 "M %5" PRId64 "M %s\n",
             i, cache[i]->seq_count, (cache[i]->seq_count * sizeof(ESL_SQ)) >> 20, cache[i]->hdr_size >> 20, cache[i]->res_size >> 20, cache[i]->filename);
#endif
  }

  return cache;

 ERROR:
  fprintf(stderr, "%s: malloc error %d - %s\n", go->argv[0], errno, strerror(errno));
  exit(1);
}


void
worker_process(ESL_GETOPTS *go)
{
  FILE            *ofp        = stdout;            /* results output file (-o)                        */
  ESL_ALPHABET    *abc;                            /* digital alphabet                                */
  ESL_STOPWATCH   *w;                              /* timer used for profiling statistics             */
  int              textw      = 0;
  int              status     = eslOK;
  int              i;
  int              fd;

  int              ncpus      = 0;
  WORKER_INFO     *info       = NULL;

  ESL_THREADS     *threadObj  = NULL;
  pthread_mutex_t  inx_mutex;
  int              current_index;

  ESL_SQCACHE    **cache      = NULL;

  QUEUE_DATA      *query      = NULL;

  /* Set processor specific flags */
  impl_Init();

  /* Initializations */
  p7_FLogsumInit();      /* we're going to use table-driven Logsum() approximations at times */

  w = esl_stopwatch_Create();
  abc = esl_alphabet_Create(eslAMINO);

  cache = cache_dbs(go, abc);

  /* start the communications with the web clients */
  fd = setup_masterside_comm(go);

  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                   esl_threads_CPUCount(&ncpus);

  threadObj = esl_threads_Create(&pipeline_thread);
  if (pthread_mutex_init(&inx_mutex, NULL) != 0) p7_Fail("mutex init failed");

  ESL_ALLOC(info, sizeof(*info) * ncpus);

  /* read query hmm/sequence */
  while ((query = read_QueryCmd(fd)) != NULL) {

    textw = (esl_opt_GetBoolean(query->opts, "--notextw")) ? 0 : esl_opt_GetInteger(query->opts, "--textw");

    esl_stopwatch_Start(w);

    if (query->hmm == NULL) {
      fprintf(ofp, "Query (%d):       %s  [L=%ld]\n", query->sock, query->seq->name, (long) query->seq->n);
    } else {
      fprintf(ofp, "Query (%d):       %s  [M=%d]\n", query->sock, query->hmm->name, query->hmm->M);
    }
    fprintf(ofp, "Database %s [%d - %d]\n", cache[query->dbx]->filename, query->sq_inx, query->sq_inx + query->sq_cnt - 1);

    /* Create processing pipeline and hit list */
    for (i = 0; i < ncpus; ++i) {
      info[i].abc   = query->abc;
      info[i].hmm   = query->hmm;
      info[i].seq   = query->seq;
      info[i].opts  = query->opts;

      info[i].th    = NULL;
      info[i].pli   = NULL;

      info[i].sq_list   = cache[query->dbx]->sq_list + query->sq_inx;
      info[i].sq_cnt    = query->sq_cnt;
      info[i].inx_mutex = inx_mutex;
      info[i].sq_inx    = &current_index;

      esl_threads_AddThread(threadObj, &info[i]);
    }

    current_index = 0;
    esl_threads_WaitForStart(threadObj);
    esl_threads_WaitForFinish(threadObj);

    esl_stopwatch_Stop(w);
#if 1
    fprintf (ofp, "   Sequences  Residues                              Elapsed\n");
    for (i = 0; i < ncpus; ++i) {
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
    /* merge the results of the search results */
    for (i = 1; i < ncpus; ++i) {
      p7_tophits_Merge(info[0].th, info[i].th);
      p7_pipeline_Merge(info[0].pli, info[i].pli);
      p7_pipeline_Destroy(info[i].pli);
      p7_tophits_Destroy(info[i].th);
    }

    p7_pli_Statistics(stdout, info[0].pli, w);
    send_results(fd, w, info);

    /* free the last of the pipeline data */
    p7_pipeline_Destroy(info->pli);
    p7_tophits_Destroy(info->th);

    free_QueueData(query);
  } /* end outer loop over query HMMs */

  pthread_mutex_destroy(&inx_mutex);
  esl_threads_Destroy(threadObj);

  for (i = 0; i < esl_opt_ArgNumber(go); ++i) {
    esl_sqfile_Free(cache[i]);
    cache[i] = NULL;
  }
  free(cache);

  free(info);

  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);

  esl_getopts_Destroy(go);
  esl_alphabet_Destroy(abc);

  return;

 ERROR:
  fprintf(stderr, "%s: malloc error %d - %s\n", go->argv[0], errno, strerror(errno));
  exit(1);
}

void
master_process(ESL_GETOPTS *go)
{
  FILE               *ofp        = stdout;    /* results output file (-o)                        */
  ESL_ALPHABET       *abc        = NULL;      /* digital alphabet                                */
  ESL_STOPWATCH      *w          = NULL;      /* timer used for profiling statistics             */
  int                 status     = eslOK;
  int                 i, n;

  int                 inx;
  int                 cnt;

  WORKER_INFO        *info       = NULL;

  int                cache_cnt   = 0;
  ESL_SQCACHE      **cache       = NULL;

  SEARCH_QUEUE       *queue      = NULL;
  QUEUE_DATA         *query      = NULL;

  CLIENTSIDE_ARGS     client_comm;
  WORKERSIDE_ARGS     worker_comm;
  WORKER_DATA        *worker;

  /* Set processor specific flags */
  impl_Init();

  /* Initializations */
  p7_FLogsumInit();     /* we're going to use table-driven Logsum() approximations at times */

  w = esl_stopwatch_Create();
  abc = esl_alphabet_Create(eslAMINO);

  cache     = cache_dbs(go, abc);
  cache_cnt = esl_opt_ArgNumber(go);

  /* initialize the search queue */
  ESL_ALLOC(queue, sizeof(SEARCH_QUEUE));
  if ((n = pthread_mutex_init(&queue->queue_mutex, NULL)) != 0) {
    errno = n;
    fprintf(stderr, "%s: mutex_init error %d - %s\n", go->argv[0], errno, strerror(errno));
    exit(1);
  }
  if ((n = pthread_cond_init(&queue->queue_cond, NULL)) != 0) {
    errno = n;
    fprintf(stderr, "%s: mutex_init error %d - %s\n", go->argv[0], errno, strerror(errno));
    exit(1);
  }

  queue->head = NULL;
  queue->tail = NULL;

  /* start the communications with the web clients */
  client_comm.pgm   = go->argv[0];
  client_comm.queue = queue;
  client_comm.cache = cache;
  client_comm.cache_size = esl_opt_ArgNumber(go);
  setup_clientside_comm(go, &client_comm);

  /* initialize the worker structure */
  if ((n = pthread_mutex_init(&worker_comm.work_mutex, NULL)) != 0) {
    errno = n;
    fprintf(stderr, "%s: mutex_init error %d - %s\n", go->argv[0], errno, strerror(errno));
    exit(1);
  }
  if ((n = pthread_cond_init(&worker_comm.start_cond, NULL)) != 0) {
    errno = n;
    fprintf(stderr, "%s: mutex_init error %d - %s\n", go->argv[0], errno, strerror(errno));
    exit(1);
  }
  if ((n = pthread_cond_init(&worker_comm.complete_cond, NULL)) != 0) {
    errno = n;
    fprintf(stderr, "%s: mutex_init error %d - %s\n", go->argv[0], errno, strerror(errno));
    exit(1);
  }
  if ((n = pthread_cond_init(&worker_comm.worker_cond, NULL)) != 0) {
    errno = n;
    fprintf(stderr, "%s: mutex_init error %d - %s\n", go->argv[0], errno, strerror(errno));
    exit(1);
  }

  worker_comm.sock_fd = -1;
  worker_comm.head    = NULL;
  worker_comm.tail    = NULL;

  worker_comm.pgm     = go->argv[0];
  setup_workerside_comm(go, &worker_comm);

  /* read query hmm/sequence */
  while ((query = read_Queue(queue)) != NULL) {

    if ((n = pthread_mutex_lock (&worker_comm.work_mutex)) != 0) {
      errno = n;
      fprintf(stderr, "%s: mutex lock error %d - %s\n", go->argv[0], errno, strerror(errno));
      exit(1);
    }

    /* wait until we have atleast one worker */
    while (worker_comm.head == NULL) {
      if ((n = pthread_cond_wait(&worker_comm.worker_cond, &worker_comm.work_mutex)) != 0) {
        errno = n;
        fprintf(stderr, "%s: mutex cond wait error %d - %s\n", go->argv[0], errno, strerror(errno));
        exit(1);
      }
    }

    /* add to all workers the new query to process */
    worker_comm.started = 0;
    worker = worker_comm.head;
    while (worker != NULL) {
      ++worker_comm.started;

      if ((worker->cmd = malloc(query->cmd->length)) == NULL) {
        fprintf(stderr, "hmmgpmd: malloc error %d - %s\n", errno, strerror(errno));
        exit(1);
      }
      memcpy(worker->cmd, query->cmd, query->cmd->length);

      worker->started    = 1;
      worker->completed  = 0;
      worker->terminated = 0;
      worker->total      = 0;

      worker            = worker->next;
    }

    /* assign each worker a portion of the database */
    inx = 0;
    i = worker_comm.started;
    cnt = cache[query->dbx]->seq_count;
    printf ("DBX: %d - %d %d\n", query->dbx, i, cnt);
    worker = worker_comm.head;
    while (worker != NULL) {
      if (worker->cmd != NULL) {
        worker->cmd->sq_inx = inx;
        worker->cmd->sq_cnt = cnt / i;
        printf ("WORKER %d: %d - %d\n", worker_comm.started - i, worker->cmd->sq_inx, worker->cmd->sq_cnt);

        inx += worker->cmd->sq_cnt;
        cnt -= worker->cmd->sq_cnt;
        --i;

        worker = worker->next;
      }
    }

    worker_comm.completed = worker_comm.started;
    worker_comm.errors    = 0;

    /* notify all the worker threads of the new query */
    if ((n = pthread_cond_broadcast(&worker_comm.start_cond)) != 0) {
      errno = n;
      fprintf(stderr, "%s: mutex cond broadcast error %d - %s\n", go->argv[0], errno, strerror(errno));
      exit(1);
    }

    if ((n = pthread_mutex_unlock (&worker_comm.work_mutex)) != 0) {
      errno = n;
      fprintf(stderr, "%s: mutex unlock error %d - %s\n", go->argv[0], errno, strerror(errno));
      exit(1);
    }

    esl_stopwatch_Start(w);

    if (query->hmm == NULL) {
      fprintf(ofp, "Query (%d):  %s  [L=%d]\n", query->sock, query->seq->name, (int)query->seq->n);
    } else {
      fprintf(ofp, "Query (%d):  %s  [M=%d]\n", query->sock, query->hmm->name, query->hmm->M);
    }

    /* Wait for all the workers to complete */
    if ((n = pthread_mutex_lock (&worker_comm.work_mutex)) != 0) {
      errno = n;
      fprintf(stderr, "%s: mutex lock error %d - %s\n", go->argv[0], errno, strerror(errno));
      exit(1);
    }

    printf("WAITING: %d\n", worker_comm.completed); fflush(stdout);
    while (worker_comm.completed > 0) {
      if ((n = pthread_cond_wait (&worker_comm.complete_cond, &worker_comm.work_mutex)) != 0) {
        errno = n;
        fprintf(stderr, "%s: cond wait error %d - %s\n", go->argv[0], errno, strerror(errno));
        exit(1);
      }
      printf("WAITING: %d\n", worker_comm.completed); fflush(stdout);
    }

    if ((n = pthread_mutex_unlock (&worker_comm.work_mutex)) != 0) {
      errno = n;
      fprintf(stderr, "%s: mutex unlock error %d - %s\n", go->argv[0], errno, strerror(errno));
      exit(1);
    }

    esl_stopwatch_Stop(w);

    if (worker_comm.errors == 0) {
      /* send the merged results to the web client */
      forward_results(query, w, &worker_comm);
      pop_Queue(queue);
    } else {
      /* TODO: after how many retries should we just kill the query???? */
      clear_results(&worker_comm);
      ++query->retry_cnt;
    }
  }

  for (i = 0; i < esl_opt_ArgNumber(go); ++i) {
    esl_sqfile_Free(cache[i]);
    cache[i] = NULL;
  }
  free(cache);

  free(info);
  free(queue);

  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);

  esl_getopts_Destroy(go);

  return;

 ERROR:
  fprintf(stderr, "%s: malloc error %d - %s\n", go->argv[0], errno, strerror(errno));
  exit(1);
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go = NULL;      /* command line processing */

  process_commandline(argc, argv, &go);    

  /* if we write to a broken socket, ignore the signal and handle the error. */
  signal(SIGPIPE, SIG_IGN);

  if (esl_opt_IsUsed(go, "--master")) master_process(go);
  if (esl_opt_IsUsed(go, "--worker")) worker_process(go);

  puts("Options --master or --worker must be specified.");

  esl_getopts_Destroy(go);

  return eslOK;
}


static QUEUE_DATA *
read_Queue(SEARCH_QUEUE *queue)
{
  int         n;
  QUEUE_DATA *data;

  if ((n = pthread_mutex_lock (&queue->queue_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex lock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }

  while (queue->head == NULL) {
    if ((n = pthread_cond_wait (&queue->queue_cond, &queue->queue_mutex)) != 0) {
      errno = n;
      fprintf(stderr, "%08X: cond wait error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }
  }

  data = queue->head;

  if ((n = pthread_mutex_unlock (&queue->queue_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }

  return data;
}

static void
pop_Queue(SEARCH_QUEUE *queue)
{
  int         n;
  QUEUE_DATA *data;

  if ((n = pthread_mutex_lock (&queue->queue_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex lock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }

  data = queue->head;
  queue->head = data->next;
  if (queue->head == NULL) queue->tail = NULL;

  if ((n = pthread_mutex_unlock (&queue->queue_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }

  free_QueueData(data);
}

static void
push_Queue(QUEUE_DATA *data, SEARCH_QUEUE *queue)
{
  int n;

  /* add the search request to the queue */
  if ((n = pthread_mutex_lock (&queue->queue_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex lock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }

  if (queue->head == NULL) {
    queue->head = data;
    queue->tail = data;
  } else {
    queue->tail->next = data;
    data->prev = queue->tail;
    queue->tail = data;
  }

  if ((n = pthread_mutex_unlock (&queue->queue_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }

  /* if anyone is waiting, wake them up */
  if ((n = pthread_cond_broadcast (&queue->queue_cond)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }
}

static void
free_QueueData(QUEUE_DATA *data)
{
  /* free the query data */
  esl_getopts_Destroy(data->opts);

  if (data->abc != NULL) esl_alphabet_Destroy(data->abc);
  if (data->hmm != NULL) p7_hmm_Destroy(data->hmm);
  if (data->seq != NULL) esl_sq_Destroy(data->seq);
  if (data->cmd != NULL) free(data->cmd);
  memset(data, 0, sizeof(*data));
  free(data);
}

static QUEUE_DATA *
read_QueryCmd(int fd)
{
  int                i;
  int                n;
  int                size;

  char              *p;
  char              *name;
  char              *desc;
  ESL_DSQ           *dsq;

  HMMD_SEARCH_CMD    hdr;
  HMMD_SEARCH_CMD   *cmd    = NULL;
  
  QUEUE_DATA        *query  = NULL;

  n = sizeof(hdr);
  if ((size = readn(fd, &hdr, n)) == -1) {
    fprintf(stderr, "hmmgpmd: read error %d - %s\n", errno, strerror(errno));
    exit(1);
  }

  if ((cmd = malloc(hdr.length)) == NULL) {
    fprintf(stderr, "hmmgpmd: malloc error %d - %s\n", errno, strerror(errno));
    exit(1);
  }

  *cmd = hdr;
  p = (char *)cmd;
  p += sizeof(hdr);

  n = hdr.length - sizeof(hdr);
  if ((size = readn(fd, p, n)) == -1) {
    fprintf(stderr, "hmmgpmd: read error %d - %s\n", errno, strerror(errno));
    exit(1);
  }

  if ((query = malloc(sizeof(QUEUE_DATA))) == NULL) {
    fprintf(stderr, "hmmgpmd: malloc error %d - %s\n", errno, strerror(errno));
    exit(1);
  }

  query->dbx       = cmd->db_inx;
  query->sq_inx    = cmd->sq_inx;
  query->sq_cnt    = cmd->sq_cnt;
  query->sock      = fd;
  query->retry     = 0;
  query->retry_cnt = 0;
  query->cmd       = NULL;
  query->next      = NULL;
  query->prev      = NULL;

  /* process search specific options */
  query->opts = esl_getopts_Create(searchOpts);
  if (query->opts == NULL)  {
    client_error(fd, eslFAIL, "hmmgpmd: internal failure creating options object");
    return 0;
  }

  if (cmd->opts_length > 1) process_searchline(fd, p, query->opts);

  query->hmm = NULL;
  query->seq = NULL;
  query->abc = esl_alphabet_Create(eslAMINO);

  /* check if we are processing a sequence or hmm */
  if (cmd->query_type == HMMD_SEQUENCE) {
    n    = cmd->query_length - 2;
    name = p + cmd->opts_length;
    desc = name + strlen(name) + 1;
    dsq  = (ESL_DSQ *) (desc + strlen(desc) + 1);
    query->seq = esl_sq_CreateDigitalFrom(query->abc, name, dsq, n, desc, NULL, NULL);
  } else {
    P7_HMM *hmm = p7_hmm_CreateShell();
    P7_HMM  thmm;

    /* allocate memory for the hmm and initialize */
    p += cmd->opts_length;
    memcpy(&thmm, p, sizeof(P7_HMM));

    hmm->flags = thmm.flags;
    p7_hmm_CreateBody(hmm, cmd->query_length, query->abc);
    p += sizeof(P7_HMM);

    /* initialize fields */
    hmm->nseq       = thmm.nseq;
    hmm->eff_nseq   = thmm.eff_nseq;
    hmm->max_length = thmm.max_length;
    hmm->checksum   = thmm.checksum;
    hmm->ctime      = NULL;
    hmm->comlog     = NULL;

    for (i = 0; i < p7_NCUTOFFS; i++) hmm->cutoff[i]  = thmm.cutoff[i];
    for (i = 0; i < p7_NEVPARAM; i++) hmm->evparam[i] = thmm.evparam[i];
    for (i = 0; i < p7_MAXABET;  i++) hmm->compo[i]   = thmm.compo[i];

    /* fill in the hmm pointers */
    n = sizeof(float) * (hmm->M + 1) * p7H_NTRANSITIONS;
    memcpy(*hmm->t, p, n);
    p += n;

    n = sizeof(float) * (hmm->M + 1) * query->abc->K;
    memcpy(*hmm->mat, p, n);
    p += n;
    memcpy(*hmm->ins, p, n);
    p += n;

    if (thmm.name != NULL) {
      hmm->name = strdup(p);
      p += strlen(hmm->name) + 1;
    }

    if (thmm.acc != NULL) {
      hmm->acc = strdup(p);
      p += strlen(hmm->name) + 1;
    }

    if (thmm.desc != NULL) {
      hmm->desc = strdup(p);
      p += strlen(hmm->name) + 1;
    }

    n = hmm->M + 2;
    if (hmm->flags & p7H_RF) {
      memcpy(hmm->rf, p, n);
      p += n;
    }

    if (hmm->flags & p7H_CS) {
      memcpy(hmm->cs, p, n);
      p += n;
    }

    if (hmm->flags & p7H_CA) {
      memcpy(hmm->ca, p, n);
      p += n;
    }

    if (hmm->flags & p7H_MAP) {
      n = sizeof(int) * (hmm->M + 1);
      memcpy(hmm->map, p, n);
      p += n;
    }

    query->hmm = hmm;
  }

  free(cmd);

  return query;
}


static void 
pipeline_thread(void *arg)
{
  int               i;
  int               count;
  int               seed;
  int               status;
  int               workeridx;
  WORKER_INFO      *info;
  ESL_THREADS      *obj;

  ESL_STOPWATCH    *w;

  P7_BUILDER       *bld      = NULL;         /* HMM construction configuration */
  P7_BG            *bg       = NULL;         /* null model                     */
  P7_PIPELINE      *pli      = NULL;         /* work pipeline                  */
  P7_TOPHITS       *th       = NULL;         /* top hit results                */
  P7_PROFILE       *gm       = NULL;         /* generic model                  */
  P7_OPROFILE      *om       = NULL;         /* optimized query profile        */

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
    status = p7_builder_SetScoreSystem(bld, 
                                       esl_opt_GetString(info->opts, "--mxfile"), 
                                       NULL, 
                                       esl_opt_GetReal(info->opts, "--popen"), 
                                       esl_opt_GetReal(info->opts, "--pextend"));
    if (status != eslOK) {
      //client_error(info->sock, status, "hmmgpmd: failed to set single query sequence score system: %s", bld->errbuf);
      fprintf(stderr, "hmmgpmd: failed to set single query sequence score system: %s", bld->errbuf);
      pthread_exit(NULL);
      return;
    }
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
    if (pthread_mutex_lock(&info->inx_mutex) != 0) p7_Fail("mutex lock failed");
    inx = *info->sq_inx;
    *info->sq_inx += BLOCK_SIZE;
    if (pthread_mutex_unlock(&info->inx_mutex) != 0) p7_Fail("mutex unlock failed");

    dbsq = info->sq_list + inx;

    count = info->sq_cnt - inx;
    if (count > BLOCK_SIZE) count = BLOCK_SIZE;
    //printf("THREAD %08x: %d %d\n", workeridx, inx, count);
#endif

#if 0
    /* grab the next block of sequences */
    if (pthread_mutex_lock(&info->inx_mutex) != 0) p7_Fail("mutex lock failed");
    inx = *info->sq_inx;
    count = info->sq_cnt - inx;
    count = count >> 7;
    if (count > 2500) count = 2500;
    if (count < 1000) count = 1000;
    *info->sq_inx += count;
    if (pthread_mutex_unlock(&info->inx_mutex) != 0) p7_Fail("mutex unlock failed");

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

  esl_stopwatch_Stop(w);
  info->elapsed = w->elapsed;

  esl_stopwatch_Destroy(w);

  esl_threads_Finished(obj, workeridx);

  pthread_exit(NULL);
  return;
}

static void
forward_results(QUEUE_DATA *query, ESL_STOPWATCH *w, WORKERSIDE_ARGS *comm)
{
  int fd;
  int i, j;
  int inx;
  int n;

  int total = 0;

  P7_HIT    **hit;
  char      **hit_data;

  WORKER_DATA        *worker;

  HMMD_SEARCH_STATS   stats;
  HMMD_SEARCH_STATUS  status;

  fd = query->sock;

  stats.nhits       = 0;
  stats.nreported   = 0;
  stats.nincluded   = 0;

  stats.nmodels     = 0;
  stats.nseqs       = 0;
  stats.nres        = 0;
  stats.nnodes      = 0;
  stats.n_past_msv  = 0;
  stats.n_past_bias = 0;
  stats.n_past_vit  = 0;
  stats.n_past_fwd  = 0;
  stats.Z           = 0;

  /* copy the search stats */
  stats.elapsed     = w->elapsed;
  stats.user        = w->user;
  stats.sys         = w->sys;

  /* lock the workers until we have merged the results */
  if ((n = pthread_mutex_lock (&comm->work_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "hmmpgmd: mutex lock error %d - %s\n", errno, strerror(errno));
    exit(1);
  }

  /* count the number of hits */
  n = 0;
  worker = comm->head;
  while (worker != NULL) {
      if (worker->completed) n += worker->stats.nhits;
      worker = worker->next;
  }

  /* allocate spaces to hold all the hits */
  hit      = malloc(sizeof(char *) * n);
  hit_data = malloc(sizeof(char *) * n);

  /* gather up the results and merge them */
  inx = 0;
  worker = comm->head;
  while (worker != NULL) {
    if (worker->completed) {
      stats.nhits       += worker->stats.nhits;
      stats.nreported   += worker->stats.nreported;
      stats.nincluded   += worker->stats.nincluded;

      stats.nmodels     += worker->stats.nmodels;
      stats.nseqs       += worker->stats.nseqs;
      stats.nres        += worker->stats.nres;
      stats.nnodes      += worker->stats.nnodes;
      stats.n_past_msv  += worker->stats.n_past_msv;
      stats.n_past_bias += worker->stats.n_past_bias;
      stats.n_past_vit  += worker->stats.n_past_vit;
      stats.n_past_fwd  += worker->stats.n_past_fwd;

      stats.domZ = worker->stats.domZ;
      if (worker->stats.Z_setby == p7_ZSETBY_NTARGETS) {
        stats.Z += (query->hmm != NULL) ? worker->stats.nmodels : worker->stats.nseqs;
      } else {
        stats.Z = worker->stats.Z;
      }

      if (worker->stats.nhits > 0) {
        n = worker->stats.nhits;
        memcpy(hit + inx, worker->hit, sizeof(char *) * n);
        memcpy(hit_data + inx, worker->hit_data, sizeof(char *) * n);

        memset(worker->hit, 0, sizeof(char *) * n);
        memset(worker->hit_data, 0, sizeof(char *) * n);

        free(worker->hit);
        free(worker->hit_data);

        worker->hit      = NULL;
        worker->hit_data = NULL;

        inx += worker->stats.nhits;
      }

      worker->completed = 0;
    }

    worker = worker->next;
  }

  if ((n = pthread_mutex_unlock (&comm->work_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "hmmpgmd: mutex unlock error %d - %s\n", errno, strerror(errno));
    exit(1);
  }

  /* send back a successful status message */
  status.status  = eslOK;
  status.err_len = 0;
  n = sizeof(status);
  total += n;
  if (writen(fd, &status, n) != n) {
    fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
    exit(1);
  }

  n = sizeof(stats);
  total += n;
  if (writen(fd, &stats, n) != n) {
    fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
    exit(1);
  }

  /* loop through the hit list sending to dest */
  for (i = 0; i < stats.nhits; ++i) {
    n = sizeof(P7_HIT);
    total += n;
    if (writen(fd, hit[i], n) != n) {
      fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
      exit(1);
    }

    if (hit_data[i] != NULL) {
      /* the hit string pointers contain the length of the strings */
      n = 0;
      if (hit[i]->name != NULL) n += (socklen_t)(hit[i]->name - (char *)NULL);
      if (hit[i]->acc  != NULL) n += (socklen_t)(hit[i]->acc  - (char *)NULL);
      if (hit[i]->desc != NULL) n += (socklen_t)(hit[i]->desc - (char *)NULL);

      total += n;
      if (writen(fd, hit_data[i], n) != n) {
        fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
        exit(1);
      }
    }

    /* send the domains */
    for (j = 0; j < hit[i]->ndom; ++j) {
      P7_DOMAIN      *dcl = &hit[i]->dcl[j];

      n = (socklen_t)sizeof(P7_DOMAIN);
      total += n;
      if (writen(fd, dcl, n) != n) {
        fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
        exit(1);
      }

      n = (socklen_t)sizeof(P7_ALIDISPLAY);
      total += n;
      if (writen(fd, dcl->ad, n) != n) {
        fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
        exit(1);
      }

      n = (socklen_t)dcl->ad->memsize;
      total += n;
      if (writen(fd, dcl->ad->mem, n) != n) {
        fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
        exit(1);
      }
    }
  }

  printf("Bytes %d sent on socket %d\n", total, fd);

  /* free all the data */
  for (i = 0; i < stats.nhits; ++i) {
    for (j = 0; j < hit[i]->ndom; ++j) {
      P7_DOMAIN *dcl = &hit[i]->dcl[j];
      free (dcl->ad);
      free (dcl->ad->mem);
    }
    free(hit[i]->dcl);
    free(hit[i]);
    
    if (hit_data[i] != NULL) free(hit_data[i]);
  }

  memset(hit, 0, sizeof(char *) * stats.nhits);
  memset(hit_data, 0, sizeof(char *) * stats.nhits);

  free(hit);
  free(hit_data);
}

static void
clear_results(WORKERSIDE_ARGS *comm)
{
  int i, j;
  int inx;
  int n;

  WORKER_DATA *worker;

  /* lock the workers until we have freed the results */
  if ((n = pthread_mutex_lock (&comm->work_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "hmmpgmd: mutex lock error %d - %s\n", errno, strerror(errno));
    exit(1);
  }

  /* free all the results */
  worker = comm->head;
  while (worker != NULL) {
    if (worker->completed) {
      if (worker->stats.nhits > 0) {
        /* free all the data */
        for (i = 0; i < worker->stats.nhits; ++i) {
          for (j = 0; j < worker->hit[i]->ndom; ++j) {
            P7_DOMAIN *dcl = &worker->hit[i]->dcl[j];
            free (dcl->ad);
            free (dcl->ad->mem);
          }
          free(worker->hit[i]->dcl);
          free(worker->hit[i]);
          
          if (worker->hit_data[i] != NULL) free(worker->hit_data[i]);
        }

        memset(worker->hit, 0, sizeof(char *) * n);
        memset(worker->hit_data, 0, sizeof(char *) * n);

        free(worker->hit);
        free(worker->hit_data);

        worker->hit      = NULL;
        worker->hit_data = NULL;

        if (worker->err_buf != NULL) {
          free(worker->err_buf);
          worker->err_buf = NULL;
        }

        inx += worker->stats.nhits;
      }

      worker->completed = 0;
    }

    /* check if the worker node has terminated */
    if (worker->terminated) {
      WORKER_DATA *temp = worker;

      worker = worker->next;

      if (temp->prev == NULL) {
        comm->head = temp->next;
      } else {
        temp->prev->next = temp->next;
      }
      if (temp->next == NULL) {
        comm->tail = temp->prev;
      } else {
        temp->next->prev = temp->prev;
      }

      temp->prev   = NULL;
      temp->next   = NULL;
      temp->parent = NULL;

      free(temp);
    } else {
      worker = worker->next;
    }
  }

  if ((n = pthread_mutex_unlock (&comm->work_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "hmmpgmd: mutex unlock error %d - %s\n", errno, strerror(errno));
    exit(1);
  }
}

static void
send_results(int fd, ESL_STOPWATCH *w, WORKER_INFO *info)
{
  int i, j;
  int n;

  int total = 0;

  P7_HIT    *hit;
  P7_DOMAIN *dcl;

  HMMD_SEARCH_STATS   stats;
  HMMD_SEARCH_STATUS  status;

  /* send back a successful status message */
  status.status  = eslOK;
  status.err_len = 0;
  n = sizeof(status);
  total += n;
  if (writen(fd, &status, n) != n) {
    fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
    exit(1);
  }

  /* copy the search stats */
  stats.elapsed     = w->elapsed;
  stats.user        = w->user;
  stats.sys         = w->sys;

  stats.nmodels     = info->pli->nmodels;
  stats.nseqs       = info->pli->nseqs;
  stats.nres        = info->pli->nres;
  stats.nnodes      = info->pli->nnodes;
  stats.n_past_msv  = info->pli->n_past_msv;
  stats.n_past_bias = info->pli->n_past_bias;
  stats.n_past_vit  = info->pli->n_past_vit;
  stats.n_past_fwd  = info->pli->n_past_fwd;

  stats.Z           = info->pli->Z;
  stats.domZ        = info->pli->domZ;
  stats.Z_setby     = info->pli->Z_setby;
  stats.domZ_setby  = info->pli->domZ_setby;

  stats.nhits       = info->th->N;
  stats.nreported   = info->th->nreported;
  stats.nincluded   = info->th->nincluded;

  n = sizeof(stats);
  total += n;
  if (writen(fd, &stats, n) != n) {
    fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
    exit(1);
  }
  printf("SENT STATS %d\n", (int)info->th->N); fflush(stdout);

  /* loop through the hit list sending to dest */
  printf("SENDING %d hits ", (int)info->th->N); fflush(stdout);
  hit = info->th->unsrt;
  for (i = 0; i < info->th->N; ++i) {
    int     l;
    char   *h;
    char   *p;

    n = sizeof(P7_HIT);
    if (hit->name != NULL) n += strlen(hit->name) + 1;
    if (hit->acc  != NULL) n += strlen(hit->acc)  + 1;
    if (hit->desc != NULL) n += strlen(hit->desc) + 1;

    if ((h = malloc(n)) == NULL) {
      fprintf(stderr, "hmmpgmd: malloc error (size %d)\n", n);
      exit(1);
    }

    memcpy(h, hit, sizeof(P7_HIT));

    /* copy the strings to the memory space after the P7_HIT structure.
     * then replace the pointers with the lenght of the strings.  this
     * will allow the client to get the hit data with two receives.  the
     * first receive will be the hit structure.  then the space for the
     * strings can be calculated and one more receive for the string data.
     */
    p = h + sizeof(P7_HIT);
    if (hit->name != NULL) {
      strcpy(p, hit->name);
      l = strlen(hit->name) + 1;
      ((P7_HIT *)h)->name = ((char *) NULL) + l;
      p += l;
    }
    if (hit->acc  != NULL) {
      strcpy(p, hit->acc);
      l = strlen(hit->acc) + 1;
      ((P7_HIT *)h)->acc = ((char *) NULL) + l;
      p += l;
    }
    if (hit->desc != NULL) {
      strcpy(p, hit->desc);
      l = strlen(hit->desc) + 1;
      ((P7_HIT *)h)->desc = ((char *) NULL) + l;
      p += l;
    }

    total += n;
    if (writen(fd, h, n) != n) {
      fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
      exit(1);
    }

    /* send the domains for this hit */
    dcl = hit->dcl;
    for (j = 0; j < hit->ndom; ++j) {
      n = sizeof(P7_DOMAIN);
      total += n;
      if (writen(fd, dcl, n) != n) {
        fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
        exit(1);
      }

      n = sizeof(P7_ALIDISPLAY);
      total += n;
      if (writen(fd, dcl->ad, n) != n) {
        fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
        exit(1);
      }

      n = dcl->ad->memsize;
      total += n;
      if (writen(fd, dcl->ad->mem, n) != n) {
        fprintf(stderr, "hmmpgmd: write (size %d) error %d - %s\n", n, errno, strerror(errno));
        exit(1);
      }

      ++dcl;
    }

    free(h);

    ++hit;
    if (((i + 1) % 10) == 0) { printf("."); fflush(stdout); }
  }
  printf("\n"); fflush(stdout);

  printf("Bytes %d sent on socket %d\n", total, fd);
  fflush(stdout);
}


static int 
setup_masterside_comm(ESL_GETOPTS *opts)
{
  int                  fd;
  struct sockaddr_in   addr;

  /* Create a reliable, stream socket using TCP */
  if ((fd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
    fprintf(stderr, "%s: socket error %d - %s\n", opts->argv[0], errno, strerror(errno));
    exit(1);
  }

  /* Construct the server address structure */
  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_port   = htons(esl_opt_GetInteger(opts, "--wport"));
  if ((inet_pton(AF_INET, esl_opt_GetString(opts, "--worker"), &addr.sin_addr)) < 0) {
    fprintf(stderr, "%s: inet_pton error %d - %s\n", opts->argv[0], errno, strerror(errno));
    exit(1);
  }

  /* Establish the connection to the echo server */
  if (connect(fd, (struct sockaddr *) &addr, sizeof(addr)) < 0) {
    fprintf(stderr, "%s: connect error %d - %s\n", opts->argv[0], errno, strerror(errno));
    exit(1);
  }

  return fd;
}


static int
clientside_loop(CLIENTSIDE_ARGS *data)
{
  int                i;
  int                status;

  char              *ptr;
  char              *buffer;
  char              *opt_str;

  int                dbx;
  int                buf_size;
  int                remaining;
  int                amount;
  int                eod;
  int                n;

  P7_HMM            *hmm     = NULL;     /* query HMM                      */
  ESL_SQ            *seq     = NULL;     /* query sequence                 */
  ESL_ALPHABET      *abc     = NULL;     /* digital alphabet               */
  ESL_GETOPTS       *opts    = NULL;     /* search specific options        */
  HMMD_SEARCH_CMD   *cmd     = NULL;     /* search cmd to send to workers  */

  SEARCH_QUEUE      *queue   = data->queue;
  ESL_SQCACHE      **cache   = data->cache;
  QUEUE_DATA        *parms;

  jmp_buf            jmp_env;

  buf_size = MAX_BUFFER;
  buffer = malloc(buf_size);
  if (buffer == NULL) {
    fprintf(stderr, "%s: malloc error\n", data->pgm);
    exit(1);
  }

  opt_str = malloc(MAX_BUFFER);
  if (opt_str == NULL) {
    fprintf(stderr, "%s: malloc error\n", data->pgm);
    exit(1);
  }

  ptr = buffer;
  remaining = buf_size;
  amount = 0;

  eod = 0;
  while (!eod) {

    /* Receive message from client */
    if ((n = read(data->sock_fd, ptr, remaining)) < 0) {
      fprintf(stderr, "%08X: recv error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    if (n == 0) return 1;

    ptr += n;
    amount += n;
    remaining -= n;
    eod = (amount > 1 && *(ptr - 2) == '/' && *(ptr - 2) == '/' );

    /* if the buffer is full, make it larger */
    if (!eod && remaining == 0) {
      buffer = realloc(buffer, buf_size * 2);
      if (buffer == NULL) {
        fprintf(stderr, "%s: realloc error\n", data->pgm);
        exit(1);
      }

      ptr = buffer + buf_size;
      remaining = buf_size;
      buf_size *= 2;
    }
  }

  /* zero terminate the buffer */
  if (remaining == 0) {
    buffer = realloc(buffer, buf_size + 1);
    if (buffer == NULL) {
      fprintf(stderr, "%s: realloc error\n", data->pgm);
      exit(1);
    }

    ptr = buffer + buf_size;
  }
  *ptr = 0;

  /* skip all leading white spaces */
  ptr = buffer;
  while (*ptr && isspace(*ptr)) ++ptr;

  /* process search specific options */
  opts = esl_getopts_Create(searchOpts);
  if (opts == NULL)  {
    client_error(data->sock_fd, eslFAIL, "Internal failure creating options object");
    free(buffer);
    return 0;
  }

  opt_str[0] = 0;
  if (*ptr == '@') {
    char *s = ++ptr;

    /* skip to the end of the line */
    while (*ptr && (*ptr != '\n' && *ptr != '\r')) ++ptr;
    *ptr++ = 0;

    /* create a commandline string with dummy program name for
     * the esl_opt_ProcessSpoof() function to parse.
     */
    strncpy(opt_str, "hmmpgmd ", MAX_BUFFER);
    strncat(opt_str, s, MAX_BUFFER);
    strncat(opt_str, "\n", MAX_BUFFER);
    opt_str[MAX_BUFFER-1] = 0;

    /* skip remaining white spaces */
    while (*ptr && isspace(*ptr)) ++ptr;
  }

  if (strncmp(ptr, "//", 2) == 0) {
    client_error(data->sock_fd, eslEFORMAT, "Missing search sequence/hmm");
    free(buffer);
    return 0;
  }

  if (!setjmp(jmp_env)) {
    dbx = 0;
    
    if (opt_str[0] != 0) process_searchline(data->sock_fd, opt_str, opts);

    /* figure out which cached database to use */
    if (esl_opt_ArgNumber(opts) == 1) {
      char *db  = esl_opt_GetArg(opts, 1);
      int   len = strlen(db);
      for (i = 0; i < data->cache_size; ++i) {
        int n = strlen(cache[i]->filename);
        if (n >= len) {
          n = n - len;
          if (strcmp(cache[i]->filename + n, db) == 0) {
            dbx = i;
            break;
          }
        }
      }
      if (i >= data->cache_size) {
        client_error_longjmp(data->sock_fd, eslENOTFOUND, &jmp_env, "Database %s was not cached", esl_opt_GetArg(opts, 1));
      }
    }

    abc = esl_alphabet_Create(eslAMINO);

    seq = NULL;
    hmm = NULL;
    if (*ptr == '>') {
      /* try to parse the input buffer as a FASTA sequence */
      seq = esl_sq_CreateDigital(abc);
      status = esl_sqio_Parse(ptr, strlen(ptr), seq, eslSQFILE_DAEMON);
      if (status != eslOK) client_error_longjmp(data->sock_fd, status, &jmp_env, "Error parsing FASTA sequence");

    } else if (strncmp(ptr, "HMM", 3) == 0) {
      P7_HMMFILE   *hfp     = NULL;

      /* try to parse the buffer as an hmm */
      status = p7_hmmfile_OpenBuffer(ptr, strlen(ptr), &hfp);
      if (status != eslOK) client_error_longjmp(data->sock_fd, status, &jmp_env, "Error opening query hmm: %s", hfp->errbuf);

      status = p7_hmmfile_Read(hfp, &abc,  &hmm);
      if (status != eslOK) client_error_longjmp(data->sock_fd, status, &jmp_env, "Error reading query hmm: %s", hfp->errbuf);

      p7_hmmfile_Close(hfp);

    } else {
      /* no idea what we are trying to parse */
      client_error_longjmp(data->sock_fd, eslEFORMAT, &jmp_env, "Unknown query sequence/hmm format");
    }
  } else {
    /* an error occured some where, so try to clean up */
    if (opts != NULL) esl_getopts_Destroy(opts);
    if (abc  != NULL) esl_alphabet_Destroy(abc);
    if (hmm  != NULL) p7_hmm_Destroy(hmm);
    if (seq  != NULL) esl_sq_Destroy(seq);

    free(buffer);

    return 0;
  }

  parms = malloc(sizeof(QUEUE_DATA));
  if (parms == NULL) {
    fprintf(stderr, "%s: malloc error\n", data->pgm);
    exit(1);
  }

  /* build the search structure that will be sent to all the workers */
  n = sizeof(HMMD_SEARCH_CMD);
  n = n + strlen(opt_str) + 1;

  if (seq != NULL) {
    n = n + strlen(seq->name) + 1;
    n = n + strlen(seq->desc) + 1;
    n = n + seq->n + 2;
  } else {
    n = n + sizeof(P7_HMM);
    n = n + sizeof(float) * (hmm->M + 1) * p7H_NTRANSITIONS;
    n = n + sizeof(float) * (hmm->M + 1) * abc->K;
    n = n + sizeof(float) * (hmm->M + 1) * abc->K;
    if (hmm->name   != NULL) n = n + strlen(hmm->name) + 1;
    if (hmm->acc    != NULL) n = n + strlen(hmm->acc)  + 1;
    if (hmm->desc   != NULL) n = n + strlen(hmm->desc) + 1;
    if (hmm->flags & p7H_RF) n = n + hmm->M + 2;
    if (hmm->flags & p7H_CS) n = n + hmm->M + 2;
    if (hmm->flags & p7H_CA) n = n + hmm->M + 2;
    if (hmm->flags & p7H_MAP) n = n + sizeof(int) * (hmm->M + 1);
  }

  cmd = malloc(n);
  if (cmd == NULL) {
    fprintf(stderr, "%s: malloc error\n", data->pgm);
    exit(1);
  }

  cmd->length       = n;
  cmd->command      = HMMD_SEARCH;
  cmd->db_inx       = dbx;
  cmd->opts_length  = strlen(opt_str) + 1;

  ptr = (char *) cmd;
  ptr += sizeof(HMMD_SEARCH_CMD);

  memcpy(ptr, opt_str, cmd->opts_length);
  ptr += cmd->opts_length;
  
  if (seq != NULL) {
    cmd->query_type   = HMMD_SEQUENCE;
    cmd->query_length = seq->n + 2;

    n = strlen(seq->name) + 1;
    memcpy(ptr, seq->name, n);
    ptr += n;

    n = strlen(seq->desc) + 1;
    memcpy(ptr, seq->desc, n);
    ptr += n;

    n = seq->n + 2;
    memcpy(ptr, seq->dsq, n);
    ptr += n;
  } else {
    cmd->query_type   = HMMD_HMM;
    cmd->query_length = hmm->M;

    n = sizeof(P7_HMM);
    memcpy(ptr, hmm, n);
    ptr += n;

    n = sizeof(float) * (hmm->M + 1) * p7H_NTRANSITIONS;
    memcpy(ptr, *hmm->t, n);
    ptr += n;

    n = sizeof(float) * (hmm->M + 1) * abc->K;
    memcpy(ptr, *hmm->mat, n);
    ptr += n;
    memcpy(ptr, *hmm->ins, n);
    ptr += n;

    if (hmm->name != NULL) {
      n = strlen(hmm->name) + 1;
      memcpy(ptr, hmm->name, n);
      ptr += n;
    }

    if (hmm->acc != NULL) {
      n = strlen(hmm->acc)  + 1;
      memcpy(ptr, hmm->name, n);
      ptr += n;
    }

    if (hmm->desc != NULL) {
      n = strlen(hmm->desc) + 1;
      memcpy(ptr, hmm->desc, n);
      ptr += n;
    }

    n = hmm->M + 2;
    if (hmm->flags & p7H_RF) {
      memcpy(ptr, hmm->rf, n);
      ptr += n;
    }

    if (hmm->flags & p7H_CS) {
      memcpy(ptr, hmm->cs, n);
      ptr += n;
    }

    if (hmm->flags & p7H_CA) {
      memcpy(ptr, hmm->ca, n);
      ptr += n;
    }

    if (hmm->flags & p7H_MAP) {
      n = sizeof(int) * (hmm->M + 1);
      memcpy(ptr, hmm->map, n);
      ptr += n;
    }
  }

  parms->hmm  = hmm;
  parms->seq  = seq;
  parms->abc  = abc;
  parms->opts = opts;
  parms->dbx  = dbx;
  parms->cmd  = cmd;

  parms->sock = data->sock_fd;
  parms->next = NULL;
  parms->prev = NULL;

  parms->retry     = 0;
  parms->retry_cnt = 0;

  if (parms->seq != NULL) {
    printf("Waiting to queue sequence %s on %d\n", parms->seq->name, parms->sock);
  } else {
    printf("Waiting to queue hmm %s on %d\n", parms->hmm->name, parms->sock);
  }

  push_Queue(parms, queue);

  if (parms->seq != NULL) {
    printf("Queued sequence %s on %d\n", parms->seq->name, parms->sock);
  } else {
    printf("Queued hmm %s on %d\n", parms->hmm->name, parms->sock);
  }

  free(buffer);
  return 0;
}

static void *
clientside_thread(void *arg)
{
  int              eof;
  CLIENTSIDE_ARGS *data = (CLIENTSIDE_ARGS *)arg;

  /* Guarantees that thread resources are deallocated upon return */
  pthread_detach(pthread_self()); 

  eof = 0;
  while (!eof) {
    eof = clientside_loop(data);
  }

  printf("Closing socket %d\n", data->sock_fd);

  close(data->sock_fd);
  free(data);

  pthread_exit(NULL);
}

static void *
client_comm_thread(void *arg)
{
  int                  n;
  int                  fd;
  pthread_t            thread_id;

  struct sockaddr_in   addr;

  CLIENTSIDE_ARGS     *targs = NULL;
  CLIENTSIDE_ARGS     *data  = (CLIENTSIDE_ARGS *)arg;
  SEARCH_QUEUE        *queue = data->queue;

  for ( ;; ) {

    /* Wait for a client to connect */
    n = sizeof(addr);
    if ((fd = accept(data->sock_fd, (struct sockaddr *)&addr, (unsigned int *)&n)) < 0) {
      fprintf(stderr, "%s: accept error %d - %s\n", data->pgm, errno, strerror(errno));
      exit(1);
    }

    printf("Handling client %s (%d)\n", inet_ntoa(addr.sin_addr), fd);

    targs = malloc(sizeof(CLIENTSIDE_ARGS));
    if (targs == NULL) {
      fprintf(stderr, "%s: malloc error %d - %s\n", data->pgm, errno, strerror(errno));
      exit(1);
    }
    targs->pgm        = data->pgm;
    targs->queue      = queue;
    targs->sock_fd    = fd;
    targs->cache      = data->cache;
    targs->cache_size = data->cache_size;

    if (pthread_create(&thread_id, NULL, clientside_thread, (void *)targs) != 0) {
      fprintf(stderr, "%s: pthread_create error %d - %s\n", data->pgm, errno, strerror(errno));
      exit(1);
    }
  }
  
  pthread_exit(NULL);
}

static void 
setup_clientside_comm(ESL_GETOPTS *opts, CLIENTSIDE_ARGS *args)
{
  int                  n;
  int                  reuse;
  int                  sock_fd;
  pthread_t            thread_id;

  struct linger        linger;
  struct sockaddr_in   addr;

  /* Create socket for incoming connections */
  if ((sock_fd = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
    fprintf(stderr, "%s: socket error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }
      
  /* incase the server went down in an ungraceful way, allow the port to be
   * reused avoiding the timeout.
   */
  reuse = 1;
  if (setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, (char *)&reuse, sizeof(reuse)) < 0) {
    fprintf(stderr, "%s: setsockopt error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }

  /* the sockets are never closed, so if the server exits, force the kernel to
   * close the socket and clear it so the server can be restarted immediately.
   */
  linger.l_onoff = 1;
  linger.l_linger = 0;
  if (setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, (char *)&linger, sizeof(linger)) < 0) {
    fprintf(stderr, "%s: setsockopt error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }

  /* Construct local address structure */
  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_ANY);
  addr.sin_port = htons(esl_opt_GetInteger(opts, "--cport"));

  /* Bind to the local address */
  if (bind(sock_fd, (struct sockaddr *) &addr, sizeof(addr)) < 0) {
    fprintf(stderr, "%s: bind error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }

  /* Mark the socket so it will listen for incoming connections */
  if (listen(sock_fd, esl_opt_GetInteger(opts, "--ccncts")) < 0) {
    fprintf(stderr, "%s: listen error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }

  args->sock_fd = sock_fd;

  if ((n = pthread_create(&thread_id, NULL, client_comm_thread, (void *)args)) != 0) {
    errno = n;
    fprintf(stderr, "%s: pthread_create error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }
}

static int
workerside_loop(WORKERSIDE_ARGS *data, WORKER_DATA *worker)
{
  int          n;
  int          i, j;
  int          size;
  int          total;

  P7_HIT      *hit     = NULL;
  P7_DOMAIN   *dcl     = NULL;

  HMMD_SEARCH_STATS  *stats;

  for ( ; ; ) {

    /* wait for the next search object */
    if ((n = pthread_mutex_lock (&data->work_mutex)) != 0) {
      errno = n;
      fprintf(stderr, "%08X: mutex lock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    /* wait for the master's signal to start the calculations */
    while (worker->cmd == NULL) {
      if ((n = pthread_cond_wait(&data->start_cond, &data->work_mutex)) != 0) {
        errno = n;
        fprintf(stderr, "%08X: mutex cond wait error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
        exit(1);
      }
    }

    if ((n = pthread_mutex_unlock (&data->work_mutex)) != 0) {
      errno = n;
      fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    n = worker->cmd->length;
    printf ("Sending data %d:\n", n);
    if (writen(worker->sock_fd, worker->cmd, n) != n) {
      fprintf(stderr, "%08X: write (size %d) error %d - %s\n", (unsigned int)pthread_self(), n, errno, strerror(errno));
      return 0;
    }

    total = 0;
    worker->total = 0;

    n = sizeof(worker->status);
    total += n;
    if ((size = readn(worker->sock_fd, &worker->status, n)) == -1) {
      fprintf(stderr, "%08X: read error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      return 0;
    }

    if (worker->status.status != eslOK) {
      n = worker->status.err_len;
      total += n; 
      worker->err_buf = malloc(n);
      if ((size = readn(worker->sock_fd, worker->err_buf, n)) == -1) {
        fprintf(stderr, "%08X: read error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
        return 0;
      }
      break;
    }

    n = sizeof(worker->stats);
    total += n;
    if ((size = readn(worker->sock_fd, &worker->stats, n)) == -1) {
      fprintf(stderr, "%08X: read error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      return 0;
    }

    stats = &worker->stats;
    worker->hit      = malloc(sizeof(char *) * stats->nhits);
    worker->hit_data = malloc(sizeof(char *) * stats->nhits);

    /* loop through the hit list sending to dest */
    for (i = 0; i < stats->nhits; ++i) {
      n = sizeof(P7_HIT);
      if ((hit = malloc(n)) == NULL) {
        fprintf(stderr, "hmmpgmd: malloc error (size %d)\n", (int)n);
        return 0;
      }
      worker->hit[i] = hit;
      total += n;
      if ((size = readn(worker->sock_fd, hit, n)) == -1) {
        fprintf(stderr, "%08X: read error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
        return 0;
      }

      /* the hit string pointers contain the length of the string including
       * the null terminator at the end.
       */
      n = 0;
      if (hit->name != NULL) n += (socklen_t)(hit->name - (char *)NULL);
      if (hit->acc  != NULL) n += (socklen_t)(hit->acc  - (char *)NULL);
      if (hit->desc != NULL) n += (socklen_t)(hit->desc - (char *)NULL);

      if (n == 0) {
        worker->hit_data[i] = NULL;
      } else {
        total += n; 
        if ((worker->hit_data[i] = malloc(n)) == NULL) {
          fprintf(stderr, "hmmpgmd: malloc error (size %d)\n", (int)n);
          return 0;
        }
        if ((size = readn(worker->sock_fd, worker->hit_data[i], n)) == -1) {
          fprintf(stderr, "%08X: read error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
          return 0;
        }
      }

      n = sizeof(P7_DOMAIN) * hit->ndom;
      if ((hit->dcl = malloc(n)) == NULL) {
        fprintf(stderr, "hmmpgmd: malloc error (size %d)\n", (int)n);
        return 0;
      }

      /* read the domains for this hit */
      dcl = hit->dcl;
      for (j = 0; j < hit->ndom; ++j) {
        char *base;
        P7_ALIDISPLAY *ad = NULL;

        n = (socklen_t)sizeof(P7_DOMAIN);
        total += n;
        if ((size = readn(worker->sock_fd, dcl, n)) == -1) {
          fprintf(stderr, "%08X: read error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
          return 0;
        }

        n = (socklen_t)sizeof(P7_ALIDISPLAY);
        total += n;
        if ((dcl->ad = malloc(n)) == NULL) {
          fprintf(stderr, "hmmpgmd: malloc error (size %d)\n", (int)n);
          return 0;
        }
        if ((size = readn(worker->sock_fd, dcl->ad, n)) == -1) {
          fprintf(stderr, "%08X: read error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
          return 0;
        }

        /* save off the original mem pointer so all the pointers can be adjusted
         * to the new block of memory.
         */
        base = dcl->ad->mem;

        n = (socklen_t)dcl->ad->memsize;
        total += n;
        if ((dcl->ad->mem = malloc(n)) == NULL) {
          fprintf(stderr, "hmmpgmd: malloc error (size %d)\n", n);
          return 0;
        }
        if ((size = readn(worker->sock_fd, dcl->ad->mem, n)) == -1) {
          fprintf(stderr, "%08X: read error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
          return 0;
        }

        /* readjust all the pointers to the new memory block */
        ad = dcl->ad;
        if (ad->rfline  != NULL) ad->rfline  = ad->rfline  - base + ad->mem;
        if (ad->csline  != NULL) ad->csline  = ad->csline  - base + ad->mem;
        if (ad->model   != NULL) ad->model   = ad->model   - base + ad->mem;
        if (ad->mline   != NULL) ad->mline   = ad->mline   - base + ad->mem;
        if (ad->aseq    != NULL) ad->aseq    = ad->aseq    - base + ad->mem;
        if (ad->ppline  != NULL) ad->ppline  = ad->ppline  - base + ad->mem;
        if (ad->hmmname != NULL) ad->hmmname = ad->hmmname - base + ad->mem;
        if (ad->hmmacc  != NULL) ad->hmmacc  = ad->hmmacc  - base + ad->mem;
        if (ad->hmmdesc != NULL) ad->hmmdesc = ad->hmmdesc - base + ad->mem;
        if (ad->sqname  != NULL) ad->sqname  = ad->sqname  - base + ad->mem;
        if (ad->sqacc   != NULL) ad->sqacc   = ad->sqacc   - base + ad->mem;
        if (ad->sqdesc  != NULL) ad->sqdesc  = ad->sqdesc  - base + ad->mem;

        ++dcl;
      }
    }

    if ((n = pthread_mutex_lock (&data->work_mutex)) != 0) {
      errno = n;
      fprintf(stderr, "%08X: mutex lock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    /* set the state of the worker to completed */
    worker->cmd       = NULL;
    worker->completed = 1;
    worker->total     = total;
    --data->completed;

    /* notify all the worker threads of the new query */
    if ((n = pthread_cond_broadcast(&data->complete_cond)) != 0) {
      errno = n;
      fprintf(stderr, "%08X: mutex cond broadcast error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    if ((n = pthread_mutex_unlock (&data->work_mutex)) != 0) {
      errno = n;
      fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    printf ("WORKER COMPLETE: received %d bytes\n", total), fflush(stdout);
  }

  return 1;
}

static void *
workerside_thread(void *arg)
{
  int              n;

  int              fd;

  WORKER_DATA      *worker  = (WORKER_DATA *)arg;
  WORKERSIDE_ARGS  *parent  = (WORKERSIDE_ARGS *)worker->parent;

  /* Guarantees that thread resources are deallocated upon return */
  pthread_detach(pthread_self()); 

  worker->next = NULL;
  worker->prev = NULL;

  /* add the worker to the new fd list */
  if ((n = pthread_mutex_lock (&parent->work_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex lock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }

  if (parent->head == NULL) {
    parent->head = worker;
    parent->tail = worker;
  } else {
    parent->tail->next = worker;
    worker->prev = parent->tail;
    parent->tail = worker;
  }

  /* notify the master process that there is atleast one worker now */
  if (parent->head == parent->tail) {
    if ((n = pthread_cond_broadcast (&parent->worker_cond)) != 0) {
      errno = n;
      fprintf(stderr, "%08X: mutex cond broadcast error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }
  }

  if ((n = pthread_mutex_unlock (&parent->work_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }

  workerside_loop(parent, worker);

  if ((n = pthread_mutex_lock (&parent->work_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex lock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }

  fd = worker->sock_fd;

  worker->sock_fd    = -1;
  worker->terminated = 1;
  --parent->completed;
  ++parent->errors;

  if ((n = pthread_mutex_unlock (&parent->work_mutex)) != 0) {
    errno = n;
    fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }

  printf("Closing socket %d\n", fd);

  close(fd);

  pthread_exit(NULL);
}

static void *
worker_comm_thread(void *arg)
{
  int                  n;
  int                  fd;
  pthread_t            thread_id;

  struct sockaddr_in   addr;

  WORKERSIDE_ARGS     *data  = (WORKERSIDE_ARGS *)arg;
  WORKER_DATA         *worker;

  for ( ;; ) {

    /* Wait for a worker to connect */
    n = sizeof(addr);
    if ((fd = accept(data->sock_fd, (struct sockaddr *)&addr, (unsigned int *)&n)) < 0) {
      fprintf(stderr, "%s: accept error %d - %s\n", data->pgm, errno, strerror(errno));
      exit(1);
    }

    printf("Handling worker %s (%d)\n", inet_ntoa(addr.sin_addr), fd);

    worker = malloc(sizeof(WORKER_DATA));
    if (worker == NULL) {
      fprintf(stderr, "%s: malloc error %d - %s\n", data->pgm, errno, strerror(errno));
      exit(1);
    }
    worker->parent     = data;
    worker->sock_fd    = fd;
    worker->next       = NULL;
    worker->prev       = NULL;

    worker->err_buf    = NULL;
    worker->hit        = NULL;
    worker->hit_data   = NULL;
    worker->total      = 0;
    worker->completed  = 0;
    worker->cmd        = NULL;

    if (pthread_create(&thread_id, NULL, workerside_thread, (void *)worker) != 0) {
      fprintf(stderr, "%s: pthread_create error %d - %s\n", data->pgm, errno, strerror(errno));
      exit(1);
    }
  }
  
  pthread_exit(NULL);
}

static void 
setup_workerside_comm(ESL_GETOPTS *opts, WORKERSIDE_ARGS *args)
{
  int                  n;
  int                  reuse;
  int                  sock_fd;
  pthread_t            thread_id;

  struct linger        linger;
  struct sockaddr_in   addr;

  /* Create socket for incoming connections */
  if ((sock_fd = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
    fprintf(stderr, "%s: socket error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }
      
  /* incase the server went down in an ungraceful way, allow the port to be
   * reused avoiding the timeout.
   */
  reuse = 1;
  if (setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, (char *)&reuse, sizeof(reuse)) < 0) {
    fprintf(stderr, "%s: setsockopt error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }

  /* the sockets are never closed, so if the server exits, force the kernel to
   * close the socket and clear it so the server can be restarted immediately.
   */
  linger.l_onoff = 1;
  linger.l_linger = 0;
  if (setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, (char *)&linger, sizeof(linger)) < 0) {
    fprintf(stderr, "%s: setsockopt error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }

  /* Construct local address structure */
  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_ANY);
  addr.sin_port = htons(esl_opt_GetInteger(opts, "--wport"));

  /* Bind to the local address */
  if (bind(sock_fd, (struct sockaddr *) &addr, sizeof(addr)) < 0) {
    fprintf(stderr, "%s: bind error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }

  /* Mark the socket so it will listen for incoming connections */
  if (listen(sock_fd, esl_opt_GetInteger(opts, "--wcncts")) < 0) {
    fprintf(stderr, "%s: listen error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }

  args->sock_fd = sock_fd;

  if ((n = pthread_create(&thread_id, NULL, worker_comm_thread, (void *)args)) != 0) {
    errno = n;
    fprintf(stderr, "%s: pthread_create error %d - %s\n", args->pgm, errno, strerror(errno));
    exit(1);
  }
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

