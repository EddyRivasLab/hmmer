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
#include <syslog.h>

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
#include "cachedb.h"

#define MAX_BUFFER 4096

typedef struct queue_data_s {
  uint32_t            srch_type;   /* type of search to preform      */
  uint32_t            query_type;  /* type of the query              */
  P7_HMM             *hmm;         /* query HMM                      */
  ESL_SQ             *seq;         /* query sequence                 */
  ESL_ALPHABET       *abc;         /* digital alphabet               */
  ESL_GETOPTS        *opts;        /* search specific options        */
  HMMD_SEARCH_CMD    *cmd;         /* workers search command         */

  int                 sock;        /* socket descriptor of client    */
  char                ip_addr[INET_ADDRSTRLEN];

  int                 retry;
  int                 retry_cnt;

  int                 dbx;         /* database index to search       */
  int                 inx;         /* sequence index to start search */
  int                 cnt;         /* number of sequences to search  */

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
  int             sock_fd;
  char            ip_addr[INET_ADDRSTRLEN];

  SEARCH_QUEUE   *queue;
  SEQ_CACHE      *seq_db;
  HMM_CACHE      *hmm_db;
} CLIENTSIDE_ARGS;

typedef struct {
  int              sock_fd;

  pthread_mutex_t  work_mutex;
  pthread_cond_t   start_cond;
  pthread_cond_t   worker_cond;
  pthread_cond_t   complete_cond;

  SEQ_CACHE       *seq_db;
  HMM_CACHE       *hmm_db;

  int              started;
  struct worker_s *head;
  struct worker_s *tail;

  int              completed;
  int              errors;
} WORKERSIDE_ARGS;

typedef struct worker_s {
  int                   sock_fd;
  char                  ip_addr[INET_ADDRSTRLEN];
  
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
  HMMER_SEQ       **sq_list;     /* list of sequences to process     */
  int               sq_cnt;      /* number of sequences              */

  P7_OPROFILE      *om_list;     /* list of profiles to process      */
  int               om_cnt;      /* number of profiles               */

  pthread_mutex_t   inx_mutex;   /* protect data                     */
  int              *inx;         /* next index to process            */

  P7_HMM           *hmm;         /* query HMM                        */
  ESL_SQ           *seq;         /* query sequence                   */
  ESL_ALPHABET     *abc;         /* digital alphabet                 */
  ESL_GETOPTS      *opts;        /* search specific options          */

  double            elapsed;     /* elapsed search time              */

  /* Structure created and populated by the individual threads.
   * The main thread is responsible for freeing up the memory.
   */
  P7_PIPELINE      *pli;         /* work pipeline                    */
  P7_TOPHITS       *th;          /* top hit results                  */
} WORKER_INFO;

static void free_QueueData(QUEUE_DATA *data);
static void pop_Queue(SEARCH_QUEUE *queue);
static void push_Queue(QUEUE_DATA *data, SEARCH_QUEUE *queue);
static QUEUE_DATA *read_Queue(SEARCH_QUEUE *queue);

static QUEUE_DATA *read_QueryCmd(int fd);

static int  setup_masterside_comm(ESL_GETOPTS *opts, SEQ_CACHE *seq_db, HMM_CACHE *hmm_db);
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
  { "--pid",        eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "write process id to file [default: /var/run/hmmpgmd.pid]",    12 },
  { "--daemon",     eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  NULL,            "run as a daemon using config file: /etc/hmmpgmd.conf",        12 },

  { "--seqdb",      eslARG_INFILE,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "protein database to cache for searches",                      12 },
  { "--hmmdb",      eslARG_INFILE,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "hmm database to cache for searches",                          12 },

  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>0", NULL,  NULL,  "--master",      "number of parallel CPU workers to use for multithreads",      12 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static ESL_OPTIONS searchOpts[] = {
  /* Control of output */
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                       2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
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
  { "--seqdb",      eslARG_INT,    NULL,  NULL, "n>0",   NULL,  NULL,  "--hmmdb",       "protein database to search",                                  12 },
  { "--hmmdb",      eslARG_INT,    NULL,  NULL, "n>0",   NULL,  NULL,  "--seqdb",       "hmm database to search",                                      12 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[options]";

static char banner[] = "search a query against a database";

#define BLOCK_SIZE 1000
static void search_thread(void *arg);
static void scan_thread(void *arg);

#define LOG_TO_STDOUT
#ifdef LOG_TO_STDOUT
void openlog(const char *ident, int option, int facility)
{
  /* do nothing */
  return;
}
void syslog(int priority, const char *format, ...)
{
  va_list ap;

  va_start(ap, format);
  vprintf(format, ap);
  va_end(ap);

  return;
}
void closelog(void)
{
  /* do nothing */
  return;
}
#endif

#define LOG_FATAL_MSG(str, err) {                                               \
    syslog(LOG_CRIT,"[%s:%d] - %s error %d - %s\n", __FILE__, __LINE__, str, err, strerror(err)); \
    exit(0); \
  }

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
  vsyslog(LOG_ERR, format, ap);

  /* send back an unsuccessful status message */
  n = sizeof(s);
  if (writen(fd, &s, n) != n) {
    syslog(LOG_ERR,"[%s:%d] - writing (%d) error %d - %s\n", __FILE__, __LINE__, fd, errno, strerror(errno));
    return;
  }
  if (writen(fd, ebuf, s.err_len) != s.err_len)  {
    syslog(LOG_ERR,"[%s:%d] - writing (%d) error %d - %s\n", __FILE__, __LINE__, fd, errno, strerror(errno));
    return;
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
write_pid(ESL_GETOPTS *go)
{
  FILE   *fp;

  pid_t   pid;

  char   *file    = NULL;
  char   *def_pid = "/var/run/hmmpgmd.pid";

  file = esl_opt_GetString(go, "--pid");
  if (file == NULL) file = def_pid;

  if ((fp = fopen(file, "w")) == NULL) {
  }

  fprintf(fp,"%ld\n", (long)getpid());
  fclose(fp);
}

static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go)
{
  int n;
  int status;
  ESL_GETOPTS *go = NULL;

  FILE *fp;
  const char *conf_file = "/etc/hmmpgmd.conf";

  if (go->argc == 1) {
  }

  if ((go = esl_getopts_Create(cmdlineOpts)) == NULL)    p7_Die("Internal failure creating options object");

  /* if there are no command line arguements, lets try and read /etc/hmmpgmd.conf
   * for any configuration data.
   */
  if (argc == 1) {
    if ((fp = fopen("/etc/hmmpgmd.conf", "r")) == NULL) {
      puts("Options --master or --worker must be specified.");
      esl_getopts_Destroy(go);
      return eslOK;
    }
    status = esl_opt_ProcessConfigfile(go, conf_file, fp);
    fclose(fp);

    if (status != eslOK) {
      printf("Failed to parse configuration file %s: %s\n",  conf_file, go->errbuf); 
      goto ERROR; 
    }
  } else {
    if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) { 
      printf("Failed to parse command line: %s\n",  go->errbuf); 
      goto ERROR; 
    }
  }

  if (esl_opt_VerifyConfig(go) != eslOK) { 
    printf("Failed to parse command line: %s\n", go->errbuf); 
    goto ERROR; 
  }
 
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
  if (n != 0) { puts("Incorrect number of command line arguments."); goto ERROR; }

  if (!esl_opt_IsUsed(go, "--seqdb") && !esl_opt_IsUsed(go, "--hmmdb")) {
    puts("At least one --seqdb or --hmmdb must be specified."); 
    goto ERROR;
  }

  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere most common options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(0);  
}

static int
process_searchline(int fd, char *cmdstr, ESL_GETOPTS *go)
{
  int status;

  if ((status = esl_opt_ProcessSpoof(go, cmdstr)) != eslOK) return status;
  if ((status = esl_opt_VerifyConfig(go))         != eslOK) return status;

  return eslOK;
}

static void
print_timings(int i, double elapsed, P7_PIPELINE *pli)
{
  char buf1[16];
  int h, m, s, hs;

  h  = (int) (elapsed / 3600.);
  m  = (int) (elapsed / 60.) - h * 60;
  s  = (int) (elapsed) - h * 3600 - m * 60;
  hs = (int) (elapsed * 100.) - h * 360000 - m * 6000 - s * 100;
  sprintf(buf1, "%02d:%02d.%02d", m,s,hs);

  fprintf (stdout, "%2d %9" PRId64 " %9" PRId64 " %7" PRId64 " %7" PRId64 " %6" PRId64 " %5" PRId64 " %s\n",
           i, pli->nseqs, pli->nres, pli->n_past_msv, pli->n_past_bias, pli->n_past_vit, pli->n_past_fwd, buf1);
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

void
worker_process(ESL_GETOPTS *go)
{
  FILE            *ofp        = stdout;            /* results output file (-o)                        */
  ESL_ALPHABET    *abc;                            /* digital alphabet                                */
  ESL_STOPWATCH   *w;                              /* timer used for profiling statistics             */
  int              status     = eslOK;
  int              i;
  int              fd;

  int              ncpus      = 0;
  WORKER_INFO     *info       = NULL;

  ESL_THREADS     *threadObj  = NULL;
  pthread_mutex_t  inx_mutex;
  int              current_index;

  SEQ_CACHE       *seq_db     = NULL;
  HMM_CACHE       *hmm_db     = NULL;
  QUEUE_DATA      *query      = NULL;

  /* Set processor specific flags */
  impl_Init();

  /* Initializations */
  p7_FLogsumInit();      /* we're going to use table-driven Logsum() approximations at times */

  w = esl_stopwatch_Create();
  abc = esl_alphabet_Create(eslAMINO);

  if (esl_opt_IsUsed(go, "--seqdb")) {
    char *name = esl_opt_GetString(go, "--seqdb");
    if ((status = cache_SeqDb(name, &seq_db)) != eslOK) p7_Fail("Failed to cache %s (%d)", name, status);
  }

  if (esl_opt_IsUsed(go, "--hmmdb")) {
    char *name = esl_opt_GetString(go, "--hmmdb");
    if ((status = cache_HmmDb(name, &hmm_db)) != eslOK) p7_Fail("Failed to cache %s (%d)", name, status);
  }

  /* start the communications with the web clients */
  fd = setup_masterside_comm(go, seq_db, hmm_db);

  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                   esl_threads_CPUCount(&ncpus);

  if (pthread_mutex_init(&inx_mutex, NULL) != 0) p7_Fail("mutex init failed");

  ESL_ALLOC(info, sizeof(*info) * ncpus);

  /* read query hmm/sequence */
  while ((query = read_QueryCmd(fd)) != NULL) {

    esl_stopwatch_Start(w);

    if (query->srch_type = HMMD_CMD_SEARCH) {
      threadObj = esl_threads_Create(&search_thread);
    } else {
      threadObj = esl_threads_Create(&scan_thread);
    }

    if (query->query_type == HMMD_SEQUENCE) {
      fprintf(ofp, "Query (%s):       %s  [L=%ld]\n", query->ip_addr, query->seq->name, (long) query->seq->n);
    } else {
      fprintf(ofp, "Query (%s):       %s  [M=%d]\n", query->ip_addr, query->hmm->name, query->hmm->M);
    }

    fprintf(ofp, "%s database %d [%d - %d]\n", 
            (query->srch_type == HMMD_CMD_SEARCH) ? "Sequence" : "HMM", 
            query->dbx, query->inx, query->inx + query->cnt - 1);

    /* Create processing pipeline and hit list */
    for (i = 0; i < ncpus; ++i) {
      info[i].abc   = query->abc;
      info[i].hmm   = query->hmm;
      info[i].seq   = query->seq;
      info[i].opts  = query->opts;

      info[i].th    = NULL;
      info[i].pli   = NULL;

      info[i].inx_mutex = inx_mutex;
      info[i].inx       = &current_index;

      if (query->srch_type == HMMD_CMD_SEARCH) {
        HMMER_SEQ **list = seq_db->db[query->dbx].list;
        info[i].sq_list   = &list[query->inx];
        info[i].sq_cnt    = query->cnt;
        info[i].om_list   = NULL;
        info[i].om_cnt    = 0;
      } else {
        info[i].sq_list   = NULL;
        info[i].sq_cnt    = 0;
        info[i].om_list   = hmm_db->list + query->inx;
        info[i].om_cnt    = query->cnt;
      }

      esl_threads_AddThread(threadObj, &info[i]);
    }

    current_index = 0;
    esl_threads_WaitForStart(threadObj);
    esl_threads_WaitForFinish(threadObj);

    esl_stopwatch_Stop(w);
#if 1
    fprintf (ofp, "   Sequences  Residues                              Elapsed\n");
    for (i = 0; i < ncpus; ++i) {
      print_timings(i, info[i].elapsed, info[i].pli);
    }
#endif
    /* merge the results of the search results */
    for (i = 1; i < ncpus; ++i) {
      p7_tophits_Merge(info[0].th, info[i].th);
      p7_pipeline_Merge(info[0].pli, info[i].pli);
      p7_pipeline_Destroy(info[i].pli);
      p7_tophits_Destroy(info[i].th);
    }

    print_timings(99, w->elapsed, info[0].pli);
    send_results(fd, w, info);

    /* free the last of the pipeline data */
    p7_pipeline_Destroy(info->pli);
    p7_tophits_Destroy(info->th);

    free_QueueData(query);

    esl_threads_Destroy(threadObj);
  } /* end outer loop over query HMMs */

  pthread_mutex_destroy(&inx_mutex);

  if (hmm_db != NULL) cache_HmmDestroy(hmm_db);
  if (seq_db != NULL) cache_SeqDestroy(seq_db);

  free(info);

  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);

  esl_getopts_Destroy(go);
  esl_alphabet_Destroy(abc);

  return;

 ERROR:
  LOG_FATAL_MSG("malloc", errno);
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

  SEQ_CACHE          *seq_db     = NULL;
  HMM_CACHE          *hmm_db     = NULL;
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

  if (esl_opt_IsUsed(go, "--seqdb")) {
    char *name = esl_opt_GetString(go, "--seqdb");
    if ((status = cache_SeqDb(name, &seq_db)) != eslOK) p7_Fail("Failed to cache %s (%d)", name, status);
  }

  if (esl_opt_IsUsed(go, "--hmmdb")) {
    char *name = esl_opt_GetString(go, "--hmmdb");
    if ((status = cache_HmmDb(name, &hmm_db)) != eslOK) p7_Fail("Failed to cache %s (%d)", name, status);
  }

  /* initialize the search queue */
  ESL_ALLOC(queue, sizeof(SEARCH_QUEUE));
  if ((n = pthread_mutex_init(&queue->queue_mutex, NULL)) != 0)  LOG_FATAL_MSG("mutex init", n);
  if ((n = pthread_cond_init(&queue->queue_cond, NULL)) != 0)    LOG_FATAL_MSG("cond init", n);

  queue->head = NULL;
  queue->tail = NULL;

  /* start the communications with the web clients */
  client_comm.queue = queue;
  client_comm.seq_db = seq_db;
  client_comm.hmm_db = hmm_db;
  setup_clientside_comm(go, &client_comm);

  /* initialize the worker structure */
  if ((n = pthread_mutex_init(&worker_comm.work_mutex, NULL)) != 0)   LOG_FATAL_MSG("mutex init", n);
  if ((n = pthread_cond_init(&worker_comm.start_cond, NULL)) != 0)    LOG_FATAL_MSG("cond init", n);
  if ((n = pthread_cond_init(&worker_comm.complete_cond, NULL)) != 0) LOG_FATAL_MSG("cond init", n);
  if ((n = pthread_cond_init(&worker_comm.worker_cond, NULL)) != 0)   LOG_FATAL_MSG("cond init", n);

  worker_comm.sock_fd = -1;
  worker_comm.head    = NULL;
  worker_comm.tail    = NULL;
  worker_comm.seq_db  = seq_db;
  worker_comm.hmm_db  = hmm_db;

  setup_workerside_comm(go, &worker_comm);

  /* read query hmm/sequence */
  while ((query = read_Queue(queue)) != NULL) {

    printf("Processing query from %s\n", query->ip_addr);

    if ((n = pthread_mutex_lock (&worker_comm.work_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

    /* wait until we have atleast one worker */
    while (worker_comm.head == NULL) {
      if ((n = pthread_cond_wait(&worker_comm.worker_cond, &worker_comm.work_mutex)) != 0) LOG_FATAL_MSG("cond wait", n);
    }

    /* add to all workers the new query to process */
    worker_comm.started = 0;
    worker = worker_comm.head;
    while (worker != NULL) {
      ++worker_comm.started;

      if ((worker->cmd = malloc(query->cmd->length)) == NULL) LOG_FATAL_MSG("malloc", errno);
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
    if (query->srch_type == HMMD_CMD_SEARCH) {
      cnt = seq_db->db[query->dbx].count;
    } else {
      cnt = hmm_db->count;
    }
    worker = worker_comm.head;
    while (worker != NULL) {
      if (worker->cmd != NULL) {
        worker->cmd->inx = inx;
        worker->cmd->cnt = cnt / i;

        inx += worker->cmd->cnt;
        cnt -= worker->cmd->cnt;
        --i;

        worker = worker->next;
      }
    }

    worker_comm.completed = worker_comm.started;
    worker_comm.errors    = 0;

    /* notify all the worker threads of the new query */
    if ((n = pthread_cond_broadcast(&worker_comm.start_cond)) != 0) LOG_FATAL_MSG("cond broadcast", n);
    if ((n = pthread_mutex_unlock (&worker_comm.work_mutex)) != 0)  LOG_FATAL_MSG("mutex unlock", n);

    esl_stopwatch_Start(w);

    /* Wait for all the workers to complete */
    if ((n = pthread_mutex_lock (&worker_comm.work_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

    while (worker_comm.completed > 0) {
      if ((n = pthread_cond_wait (&worker_comm.complete_cond, &worker_comm.work_mutex)) != 0) LOG_FATAL_MSG("cond wait", n);
    }

    if ((n = pthread_mutex_unlock (&worker_comm.work_mutex)) != 0) LOG_FATAL_MSG("mutex unlock", n);

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

  esl_stopwatch_Destroy(w);

  if (hmm_db != NULL) cache_HmmDestroy(hmm_db);
  if (seq_db != NULL) cache_SeqDestroy(seq_db);

  if (ofp != stdout) fclose(ofp);

  esl_getopts_Destroy(go);

  free(info);
  free(queue);

  return;

 ERROR:
  LOG_FATAL_MSG("malloc", errno);
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go = NULL;      /* command line processing */

  process_commandline(argc, argv, &go);

  /* if we write to a broken socket, ignore the signal and handle the error. */
  signal(SIGPIPE, SIG_IGN);

  /* check if we need to write out our pid */
  if (esl_opt_IsOn(go, "--pid")) write_pid(go);

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

  if ((n = pthread_mutex_lock (&queue->queue_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

  while (queue->head == NULL) {
    if ((n = pthread_cond_wait (&queue->queue_cond, &queue->queue_mutex)) != 0) LOG_FATAL_MSG("cond wait", n);
  }

  data = queue->head;

  if ((n = pthread_mutex_unlock (&queue->queue_mutex)) != 0) LOG_FATAL_MSG("mutex unlock", n);

  return data;
}

static void
pop_Queue(SEARCH_QUEUE *queue)
{
  int         n;
  QUEUE_DATA *data;

  if ((n = pthread_mutex_lock (&queue->queue_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

  data = queue->head;
  queue->head = data->next;
  if (queue->head == NULL) queue->tail = NULL;

  if ((n = pthread_mutex_unlock (&queue->queue_mutex)) != 0) LOG_FATAL_MSG("mutex unlock", n);

  free_QueueData(data);
}

static void
push_Queue(QUEUE_DATA *data, SEARCH_QUEUE *queue)
{
  int n;

  /* add the search request to the queue */
  if ((n = pthread_mutex_lock (&queue->queue_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

  if (queue->head == NULL) {
    queue->head = data;
    queue->tail = data;
  } else {
    queue->tail->next = data;
    data->prev = queue->tail;
    queue->tail = data;
  }

  /* if anyone is waiting, wake them up */
  if ((n = pthread_cond_broadcast (&queue->queue_cond)) != 0) LOG_FATAL_MSG("cond broadcast", n);
  if ((n = pthread_mutex_unlock (&queue->queue_mutex)) != 0)  LOG_FATAL_MSG("mutex unlock", n);
}

static void
remove_Queue(int fd, SEARCH_QUEUE *queue)
{
  int n;
  QUEUE_DATA *data   = NULL;
  QUEUE_DATA *next   = NULL;
  QUEUE_DATA *list   = NULL;

  if ((n = pthread_mutex_lock (&queue->queue_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

  /* skip over the first queued object since it is being searched */
  data = queue->head;
  if (data != NULL) data = data->next;

  /* remove all search request from the queue for this socket */
  while (data != NULL) {
    next = data->next;
    if (data->sock == fd) {
      data->prev->next = data->next;
      if (data->next == NULL) {
        queue->tail = data->prev;
      } else {
        data->next->prev = data->prev;
      }
      data->next = list;
      list = data;
    }
    data = next;
  }

  /* unlock the queue */
  if ((n = pthread_mutex_unlock (&queue->queue_mutex)) != 0)  LOG_FATAL_MSG("mutex unlock", n);

  /* free all the removed queries */
  while (list != NULL) {
    data = list;
    list = list->next;
    free_QueueData(data);
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
  if ((size = readn(fd, &hdr, n)) == -1)  LOG_FATAL_MSG("read", errno);

  if ((cmd = malloc(hdr.length)) == NULL) LOG_FATAL_MSG("malloc", errno);

  *cmd = hdr;
  p = (char *)cmd;
  p += sizeof(hdr);

  n = hdr.length - sizeof(hdr);
  if ((size = readn(fd, p, n)) == -1)               LOG_FATAL_MSG("read", errno);

  if ((query = malloc(sizeof(QUEUE_DATA))) == NULL) LOG_FATAL_MSG("malloc", errno);

  query->srch_type  = cmd->command;
  query->query_type = cmd->query_type;
  query->dbx        = cmd->db_inx;
  query->inx        = cmd->inx;
  query->cnt        = cmd->cnt;
  query->sock       = fd;
  query->retry      = 0;
  query->retry_cnt  = 0;
  query->cmd        = NULL;
  query->next       = NULL;
  query->prev       = NULL;

  /* process search specific options */
  query->opts = esl_getopts_Create(searchOpts);
  if (query->opts == NULL)  LOG_FATAL_MSG("esl_getopts_Create", EINVAL);
  process_searchline(fd, p, query->opts);

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
    P7_HMM  thmm;
    P7_HMM *hmm = p7_hmm_CreateShell();

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
search_thread(void *arg)
{
  int               i;
  int               count;
  int               seed;
  int               status;
  int               workeridx;
  WORKER_INFO      *info;
  ESL_THREADS      *obj;

  ESL_SQ            dbsq;
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

  /* set up the dummy description and accession fields */
  dbsq.desc = "";
  dbsq.acc  = "";

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
    p7_ProfileConfig(info->hmm, bg, gm, 100, p7_LOCAL);
    p7_oprofile_Convert(gm, om);
  }

  /* Create processing pipeline and hit list */
  th  = p7_tophits_Create(); 
  pli = p7_pipeline_Create(info->opts, om->M, 100, FALSE, p7_SEARCH_SEQS);
  p7_pli_NewModel(pli, om, bg);

  /* loop until all sequences have been processed */
  count = 1;
  while (count > 0) {
    int          inx;
    HMMER_SEQ  **sq;

    /* grab the next block of sequences */
    if (pthread_mutex_lock(&info->inx_mutex) != 0) p7_Fail("mutex lock failed");
    inx = *info->inx;
    *info->inx += BLOCK_SIZE;
    if (pthread_mutex_unlock(&info->inx_mutex) != 0) p7_Fail("mutex unlock failed");

    sq = info->sq_list + inx;

    count = info->sq_cnt - inx;
    if (count > BLOCK_SIZE) count = BLOCK_SIZE;
    //printf("THREAD %08x: %d %d\n", workeridx, inx, count);

    /* Main loop: */
    for (i = 0; i < count; ++i, ++sq) {

      dbsq.name  = (*sq)->name;
      dbsq.dsq   = (*sq)->dsq;
      dbsq.n     = (*sq)->n;
      dbsq.idx   = (*sq)->idx;

      p7_pli_NewSeq(pli, &dbsq);
      p7_bg_SetLength(bg, dbsq.n);
      p7_oprofile_ReconfigLength(om, dbsq.n);

      p7_Pipeline(pli, om, bg, &dbsq, th);

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
scan_thread(void *arg)
{
  int               i;
  int               count;
  int               workeridx;
  WORKER_INFO      *info;
  ESL_THREADS      *obj;

  ESL_STOPWATCH    *w;

  P7_BG            *bg       = NULL;         /* null model                     */
  P7_PIPELINE      *pli      = NULL;         /* work pipeline                  */
  P7_TOPHITS       *th       = NULL;         /* top hit results                */

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  w = esl_stopwatch_Create();
  esl_stopwatch_Start(w);

  /* Convert to an optimized model */
  bg = p7_bg_Create(info->abc);

  /* Create processing pipeline and hit list */
  th  = p7_tophits_Create(); 
  pli = p7_pipeline_Create(info->opts, 100, 100, FALSE, p7_SCAN_MODELS);

  p7_pli_NewSeq(pli, info->seq);
  p7_bg_SetLength(bg, info->seq->n);

  /* loop until all sequences have been processed */
  count = 1;
  while (count > 0) {
    int           inx;
    P7_OPROFILE  *om;

    /* grab the next block of sequences */
    if (pthread_mutex_lock(&info->inx_mutex) != 0) p7_Fail("mutex lock failed");
    inx = *info->inx;
    *info->inx += BLOCK_SIZE;
    if (pthread_mutex_unlock(&info->inx_mutex) != 0) p7_Fail("mutex unlock failed");

    om = info->om_list + inx;

    count = info->om_cnt - inx;
    if (count > BLOCK_SIZE) count = BLOCK_SIZE;
    //printf("THREAD %08x: %d %d\n", workeridx, inx, count);

    /* Main loop: */
    for (i = 0; i < count; ++i, ++om) {
      p7_pli_NewModel(pli, om, bg);
      p7_oprofile_ReconfigLength(om, info->seq->n);
	      
      p7_Pipeline(pli, om, bg, info->seq, th);
	      
      p7_oprofile_Destroy(om);
      p7_pipeline_Reuse(pli);
    }
  }

  /* make available the pipeline objects to the main thread */
  info->th = th;
  info->pli = pli;

  /* clean up */
  p7_bg_Destroy(bg);

  esl_stopwatch_Stop(w);
  info->elapsed = w->elapsed;

  esl_stopwatch_Destroy(w);

  esl_threads_Finished(obj, workeridx);

  pthread_exit(NULL);
  return;
}

typedef struct {
  P7_HIT         *hit;
  char           *data;
} HIT_LIST;

static int
hit_sorter(const void *p1, const void *p2)
{
  int cmp;

  const P7_HIT *h1 = ((HIT_LIST *) p1)->hit;
  const P7_HIT *h2 = ((HIT_LIST *) p2)->hit;

  cmp  = (h1->sortkey < h2->sortkey);
  cmp -= (h1->sortkey > h2->sortkey);

  return cmp;
}

static void
forward_results(QUEUE_DATA *query, ESL_STOPWATCH *w, WORKERSIDE_ARGS *comm)
{
  int fd;
  int i, j;
  int inx;
  int n;

  int total = 0;

  HIT_LIST           *list;
  HIT_LIST           *ptr;

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
  if ((n = pthread_mutex_lock (&comm->work_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

  /* count the number of hits */
  n = 0;
  worker = comm->head;
  while (worker != NULL) {
      if (worker->completed) n += worker->stats.nhits;
      worker = worker->next;
  }

  /* allocate spaces to hold all the hits */
  if ((list = malloc(sizeof(HIT_LIST) * n)) == NULL) LOG_FATAL_MSG("malloc", errno);
  ptr = list;

  /* gather up the results and merge them */
  inx = 0;
  worker = comm->head;
  while (worker != NULL) {
    if (worker->completed) {
      stats.nmodels      = worker->stats.nmodels;
      stats.nseqs       += worker->stats.nseqs;
      stats.nres        += worker->stats.nres;
      stats.nnodes       = worker->stats.nnodes;

      stats.nhits       += worker->stats.nhits;
      stats.nreported   += worker->stats.nreported;
      stats.nincluded   += worker->stats.nincluded;

      stats.n_past_msv  += worker->stats.n_past_msv;
      stats.n_past_bias += worker->stats.n_past_bias;
      stats.n_past_vit  += worker->stats.n_past_vit;
      stats.n_past_fwd  += worker->stats.n_past_fwd;

      stats.Z_setby      = worker->stats.Z_setby;
      stats.domZ_setby   = worker->stats.domZ_setby;
      stats.domZ         = worker->stats.domZ;
      stats.Z            = worker->stats.Z;

      if (worker->stats.nhits > 0) {
        for (i = 0; i < worker->stats.nhits; ++i) {
          ptr->hit  = worker->hit[i];
          ptr->data = worker->hit_data[i];
          ++ptr;
        }
        memset(worker->hit, 0, sizeof(char *) * worker->stats.nhits);
        memset(worker->hit_data, 0, sizeof(char *) * worker->stats.nhits);

        free(worker->hit);
        free(worker->hit_data);

        worker->hit      = NULL;
        worker->hit_data = NULL;
      }

      worker->completed = 0;
    }

    worker = worker->next;
  }

  if ((n = pthread_mutex_unlock (&comm->work_mutex)) != 0) LOG_FATAL_MSG("mutex unlock", n);

  if (comm->seq_db != NULL) {
    stats.nseqs = comm->seq_db->db[query->dbx].K;
  }
    
  if (stats.Z_setby == p7_ZSETBY_NTARGETS) {
    stats.Z = (query->srch_type == HMMD_CMD_SEARCH) ? stats.nseqs : stats.nmodels;
  }

  /* sort the hits and apply score and E-value thresholds */
  if (stats.nhits > 0) {
    P7_PIPELINE     *pli;
    P7_TOPHITS       th;

    th.unsrt     = NULL;
    th.N         = stats.nhits;
    th.nreported = 0;
    th.nincluded = 0;
    th.is_sorted = 0;
      
    pli = p7_pipeline_Create(query->opts, 100, 100, FALSE, p7_SEARCH_SEQS);
    pli->nmodels     = stats.nmodels;
    pli->nseqs       = stats.nseqs;
    pli->nres        = stats.nres;
    pli->nnodes      = stats.nnodes;
    pli->n_past_msv  = stats.n_past_msv;
    pli->n_past_bias = stats.n_past_bias;
    pli->n_past_vit  = stats.n_past_vit;
    pli->n_past_fwd  = stats.n_past_fwd;

    pli->Z           = stats.Z;
    pli->domZ        = stats.domZ;
    pli->Z_setby     = stats.Z_setby;
    pli->domZ_setby  = stats.domZ_setby;

    if ((th.hit = malloc(sizeof(P7_HIT *) * stats.nhits)) == NULL) LOG_FATAL_MSG("malloc", errno);

    qsort(list, stats.nhits, sizeof(HIT_LIST), hit_sorter);

    for (i = 0; i < th.N; ++i) th.hit[i] = list[i].hit;
    p7_tophits_Threshold(&th, pli);

    stats.nreported = th.nreported;
    stats.nincluded = th.nreported;
    stats.domZ      = pli->domZ;
    stats.Z         = pli->Z;

    memset(th.hit, 0, sizeof(P7_HIT *) * th.N);
    free(th.hit);
  }

  /* send back a successful status message */
  status.status  = eslOK;
  status.err_len = 0;
  n = sizeof(status);
  total += n;
  if (writen(fd, &status, n) != n) {
    syslog(LOG_ERR,"[%s:%d] - writing %s error %d - %s\n", __FILE__, __LINE__, query->ip_addr, errno, strerror(errno));
    goto CLEAR;
  }

  n = sizeof(stats);
  total += n;
  if (writen(fd, &stats, n) != n) {
    syslog(LOG_ERR,"[%s:%d] - writing %s error %d - %s\n", __FILE__, __LINE__, query->ip_addr, errno, strerror(errno));
    goto CLEAR;
  }

  /* loop through the hit list sending to dest */
  for (i = 0; i < stats.nhits; ++i) {
    P7_HIT *hit = list[i].hit;

    n = sizeof(P7_HIT);
    total += n;
    if (writen(fd, hit, n) != n) {
      syslog(LOG_ERR,"[%s:%d] - writing %s error %d - %s\n", __FILE__, __LINE__, query->ip_addr, errno, strerror(errno));
      goto CLEAR;
    }

    if (list[i].data != NULL) {
      /* the hit string pointers contain the length of the strings */
      n = 0;
      if (hit->name != NULL) n += (socklen_t)(hit->name - (char *)NULL);
      if (hit->acc  != NULL) n += (socklen_t)(hit->acc  - (char *)NULL);
      if (hit->desc != NULL) n += (socklen_t)(hit->desc - (char *)NULL);

      total += n;
      if (writen(fd, list[i].data, n) != n) {
        syslog(LOG_ERR,"[%s:%d] - writing %s error %d - %s\n", __FILE__, __LINE__, query->ip_addr, errno, strerror(errno));
        goto CLEAR;
      }
    }

    /* send the domains */
    for (j = 0; j < hit->ndom; ++j) {
      P7_DOMAIN      *dcl = &hit->dcl[j];

      n = (socklen_t)sizeof(P7_DOMAIN);
      total += n;
      if (writen(fd, dcl, n) != n) {
        syslog(LOG_ERR,"[%s:%d] - writing %s error %d - %s\n", __FILE__, __LINE__, query->ip_addr, errno, strerror(errno));
        goto CLEAR;
      }

      n = (socklen_t)sizeof(P7_ALIDISPLAY);
      total += n;
      if (writen(fd, dcl->ad, n) != n) {
        syslog(LOG_ERR,"[%s:%d] - writing %s error %d - %s\n", __FILE__, __LINE__, query->ip_addr, errno, strerror(errno));
        goto CLEAR;
      }

      n = (socklen_t)dcl->ad->memsize;
      total += n;
      if (writen(fd, dcl->ad->mem, n) != n) {
        syslog(LOG_ERR,"[%s:%d] - writing %s error %d - %s\n", __FILE__, __LINE__, query->ip_addr, errno, strerror(errno));
        goto CLEAR;
      }
    }
  }

  printf("Results for %s (%d) sent %d bytes\n", query->ip_addr, fd, total);
  printf("Hits:%"PRId64 "  reported:%"PRId64 "  included:%"PRId64 "\n", stats.nhits, stats.nreported, stats.nincluded);

  /* free all the data */
 CLEAR:
  for (i = 0; i < stats.nhits; ++i) {
    P7_HIT *hit = list[i].hit;
    for (j = 0; j < hit->ndom; ++j) {
      P7_DOMAIN *dcl = &hit->dcl[j];
      free (dcl->ad);
      free (dcl->ad->mem);
    }
    free(hit->dcl);
    free(hit);
    
    if (list[i].data != NULL) free(list[i].data);
  }

  memset(list, 0, sizeof(HIT_LIST) * stats.nhits);
  free(list);
}

static void
clear_results(WORKERSIDE_ARGS *comm)
{
  int i, j;
  int inx;
  int n;

  WORKER_DATA *worker;

  /* lock the workers until we have freed the results */
  if ((n = pthread_mutex_lock (&comm->work_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

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

  if ((n = pthread_mutex_unlock (&comm->work_mutex)) != 0)  LOG_FATAL_MSG("mutex unlock", n);
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
  if (writen(fd, &status, n) != n) LOG_FATAL_MSG("write", errno);

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
  if (writen(fd, &stats, n) != n) LOG_FATAL_MSG("write", errno);
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

    if ((h = malloc(n)) == NULL) LOG_FATAL_MSG("malloc", errno);
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
    if (writen(fd, h, n) != n) LOG_FATAL_MSG("write", errno);

    /* send the domains for this hit */
    dcl = hit->dcl;
    for (j = 0; j < hit->ndom; ++j) {
      n = sizeof(P7_DOMAIN);
      total += n;
      if (writen(fd, dcl, n) != n) LOG_FATAL_MSG("write", errno);

      n = sizeof(P7_ALIDISPLAY);
      total += n;
      if (writen(fd, dcl->ad, n) != n) LOG_FATAL_MSG("write", errno);

      n = dcl->ad->memsize;
      total += n;
      if (writen(fd, dcl->ad->mem, n) != n) LOG_FATAL_MSG("write", errno);

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
setup_masterside_comm(ESL_GETOPTS *opts, SEQ_CACHE *seq_db, HMM_CACHE *hmm_db)
{
  int                  n;
  int                  fd;
  struct sockaddr_in   addr;

  HMMD_INIT_CMD        cmd;

  /* create a reliable, stream socket using TCP */
  if ((fd = socket(AF_INET, SOCK_STREAM, 0)) < 0) LOG_FATAL_MSG("socket", errno);

  /* construct the server address structure */
  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_port   = htons(esl_opt_GetInteger(opts, "--wport"));
  if ((inet_pton(AF_INET, esl_opt_GetString(opts, "--worker"), &addr.sin_addr)) < 0) LOG_FATAL_MSG("inet pton", errno);

  /* establish the connection to the master server */
  if (connect(fd, (struct sockaddr *) &addr, sizeof(addr)) < 0) LOG_FATAL_MSG("connect", errno);

  /* read the database information */
  n = sizeof(cmd);
  if (readn(fd, &cmd, n) == -1) {
      syslog(LOG_ERR,"[%s:%d] - reading %s error %d - %s\n", __FILE__, __LINE__, inet_ntoa(addr.sin_addr), errno, strerror(errno));
      LOG_FATAL_MSG("reading init command", errno);
  }

  if (cmd.command != HMMD_CMD_INIT || cmd.length != sizeof(cmd)) {
      syslog(LOG_ERR,"[%s:%d] - malformed init cmd from %s %d %d\n", __FILE__, __LINE__, inet_ntoa(addr.sin_addr), cmd.command, cmd.length);
      LOG_FATAL_MSG("malformed message", 0);
  }

  if (seq_db != NULL) {
    cmd.sid[MAX_INIT_DESC-1] = 0;
    if (strcmp (cmd.sid, seq_db->id) != 0 || cmd.db_cnt != seq_db->db_cnt || cmd.seq_cnt != seq_db->count) {
      syslog(LOG_ERR,"[%s:%d] - database integrity error %s - %s\n", __FILE__, __LINE__, cmd.sid, seq_db->id);
      LOG_FATAL_MSG("database integrity error", 0);
    }
  }

  return fd;    
}


static int
clientside_loop(CLIENTSIDE_ARGS *data)
{
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
  SEQ_CACHE         *seq_db  = data->seq_db;
  QUEUE_DATA        *parms;

  jmp_buf            jmp_env;

  buf_size = MAX_BUFFER;
  if ((buffer  = malloc(buf_size))   == NULL) LOG_FATAL_MSG("malloc", errno);
  if ((opt_str = malloc(MAX_BUFFER)) == NULL) LOG_FATAL_MSG("malloc", errno);

  ptr = buffer;
  remaining = buf_size;
  amount = 0;

  eod = 0;
  while (!eod) {

    /* Receive message from client */
    if ((n = read(data->sock_fd, ptr, remaining)) < 0) {
      syslog(LOG_ERR,"[%s:%d] - reading %s error %d - %s\n", __FILE__, __LINE__, data->ip_addr, errno, strerror(errno));
      return 1;
    }

    if (n == 0) return 1;

    ptr += n;
    amount += n;
    remaining -= n;
    eod = (amount > 1 && *(ptr - 2) == '/' && *(ptr - 2) == '/' );

    /* if the buffer is full, make it larger */
    if (!eod && remaining == 0) {
      if ((buffer = realloc(buffer, buf_size * 2)) == NULL) LOG_FATAL_MSG("realloc", errno);
      ptr = buffer + buf_size;
      remaining = buf_size;
      buf_size *= 2;
    }
  }

  /* zero terminate the buffer */
  if (remaining == 0) {
    if ((buffer = realloc(buffer, buf_size + 1)) == NULL) LOG_FATAL_MSG("realloc", errno);
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
  if (*ptr != '@') {
    client_error(data->sock_fd, eslEFORMAT, "Missing options string");
    free(buffer);
    return 0;
  } else {
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
    
    status = process_searchline(data->sock_fd, opt_str, opts);
    if (status != eslOK) {
      client_error_longjmp(data->sock_fd, status, &jmp_env, "Failed to parse options string: %s", opts->errbuf);
    }

    /* the options string can handle an optional database */
    if (esl_opt_ArgNumber(opts) > 0) {
      client_error_longjmp(data->sock_fd, status, &jmp_env, "Incorrect number of command line arguments.");
    }

    if (esl_opt_IsUsed(opts, "--seqdb")) {
      dbx = esl_opt_GetInteger(opts, "--seqdb");
      if (dbx < 1 || dbx > seq_db->db_cnt) {
        client_error_longjmp(data->sock_fd, eslEINVAL, &jmp_env, "Database out of range (1 - %d).", seq_db->db_cnt);
      }
    } else if (esl_opt_IsUsed(opts, "--hmmdb")) {
      dbx = esl_opt_GetInteger(opts, "--hmmdb");
      if (dbx != 1) {
        client_error_longjmp(data->sock_fd, eslEINVAL, &jmp_env, "Database out of range (1 - 1).");
      }
    } else {
      client_error_longjmp(data->sock_fd, eslEINVAL, &jmp_env, "No search database specified, --seqdb or --hmmdb.");
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

      if (esl_opt_IsUsed(opts, "--hmmdb")) {
        client_error_longjmp(data->sock_fd, status, &jmp_env, "A HMM cannot be used to search a hmm database");
      }

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

  if ((parms = malloc(sizeof(QUEUE_DATA))) == NULL) LOG_FATAL_MSG("malloc", errno);

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

  if ((cmd = malloc(n)) == NULL) LOG_FATAL_MSG("malloc", errno);

  cmd->length       = n;
  cmd->command      = (esl_opt_IsUsed(opts, "--seqdb")) ? HMMD_CMD_SEARCH : HMMD_CMD_SCAN;
  cmd->db_inx       = dbx - 1;   /* the program indexes databases 0 .. n-1 */
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
  parms->dbx  = dbx - 1;
  parms->cmd  = cmd;

  strcpy(parms->ip_addr, data->ip_addr);
  parms->sock = data->sock_fd;
  parms->next = NULL;
  parms->prev = NULL;

  parms->retry      = 0;
  parms->retry_cnt  = 0;

  parms->srch_type  = cmd->command;
  parms->query_type = (seq != NULL) ? HMMD_SEQUENCE : HMMD_HMM;

  push_Queue(parms, queue);

  if (parms->seq != NULL) {
    printf("Queued sequence %s from %s (%d)\n", parms->seq->name, parms->ip_addr, parms->sock);
  } else {
    printf("Queued hmm %s from %s (%d)\n", parms->hmm->name, parms->ip_addr, parms->sock);
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

  remove_Queue(data->sock_fd, data->queue);

  printf("Closing %s (%d)\n", data->ip_addr, data->sock_fd);

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
    if ((fd = accept(data->sock_fd, (struct sockaddr *)&addr, (unsigned int *)&n)) < 0) LOG_FATAL_MSG("accept", errno);

    if ((targs = malloc(sizeof(CLIENTSIDE_ARGS))) == NULL) LOG_FATAL_MSG("malloc", errno);
    targs->queue      = queue;
    targs->sock_fd    = fd;
    targs->seq_db     = data->seq_db;
    targs->hmm_db     = data->hmm_db;

    strncpy(targs->ip_addr, inet_ntoa(addr.sin_addr), INET_ADDRSTRLEN);
    targs->ip_addr[INET_ADDRSTRLEN-1] = 0;

    if ((n = pthread_create(&thread_id, NULL, clientside_thread, targs)) != 0) LOG_FATAL_MSG("thread create", n);
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
  if ((sock_fd = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) LOG_FATAL_MSG("socket", errno);
      
  /* incase the server went down in an ungraceful way, allow the port to be
   * reused avoiding the timeout.
   */
  reuse = 1;
  if (setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, (char *)&reuse, sizeof(reuse)) < 0) LOG_FATAL_MSG("setsockopt", errno);

  /* the sockets are never closed, so if the server exits, force the kernel to
   * close the socket and clear it so the server can be restarted immediately.
   */
  linger.l_onoff = 1;
  linger.l_linger = 0;
  if (setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, (char *)&linger, sizeof(linger)) < 0) LOG_FATAL_MSG("setsockopt", errno);

  /* Construct local address structure */
  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_ANY);
  addr.sin_port = htons(esl_opt_GetInteger(opts, "--cport"));

  /* Bind to the local address */
  if (bind(sock_fd, (struct sockaddr *) &addr, sizeof(addr)) < 0) LOG_FATAL_MSG("bind", errno);

  /* Mark the socket so it will listen for incoming connections */
  if (listen(sock_fd, esl_opt_GetInteger(opts, "--ccncts")) < 0) LOG_FATAL_MSG("listen", errno);
  args->sock_fd = sock_fd;

  if ((n = pthread_create(&thread_id, NULL, client_comm_thread, (void *)args)) != 0) LOG_FATAL_MSG("socket", n);
}

static void
workerside_loop(WORKERSIDE_ARGS *data, WORKER_DATA *worker)
{
  int          n;
  int          i, j;
  int          size;
  int          total;

  P7_HIT      *hit     = NULL;
  P7_DOMAIN   *dcl     = NULL;

  ESL_STOPWATCH   *w;

  HMMD_SEARCH_STATS  *stats;

  w = esl_stopwatch_Create();

  for ( ; ; ) {

    /* wait for the next search object */
    if ((n = pthread_mutex_lock (&data->work_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

    /* wait for the master's signal to start the calculations */
    while (worker->cmd == NULL) {
      if ((n = pthread_cond_wait(&data->start_cond, &data->work_mutex)) != 0) LOG_FATAL_MSG("cond wait", n);
    }

    if ((n = pthread_mutex_unlock (&data->work_mutex)) != 0) LOG_FATAL_MSG("mutex unlock", n);

    esl_stopwatch_Start(w);

    n = worker->cmd->length;
    if (writen(worker->sock_fd, worker->cmd, n) != n) {
      syslog(LOG_ERR,"[%s:%d] - writing %s error %d - %s\n", __FILE__, __LINE__, worker->ip_addr, errno, strerror(errno));
      break;
    }

    total = 0;
    worker->total = 0;

    n = sizeof(worker->status);
    total += n;
    if ((size = readn(worker->sock_fd, &worker->status, n)) == -1) {
      syslog(LOG_ERR,"[%s:%d] - reading %s error %d - %s\n", __FILE__, __LINE__, worker->ip_addr, errno, strerror(errno));
      break;
    }

    if (worker->status.status != eslOK) {
      n = worker->status.err_len;
      total += n; 
      if ((worker->err_buf = malloc(n)) == NULL) LOG_FATAL_MSG("malloc", errno);
      worker->err_buf[0] = 0;
      if ((size = readn(worker->sock_fd, worker->err_buf, n)) == -1) {
        syslog(LOG_ERR,"[%s:%d] - reading %s error %d - %s\n", __FILE__, __LINE__, worker->ip_addr, errno, strerror(errno));
        break;
      }
    } else {

      n = sizeof(worker->stats);
      total += n;
      if ((size = readn(worker->sock_fd, &worker->stats, n)) == -1) {
        syslog(LOG_ERR,"[%s:%d] - reading %s error %d - %s\n", __FILE__, __LINE__, worker->ip_addr, errno, strerror(errno));
        break;
      }

      stats = &worker->stats;
      if ((worker->hit      = malloc(sizeof(char *) * stats->nhits)) == NULL) LOG_FATAL_MSG("malloc", errno);
      if ((worker->hit_data = malloc(sizeof(char *) * stats->nhits)) == NULL) LOG_FATAL_MSG("malloc", errno);

      /* loop through the hit list sending to dest */
      for (i = 0; i < stats->nhits; ++i) {
        n = sizeof(P7_HIT);
        if ((hit = malloc(n)) == NULL) LOG_FATAL_MSG("malloc", errno);
        worker->hit[i] = hit;
        total += n;
        if ((size = readn(worker->sock_fd, hit, n)) == -1) {
          syslog(LOG_ERR,"[%s:%d] - reading %s error %d - %s\n", __FILE__, __LINE__, worker->ip_addr, errno, strerror(errno));
          break;
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
          if ((worker->hit_data[i] = malloc(n)) == NULL) LOG_FATAL_MSG("malloc", errno);
          if ((size = readn(worker->sock_fd, worker->hit_data[i], n)) == -1) {
            syslog(LOG_ERR,"[%s:%d] - reading %s error %d - %s\n", __FILE__, __LINE__, worker->ip_addr, errno, strerror(errno));
            break;
          }
        }

        n = sizeof(P7_DOMAIN) * hit->ndom;
        if ((hit->dcl = malloc(n)) == NULL) LOG_FATAL_MSG("malloc", errno);

        /* read the domains for this hit */
        dcl = hit->dcl;
        for (j = 0; j < hit->ndom; ++j) {
          char *base;
          P7_ALIDISPLAY *ad = NULL;

          n = (socklen_t)sizeof(P7_DOMAIN);
          total += n;
          if ((size = readn(worker->sock_fd, dcl, n)) == -1) {
            syslog(LOG_ERR,"[%s:%d] - reading %s error %d - %s\n", __FILE__, __LINE__, worker->ip_addr, errno, strerror(errno));
            break;
          }

          n = (socklen_t)sizeof(P7_ALIDISPLAY);
          total += n;
          if ((dcl->ad = malloc(n)) == NULL) LOG_FATAL_MSG("malloc", errno);
          if ((size = readn(worker->sock_fd, dcl->ad, n)) == -1) {
            syslog(LOG_ERR,"[%s:%d] - reading %s error %d - %s\n", __FILE__, __LINE__, worker->ip_addr, errno, strerror(errno));
            break;
          }

          /* save off the original mem pointer so all the pointers can be adjusted
           * to the new block of memory.
           */
          base = dcl->ad->mem;

          n = (socklen_t)dcl->ad->memsize;
          total += n;
          if ((dcl->ad->mem = malloc(n)) == NULL) LOG_FATAL_MSG("malloc", errno);
          if ((size = readn(worker->sock_fd, dcl->ad->mem, n)) == -1) {
            syslog(LOG_ERR,"[%s:%d] - reading %s error %d - %s\n", __FILE__, __LINE__, worker->ip_addr, errno, strerror(errno));
            break;
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
    }

    esl_stopwatch_Stop(w);

    if ((n = pthread_mutex_lock (&data->work_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

    /* set the state of the worker to completed */
    worker->cmd       = NULL;
    worker->completed = 1;
    worker->total     = total;
    --data->completed;

    /* notify all the worker threads of the new query */
    if ((n = pthread_cond_broadcast(&data->complete_cond)) != 0) LOG_FATAL_MSG("cond broadcast", n);

    if ((n = pthread_mutex_unlock (&data->work_mutex)) != 0) LOG_FATAL_MSG("mutex unlock", n);

    printf ("WORKER %s COMPLETED: %.2f sec received %d bytes\n", worker->ip_addr, w->elapsed, total), fflush(stdout);
  }

  esl_stopwatch_Destroy(w);

  return;
}

static void *
workerside_thread(void *arg)
{
  int               n;
  int               fd;

  HMMD_INIT_CMD     cmd;

  WORKER_DATA      *worker  = (WORKER_DATA *)arg;
  WORKERSIDE_ARGS  *parent  = (WORKERSIDE_ARGS *)worker->parent;

  /* Guarantees that thread resources are deallocated upon return */
  pthread_detach(pthread_self()); 

  printf("Handling worker %s (%d)\n", worker->ip_addr, worker->sock_fd);

  memset(&cmd, 0, sizeof(cmd));

  n = sizeof(cmd);

  cmd.length  = n;
  cmd.command = HMMD_CMD_INIT;
  cmd.db_cnt  = parent->seq_db->db_cnt;
  cmd.seq_cnt = parent->seq_db->count;

  strncpy(cmd.sid, parent->seq_db->id, sizeof(cmd.sid));
  cmd.sid[n-1]  = 0;

  if (writen(worker->sock_fd, &cmd, n) != n) {
    syslog(LOG_ERR,"[%s:%d] - writing (%d) error %d - %s\n", __FILE__, __LINE__, fd, errno, strerror(errno));
    printf("Closing worker %s (%d)\n", worker->ip_addr, fd);

    close(fd);

    pthread_exit(NULL);
  }

  worker->next = NULL;
  worker->prev = NULL;

  /* add the worker to the new fd list */
  if ((n = pthread_mutex_lock (&parent->work_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

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
    if ((n = pthread_cond_broadcast (&parent->worker_cond)) != 0) LOG_FATAL_MSG("cond broadcast", n);
  }

  if ((n = pthread_mutex_unlock (&parent->work_mutex)) != 0)  LOG_FATAL_MSG("mutex unlock", n);

  workerside_loop(parent, worker);

  if ((n = pthread_mutex_lock (&parent->work_mutex)) != 0) LOG_FATAL_MSG("mutex lock", n);

  fd = worker->sock_fd;

  worker->sock_fd    = -1;
  worker->terminated = 1;
  --parent->completed;
  ++parent->errors;

  if ((n = pthread_mutex_unlock (&parent->work_mutex)) != 0) LOG_FATAL_MSG("mutex unlock", n);

  printf("Closing worker %s (%d)\n", worker->ip_addr, fd);

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
    if ((fd = accept(data->sock_fd, (struct sockaddr *)&addr, (unsigned int *)&n)) < 0) LOG_FATAL_MSG("accept", errno);

    if ((worker = malloc(sizeof(WORKER_DATA))) == NULL) LOG_FATAL_MSG("thread create", errno);
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

    strncpy(worker->ip_addr, inet_ntoa(addr.sin_addr), INET_ADDRSTRLEN);
    worker->ip_addr[INET_ADDRSTRLEN-1] = 0;

    if ((n = pthread_create(&thread_id, NULL, workerside_thread, worker)) != 0) LOG_FATAL_MSG("thread create", n);
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
  if ((sock_fd = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) LOG_FATAL_MSG("socket", errno);
      
  /* incase the server went down in an ungraceful way, allow the port to be
   * reused avoiding the timeout.
   */
  reuse = 1;
  if (setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, (char *)&reuse, sizeof(reuse)) < 0) LOG_FATAL_MSG("setsockopt", errno);

  /* the sockets are never closed, so if the server exits, force the kernel to
   * close the socket and clear it so the server can be restarted immediately.
   */
  linger.l_onoff = 1;
  linger.l_linger = 0;
  if (setsockopt(sock_fd, SOL_SOCKET, SO_REUSEADDR, (char *)&linger, sizeof(linger)) < 0) LOG_FATAL_MSG("setsockopt", errno);

  /* Construct local address structure */
  memset(&addr, 0, sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_ANY);
  addr.sin_port = htons(esl_opt_GetInteger(opts, "--wport"));

  /* Bind to the local address */
  if (bind(sock_fd, (struct sockaddr *) &addr, sizeof(addr)) < 0) LOG_FATAL_MSG("bind", errno);

  /* Mark the socket so it will listen for incoming connections */
  if (listen(sock_fd, esl_opt_GetInteger(opts, "--wcncts")) < 0) LOG_FATAL_MSG("listen", errno);

  args->sock_fd = sock_fd;

  if ((n = pthread_create(&thread_id, NULL, worker_comm_thread, (void *)args)) != 0) LOG_FATAL_MSG("thread create", n);
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/

