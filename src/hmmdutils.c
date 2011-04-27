/* hmmpgmd: hmmer deamon searchs against a sequence database.
 * 
 * MSF, Thu Aug 12, 2010 [Janelia]
 * SVN $Id$
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
#include <assert.h>

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

#define MAX_WORKERS  64
#define MAX_BUFFER   4096

#define CONF_FILE "/etc/hmmpgmd.conf"

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define STAGEOPTS   "--F1,--F2,--F3"

static ESL_OPTIONS searchOpts[] = {
  /* Control of output */
  { "--acc",        eslARG_NONE,      FALSE, NULL, NULL,      NULL,  NULL, NULL,        "prefer accessions over names in output",                       2 },
  { "--noali",      eslARG_NONE,      FALSE, NULL, NULL,      NULL,  NULL, NULL,        "don't output alignments, so output is smaller",                2 },
  /* Control of scoring system */
  { "--popen",      eslARG_REAL,     "0.02", NULL, "0<=x<0.5",NULL,  NULL, NULL,        "gap open probability",                                         3 },
  { "--pextend",    eslARG_REAL,      "0.4", NULL, "0<=x<1",  NULL,  NULL, NULL,        "gap extend probability",                                       3 },
  { "--mxfile",   eslARG_STRING, "BLOSUM62", NULL, NULL,      NULL,  NULL, NULL,        "substitution score matrix [default: BLOSUM62]",                3 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,     "10.0", NULL, "x>0",     NULL,  NULL, REPOPTS,     "report sequences <= this E-value threshold in output",         4 },
  { "-T",           eslARG_REAL,      FALSE, NULL, NULL,      NULL,  NULL, REPOPTS,     "report sequences >= this score threshold in output",           4 },
  { "--domE",       eslARG_REAL,     "10.0", NULL, "x>0",     NULL,  NULL, DOMREPOPTS,  "report domains <= this E-value threshold in output",           4 },
  { "--domT",       eslARG_REAL,      FALSE, NULL, NULL,      NULL,  NULL, DOMREPOPTS,  "report domains >= this score cutoff in output",                4 },
  /* Control of inclusion (significance) thresholds */
  { "--incE",       eslARG_REAL,     "0.01", NULL, "x>0",     NULL,  NULL, INCOPTS,     "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",       eslARG_REAL,      FALSE, NULL, NULL,      NULL,  NULL, INCOPTS,     "consider sequences >= this score threshold as significant",    5 },
  { "--incdomE",    eslARG_REAL,     "0.01", NULL, "x>0",     NULL,  NULL, INCDOMOPTS,  "consider domains <= this E-value threshold as significant",    5 },
  { "--incdomT",    eslARG_REAL,      FALSE, NULL, NULL,      NULL,  NULL, INCDOMOPTS,  "consider domains >= this score threshold as significant",      5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,      FALSE, NULL, NULL,      NULL,  NULL, THRESHOPTS,  "use profile's GA gathering cutoffs to set all thresholding",   6 },
  { "--cut_nc",     eslARG_NONE,      FALSE, NULL, NULL,      NULL,  NULL, THRESHOPTS,  "use profile's NC noise cutoffs to set all thresholding",       6 },
  { "--cut_tc",     eslARG_NONE,      FALSE, NULL, NULL,      NULL,  NULL, THRESHOPTS,  "use profile's TC trusted cutoffs to set all thresholding",     6 },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,      FALSE, NULL, NULL,      NULL,  NULL, STAGEOPTS,   "Turn all heuristic filters off (less speed, more power)",      7 },
  { "--F1",         eslARG_REAL,     "0.02", NULL, NULL,      NULL,  NULL, "--max",     "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F2",         eslARG_REAL,     "1e-3", NULL, NULL,      NULL,  NULL, "--max",     "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F3",         eslARG_REAL,     "1e-5", NULL, NULL,      NULL,  NULL, "--max",     "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--nobias",     eslARG_NONE,       NULL, NULL, NULL,      NULL,  NULL, "--max",     "turn off composition bias filter",                             7 },
  /* Control of E-value calibration */
  { "--EmL",        eslARG_INT,       "200", NULL,"n>0",      NULL,  NULL, NULL,        "length of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EmN",        eslARG_INT,       "200", NULL,"n>0",      NULL,  NULL, NULL,        "number of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EvL",        eslARG_INT,       "200", NULL,"n>0",      NULL,  NULL, NULL,        "length of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EvN",        eslARG_INT,       "200", NULL,"n>0",      NULL,  NULL, NULL,        "number of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EfL",        eslARG_INT,       "100", NULL,"n>0",      NULL,  NULL, NULL,        "length of sequences for Forward exp tail tau fit",            11 },   
  { "--EfN",        eslARG_INT,       "200", NULL,"n>0",      NULL,  NULL, NULL,        "number of sequences for Forward exp tail tau fit",            11 },   
  { "--Eft",        eslARG_REAL,     "0.04", NULL,"0<x<1",    NULL,  NULL, NULL,        "tail mass for Forward exponential tail tau fit",              11 },   
  /* Other options */
  { "--seed",       eslARG_INT,        "42", NULL, "n>=0",    NULL,  NULL, NULL,        "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--nonull2",    eslARG_NONE,       NULL, NULL, NULL,      NULL,  NULL, NULL,        "turn off biased composition score corrections",               12 },
  { "-Z",           eslARG_REAL,      FALSE, NULL, "x>0",     NULL,  NULL, NULL,        "set # of comparisons done, for E-value calculation",          12 },
  { "--domZ",       eslARG_REAL,      FALSE, NULL, "x>0",     NULL,  NULL, NULL,        "set # of significant seqs, for domain E-value calculation",   12 },
  { "--seqdb",      eslARG_INT,        NULL, NULL, "n>0",     NULL,  NULL, "--hmmdb",   "protein database to search",                                  12 },
  { "--hmmdb",      eslARG_INT,        NULL, NULL, "n>0",     NULL,  NULL, "--seqdb",   "hmm database to search",                                      12 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

size_t
writen(int fd, const void *vptr, size_t n)
{
  ssize_t     remaining;
  ssize_t     outn;
  const char *ptr;

  ptr = vptr;
  remaining = n;
  while (remaining > 0) {
    if ((outn = write(fd, ptr, remaining)) <= 0) {
      if (outn < 0 && errno == EINTR) {
        outn = 0;
      } else {
        return -1;
      }
    }

    remaining -= outn;
    ptr += outn;
  }

  return n;
}

size_t
readn(int fd, void *vptr, size_t n)
{
  size_t      remaining;
  size_t      bytes;
  char       *ptr;

  ptr = vptr;
  remaining = n;
  while (remaining > 0) {
    if ((bytes = read(fd, ptr, remaining)) <= 0) {
      if (errno == EINTR) {
        bytes = 0;
      } else {
        return -1;
      }
    }
    
    remaining -= bytes;
    ptr += bytes;
  }

  return n - remaining;
}

#define LOG_TO_STDOUT
#ifdef LOG_TO_STDOUT
void p7_openlog(const char *ident, int option, int facility)
{
  /* do nothing */
  return;
}

void p7_syslog(int priority, const char *format, ...)
{
  va_list ap;

  printf("\n*** ERROR ***\n");

  va_start(ap, format);
  vprintf(format, ap);
  va_end(ap);

  printf("\n");
  fflush(stdout);

  return;
}
void p7_closelog(void)
{
  /* do nothing */
  return;
}
#endif

int
process_searchopts(int fd, char *cmdstr, ESL_GETOPTS **ret_opts)
{
  int status;

  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(searchOpts))       == NULL)  return eslEMEM;
  if ((status = esl_opt_ProcessSpoof(go, cmdstr)) != eslOK) return status;
  if ((status = esl_opt_VerifyConfig(go))         != eslOK) return status;

  *ret_opts = go;
  return eslOK;
}

void
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
