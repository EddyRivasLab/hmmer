#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <sys/socket.h>
#ifdef HAVE_NETINET_IN_H
#include <netinet/in.h>     /* On FreeBSD, you need netinet/in.h for struct sockaddr_in            */
#endif                      /* On OpenBSD, netinet/in.h is required for (must precede) arpa/inet.h */
#include <arpa/inet.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <net/if.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"
#include "errno.h"

#include "hmmer.h"
#include "hmmpgmd.h"

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define MAX_READ_LEN     4096

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                         1 },
  /* Interface with server */
  { "-s",           eslARG_STRING,   "localhost", NULL, NULL,    NULL,  NULL,  NULL,            "name of the machine running hmmserver (default localhost)",                         8 },
  { "-p",           eslARG_INT,   "51371", NULL, NULL,    NULL,  NULL,  NULL,            "number of the port that hmmserver is listening on (default 51371)",                         8 },
  { "--db",           eslARG_STRING,   "1", NULL, NULL,    NULL,  NULL,  NULL,            "number of the database to search (default 1)",                         8 },
      /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                        2 },
  { "-A",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save multiple alignment of all hits to file <f>",              2 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-sequence hits to file <f>",        2 },
  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-domain hits to file <f>",          2 },
  { "--pfamtblout", eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save table of hits and domains to file, in Pfam format <f>",   2 },
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                       2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                         2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                     2 },
/* Control of scoring system */
  { "--popen",      eslARG_REAL,       "0.02", NULL, "0<=x<0.5",NULL,  NULL,  NULL,              "gap open probability",                                         3 },
  { "--pextend",    eslARG_REAL,        "0.4", NULL, "0<=x<1",  NULL,  NULL,  NULL,              "gap extend probability",                                       3 },
  { "--mx",         eslARG_STRING, "BLOSUM62", NULL, NULL,      NULL,  NULL,  "--mxfile",        "substitution score matrix choice (of some built-in matrices)", 3 },
  { "--mxfile",     eslARG_INFILE,       NULL, NULL, NULL,      NULL,  NULL,  "--mx",            "read substitution score matrix from file <f>",                 3 },
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
  { "--EmL",        eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,              "length of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EmN",        eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,              "number of sequences for MSV Gumbel mu fit",                   11 },   
  { "--EvL",        eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,              "length of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EvN",        eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,              "number of sequences for Viterbi Gumbel mu fit",               11 },   
  { "--EfL",        eslARG_INT,         "100", NULL,"n>0",      NULL,  NULL,  NULL,              "length of sequences for Forward exp tail tau fit",            11 },   
  { "--EfN",        eslARG_INT,         "200", NULL,"n>0",      NULL,  NULL,  NULL,              "number of sequences for Forward exp tail tau fit",            11 },   
  { "--Eft",        eslARG_REAL,       "0.04", NULL,"0<x<1",    NULL,  NULL,  NULL,              "tail mass for Forward exponential tail tau fit",              11 },    
/* Other options */
  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,    NULL,  NULL,  NULL,            "turn off biased composition score corrections",               12 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of comparisons done, for E-value calculation",          12 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of significant seqs, for domain E-value calculation",   12 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--qformat",    eslARG_STRING,      NULL, NULL, NULL,      NULL,  NULL,  NULL,              "assert query <seqfile> is in format <s>: no autodetection",   12 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = " <queryfile> <database number to search>";
static char banner[] = "search profile(s) against a sequence database";

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

  char             *firstseq_key;     /* name of the first sequence in the restricted db range */
  int              n_targetseq;       /* number of sequences in the restricted range */
};



static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_query, char **ret_db)
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

      if (puts("\nOptions directing output:")                                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      if (puts("\nOptions controlling reporting thresholds:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      if (puts("\nOptions controlling inclusion (significance) thresholds:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      if (puts("\nOptions controlling model-specific thresholding:")         < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 

      if (puts("\nOptions controlling acceleration heuristics:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      if (puts("\nOther expert options:")                                    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                  != 1)     { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_query= esl_opt_GetArg(go, 1)) == NULL)  { if (puts("Failed to get <queryfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  *ret_db = esl_opt_GetString (go, "--db");
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
output_header(FILE *ofp, const ESL_GETOPTS *go, char *hmmfile, char *seqfile)
{
  p7_banner(ofp, go->argv[0], banner);
  
  if (fprintf(ofp, "# query HMM file:                  %s\n", hmmfile)                                                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# target sequence database:        %s\n", seqfile)                                                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-o")           && fprintf(ofp, "# output directed to file:         %s\n",             esl_opt_GetString(go, "-o"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-A")           && fprintf(ofp, "# MSA of all hits saved to file:   %s\n",             esl_opt_GetString(go, "-A"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--tblout")     && fprintf(ofp, "# per-seq hits tabular output:     %s\n",             esl_opt_GetString(go, "--tblout"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domtblout")  && fprintf(ofp, "# per-dom hits tabular output:     %s\n",             esl_opt_GetString(go, "--domtblout"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--pfamtblout") && fprintf(ofp, "# pfam-style tabular hit output:   %s\n",             esl_opt_GetString(go, "--pfamtblout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--acc")        && fprintf(ofp, "# prefer accessions over names:    yes\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--noali")      && fprintf(ofp, "# show alignments in output:       no\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notextw")    && fprintf(ofp, "# max ASCII text line length:      unlimited\n")                                             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--textw")      && fprintf(ofp, "# max ASCII text line length:      %d\n",             esl_opt_GetInteger(go, "--textw"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-E")           && fprintf(ofp, "# sequence reporting threshold:    E-value <= %g\n",  esl_opt_GetReal(go, "-E"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-T")           && fprintf(ofp, "# sequence reporting threshold:    score >= %g\n",    esl_opt_GetReal(go, "-T"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domE")       && fprintf(ofp, "# domain reporting threshold:      E-value <= %g\n",  esl_opt_GetReal(go, "--domE"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domT")       && fprintf(ofp, "# domain reporting threshold:      score >= %g\n",    esl_opt_GetReal(go, "--domT"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incE")       && fprintf(ofp, "# sequence inclusion threshold:    E-value <= %g\n",  esl_opt_GetReal(go, "--incE"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incT")       && fprintf(ofp, "# sequence inclusion threshold:    score >= %g\n",    esl_opt_GetReal(go, "--incT"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incdomE")    && fprintf(ofp, "# domain inclusion threshold:      E-value <= %g\n",  esl_opt_GetReal(go, "--incdomE"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incdomT")    && fprintf(ofp, "# domain inclusion threshold:      score >= %g\n",    esl_opt_GetReal(go, "--incdomT"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_ga")     && fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--cut_nc")     && fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--cut_tc")     && fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--max")        && fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n")                        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F1")         && fprintf(ofp, "# MSV filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F1"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F2")         && fprintf(ofp, "# Vit filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F2"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F3")         && fprintf(ofp, "# Fwd filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F3"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nobias")     && fprintf(ofp, "# biased composition HMM filter:   off\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--restrictdb_stkey") && fprintf(ofp, "# Restrict db to start at seq key: %s\n",            esl_opt_GetString(go, "--restrictdb_stkey"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--restrictdb_n")     && fprintf(ofp, "# Restrict db to # target seqs:    %d\n",            esl_opt_GetInteger(go, "--restrictdb_n")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--ssifile")          && fprintf(ofp, "# Override ssi file to:            %s\n",            esl_opt_GetString(go, "--ssifile"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--nonull2")    && fprintf(ofp, "# null2 bias corrections:          off\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-Z")           && fprintf(ofp, "# sequence search space set to:    %.0f\n",           esl_opt_GetReal(go, "-Z"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domZ")       && fprintf(ofp, "# domain search space set to:      %.0f\n",           esl_opt_GetReal(go, "--domZ"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0 && fprintf(ofp, "# random number seed:              one-time arbitrary\n")                               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if (                               fprintf(ofp, "# random number seed set to:       %d\n",             esl_opt_GetInteger(go, "--seed"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--tformat")    && fprintf(ofp, "# targ <seqfile> format asserted:  %s\n",             esl_opt_GetString(go, "--tformat"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go       = NULL;	       
  int              status   = eslOK;
  char *buffer=NULL;
  char *cmd=NULL;
  uint8_t *buf;
  uint32_t buf_offset, hits_start;
  impl_Init();                  /* processor specific initialization */
  char *queryfile, *search_db;
  process_commandline(argc, argv, &go, &queryfile, &search_db);
  char *optsstring=NULL;
  optsstring =esl_getopts_CreateOptsLine(go); 
  ESL_ALLOC(buffer, MAX_READ_LEN);
  ESL_ALLOC(cmd, MAX_READ_LEN);
  int rem = MAX_READ_LEN;
  int cmdlen = MAX_READ_LEN;
  int eod = 0; 
  int i;
  ESL_STOPWATCH   *w       = NULL;
  P7_PIPELINE     *pli     = NULL;
  P7_TOPHITS      *th      = NULL;
  FILE *qf = fopen(queryfile, "r"); 
  int                  sock, n;
  unsigned short       serv_port;
  HMMD_SEARCH_STATS   *stats;
  HMMD_SEARCH_STATUS   sstatus;
  struct sockaddr_in   serv_addr;
  P7_HMM          *hmm     = NULL;        /* query HMM                       */
  ESL_SQ          *sq      = NULL;        /* one target sequence (digital)   */
  ESL_ALPHABET    *abc     = NULL;        /* digital alphabet                */
  //set up the sockets connection to the server

  // step 1: networking voodoo to get the IP address of the hostname for the server 
  struct addrinfo hints, *info;
  memset(&hints, 0, sizeof(hints));
  hints.ai_family = AF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_flags = AI_CANONNAME;
  getaddrinfo(esl_opt_GetString(go, "-s"), NULL, &hints, &info); //info->ai_addr now should have the sockaddr data structure we want
  char server_ip[256]; //should be longer than any possible ip address
  struct sockaddr_in *addr_temp;
  addr_temp = (struct sockaddr_in *) info->ai_addr;
  if (strlen(inet_ntoa(addr_temp->sin_addr)) > 255){
    p7_Fail("IP address of server %s appears to be more than 255 characters long, something has gone horribly wrong.\n",esl_opt_GetString(go, "-s"));
  }

  strcpy(server_ip, inet_ntoa(addr_temp->sin_addr));
  printf("Trying to connect to server at IP %s\n", server_ip);

  /* Create a reliable, stream socket using TCP */
  if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
    p7_Fail("[%s:%d] socket error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
    exit(1);
  }

  serv_port = esl_opt_GetInteger(go, "-p");
  /* Construct the server address structure */
  memset(&serv_addr, 0, sizeof(serv_addr));
  serv_addr.sin_family      = AF_INET;
  serv_addr.sin_addr.s_addr = inet_addr(server_ip);
  serv_addr.sin_port        = htons(serv_port);
  if ((inet_pton(AF_INET, server_ip, &serv_addr.sin_addr)) < 0) {
    p7_Fail("[%s:%d] inet_pton error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
  }

  /* Establish the connection to the server */
  if (connect(sock, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
    p7_Fail("[%s:%d] connect error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
  }

  if(qf == NULL){
    p7_Fail("Unable to open query file %s\n", queryfile);
  }
  // Start building the command string to send to the server
  // Can skip the length checks on cmd, because we know it has enough space for the opening command chars
  if(esl_opt_IsDefault(go, "--db")){//need to construct default database specifier
    strcpy(cmd, "@--db 1 ");
    rem -= 8;
  }
  else{//just need the starting @ symbol
    strcpy(cmd, "@");
    rem -=1;
    strcat(cmd, search_db);
    rem -= strlen(search_db);
    strcat(cmd, " ");
    rem -= 1;
  }
  

  while(strlen(optsstring)+1 >rem){
    cmd = realloc(cmd, 2*cmdlen); 
    if(cmd ==NULL){
      p7_Fail("[%s:%d] realloc error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
    }
    rem += cmdlen;
    cmdlen *=2;
  }
  
  strcat(cmd, optsstring);
  strcat(cmd, "\n");

  while (!eod) {
    if (fgets(buffer, MAX_READ_LEN, qf) != NULL) {
        n = strlen(buffer);
        if (n >= rem) {
          if ((cmd = realloc(cmd, cmdlen * 2)) == NULL) {
            p7_Fail("[%s:%d] realloc error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
          }
          rem += cmdlen;
          cmdlen *= 2;
        }
        strcat(cmd, buffer);
        rem -= n;
        cmdlen += n;
        eod = (strncmp(buffer, "//", 2) == 0);
    } 
    else {
        eod = 2;
    }
  }
  if(eod ==2){
    // We reached the end of the query object without finding a terminating '//' pair.  This is expected when the query object is a FASTA-format sequence
    if(rem < 3){
      if ((cmd = realloc(cmd, cmdlen +3)) == NULL) {
        p7_Fail("[%s:%d] realloc error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
      }
    }
    rem +=3;
    cmdlen +=3; 
    strcat(cmd, "//\n");
  }

  printf("About to send this command string:\n%s",cmd);
  n = strlen(cmd);
  printf ("Sending data %d\n", n);
  if (writen(sock, cmd, n) != n) {
    p7_Fail("[%s:%d] write (size %" PRIu64 ") error %d - %s\n", __FILE__, __LINE__, n, errno, strerror(errno));
    exit(1);
  }

  // Get the status structure back from the server
  buf = malloc(HMMD_SEARCH_STATUS_SERIAL_SIZE);
  buf_offset = 0;
  n = HMMD_SEARCH_STATUS_SERIAL_SIZE;
  int size;
  if(buf == NULL){
    p7_Fail("Unable to allocate memory for search status structure\n");
  }

  if ((size = readn(sock, buf, n)) == -1) {
    p7_Fail("[%s:%d] read error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
  }

  if(hmmd_search_status_Deserialize(buf, &buf_offset, &sstatus) != eslOK){
    p7_Fail("Unable to deserialize search status object \n");
  }

  if (sstatus.status != eslOK) {
    char *ebuf;
    n = sstatus.msg_size;
    ebuf = malloc(n);
    if ((size = readn(sock, ebuf, n)) == -1) {
      p7_Fail("[%s:%d] read error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
    }

    if (abc) esl_alphabet_Destroy(abc);
    if (hmm) p7_hmm_Destroy(hmm);
    if (sq)  esl_sq_Destroy(sq);
    free(ebuf);
    p7_Fail("ERROR (%d): %s\n", sstatus.status, ebuf);
  }

  free(buf); // clear this out 
  buf_offset = 0; // reset to beginning for next serialized object
  n = sstatus.msg_size;

  if ((buf = malloc(n)) == NULL) {
    p7_Fail("[%s:%d] malloc error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
  }
  // Grab the serialized search results
  if ((size = readn(sock, buf, n)) == -1) {
    p7_Fail("[%s:%d] read error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
  }

  if ((stats = malloc(sizeof(HMMD_SEARCH_STATS))) == NULL) {
    p7_Fail("[%s:%d] malloc error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
  }
  stats->hit_offsets = NULL; // force allocation of memory for this in _Deserialize
  if(p7_hmmd_search_stats_Deserialize(buf, &buf_offset, stats) != eslOK){
    p7_Fail("Unable to deserialize search stats object \n");
  }


  if(sstatus.type == HMMD_CMD_SEARCH){
    // Create the structures we'll deserialize the hits into
    pli = p7_pipeline_Create(go, 100, 100, FALSE, p7_SEARCH_SEQS);
    printf("Search status received\n");
  }
  else if(sstatus.type == HMMD_CMD_SCAN){
    // Create the structures we'll deserialize the hits into
    pli = p7_pipeline_Create(go, 100, 100, FALSE, p7_SCAN_MODELS);
    printf("Scan status received\n");
  }
  else p7_Fail("Illegal search type of %d found in HMMD_SEARCH_STATUS\n", sstatus.type);

  if(pli == NULL){
    p7_Fail("Unable to create pipeline data structure in hmmclient\n");
  }
        /* copy the search stats */
  w->elapsed       = stats->elapsed;
  w->user          = stats->user;
  w->sys           = stats->sys;

  pli->nmodels     = stats->nmodels;
  pli->nnodes      = stats->nnodes;
  pli->nseqs       = stats->nseqs;  
  pli->nres        = stats->nres;
  pli->n_past_msv  = stats->n_past_msv;
  pli->n_past_bias = stats->n_past_bias;
  pli->n_past_vit  = stats->n_past_vit;
  pli->n_past_fwd  = stats->n_past_fwd;

  pli->Z           = stats->Z;
  pli->domZ        = stats->domZ;
  pli->Z_setby     = stats->Z_setby;
  pli->domZ_setby  = stats->domZ_setby;

  th = p7_tophits_Create(); 

  free(th->unsrt); // free these because p7_tophits_Create() allocates a default amount of memory for them, and we're going to allocate the exact right amount next
  free(th->hit);

  th->N         = stats->nhits;
  if ((th->unsrt = malloc(stats-> nhits *sizeof(P7_HIT))) == NULL) {
    p7_Fail("[%s:%d] malloc error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
  }
  th->nreported = stats->nreported;
  th->nincluded = stats->nincluded;
  th->is_sorted_by_seqidx  = FALSE;
  th->is_sorted_by_sortkey = TRUE;

  if ((th->hit = malloc(sizeof(void *) * stats->nhits)) == NULL) {
    p7_Fail("[%s:%d] malloc error %d - %s\n", __FILE__, __LINE__, errno, strerror(errno));
  }
  hits_start = buf_offset;
  // deserialize the hits
  for (i = 0; i < stats->nhits; ++i) {
    // set all internal pointers of the hit to NULL before deserializing into it
    th->unsrt[i].name = NULL;
    th->unsrt[i].acc = NULL;
    th->unsrt[i].desc = NULL;
    th->unsrt[i].dcl = NULL;
    if((buf_offset -hits_start) != stats->hit_offsets[i]){
      printf("Hit offset %d did not match expected.  Found %d, expected %" PRIu64 "\n", i, (buf_offset-hits_start), stats->hit_offsets[i]);
    }
    if(p7_hit_Deserialize(buf, &buf_offset, &(th->unsrt[i])) != eslOK){
            p7_Fail("Unable to deserialize hit %d\n", i);
    }
    th->hit[i] = &(th->unsrt[i]);  
  }




  esl_getopts_Destroy(go);
  if(buffer) free(buffer);
  if(cmd) free(cmd);
  if(optsstring) free(optsstring);
  if (abc) esl_alphabet_Destroy(abc);
  if (hmm) p7_hmm_Destroy(hmm);
  if (sq)  esl_sq_Destroy(sq);

  return status;
ERROR:  
    p7_Fail("Unable to allocate memory in hmmclient\n");
}
