/* hmmscant: search sequence(s) against a profile HMM database
 *
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

/* for hmmscant */
#include "esl_gencode.h"

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  ESL_SQ                *ntqsq;  /* query or target sequence; this is a DNA sequence in the case of hmmscant */
  P7_BG                 *bg;	 /* null model                                                               */
  P7_PIPELINE           *pli;    /* work pipeline                                                            */
  P7_TOPHITS            *th;     /* top hit results                                                          */
  ESL_GENCODE           *gcode;  /* used for translating ORFs                                                */
  ESL_GENCODE_WORKSTATE *wrk;    /* maintain state of nucleotide sequence in the midst of processing ORFs    */

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
  /* name           type            default  env          range    toggles reqs   incomp             help                                                          docgroup*/
  { "-h",           eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  NULL,             "show brief help on version and usage",                          1  },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL,    NULL,        NULL,     NULL,  NULL,  NULL,             "direct output to file <f>, not stdout",                         2  },
  { "--tblout",     eslARG_OUTFILE, NULL,    NULL,        NULL,     NULL,  NULL,  NULL,             "save parseable table of per-sequence hits to file <f>",         2  },
  { "--domtblout",  eslARG_OUTFILE, NULL,    NULL,        NULL,     NULL,  NULL,  NULL,             "save parseable table of per-domain hits to file <f>",           2  },
  { "--acc",        eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  NULL,             "prefer accessions over names in output",                        2  },
  { "--noali",      eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  NULL,             "don't output alignments, so output is smaller",                 2  },
  { "--notextw",    eslARG_NONE,    NULL,    NULL,        NULL,     NULL,  NULL,  "--textw",        "unlimit ASCII text output line width",                          2  },
  { "--textw",      eslARG_INT,     "120",   NULL,        "n>=120", NULL,  NULL,  "--notextw",      "set max width of ASCII text output lines",                      2  },
  { "--notrans",    eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  NULL,             "don't show the translated DNA sequence in domain alignment",    2  }, /*for hmmscant */
  { "--vertcodon",  eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  NULL,             "show the DNA vertically in domain alignment",                   2  }, /*for hmmscant */
  /* Translation Options */
  { "-c",           eslARG_INT,     "1",     NULL,        NULL,     NULL,  NULL,  NULL,             "use alt genetic code of NCBI transl table <n>",                 15 },
  { "-l",           eslARG_INT,     "20",    NULL,        NULL,     NULL,  NULL,  NULL,             "minimum ORF length",                                            15 },
  { "-m",           eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  "-M",             "ORFs must initiate with AUG only",                              15 },
  { "-M",           eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  "-m",             "ORFs must start with allowed initiation codon",                 15 },
  { "--watson",     eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  NULL,             "only translate top strand",                                     15 },
  { "--crick",      eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  NULL,             "only translate bottom strand",                                  15 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,    "10.0",  NULL,        "x>0",    NULL,  NULL,  REPOPTS,          "report models <= this E-value threshold in output",             4  },
  { "-T",           eslARG_REAL,    FALSE,   NULL,        NULL,     NULL,  NULL,  REPOPTS,          "report models >= this score threshold in output",               4  },
  { "--domE",       eslARG_REAL,    "10.0",  NULL,        "x>0",    NULL,  NULL,  DOMREPOPTS,       "report domains <= this E-value threshold in output",            4  },
  { "--domT",       eslARG_REAL,    FALSE,   NULL,        NULL,     NULL,  NULL,  DOMREPOPTS,       "report domains >= this score cutoff in output",                 4  },
  /* Control of inclusion (significance) thresholds: */
  { "--incE",       eslARG_REAL,    "0.01",  NULL,        "x>0",    NULL,  NULL,  INCOPTS,          "consider models <= this E-value threshold as significant",      5  },
  { "--incT",       eslARG_REAL,    FALSE,   NULL,        NULL,     NULL,  NULL,  INCOPTS,          "consider models >= this score threshold as significant",        5  },
  { "--incdomE",    eslARG_REAL,    "0.01",  NULL,        "x>0",    NULL,  NULL,  INCDOMOPTS,       "consider domains <= this E-value threshold as significant",     5  },
  { "--incdomT",    eslARG_REAL,    FALSE,   NULL,        NULL,     NULL,  NULL,  INCDOMOPTS,       "consider domains >= this score threshold as significant",       5  },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  THRESHOPTS,       "use profile's GA gathering cutoffs to set all thresholding",    6  },
  { "--cut_nc",     eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  THRESHOPTS,       "use profile's NC noise cutoffs to set all thresholding",        6  },
  { "--cut_tc",     eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  THRESHOPTS,       "use profile's TC trusted cutoffs to set all thresholding",      6  },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,    FALSE,   NULL,        NULL,     NULL,  NULL,  "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",       7  },
  { "--F1",         eslARG_REAL,    "0.02",  NULL,        NULL,     NULL,  NULL,  "--max",          "MSV threshold: promote hits w/ P <= F1",                        7  },
  { "--F2",         eslARG_REAL,    "1e-3",  NULL,        NULL,     NULL,  NULL,  "--max",          "Vit threshold: promote hits w/ P <= F2",                        7  },
  { "--F3",         eslARG_REAL,    "1e-5",  NULL,        NULL,     NULL,  NULL,  "--max",          "Fwd threshold: promote hits w/ P <= F3",                        7  },
  { "--nobias",     eslARG_NONE,    NULL,    NULL,        NULL,     NULL,  NULL,  "--max",          "turn off composition bias filter",                              7  },
  /* Other options */
  { "--nonull2",    eslARG_NONE,    NULL,    NULL,        NULL,     NULL,  NULL,  NULL,             "turn off biased composition score corrections",                 12 },
  { "-Z",           eslARG_REAL,    FALSE,   NULL,        "x>0",    NULL,  NULL,  NULL,             "set # of comparisons done, for E-value calculation",            12 },
  { "--domZ",       eslARG_REAL,    FALSE,   NULL,        "x>0",    NULL,  NULL,  NULL,             "set # of significant seqs, for domain E-value calculation",     12 },
  { "--seed",       eslARG_INT,     "42",    NULL,        "n>=0",   NULL,  NULL,  NULL,             "set RNG seed to <n> (if 0: one-time arbitrary seed)",           12 },
  { "--qformat",    eslARG_STRING,  NULL,    NULL,        NULL,     NULL,  NULL,  NULL,             "assert input <seqfile> is in format <s>: no autodetection",     12 },
#ifdef HMMER_THREADS
  { "--cpu",        eslARG_INT,     p7_NCPU, "HMMER_NCPU","n>=0",   NULL,  NULL,  CPUOPTS,          "number of parallel CPU workers to use for multithreads",        12 },

#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads).
 */
struct cfg_s {
  char            *seqfile;           /* query sequence file                             */
  char            *hmmfile;           /* database HMM file                               */
};

static char usage[]  = "[-options] <hmmdb> <seqfile>";
static char banner[] = "search DNA sequence(s) against a protein profile database";

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, P7_HMMFILE *hfp);
#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ID_LENGTH_LIST *id_length_list, int64_t seqidx, int64_t seqL, P7_HMMFILE *hfp);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/


/* process_commandline()
 * 
 * Processes the commandline, filling in fields in <cfg> and creating and returning
 * an <ESL_GETOPTS> options structure. The help page (hmmscant -h) is formatted
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

      if (puts("\nTranslation options:")                                     < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 15, 2, 80); 

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

      if (puts("\nAvailable NCBI genetic code tables (for -c <id>):")        < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_gencode_DumpAltCodeTable(stdout);

      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)      { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { if (puts("Failed to get <hmmdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { if (puts("Failed to get <seqfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_hmmfile, "-") == 0) 
    { if (puts("hmmscan cannot read <hmm database> from stdin stream, because it must have hmmpress'ed auxfiles") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");   goto FAILURE;  }

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
  if (esl_opt_IsUsed(go, "--domtblout") && fprintf(ofp, "# per-dom hits tabular output:     %s\n",            esl_opt_GetString(go, "--domtblout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--acc")       && fprintf(ofp, "# prefer accessions over names:    yes\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--noali")     && fprintf(ofp, "# show alignments in output:       no\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notextw")   && fprintf(ofp, "# max ASCII text line length:      unlimited\n")                                           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--textw")     && fprintf(ofp, "# max ASCII text line length:      %d\n",            esl_opt_GetInteger(go, "--textw"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");  
  if (esl_opt_IsUsed(go, "--notrans")   && fprintf(ofp, "# show translated DNA sequence:    no\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--vertcodon") && fprintf(ofp, "# show DNA codon vertically:       no\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-E")          && fprintf(ofp, "# profile reporting threshold:     E-value <= %g\n", esl_opt_GetReal(go, "-E"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-T")          && fprintf(ofp, "# profile reporting threshold:     score >= %g\n",   esl_opt_GetReal(go, "-T"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domE")      && fprintf(ofp, "# domain reporting threshold:      E-value <= %g\n", esl_opt_GetReal(go, "--domE"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domT")      && fprintf(ofp, "# domain reporting threshold:      score >= %g\n",   esl_opt_GetReal(go, "--domT"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incE")      && fprintf(ofp, "# profile inclusion threshold:     E-value <= %g\n", esl_opt_GetReal(go, "--incE"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incT")      && fprintf(ofp, "# profile inclusion threshold:     score >= %g\n",   esl_opt_GetReal(go, "--incT"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incdomE")   && fprintf(ofp, "# domain inclusion threshold:      E-value <= %g\n", esl_opt_GetReal(go, "--incdomE"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incdomT")   && fprintf(ofp, "# domain inclusion threshold:      score >= %g\n",   esl_opt_GetReal(go, "--incdomT"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_ga")    && fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_nc")    && fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_tc")    && fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--max")       && fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n")                      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F1")        && fprintf(ofp, "# MSV filter P threshold:       <= %g\n",            esl_opt_GetReal(go, "--F1"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F2")        && fprintf(ofp, "# Vit filter P threshold:       <= %g\n",            esl_opt_GetReal(go, "--F2"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F3")        && fprintf(ofp, "# Fwd filter P threshold:       <= %g\n",            esl_opt_GetReal(go, "--F3"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nobias")    && fprintf(ofp, "# biased composition HMM filter:   off\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nonull2")   && fprintf(ofp, "# null2 bias corrections:          off\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-Z")          && fprintf(ofp, "# sequence search space set to:    %.0f\n",          esl_opt_GetReal(go, "-Z"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domZ")      && fprintf(ofp, "# domain search space set to:      %.0f\n",          esl_opt_GetReal(go, "--domZ"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--seed")) 
  {
    if (esl_opt_GetInteger(go, "--seed")==0 && fprintf(ofp, "# random number seed:              one-time arbitrary\n")                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if (                                  fprintf(ofp, "# random number seed set to:       %d\n",        esl_opt_GetInteger(go, "--seed"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--qformat")   && fprintf(ofp, "# input seqfile format asserted:   %s\n",            esl_opt_GetString(go, "--qformat"))   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu")       && fprintf(ofp, "# number of worker threads:        %d\n",            esl_opt_GetInteger(go, "--cpu"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");  
#endif
  if (esl_opt_IsUsed(go, "-c")          && fprintf(ofp, "# use alt genetic code of NCBI transl table: %d\n", esl_opt_GetInteger(go, "-c"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-l")          && fprintf(ofp, "# minimum ORF length: %d\n",                        esl_opt_GetInteger(go, "-l"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-m")          && fprintf(ofp, "# ORFs must initiate with AUG only:    yes\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-M")          && fprintf(ofp, "# ORFs must start with allowed initiation codon:    yes\n")                               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--watson")    && fprintf(ofp, "# only translate top strand:    yes\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--crick")     && fprintf(ofp, "# only translate bottom strand:    yes\n")                                                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (                                     fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go  = NULL;	
  struct cfg_s     cfg;         
  int              status   = eslOK;

  impl_Init();			/* processor-specific initialization */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */

  /* Initialize what we can in the config structure (without knowing the alphabet yet) */
  cfg.hmmfile    = NULL;
  cfg.seqfile    = NULL;

  process_commandline(argc, argv, &go, &cfg.hmmfile, &cfg.seqfile);    


  status = serial_master(go, &cfg);

  esl_getopts_Destroy(go);
  return status;
}


/* translate_sequence()
 * For input DNA sequence, add all ORFs (6 frames) to wrk block
 */
static int
translate_sequence(ESL_GENCODE *gcode, ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq)
{
      if (wrk->do_watson) {
        esl_gencode_ProcessStart(gcode, wrk, sq);
        esl_gencode_ProcessPiece(gcode, wrk, sq);
        esl_gencode_ProcessEnd(wrk, sq);
      }
      if (wrk->do_crick) {
        esl_sq_ReverseComplement(sq);
        esl_gencode_ProcessStart(gcode, wrk, sq);
        esl_gencode_ProcessPiece(gcode, wrk, sq);
        esl_gencode_ProcessEnd(wrk, sq);
        esl_sq_ReverseComplement(sq);
      }

  return eslOK;
}


/* serial_master()
 * The serial version of hmmscant.
 * For each query HMM in <hmmdb> search the database for hits.
 * 
 * A master can only return if it's successful. All errors are handled
 * immediately and fatally with p7_Fail().
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;	         /* output file for results (default stdout)                */
  FILE            *tblfp    = NULL;		 /* output stream for tabular per-seq (--tblout)            */
  FILE            *domtblfp = NULL;	  	 /* output stream for tabular per-seq (--domtblout)         */
  int              seqfmt   = eslSQFILE_UNKNOWN; /* format of seqfile                                       */
  ESL_SQFILE      *sqfp     = NULL;              /* open seqfile                                            */
  P7_HMMFILE      *hfp      = NULL;		 /* open HMM database file                                  */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                                       */
  P7_OPROFILE     *om       = NULL;		 /* target profile                                          */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                                  */
  int              nquery   = 0;
  int              textw;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;
  uint64_t         prev_char_cnt = 0;
  
  /* used to keep track of the lengths of the sequences that are processed */
  ID_LENGTH_LIST  *id_length_list = NULL;

  int              ncpus    = 0;
  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  P7_OM_BLOCK     *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char             errbuf[eslERRBUFSIZE];

  /* hmmscant */
  P7_TOPHITS      *tophits_accumulator       = NULL; /* to hold the top hits information from all 6 frame translations     */
  P7_PIPELINE     *pipelinehits_accumulator  = NULL; /* to hold the pipeline hit information from all 6 frame translations */

  ESL_ALPHABET    *abcDNA                    = NULL; /* DNA sequence alphabet                                              */
  ESL_ALPHABET    *abcAMINO                  = NULL; /* DNA sequence alphabet                                              */

  ESL_SQ          *qsqDNA                    = NULL; /* DNA query sequence                                                 */

  ESL_GENCODE     *gcode                     = NULL;
  ESL_GENCODE_WORKSTATE *wrk                 = NULL;
  /* end hmmscant */

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* If caller declared an input format, decode it */
  if (esl_opt_IsOn(go, "--qformat")) {
    seqfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (seqfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }

  /*the query sequence will be DNA but will be translated to amino acids */
  abcDNA   = esl_alphabet_Create(eslDNA); 
  abcAMINO = esl_alphabet_Create(eslAMINO); 
  qsqDNA   = esl_sq_CreateDigital(abcDNA);


  /* Open the target profile database to get the sequence alphabet */
  status = p7_hmmfile_OpenE(cfg->hmmfile, p7_HMMDBENV, &hfp, errbuf);
  if      (status == eslENOTFOUND)  p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)    p7_Fail("File format problem, trying to open HMM file %s.\n%s\n",                  cfg->hmmfile, errbuf);
  else if (status != eslOK)         p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  
  if (! hfp->is_pressed)            p7_Fail("Failed to open binary auxfiles for %s: use hmmpress first\n",             hfp->fname);

  hstatus = p7_oprofile_ReadMSV(hfp, &abc, &om);
  if      (hstatus == eslEFORMAT)   p7_Fail("bad format, binary auxfiles, %s:\n%s",     cfg->hmmfile, hfp->errbuf);
  else if (hstatus == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets", cfg->hmmfile);
  else if (hstatus != eslOK)        p7_Fail("Unexpected error in reading HMMs from %s", cfg->hmmfile); 

  if (abc->type != eslAMINO)        p7_Fail("hmmscant only supports amino acid HMMs; %s uses a different alphabet", cfg->hmmfile);
  
  p7_oprofile_Destroy(om);
  p7_hmmfile_Close(hfp);

  /* Open the query sequence database */
  status = esl_sqfile_OpenDigital(abcDNA, cfg->seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND)  p7_Fail("Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)    p7_Fail("Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)     p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)         p7_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);


 
  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"),          "w")) == NULL)  esl_fatal("Failed to open output file %s for writing\n",                 esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }
  if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblout")); }

  output_header(ofp, go, cfg->hmmfile, cfg->seqfile);


  /* Set up the genetic code. Default = NCBI 1, the standard code; allow ORFs to start at any Amino Acid */
  gcode = esl_gencode_Create(abcDNA, abcAMINO);
  esl_gencode_Set(gcode, esl_opt_GetInteger(go, "-c"));  // default = 1, the standard genetic code

  if      (esl_opt_GetBoolean(go, "-m"))   esl_gencode_SetInitiatorOnlyAUG(gcode);
  else if (! esl_opt_GetBoolean(go, "-M")) esl_gencode_SetInitiatorAny(gcode);      // note this is the default, if neither -m or -M are set


#ifdef HMMER_THREADS
  /* initialize thread data */
  ncpus = ESL_MIN( esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());


  if (ncpus > 0)
    {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue     = esl_workqueue_Create(ncpus * 2);
    }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  for (i = 0; i < infocnt; ++i)
    {
      info[i].bg    = p7_bg_Create(abc);
#ifdef HMMER_THREADS
      info[i].queue = queue;
#endif
    }

#ifdef HMMER_THREADS
  for (i = 0; i < ncpus * 2; ++i)
    {
      block = p7_oprofile_CreateBlock(BLOCK_SIZE);
      if (block == NULL)    esl_fatal("Failed to allocate sequence block");

      status = esl_workqueue_Init(queue, block);
      if (status != eslOK)  esl_fatal("Failed to add block to work queue");
    }
#endif

  /* Outside loop: over each query sequence in <seqfile>. */
  while ((sstatus = esl_sqio_Read(sqfp, qsqDNA)) == eslOK)
  {
     if (qsqDNA->n < 3) continue; /* do not process sequence of less than 1 codon */

     nquery++;
     esl_stopwatch_Start(w);	                          


     /* Create processing pipeline and hit list accumulators */
     tophits_accumulator      = p7_tophits_Create(); 
     pipelinehits_accumulator = p7_pipeline_Create(go, 100, 100, FALSE, p7_SCAN_MODELS);
     
     pipelinehits_accumulator->nseqs         = 0;
     pipelinehits_accumulator->nres          = 0;
     pipelinehits_accumulator->is_translated = TRUE;

     if (fprintf(ofp, "Query:       %s  [L=%ld]\n", qsqDNA->name, (long) qsqDNA->n) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
     if (qsqDNA->acc[0]  != 0 && fprintf(ofp, "Accession:   %s\n", qsqDNA->acc)     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
     if (qsqDNA->desc[0] != 0 && fprintf(ofp, "Description: %s\n", qsqDNA->desc)    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");


     /* Open the target profile database */
     status = p7_hmmfile_OpenE(cfg->hmmfile, p7_HMMDBENV, &hfp, NULL);
     if (status != eslOK)        p7_Fail("Unexpected error %d in opening hmm file %s.\n",           status, cfg->hmmfile);
  
#ifdef HMMER_THREADS
     /* if we are threaded, create a lock to prevent multiple readers */
     if (ncpus > 0)
     {
         status = p7_hmmfile_CreateLock(hfp);
         if (status != eslOK) p7_Fail("Unexpected error %d creating lock\n", status);
     }
#endif


     for (i = 0; i < infocnt; ++i)
     {
       /* Create processing pipeline and hit list */

       info[i].gcode              = gcode;
       info[i].wrk                = esl_gencode_WorkstateCreate(go, gcode);
       info[i].wrk->orf_block     = esl_sq_CreateDigitalBlock(3, abcAMINO);
       info[i].th                 = p7_tophits_Create();
       info[i].pli                = p7_pipeline_Create(go, 100, 100, FALSE, p7_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
       info[i].pli->hfp           = hfp;  /* for two-stage input, pipeline needs <hfp> */
       info[i].pli->is_translated = TRUE;
       info[i].ntqsq              = qsqDNA;
       info[i].ntqsq->prev_n      = prev_char_cnt;

#ifdef HMMER_THREADS
       if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
     }

     /* establish the id_lengths data structutre */
      id_length_list = init_id_length(1000);

#ifdef HMMER_THREADS
     if (ncpus > 0)
         hstatus = thread_loop(threadObj, queue, id_length_list, qsqDNA->idx, qsqDNA->L, hfp);
     else
#endif
         hstatus = serial_loop(info, id_length_list, hfp);

     switch(hstatus)
     {
       case eslEFORMAT:   p7_Fail("bad file format in HMM file %s",             cfg->hmmfile);	  break;
       case eslEINCOMPAT: p7_Fail("HMM file %s contains different alphabets",   cfg->hmmfile);	  break;
       case eslEOF: 	  /* do nothing */                                                 	  break;
       default: 	  p7_Fail("Unexpected error in reading HMMs from %s",   cfg->hmmfile);
     }

     
     /* merge the results of the search results */
     for (i = 0; i < infocnt; ++i)
     {
       p7_tophits_Merge(tophits_accumulator, info[i].th);
       p7_pipeline_Merge(pipelinehits_accumulator, info[i].pli);
       
       pipelinehits_accumulator->nseqs += info[i].pli->nseqs;
       pipelinehits_accumulator->nres  += info[i].pli->nres;

       p7_pipeline_Destroy(info[i].pli);
       p7_tophits_Destroy(info[i].th);
     }

     p7_hmmfile_Close(hfp);
     
     /* Sort and remove duplicates */
     p7_tophits_SortBySeqidxAndAlipos(tophits_accumulator);
     assign_Lengths(tophits_accumulator, id_length_list);
     p7_tophits_RemoveDuplicates(tophits_accumulator, pipelinehits_accumulator->use_bit_cutoffs);

     p7_tophits_SortBySortkey(tophits_accumulator);
     p7_tophits_Threshold(tophits_accumulator, pipelinehits_accumulator);


     /* Print results */	 
     p7_tophits_Targets(ofp, tophits_accumulator, pipelinehits_accumulator, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
     p7_tophits_Domains(ofp, tophits_accumulator, pipelinehits_accumulator, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

     if (tblfp)     p7_tophits_TabularTargets(tblfp,    qsqDNA->name, qsqDNA->acc, tophits_accumulator, pipelinehits_accumulator, (nquery == 1));
     if (domtblfp)  p7_tophits_TabularDomains(domtblfp, qsqDNA->name, qsqDNA->acc, tophits_accumulator, pipelinehits_accumulator, (nquery == 1));

     esl_stopwatch_Stop(w);
     p7_pli_Statistics(ofp, pipelinehits_accumulator, w);
     if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
     fflush(ofp);

     p7_pipeline_Destroy(pipelinehits_accumulator);
     p7_tophits_Destroy(tophits_accumulator);
     destroy_id_length(id_length_list);     

     prev_char_cnt += qsqDNA->n;
     esl_sq_Reuse(qsqDNA);

     for (i = 0; i < infocnt; ++i)
       esl_sq_ReuseBlock(info[i].wrk->orf_block);
     

  } /*end hmmscant while loop */

  if      (sstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					    sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (sstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    sstatus, sqfp->filename);

  /* Terminate outputs - any last words? */
  if (tblfp)     p7_tophits_TabularTail(tblfp,    "hmmscant", p7_SCAN_MODELS, cfg->seqfile, cfg->hmmfile, go);
  if (domtblfp)  p7_tophits_TabularTail(domtblfp, "hmmscant", p7_SCAN_MODELS, cfg->seqfile, cfg->hmmfile, go);
  if (ofp)       { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }


  /* Cleanup - prepare for successful exit */
  for (i = 0; i < infocnt; ++i) {
    p7_bg_Destroy(info[i].bg);
    esl_gencode_WorkstateDestroy(info[i].wrk);
  }


#ifdef HMMER_THREADS
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

  esl_gencode_WorkstateDestroy(wrk);
  esl_gencode_Destroy(gcode);

  esl_sq_Destroy(qsqDNA); 
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_alphabet_Destroy(abcDNA); 
  esl_alphabet_Destroy(abcAMINO); 
  esl_sqfile_Close(sqfp);

  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);
  return eslOK;

 ERROR:
  return status;
}



static int
serial_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, P7_HMMFILE *hfp)
{
  int            status;

  int            k;
  P7_OPROFILE   *om;
  ESL_ALPHABET  *abc = NULL;
  ESL_SQ        *qsq_aa = NULL;                  /* query sequence, amino acid                                 */
  ESL_SQ        *qsqDNATxt = esl_sq_Create();    /* DNA query sequence that will be in text mode for printing */

  /* copy and convert the DNA sequence to text so we can print it in the domain alignment display */
  esl_sq_Copy(info->ntqsq, qsqDNATxt);
  qsqDNATxt->abc = info->ntqsq->abc;

  add_id_length(id_length_list, info->ntqsq->idx, info->ntqsq->L);
  /* translate DNA sequence to 6 frame ORFs */
  translate_sequence(info->gcode, info->wrk, info->ntqsq);
  for (k = 0; k < info->wrk->orf_block->count; ++k)
  {
      qsq_aa = &(info->wrk->orf_block->list[k]);
      p7_pli_NewSeq(info->pli, qsq_aa);
  }

  /* Main loop: */
  while ((status = p7_oprofile_ReadMSV(hfp, &abc, &om)) == eslOK)
  {
      p7_pli_NewModel(info->pli, om, info->bg);

      /*process each 6 frame translated sequence */
      for (k = 0; k < info->wrk->orf_block->count; ++k)
      {
          qsq_aa = &(info->wrk->orf_block->list[k]);

          /* use the name, accession, and description from the DNA sequence and not
           * from the ORF which is generated by gencode and only for internal use */
          qsq_aa->idx = info->ntqsq->prev_n + k;
          sprintf(qsq_aa->orfid, "orf%" PRId64 "", qsq_aa->idx);

          if ((status = esl_sq_SetName     (qsq_aa, info->ntqsq->name))   != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence name failed");
          if ((status = esl_sq_SetAccession(qsq_aa, info->ntqsq->acc))    != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence accession failed");
          if ((status = esl_sq_SetDesc     (qsq_aa, info->ntqsq->desc))   != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence description failed");

          qsq_aa->idx = info->ntqsq->idx;

          p7_bg_SetLength(info->bg, qsq_aa->n);

          p7_oprofile_ReconfigLength(om, qsq_aa->n);

          p7_Pipeline(info->pli, om, info->bg, qsq_aa, qsqDNATxt, info->th, NULL);
          p7_pipeline_Reuse(info->pli);

      }
      p7_oprofile_Destroy(om);

  }

  esl_sq_Destroy(qsqDNATxt);
  esl_alphabet_Destroy(abc);

  return status;
}

#ifdef HMMER_THREADS
static int
thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ID_LENGTH_LIST *id_length_list, int64_t seqidx, int64_t seqL, P7_HMMFILE *hfp)
{
  int  status   = eslOK;
  int  sstatus  = eslOK;
  int  eofCount = 0;
  P7_OM_BLOCK   *block;
  ESL_ALPHABET  *abc = NULL;
  void          *newBlock;

  add_id_length(id_length_list, seqidx, seqL);  
  
  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
      
  /* Main loop: */
  while (sstatus == eslOK)
  {
      block = (P7_OM_BLOCK *) newBlock;
      sstatus = p7_oprofile_ReadBlockMSV(hfp, &abc, block);
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
  int i, k;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  P7_OM_BLOCK   *block;
  void          *newBlock;
  ESL_SQ        *qsq_aa = NULL;       /* query sequence, amino acid                                */
  ESL_SQ        *qsqDNATxt = NULL;    /* DNA query sequence that will be in text mode for printing */


  impl_Init();

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  /* loop until all profile blocks have been processed */
  block = (P7_OM_BLOCK *) newBlock;
  while (block->count > 0)
  {

      /* Main loop over hmms */
    for (i = 0; i < block->count; ++i)
      {

      if (qsqDNATxt == NULL) {
        /* copy digital sequence to avoid race condition during esl_translate's revcomp;
         * the copy will then be textized, for use in the domain alignemnt display
         */
        qsqDNATxt = esl_sq_CreateDigitalFrom(info->ntqsq->abc,info->ntqsq->name,info->ntqsq->dsq,
                                  info->ntqsq->n,info->ntqsq->desc,info->ntqsq->acc,info->ntqsq->ss);


        /* translate DNA sequence to 6 frame ORFs  */
        translate_sequence(info->gcode, info->wrk, qsqDNATxt);
        for (k = 0; k < info->wrk->orf_block->count; ++k)
        {
            qsq_aa = &(info->wrk->orf_block->list[k]);
            p7_pli_NewSeq(info->pli, qsq_aa);
        }

        /* convert the DNA sequence to text so we can print it in the domain alignment display */
        esl_sq_Textize(qsqDNATxt);
        qsqDNATxt->abc = info->ntqsq->abc; /* add the alphabet back, since it's used in the pipeline */
      }

	  P7_OPROFILE *om = block->list[i];
	  p7_pli_NewModel(info->pli, om, info->bg);

	  /*process each 6 frame translated sequence */
      for (k = 0; k < info->wrk->orf_block->count; ++k)
        {
          qsq_aa = &(info->wrk->orf_block->list[k]);

          /* use the name, accession, and description from the DNA sequence and not 
           * from the ORF which is generated by gencode and only for internal use */
          qsq_aa->idx = info->ntqsq->prev_n + k;
          sprintf(qsq_aa->orfid, "orf%" PRId64 "", qsq_aa->idx);

          if ((status = esl_sq_SetName     (qsq_aa, info->ntqsq->name))   != eslOK)  esl_fatal("Set query sequence name failed");
          if ((status = esl_sq_SetAccession(qsq_aa, info->ntqsq->acc))    != eslOK)  esl_fatal("Set query sequence accession failed");
          if ((status = esl_sq_SetDesc     (qsq_aa, info->ntqsq->desc))   != eslOK)  esl_fatal("Set query sequence description failed");
          qsq_aa->idx = info->ntqsq->idx; 
          //printf("name %s ntqsq idx %d qsqDNAtxt idx %d\n", info->ntqsq->name, info->ntqsq->idx, qsqDNATxt->idx);
          p7_bg_SetLength(info->bg, qsq_aa->n);

          p7_oprofile_ReconfigLength(om, qsq_aa->n);

          p7_Pipeline(info->pli, om, info->bg, qsq_aa, qsqDNATxt, info->th, NULL);
          p7_pipeline_Reuse(info->pli);

        }

      p7_oprofile_Destroy(om);

      block->list[i] = NULL;
      }

    status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
    if (status != eslOK) esl_fatal("Work queue worker failed");

    block = (P7_OM_BLOCK *) newBlock;
  }

  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_sq_Destroy(qsqDNATxt);

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

  int i,d;
  int j = 0;
  for (i=0; i<th->N; i++) {
    while (th->hit[i]->seqidx != id_length_list->id_lengths[j].id) { j++; }
    for(d=0; d<th->hit[i]->ndom; d++) {
      th->hit[i]->dcl[d].ad->L = id_length_list->id_lengths[j].length;
    }
  }
  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

