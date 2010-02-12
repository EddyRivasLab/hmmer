/* jackhmmer: iterative search of a protein sequence against a protein database
 * 
 * SRE, Thu Dec 11 15:20:27 2008 [Janelia] [maestro Bear McCreary, conducting]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_keyhash.h"
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

static void checkpoint_hmm(int nquery, P7_HMM *hmm,  char *basename, int iteration);
static void checkpoint_msa(int nquery, ESL_MSA *msa, char *basename, int iteration);

#define CONOPTS "--fast,--hand"                            /* Exclusive options for model construction                    */
#define EFFOPTS "--eent,--eclust,--eset,--enone"           /* Exclusive options for effective sequence number calculation */
#define WGTOPTS "--wgsc,--wblosum,--wpb,--wnone,--wgiven"  /* Exclusive options for relative weighting                    */

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles     reqs   incomp                             help                                                  docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,  NULL,                          "show brief help on version and usage",                         1 },
  { "-N",           eslARG_INT,      "5", NULL, NULL,      NULL,    NULL,  NULL,                          "set maximum number of iterations to <n>",                      1 },
/* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,    NULL,  NULL,                          "direct output to file <f>, not stdout",                        10 },
  { "-A",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,    NULL,  NULL,                          "save multiple alignment of hits to file <s>",                  10 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,    NULL,  NULL,                          "save parseable table of per-sequence hits to file <s>",        10 },
  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,      NULL,    NULL,  NULL,                          "save parseable table of per-domain hits to file <s>",          10 },
  { "--chkhmm",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,    NULL,  NULL,                          "save HMM checkpoints to files <s>-<iteration>.hmm",            10 },
  { "--chkali",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,    NULL,  NULL,                          "save alignment checkpoints to files <s>-<iteration>.sto",      10 },
/* Control of scoring system */
  { "--popen",      eslARG_REAL,  "0.02", NULL, "0<=x<0.5",NULL,    NULL,  NULL,                          "gap open probability",                                         2 },
  { "--pextend",    eslARG_REAL,   "0.4", NULL, "0<=x<1",  NULL,    NULL,  NULL,                          "gap extend probability",                                       2 },
  { "--mxfile",     eslARG_INFILE,  NULL, NULL, NULL,      NULL,    NULL,  NULL,                          "substitution score matrix [default: BLOSUM62]",                2 },
/* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",     NULL,    NULL,  "--cut_ga,--cut_nc,--cut_tc",  "report sequences <= this E-value threshold in output",         3 },
  { "-T",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,    NULL,  "--cut_ga,--cut_nc,--cut_tc",  "report sequences >= this score threshold in output",           3 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,    NULL,  NULL,                          "set # of comparisons done, for E-value calculation",           3 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",     NULL,    NULL,  "--cut_ga,--cut_nc,--cut_tc",  "report domains <= this E-value threshold in output",           3 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,    NULL,  "--cut_ga,--cut_nc,--cut_tc",  "report domains >= this score cutoff in output",                3 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,    NULL,  NULL,                          "set # of significant seqs, for domain E-value calculation",    3 },
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,  "-E,-T,--domE,--domT",         "use profile's GA gathering cutoffs to set -T, --domT",         99 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,  "-E,-T,--domE,--domT",         "use profile's NC noise cutoffs to set -T, --domT",             99 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL,  "-E,-T,--domE,--domT",         "use profile's TC trusted cutoffs to set -T, --domT",           99 },
/* Control of inclusion thresholds */
  { "--incE",       eslARG_REAL, "0.001", NULL, "x>0",     NULL,  NULL, "--inc_ga,--inc_nc,--inc_tc",        "include sequences <= this E-value threshold in output ali",    4 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL, "--inc_ga,--inc_nc,--inc_tc",        "include sequences >= this score threshold in output ali",      4 },
  { "--incdomE",    eslARG_REAL, "0.001", NULL, "x>0",     NULL,  NULL, "--inc_ga,--inc_nc,--inc_tc",        "include domains <= this E-value threshold in output ali",      4 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL, "--inc_ga,--inc_nc,--inc_tc",        "include domains >= this score threshold in output ali",        4 },
  { "--inc_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--incE,--incT,--incdomE,--incdomT", "use profile's GA gathering cutoffs to set --incT, --incdomT",  99 },
  { "--inc_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--incE,--incT,--incdomE,--incdomT", "use profile's NC noise cutoffs to set --incT, --incdomT",      99 },
  { "--inc_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--incE,--incT,--incdomE,--incdomT", "use profile's TC trusted cutoffs to set --incT, --incdomT",    99 },
/* Control of filter pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,    NULL, "--F1,--F2,--F3",               "Turn all heuristic filters off (less speed, more power)",      5 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,      NULL,    NULL, "--max",                        "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             5 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,      NULL,    NULL, "--max",                        "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             5 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,      NULL,    NULL, "--max",                        "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             5 },
  { "--nobias",     eslARG_NONE,    NULL, NULL, NULL,      NULL,    NULL, "--max",                        "turn off composition bias filter",                             5 },
  { "--nonull2",    eslARG_NONE,    NULL, NULL, NULL,      NULL,    NULL,    NULL,                        "turn off biased composition score corrections",                5 },
/* Alternate model construction strategies */
  { "--fast",       eslARG_NONE,   FALSE, NULL, NULL,   CONOPTS,    NULL,    NULL, "assign cols w/ >= symfrac residues as consensus",       6 },
  { "--hand",       eslARG_NONE,"default",NULL, NULL,   CONOPTS,    NULL,    NULL, "manual construction (requires reference annotation)",   6 },
  { "--symfrac",    eslARG_REAL,   "0.5", NULL, "0<=x<=1", NULL,"--fast",    NULL, "sets sym fraction controlling --fast construction",     6 },
  { "--fragthresh", eslARG_REAL,   "0.5", NULL, "0<=x<=1", NULL,    NULL,    NULL, "if L < x<L>, tag sequence as a fragment",               6 },
/* Alternate relative sequence weighting strategies */
  { "--wpb",        eslARG_NONE,"default",NULL, NULL,   WGTOPTS,    NULL,    NULL, "Henikoff position-based weights",                      7 },
  { "--wgsc",       eslARG_NONE,    NULL, NULL, NULL,   WGTOPTS,    NULL,    NULL, "Gerstein/Sonnhammer/Chothia tree weights",             7 },
  { "--wblosum",    eslARG_NONE,    NULL, NULL, NULL,   WGTOPTS,    NULL,    NULL, "Henikoff simple filter weights",                       7 },
  { "--wnone",      eslARG_NONE,    NULL, NULL, NULL,   WGTOPTS,    NULL,    NULL, "don't do any relative weighting; set all to 1",        7 },
  { "--wgiven",     eslARG_NONE,    NULL, NULL, NULL,   WGTOPTS,    NULL,    NULL, "use weights as given in MSA file",                    99 }, /* no-op in jackhmmer */
  { "--wid",        eslARG_REAL,  "0.62", NULL,"0<=x<=1",  NULL,"--wblosum", NULL, "for --wblosum: set identity cutoff",                   7 },
/* Alternate effective sequence weighting strategies */
  { "--eent",       eslARG_NONE,"default",NULL, NULL,   EFFOPTS,    NULL,    NULL, "adjust eff seq # to achieve relative entropy target",  8 },
  { "--eclust",     eslARG_NONE,   FALSE, NULL, NULL,   EFFOPTS,    NULL,    NULL, "eff seq # is # of single linkage clusters",            8 },
  { "--enone",      eslARG_NONE,   FALSE, NULL, NULL,   EFFOPTS,    NULL,    NULL, "no effective seq # weighting: just use nseq",          8 },
  { "--eset",       eslARG_REAL,    NULL, NULL, NULL,   EFFOPTS,    NULL,    NULL, "set eff seq # for all models to <x>",                  8 },
  { "--ere",        eslARG_REAL,    NULL, NULL,"x>0",      NULL, "--eent",   NULL, "for --eent: set minimum rel entropy/position to <x>",  8 },
  { "--esigma",     eslARG_REAL,  "45.0", NULL,"x>0",      NULL, "--eent",   NULL, "for --eent: set sigma param to <x>",                   8 },
  { "--eid",        eslARG_REAL,  "0.62", NULL,"0<=x<=1",  NULL,"--eclust",  NULL, "for --eclust: set fractional identity cutoff to <x>",  8 },
/* Control of E-value calibration */
  { "--EmL",         eslARG_INT,   "200", NULL,"n>0",      NULL,    NULL,    NULL, "length of sequences for MSV Gumbel mu fit",                    9 },   
  { "--EmN",         eslARG_INT,   "200", NULL,"n>0",      NULL,    NULL,    NULL, "number of sequences for MSV Gumbel mu fit",                    9 },   
  { "--EvL",         eslARG_INT,   "200", NULL,"n>0",      NULL,    NULL,    NULL, "length of sequences for Viterbi Gumbel mu fit",                9 },   
  { "--EvN",         eslARG_INT,   "200", NULL,"n>0",      NULL,    NULL,    NULL, "number of sequences for Viterbi Gumbel mu fit",                9 },   
  { "--EfL",         eslARG_INT,   "100", NULL,"n>0",      NULL,    NULL,    NULL, "length of sequences for Forward exp tail tau fit",             9 },   
  { "--EfN",         eslARG_INT,   "200", NULL,"n>0",      NULL,    NULL,    NULL, "number of sequences for Forward exp tail tau fit",             9 },   
  { "--Eft",         eslARG_REAL, "0.04", NULL,"0<x<1",    NULL,    NULL,    NULL, "tail mass for Forward exponential tail tau fit",               9 },   
/* Other options */
  { "--acc",        eslARG_NONE,  FALSE,  NULL, NULL,      NULL,    NULL,    NULL,   "output target accessions instead of names if possible",     11 },
  { "--seed",        eslARG_INT,    "42", NULL, "n>=0",    NULL,    NULL,    NULL,   "set RNG seed to <n> (if 0: one-time arbitrary seed)",       11 },
  { "--textw",       eslARG_INT,   "120", NULL, "n>=120",  NULL,    NULL,"--notextw","set max width of ASCII text output lines",                  11 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,      NULL,    NULL,"--textw",  "unlimit ASCII text output line width",                      11 },
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,      NULL,    NULL,   NULL,    "assert query <seqfile> is in format <s>: no autodetection", 11 },
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,      NULL,    NULL,   NULL,    "assert target <seqdb> is in format <s>>: no autodetection", 11 },
#ifdef HMMER_THREADS
  { "--cpu",        eslARG_INT,     NULL,"HMMER_NCPU", "n>0",      NULL,    NULL,   NULL,    "number of worker threads",                                  11 },
#endif
#ifdef HAVE_MPI
  // { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "arrest after start: for debugging MPI under gdb",          4 },  
  // { "--mpi",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "run as an MPI parallel program",                           4 },
#endif 
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <query seqfile> <target seqdb>";
static char banner[] = "iteratively search a protein sequence against a protein database";

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
      esl_opt_DisplayHelp(stdout, go, 1, 2, 120); /* 1= group; 2 = indentation; 120=textwidth*/

      puts("\nOptions directing output:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80); 

      puts("\nOptions controlling scoring system in iteration one:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 120); 

      puts("\nOptions controlling significance thresholds for reporting:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 120); 

      puts("\nOptions controlling significance thresholds for inclusion in next round:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 120); 

      puts("\nOptions controlling acceleration heuristics:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 120); 

      puts("\nOptions controlling model construction after first iteration:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 120); 

      puts("\nOptions controlling relative weights in models after first iteration:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 120); 

      puts("\nOptions controlling effective seq number in models after first iteration:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 120); 

      puts("\nOptions controlling E value calibration:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 120); 

      puts("\nOther expert options:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 120); 
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
  if (esl_opt_IsUsed(go, "-N"))          fprintf(ofp, "# maximum iterations set to:       %d\n",      esl_opt_GetInteger(go, "-N"));
  if (esl_opt_IsUsed(go, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",      esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "-A"))          fprintf(ofp, "# MSA of hits saved to file:       %s\n",      esl_opt_GetString(go, "-A"));
  if (esl_opt_IsUsed(go, "--tblout"))    fprintf(ofp, "# per-seq hits tabular output:     %s\n",      esl_opt_GetString(go, "--tblout"));
  if (esl_opt_IsUsed(go, "--domtblout")) fprintf(ofp, "# per-dom hits tabular output:     %s\n",      esl_opt_GetString(go, "--domtblout"));
  if (esl_opt_IsUsed(go, "--chkhmm"))    fprintf(ofp, "# HMM checkpoint files output:     %s-<i>.hmm\n", esl_opt_GetString(go, "--chkhmm"));
  if (esl_opt_IsUsed(go, "--chkali"))    fprintf(ofp, "# MSA checkpoint files output:     %s-<i>.sto\n", esl_opt_GetString(go, "--chkali"));
  if (esl_opt_IsUsed(go, "--popen"))     fprintf(ofp, "# gap open probability:            %f\n",      esl_opt_GetReal  (go, "--popen"));
  if (esl_opt_IsUsed(go, "--pextend"))   fprintf(ofp, "# gap extend probability:          %f\n",      esl_opt_GetReal  (go, "--pextend"));
  if (esl_opt_IsUsed(go, "--mxfile"))    fprintf(ofp, "# subst score matrix:              %s\n",      esl_opt_GetString(go, "--mxfile"));
  if (esl_opt_IsUsed(go, "-E"))          fprintf(ofp, "# sequence E-value threshold:   <= %g\n",      esl_opt_GetReal(go, "-E"));
  if (esl_opt_IsUsed(go, "-T"))          fprintf(ofp, "# sequence bit score threshold: <= %g\n",      esl_opt_GetReal(go, "-T"));
  if (esl_opt_IsUsed(go, "-Z"))          fprintf(ofp, "# sequence search space set to:    %.0f\n",    esl_opt_GetReal(go, "-Z"));
  if (esl_opt_IsUsed(go, "--domE"))      fprintf(ofp, "# domain E-value threshold:     <= %g\n",      esl_opt_GetReal(go, "--domE"));
  if (esl_opt_IsUsed(go, "--domT"))      fprintf(ofp, "# domain bit score threshold:   <= %g\n",      esl_opt_GetReal(go, "--domT"));
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
  if (esl_opt_IsUsed(go, "--fast"))      fprintf(ofp, "# model architecture construction: fast/heuristic\n");
  if (esl_opt_IsUsed(go, "--hand"))      fprintf(ofp, "# model architecture construction: hand-specified by RF annotation\n");
  if (esl_opt_IsUsed(go, "--symfrac"))   fprintf(ofp, "# sym frac for model structure:    %.3f\n", esl_opt_GetReal(go, "--symfrac"));
  if (esl_opt_IsUsed(go, "--fragthresh"))fprintf(ofp, "# define fragments if < xL    :    %.3f\n", esl_opt_GetReal(go, "--fragthresh"));
  if (esl_opt_IsUsed(go, "--wpb"))       fprintf(ofp, "# relative weighting scheme:       Henikoff PB\n");
  if (esl_opt_IsUsed(go, "--wgsc"))      fprintf(ofp, "# relative weighting scheme:       G/S/C\n");
  if (esl_opt_IsUsed(go, "--wblosum"))   fprintf(ofp, "# relative weighting scheme:       BLOSUM filter\n");
  if (esl_opt_IsUsed(go, "--wnone"))     fprintf(ofp, "# relative weighting scheme:       none\n");
  if (esl_opt_IsUsed(go, "--wid"))       fprintf(ofp, "# frac id cutoff for BLOSUM wgts:  %f\n",   esl_opt_GetReal(go, "--wid"));
  if (esl_opt_IsUsed(go, "--eent"))      fprintf(ofp, "# effective seq number scheme:     entropy weighting\n");
  if (esl_opt_IsUsed(go, "--eclust"))    fprintf(ofp, "# effective seq number scheme:     single linkage clusters\n");
  if (esl_opt_IsUsed(go, "--enone"))     fprintf(ofp, "# effective seq number scheme:     none\n");
  if (esl_opt_IsUsed(go, "--eset"))      fprintf(ofp, "# effective seq number:            set to %f\n", esl_opt_GetReal(go, "--eset"));
  if (esl_opt_IsUsed(go, "--ere") )      fprintf(ofp, "# minimum rel entropy target:      %f bits\n",   esl_opt_GetReal(go, "--ere"));
  if (esl_opt_IsUsed(go, "--esigma") )   fprintf(ofp, "# entropy target sigma parameter:  %f bits\n",   esl_opt_GetReal(go, "--esigma"));
  if (esl_opt_IsUsed(go, "--eid") )      fprintf(ofp, "# frac id cutoff for --eclust:     %f\n",        esl_opt_GetReal(go, "--eid"));
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
  ESL_GETOPTS     *go       = NULL;               /* application configuration options               */
  FILE            *ofp      = stdout;             /* output file for results (default stdout)        */
  FILE            *afp      = NULL;               /* alignment output file (-A option)               */
  FILE            *tblfp    = NULL;		  /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;		  /* output stream for tabular per-seq (--domtblout) */
  char            *qfile    = NULL;               /* file to read query sequence from                */
  char            *dbfile   = NULL;               /* file to read sequence(s) from                   */
  int              qformat  = eslSQFILE_UNKNOWN;  /* format of qfile                                 */
  int              dbformat = eslSQFILE_UNKNOWN;  /* format of dbfile                                */
  ESL_SQFILE      *qfp      = NULL;		  /* open qfile                                      */
  ESL_SQFILE      *dbfp     = NULL;               /* open dbfile                                     */
  ESL_ALPHABET    *abc      = NULL;               /* sequence alphabet                               */
  P7_BUILDER      *bld      = NULL;               /* HMM construction configuration                  */
  ESL_SQ          *qsq      = NULL;               /* query sequence                                  */
  ESL_SQ          *dbsq     = NULL;               /* target sequence                                 */
  ESL_KEYHASH     *kh       = NULL;		  /* hash of previous top hits' ranks                */
  ESL_STOPWATCH   *w        = NULL;               /* for timing                                      */
  int              nquery   = 0;
  int              textw;
  int              iteration;
  int              maxiterations;
  int              nnew_targets;
  int              prv_msa_nseq;
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
  abc           = esl_alphabet_Create(eslAMINO);
  w             = esl_stopwatch_Create();
  kh            = esl_keyhash_Create();
  maxiterations = esl_opt_GetInteger(go, "-N");
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

  /* Initialize builder configuration */
  bld = p7_builder_Create(go, abc);
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
  
  if (! esl_sqfile_IsRewindable(dbfp)) 
    esl_fatal("Target sequence file %s isn't rewindable; jackhmmer requires that it is", dbfile);

  /* Open the query sequence file  */
  status = esl_sqfile_OpenDigital(abc, qfile, qformat, NULL, &qfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open sequence file %s for reading\n",      qfile);
  else if (status == eslEFORMAT)   esl_fatal("Sequence file %s is empty or misformatted\n",        qfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        esl_fatal ("Unexpected error %d opening sequence file %s\n", status, qfile);
  qsq = esl_sq_CreateDigital(abc);

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

  /* Ready to begin */
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

  /* Outer loop over sequence queries, if more than one */
  while ((qstatus = esl_sqio_Read(qfp, qsq)) == eslOK)
    {
      P7_HMM          *hmm     = NULL;	     /* HMM - only needed if checkpointed        */
      P7_HMM         **ret_hmm = NULL;	     /* HMM - only needed if checkpointed        */
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */
      P7_TRACE        *qtr     = NULL;       /* faux trace for query sequence            */
      ESL_MSA         *msa     = NULL;       /* multiple alignment of included hits      */
      
      if (esl_opt_IsOn(go, "--chkhmm")) ret_hmm = &hmm;

      nquery++;
      if (qsq->n == 0) continue; /* skip zero length queries as if they aren't even present. */

      fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n);
      if (qsq->acc[0]  != '\0') fprintf(ofp, "Accession:   %s\n", qsq->acc);
      if (qsq->desc[0] != '\0') fprintf(ofp, "Description: %s\n", qsq->desc);  
      fprintf(ofp, "\n");

      for (iteration = 1; iteration <= maxiterations; iteration++)
	{       /* We enter each iteration with an optimized profile. */
	  esl_stopwatch_Start(w);

	  if (om  != NULL)       p7_oprofile_Destroy(om);

	  if (info->pli != NULL) p7_pipeline_Destroy(info->pli);
	  if (info->th  != NULL) p7_tophits_Destroy(info->th);
	  if (info->om  != NULL) p7_oprofile_Destroy(info->om);

 	  /* Create the search model: from query alone (round 1) or from MSA (round 2+) */
	  if (msa == NULL)	/* round 1 */
	    {
	      p7_SingleBuilder(bld, qsq, info->bg, ret_hmm, &qtr, NULL, &om); /* bypass HMM - only need model */

	      prv_msa_nseq = 1;
	    }
	  else 
	    {
	      /* Throw away old model. Build new one. */
	      status = p7_Builder(bld, msa, info->bg, ret_hmm, NULL, NULL, &om, NULL);
	      if      (status == eslENORESULT) esl_fatal("Failed to construct new model from iteration %d results:\n%s", iteration, bld->errbuf);
	      else if (status == eslEFORMAT)   esl_fatal("Failed to construct new model from iteration %d results:\n%s", iteration, bld->errbuf);
	      else if (status != eslOK)        esl_fatal("Unexpected error constructing new model at iteration %d:",     iteration);

	      fprintf(ofp, "@@\n");
	      fprintf(ofp, "@@ Round:                  %d\n", iteration);
	      fprintf(ofp, "@@ Included in MSA:        %d subsequences (query + %d subseqs from %d targets)\n", 
		      msa->nseq, msa->nseq-1, kh->nkeys);
	      fprintf(ofp, "@@ Model size:             %d positions\n", om->M);
	      fprintf(ofp, "@@\n\n");
	  
	      prv_msa_nseq = msa->nseq;
	      esl_msa_Destroy(msa);
	    }

	  /* HMM checkpoint output */
	  if (esl_opt_IsOn(go, "--chkhmm")) {
	    checkpoint_hmm(nquery, hmm, esl_opt_GetString(go, "--chkhmm"), iteration);
	    p7_hmm_Destroy(hmm);
	    hmm = NULL;
	  }

	  /* Create new processing pipeline and top hits list; destroy old. (TODO: reuse rather than recreate) */
	  for (i = 0; i < ncpus; ++i)
	    {
	      info[i].th  = p7_tophits_Create(); 
	      info[i].om  = p7_oprofile_Clone(om);
	      info[i].pli = p7_pipeline_Create(go, om->M, 400, p7_SEARCH_SEQS); /* 400 is a dummy length for now */
	      p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);

#ifdef HMMER_THREADS
	      esl_threads_AddThread(threadObj, &info[i]);
#endif
	    }

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

	  /* Print the results. */
	  p7_tophits_Sort(info->th);
	  p7_tophits_Threshold(info->th, info->pli);
	  p7_tophits_CompareRanking(info->th, kh, &nnew_targets);
	  p7_tophits_Targets(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
	  p7_tophits_Domains(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");

	  /* Create alignment of the top hits */
	  p7_tophits_Alignment(info->th, abc, &qsq, &qtr, 1, p7_ALL_CONSENSUS_COLS, &msa);
	  esl_msa_Digitize(abc,msa,NULL);
	  esl_msa_FormatName(msa, "%s-i%d", qsq->name, iteration);

	  /* Optional checkpointing */
	  if (esl_opt_IsOn(go, "--chkali")) checkpoint_msa(nquery, msa, esl_opt_GetString(go, "--chkali"), iteration);

	  esl_stopwatch_Stop(w);
	  p7_pli_Statistics(ofp, info->pli, w);
	  fprintf(ofp, "\n");

	  /* Convergence test */
	  fprintf(ofp, "@@ New targets included:   %d\n", nnew_targets);
	  fprintf(ofp, "@@ New alignment includes: %d subseqs (was %d), including original query\n", 
		  msa->nseq, prv_msa_nseq);
	  if (nnew_targets == 0 && msa->nseq <= prv_msa_nseq) 
	    { 
	      fprintf(ofp, "@@\n");
	      fprintf(ofp, "@@ CONVERGED (in %d rounds). \n", iteration);
	      fprintf(ofp, "@@\n\n");
	      break;
	    }
	  else if (iteration < maxiterations)
	    { fprintf(ofp, "@@ Continuing to next round.\n\n"); }

	  esl_sqfile_Position(dbfp, 0);
	} /* end iteration loop */

      /* Because we destroy/create the hitlist, om, pipeline, and msa above, rather than create/destroy,
       * the results of the last iteration have carried through to us now, and we can output
       * whatever final results we care to.
       */
      if (tblfp)    p7_tophits_TabularTargets(tblfp,    qsq->name, qsq->acc, info->th, info->pli, (nquery == 1));
      if (domtblfp) p7_tophits_TabularDomains(domtblfp, qsq->name, qsq->acc, info->th, info->pli, (nquery == 1));
      if (afp) 
	{
	  if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
	  else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);

	  fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
	}
      fprintf(ofp, "//\n");

      p7_pipeline_Destroy(info->pli);
      p7_tophits_Destroy(info->th);
      p7_oprofile_Destroy(info->om);

      esl_msa_Destroy(msa);
      p7_oprofile_Destroy(om);
      p7_trace_Destroy(qtr);
      esl_sq_Reuse(qsq);
      esl_keyhash_Reuse(kh);
      esl_sqfile_Position(dbfp, 0);
    }
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
  while (esl_workqueue_Remove(queue, (void **) &block) == eslOK)
    {
      esl_sq_DestroyBlock(block);
    }
  esl_workqueue_Destroy(queue);
  esl_threads_Destroy(threadObj);
#endif

  free(info);

  esl_keyhash_Destroy(kh);
  esl_sqfile_Close(qfp);
  esl_sqfile_Close(dbfp);
  esl_sq_Destroy(dbsq);
  esl_sq_Destroy(qsq);  
  esl_stopwatch_Destroy(w);
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



/* checkpoint_hmm()
 * Incept:    SRE, Tue Jun  2 10:20:08 2009 [Janelia]
 *
 * Purpose:   Save <hmm> to a file <basename>-<iteration>.hmm.
 *            If <nquery == 1>, start a new checkpoint file;
 *            for <nquery > 1>, append to existing one.
 */
static void
checkpoint_hmm(int nquery, P7_HMM *hmm, char *basename, int iteration)
{
  FILE *fp         = NULL;
  char *filename   = NULL;

  esl_sprintf(&filename, "%s-%d.hmm", basename, iteration);
  if (nquery == 1) { if ((fp = fopen(filename, "w")) == NULL) esl_fatal("Failed to open HMM checkpoint file %s for writing\n", filename); }
  else             { if ((fp = fopen(filename, "a")) == NULL) esl_fatal("Failed to open HMM checkpoint file %s for append\n",  filename); }
  p7_hmmfile_WriteASCII(fp, -1, hmm);
  
  fclose(fp);
  free(filename);
  return;
}


/* checkpoint_msa()
 * Incept:    SRE, Tue Jun  2 10:35:38 2009 [Janelia]
 *
 * Purpose:   Save <msa> to a file <basename>-<iteration>.sto.
 *            If <nquery == 1>, start a new checkpoint file;
 *            for <nquery > 1>, append to existing one.
 */
static void
checkpoint_msa(int nquery, ESL_MSA *msa, char *basename, int iteration)
{
  FILE *fp         = NULL;
  char *filename   = NULL;

  esl_sprintf(&filename, "%s-%d.sto", basename, iteration);
  if (nquery == 1) { if ((fp = fopen(filename, "w")) == NULL) esl_fatal("Failed to open MSA checkpoint file %s for writing\n", filename); }
  else             { if ((fp = fopen(filename, "a")) == NULL) esl_fatal("Failed to open MSA checkpoint file %s for append\n",  filename); }
  esl_msa_Write(fp, msa, eslMSAFILE_PFAM);
  
  fclose(fp);
  free(filename);
  return;

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
