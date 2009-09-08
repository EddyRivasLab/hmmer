/* hmmscan: search sequence(s) against a profile HMM database
 * 
 * SRE, Mon Oct 20 08:28:05 2008 [Janelia]
 * SVN $Id$
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
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif

#include "hmmer.h"

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif
  ESL_SQ           *qsq;
  P7_BG            *bg;
  P7_PIPELINE      *pli;
  P7_TOPHITS       *th;
} WORKER_INFO;

#ifdef HMMER_THREADS
#define BLOCK_SIZE 500

static int  threadedLoop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, P7_HMMFILE *hfp);
static void pipelineThread(void *arg);

#else
static int serialLoop(WORKER_INFO *info, P7_HMMFILE *hfp);
#endif

static ESL_OPTIONS options[] = {
  /* name           type          default  env  range toggles  reqs   incomp                         help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL,  NULL,                                   "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL,                                "direct output to file <f>, not stdout",                        7 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL,                                "save parseable table of per-sequence hits to file <s>",        7 },
  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL,                                "save parseable table of per-domain hits to file <s>",          7 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",        "report sequences <= this E-value threshold in output",         2 },
  { "-T",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",        "report sequences >= this score threshold in output",           2 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                                "set # of comparisons done, for E-value calculation",           2 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",        "report domains <= this E-value threshold in output",           2 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--cut_ga,--cut_nc,--cut_tc",        "report domains >= this score cutoff in output",                2 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  NULL,                                "set # of significant seqs, for domain E-value calculation",    2 },
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "-E,-T,--domE,--domT",               "use profile's GA gathering cutoffs to set -T, --domT",         2 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "-E,-T,--domE,--domT",               "use profile's NC noise cutoffs to set -T, --domT",             2 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "-E,-T,--domE,--domT",               "use profile's TC trusted cutoffs to set -T, --domT",           2 },
  /* Control of inclusion thresholds: in hmmscan, these make no sense, so hide them (group 99) but they must be present for the P7_PIPELINE */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",     NULL,  NULL,  "--inc_ga,--inc_nc,--inc_tc",        "include sequences <= this E-value threshold in output ali",    99 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--inc_ga,--inc_nc,--inc_tc",        "include sequences >= this score threshold in output ali",      99 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",     NULL,  NULL,  "--inc_ga,--inc_nc,--inc_tc",        "include domains <= this E-value threshold in output ali",      99 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  NULL,  "--inc_ga,--inc_nc,--inc_tc",        "include domains >= this score threshold in output ali",        99 },
  { "--inc_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--incE,--incT,--incdomE,--incdomT", "use profile's GA gathering cutoffs to set --incT, --incdomT",  99 },
  { "--inc_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--incE,--incT,--incdomE,--incdomT", "use profile's NC noise cutoffs to set --incT, --incdomT",      99 },
  { "--inc_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  "--incE,--incT,--incdomE,--incdomT", "use profile's TC trusted cutoffs to set --incT, --incdomT",    99 },
  /* Control of filter pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, "--F1,--F2,--F3",   "Turn all heuristic filters off (less speed, more power)",      4 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,   NULL,  NULL, "--max",            "MSV threshold: promote hits w/ P <= F1",                       4 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,   NULL,  NULL, "--max",            "Vit threshold: promote hits w/ P <= F2",                       4 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,   NULL,  NULL, "--max",            "Fwd threshold: promote hits w/ P <= F3",                       4 },
  { "--nobias",     eslARG_NONE,    NULL, NULL, NULL,   NULL,  NULL, "--max",            "turn off composition bias filter",                             4 },
  { "--nonull2",    eslARG_NONE,    NULL, NULL, NULL,   NULL,  NULL,    NULL,            "turn off biased composition score corrections",                4 },
  /* Other options */
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",    NULL,  NULL,    NULL,         "set RNG seed to <n> (if 0: one-time arbitrary seed)",          6 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",  NULL,  NULL,  "--notextw",    "set max width of ASCII text output lines",                     6 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,      NULL,  NULL,  "--textw",      "unlimit ASCII text output line width",                         6 },
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,   NULL,   NULL,   NULL,            "assert input <seqfile> is in format <s>: no autodetection",    6 },
#ifdef HMMER_THREADS
  { "--cpu",        eslARG_INT,     NULL,"HMMER_NCPU", "n>0",    NULL,   NULL,   NULL,            "number of worker threads",                                     6 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmm database> <query seqfile>";
static char banner[] = "search sequence(s) against a profile HMM database";



/* process_commandline()
 * 
 * Processes the commandline, filling in fields in <cfg> and creating and returning
 * an <ESL_GETOPTS> options structure. The help page (hmmsearch -h) is formatted
 * here.
 */
static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_hmmfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     p7_Die("Internal failure creating options object");
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere most common options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

      puts("\nOptions directing output:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      puts("\nOptions controlling reporting thresholds:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      puts("\nOptions controlling acceleration heuristics:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      puts("\nOther expert options:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)      { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <hmmfile> argument on command line"); goto ERROR; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }
  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere most common options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}


static int
output_header(FILE *ofp, ESL_GETOPTS *go, char *hmmfile, char *seqfile)
{
  p7_banner(ofp, go->argv[0], banner);
  
  fprintf(ofp, "# query sequence file:             %s\n", seqfile);
  fprintf(ofp, "# target HMM database:             %s\n", hmmfile);
  if (esl_opt_IsUsed(go, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",      esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "--tblout"))    fprintf(ofp, "# per-seq hits tabular output:     %s\n",      esl_opt_GetString(go, "--tblout"));
  if (esl_opt_IsUsed(go, "--domtblout")) fprintf(ofp, "# per-dom hits tabular output:     %s\n",      esl_opt_GetString(go, "--domtblout"));
  if (esl_opt_IsUsed(go, "-E"))          fprintf(ofp, "# sequence E-value threshold:   <= %g\n",      esl_opt_GetReal(go, "-E"));
  if (esl_opt_IsUsed(go, "-T"))          fprintf(ofp, "# sequence bit score threshold: >= %g\n",      esl_opt_GetReal(go, "-T"));
  if (esl_opt_IsUsed(go, "-Z"))          fprintf(ofp, "# sequence search space set to:    %.0f\n",    esl_opt_GetReal(go, "-Z"));
  if (esl_opt_IsUsed(go, "--domE"))      fprintf(ofp, "# domain E-value threshold:     <= %g\n",      esl_opt_GetReal(go, "--domE"));
  if (esl_opt_IsUsed(go, "--domT"))      fprintf(ofp, "# domain bit score threshold:   >= %g\n",      esl_opt_GetReal(go, "--domT"));
  if (esl_opt_IsUsed(go, "--domZ"))      fprintf(ofp, "# domain search space set to:      %.0f\n",    esl_opt_GetReal(go, "--domZ"));
  if (esl_opt_IsUsed(go, "--cut_ga"))    fprintf(ofp, "# set reporting thresholds to:     GA cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--cut_nc"))    fprintf(ofp, "# set reporting thresholds to:     NC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--cut_tc"))    fprintf(ofp, "# set reporting thresholds to:     TC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--incE"))      fprintf(ofp, "# seq inclusion E-val thresh:   <= %g\n",      esl_opt_GetReal(go, "--incE"));
  if (esl_opt_IsUsed(go, "--incT"))      fprintf(ofp, "# seq inclusion score thresh:   >= %g\n",      esl_opt_GetReal(go, "--incT"));
  if (esl_opt_IsUsed(go, "--incdomE"))   fprintf(ofp, "# dom inclusion E-val thresh:   <= %g\n",      esl_opt_GetReal(go, "--incdomE"));
  if (esl_opt_IsUsed(go, "--incdomT"))   fprintf(ofp, "# dom inclusion score thresh:   >= %g\n",      esl_opt_GetReal(go, "--incdomT"));
  if (esl_opt_IsUsed(go, "--inc_ga"))    fprintf(ofp, "# set inclusion thresholds to:     GA cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--inc_nc"))    fprintf(ofp, "# set inclusion thresholds to:     NC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--inc_tc"))    fprintf(ofp, "# set inclusion thresholds to:     TC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--max"))       fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n");
  if (esl_opt_IsUsed(go, "--F1"))        fprintf(ofp, "# MSV filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F1"));
  if (esl_opt_IsUsed(go, "--F2"))        fprintf(ofp, "# Vit filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F2"));
  if (esl_opt_IsUsed(go, "--F3"))        fprintf(ofp, "# Fwd filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F3"));
  if (esl_opt_IsUsed(go, "--nonull2"))   fprintf(ofp, "# null2 bias corrections:          off\n");
  if (esl_opt_IsUsed(go, "--nobias"))    fprintf(ofp, "# biased composition HMM filter:   off\n");
  if (esl_opt_IsUsed(go, "--nonull2"))   fprintf(ofp, "# null2 bias corrections:          off\n");
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed")==0)fprintf(ofp, "# random number seed:              one-time arbitrary\n");
    else                                    fprintf(ofp, "# random number seed set to:       %d\n", esl_opt_GetInteger(go, "--seed"));
  }
  if (esl_opt_IsUsed(go, "--textw"))     fprintf(ofp, "# max ASCII text line length:      %d\n",     esl_opt_GetInteger(go, "--textw"));
  if (esl_opt_IsUsed(go, "--notextw"))   fprintf(ofp, "# max ASCII text line length:      unlimited\n");
  if (esl_opt_IsUsed(go, "--qformat"))   fprintf(ofp, "# input seqfile format asserted:   %s\n", esl_opt_GetString(go, "--qformat"));
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}



int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go	    = NULL;              /* command line processing                         */
  FILE            *ofp      = stdout;	         /* output file for results (default stdout)        */
  FILE            *tblfp    = NULL;		 /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;	  	 /* output stream for tabular per-seq (--domtblout) */
  char            *seqfile  = NULL;              /* file to read query sequence(s) from             */
  int              seqfmt   = eslSQFILE_UNKNOWN; /* format of seqfile                               */
  ESL_SQFILE      *sqfp     = NULL;              /* open seqfile                                    */
  P7_HMMFILE      *hfp      = NULL;		 /* open HMM database file                          */
  char            *hmmfile  = NULL;              /* file to read HMM(s) from                        */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                               */
  P7_OPROFILE     *om       = NULL;		 /* target profile                                  */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                          */
  ESL_SQ          *qsq      = NULL;		 /* query sequence                                  */
  int              nquery   = 0;
  int              textw;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;

  int              ncpus    = 1;

  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  P7_OM_BLOCK     *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif

  /* Initializations */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */
  process_commandline(argc, argv, &go, &hmmfile, &seqfile);    
  w = esl_stopwatch_Create();
  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* If caller declared an input format, decode it */
  if (esl_opt_IsOn(go, "--qformat")) {
    seqfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (seqfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }

  /* Open the target profile database to get the sequence alphabet */
  status = p7_hmmfile_Open(hmmfile, p7_HMMDBENV, &hfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open hmm file %s for reading.\n",                       hmmfile);
  else if (status == eslEFORMAT)   p7_Fail("Unrecognized format, trying to open hmm file %s for reading.\n",  hmmfile);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening hmm file %s.\n",           status, hmmfile);  
  if (! hfp->is_pressed)           p7_Fail("Failed to open binary dbs for HMM file %s: use hmmpress first\n", hmmfile);
  
  hstatus = p7_oprofile_ReadMSV(hfp, &abc, &om);
  if      (hstatus == eslEFORMAT)   p7_Fail("bad file format in HMM file %s",             hmmfile);
  else if (hstatus == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets",   hmmfile);
  else if (hstatus != eslOK)        p7_Fail("Unexpected error in reading HMMs from %s",   hmmfile); 

  p7_oprofile_Destroy(om);
  p7_hmmfile_Close(hfp);

  /* Open the query sequence database */
  status = esl_sqfile_OpenDigital(abc, seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",      seqfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",        seqfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, seqfile);
  qsq = esl_sq_CreateDigital(abc);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"),          "w")) == NULL)  esl_fatal("Failed to open output file %s for writing\n",                 esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp")); }
  if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblfp")); }
 
  output_header(ofp, go, hmmfile, seqfile);

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

  for (i = 0; i < ncpus; ++i)
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
  while ( (sstatus = esl_sqio_Read(sqfp, qsq)) == eslOK)
    {
      nquery++;
      esl_stopwatch_Start(w);	                          

      /* Open the target profile database */
      status = p7_hmmfile_Open(hmmfile, p7_HMMDBENV, &hfp);
      if (status != eslOK)        p7_Fail("Unexpected error %d in opening hmm file %s.\n",           status, hmmfile);  
  
#ifdef HMMER_THREADS
      /* if we are threaded, create a lock to prevent multiple readers */
      status = p7_hmmfile_CreateLock(hfp);
      if (status != eslOK) p7_Fail("Unexpected error %d createing lock\n", status);
#endif

      fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n);
      if (qsq->acc[0]  != 0) fprintf(ofp, "Accession:   %s\n", qsq->acc);
      if (qsq->desc[0] != 0) fprintf(ofp, "Description: %s\n", qsq->desc);

      for (i = 0; i < ncpus; ++i)
	{
	  /* Create processing pipeline and hit list */
	  info[i].th  = p7_tophits_Create(); 
	  info[i].pli = p7_pipeline_Create(go, 100, 100, p7_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
	  info[i].pli->hfp = hfp;  /* for two-stage input, pipeline needs <hfp> */

	  p7_pli_NewSeq(info[i].pli, qsq);
	  p7_bg_SetLength(info[i].bg, qsq->n);
	  info[i].qsq = qsq;

#ifdef HMMER_THREADS
	  esl_threads_AddThread(threadObj, &info[i]);
#endif
	}

#ifdef HMMER_THREADS
      hstatus = threadedLoop(threadObj, queue, hfp);
#else
      hstatus = serialLoop(info, hfp);
#endif
      switch(hstatus)
	{
	case eslEFORMAT:
	  p7_Fail("bad file format in HMM file %s",             hmmfile);
	  break;
	case eslEINCOMPAT:
	  p7_Fail("HMM file %s contains different alphabets",   hmmfile);
	  break;
	case eslEOF:
	  /* do nothing */
	  break;
	default:
	  p7_Fail("Unexpected error in reading HMMs from %s",   hmmfile); 
	}

      /* merge the results of the search results */
      for (i = 1; i < ncpus; ++i)
	{
	  p7_tophits_Merge(info[0].th, info[i].th);
	  p7_pipeline_Merge(info[0].pli, info[i].pli);

	  p7_pipeline_Destroy(info[i].pli);
	  p7_tophits_Destroy(info[i].th);
	}

      /* Print results */
      p7_tophits_Sort(info->th);
      p7_tophits_Threshold(info->th, info->pli);
      p7_tophits_Targets(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
      p7_tophits_Domains(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");

      if (tblfp)    p7_tophits_TabularTargets(tblfp,    qsq->name, info->th, info->pli, (nquery == 1));
      if (domtblfp) p7_tophits_TabularDomains(domtblfp, qsq->name, info->th, info->pli, (nquery == 1));

      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, info->pli, w);
      fprintf(ofp, "//\n");

      p7_hmmfile_Close(hfp);
      p7_pipeline_Destroy(info->pli);
      p7_tophits_Destroy(info->th);
      esl_sq_Reuse(qsq);
    }
  if      (sstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s line %" PRId64 "):\n%s\n",
					    sqfp->filename, sqfp->linenumber, sqfp->errbuf);     
  else if (sstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    sstatus, sqfp->filename);

  for (i = 0; i < ncpus; ++i)
    {
      p7_bg_Destroy(info[i].bg);
    }

#ifdef HMMER_THREADS
  esl_workqueue_Reset(queue);
  while ((block = esl_workqueue_Remove(queue)) != NULL)
    {
      p7_oprofile_DestroyBlock(block);
    }
  esl_workqueue_Destroy(queue);
  esl_threads_Destroy(threadObj);
#endif

  free(info);

  esl_sq_Destroy(qsq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);
  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);
  esl_getopts_Destroy(go);
  return eslOK;

 ERROR:
  return eslFAIL;
}

#ifdef HMMER_THREADS
static int
threadedLoop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, P7_HMMFILE *hfp)
{
  int  status   = eslOK;
  int  sstatus  = eslOK;
  int  eofCount = 0;
  P7_OM_BLOCK   *block;
  ESL_ALPHABET  *abc = NULL;
  void          *newBlock;

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
pipelineThread(void *arg)
{
  int i;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  P7_OM_BLOCK   *block;
  void          *newBlock;
  
  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  /* loop until all blocks have been processed */
  block = (P7_OM_BLOCK *) newBlock;
  while (block->count > 0)
    {
      /* Main loop: */
      for (i = 0; i < block->count; ++i)
	{
	  P7_OPROFILE *om = block->list[i];

	  p7_pli_NewModel(info->pli, om, info->bg);
	  p7_oprofile_ReconfigLength(om, info->qsq->n);

	  p7_Pipeline(info->pli, om, info->bg, info->qsq, info->th);

	  p7_oprofile_Destroy(om);
	  p7_pipeline_Reuse(info->pli);

	  block->list[i] = NULL;
	}

      status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      block = (P7_OM_BLOCK *) newBlock;
    }

  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
  return;
}
#else
static int
serialLoop(WORKER_INFO *info, P7_HMMFILE *hfp)
{
  int            status;

  P7_OPROFILE   *om;
  ESL_ALPHABET  *abc = NULL;

  /* Main loop: */
  while ((status = p7_oprofile_ReadMSV(hfp, &abc, &om)) == eslOK)
    {
      p7_pli_NewModel(info->pli, om, info->bg);
      p7_oprofile_ReconfigLength(om, info->qsq->n);

      p7_Pipeline(info->pli, om, info->bg, info->qsq, info->th);

      p7_oprofile_Destroy(om);
      p7_pipeline_Reuse(info->pli);
    }

  return status;
}
#endif


/*****************************************************************
 * @LICENSE@
 *****************************************************************/

