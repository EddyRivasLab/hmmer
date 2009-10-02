/* hmmsearch: search profile HMM(s) against a sequence database.
 * 
 * SRE, Thu Dec 20 07:07:25 2007 [Janelia]
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

#ifdef HAVE_MPI
#include "mpi.h"
#include "esl_mpi.h"
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
  /* name           type      default  env  range     toggles  reqs   incomp  help   docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL, "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL, "direct output to file <f>, not stdout",                        7 },
  { "-A",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL, "save multiple alignment of all hits to file <s>",              7 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL, "save parseable table of per-sequence hits to file <s>",        7 },
  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL, "save parseable table of per-domain hits to file <s>",          7 },
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
  /* Control of inclusion thresholds */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include sequences <= this E-value threshold in output ali",    3 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include sequences >= this score threshold in output ali",      3 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include domains <= this E-value threshold in output ali",      3 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, "x>0",     NULL,  "-A",  "--inc_ga,--inc_nc,--inc_tc",        "include domains >= this score threshold in output ali",        3 },
  { "--inc_ga",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  "-A",  "--incE,--incT,--incdomE,--incdomT", "use profile's GA gathering cutoffs to set --incT, --incdomT",  3 },
  { "--inc_nc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  "-A",  "--incE,--incT,--incdomE,--incdomT", "use profile's NC noise cutoffs to set --incT, --incdomT",      3 },
  { "--inc_tc",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,  "-A",  "--incE,--incT,--incdomE,--incdomT", "use profile's TC trusted cutoffs to set --incT, --incdomT",    3 },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",      4 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,      NULL,  NULL, "--max",          "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             4 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,      NULL,  NULL, "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             4 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,      NULL,  NULL, "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             4 },
  { "--nobias",     eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL, "--max",          "turn off composition bias filter",                             4 },
  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,      NULL,  NULL,    NULL,          "turn off biased composition score corrections",                4 },
/* Other options */
  { "--acc",        eslARG_NONE,  FALSE,  NULL, NULL,      NULL,  NULL,    NULL,          "output target accessions instead of names if possible",        6 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",    NULL,  NULL,    NULL,          "set RNG seed to <n> (if 0: one-time arbitrary seed)",          6 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",  NULL,  NULL,  "--notextw",     "set max width of ASCII text output lines",                     6 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,      NULL,  NULL,  "--textw",       "unlimit ASCII text output line width",                         6 },
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,      NULL,  NULL,   NULL,           "assert target <seqfile> is in format <s>>: no autodetection",  6 },
#ifdef HMMER_THREADS
  { "--cpu",        eslARG_INT,     NULL,"HMMER_NCPU", "n>0", NULL,  NULL,  NULL,            "number of worker threads",                                     6 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[options] <query hmmfile> <target seqfile>";
static char banner[] = "search profile HMM(s) against a sequence database";



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

      puts("\nOptions controlling significance thresholds for reporting:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      puts("\nOptions controlling significance thresholds for inclusion in output alignment:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 

      puts("\nOptions controlling acceleration heuristics:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      puts("\nOther expert options:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                  != 2)     { puts("Incorrect number of command line arguments.");      goto ERROR; }
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
  
  fprintf(ofp, "# query HMM file:                  %s\n", hmmfile);
  fprintf(ofp, "# target sequence database:        %s\n", seqfile);
  if (esl_opt_IsUsed(go, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",      esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "-A"))          fprintf(ofp, "# MSA of all hits saved to file:   %s\n",      esl_opt_GetString(go, "-A"));
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
  if (esl_opt_IsUsed(go, "--nobias"))    fprintf(ofp, "# biased composition HMM filter:   off\n");
  if (esl_opt_IsUsed(go, "--nonull2"))   fprintf(ofp, "# null2 bias corrections:          off\n");
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0) fprintf(ofp, "# random number seed:              one-time arbitrary\n");
    else                                       fprintf(ofp, "# random number seed set to:       %d\n", esl_opt_GetInteger(go, "--seed"));
  }
  if (esl_opt_IsUsed(go, "--textw"))     fprintf(ofp, "# max ASCII text line length:      %d\n", esl_opt_GetInteger(go, "--textw"));
  if (esl_opt_IsUsed(go, "--notextw"))   fprintf(ofp, "# max ASCII text line length:      unlimited\n");
  if (esl_opt_IsUsed(go, "--tformat"))   fprintf(ofp, "# targ <seqfile> format asserted:  %s\n", esl_opt_GetString(go, "--tformat"));
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go	    = NULL;              /* command line processing                         */
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;		 /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;		 /* output stream for tabular per-seq (--domtblout) */
  char            *hmmfile  = NULL;              /* query HMM file                                  */
  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  char            *dbfile   = NULL;              /* target sequence database file                   */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  ESL_STOPWATCH   *w        = NULL;              /* for clocking performance                        */
  int              nquery   = 0;
  int              textw;
  int              status   = eslOK;
  int              hstatus  = eslOK;
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
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */
  process_commandline(argc, argv, &go, &hmmfile, &dbfile);    
  w = esl_stopwatch_Create();
  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  if (esl_opt_IsOn(go, "--tformat")) {
    dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Open the target sequence database */
  status = esl_sqfile_Open(dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          dbfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            dbfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, dbfile);  

  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open hmm file %s for reading.\n",                      hmmfile);
  else if (status == eslEFORMAT)   p7_Fail("Unrecognized format, trying to open hmm file %s for reading.\n", hmmfile);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening hmm file %s.\n", status,          hmmfile);  

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp")); }
  if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblfp")); }

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

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, hmmfile, dbfile);

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
    }

  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hstatus == eslOK) 
    {
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */

      nquery++;
      esl_stopwatch_Start(w);

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1)
	{
	  if (! esl_sqfile_IsRewindable(dbfp)) 
	    esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", dbfile);
	  esl_sqfile_Position(dbfp, 0);
	}

      fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M);
      if (hmm->acc  != NULL) fprintf(ofp, "Accession:   %s\n", hmm->acc);
      if (hmm->desc != NULL) fprintf(ofp, "Description: %s\n", hmm->desc);

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, info->bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
      p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */

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

      if (tblfp)    p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, info->th, info->pli, (nquery == 1));
      if (domtblfp) p7_tophits_TabularDomains(domtblfp, hmm->name, hmm->acc, info->th, info->pli, (nquery == 1));
  
      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, info->pli, w);
      fprintf(ofp, "//\n");

      /* Output the results in an MSA (-A option) */
      if (afp) {
	ESL_MSA *msa = NULL;

	if (p7_tophits_Alignment(info->th, abc, NULL, NULL, 0, p7_DEFAULT, &msa) == eslOK)
	  {
	    if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
	    else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);
	  
	    fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
	  } 
	else fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n");
	  
	esl_msa_Destroy(msa);
      }

      p7_pipeline_Destroy(info->pli);
      p7_tophits_Destroy(info->th);
      p7_oprofile_Destroy(info->om);
      p7_oprofile_Destroy(om);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);

      hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
    } /* end outer loop over query HMMs */

  switch(hstatus)
    {
    case eslEOD:
      p7_Fail("read failed, HMM file %s may be truncated?", hmmfile);
      break;
    case eslEFORMAT:
      p7_Fail("bad file format in HMM file %s", hmmfile);
      break;
    case eslEINCOMPAT:
      p7_Fail("HMM file %s contains different alphabets", hmmfile);
      break;
    case eslEOF:
      /* do nothing */
      break;
    default:
      p7_Fail("Unexpected error (%d) in reading HMMs from %s", hstatus, hmmfile);
    }

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

  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);
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

