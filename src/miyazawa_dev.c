/* miyazawa: partition function on exponential scores
 * 
 * In this file, I remove all I don't need for a simple
 * implementation of Miyazawa's partition function
 *
 * SC, Fri Feb 12 11:47:01 EST 2010 [Janelia] [Come And Find Me]
 * SVN $Id: phmmer.c 3143 2010-01-30 15:59:57Z eddys $
 */

/*
 * gcc -o miyazawa_dev miyazawa_dev.c -g -W -Wall -Wstrict-prototypes -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -ansi -pedantic -std=c99 -O2 -I. -L. -lhmmer -I./impl -L./impl -lhmmerimpl -I../easel -L../easel -leasel -lm
 *
 * Note the linking order of the libraries (they call each other from left to right)
 *
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

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles   reqs   incomp                             help                                       docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,              "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,      NULL,  NULL,  NULL,              "direct output to file <f>, not stdout",                        2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,      NULL,  NULL, "--textw",          "unlimit ASCII text output line width",                         2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",  NULL,  NULL, "--notextw",        "set max width of ASCII text output lines",                     2 },
/* Control of scoring system */
  { "--popen",      eslARG_REAL,  "0.02", NULL, "0<=x<0.5",NULL,  NULL,  NULL,              "gap open probability",                                         3 },
  { "--pextend",    eslARG_REAL,   "0.4", NULL, "0<=x<1",  NULL,  NULL,  NULL,              "gap extend probability",                                       3 },
  { "--mxfile",     eslARG_INFILE,  NULL, NULL, NULL,      NULL,  NULL,  NULL,              "substitution score matrix [default: BLOSUM62]",                3 },
  { "--lambda",     eslARG_REAL, "0.6931", NULL, NULL, NULL,  NULL,  NULL,                 "base of the score matrix [default: BLOSUM62]",                  3 }, /* WHAT IS THE BASE OF BLOSUM62, BITS OR HALF-BITS? */
  /* other options */
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,      NULL,  NULL,  NULL,              "assert query <seqfile> is in format <s>: no autodetection",   12 },
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,      NULL,  NULL,  NULL,              "assert target <seqdb> is in format <s>>: no autodetection",   12 },
  /* Control of reporting thresholds */
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <query seqfile> <target seqdb>";
static char banner[] = "partition function";

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * We use it to provide configuration options to serial_master. serial_loop gets the
 * WORKER_INFO object. Not sure I need to do it this way.
 *
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char            *qfile;             /* query sequence file                             */
  char            *dbfile;            /* database file                                   */

  int              do_mpi;            /* TRUE if we're doing MPI parallelization         */
  int              nproc;             /* how many MPI processes, total                   */
  int              my_rank;           /* who am I, in 0..nproc-1                         */
};

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (ESL_ALPHABET *abc, double popen, double pextend, double lambda, char* mxfile, ESL_SQFILE *dbfp);        /* WORKER_INFO SHOULD HAVE SEQUENCE AND SCORE SYSTEM */

/* process_commandline()
 * Take argc, argv, and options; parse the command line;
 * display help/usage info.
 */
static void 
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_qfile, char **ret_dbfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     esl_fatal("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n",  go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n",  go->errbuf); goto ERROR; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/

      puts("\noptions directing output:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      puts("\noptions controlling scoring system:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 

      puts("\noptions controlling reporting thresholds:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      puts("\noptions controlling inclusion (significance) thresholds:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      puts("\noptions controlling acceleration heuristics:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      puts("\noptions controlling E value calibration:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80); 

      puts("\nother expert options:");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
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

/* HEADER */

/* We call in serial_master just before starting the computation */

/* Only prints options set to a non-default value */

static int
output_header(FILE *ofp, ESL_GETOPTS *go, char *qfile, char *dbfile)
{
  p7_banner(ofp, go->argv[0], banner);  /* It includes some standard hmmer message along the banner defined here */
  
  fprintf(ofp, "# query sequence file:             %s\n", qfile);
  fprintf(ofp, "# target sequence database:        %s\n", dbfile);
  if (esl_opt_IsUsed(go, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",      esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "--notextw"))   fprintf(ofp, "# max ASCII text line length:      unlimited\n");
  if (esl_opt_IsUsed(go, "--textw"))     fprintf(ofp, "# max ASCII text line length:      %d\n",             esl_opt_GetInteger(go, "--textw"));  
  if (esl_opt_IsUsed(go, "--popen"))     fprintf(ofp, "# gap open probability:            %f\n",             esl_opt_GetReal  (go, "--popen"));
  if (esl_opt_IsUsed(go, "--pextend"))   fprintf(ofp, "# gap extend probability:          %f\n",             esl_opt_GetReal  (go, "--pextend"));
  if (esl_opt_IsUsed(go, "--mxfile"))    fprintf(ofp, "# subst score matrix:              %s\n",             esl_opt_GetString(go, "--mxfile"));
  if (esl_opt_IsUsed(go, "--lambda"))    fprintf(ofp, "# lambda:              %f\n",             esl_opt_GetString(go, "--lambda"));
  if (esl_opt_IsUsed(go, "--qformat"))   fprintf(ofp, "# query <seqfile> format asserted: %s\n",     esl_opt_GetString(go, "--qformat"));
  if (esl_opt_IsUsed(go, "--tformat"))   fprintf(ofp, "# target <seqdb> format asserted:  %s\n",     esl_opt_GetString(go, "--tformat"));

  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}

/********/
/* MAIN */
/********/

int
main(int argc, char **argv)
{
  int              status   = eslOK;

  ESL_GETOPTS     *go  = NULL;	/* command line processing                 */
  struct cfg_s     cfg;         /* configuration data                      */

  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet) 
   */
  cfg.qfile      = NULL;
  cfg.dbfile     = NULL;

  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */

  /* PROCESS COMMAND-LINE  */

  process_commandline(argc, argv, &go, &cfg.qfile, &cfg.dbfile);    

  /* STATUS IS SERIAL */

  status = serial_master(go, &cfg);

  esl_getopts_Destroy(go);

  return status;
}

/* serial_master()
 * The serial version of hmmsearch.
 * For each query HMM in <hmmfile> search the database for hits.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with p7_Fail().
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;             /* output file for results (default stdout)         */
  int              qformat  = eslSQFILE_UNKNOWN;  /* format of qfile (it will autodetect the format unless it is asserted */
  ESL_SQFILE      *qfp      = NULL;		 		  /* open qfile                                       */
  ESL_SQ          *qsq      = NULL;               /* query sequence                                   */
  int              dbformat = eslSQFILE_UNKNOWN;  /* format of dbfile                                 */
  ESL_SQFILE      *dbfp     = NULL;               /* open dbfile                                      */
  ESL_ALPHABET    *abc      = NULL;               /* sequence alphabet                                */
  ESL_STOPWATCH   *w        = NULL;               /* for timing                                       */
  int              nquery   = 0;
  int              seed;
  int              textw;                        /* set max width of ASCII text output lines */
  int              status   = eslOK;
  int              qstatus  = eslOK;             /* WHAT IS Q? */
  int              sstatus  = eslOK;             /* WHAT IS S? */
  int              i;

  double           popen;
  double           pextend;
  double           lambda;
  char             *mxfile;

  int              ncpus    = 0;

  int              infocnt  = 0;                 /* WHAT IS THIS? */
 // WORKER_INFO     *info     = NULL;

  /* Initializations */
  abc     = esl_alphabet_Create(eslAMINO);      /* returns a pointer to a ESL_ALPHABET object */
  w       = esl_stopwatch_Create();
  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* START WATCH */
  esl_stopwatch_Start(w);

  /* INPUT FORMATS */
  /* Query */

  if (esl_opt_IsOn(go, "--qformat")) {
    qformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (qformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }

  /* Target */

  if (esl_opt_IsOn(go, "--tformat")) {
    dbformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat")); /* target format */
  }

  /* SCORE SYSTEM */

  popen = esl_opt_GetReal(go, "--popen");
  pextend = esl_opt_GetReal(go, "--pextend");
  lambda = esl_opt_GetReal(go, "--lambda");
  mxfile = esl_opt_GetString(go, "--mxfile");

  /* OPEN OUTPUT FILES */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"),          "w")) == NULL)  esl_fatal("Failed to open output file %s for writing\n",                 esl_opt_GetString(go, "-o")); } 
  /*
  if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"),          "w")) == NULL)  esl_fatal("Failed to open alignment output file %s for writing\n",       esl_opt_GetString(go, "-A")); } 
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp")); }
  if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblfp")); }
   */

  /* OPEN TARGET FILE (autodetecting format unless given) */ /* contrary to documentation, if I open digital, I need to read digital. Otherwise, it compiles but get no sequence */
  status =  esl_sqfile_OpenDigital(abc, cfg->dbfile, dbformat, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open target sequence database %s for reading\n",      cfg->dbfile);
  else if (status == eslEFORMAT)   esl_fatal("Target sequence database file %s is empty or misformatted\n",   cfg->dbfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        esl_fatal("Unexpected error %d opening target sequence database file %s\n", status, cfg->dbfile);

  /* OPEN QUERY FILE (autodetecting format unless given) */
  status = esl_sqfile_OpenDigital(abc, cfg->qfile, qformat, NULL, &qfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open sequence file %s for reading\n",      cfg->qfile);
  else if (status == eslEFORMAT)   esl_fatal("Sequence file %s is empty or misformatted\n",        cfg->qfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        esl_fatal ("Unexpected error %d opening sequence file %s\n", status, cfg->qfile);

  /* CREATE DIGITAL QUERY SEQUENCE */
  qsq  = esl_sq_CreateDigital(abc);

//  infocnt = (ncpus == 0) ? 1 : ncpus;        /* not sure what this is doing here! */
//  ESL_ALLOC(info, sizeof(*info) * infocnt);  /* not sure what this is doing here! */

  /* HEADER (INPUT FILES AND OTHER OPTIONS) */
  output_header(ofp, go, cfg->qfile, cfg->dbfile);

  /*
  for (i = 0; i < infocnt; ++i)
    {
      info[i].pli   = NULL;
      info[i].th    = NULL;
      info[i].om    = NULL;
      info[i].bg    = p7_bg_Create(abc);
    }
*/

  /* OUTER LOOP OVER SEQUENCE QUERIES */
  while ((qstatus = esl_sqio_Read(qfp, qsq)) == eslOK) /* qsq is set for digital, but I am not sure at which point the conversion is made */
    {
      nquery++;
      if (qsq->n == 0) continue; /* skip zero length seqs as if they aren't even present */

      /* START WATCH */
      esl_stopwatch_Start(w);

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1)
	{
	  if (! esl_sqfile_IsRewindable(dbfp)) esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile); /* NOT SURE WHAT IS GOING ON HERE!!! */
	  esl_sqfile_Position(dbfp, 0);
	}

      fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n);
      if (qsq->acc[0]  != '\0') fprintf(ofp, "Accession:   %s\n", qsq->acc);
      if (qsq->desc[0] != '\0') fprintf(ofp, "Description: %s\n", qsq->desc);  

      /*
      for (i = 0; i < infocnt; ++i)
	{*/
	  /* Create processing pipeline and hit list */
	 /* info[i].th  = p7_tophits_Create();
	  info[i].om  = p7_oprofile_Clone(om);
	  info[i].pli = p7_pipeline_Create(go, om->M, 100, p7_SEARCH_SEQS); */ /* L_hint = 100 is just a dummy for now */
	  /*p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);
	}*/

      /* INNER LOOP OVER SEQUENCE TARGETS */

      sstatus = serial_loop(abc, popen, pextend, lambda, mxfile, dbfp);

      switch(sstatus)
	{
	case eslEFORMAT: 
	  esl_fatal("Parse failed (sequence file %s):\n%s\n",
		    dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
	  break;
	case eslEOF:
	  /* do nothing */
	  break;
	default:
	  esl_fatal("Unexpected error %d reading sequence file %s",
		    sstatus, dbfp->filename);
	}


      /* merge the results of the search results */
/*      for (i = 1; i < infocnt; ++i)
	{
	  p7_tophits_Merge(info[0].th, info[i].th);
	  p7_pipeline_Merge(info[0].pli, info[i].pli);

	  p7_pipeline_Destroy(info[i].pli);
	  p7_tophits_Destroy(info[i].th);
	  p7_oprofile_Destroy(info[i].om);
	}
*/
      /* Print the results.  */
/*      p7_tophits_Sort(info->th);
      p7_tophits_Threshold(info->th, info->pli);
      p7_tophits_Targets(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
      p7_tophits_Domains(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
  
      if (tblfp)    p7_tophits_TabularTargets(tblfp,    qsq->name, qsq->acc, info->th, info->pli, (nquery == 1));
      if (domtblfp) p7_tophits_TabularDomains(domtblfp, qsq->name, qsq->acc, info->th, info->pli, (nquery == 1));
*/
	  /* STOP WATCH */
      esl_stopwatch_Stop(w);
/*      p7_pli_Statistics(ofp, info->pli, w); */
      fprintf(ofp, "//\n");

/*      p7_tophits_Destroy(info->th);
      p7_pipeline_Destroy(info->pli);
      p7_oprofile_Destroy(info->om);
      p7_oprofile_Destroy(om);
*/
      esl_sq_Reuse(qsq);
    } /* end outer loop over query sequences */

  if      (qstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					    qfp->filename, esl_sqfile_GetErrorBuf(qfp));
  else if (qstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    qstatus, qfp->filename);

/*  for (i = 0; i < infocnt; ++i)
    {
      p7_bg_Destroy(info[i].bg);
    }
*/

//  free(info);

  esl_sqfile_Close(dbfp);
  esl_sqfile_Close(qfp);
  esl_stopwatch_Destroy(w);
  esl_sq_Destroy(qsq);
  esl_alphabet_Destroy(abc);

  if (ofp      != stdout) fclose(ofp);
  return eslOK;

 ERROR:
  return eslFAIL;
}

/************************************/
/* OUTER LOOP OVER SEQUENCE TARGETS */
/************************************/

static int
serial_loop(ESL_ALPHABET *abc, double popen, double pextend, double lambda, char *mxfile, ESL_SQFILE *dbfp)
{
  int      sstatus;
  ESL_SQ   *dbsq     = NULL;   /* one target sequence (digital)  */

  dbsq = esl_sq_CreateDigital(abc);  /* Going digital here! deal with it! */

  printf("popen: %f\n", popen);
  printf("pextend: %f\n", pextend);
  printf("lambda: %f\n", lambda);
  printf("mxfile: %s\n", mxfile);

  /* Main loop: */
  while ((sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
    {

	  /*
	   * HERE I WANT THE ZCORE COMPUTATION -need to call a zcore() function-
	   *
	   * zscore()
	   *
	   */

	  /*
      p7_pli_NewSeq(info->pli, dbsq);
      p7_bg_SetLength(info->bg, dbsq->n);
      p7_oprofile_ReconfigLength(info->om, dbsq->n);
      
      p7_Pipeline(info->pli, info->om, info->bg, dbsq, info->th);
	  
      esl_sq_Reuse(dbsq);
      p7_pipeline_Reuse(info->pli);
      */
    }

  esl_sq_Destroy(dbsq);

  return sstatus;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
