/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* hmmpfam.c
 * SRE, Mon Aug 25 17:03:14 1997 [Denver] 
 *
 * Search a single sequence against an HMM database.
 * Conditionally includes PVM parallelization when HMMER_PVM is defined
 *    at compile time; hmmpfam --pvm runs the PVM version.
 *    
 * RCS $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#ifdef HMMER_THREADS
#include <pthread.h>
#endif
#ifdef HMMER_PVM
#include <pvm3.h>
#endif

#include "squid.h"		/* general sequence analysis library    */
#include "config.h"		/* compile-time configuration constants */
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "version.h"		/* version info                         */

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static char banner[] = "hmmpfam - search a single seq against HMM database";

static char usage[]  = "\
Usage: hmmpfam [-options] <hmm database> <sequence file>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
   -n        : nucleic acid models/sequence (default protein)\n\
   -A <n>    : sets alignment output limit to <n> best domain alignments\n\
   -E <x>    : sets E value cutoff (globE) to <x>\n\
   -T <x>    : sets T bit threshold (globT) to <x>\n\
   -Z <n>    : sets Z (# seqs) for E-value calculation\n\
";

static char experts[] = "\
   --cpu <n> : run <n> threads in parallel (if threaded)\n\
   --domE <x>: sets domain Eval cutoff (2nd threshold) to <x>\n\
   --domT <x>: sets domain T bit thresh (2nd threshold) to <x>\n\
   --forward : use the full Forward() algorithm instead of Viterbi\n\
   --null2   : turn OFF the post hoc second null model\n\
   --pvm     : run on a PVM (Parallel Virtual Machine) cluster\n\
   --xnu     : turn ON XNU filtering of query protein sequence\n\
\n";


static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-n",        TRUE,  sqdARG_NONE },
  { "-A",        TRUE,  sqdARG_INT  },  
  { "-E",        TRUE,  sqdARG_FLOAT},  
  { "-T",        TRUE,  sqdARG_FLOAT},  
  { "-Z",        TRUE,  sqdARG_INT },
  { "--cpu",     FALSE, sqdARG_INT },
  { "--domE",    FALSE, sqdARG_FLOAT},
  { "--domT",    FALSE, sqdARG_FLOAT},
  { "--forward", FALSE, sqdARG_NONE },
  { "--null2",   FALSE, sqdARG_NONE },  
  { "--pvm",     FALSE, sqdARG_NONE },  
  { "--xnu",     FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

#ifdef HMMER_THREADS
/* POSIX threads version:
 * the threads share a workpool_s structure amongst themselves,
 * for obtaining locks on input HMM file and output histogram and
 * tophits structures.
 */
struct workpool_s {
  /* Shared configuration resources which don't change:
   */
  char  *hmmfile;               /* name of HMM file                */
  char  *dsq;                   /* digitized query sequence        */
  char  *seqname;               /* sequence name                   */
  int    L;			/* length of dsq                   */
  int    do_forward;		/* TRUE to score using Forward     */
  int    do_null2;		/* TRUE to apply null2 ad hoc correction */
  float  globT;                 /* global score threshold          */
  double globE;			/* global E-value threshold        */
  int    Z;                     /* effective # of seqs in database */
  
  /* Shared (mutex-protected) input resources:
   */
  HMMFILE *hmmfp;               /* ptr to open HMM file           */
  int nhmm;			/* number of HMMs searched so far */
  pthread_mutex_t input_lock;   /* mutex for locking input        */

  /* Shared (mutex-protected) output resources:
   */
  struct tophit_s *ghit;        /* per-sequence top hits */
  struct tophit_s *dhit;        /* per-domain top hits */
  pthread_mutex_t output_lock;  /* mutex for locking output */

  /* Thread pool information
   */
  pthread_t *thread;            /* our pool of threads */
  int        num_threads;       /* number of threads   */
};

static struct workpool_s *workpool_start(char *hmmfile, HMMFILE *hmmfp, char *dsq, char *seqname, int L,
					 int do_forward, int do_null2, float globT, double globE, int Z, 
					 struct tophit_s *ghit, struct tophit_s *dhit, 
					 int num_threads);
static void  workpool_stop(struct workpool_s *wpool);
static void  workpool_free(struct workpool_s *wpool);
static void *worker_thread(void *ptr);
#endif /* HMMER_THREADS */


#ifdef HMMER_PVM
static void main_loop_pvm(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo, 
			  float globT, double globE, int Z, int do_xnu, int do_forward, int do_null2,
			  struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm);
static void init_slaves(int *slave_tid, int nslaves, char *hmmfile, char *seq, int L,
			float globT, double globE, int Z, 
			int do_xnu, int do_forward, int do_null2);
#endif
static void main_loop_serial(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo, 
			     float globT, double globE, int Z, int do_xnu, int do_forward, int do_null2,
			     int num_threads,
			     struct tophit_s *ghit, struct tophit_s *dhit, int *nhmm);
static void record_domains(struct tophit_s *h, 
			   struct plan7_s *hmm, char *dsq, int L, char *seqname,
			   struct p7trace_s *tr, double whole_pval, float whole_sc,
			   int do_null2);

int
main(int argc, char **argv) 
{
  char            *hmmfile;	/* file to read HMMs from                  */
  HMMFILE         *hmmfp;       /* opened hmmfile for reading              */
  char            *seqfile;     /* file to read target sequence from       */ 
  SQFILE          *sqfp;        /* opened seqfile for reading              */
  int              format;	/* format of seqfile                       */
  char              *seq;	/* target sequence                         */
  SQINFO             sqinfo;	/* optional info for seq                   */
  struct fancyali_s *ali;	/* an alignment for display                */
  struct tophit_s   *ghit;      /* list of top hits and alignments for seq  */
  struct tophit_s   *dhit;	/* list of top hits/alignments for domains  */

  float   sc;			/* log-odds score in bits                  */
  double  pvalue;		/* pvalue of an HMM score                  */
  double  evalue;		/* evalue of an HMM score                  */
  double  motherp;		/* pvalue of a whole seq HMM score         */
  float   mothersc;		/* score of a whole seq parent of domain   */
  int     sqfrom, sqto;		/* coordinates in sequence                 */
  int     hmmfrom, hmmto;	/* coordinate in HMM                       */
  char   *name, *desc;          /* hit HMM name and description            */
  int     hmmlen;		/* length of HMM hit                       */
  int     nhmm;			/* number of HMMs searched                 */
  int     domidx;		/* number of this domain                   */
  int     ndom;			/* total # of domains in this seq          */
  int     namewidth;		/* max width of sequence name              */

  int    Alimit;		/* A parameter limiting output alignments   */
  float  globT;			/* T parameter: keep only hits > globT bits */
  double globE;			/* E parameter: keep hits < globE E-value   */
  float  domT;			/* T parameter for individual domains       */
  double domE;			/* E parameter for individual domains       */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  int   do_forward;		/* TRUE to use Forward() not Viterbi()      */
  int   do_nucleic;		/* TRUE to do DNA/RNA instead of protein    */
  int   do_null2;		/* TRUE to adjust scores with null model #2 */
  int   do_pvm;			/* TRUE to run on PVM                       */
  int   do_xnu;			/* TRUE to do XNU filtering                 */
  int   Z;			/* nseq to base E value calculation on      */
  int   i;

  int   num_threads;            /* number of worker threads */   

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  do_forward  = FALSE;
  do_nucleic  = FALSE;
  do_null2    = TRUE;
  do_pvm      = FALSE;
  do_xnu      = FALSE;
  Z           = 59021;		/* default: nseq in Swissprot34     */
  
  Alimit      = INT_MAX;	/* no limit on alignment output     */
  globE       = 10.0;		/* use a reasonable Eval threshold; */
  globT       = -FLT_MAX;	/*   but no bit threshold,          */
  domT        = -FLT_MAX;	/*   no domain bit threshold,       */
  domE        = FLT_MAX;        /*   and no domain Eval threshold.  */
#ifdef HMMER_THREADS
  num_threads = ThreadNumber(); /* only matters if we're threaded */
#else
  num_threads = 0;
#endif 

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n")        == 0) do_nucleic = TRUE; 
    else if (strcmp(optname, "-A")        == 0) Alimit     = atoi(optarg);  
    else if (strcmp(optname, "-E")        == 0) globE      = atof(optarg);
    else if (strcmp(optname, "-T")        == 0) globT      = atof(optarg);
    else if (strcmp(optname, "-Z")        == 0) Z          = atoi(optarg);
    else if (strcmp(optname, "--cpu")     == 0) num_threads= atoi(optarg);
    else if (strcmp(optname, "--domE")    == 0) domE       = atof(optarg);
    else if (strcmp(optname, "--domT")    == 0) domT       = atof(optarg);
    else if (strcmp(optname, "--forward") == 0) do_forward = TRUE;
    else if (strcmp(optname, "--null2")   == 0) do_null2   = FALSE;
    else if (strcmp(optname, "--pvm")     == 0) do_pvm     = TRUE;
    else if (strcmp(optname, "--xnu")     == 0) do_xnu     = TRUE;
    else if (strcmp(optname, "-h")      == 0) {
      Banner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }
  if (argc - optind != 2)
    Die("Incorrect number of arguments.\n%s\n", usage);

  hmmfile = argv[optind++];
  seqfile = argv[optind++]; 

#ifndef HMMER_PVM
  if (do_pvm) Die("PVM support is not compiled into HMMER; --pvm doesn't work.");
#endif
#ifndef HMMER_THREADS
  if (num_threads) Die("Posix threads support is not compiled into HMMER; --cpu doesn't have any effect");
#endif

  /*********************************************** 
   * Open sequence database (must be in curr directory);
   * get target sequence.
   ***********************************************/

  if (do_nucleic) SetAlphabet(hmmNUCLEIC);
  else            SetAlphabet(hmmAMINO);

  if (! SeqfileFormat(seqfile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: 
      Die("Sequence file %s could not be opened for reading", seqfile); 
      break;
    case SQERR_FORMAT: 
    default:           
      Die("Failed to determine format of sequence file %s", seqfile);
    }
  if ((sqfp = SeqfileOpen(seqfile, format, NULL)) == NULL)
    Die("Failed to open sequence file %s\n%s\n", seqfile, usage);

  /*********************************************** 
   * Open HMM database (might be in HMMERDB or current directory)
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM database %s\n%s", hmmfile, usage);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  Banner(stdout, banner);
  printf(   "HMM file:                 %s\n", hmmfile);
  printf(   "Sequence file:            %s\n", seqfile);
  if (do_pvm)
    printf( "PVM:                      ACTIVE\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  /*********************************************** 
   * Search each HMM against each sequence
   ***********************************************/

  while (ReadSeq(sqfp, format, &seq, &sqinfo)) 
    {
      ghit = AllocTophits(20);   /* keeps full seq scores */
      dhit = AllocTophits(20);   /* keeps domain scores   */

      /* 1. Search sequence against all HMMs.
       *    Significant scores+alignments accumulate in ghit, dhit.
       */
      if (!do_pvm)
	main_loop_serial(hmmfile, hmmfp, seq, &sqinfo, 
			 globT, globE, Z, do_xnu, do_forward, do_null2, num_threads,
			 ghit, dhit, &nhmm);
#ifdef HMMER_PVM
      else if (do_pvm)
	{
	  SQD_DPRINTF1(("Entering PVM main loop\n"));
	  main_loop_pvm(hmmfile, hmmfp, seq, &sqinfo, 
			globT, globE, Z, do_xnu, do_forward, do_null2,
			ghit, dhit, &nhmm);
	}
#endif
      else Die("wait. that can't happen. I didn't do anything.");

      /* 2. (Done searching all HMMs for this query seq; start output)
       *    Report the overall sequence hits, sorted by significance.
       */
      printf("Query:  %s  %s\n", sqinfo.name, 
	     sqinfo.flags & SQINFO_DESC ? sqinfo.desc : "");
      FullSortTophits(ghit);
      namewidth =  MAX(8, TophitsMaxName(ghit));

      printf("\nScores for sequence family classification (score includes all domains):\n");
      printf("%-*s %-*s %7s %10s %3s\n", namewidth, "Sequence", 52-namewidth, "Description", "Score", "E-value", " N ");
      printf("%-*s %-*s %7s %10s %3s\n", namewidth, "--------", 52-namewidth, "-----------", "-----", "-------", "---");
      for (i = 0; i < ghit->num; i++)
	{
	  char *safedesc;
	  GetRankedHit(ghit, i, 
		       &pvalue, &sc, NULL, NULL,
		       &name, &desc,
		       NULL, NULL, NULL,           /* seq positions */
		       NULL, NULL, NULL,           /* HMM positions */
		       NULL, &ndom,                /* domain info   */
		       NULL);	                   /* alignment info*/

	  evalue = pvalue * (double) Z;

	  /* safedesc is a workaround for an apparent Linux printf()
	   * bug with the *.*s format. dbmalloc crashes with a memchr() ptr out of bounds
	   * flaw if the malloc'ed space for desc is short. The workaround
           * is to make sure the ptr for *.* has a big malloc space.
	   */
	  if (desc != NULL && strlen(desc) < 80) 
	    {
	      safedesc = MallocOrDie(sizeof(char) * 80);
	      strcpy(safedesc, desc);
	    }
	  else safedesc = Strdup(desc);

	  if (evalue < globE && sc >= globT) 
	    printf("%-*s %-*.*s %7.1f %10.2g %3d\n", 
		   namewidth, name, 
		   52-namewidth, 52-namewidth, safedesc != NULL ? safedesc : "",
		   sc, evalue, ndom);
	  else if (evalue >= globE)
	    {
	      if (i > 0) printf("\t[no more scores below E threshold]\n");
	      break;
	    }
	  else if (sc <= globT)
	    {
	      if (i > 0) printf("\t[no more scores above T threshold]");
	      break;
	    }
	  free(safedesc);
	}
      if (i == 0) printf("\t[no hits above thresholds]\n");

      /* 3. Report domain hits (sorted on sqto coordinate)
       */
      FullSortTophits(dhit);
      namewidth = MAX(8, TophitsMaxName(dhit));

      printf("\nParsed for domains:\n");
      printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
	     namewidth, "Sequence", "Domain ", "seq-f", "seq-t", "hmm-f", "hmm-t", "score", "E-value");
      printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
	     namewidth, "--------", "-------", "-----", "-----", "-----", "-----", "-----", "-------");
      
      for (i = 0; i < dhit->num; i++)
	{
	  GetRankedHit(dhit, i, 
		       &pvalue, &sc, &motherp, &mothersc,
		       &name, NULL,
		       &sqfrom, &sqto, NULL, 
		       &hmmfrom, &hmmto, &hmmlen, 
		       &domidx, &ndom,
		       NULL);
	  evalue = pvalue * (double) Z;
	  
	  if (motherp * (double) Z >= globE || mothersc <= globT) 
	    continue;
	  else if (evalue < domE && sc > domT)
	    printf("%-*s %3d/%-3d %5d %5d %c%c %5d %5d %c%c %7.1f %8.2g\n",
		   namewidth, name, 
		   domidx, ndom,
		   sqfrom, sqto, 
		   sqfrom == 1 ? '[' : '.', sqto == sqinfo.len ? ']' : '.',
		   hmmfrom, hmmto,
		   hmmfrom == 1 ? '[':'.',  hmmto == hmmlen ? ']' : '.',
		   sc, evalue);
	  else if (evalue >= domE) {
	    if (i > 0) printf("\t[no more scores below domE threshold]\n");
	    break;
	  }
	  else if (sc <= domT) {
	    if (i > 0) printf("\t[no more scores above domT threshold]\n");
	    break;
	  }
	}
      if (i == 0) printf("\t[no hits above thresholds]\n");      


      /* 3. Alignment output, also by domain.
       *    dhits is already sorted and namewidth is set, from above code.
       *    Number of displayed alignments is limited by Alimit parameter;
       *    also by domE (evalue threshold), domT (score theshold).
       */
      if (Alimit != 0)
	{
	  printf("\nAlignments of top-scoring domains:\n");
	  for (i = 0; i < dhit->num; i++)
	    {
	      if (i == Alimit) break; /* limit to Alimit output alignments */
	      GetRankedHit(dhit, i, 
			   &pvalue, &sc, &motherp, &mothersc,
			   &name, NULL,
			   &sqfrom, &sqto, NULL,              /* seq position info  */
			   &hmmfrom, &hmmto, &hmmlen,         /* HMM position info  */
			   &domidx, &ndom,                    /* domain info        */
			   &ali);	                      /* alignment info     */
	      evalue = pvalue * (double) Z;

	      if (motherp * (double) Z >= globE || mothersc <= globT) 
		continue;
	      else if (evalue < domE && sc > domT) 
		{
		  printf("%s: domain %d of %d, from %d to %d: score %.1f, E = %.2g\n", 
			 name, domidx, ndom, sqfrom, sqto, sc, evalue);
		  PrintFancyAli(stdout, ali);
		}
	      else if (evalue >= domE) {
		if (i > 0) printf("\t[no more alignments below E threshold]\n");
		break;
	      }
	      else if (sc <= domT) {
		if (i > 0) printf("\t[no more alignments above T threshold]\n");
		break;
	      }
	    }
	  if (i == 0)      printf("\t[no hits above thresholds]\n");
	  if (i == Alimit) printf("\t[output cut off at A = %d top alignments]\n", Alimit);
	}


      printf("//\n");
      FreeSequence(seq, &sqinfo); 
      FreeTophits(ghit);
      FreeTophits(dhit);

      HMMFileRewind(hmmfp);
    }

  /*********************************************** 
   * Clean-up and exit.
   ***********************************************/
  SeqfileClose(sqfp);
  HMMFileClose(hmmfp);
  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) malloc_list(2, histid1, histid2);
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  return 0;
}


/* Function: main_loop_serial()
 * Date:     SRE, Fri Aug  7 13:46:48 1998 [St. Louis]
 *
 * Purpose:  Search a sequence against an HMM database;
 *           main loop for the serial (non-PVM, non-threads)
 *           version. 
 *
 * Args:     hmmfile - name of HMM file
 *           hmmfp   - open HMM file (and at start of file)
 *           dsq     - digitized sequence
 *           sqinfo  - ptr to SQINFO optional info for dsq
 *           globT   - bit threshold for significant global score
 *           globE   - E-value threshold for significant global sc
 *           Z       - effective number of seqs to calc E with
 *           do_xnu     - TRUE to apply XNU filter to sequence
 *           do_forward - TRUE to use Forward() scores
 *           do_null2   - TRUE to adjust scores w/ ad hoc null2 model
 *           num_threads- number of threads, if threaded
 *           ghit    - global hits list
 *           dhit    - domain hits list
 *           ret_nhmm    - number of HMMs searched.
 *
 * Returns:  (void)
 */
static void
main_loop_serial(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo, 
		 float globT, double globE, int Z, int do_xnu, int do_forward, int do_null2,
		 int num_threads,
		 struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm) 
{
  char              *dsq;       /* digitized sequence                      */
  int     nhmm;			/* number of HMMs searched                 */
#ifdef HMMER_THREADS
  struct workpool_s *wpool;     /* pool of worker threads  */
#else
  struct plan7_s    *hmm;       /* current HMM to search with              */ 
  struct p7trace_s  *tr;	/* traceback of alignment                  */
  float   sc;                   /* an alignment score                      */ 
  double  pvalue;		/* pvalue of an HMM score                  */
#endif

  /* Prepare sequence.
   */
  dsq = DigitizeSequence(seq, sqinfo->len);
  if (do_xnu) XNU(dsq, sqinfo->len);

#ifdef HMMER_THREADS
  wpool = workpool_start(hmmfile, hmmfp, dsq, sqinfo->name, sqinfo->len, 
			 do_forward, do_null2, globT, globE, Z, 
			 ghit, dhit, num_threads);
  workpool_stop(wpool);
  nhmm = wpool->nhmm;
  workpool_free(wpool);
#else /* unthreaded code: */

  nhmm = 0;
  while (HMMFileRead(hmmfp, &hmm)) {
    if (hmm == NULL) 
      Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);
    P7Logoddsify(hmm, !(do_forward));

#ifdef MEMDEBUG
    printf("working on HMM %s, len %d.\n", hmm->name, hmm->M);
    printf("current allocation = %ld\n", malloc_inuse(NULL));
#endif

    /* Score sequence, do alignment (Viterbi), recover trace
     */
    if (P7ViterbiSize(sqinfo->len, hmm->M) <= RAMLIMIT)
      sc = P7Viterbi(dsq, sqinfo->len, hmm, &tr);
    else
      sc = P7SmallViterbi(dsq, sqinfo->len, hmm, &tr);
	
    if (do_forward) sc  = P7Forward(dsq, sqinfo->len, hmm, NULL);
    if (do_null2)   sc -= TraceScoreCorrection(hmm, tr, dsq);

    /* Store scores/pvalue for each HMM aligned to this sequence, overall
     */
    pvalue = PValue(hmm, sc);
    if (sc > globT && pvalue * (float) Z < globE) 
      { 
	RegisterHit(ghit, sc, pvalue, sc,
		    0., 0.,	                  /* no mother seq */
		    hmm->name, hmm->desc, 
		    0,0,0,                	  /* seq positions  */
		    0,0,0,	                  /* HMM positions  */
		    0, TraceDomainNumber(tr), /* domain info    */
		    NULL);	                  /* alignment info */
	/* 1c. Store scores/evalue/position/alignment for each HMM
	 *    aligned to domains of this sequence; UNFINISHED
	 */
	record_domains(dhit, hmm, dsq, sqinfo->len, sqinfo->name, tr, pvalue, sc, do_null2);
      }
    P7FreeTrace(tr);
    FreePlan7(hmm);
    nhmm++;
  }
#endif /* threaded vs. unthreaded code */

  free(dsq);
  *ret_nhmm = nhmm;
  return;
}

/* Function: record_domains()
 * Date:     SRE, Tue Oct 28 14:20:56 1997 [Cambridge UK]
 * 
 * Purpose:  Decompose a trace, and register scores, P-values, alignments,
 *           etc. for individual domains in a hitlist. Almost a duplicate
 *           of the eponymous function in hmmsearch, but we sort by
 *           start position of the domain. 
 *           
 * Return:   (void)          
 */
static void
record_domains(struct tophit_s *h, 
	       struct plan7_s *hmm, char *dsq, int L, char *seqname,
	       struct p7trace_s *tr, double whole_pval, float whole_sc,
	       int do_null2)
{
  struct p7trace_s **tarr;      /* array of per-domain traces */
  struct fancyali_s *ali;       /* alignment of a domain      */ 
  int ntr;			/* number of domain traces    */
  int idx;			/* index for traces           */
  int k1, k2;			/* start, stop coord in model */
  int i1, i2;			/* start, stop in sequence    */
  float  score;
  double pvalue;

  TraceDecompose(tr, &tarr, &ntr);
  if (ntr == 0) Die("TraceDecompose() screwup"); /* "can't happen" (!) */

  for (idx = 0; idx < ntr; idx++)
    {
      /* Get the score and bounds of the match.
       */
      score  = P7TraceScore(hmm, dsq, tarr[idx]);
      if (do_null2) 
	score -= TraceScoreCorrection(hmm, tarr[idx], dsq);
      pvalue = PValue(hmm, score); 
      TraceSimpleBounds(tarr[idx], &i1, &i2, &k1, &k2);
      
      /* Record the match
       */
      ali = CreateFancyAli(tarr[idx], hmm, dsq, seqname);
      RegisterHit(h, -1.*(double)i2,           /* sortkey= -(start) (so 1 comes first) */
		  pvalue, score, whole_pval, whole_sc,
		  hmm->name, hmm->desc, 
		  i1,i2, L, 
		  k1,k2, hmm->M, 
		  idx+1,ntr,
		  ali);
    }
  
  for (idx = 0; idx < ntr; idx++)
    P7FreeTrace(tarr[idx]);
  free(tarr);
  return;
}

#ifdef HMMER_PVM
/*****************************************************************
 * PVM specific functions
 ****************************************************************/ 

/* Function: main_loop_pvm()
 * Date:     SRE, Fri Aug  7 13:58:34 1998 [St. Louis]
 *
 * Purpose:  Search a sequence against an HMM database;
 *           main loop for the PVM version.
 *
 * Args:     hmmfile - name of HMM file
 *           hmmfp   - open HMM file (and at start of file)
 *           seq     - sequence to search against
 *           sqinfo  - ptr to SQINFO optional info for dsq
 *           globT   - bit threshold for significant global score
 *           globE   - E-value threshold for significant global sc
 *           Z       - effective number of seqs to calc E with
 *           do_xnu     - TRUE to apply XNU filter to sequence
 *           do_forward - TRUE to use Forward() scores
 *           do_null2   - TRUE to adjust scores w/ ad hoc null2 model
 *           ghit    - global hits list
 *           dhit    - domain hits list
 *           nhmm    - number of HMMs searched.
 *
 * Returns:  (void)
 */
static void
main_loop_pvm(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo, 
	      float globT, double globE, int Z, int do_xnu, int do_forward, int do_null2,
	      struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm)
{
  struct plan7_s   *hmm;        /* HMM that was searched with */
  struct p7trace_s *tr;         /* a traceback structure */
  char  *dsq;                   /* digitized sequence */ 
  float  sc;			/* score of an HMM match */
  int    master_tid;		/* master's ID */
  int   *slave_tid;             /* array of slave IDs */
  int   *hmmlist;               /* array of hmm indexes being worked on by slaves */
  int    nslaves;	       	/* number of slaves in virtual machine   */
  int    nhmm;			/* number of HMMs searched */
  int    slaveidx;		/* index of a slave wanting work         */
  int    slave, msg;
  int    sent_trace;		/* TRUE if slave sent us a trace */
  char   slavename[32];		/* name of HMM that slave actually did */
  double pvalue;		/* pvalue of HMM score */

  /* Sanity checks.
   */
  if (hmmfp->gsi == NULL)
    Die("HMM file %s needs a GSI index to use PVM. See: hmmindex.", hmmfile);

  /* Prepare sequence.
   */
  dsq = DigitizeSequence(seq, sqinfo->len);
  if (do_xnu) XNU(dsq, sqinfo->len);

  /* Initialize PVM
   */
  SQD_DPRINTF1(("Requesting master TID...\n"));
  master_tid = pvm_mytid();
  SQD_DPRINTF1(("Spawning slaves...\n"));
  PVMSpawnSlaves("hmmpfam-slave", &slave_tid, &nslaves);
  hmmlist   = MallocOrDie(sizeof(int) * nslaves);
  SQD_DPRINTF1(("Spawned a total of %d slaves...\n", nslaves));

  /* Initialize the slaves
   */
  init_slaves(slave_tid, nslaves, hmmfile, seq, sqinfo->len,
	      globT, globE, Z, do_xnu, do_forward, do_null2);

  /* Load the slaves.
   * For efficiency reasons, we don't want the master to
   * load HMMs from disk until she absolutely needs them.
   */
  for (nhmm = 0; nhmm < nslaves && nhmm < hmmfp->gsi->recnum; nhmm++) {
    pvm_initsend(PvmDataDefault);
    pvm_pkint(&nhmm, 1, 1);	/* side effect: also tells him what number he is. */
    pvm_send(slave_tid[nhmm], HMMPVM_WORK);
    hmmlist[nhmm] = nhmm;
  }
  SQD_DPRINTF1(("%d slaves are loaded\n", nhmm));

  
  /* Receive/send loop
   */
  for (; nhmm < hmmfp->gsi->recnum; nhmm++)
    {
				/* check slaves before blocking */
      PVMCheckSlaves(slave_tid, nslaves);
				/* receive output */
      SQD_DPRINTF1(("Waiting for a slave to give me output...\n"));
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&slaveidx, 1, 1);     /* # of slave who's sending us stuff */
      pvm_upkstr(slavename);           /* name of HMM that slave did */
      pvm_upkfloat(&sc, 1, 1);         /* score   */
      pvm_upkdouble(&pvalue, 1, 1);	   /* P-value */
      pvm_upkint(&sent_trace, 1, 1);   /* TRUE if trace is coming */
      tr = (sent_trace) ? PVMUnpackTrace() : NULL;
      SQD_DPRINTF1(("Slave %d finished %s for me...\n", slaveidx, slavename));

				/* send new work */
      pvm_initsend(PvmDataDefault);
      pvm_pkint(&nhmm, 1, 1);
      pvm_send(slave_tid[slaveidx], HMMPVM_WORK);
      SQD_DPRINTF1(("Assigned %d -> slave %d\n", nhmm, slaveidx));

				/* process output */
      /* 1b. Store scores/pvalue for each HMM aligned to this sequence, overall
       */
      SQD_DPRINTF1(("%15s : %2d : %f\n", slavename, slaveidx, sc));
      if (sent_trace)
	{ 
				/* now load the HMM, because the hit is significant */
	  HMMFilePositionByIndex(hmmfp, hmmlist[slaveidx]);
	  if (!HMMFileRead(hmmfp, &hmm))
	    { pvm_exit(); Die("Unexpected failure to read HMM file %s", hmmfile); }
	  if (hmm == NULL) 
	    { pvm_exit(); Die("HMM file %s may be corrupt; parse failed", hmmfile); }
	  P7Logoddsify(hmm, TRUE);
	
	  RegisterHit(ghit, sc, pvalue, sc,
		      0., 0.,	                  /* no mother seq */
		      hmm->name, hmm->desc, 
		      0,0,0,                	  /* seq positions  */
		      0,0,0,	                  /* HMM positions  */
		      0, TraceDomainNumber(tr), /* domain info    */
		      NULL);	                  /* alignment info */
	  record_domains(dhit, hmm, dsq, sqinfo->len, sqinfo->name,
			 tr, pvalue, sc, do_null2);
	  FreePlan7(hmm);
	  P7FreeTrace(tr);
	}
      hmmlist[slaveidx] = nhmm;
    }

  /* Collect the output. all n slaves are still working, so wait for them.
   */
  for (slave = 0; slave < nslaves && slave < nhmm; slave++)
    {
				/* don't check slaves (they're exiting normally);
				   window of vulnerability here to slave crashes */
				/* receive output */
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&slaveidx, 1, 1);     /* slave who's sending us stuff */
      pvm_upkstr(slavename);
      pvm_upkfloat(&sc, 1, 1);         /* one score */
      pvm_upkdouble(&pvalue, 1, 1);	   /* P-value */
      pvm_upkint(&sent_trace, 1, 1);   /* TRUE if trace is coming */
      tr = (sent_trace) ? PVMUnpackTrace() : NULL;

				/* process output */
      SQD_DPRINTF1(("%15s : %2d : %f\n", slavename, slaveidx, sc));
      if (sent_trace)
	{ 
	  /* now load the HMM, because the hit is significant */
	  HMMFilePositionByIndex(hmmfp, hmmlist[slaveidx]);
	  if (!HMMFileRead(hmmfp, &hmm))
	    { pvm_exit(); Die("Unexpected failure to read HMM file %s", hmmfile);}
	  if (hmm == NULL) 
	    { pvm_exit(); Die("HMM file %s may be corrupt; parse failed", hmmfile); }
	  P7Logoddsify(hmm, TRUE);
	  
	  RegisterHit(ghit, sc, pvalue, sc,
		      0., 0.,	                  /* no mother seq */
		      hmm->name, hmm->desc, 
		      0,0,0,                	  /* seq positions  */
		      0,0,0,	                  /* HMM positions  */
		      0, TraceDomainNumber(tr), /* domain info    */
		      NULL);	                  /* alignment info */
	  record_domains(dhit, hmm, dsq, sqinfo->len, sqinfo->name, 
			 tr, pvalue, sc, do_null2);
	  FreePlan7(hmm);
	  P7FreeTrace(tr);
	}
				/* send cleanup/shutdown flag */
      pvm_initsend(PvmDataDefault);
      msg = -1;
      pvm_pkint(&msg, 1, 1);
      pvm_send(slave_tid[slaveidx], HMMPVM_WORK);
    }

  /* Cleanup; quit the VM; and return
   */
  free(slave_tid);
  free(hmmlist);
  free(dsq);
  pvm_exit();
  *ret_nhmm = nhmm;
  return;
}



/* Function: init_slaves()
 * Date:     SRE, Thu Aug 13 16:11:09 1998 [St. Louis]
 *
 * Purpose:  Send initialization information to slaves.
 *           Wait until slaves reply "ok" before going ahead.
 *
 * Args:     slave_tid- slave tid array to broadcast to
 *           nslaves  - number of slaves
 *           hmmfile  - name of hmm file for slaves to open
 *           seq      - sequence to search 
 *           L        - length of sequence
 *           globT    - bit threshold
 *           globE    - Evalue threshold
 *           Z        - effective # seqs to use in calc of Eval
 *           do_xnu   - TRUE to use XNU on sequence
 *           do_forward - TRUE to use Forward() to score
 *           do_null2 - TRUE to apply null2 bias filter
 *
 * Returns:  void
 */
static void
init_slaves(int *slave_tid, int nslaves, char *hmmfile, char *seq, int L,
	    float globT, double globE, int Z, 
	    int do_xnu, int do_forward, int do_null2)
{
  int arglen;
  int i;
  int initflag;

  /* Send init info
   */
  SQD_DPRINTF1(("Broadcasting to %d slaves...\n", nslaves));
  pvm_initsend(PvmDataDefault);
  arglen = strlen(hmmfile);
  pvm_pkint(&arglen, 1, 1);
  pvm_pkstr(hmmfile);
  pvm_pkint(&L, 1, 1);
  pvm_pkstr(seq);
  pvm_pkfloat(&globT, 1, 1);
  pvm_pkdouble(&globE, 1, 1);
  pvm_pkint(&Z, 1, 1);
  pvm_pkint(&do_xnu, 1, 1);
  pvm_pkint(&do_forward, 1, 1);
  pvm_pkint(&do_null2, 1, 1);
  pvm_pkint(&Alphabet_type, 1, 1);
  pvm_mcast(slave_tid, nslaves, HMMPVM_INIT);
  SQD_DPRINTF1(("Slaves should be ready...\n"));

  /* Wait for slaves to say ok; else they send error codes
   */
  for (i = 0; i < nslaves; i++)
    {
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&initflag, 1, 1);
      if (initflag != HMMPVM_OK)
	{
	  PVMKillSlaves(slave_tid, nslaves);
	  pvm_exit();
	  switch (initflag) {
	  case HMMPVM_NO_HMMFILE: 
	    Die("One or more PVM slaves couldn't open hmm file. Check installation of %s", hmmfile);
	  case HMMPVM_NO_INDEX:
	    Die("One or more PVM slaves couldn't open GSI index. Check installation of %s", hmmfile);
	  default:
	    Die("Unknown error code. A slave is confused.");
	  }
	}
    }
  return;
}
#endif /*HMMER_PVM*/


#ifdef HMMER_THREADS
/*****************************************************************
 * POSIX threads implementation.
 * 
 * API:
 *      workpool_start()   (makes a workpool_s structure. Starts calculations.)
 *      workpool_stop()    (waits for threads to finish.)
 *      workpool_free()    (destroys the structure)
 *      
 * Threads:
 *      worker_thread()    (the actual parallelized worker thread).
 *****************************************************************/

/* Function: workpool_start()
 * Date:     SRE, Mon Sep 28 11:10:58 1998 [St. Louis]
 *
 * Purpose:  Initialize a workpool_s structure, and return it.
 *
 * Args:     hmmfile    - name of HMM file
 *           hmmfp      - open HMM file, at start
 *           dsq        - ptr to sequence to search
 *           seqname    - ptr to name of dsq
 *           L          - length of dsq
 *           do_forward - TRUE to score using Forward
 *           do_null2   - TRUE to apply null2 ad hoc correction
 *           globT      - per-sequence score threshold
 *           globE      - per-sequence E-value threshold
 *           Z          - effective # of seqs in database             
 *           ghit       - per-seq hit list
 *           dhit       - per-domain hit list             
 *           num_threads- number of worker threads to run.
 *
 * Returns:  ptr to struct workpool_s.
 *           Caller must wait for threads to finish with workpool_stop(),
 *           then free the structure with workpool_free().
 */
static struct workpool_s *
workpool_start(char *hmmfile, HMMFILE *hmmfp, char *dsq, char *seqname, int L,
	       int do_forward, int do_null2, float globT, double globE, int Z, 
	       struct tophit_s *ghit, struct tophit_s *dhit, 
	       int num_threads)
{
  struct workpool_s *wpool;
  pthread_attr_t    attr;
  int i;
  int rtn;

  wpool             = MallocOrDie(sizeof(struct workpool_s));
  wpool->thread     = MallocOrDie(num_threads * sizeof(pthread_t));
  wpool->hmmfile    = hmmfile;
  wpool->dsq        = dsq;
  wpool->L          = L;
  wpool->seqname    = seqname;
  wpool->do_forward = do_forward;
  wpool->do_null2   = do_null2;
  wpool->globT      = globT;
  wpool->globE      = globE;
  wpool->Z          = Z;

  wpool->hmmfp      = hmmfp;
  wpool->nhmm       = 0;
  if ((rtn = pthread_mutex_init(&(wpool->input_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));

  wpool->ghit       = ghit;
  wpool->dhit       = dhit;
  if ((rtn = pthread_mutex_init(&(wpool->output_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));

  wpool->num_threads= num_threads;

  /* Create slave threads
   */
  pthread_attr_init(&attr);
#ifdef HAVE_PTHREAD_SETCONCURRENCY
  pthread_setconcurrency(num_threads+1);
#endif
  for (i = 0; i < num_threads; i++)
    if ((rtn = pthread_create(&(wpool->thread[i]), &attr,
			      worker_thread , (void *) wpool)) != 0)
      Die("Failed to create thread %d; return code %d\n", i, rtn);

  pthread_attr_destroy(&attr);
  return wpool;
}

/* Function: workpool_stop()
 * Date:     SRE, Thu Jul 16 11:20:16 1998 [St. Louis]
 *
 * Purpose:  Waits for threads in a workpool to finish.
 *
 * Args:     wpool -- ptr to the workpool structure
 *
 * Returns:  (void)
 */
static void
workpool_stop(struct workpool_s *wpool)
{
  int i;
				/* wait for threads to stop */
  for (i = 0; i < wpool->num_threads; i++)
    if (pthread_join(wpool->thread[i],NULL) != 0)
      Die("pthread_join failed");
  return;
}

/* Function: workpool_free()
 * Date:     SRE, Thu Jul 16 11:26:27 1998 [St. Louis]
 *
 * Purpose:  Free a workpool_s structure, after the threads
 *           have finished.
 *
 * Args:     wpool -- ptr to the workpool.
 *
 * Returns:  (void)
 */
static void
workpool_free(struct workpool_s *wpool)
{
  free(wpool->thread);
  free(wpool);
  return;
}


/* Function: worker_thread()
 * Date:     SRE, Mon Sep 28 10:48:29 1998 [St. Louis]
 *
 * Purpose:  The procedure executed by the worker threads.
 *
 * Args:     ptr  - (void *) that is recast to a pointer to
 *                  the workpool.
 *
 * Returns:  (void *)
 */
void *
worker_thread(void *ptr)
{
  struct workpool_s *wpool;     /* our working threads structure   */
  struct plan7_s    *hmm;       /* an HMM to search with           */
  struct p7trace_s  *tr;        /* traceback from an alignment     */
  float  sc;			/* score of an alignment           */
  int    rtn;			/* a return code from pthreads lib */
  double pvalue;		/* P-value of score                */

  wpool = (struct workpool_s *) ptr;
  for (;;) {

    /* 1. acquire lock on HMM input, and get
     *    the next HMM to work on.
     */
				/* acquire a lock */
    if ((rtn = pthread_mutex_lock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    wpool->nhmm++;
    
    if (! HMMFileRead(wpool->hmmfp, &hmm)) 
      {	/* we're done. release lock, exit thread */
	if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
	  Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
	pthread_exit(NULL);
      }
    SQD_DPRINTF1(("a thread is working on %s\n", hmm->name));
				/* release the lock */
    if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

    if (hmm == NULL)
      Die("HMM file %s may be corrupt or in incorrect format; parse failed", wpool->hmmfile);
    P7Logoddsify(hmm, !(wpool->do_forward));

    /* 2. We have an HMM in score form.
     *    Score the sequence.
     */
    if (P7ViterbiSize(wpool->L, hmm->M) <= RAMLIMIT)
      sc = P7Viterbi(wpool->dsq, wpool->L, hmm, &tr);
    else
      sc = P7SmallViterbi(wpool->dsq, wpool->L, hmm, &tr);
    
    if (wpool->do_forward) sc  = P7Forward(wpool->dsq, wpool->L, hmm, NULL);
    if (wpool->do_null2)   sc -= TraceScoreCorrection(hmm, tr, wpool->dsq);

    
    /* 3. Save the output in tophits structures, after acquiring a lock
     */
    if ((rtn = pthread_mutex_lock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    SQD_DPRINTF1(("model %s scores %f\n", hmm->name, sc));

    pvalue = PValue(hmm, sc);
    if (sc > wpool->globT && pvalue * (float) wpool->Z < wpool->globE) 
      { 
	RegisterHit(wpool->ghit, sc, pvalue, sc,
		    0., 0.,	                  /* no mother seq */
		    hmm->name, hmm->desc, 
		    0,0,0,                	  /* seq positions  */
		    0,0,0,	                  /* HMM positions  */
		    0, TraceDomainNumber(tr), /* domain info    */
		    NULL);	                  /* alignment info */
	record_domains(wpool->dhit, hmm, wpool->dsq, wpool->L, wpool->seqname, 
		       tr, pvalue, sc, wpool->do_null2);
      }
    if ((rtn = pthread_mutex_unlock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

    P7FreeTrace(tr);
    FreePlan7(hmm);

  } /* end 'infinite' loop over HMMs in this thread */
}

#endif /* HMMER_THREADS */
