/************************************************************
 * @LICENSE@
 ************************************************************/

/* hmmsearch.c
 * SRE, Tue Jan  7 17:19:20 1997 [St. Louis]
 *
 * Search a sequence database with a profile HMM.
 * Conditionally includes PVM parallelization when HMMER_PVM is defined
 *    at compile time; hmmsearch --pvm runs the PVM version.
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

static char banner[] = "hmmsearch - search a sequence database with a profile HMM";

static char usage[]  = "\
Usage: hmmsearch [-options] <hmmfile> <sequence file or database>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
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
   --pvm     : run on a Parallel Virtual Machine (PVM)\n\
   --xnu     : turn ON XNU filtering of target protein sequences\n\
";

static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-A",        TRUE,  sqdARG_INT  },  
  { "-E",        TRUE,  sqdARG_FLOAT},  
  { "-T",        TRUE,  sqdARG_FLOAT},  
  { "-Z",        TRUE,  sqdARG_INT  },
  { "--cpu",     FALSE, sqdARG_INT  },
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
  struct plan7_s *hmm;		/* HMM to search with              */
  int    do_xnu;		/* TRUE to apply XNU filter  */
  int    do_forward;		/* TRUE to score using Forward     */
  int    do_null2;		/* TRUE to apply null2 ad hoc correction */
  float  globT;                 /* global score threshold          */
  double globE;			/* global E-value threshold        */
  int    Z;                     /* effective # of seqs in database */
  
  /* Shared (mutex-protected) input resources:
   */
  SQFILE *sqfp;                 /* ptr to open sequence file      */
  int nseq;			/* number of seqs searched so far */
  pthread_mutex_t input_lock;   /* mutex for locking input        */

  /* Shared (mutex-protected) output resources:
   */
  struct tophit_s *ghit;        /* per-sequence top hits */
  struct tophit_s *dhit;        /* per-domain top hits */
  struct histogram_s *hist;     /* histogram of scores */
  pthread_mutex_t output_lock;  /* mutex for locking output */

  /* Thread pool information
   */
  pthread_t *thread;            /* our pool of threads */
  int        num_threads;       /* number of threads   */
};
static struct workpool_s *workpool_start(struct plan7_s *hmm, SQFILE *sqfp, 
					 int do_xnu, int do_forward, int do_null2, float globT, double globE, int Z, 
					 struct tophit_s *ghit, struct tophit_s *dhit, 
					 struct histogram_s *hist, int num_threads);
static void  workpool_stop(struct workpool_s *wpool);
static void  workpool_free(struct workpool_s *wpool);
static void *worker_thread(void *ptr);
#endif /* HMMER_THREADS */

static void record_domains(struct tophit_s *h, 
			   struct plan7_s *hmm, char *dsq, 
			   char *sqname, char *sqdesc, int L,
			   struct p7trace_s *tr, double whole_pval, float whole_sc,
			   int do_null2);
static void main_loop_serial(struct plan7_s *hmm, SQFILE *sqfp, 
			     float globT, double globE, int Z, int do_forward,
			     int do_null2, int do_xnu, int num_threads,
			     struct histogram_s *histogram, struct tophit_s *ghit, 
			     struct tophit_s *dhit, int *ret_nseq);
#ifdef HMMER_PVM
static void main_loop_pvm(struct plan7_s *hmm, SQFILE *sqfp, 
			  float globT, double globE, int Z, int do_forward,
			  int do_null2, int do_xnu, struct histogram_s *histogram, 
			  struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nseq);
#endif


int
main(int argc, char **argv) 
{
  char    *hmmfile;	        /* file to read HMM(s) from                */
  HMMFILE *hmmfp;               /* opened hmmfile for reading              */
  char    *seqfile;             /* file to read target sequence(s) from    */ 
  SQFILE   *sqfp;               /* opened seqfile for reading              */
  int       format;	        /* format of seqfile                       */
  int       i; 
  struct plan7_s  *hmm;         /* HMM to search with                      */ 
  struct histogram_s *histogram;/* histogram of all scores                 */
  struct fancyali_s *ali;       /* displayed alignment info                */ 
  struct tophit_s   *ghit;      /* list of top hits for whole sequences    */
  struct tophit_s   *dhit;	/* list of top hits for domains            */

  float     sc;	        	/* score of an HMM search                  */
  double  pvalue;		/* pvalue of an HMM score                  */
  double  evalue;		/* evalue of an HMM score                  */
  double  motherp;		/* pvalue of a whole seq HMM score         */
  float   mothersc;		/* score of a whole seq parent of domain   */
  int     sqfrom, sqto;		/* coordinates in sequence                 */
  int     hmmfrom, hmmto;	/* coordinate in HMM                       */
  char   *name, *desc;          /* hit sequence name and description       */
  int     sqlen;		/* length of seq that was hit              */
  int     nseq;			/* number of sequences searched            */
  int     Z;			/* # of seqs for purposes of E-val calc    */
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
  int   do_null2;		/* TRUE to adjust scores with null model #2 */
  int   do_forward;		/* TRUE to use Forward() not Viterbi()      */
  int   do_xnu;			/* TRUE to filter sequences thru XNU        */
  int   do_pvm;			/* TRUE to run on Parallel Virtual Machine  */
  int   num_threads;		/* number of worker threads                 */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  do_forward  = FALSE;
  do_null2    = TRUE;
  do_xnu      = FALSE;
  do_pvm      = FALSE;  
  Z           = 0;

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
    if      (strcmp(optname, "-A") == 0)        Alimit     = atoi(optarg);  
    else if (strcmp(optname, "-E") == 0)        globE      = atof(optarg);
    else if (strcmp(optname, "-T") == 0)        globT      = atof(optarg);
    else if (strcmp(optname, "-Z") == 0)        Z          = atoi(optarg);
    else if (strcmp(optname, "--cpu")     == 0) num_threads= atoi(optarg);
    else if (strcmp(optname, "--domE")    == 0) domE       = atof(optarg);
    else if (strcmp(optname, "--domT")    == 0) domT       = atof(optarg);
    else if (strcmp(optname, "--forward") == 0) do_forward = TRUE;
    else if (strcmp(optname, "--null2")   == 0) do_null2   = FALSE;
    else if (strcmp(optname, "--pvm")     == 0) do_pvm     = TRUE;
    else if (strcmp(optname, "--xnu")     == 0) do_xnu     = TRUE;
    else if (strcmp(optname, "-h") == 0) {
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
  if (do_pvm) Die("PVM support is not compiled into your HMMER software; --pvm doesn't work.");
#endif
#ifndef HMMER_THREADS
  if (num_threads) Die("Posix threads support is not compiled into HMMER; --cpu doesn't have any effect");
#endif

  /*********************************************** 
   * Open sequence database (might be in BLASTDB or current directory)
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, "BLASTDB"))
    switch (squid_errno) {
    case SQERR_NOFILE: 
      Die("Sequence file %s could not be opened for reading", seqfile); break;
    case SQERR_FORMAT: 
    default:           
      Die("Failed to determine format of sequence file %s", seqfile);
    }
  if ((sqfp = SeqfileOpen(seqfile, format, "BLASTDB")) == NULL)
    Die("Failed to open sequence database file %s\n%s\n", seqfile, usage);

  /*********************************************** 
   * Open HMM file (might be in HMMERDB or current directory).
   * Read a single HMM from it. (Config HMM, if necessary).
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  if (hmm == NULL) 
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  P7Logoddsify(hmm, !do_forward);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  Banner(stdout, banner);
  printf(   "HMM file:                 %s [%s]\n", hmmfile, hmm->name);
  printf(   "Sequence database:        %s\n", seqfile); 
  if (do_pvm)
    printf( "PVM:                      ACTIVE\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  /*********************************************** 
   * Search HMM against each sequence
   ***********************************************/

				/* set up structures for storing output  */
  histogram = AllocHistogram(-200, 200, 100);  /* keeps full histogram */
  ghit      = AllocTophits(200);         /* per-seq hits: 200=lumpsize */
  dhit      = AllocTophits(200);         /* domain hits:  200=lumpsize */

  if (! do_pvm)
    main_loop_serial(hmm, sqfp, 
		     globT, globE, Z, do_forward, do_null2, do_xnu, num_threads,
		     histogram, ghit, dhit, &nseq);
#ifdef HMMER_PVM
  else
    main_loop_pvm(hmm, sqfp, 		     
		  globT, globE, Z, do_forward, do_null2, do_xnu, 
		  histogram, ghit, dhit, &nseq);
#endif

  /*********************************************** 
   * Process hit lists, produce text output
   ***********************************************/

  /* Set the theoretical EVD curve in our histogram using 
   * calibration in the HMM, if available. 
   */
  if (hmm->flags & PLAN7_STATS)
    ExtremeValueSetHistogram(histogram, hmm->mu, hmm->lambda, 
			     histogram->lowscore, histogram->highscore, 0);
  if (!Z) Z = nseq;		/* set Z for good now that we're done. */

  /* Format and report our output 
   */
  /* 1. Report overall sequence hits (sorted on E-value) */
  printf("\nQuery HMM:  %s  %s\n", 
	 hmm->name, hmm->desc != NULL ? hmm->desc : "");
  if (hmm->flags & PLAN7_STATS)
    printf("  [HMM has been calibrated; E-values are empirical estimates]\n");
  else
    printf("  [No calibration for HMM; E-values are upper bounds]\n");

  FullSortTophits(ghit);
  namewidth = MAX(8, TophitsMaxName(ghit));

  printf("\nScores for complete sequences (score includes all domains):\n");
  printf("%-*s %-*s %7s %10s %3s\n", namewidth, "Sequence", 52-namewidth, "Description", "Score", "E-value", " N ");
  printf("%-*s %-*s %7s %10s %3s\n", namewidth, "--------", 52-namewidth, "-----------", "-----", "-------", "---");
  for (i = 0; i < ghit->num; i++)
    {
      char *safedesc;
      GetRankedHit(ghit, i, 
		   &pvalue, &sc, NULL, NULL,
		   &name, &desc,
		   NULL, NULL, NULL,               /* sequence positions */
		   NULL, NULL, NULL,               /* HMM positions      */
		   NULL, &ndom,	                   /* domain info        */
		   NULL);	                   /* alignment info     */
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




  /* 2. Report domain hits (also sorted on E-value) */
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
		   &sqfrom, &sqto, &sqlen,            /* seq position info  */
		   &hmmfrom, &hmmto, NULL,            /* HMM position info  */
		   &domidx, &ndom,                    /* domain info        */
		   NULL);	                      /* alignment info     */
      evalue = pvalue * (double) Z;

      if (motherp * (double) Z >= globE || mothersc <= globT) 
	continue;
      else if (evalue < domE && sc > domT)
	printf("%-*s %3d/%-3d %5d %5d %c%c %5d %5d %c%c %7.1f %8.2g\n",
	       namewidth, name, 
	       domidx, ndom,
	       sqfrom, sqto, 
	       sqfrom == 1 ? '[' : '.', sqto == sqlen ? ']' : '.',
	       hmmfrom, hmmto,
	       hmmfrom == 1 ? '[':'.', hmmto == hmm->M ? ']' : '.',
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
  if (i == 0) printf("\t[no hits above thresholds\n");



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
		       &sqfrom, &sqto, &sqlen,            /* seq position info  */
		       &hmmfrom, &hmmto, NULL,            /* HMM position info  */
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
	    if (i > 0) printf("\t[no more alignments below domE threshold\n");
	    break;
	  }
	  else if (sc <= domT) {
	    if (i > 0) printf("\t[no more alignments above domT threshold\n");
	    break;
	  }
	}
      if (i == 0)      printf("\t[no hits above thresholds\n");
      if (i == Alimit) printf("\t[output cut off at A = %d top alignments]\n", Alimit);
    }

  /* 4. Histogram output */
  printf("\nHistogram of all scores:\n");
  PrintASCIIHistogram(stdout, histogram);

  /* 5. Tophits summaries, while developing...
   */
  printf("\nTotal sequences searched: %d\n", nseq);
  printf("\nWhole sequence top hits:\n");
  TophitsReport(ghit, globE, nseq);
  printf("\nDomain top hits:\n");
  TophitsReport(dhit, domE, nseq);

  /*********************************************** 
   * Clean-up and exit.
   ***********************************************/

  FreeHistogram(histogram);
  HMMFileClose(hmmfp);
  SeqfileClose(sqfp);
  FreeTophits(ghit);
  FreeTophits(dhit);
  FreePlan7(hmm);
  SqdClean();

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) malloc_list(2, histid1, histid2);
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  return 0;
}


/* Function: main_loop_serial()
 * Date:     SRE, Wed Sep 23 10:20:49 1998 [St. Louis]
 *
 * Purpose:  Search an HMM against a sequence database.
 *           main loop for the serial (non-PVM, non-threads)
 *           version.
 *           
 *           In:   HMM and open sqfile, plus options
 *           Out:  histogram, global hits list, domain hits list, nseq.
 *
 * Args:     hmm        - the HMM to search with. 
 *           sqfp       - open SQFILE for sequence database
 *           globT      - bit score significance threshold 
 *           globE      - E value significance threshold            
 *           Z          - 0, or forced number of seqs for E-value calc's
 *           do_forward - TRUE to score using Forward()        
 *           do_null2   - TRUE to use ad hoc null2 score correction
 *           do_xnu     - TRUE to apply XNU mask
 *           num_threads- number of worker threads to start, or 0
 *           histogram  - RETURN: score histogram
 *           ghit       - RETURN: ranked global scores
 *           dhit       - RETURN: ranked domain scores
 *           ret_nseq   - RETURN: actual number of seqs searched
 *           
 * Returns:  (void)
 */
static void
main_loop_serial(struct plan7_s *hmm, SQFILE *sqfp, 
		 float globT, double globE, int Z, int do_forward,
		 int do_null2, int do_xnu, int num_threads,
		 struct histogram_s *histogram, 
		 struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nseq)
{
#ifdef HMMER_THREADS
  struct workpool_s *wpool;	/* pool of worker threads                  */
#else
  struct p7trace_s *tr;         /* traceback                               */
  char   *seq;                  /* target sequence                         */
  char   *dsq;		        /* digitized target sequence               */
  SQINFO sqinfo;		/* optional info for seq                   */
  float  sc;	        	/* score of an HMM search                  */
  double pvalue;		/* pvalue of an HMM score                  */
  double evalue;		/* evalue of an HMM score                  */
#endif
  int    nseq;			/* number of sequences searched            */
 
#ifdef HMMER_THREADS
  wpool = workpool_start(hmm, sqfp, do_xnu, do_forward, do_null2, globT, globE, Z, 
			 ghit, dhit, histogram, num_threads);
  workpool_stop(wpool);
  nseq = wpool->nseq;
  workpool_free(wpool);

#else /* unthreaded code: */
  nseq = 0;
  while (ReadSeq(sqfp, sqfp->format, &seq, &sqinfo))
    {
				/* silently skip len 0 seqs */
      if (sqinfo.len == 0) continue;

      nseq++;
      dsq = DigitizeSequence(seq, sqinfo.len);
      
      if (do_xnu) XNU(dsq, sqinfo.len);
      
      /* 1. Recover a trace by Viterbi.
       */
      if (P7ViterbiSize(sqinfo.len, hmm->M) <= RAMLIMIT)
	sc = P7Viterbi(dsq, sqinfo.len, hmm, &tr);
      else
	sc = P7SmallViterbi(dsq, sqinfo.len, hmm, &tr);

      /* 2. If we're using Forward scores, do another DP
       *    to get it; else, we already have a Viterbi score
       *    in sc.
       */
      if (do_forward) sc  = P7Forward(dsq, sqinfo.len, hmm, NULL);
      if (do_null2)   sc -= TraceScoreCorrection(hmm, tr, dsq);

#if DEBUGLEVEL >= 2
      P7PrintTrace(stdout, tr, hmm, dsq); 
#endif

      /* 2. Store score/pvalue for global alignment; will sort on score. 
       *    Keep all domains in a significant sequence hit.
       *    We can only make a lower bound estimate of E-value since
       *    we don't know the final value of nseq yet. 
       */
      pvalue = PValue(hmm, sc);
      evalue = Z ? (double) Z * pvalue : (double) nseq * pvalue;
      if (sc > globT && evalue < globE) 
	{
	  RegisterHit(ghit, sc, pvalue, sc, 
		      0., 0.,	                /* no mother seq */
		      sqinfo.name, 
		      sqinfo.flags & SQINFO_DESC ? sqinfo.desc : NULL, 
		      0,0,0,                	/* seq positions  */
		      0,0,0,	                /* HMM positions  */
		      0, TraceDomainNumber(tr), /* domain info    */
		      NULL);	                /* alignment info */

	  record_domains(dhit, hmm, dsq, sqinfo.name, sqinfo.desc, sqinfo.len, 
			 tr, pvalue, sc, do_null2); 
	}
      AddToHistogram(histogram, sc);

      FreeSequence(seq, &sqinfo); 
      P7FreeTrace(tr);
      free(dsq);
    }
#endif

  *ret_nseq = nseq;
  return;
}


/* Function: record_domains()
 * Date:     SRE, Tue Nov  4 11:25:14 1997 [St. Louis]
 * 
 * Purpose:  Decompose a trace, and register scores, P-values, alignments,
 *           etc. for individual domains in a hitlist. 
 *           
 * Args:     hmm    - the HMM structure
 *           dsq    - digitized sequence 1..L
 *           sqname - name of sequence
 *           sqdesc - description of sequence
 *           L      - length of sequence         
 *           tr     - traceback of the whole sequence aligned to HMM
 *           whole_pval - P-value of complete alignment
 *           whole_sc   - score of complete alignment (bits)
 *           do_null2   - TRUE to use post hoc null model correction 
 *           
 * Return:   (void)          
 */
static void
record_domains(struct tophit_s *h, 
	       struct plan7_s *hmm, char *dsq, 
	       char *sqname, char *sqdesc, int L,
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

  SQD_DPRINTF1(("dsq is %d%d%d%d\n", dsq[0], dsq[1], dsq[2], dsq[3]));
  for (idx = 0; idx < ntr; idx++)
    {
      /* Get the score and bounds of the match.
       */
      score  = P7TraceScore(hmm, dsq, tarr[idx]);
      if (do_null2) 
	score -= TraceScoreCorrection(hmm, tarr[idx], dsq);
      pvalue = PValue(hmm, score); 
      TraceSimpleBounds(tarr[idx], &i1, &i2, &k1, &k2);
      
#if DEBUGLEVEL >= 2
      P7PrintTrace(stdout, tarr[idx], hmm, dsq); 
#endif

      /* Record the match. Use score as the sort key.
       */
      ali = CreateFancyAli(tarr[idx], hmm, dsq, sqname);
      RegisterHit(h, score, pvalue, score, whole_pval, whole_sc,
		  sqname, 
		  sqdesc,
		  i1,i2, L, 
		  k1,k2, hmm->M, 
		  idx+1, ntr,
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
 * Date:     SRE, Wed Sep 23 10:36:44 1998 [St. Louis]
 *
 * Purpose:  Search an HMM against a sequence database.
 *           main loop for the PVM version.
 *           
 *           In:   HMM and open sqfile, plus options
 *           Out:  histogram, global hits list, domain hits list, nseq.
 *
 * Args:     hmm        - the HMM to search with. scoring form.
 *           sqfp       - open SQFILE for sequence database
 *           globT      - bit score significance threshold 
 *           globE      - E value significance threshold            
 *           Z          - 0, or forced number of seqs for E-value calc's
 *           do_forward - TRUE to score using Forward()        
 *           do_null2   - TRUE to use ad hoc null2 score correction
 *           do_xnu     - TRUE to apply XNU mask
 *           histogram  - RETURN: score histogram
 *           ghit       - RETURN: ranked global scores
 *           dhit       - RETURN: ranked domain scores
 *           ret_nseq   - RETURN: actual number of seqs searched
 *           
 * Returns:  (void)
 */
static void
main_loop_pvm(struct plan7_s *hmm, SQFILE *sqfp, 
	      float globT, double globE, int Z, int do_forward,
	      int do_null2, int do_xnu, struct histogram_s *histogram, 
	      struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nseq)
{
  char *seq;                    /* target sequence */
  char *dsq;                    /* digitized target seq */
  SQINFO sqinfo;                /* optional info about target seq */
  int   master_tid;		/* master's (my) PVM TID */
  int  *slave_tid;              /* array of slave TID's  */
  int   nslaves;		/* number of slaves      */
  int   code;			/* status code rec'd from a slave */
  int   nseq;			/* number of sequences searched */
  int   sent_trace;		/* TRUE if slave gave us a trace */
  char **dsqlist;               /* remember what seqs slaves are doing */
  char **namelist;              /* remember what seq names slaves are doing */
  char **desclist;              /* remember what seq desc's slaves are doing */
  int   *lenlist;               /* remember lengths of seqs slaves are doing */
  int    slaveidx;		/* counter for slaves */
  float  sc;			/* score of an alignment */
  double pvalue;		/* P-value of a score of an alignment */
  struct p7trace_s *tr;         /* Viterbi traceback of an alignment */
  int    i;			/* generic counter */

  /* Initialize PVM.
   */
  SQD_DPRINTF1(("Requesting master TID...\n"));
  master_tid = pvm_mytid();
  SQD_DPRINTF1(("Spawning slaves...\n"));
  PVMSpawnSlaves("hmmsearch-pvm", &slave_tid, &nslaves);
  SQD_DPRINTF1(("Spawned a total of %d slaves...\n", nslaves));
 
  /* Initialize the slaves by broadcast.
   */
  SQD_DPRINTF1(("Broadcasting to %d slaves...\n", nslaves));
  pvm_initsend(PvmDataDefault);
  pvm_pkfloat(&globT, 1, 1);
  pvm_pkdouble(&globE, 1, 1);
  pvm_pkint(&Z,          1, 1);
  pvm_pkint(&do_forward, 1, 1);
  pvm_pkint(&do_null2,   1, 1);
  pvm_pkint(&Alphabet_type, 1, 1);
  PVMPackHMM(hmm);
  pvm_mcast(slave_tid, nslaves, HMMPVM_INIT);
  SQD_DPRINTF1(("Slaves should be ready...\n"));

  /* Confirm slaves' OK status.
   */
  for (slaveidx = 0; slaveidx < nslaves; slaveidx++)
    {
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&code, 1, 1);
      if (code != HMMPVM_OK)
	{
	  PVMKillSlaves(slave_tid, nslaves);
	  pvm_exit();
	  switch (code) {
	  case HMMPVM_BAD_INIT: 
	    Die("One or more PVM slaves didn't initialize properly.");
	  default:
	    Die("Unknown init error code. One or more slaves is confused.");
	  }
	}
    }
  SQD_DPRINTF1(("Slaves confirm that they're ok...\n"));
  
  /* Alloc arrays for remembering what seq each
   * slave was working on.
   */
  namelist = MallocOrDie(sizeof(char *) * nslaves);
  desclist = MallocOrDie(sizeof(char *) * nslaves);
  dsqlist  = MallocOrDie(sizeof(char *) * nslaves);
  lenlist  = MallocOrDie(sizeof(int) * nslaves);

  /* Load the slaves.
   * Give them all a sequence number and a digitized sequence
   * to work on.
   * A side effect of the seq number is that we assign each slave
   * a number from 0..nslaves-1.
   */
  for (nseq = 0; nseq < nslaves; nseq++)
    {
      if (! ReadSeq(sqfp, sqfp->format, &seq, &sqinfo)) break;
      if (sqinfo.len == 0) { nseq--; continue; }

      dsq = DigitizeSequence(seq, sqinfo.len);
      if (do_xnu) XNU(dsq, sqinfo.len);

      pvm_initsend(PvmDataDefault);
      pvm_pkint(&nseq, 1, 1);
      pvm_pkint(&(sqinfo.len), 1, 1);
      pvm_pkbyte(dsq, sqinfo.len+2, 1);
      pvm_send(slave_tid[nseq], HMMPVM_WORK);
      SQD_DPRINTF1(("sent a dsq : %d bytes\n", sqinfo.len+2));

      namelist[nseq] = Strdup(sqinfo.name);
      desclist[nseq] = (sqinfo.flags & SQINFO_DESC) ? Strdup(sqinfo.desc) : NULL;
      lenlist[nseq]  = sqinfo.len;
      dsqlist[nseq]  = dsq;

      FreeSequence(seq, &sqinfo);
    }
  SQD_DPRINTF1(("%d slaves are loaded\n", nseq));
  
  /* main receive/send loop
   */
  while (ReadSeq(sqfp, sqfp->format, &seq, &sqinfo)) 
    {
      if (sqinfo.len == 0) { continue; }
      nseq++;
				/* check slaves before blocking */
      PVMCheckSlaves(slave_tid, nslaves);

				/* receive output */
      SQD_DPRINTF1(("Waiting for a slave to give me output...\n"));
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&slaveidx, 1, 1);     /* # of slave who's sending us stuff */
      pvm_upkfloat(&sc, 1, 1);         /* score   */
      pvm_upkdouble(&pvalue, 1, 1);    /* P-value */
      pvm_upkint(&sent_trace, 1, 1);   /* TRUE if trace is coming */
      tr = (sent_trace) ? PVMUnpackTrace() : NULL;
      SQD_DPRINTF1(("Slave %d finished %s for me...\n", slaveidx, namelist[slaveidx]));

				/* send new work */
      dsq = DigitizeSequence(seq, sqinfo.len);
      if (do_xnu) XNU(dsq, sqinfo.len);

      pvm_initsend(PvmDataDefault);
      pvm_pkint(&nseq, 1, 1);
      pvm_pkint(&(sqinfo.len), 1, 1);
      pvm_pkbyte(dsq, sqinfo.len+2, 1);
      pvm_send(slave_tid[slaveidx], HMMPVM_WORK);
      
				/* process output */
      if (sent_trace)
	{
	  RegisterHit(ghit, sc, pvalue, sc, 
		      0., 0.,	                /* no mother seq */
		      namelist[slaveidx],
		      desclist[slaveidx],
		      0,0,0,                	/* seq positions  */
		      0,0,0,	                /* HMM positions  */
		      0, TraceDomainNumber(tr), /* domain info    */
		      NULL);	                /* alignment info */
	  record_domains(dhit, hmm, dsqlist[slaveidx], 
			 namelist[slaveidx], desclist[slaveidx], lenlist[slaveidx],
			 tr, pvalue, sc, do_null2); 
	  P7FreeTrace(tr);
	}
      AddToHistogram(histogram, sc);

				/* record seq info for seq we just sent */
      free(namelist[slaveidx]);
      if (desclist[slaveidx] != NULL) free(desclist[slaveidx]);
      free(dsqlist[slaveidx]);

      dsqlist[slaveidx]  = dsq;
      namelist[slaveidx] = Strdup(sqinfo.name);
      desclist[slaveidx] = (sqinfo.flags & SQINFO_DESC) ? Strdup(sqinfo.desc) : NULL;
      lenlist[slaveidx]  = sqinfo.len;

      FreeSequence(seq, &sqinfo); 
    }
  SQD_DPRINTF1(("End of receive/send loop\n"));

  /* Collect the output. All n slaves are still working.
   */
  for (i = 0; i < nslaves && i < nseq; i++)
    {
				/* don't check slaves (they're exiting normally);
				   window of vulnerability here to slave crashes */
				/* receive output */
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&slaveidx, 1, 1);     /* # of slave who's sending us stuff */
      pvm_upkfloat(&sc, 1, 1);         /* score   */
      pvm_upkdouble(&pvalue, 1, 1);    /* P-value */
      pvm_upkint(&sent_trace, 1, 1);   /* TRUE if trace is coming */
      tr = (sent_trace) ? PVMUnpackTrace() : NULL;
      SQD_DPRINTF1(("Slave %d finished %s for me...\n", slaveidx, namelist[slaveidx]));

      			/* process output */
      if (sent_trace)
	{
	  RegisterHit(ghit, sc, pvalue, sc, 
		      0., 0.,	                /* no mother seq */
		      namelist[slaveidx],
		      desclist[slaveidx],
		      0,0,0,                	/* seq positions  */
		      0,0,0,	                /* HMM positions  */
		      0, TraceDomainNumber(tr), /* domain info    */
		      NULL);	                /* alignment info */
	  record_domains(dhit, hmm, dsqlist[slaveidx], 
			 namelist[slaveidx], desclist[slaveidx], lenlist[slaveidx],
			 tr, pvalue, sc, do_null2); 
	  P7FreeTrace(tr);
	}
      AddToHistogram(histogram, sc);

				/* free seq info */
      free(namelist[slaveidx]);
      if (desclist[slaveidx] != NULL) free(desclist[slaveidx]);
      free(dsqlist[slaveidx]);

      				/* send cleanup/shutdown flag to slave */
      pvm_initsend(PvmDataDefault);
      code = -1;
      pvm_pkint(&code, 1, 1);
      pvm_send(slave_tid[slaveidx], HMMPVM_WORK);
    }

  
  /* Cleanup; quit the VM; and return
   */
  free(slave_tid);
  free(dsqlist);
  free(namelist);
  free(desclist);
  free(lenlist);
  pvm_exit();
  *ret_nseq = nseq;
  return;
}
#endif /* HMMER_PVM */

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
 * Date:     SRE, Mon Oct  5 16:44:53 1998
 *
 * Purpose:  Initialize a workpool_s structure, and return it.
 *
 * Args:     sqfp       - open sequence file, at start
 *           do_xnu     - TRUE to apply XNU filter
 *           do_forward - TRUE to score using Forward
 *           do_null2   - TRUE to apply null2 ad hoc correction
 *           globT      - per-sequence score threshold
 *           globE      - per-sequence E-value threshold
 *           Z          - effective # of seqs in database             
 *           ghit       - per-seq hit list
 *           dhit       - per-domain hit list             
 *           hist       - histogram (alloced but empty)
 *           num_threads- number of worker threads to run.
 *
 * Returns:  ptr to struct workpool_s.
 *           Caller must wait for threads to finish with workpool_stop(),
 *           then free the structure with workpool_free().
 */
static struct workpool_s *
workpool_start(struct plan7_s *hmm, SQFILE *sqfp, int do_xnu,
	       int do_forward, int do_null2, float globT, double globE, int Z, 
	       struct tophit_s *ghit, struct tophit_s *dhit, 
	       struct histogram_s *hist, int num_threads)
{
  struct workpool_s *wpool;
  pthread_attr_t    attr;
  int i;
  int rtn;

  wpool             = MallocOrDie(sizeof(struct workpool_s));
  wpool->thread     = MallocOrDie(num_threads * sizeof(pthread_t));
  wpool->hmm        = hmm;

  wpool->do_xnu     = do_xnu;
  wpool->do_forward = do_forward;
  wpool->do_null2   = do_null2;
  wpool->globT      = globT;
  wpool->globE      = globE;
  wpool->Z          = Z;

  wpool->sqfp       = sqfp;
  wpool->nseq       = 0;
  if ((rtn = pthread_mutex_init(&(wpool->input_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));

  wpool->ghit       = ghit;
  wpool->dhit       = dhit;
  wpool->hist       = hist;
  if ((rtn = pthread_mutex_init(&(wpool->output_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));

  wpool->num_threads= num_threads;

  /* Create slave threads
   * IRIX note: IRIX pthreads isn't so great at determining the
   * right number of "execution vehicles". Recommended fix
   * is to provide a hint using pthread_setconcurrency(), but
   * this call is not portable, hence the ifdef.
   */
  pthread_attr_init(&attr);
#ifdef HAVE_PTHREAD_SETCONCURRENCY
  pthread_setconcurrency(num_threads+1);
#endif
  /* pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); */
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
  char  *seq;                   /* target sequence                 */
  SQINFO sqinfo;		/* information assoc w/ seq        */
  char  *dsq;                   /* digitized sequence              */
  struct p7trace_s  *tr;        /* traceback from an alignment     */
  float  sc;			/* score of an alignment           */
  int    rtn;			/* a return code from pthreads lib */
  double pvalue;		/* P-value of score                */
  double evalue;		/* E-value of score                */

  wpool = (struct workpool_s *) ptr;
  for (;;) {

    /* 1. acquire lock on sequence input, and get
     *    the next seq to work on.
     */
				/* acquire a lock */
    if ((rtn = pthread_mutex_lock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    if (! ReadSeq(wpool->sqfp, wpool->sqfp->format, &seq, &sqinfo))
      {	/* we're done. release lock, exit thread */
	if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
	  Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
	pthread_exit(NULL);
      }
    SQD_DPRINTF1(("a thread is working on %s\n", sqinfo->name));
				/* release the lock */
    if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

    if (sqinfo.len == 0) continue; /* silent skip of len=0 seqs (wormpep!?!) */

    wpool->nseq++;
    dsq = DigitizeSequence(seq, sqinfo.len);
    if (wpool->do_xnu) XNU(dsq, sqinfo.len);
      
    /* 1. Recover a trace by Viterbi.
     */
    if (P7ViterbiSize(sqinfo.len, wpool->hmm->M) <= RAMLIMIT)
      sc = P7Viterbi(dsq, sqinfo.len, wpool->hmm, &tr);
    else
      sc = P7SmallViterbi(dsq, sqinfo.len, wpool->hmm, &tr);

    /* 2. If we're using Forward scores, do another DP
     *    to get it; else, we already have a Viterbi score
     *    in sc.
     */
    if (wpool->do_forward) sc  = P7Forward(dsq, sqinfo.len, wpool->hmm, NULL);
    if (wpool->do_null2)   sc -= TraceScoreCorrection(wpool->hmm, tr, dsq);

    /* 3. Save the output in tophits and histogram structures, after acquiring a lock
     */
    if ((rtn = pthread_mutex_lock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    SQD_DPRINTF1(("seq %s scores %f\n", sqinfo.name, sc));

    pvalue = PValue(wpool->hmm, sc);
    evalue = wpool->Z ? (double) wpool->Z * pvalue : (double) wpool->nseq * pvalue;
 
    if (sc > wpool->globT && evalue < wpool->globE) 
      { 
	RegisterHit(wpool->ghit, sc, pvalue, sc,
		    0., 0.,	                  /* no mother seq */
		    sqinfo.name, 
		    sqinfo.flags & SQINFO_DESC ? sqinfo.desc : NULL, 
		    0,0,0,                	  /* seq positions  */
		    0,0,0,	                  /* HMM positions  */
		    0, TraceDomainNumber(tr), /* domain info    */
		    NULL);	                  /* alignment info */
	record_domains(wpool->dhit, wpool->hmm, dsq, 
		       sqinfo.name, sqinfo.desc, sqinfo.len,
		       tr, pvalue, sc, wpool->do_null2);
      }
    AddToHistogram(wpool->hist, sc);
    if ((rtn = pthread_mutex_unlock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

    P7FreeTrace(tr);
    FreeSequence(seq, &sqinfo);
    free(dsq);
  } /* end 'infinite' loop over seqs in this thread */
}

#endif /* HMMER_THREADS */

