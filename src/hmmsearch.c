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
 * CVS $Id$
 */

#include "config.h"		/* compile-time configuration constants */
#include "squidconf.h"

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
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */

static char banner[] = "hmmsearch - search a sequence database with a profile HMM";

static char usage[]  = "\
Usage: hmmsearch [-options] <hmmfile> <sequence file or database>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
   -A <n>    : sets alignment output limit to <n> best domain alignments\n\
   -E <x>    : sets E value cutoff (globE) to <= x\n\
   -T <x>    : sets T bit threshold (globT) to >= x\n\
   -Z <n>    : sets Z (# seqs) for E-value calculation\n\
";

static char experts[] = "\
   --compat       : make best effort to use last version's output style\n\
   --cpu <n>      : run <n> threads in parallel (if threaded)\n\
   --cut_ga       : use Pfam GA gathering threshold cutoffs\n\
   --cut_nc       : use Pfam NC noise threshold cutoffs\n\
   --cut_tc       : use Pfam TC trusted threshold cutoffs\n\
   --domE <x>     : sets domain Eval cutoff (2nd threshold) to <= x\n\
   --domT <x>     : sets domain T bit thresh (2nd threshold) to >= x\n\
   --forward      : use the full Forward() algorithm instead of Viterbi\n\
   --informat <s> : sequence file is in format <s>\n\
   --null2        : turn OFF the post hoc second null model\n\
   --pvm          : run on a Parallel Virtual Machine (PVM)\n\
   --xnu          : turn ON XNU filtering of target protein sequences\n\
";

static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-A",        TRUE,  sqdARG_INT  },  
  { "-E",        TRUE,  sqdARG_FLOAT},  
  { "-T",        TRUE,  sqdARG_FLOAT},  
  { "-Z",        TRUE,  sqdARG_INT  },
  { "--compat",  FALSE, sqdARG_NONE },
  { "--cpu",     FALSE, sqdARG_INT  },
  { "--cut_ga",  FALSE, sqdARG_NONE },
  { "--cut_nc",  FALSE, sqdARG_NONE },
  { "--cut_tc",  FALSE, sqdARG_NONE },
  { "--domE",    FALSE, sqdARG_FLOAT},
  { "--domT",    FALSE, sqdARG_FLOAT},
  { "--forward", FALSE, sqdARG_NONE },  
  { "--informat",FALSE, sqdARG_STRING},
  { "--null2",   FALSE, sqdARG_NONE },
  { "--pvm",     FALSE, sqdARG_NONE },
  { "--xnu",     FALSE, sqdARG_NONE },

};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


static void main_loop_serial(struct plan7_s *hmm, SQFILE *sqfp, struct threshold_s *thresh,
			     int do_forward, int do_null2, int do_xnu, 
			     struct histogram_s *histogram, struct tophit_s *ghit, 
			     struct tophit_s *dhit, int *ret_nseq);
static void main_loop_pvm(struct plan7_s *hmm, SQFILE *sqfp, struct threshold_s *thresh, 
			  int do_forward, int do_null2, int do_xnu, 
			  struct histogram_s *histogram, struct tophit_s *ghit, 
			  struct tophit_s *dhit, int *ret_nseq);
static void main_loop_threaded(struct plan7_s *hmm, SQFILE *sqfp, struct threshold_s *thresh,
			       int do_forward, int do_null2, int do_xnu, int num_threads,
			       struct histogram_s *histogram, struct tophit_s *ghit, 
			       struct tophit_s *dhit, int *ret_nseq);

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
  struct threshold_s *thresh;   /* score/evalue threshold info     */
  
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
					 int do_xnu, int do_forward, int do_null2, 
					 struct threshold_s *thresh,
					 struct tophit_s *ghit, struct tophit_s *dhit, 
					 struct histogram_s *hist, int num_threads);
static void  workpool_stop(struct workpool_s *wpool);
static void  workpool_free(struct workpool_s *wpool);
static void *worker_thread(void *ptr);
#endif /* HMMER_THREADS */



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
  int     descwidth;		/* max width of description */
  int     nreported;		/* # of hits reported in a list            */

  int    Alimit;		/* A parameter limiting output alignments   */
  struct threshold_s thresh;    /* contains all threshold (cutoff) info     */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  int   do_null2;		/* TRUE to adjust scores with null model #2 */
  int   do_forward;		/* TRUE to use Forward() not Viterbi()      */
  int   do_xnu;			/* TRUE to filter sequences thru XNU        */
  int   do_pvm;			/* TRUE to run on Parallel Virtual Machine  */
  int   be_backwards;		/* TRUE to be backwards-compatible in output*/
  int   num_threads;		/* number of worker threads                 */
  int   threads_support;	/* TRUE if threads support compiled in      */
  int   pvm_support;		/* TRUE if PVM support compiled in          */

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  format      = SQFILE_UNKNOWN;	/* default: autodetect seq file format  */
  do_forward  = FALSE;
  do_null2    = TRUE;
  do_xnu      = FALSE;
  do_pvm      = FALSE;  
  Z           = 0;
  be_backwards= FALSE; 

  pvm_support     = FALSE;
  threads_support = FALSE;
  num_threads     = 0;
#ifdef HMMER_THREADS
  num_threads     = ThreadNumber(); 
  threads_support = TRUE;
#endif
#ifdef HMMER_PVM
  pvm_support     = TRUE;
#endif

  Alimit         = INT_MAX;	/* no limit on alignment output       */
  thresh.globE   = 10.0;	/* use a reasonable Eval threshold;   */
  thresh.globT   = -FLT_MAX;	/*   but no bit threshold,            */
  thresh.domT    = -FLT_MAX;	/*   no domain bit threshold,         */
  thresh.domE    = FLT_MAX;     /*   and no domain Eval threshold.    */
  thresh.autocut = CUT_NONE;	/*   and no Pfam cutoffs used         */
  thresh.Z       = 0;           /* Z not preset; use actual # of seqs */

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-A") == 0)        Alimit         = atoi(optarg); 
    else if (strcmp(optname, "-E") == 0)        thresh.globE   = atof(optarg);
    else if (strcmp(optname, "-T") == 0)        thresh.globT   = atof(optarg);
    else if (strcmp(optname, "-Z") == 0)        thresh.Z       = atoi(optarg);
    else if (strcmp(optname, "--compat")  == 0) be_backwards   = TRUE;
    else if (strcmp(optname, "--cpu")     == 0) num_threads    = atoi(optarg);
    else if (strcmp(optname, "--cut_ga")  == 0) thresh.autocut = CUT_GA;
    else if (strcmp(optname, "--cut_nc")  == 0) thresh.autocut = CUT_NC;
    else if (strcmp(optname, "--cut_tc")  == 0) thresh.autocut = CUT_TC;
    else if (strcmp(optname, "--domE")    == 0) thresh.domE    = atof(optarg);
    else if (strcmp(optname, "--domT")    == 0) thresh.domT    = atof(optarg);
    else if (strcmp(optname, "--forward") == 0) do_forward     = TRUE;
    else if (strcmp(optname, "--null2")   == 0) do_null2       = FALSE;
    else if (strcmp(optname, "--pvm")     == 0) do_pvm         = TRUE;
    else if (strcmp(optname, "--xnu")     == 0) do_xnu         = TRUE;
    else if (strcmp(optname, "--informat") == 0) {
      format = String2SeqfileFormat(optarg);
      if (format == SQFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
    }
    else if (strcmp(optname, "-h") == 0) {
      HMMERBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }
  if (argc - optind != 2)
    Die("Incorrect number of arguments.\n%s\n", usage);

  hmmfile = argv[optind++];
  seqfile = argv[optind++]; 
  
  if (do_pvm && ! pvm_support) 
    Die("PVM support is not compiled into your HMMER software; --pvm doesn't work.");
  if (num_threads && ! threads_support)
    Die("POSIX threads support is not compiled into HMMER; --cpu doesn't have any effect");

  /* Try to work around inability to autodetect from a pipe or .gz:
   * assume FASTA format
   */
  if (format == SQFILE_UNKNOWN &&
      (Strparse("^.*\\.gz$", seqfile, 0) || strcmp(seqfile, "-") == 0))
    format = SQFILE_FASTA;

  /*********************************************** 
   * Open sequence database (might be in BLASTDB or current directory)
   ***********************************************/

  if ((sqfp = SeqfileOpen(seqfile, format, "BLASTDB")) == NULL)
    Die("Failed to open sequence database file %s\n%s\n", seqfile, usage);

  /*********************************************** 
   * Open HMM file (might be in HMMERDB or current directory).
   * Read a single HMM from it. (Config HMM, if necessary).
   * Alphabet globals are set by reading the HMM.
   ***********************************************/

  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  if (hmm == NULL) 
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  P7Logoddsify(hmm, !do_forward);

  if (do_xnu && Alphabet_type == hmmNUCLEIC) 
    Die("The HMM is a DNA model, and you can't use the --xnu filter on DNA data");

  /*****************************************************************
   * Set up optional Pfam score thresholds. 
   * Can do this before starting any searches, since we'll only use 1 HMM.
   *****************************************************************/ 

  if (! SetAutocuts(&thresh, hmm)) 
    Die("HMM %s did not contain the GA, TC, or NC cutoffs you needed",
	hmm->name);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  HMMERBanner(stdout, banner);
  printf(   "HMM file:                   %s [%s]\n", hmmfile, hmm->name);
  printf(   "Sequence database:          %s\n", seqfile); 
  if (do_pvm)
    printf( "PVM:                        ACTIVE\n");
  printf(   "per-sequence score cutoff:  ");
  if (thresh.globT == -FLT_MAX) printf("[none]\n");
  else  {
    printf(">= %.1f", thresh.globT);
    if      (thresh.autocut == CUT_GA) printf(" [GA1]\n");
    else if (thresh.autocut == CUT_NC) printf(" [NC1]\n");
    else if (thresh.autocut == CUT_TC) printf(" [TC1]\n");
    else                               printf("\n");
  }
  printf(   "per-domain score cutoff:    ");
  if (thresh.domT == -FLT_MAX) printf("[none]\n");
  else  {
    printf(">= %.1f", thresh.domT);
    if      (thresh.autocut == CUT_GA) printf(" [GA2]\n");
    else if (thresh.autocut == CUT_NC) printf(" [NC2]\n");
    else if (thresh.autocut == CUT_TC) printf(" [TC2]\n");
    else                               printf("\n");
  }
  printf(   "per-sequence Eval cutoff:   ");
  if (thresh.globE == FLT_MAX) printf("[none]\n");
  else                  printf("<= %-10.2g\n", thresh.globE);
    
  printf(   "per-domain Eval cutoff:     ");
  if (thresh.domE == FLT_MAX) printf("[none]\n");
  else                 printf("<= %10.2g\n", thresh.domE);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  /*********************************************** 
   * Search HMM against each sequence
   ***********************************************/

				/* set up structures for storing output  */
  histogram = AllocHistogram(-200, 200, 100);  /* keeps full histogram */
  ghit      = AllocTophits(200);         /* per-seq hits: 200=lumpsize */
  dhit      = AllocTophits(200);         /* domain hits:  200=lumpsize */

  if (pvm_support && do_pvm)
    main_loop_pvm(hmm, sqfp, &thresh, do_forward, do_null2, do_xnu, 
		  histogram, ghit, dhit, &nseq);
  else if (threads_support && num_threads)
    main_loop_threaded(hmm, sqfp, &thresh, do_forward, do_null2, do_xnu, num_threads,
		       histogram, ghit, dhit, &nseq);    
  else
    main_loop_serial(hmm, sqfp, &thresh, do_forward, do_null2, do_xnu, 
		     histogram, ghit, dhit, &nseq);

  /*********************************************** 
   * Process hit lists, produce text output
   ***********************************************/

  /* Set the theoretical EVD curve in our histogram using 
   * calibration in the HMM, if available. 
   */
  if (hmm->flags & PLAN7_STATS)
    ExtremeValueSetHistogram(histogram, hmm->mu, hmm->lambda, 
			     histogram->lowscore, histogram->highscore, 0);
  if (!thresh.Z) thresh.Z = nseq;		/* set Z for good now that we're done. */

  /* Format and report our output 
   */
  /* 1. Report overall sequence hits (sorted on E-value) */
  if (be_backwards) 
    {
      printf("\nQuery HMM: %s|%s|%s\n", 
	     hmm->name, 
	     hmm->flags & PLAN7_ACC  ? hmm->acc  : "",
	     hmm->flags & PLAN7_DESC ? hmm->desc : "");
    }
  else 
    {
      printf("\nQuery HMM:   %s\n", hmm->name);
      printf("Accession:   %s\n", hmm->flags & PLAN7_ACC  ? hmm->acc  : "[none]");
      printf("Description: %s\n", hmm->flags & PLAN7_DESC ? hmm->desc : "[none]");
    }

  if (hmm->flags & PLAN7_STATS)
    printf("  [HMM has been calibrated; E-values are empirical estimates]\n");
  else
    printf("  [No calibration for HMM; E-values are upper bounds]\n");

  FullSortTophits(ghit);
  namewidth = MAX(8, TophitsMaxName(ghit)); /* cannot truncate name. */
  descwidth = MAX(52-namewidth, 11);/* may truncate desc, but need strlen("Description") */

  printf("\nScores for complete sequences (score includes all domains):\n");
  printf("%-*s %-*s %7s %10s %3s\n", namewidth, "Sequence", descwidth, "Description", "Score", "E-value", " N ");
  printf("%-*s %-*s %7s %10s %3s\n", namewidth, "--------", descwidth, "-----------", "-----", "-------", "---");
  for (i = 0, nreported = 0; i < ghit->num; i++)
    {
      char *safedesc;
      GetRankedHit(ghit, i, 
		   &pvalue, &sc, NULL, NULL,
		   &name, NULL, &desc,
		   NULL, NULL, NULL,               /* sequence positions */
		   NULL, NULL, NULL,               /* HMM positions      */
		   NULL, &ndom,	                   /* domain info        */
		   NULL);	                   /* alignment info     */
      evalue = pvalue * (double) thresh.Z;

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

      if (evalue <= thresh.globE && sc >= thresh.globT) {
	printf("%-*s %-*.*s %7.1f %10.2g %3d\n", 
	       namewidth, name, 
	       descwidth, descwidth, safedesc != NULL ? safedesc : "",
	       sc, evalue, ndom);
	nreported++;
      }
      free(safedesc);
    }
  if (nreported == 0) printf("\t[no hits above thresholds]\n");


  /* 2. Report domain hits (also sorted on E-value) */
  FullSortTophits(dhit);
  namewidth = MAX(8, TophitsMaxName(dhit));

  printf("\nParsed for domains:\n");
  printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
	 namewidth, "Sequence", "Domain ", "seq-f", "seq-t", "hmm-f", "hmm-t", "score", "E-value");
  printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
	 namewidth, "--------", "-------", "-----", "-----", "-----", "-----", "-----", "-------");
      
  for (i = 0, nreported = 0; i < dhit->num; i++)
    {
      GetRankedHit(dhit, i, 
		   &pvalue, &sc, &motherp, &mothersc,
		   &name, NULL, NULL,
		   &sqfrom, &sqto, &sqlen,            /* seq position info  */
		   &hmmfrom, &hmmto, NULL,            /* HMM position info  */
		   &domidx, &ndom,                    /* domain info        */
		   NULL);	                      /* alignment info     */
      evalue = pvalue * (double) thresh.Z;

      if (motherp * (double) thresh.Z > thresh.globE || mothersc < thresh.globT) 
	continue;
      else if (evalue <= thresh.domE && sc >= thresh.domT) {
	printf("%-*s %3d/%-3d %5d %5d %c%c %5d %5d %c%c %7.1f %8.2g\n",
	       namewidth, name, 
	       domidx, ndom,
	       sqfrom, sqto, 
	       sqfrom == 1 ? '[' : '.', sqto == sqlen ? ']' : '.',
	       hmmfrom, hmmto,
	       hmmfrom == 1 ? '[':'.', hmmto == hmm->M ? ']' : '.',
	       sc, evalue);
	nreported++;
      }
    }
  if (nreported == 0) printf("\t[no hits above thresholds]\n");


  /* 3. Alignment output, also by domain.
   *    dhits is already sorted and namewidth is set, from above code.
   *    Number of displayed alignments is limited by Alimit parameter;
   *    also by domE (evalue threshold), domT (score theshold).
   */
  if (Alimit != 0)
    {
      printf("\nAlignments of top-scoring domains:\n");
      for (i = 0, nreported = 0; i < dhit->num; i++)
	{
	  if (nreported == Alimit) break; /* limit to Alimit output alignments */
	  GetRankedHit(dhit, i, 
		       &pvalue, &sc, &motherp, &mothersc,
		       &name, NULL, NULL,
		       &sqfrom, &sqto, &sqlen,            /* seq position info  */
		       &hmmfrom, &hmmto, NULL,            /* HMM position info  */
		       &domidx, &ndom,                    /* domain info        */
		       &ali);	                      /* alignment info     */
	  evalue = pvalue * (double) thresh.Z;

	  if (motherp * (double) thresh.Z > thresh.globE || mothersc < thresh.globT) 
	    continue;
	  else if (evalue <= thresh.domE && sc >= thresh.domT) 
	    {
	      printf("%s: domain %d of %d, from %d to %d: score %.1f, E = %.2g\n", 
		     name, domidx, ndom, sqfrom, sqto, sc, evalue);
	      PrintFancyAli(stdout, ali);
	      nreported++;
	    }
	}
      if (nreported == 0)      printf("\t[no hits above thresholds]\n");
      if (nreported == Alimit) printf("\t[output cut off at A = %d top alignments]\n", Alimit);
    }

  /* 4. Histogram output */
  printf("\nHistogram of all scores:\n");
  PrintASCIIHistogram(stdout, histogram);

  /* 5. Tophits summaries, while developing...
   */
  printf("\nTotal sequences searched: %d\n", nseq);
  printf("\nWhole sequence top hits:\n");
  TophitsReport(ghit, thresh.globE, nseq);
  printf("\nDomain top hits:\n");
  TophitsReport(dhit, thresh.domE, nseq);

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
 *           thresh     - score/evalue threshold info
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
main_loop_serial(struct plan7_s *hmm, SQFILE *sqfp, struct threshold_s *thresh, int do_forward,
		 int do_null2, int do_xnu,
		 struct histogram_s *histogram, 
		 struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nseq)
{
  struct dpmatrix_s *mx;        /* DP matrix, growable                     */
  struct p7trace_s *tr;         /* traceback                               */
  char   *seq;                  /* target sequence                         */
  unsigned char   *dsq;	        /* digitized target sequence               */
  SQINFO sqinfo;		/* optional info for seq                   */
  float  sc;	        	/* score of an HMM search                  */
  double pvalue;		/* pvalue of an HMM score                  */
  double evalue;		/* evalue of an HMM score                  */
  int    nseq;			/* number of sequences searched            */
 
  /* Create a DP matrix; initially only two rows big, but growable;
   * we overalloc by 25 rows (L dimension) when we grow; not growable
   * in model dimension, since we know the hmm size
   */
  mx = CreatePlan7Matrix(1, hmm->M, 25, 0); 

  nseq = 0;
  while (ReadSeq(sqfp, sqfp->format, &seq, &sqinfo))
    {
      /* Silently skip length 0 seqs. 
       * What, you think this doesn't occur? Welcome to genomics, 
       * young grasshopper.
       */
      if (sqinfo.len == 0) continue;

      nseq++;
      dsq = DigitizeSequence(seq, sqinfo.len);
      
      if (do_xnu && Alphabet_type == hmmAMINO) XNU(dsq, sqinfo.len);
      
      /* 1. Recover a trace by Viterbi.
       *    In extreme cases, the alignment may be literally impossible;
       *    in which case, the score comes out ridiculously small (but not
       *    necessarily <= -INFTY, because we're not terribly careful
       *    about underflow issues), and tr will be returned as NULL.
       */
      if (P7ViterbiSpaceOK(sqinfo.len, hmm->M, mx))
	sc = P7Viterbi(dsq, sqinfo.len, hmm, mx, &tr);
      else
	sc = P7SmallViterbi(dsq, sqinfo.len, hmm, mx, &tr);

      /* 2. If we're using Forward scores, calculate the
       *    whole sequence score; this overrides anything
       *    PostprocessSignificantHit() is going to do to the per-seq score.
       */
      if (do_forward) {
	sc  = P7Forward(dsq, sqinfo.len, hmm, NULL);
	if (do_null2)   sc -= TraceScoreCorrection(hmm, tr, dsq); 
      }

#if DEBUGLEVEL >= 2
      P7PrintTrace(stdout, tr, hmm, dsq); 
#endif

      /* 2. Store score/pvalue for global alignment; will sort on score,
       *    which in hmmsearch is monotonic with E-value. 
       *    Keep all domains in a significant sequence hit.
       *    We can only make a lower bound estimate of E-value since
       *    we don't know the final value of nseq yet, so the list
       *    of hits we keep in memory is >= the list we actually
       *    output. 
       */
      pvalue = PValue(hmm, sc);
      evalue = thresh->Z ? (double) thresh->Z * pvalue : (double) nseq * pvalue;
      if (sc >= thresh->globT && evalue <= thresh->globE) 
	{
	  sc = PostprocessSignificantHit(ghit, dhit, 
					 tr, hmm, dsq, sqinfo.len,
					 sqinfo.name, 
					 sqinfo.flags & SQINFO_ACC  ? sqinfo.acc  : NULL, 
					 sqinfo.flags & SQINFO_DESC ? sqinfo.desc : NULL, 
					 do_forward, sc,
					 do_null2,
					 thresh,
					 FALSE); /* FALSE-> not hmmpfam mode, hmmsearch mode */
	}
      SQD_DPRINTF2(("AddToHistogram: %s\t%f\n", sqinfo.name, sc));
      AddToHistogram(histogram, sc);
      FreeSequence(seq, &sqinfo); 
      P7FreeTrace(tr);
      free(dsq);
    }

  FreePlan7Matrix(mx);
  *ret_nseq = nseq;
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
 *           thresh     - score/evalue threshold information
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
main_loop_pvm(struct plan7_s *hmm, SQFILE *sqfp, struct threshold_s *thresh, int do_forward,
	      int do_null2, int do_xnu, struct histogram_s *histogram, 
	      struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nseq)
{
  char *seq;                    /* target sequence */
  unsigned char *dsq;           /* digitized target seq */
  SQINFO sqinfo;                /* optional info about target seq */
  int   master_tid;		/* master's (my) PVM TID */
  int  *slave_tid;              /* array of slave TID's  */
  int   nslaves;		/* number of slaves      */
  int   code;			/* status code rec'd from a slave */
  int   nseq;			/* number of sequences searched */
  int   sent_trace;		/* TRUE if slave gave us a trace */
  unsigned char **dsqlist;      /* remember what seqs slaves are doing */
  char **namelist;              /* remember what seq names slaves are doing */
  char **acclist ;              /* remember what seq accessions slaves are doing */
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
#if DEBUGLEVEL >= 1
  pvm_catchout(stderr);		/* catch output for debugging */
#endif
  SQD_DPRINTF1(("Spawning slaves...\n"));
  PVMSpawnSlaves("hmmsearch-pvm", &slave_tid, &nslaves);
  SQD_DPRINTF1(("Spawned a total of %d slaves...\n", nslaves));
 
  /* Initialize the slaves by broadcast.
   */
  SQD_DPRINTF1(("Broadcasting to %d slaves...\n", nslaves));
  pvm_initsend(PvmDataDefault);
  pvm_pkfloat(&(thresh->globT), 1, 1);
  pvm_pkdouble(&(thresh->globE), 1, 1);
  pvm_pkint(&(thresh->Z),          1, 1);
  pvm_pkint(&do_forward, 1, 1);
  pvm_pkint(&do_null2,   1, 1);
  pvm_pkint(&Alphabet_type, 1, 1);
  PVMPackHMM(hmm);
  pvm_mcast(slave_tid, nslaves, HMMPVM_INIT);
  SQD_DPRINTF1(("Slaves should be ready...\n"));

  /* Confirm slaves' OK status.
   */
  PVMConfirmSlaves(slave_tid, nslaves);
  SQD_DPRINTF1(("Slaves confirm that they're ok...\n"));
  
  /* Alloc arrays for remembering what seq each
   * slave was working on.
   */
  namelist = MallocOrDie(sizeof(char *) * nslaves);
  acclist  = MallocOrDie(sizeof(char *) * nslaves);
  desclist = MallocOrDie(sizeof(char *) * nslaves);
  dsqlist  = MallocOrDie(sizeof(unsigned char *) * nslaves);
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
      if (do_xnu && Alphabet_type == hmmAMINO) XNU(dsq, sqinfo.len);

      pvm_initsend(PvmDataDefault);
      pvm_pkint(&nseq, 1, 1);
      pvm_pkint(&(sqinfo.len), 1, 1);
      pvm_pkbyte(dsq, sqinfo.len+2, 1);
      pvm_send(slave_tid[nseq], HMMPVM_WORK);
      SQD_DPRINTF1(("sent a dsq : %d bytes\n", sqinfo.len+2));

      namelist[nseq] = Strdup(sqinfo.name);
      acclist[nseq]  = (sqinfo.flags & SQINFO_ACC)  ? Strdup(sqinfo.acc)  : NULL;
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
	  sc = PostprocessSignificantHit(ghit, dhit, 
					 tr, hmm, dsqlist[slaveidx], lenlist[slaveidx],
					 namelist[slaveidx], acclist[slaveidx], desclist[slaveidx],
					 do_forward, sc,
					 do_null2,
					 thresh,
					 FALSE); /* FALSE-> not hmmpfam mode, hmmsearch mode */
	  P7FreeTrace(tr);
	}
      AddToHistogram(histogram, sc);

				/* record seq info for seq we just sent */
      free(namelist[slaveidx]);
      if (acclist[slaveidx]  != NULL) free(acclist[slaveidx]);
      if (desclist[slaveidx] != NULL) free(desclist[slaveidx]);
      free(dsqlist[slaveidx]);

      dsqlist[slaveidx]  = dsq;
      namelist[slaveidx] = Strdup(sqinfo.name);
      acclist[slaveidx]  = (sqinfo.flags & SQINFO_ACC)  ? Strdup(sqinfo.acc)  : NULL;
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
	  sc = PostprocessSignificantHit(ghit, dhit, 
					 tr, hmm, dsqlist[slaveidx], lenlist[slaveidx],
					 namelist[slaveidx], acclist[slaveidx], desclist[slaveidx],
					 do_forward, sc,
					 do_null2,
					 thresh,
					 FALSE); /* FALSE-> not hmmpfam mode, hmmsearch mode */
	  P7FreeTrace(tr);
	}
      SQD_DPRINTF2(("AddToHistogram: %s\t%f\n", namelist[slaveidx], sc));
      AddToHistogram(histogram, sc);

				/* free seq info */
      free(namelist[slaveidx]);
      if (acclist[slaveidx]  != NULL) free(acclist[slaveidx]);
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
  free(acclist);
  free(desclist);
  free(lenlist);
  pvm_exit();
  *ret_nseq = nseq;
  return;
}
#else /* HMMER_PVM off, no PVM support, dummy function: */
static void
main_loop_pvm(struct plan7_s *hmm, SQFILE *sqfp, struct threshold_s *thresh, int do_forward,
	      int do_null2, int do_xnu, struct histogram_s *histogram, 
	      struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nseq)
{
  Die("No PVM support");
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
/* Function: main_loop_threaded
 * Date:     SRE, Wed Feb 27 11:38:11 2002 [St. Louis]
 *
 * Purpose:  Search an HMM against a sequence database.
 *           main loop for the threaded version.
 *           
 *           In:   HMM and open sqfile, plus options
 *           Out:  histogram, global hits list, domain hits list, nseq.
 *
 * Args:     hmm        - the HMM to search with. 
 *           sqfp       - open SQFILE for sequence database
 *           thresh     - score/evalue threshold info
 *           do_forward - TRUE to score using Forward()        
 *           do_null2   - TRUE to use ad hoc null2 score correction
 *           do_xnu     - TRUE to apply XNU mask
 *           num_threads- number of worker threads to start, >=1
 *           histogram  - RETURN: score histogram
 *           ghit       - RETURN: ranked global scores
 *           dhit       - RETURN: ranked domain scores
 *           ret_nseq   - RETURN: actual number of seqs searched
 *           
 * Returns:  (void)
 */
static void
main_loop_threaded(struct plan7_s *hmm, SQFILE *sqfp, struct threshold_s *thresh, int do_forward,
		   int do_null2, int do_xnu, int num_threads,
		   struct histogram_s *histogram, 
		   struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nseq)
{
  struct workpool_s *wpool;	/* pool of worker threads                  */
  int    nseq;			/* number of sequences searched            */
 
  wpool = workpool_start(hmm, sqfp, do_xnu, do_forward, do_null2, thresh,
			 ghit, dhit, histogram, num_threads);
  workpool_stop(wpool);
  nseq = wpool->nseq;
  workpool_free(wpool);

  *ret_nseq = nseq;
  return;
}
/* Function: workpool_start()
 * Date:     SRE, Mon Oct  5 16:44:53 1998
 *
 * Purpose:  Initialize a workpool_s structure, and return it.
 *
 * Args:     sqfp       - open sequence file, at start
 *           do_xnu     - TRUE to apply XNU filter
 *           do_forward - TRUE to score using Forward
 *           do_null2   - TRUE to apply null2 ad hoc correction
 *           thresh     - score/evalue threshold info
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
	       int do_forward, int do_null2, struct threshold_s *thresh,
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
  wpool->thresh     = thresh;

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

  /* Create slave threads. See comments in hmmcalibrate.c at this
   * step, regarding concurrency, system scope, and portability
   * amongst various UNIX implementations of pthreads.
   */
  pthread_attr_init(&attr);
#ifndef __sgi
#ifdef HAVE_PTHREAD_ATTR_SETSCOPE
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
#endif
#endif
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
  unsigned char     *dsq;       /* digitized sequence              */
  struct dpmatrix_s *mx;        /* growable DP matrix              */
  struct p7trace_s  *tr;        /* traceback from an alignment     */
  float  sc;			/* score of an alignment           */
  int    rtn;			/* a return code from pthreads lib */
  double pvalue;		/* P-value of score                */
  double evalue;		/* E-value of score                */

  wpool = (struct workpool_s *) ptr;

  /* Init with a small DP matrix; we'll grow in the sequence dimension
   * overalloc'ing by 25 rows (residues).
   */
  mx = CreatePlan7Matrix(1, wpool->hmm->M, 25, 0);
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
	FreePlan7Matrix(mx);
	pthread_exit(NULL);
      }
    SQD_DPRINTF1(("a thread is working on %s\n", sqinfo.name));
    wpool->nseq++;
				/* release the lock */
    if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

    if (sqinfo.len == 0) continue; /* silent skip of len=0 seqs (wormpep!?!) */

    dsq = DigitizeSequence(seq, sqinfo.len);
    if (wpool->do_xnu) XNU(dsq, sqinfo.len);
      
    /* 1. Recover a trace by Viterbi.
     */
    if (P7ViterbiSpaceOK(sqinfo.len, wpool->hmm->M, mx))
      sc = P7Viterbi(dsq, sqinfo.len, wpool->hmm, mx, &tr);
    else
      sc = P7SmallViterbi(dsq, sqinfo.len, wpool->hmm, mx, &tr);

    /* 2. If we're using Forward scores, do another DP
     *    to get it; else, we already have a Viterbi score
     *    in sc.
     */
    if (wpool->do_forward) {
      sc  = P7Forward(dsq, sqinfo.len, wpool->hmm, NULL);
      if (wpool->do_null2) sc -= TraceScoreCorrection(wpool->hmm, tr, dsq);
    }

    /* 3. Save the output in tophits and histogram structures, after acquiring a lock
     */
    if ((rtn = pthread_mutex_lock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    SQD_DPRINTF1(("seq %s scores %f\n", sqinfo.name, sc));

    pvalue = PValue(wpool->hmm, sc);
    evalue = wpool->thresh->Z ? (double) wpool->thresh->Z * pvalue : (double) wpool->nseq * pvalue;
 
    if (sc >= wpool->thresh->globT && evalue <= wpool->thresh->globE) 
      { 
	sc = PostprocessSignificantHit(wpool->ghit, wpool->dhit, 
				       tr, wpool->hmm, dsq, sqinfo.len,
				       sqinfo.name, 
				       sqinfo.flags & SQINFO_ACC  ? sqinfo.acc  : NULL, 
				       sqinfo.flags & SQINFO_DESC ? sqinfo.desc : NULL, 
				       wpool->do_forward, sc,
				       wpool->do_null2,
				       wpool->thresh,
				       FALSE); /* FALSE-> not hmmpfam mode, hmmsearch mode */
      }
    SQD_DPRINTF2(("AddToHistogram: %s\t%f\n", sqinfo.name, sc));
    AddToHistogram(wpool->hist, sc);
    if ((rtn = pthread_mutex_unlock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

    P7FreeTrace(tr);
    FreeSequence(seq, &sqinfo);
    free(dsq);
  } /* end 'infinite' loop over seqs in this thread */
}
#else /*HMMER_THREADS off; no threads support; dummy stub: */
static void
main_loop_threaded(struct plan7_s *hmm, SQFILE *sqfp, struct threshold_s *thresh, int do_forward,
		   int do_null2, int do_xnu, int num_threads,
		   struct histogram_s *histogram, 
		   struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nseq)
{
  Die("No threads support");
}
#endif /* HMMER_THREADS */

