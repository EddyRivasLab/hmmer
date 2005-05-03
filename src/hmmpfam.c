/* hmmpfam.c
 * Search a single sequence against an HMM database.
 * Conditionally includes PVM parallelization when HMMER_PVM is defined
 *    at compile time; hmmpfam --pvm runs the PVM version.
 *    
 * SRE, Mon Aug 25 17:03:14 1997 [Denver] 
 * SVN $Id$
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

#include "plan7.h"              /* plan 7 profile HMM structure         */
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */

static char banner[] = "hmmpfam - search one or more sequences against HMM database";

static char usage[]  = "\
Usage: hmmpfam [-options] <hmm database> <sequence file or database>\n\
  Available options are:\n\
   -h        : help; print brief help on version and usage\n\
   -n        : nucleic acid models/sequence (default protein)\n\
   -A <n>    : sets alignment output limit to <n> best domain alignments\n\
   -E <x>    : sets E value cutoff (globE) to <x>; default 10\n\
   -T <x>    : sets T bit threshold (globT) to <x>; no threshold by default\n\
   -Z <n>    : sets Z (# models) for E-value calculation\n\
";

static char experts[] = "\
   --acc         : use HMM accession numbers instead of names in output\n\
   --compat      : make best effort to use last version's output style\n\
   --cpu <n>     : run <n> threads in parallel (if threaded)\n\
   --cut_ga      : use Pfam GA gathering threshold cutoffs\n\
   --cut_nc      : use Pfam NC noise threshold cutoffs\n\
   --cut_tc      : use Pfam TC trusted threshold cutoffs\n\
   --domE <x>    : sets domain Eval cutoff (2nd threshold) to <x>\n\
   --domT <x>    : sets domain T bit thresh (2nd threshold) to <x>\n\
   --forward     : use the full Forward() algorithm instead of Viterbi\n\
   --informat <s>: sequence file is in format <s>\n\
   --null2       : turn OFF the post hoc second null model\n\
   --pvm         : run on a PVM (Parallel Virtual Machine) cluster\n\
   --xnu         : turn ON XNU filtering of query protein sequence\n\
\n";


static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-n",        TRUE,  sqdARG_NONE },
  { "-A",        TRUE,  sqdARG_INT  },  
  { "-E",        TRUE,  sqdARG_FLOAT},  
  { "-T",        TRUE,  sqdARG_FLOAT},  
  { "-Z",        TRUE,  sqdARG_INT  },
  { "--acc",     FALSE, sqdARG_NONE },
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

static void main_loop_serial(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo, 
			     struct threshold_s *thresh, int do_xnu, int do_forward, int do_null2,
			     struct tophit_s *ghit, struct tophit_s *dhit, int *nhmm);
static void main_loop_pvm(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo, 
			  struct threshold_s *thresh, int do_xnu, int do_forward, int do_null2,
			  struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm);
static void main_loop_threaded(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo, 
			       struct threshold_s *thresh, int do_xnu, 
			       int do_forward, int do_null2, int num_threads,
			       struct tophit_s *ghit, struct tophit_s *dhit, int *nhmm);

#ifdef HMMER_THREADS
/* POSIX threads version:
 * the threads share a workpool_s structure amongst themselves,
 * for obtaining locks on input HMM file and output histogram and
 * tophits structures.
 */
struct workpool_s {
  /* Shared configuration resources that don't change:
   */
  char          *hmmfile;       /* name of HMM file                */
  unsigned char *dsq;           /* digitized query sequence        */
  char  *seqname;               /* sequence name                   */
  int    L;			/* length of dsq                   */
  int    do_forward;		/* TRUE to score using Forward     */
  int    do_null2;		/* TRUE to apply null2 correction  */
  struct threshold_s *thresh;   /* score/evalue cutoff information */ 
  
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

static struct workpool_s *workpool_start(char *hmmfile, HMMFILE *hmmfp, 
					 unsigned char *dsq, char *seqname, int L,
					 int do_forward, int do_null2, 
					 struct threshold_s *thresh,
					 struct tophit_s *ghit, struct tophit_s *dhit, 
					 int num_threads);
static void  workpool_stop(struct workpool_s *wpool);
static void  workpool_free(struct workpool_s *wpool);
static void *worker_thread(void *ptr);
#endif /* HMMER_THREADS */




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
  char   *name, *acc, *desc;    /* hit HMM name, accession, description            */
  int     hmmlen;		/* length of HMM hit                       */
  int     nhmm;			/* number of HMMs searched                 */
  int     domidx;		/* number of this domain                   */
  int     ndom;			/* total # of domains in this seq          */
  int     namewidth;		/* max width of printed HMM name           */
  int     descwidth;		/* max width of printed description        */

  int    Alimit;		/* A parameter limiting output alignments   */
  struct threshold_s thresh;    /* contains all threshold (cutoff) info     */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  int   do_forward;		/* TRUE to use Forward() not Viterbi()      */
  int   do_nucleic;		/* TRUE to do DNA/RNA instead of protein    */
  int   do_null2;		/* TRUE to adjust scores with null model #2 */
  int   do_pvm;			/* TRUE to run on PVM                       */
  int   do_xnu;			/* TRUE to do XNU filtering                 */
  int   be_backwards;		/* TRUE to be backwards-compatible in output*/
  int   show_acc;		/* TRUE to sub HMM accessions for names     */
  int   i;
  int   nreported;

  int   num_threads;            /* number of worker threads */   
  int   threads_support;	/* TRUE if threads support compiled in */
  int   pvm_support;		/* TRUE if PVM support compiled in     */

  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  format      = SQFILE_UNKNOWN;	/* default: autodetect format w/ Babelfish */
  do_forward  = FALSE;
  do_nucleic  = FALSE;
  do_null2    = TRUE;
  do_pvm      = FALSE;
  do_xnu      = FALSE;
  be_backwards= FALSE; 
  show_acc    = FALSE;
  
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

  Alimit         = INT_MAX;	/* no limit on alignment output     */
  thresh.globE   = 10.0;	/* use a reasonable Eval threshold; */
  thresh.globT   = -FLT_MAX;	/*   but no bit threshold,          */
  thresh.domT    = -FLT_MAX;	/*   no domain bit threshold,       */
  thresh.domE    = FLT_MAX;     /*   and no domain Eval threshold.  */
  thresh.autocut = CUT_NONE;	/*   and no Pfam cutoffs used.      */
  thresh.Z       = 0;		/* Z not preset, so determined by # of HMMs */

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n")        == 0) do_nucleic     = TRUE; 
    else if (strcmp(optname, "-A")        == 0) Alimit         = atoi(optarg);  
    else if (strcmp(optname, "-E")        == 0) thresh.globE   = atof(optarg);
    else if (strcmp(optname, "-T")        == 0) thresh.globT   = atof(optarg);
    else if (strcmp(optname, "-Z")        == 0) thresh.Z       = atoi(optarg);
    else if (strcmp(optname, "--acc")     == 0) show_acc       = TRUE;
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
    else if (strcmp(optname, "-h")      == 0) {
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
    Die("PVM support is not compiled into HMMER; --pvm doesn't work.");
  if (num_threads && ! threads_support)
    Die("Posix threads support is not compiled into HMMER; --cpu doesn't have any effect");

  /* Try to work around inability to autodetect from a pipe or .gz:
   * assume FASTA format
   */
  if (format == SQFILE_UNKNOWN &&
      (Strparse("^.*\\.gz$", seqfile, 0) || strcmp(seqfile, "-") == 0))
    format = SQFILE_FASTA;

  /*********************************************** 
   * Open sequence database (must be in curr directory);
   * get target sequence.
   ***********************************************/

  if (do_nucleic) SetAlphabet(hmmNUCLEIC);
  else            SetAlphabet(hmmAMINO);

  if (do_nucleic && do_xnu) 
    Die("You can't use -n and --xnu together: I can't xnu DNA data.");

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

  HMMERBanner(stdout, banner);
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
      if (pvm_support && do_pvm)
	main_loop_pvm(hmmfile, hmmfp, seq, &sqinfo, 
		      &thresh, do_xnu, do_forward, do_null2,
		      ghit, dhit, &nhmm);
      else if (threads_support && num_threads)
	main_loop_threaded(hmmfile, hmmfp, seq, &sqinfo, 
			   &thresh, do_xnu, do_forward, do_null2, num_threads,
			   ghit, dhit, &nhmm);
      else 
	main_loop_serial(hmmfile, hmmfp, seq, &sqinfo, 
			 &thresh, do_xnu, do_forward, do_null2, 
			 ghit, dhit, &nhmm);

				/* set Z for good now that we're done */
      if (!thresh.Z) thresh.Z = nhmm;	

      /* 2. (Done searching all HMMs for this query seq; start output)
       *    Report the overall sequence hits, sorted by significance.
       */
      if (be_backwards)
	{
	  printf("Query:  %s  %s\n", sqinfo.name, 
		 sqinfo.flags & SQINFO_DESC ? sqinfo.desc : "");
	}
      else
	{
	  printf("\nQuery sequence: %s\n", sqinfo.name);
	  printf("Accession:      %s\n", sqinfo.flags &SQINFO_ACC ? sqinfo.acc  : "[none]");
	  printf("Description:    %s\n", sqinfo.flags &SQINFO_DESC? sqinfo.desc : "[none]");
	}
      /* We'll now sort the global hit list by evalue... 
       * (not score! that was bug #12. in hmmpfam, score and evalue are not 
       *  monotonic.)
       */
      FullSortTophits(ghit);	
      namewidth = MAX(8, TophitsMaxName(ghit)); /* must print whole name, no truncation  */
      descwidth = MAX(52-namewidth, 11);        /* may truncate desc, but avoid neg len! */

      printf("\nScores for sequence family classification (score includes all domains):\n");
      printf("%-*s %-*s %7s %10s %3s\n", namewidth, "Model", descwidth, "Description", "Score", "E-value", " N ");
      printf("%-*s %-*s %7s %10s %3s\n", namewidth, "--------", descwidth, "-----------", "-----", "-------", "---");
      for (i = 0, nreported = 0; i < ghit->num; i++)
	{
	  char *safedesc;
	  GetRankedHit(ghit, i, 
		       &pvalue, &sc, NULL, NULL,
		       &name, &acc, &desc,
		       NULL, NULL, NULL,           /* seq positions */
		       NULL, NULL, NULL,           /* HMM positions */
		       NULL, &ndom,                /* domain info   */
		       NULL);	                   /* alignment info*/

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

	  /* sneaky trick warning:
	   * if we're using dynamic Pfam score cutoffs (GA, TC, NC),
	   * then the list of hits is already correct and does not
	   * need any score cutoffs. Unset the thresholds. They'll
           * be reset in the main_loop if we still have sequences
           * to process.
	   */
	  if (thresh.autocut != CUT_NONE) {
	    thresh.globE = thresh.domE = FLT_MAX;
	    thresh.globT = thresh.domT = -FLT_MAX;
	  }

	  if (evalue <= thresh.globE && sc >= thresh.globT) 
	    {
	      printf("%-*s %-*.*s %7.1f %10.2g %3d\n", 
		     namewidth, 
		     (show_acc && acc != NULL) ?  acc : name,
		     descwidth, descwidth, safedesc != NULL ? safedesc : "",
		     sc, evalue, ndom);
	      nreported++;
	    }
	  free(safedesc);
	}
      if (nreported == 0) printf("\t[no hits above thresholds]\n");

      /* 3. Report domain hits (sorted on sqto coordinate)
       */
      FullSortTophits(dhit);
      namewidth = MAX(8, TophitsMaxName(dhit)); /* must print whole name, no truncation  */

      printf("\nParsed for domains:\n");
      printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
	     namewidth, "Model", "Domain ", "seq-f", "seq-t", "hmm-f", "hmm-t", "score", "E-value");
      printf("%-*s %7s %5s %5s    %5s %5s    %7s %8s\n",
	     namewidth, "--------", "-------", "-----", "-----", "-----", "-----", "-----", "-------");
      
      for (i = 0, nreported = 0; i < dhit->num; i++)
	{
	  GetRankedHit(dhit, i, 
		       &pvalue, &sc, &motherp, &mothersc,
		       &name, &acc, NULL,
		       &sqfrom, &sqto, NULL, 
		       &hmmfrom, &hmmto, &hmmlen, 
		       &domidx, &ndom,
		       NULL);
	  evalue = pvalue * (double) thresh.Z;
	  
				/* Does the "mother" (complete) sequence satisfy global thresholds? */
	  if (motherp * (double)thresh. Z > thresh.globE || mothersc < thresh.globT) 
	    continue;
	  else if (evalue <= thresh.domE && sc >= thresh.domT) {
	    printf("%-*s %3d/%-3d %5d %5d %c%c %5d %5d %c%c %7.1f %8.2g\n",
		   namewidth, 
		   (show_acc && acc != NULL) ?  acc : name,
		   domidx, ndom,
		   sqfrom, sqto, 
		   sqfrom == 1 ? '[' : '.', sqto == sqinfo.len ? ']' : '.',
		   hmmfrom, hmmto,
		   hmmfrom == 1 ? '[':'.',  hmmto == hmmlen ? ']' : '.',
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
			   &name, &acc, NULL,
			   &sqfrom, &sqto, NULL,              /* seq position info  */
			   &hmmfrom, &hmmto, &hmmlen,         /* HMM position info  */
			   &domidx, &ndom,                    /* domain info        */
			   &ali);	                      /* alignment info     */
	      evalue = pvalue * (double) thresh.Z;

	      if (motherp * (double) thresh.Z > thresh.globE || mothersc < thresh.globT) 
		continue;
	      else if (evalue <= thresh.domE && sc >= thresh.domT) 
		{
		  printf("%s: domain %d of %d, from %d to %d: score %.1f, E = %.2g\n", 
			 (show_acc && acc != NULL) ?  acc : name,
			 domidx, ndom, sqfrom, sqto, sc, evalue);
		  PrintFancyAli(stdout, ali);
		  nreported++;
		}
	    }
	  if (nreported == 0)      printf("\t[no hits above thresholds]\n");
	  if (nreported == Alimit) printf("\t[output cut off at A = %d top alignments]\n", Alimit);
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

  return 0;
}


/* Function: main_loop_serial()
 * Date:     SRE, Fri Aug  7 13:46:48 1998 [St. Louis]
 *
 * Purpose:  Search a sequence against an HMM database;
 *           main loop for the serial (non-PVM, non-threads)
 *           version. 
 *           
 *           On return, ghit and dhit contain info for all hits
 *           that satisfy the set thresholds. If an evalue
 *           cutoff is used at all, the lists will be overestimated --
 *           because the evalue will be underestimated until
 *           we know the final Z. (Thus the main program must recheck
 *           thresholds before printing any results.) If only
 *           score cutoffs are used, then the lists are correct,
 *           and may be printed exactly as they come (after 
 *           appropriate sorting, anyway). This is especially
 *           important for dynamic thresholding using Pfam
 *           score cutoffs -- the main caller cannot afford to
 *           rescan the HMM file just to get the GA/TC/NC cutoffs
 *           back out for each HMM, and neither do I want to
 *           burn the space to store them as I make a pass thru
 *           Pfam.
 *
 * Args:     hmmfile - name of HMM file
 *           hmmfp   - open HMM file (and at start of file)
 *           dsq     - digitized sequence
 *           sqinfo  - ptr to SQINFO optional info for dsq
 *           thresh  - score/evalue threshold information
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
		 struct threshold_s *thresh, int do_xnu, int do_forward, int do_null2,
		 struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm) 
{
  unsigned char     *dsq;       /* digitized sequence                      */
  int                nhmm;	/* number of HMMs searched                 */
  struct plan7_s    *hmm;       /* current HMM to search with              */ 
  struct p7trace_s  *tr;	/* traceback of alignment                  */
  struct dpmatrix_s *mx;        /* growable DP matrix                      */
  float   sc;                   /* an alignment score                      */ 
  double  pvalue;		/* pvalue of an HMM score                  */
  double  evalue;		/* evalue of an HMM score                  */
    
  tr = NULL; 
  
  /* Prepare sequence.
   */
  SQD_DPRINTF1(("main_loop_serial:\n"));

  dsq = DigitizeSequence(seq, sqinfo->len);
  if (do_xnu && Alphabet_type == hmmAMINO) XNU(dsq, sqinfo->len);

  /* 
   * We'll create for at least N=300xM=300, and thus consume at least 1 MB,
   * regardless of RAMLIMIT -- this helps us avoid reallocating some weird
   * asymmetric matrices.
   * 
   * We're growable in both M and N, because inside of P7SmallViterbi,
   * we're going to be calling P7Viterbi on subpieces that vary in size,
   * and for different models.
   */
  mx = CreatePlan7Matrix(300, 300, 25, 25);

  nhmm = 0;
  while (HMMFileRead(hmmfp, &hmm)) {
    if (hmm == NULL) 
      Die("HMM file %s may be corrupt or in incorrect format; parse failed", hmmfile);
    SQD_DPRINTF1(("   ... working on HMM %s\n", hmm->name));
    SQD_DPRINTF1(("   ... mx is now M=%d by N=%d\n", mx->maxM, mx->maxN));

    if (! SetAutocuts(thresh, hmm)) 
      Die("HMM %s did not contain the GA, TC, or NC cutoffs you needed",
	  hmm->name);

    /* Score sequence, do alignment (Viterbi), recover trace
     */

#ifdef ALTIVEC 
    
    /* By default we call an Altivec routine that doesn't save the full
     * trace. This also means the memory usage is just proportional to the
     * model (hmm->M), so we don't need to check if space is OK unless
     * we need the trace.   
     */   
    
    if(do_forward && do_null2)
    {
        /* Need the trace - first check space */
        if (P7ViterbiSpaceOK(sqinfo->len, hmm->M, mx))
        {
            /* Slower altivec version */
            sc = P7Viterbi(dsq, sqinfo->len, hmm, mx, &tr);
        }
        else
        {
            /* Low-memory C version */
            sc = P7SmallViterbi(dsq, sqinfo->len, hmm, mx, &tr);
        }
    }
    else
    {
        /* Fastest altivec version */
        sc = P7ViterbiNoTrace(dsq, sqinfo->len, hmm, mx);
        tr = NULL;
    }

#else /* not altivec */
    
    if (P7ViterbiSpaceOK(sqinfo->len, hmm->M, mx))
    {
        SQD_DPRINTF1(("   ... using normal P7Viterbi(); size ~%d MB\n", 
                      P7ViterbiSize(sqinfo->len, hmm->M)));
        sc = P7Viterbi(dsq, sqinfo->len, hmm, mx, &tr);
    }
    else
    {
        SQD_DPRINTF1(("   ... using P7SmallViterbi(); size ~%d MB\n",
                      P7ViterbiSize(sqinfo->len, hmm->M)));
        sc = P7SmallViterbi(dsq, sqinfo->len, hmm, mx, &tr);
    }

#endif
    
    /* Implement do_forward; we'll override the whole_sc with a P7Forward()
     * calculation.
     * HMMER is so trace- (alignment-) dependent that this gets a bit hacky.
     * Some important implications: 
     *   1) if --do_forward is selected, the domain (Viterbi) scores do not
     *      necessarily add up to the whole sequence (Forward) score.
     *   2) The implementation of null2 for a Forward score is undefined,
     *      since the null2 correction is trace-dependent. As a total hack, 
     *      we use a null2 correction derived from the whole trace
     *      (which was the behavior of HMMER 2.1.1 and earlier, anyway).
     *      This could put the sum of domain scores and whole seq score even
     *      further in disagreement. 
     * 
     * Note that you can't move the Forward calculation into 
     * PostprocessSignificantHit(). The Forward score will exceed the
     * Viterbi score, so you can't apply thresholds until you
     * know the Forward score. Also, since PostprocessSignificantHit()
     * is wrapped by a mutex in the threaded implementation,
     * you'd destroy all useful parallelism if PostprocessSignificantHit()
     * did anything compute intensive. 
     */
    if (do_forward) 
    {
        
        sc = P7Forward(dsq, sqinfo->len, hmm, NULL);
        if (do_null2) 
        {
            /* We need the trace - recalculate it if we didn't already do it */
#ifdef ALTIVEC
            if(tr == NULL)
            {
                if (P7ViterbiSpaceOK(sqinfo->len, hmm->M, mx))
                {
                    /* Slower altivec version */
                    P7Viterbi(dsq, sqinfo->len, hmm, mx, &tr);
                }
                else
                {
                    /* Low-memory C version */
                    P7SmallViterbi(dsq, sqinfo->len, hmm, mx, &tr);
                }
            }                
#endif            
            sc -= TraceScoreCorrection(hmm, tr, dsq);
        }
	}
    /* Store scores/pvalue for each HMM aligned to this sequence, overall
     */
    pvalue = PValue(hmm, sc);
    evalue = thresh->Z ? (double) thresh->Z * pvalue : (double) nhmm * pvalue;
    if (sc >= thresh->globT && evalue <= thresh->globE) {
        /* Recalculate trace if we used altivec */
#ifdef ALTIVEC
        if(tr == NULL)
        {
            if (P7ViterbiSpaceOK(sqinfo->len, hmm->M, mx))
            {
                /* Slower altivec version */
                P7Viterbi(dsq, sqinfo->len, hmm, mx, &tr);
            }
            else
            {
                /* Low-memory C version */
                P7SmallViterbi(dsq, sqinfo->len, hmm, mx, &tr);
            }
        }
#endif
        sc = PostprocessSignificantHit(ghit, dhit, 
                                       tr, hmm, dsq, sqinfo->len, 
                                       sqinfo->name, NULL, NULL, /* won't need acc or desc even if we have 'em */
                                       do_forward, sc,
                                       do_null2,
                                       thresh,
                                       TRUE); 
        /* TRUE -> hmmpfam mode */
    }
    if(tr != NULL)
        P7FreeTrace(tr);
    
    tr = NULL;
    
    FreePlan7(hmm);
    nhmm++;
  }

  FreePlan7Matrix(mx);
  free(dsq);
  *ret_nhmm = nhmm;
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
 *           thresh  - score/evalue threshold settings
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
	      struct threshold_s *thresh, int do_xnu, int do_forward, int do_null2,
	      struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm)
{
  struct plan7_s   *hmm;        /* HMM that was searched with */
  struct p7trace_s *tr;         /* a traceback structure */
  unsigned char    *dsq;        /* digitized sequence */ 
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
  int    arglen;

  /* Sanity checks.
   */
  if (hmmfp->ssi == NULL)
    Die("HMM file %s needs an SSI index to use PVM. See: hmmindex.", hmmfile);

  /* Prepare sequence.
   */
  dsq = DigitizeSequence(seq, sqinfo->len);
  if (do_xnu && Alphabet_type == hmmAMINO) XNU(dsq, sqinfo->len);

  /* Initialize PVM
   */
  master_tid = pvm_mytid();
#if DEBUGLEVEL >= 1  
  pvm_catchout(stderr);		/* catch output for debugging */
#endif
  SQD_DPRINTF1(("Spawning slaves...\n"));
  PVMSpawnSlaves("hmmpfam-pvm", &slave_tid, &nslaves);
  hmmlist   = MallocOrDie(sizeof(int) * nslaves);
  SQD_DPRINTF1(("Spawned a total of %d slaves...\n", nslaves));

  /* Initialize the slaves
   */
  SQD_DPRINTF1(("Broadcasting to %d slaves...\n", nslaves));
  pvm_initsend(PvmDataDefault);
  arglen = strlen(hmmfile);
  pvm_pkint(&arglen, 1, 1);
  pvm_pkstr(hmmfile);
  pvm_pkint(&(sqinfo->len), 1, 1);
  pvm_pkstr(seq);
  pvm_pkfloat(&(thresh->globT), 1, 1);
  pvm_pkdouble(&(thresh->globE), 1, 1);
  pvm_pkint(&(thresh->Z), 1, 1);
  pvm_pkint((int *)&(thresh->autocut), 1, 1);
  pvm_pkint(&do_xnu, 1, 1);
  pvm_pkint(&do_forward, 1, 1);
  pvm_pkint(&do_null2, 1, 1);
  pvm_pkint(&Alphabet_type, 1, 1);
  pvm_mcast(slave_tid, nslaves, HMMPVM_INIT);
  SQD_DPRINTF1(("Slaves should be ready...\n"));
				/* get their OK codes. */
  PVMConfirmSlaves(slave_tid, nslaves);
  SQD_DPRINTF1(("Slaves confirm that they're ok...\n"));

  /* Load the slaves.
   * For efficiency reasons, we don't want the master to
   * load HMMs from disk until she absolutely needs them.
   */
  for (nhmm = 0; nhmm < nslaves && nhmm < hmmfp->ssi->nprimary; nhmm++) {
    pvm_initsend(PvmDataDefault);
    pvm_pkint(&nhmm, 1, 1);	/* side effect: also tells him what number he is. */
    pvm_send(slave_tid[nhmm], HMMPVM_WORK);
    hmmlist[nhmm] = nhmm;
  }
  SQD_DPRINTF1(("%d slaves are loaded\n", nhmm));

  
  /* Receive/send loop
   */
  for (; nhmm < hmmfp->ssi->nprimary; nhmm++)
    {
				/* check slaves before blocking */
      PVMCheckSlaves(slave_tid, nslaves);
				/* receive output */
      SQD_DPRINTF1(("Waiting for a slave to give me output...\n"));
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&slaveidx, 1, 1);     /* # of slave who's sending us stuff */
      pvm_upkstr(slavename);           /* name of HMM that slave did */
      pvm_upkfloat(&sc, 1, 1);         /* score   */
      pvm_upkdouble(&pvalue, 1, 1);    /* P-value */
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
	  if (! SetAutocuts(thresh, hmm))
	    Die("HMM %s did not contain your GA, NC, or TC cutoffs", hmm->name);
	
	  sc = PostprocessSignificantHit(ghit, dhit, 
					 tr, hmm, dsq, sqinfo->len, 
					 sqinfo->name, 
					 sqinfo->flags & SQINFO_ACC ? sqinfo->acc : NULL,
					 sqinfo->flags & SQINFO_DESC ? sqinfo->desc : NULL,
					 do_forward, sc,
					 do_null2,
					 thresh,
					 TRUE); /* TRUE -> hmmpfam mode */

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
	  if (! SetAutocuts(thresh, hmm))
	    Die("HMM %s did not contain your GA, NC, or TC cutoffs", hmm->name);
	  
	  sc = PostprocessSignificantHit(ghit, dhit, 
					 tr, hmm, dsq, sqinfo->len, 
					 sqinfo->name, NULL, NULL, /* won't need acc or desc even if we have 'em */
					 do_forward, sc,
					 do_null2,
					 thresh,
					 TRUE); /* TRUE -> hmmpfam mode */

	  FreePlan7(hmm);
	  P7FreeTrace(tr);
	}
    }

  /* Shut down all the slaves. (xref bug #15/STL5 p.66).
   */
  for (slave = 0; slave < nslaves; slave++)
    {
      pvm_initsend(PvmDataDefault);
      msg = -1;			/* the special "shut down now" flag */
      pvm_pkint(&msg, 1, 1);
      pvm_send(slave_tid[slave], HMMPVM_WORK);
    }

  /* Cleanup; quit the PVM; and return
   */
  free(slave_tid);
  free(hmmlist);
  free(dsq);
  pvm_exit();
  *ret_nhmm = nhmm;
  return;
}
#else /* HMMER_PVM off; insert dummy stub function: */
static void
main_loop_pvm(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo, 
	      struct threshold_s *thresh, int do_xnu, int do_forward, int do_null2,
	      struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm)
{
  Die("no PVM support");
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
/* Function: main_loop_threaded()
 * Date:     SRE, Wed Feb 27 11:10:13 2002 [St. Louis]
 *
 * Purpose:  Search a sequence against an HMM database;
 *           main loop for the threaded version. 
 *           
 *           On return, ghit and dhit contain info for all hits
 *           that satisfy the set thresholds. If an evalue
 *           cutoff is used at all, the lists will be overestimated --
 *           because the evalue will be underestimated until
 *           we know the final Z. (Thus the main program must recheck
 *           thresholds before printing any results.) If only
 *           score cutoffs are used, then the lists are correct,
 *           and may be printed exactly as they come (after 
 *           appropriate sorting, anyway). This is especially
 *           important for dynamic thresholding using Pfam
 *           score cutoffs -- the main caller cannot afford to
 *           rescan the HMM file just to get the GA/TC/NC cutoffs
 *           back out for each HMM, and neither do I want to
 *           burn the space to store them as I make a pass thru
 *           Pfam.
 *
 * Args:     hmmfile - name of HMM file
 *           hmmfp   - open HMM file (and at start of file)
 *           dsq     - digitized sequence
 *           sqinfo  - ptr to SQINFO optional info for dsq
 *           thresh  - score/evalue threshold information
 *           do_xnu     - TRUE to apply XNU filter to sequence
 *           do_forward - TRUE to use Forward() scores
 *           do_null2   - TRUE to adjust scores w/ ad hoc null2 model
 *           num_threads- number of threads, >= 1
 *           ghit    - global hits list
 *           dhit    - domain hits list
 *           ret_nhmm    - number of HMMs searched.
 *
 * Returns:  (void)
 */
static void
main_loop_threaded(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo, 
		   struct threshold_s *thresh, int do_xnu, int do_forward, int do_null2,
		   int num_threads,
		   struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm) 
{
  unsigned char     *dsq;       /* digitized sequence      */
  int                nhmm;	/* number of HMMs searched */
  struct workpool_s *wpool;     /* pool of worker threads  */

  /* Prepare sequence.
   */
  dsq = DigitizeSequence(seq, sqinfo->len);
  if (do_xnu && Alphabet_type == hmmAMINO) XNU(dsq, sqinfo->len);

  wpool = workpool_start(hmmfile, hmmfp, dsq, sqinfo->name, sqinfo->len, 
			 do_forward, do_null2, thresh,
			 ghit, dhit, num_threads);
  workpool_stop(wpool);
  nhmm = wpool->nhmm;
  workpool_free(wpool);

  free(dsq);
  *ret_nhmm = nhmm;
  return;
}
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
 *           threshold  - evalue/score threshold settings
 *           ghit       - per-seq hit list
 *           dhit       - per-domain hit list             
 *           num_threads- number of worker threads to run.
 *
 * Returns:  ptr to struct workpool_s.
 *           Caller must wait for threads to finish with workpool_stop(),
 *           then free the structure with workpool_free().
 */
static struct workpool_s *
workpool_start(char *hmmfile, HMMFILE *hmmfp, unsigned char *dsq, char *seqname, int L,
	       int do_forward, int do_null2, struct threshold_s *thresh,
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
  wpool->thresh     = thresh;

  wpool->hmmfp      = hmmfp;
  wpool->nhmm       = 0;
  if ((rtn = pthread_mutex_init(&(wpool->input_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));

  wpool->ghit       = ghit;
  wpool->dhit       = dhit;
  if ((rtn = pthread_mutex_init(&(wpool->output_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));

  wpool->num_threads= num_threads;

  /* Create slave threads. See comments in hmmcalibrate.c at 
   * this step regarding concurrency and system scope.
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
  double evalue;		/* E-value of a score              */
  struct threshold_s thresh;	/* a local copy of thresholds      */
  struct dpmatrix_s *mx;        /* growable DP matrix              */

  tr = NULL;
  
  wpool = (struct workpool_s *) ptr;
  /* Because we might dynamically change the thresholds using
   * Pfam GA/NC/TC cutoffs, we make a local copy of the threshold
   * structure in this thread.
   */ 
  thresh.globT   = wpool->thresh->globT;
  thresh.globE   = wpool->thresh->globE;
  thresh.domT    = wpool->thresh->domT;
  thresh.domE    = wpool->thresh->domE;
  thresh.autocut = wpool->thresh->autocut;
  thresh.Z       = wpool->thresh->Z;
  
  /* 
   * We'll create for at least N=300xM=300, and thus consume at least 1 MB,
   * regardless of RAMLIMIT -- this helps us avoid reallocating some weird
   * asymmetric matrices.
   * 
   * We're growable in both M and N, because inside of P7SmallViterbi,
   * we're going to be calling P7Viterbi on subpieces that vary in size,
   * and for different models.
   */
  mx = CreatePlan7Matrix(300, 300, 25, 25);

  for (;;) {

    /* 1. acquire lock on HMM input, and get
     *    the next HMM to work on.
     */
				/* acquire a lock */
    if ((rtn = pthread_mutex_lock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    
    if (! HMMFileRead(wpool->hmmfp, &hmm)) 
      {	/* we're done. release lock, exit thread */
	if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
	  Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
	FreePlan7Matrix(mx);
	pthread_exit(NULL);
      }
    wpool->nhmm++;
    SQD_DPRINTF1(("a thread is working on %s\n", hmm->name));
				/* release the lock */
    if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

    if (hmm == NULL)
      Die("HMM file %s may be corrupt or in incorrect format; parse failed", wpool->hmmfile);

    if (!SetAutocuts(&thresh, hmm)) 
      Die("HMM %s did not have the right GA, NC, or TC cutoffs", hmm->name);

    /* 2. We have an HMM in score form.
     *    Score the sequence.
     */
    
#ifdef ALTIVEC 
    
    /* By default we call an Altivec routine that doesn't save the full
     * trace. This also means the memory usage is just proportional to the
     * model (hmm->M), so we don't need to check if space is OK unless
     * we need the trace.   
     */   
    
    if(wpool->do_forward && wpool->do_null2)
    {
        /* Need the trace - first check space */
        if (P7ViterbiSpaceOK(wpool->L, hmm->M, mx))
        {
            /* Slower altivec version */
            sc = P7Viterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
        }
        else
        {
            /* Low-memory C version */
            sc = P7SmallViterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
        }
    }
    else
    {
        /* Fastest altivec version */
        sc = P7ViterbiNoTrace(wpool->dsq, wpool->L, hmm, mx);
        tr = NULL;
    }

#else /* not altivec */
    
    if (P7ViterbiSpaceOK(wpool->L, hmm->M, mx))
    {
        sc = P7Viterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
    }
    else
    {
        sc = P7SmallViterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
    }

#endif
    
    /* The Forward score override (see comments in serial vers)
     */
    if (wpool->do_forward) {
      sc  = P7Forward(wpool->dsq, wpool->L, hmm, NULL);
      if (wpool->do_null2)  
      {
#ifdef ALTIVEC
          if(tr == NULL)
          {
              if (P7ViterbiSpaceOK(wpool->L, hmm->M, mx))
              {
                  /* Slower altivec version */
                  sc = P7Viterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
              }
              else
              {
                  /* Low-memory C version */
                  sc = P7SmallViterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
              }            
          }
#endif
          sc -= TraceScoreCorrection(hmm, tr, wpool->dsq);
      }
    }
    
    /* 3. Save the output in tophits structures, after acquiring a lock
     */
    if ((rtn = pthread_mutex_lock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    SQD_DPRINTF1(("model %s scores %f\n", hmm->name, sc));
    
    pvalue = PValue(hmm, sc);
    evalue = thresh.Z ? (double) thresh.Z * pvalue : (double) wpool->nhmm * pvalue;
    if (sc >= thresh.globT && evalue <= thresh.globE) 
    { 
#ifdef ALTIVEC
        if(tr == NULL)
        {
            if (P7ViterbiSpaceOK(wpool->L, hmm->M, mx))
            {
                /* Slower altivec version */
                sc = P7Viterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
            }
            else
            {
                /* Low-memory C version */
                sc = P7SmallViterbi(wpool->dsq, wpool->L, hmm, mx, &tr);
            }            
        }
#endif        
        sc = PostprocessSignificantHit(wpool->ghit, wpool->dhit, 
                                       tr, hmm, wpool->dsq, wpool->L, 
                                       wpool->seqname, 
                                       NULL, NULL, /* won't need seq's acc or desc */
                                       wpool->do_forward, sc,
                                       wpool->do_null2,
                                       &thresh,
                                       TRUE); 
        /* TRUE -> hmmpfam mode */
    }
    if ((rtn = pthread_mutex_unlock(&(wpool->output_lock))) != 0)
        Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
    
    if(tr != NULL)
        P7FreeTrace(tr);
    
    tr = NULL;
    
    FreePlan7(hmm);

  } /* end 'infinite' loop over HMMs in this thread */
}
#else /* HMMER_THREADS off; put in a dummy stub: */
static void
main_loop_threaded(char *hmmfile, HMMFILE *hmmfp, char *seq, SQINFO *sqinfo, 
		   struct threshold_s *thresh, int do_xnu, int do_forward, int do_null2,
		   int num_threads,
		   struct tophit_s *ghit, struct tophit_s *dhit, int *ret_nhmm) 
{
  Die("no threads support");
}
#endif /* HMMER_THREADS */

/************************************************************
 * @LICENSE@
 ************************************************************/

