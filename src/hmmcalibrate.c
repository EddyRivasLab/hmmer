/* hmmcalibrate.c
 * Score an HMM against random sequence data to set the statistical
 * parameters for E-value determination.
 * 
 * SRE, Fri Oct 31 09:25:21 1997 [St. Louis]
 * SVN $Id$
 */

#include "config.h"		/* compile-time configuration constants */
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#ifdef HMMER_THREADS
#include <pthread.h>
#endif
#ifdef HMMER_PVM
#include <pvm3.h>
#endif

#include "squid.h"		/* general sequence analysis library    */
#include "stopwatch.h"		/* process timings                      */

#include "plan7.h"		/* plan 7 profile HMM structure         */
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */


static char banner[] = "hmmcalibrate -- calibrate HMM search statistics";

static char usage[] = "\
Usage: hmmcalibrate [-options] <hmmfile>\n\
Available options are:\n\
  -h             : print short usage and version info, then exit\n\
";

static char experts[] = "\
  --cpu <n>      : run <n> threads in parallel (if threaded)\n\
  --dry          : show output, but don't save params in the HMM(s)\n\
  --fixed <n>    : fix random sequence length at <n>\n\
  --histfile <f> : save fitted histogram(s) to file <f>\n\
  --mult <x>     : set length multiplier to <x> [3.0]\n\
  --num <n>      : set number of sampled seqs to <n> [5000]\n\
  --phist <f>    : save predicted histogram(s) to file <f>\n\
  --pvm          : run on a Parallel Virtual Machine (PVM)\n\
  --sctbl <f>    : save score/Evalue table to file <f>\n\
  --seed <n>     : set random seed to <n> [time()]\n\
";

static struct opt_s OPTIONS[] = {
   { "-h",         TRUE,  sqdARG_NONE  },
   { "--cpu",      FALSE, sqdARG_INT },
   { "--dry",      FALSE, sqdARG_NONE },
   { "--fixed",    FALSE, sqdARG_INT   },
   { "--histfile", FALSE, sqdARG_STRING },
   { "--mult",     FALSE, sqdARG_FLOAT },
   { "--num",      FALSE, sqdARG_INT   },
   { "--phist",    FALSE, sqdARG_STRING },
   { "--pvm",      FALSE, sqdARG_NONE  },
   { "--sctbl",    FALSE, sqdARG_STRING },
   { "--seed",     FALSE, sqdARG_INT}, 
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


static void   imean_sd(int *v, int N, double *ret_mean, double *ret_sd);
static void   dmean_sd(double *v, int N, double *ret_mean, double *ret_sd);
static double dmode(double *v, int N, double window);
static void sort_scores(double *sc, int N);
static void save_fitted_histogram(FILE *hfp, char *name, double *sc, int N, double mu, double lambda);
static void save_predicted_histogram(FILE *hfp, struct plan7_s *hmm, double *sc, int N, int L);
static void save_score_list(FILE *sfp, double *sc, int N, int L, 
			    struct plan7_s *hmm, double mu, double lambda);


static void main_loop_serial(struct plan7_s *hmm, int seed, int N, int L,
			     double *sc, int *alen);

#ifdef HMMER_THREADS
/* A structure of this type is shared by worker threads in the POSIX
 * threads parallel version.
 */
struct workpool_s {
  /* Static configuration:
   */
  struct plan7_s  *hmm;	     /* ptr to single HMM to search with    */
  int    N;		     /* number of random seqs to do         */
  int    L;		     /* length of random sequences          */
  float *randomseq;          /* 0..Alphabet_size-1 i.i.d. probs     */

  /* Shared (mutex-protected) input:
   */
  int    nseq;		     /* current number of seqs searched     */

  /* Shared (mutex-protected) output:
   */
  double        *sc;		/* unsorted score array */
  int           *alen;		/* opt ali length array  */
  Stopwatch_t    watch;		/* Timings accumulated for threads */

  /* Thread pool information:
   */
  pthread_t      *thread;       /* our pool of threads */
  int             num_threads;  /* number of threads   */
  pthread_mutex_t input_lock;	/* a mutex protecting input fields */
  pthread_mutex_t output_lock;  /* a mutex protecting output fields */
};
static void main_loop_threaded(struct plan7_s *hmm, int seed, int N, int L,
			       int nthreads,
			       double *sc, int *alen, Stopwatch_t *twatch);
static struct workpool_s *workpool_start(struct plan7_s *hmm, int N, int L,
					 float *randomseq, 
					 double *sc, int *alen, int num_threads);
static void  workpool_stop(struct workpool_s *wpool);
static void  workpool_free(struct workpool_s *wpool);
static void *worker_thread(void *ptr);
#endif /* HMMER_THREADS */

#ifdef HMMER_PVM
static void main_loop_pvm(struct plan7_s *hmm, int seed, int N, int L,
			  int lumpsize,
			  double *sc, int *alen, 
			  Stopwatch_t *extrawatch, int *ret_nslaves);
#endif /* HMMER_PVM */


int
main(int argc, char **argv)
{
  char    *hmmfile;             /* HMM file to open                */
  char    *tmpfile;             /* temporary calibrated HMM file   */
  HMMFILE *hmmfp;               /* opened hmm file pointer         */
  FILE    *outfp;               /* for writing HMM(s) into tmpfile */
  char    *mode;                /* write mode, "w" or "wb"         */
  struct plan7_s *hmm;          /* the hidden Markov model         */
  int     idx;			/* counter over sequences          */

  /* Controls over the random sequence simulation: */
  int     N;			/* number of random seqs to sample */
  float   mult;			/* default: L= mult * hmm->M       */
  int     fixedL;		/* if >0, L is fixed to this       */
  int     seed;			/* random number seed              */
  double *sc;			/* array of <nsample> scores       */
  int    *alen;			/* array of <nsample> ali lengths  */

  /* For each HMM, we obtain these results: */
  int    *L;		        /* random sequence lengths per HMM */
  double *mu;			/* array of EVD mu's for HMMs      */
  double *lambda;		/* array of EVD lambda's for HMMs  */
  double *kappa;		/* array of edge correction kappa's*/
  double *sigma;		/* array of edge correction sigma's*/
  int     nhmm;			/* number of HMMs calibrated       */
  int     halloc;		/* number of HMMs allocated for    */
  double  mean, median, sd;	/* mean, median, sd for scores     */
  double  scmode;

  /* Control over optional output files:
   */
  char   *histfile;             /* histogram save file             */
  FILE   *hfp;                  /* open file pointer for histfile  */
  char   *phistfile;		/* predicted histogram save file   */
  FILE   *pfp;			/* open file ptr for phistfile     */
  char   *sctblfile;		/* score/Evalue table              */
  FILE   *sfp;			/* open file ptr for sctblfile     */

  /* Control over PVM parallelization:
   */
  int     do_pvm;		/* TRUE to use PVM                 */
  int     pvm_lumpsize;		/* # of seqs to do per PVM slave exchange */
  int     pvm_nslaves;		/* number of slaves used in the PVM */

  /* Parsing the command line:
   */
  char *optname;		/* name of option found by Getopt() */
  char *optarg;			/* argument found by Getopt()       */
  int   optind;		        /* index in argv[]                  */
  int   do_dryrun;		/* TRUE to not save params          */

  /* Control over POSIX thread parallelization
   */
  int   num_threads;            /* number of worker threads */   

  /* Misc timing and OS control
   */
  Stopwatch_t stopwatch;	/* main stopwatch for process       */
  Stopwatch_t extrawatch;	/* stopwatch for threads/PVM slaves */
  sigset_t blocksigs;		/* list of signals to protect from */

  StopwatchStart(&stopwatch);
  StopwatchZero(&extrawatch);

  /***********************************************
   * Parse the command line
   ***********************************************/

  do_dryrun    = FALSE;
  N            = 5000;
  mult         = 3.0;		/* default: make seqs 3.0 * HMM length */
  fixedL       = 0;
  seed         = (int) time ((time_t *) NULL);
  histfile     = NULL;
  phistfile    = NULL;
  sctblfile    = NULL;
  do_pvm       = FALSE;
  pvm_lumpsize = 20;		/* 20 seqs/PVM exchange: sets granularity */
  pvm_nslaves  = 0;	/* solely to silence compiler warnings: PVM main sets this */
  num_threads  = ThreadNumber(); /* 0 if unthreaded; else >=1 default */

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "--cpu")      == 0) num_threads  = atoi(optarg);
      else if (strcmp(optname, "--dry")      == 0) do_dryrun    = TRUE;
      else if (strcmp(optname, "--fixed")    == 0) fixedL       = atoi(optarg);
      else if (strcmp(optname, "--histfile") == 0) histfile     = optarg;
      else if (strcmp(optname, "--mult")     == 0) mult         = atof(optarg); 
      else if (strcmp(optname, "--num")      == 0) N            = atoi(optarg); 
      else if (strcmp(optname, "--phist")    == 0) phistfile    = optarg;
      else if (strcmp(optname, "--pvm")      == 0) do_pvm       = TRUE;
      else if (strcmp(optname, "--sctbl")    == 0) sctblfile    = optarg;
      else if (strcmp(optname, "--seed")     == 0) seed         = atoi(optarg);
      else if (strcmp(optname, "-h") == 0)
	{
	  HMMERBanner(stdout, banner);
	  puts(usage);
	  puts(experts);
	  exit(0);
	}
    }

  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  hmmfile = argv[optind++];

#ifndef HMMER_PVM
  if (do_pvm) Die("PVM support is not compiled into HMMER; --pvm doesn't work.");
#endif
#ifndef HMMER_THREADS
  if (num_threads) Die("Posix threads support is not compiled into HMMER; --cpu doesn't have any effect");
#endif

  /***********************************************
   * Open our i/o file pointers, make sure all is well
   ***********************************************/

				/* HMM file */
  if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
    Die("failed to open HMM file %s for reading.", hmmfile);

				/* histogram files, fitted and predicted */
  hfp = NULL;
  if (histfile != NULL) {
    if ((hfp = fopen(histfile, "w")) == NULL)
      Die("Failed to open histogram save file %s for writing\n", histfile);
  }
  pfp = NULL;
  if (phistfile != NULL) {
    if ((pfp = fopen(phistfile, "w")) == NULL)
      Die("Failed to open predicted histogram file %s for writing\n", phistfile);
  }
  sfp = NULL;
  if (sctblfile != NULL) {
    if ((sfp = fopen(sctblfile, "w")) == NULL)
      Die("Failed to open score table file 5s for writing\n", sctblfile);
  }

  /* Generate calibrated HMM(s) in a tmp file in the current
   * directory. When we're finished, we delete the original
   * HMM file and rename() this one. That way, the worst
   * effect of a catastrophic failure should be that we
   * leave a tmp file lying around, but the original HMM
   * file remains uncorrupted. tmpnam() doesn't work portably here,
   * because it'll put the file in /tmp and we won't
   * necessarily be able to rename() it from there.
   */
  tmpfile = MallocOrDie(strlen(hmmfile) + 5);
  strcpy(tmpfile, hmmfile);
  strcat(tmpfile, ".xxx");	/* could be more inventive here... */
  if (FileExists(tmpfile))
    Die("temporary file %s already exists; please delete it first", tmpfile);
  if (hmmfp->is_binary) mode = "wb";
  else                  mode = "w"; 

  /*********************************************** 
   * Show the banner
   ***********************************************/

  HMMERBanner(stdout, banner);
  printf("HMM file:                 %s\n", hmmfile);
  if (fixedL) 
    printf("Length fixed to:          %d\n", fixedL);
  else {
    printf("Length multiplier:        %.0f\n", mult);
  }
  printf("Number of samples:        %d\n", N);
  printf("random seed:              %d\n", seed);
  if (histfile != NULL) 
    printf("histogram(s) saved to:    %s\n", histfile);
  if (phistfile != NULL) 
    printf("predicted histogram(s) to: %s\n", phistfile);

  if (do_pvm)
    printf("PVM:                      ACTIVE\n");
  else if (num_threads > 0)
    printf("POSIX threads:            %d\n", num_threads);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");


  /* Initial allocations for results per HMM;
   * we'll resize these arrays dynamically as we read more HMMs.
   */
  halloc = 128;
  L      = MallocOrDie(sizeof(int)    * halloc); 
  mu     = MallocOrDie(sizeof(double) * halloc); 
  lambda = MallocOrDie(sizeof(double) * halloc); 
  kappa  = MallocOrDie(sizeof(double) * halloc); 
  sigma  = MallocOrDie(sizeof(double) * halloc); 
  nhmm   = 0;

  /* Allocation for results per simulation, for one HMM.
   */
  sc     = MallocOrDie(sizeof(double) * N);
  alen   = MallocOrDie(sizeof(int)    * N);


  /***********************************************
   * Read the HMMs one at a time, and send them off
   * to one of the main loops to collect sc, alen arrays.
   ***********************************************/

  while (HMMFileRead(hmmfp, &hmm))
    {
      if (hmm == NULL)
	Die("HMM file may be corrupt or in incorrect format; parse failed");

      /* Configure for SW calibration.
       */
      P7Config(hmm, P7_SW_MODE); 

      /* Determine what length the random seqs should be.
       */
      if (fixedL) L[nhmm] = fixedL;
      else        L[nhmm] = (int) (mult * (float) hmm->M);

      /* Use one of the engines to collect sc, alen arrays
       * for N random seqs of length L. This is the compute
       * intensive part.
       */
      if (! do_pvm && num_threads == 0)
	main_loop_serial(hmm, seed, N, L[nhmm], sc, alen);
#ifdef HMMER_PVM
      else if (do_pvm) {
	main_loop_pvm(hmm, seed, N, L[nhmm], pvm_lumpsize, 
		      sc, alen, &extrawatch, &pvm_nslaves);
      }
#endif 
#ifdef HMMER_THREADS
      else if (num_threads > 0)
	main_loop_threaded(hmm, seed, N, L[nhmm], num_threads, 
			   sc, alen, &extrawatch);
#endif
      else Die("That can't happen. I didn't do anything.");

      /* Fit an EVD to the scores to get mu, lambda params.
       */
      EVDMaxLikelyFit(sc, NULL, N, &(mu[nhmm]), &(lambda[nhmm]));

      /* Calculate mean ali length (kappa) and std. dev of ali lengths (sigma).
       */
      imean_sd(alen, N, &(kappa[nhmm]), &(sigma[nhmm]));

      sort_scores(sc, N);
      dmean_sd(sc, N, &mean, &sd);
      median = sc[N/2];
      scmode = dmode(sc, N, 2.0); /* mode in 2.0 bit windows */

      /* Output
       */
      printf("HMM    : %s\n",   hmm->name);
      printf("L      : %d\n",   L[nhmm]);
      printf("mean   : %12f\n", mean);
      printf("stddev : %12f\n", sd);
      printf("median : %12f\n", median);
      printf("scmode : %12f\n", scmode);
      printf("mu     : %12f\n", mu[nhmm]);
      printf("lambda : %12f\n", lambda[nhmm]);
      printf("kappa  : %12f\n", kappa[nhmm]);
      printf("sigma  : %12f\n", sigma[nhmm]);
      printf("//\n");

      if (hfp != NULL) 	save_fitted_histogram(hfp, hmm->name, sc, N, mu[nhmm], lambda[nhmm]);
      if (pfp != NULL)	save_predicted_histogram(pfp, hmm, sc, N, L[nhmm]); 
      if (sfp != NULL)	save_score_list(sfp, sc, N, L[nhmm], hmm, mu[nhmm], lambda[nhmm]);

      /* Reallocation, if needed.
       */
      nhmm++;
      if (nhmm == halloc) {
	halloc *= 2;		/* realloc by doubling */
	L      = ReallocOrDie(L,      sizeof(int)    * halloc);
	mu     = ReallocOrDie(mu,     sizeof(double) * halloc);
	lambda = ReallocOrDie(lambda, sizeof(double) * halloc);
	kappa  = ReallocOrDie(kappa,  sizeof(double) * halloc);
	sigma  = ReallocOrDie(sigma,  sizeof(double) * halloc);
      }      

      FreePlan7(hmm);
    }
  SQD_DPRINTF1(("Main body believes it has calibrations for %d HMMs\n", nhmm));

  /*****************************************************************
   * Rewind the HMM file for a second pass.
   * Write a temporary HMM file with new mu, lambda values in it
   *****************************************************************/
  if (do_dryrun) 
    {
      HMMFileClose(hmmfp);
    }
  else 
    {
      HMMFileRewind(hmmfp);
      if (FileExists(tmpfile))
	Die("Ouch. Temporary file %s appeared during the run.", tmpfile);
      if ((outfp = fopen(tmpfile, mode)) == NULL)
	Die("Ouch. Temporary file %s couldn't be opened for writing.", tmpfile); 
      
      for (idx = 0; idx < nhmm; idx++)
	{
	  /* Sanity checks 
	   */
	  if (!HMMFileRead(hmmfp, &hmm))
	    Die("Ran out of HMMs too early in pass 2");
	  if (hmm == NULL) 
	    Die("HMM file %s was corrupted? Parse failed in pass 2", hmmfile);

	  /* Put results in HMM
	   */
	  hmm->mu     = mu[idx];
	  hmm->lambda = lambda[idx];
	  hmm->kappa  = kappa[idx];
	  hmm->sigma  = sigma[idx];
	  hmm->Lbase  = L[idx];
	  hmm->flags |= PLAN7_STATS;
	  Plan7ComlogAppend(hmm, argc, argv);

	  /* Save HMM to tmpfile
	   */
	  if (hmmfp->is_binary) WriteBinHMM(outfp, hmm);
	  else                  WriteAscHMM(outfp, hmm); 

	  FreePlan7(hmm);
	}
  
      /*****************************************************************
       * Now, carefully remove original file and replace it
       * with the tmpfile. Note the protection from signals;
       * we wouldn't want a user to ctrl-C just as we've deleted
       * their HMM file but before the new one is moved.
       *****************************************************************/

      HMMFileClose(hmmfp);
      if (fclose(outfp)   != 0)                            PANIC;
      if (sigemptyset(&blocksigs) != 0)                    PANIC;
      if (sigaddset(&blocksigs, SIGINT) != 0)              PANIC;
      if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   PANIC;
      if (remove(hmmfile) != 0)                            PANIC;
      if (rename(tmpfile, hmmfile) != 0)                   PANIC;
      if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) PANIC;
      free(tmpfile);
    } /* end of ! do_dryrun section */

  /***********************************************
   * Exit
   ***********************************************/

  StopwatchStop(&stopwatch);
#ifdef HMMER_PVM
  if (do_pvm > 0) {
    printf("PVM processors used: %d\n", pvm_nslaves);
    StopwatchInclude(&stopwatch, &extrawatch);
  }
#endif
#ifdef PTHREAD_TIMES_HACK
  else if (num_threads > 0) StopwatchInclude(&stopwatch, &extrawatch);
#endif

  /*  StopwatchDisplay(stdout, "CPU Time: ", &stopwatch); */

  free(sc);
  free(alen);
  free(L);
  free(mu);
  free(lambda);
  free(kappa);
  free(sigma);
  if (hfp != NULL) fclose(hfp);
  if (pfp != NULL) fclose(pfp);
  SqdClean();
  return 0;
}

/* Function: main_loop_serial()
 * Date:     SRE, Tue Aug 18 16:18:28 1998 [St. Louis]
 *
 * Purpose:  Given an HMM and parameters for synthesizing random
 *           sequences; return a histogram of scores.
 *           (Serial version)  
 *
 * Args:     hmm      - an HMM to calibrate.
 *           seed     - random number seed
 *           N        - number of seqs to synthesize
 *           L        - mean length of random sequence
 *           sc       - RETURN: unsorted array of <N> scores
 *           alen     - RETURN: array of <N> optimal ali lengths
 *
 * Returns:  (void)
 */
static void
main_loop_serial(struct plan7_s *hmm, int seed, int N, int L,
		 double *sc, int *alen)
{
  struct p7trace_s   *tr;
  struct dpmatrix_s  *mx;
  float  randomseq[MAXABET];
  float  p1;
  char           *seq;
  unsigned char  *dsq;
  float  score;
  int    idx;
  
  /* Initialize.
   * We assume we've already set the alphabet (safe, because
   * HMM input sets the alphabet).
   */
  sre_srandom(seed);
  if (! (hmm->flags & PLAN7_HASBITS)) Die("oops: that model isn't configured");
  P7DefaultNullModel(randomseq, &p1);
  mx = CreatePlan7Matrix(L, hmm->M, 0, 0);

  for (idx = 0; idx < N; idx++)
  {
      seq = RandomSequence(Alphabet, randomseq, Alphabet_size, L);
      dsq = DigitizeSequence(seq, L);
      
      if (P7ViterbiSpaceOK(L, hmm->M, mx))
          score = P7Viterbi(dsq, L, hmm, mx, &tr);
      else
          score = P7SmallViterbi(dsq, L, hmm, mx, &tr);

      sc[idx] = (double) score;
      Trace_GetAlignmentBounds(tr, 1, NULL, NULL, NULL, NULL, &(alen[idx]));

      P7FreeTrace(tr);
      free(dsq); 
      free(seq);
    }

  FreePlan7Matrix(mx);
  return;
}


#ifdef HMMER_THREADS
/* Function: main_loop_threaded()
 * Date:     SRE, Wed Dec  1 12:43:09 1999 [St. Louis]
 *
 * Purpose:  Given an HMM and parameters for synthesizing random
 *           sequences; return a histogram of scores.
 *           (Threaded version.)  
 *
 * Args:     hmm      - an HMM to calibrate.
 *           seed     - random number seed
 *           N        - number of seqs to synthesize
 *           L        - sequence length
 *           nthreads - number of threads to start
 *           sc       - RETURN: unsorted array of <N> scores
 *           alen     - RETURN: array of <N> optimal ali lengths
 *           twatch   - RETURN: accumulation of thread times
 *
 * Returns:  (void)
 */
static void
main_loop_threaded(struct plan7_s *hmm, int seed, int N, int L,
		   int nthreads,
		   double *sc, int *alen, Stopwatch_t *twatch)
{
  float  randomseq[MAXABET];
  float  p1;
  struct workpool_s *wpool;     /* pool of worker threads  */
  
  /* Initialize.
   * We assume we've already set the alphabet (safe, because
   * HMM input sets the alphabet).
   */
  sre_srandom(seed);
  P7Logoddsify(hmm, TRUE);
  P7DefaultNullModel(randomseq, &p1);

  wpool = workpool_start(hmm, N, L, randomseq, sc, alen, nthreads);
  workpool_stop(wpool);

  StopwatchInclude(twatch, &(wpool->watch));

  workpool_free(wpool);
  return;
}

/*****************************************************************
 * POSIX threads implementation.
 * API:
 *      workpool_start()   (makes a workpool_s structure. Starts calculations.)
 *      workpool_stop()    (waits for threads to finish.)
 *      [process histogram]
 *      workpool_free()    (destroys the structure)
 *      
 * Threads:
 *      worker_thread()    (the actual parallelized worker thread).
 *****************************************************************/

/* Function: workpool_start()
 * Date:     SRE, Thu Jul 16 11:09:05 1998 [St. Louis]
 *
 * Purpose:  Initialize a workpool_s structure, and return it.
 *
 * Args:     hmm      - the HMM to calibrate
 *           fixedlen - 0, or a fixed length for seqs (bypass of Gaussian)
 *           lenmean  - mean sequence length 
 *           lensd    - std. dev. for sequence length
 *           randomseq- i.i.d. frequencies for residues, 0..Alphabet_size-1
 *           nsample  - how many seqs to calibrate on
 *           hist     - histogram structure for storing results
 *           sc       - unsorted array of score results 0..nsample-1
 *           alen     - RETURN: array of <nsample> optimal ali lengths
 *           num_threads - how many processors to run on
 *
 * Returns:  ptr to struct workpool_s.
 *           Caller must wait for threads to finish with workpool_stop(),
 *           then free the structure with workpool_free().
 */
static struct workpool_s *
workpool_start(struct plan7_s *hmm, int N, int L, float *randomseq, 
	       double *sc, int *alen, int num_threads)
{
  struct workpool_s *wpool;
  pthread_attr_t    attr;
  int i;
  int rtn;

  wpool         = MallocOrDie(sizeof(struct workpool_s));
  wpool->thread = MallocOrDie(num_threads * sizeof(pthread_t));
  wpool->hmm        = hmm;
  wpool->N          = N;
  wpool->L          = L;
  wpool->randomseq  = randomseq;
  
  wpool->nseq       = 0;
  wpool->sc         = sc;
  wpool->alen       = alen;
  wpool->num_threads= num_threads;

  StopwatchZero(&(wpool->watch));
  
  if ((rtn = pthread_mutex_init(&(wpool->input_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));
  if ((rtn = pthread_mutex_init(&(wpool->output_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));

  /* Create slave threads.
   * Note the crazy machinations we have to go through to achieve concurrency.
   * You'd think that POSIX threads were portable... ha.
   * On IRIX 6.5, system scope threads are only available to root, or if
   *   /etc/capability has been configured specially, so to avoid strange
   *   permissions errors we can't set PTHREAD_SCOPE_SYSTEM for IRIX.
   * On IRIX pre-6.5, we can't get good concurrency, period. As of 6.5,
   *   SGI provides the nonportable pthread_setconcurrency() call.
   * On FreeBSD (3.0 snapshots), the pthread_attr_setscope() call isn't
   *   even provided, apparently on grounds of "if it doesn't do anything,
   *   why provide it?" Hello? POSIX compliance, perhaps?
   * On Sun Solaris, we need to set system scope to achieve concurrency.
   * Linux and DEC Digital UNIX seem to work fine in either process scope
   *   or system scope, without a pthread_setconcurrency call.
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
 * Date:     SRE, Thu Jul 16 10:41:02 1998 [St. Louis]
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
  struct plan7_s    *hmm;
  struct p7trace_s  *tr;
  struct dpmatrix_s *mx;
  struct workpool_s *wpool;
  char              *seq;
  unsigned char     *dsq;
  float       sc;
  int         rtn;
  Stopwatch_t thread_watch;

  StopwatchStart(&thread_watch);
  wpool = (struct workpool_s *) ptr;
  hmm   = wpool->hmm;
  mx    = CreatePlan7Matrix(wpool->L, hmm->M, 0, 0);
  for (;;)
    {
      /* 1. Synthesize a random sequence. 
       *    The input sequence number is a shared resource,
       *    and sre_random() isn't thread-safe, so protect
       *    the whole section with mutex.
       */
				/* acquire a lock */
      if ((rtn = pthread_mutex_lock(&(wpool->input_lock))) != 0)
	Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
				/* generate a sequence */
      wpool->nseq++;
      if (wpool->nseq > wpool->nsample) 
	{ /* we're done; release input lock, break loop */
	  if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
	    Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
	  break;
	}
      seq = RandomSequence(Alphabet, wpool->randomseq, Alphabet_size, wpool->L);

				/* release the lock */
      if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
	Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

      /* 2. Score the sequence against the model.
       */
      dsq = DigitizeSequence(seq, len);
      if (P7ViterbiSpaceOK(len, hmm->M, mx))
          sc = P7Viterbi(dsq, len, hmm, mx, &tr);
      else
          sc = P7SmallViterbi(dsq, len, hmm, mx, &tr);
      
      free(dsq); 
      free(seq);
      
      /* 3. Save the output; hist and max_score are shared,
       *    so protect this section with the output mutex.
       */
      /* acquire lock on the output queue */
      if ((rtn = pthread_mutex_lock(&(wpool->output_lock))) != 0)
          Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
      /* save output */
      wpool->sc[wpool->nseq-1] = (double) sc;
      Trace_GetAlignmentBounds(tr, 1, NULL, NULL, NULL, NULL,
			       &(wpool->alen[wpool->nseq-1]));
      /* release our lock */
      if ((rtn = pthread_mutex_unlock(&(wpool->output_lock))) != 0)
          Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
      P7FreeTrace(tr);
    }

  StopwatchStop(&thread_watch);
  /* acquire lock on the output queue */
  if ((rtn = pthread_mutex_lock(&(wpool->output_lock))) != 0)
    Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
				/* accumulate cpu time into main stopwatch */
  StopwatchInclude(&(wpool->watch), &thread_watch);
    				/* release our lock */
  if ((rtn = pthread_mutex_unlock(&(wpool->output_lock))) != 0)
    Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

  FreePlan7Matrix(mx);
  pthread_exit(NULL);
  return NULL; /* solely to silence compiler warnings */
}
#endif /* HMMER_THREADS */



#ifdef HMMER_PVM
/* Function: main_loop_pvm()
 * Date:     SRE, Wed Aug 19 13:59:54 1998 [St. Louis]
 *
 * Purpose:  Given an HMM and parameters for synthesizing random
 *           sequences; return a histogram of scores.
 *           (PVM version)  
 *
 * Args:     hmm     - an HMM to calibrate.
 *           seed    - random number seed
 *           N       - number of seqs to synthesize
 *           L       - length of random sequences
 *           lumpsize- # of seqs per slave exchange; controls granularity
 *           sc         - RETURN: unsorted array of <nsample> scores
 *           alen     - RETURN: array of <nsample> optimal ali lengths
 *           extrawatch - RETURN: total CPU time spend in slaves.
 *           ret_nslaves- RETURN: number of PVM slaves run.
 *
 * Returns:  (void)
 *           hist is alloc'ed here, and must be free'd by caller.
 */
static void
main_loop_pvm(struct plan7_s *hmm, int seed, int N, int L, int lumpsize,
	      float *sc, int *alen, Stopwatch_t *extrawatch, int *ret_nslaves)
{
  int                 master_tid;
  int                *slave_tid;
  int                 nslaves;
  int                 nsent;	/* # of seqs we've asked for so far       */
  int                 ndone;	/* # of seqs we've got results for so far */
  int		      packet;	/* # of seqs to have a slave do           */
  int                 slaveidx;	/* id of a slave */
  float              *tmpsc;    /* scores returned by a slave */
  int                *tmplen;   /* ali lengths returned by a slave */
  Stopwatch_t         slavewatch;
  int                 i;
  
  StopwatchZero(extrawatch);

  /* Initialize PVM
   */
  if ((master_tid = pvm_mytid()) < 0)
    Die("pvmd not responding -- do you have PVM running?");
#if DEBUGLEVEL >= 1
  pvm_catchout(stderr);		/* catch output for debugging */
#endif
  PVMSpawnSlaves("hmmcalibrate-pvm", &slave_tid, &nslaves);

  /* Initialize the slaves
   */
  pvm_initsend(PvmDataDefault);
  pvm_pkint(  &L,             1, 1);
  pvm_pkint(  &Alphabet_type, 1, 1);
  pvm_pkint(  &seed,          1, 1);
  if (! PVMPackHMM(hmm)) Die("Failed to pack the HMM");
  pvm_mcast(slave_tid, nslaves, HMMPVM_INIT);
  SQD_DPRINTF1(("Initialized %d slaves\n", nslaves));

  /* Confirm slaves' OK status.
   */
  PVMConfirmSlaves(slave_tid, nslaves);
  SQD_DPRINTF1(("Slaves confirm that they're ok...\n"));
 
  /* Load the slaves
   */
  nsent = ndone = 0;
  for (slaveidx = 0; slaveidx < nslaves; slaveidx++)
    {
      packet    = (nsample - nsent > lumpsize ? lumpsize : nsample - nsent);

      pvm_initsend(PvmDataDefault);
      pvm_pkint(&packet,    1, 1);
      pvm_pkint(&slaveidx,  1, 1);
      pvm_send(slave_tid[slaveidx], HMMPVM_WORK);
      nsent += packet;
    }
  SQD_DPRINTF1(("Loaded %d slaves\n", nslaves));

  /* Receive/send loop
   */
  tmpsc  = MallocOrDie(sizeof(float) * lumpsize);
  tmplen = MallocOrDie(sizeof(int)   * lumpsize);
  while (nsent < nsample)
    {
				/* integrity check of slaves */
      PVMCheckSlaves(slave_tid, nslaves);

				/* receive results */
      SQD_DPRINTF2(("Waiting for results...\n"));
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&slaveidx,   1, 1);
      pvm_upkint(&packet,     1, 1);
      pvm_upkfloat(tmpsc, packet, 1);
      pvm_upkfloat(tmplen, packet, 1);
      SQD_DPRINTF2(("Got results.\n"));

				/* store results */
      for (i = 0; i < packet; i++) {
	sc[ndone]   = (double) tmpsc[i];
	alen[ndone] = tmplen[i];
	ndone++;
      }
				/* send new work */
      packet    = (nsample - nsent > lumpsize ? lumpsize : nsample - nsent);

      pvm_initsend(PvmDataDefault);
      pvm_pkint(&packet,    1, 1);
      pvm_pkint(&slaveidx,  1, 1);
      pvm_send(slave_tid[slaveidx], HMMPVM_WORK);
      SQD_DPRINTF2(("Told slave %d to do %d more seqs.\n", slaveidx, packet));
      nsent += packet;
    }
      
  /* Wait for the last output to come in.
   */
  while (ndone < nsample)
    {
				/* integrity check of slaves */
      PVMCheckSlaves(slave_tid, nslaves);

				/* receive results */
      SQD_DPRINTF1(("Waiting for final results...\n"));
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&slaveidx, 1, 1);
      pvm_upkint(&packet,   1, 1);
      pvm_upkfloat(tmpsc, packet, 1);
      pvm_upkfloat(tmplen, packet, 1);
      SQD_DPRINTF2(("Got some final results.\n"));
				/* store results */
      for (i = 0; i < packet; i++) {
	sc[ndone]   = (double) tmpsc[i];
	alen[ndone] = tmplen[i];
	ndone++;
      }
    }

  /* Shut down the slaves: send -1,-1,-1.
   */
  pvm_initsend(PvmDataDefault);
  packet = -1;
  pvm_pkint(&packet, 1, 1);
  pvm_pkint(&packet, 1, 1);
  pvm_pkint(&packet, 1, 1);
  pvm_mcast(slave_tid, nslaves, HMMPVM_WORK);

  /* Collect stopwatch results; quit the VM; return.
   */
  for (i = 0; i < nslaves; i++)
    {
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&slaveidx, 1, 1);
      StopwatchPVMUnpack(&slavewatch);

      SQD_DPRINTF1(("Slave %d finished; says it used %.2f cpu, %.2f sys\n",
		    slaveidx, slavewatch.user, slavewatch.sys));

      StopwatchInclude(extrawatch, &slavewatch);
    }

  free(slave_tid);
  free(tmpsc);
  free(tmplen);
  pvm_exit();
  *ret_nslaves = nslaves;
  return;
}
#endif /* HMMER_PVM */



/* imean-sd()
 * 
 * Calculate the mean and sd of an array of integer values.
 */
static void
imean_sd(int *v, int N, double *ret_mean, double *ret_sd)
{
  int    i;
  double sum = 0.;
  double sqsum = 0.;
  double mean;
  double var;
  double sd;

  for (i = 0; i < N; i++)
    {
      sum += (double) v[i];
      sqsum += (double) v[i] * (double) v[i];
    }
  mean = sum / (double) N;
  var  = (sqsum - (sum*sum / (double) N)) / ((double) N - 1.);
  sd   = sqrt(var);
  
  if (ret_mean != NULL) *ret_mean = mean;
  if (ret_sd   != NULL) *ret_sd   = sd;
}

/* dmean-sd()
 * 
 * Calculate the mean and sd of an array of double values.
 */
static void
dmean_sd(double *v, int N, double *ret_mean, double *ret_sd)
{
  int    i;
  double sum = 0.;
  double sqsum = 0.;
  double mean;
  double var;
  double sd;

  for (i = 0; i < N; i++)
    {
      sum += v[i];
      sqsum += v[i] * v[i];
    }
  mean = sum / N;
  var  = (sqsum - (sum*sum / N)) / (N - 1.);
  sd   = sqrt(var);
  
  if (ret_mean != NULL) *ret_mean = mean;
  if (ret_sd   != NULL) *ret_sd   = sd;
}



/* dmode()
 * This could surely be optimized.
 */
static double
dmode(double *v, int N, double window)
{
  int i,j;
  int n;
  int nmax;
  double mode;

  nmax = -1;
  mode = 0.;
  for (i = 0; i < N; i++)
    {
      n = 1;
      for (j = i-1; j >= 0 && v[j]-v[i] < window/2.; j--) n++;
      for (j = i+1; j <  N && v[i]-v[j] < window/2.; j++) n++;
      if (n > nmax) { nmax = n; mode = v[i]; }
    }
  return mode;
}

static int
cmp_scores(const void *sp1, const void *sp2)
{
  double sc1;
  double sc2; 
  sc1 = * (double *) sp1;
  sc2 = * (double *) sp2;
  if (sc1 <  sc2) return 1;
  if (sc1 >  sc2) return -1;
  return 0;
}
static void
sort_scores(double *sc, int N)
{
  qsort((void *) sc, N, sizeof(double), cmp_scores);
}

static void
save_fitted_histogram(FILE *hfp, char *name, double *sc, int N, double mu, double lambda)
{
  struct histogram_s *h;     
  int i;

  fprintf(hfp, "HMM: %s\n", name);
  h = AllocHistogram(-200, 200, 100);
  for (i = 0; i < N; i++) AddToHistogram(h, sc[i]);
  ExtremeValueSetHistogram(h, mu, lambda, h->lowscore, h->highscore, 2);
  
  PrintASCIIHistogram(hfp, h);
  FreeHistogram(h);
  fprintf(hfp, "//\n");
}


static void
save_predicted_histogram(FILE *hfp, struct plan7_s *hmm, double *sc, int N, int L)
{
  struct histogram_s *h;     
  double pmu;
  double L1, L2;
  int    i;
  
  fprintf(hfp, "HMM: %s\n", hmm->name);
  if (! (hmm->flags & PLAN7_STATS)) 
    fprintf(hfp, "[No previous EVD parameters set in this model.]\n");
  else
    {
      L1 = EdgeCorrection((double) L,          hmm->kappa, hmm->sigma);
      L2 = EdgeCorrection((double) hmm->Lbase, hmm->kappa, hmm->sigma);
      pmu = hmm->mu + log(L1/L2) / hmm->lambda;

      h = AllocHistogram(-200, 200, 100);
      for (i = 0; i < N; i++) AddToHistogram(h, sc[i]);
      ExtremeValueSetHistogram(h, pmu, hmm->lambda, h->lowscore, h->highscore, 0);
  
      PrintASCIIHistogram(hfp, h);
      FreeHistogram(h);
    }
  fprintf(hfp, "//\n");
}

/* save_score_list()
 * 
 * Saves a file of sorted scores w/ four fields:
 *  <#>  <sc>  <fitted E-val>  <predicted E-val>
 *  
 * where the <fitted E-val> uses mu, lambda that hmmcalibrate just fitted; 
 * and the <predicted E-val> uses the previous statistics stored
 * in the HMM.
 */
static void
save_score_list(FILE *sfp, double *sc, int N, int L, 
		struct plan7_s *hmm, double mu, double lambda)
{
  int i;
  double Epredicted, Efitted;

  fprintf(sfp, "HMM: %s\n", hmm->name);
  for (i = 0; i < N; i++) 
    {
      Efitted    = ExtremeValueE(sc[i], mu, lambda, N);
      Epredicted = LPValue(hmm, L, sc[i]) * (double) N;
      fprintf(sfp, "%-6d %12.2f %16.2f %16.2f\n", 
	      i+1, sc[i], Efitted, Epredicted);
    }
  fprintf(sfp, "//\n");
}

/************************************************************
 * @LICENSE@
 ************************************************************/

