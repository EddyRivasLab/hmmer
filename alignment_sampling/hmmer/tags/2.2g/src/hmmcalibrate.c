/************************************************************
 * @LICENSE@
 ************************************************************/

/* hmmcalibrate.c
 * SRE, Fri Oct 31 09:25:21 1997 [St. Louis]
 * 
 * Score an HMM against random sequence data sets;
 * set histogram fitting parameters.
 * 
 * CVS $Id$
 */

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
#include "config.h"		/* compile-time configuration constants */
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "version.h"		/* release version info                 */
#include "stopwatch.h"		/* process timings                      */

static char banner[] = "hmmcalibrate -- calibrate HMM search statistics";

static char usage[] = "\
Usage: hmmcalibrate [-options] <hmmfile>\n\
Available options are:\n\
  -h             : print short usage and version info, then exit\n\
";

static char experts[] = "\
  --cpu <n>      : run <n> threads in parallel (if threaded)\n\
  --fixed <n>    : fix random sequence length at <n>\n\
  --histfile <f> : save histogram(s) to file <f>\n\
  --mean <x>     : set random seq length mean at <x> [350]\n\
  --num <n>      : set number of sampled seqs to <n> [5000]\n\
  --pvm          : run on a Parallel Virtual Machine (PVM)\n\
  --sd <x>       : set random seq length std. dev to <x> [350]\n\
  --seed <n>     : set random seed to <n> [time()]\n\
";

static struct opt_s OPTIONS[] = {
   { "-h",         TRUE,  sqdARG_NONE  },
   { "--cpu",      FALSE, sqdARG_INT },
   { "--fixed",    FALSE, sqdARG_INT   },
   { "--histfile", FALSE, sqdARG_STRING },
   { "--mean",     FALSE, sqdARG_FLOAT },
   { "--num",      FALSE, sqdARG_INT   },
   { "--pvm",      FALSE, sqdARG_NONE  },
   { "--sd",       FALSE, sqdARG_FLOAT },   
   { "--seed",     FALSE, sqdARG_INT}, 
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


static void main_loop_serial(struct plan7_s *hmm, int seed, int nsample,
			     float lenmean, float lensd, int fixedlen,
			     struct histogram_s **ret_hist, float *ret_max);

#ifdef HMMER_THREADS
/* A structure of this type is shared by worker threads in the POSIX
 * threads parallel version.
 */
struct workpool_s {
  /* Static configuration:
   */
  struct plan7_s  *hmm;		/* ptr to single HMM to search with    */
  int    fixedlen;		/* if >0, fix random seq len to this   */
  float  lenmean;		/* mean of Gaussian for random seq len */
  float  lensd;			/* s.d. of Gaussian for random seq len */
  float *randomseq;             /* 0..Alphabet_size-1 i.i.d. probs     */
  int    nsample;		/* number of random seqs to do         */

  /* Shared (mutex-protected) input:
   */
  int    nseq;			/* current number of seqs searched     */

  /* Shared (mutex-protected) output:
   */
  struct histogram_s *hist;     /* histogram          */
  float          max_score;     /* maximum score seen */
  Stopwatch_t    watch;		/* Timings accumulated for threads */

  /* Thread pool information:
   */
  pthread_t      *thread;       /* our pool of threads */
  int             num_threads;  /* number of threads   */
  pthread_mutex_t input_lock;	/* a mutex protecting input fields */
  pthread_mutex_t output_lock;  /* a mutex protecting output fields */
};
static void main_loop_threaded(struct plan7_s *hmm, int seed, int nsample, 
			       float lenmean, float lensd, int fixedlen,
			       int nthreads,
			       struct histogram_s **ret_hist, float *ret_max,
			       Stopwatch_t *twatch);
static struct workpool_s *workpool_start(struct plan7_s *hmm, 
				 float lenmean, float lensd, int fixedlen,
				 float *randomseq, int nsample, 
				 struct histogram_s *hist, 
				 int num_threads);
static void  workpool_stop(struct workpool_s *wpool);
static void  workpool_free(struct workpool_s *wpool);
static void *worker_thread(void *ptr);
#endif /* HMMER_THREADS */

#ifdef HMMER_PVM
static void main_loop_pvm(struct plan7_s *hmm, int seed, int nsample, 
			  int lumpsize,
			  float lenmean, float lensd, int fixedlen,
			  struct histogram_s **ret_hist, float *ret_max, 
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
  sigset_t blocksigs;		/* list of signals to protect from */
  int     nhmm;			/* number of HMMs calibrated       */

  struct histogram_s *hist;     /* a resulting histogram           */
  float   max;			/* maximum score from an HMM       */
  char   *histfile;             /* histogram save file             */
  FILE   *hfp;                  /* open file pointer for histfile  */

  Stopwatch_t stopwatch;	/* main stopwatch for process       */
  Stopwatch_t extrawatch;	/* stopwatch for threads/PVM slaves */

  float  *mu;			/* array of EVD mu's for HMMs      */
  float  *lambda;		/* array of EVD lambda's for HMMs  */
  int     mu_lumpsize;		/* allocation lumpsize for mu, lambda */

  int     nsample;		/* number of random seqs to sample */
  int     seed;			/* random number seed              */
  int     fixedlen;		/* fixed length, or 0 if unused    */
  float   lenmean;		/* mean of length distribution     */
  float   lensd;		/* std dev of length distribution  */
  int     do_pvm;		/* TRUE to use PVM                 */
  int     pvm_lumpsize;		/* # of seqs to do per PVM slave exchange */
  int     pvm_nslaves;		/* number of slaves used in the PVM */


  char *optname;		/* name of option found by Getopt() */
  char *optarg;			/* argument found by Getopt()       */
  int   optind;		        /* index in argv[]                  */

  int   num_threads;            /* number of worker threads */   


  /***********************************************
   * Parse the command line
   ***********************************************/
  StopwatchStart(&stopwatch);
  StopwatchZero(&extrawatch);

  nsample      = 5000;
  fixedlen     = 0;
  lenmean      = 325.;
  lensd        = 200.;
  seed         = (int) time ((time_t *) NULL);
  histfile     = NULL;
  do_pvm       = FALSE;
  pvm_lumpsize = 20;		/* 20 seqs/PVM exchange: sets granularity */
  mu_lumpsize  = 100;
#ifdef HMMER_THREADS
  num_threads  = ThreadNumber(); /* only matters if we're threaded */
#else
  num_threads  = 0;
#endif

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))
    {
      if      (strcmp(optname, "--cpu")      == 0) num_threads  = atoi(optarg);
      else if (strcmp(optname, "--fixed")    == 0) fixedlen = atoi(optarg);
      else if (strcmp(optname, "--histfile") == 0) histfile = optarg;
      else if (strcmp(optname, "--mean")     == 0) lenmean  = atof(optarg); 
      else if (strcmp(optname, "--num")      == 0) nsample  = atoi(optarg); 
      else if (strcmp(optname, "--pvm")      == 0) do_pvm   = TRUE;
      else if (strcmp(optname, "--sd")       == 0) lensd    = atof(optarg); 
      else if (strcmp(optname, "--seed")     == 0) seed     = atoi(optarg);
      else if (strcmp(optname, "-h") == 0)
	{
	  Banner(stdout, banner);
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

				/* histogram file */
  hfp = NULL;
  if (histfile != NULL) {
    if ((hfp = fopen(histfile, "w")) == NULL)
      Die("Failed to open histogram save file %s for writing\n", histfile);
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

  Banner(stdout, banner);
  printf("HMM file:                 %s\n", hmmfile);
  if (fixedlen) 
    printf("Length fixed to:          %d\n", fixedlen);
  else {
    printf("Length distribution mean: %.0f\n", lenmean);
    printf("Length distribution s.d.: %.0f\n", lensd);
  }
  printf("Number of samples:        %d\n", nsample);
  printf("random seed:              %d\n", seed);
  printf("histogram(s) saved to:    %s\n",
	 histfile != NULL ? histfile : "[not saved]");
  if (do_pvm)
    printf("PVM:                      ACTIVE\n");
  else if (num_threads > 0)
    printf("POSIX threads:            %d\n", num_threads);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  /***********************************************
   * Read the HMMs one at a time, and send them off
   * in probability form to one of the main loops.
   * The main loop functions are responsible for 
   * synthesizing random sequences and returning
   * a score histogram for each HMM.
   ***********************************************/

  nhmm = 0;
  mu     = MallocOrDie(sizeof(float) * mu_lumpsize);
  lambda = MallocOrDie(sizeof(float) * mu_lumpsize);

  while (HMMFileRead(hmmfp, &hmm))
    {
      if (hmm == NULL)
	Die("HMM file may be corrupt or in incorrect format; parse failed");

      if (! do_pvm && num_threads == 0)
	main_loop_serial(hmm, seed, nsample, lenmean, lensd, fixedlen, 
			 &hist, &max);
#ifdef HMMER_PVM
      else if (do_pvm) {
	pvm_nslaves = 0;	/* solely to silence compiler warnings */
	main_loop_pvm(hmm, seed, nsample, pvm_lumpsize, 
		      lenmean, lensd, fixedlen, 
		      &hist, &max, &extrawatch, &pvm_nslaves);
      }
#endif 
#ifdef HMMER_THREADS
      else if (num_threads > 0)
	main_loop_threaded(hmm, seed, nsample, lenmean, lensd, fixedlen,
			   num_threads, &hist, &max, &extrawatch);
#endif
      else 
	Die("wait. that can't happen. I didn't do anything.");


      /* Fit an EVD to the observed histogram.
       * The TRUE left-censors and fits only the right slope of the histogram.
       * The 9999. is an arbitrary high number that means we won't trim
       * outliers on the right.
       */
      if (! ExtremeValueFitHistogram(hist, TRUE, 9999.))
	Die("fit failed; -n may be set too small?\n");
      
      mu[nhmm]     = hist->param[EVD_MU];
      lambda[nhmm] = hist->param[EVD_LAMBDA];
      nhmm++;
      if (nhmm % 100 == 0) {
	mu     = ReallocOrDie(mu,     sizeof(float) * (nhmm+mu_lumpsize));
	lambda = ReallocOrDie(lambda, sizeof(float) * (nhmm+mu_lumpsize));
      }      

      /* Output
       */
      printf("HMM    : %s\n",   hmm->name);
      printf("mu     : %12f\n", hist->param[EVD_MU]);
      printf("lambda : %12f\n", hist->param[EVD_LAMBDA]);
      printf("max    : %12f\n", max);
      printf("//\n");

      if (hfp != NULL) 
	{
	  fprintf(hfp, "HMM: %s\n", hmm->name);
	  PrintASCIIHistogram(hfp, hist);
	  fprintf(hfp, "//\n");
	}

      FreeHistogram(hist);
    }
  SQD_DPRINTF1(("Main body believes it has calibrations for %d HMMs\n", nhmm));

  /*****************************************************************
   * Rewind the HMM file for a second pass.
   * Write a temporary HMM file with new mu, lambda values in it
   *****************************************************************/

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
  if (fclose(outfp)   != 0) PANIC;

  if (sigemptyset(&blocksigs) != 0) PANIC;
  if (sigaddset(&blocksigs, SIGINT) != 0) PANIC;
  if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   PANIC;
  if (remove(hmmfile) != 0)                            PANIC;
  if (rename(tmpfile, hmmfile) != 0)                   PANIC;
  if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) PANIC;

  /***********************************************
   * Exit
   ***********************************************/

  StopwatchStop(&stopwatch);
  if (do_pvm > 0) {
    printf("PVM processors used: %d\n", pvm_nslaves);
    StopwatchInclude(&stopwatch, &extrawatch);
  }
#ifdef PTHREAD_TIMES_HACK
  else if (num_threads > 0) StopwatchInclude(&stopwatch, &extrawatch);
#endif

  /*  StopwatchDisplay(stdout, "CPU Time: ", &stopwatch); */

  free(mu);
  free(lambda);
  free(tmpfile);
  if (hfp != NULL) fclose(hfp);
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
 *           nsample  - number of seqs to synthesize
 *           lenmean  - mean length of random sequence
 *           lensd    - std dev of random seq length
 *           fixedlen - if nonzero, override lenmean, always this len
 *           ret_hist - RETURN: the score histogram 
 *           ret_max  - RETURN: highest score seen in simulation
 *
 * Returns:  (void)
 *           hist is alloc'ed here, and must be free'd by caller.
 */
static void
main_loop_serial(struct plan7_s *hmm, int seed, int nsample, 
		 float lenmean, float lensd, int fixedlen,
		 struct histogram_s **ret_hist, float *ret_max)
{
  struct histogram_s *hist;
  float  randomseq[MAXABET];
  float  p1;
  float  max;
  char  *seq;
  char  *dsq;
  float  score;
  int    sqlen;
  int    idx;
  
  /* Initialize.
   * We assume we've already set the alphabet (safe, because
   * HMM input sets the alphabet).
   */
  sre_srandom(seed);
  P7Logoddsify(hmm, TRUE);
  P7DefaultNullModel(randomseq, &p1);
  hist = AllocHistogram(-200, 200, 100);
  max = -FLT_MAX;

  for (idx = 0; idx < nsample; idx++)
    {
				/* choose length of random sequence */
      if (fixedlen) sqlen = fixedlen;
      else do sqlen = (int) Gaussrandom(lenmean, lensd); while (sqlen < 1);
				/* generate it */
      seq = RandomSequence(Alphabet, randomseq, Alphabet_size, sqlen);
      dsq = DigitizeSequence(seq, sqlen);

      if (P7ViterbiSize(sqlen, hmm->M) <= RAMLIMIT)
	score = P7Viterbi(dsq, sqlen, hmm, NULL);
      else
	score = P7SmallViterbi(dsq, sqlen, hmm, NULL);

      AddToHistogram(hist, score);
      if (score > max) max = score;

      free(dsq); 
      free(seq);
    }

  *ret_hist   = hist;
  *ret_max    = max;
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
 *           nsample  - number of seqs to synthesize
 *           lenmean  - mean length of random sequence
 *           lensd    - std dev of random seq length
 *           fixedlen - if nonzero, override lenmean, always this len
 *           nthreads - number of threads to start
 *           ret_hist - RETURN: the score histogram 
 *           ret_max  - RETURN: highest score seen in simulation
 *           twatch   - RETURN: accumulation of thread times
 *
 * Returns:  (void)
 *           hist is alloc'ed here, and must be free'd by caller.
 */
static void
main_loop_threaded(struct plan7_s *hmm, int seed, int nsample, 
		   float lenmean, float lensd, int fixedlen,
		   int nthreads,
		   struct histogram_s **ret_hist, float *ret_max,
		   Stopwatch_t *twatch)
{
  struct histogram_s *hist;
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
  hist = AllocHistogram(-200, 200, 100);

  wpool = workpool_start(hmm, lenmean, lensd, fixedlen, randomseq, nsample,
			 hist, nthreads);
  workpool_stop(wpool);

  *ret_hist = hist;
  *ret_max  = wpool->max_score;
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
 *           num_threads - how many processors to run on
 *
 * Returns:  ptr to struct workpool_s.
 *           Caller must wait for threads to finish with workpool_stop(),
 *           then free the structure with workpool_free().
 */
static struct workpool_s *
workpool_start(struct plan7_s *hmm, float lenmean, float lensd, int fixedlen,
	       float *randomseq, int nsample, struct histogram_s *hist, 
	       int num_threads)
{
  struct workpool_s *wpool;
  pthread_attr_t    attr;
  int i;
  int rtn;

  wpool         = MallocOrDie(sizeof(struct workpool_s));
  wpool->thread = MallocOrDie(num_threads * sizeof(pthread_t));
  wpool->hmm        = hmm;
  wpool->fixedlen   = fixedlen;
  wpool->lenmean    = lenmean;
  wpool->lensd      = lensd;
  wpool->randomseq  = randomseq;
  wpool->nsample    = nsample;
  
  wpool->nseq       = 0;
  wpool->hist       = hist;
  wpool->max_score  = -FLT_MAX;
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
  struct workpool_s *wpool;
  char       *seq;
  char       *dsq;
  int         len;
  float       sc;
  int         rtn;
  Stopwatch_t thread_watch;

  StopwatchStart(&thread_watch);
  wpool = (struct workpool_s *) ptr;
  hmm   = wpool->hmm;
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
      if (wpool->fixedlen) len = wpool->fixedlen;
      else do len = (int) Gaussrandom(wpool->lenmean, wpool->lensd); while (len < 1);
      seq = RandomSequence(Alphabet, wpool->randomseq, Alphabet_size, len);

				/* release the lock */
      if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
	Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

      /* 2. Score the sequence against the model.
       */
      dsq = DigitizeSequence(seq, len);
      
      if (P7ViterbiSize(len, hmm->M) <= RAMLIMIT)
	sc = P7Viterbi(dsq, len, hmm, NULL);
      else
	sc = P7SmallViterbi(dsq, len, hmm, NULL);
      free(dsq); 
      free(seq);
      
      /* 3. Save the output; hist and max_score are shared,
       *    so protect this section with the output mutex.
       */
				/* acquire lock on the output queue */
      if ((rtn = pthread_mutex_lock(&(wpool->output_lock))) != 0)
	Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
				/* save output */
      AddToHistogram(wpool->hist, sc);
      if (sc > wpool->max_score) wpool->max_score = sc;
    				/* release our lock */
      if ((rtn = pthread_mutex_unlock(&(wpool->output_lock))) != 0)
	Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
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
 *           nsample - number of seqs to synthesize
 *           lumpsize- # of seqs per slave exchange; controls granularity
 *           lenmean - mean length of random sequence
 *           lensd   - std dev of random seq length
 *           fixedlen- if nonzero, override lenmean, always this len
 *           hist       - RETURN: the score histogram 
 *           ret_max    - RETURN: highest score seen in simulation
 *           extrawatch - RETURN: total CPU time spend in slaves.
 *           ret_nslaves- RETURN: number of PVM slaves run.
 *
 * Returns:  (void)
 *           hist is alloc'ed here, and must be free'd by caller.
 */
static void
main_loop_pvm(struct plan7_s *hmm, int seed, int nsample, int lumpsize,
	      float lenmean, float lensd, int fixedlen,
	      struct histogram_s **ret_hist, float *ret_max, 
	      Stopwatch_t *extrawatch, int *ret_nslaves)
{
  struct histogram_s *hist;
  int                 master_tid;
  int                *slave_tid;
  int                 nslaves;
  int                 nsent;	/* # of seqs we've asked for so far       */
  int                 ndone;	/* # of seqs we've got results for so far */
  int		      packet;	/* # of seqs to have a slave do           */
  float               max;
  int                 slaveidx;	/* id of a slave */
  float              *sc;        /* scores returned by a slave */
  Stopwatch_t         slavewatch;
  int                 i;
  
  StopwatchZero(extrawatch);
  hist = AllocHistogram(-200, 200, 100);
  max  = -FLT_MAX;

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
  pvm_pkfloat(&lenmean,       1, 1);
  pvm_pkfloat(&lensd,         1, 1);
  pvm_pkint(  &fixedlen,      1, 1);
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
  sc = MallocOrDie(sizeof(float) * lumpsize);
  while (nsent < nsample)
    {
				/* integrity check of slaves */
      PVMCheckSlaves(slave_tid, nslaves);

				/* receive results */
      SQD_DPRINTF2(("Waiting for results...\n"));
      pvm_recv(-1, HMMPVM_RESULTS);
      pvm_upkint(&slaveidx,   1, 1);
      pvm_upkint(&packet,     1, 1);
      pvm_upkfloat(sc,   packet, 1);
      SQD_DPRINTF2(("Got results.\n"));
      ndone += packet;

				/* store results */
      for (i = 0; i < packet; i++) {
	AddToHistogram(hist, sc[i]);
	if (sc[i] > max) max = sc[i];
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
      pvm_upkfloat(sc, packet, 1);
      SQD_DPRINTF2(("Got some final results.\n"));
      ndone += packet;
				/* store results */
      for (i = 0; i < packet; i++) {
	AddToHistogram(hist, sc[i]);
	if (sc[i] > max) max = sc[i];
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
  free(sc);
  pvm_exit();
  *ret_hist    = hist;
  *ret_max     = max;
  *ret_nslaves = nslaves;
  return;
}
#endif /* HMMER_PVM */



