/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* threads.c
 * SRE, Fri Jul 10 10:05:44 1998
 * 
 * Pthreads code shared by hmmsearch, hmmcalibrate, and hmmpfam
 * to coarse-grain parallelize on platforms capable of POSIX
 * threads.
 * 
 *****************************************************************
 * API to the pthreads implementation:
 * 
 *  vpool = VpoolInit(FALSE, FALSE, num_threads, num_threads);
 *
 *  [for all (digitized) sequences, HMMs, or whatever:]
 *      VpoolAddWork(vpool, hmm, dsq, sqinfo, sqinfo.len);
 *      else if [no work to add] VpoolShutdown(vpool);
 *  
 *      while (VpoolGetResults(vpool, &hmm, &dsq, &sqinfo, &len, &score, &tr)) 
 *           [process results: add to histograms, etc.]
 *           [free any data that is being looped over and overwritten]
 *  
 *  VpoolDestroy(vpool);
 *           
 *****************************************************************     
 * Overview of the internals of the pthreads implementation:
 * 
 * All the business is in struct vpool_s, including the
 * thread handles and the input and output queues.
 * The input/output queues are a double buffered system.
 *                    
 *     VpoolInit()       - allocates the vpool_s, starts threads.
 *     VpoolShutdown()   - stops the threads, fills up the output queue.
 *     VpoolDestroy()    - frees the vpool_s after output queue has been emptied.
 *     VpoolThread()     - what the threads run: takes seq off input queue,
 *                         runs P7Viterbi() on it, puts results on output
 *                         queue.
 *     VpoolAddWork()    - put one new sequence on the input queue.
 *     VpoolGetResults() - take one result off the output queue.
 *                       
 * 
 * RCS $Id$
 */

#ifdef HMMER_THREADS		/* conditional inclusion of the entire file */

#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <sched.h>		/* apparently SGI IRIX pthread.h is not complete */
#include <assert.h>

#include "structs.h"
#include "funcs.h"
#include "squid.h"
#include "sqfuncs.h"


/* Function: ThreadNumber()
 * Date:     SRE, Sat Jul 11 11:03:50 1998 [St. Louis]
 *
 * Purpose:  Recommend how many threads to use. 
 *
 *             - if we can determine the number of processors
 *               on the machine by SQD_NPROC, use that. This
 *               should succeed for SGI IRIX, Digital UNIX, and 
 *               Sun Solaris platforms.
 *             - if not, assume two processors. We're probably
 *               on a FreeBSD or Linux box, and odds are that its
 *               a dualprocessor.
 *             - if HMMER_NCPU is defined in config.h, use that
 *               number instead; allows Linux or FreeBSD machines
 *               to compile code for a quadprocessor, for instance.
 *               That define can be overridden at compile
 *               time by a -DHMMER_NCPU=x, where x is the
 *               number of threads..
 *             - if HMMER_NCPU is defined in the environment,
 *               use that number, overriding all others.
 *
 *           Typically, we'll set the default number of
 *           threads with ThreadNumber() but allow it
 *           to be overridden at the command line with --cpu.    
 *           
 *           Summarizing priority:
 *                  --ncpu <x> option
 *                  environment variable, setenv HMMER_NCPU x
 *                  compile-time, MDEFS=HMMER_NCPU=x
 *                  compile-time, config.h definition of HMMER_NCPU
 *                  SQD_NPROC, or 2 if SQD_NPROC doesn't work.
 *
 * Args:     void
 *
 * Returns:  >= 1, recommended number of threads
 */
int
ThreadNumber(void)
{
  int   num;
  char *env;

  num = SQD_NPROC;		/* SGI, Sun, Digital: get # of available CPUs */
  if (num == -1) num = 2;	/* Linux, FreeBSD: assume dualprocessor       */
#ifdef HMMER_NCPU	
  num = HMMER_NCPU;		/* allow config.h to override; usually we don't */
#endif
				/* allow environment variable to override */
  if ((env = getenv("HMMER_NCPU")) != NULL)
    num = atoi(env);		
  if (num <= 0) num = 1;	/* silent sanity check */
  SQD_DPRINTF1(("ThreadNumber(): setting number of threads to %d\n", num));
  return num;
}



/* Function: VpoolInit()
 * Date:     SRE, Fri Jul 10 11:04:16 1998 [St. Louis]
 *
 * Purpose:  Allocates and initializes the worker thread pool.
 *
 * Args:     do_forward  - TRUE to use P7Forward() to calculate scores (us. FALSE)
 *           do_null     - TRUE to correct scores with null2 model (us. TRUE)
 *           num_threads - number of threads (us. # of CPUs)
 *           max_queue   - size of input queue in double buffer (us. == num_threads)
 *
 * Returns:  ptr to allocated/initialized vpool_s structure.
 *           Caller (boss thread) must free with VpoolDestroy().
 */
struct vpool_s *
VpoolInit(int do_forward, int do_null, int num_threads, int max_queue)
{
  struct vpool_s *vpool;
  int i;
  int rtn;

  vpool = MallocOrDie(sizeof(struct vpool_s));
  vpool->do_forward       = do_forward;
  vpool->do_null          = do_null;
  vpool->num_threads      = num_threads;
  vpool->max_input_queue  = max_queue;
  vpool->max_output_queue = max_queue + num_threads;

  vpool->thread  = MallocOrDie(sizeof(pthread_t)          * num_threads);
  vpool->nin     = 0;
  vpool->hmm     = MallocOrDie(sizeof(struct plan7_s *)   * vpool->max_input_queue);
  vpool->dsq     = MallocOrDie(sizeof(char *)             * vpool->max_input_queue);
  vpool->sqinfo  = MallocOrDie(sizeof(SQINFO *)           * vpool->max_input_queue);
  vpool->len     = MallocOrDie(sizeof(int)                * vpool->max_input_queue);
  vpool->nout    = 0;
  vpool->ohmm    = MallocOrDie(sizeof(struct plan7_s *)   * vpool->max_output_queue);
  vpool->odsq    = MallocOrDie(sizeof(char *)             * vpool->max_output_queue);
  vpool->osqinfo = MallocOrDie(sizeof(SQINFO *)           * vpool->max_output_queue);
  vpool->olen    = MallocOrDie(sizeof(int)                * vpool->max_output_queue);
  vpool->score   = MallocOrDie(sizeof(float)              * vpool->max_output_queue);
  vpool->tr      = MallocOrDie(sizeof(struct p7trace_s * )* vpool->max_output_queue);

  vpool->shutdown = 0;

  if ((rtn = pthread_mutex_init(&(vpool->input_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));
  if ((rtn = pthread_mutex_init(&(vpool->output_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));
  if ((rtn = pthread_cond_init(&(vpool->input_ready), NULL)) != 0)
    Die("pthread_cond_init FAILED; %s\n", strerror(rtn));
  if ((rtn = pthread_cond_init(&(vpool->output_ready), NULL)) != 0)
    Die("pthread_cond_init FAILED; %s\n", strerror(rtn));

  /* Create slave threads
   */
  for (i = 0; i < num_threads; i++)
    if ((rtn = pthread_create(&(vpool->thread[i]), NULL,
			      VpoolThread , (void *) vpool)) != 0)
      Die("Failed to create thread %d; return code %d\n", i, rtn);

  return vpool;
}



/* Function: VpoolShutdown()
 * Date:     SRE, Fri Jul 10 13:57:08 1998 [St. Louis]
 *
 * Purpose:  Orders the worker threads to shutdown, then
 *           waits for them to do so. 
 *           
 *           This is called when the boss runs out of input.
 *           After VpoolShutdown(), the boss processes the
 *           output queue one last time, then can call VpoolDestroy().
 *           
 *           Output queue must be flushed recently; no calls to VpoolAddWork()
 *           between that time and VpoolShutdown(), else may overflow the
 *           output queue.
 *
 * Args:     vpool  - thread pool to shut down
 *
 * Returns:  void
 */
void
VpoolShutdown(struct vpool_s *vpool)
{
  int i;

			/* young Zaphod plays it safe, and acquires a lock */
  if (pthread_mutex_lock(&(vpool->input_lock)) != 0)
    Die("pthread_mutex_lock failed");
			/* toggle the shutdown flag */
  SQD_DPRINTF2(("VpoolShutdown() says, time to die\n"));
  vpool->shutdown = 1;
			/* poke sleeping workers to make them look at shutdown */
  SQD_DPRINTF2(("VpoolShutdown() says, this means you, slave\n"));
  if (pthread_cond_broadcast(&(vpool->input_ready)) != 0) 
    Die("pthread_cond_broadcast failed");
			/* release our lock so workers can get remaining input */
  if (pthread_mutex_unlock(&(vpool->input_lock)) != 0)
    Die("pthread_mutex_unlock failed");
			/* wait for our workers to finish */
  SQD_DPRINTF2(("VpoolShutdown() says everybody out of the pool!\n"));
  for (i = 0; i < vpool->num_threads; i++)
    if (pthread_join(vpool->thread[i],NULL) != 0)
      Die("pthread_join failed");
}


/* Function: VpoolDestroy()
 * Date:     SRE, Fri Jul 10 13:24:53 1998 [St. Louis]
 *
 * Purpose:  Frees memory associated with a vpool_s.
 *           Caller must be sure that i/o queues have been flushed,
 *           and threads have terminated.
 *
 * Args:     vpool - thread pool to destroy.
 *
 * Returns:  (void)
 */
void
VpoolDestroy(struct vpool_s *vpool)
{
  SQD_DASSERT1((vpool->nin != 0));  /* sanity: don't die before input is consumed */
  SQD_DASSERT1((vpool->nout != 0)); /* sanity: don't die before output is consumed */

  free(vpool->thread);
  free(vpool->hmm);
  free(vpool->dsq);
  free(vpool->sqinfo);
  free(vpool->len);
  free(vpool->ohmm);
  free(vpool->odsq);
  free(vpool->osqinfo);
  free(vpool->olen);
  free(vpool->score);
  free(vpool->tr);
  free(vpool);
}

           


/* Function: VpoolThread()
 * Date:     SRE, Fri Jul 10 11:17:29 1998 [St. Louis]
 *
 * Purpose:  The procedure executed by the worker threads.
 *           Gets seq/HMM from input queue (waits for input
 *           or shutdown signal if input isn't immediately
 *           available); runs alignment algorithms to
 *           get score/trace; puts results on output queue.
 *
 * Args:     vpool   - ptr to the thread pool, with i/o queue
 *
 * Returns:  (void)
 */
void *
VpoolThread(void *ptr)
{
  struct plan7_s   *hmm;
  struct p7trace_s *tr;
  struct vpool_s   *vpool;
  char   *dsq;
  int     len;
  SQINFO *sqinfo;
  float   sc;
  int     rtn;

  /* The only way that the worker thread terminates is when it
   * sees a shutdown flag from the boss (which tells it that
   * the boss has no more data for it), and the data in the
   * input queue have run out.
   */
  vpool = (struct vpool_s *) ptr;
  for (;;) {

    /*****************************************************************
     * Get the next sequence to work on.
     * If none are available, wait for one to appear from the boss.
     *    (boss signals a condition input_ready)
     *****************************************************************/
                                 /* acquire lock on the input queue. */
    SQD_DPRINTF2(("VpoolThread(): thread %x wants the input lock!\n",  
		      (unsigned int) pthread_self()));
    if ((rtn = pthread_mutex_lock(&(vpool->input_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    SQD_DPRINTF2(("VpoolThread(): thread %x gets the input lock.\n", 
		      (unsigned int) pthread_self()));

				/* wait for input seqs to be available */
    while ((vpool->nin == 0) && (! vpool->shutdown)) 
      {
	SQD_DPRINTF2(("VpoolThread(): thread %x blocking, waiting for input\n", 
		      (unsigned int) pthread_self()));
	if ((rtn = pthread_cond_wait(&(vpool->input_ready), &(vpool->input_lock))) != 0)
	  Die("pthread_cond_wait failure: %s\n", strerror(rtn));
	SQD_DPRINTF2(("VpoolThread(): thread %x sees a signal, and gets the input lock.\n", 
		      (unsigned int) pthread_self()));
      }
    SQD_DPRINTF2(("VpoolThread(): thread %x taking input from the queue, going to work\n", 
		  (unsigned int) pthread_self()));

				/* check for shutdown request, and lack of input data */
    if (vpool->nin == 0 && vpool->shutdown) {
      SQD_DPRINTF2(("VPoolThread(): thread %x shutting down normally\n", 
		      (unsigned int) pthread_self()));
      if ((rtn = pthread_mutex_unlock(&(vpool->input_lock))) != 0)
	Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
      pthread_exit(NULL);
    }
				/* get HMM and sequence off the queue */
    vpool->nin--;
    hmm    = vpool->hmm[vpool->nin];
    dsq    = vpool->dsq[vpool->nin];
    len    = vpool->len[vpool->nin];
    sqinfo = vpool->sqinfo[vpool->nin];

				/* release our lock on the input queue */
    if ((rtn = pthread_mutex_unlock(&(vpool->input_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));


    /*****************************************************************
     * Do the work; recover a score and a trace.
     *****************************************************************/

				/* recover Viterbi trace */
    if (P7ViterbiSize(len, hmm->M) <= RAMLIMIT)
      sc = P7Viterbi(dsq, len, hmm, &tr);
    else
      sc = P7SmallViterbi(dsq, len, hmm, &tr);

				/* maybe rescore with Forward */
    if (vpool->do_forward) 
      sc = P7Forward(dsq, len, hmm, NULL);
				/* maybe correct score with null2 */
    if (vpool->do_null)
      sc -= TraceScoreCorrection(hmm, tr, dsq);
    
    /*****************************************************************
     * Put the results on the output queue
     * In theory, output queue cannot overflow: in the boss, each
     *   time we add a sequence to the input queue, we check the
     *   whole output queue.
     *****************************************************************/

				/* acquire lock on the output queue */
    if ((rtn = pthread_mutex_lock(&(vpool->output_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));

				/* sanity check on the output queue */
    SQD_DASSERT1((vpool->nout < vpool->max_output_queue));
    SQD_DPRINTF2(("VpoolThread(): thread %x is done, putting output on the queue\n", 
		  (unsigned int) pthread_self()));

    				/* put results on the output queue */
    vpool->ohmm[vpool->nout]    = hmm;
    vpool->odsq[vpool->nout]    = dsq;
    vpool->olen[vpool->nout]    = len;
    vpool->osqinfo[vpool->nout] = sqinfo;
    vpool->tr[vpool->nout]      = tr;
    vpool->score[vpool->nout]   = sc;
    vpool->nout++;
				/* tell the boss we're done */
    if (pthread_cond_signal(&(vpool->output_ready)) != 0)
      Die("pthread_cond_signal failed");
				/* release our lock */
    if ((rtn = pthread_mutex_unlock(&(vpool->output_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
  }
}


/* Function: VpoolAddWork()
 * Date:     SRE, Fri Jul 10 12:10:03 1998 [St. Louis]
 *
 * Purpose:  Add hmm/sequence data to the input queue of 
 *           the thread pool.
 *
 * Args:     vpool    - thread pool
 *           hmm      - model to align to
 *           dsq      - digitized sequence
 *           sqinfo   - SQINFO structure w/ seq name (may be NULL)
 *           len      - length of dsq
 *
 * Returns:  void
 */
void
VpoolAddWork(struct vpool_s *vpool, struct plan7_s *hmm, char *dsq, SQINFO *sqinfo, int len)
{
				/* acquire lock on input queue */
  SQD_DPRINTF2(("VpoolAddWork(): give me the input lock!\n"));
  if (pthread_mutex_lock(&(vpool->input_lock)) != 0)
    Die("pthread_mutex_lock failed");
  SQD_DPRINTF2(("VpoolAddWork(): thanks, I've got the input lock.\n"));

				/* sanity check: queue is never full */
  SQD_DASSERT1((vpool->nin < vpool->max_input_queue));
  SQD_DPRINTF2(("VpoolAddWork(): putting work on the input queue\n"));

				/* add work to input queue */
  vpool->hmm[vpool->nin]    = hmm;
  vpool->dsq[vpool->nin]    = dsq;
  vpool->sqinfo[vpool->nin] = sqinfo;
  vpool->len[vpool->nin]    = len;
  vpool->nin++;
				/* wake up slaves */
  SQD_DPRINTF2(("VpoolAddWork(): wake up, slaves!\n"));
  if (pthread_cond_signal(&(vpool->input_ready)) != 0)
    Die("pthread_cond_signal failed");
				/* release lock on the input queue */
  if (pthread_mutex_unlock(&(vpool->input_lock)) != 0)
    Die("pthread_mutex_unlock failed");

  return;
}


/* Function: VpoolGetResults()
 * Date:     SRE, Fri Jul 10 13:34:33 1998 [St. Louis]
 *
 * Purpose:  Pop results off the output queue, if any.
 *           Since the logic in our threads implementation
 *           requires completely flushing the output
 *           queue, usually we do this in a while(VpoolGetResults)
 *           loop.
 *
 * Args:     vpool      - thread pool
 *           ret_hmm    - HMM 
 *           ret_dsq    - digitized sequence 
 *           ret_sqinfo - optional SQINFO information (w/ seq name) or NULL
 *           ret_len    - length of dsq
 *           ret_score  - log-odds score of sequence, in bits
 *           ret_tr     - traceback for sequence
 *
 * Returns:  1 if results are returned. 0 if we're shutting down and output queue is empty.
 */
int
VpoolGetResults(struct vpool_s *vpool, struct plan7_s **ret_hmm,
		char **ret_dsq, SQINFO **ret_sqinfo, int *ret_len,
		float *ret_score, struct p7trace_s **ret_tr)
{
				/* acquire a lock on the output queue */
  if (pthread_mutex_lock(&(vpool->output_lock)) != 0)
    Die("pthread_mutex_lock failed");

  if (vpool->shutdown && vpool->nout == 0) return 0;

				/* block and wait until output is in queue */
  while (vpool->nout == 0) 
    {	
      SQD_DPRINTF2(("VpoolGetResults(): no output ready; twiddling my thumbs\n"));
      if (pthread_cond_wait(&(vpool->output_ready), &(vpool->output_lock)) != 0)
	Die("pthread_cond_wait failed");
    }
				/* pop results off the queue */
  vpool->nout--;
  if (ret_hmm    != NULL) *ret_hmm    = vpool->ohmm[vpool->nout];
  if (ret_dsq    != NULL) *ret_dsq    = vpool->odsq[vpool->nout];
  if (ret_sqinfo != NULL) *ret_sqinfo = vpool->osqinfo[vpool->nout];
  if (ret_len    != NULL) *ret_len    = vpool->olen[vpool->nout];
  if (ret_score  != NULL) *ret_score  = vpool->score[vpool->nout];
  if (ret_tr     != NULL) *ret_tr     = vpool->tr[vpool->nout];
  else                    P7FreeTrace(vpool->tr[vpool->nout]);

  SQD_DPRINTF2(("VpoolGetResults(): took some output off the queue\n"));

				/* release lock on the output queue */
  if (pthread_mutex_unlock(&(vpool->output_lock)) != 0)
    Die("pthread_mutex_unlock failed");

  return 1;
}


/* Function: VpoolPrintInputQueue()
 * Date:     SRE, Fri Jul 10 16:14:33 1998 [St. Louis]
 *
 * Purpose:  Debugging.
 *           Print the input queue.
 *
 * Args:     vpool  - the thread pool
 *
 * Returns:  void
 */
void
VpoolPrintInputQueue(struct vpool_s *vpool)
{
  int i;
				/* acquire lock on input queue */
  if (pthread_mutex_lock(&(vpool->input_lock)) != 0)
    Die("pthread_mutex_lock failed");

  printf("Input queue has %d sequences:\n", vpool->nin);
  for (i = 0; i < vpool->nin; i++)
    printf("%d: %s\n", i, DedigitizeSequence(vpool->dsq[i], 10));
  fflush(stdout);

  if (pthread_mutex_unlock(&(vpool->input_lock)) != 0)
    Die("pthread_mutex_unlock failed");
}

/* Function: VpoolPrintOutputQueue()
 * Date:     SRE, Fri Jul 10 16:22:51 1998 [St. Louis]
 *
 * Purpose:  Debugging.
 *           Print the output queue.
 *
 * Args:     vpool  - the thread pool
 *
 * Returns:  void
 */
void
VpoolPrintOutputQueue(struct vpool_s *vpool)
{
  int i;
				/* acquire lock on input queue */
  if (pthread_mutex_lock(&(vpool->output_lock)) != 0)
    Die("pthread_mutex_lock failed");

  printf("Output queue has %d sequences:\n", vpool->nout);
  for (i = 0; i < vpool->nout; i++)
    printf("%d: %s\n", i, DedigitizeSequence(vpool->odsq[i], 10));
  fflush(stdout);

  if (pthread_mutex_unlock(&(vpool->output_lock)) != 0)
    Die("pthread_mutex_unlock failed");
}


#endif /*HMMER_THREADS*/
