/************************************************************
 * @LICENSE@
 ************************************************************/

/* threads.c
 * SRE, Fri Jul 10 10:05:44 1998
 * 
 * Pthreads code shared by hmmsearch, hmmcalibrate, and hmmpfam
 * to coarse-grain parallelize on platforms capable of POSIX
 * threads. Most of the threads code, however, is in the respective
 * main's, i.e. hmmsearch.c, hmmpfam.c, hmmcalibrate.c
 * 
 * RCS $Id$
 */

#include "config.h"
#include "squidconf.h"

#ifdef HMMER_THREADS		/* conditional inclusion of the entire file */
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <pthread.h>

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

#endif /*HMMER_THREADS*/
