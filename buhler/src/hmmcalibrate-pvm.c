/* hmmcalibrate-pvm.c
 * PVM slave for hmmcalibrate.
 * Redesigned for better parallelization: SRE, Wed Dec  1 09:48:58 1999
 *
 * Design:
 *   Initialization: 
 *       receive parameters of random sequence synthesis, and an HMM.
 *       send an OK signal to the master.
 *   
 *   Main loop:
 *       receive work packet: # of seqs to make
 *       Synthesize and score # seqs
 *       send results: # raw scores.
 *   
 *   Termination: 
 *       master sends a shutdown signal instead of a work packet.
 * 
 * SRE, Tue Aug 18 15:19:28 1998
 * SVN $Id$
 */

#include "config.h"		/* compile-time configuration constants */
#include "squidconf.h"

#ifdef HMMER_PVM
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <pvm3.h>

#include "squid.h"		/* general sequence analysis library    */
#include "stopwatch.h"		/* CPU timing routines                  */

#include "plan7.h"		/* plan7 profile HMM structure          */
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */

static void leave_pvm(void);

int 
main(void)
{
  int      master_tid;		/* PVM TID of our master */
  int      slaveidx;		/* my slave index (0..nslaves-1) */
  struct plan7_s *hmm;		/* HMM to calibrate, sent from master */
  struct p7trace_s *tr;         /* traceback from an alignment */
  struct dpmatrix_s *mx;        /* growable DP matrix */
  char              *seq;	/* synthetic random sequence */
  unsigned char     *dsq;	/* digitized seq */
  float   *sc;			/* scores of seqs */
  int     *alen;		/* lengths of alignments */
  int      seed;		/* random number seed */
  int      N;			/* number of seqs to sample */
  int      L;			/* length of random sequences */
  float    randomseq[MAXABET];	/* iid frequencies of residues */
  float    p1;
  int      alphatype;		/* alphabet type, hmmAMINO or hmmNUCLEIC    */
  int      idx;
  int      code;
  Stopwatch_t stopwatch;       /* CPU timings */

  /* Register leave_pvm() cleanup function so any exit() call
   * first calls pvm_exit().
   */
  if (atexit(leave_pvm) != 0) { 
    pvm_exit(); Die("slave couldn't register leave_pvm()"); 
  }

  /*****************************************************************
   * initialization.
   * Master broadcasts the problem to us: 
   *    an HMM;
   *    parameters of the HMM calibration.
   * We send back:
   *    an OK flag, and our version number, for some sanity checking.  
   ******************************************************************/
  
  StopwatchStart(&stopwatch);

  master_tid = pvm_parent();	/* who's our master? */

  pvm_recv(master_tid, HMMPVM_INIT);
  pvm_upkint(&L,          1, 1); /* if non-zero, override lenmean */
  pvm_upkint(&alphatype,  1, 1); /* alphabet type, hmmAMINO or hmmNUCLEIC */
  pvm_upkint(&seed,       1, 1); /* random number seed */
  SetAlphabet(alphatype);	 /* must set alphabet before reading HMM! */
  hmm = PVMUnpackHMM();
  if (hmm == NULL) Die("oh no, the HMM never arrived");
  if (! (hmm->flags & PLAN7_HASBITS)) Die("Oops, that model isn't configured");
  P7DefaultNullModel(randomseq, &p1);
  mx = CreateDPMatrix(L, hmm->M, 0, 0);

  /* tell the master we're OK and ready to go (or not)
   */
  code = HMMPVM_OK;
  pvm_initsend(PvmDataDefault);
  pvm_pkint(&code, 1, 1);	
  PVMPackString(PACKAGE_VERSION);
  pvm_send(master_tid, HMMPVM_RESULTS);

  /*****************************************************************
   * Main loop.
   * Receive: a number of sequences we're supposed to do.
   *          If we receive a 0, we have no work, so wait for shutdown;
   *          if we receive a -1, shut down.
   *****************************************************************/ 
  slaveidx = -1;
  for (;;) 
    {
      pvm_recv(master_tid, HMMPVM_WORK);
      pvm_upkint(&nsample,  1, 1);
      pvm_upkint(&idx,      1, 1);

      if (nsample == 0)  continue;  /* go into stasis */
      if (nsample == -1) break;	    /* shut down      */

      if (slaveidx == -1) {	/* first time: set id, seed sre_random  */
	slaveidx = idx;
	sre_srandom(seed+idx);	/* unique seed in current PVM   */
      }

      sc   = MallocOrDie(sizeof(float) * nsample);
      alen = MallocOrDie(sizeof(int)   * nsample);
      for (idx = 0; idx < nsample; idx++)
	{
	  seq = RandomSequence(Alphabet, randomseq, Alphabet_size, L);
	  dsq = DigitizeSequence(seq, L);
	  SQD_DPRINTF2(("slave %d seq: %d : %20.20s...\n", slaveidx, L, seq));

	  if (P7ViterbiSpaceOK(L, hmm->M, mx))
	    sc[idx] = Viterbi(dsq, L, hmm, mx, &tr);
	  else
	    sc[idx] = P7SmallViterbi(dsq, L, hmm, mx, &tr);
	  TraceGetAlignmentBounds(tr, 1, NULL, NULL, NULL, NULL, &(alen[idx]));
	  
	  P7FreeTrace(tr);
	  free(seq);
	  free(dsq);
	}

      /* Return output to master, some of which is sanity checking.
       *   1. our slave index.
       *   2. how many seqs we simulated.
       *   3. the array of scores we got, so the master can stuff
       *      them into a histogram.
       *   4. array of alignment lengths, for calculating edge corrections.
       */
      pvm_initsend(PvmDataDefault);
      pvm_pkint(&slaveidx, 1, 1);
      pvm_pkint(&nsample,  1, 1);
      pvm_pkfloat(sc,   nsample, 1); 
      pvm_pkfloat(alen, nsample, 1); 
      pvm_send(master_tid, HMMPVM_RESULTS);

      /* cleanup
       */
      free(alen);
      free(sc);
    }

  /*********************************************** 
   * Cleanup, return.
   ***********************************************/
  
  FreePlan7(hmm);
  FreeDPMatrix(mx);
  StopwatchStop(&stopwatch);

  /* tell the master we heard his shutdown signal, and
   * give him our CPU times; then exit.
   */
  pvm_initsend(PvmDataDefault);
  pvm_pkint(&slaveidx, 1, 1);	
  StopwatchPVMPack(&stopwatch);
  pvm_send(master_tid, HMMPVM_RESULTS);  

  return 0;	/* pvm_exit() is called by atexit() registration. */
}

/* Function: leave_pvm()
 * 
 * Purpose:  Cleanup function, to deal with crashes. We register
 *           this function using atexit() so it gets called before
 *           the slave dies.
 */
void leave_pvm(void)
{
  SQD_DPRINTF1(("slave leaving PVM.\n"));
  pvm_exit();
}

#else /* if HMMER_PVM not defined: include a dummy */

#include <stdio.h>
int main(void)
{
  printf("hmmcalibrate-pvm disabled. PVM support was not compiled into HMMER.\n");
  exit(0);
} 

#endif

/************************************************************
 * @LICENSE@
 ************************************************************/

