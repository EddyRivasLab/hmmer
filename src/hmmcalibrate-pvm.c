/************************************************************
 * @LICENSE@
 ************************************************************/

#ifdef HMMER_PVM

/* hmmcalibrate-slave.c
 * SRE, Tue Aug 18 15:19:28 1998
 * 
 * PVM slave for hmmcalibrate.
 * RCS $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <pvm3.h>

#include "version.h"
#include "structs.h"		/* data structures, macros, #define's   */
#include "config.h"		/* compile-time configuration constants */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "squid.h"		/* general sequence analysis library    */

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static void leave_pvm(void);

int 
main(void)
{
  int      master_tid;		/* PVM TID of our master */
  int      slaveidx;		/* my slave index (0..nslaves-1) */
  struct plan7_s *hmm;		/* HMM to calibrate, sent from master */
  struct histogram_s *hist;     /* score histogram */
  int      hmmidx;		/* index of this HMM */
  char    *seq;			/* synthetic random sequence */
  char    *dsq;			/* digitized seq */
  int      len;			/* length of seq */
  float    sc;			/* score of seq aligned to HMM */
  float    max;			/* maximum score seen in sample */
  int      seed;		/* random number seed */
  int      nsample;		/* number of seqs to sample */
  int      fixedlen;		/* if nonzero, fixed length of seq */
  float    lenmean;		/* Gaussian mean length of seq */
  float    lensd;		/* Gaussian length std. dev. for seq */
  int      fitok;		/* TRUE if EVD fit was OK */
  float    randomseq[MAXABET];	/* iid frequencies of residues */
  float    p1;
  int      alphatype;		/* alphabet type, hmmAMINO or hmmNUCLEIC    */
  int      idx;
  int      code;

  /* Register leave_pvm() cleanup function so any exit() call
   * first calls pvm_exit().
   */
  if (atexit(leave_pvm) != 0) { pvm_exit(); Die("slave couldn't register leave_pvm()"); }

  /*****************************************************************
   * initialization.
   * Master broadcasts the problem to us: parameters of the
   * HMM calibration.  
   ******************************************************************/

  master_tid = pvm_parent();	/* who's our master? */

  pvm_recv(master_tid, HMMPVM_INIT);
  pvm_upkint(&nsample,  1, 1);
  pvm_upkint(&fixedlen, 1, 1);
  pvm_upkfloat(&lenmean,  1, 1);
  pvm_upkfloat(&lensd,    1, 1);

  /* tell the master we're OK and ready to go (or not)
   */
  code = HMMPVM_OK;
  pvm_initsend(PvmDataDefault);
  pvm_pkint(&code, 1, 1);	
  PVMPackString(RELEASE);
  pvm_send(master_tid, HMMPVM_RESULTS);

  /*****************************************************************
   * Main loop.
   * Receive a random number seed, then an HMM to search against.
   * If we receive a -1 seed, we shut down. 
   *****************************************************************/ 
  
  slaveidx = -1;
  for (;;) 
    {
      pvm_recv(master_tid, HMMPVM_WORK);
      pvm_upkint(&seed, 1, 1);
      if (seed == -1) break;	/* shutdown signal */
      pvm_upkint(&hmmidx, 1, 1);
      pvm_upkint(&alphatype,1, 1);
      SetAlphabet(alphatype);
      hmm = PVMUnpackHMM();
      if (hmm == NULL) Die("oh no, the HMM never arrived");

      if (slaveidx == -1) slaveidx = hmmidx; 
      P7DefaultNullModel(randomseq, &p1);

      sre_srandom(seed);
      P7Logoddsify(hmm, TRUE);
      hist = AllocHistogram(-200, 200, 100);
      max  = -FLT_MAX;

      for (idx = 0; idx < nsample; idx++)
	{
  				/* choose length of random sequence */
	  if (fixedlen) len = fixedlen;
	  else do len = (int) Gaussrandom(lenmean, lensd); while (len < 1);
				/* generate it */
	  seq = RandomSequence(Alphabet, randomseq, Alphabet_size, len);
	  dsq = DigitizeSequence(seq, len);

	  if (P7ViterbiSize(len, hmm->M) <= RAMLIMIT)
	    sc = P7Viterbi(dsq, len, hmm, NULL);
	  else
	    sc = P7SmallViterbi(dsq, len, hmm, NULL);

	  AddToHistogram(hist, sc);
	  if (sc > max) max = sc;
	  
	  free(seq);
	  free(dsq);
	}

      /* Fit an EVD to the observed histogram.
       * The TRUE left-censors and fits only the right slope of the histogram.
       * The 9999. is an arbitrary high number that means we won't trim outliers
       * on the right.
       */
      fitok = ExtremeValueFitHistogram(hist, TRUE, 9999.);

      /* Return output to master.
       * Currently we don't send the histogram back, but we could.
       */
      pvm_initsend(PvmDataDefault);
      pvm_pkint(&slaveidx, 1, 1);
      pvm_pkint(&hmmidx, 1, 1);	
      PVMPackString(hmm->name);
      pvm_pkint(&fitok,  1, 1);
      pvm_pkfloat(&(hist->param[EVD_MU]), 1, 1);
      pvm_pkfloat(&(hist->param[EVD_LAMBDA]), 1, 1);
      pvm_pkfloat(&max, 1, 1);
      pvm_send(master_tid, HMMPVM_RESULTS);

      /* cleanup
       */
      FreeHistogram(hist);
      FreePlan7(hmm);
    }

  /*********************************************** 
   * Cleanup, return.
   ***********************************************/

  return 0;			/* pvm_exit() is called by atexit() registration. */
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
  printf("hmmcalibrate-slave is disabled. PVM support was not compiled into HMMER.\n");
  exit(0);
} 

#endif
