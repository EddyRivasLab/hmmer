/************************************************************
 * @LICENSE@
 ************************************************************/

#ifdef HMMER_PVM

/* hmmslave-pvm.c
 * SRE, Sun Jul 12 17:15:36 1998
 * 
 * PVM slave for hmmpfam-pvm and hmmsearch-pvm.
 * RCS $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <pvm3.h>

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
  struct p7trace_s *tr;         /* traceback of an alignment               */
  int      master_tid;		/* PVM TID of our master */
  char    *hmmfile;	        /* file to read HMM(s) from                */
  HMMFILE *hmmfp;               /* opened hmmfile for reading              */
  struct plan7_s *hmm;
  char    *seq;
  char    *dsq;
  int      len;
  int      nhmm;		/* number of HMM to work on                */
  float    sc;
  int      my_idx = -1;		/* my index, 0..nslaves-1 */
  float    globT;		/* T parameter: keep only hits > globT bits */
  double   globE;		/* E parameter: keep hits < globE E-value   */
  double   pvalue;		/* Z*pvalue = Evalue                        */
  int      Z;			/* nseq to base E value calculation on      */
  int      send_trace;		/* TRUE if score is significant             */
  int      do_xnu;		/* TRUE to do XNU filter on seq             */
  int      do_forward;		/* TRUE to use Forward() scores not Viterbi */
  int      do_null2;		/* TRUE to correct scores w/ ad hoc null2   */
  int      alphatype;		/* alphabet type, hmmAMINO or hmmNUCLEIC    */
  int      code;		/* return code after initialization         */

  
  /* Register leave_pvm() cleanup function so any exit() call
   * first calls pvm_exit().
   */
  if (atexit(leave_pvm) != 0) { pvm_exit(); Die("slave couldn't register leave_pvm()"); }

  /*****************************************************************
   * initialization.
   * Master broadcasts to us: 
   *     1) len of HMM file name        (int)
   *     2) name of HMM file            (string)
   *     3) length of sequence string   (int) 
   *     4) sequence                    (string)
   *     5) globT threshold
   *     6) globE threshold
   *     7) Z 
   *     8) do_xnu flag
   *     9) do_forward flag
   *    10) do_null2 flag
   *    11) alphabet type
   * We receive the broadcast and open the files.    
   ******************************************************************/

  master_tid = pvm_parent();	/* who's our master? */

  pvm_recv(master_tid, HMMPVM_INIT);
  pvm_upkint(&len, 1, 1);
  hmmfile = MallocOrDie(sizeof(char *) * (len+1));
  pvm_upkstr(hmmfile);
  pvm_upkint(&len, 1, 1);
  seq = MallocOrDie(sizeof(char *) * (len+1));
  pvm_upkstr(seq);
  pvm_upkfloat(&globT, 1, 1);
  pvm_upkdouble(&globE, 1, 1);
  pvm_upkint(&Z, 1, 1);
  pvm_upkint(&do_xnu, 1, 1);
  pvm_upkint(&do_forward, 1, 1);
  pvm_upkint(&do_null2, 1, 1);
  pvm_upkint(&alphatype, 1, 1);

  SetAlphabet(alphatype);
				/* Open HMM file (maybe in HMMERDB) */
  code = HMMPVM_OK;
  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    code = HMMPVM_NO_HMMFILE;
  else if (hmmfp->gsi == NULL)
    code = HMMPVM_NO_INDEX;
  
  /* report our status.
   */
  pvm_initsend(PvmDataDefault);
  pvm_pkint(&code, 1, 1);	
  pvm_send(master_tid, HMMPVM_RESULTS);

  dsq = DigitizeSequence(seq, len);
  if (do_xnu) XNU(dsq, len);

  /*****************************************************************
   * Main loop.
   * Receive an integer 0..nhmm-1 for which HMM to search against.
   * If we receive a -1, we shut down. 
   *****************************************************************/ 
  
  for (;;) 
    {
      pvm_recv(master_tid, HMMPVM_WORK);
      pvm_upkint(&nhmm, 1, 1);
      if (my_idx < 0) my_idx = nhmm; /* first time thru, remember what index we are. */

      if (nhmm == -1) break;	/* shutdown signal */

      /* move to our assigned HMM in the HMM file, and read it
       */
      HMMFilePositionByIndex(hmmfp, nhmm);
      if (! HMMFileRead(hmmfp, &hmm)) Die("unexpected end of HMM file"); 
      if (hmm == NULL)                Die("unexpected failure to parse HMM file"); 
      P7Logoddsify(hmm, TRUE);
      
      /* Score sequence, do alignment (Viterbi), recover trace
       */
      if (P7ViterbiSize(len, hmm->M) <= RAMLIMIT)
	{
	  SQD_DPRINTF1(("P7Viterbi(): Estimated size %d Mb\n", P7ViterbiSize(len, hmm->M)));
	  sc = P7Viterbi(dsq, len, hmm, &tr);
	}
      else
	{
	  SQD_DPRINTF1(("P7SmallViterbi() called; %d Mb > %d\n", P7ViterbiSize(len, hmm->M), RAMLIMIT));
	  sc = P7SmallViterbi(dsq, len, hmm, &tr);
	}

      if (do_forward) sc  = P7Forward(dsq, len, hmm, NULL);
      if (do_null2)   sc -= TraceScoreCorrection(hmm, tr, dsq);
	
      pvalue = PValue(hmm, sc);
      send_trace = (sc > globT && pvalue * (float) Z < globE) ? 1 : 0;

      /* return output
       */
      pvm_initsend(PvmDataDefault);
      pvm_pkint(&my_idx, 1, 1);	/* tell master who we are */
      pvm_pkstr(hmm->name);	/* double check that we did the right thing */
      pvm_pkfloat(&sc, 1, 1);
      pvm_pkdouble(&pvalue, 1, 1);
      pvm_pkint(&send_trace, 1, 1); /* flag for whether a trace structure is coming */
      if (send_trace) PVMPackTrace(tr);
      pvm_send(master_tid, HMMPVM_RESULTS);

      /* cleanup
       */
      FreePlan7(hmm);
      P7FreeTrace(tr);
    }

  /*********************************************** 
   * Cleanup, return.
   ***********************************************/

  HMMFileClose(hmmfp);
  free(seq);
  free(dsq);
  free(hmmfile);
  return 0;
}


/* Function: leave_pvm()
 * 
 * Purpose:  Cleanup function, to deal with crashes. We register
 *           this function using atexit() so it gets called before
 *           the slave dies.
 */
static void leave_pvm(void)
{
  SQD_DPRINTF1(("slave leaving PVM.\n"));
  pvm_exit();
}



#else /* if HMMER_PVM not defined: include a dummy */

#include <stdio.h>
int main(void)
{
  printf("hmmpfam-slave is disabled. PVM support was not compiled into HMMER.\n");
  exit(0);
} 

#endif

