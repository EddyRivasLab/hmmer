/* hmmsearch-pvm.c
 * PVM slave for hmmsearch.
 * 
 * SRE, Wed Sep 23 09:30:53 1998
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

#include "plan7.h"		/* plan 7 profile HMM structure         */
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */

static void leave_pvm(void);

int 
main(void)
{
  struct plan7_s *hmm;          /* HMM to search with */
  struct dpmatrix_s *mx;        /* growable DP matrix */
  struct p7trace_s *tr;         /* trace structure for a Viterbi alignment */
  int master_tid;		/* PVM TID of our master */
  int alphatype;		/* alphabet type */
  int code;			/* status code for whether we're ok */
  int my_idx;			/* my slave index: 0..nslaves-1, master assigns */
  int L;			/* length of sequence */
  unsigned char *dsq;           /* digitized sequence 1..L */
  float  sc;			/* log odds score for seq + HMM */
  double pvalue;		/* P-value of sc */
  double evalue;		/* bounded E-value of sc (we don't know nseq yet) */
  int    do_forward;		/* TRUE to score using Forward() */
  int    do_null2;		/* TRUE to use null2 ad hoc correction */
  float  globT;			/* T parameter: keep only hits > globT bits */
  double globE;			/* E parameter: keep hits < globE E-value   */
  int    Z;			/* nseq to base E value calculation on      */
  int    nseq;			/* actual nseq so far (master keeps updating this) */
  int    send_trace;		/* TRUE if sc looks significant and we return tr */
  
  tr = NULL;
  
  /* Register leave_pvm() cleanup function so any exit() call
   * first calls pvm_exit().
   */
  if (atexit(leave_pvm) != 0) { pvm_exit(); Die("slave couldn't register leave_pvm()"); }

  /*****************************************************************
   * Initialization.
   * Master broadcasts the problem to us:
   * globT, globE, Z, do_forward, do_null2, alphabet type, HMM, 
   ******************************************************************/

  master_tid = pvm_parent();	/* who's our master? */
  my_idx     = -1;

  /* wait for a HMMPVM_INIT message, and unpack it;
   * get options, set alphabet type, get HMM.
   */
  pvm_recv(master_tid, HMMPVM_INIT);
  pvm_upkfloat(&globT, 1, 1);
  pvm_upkdouble(&globE, 1, 1);
  pvm_upkint(&Z,          1, 1);
  pvm_upkint(&do_forward, 1, 1);
  pvm_upkint(&do_null2,   1, 1);
  pvm_upkint(&alphatype,  1, 1);
  SetAlphabet(alphatype);
  hmm = PVMUnpackHMM();

  mx = CreatePlan7Matrix(1, hmm->M, 25, 0);
  if (! (hmm->flags & PLAN7_HASBITS)) Die("no scores in the model");

  /* tell the master we're OK and ready to go (or not)
   */
  code = HMMPVM_OK;
  if (hmm == NULL) code = HMMPVM_BAD_INIT;
  pvm_initsend(PvmDataDefault);
  pvm_pkint(&code, 1, 1);	
  PVMPackString(PACKAGE_VERSION);
  pvm_send(master_tid, HMMPVM_RESULTS);

  /*****************************************************************
   * Main loop.
   * Receive a digitized sequence to search against.
   *****************************************************************/ 
 

  for (;;)
    {
      SQD_DPRINTF1(("Slave about to do a blocking receive, waiting for input.\n"));
      pvm_recv(master_tid, HMMPVM_WORK);
      pvm_upkint(&nseq, 1, 1);
      if (nseq == -1) break;	/* shutdown signal */
      if (my_idx == -1) my_idx = nseq;
      pvm_upkint(&L,    1, 1);
      SQD_DPRINTF1(("Slave received nseq=%d L=%d my_idx=%d\n", nseq, L, my_idx));
      dsq = MallocOrDie(sizeof(unsigned char) * (L + 2));
      pvm_upkbyte((char *) dsq, L+2, 1); /* char * cast is dubious practice */
      SQD_DPRINTF1(("Slave unpacked a seq of %d bytes; beginning processing\n", L+2));
      
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
          if (P7ViterbiSpaceOK(sqinfo.len, hmm->M, mx))
          {
              sc = P7Viterbi(dsq, L, hmm, mx, &tr);
          }
          else
          {
              /* Low-memory version */
              sc = P7SmallViterbi(dsq, L, hmm, mx, &tr);
          }
      }
      else
      {
          /* When using Altivec, the new dynamic programming process is so 
          * effective that it is constrained by the memory bandwidth used
          * for loading/storing values in the dynamic programming matrix.
          *
          * In 99% of all cases the score is anyway so low that we are not
          * interested in the alignment, and thus don't NEED the trace.
          * which is the only reason to store the full programming matrix.
          *
          * Do the programming without a trace first, and recalculate the
          * trace further down if it turns out we need it.
          */          
          sc = P7ViterbiNoTrace(dsq, L, hmm, mx);
          tr = NULL;
      }
      
#else /* not altivec */
      
      if (P7ViterbiSpaceOK(L, hmm->M, mx))
      {
          SQD_DPRINTF1(("Slave doing Viterbi after estimating %d MB\n", (P7ViterbiSize(L, hmm->M))));
          sc = P7Viterbi(dsq, L, hmm, mx, &tr);
      }
      else
      {
          SQD_DPRINTF1(("Slave going small after estimating %d MB\n", (P7ViterbiSize(L, hmm->M))));
          sc = P7SmallViterbi(dsq, L, hmm, mx, &tr);
      }

#endif
      
      if (do_forward) 
      {
          sc  = P7Forward(dsq, L, hmm, NULL);
          if (do_null2) 
          {
              /* We need the trace - recalculate it if we didn't already do it */
#ifdef ALTIVEC
              if(tr == NULL)
              {
                  if (P7ViterbiSpaceOK(L, hmm->M, mx))
                  {
                      P7Viterbi(dsq, L, hmm, mx, &tr);
                  }
                  else
                  {
                      P7SmallViterbi(dsq, L, hmm, mx, &tr);
                  }
              }
#endif
              sc -= TraceScoreCorrection(hmm, tr, dsq);
          }
      }
	
      pvalue = PValue(hmm, sc);
      evalue = Z ? (double) Z * pvalue : (double) nseq * pvalue;
      send_trace = (sc >= globT && evalue <= globE) ? 1 : 0;
     
      /* return output
       */
      SQD_DPRINTF1(("Slave has a result (sc = %.1f); sending back to master\n", sc));
      pvm_initsend(PvmDataDefault);
      pvm_pkint   (&my_idx,  1, 1);   
      pvm_pkfloat (&sc,      1, 1);
      pvm_pkdouble(&pvalue,  1, 1);
      pvm_pkint(&send_trace, 1, 1); /* flag for whether a trace structure is coming */
      if (send_trace)
      {
          /* We need the trace - recalculate it if we didn't already do it */
#ifdef ALTIVEC
          if(tr == NULL)
          {
              if (P7ViterbiSpaceOK(L, hmm->M, mx))
              {
                  P7Viterbi(dsq, L, hmm, mx, &tr);
              }
              else
              {
                  P7SmallViterbi(dsq, L, hmm, mx, &tr);
              }              
          }
#endif
          PVMPackTrace(tr);
      }
      pvm_send(master_tid, HMMPVM_RESULTS);
      
      /* cleanup
          */
      free(dsq);
      
      if(tr != NULL)
          P7FreeTrace(tr);
      
      tr = NULL;
    }

  /*********************************************** 
   * Cleanup, return.
   ***********************************************/

  SQD_DPRINTF1(("Slave is done; performing a normal exit.\n"));
  FreePlan7Matrix(mx);
  FreePlan7(hmm);
  exit(0);			/* pvm_exit() gets called by atexit() registration. */
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
  printf("hmmsearch-pvm is disabled. PVM support was not compiled into HMMER.\n");
  exit(0);
} 

#endif

/************************************************************
 * @LICENSE@
 ************************************************************/

