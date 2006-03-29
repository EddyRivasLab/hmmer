/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* emit.c
 * SRE, Sun Mar  8 12:26:58 1998
 * RCS $Id$
 * 
 * Generation of sequences/traces from an HMM.
 */

#include "structs.h"
#include "config.h"
#include "funcs.h"
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: EmitSequence()
 * Date:     SRE, Sun Mar  8 12:28:03 1998 [St. Louis]
 *
 * Purpose:  Given a model, sample a sequence and/or traceback.
 *
 * Args:     hmm     - the model
 *           ret_dsq - RETURN: generated digitized sequence (pass NULL if unwanted)
 *           ret_L   - RETURN: length of generated sequence 
 *           ret_tr  - RETURN: generated trace (pass NULL if unwanted)
 *
 * Returns:  void
 */
void
EmitSequence(struct plan7_s *hmm, char **ret_dsq, int *ret_L, struct p7trace_s **ret_tr)
{
  struct p7trace_s *tr;
  enum   p7stype    type;	/* current state type */
  int   k;			/* current node index */
  char *dsq;                    /* generated sequence, digitized */
  int   L;			/* length of sequence */
  int   alloc_tlen;		/* allocated space for traceback */
  int   alloc_L;		/* allocated space for sequence  */
  int   tpos;			/* position in traceback */
  int   sym;			/* a generated symbol index */
  float t[4];			/* little array for choosing M transition from */
  
  /* Initialize; allocations
   */
  P7AllocTrace(64, &tr);
  alloc_tlen = 64;
  dsq = MallocOrDie(sizeof(char) * 64);
  alloc_L = 64;

  TraceSet(tr, 0, STS, 0, 0);
  TraceSet(tr, 1, STN, 0, 0);
  dsq[0] = (char) Alphabet_iupac;
  L      = 1;
  k      = 0;
  type   = STN;
  tpos   = 2;

  while (type != STT) 
    {
      /* Deal with state transition
       */
      switch (type) {
      case STB:	type = STM; k = FChoose(hmm->begin+1, hmm->M) + 1; break;
      case STI:	type = (FChoose(hmm->t[k]+TIM, 2) == 0)    ? STM : STI; if (type == STM) k++; break;
      case STN: type = (FChoose(hmm->xt[XTN], 2)  == LOOP) ? STN : STB; k = 0; break;
      case STE:	type = (FChoose(hmm->xt[XTE], 2)  == LOOP) ? STJ : STC; k = 0; break;
      case STC:	type = (FChoose(hmm->xt[XTC], 2)  == LOOP) ? STC : STT; k = 0; break;
      case STJ:	type = (FChoose(hmm->xt[XTJ], 2)  == LOOP) ? STJ : STB; k = 0; break;

      case STD:	
	if (k < hmm->M) {
	  type = (FChoose(hmm->t[k]+TDM, 2) == 0) ? STM : STD; 
	  k++;   
	} else {
	  type = STE;
	  k = 0;
	}
	break;

      case STM:
	if (k < hmm->M) {
	  FCopy(t, hmm->t[k], 3);
	  t[3] = hmm->end[k];
	  switch (FChoose(t,4)) {
	  case 0: k++;  type = STM; break;
	  case 1:       type = STI; break;
	  case 2: k++;  type = STD; break;
	  case 3: k=0;  type = STE; break;
	  default: Die("never happens");
	  }
	} else {
	  k    = 0;
	  type = STE;
	}
	break;

      case STT:
      case STBOGUS:
      default:
	Die("can't happen.");
      }
  
      /* Choose a symbol emission, if necessary
       */
      sym = -1;
      if      (type == STM) sym = FChoose(hmm->mat[k], Alphabet_size);
      else if (type == STI) sym = FChoose(hmm->ins[k], Alphabet_size); 
      else if ((type == STN && tr->statetype[tpos-1] == STN) ||
	       (type == STC && tr->statetype[tpos-1] == STC) ||
	       (type == STJ && tr->statetype[tpos-1] == STJ))
	sym = FChoose(hmm->null, Alphabet_size);
	
      /* Add to the traceback; deal with realloc if necessary
       */
      TraceSet(tr, tpos, type, k, (sym != -1) ? L : 0);
      tpos++;
      if (tpos == alloc_tlen) {
	alloc_tlen += 64; 
	P7ReallocTrace(tr, alloc_tlen);
      }

      /* Add to the digitized seq; deal with realloc, if necessary
       */
      if (sym != -1) {
	dsq[L] = (char) sym;
	L++;
	if (L+1 == alloc_L) {	/* L+1 leaves room for sentinel byte + \0 */
	  alloc_L += 64;
	  dsq = ReallocOrDie(dsq, sizeof(char) * alloc_L);
	}
      }
    }
  
  /* Finish off the trace
   */ 
  tr->tlen = tpos;

  /* Finish off the dsq with sentinel byte and null terminator.
   * Emitted Sequence length is L-1.
   */
  dsq[L]   = (char) Alphabet_iupac;
  dsq[L+1] = '\0';
  L--;

  /* Return
   */
  if (ret_dsq != NULL) *ret_dsq = dsq; else free(dsq);
  if (ret_L   != NULL) *ret_L   = L;
  if (ret_tr  != NULL) *ret_tr  = tr;  else P7FreeTrace(tr);
  return;
}

