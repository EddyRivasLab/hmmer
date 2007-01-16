/* Traceback structure
 * 
 * SRE, Tue Jan  2 2007 [Casa de Gatos] 
 * SVN $Id$
 */
#ifndef P7_TRACEH_INCLUDED
#define P7_TRACEH_INCLUDED

#include "p7_config.h"
#include "p7_hmm.h"

#include "easel.h"

/* P7_TRACE.
 * Traceback structure for alignments of model to sequence.
 *
 * A traceback only makes sense in a triplet (tr, hmm, ax|dsq),
 * for a given HMM (with nodes 1..M) and a given digital sequence 
 * (with positions 1..L).
 * 
 * A traceback is always relative to the search form of the model,
 * so they always include the S,N,B... E,C,T states. Element 0 always 
 * a p7_STS. Element N-1 always a p7_STT.
 *
 * A traceback may be relative to an aligned sequence or an unaligned
 * sequence in digital mode.
 */
typedef struct {
  int   N;		/* length of traceback                       */
  int   nalloc;		/* allocated length of traceback             */
  char *st;		/* state type code                   [0..N-1]*/
  int  *k;		/* node index; 1..M if M,D,I; else 0 [0..N-1]*/
  int  *i;		/* position in dsq, 1..L; 0 if none  [0..N-1]*/
} P7_TRACE;

extern int  p7_trace_Create(int N, P7_TRACE **ret_tr);
extern int  p7_trace_Reuse(P7_TRACE *tr);
extern int  p7_trace_Expand(P7_TRACE *tr);
extern int  p7_trace_ExpandTo(P7_TRACE *tr, int N);
extern void p7_trace_Destroy(P7_TRACE *tr);
extern void p7_trace_DestroyArray(P7_TRACE **tr, int N);
extern int  p7_trace_Validate(P7_TRACE *tr, ESL_ALPHABET *abc, ESL_DSQ *sq);
extern int  p7_trace_Dump(FILE *fp, P7_TRACE *tr, void *gm, ESL_DSQ *dsq);

extern int  p7_trace_Append(P7_TRACE *tr, char st, int k, int i);
extern int  p7_trace_Reverse(P7_TRACE *tr);
extern int  p7_trace_Count(P7_HMM *hmm, ESL_DSQ *dsq, float wt, P7_TRACE *tr);

#endif /*P7_TRACEH_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
