/* p7_trace.h
 * The traceback structure, P7_TRACE:
 * alignment of a profile to a target sequence.
 * 
 * SVN $Id$
 * SRE, Tue Apr 11 16:32:50 2006
 */
#ifndef P7_TRACE_INCLUDED
#define P7_TRACE_INCLUDED

#include <stdio.h>		
#include "p7_hmm.h"		
#include "p7_profile.h"		

/* P7_TRACE.
 * Traceback structure for alignments of model to sequence.
 * Element 0 always p7_STS. Element N-1 always p7_STT.
 */
typedef struct {
  int   N;		/* length of traceback                       */
  int   nalloc;		/* allocated length of traceback             */
  char *st;		/* state type code                   [0..N-1]*/
  int  *k;		/* node index; 1..M if M,D,I; else 0 [0..N-1]*/
  int  *i;		/* position in dsq, 1..L; 0 if none  [0..N-1]*/
} P7_TRACE;

extern int  p7_trace_Create(int N, P7_TRACE **ret_tr);
extern int  p7_trace_Expand(P7_TRACE *tr);
extern int  p7_trace_ExpandTo(P7_TRACE *tr, int N);
extern void p7_trace_Destroy(P7_TRACE *tr);
extern int  p7_trace_Dump(FILE *fp, P7_TRACE *tr, P7_PROFILE *gm, char *dsq);
extern int  p7_trace_Append(P7_TRACE *tr, char st, int k, int i);
extern int  p7_trace_Reverse(P7_TRACE *tr);
extern int  p7_trace_Count(P7_HMM *hmm, char *dsq, float wt, P7_TRACE *tr, int mode);
extern int  p7_trace_Score(P7_PROFILE *gm, char *dsq, P7_TRACE *tr, int *ret_sc);
extern int  p7_trace_DomainCount(P7_TRACE *tr);
extern int  p7_trace_Decompose(P7_TRACE *otr, P7_TRACE ***ret_tr, int *ret_ntr);
extern int  p7_trace_GetDomainCoords(P7_TRACE *tr, int which,
				     int *ret_i1, int *ret_i2,
				     int *ret_k1, int *ret_k2,
				     int *ret_avlen);

#endif /*P7_TRACE_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
