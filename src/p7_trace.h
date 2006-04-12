/* p7_trace.h
 * The traceback structure, P7_TRACE:
 * alignment of a profile to a target sequence.
 * 
 * SVN $Id$
 * SRE, Tue Apr 11 16:32:50 2006
 */
#ifndef P7_TRACE_INCLUDED
#define P7_TRACE_INCLUDED

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

/* Plan 7 model state types used in traceback structure */
#define p7_STBOGUS 0
#define p7_STM     1
#define p7_STD     2
#define p7_STI     3
#define p7_STS     4
#define p7_STN     5
#define p7_STB     6
#define p7_STE     7
#define p7_STC     8
#define p7_STT     9
#define p7_STJ     10     


extern int  p7_trace_Create(int N, P7_TRACE **ret_tr);
extern int  p7_trace_Expand(P7_TRACE *tr);
extern int  p7_trace_ExpandTo(P7_TRACE *tr, int N);
extern void p7_trace_Destroy(P7_TRACE *tr);
extern void p7_trace_Dump(FILE *fp, P7_TRACE *tr, P7_PROFILE *gm, char *dsq);
extern int  p7_trace_Append(P7_TRACE *tr, char st, int k, int i);
extern int  p7_trace_Reverse(struct p7trace_s *tr);

#endif /*P7_TRACE_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
