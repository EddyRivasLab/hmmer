/* P7_ANCHORS: arrays of i0,k0 anchor points that define domain locations.
 * 
 * Also defines P7_ANCHOR, singular. P7_ANCHOR is a single i0,k0
 * anchor point, whereas P7_ANCHORS is a reusable, memory managed
 * container for an array of anchors. Sometimes we use a naked array
 * of <P7_ANCHOR> structures without the wrapper, so that a caller can
 * either be passing <anch->arr> from <P7_ANCHORS>, or its own
 * allocation of a <P7_ANCHOR *> array. We have a nasty tendency to
 * call both the array and the container <anch> in the code, so beware
 * the distinction between the singular and plural.
 */
#ifndef p7ANCHORS_INCLUDED
#define p7ANCHORS_INCLUDED

#include "p7_config.h"

#include "base/p7_trace.h"

#include "esl_random.h"

/* P7_ANCHOR
 *    Is just an integer pair, usually used in an array.
 *    Doesn't need Create/Destroy; arrays of these are just
 *    allocated normally, or as part of P7_ANCHORS.
 */
typedef struct {
  int32_t i0;
  int32_t k0;
} P7_ANCHOR;


/* P7_ANCHORS
 *    Is a memory-managed wrapper around a P7_ANCHOR array,
 *    allowing arrays of start/end positions (as in domain
 *    definition code) to be rapidly created and reused.
 *
 *    Anchors are indexed 1..D, with a[0] and a[D+1] used as
 *    sentinels:
 *      a[0].i0,k0   = (0,0)
 *      a[D+1].i0,k0 = (L+1,M+1)
 *
 *    Therefore, for D domains, anch->a[] is always
 *    allocated for at least D+2 <P7_ANCHOR> structures,
 *    (0) 1..D (D+1).
 */
typedef struct {
  P7_ANCHOR *a;	   	 // array of anchor (i0,k0) coords: (0) 1..D (D+1)   
  int32_t    D;		 // number of domains in <a>, exclusive of sentinels 

  int32_t    nalloc;	 // current allocation size for <a>; >= D+2 because of sentinels
  int32_t    nredline;	 // a Reuse() pulls alloc back down to this                     
} P7_ANCHORS;


/* 1. P7_ANCHOR arrays, either naked or as part of P7_ANCHORS */
extern int p7_anchor_GetSentinels(P7_ANCHOR *anch, int D, int *ret_L, int *ret_M);
extern int p7_anchor_SetSentinels(P7_ANCHOR *anch, int D, int L,      int M);

/* 2. Standard object stuff for P7_ANCHORS */
extern P7_ANCHORS *p7_anchors_Create  (void);
extern int         p7_anchors_Resize  (P7_ANCHORS *anch, int D);
extern int         p7_anchors_Copy    (const P7_ANCHORS *src, P7_ANCHORS *dst);
extern int         p7_anchors_Catenate(P7_ANCHORS *anch, int D0, P7_ANCHOR *arr, int Dg);
extern int         p7_anchors_Reuse   (P7_ANCHORS *anch);
extern void        p7_anchors_Destroy (P7_ANCHORS *anch);

/* 3. Debugging and development tools */
extern int  p7_anchors_Dump(FILE *fp, const P7_ANCHORS *anch);
extern int  p7_anchors_DumpOneLine(FILE *ofp, const P7_ANCHORS *anch);
extern int  p7_anchors_Compare(P7_ANCHORS *anch1, P7_ANCHORS *anch2);
extern int  p7_anchors_Validate(P7_ANCHORS *anch, int L, int M, char *errbuf);
extern int  p7_anchors_Sample(ESL_RANDOMNESS *rng, int L, int M, int maxD, P7_ANCHORS *anch);
extern int  p7_anchors_SampleFromTrace(P7_ANCHORS *anch, ESL_RANDOMNESS *rng, const P7_TRACE *tr);

#endif /* p7ANCHORS_INCLUDED */

  
  


