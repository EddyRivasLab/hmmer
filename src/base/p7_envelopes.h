#ifndef P7_ENVELOPES_INCLUDED
#define P7_ENVELOPES_INCLUDED

#include "p7_config.h"


/* P7_ENVELOPE
 *    Contains information about a domain's "envelope",
 *    as calculated from an anchor set. Doesn't need
 *    Create/Destroy; arrays of P7_ENVELOPE are just
 *    allocated normally. 
 * 
 *    48 bytes per domain.
 */
typedef struct {
  int32_t  i0, k0;	/* anchor for this domain                   */
  int32_t  oea, oeb;	/* outer envelope start, stop on seq (1..L) */
  int32_t  ia, ib;	/* envelope start, end on sequence   (1..L) */
  int32_t  alia, alib;	/* ali start, end on sequence        (1..L) */
  int32_t  ka, kb;	/* ali start, end on model           (1..M) */
  float    env_sc;	/* envelope score (nats)                    */
  uint32_t flags;
} P7_ENVELOPE;
  
#define p7E_ENVSC_APPROX (1<<0)
#define p7E_IS_GLOCAL    (1<<1)
#define p7E_IS_REPORTED  (1<<2)
#define p7E_IS_INCLUDED  (1<<3)
  
  

/* P7_ENVELOPES
 *    A memory-managed wrapper around a P7_ENVELOPE
 *    array, allowing arrays of envelopes to be 
 *    reused efficiently.
 */
typedef struct {
  P7_ENVELOPE *arr;		/* array of envelope structures       */
  int32_t      n;		/* number of valid envelopes in <arr> */

  int32_t      L;		/* length of sequence that envelopes are in */
  int32_t      M;		/* length of model that envelopes are for   */

  int32_t      nalloc;		/* current allocation size for <arr> */
  int32_t      nredline;	/* _Reuse() pulls alloc back down to this */
} P7_ENVELOPES;



extern P7_ENVELOPES *p7_envelopes_Create (int32_t nalloc, int32_t nredline);
extern int           p7_envelopes_GrowTo (P7_ENVELOPES *envs, int32_t nalloc);
extern int           p7_envelopes_Reuse  (P7_ENVELOPES *envs);
extern void          p7_envelopes_Destroy(P7_ENVELOPES *envs);

extern int p7_envelopes_Dump(FILE *ofp, P7_ENVELOPES *env);


#endif /* P7_ENVELOPES_INCLUDED */
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
