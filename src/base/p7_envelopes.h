/* P7_ENVELOPES:  information about domain locations on a target sequence
 * 
 * Also defines P7_ENVELOPE, singular. P7_ENVELOPE contains
 * information about one domain, whereas P7_ENVELOPES is a reusable,
 * memory managed container for an array of envelopes. Sometimes we
 * use a naked array of <P7_ENVELOPE> structures without the wrapper,
 * so a caller can either be passing <env->arr> from <P7_ENVELOPES> or
 * its own allocation of a <P7_ENVELOPE *> array. We have an evil
 * tendency to call both the array and the container <env> in the code,
 * because the meaning is unambiguous if you're alert to the context;
 * so beware the distinction between the singular and plural forms.
 */
#ifndef p7ENVELOPES_INCLUDED
#define p7ENVELOPES_INCLUDED

#include "p7_config.h"

#include "easel.h"       // Needed for int32_t

/* P7_ENVELOPE
 *    Contains information about a domain's "envelope",
 *    as calculated from an anchor set. Doesn't need
 *    Create/Destroy; arrays of P7_ENVELOPE are just
 *    allocated normally. 
 *    
 *    Arrays of P7_ENVELOPE are always indexed 1..D, with sentinels at
 *    0 and D+1. Sentinel [0] has all coords as 0, flags 0, and env_sc
 *    0.0. Sentinel [D+1] has all i coords as L+1, all k coords as M+1,
 *    flags 0, and env_sc 0.0.
 */
typedef struct {
  int32_t  i0,   k0;	// anchor for this domain            (1..L, 1..M)
  int32_t  oea,  oeb;	// outer envelope start, stop on seq (1..L) 
  int32_t  ia,   ib;	// envelope start, end on sequence   (1..L) 
  int32_t  ka,   kb;	// ali start, end on model           (1..M) 
  float    env_sc;	// envelope score (nats)                    
  uint32_t flags;       // p7E_ENVSC_APPROX | p7E_IS_GLOCAL | p7E_IS_REPORTED | p7E_IS_INCLUDED
} P7_ENVELOPE;
  
#define p7E_ENVSC_APPROX (1<<0)
#define p7E_IS_GLOCAL    (1<<1)
#define p7E_IS_REPORTED  (1<<2)
#define p7E_IS_INCLUDED  (1<<3)
  
  

/* P7_ENVELOPES
 *    A memory-managed wrapper around a P7_ENVELOPE
 *    array, allowing arrays of envelopes to be 
 *    reused efficiently.
 *    
 *    arr[] is indexed 1..D, with arr[0] and arr[D+1]
 *    as sentinels; see note above on P7_ENVELOPES.
 */
typedef struct {
  P7_ENVELOPE *arr;		// array of envelope structures       
  int32_t      D;		// number of valid envelopes in <arr> 

  int32_t      L;		// length of sequence that envelopes are in 
  int32_t      M;		// length of model that envelopes are for   

  int32_t      nalloc;		// current allocation size for <arr>;  >= D+2 because of sentinels
  int32_t      nredline;	// _Reuse() pulls alloc back down to this 
} P7_ENVELOPES;

extern int           p7_envelope_SetSentinels(P7_ENVELOPE *env, int D, int L, int M);

extern P7_ENVELOPES *p7_envelopes_Create (void);
extern int           p7_envelopes_Reinit (P7_ENVELOPES *envs, int D);
extern int           p7_envelopes_Reuse  (P7_ENVELOPES *envs);
extern void          p7_envelopes_Destroy(P7_ENVELOPES *envs);

extern int p7_envelopes_Dump(FILE *ofp, P7_ENVELOPES *env);


#endif /*p7ENVELOPES_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
