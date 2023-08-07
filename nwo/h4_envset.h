#ifndef h4ENVSET_INCLUDED
#define h4ENVSET_INCLUDED

#include <h4_config.h>

#include "h4_anchorset.h"

/* H4_ENVELOPE
 *    Contains information about a domain's "envelope".
 */
typedef struct {
  int   i0,   k0;	// anchor for this domain            (1..L, 1..M)
  int   ia,   ib;	// envelope start, end on sequence   (1..L) 
  int   ka,   kb;	// ali start, end on model           (1..M) 
  int   oa,   ob;	// outer envelope start, stop on seq (1..L) 
  float env_sc;	        // envelope raw score s_d, nats                    
  float null2_sc;       // domain null2 score r_d, nats
  
  uint32_t flags;       // h4E_ENVSC_APPROX | h4E_IS_GLOCAL 
} H4_ENVELOPE;
  
#define h4E_ENVSC_APPROX (1<<0)
#define h4E_IS_GLOCAL    (1<<1)


/* H4_ENVSET
 *    A memory-managed wrapper around an array of H4_ENVELOPE structures,
 *    one per domain defined by the anchor set.
 *    
 *    Envelopes are indexed 1..D, with sentinels at 0 and D+1.
 *
 *    Sentinels are set the same as in an H4_ANCHORSET:
 *         e[0]   (i',k') = (0, M+1)  
 *         e[D+1] (i',k') = (L+1, 0)
 *    for all coords i' = {oa,ia,i0,ib,ob} on the sequence and
 *    k' = {ka,k0,kb} on the profile.
 *
 *    AEC traceback depends in the ia[D+1] sentinel.
 */
typedef struct {
  H4_ENVELOPE *e;	// array of envelope structures       
  int          D;	// number of valid envelopes in <e> 

  int  L;		// length of sequence that envelopes are in 
  int  M;		// length of model that envelopes are for   

  int  nalloc;		// current allocation size for <e>;  >= D+2 because of sentinels
  int  nredline;	// _Reuse() pulls alloc back down to this 
} H4_ENVSET;


extern H4_ENVSET *h4_envset_Create(int D, int L, int M);
extern int        h4_envset_Resize(H4_ENVSET *env, int D);
extern int        h4_envset_CopyFromAnchorset(const H4_ANCHORSET *anch, H4_ENVSET *env);
extern int        h4_envset_Dump(FILE *ofp, const H4_ENVSET *env);
extern void       h4_envset_Destroy(H4_ENVSET *env);

#endif
