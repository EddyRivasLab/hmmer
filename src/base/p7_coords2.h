#ifndef P7_COORDS2_INCLUDED
#define P7_COORDS2_INCLUDED

#include "p7_config.h"

#include "base/p7_trace.h"

#include "esl_random.h"

/* P7_COORD2 
 *    Is just an integer pair, usually used in an array.
 *    Doesn't need Create/Destroy; arrays of these are just
 *    allocated normally.
 */
typedef struct {
  int32_t n1;
  int32_t n2;
} P7_COORD2;


/* P7_COORDS2
 *    Is a memory-managed wrapper around a P7_COORD2 array,
 *    allowing arrays of start/end positions (as in domain
 *    definition code) to be rapidly created and reused.
 *
 *    Watch out for the distinction between P7_COORD2, P7_COORDS2;
 *    perhaps we should make more distinctive names.
 */
typedef struct {
  P7_COORD2 *arr;		/* array of coord pairs                  */
  int32_t    n;			/* number of coord pairs in <arr>        */

  int32_t    dim1;		/* max value of coordinate 1, arr[].n1   */
  int32_t    dim2;		/* max value of coordinate 2, arr[].n2   */

  int32_t    nalloc;		/* current allocation size for <arr>     */
  int32_t    nredline;		/* Reuse() pulls alloc back down to this */
} P7_COORDS2;


/* P7_COORD2_HASH
 *    Experimental, for now.
 *    A hash table for storing P7_COORD2 arrays.
 *    To support ongoing experiments with MPL domain definition.
 *    Patterned on ESL_KEYHASH implementation.
 *    Pointerless, to facilitate copying.
 *
 * Example: 
 *   suppose we've stored 2 different domain definitions for seq of length 100.
 *   1. Two domains, 42..61 72..90
 *   2. Three domains, 23..37 43..62 71..89
 * then:
 *   cmem[]  = 2 42 61 72 90 3 23 37 43 62 71 89
 *   cn      = 12
 *   calloc >= 12
 *   
 *   key_offset[] = 0 5
 *   nkeys        = 2
 *   kalloc       >= 2
 *   
 *   hashtable    = { -1 ... 0 ... -1 ... 1 ... -1 }
 *   nxt[]        = -1 -1
 *   
 * Note: structure is suitable for becoming more general, a hash
 *  of any integer array.
 *  
 */
typedef struct {
  int32_t  *hashtable;	  /* [0..hashsize-1]; head node for chains; index of 1st elem, 0..nkeys-1  */
  int32_t   hashsize;	  /* size of the hashtable (# of buckets). Must be a power of 2 */

  int32_t  *key_offset;	  /* [k=0..nkeys-1]; key[k] data starts at cmem + key_offset[k] */
  int32_t  *nxt;  	  /* [k=0..nkeys-1]; nxt[k] = next elem in chain, or -1 for end */
  int32_t   nkeys;     	  /* number of keys stored */
  int32_t   kalloc;	  /* number of keys allocated for */

  int32_t  *cmem;	  /* memory for storing coord2 data */
  int32_t   calloc;	  /* current allocated size of <cmem>, in # of ints */
  int32_t   cn;		  /* current used size of coord2 data */
} P7_COORD2_HASH;

extern P7_COORDS2 *p7_coords2_Create      (int32_t nalloc, int32_t nredline);
extern int         p7_coords2_Grow        (P7_COORDS2 *c2);
extern int         p7_coords2_GrowTo      (P7_COORDS2 *c2, int32_t nalloc);
extern int         p7_coords2_Copy        (const P7_COORDS2 *src, P7_COORDS2 *dst);
extern int         p7_coords2_SetFromTrace(P7_COORDS2 *c2, const P7_TRACE *tr);
extern int         p7_coords2_Reuse       (P7_COORDS2 *c2);
extern void        p7_coords2_Destroy     (P7_COORDS2 *c2);

extern P7_COORD2_HASH *p7_coord2_hash_Create (int32_t init_hashsize, int32_t init_nkeyalloc, int32_t init_calloc);
extern size_t          p7_coord2_hash_Sizeof (const P7_COORD2_HASH *ch);
extern void            p7_coord2_hash_Destroy(P7_COORD2_HASH *ch);
extern int             p7_coord2_hash_Store  (P7_COORD2_HASH *ch, const P7_COORD2 *seg, int32_t nseg, int32_t *opt_index);
extern int             p7_coord2_hash_Get    (const P7_COORD2_HASH *ch, int32_t L, int32_t keyidx, P7_COORDS2 *c2);
extern int             p7_coord2_hash_Dump   (FILE *ofp, const P7_COORD2_HASH *ch);

extern int p7_coords2_Sample(ESL_RANDOMNESS *rng, P7_COORDS2 *c2, int32_t maxseg, int32_t L, int32_t **byp_wrk);

#endif /* P7_COORDS2_INCLUDED */

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

  
  


