/* H4_ANCHORHASH: an auxiliary data structure that the MPAS (most probable anchor
 * set) algorithm uses, to know quickly when it's sampled an anchor set that it's
 * already scored.
 *
 * Patterned on ESL_KEYHASH implementation.
 * Pointerless, to facilitate copying.
 *
 * Each "element" we store is one anchor array, of 2D+1 integers:
 *     D [i0(1),k0(1)]..[i0(D),k0(D)]
 * Each element is associated with a unique integer key, 0..nkeys-1.
 *
 * <amem> consists of packed data for all stored elements. 
 * <an> is the total number of integers stored in amem now.
 * 
 * Data element <k=0..nkeys-1> starts at amem+key_offset[k].  We use offsets, rather
 * than pointers, because we may need to reallocate <amem>, which would invalidate
 * pointers.
 * 
 * To read the data from element <k>, you first read its <D> at
 * *(amem+key_offset[k]), then read <D> integer (i0,k0) pairs.
 * 
 * Anchor sets come with leading (d=0) and trailing (d=D+1) sentinels, set to (0,M+1)
 * and (L+1,0). We assume that all anchor sets are for the same profile/sequence
 * comparison, so the anchorhash just stores L,M itself instead of storing the
 * trailing sentinels for every anchorset.
 * 
 * Each anchorset is hashed, generating a hash value from 0..hashsize-1, and these
 * integer keys are stored in a hashtable of linked lists.  hashtable[z] is the head
 * node of the list for hash value z.  It's either -1 if the list is empty, or a key,
 * 0..nkeys-1.  The linked list is implemented by having <nxt> values for each key;
 * also -1 (no more keys in list) or a key 0..nkeys-1.  Again, implementing this in
 * terms of keys, not pointers, is important because we will need to reallocate
 * memory sometimes.
 * 
 * Example: 
 *   suppose we've stored 2 different anchor sets for seq of length 100
 *   and a model of length 85.
 *   1. Two domains,  (42,61) and (71,53)
 *   2. Three domains,(23,37), (43,62), and (70,52)
 * then:
 *   amem[]  = 2 42 61 71 53 3 23 37 43 62 70 52
 *   an      = 12
 *   aalloc >= 12
 *   
 *   key_offset[] = { 0, 5 }
 *   nkeys        = 2
 *   kalloc       >= 2
 *   
 *   hashtable    = { -1 ... 0 ... -1 ... 1 ... -1 }
 *   nxt[]        = -1 -1
 *   
 *   L            = 100
 *   M            = 85
 *
 * Note: structure is suitable for becoming more general, a hash
 *  of any integer array.
 *
 * We also keep <key_count[k]>, the number of times that <_Store()>
 * has been called for each different key <k>. This allows debug/test
 * code to compare the observed frequency of an anchor set with its
 * calculated probability. Production code doesn't need it, but
 * overhead is negligible.
 *  
 */
#ifndef h4ANCHORHASH_INCLUDED
#define h4ANCHORHASH_INCLUDED

#include "h4_config.h"

#include "h4_anchorset.h"

typedef struct {
  int32_t   L;            // trailing sentinel i0(D+1) of all anchor sets = L+1
  int32_t   M;            //           ... and k0(0) = M+1

  int32_t  *hashtable;	  // [0..hashsize-1]; head node for chains; index of 1st elem, 0..nkeys-1  
  int32_t   hashsize;	  // size of the hashtable (# of buckets). Must be a power of 2.

  int32_t  *key_offset;	  // [k=0..nkeys-1]; key[k]'s data element starts at amem + key_offset[k] 
  int32_t  *key_count;    // how many times key[k] has been _Store()'d     
  int32_t  *nxt;  	  // [k=0..nkeys-1]; nxt[k] = next elem in a hashtable[z] chain, or -1 for end 
  int32_t   nkeys;     	  // number of keys/data elements stored 
  int32_t   kalloc;	  // number of keys allocated for 

  int32_t  *amem;	  // memory for storing anchor array data "elements" 
  int32_t   aalloc;	  // current allocated size of <amem>, in total # of ints (each element of D anchors uses 2D+1) 
  int32_t   an;		  // current used size of <amem> data, in total # of ints 
} H4_ANCHORHASH;

extern H4_ANCHORHASH *h4_anchorhash_Create (void);
extern size_t         h4_anchorhash_Sizeof (const H4_ANCHORHASH *ah);
extern int            h4_anchorhash_Reuse  (H4_ANCHORHASH *ah);
extern void           h4_anchorhash_Destroy(H4_ANCHORHASH *ah);

extern int            h4_anchorhash_Store  (H4_ANCHORHASH *ah, const H4_ANCHORSET *anch, int D0, int32_t *opt_index);
extern int            h4_anchorhash_Get    (const H4_ANCHORHASH *ah, int32_t keyidx, int D0, H4_ANCHORSET *anch);

extern int            h4_anchorhash_Dump   (FILE *ofp, const H4_ANCHORHASH *ah);

#endif /*h4ANCHORHASH_INCLUDED*/

  
