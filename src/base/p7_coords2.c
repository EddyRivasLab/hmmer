/* P7_COORDS2 maintains a resizeable array of coordinate pairs.
 * 
 * The coord pair might be a start/end (i,j) pair on one thing (domain
 * locations, for example), or an (i,k) correspondence between two
 * things (as in anchor points for profile/seq comparison).
 * 
 * Contents:
 *   1. P7_COORDS2: domain or segment start/end coordinates object.
 *   2. Debugging and development, for P7_COORDS2
 *   3. P7_COORDS2_HASH: hash table for storing alternative <P7_COORDS2> data.
 *   4. Debugging and development, for P7_COORDS2_HASH
 *   x. Unit tests.
 *   x. Example driver
 *   x. Copyright and license information.
 */
#include "p7_config.h"

#include <stdlib.h>
#include <limits.h> 		/* INT_MAX */

#include "easel.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "base/p7_coords2.h"

/*****************************************************************
 * 1. The P7_COORDS2 object.
 *****************************************************************/

/* Function:  p7_coords2_Create()
 * Synopsis:  Create a new <P7_COORDS2>
 *
 * Purpose:   Returns a pointer to a newly created <P7_COORDS2> 
 *            object. 
 *            
 *            Caller will typically pass 0's for both arguments, which
 *            means to use default initial allocation sizes: thus,
 *            <p7_coords2_Create(0,0)> is typical.
 *            
 *            To customize allocation sizes, use nonzero arguments for
 *            <nalloc> and/or <nredline>.  <nalloc> is the initially
 *            allocated number of coord pairs; default is
 *            8. <nredline> is the maximum allocation retained after
 *            the object is Reuse()'d; default is 64.
 *            
 * Args:      nalloc   - initial allocation for # of coord pairs; or 0 to use default
 *            nredline - max allocation retained after Reuse();   or 0 to use default           
 *
 * Returns:   pointer to the new <P7_COORDS2>
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_COORDS2 *
p7_coords2_Create(int32_t nalloc, int32_t nredline)
{
  P7_COORDS2 *c2 = NULL;
  int         status;

  ESL_ALLOC(c2, sizeof(P7_COORDS2));
  c2->arr  = NULL;
  c2->n    = 0;

  c2->nalloc   = (nalloc   > 0 ? nalloc   : 8);
  c2->nredline = (nredline > 0 ? nredline : 64);

  ESL_ALLOC(c2->arr, sizeof(P7_COORD2) * (c2->nalloc));

  return c2;

 ERROR:
  p7_coords2_Destroy(c2);
  return NULL;
}

/* Function:  p7_coords2_Grow()
 * Synopsis:  Increase allocation for coord pairs, if needed.
 *
 * Purpose:   Check if there's enough space in <c2> to hold
 *            a new coord pair. If not, increase the allocation
 *            in <c2> by doubling it.
 *
 * Args:      c2  : coord pair array
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_coords2_Grow(P7_COORDS2 *c2)
{
  int status;

  if (c2->n < c2->nalloc) return eslOK;

  ESL_REALLOC(c2->arr, sizeof(P7_COORD2) * c2->nalloc * 2);
  c2->nalloc = c2->nalloc * 2;
  return eslOK;

 ERROR:
  return status;
}

int
p7_coords2_GrowTo(P7_COORDS2 *c2, int32_t nalloc)
{
  int status;

  if (c2->nalloc >= nalloc) return eslOK;

  ESL_REALLOC(c2->arr, sizeof(P7_COORD2) * nalloc);
  c2->nalloc = nalloc;
  return eslOK;

 ERROR:
  return status;
}

int
p7_coords2_Copy(const P7_COORDS2 *src, P7_COORDS2 *dst)
{
  int32_t d;
  int     status;

  if ((status = p7_coords2_GrowTo(dst, src->n)) != eslOK) goto ERROR;
  
  for (d = 0; d < src->n; d++)
    {
      dst->arr[d].n1 = src->arr[d].n1;
      dst->arr[d].n2 = src->arr[d].n2;
    }
  dst->n    = src->n;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_coords2_SetFromTrace()
 * Synopsis:  Convert domain coords from indexed trace to <P7_COORDS2>
 *
 * Purpose:   Given an indexed trace <tr>, convert the domain
 *            start/end coords on the sequence into the 
 *            <P7_COORDS2> object <c2>, using <arr[].n1> for
 *            start and <.n2> for end.
 *
 *            If needed, <c2> will be reallocated to fit the
 *            number of domains in the <tr>.
 *
 *            Note that the trace must be indexed first by the
 *            caller (by calling <p7_trace_Index()>).
 *
 * Args:      c2 : P7_COORDS2 object to copy domain coords into
 *            tr : P7_TRACE object to copy them from
 * 
 * Returns:   <eslOK> on sucess.
 *
 * Throws:    <eslEMEM> on allocation failure. Now the state of <c2>
 *            is undefined.
 *
 * Note:      Dependence on P7_TRACE might not be appropriate for something
 *            in base/. Someday we might need to move this elsewhere.
 */
int
p7_coords2_SetFromTrace(P7_COORDS2 *c2, const P7_TRACE *tr)
{
  int32_t d;
  int     status;

  if ((status = p7_coords2_GrowTo(c2, tr->ndom)) != eslOK) goto ERROR;

  for (d = 0; d < tr->ndom; d++)
    {
      c2->arr[d].n1 = tr->sqfrom[d];
      c2->arr[d].n2 = tr->sqto[d];
    }

  c2->n    = tr->ndom;
  return eslOK;

 ERROR:
  return status;
}

int
p7_coords2_Reuse(P7_COORDS2 *c2)
{
  int status;

  if (c2->nalloc > c2->nredline) {
    ESL_REALLOC(c2->arr, sizeof(P7_COORD2) * c2->nredline);
    c2->nalloc = c2->nredline;
  }
  c2->n    = 0;
  return eslOK;

 ERROR:
  return status;
}

void
p7_coords2_Destroy(P7_COORDS2 *c2)
{
  if (c2) {
    if (c2->arr) free(c2->arr);
    free(c2);
  }
  return;
}

/*****************************************************************
 * 2. Debugging and development tools, P7_COORDS2
 *****************************************************************/

/* Sample random domain segment positions, start/end pairs, sorted and nonoverlapping.
 */
int
p7_coords2_Sample(ESL_RANDOMNESS *rng, P7_COORDS2 *c2, int32_t maxseg, int32_t L, int32_t **byp_wrk)
{
  int32_t *wrk  = NULL;
  int32_t  nseg = 1 + esl_rnd_Roll(rng, maxseg); /* 1..maxseg */
  int32_t  i;
  int      status;

  /* Using the bypass idiom, make sure we have a workspace for <L> coords */
  if      (esl_byp_IsInternal(byp_wrk) ) ESL_ALLOC(wrk, sizeof(int32_t) * L);
  else if (esl_byp_IsReturned(byp_wrk) ) ESL_ALLOC(wrk, sizeof(int32_t) * L);
  else if (esl_byp_IsProvided(byp_wrk) ) { wrk = *byp_wrk; ESL_REALLOC(wrk, sizeof(int32_t) * L); }
			      
  /* We put the numbers 1..L into the workspace <wrk>; shuffle them;
   * then sort the top nseg*2 of them. This gives us <nseg>
   * nonoverlapping start/end coords, in order.
   */
  for (i = 0; i < L; i++) wrk[i] = i+1;
  esl_vec_IShuffle(rng, wrk, L);
  esl_vec_ISortIncreasing(wrk, nseg*2);

  /* Store those randomized coords now in the data structure. */
  p7_coords2_GrowTo(c2, nseg);
  c2->n    = nseg;
  for (i = 0; i < nseg; i++)
    {
      c2->arr[i].n1 = wrk[i*2];
      c2->arr[i].n2 = wrk[i*2+1];
    }
  
  /* Using the bypass idiom, recycle workspace, if we're supposed to */
  if      (esl_byp_IsInternal(byp_wrk)) free(wrk);
  else if (esl_byp_IsReturned(byp_wrk)) *byp_wrk = wrk;
  else if (esl_byp_IsProvided(byp_wrk)) *byp_wrk = wrk;
  return eslOK;

 ERROR:
  if (esl_byp_IsInternal(byp_wrk) && wrk) free(wrk);
  return status;
}


/*****************************************************************
 * 3. The P7_COORDS2_HASH object
 *****************************************************************/
static uint32_t p7_coords2_hash_function    (const P7_COORD2 *seg, int32_t nseg, uint32_t hashsize);
static uint32_t p7_coords2_hash_function_alt(int32_t *keydata, int32_t hashsize);
static int      p7_coords2_hash_compare     (const P7_COORD2 *seg, int32_t nseg,  int32_t *keydata);
static int      p7_coords2_hash_upsize      (P7_COORDS2_HASH *ch);

/* Function:  p7_coords2_hash_Create()
 * Synopsis:  Create a <P7_COORDS2_HASH>
 *
 * Purpose:   Allocate and initialize a <P7_COORDS2_HASH> hash table for storing
 *            lots of coord2 arrays (i.e. domain annotations).
 * 
 *            The <init_*> arguments let you set non-default initial
 *            allocation sizes. To use the default for any of these,
 *            pass a 0 value. Defaults are 128 for the initial 
 *            hashtable size <init_hashsize>; 128 for the initial
 *            allocation for number of keys to be stored <init_nkeyalloc>;
 *            and 2048 for the initial allocation for the number
 *            of integers to be stored in key data. 
 *            
 *            In general the initialization defaults should be
 *            fine. All three are grown automatically as needed, as
 *            you add keys to the hash.
 *            
 *            "key data" means <n> <start>/<end> pairs, plus <n>
 *            itself: it takes 2n+1 integers to store a <P7_COORD2>
 *            array of length <n>.
 *            
 *            <hashsize> must be a power of 2; remember that if you
 *            pass a non-default value.
 *            
 * Args:      init_hashsize : initial hashtable size. Power of 2; >0.
 *            init_keyalloc : initial allocation for # keys. >0.
 *            init_calloc   : initial allocation for key data. >0.
 *
 * Returns:   pointer to the new <P7_COORDS2_HASH> object on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_COORDS2_HASH *
p7_coords2_hash_Create(int32_t init_hashsize, int32_t init_nkeyalloc, int32_t init_calloc)
{
  P7_COORDS2_HASH *ch = NULL;
  int32_t          i;
  int              status;

  ESL_DASSERT1(( init_hashsize == 0 || (init_hashsize && ((init_hashsize & (init_hashsize-1)) == 0)))); /* hashsize is a power of 2 (bitshifting trickery) */
  
  ESL_ALLOC(ch, sizeof(P7_COORDS2_HASH));
  ch->hashtable  = NULL;
  ch->key_offset = NULL;
  ch->nxt        = NULL;
  ch->cmem       = NULL;

  ch->nkeys      = 0;
  ch->cn         = 0;

  ch->hashsize   = (init_hashsize  > 0 ? init_hashsize  : 128);
  ch->kalloc     = (init_nkeyalloc > 0 ? init_nkeyalloc : 128);
  ch->calloc     = (init_calloc    > 0 ? init_calloc    : 2048);
  
  ESL_ALLOC(ch->hashtable, sizeof(int32_t) * ch->hashsize);
  for (i = 0; i < ch->hashsize; i++) ch->hashtable[i] = -1;

  ESL_ALLOC(ch->key_offset, sizeof(int32_t) * ch->kalloc);
  ESL_ALLOC(ch->nxt,        sizeof(int32_t) * ch->kalloc);
  ESL_ALLOC(ch->cmem,       sizeof(int32_t) * ch->calloc);
  return ch;
  
 ERROR:
  p7_coords2_hash_Destroy(ch);
  return NULL;
}

size_t
p7_coords2_hash_Sizeof(const P7_COORDS2_HASH *ch)
{
  size_t n = 0;

  n += sizeof(P7_COORDS2_HASH);
  n += sizeof(int32_t) * ch->hashsize;	 /* hashtable */
  n += sizeof(int32_t) * ch->kalloc * 2; /* key_offset, nxt */
  n += sizeof(int32_t) * ch->calloc;	 /* cmem */
  return n;
}


/* Function:  p7_coords2_hash_Reuse()
 * Synopsis:  Reuse a <P7_COORDS2>
 *
 * Purpose:   Clear a <P7_COORDS2_HASH> hash table for reuse.
 *
 *            If any allocations are overly large, drop them
 *            back to 'redline' values. Default redlines
 *            are 1024 keys (i.e. different coord pair arrays),
 *            1024 hash values, and 16384 total integers of
 *            raw data. Redlines are all 8x the default
 *            initial allocations.
 *
 * Args:      ch :  hash table to reuse
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            (But any reallocations here are shrinkages, so I don't
 *            believe they can fail.)
 */
int
p7_coords2_hash_Reuse(P7_COORDS2_HASH *ch)
{
  int hashsize_redline = 1024;
  int kalloc_redline   = 1024;
  int calloc_redline   = 16384;
  int i;
  int status;

  if (ch->hashsize > hashsize_redline)
    {
      ESL_REALLOC(ch->hashtable, sizeof(int32_t) * hashsize_redline);
      ch->hashsize = hashsize_redline;
    }
  if (ch->kalloc > kalloc_redline)
    { 
      ESL_REALLOC(ch->nxt,        sizeof(int32_t) * kalloc_redline);
      ESL_REALLOC(ch->key_offset, sizeof(int32_t) * kalloc_redline);
      ch->kalloc = kalloc_redline;
    }
  if (ch->calloc > calloc_redline)
    {
      ESL_REALLOC(ch->cmem, sizeof(int32_t) * ch->calloc);
      ch->calloc = calloc_redline;
    }

  for (i = 0; i < ch->hashsize; i++) ch->hashtable[i] = -1;
  ch->nkeys = 0;
  ch->cn    = 0;
  return eslOK;

 ERROR:
  return status;
}



/* Function:  p7_coords2_hash_Destroy()
 * Synopsis:  Destroys a <P7_COORDS2_HASH> hash table.
 */
void
p7_coords2_hash_Destroy(P7_COORDS2_HASH *ch)
{
  if (ch)
    {
      if (ch->hashtable)  free(ch->hashtable);
      if (ch->key_offset) free(ch->key_offset);
      if (ch->nxt)        free(ch->nxt);
      if (ch->cmem)       free(ch->cmem);
      free(ch);
    }
  return;
}


/* Function:  p7_coords2_hash_Store()
 * Synopsis:  Store a <P7_COORDS2> array and get a key index for it.
 *
 * Purpose:   In the hash table <ch>, store the array of coordinate
 *            pairs in <c2>.  Associate it with a unique key index,
 *            counting from 0. This index lets us map the hashed data
 *            to integer-based C arrays. Return the index through <opt_index>.
 *            
 *            If an identical array of paired coords has already been
 *            stored, then set <*opt_index> to the index of where the
 *            data were already stored, and return <eslEDUP>
 *
 * Args:      ch         : hash table holding different arrays of coord pairs
 *            c2         : new array of coord pairs to try to store
 *            opt_index  : optRETURN: index of stored data
 *            
 * Returns:   <eslOK> if <seg>/<nseg> is new; the data are stored, 
 *            and <opt_index>, if requested, is set to the lookup 
 *            key index for the stored data.
 *            
 *            <eslEDUP> if <seg>/<nseg> has already been stored before;
 *            <opt_index>, if requested, is set to the lookup key
 *            index of the previously stored data.
 *
 * Throws:    <eslEMEM> on allocation failure. 
 */
int
p7_coords2_hash_Store(P7_COORDS2_HASH *ch, const P7_COORDS2 *c2, int32_t *opt_index)
{
  uint32_t  val = p7_coords2_hash_function(c2->arr, c2->n, ch->hashsize);
  int32_t  *ptr;
  int32_t   idx;
  int32_t   d;
  int       status;
  
  /* Was this key already stored? */
  for (idx = ch->hashtable[val]; idx != -1; idx = ch->nxt[idx])
    {
      if (p7_coords2_hash_compare(c2->arr, c2->n, ch->cmem + ch->key_offset[idx]) == eslOK)
	{
	  if (opt_index) *opt_index = idx;
	  return eslEDUP;
	}
    }

  /* Reallocate key memory if needed */
  if (ch->nkeys == ch->kalloc)
    {
      ESL_REALLOC(ch->key_offset, sizeof(int32_t) * ch->kalloc * 2);
      ESL_REALLOC(ch->nxt,        sizeof(int32_t) * ch->kalloc * 2);
      ch->kalloc *= 2;
    }

  /* Reallocate key data memory if needed */
  while (ch->cn + 2 * c2->n + 1 > ch->calloc)
    {
      ESL_REALLOC(ch->cmem, sizeof(int32_t) * ch->calloc * 2);
      ch->calloc *= 2;
    }

  /* Copy the key, assign its index */
  idx                 = ch->nkeys;
  ch->key_offset[idx] = ch->cn;
  ch->cn             += 2 * c2->n + 1;
  ch->nkeys++;

  ptr  = ch->cmem + ch->key_offset[idx];
  *ptr = c2->n;
  for (d = 0; d < c2->n; d++) 
    {
      ptr++; *ptr = c2->arr[d].n1;
      ptr++; *ptr = c2->arr[d].n2;
    }

  /* Insert new element at head of the approp chain in hashtable */
  ch->nxt[idx]       = ch->hashtable[val];
  ch->hashtable[val] = idx;

  /* Time to upsize? If we're 3x saturated, expand the hash table */
  if (ch->nkeys > 3 * ch->hashsize)
    if ((status = p7_coords2_hash_upsize(ch)) != eslOK) goto ERROR;

  if (opt_index) *opt_index = idx;
  return eslOK;

 ERROR:
  if (opt_index) *opt_index = -1;
  return status;

}

/* Function:  p7_coords2_hash_Get()
 * Synopsis:  Get a set of coordinate pairs back from the hash.
 *
 * Purpose:   From hash <ch>, retrieve coord pair array <keyidx>,
 *            putting it in <c2>. <c2> is reallocated if necessary.
 *            
 * Args:      ch      : hash storage for alternative annotations
 *            keyidx  : which annotation to get [0..ch->nkeys-1]
 *            c2      : RETURN: coord pair array
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on failure, and the state of <c2> is undefined.
 */
int
p7_coords2_hash_Get(const P7_COORDS2_HASH *ch, int32_t keyidx, P7_COORDS2 *c2)
{
  int32_t *ptr  = ch->cmem + ch->key_offset[keyidx];
  int32_t  n    = *ptr;
  int32_t  d;
  int      status;

  if ((status = p7_coords2_GrowTo(c2, n)) != eslOK) goto ERROR;

  for (d = 0; d < n; d++)
    {
      ptr++; c2->arr[d].n1 = *ptr;
      ptr++; c2->arr[d].n2 = *ptr;
    }
  c2->n    = n;
  return eslOK;

 ERROR:
  return status;
}
  



/* p7_coords2_hash_function()
 *   
 * Given <arr>/<n> data, and the current <hashsize>;
 * calculate and return a hash function on that data, in range <0..hashsize-1>.
 */
static uint32_t
p7_coords2_hash_function(const P7_COORD2 *arr, int32_t n, uint32_t hashsize)
{
  uint32_t hashval = 0;
  int32_t  d;

  hashval = (hashval * 33 + n) % hashsize;
  for (d = 0; d < n; d++)
    {
      hashval = (hashval * 33 + arr[d].n1) % hashsize;
      hashval = (hashval * 33 + arr[d].n2)   % hashsize;
    }
  return hashval;
}

/* p7_coords2_hash_function_alt()
 * 
 * Exactly the same as above (indeed, MUST be the same hash function),
 * but it works on the internal stored version of the <arr>/<n> data,
 * given a pointer to the start of the data.
 * 
 * The data are stored in an array of 2n+1 integers:
 *    <n> <s1> <e1> ... <sn> <en>
 */
static uint32_t
p7_coords2_hash_function_alt(int32_t *keydata, int32_t hashsize)
{
  uint32_t hashval = 0;
  int32_t  d;
  int32_t  n = keydata[0];

  hashval = (hashval * 33 + n) % hashsize;
  for (d = 0; d < n; d++)
    {
      keydata++; hashval = (hashval * 33 + *keydata) % hashsize;
      keydata++; hashval = (hashval * 33 + *keydata) % hashsize;
    }
  return hashval;
}

/* p7_coords2_hash_compare()
 * 
 * Compare a <seg>/<nseg> array (a P7_COORD2 array) to 
 * data from a P7_COORD2 array that's already been stored,
 * starting at <keydata>.
 * 
 * Return <eslOK> if the two are identical; 
 * return <eslFAIL> if not.
 */
static int
p7_coords2_hash_compare(const P7_COORD2 *arr, int n, int32_t *keydata)
{                  /* <keydata> = [ <n> <s1> <e1> ... <sn> <en> */
  int d;
  if (n != *keydata) return eslFAIL;
  for (d = 0; d < n; d++) 
    {
      keydata++; if (*keydata != arr[d].n1) return eslFAIL; 
      keydata++; if (*keydata != arr[d].n2) return eslFAIL; 
    }
  return eslOK;
}
      
/* p7_coords2_hash_upsize()
 * 
 * Increase the hash table size in <ch>, because it's getting
 * too full. This requires recalculating the hash functions for
 * all the previously stored keys, and re-storing them.
 *
 * Throws: <eslEMEM> on allocation failure.
 */
int
p7_coords2_hash_upsize(P7_COORDS2_HASH *ch)
{
  uint32_t val;
  int32_t  i;
  int      status;

  /* 28, because we're going to upsize in steps of 8x, 2^3, so need <2^(31-3) */
  if (ch->hashsize >= (1<<28)) return eslOK; /* quasi-success: don't grow any more */

  /* The catch: upsizing table changes all hash functions, so all
   * keys have to be re-hashed and re-stored. But they can stay
   * where they are in the data storage array.
   */
  ESL_REALLOC(ch->hashtable, sizeof(int32_t) * (ch->hashsize << 3));
  ch->hashsize = ch->hashsize << 3; /* x8 */
  for (i = 0; i < ch->hashsize; i++) 
    ch->hashtable[i] = -1;

  for (i = 0; i < ch->nkeys; i++)
    {
      val        = p7_coords2_hash_function_alt(ch->cmem + ch->key_offset[i], ch->hashsize);
      ch->nxt[i] = ch->hashtable[val];
      ch->hashtable[val] = i;
    }
  return eslOK;

 ERROR:
  return eslEMEM;
}

/*****************************************************************
 * 4. Debugging and development tools, P7_COORDS2_HASH
 *****************************************************************/

int
p7_coords2_hash_Dump(FILE *ofp, const P7_COORDS2_HASH *ch)
{
  int32_t nempty  = 0;
  int32_t maxkeys = -1;
  int32_t minkeys = INT32_MAX;
  int32_t h;
  int32_t idx;
  int32_t n;

  for (h = 0; h < ch->hashsize; h++)
    {
      for (n = 0, idx = ch->hashtable[h]; idx != -1; idx = ch->nxt[idx]) n++;

      if (n == 0) nempty++;
      if (n > maxkeys) maxkeys = n;
      if (n < minkeys) minkeys = n;
    }

  fprintf(ofp, "Total keys:             %d\n", ch->nkeys);
  fprintf(ofp, "Hash table size:        %d\n", ch->hashsize);
  fprintf(ofp, "Average occupancy:      %.2f\n", (float) ch->nkeys /(float) ch->hashsize);
  fprintf(ofp, "Unoccupied slots:       %d\n", nempty);
  fprintf(ofp, "Most in one slot:       %d\n", maxkeys);
  fprintf(ofp, "Least in one slot:      %d\n", minkeys);
  fprintf(ofp, "Keys allocated for:     %d\n", ch->kalloc);
  fprintf(ofp, "Key data space alloc:   %d\n", ch->calloc);
  fprintf(ofp, "Key data space used:    %d\n", ch->cn);
  fprintf(ofp, "Total obj size, bytes:  %d\n", (int) p7_coords2_hash_Sizeof(ch));

  return eslOK;
}


/*****************************************************************
 * x. Unit tests.
 *****************************************************************/


/*****************************************************************
 * x. Example driver
 *****************************************************************/
#ifdef p7COORDS2_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles   reqs   incomp   help                                    docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,   "show brief help on version and usage",        0 },
  { "-s",           eslARG_INT,      "0", NULL, NULL,      NULL,  NULL,  NULL,   "set random number seed to <n>",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "example driver for P7_COORDS2, P7_COORDS2_HASH";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng      = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  P7_COORDS2     *c2       = p7_coords2_Create(0, 0);
  P7_COORDS2_HASH *hash     = p7_coords2_hash_Create(0, 0, 0);
  int             L        = 20;
  int             maxseg   = 1;
  int             nsamples = 1000;
  int32_t        *wrk      = NULL;
  int32_t         keyidx;
  int             i;

  for (i = 0; i < nsamples; i++)
    {
      p7_coords2_Sample(rng, c2, maxseg, L, &wrk);

      p7_coords2_hash_Store(hash, c2->seg, c2->nseg, &keyidx);

      p7_coords2_Reuse(c2);
    }
  
  p7_coords2_hash_Dump(stdout, hash);

  if (wrk) free(wrk);
  p7_coords2_hash_Destroy(hash);
  p7_coords2_Destroy(c2);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7COORDS2_EXAMPLE*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

