/* P7_COORDS2 maintains a resizeable array of start/end coord pairs. 
 * Used for domain segment coords in a target sequence, for example.
 * 
 * Contents:
 *   1. P7_COORDS2: domain or segment start/end coordinates object.
 *   2. Debugging and development, for P7_COORDS2
 *   3. P7_COORD2_HASH: hash table for storing alternative <P7_COORDS2> data.
 *   4. Debugging and development, for P7_COORD2_HASH
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

/* Function:  
 * Synopsis:  
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Note:      Need to document redline, and convention
 *            for default allocation size, in codestyle.
 */

P7_COORDS2 *
p7_coords2_Create(int nalloc, int nredline)
{
  P7_COORDS2 *c2 = NULL;
  int         status;

  ESL_ALLOC(c2, sizeof(P7_COORDS2));
  c2->seg  = NULL;
  c2->nseg = 0;
  c2->L    = 0;  

  c2->nalloc   = (nalloc   > 0 ? nalloc   : 8);
  c2->nredline = (nredline > 0 ? nredline : 64);

  ESL_ALLOC(c2->seg, sizeof(P7_COORD2) * (c2->nalloc));

  return c2;

 ERROR:
  p7_coords2_Destroy(c2);
  return NULL;
}

int
p7_coords2_Grow(P7_COORDS2 *c2)
{
  int status;

  if (c2->nseg < c2->nalloc) return eslOK;

  ESL_REALLOC(c2->seg, sizeof(P7_COORD2) * c2->nalloc * 2);
  c2->nalloc = c2->nalloc * 2;
  return eslOK;

 ERROR:
  return status;
}

int
p7_coords2_GrowTo(P7_COORDS2 *c2, int nalloc)
{
  int status;

  if (c2->nalloc >= nalloc) return eslOK;

  ESL_REALLOC(c2->seg, sizeof(P7_COORD2) * nalloc);
  c2->nalloc = nalloc;
  return eslOK;

 ERROR:
  return status;
}

int
p7_coords2_Copy(const P7_COORDS2 *src, P7_COORDS2 *dst)
{
  int d;
  int status;

  if ((status = p7_coords2_GrowTo(dst, src->nseg)) != eslOK) goto ERROR;
  
  for (d = 0; d < src->nseg; d++)
    {
      dst->seg[d].start = src->seg[d].start;
      dst->seg[d].end   = src->seg[d].end;
    }
  dst->nseg = src->nseg;
  dst->L    = src->L;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_coords2_SetFromTrace()
 * Synopsis:  Convert domain coords from indexed trace to <P7_COORDS2>
 *
 * Purpose:   Given an indexed trace <tr>, convert the domain
 *            start/end coords on the sequence into the 
 *            <P7_COORDS2> object <c2>. 
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
  int d;
  int status;

  if ((status = p7_coords2_GrowTo(c2, tr->ndom)) != eslOK) goto ERROR;

  for (d = 0; d < tr->ndom; d++)
    {
      c2->seg[d].start = tr->sqfrom[d];
      c2->seg[d].end   = tr->sqto[d];
    }

  c2->nseg = tr->ndom;
  c2->L    = tr->L;
  return eslOK;

 ERROR:
  return status;
}

int
p7_coords2_Reuse(P7_COORDS2 *c2)
{
  int status;

  if (c2->nalloc > c2->nredline) 
    {
      ESL_REALLOC(c2->seg, sizeof(P7_COORD2) * c2->nredline);
      c2->nalloc = c2->nredline;
    }

  c2->nseg = 0;
  c2->L    = 0;
  return eslOK;

 ERROR:
  return status;
}

void
p7_coords2_Destroy(P7_COORDS2 *c2)
{
  if (c2) 
    {
      if (c2->seg) free(c2->seg);
      free(c2);
    }
  return;
}

/*****************************************************************
 * 2. Debugging and development tools, P7_COORDS2
 *****************************************************************/

int
p7_coords2_Sample(ESL_RANDOMNESS *rng, P7_COORDS2 *c2, int maxseg, int L, int **byp_wrk)
{
  int *wrk  = NULL;
  int  nseg = 1 + esl_rnd_Roll(rng, maxseg); /* 1..maxseg */
  int  i;
  int  status;

  /* Using the bypass idiom, make sure we have a workspace for <L> coords */
  if      (esl_byp_IsInternal(byp_wrk) ) ESL_ALLOC(wrk, sizeof(int) * L);
  else if (esl_byp_IsReturned(byp_wrk) ) ESL_ALLOC(wrk, sizeof(int) * L);
  else if (esl_byp_IsProvided(byp_wrk) ) { wrk = *byp_wrk; ESL_REALLOC(wrk, sizeof(int) * L); }
			      
  /* We put the numbers 1..L into the workspace <wrk>; shuffle them;
   * then sort the top nseg*2 of them. This gives us <nseg>
   * nonoverlapping start/end coords, in order.
   */
  for (i = 0; i < L; i++) wrk[i] = i+1;
  esl_vec_IShuffle(rng, wrk, L);
  esl_vec_ISortIncreasing(wrk, nseg*2);

  /* Store those randomized coords now in the data structure. */
  p7_coords2_GrowTo(c2, nseg);
  c2->L    = L;
  c2->nseg = nseg;
  for (i = 0; i < nseg; i++)
    {
      c2->seg[i].start = wrk[i*2];
      c2->seg[i].end   = wrk[i*2+1];
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
 * 3. The P7_COORD2_HASH object
 *****************************************************************/
static uint32_t p7_coord2_hash_function    (const P7_COORD2 *seg, int nseg, uint32_t hashsize);
static uint32_t p7_coord2_hash_function_alt(int32_t *keydata, int32_t hashsize);
static int      p7_coord2_hash_compare     (const P7_COORD2 *seg, int nseg,  int32_t *keydata);
static int      p7_coord2_hash_upsize      (P7_COORD2_HASH *ch);

/* Function:  p7_coord2_hash_Create()
 * Synopsis:  Create a <P7_COORD2_HASH>
 *
 * Purpose:   Allocate and initialize a <P7_COORD2_HASH> hash table for storing
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
 * Returns:   pointer to the new <P7_COORD2_HASH> object on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_COORD2_HASH *
p7_coord2_hash_Create(int32_t init_hashsize, int32_t init_nkeyalloc, int32_t init_calloc)
{
  P7_COORD2_HASH *ch = NULL;
  int             i;
  int             status;

  ESL_DASSERT1(( init_hashsize == 0 || (init_hashsize && ((init_hashsize & (init_hashsize-1)) == 0)))); /* hashsize is a power of 2 (bitshifting trickery) */
  
  ESL_ALLOC(ch, sizeof(P7_COORD2_HASH));
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
  p7_coord2_hash_Destroy(ch);
  return NULL;
}

size_t
p7_coord2_hash_Sizeof(const P7_COORD2_HASH *ch)
{
  size_t n = 0;

  n += sizeof(P7_COORD2_HASH);
  n += sizeof(int32_t) * ch->hashsize;	 /* hashtable */
  n += sizeof(int32_t) * ch->kalloc * 2; /* key_offset, nxt */
  n += sizeof(int32_t) * ch->calloc;	 /* cmem */
  return n;
}



/* Function:  p7_coord2_hash_Destroy()
 * Synopsis:  Destroys a <P7_COORD2_HASH> hash table.
 */
void
p7_coord2_hash_Destroy(P7_COORD2_HASH *ch)
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


/* Function:  p7_coord2_hash_Store()
 * Synopsis:  Store a <P7_COORD2> array and get a key index for it.
 *
 * Purpose:   In the hash table <ch>, store <P7_COORD2> array <seg> containing
 *            <nseg> start/end pairs. Associate it with a unique key index, 
 *            counting from 0. This index lets us map the hashed data to
 *            integer-based C arrays. 
 *            
 *            If an identical <seg>,<nseg> array has already been
 *            stored, do nothing, and return <eslEDUP>.
 *            
 *            Optionally, return the index through <opt_index>.
 *
 * Args:      ch
 *            seg
 *            nseg
 *            opt_index 
 *            
 * Returns:   <eslOK> if <seg>/<nseg> is new; the data are stored, 
 *            and <opt_index>, if requested, is set to the lookup 
 *            key index for the stored data.
 *            
 *            <eslEDUP> if <seg>/<nseg> has already been stored before;
 *            <opt_index>, if requested, is set to the lookup key
 *            index of the previously stored data.
 *
 * Throws:    <eslEMEM> on allocation failure; <opt_index>, if requested,
 *            is -1.
 */
int
p7_coord2_hash_Store(P7_COORD2_HASH *ch, const P7_COORD2 *seg, int nseg, int *opt_index)
{
  uint32_t  val = p7_coord2_hash_function(seg, nseg, ch->hashsize);
  int32_t  *ptr;
  int       idx;
  int       d;
  int       status;
  
  /* Was this key already stored? */
  for (idx = ch->hashtable[val]; idx != -1; idx = ch->nxt[idx])
    {
      if (p7_coord2_hash_compare(seg, nseg, ch->cmem + ch->key_offset[idx]) == eslOK)
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
  while (ch->cn + 2*nseg + 1 > ch->calloc)
    {
      ESL_REALLOC(ch->cmem, sizeof(int32_t) * ch->calloc * 2);
      ch->calloc *= 2;
    }

  /* Copy the key, assign its index */
  idx                 = ch->nkeys;
  ch->key_offset[idx] = ch->cn;
  ch->cn             += 2*nseg + 1;
  ch->nkeys++;

  ptr  = ch->cmem + ch->key_offset[idx];
  *ptr = nseg;
  for (d = 0; d < nseg; d++) 
    {
      ptr++; *ptr = seg[d].start;
      ptr++; *ptr = seg[d].end;
    }

  /* Insert new element at head of the approp chain in hashtable */
  ch->nxt[idx]       = ch->hashtable[val];
  ch->hashtable[val] = idx;

  /* Time to upsize? If we're 3x saturated, expand the hash table */
  if (ch->nkeys > 3 * ch->hashsize)
    if ((status = p7_coord2_hash_upsize(ch)) != eslOK) goto ERROR;

  if (opt_index) *opt_index = idx;
  return eslOK;

 ERROR:
  if (opt_index) *opt_index = -1;
  return status;

}

/* Function:  p7_coord2_hash_Get()
 * Synopsis:  Get a set of segment coords back from the hash.
 *
 * Purpose:   From hash <ch>, retrieve annotation number <keyidx>,
 *            for a sequence of length <L>. Store this segment
 *            set in <c2>. <c2> is reallocated if necessary.
 *            
 * Args:      ch      : hash storage for alternative annotations
 *            L       : length of the annotated sequence
 *            keyidx  : which annotation to get [0..ch->nkeys-1]
 *            c2      : RETURN: domain segment coords
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on failure, and the state of <c2> is undefined.
 */
int
p7_coord2_hash_Get(const P7_COORD2_HASH *ch, int L, int keyidx, P7_COORDS2 *c2)
{
  int32_t *ptr  = ch->cmem + ch->key_offset[keyidx];
  int      nseg = *ptr;
  int      d;
  int      status;

  if ((status = p7_coords2_GrowTo(c2, nseg)) != eslOK) goto ERROR;

  for (d = 0; d < nseg; d++)
    {
      ptr++; c2->seg[d].start = *ptr;
      ptr++; c2->seg[d].end   = *ptr;
    }
  c2->nseg = nseg;
  c2->L    = L;
  return eslOK;

 ERROR:
  return status;
}
  



/* p7_coord2_hash_function()
 *   
 * Given <seg>/<nseg> data, and the current <hashsize>;
 * calculate and return a hash function on that data, in range <0..hashsize-1>.
 */
static uint32_t
p7_coord2_hash_function(const P7_COORD2 *seg, int nseg, uint32_t hashsize)
{
  uint32_t hashval = 0;
  int      d;

  hashval = (hashval * 33 + nseg) % hashsize;
  for (d = 0; d < nseg; d++)
    {
      hashval = (hashval * 33 + seg[d].start) % hashsize;
      hashval = (hashval * 33 + seg[d].end)   % hashsize;
    }
  return hashval;
}

/* p7_coord2_hash_function_alt()
 * 
 * Exactly the same as above (indeed, MUST be the same hash function),
 * but it works on the internal stored version of the <seg>/<nseg> data,
 * given a pointer to the start of the data.
 * 
 * The data are stored in an array of 2n+1 integers:
 *    <n> <s1> <e1> ... <sn> <en>
 */
static uint32_t
p7_coord2_hash_function_alt(int32_t *keydata, int32_t hashsize)
{
  uint32_t hashval = 0;
  int      d;
  int      nseg    = (int) keydata[0];

  hashval = (hashval * 33 + nseg) % hashsize;
  for (d = 0; d < nseg; d++)
    {
      keydata++; hashval = (hashval * 33 + *keydata) % hashsize;
      keydata++; hashval = (hashval * 33 + *keydata) % hashsize;
    }
  return hashval;
}

/* p7_coord2_hash_compare()
 * 
 * Compare a <seg>/<nseg> array (a P7_COORD2 array) to 
 * data from a P7_COORD2 array that's already been stored,
 * starting at <keydata>.
 * 
 * Return <eslOK> if the two are identical; 
 * return <eslFAIL> if not.
 */
static int
p7_coord2_hash_compare(const P7_COORD2 *seg, int nseg, int32_t *keydata)
{                  /* <keydata> = [ <n> <s1> <e1> ... <sn> <en> */
  int d;
  if (nseg != (int) *keydata) return eslFAIL;
  for (d = 0; d < nseg; d++) 
    {
      keydata++; if ((int) *keydata != seg[d].start) return eslFAIL; 
      keydata++; if ((int) *keydata != seg[d].end)   return eslFAIL; 
    }
  return eslOK;
}
      
/* p7_coord2_hash_upsize()
 * 
 * Increase the hash table size in <ch>, because it's getting
 * too full. This requires recalculating the hash functions for
 * all the previously stored keys, and re-storing them.
 *
 * Throws: <eslEMEM> on allocation failure.
 */
int
p7_coord2_hash_upsize(P7_COORD2_HASH *ch)
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
      val        = p7_coord2_hash_function_alt(ch->cmem + ch->key_offset[i], ch->hashsize);
      ch->nxt[i] = ch->hashtable[val];
      ch->hashtable[val] = i;
    }
  return eslOK;

 ERROR:
  return eslEMEM;
}

/*****************************************************************
 * 4. Debugging and development tools, P7_COORD2_HASH
 *****************************************************************/

int
p7_coord2_hash_Dump(FILE *ofp, const P7_COORD2_HASH *ch)
{
  int nempty  = 0;
  int maxkeys = -1;
  int minkeys = INT_MAX;
  int h;
  int idx;
  int n;

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
  fprintf(ofp, "Total obj size, bytes:  %d\n", (int) p7_coord2_hash_Sizeof(ch));

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
  P7_COORD2_HASH *hash     = p7_coord2_hash_Create(0, 0, 0);
  int             L        = 20;
  int             maxseg   = 1;
  int             nsamples = 1000;
  int            *wrk      = NULL;
  int             keyidx;
  int             i;

  for (i = 0; i < nsamples; i++)
    {
      p7_coords2_Sample(rng, c2, maxseg, L, &wrk);

      p7_coord2_hash_Store(hash, c2->seg, c2->nseg, &keyidx);

      p7_coords2_Reuse(c2);
    }
  
  p7_coord2_hash_Dump(stdout, hash);

  if (wrk) free(wrk);
  p7_coord2_hash_Destroy(hash);
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

