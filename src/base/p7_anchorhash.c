/* P7_ANCHORHASH implements an auxiliary data structure that the
 * MPAS (most probable anchor set) algorithm uses.
 * 
 * Contents:
 *   1. P7_ANCHORHASH object
 *   2. Debugging and development tools
 *   3. Unit tests
 *   4. Test driver
 *   5. Copyright and license information
 */
#include "p7_config.h"

#include <stdlib.h>
#include <limits.h>  // for INT32_MAX

#include "easel.h"

#include "base/p7_anchors.h"
#include "base/p7_anchorhash.h"


/*****************************************************************
 * 1. The P7_ANCHORHASH object
 *****************************************************************/
static uint32_t anchorhash_function    (const P7_ANCHOR *arr, int32_t D, uint32_t hashsize);
static uint32_t anchorhash_function_alt(int32_t *keydata, int32_t hashsize);
static int      anchorhash_compare     (const P7_ANCHOR *arr, int32_t D,  int32_t *keydata);
static int      anchorhash_upsize      (P7_ANCHORHASH *ah);

/* Function:  p7_anchorhash_Create()
 * Synopsis:  Create a <P7_ANCHORHASH>
 *
 * Purpose:   Allocate and initialize a <P7_ANCHORHASH> hash table for storing
 *            lots of anchor arrays.
 *
 * Returns:   pointer to the new <P7_ANCHORHASH> object on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_ANCHORHASH *
p7_anchorhash_Create(void)
{
  P7_ANCHORHASH *ah                = NULL;
  int32_t        default_hashsize  = 128;  // Size of the hash table. Must be a power of 2; grows by doubling.
  int32_t        default_nkeyalloc = 128;  // Number of keys (anchor sets)
  int32_t        default_aalloc    = 2048; // Total number of integers in data elements (2D+1 for each anchor set of D anchors)
  int32_t        i;
  int            status;

  ESL_ALLOC(ah, sizeof(P7_ANCHORHASH));
  ah->hashtable  = NULL;
  ah->key_offset = NULL;
  ah->nxt        = NULL;
  ah->amem       = NULL;

  ah->L          = 0;
  ah->M          = 0;
  ah->nkeys      = 0;
  ah->an         = 0;

  ah->hashsize   = default_hashsize;
  ah->kalloc     = default_nkeyalloc;
  ah->aalloc     = default_aalloc;
  
  ESL_ALLOC(ah->hashtable, sizeof(int32_t) * ah->hashsize);
  for (i = 0; i < ah->hashsize; i++) ah->hashtable[i] = -1;

  ESL_ALLOC(ah->key_offset, sizeof(int32_t) * ah->kalloc);
  ESL_ALLOC(ah->nxt,        sizeof(int32_t) * ah->kalloc);
  ESL_ALLOC(ah->amem,       sizeof(int32_t) * ah->aalloc);
  return ah;
  
 ERROR:
  p7_anchorhash_Destroy(ah);
  return NULL;
}

size_t
p7_anchorhash_Sizeof(const P7_ANCHORHASH *ah)
{
  size_t n = 0;

  n += sizeof(P7_ANCHORHASH);
  n += sizeof(int32_t) * ah->hashsize;	 // hashtable 
  n += sizeof(int32_t) * ah->kalloc * 2; // key_offset, nxt 
  n += sizeof(int32_t) * ah->aalloc;	 // amem 
  return n;
}


/* Function:  p7_anchorhash_Reuse()
 * Synopsis:  Reuse a <P7_ANCHORHASH>
 *
 * Purpose:   Clear a <P7_ANCHORHASH> hash table for reuse.
 *
 *            If any allocations are overly large, drop them
 *            back to 'redline' values. Default redlines
 *            are 1024 keys (i.e. different anchor sets)
 *            1024 hash values, and 16384 total integers of
 *            raw data. Redlines are all 8x the default
 *            initial allocations.
 *
 * Args:      ah :  hash table to reuse
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            (But any reallocations here are shrinkages, so I don't
 *            believe they can fail.)
 */
int
p7_anchorhash_Reuse(P7_ANCHORHASH *ah)
{
  int hashsize_redline = 1024;
  int kalloc_redline   = 1024;
  int aalloc_redline   = 16384;
  int i;
  int status;

  if (ah->hashsize > hashsize_redline)
    {
      ESL_REALLOC(ah->hashtable, sizeof(int32_t) * hashsize_redline);
      ah->hashsize = hashsize_redline;
    }
  if (ah->kalloc > kalloc_redline)
    { 
      ESL_REALLOC(ah->nxt,        sizeof(int32_t) * kalloc_redline);
      ESL_REALLOC(ah->key_offset, sizeof(int32_t) * kalloc_redline);
      ah->kalloc = kalloc_redline;
    }
  if (ah->aalloc > aalloc_redline)
    {
      ESL_REALLOC(ah->amem, sizeof(int32_t) * ah->aalloc);
      ah->aalloc = aalloc_redline;
    }

  for (i = 0; i < ah->hashsize; i++) ah->hashtable[i] = -1;
  ah->L     = 0;
  ah->M     = 0;
  ah->nkeys = 0;
  ah->an    = 0;
  return eslOK;

 ERROR:
  return status;
}



/* Function:  p7_anchorhash_Destroy()
 * Synopsis:  Destroys a <P7_ANCHORHASH> hash table.
 */
void
p7_anchorhash_Destroy(P7_ANCHORHASH *ah)
{
  if (ah)
    {
      if (ah->hashtable)  free(ah->hashtable);
      if (ah->key_offset) free(ah->key_offset);
      if (ah->nxt)        free(ah->nxt);
      if (ah->amem)       free(ah->amem);
      free(ah);
    }
  return;
}


/* Function:  p7_anchorhash_Store()
 * Synopsis:  Store a <P7_ANCHORS> array and get a key index for it.
 *
 * Purpose:   Try to store anchor set <anch> in hash table <ah>.
 *            Associate it with a unique key index, counting from
 *            0. This index lets us map the hashed data to
 *            integer-based C arrays. Return the index through
 *            <opt_index>.
 *            
 *            If an identical anchor set is already stored in <ah>,
 *            set <*opt_index> to the key for that anchor set, and
 *            return <eslEDUP>.
 *
 * Args:      ah         : hash table holding different anchor sets
 *            anch       : new anchor set to try to store
 *            opt_index  : optRETURN: index of stored data
 *            
 * Returns:   <eslOK> if <anch>/<D> is new; the anchor set data are stored, 
 *            and <opt_index>, if requested, is set to the lookup 
 *            key index for the stored data.
 *            
 *            <eslEDUP> if this anchor set has already been stored before;
 *            <opt_index>, if requested, is set to the lookup key
 *            index of the previously stored data.
 *
 * Throws:    <eslEMEM> on allocation failure. 
 */
int
p7_anchorhash_Store(P7_ANCHORHASH *ah, const P7_ANCHORS *anch, int32_t *opt_index)
{
  uint32_t  val = anchorhash_function(anch->a, anch->D, ah->hashsize);
  int32_t  *ptr;
  int32_t   idx;
  int32_t   d;
  int       status;
  
  /* Was this key already stored? */
  for (idx = ah->hashtable[val]; idx != -1; idx = ah->nxt[idx])
    {
      if (anchorhash_compare(anch->a, anch->D, ah->amem + ah->key_offset[idx]) == eslOK)
	{
	  if (opt_index) *opt_index = idx;
	  return eslEDUP;
	}
    }

  /* Reallocate key memory if needed */
  if (ah->nkeys == ah->kalloc)
    {
      ESL_REALLOC(ah->key_offset, sizeof(int32_t) * ah->kalloc * 2);
      ESL_REALLOC(ah->nxt,        sizeof(int32_t) * ah->kalloc * 2);
      ah->kalloc *= 2;
    }

  /* Reallocate key data memory if needed (by doubling) */
  while (ah->an + 2 * anch->D + 1 > ah->aalloc)
    {
      ESL_REALLOC(ah->amem, sizeof(int32_t) * ah->aalloc * 2);
      ah->aalloc *= 2;
    }

  /* Copy the key, assign its index */
  idx                 = ah->nkeys;
  ah->key_offset[idx] = ah->an;
  ah->an             += 2 * anch->D + 1;
  ah->nkeys++;

  ptr  = ah->amem + ah->key_offset[idx];
  *ptr = anch->D;
  for (d = 1; d <= anch->D; d++) 
    {
      ptr++; *ptr = anch->a[d].i0;
      ptr++; *ptr = anch->a[d].k0;
    }
  
  /* Trailing sentinel */
  if (ah->nkeys == 1)    // first one?
    {
      ah->L = anch->a[anch->D+1].i0 - 1;  // L+1 -> stored as L
      ah->M = anch->a[0].k0 - 1;          // M+1 -> stored as M
    }
  ESL_DASSERT1(( anch->a[anch->D+1].i0 = ah->L+1 ));
  ESL_DASSERT1(( anch->a[0].k0         = ah->M+1 ));

  /* Insert new element at head of the approp chain in hashtable */
  ah->nxt[idx]       = ah->hashtable[val];
  ah->hashtable[val] = idx;

  /* Time to upsize? If we're 3x saturated, expand the hash table */
  if (ah->nkeys > 3 * ah->hashsize)
    if ((status = anchorhash_upsize(ah)) != eslOK) goto ERROR;

  if (opt_index) *opt_index = idx;
  return eslOK;

 ERROR:
  if (opt_index) *opt_index = -1;
  return status;

}

/* Function:  p7_anchorhash_Get()
 * Synopsis:  Get an anchor set back from the hash.
 *
 * Purpose:   Retrieve anchor set <keyidx> from hash <ah>, and put
 *            it in <anch>. <anch> is reallocated if necessary.
 *            
 * Args:      ah      : hash storage for alternative annotations
 *            keyidx  : which annotation to get [0..ah->nkeys-1]
 *            anch    : RETURN: anchor set
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on failure, and the state of <anch> is undefined.
 */
int
p7_anchorhash_Get(const P7_ANCHORHASH *ah, int32_t keyidx, P7_ANCHORS *anch)
{
  int32_t *ptr  = ah->amem + ah->key_offset[keyidx];
  int32_t  D    = *ptr;
  int32_t  d;
  int      status;

  if ((status = p7_anchors_Reinit(anch, D)) != eslOK) goto ERROR;

  for (d = 1; d <= D; d++)
    {
      ptr++; anch->a[d].i0 = *ptr;
      ptr++; anch->a[d].k0 = *ptr;
    }

  p7_anchor_SetSentinels(anch->a, D, ah->L, ah->M);
  anch->D = D;
  return eslOK;

 ERROR:
  return status;
}
  



/* p7_anchorhash_function()
 *   
 * Given <arr>/<n> data, and the current <hashsize>;
 * calculate and return a hash function on that data, in range <0..hashsize-1>.
 */
static uint32_t
anchorhash_function(const P7_ANCHOR *arr, int32_t D, uint32_t hashsize)
{
  uint32_t hashval = 0;
  int32_t  d;

  hashval = (hashval * 33 + D) % hashsize;
  for (d = 1; d <= D; d++)
    {
      hashval = (hashval * 33 + arr[d].i0) % hashsize;
      hashval = (hashval * 33 + arr[d].k0) % hashsize;
    }
  return hashval;
}

/* anchorhash_function_alt()
 * 
 * Exactly the same as above (indeed, MUST be the same hash function),
 * but it works on the internal stored version of the <arr>/<n> data,
 * given a pointer to the start of the data.
 * 
 * The data are stored in an array of 2D+1 integers:
 *    <D> <(i0,k0)[1]> ... <(i0,k0)[2]>
 */
static uint32_t
anchorhash_function_alt(int32_t *keydata, int32_t hashsize)
{
  uint32_t hashval = 0;
  int32_t  d;
  int32_t  D = keydata[0];

  hashval = (hashval * 33 + D) % hashsize;
  for (d = 1; d <= D; d++)
    {
      keydata++; hashval = (hashval * 33 + *keydata) % hashsize;
      keydata++; hashval = (hashval * 33 + *keydata) % hashsize;
    }
  return hashval;
}

/* anchorhash_compare()
 * 
 * Compare a <arr>/<D> anchor set (a P7_ANCHOR array) to 
 * data from a P7_ANCHOR array that's already been stored,
 * starting at <keydata>.
 * 
 * Return <eslOK> if the two are identical; 
 * return <eslFAIL> if not.
 */
static int
anchorhash_compare(const P7_ANCHOR *arr, int D, int32_t *keydata)
{                  /* <keydata> = [ <D> <i0 k0> ... <i0 k0> */
  int d;
  if (D != *keydata) return eslFAIL;
  for (d = 1; d <= D; d++) 
    {
      keydata++; if (*keydata != arr[d].i0) return eslFAIL; 
      keydata++; if (*keydata != arr[d].k0) return eslFAIL; 
    }
  return eslOK;
}
      
/* anchorhash_upsize()
 * 
 * Increase the hash table size in <ah>, because it's getting
 * too full. This requires recalculating the hash functions for
 * all the previously stored keys, and re-storing them.
 *
 * Throws: <eslEMEM> on allocation failure.
 */
int
anchorhash_upsize(P7_ANCHORHASH *ah)
{
  uint32_t val;
  int32_t  i;
  int      status;

  /* 28, because we're going to upsize in steps of 8x, 2^3, so need <2^(31-3) */
  if (ah->hashsize >= (1<<28)) return eslOK; /* quasi-success: don't grow any more */

  /* The catch: upsizing table changes all hash functions, so all
   * keys have to be re-hashed and re-stored. But they can stay
   * where they are in the data storage array.
   */
  ESL_REALLOC(ah->hashtable, sizeof(int32_t) * (ah->hashsize << 3));
  ah->hashsize = ah->hashsize << 3; /* x8 */
  for (i = 0; i < ah->hashsize; i++) 
    ah->hashtable[i] = -1;

  for (i = 0; i < ah->nkeys; i++)
    {
      val        = anchorhash_function_alt(ah->amem + ah->key_offset[i], ah->hashsize);
      ah->nxt[i] = ah->hashtable[val];
      ah->hashtable[val] = i;
    }
  return eslOK;

 ERROR:
  return eslEMEM;
}

/*****************************************************************
 * 2. Debugging and development tools
 *****************************************************************/

int
p7_anchorhash_Dump(FILE *ofp, const P7_ANCHORHASH *ah)
{
  int32_t nempty  = 0;
  int32_t maxkeys = -1;
  int32_t minkeys = INT32_MAX;
  int32_t h;
  int32_t idx;
  int32_t n;

  for (h = 0; h < ah->hashsize; h++)
    {
      for (n = 0, idx = ah->hashtable[h]; idx != -1; idx = ah->nxt[idx]) n++;

      if (n == 0) nempty++;
      if (n > maxkeys) maxkeys = n;
      if (n < minkeys) minkeys = n;
    }

  fprintf(ofp, "Total keys:             %d\n", ah->nkeys);
  fprintf(ofp, "Hash table size:        %d\n", ah->hashsize);
  fprintf(ofp, "Average occupancy:      %.2f\n", (float) ah->nkeys /(float) ah->hashsize);
  fprintf(ofp, "Unoccupied slots:       %d\n", nempty);
  fprintf(ofp, "Most in one slot:       %d\n", maxkeys);
  fprintf(ofp, "Least in one slot:      %d\n", minkeys);
  fprintf(ofp, "Keys allocated for:     %d\n", ah->kalloc);
  fprintf(ofp, "Key data space alloc:   %d\n", ah->aalloc);
  fprintf(ofp, "Key data space used:    %d\n", ah->an);
  fprintf(ofp, "Total obj size, bytes:  %d\n", (int) p7_anchorhash_Sizeof(ah));

  return eslOK;
}


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7ANCHORHASH_TESTDRIVE

#include "hmmer.h"

static void
utest_sampling(ESL_RANDOMNESS *rng)
{
  char           msg[]    = "p7_anchorhash.c :: sampling unit test failed";
  int            nsamples = 1000;
  int            M        = 100;
  int            L        = 200;
  int            maxD     = 50;  
  P7_ANCHORHASH *ah       = p7_anchorhash_Create();
  P7_ANCHORS    *anch     = p7_anchors_Create();
  P7_ANCHORS   **aarr     = malloc(sizeof(P7_ANCHORS *) * nsamples);
  int32_t       *keys     = malloc(sizeof(int32_t) * nsamples);
  int32_t        keyidx;        
  int32_t        nk;
  int            iteration;
  int            s;
  char           errmsg[eslERRBUFSIZE];
  int            status;

  for (iteration = 0; iteration < 2; iteration++)  // do it twice: this tests _Reuse()
    {
      for (s = 0; s < nsamples; s++)
	aarr[s] = p7_anchors_Create();

      nk = 0;
      for (s = 0; s < nsamples; s++)
	{
	  if ( p7_anchors_Sample(rng, L, M, maxD, anch) != eslOK) esl_fatal(msg);
	  if ( p7_anchors_Copy(anch, aarr[s])           != eslOK) esl_fatal(msg);

	  status = p7_anchorhash_Store(ah, anch, &keyidx);
	  keys[s] = keyidx;
	  if      (status == eslOK)   { if (keyidx != nk) esl_fatal(msg); nk++; }
	  else if (status == eslEDUP) { if (keyidx >= nk) esl_fatal(msg);       }
	  else                        {                   esl_fatal(msg);       }

	  p7_anchors_Reuse(anch);
	}

      for (s = 0; s < nsamples; s++)
	{
	  if ( p7_anchorhash_Get(ah, keys[s], anch) != eslOK) esl_fatal(msg);

	  //p7_anchors_Dump(stdout, anch);

	  if ( p7_anchors_Validate(anch, errmsg)    != eslOK) esl_fatal("%s:\n  %s", msg, errmsg);
	  if ( p7_anchors_Compare(anch, aarr[s])    != eslOK) esl_fatal(msg);
	  if ( p7_anchors_Reuse(anch)               != eslOK) esl_fatal(msg);
	}

      for (s = 0; s < nsamples; s++)
	p7_anchors_Destroy(aarr[s]);

      p7_anchorhash_Reuse(ah);  
    }

  p7_anchorhash_Destroy(ah);
  p7_anchors_Destroy(anch);
  free(aarr);
  free(keys);
}


#endif /*p7ANCHORHASH_TESTDRIVE*/

/*****************************************************************
 * 4. Test driver
 *****************************************************************/

#ifdef p7ANCHORHASH_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for p7_anchorhash.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_sampling(rng);

  fprintf(stderr, "#  status = ok\n");

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  exit(0); 
}

#endif /*p7ANCHORS_TESTDRIVE*/
/*-------------------- end of test driver ---------------------*/

/*****************************************************************
 * 5. Example driver
 *****************************************************************/
#ifdef p7ANCHORHASH_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles   reqs   incomp   help                                    docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,   "show brief help on version and usage",        0 },
  { "-s",           eslARG_INT,      "0", NULL, NULL,      NULL,  NULL,  NULL,   "set random number seed to <n>",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "example driver for P7_ANCHORHASH";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng      = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  P7_ANCHORS     *anch     = p7_anchors_Create();
  P7_ANCHORHASH  *ah       = p7_anchorhash_Create();
  int             L        = 400;
  int             M        = 200;
  int             maxD     = 50;
  int             nsamples = 1000;
  int32_t         keyidx;
  int             s;

  for (s = 0; s < nsamples; s++)
    {
      p7_anchors_Sample(rng, L, M, maxD, anch);
      p7_anchorhash_Store(ah, anch, &keyidx);
      p7_anchors_Reuse(anch);
    }
  
  p7_anchorhash_Dump(stdout, ah);

  p7_anchorhash_Destroy(ah);
  p7_anchors_Destroy(anch);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7ANCHORHASH_EXAMPLE*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

