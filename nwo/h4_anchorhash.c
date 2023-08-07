/* H4_ANCHORHASH implements an auxiliary data structure that the MPAS
 * (most probable anchor set) algorithm uses to store anchor sets it's
 * already seen.
 * 
 * Needs to deal with both complete anchor sets, and suffixes of
 * partial anchor sets. The segmental divide & conquer version of the
 * MPAS algorithm solves the optimal anchor set sequentially, one
 * segment at a time.
 * 
 * Contents:
 *   1. H4_ANCHORHASH object
 *   2. Internal (static) functions
 *   3. Debugging and development tools
 *   4. Unit tests
 *   5. Test driver
 */
#include <h4_config.h>

#include <stdlib.h>
#include <limits.h>  // for INT32_MAX

#include "easel.h"

#include "h4_anchorset.h"
#include "h4_anchorhash.h"

static uint32_t anchorhash_function    (const H4_ANCHOR *arr, int32_t D, uint32_t hashsize);
static uint32_t anchorhash_function_alt(int32_t *keydata, int32_t hashsize);
static int      anchorhash_compare     (const H4_ANCHOR *arr, int32_t D, int32_t *keydata);
static int      anchorhash_upsize      (H4_ANCHORHASH *ah);


/*****************************************************************
 * 1. The H4_ANCHORHASH object
 *****************************************************************/

/* Function:  h4_anchorhash_Create()
 * Synopsis:  Create a <H4_ANCHORHASH>
 *
 * Purpose:   Allocate and initialize a <H4_ANCHORHASH> hash table for storing
 *            lots of anchor arrays.
 *
 * Returns:   pointer to the new <H4_ANCHORHASH> object on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
H4_ANCHORHASH *
h4_anchorhash_Create(void)
{
  H4_ANCHORHASH *ah                = NULL;
  int32_t        default_hashsize  = 128;  // Size of the hash table. Must be a power of 2; grows by doubling.
  int32_t        default_nkeyalloc = 128;  // Number of keys (anchor sets)
  int32_t        default_aalloc    = 2048; // Total number of integers in data elements (2D+1 for each anchor set of D anchors)
  int32_t        i;
  int            status;

  ESL_ALLOC(ah, sizeof(H4_ANCHORHASH));
  ah->hashtable  = NULL;
  ah->key_offset = NULL;
  ah->key_count  = NULL;
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
  ESL_ALLOC(ah->key_count,  sizeof(int32_t) * ah->kalloc);
  ESL_ALLOC(ah->nxt,        sizeof(int32_t) * ah->kalloc);
  ESL_ALLOC(ah->amem,       sizeof(int32_t) * ah->aalloc);
  return ah;
  
 ERROR:
  h4_anchorhash_Destroy(ah);
  return NULL;
}

/* Function:  h4_anchorhash_Sizeof()
 * Synopsis:  Returns the allocated size (in bytes) of an H4_ANCHORHASH
 * Incept:    SRE, Thu 18 Feb 2021
 */
size_t
h4_anchorhash_Sizeof(const H4_ANCHORHASH *ah)
{
  size_t n = 0;

  n += sizeof(H4_ANCHORHASH);
  n += sizeof(int32_t) * ah->hashsize;	 // hashtable 
  n += sizeof(int32_t) * ah->kalloc * 3; // key_offset, key_count, nxt 
  n += sizeof(int32_t) * ah->aalloc;	 // amem 
  return n;
}


/* Function:  h4_anchorhash_Reuse()
 * Synopsis:  Reuse a <H4_ANCHORHASH>
 *
 * Purpose:   Clear a <H4_ANCHORHASH> hash table for reuse.
 *
 *            If any allocations are overly large, drop them back to 'redline'
 *            values. Default redlines are 1024 keys (i.e. different anchor sets)
 *            1024 hash values, and 16384 total integers of raw data. These redlines
 *            are all 8x the default initial allocations.
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
h4_anchorhash_Reuse(H4_ANCHORHASH *ah)
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
      ESL_REALLOC(ah->key_count,  sizeof(int32_t) * kalloc_redline);
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


/* Function:  h4_anchorhash_Destroy()
 * Synopsis:  Destroys a <H4_ANCHORHASH> hash table.
 */
void
h4_anchorhash_Destroy(H4_ANCHORHASH *ah)
{
  if (ah)
    {
      if (ah->hashtable)  free(ah->hashtable);
      if (ah->key_offset) free(ah->key_offset);
      if (ah->key_count)  free(ah->key_count);
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
 *            <D0> allows us to store suffixes, supporting the
 *            segmental divide and conquer version of the MPAS
 *            algorithm. Do not store the first <D0> anchors; only
 *            store <D0+1..D>. To store the complete anchor set, pass
 *            <D0=0>.
 *
 *            If an identical anchor set is already stored in <ah>,
 *            set <*opt_index> to the key for that anchor set, and
 *            return <eslEDUP>.
 *            
 *            Increment <ah->key_count[]> counter every time we call
 *            <_Store()> on a given anchorset suffix (not counting
 *            D0). This collects the observed frequency of sampling
 *            the anchorset suffix, which we can compare to its
 *            calculated probability.
 *
 * Args:      ah         : hash table holding different anchor sets
 *            anch       : new anchor set to try to store
 *            D0         : ignore first <D0> anchors, store <D0+1..D> (0 = store all)
 *            opt_index  : optRETURN: index of stored data
 *            
 * Returns:   <eslOK> if <anch> is new; the anchor set data are stored, 
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
h4_anchorhash_Store(H4_ANCHORHASH *ah, const H4_ANCHORSET *anch, int D0, int32_t *opt_index)
{
  uint32_t  val = anchorhash_function(anch->a + D0, anch->D - D0, ah->hashsize);
  int32_t  *ptr;
  int32_t   idx;
  int32_t   d;
  int       status;
  
  /* Was this key already stored? */
  for (idx = ah->hashtable[val]; idx != -1; idx = ah->nxt[idx])
    {
      if (anchorhash_compare(anch->a + D0, anch->D - D0, ah->amem + ah->key_offset[idx]) == eslOK)
	{
	  ah->key_count[idx]++;
	  if (opt_index) *opt_index = idx;
	  return eslEDUP;
	}
    }

  /* Reallocate key memory if needed */
  if (ah->nkeys == ah->kalloc)
    {
      ESL_REALLOC(ah->key_offset, sizeof(int32_t) * ah->kalloc * 2);
      ESL_REALLOC(ah->key_count,  sizeof(int32_t) * ah->kalloc * 2);
      ESL_REALLOC(ah->nxt,        sizeof(int32_t) * ah->kalloc * 2);
      ah->kalloc *= 2;
    }

  /* Reallocate key data memory if needed (by doubling) */
  while (ah->an + 2 * (anch->D - D0) + 1 > ah->aalloc)
    {
      ESL_REALLOC(ah->amem, sizeof(int32_t) * ah->aalloc * 2);
      ah->aalloc *= 2;
    }

  /* Copy the key, assign its index */
  idx                 = ah->nkeys;
  ah->key_offset[idx] = ah->an;
  ah->key_count[idx]  = 1;                      // Not ++. This is an initialization.
  ah->an             += 2 * (anch->D - D0) + 1;
  ah->nkeys++;

  ptr  = ah->amem + ah->key_offset[idx];
  *ptr = anch->D - D0;
  for (d = D0 + 1; d <= anch->D; d++) 
    {
      ptr++; *ptr = anch->a[d].i0;
      ptr++; *ptr = anch->a[d].k0;
    }
  
  /* anchorhash needs to remember L,M so when caller asks
   * to _Get() an anchor set, anchorhash can set the sentinels
   * correctly. Fortunately even when we're only storing a 
   * suffix of <anch>, we still get the whole <anch> object,
   * which has valid sentinels, so we can deduce from them
   * what L,M are.
   */
  if (ah->nkeys == 1) 
    h4_anchorset_GetSentinels(anch, &(ah->L), &(ah->M));

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

/* Function:  h4_anchorhash_Get()
 * Synopsis:  Get an anchor set back from the hash.
 *
 * Purpose:   Retrieve anchor set <keyidx> from hash <ah>, and put
 *            it in <anch>. <anch> is reallocated if necessary.
 *            
 *            <D0> supports the segmental divide and conquer version
 *            of the MPAS algorithm, where we're storing and
 *            retrieving suffixes of a growing anchor set, rather than
 *            complete anchor sets. <D0> says to keep the first <D0>
 *            anchors in a growing <anch>, and appending the retrieved
 *            anchorset starting at <D0+1>. To get a complete anchor
 *            set, pass <D0=0>.
 *            
 * Args:      ah      : hash storage for alternative annotations
 *            keyidx  : which annotation to get [0..ah->nkeys-1]
 *            D0      : keep 1..D0 in <anch>; append starting at D0+1; 0=get complete anchor set
 *            anch    : RETURN: anchor set
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on failure, and the state of <anch> is undefined.
 */
int
h4_anchorhash_Get(const H4_ANCHORHASH *ah, int32_t keyidx, int D0, H4_ANCHORSET *anch)
{
  int32_t *ptr  = ah->amem + ah->key_offset[keyidx];
  int32_t  D    = *ptr + D0;                          // deal w/ case where *ptr is Dg, a suffix, to be added to D0
  int32_t  d;
  int      status;

  if ((status = h4_anchorset_GrowFor(anch, D)) != eslOK) goto ERROR;

  /* Get the data and append it to <anch> */
  for (d = D0+1; d <= D; d++)
    {
      ptr++; anch->a[d].i0 = *ptr;
      ptr++; anch->a[d].k0 = *ptr;
    }

  anch->D = D;
  h4_anchorset_SetSentinels(anch, ah->L, ah->M);
  return eslOK;

 ERROR:
  return status;
}
  
/*****************************************************************
 * 2. Internal (static) functions
 *****************************************************************/


/* anchorhash_function()
 *   
 * Given <arr>/<n> data, and the current <hashsize>;
 * calculate and return a hash function on that data, in range <0..hashsize-1>.
 * 
 * Does not depend on, nor access, sentinels; so may be safely called on
 * a subsequence of an anchor set. For example, when segmental divide and
 * conquer MPAS stores/retrieves a suffix D0+1..D of an anchor set,
 * we pass <arr = anch->a+D0> and <D = anch->D-D0>.
 */
static uint32_t
anchorhash_function(const H4_ANCHOR *arr, int32_t D, uint32_t hashsize)
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
 * The data are stored in a serialized array of 2D+1 integers:
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
 * Compare a <arr>/<D> anchor set (a H4_ANCHOR array) to 
 * data from a H4_ANCHOR array that's already been stored,
 * starting at <keydata>.
 *
 * When comparing a suffix of an anchor set, <arr> is <anch->a + D0>,
 * and D is really D-D0,
 * and we start accessing at <anch->a + D0 + 1>.  This internal access
 * of a subsequence of <anch->a> array is fine because
 * <anchorhash_compare> does not depend on the sentinel values at
 * <0,D+1>.
 * 
 * Return <eslOK> if the two are identical; 
 * return <eslFAIL> if not.
 */
static int
anchorhash_compare(const H4_ANCHOR *arr, int D, int32_t *keydata)
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
anchorhash_upsize(H4_ANCHORHASH *ah)
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
 * 3. Debugging and development tools
 *****************************************************************/

/* Function:  h4_anchorhash_Dump()
 * Synopsis:  Dump contents of an H4_ANCHORHASH to a stream.
 * Incept:    SRE, Fri 19 Feb 2021
 */
int
h4_anchorhash_Dump(FILE *ofp, const H4_ANCHORHASH *ah)
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
  fprintf(ofp, "Total obj size, bytes:  %d\n", (int) h4_anchorhash_Sizeof(ah));

  return eslOK;
}


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef h4ANCHORHASH_TESTDRIVE

#include "esl_random.h"

static void
utest_sampling(ESL_RANDOMNESS *rng)
{
  char           msg[]    = "h4_anchorhash :: sampling unit test failed";
  int            nsamples = 1000;
  int            M        = 100;
  int            L        = 200;
  int            maxD     = 50;  
  H4_ANCHORHASH *ah       = h4_anchorhash_Create();
  H4_ANCHORSET  *anch     = h4_anchorset_Create( (esl_rnd_Roll(rng, 2) ? maxD : 0), L, M);  // test reallocation sometimes
  H4_ANCHORSET **aarr     = malloc(sizeof(H4_ANCHORSET *) * nsamples);
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
	aarr[s] = h4_anchorset_Create( (esl_rnd_Roll(rng, 2) ? maxD : 0), L, M);

      nk = 0;
      for (s = 0; s < nsamples; s++)
	{
	  if ( h4_anchorset_Sample(rng, L, M, maxD, anch) != eslOK) esl_fatal(msg);
	  if ( h4_anchorset_Copy(anch, aarr[s])           != eslOK) esl_fatal(msg);

	  status = h4_anchorhash_Store(ah, anch, 0, &keyidx);
	  keys[s] = keyidx;
	  if      (status == eslOK)   { if (keyidx != nk) esl_fatal(msg); nk++; }
	  else if (status == eslEDUP) { if (keyidx >= nk) esl_fatal(msg);       }
	  else                        {                   esl_fatal(msg);       }

	  h4_anchorset_Reuse(anch);
	}

      for (s = 0; s < nsamples; s++)
	{
	  if ( h4_anchorhash_Get(ah, keys[s], 0, anch) != eslOK) esl_fatal(msg);

	  //p7_anchors_Dump(stdout, anch);

	  if ( h4_anchorset_Validate(anch, errmsg) != eslOK) esl_fatal("%s:\n  %s", msg, errmsg);
	  if ( h4_anchorset_Compare(anch, aarr[s]) != eslOK) esl_fatal(msg);
	  if ( h4_anchorset_Reuse(anch)            != eslOK) esl_fatal(msg);
	}

      for (s = 0; s < nsamples; s++)
	h4_anchorset_Destroy(aarr[s]);

      h4_anchorhash_Reuse(ah);  
    }

  h4_anchorhash_Destroy(ah);
  h4_anchorset_Destroy(anch);
  free(aarr);
  free(keys);
}
#endif /*h4ANCHORHASH_TESTDRIVE*/

/*****************************************************************
 * 5. Test driver
 *****************************************************************/

#ifdef h4ANCHORHASH_TESTDRIVE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",         0 },
  { "--version", eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for h4_anchorhash";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = h4_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_sampling(rng);

  fprintf(stderr, "#  status = ok\n");

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  exit(0); 
}
#endif //h4ANCHORHASH_TESTDRIVE
/*-------------------- end of test driver ---------------------*/


/*****************************************************************
 * 6. Example driver
 *****************************************************************/
#ifdef h4ANCHORHASH_EXAMPLE

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles   reqs   incomp   help                                    docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL,  NULL,   "show brief help on version and usage",        0 },
  { "-s",           eslARG_INT,      "0", NULL, NULL,      NULL,  NULL,  NULL,   "set random number seed to <n>",               0 },
  { "--version",    eslARG_NONE,    NULL, NULL, NULL,      NULL,  NULL,  NULL,   "show HMMER version number",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "example driver for H4_ANCHORHASH";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go       = h4_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng      = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  H4_ANCHORSET   *anch     = NULL;
  H4_ANCHORHASH  *ah       = h4_anchorhash_Create();
  int             L        = 400;
  int             M        = 200;
  int             maxD     = 50;
  int             nsamples = 1000;
  int32_t         keyidx;
  int             s;

  anch = h4_anchorset_Create(maxD, L, M);
  for (s = 0; s < nsamples; s++)
    {
      h4_anchorset_Sample(rng, L, M, maxD, anch);
      h4_anchorhash_Store(ah, anch, 0, &keyidx);
      h4_anchorset_Reuse(anch);
    }
  
  h4_anchorhash_Dump(stdout, ah);

  h4_anchorhash_Destroy(ah);
  h4_anchorset_Destroy(anch);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*h4ANCHORHASH_EXAMPLE*/



