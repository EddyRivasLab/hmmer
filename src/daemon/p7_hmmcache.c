/* A cached profile database. Used by the hmmpgmd daemon.
 * 
 * Contents:
 *   1. P7_HMMCACHE : a daemon's cached profile database.
 *   2. Benchmark driver.
 *   3. License and copyright information
 */
#include "p7_config.h"

#include <stdlib.h>
#include <string.h>

#include "easel.h"

#include "base/p7_hmmfile.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/io.h"

#include "daemon/p7_hmmcache.h"

#include "search/modelconfig.h"

/*****************************************************************
 * 1. P7_HMMCACHE: a daemon's cached profile database
 *****************************************************************/ 

/* Function:  p7_hmmcache_Open()
 * Synopsis:  Cache a profile database.
 *
 * Purpose:   Open <hmmfile> and read all of its contents, creating
 *            a cached profile database in memory. Return a ptr to the 
 *            cached profile database in <*ret_cache>. 
 *            
 *            Caller may optionally provide an <errbuf> ptr to
 *            at least <eslERRBUFSIZE> bytes, to capture an 
 *            informative error message on failure. 
 *            
 * Args:      hmmfile   - (base) name of profile file to open
 *            ret_cache - RETURN: cached profile database
 *            errbuf    - optRETURN: error message for a failure
 *
 * Returns:   <eslOK> on success. <*ret_cache> points to the
 *            cached db. <errbuf> is unchanged.
 *            
 *            Failure codes:
 *            <eslENOTFOUND> : <hmmfile> couldn't be opened for reading
 *            <eslEFORMAT>   : <hmmfile> isn't in recognized HMMER file format
 *            <eslEINCOMPAT> : profiles in <hmmfile> have different alphabets
 *
 *            On any failure, <*ret_cache> is <NULL> and <errbuf> contains
 *            an informative error message for the user.
 *
 * Throws:    <eslEMEM> : memory allocation error.
 */
int
p7_hmmcache_Open(char *hmmfile, P7_HMMCACHE **ret_cache, char *errbuf)
{
  P7_HMMCACHE *cache    = NULL;
  P7_HMMFILE  *hfp      = NULL;     
  P7_HMM      *hmm      = NULL;
  P7_BG       *bg       = NULL;
  P7_PROFILE  *gm       = NULL;
  P7_OPROFILE *om       = NULL;     
  int          status;
  P7_HARDWARE *hw;
  if ((hw = p7_hardware_Create ()) == NULL)  p7_Fail("Couldn't get HW information data structure"); 
  if (errbuf) errbuf[0] = '\0';

  ESL_ALLOC(cache, sizeof(P7_HMMCACHE));
  cache->name      = NULL;
  cache->abc       = NULL;
  cache->omlist    = NULL;
  cache->gmlist    = NULL;
  cache->lalloc    = 4096;	/* allocation chunk size for <list> of ptrs  */
  cache->n         = 0;

  if ( ( status = esl_strdup(hmmfile, -1, &cache->name) != eslOK)) goto ERROR; 
  ESL_ALLOC(cache->omlist, sizeof(P7_OPROFILE *) * cache->lalloc);
  ESL_ALLOC(cache->gmlist, sizeof(P7_PROFILE *)  * cache->lalloc);

  if ( (status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf)) != eslOK) goto ERROR;  // eslENOTFOUND | eslEFORMAT; <errbuf> 

  while ((status = p7_hmmfile_Read(hfp, &(cache->abc), &hmm)) != eslEOF)  // eslEFORMAT | eslEINCOMPAT; <errbuf>
    {
      if (status != eslOK) ESL_XFAIL(status, errbuf, "%s", hfp->errbuf); 

      if (!bg && (bg = p7_bg_Create(cache->abc)) == NULL)  { status = eslEMEM; goto ERROR; }

      if ( (    gm = p7_profile_Create(hmm->M, cache->abc)) == NULL)  { status = eslEMEM; goto ERROR; }
      if ( (status = p7_profile_Config(gm, hmm, bg)) != eslOK) goto ERROR;
 
      if ( (status = p7_oprofile_ReadMSV (hfp, &(cache->abc), &om, hw->simd)) != eslOK || /* eslEFORMAT: hfp->errbuf | eslEINCOMPAT | eslEOF */
	   (status = p7_oprofile_ReadRest(hfp, om))                 != eslOK)   /* eslEFORMAT: hfp->errbuf */
	{
	  if (status == eslEOF) ESL_XFAIL(eslEFORMAT, errbuf, "Premature EOF in vectorized profile files");
	  else                  goto ERROR;
	}

      ESL_DASSERT1(( strcmp(gm->name, om->name) == 0 ));

      if (cache->n >= cache->lalloc) {
	ESL_REALLOC(cache->gmlist, sizeof(P7_PROFILE  *) * cache->lalloc * 2);
	ESL_REALLOC(cache->omlist, sizeof(P7_OPROFILE *) * cache->lalloc * 2);
	cache->lalloc *= 2;
      }

      cache->omlist[cache->n] = om;
      cache->gmlist[cache->n] = gm;
      cache->n++;

      om = NULL;
      gm = NULL;
      p7_hmm_Destroy(hmm);
    }

  //printf("\nfinal:: %d  memory %" PRId64 "\n", inx, total_mem);
  p7_hmmfile_Close(hfp);
  p7_bg_Destroy(bg);
  *ret_cache = cache;
  return eslOK;

 ERROR:
  if (cache) p7_hmmcache_Close(cache);
  if (om)    p7_oprofile_Destroy(om);
  if (gm)    p7_profile_Destroy(gm);
  if (hmm)   p7_hmm_Destroy(hmm);
  if (bg)    p7_bg_Destroy(bg);
  if (hfp)   p7_hmmfile_Close(hfp);
  return status;
}


/* Function:  p7_hmmcache_Sizeof()
 * Synopsis:  Returns total size of a profile cache, in bytes.
 * 
 * Purpose:   Calculate and return the size of a profile cache,
 *            in bytes. 
 *            
 *            The cache contains both standard and vectorized
 *            profiles. Very roughly, for a total number of consensus
 *            positions M, we consume M*560 bytes; 276 in a standard
 *            profile, and 284 in a vectorized one. This is a lot, and
 *            a target of future optimization.
 */
size_t
p7_hmmcache_Sizeof(P7_HMMCACHE *cache)
{
  size_t n = sizeof(P7_HMMCACHE);
  int    i;

  n += sizeof(char) * (strlen(cache->name) + 1);
  n += esl_alphabet_Sizeof(cache->abc);
  n += sizeof(P7_OPROFILE *) * cache->lalloc;     /* cache->omlist */
  n += sizeof(P7_PROFILE *)  * cache->lalloc;     /* cache->gmlist */

  for (i = 0; i < cache->n; i++) 
    {
      n += p7_oprofile_Sizeof(cache->omlist[i]);
      n += p7_profile_Sizeof (cache->gmlist[i]);
    }
  return n;
}
  

/* Function:  p7_hmmcache_SetNumericNames()
 * Synopsis:  Rename each profile in cache with a numeric name.
 *
 * Purpose:   Rename every profile in profile cache <cache> 
 *            with a numeric code, starting from "000000001".
 *
 *            The code is nine digits long, left padded with
 *            0's.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_hmmcache_SetNumericNames(P7_HMMCACHE *cache)
{
  int          namelen = 9;	/* 9 digit numeric code: 000000001, 000000002... */
  P7_PROFILE  *gm;
  P7_OPROFILE *om;
  int          i;
  int          status;

  for (i = 0; i < cache->n; i++)
    {
      gm = cache->gmlist[i];
      if (gm->name) free(gm->name);
      if (( status = esl_sprintf(&(gm->name), "%0*d", namelen, i+1)) != eslOK) return status;

      om = cache->omlist[i];
      if (om->name) free(om->name);
      if (( status = esl_sprintf(&(om->name), "%0*d", namelen, i+1)) != eslOK) return status;
    }
  return eslOK;
}


/* Function:  p7_hmmcache_Close()
 * Synopsis:  Free a profile cache.
 */
void
p7_hmmcache_Close(P7_HMMCACHE *cache)
{
  int i;

  if (cache)
    {
      if (cache->name) free(cache->name);
      if (cache->abc)  esl_alphabet_Destroy(cache->abc);
      for (i = 0; i < cache->n; i++)
	{
	  if (cache->gmlist) p7_profile_Destroy (cache->gmlist[i]);
	  if (cache->omlist) p7_oprofile_Destroy(cache->omlist[i]);
	}
      if (cache->gmlist) free(cache->gmlist);
      if (cache->omlist) free(cache->omlist);
      free(cache);
    }
}

/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/
#ifdef p7HMMCACHE_BENCHMARK

#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "daemon/p7_hmmcache.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <HMM file>";
static char banner[] = "benchmark driver for profile database cache";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_STOPWATCH *w       = esl_stopwatch_Create();
  char          *hmmfile = esl_opt_GetArg(go, 1);
  P7_HMMCACHE   *hcache  = NULL;
  char           errbuf[eslERRBUFSIZE];
  size_t         tot_mem;
  int            status;

  esl_stopwatch_Start(w);

  status = p7_hmmcache_Open(hmmfile, &hcache, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("Failed to read %s\n  %s\n",           hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("Failed to parse %s\n  %s\n",          hmmfile, errbuf);
  else if (status == eslEINCOMPAT) p7_Fail("Mixed profile types in %s\n  %s\n",   hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Failed to cache %s: error code %d\n", hmmfile, status);

  p7_hmmcache_SetNumericNames(hcache);
  tot_mem = p7_hmmcache_Sizeof(hcache);

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("models     = %d\n",          hcache->n);
  printf("tot memory = %" PRIu64 "\n", (uint64_t) tot_mem);
  
  p7_hmmcache_Close(hcache);
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return 0;
}
#endif /*p7HMMCACHE_BENCHMARK*/
/*--------------- end, benchmark driver -------------------------*/



/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
