/* A cached profile database. Used by the hmmpgmd daemon.
 */
#ifndef P7_HMMCACHE_INCLUDED
#define P7_HMMCACHE_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"

#include "base/p7_profile.h"
#include "dp_vector/p7_oprofile.h"

typedef struct {
  char               *name;        /* name of the hmm database              */
  ESL_ALPHABET       *abc;         /* alphabet for database                 */

  P7_OPROFILE       **omlist;      /* list of vectorized profiles [0..n-1]  */
  P7_PROFILE        **gmlist;      /* list of standard profiles [0..n-1]    */

  uint32_t            lalloc;	   /* allocated length of <list>            */
  uint32_t            n;           /* number of entries in <list>           */
} P7_HMMCACHE;

extern int    p7_hmmcache_Open (char *hmmfile, P7_HMMCACHE **ret_cache, char *errbuf);
extern size_t p7_hmmcache_Sizeof         (P7_HMMCACHE *cache);
extern int    p7_hmmcache_SetNumericNames(P7_HMMCACHE *cache);
extern void   p7_hmmcache_Close          (P7_HMMCACHE *cache);

#endif /*P7_HMMCACHE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
