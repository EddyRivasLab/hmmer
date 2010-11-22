#ifndef P7_CACHEDB_INCLUDED
#define P7_CACHEDB_INCLUDED

typedef struct {
  char    *name;                   /* name; ("\0" if no name)               */
  ESL_DSQ *dsq;                    /* digitized sequence [1..n]             */
  int64_t  n;                      /* length of dsq                         */
  int64_t  idx;	                   /* ctr for this seq                      */
  uint64_t db_key;                 /* flag for included databases           */
} HMMER_SEQ;

typedef struct {
  uint32_t            count;       /* number of entries                     */
  uint32_t            K;           /* original number of entries            */
  HMMER_SEQ         **list;        /* list of sequences [0 .. count-1]      */
} SEQ_DB;

typedef struct {
  char               *name;        /* name of the seq database              */
  char               *id;          /* unique identifier string              */
  uint32_t            db_cnt;      /* number of sub databases               */
  SEQ_DB             *db;          /* list of databases [0 .. db_cnt-1]     */

  ESL_ALPHABET       *abc;         /* alphabet for database                 */

  uint32_t            count;       /* total number of sequences             */
  HMMER_SEQ          *list;        /* complete list of sequences (count)    */
  void               *residue_mem; /* memory holding the residues           */
  char               *header_mem;  /* memory holding the header strings     */

  uint64_t            res_size;    /* size of residue memory allocation     */
  uint64_t            hdr_size;    /* size of header memory allocation      */
} SEQ_CACHE;

typedef struct {
  char               *name;        /* name of the hmm database              */
  char               *id;          /* unique identifier string              */
  ESL_ALPHABET       *abc;         /* alphabet for database                 */

  uint32_t            count;       /* number of entries                     */
  P7_OPROFILE       **list;        /* list of profiles [0 .. count-1]       */
} HMM_CACHE;


int  cache_HmmDb(char *hmmfile, HMM_CACHE **ret_cache);
void cache_HmmDestroy(HMM_CACHE *cache);

int  cache_SeqDb(char *seqfile, SEQ_CACHE **ret_cache);
void cache_SeqDestroy(SEQ_CACHE *cache);

#endif /*P7_CACHEDB_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
