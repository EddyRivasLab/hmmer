/* P7_TOPHITS: ranking lists of top-scoring hits
 */
#ifndef p7TOPHITS_INCLUDED
#define p7TOPHITS_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_keyhash.h"
#include "esl_msa.h"
#include "esl_sq.h"

#include "base/p7_domain.h"
#include "base/p7_alidisplay.h"

#define p7_HITFLAGS_DEFAULT 0
#define p7_IS_INCLUDED      (1<<0)
#define p7_IS_REPORTED      (1<<1)
#define p7_IS_NEW           (1<<2)
#define p7_IS_DROPPED       (1<<3)
#define p7_IS_DUPLICATE     (1<<4)


/* Structure: P7_HIT
 * 
 * Info about a high-scoring database hit, kept so we can output a
 * sorted list of high hits at the end.
 *
 * sqfrom and sqto are the coordinates that will be shown in the
 * results, not coords in arrays... therefore, reverse complements
 * have sqfrom > sqto
 */
typedef struct p7_hit_s {
  char   *name;		/* name of the target               (mandatory)           */
  char   *acc;		/* accession of the target          (optional; else NULL) */
  char   *desc;		/* description of the target        (optional; else NULL) */
  int    window_length; /* used in e-value computation, when splitting long seqs  */
  double sortkey;	/* number to sort by; big is better                       */

  float  score;		/* bit score of the sequence (all domains, w/ correction) */
  float  pre_score;	/* bit score of sequence before null2 correction          */
  float  sum_score;	/* bit score reconstructed from sum of domain envelopes   */

  double lnP;		/* log(P-value) of the score               */
  double pre_lnP;	/* log(P-value) of the pre_score           */
  double sum_lnP;	/* log(P-value) of the sum_score           */

  int    ndom;		/* total # of domains identified in this seq              */
  int    noverlaps;	/* # of domain envelopes that overlap with a previous one */
  float  nexpected;	/* expected # of domains, by posterior decoding           */

  uint32_t flags;      	/* p7_IS_REPORTED | p7_IS_INCLUDED | p7_IS_NEW | p7_IS_DROPPED */
  int      nreported;	/* # of domains satisfying reporting thresholding              */
  int      nincluded;	/* # of domains satisfying inclusion thresholding              */
  int      best_domain;	/* index of best-scoring domain in dcl                         */

  int64_t  seqidx;       /* unique identifier to track the database sequence from which this hit came*/
  int64_t  subseq_start; /* used to track which subseq of full len target this hit came from, for purposes of removing duplicates */

  P7_DOMAIN *dcl;	/* domain coordinate list and alignment display */
  esl_pos_t  offset;	/* used in socket communications, in serialized communication: offset of P7_DOMAIN msg for this P7_HIT */
} P7_HIT;


/* Structure: P7_TOPHITS
 * merging when we prepare to output results. "hit" list is NULL and
 * unavailable until after we do a sort.  
 */
typedef struct p7_tophits_s {
  P7_HIT **hit;                  /* sorted pointer array                     */
  P7_HIT  *unsrt;                /* unsorted data storage                    */
  int      Nalloc;	         /* current allocation size                  */
  int      N;		         /* number of hits in list now               */
  int      nreported;	         /* number of hits that are reportable       */
  int      nincluded;	         /* number of hits that are includable       */
  int      is_sorted_by_sortkey; /* TRUE when hits sorted by sortkey and th->hit valid for all N hits */
  int      is_sorted_by_seqidx;  /* TRUE when hits sorted by seq_idx, position, and th->hit valid for all N hits */
} P7_TOPHITS;

#define p7_TOPHITS_DEFAULT_INIT_ALLOC 100


/* 1. The P7_TOPHITS object */
extern P7_TOPHITS *p7_tophits_Create(int init_hit_alloc);
extern int         p7_tophits_Grow(P7_TOPHITS *th);
extern int         p7_tophits_CreateNextHit(P7_TOPHITS *th, P7_HIT **ret_hit);
extern int         p7_tophits_SortBySortkey(P7_TOPHITS *th);
extern int         p7_tophits_SortBySeqidxAndAlipos(P7_TOPHITS *th);
extern int         p7_tophits_SortByModelnameAndAlipos(P7_TOPHITS *th);
extern int         p7_tophits_Merge(P7_TOPHITS *th1, P7_TOPHITS *th2);
extern int         p7_tophits_GetMaxPositionLength(P7_TOPHITS *th);
extern int         p7_tophits_GetMaxNameLength(P7_TOPHITS *th);
extern int         p7_tophits_GetMaxAccessionLength(P7_TOPHITS *th);
extern int         p7_tophits_GetMaxShownLength(P7_TOPHITS *th);
extern int         p7_tophits_Reuse(P7_TOPHITS *th);
extern void        p7_tophits_Destroy(P7_TOPHITS *th);

/* 2. The P7_HIT object array in P7_TOPHITS */
extern P7_HIT *p7_hit_Create(int nhit_alloc);
extern int     p7_hit_Grow(P7_HIT **hitp, int oldalloc, int newalloc);
extern void    p7_hit_Destroy(P7_HIT *hits, int nhits);

/* 3. Debugging and development tools */
extern int p7_tophits_TestSample(ESL_RANDOMNESS *rng, P7_TOPHITS **ret_th);
extern int p7_tophits_Validate(const P7_TOPHITS *th,  char *errbuf);
extern int p7_tophits_Compare (const P7_TOPHITS *th1, const P7_TOPHITS *th2, float tol);
extern int p7_hit_TestSample(ESL_RANDOMNESS *rng, P7_HIT *hit);
extern int p7_hit_Validate(const P7_HIT *hit, char *errbuf);
extern int p7_hit_Compare(const P7_HIT *h1, const P7_HIT *h2, float tol);


#endif /*p7TOPHITS_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
