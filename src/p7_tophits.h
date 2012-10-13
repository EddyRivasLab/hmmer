/* P7_TOPHITS: ranking lists of top-scoring hits
 */
#ifndef P7_TOPHITS_INCLUDED
#define P7_TOPHITS_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_keyhash.h"
#include "esl_msa.h"
#include "esl_sq.h"

#include "base/p7_domain.h"
#include "base/p7_alidisplay.h"
#include "p7_pipeline.h"

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
 * sqfrom and sqto are the coordinates that will be shown
 * in the results, not coords in arrays... therefore, reverse
 * complements have sqfrom > sqto
 */
typedef struct p7_hit_s {
  char   *name;			/* name of the target               (mandatory)           */
  char   *acc;			/* accession of the target          (optional; else NULL) */
  char   *desc;			/* description of the target        (optional; else NULL) */
  int    window_length;         /* for later use in e-value computation, when splitting long sequences */
  double sortkey;		/* number to sort by; big is better                       */

  float  score;			/* bit score of the sequence (all domains, w/ correction) */
  float  pre_score;		/* bit score of sequence before null2 correction          */
  float  sum_score;		/* bit score reconstructed from sum of domain envelopes   */

  double lnP;		        /* log(P-value) of the score               */
  double pre_lnP;		/* log(P-value) of the pre_score           */
  double sum_lnP;		/* log(P-value) of the sum_score           */

  float  nexpected;     /* posterior expected number of domains in the sequence (from posterior arrays) */
  int    nregions;	/* number of regions evaluated */
  int    nclustered;	/* number of regions evaluated by clustering ensemble of tracebacks */
  int    noverlaps;	/* number of envelopes defined in ensemble clustering that overlap w/ prev envelope */
  int    nenvelopes;	/* number of envelopes handed over for domain definition, null2, alignment, and scoring. */
  int    ndom;		/* total # of domains identified in this seq   */

  uint32_t flags;      	/* p7_IS_REPORTED | p7_IS_INCLUDED | p7_IS_NEW | p7_IS_DROPPED */
  int      nreported;	/* # of domains satisfying reporting thresholding  */
  int      nincluded;	/* # of domains satisfying inclusion thresholding */
  int      best_domain;	/* index of best-scoring domain in dcl */

  int64_t  seqidx;          /*unique identifier to track the database sequence from which this hit came*/
  int64_t  subseq_start; /*used to track which subsequence of a full_length target this hit came from, for purposes of removing duplicates */

  P7_DOMAIN *dcl;	/* domain coordinate list and alignment display */
  esl_pos_t  offset;	/* used in socket communications, in serialized communication: offset of P7_DOMAIN msg for this P7_HIT */
} P7_HIT;


/* Structure: P7_TOPHITS
 * merging when we prepare to output results. "hit" list is NULL and
 * unavailable until after we do a sort.  
 */
typedef struct p7_tophits_s {
  P7_HIT **hit;         /* sorted pointer array                     */
  P7_HIT  *unsrt;	/* unsorted data storage                    */
  uint64_t Nalloc;	/* current allocation size                  */
  uint64_t N;		/* number of hits in list now               */
  uint64_t nreported;	/* number of hits that are reportable       */
  uint64_t nincluded;	/* number of hits that are includable       */
  int      is_sorted_by_sortkey; /* TRUE when hits sorted by sortkey and th->hit valid for all N hits */
  int      is_sorted_by_seqidx; /* TRUE when hits sorted by seq_idx, position, and th->hit valid for all N hits */
} P7_TOPHITS;



extern P7_TOPHITS *p7_tophits_Create(void);
extern int         p7_tophits_Grow(P7_TOPHITS *h);
extern int         p7_tophits_CreateNextHit(P7_TOPHITS *h, P7_HIT **ret_hit);
extern int         p7_tophits_Add(P7_TOPHITS *h,
				  char *name, char *acc, char *desc, 
				  double sortkey, 
				  float score,    double lnP, 
				  float mothersc, double mother_lnP,
				  int sqfrom, int sqto, int sqlen,
				  int hmmfrom, int hmmto, int hmmlen, 
				  int domidx, int ndom,
				  P7_ALIDISPLAY *ali);
extern int         p7_tophits_SortBySortkey(P7_TOPHITS *h);
extern int         p7_tophits_SortBySeqidxAndAlipos(P7_TOPHITS *h);
extern int         p7_tophits_SortByModelnameAndAlipos(P7_TOPHITS *h);

extern int         p7_tophits_Merge(P7_TOPHITS *h1, P7_TOPHITS *h2);
extern int         p7_tophits_GetMaxPositionLength(P7_TOPHITS *h);
extern int         p7_tophits_GetMaxNameLength(P7_TOPHITS *h);
extern int         p7_tophits_GetMaxAccessionLength(P7_TOPHITS *h);
extern int         p7_tophits_GetMaxShownLength(P7_TOPHITS *h);
extern void        p7_tophits_Destroy(P7_TOPHITS *h);

extern int p7_tophits_ComputeNhmmerEvalues(P7_TOPHITS *th, double N, int W);
extern int p7_tophits_RemoveDuplicates(P7_TOPHITS *th, int using_bit_cutoffs);
extern int p7_tophits_Threshold(P7_TOPHITS *th, P7_PIPELINE *pli);
extern int p7_tophits_CompareRanking(P7_TOPHITS *th, ESL_KEYHASH *kh, int *opt_nnew);
extern int p7_tophits_Targets(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw);
extern int p7_tophits_Domains(FILE *ofp, P7_TOPHITS *th, P7_PIPELINE *pli, int textw);
extern int p7_tophits_Alignment(const P7_TOPHITS *th, const ESL_ALPHABET *abc, 
				ESL_SQ **inc_sqarr, P7_TRACE **inc_trarr, int inc_n, int optflags,
				ESL_MSA **ret_msa);
extern int p7_tophits_TabularTargets(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header);
extern int p7_tophits_TabularDomains(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header);
extern int p7_tophits_TabularXfam(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli);
extern int p7_tophits_TabularTail(FILE *ofp, const char *progname, enum p7_pipemodes_e pipemode, 
				  const char *qfile, const char *tfile, const ESL_GETOPTS *go);
extern int p7_tophits_LongInserts(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int min_length);
extern int p7_tophits_AliScores(FILE *ofp, char *qname, P7_TOPHITS *th );

#endif /*P7_TOPHITS_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
