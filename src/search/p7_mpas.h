/* Optional parameters for controlling, and statistics for monitoring,
 * the most probable anchor set (MPAS) algorithm implementations --
 * both reference and production versions.
 * 
 */
#ifndef p7MPAS_INCLUDED
#define p7MPAS_INCLUDED
#include "p7_config.h"

#include <stdio.h>

#include "base/p7_anchors.h"
#include "base/p7_trace.h"


typedef struct {
  int     max_iterations;   // limit on # of sampled paths in MPAS (typically ~1000)
  float   loss_threshold;   // max acceptable probability of missing better solution (typical: ~0.001)
  int     be_verbose;       // TRUE to dump additional debug/devel info as it runs
} P7_MPAS_PARAMS;


typedef struct {
  int   has_part1;	        /* TRUE once we set the stuff below:             */
  int   tot_iterations;		/* # of samples we took                          */
  int   tot_asc_calculations;	/* # of non-dup samples, where we did ASC        */
  float vsc;			/* Viterbi score (raw; nats)                     */
  float fsc;			/* Forward score (raw; nats)                     */
  float vit_asc;		/* ASC score of Viterbi annotation               */
  float vit_ascprob;		/* probability of the Viterbi annotation         */
  float best_asc;		/* anchor-constrained score of solution          */
  float best_ascprob;		/* probability of the asc: exp(best_asc - fwdsc) */
  int   nsamples_in_best;	/* number of samples we saw for the solution     */
  int   best_is_viterbi;	/* TRUE if solution is Viterbi annotation        */
  int   solution_not_found;	/* TRUE if we find no solution                   */

  int   has_part2;		/* TRUE once we set the remaining stuff below:   */
  int   anch_outside;		/* # of anchors that fall outside any domain in Viterbi trace */
  int   anch_unique;		/* # of anchors that map 1:1 to a domain in Viterbi trace       */
  int   anch_multiple;		/* # of anchors that fall in same Viterbi domain with other anchor(s) */
  int   dom_zero;		/* # of Viterbi domains with no anchors */
  int   dom_one;		/* # of Viterbi domains that map 1:1 to anchor ( == anch_unique) */
  int   dom_multiple;		/* # of Viterbi domains with >1 anchor  */
} P7_MPAS_STATS;


extern int p7_mpas_stats_Init(P7_MPAS_STATS *stats);
extern int p7_mpas_stats_Dump(FILE *ofp, P7_MPAS_STATS *stats);
extern int p7_mpas_stats_CompareAS2Trace(P7_MPAS_STATS *stats, const P7_ANCHORS *anch, const P7_TRACE *tr);

#endif /*p7MPAS_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
