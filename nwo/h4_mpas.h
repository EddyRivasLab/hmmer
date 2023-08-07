/* H4_MPAS_PARAMS, H4_MPAS_STATS: optional structures for setting control parameters,
 * and collecting performance statistics for the most probable anchor set (MPAS)
 * algorithm implementations -- both reference and production versions.
 * 
 */
#ifndef h4MPAS_INCLUDED
#define h4MPAS_INCLUDED
#include <h4_config.h>

#include <stdio.h>

#include "h4_anchorset.h"
#include "h4_path.h"


typedef struct {
  float   loss_threshold;   // max acceptable probability of missing better solution (typical: ~0.001)
  int     max_iterations;   // limit on # of sampled paths in MPAS (typically ~1000)
  int     nmax_sampling;    // TRUE to sample all n_max, to see if better sol'n is found late
  int     be_verbose;       // TRUE to dump additional debug/devel info as it runs
} H4_MPAS_PARAMS;


typedef struct {
  int   has_part1;	        // TRUE once we set the stuff below:
  int   tot_iterations;		// # of samples we took
  int   tot_asc_calculations;	// # of non-dup samples, where we did ASC
  float vsc;			// Viterbi score (raw; nats)
  float fsc;			// Forward score (raw; nats)
  float vit_asc;		// ASC score of Viterbi annotation
  float vit_ascprob;		// probability of the Viterbi annotation
  float best_asc;		// anchor-constrained score of solution
  float best_ascprob;		// probability of the asc: exp(best_asc - fwdsc)
  float tot_prob;		// \sum_A P(A | x,M) for all sampled anchorsets
  int   nsamples_in_best;	// number of samples we saw for the solution
  int   best_is_viterbi;	// TRUE if solution is Viterbi annotation
  int   late_solution;		// TRUE if better sol'n found after terminations
  int   solution_not_found;	// TRUE if we reach n_max without finding sol'n

  int   has_part2;		// TRUE once we set the remaining stuff below:
  int   anch_outside;		// # of anchors that fall outside any domain in Viterbi path
  int   anch_unique;		// # of anchors that map 1:1 to a domain in Viterbi path
  int   anch_multiple;		// # of anchors that fall in same Viterbi domain with other anchor(s)
  int   dom_zero;		// # of Viterbi domains with no anchors
  int   dom_one;		// # of Viterbi domains that map 1:1 to anchor ( == anch_unique)
  int   dom_multiple;		// # of Viterbi domains with >1 anchor
} H4_MPAS_STATS;


extern H4_MPAS_PARAMS *h4_mpas_params_Create(void);
extern void            h4_mpas_params_Destroy(H4_MPAS_PARAMS *params);

extern H4_MPAS_STATS  *h4_mpas_stats_Create(void);
extern int             h4_mpas_stats_Reuse(H4_MPAS_STATS *stats);
extern int             h4_mpas_stats_Dump(FILE *ofp, H4_MPAS_STATS *stats);
extern void            h4_mpas_stats_Destroy(H4_MPAS_STATS *stats);

#endif /*h4MPAS_INCLUDED*/
