/* H4_MPAS_PARAMS, H4_MPAS_STATS structures for control parameters and statistics
 * collection for the most probable anchor set (MPAS) algorithm implementations.
 *
 * Contents:
 *    1. H4_MPAS_PARAMS: control parameters for MPAS algorithm
 *    2. H4_MPAS_STATS:  MPAS statistics collection
 *
 *    
 * See also:
 *    reference_mpas.* : reference implementation of MPAS, used for testing
 *    sparse_mpas.*    : sparse DP implementation of MPAS, used in production code
 */
#include <h4_config.h>

#include <stdio.h>

#include "easel.h"

#include "h4_mpas.h"

/*****************************************************************
 * 1. H4_MPAS_PARAMS: control parameters for MPAS algorithm
 *****************************************************************/

H4_MPAS_PARAMS *
h4_mpas_params_Create(void)
{
  H4_MPAS_PARAMS *params = NULL;
  int             status;

  ESL_ALLOC(params, sizeof(H4_MPAS_PARAMS));
  params->max_iterations = h4_MPAS_MAX_ITERATIONS;  // limit MPAS algorithm to <max_iterations> stochastic traces
  params->loss_threshold = h4_MPAS_LOSS_THRESHOLD;  // stoppage criterion: probability that better AS exists, hasn't been found yet
  params->nmax_sampling  = h4_MPAS_NMAX_SAMPLING;   // if TRUE, take all <max_iterations> samples, don't apply other stoppage tests
  params->be_verbose     = h4_MPAS_BE_VERBOSE;      // if TRUE, MPAS procedure printf's internal info for debugging
  return params;

 ERROR:
  h4_mpas_params_Destroy(params);
  return NULL;
}

void
h4_mpas_params_Destroy(H4_MPAS_PARAMS *params)
{
  free(params);
}
      
/*****************************************************************
 * 2. H4_MPAS_STATS: MPAS statistics collection 
 *****************************************************************/

H4_MPAS_STATS *
h4_mpas_stats_Create(void)
{
  H4_MPAS_STATS *stats = NULL;
  int            status;

  ESL_ALLOC(stats, sizeof(H4_MPAS_STATS));
  h4_mpas_stats_Reuse(stats);
  return stats;

 ERROR:
  h4_mpas_stats_Destroy(stats);
  return NULL;
}

int
h4_mpas_stats_Reuse(H4_MPAS_STATS *stats)
{
  stats->has_part1            = FALSE;
  stats->tot_iterations       = 0;     // 0 is a possible answer, if Viterbi path immediately gave provably best AS
  stats->tot_asc_calculations = 0;
  stats->vsc                  = 0.0;
  stats->fsc                  = 0.0;
  stats->vit_asc              = 0.0;
  stats->vit_ascprob          = 0.0;
  stats->best_asc             = 0.0;
  stats->best_ascprob         = 0.0;
  stats->tot_prob             = 0.0;
  stats->nsamples_in_best     = 0;      // only counts suboptimal path samples; Viterbi path doesn't count
  stats->best_is_viterbi      = TRUE;   // best is Viterbi, until proven otherwise by MPAS algorithm
  stats->late_solution        = FALSE;
  stats->solution_not_found   = FALSE;

  stats->has_part2            = FALSE;
  stats->anch_outside         = 0;
  stats->anch_unique          = 0;
  stats->anch_multiple        = 0;
  stats->dom_zero             = 0;
  stats->dom_one              = 0;
  stats->dom_multiple         = 0;

  return eslOK;
}


int
h4_mpas_stats_Dump(FILE *ofp, H4_MPAS_STATS *stats)
{
  fprintf(ofp, "# Stats on the ASC solution:\n");

  if (stats->has_part1)
    {
      fprintf(ofp, "# tot_iterations       = %d\n",    stats->tot_iterations);
      fprintf(ofp, "# tot_asc_calculations = %d\n",    stats->tot_asc_calculations);
      fprintf(ofp, "# Viterbi score (bits) = %.2f\n",  stats->vsc);
      fprintf(ofp, "# Forward score (bits) = %.2f\n",  stats->fsc);
      fprintf(ofp, "# vit_asc              = %.2f\n",  stats->vit_asc);
      fprintf(ofp, "# vit_ascprob          = %6.4f\n", stats->vit_ascprob);
      fprintf(ofp, "# best_asc             = %.2f\n",  stats->best_asc);
      fprintf(ofp, "# best_ascprob         = %6.4f\n", stats->best_ascprob);
      fprintf(ofp, "# tot_prob             = %6.4f\n", stats->tot_prob);
      fprintf(ofp, "# nsamples_in_best     = %d\n",    stats->nsamples_in_best);
      fprintf(ofp, "# best_is_viterbi      = %s\n",    (stats->best_is_viterbi    ? "TRUE" : "FALSE"));
      fprintf(ofp, "# late_solution        = %s\n",    (stats->late_solution      ? "TRUE" : "FALSE"));
      fprintf(ofp, "# solution_not_found   = %s\n",    (stats->solution_not_found ? "TRUE" : "FALSE"));
    }

  if (stats->has_part2) 
    {
      fprintf(ofp, "# anch_outside         = %d\n",    stats->anch_outside);
      fprintf(ofp, "# anch_unique          = %d\n",    stats->anch_unique);
      fprintf(ofp, "# anch_multiple        = %d\n",    stats->anch_multiple);
      fprintf(ofp, "# dom_zero             = %d\n",    stats->dom_zero);
      fprintf(ofp, "# dom_one              = %d\n",    stats->dom_one);
      fprintf(ofp, "# dom_multiple         = %d\n",    stats->dom_multiple);
    }
  return eslOK;
}


void
h4_mpas_stats_Destroy(H4_MPAS_STATS *stats)
{
  free(stats);
}


