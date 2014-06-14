/* Support for the MPAS (most probable anchor set) algorithm,
 * especially development, testing, and debugging; the P7_MPAS_STATS
 * and P7_MPAS_PARAMS structures.
 * 
 * MPAS algorithm is implemented twice, once in "reference" (a
 * development testbed) and once in "sparse" (production version).
 * These support routines are shared between the implementations,
 * which is why they're here in one place, rather than with the MPAS
 * implementations.
 * 
 */
#include "p7_config.h"

#include "base/p7_anchors.h"
#include "base/p7_trace.h"

#include "search/p7_mpas.h"

int
p7_mpas_stats_Init(P7_MPAS_STATS *stats)
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
  stats->nsamples_in_best     = 0;      // only counts suboptimal path samples; Viterbi path doesn't count
  stats->best_is_viterbi      = TRUE;   // best is Viterbi, until proven otherwise by MPAS algorithm
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
p7_mpas_stats_Dump(FILE *ofp, P7_MPAS_STATS *stats)
{
  fprintf(ofp, "# Stats on the ASC solution:\n");

  if (stats->has_part1)
    {
      fprintf(ofp, "# tot_iterations       = %d\n",    stats->tot_iterations);
      fprintf(ofp, "# tot_asc_calculations = %d\n",    stats->tot_asc_calculations);
      fprintf(ofp, "# Viterbi score (nats) = %.2f\n",  stats->vsc);
      fprintf(ofp, "# Forward score (nats) = %.2f\n",  stats->fsc);
      fprintf(ofp, "# vit_asc              = %.2f\n",  stats->vit_asc);
      fprintf(ofp, "# vit_ascprob          = %6.4f\n", stats->vit_ascprob);
      fprintf(ofp, "# best_asc             = %.2f\n",  stats->best_asc);
      fprintf(ofp, "# best_ascprob         = %6.4f\n", stats->best_ascprob);
      fprintf(ofp, "# nsamples_in_best     = %d\n",    stats->nsamples_in_best);
      fprintf(ofp, "# best_is_viterbi      = %s\n",    (stats->best_is_viterbi    ? "TRUE" : "FALSE"));
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

int
p7_mpas_stats_CompareAS2Trace(P7_MPAS_STATS *stats, const P7_ANCHORS *anch, const P7_TRACE *tr)
{
  int ad;
  int td              = 0;
  int anch_in_this_td = 0;

  stats->anch_outside    = 0;
  stats->anch_unique     = 0;
  stats->anch_multiple   = 0;

  stats->dom_zero        = 0;
  stats->dom_one         = 0;
  stats->dom_multiple    = 0;

  /* For n domains in tr:
   *   they can either be hit 0 times, 1 time, or 2+ times by anchors.
   * For m anchors in anch:
   *   they can either fall outside any domain, uniquely in a domain, or multiply in a domain.
   *   
   * Watch out: ad (in anchor set) is 1..D; td (in trace) is 0..D-1.  
   */
  for (ad = 1; ad <= anch->D; ad++)
    {
      if   (anch->a[ad].i0 < tr->sqfrom[td] || td == tr->ndom) 
	stats->anch_outside++;
      else if (anch->a[ad].i0 >= tr->sqfrom[td] && anch->a[ad].i0 <= tr->sqto[td])
	anch_in_this_td++;
      else 
	{
	  /* we have to advance <td>, and try again */
	  if      (anch_in_this_td == 0) { stats->dom_zero++; }
	  else if (anch_in_this_td == 1) { stats->anch_unique++; stats->dom_one++; }
	  else if (anch_in_this_td > 1)  { stats->anch_multiple += anch_in_this_td; stats->dom_multiple++; }
	  anch_in_this_td = 0;
	  td++;
	  ad--;			/* forces reevaluation of <ad> when we go back around; a bit hacky! */
	}
    }
  
  /* we're out of anchors. If td == tr->ndom, we also know we
   * handled what happened with anchors in the last <td>. But if
   * td == tr->ndom-1, we haven't yet resolved what happened with final <td> yet,
   * and if td is even smaller, we have some dom_zero's to count.
   */
  for (; td < tr->ndom; td++)
    {
      if      (anch_in_this_td == 0) { stats->dom_zero++; }
      else if (anch_in_this_td == 1) { stats->anch_unique++; stats->dom_one++; }
      else if (anch_in_this_td > 1)  { stats->anch_multiple += anch_in_this_td; stats->dom_multiple++; }
      anch_in_this_td = 0;
    }

  ESL_DASSERT1(( stats->dom_zero     + stats->dom_one     + stats->dom_multiple  == tr->ndom ));
  ESL_DASSERT1(( stats->anch_outside + stats->anch_unique + stats->anch_multiple == anch->D  ));

  stats->has_part2 = TRUE;
  return eslOK;
}

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
