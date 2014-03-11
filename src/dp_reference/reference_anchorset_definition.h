#ifndef p7REFERENCE_ASC_INCLUDED

#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "base/p7_profile.h"
#include "base/p7_coords2.h"

typedef struct {
  int   tot_iterations;		/* # of samples we took                          */
  int   tot_asc_calculations;	/* # of non-dup samples, where we did ASC        */
  int   n_to_find;		/* # of iterations it took to find the solution  */
  float vsc;			/* Viterbi score (raw; nats)                     */
  float fsc;			/* Forward score (raw; nats)                     */
  float best_asc;		/* anchor-constrained score of solution          */
  float vit_asc;		/* ASC score of Viterbi annotation               */
  float vit_ascprob;		/* probability of the Viterbi annotation         */
  float best_ascprob;		/* probability of the asc: exp(best_asc - fwdsc) */
  int   nsamples_in_best;	/* number of samples we saw for the solution     */
  int   best_is_viterbi;	/* TRUE if solution is Viterbi annotation        */
  int   best_found_late;	/* TRUE if a better solution was found late      */
  int   solution_not_found;	/* TRUE if we find no solution                   */

  int   anch_outside;		/* # of anchors that fall outside any domain in Viterbi trace */
  int   anch_unique;		/* # of anchors that map 1:1 to a domain in Viterbi trace       */
  int   anch_multiple;		/* # of anchors that fall in same Viterbi domain with other anchor(s) */

  int   dom_zero;		/* # of Viterbi domains with no anchors */
  int   dom_one;		/* # of Viterbi domains that map 1:1 to anchor ( == anch_unique) */
  int   dom_multiple;		/* # of Viterbi domains with >1 anchor  */

} P7_XSTATS_ASC;

extern int p7_ReferenceASCSearch(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_COORDS2 *anch, P7_XSTATS_ASC *stats);

#endif /*p7REFERENCE_ASC_INCLUDED*/
