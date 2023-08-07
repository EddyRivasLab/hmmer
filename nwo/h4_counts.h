/* H4_COUNTS: profile HMM in counts-collection form.
 * 
 * Profile structure used to collect counts, where we use double
 * instead of float, and we only need the probability model part of
 * the H4_PROFILE structure. Deep alignments (and simulations) can
 * exceed dynamic range of float - especially on II transitions.
 *                                 
 * Essentially just the probability model part of an H4_PROFILE, in
 * double precision. Shares transition indices with H4_PROFILE.
 */
#ifndef h4COUNTS_INCLUDED
#define h4COUNTS_INCLUDED
#include <h4_config.h>

#include "esl_alphabet.h"

#include "h4_profile.h"   // to get h4_TMM, etc.: transition indices are same as in H4_PROFILE

typedef struct {
  int    M;
  double **t;
  double **e;

  const ESL_ALPHABET *abc;
} H4_COUNTS;

extern H4_COUNTS *h4_counts_Create(const ESL_ALPHABET *abc, int M);
extern void       h4_counts_Destroy(H4_COUNTS *ctm);

extern int        h4_counts_Dump(FILE *fp, H4_COUNTS *ctm);

#endif /* h4COUNTS_INCLUDED */

