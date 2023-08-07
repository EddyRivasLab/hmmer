/* This contains a compact representation of 8-bit bias-shifted scores for use in
 * diagonal recovery (standard [MS]SV) and extension (standard and FM-[MS]SV),
 * along with MAXL-associated prefix- and suffix-lengths, and optimal extensions
 * for FM-MSV.
 */
#ifndef p7SCOREDATA_INCLUDED
#define p7SCOREDATA_INCLUDED

#include <p7_config.h>

#include "p7_hmmwindow.h"
#include "dp_vector/p7_oprofile.h"

typedef struct p7_scoredata_s {
  int      M;
  uint8_t    *msv_scores;  //implicit (M+1)*K matrix, where M = # states, and K = # characters in alphabet
  uint8_t   **opt_ext_fwd;
  uint8_t   **opt_ext_rev;
  float      *prefix_lengths;
  float      *suffix_lengths;
  float      *fwd_scores;
  float     **fwd_transitions;
} P7_SCOREDATA;


extern P7_SCOREDATA   *p7_hmm_ScoreDataCreate(P7_OPROFILE *om, int do_opt_ext);
extern P7_SCOREDATA   *p7_hmm_ScoreDataClone(P7_SCOREDATA *src, int K);
extern int            p7_hmm_ScoreDataComputeRest(P7_OPROFILE *om, P7_SCOREDATA *data );
extern void           p7_hmm_ScoreDataDestroy( P7_SCOREDATA *data );
extern int            p7_hmm_initWindows (P7_HMM_WINDOWLIST *list);
extern P7_HMM_WINDOW *p7_hmm_newWindow (P7_HMM_WINDOWLIST *list, uint32_t id, uint32_t pos, uint32_t fm_pos, uint16_t k, uint32_t length, float score, uint8_t complementarity);

#endif /*p7SCOREDATA_INCLUDED*/
