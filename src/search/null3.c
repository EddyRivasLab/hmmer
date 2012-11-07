/* "null3" model: biased composition correction
 * 
 * Contents:
 *   1. Null3 estimation algorithm.
 *   2. Copyright and license information.
 *
 * Approach is based heavily on the null3 approach used in Infernal,
 * and described in its user guide, specifically based on
 * ScoreCorrectionNull3CompUnknown()
 *
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "base/p7_bg.h"
#include "base/p7_trace.h"
#include "misc/logsum.h"

/*****************************************************************
 * 1. Null3 estimation algorithm.
 *****************************************************************/


/* Function: p7_null3_score()
 *
 * Purpose:  Calculate a correction (in log_2 odds) to be applied
 *           to a sequence, using a null model based on the
 *           composition of the target sequence.
 *           The null model is constructed /post hoc/ as the
 *           distribution of the target sequence; if the target
 *           sequence is 40% A, 5% C, 5% G, 40% T, then the null
 *           model is (0.4, 0.05, 0.05, 0.4). This function is
 *           based heavily on Infernal's ScoreCorrectionNull3(),
 *           with two important changes:
 *            - it leaves the log2 conversion from NATS to BITS
 *              for the calling function.
 *            - it doesn't include the omega score modifier
 *              (based on prior probability of using the null3
 *              model), again leaving this to the calling function.
 *
 * Args:     abc   - alphabet for hit (only used to get alphabet size)
 *           dsq   - the sequence the hit resides in
 *           tr   - trace of the alignment, used to find the match states
 *                  (non-match chars are ignored in computing freq, not used if NULL)
 *           start - start position of hit in dsq
 *           stop  - end  position of hit in dsq
 *           bg    - background, used for the default null model's emission freq
 *           ret_sc - RETURN: the correction to the score (in NATS);
 *                   caller subtracts this from hit score to get
 *                   corrected score.
 * Return:   void, ret_sc: the log-odds score correction (in NATS).
 */
void
p7_null3_score(const ESL_ALPHABET *abc, const ESL_DSQ *dsq, P7_TRACE *tr, int start, int stop, P7_BG *bg, float *ret_sc)
{
  float score = 0.;
  int status;
  int i;
  float *freq;
  int dir;
  int tr_pos;

  ESL_ALLOC(freq, sizeof(float) * abc->K);
  esl_vec_FSet(freq, abc->K, 0.0);

  /* contract check */
  if(abc == NULL) esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "p7_null3_score() alphabet is NULL.%s\n", "");
  if(dsq == NULL) esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "p7_null3_score() dsq alphabet is NULL.%s\n", "");
  if(abc->type != eslRNA && abc->type != eslDNA) esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "p7_null3_score() expects alphabet of RNA or DNA.%s\n", "");

  dir = start < stop ? 1 : -1;

  if (tr != NULL) {
    /* skip the parts of the trace that precede the first match state */
    tr_pos = 2;
    i = start;
    while (tr->st[tr_pos] != p7T_MG && tr->st[tr_pos] != p7T_ML) {
      if (tr->st[tr_pos] == p7T_N)
        i += dir;
      tr_pos++;
    }

    /* tally frequencies from characters hitting match state*/
    while (tr->st[tr_pos] != p7T_E) {
      if (tr->st[tr_pos] == p7T_MG || tr->st[tr_pos] == p7T_ML) {
        if(esl_abc_XIsGap(abc, dsq[i])) esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "in p7_null3_score(), res %d is a gap!%s\n", "");
        esl_abc_FCount(abc, freq, dsq[i], 1.);
      }
      if (tr->st[tr_pos] != p7T_DG && tr->st[tr_pos] != p7T_DL )
        i += dir;
      tr_pos++;
    }
  } else {
    /* tally frequencies from the full envelope */
    for (i=ESL_MIN(start,stop); i <= ESL_MAX(start,stop); i++)
    {
      if(esl_abc_XIsGap(abc, dsq[i])) esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "in p7_null3_score(), res %d is a gap!%s\n", "");
      esl_abc_FCount(abc, freq, dsq[i], 1.);
    }
  }

  esl_vec_FNorm(freq, abc->K);


  /* now compute score modifier (nats) - note: even with tr!=NULL, this includes the unmatched characters*/
  for (i = 0; i < abc->K; i++)
    score += freq[i]==0 ? 0.0 : esl_logf( freq[i]/bg->f[i] ) * freq[i] * ( (stop-start)*dir +1) ;

  /* Return the correction to the bit score. */
  score = p7_FLogsum(0., score);
  *ret_sc = score;

  return;

 ERROR:
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "p7_null3_score() memory allocation error.%s\n", "");
  return; /* never reached */

}





/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id: null3.c 3474 2011-02-17 13:25:32Z wheelert $
 * SVN $URL$
 *****************************************************************/

