/* P7_OPROFILE, continued: more debugging and development utilities.
 * 
 * This is split off from p7_oprofile.c just to streamline
 * dependencies. p7_oprofile itself has few dependencies, but
 * p7_oprofile_Sample(), which we use for testing, depends on more of
 * the code.
 */
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"

#include "base/p7_bg.h"
#include "base/p7_hmm.h"
#include "base/p7_profile.h"
#include "build/modelsample.h"
#include "search/modelconfig.h"
#include "dp_vector/p7_oprofile.h"

/* Function:  p7_oprofile_Sample()
 * Synopsis:  Sample a random profile.
 *
 * Purpose:   Sample a random profile of <M> nodes for alphabet <abc>,
 *            using <r> as the source of random numbers. Parameterize
 *            it for generation of target sequences of mean length
 *            <L>. Calculate its log-odds scores using background
 *            model <bg>.
 *            
 *            Caller may optionally obtain the corresponding hmm by
 *            passing a non-<NULL> <opt_hmm>, and/or the corresponding
 *            profile by passing a non-<NULL> <opt_gm>. If the <gm> is
 *            obtained, it is configured for local-only mode and for a
 *            target length of <L>, so that its scores will match the
 *            <om> (as closely as roundoff allows).
 *            
 * Args:      r       - random number generator
 *            abc     - emission alphabet 
 *            bg      - background frequency model
 *            M       - size of sampled profile, in nodes
 *            L       - configured target seq mean length
 *            opt_hmm - optRETURN: sampled HMM
 *            opt_gm  - optRETURN: sampled normal profile, (local,L) mode
 *            opt_om  - optRETURN: optimized profile, length config'ed to L
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
		   P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **opt_om)
{
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_OPROFILE    *om   = NULL;
  int             status;

  if ((gm = p7_profile_Create (M, abc)) == NULL)  { status = eslEMEM; goto ERROR; }
  if ((om = p7_oprofile_Create(M, abc)) == NULL)  { status = eslEMEM; goto ERROR; }

  if ((status = p7_modelsample(r, M, abc, &hmm))        != eslOK) goto ERROR;
  if ((status = p7_profile_ConfigLocal(gm, hmm, bg, L)) != eslOK) goto ERROR;
  if ((status = p7_oprofile_Convert(gm, om))            != eslOK) goto ERROR;
  if ((status = p7_oprofile_ReconfigLength(om, L))      != eslOK) goto ERROR;

  if (opt_hmm) *opt_hmm = hmm; else p7_hmm_Destroy(hmm);
  if (opt_gm)  *opt_gm  = gm;  else p7_profile_Destroy(gm);
  if (opt_om)  *opt_om  = om;  else p7_oprofile_Destroy(om);
  return eslOK;

 ERROR:
  if (hmm) p7_hmm_Destroy(hmm); 
  if (gm)  p7_profile_Destroy(gm);
  if (om)  p7_oprofile_Destroy(om);
  if (opt_hmm) *opt_hmm = NULL;
  if (opt_gm)  *opt_gm  = NULL; 
  if (opt_om)  *opt_om  = NULL; 
  return status;
}



