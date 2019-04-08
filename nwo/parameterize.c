#include "h4_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "h4_prior.h"
#include "h4_profile.h"


/* Function:  h4_Parameterize()
 * Synopsis:  Convert counts to probability parameters in a profile HMM
 * Incept:    SRE, Sun 24 Jun 2018 [World Cup, Poland v. Colombia]
 *
 * Purpose:   Given an <hmm> containing observed counts, and prior
 *            parameters <pri>, convert <hmm> from counts to
 *            probabilities (mean posterior parameter estimates).
 *
 *            Sets model parameter fixed boundary conditions
 *            correctly.  This allows entropy weighting to blindly
 *            rescale all values in a model (almost of which are
 *            counts, but some are 1.0 probabilities in fixed boundary
 *            conditions that aren't observed counts) before 
 *            calling <h4_Parameterize()> on the reweighted counts.
 *            
 * Returns:   <eslOK> on success, and <hmm> is converted from 
 *            counts to mean posterior probability parameters.           
 */
int
h4_Parameterize(H4_PROFILE *hmm, const H4_PRIOR *pri)
{
  double c[h4_MAXABET];
  double p[h4_MAXABET];
  int    k;
  
  /* special case of no prior. */
  if (! pri) return h4_profile_Renormalize(hmm);

  /* match emissions, 1..M
   */
  for (k = 1; k <= hmm->M; k++) {
    esl_vec_F2D(hmm->e[k], hmm->abc->K, c);
    esl_mixdchlet_MPParameters(pri->em, c, p);
    esl_vec_D2F(p, hmm->abc->K, hmm->e[k]);
  }

  /* match transitions */
  for (k = 0; k < hmm->M; k++)
    {
      esl_vec_F2D(hmm->t[k], 3, c);
      esl_mixdchlet_MPParameters(pri->tm, c, p);
      esl_vec_D2F(p, 3, hmm->t[k]);
    }
  hmm->t[0][h4_TMI] = 0.0;      // [0] is G->{M1,D1}. Above loop parameterized it with usual match
  esl_vec_FNorm(hmm->t[0], 3);  //    transition prior. Set the M->I transition to 0, and renormalize.

  /* insert transitions */
  for (k = 1; k < hmm->M; k++)
    {
      esl_vec_F2D(hmm->t[k]+3, 3, c);
      esl_mixdchlet_MPParameters(pri->ti, c, p);
      esl_vec_D2F(p, 3, hmm->t[k]+3);
    }

  /* delete transitions */
  for (k = 1; k < hmm->M; k++) {
    esl_vec_F2D(hmm->t[k]+6, 3, c);
    esl_mixdchlet_MPParameters(pri->td, c, p);
    esl_vec_D2F(p, 3, hmm->t[k]+6);
  }

  h4_profile_SetConventions(hmm);
  return eslOK;
}
