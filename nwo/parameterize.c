#include "h4_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "h4_counts.h"
#include "h4_prior.h"
#include "h4_profile.h"


/* Function:  h4_parameterize()
 * Synopsis:  Convert counts to probability parameters in a profile HMM
 * Incept:    SRE, Sun 24 Jun 2018 [World Cup, Poland v. Colombia]
 *
 * Purpose:   Given counts-collection model <ctm> and prior parameters
 *            <pri>, parameterize <hmm> with mean posterior parameter
 *            estimates. Caller provides an allocated <hmm> for the same
 *            size <M> and alphabet <abc> as <ctm>.
 *
 *            Sets model parameter fixed boundary conditions
 *            correctly.  This allows entropy weighting to blindly
 *            rescale all values in a model (almost of which are
 *            counts, but some are 1.0 probabilities in fixed boundary
 *            conditions that aren't observed counts) before 
 *            calling <h4_parameterize()> on the reweighted counts.
 *            
 *            <pri> can be NULL, in which case counts are simply
 *            renormalized to set <hmm> parameters.
 *            
 * Returns:   <eslOK> on success; <hmm> has its probability parameters
 *            set.
 *            
 * Note:      Counts in <ctm> are collected in doubles; parameters
 *            in <hmm> are floats.
 *            
 *            Deliberately doesn't configure the profile with
 *            scores. For some uses, we only need probability params
 *            and we're in an expensive critical path - for example,
 *            in the entropy weighting iterative optimization.
 */
int
h4_parameterize(const H4_COUNTS *ctm, const H4_PRIOR *pri, H4_PROFILE *hmm)
{
  double p[h4_MAXABET];
  int    M = ctm->M;
  int    K = ctm->abc->K;
  int    k;

  ESL_DASSERT1(( h4_MAXABET  >= h4_NT ));  // using single tmp space p[] for both t,e assumes this
  ESL_DASSERT1(( ctm->M      == hmm->M ));
  ESL_DASSERT1(( ctm->abc->K == hmm->abc->K ));
  
  if (pri)
    {
      esl_mixdchlet_MPParameters(pri->tm, ctm->t[0], p);
      p[h4_TMI] = 0.;                 // [0] is G->{M1,D1}. Above loop parameterized it with usual match
      esl_vec_DNorm(p, 3);            //    transition prior. Set the M->I transition to 0, and renormalize.
      esl_vec_D2F(p, 3, hmm->t[0]);
      for (k = 1; k < hmm->M; k++)
	{
	  esl_mixdchlet_MPParameters(pri->em, ctm->e[k],   p);   esl_vec_D2F(p, K, hmm->e[k]);
	  esl_mixdchlet_MPParameters(pri->tm, ctm->t[k],   p);   esl_vec_D2F(p, 3, hmm->t[k]);   // M
	  esl_mixdchlet_MPParameters(pri->ti, ctm->t[k]+3, p);   esl_vec_D2F(p, 3, hmm->t[k]+3); // I
	  esl_mixdchlet_MPParameters(pri->td, ctm->t[k]+6, p);   esl_vec_D2F(p, 3, hmm->t[k]+6); // D
	}
      esl_mixdchlet_MPParameters(pri->em, ctm->e[M], p);         esl_vec_D2F(p, K, hmm->e[k]);
    }
  else
    {
      esl_vec_DCopy(ctm->t[0], 3, p); esl_vec_DNorm(p, 3); esl_vec_D2F(p, 3, hmm->t[0]);
      for (k = 1; k < hmm->M; k++)
	{
	  esl_vec_DCopy(ctm->e[k],   K, p); esl_vec_DNorm(p, K); esl_vec_D2F(p, K, hmm->e[k]);
	  esl_vec_DCopy(ctm->t[k],   3, p); esl_vec_DNorm(p, 3); esl_vec_D2F(p, 3, hmm->t[k]);
	  esl_vec_DCopy(ctm->t[k]+3, 3, p); esl_vec_DNorm(p, 3); esl_vec_D2F(p, 3, hmm->t[k]+3);
	  esl_vec_DCopy(ctm->t[k]+6, 3, p); esl_vec_DNorm(p, 3); esl_vec_D2F(p, 3, hmm->t[k]+6);
	}
      esl_vec_DCopy(ctm->e[M], K, p); esl_vec_DNorm(p, K); esl_vec_D2F(p, K, hmm->e[M]);
    }

  h4_profile_SetConventions(hmm);
  hmm->flags |= h4_HASPROBS;
  hmm->flags &= (~h4_HASBITS);
  return eslOK;
}
