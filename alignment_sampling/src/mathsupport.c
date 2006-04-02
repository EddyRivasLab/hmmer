/* mathsupport.c
 * Miscellaneous mathematical functions.
 * General functions are in the SQUID library sre_math.c.
 * These functions are too HMM-specific to warrant being in the
 * SQUID library.
 * 
 * SRE, Mon Nov 11 15:07:33 1996
 * SVN $Id: mathsupport.c 1388 2005-05-16 12:27:46Z eddy $
 */


#include "config.h"
#include "squidconf.h"

#include <math.h>
#include <float.h>
#ifdef HMMER_THREADS
#include <pthread.h>
#endif

#include "squid.h"
#include "dirichlet.h"		/* Gammln() is in dirichlet module */

#include "plan7.h"
#include "funcs.h"
#include "structs.h"


/* Function: Prob2Score()
 * 
 * Purpose:  Convert a probability to a scaled integer log_2 odds score. 
 *           Round to nearest integer (i.e. note use of +0.5 and floor())
 *           Return the score. 
 */
int
Prob2Score(float p, float null)
{
  if   (p == 0.0) return -INFTY;
  else            return MAX((int) floor(0.5 + INTSCALE * sreLOG2(p/null)),
			     -INFTY);
}

/* Function:  LL2Score()
 * Incept:    SRE, Mon May  2 08:19:36 2005 [St. Louis]
 *
 * Purpose:   Convert a log likelihood to a scaled integer log_2 odds score,
 *            rounded to nearest integer, given a <null> probability; 
 *            return the score. 
 *            
 *            Note that <ll> is a log(prob), but <null> is a probability.
 * 
 */
int
LL2Score(float ll, float null)
{
  int sc;
  sc = (int) floor(0.5 + INTSCALE * 1.44269504 * (ll - log(null)));
  if (sc < -INFTY) sc = -INFTY;
  return sc;
}


/* Function: Score2Prob()
 * 
 * Purpose:  Convert an integer log_2 odds score back to a probability;
 *           needs the null model probability, if any, to do the conversion.
 */
float 
Score2Prob(int sc, float null)
{
  if (sc == -INFTY) return 0.;
  else              return (null * sreEXP2((float) sc / INTSCALE));
}


/* Function: Scorify()
 * 
 * Purpose:  Convert a scaled integer log-odds score to a floating
 *           point score for output. (could be a macro but who cares.)
 */
float 
Scorify(int sc)
{
  return ((float) sc / INTSCALE);
}


/* Function: PValue()
 * Date:     SRE, Mon Oct 27 12:21:02 1997 [Sanger Centre, UK]
 * 
 * Purpose:  Convert an HMM score to a P-value.
 *           We know P(S>x) is bounded by 1 / (1 + exp_2^x) for a bit score of x.
 *           We can also use EVD parameters for a tighter bound if we have
 *           them available.
 *           
 * Args:     hmm - model structure, contains EVD parameters
 *           sc  - score in bits
 *           
 * Returns:  P value for score significance.          
 */
double
PValue(struct plan7_s *hmm, double sc)
{
  double pval;
  double pval2;
				/* the bound from Bayes */
  if      (sc >= sreLOG2(DBL_MAX))       pval = 0.0;
  else if (sc <= -1. * sreLOG2(DBL_MAX)) pval = 1.0;
  else                        pval = 1. / (1.+sreEXP2(sc));

				/* try for a better estimate from EVD fit */
  if (hmm != NULL && (hmm->flags & PLAN7_STATS))
    {		
      pval2 = ExtremeValueP(sc, hmm->mu, hmm->lambda);
      if (pval2 < pval) pval = pval2;
    }
  return pval;
}



/* Function: LogSum()
 * 
 * Purpose:  Returns the log of the sum of two log probabilities.
 *           log(exp(p1)+exp(p2)) = p1 + log(1 + exp(p2-p1)) for p1 > p2
 *           Note that this is in natural log space, not log_2.
 */
float 
LogSum(float p1, float p2)
{
  if (p1 > p2)
    return (p1-p2 > 50.) ? p1 : p1 + log(1. + exp(p2-p1));
  else
    return (p2-p1 > 50.) ? p2 : p2 + log(1. + exp(p1-p2));
}


/* Function: ILogsum()
 * 
 * Purpose:  Return the scaled integer log probability of
 *           the sum of two probabilities p1 and p2, where
 *           p1 and p2 are also given as scaled log probabilities.
 *         
 *           log(exp(p1)+exp(p2)) = p1 + log(1 + exp(p2-p1)) for p1 > p2
 *           
 *           For speed, builds a lookup table the first time it's called.
 *           LOGSUM_TBL is set to 20000 by default, in config.h.
 *
 *           Because of the one-time initialization, we have to
 *           be careful in a multithreaded implementation... hence
 *           the use of pthread_once(), which forces us to put
 *           the initialization routine and the lookup table outside
 *           ILogsum(). (Thanks to Henry Gabb at Intel for pointing
 *           out this problem.)
 *           
 * Args:     p1,p2 -- scaled integer log_2 probabilities to be summed
 *                    in probability space.
 *                    
 * Return:   scaled integer log_2 probability of the sum.
 */
static int ilogsum_lookup[LOGSUM_TBL];
static void 
init_ilogsum(void)
{
  int i;
  for (i = 0; i < LOGSUM_TBL; i++) 
    ilogsum_lookup[i] = (int) (INTSCALE * 1.44269504 * 
	   (log(1.+exp(0.69314718 * (float) -i/INTSCALE))));
}
int 
ILogsum(int p1, int p2)
{
  int    diff;
#ifdef HMMER_THREADS
  static pthread_once_t firsttime = PTHREAD_ONCE_INIT;
  pthread_once(&firsttime, init_ilogsum);
#else
  static int firsttime = 1;
  if (firsttime) { init_ilogsum(); firsttime = 0; }
#endif

  diff = p1-p2;
  if      (diff >=  LOGSUM_TBL) return p1;
  else if (diff <= -LOGSUM_TBL) return p2;
  else if (diff > 0)            return p1 + ilogsum_lookup[diff];
  else                          return p2 + ilogsum_lookup[-diff];
} 

/* Function: LogNorm()
 * 
 * Purpose:  Normalize a vector of log likelihoods, changing it
 *           to a probability vector. Be careful of overflowing exp().
 *           Implementation adapted from Graeme Mitchison.
 *           [deprecated: use vectorops.c::FLogNorm(), DLogNorm().]
 *
 * Args:     vec - vector destined to become log probabilities
 *           n   - length of vec 
 */
void
LogNorm(float *vec, int n)
{
  int   x;
  float max   = -1.0e30;
  float denom = 0.;

  for (x = 0; x < n; x++)
    if (vec[x] > max) max = vec[x];
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      denom += exp(vec[x] - max);
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      vec[x] = exp(vec[x] - max) / denom;
    else
      vec[x] = 0.0;
}
  

/* Function: Logp_cvec()
 *
 * Purpose:  Calculates ln P(cvec|dirichlet), the log probability of a 
 *           count vector given a Dirichlet distribution. Adapted 
 *           from an implementation by Graeme Mitchison.
 *           
 * Args:     cvec  - count vector
 *           n     - length of cvec
 *           alpha - Dirichlet alpha terms
 *           
 * Return:   log P(cvec|dirichlet)                 
 */
float
Logp_cvec(float *cvec, int n, float *alpha)
{
  float lnp;                   /* log likelihood of P(cvec | Dirichlet) */
  float sum1, sum2, sum3;
  int   x;

  sum1 = sum2 = sum3 = lnp = 0.0;
  for (x = 0; x < n; x++)
    {
      sum1 += cvec[x] + alpha[x];
      sum2 += alpha[x];
      sum3 += cvec[x];
      lnp  += Gammln(alpha[x] + cvec[x]);
      lnp  -= Gammln(cvec[x] + 1.);
      lnp  -= Gammln(alpha[x]);
    }
  lnp -= Gammln(sum1);
  lnp += Gammln(sum2);
  lnp += Gammln(sum3 + 1.);
  return lnp;
}


/* Function: P_PvecGivenDirichlet()
 * 
 * Purpose:  Calculate the log probability of a probability
 *           vector given a single Dirichlet component, alpha.
 *           Follows Sjolander (1996) appendix, lemma 2.
 *           
 * Return:   log P(p | alpha)
 */          
float
P_PvecGivenDirichlet(float *p, int n, float *alpha)
{
  float sum;		        /* for Gammln(|alpha|) in Z     */
  float logp;			/* RETURN: log P(p|alpha)       */
  int x;

  sum = logp = 0.0;
  for (x = 0; x < n; x++)
    if (p[x] > 0.0)		/* any param that is == 0.0 doesn't exist */
      {
	logp += (alpha[x]-1.0) * log(p[x]);
	logp -= Gammln(alpha[x]);
	sum  += alpha[x];
      }
  logp += Gammln(sum);
  return logp;
}



/************************************************************
 * @LICENSE@
 ************************************************************/


