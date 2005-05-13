/* mathsupport.c
 * Miscellaneous mathematical functions.
 * General functions are in the SQUID library sre_math.c.
 * These functions are too HMM-specific to warrant being in the
 * SQUID library.
 * 
 * SRE, Mon Nov 11 15:07:33 1996
 * SVN $Id$
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
  else            return (int) floor(0.5 + INTSCALE * sreLOG2(p/null));
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


/* Function:  LPValue()
 * Incept:    SRE, Fri May 13 14:15:49 2005 [St. Louis]
 *
 * Purpose:   Convert an HMM Viterbi score to a P-value.
 *            We know P(S>=x) is bounded by 1 / (1 + exp_2^x) for a bit
 *            score of x. We can also use EVD parameters for a tighter
 *            bound if we have them available.
 *            
 *            LPValue(), unlike PValue(), takes the length of the
 *            target sequence L into account, assuming
 *            Karlin/Altschul statistics.
 *
 * Args:      hmm - the HMM contains Lbase, mu, lambda, kappa, sigma params.
 *            L   - target seq length in residues
 *            sc  - bitscore to evaluate      
 *
 * Returns:   pvalue, P(S>=x).
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      Emails from Bill Bruno 2004-2005, and STL9/87. Bill
 *            contributed the original version of LPValue().
 */
double
LPValue(struct plan7_s *hmm, double L, double sc)
{
  double L1, L2;
  double mu2;
  double pval;
  double pval2;
				/* the bound from Bayes */
  if      (sc >= sreLOG2(DBL_MAX))       pval = 0.0;
  else if (sc <= -1. * sreLOG2(DBL_MAX)) pval = 1.0;
  else    pval = 1. / (1.+sreEXP2(sc));
				/* try for a better estimate from EVD fit */
  if (hmm != NULL && (hmm->flags & PLAN7_STATS))
    {		
      L1  = EdgeCorrection(L,          hmm->kappa, hmm->sigma);
      L2  = EdgeCorrection(hmm->Lbase, hmm->kappa, hmm->sigma);
      mu2 = hmm->mu + sreLOG2(L1/L2) / hmm->lambda;

      pval2 = ExtremeValueP(sc, mu2, hmm->lambda);
      if (pval2 < pval) pval = pval2;
    }

  return pval;
}


/* Function:  EdgeCorrection()
 * Incept:    SRE, Fri May 13 13:08:18 2005 [St. Louis]
 *
 * Purpose:   Edge effect correction: calculate an "effective"
 *            sequence length L', given the real length <L> and
 *            a length distribution for optimal alignments.
 *            This length distribution P(x) is assumed to be a
 *            left-censored Gaussian with mean <kappa> and
 *            std deviation <sigma>, x > 0.
 *            
 * Args:      N     - target sequence length in residues (cast to double)
 *            kappa - mean alignment length
 *            sigma - std dev on alignment length
 *
 * Returns:   N', the effective sequence length.
 *
 * Discussion: BLAST's edge correction is L-l(s); subtract the
 *            alignment length from L, where l(s) is the expected
 *            length of an alignment that scores s. l(s), in turn,
 *            is assumed to be a linear function of s: l(s) =
 *            alpha * s + beta. [Altschul01] discusses this
 *            correction, and how alpha can alternatively be used to
 *            calculate an edge correction on lambda, instead of an
 *            edge correction on L. 
 *            
 *            Bill Bruno points out that this can break down for small
 *            target sequence lengths L, especially for L<l(s), where
 *            you'd calculate a negative effective seq length (which
 *            then gives you a negative probability for the P-value
 *            from K/A statistics).
 *            
 *            Bill suggested instead calculating the expected effective
 *            seq length, integrating over all possible alignment
 *            lengths, rather than assuming a single fixed length l(s).
 *            
 *            We assume that alignment lengths l(s) are distributed as a
 *            left-censored Gaussian, l(s) > 0, with mean kappa and std
 *            deviation sigma.
 *            
 *            Then, N' =   \int_0^N (N-x) Prob(x | k,s) dx
 *                        ----------------------------------
 *                       1 - \int_{-\infty}^0 Prob(x | k,s) dx
 *                       
 *            This function implements the solution of this equation.
 *            See STL9/92-93 for the derivation.
 *            
 * Portability: Requires erf() -- double erf(double x) -- the error
 *            function of x. erf() is in ISO C99 spec, but is not ANSI
 *            C or POSIX. 
 *            
 * Xref:      STL9/92-93;
 *            notebook/0510-bruno-lengthcorrection;
 *            [Altschul01];
 *            emails from Bill Bruno, 2004-2005.
 */
double
EdgeCorrection(double L, double kappa, double sigma)
{
  double term1, term2, term3;
  double z1, z2;
  double Lp;

  z1 = (N-kappa) / sigma;
  z2 = kappa / sigma;

  term1 = (N-kappa) * (erf(z1/sqrt(2.)) + erf(z2/sqrt(2.)));
  term2 = sigma * sqrt(2.) / sqrt(SQDCONST_PI) * (exp(-0.5 * z1 * z1) - exp(-0.5 * z2 * z2));
  term3 = 1. + erf(z2/sqrt(2.));

  Lp = (term1 + term2)/ term3;
  return Lp;
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


