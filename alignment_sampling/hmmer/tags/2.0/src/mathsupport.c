/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/


/* mathsupport.c
 * SRE, Mon Nov 11 15:07:33 1996
 * 
 * Miscellaneous mathematical functions.
 * General functions are in the SQUID library sre_math.c.
 * These functions are too HMM-specific to warrant being in the
 * SQUID library.
 * 
 */


#include <math.h>
#include "funcs.h"
#include "config.h"
#include "structs.h"
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

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
PValue(struct plan7_s *hmm, float sc)
{
  double pval;
  double pval2;
				/* the bound from Bayes */
  pval = 1. / (1.+sreEXP2(sc));
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
 */
float 
LogSum(float p1, float p2)
{
  if (p1 > p2)
    return p1 + log(1. + exp(p2-p1));
  else
    return p2 + log(1. + exp(p1-p2));
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
 *           
 * Args:     p1,p2 -- scaled integer log_2 probabilities to be summed
 *                    in probability space.
 *                    
 * Return:   scaled integer log_2 probability of the sum.
 */
int 
ILogsum(int p1, int p2)
{
  static int firsttime = 1;
  static int lookup[LOGSUM_TBL];
  int    diff;
  int    i;

  if (firsttime) {
    for (i = 0; i < LOGSUM_TBL; i++) 
      lookup[i] = (int) (INTSCALE * 1.44269504 * 
			 (log(1.+exp(0.69314718 * (float) -i/INTSCALE))));
    firsttime = 0;
  }

  diff = p1-p2;
  if      (diff >=  LOGSUM_TBL) return p1;
  else if (diff <= -LOGSUM_TBL) return p2;
  else if (diff > 0)            return p1 + lookup[diff];
  else                          return p2 + lookup[-diff];
} 


/* Function: LogNorm()
 * 
 * Purpose:  Normalize a vector of log likelihoods, changing it
 *           to a probability vector. Be careful of overflowing exp().
 *           Implementation adapted from Graeme Mitchison.
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

/* Function: SampleDirichlet()
 * 
 * Purpose:  Given a Dirichlet distribution defined by
 *           a vector of n alpha terms, sample of probability
 *           distribution of dimension n.
 *           
 *           This code was derived from source provided 
 *           by Betty Lazareva, from Gary Churchill's group.
 *           
 * Args:     alpha - vector of Dirichlet alphas components
 *           n     - number of components   
 *           ret_p - RETURN: sampled probability vector.
 *           
 * Return:   (void)
 *           ret_p, an n-dimensional array alloced by the caller,
 *           is filled.
 */
void
SampleDirichlet(float *alpha, int n, float *p)
{
  int    x;
  
  for (x = 0; x < n; x++)
    p[x] = SampleGamma(alpha[x]);
  FNorm(p, n);
}
  

/* Function: SampleGamma()
 * 
 * Purpose:  Return a random deviate distributed as Gamma(alpha, 1.0).
 *           Uses two different accept/reject algorithms, one
 *           for 0<alpha<1, the other for 1<=alpha. 
 *           
 *           Code modified from source provided by Betty Lazareva
 *           and Gary Churchill.
 *            
 * Args:     alpha - order of gamma function
 *           
 * Return:   the gamma-distributed deviate.
 */                       
float
SampleGamma(float alpha)
{
  float U,V,X,W,lambda;

  if (alpha >= 1.0) 
    {
      /*CONSTCOND*/ while (1)
	{
	  lambda = sqrt(2.0*alpha -1.0);
	  U = sre_random();
	  V = U/(1-U);
	  X = alpha * pow(V, 1/lambda);
	  W = .25*exp(-X+alpha)*pow(V,1.0+alpha/lambda)*pow(1.0+1.0/V, 2.0);
	  if (sre_random() <= W)
	    return X;
	}
    }
  else if (alpha > 0.0)
    {
      /*CONSTCOND*/ while (1) 
	{
	  U = sre_random();
	  V = U*(1+ alpha/exp(1.0));
	  if (V > 1.0)
	    {
	      X = -log( (1-V+alpha/exp(1.0))/alpha);
	      if (sre_random() <= pow(X, alpha-1.0))
		return X;
	    }
	  else
	    {
	      X = pow(V,1.0/alpha);
	      if (sre_random() <= exp(-X))
		return X;
	    }
	}
    }
  Die("Invalid argument alpha < 0.0 to SampleGamma()");
  /*NOTREACHED*/
  return 0.0;
}

/* Function: SampleCountvector()
 * 
 * Purpose:  Given a probability vector p of dimensionality
 *           n, sample c counts and store them in cvec.
 *           cvec is n-dimensional and is alloced by the caller.
 */
void
SampleCountvector(float *p, int n, int c, float *cvec)
{
  int i;

  FSet(cvec, n, 0.0);
  for (i = 0; i < c; i++)
    cvec[FChoose(p,n)] += 1.0;
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


