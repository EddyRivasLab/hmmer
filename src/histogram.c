/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1997 Sean R. Eddy
 *
 *   This source code is distributed under the terms of the
 *   GNU General Public License. See the files COPYING and
 *   GNULICENSE for details.
 *
 ************************************************************/

/* histogram.c
 * SRE, Sat Jan 20 16:16:17 1996
 * 
 * Accumulation, printing, and fitting of score histograms
 * from database searches.
 ************************************************************
 * Basic API:
 * 
 * struct histogram_s *h;
 * 
 * h = AllocHistogram(min_hint, max_hint, lumpsize);
 * 
 * while (getting scores x) AddToHistogram(h, x);
 * 
 * ExtremeValueFitHistogram(h, high_hint);   
 * PrintASCIIHistogram(fp, h);   
 * FreeHistogram(h);
 */

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include "squid.h"
#include "config.h"
#include "structs.h"
#include "funcs.h"

#include <assert.h>

static void  evd_bruteforce(struct histogram_s *h, int lowbound, int highbound);
static float evd_goodness(struct histogram_s *h, int lowbound, int highbound,
			  float mu, float lambda);


/* Function: AllocHistogram()
 * 
 * Purpose:  Allocate and return a histogram structure.
 *           min and max are your best guess. They need
 *           not be absolutely correct; the histogram
 *           will expand dynamically to accomodate scores
 *           that exceed these suggested bounds. The amount
 *           that the histogram grows by is set by "lumpsize".
 * 
 * Args:     min:      minimum score (integer)
 *           max:      maximum score (integer)
 *           lumpsize: when reallocating histogram, pad the reallocation
 *                     by this much (saves excessive reallocation)
 */
struct histogram_s *
AllocHistogram(int min, int max, int lumpsize)
{
  struct histogram_s *h;
  int            newsize;
  int            i;

  newsize = max - min + 1;

  h = (struct histogram_s *) MallocOrDie(sizeof(struct histogram_s));
  h->min       = min;
  h->max       = max;
  h->total     = 0;
  h->lowscore  = INT_MAX;
  h->highscore = INT_MIN;
  h->lumpsize  = lumpsize;
  h->histogram = (int *) MallocOrDie (sizeof(int) * newsize);
  for (i = 0; i < newsize; i++) h->histogram[i] = 0;

  h->expect    = NULL;
  h->fit_type  = HISTFIT_NONE;

  return h;
}


/* Function: FreeHistogram()
 * 
 * Purpose:  free a histogram structure.
 */
void
FreeHistogram(struct histogram_s *h)
{
  free(h->histogram);
  if (h->expect != NULL) free(h->expect);
  free(h);
} 

/* Function: UnfitHistogram()
 * 
 * Purpose:  Free only the theoretical fit part of a histogram.
 */
void
UnfitHistogram(struct histogram_s *h)
{
  if (h->expect != NULL) free(h->expect);
  h->expect   = NULL;
  h->fit_type = HISTFIT_NONE;
}


/* Function: AddToHistogram()
 * 
 * Purpose:  Bump the appropriate counter in a histogram
 *           structure, given a score. The score is
 *           rounded off from float precision to the
 *           next lower integer.
 */
void
AddToHistogram(struct histogram_s *h, float sc)
{
  int score;
  int moveby;
  int prevsize;
  int newsize;
  int i;

  /* Adding to a histogram conflicts with existing fit:
   * prohibit this.
   */
  if (h->fit_type != HISTFIT_NONE)
    Die("AddToHistogram(): Can't add to a fitted histogram\n");
  

  /* histogram bins are defined as:  score >= bin value, < bin+1 
   * -1.9 -> -2    -0.4 -> -1    1.9 -> 1
   * -2.1 -> -3     0.4 -> 0     2.1 -> 2
   */
  score = (int) floor(sc);

  /* Check to see if we must reallocate the histogram.
   */
  if (score < h->min)
    {
      prevsize = h->max - h->min + 1;
      moveby   = (h->min - score) + h->lumpsize;
      newsize  = prevsize + moveby;
      h->min  -= moveby;

      h->histogram = (int *) ReallocOrDie(h->histogram, sizeof(int) * newsize);
      memmove(h->histogram+moveby, h->histogram, sizeof(int) * prevsize);
      for (i = 0; i < moveby; i++)
	h->histogram[i] = 0;
    }
  else if (score > h->max)
    {
      prevsize = h->max - h->min + 1;
      h->max   = h->lumpsize + score;
      newsize  = h->max - h->min + 1;

      h->histogram = (int *) ReallocOrDie(h->histogram, sizeof(int) * newsize);
      for (i = prevsize; i < newsize; i++)
	h->histogram[i] = 0;
    }

  /* Bump the correct bin.
   * The bin number is score - h->min
   */

  h->histogram[score - h->min]++;
  h->total++;
  if (score < h->lowscore) h->lowscore   = score;
  if (score > h->highscore) h->highscore = score;

#if DEBUGLEVEL >= DEBUG_AGGRESSIVE
  fprintf(stderr, "AddToHistogram(): added %.1f; rounded to %d; in bin %d (%d-%d)\n",
	  sc, score, score-h->min, h->min, h->max);
#endif
  return;
}


/* Function: PrintASCIIHistogram()
 * 
 * Purpose:  Print a "prettified" histogram to a file pointer.
 *           Deliberately a clone of Bill Pearson's excellent FASTA output.
 * 
 * Args:     fp     - open file to print to (stdout works)
 *           h      - histogram to print
 */           
void
PrintASCIIHistogram(FILE *fp, struct histogram_s *h)
{
  int units;
  int maxbar;
  int num;
  int i, idx;
  char buffer[81];		/* output line buffer */
  int  pos;			/* position in output line buffer */

  /* Find out how we'll scale the histogram.
   * We have 59 characters to play with on a
   * standard 80-column terminal display:
   * leading "%5d %6d %6d|" occupies 20 chars.
   */
  maxbar = 0;
  for (i = h->lowscore - h->min; i <= h->highscore - h->min; i++)
    if (h->histogram[i] > maxbar) maxbar = h->histogram[i];
  units = ((maxbar-1)/ 59) + 1;

  /* Print the histogram
   */
  fprintf(fp, "%5s %6s %6s  (one = represents %d sequences)\n", 
	  "score", "obs", "exp", units);
  buffer[80] = '\0';
  buffer[79] = '\n';
  for (i = h->lowscore; i <= h->highscore; i++)
    {
      memset(buffer, ' ', 79 * sizeof(char));
      idx = i - h->min;

      if (h->fit_type != HISTFIT_NONE) 
	sprintf(buffer, "%5d %6d %6d|", 
		i, h->histogram[idx], (int) h->expect[idx]);
      else
	sprintf(buffer, "%5d %6d %6s|", 
		i, h->histogram[idx], "-");
      buffer[20] = ' ';		/* sprintf writes a null char */

      /* Mark the histogram bar for observed hits
       */ 
      if (h->histogram[idx] > 0) {
	num = 1 + (h->histogram[idx]-1) / units;
	for (pos = 20; num > 0; num--)  buffer[pos++] = '=';
      }
	  
      /* Mark the theoretically expected value
       */

      if (h->fit_type != HISTFIT_NONE && (int) h->expect[idx] > 0)
	{
	  pos = 20 + (int)(h->expect[idx]-1) / units;
	  if (pos >= 78) pos = 78; /* be careful of buffer bounds */
	  buffer[pos] = '*';
	}

      /* Print the line
       */
      fputs(buffer, fp);
    }

  /* Print details about the statistics
   */
  switch (h->fit_type) {
  case HISTFIT_NONE:
    fprintf(fp, "\n\n%% No statistical fit available\n");
    break;
    
  case HISTFIT_EVD:
    fprintf(fp, "\n\n%% Statistical details of theoretical EVD fit:\n");
    fprintf(fp, "              mu = %10.4f\n", h->param[EVD_MU]);
    fprintf(fp, "          lambda = %10.4f\n", h->param[EVD_LAMBDA]);
    fprintf(fp, "chi-sq statistic = %10.4f\n", h->chisq);
    fprintf(fp, "  P(chi-square)  = %10.4g\n", h->chip);
    break;

  case HISTFIT_GAUSSIAN:
    fprintf(fp, "\n\n%% Statistical details of theoretical Gaussian fit:\n");
    fprintf(fp, "            mean = %10.4f\n", h->param[GAUSS_MEAN]);
    fprintf(fp, "              sd = %10.4f\n", h->param[GAUSS_SD]);
    fprintf(fp, "chi-sq statistic = %10.4f\n", h->chisq);
    fprintf(fp, "  P(chi-square)  = %10.4g\n", h->chip);
    break;
  }    
  return;
}
  

/* Function: ExtremeValueFitHistogram()
 * 
 * Purpose:  Fit a score histogram to the extreme value 
 *           distribution. Set the parameters lambda
 *           and mu in the histogram structure. Calculate
 *           a chi-square test as a measure of goodness of fit. 
 *           
 * Methods:  Uses a linear regression fitting method [Collins88].
 *           Fits only the upper side of the distribution, starting
 *           from the noise peak, since the lower outliers are difficult
 *           to remove and screw up the EVD fit.
 *           Upper outliers are removed using method described by [Mott92].
 *           
 *           A Bayesian fit would be more consistent, religiously
 *           speaking. Likewise, we might consider a Kullback log-ratio
 *           test instead of the traditional chi-square test.
 *
 * Args:     h         - histogram to fit
 *           high_hint - score cutoff; above this are `real' hits that aren't fit
 *           
 * Return:   1 if fit is judged to be valid.
 *           else 0 if fit is invalid (too few seqs.)
 */
int
ExtremeValueFitHistogram(struct histogram_s *h, float high_hint) 
{
  float *d;            /* distribution P(S < x)          */
  float *x;            /* x-axis of P(S<x) for Linefit() */
  float *v;            /* made-up variances              */
  int    hsize;
  int    sum;
  int    max;
  int    sc, idx;		/* loop indices for score or score-h->min*/
  int    lowbound, highbound;
  int    new_highbound;
  float  slope, intercept;	/* m,b fit from WeightedLinefit()        */
  int    nbins;			/* number of bins counted towards chi-sq */
  float  delta;			/* obs - exp difference                  */
  int    offset;
  int    max_iterations = 100;

  /* Clear any previous fitting from the histogram.
   */
  UnfitHistogram(h);

  /* Determine if we have enough hits to fit the histogram;
   * arbitrarily require 1000.
   */
  if (h->total < 1000) { h->fit_type = HISTFIT_NONE; return 0; }

  /* Calculate P(S < x) distribution from histogram.
   * distribution d runs from min..max with indices 0..max-min
   *     i.e. score - min = index into d, x, v, histogram, and expect
   */
  hsize = h->max - h->min + 1;
  d         = (float *) MallocOrDie(sizeof(float) * hsize);
  x         = (float *) MallocOrDie(sizeof(float) * hsize);
  v         = (float *) MallocOrDie(sizeof(float) * hsize);
  h->expect = (float *) MallocOrDie(sizeof(float) * hsize);
  for (idx = 0; idx < hsize; idx++)
    d[idx] = x[idx] = v[idx] = h->expect[idx] = 0.;

  sum = 0;
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    {
      d[sc - h->min] = (float) sum / (float) h->total;
      x[sc - h->min] = (float) sc;
      v[sc - h->min] = 1. / (float) (1.+h->histogram[sc-h->min]); /* WAG */
      sum  += h->histogram[sc - h->min];
    }

  /* Do a linear regression fit to the log[-log(P(S<x))] "line".
   * Leave off the first point, which is P(S<x)=0, log of which is undefined.
   * Linefit() does y = mx + b
   * we have log[-log(1-P(S>x))]  = -lambda * x + lambda * mu
   * so lambda = -m  and mu = b/lambda
   * 
   * Only part of the line (between lowbound and highbound) is
   * fitted.
   */
#if DEBUGLEVEL >= DEBUG_AGGRESSIVE
  fprintf(stderr, "## Here is the distribution P(S<x) observed\n");
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    fprintf(stderr, "%f %f\n", x[sc-h->min], d[sc-h->min]);
#endif

				/* convert y axis to log[-log(P(S<x))]  */
  for (sc = h->lowscore + 1; sc <= h->highscore; sc++)
    d[sc - h->min] = log(-1. * log(d[sc - h->min]));
#if DEBUGLEVEL >= DEBUG_AGGRESSIVE
  fprintf(stderr, "## Here is the line y = log[-log P(S<x)], observed\n");
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    fprintf(stderr, "%f %f\n", x[sc-h->min], d[sc-h->min]);
#endif

  /* find lowbound as peak of distribution 
   */
  lowbound  = 0;
  max       = 0;
  for (sc = h->lowscore + 1; sc <= h->highscore; sc++)
    if (h->histogram[sc - h->min] > max) 
      { 
	max      = h->histogram[sc - h->min];
	lowbound = sc;
      }

  /* initial highbound is the lower of high_hint or right edge 
   */
  if ((int) floor(high_hint) < h->highscore) 
    highbound = (int) floor(high_hint);
  else
    highbound = h->highscore;

  while (max_iterations--)	/* while (1) might risk an infinite loop */
    {
      offset = lowbound - h->min;
      hsize  = highbound - lowbound + 1;

		/* Check that we have enough sequences within bounds */
      sum = 0;
      for (sc = lowbound; sc <= highbound; sc++)
	sum += h->histogram[sc - h->min];
      if (sum < 1000) { h->fit_type = HISTFIT_NONE; return 0; }

		/* call the line fit and solve for lambda, mu */
      WeightedLinefit(x+offset,d+offset,v+offset, hsize, &slope, &intercept);
      h->param[EVD_LAMBDA] = slope * -1.;
      h->param[EVD_MU]     = intercept / h->param[EVD_LAMBDA];

				/* calculate new high bound [Mott92] */
      new_highbound = (int) ceil(h->param[EVD_MU] - (log(-1.*log((float)h->total/((float)h->total+1.))))/(h->param[EVD_LAMBDA]));
      if (new_highbound > h->highscore) new_highbound = h->highscore;

				/* check for convergence */
      if (new_highbound == highbound) break;
      highbound = new_highbound;
    }

#if DEBUGLEVEL >= DEBUG_AGGRESSIVE
  fprintf(stderr, "## Here is the line y = -lambda*x + lambda*mu, fitted\n");
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    fprintf(stderr, "%f %f\n", x[sc-h->min], 
	    -1 *h->param[EVD_LAMBDA]*((float)sc - h->param[EVD_MU]));
  evd_goodness(h, lowbound, highbound, h->param[EVD_MU], h->param[EVD_LAMBDA]);
  printf("Before the amoeba: %f %f %f\n", 
	 h->param[EVD_MU], h->param[EVD_LAMBDA], h->chisq);

#endif

  /* Optimize the fit by a brute force grid search.
   * This code is commented out in production code;
   * only used during EVD fit testing to produce highly
   * optimized fits. Since it doesn't appear in production
   * code, there's no reason to use fancier optimization algorithms.
   */
  /* evd_bruteforce(h, lowbound, highbound);  */
  
  /* Calculate the expected values for the histogram.
   */
  h->fit_type = HISTFIT_EVD;
  for (sc = h->min; sc <= h->max; sc++)
    h->expect[sc - h->min] =
      ExtremeValueE((float)(sc), h->param[EVD_MU], h->param[EVD_LAMBDA], 
		    h->total) -
      ExtremeValueE((float)(sc+1), h->param[EVD_MU], h->param[EVD_LAMBDA],
		    h->total);

  /* Calculate the goodness-of-fit (within region that was fitted)
   */
  h->chisq = 0.;
  nbins    = 0;
  for (sc = lowbound; sc <= highbound; sc++)
    if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
      {
	delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
	h->chisq += delta * delta / h->expect[sc-h->min];
	nbins++;
      }
  /* Note on degrees of freedom: there are only *two* constraints,
   * mu and lambda. Because we fit only between the peak and the
   * right edge, instead of to the whole histogram, the normalization
   * of the histogram does not count as a constraint.
   */
  if (nbins > 2)
    h->chip = (float) IncompleteGamma((double)(nbins-2)/2., 
				      (double) h->chisq/2.);
  else
    h->chip = 0.;		

  free(x);
  free(d);
  free(v);
  return 1;
}

    
/* Function: ExtremeValueSetHistogram()
 * 
 * Purpose:  Instead of fitting the histogram to an EVD,
 *           simply set the EVD parameters from an external source.
 */
void
ExtremeValueSetHistogram(struct histogram_s *h, float mu, float lambda)
{
  int   sc;
  int   hsize, idx;
  int   nbins;
  float delta;

  UnfitHistogram(h);
  h->fit_type          = HISTFIT_EVD;
  h->param[EVD_LAMBDA] = lambda;
  h->param[EVD_MU]     = mu;

  hsize     = h->max - h->min + 1;
  h->expect = (float *) MallocOrDie(sizeof(float) * hsize);
  for (idx = 0; idx < hsize; idx++)
    h->expect[idx] = 0.;

  /* Calculate the expected values for the histogram.
   */
  for (sc = h->min; sc <= h->max; sc++)
    h->expect[sc - h->min] =
      ExtremeValueE((float)(sc), h->param[EVD_MU], h->param[EVD_LAMBDA], 
		    h->total) -
      ExtremeValueE((float)(sc+1), h->param[EVD_MU], h->param[EVD_LAMBDA],
		    h->total);
  
  /* Calculate the goodness-of-fit (within whole region)
   */
  h->chisq = 0.;
  nbins    = 0;
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
      {
	delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
	h->chisq += delta * delta / h->expect[sc-h->min];
	nbins++;
      }
  /* Since we fit the whole histogram, there is one constraint on chi-square:
   * the normalization to h->total.
   */
  if (nbins > 1)
    h->chip = (float) IncompleteGamma((double)(nbins-1)/2., 
				      (double) h->chisq/2.);
  else
    h->chip = 0.;		
}



/* Function: GaussianFitHistogram()
 * 
 * Purpose:  Fit a score histogram to a Gaussian distribution.
 *           Set the parameters mean and sd in the histogram
 *           structure, as well as a chi-squared test for
 *           goodness of fit.
 *
 * Args:     h         - histogram to fit
 *           high_hint - score cutoff; above this are `real' hits that aren't fit
 *           
 * Return:   1 if fit is judged to be valid.
 *           else 0 if fit is invalid (too few seqs.)           
 */
int
GaussianFitHistogram(struct histogram_s *h, float high_hint)
{
  float sum;
  float sqsum;
  float delta;
  int   sc;
  int   nbins;
  int   hsize, idx;
  
  /* Clear any previous fitting from the histogram.
   */
  UnfitHistogram(h);

  /* Determine if we have enough hits to fit the histogram;
   * arbitrarily require 1000.
   */
  if (h->total < 1000) { h->fit_type = HISTFIT_NONE; return 0; }

  /* Simplest algorithm for mean and sd;
   * no outlier detection yet (not even using high_hint)
   * 
   * Magic 0.5 correction is because our histogram is for
   * scores between x and x+1; we estimate the expectation
   * (roughly) as x + 0.5. 
   */
  sum = sqsum = 0.;
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    {
      delta  = (float) sc + 0.5;
      sum   += (float) h->histogram[sc-h->min] * delta;
      sqsum += (float) h->histogram[sc-h->min] * delta * delta;
    }
  h->fit_type          = HISTFIT_GAUSSIAN;
  h->param[GAUSS_MEAN] = sum / (float) h->total;
  h->param[GAUSS_SD]   = sqrt((sqsum - (sum*sum/(float)h->total)) / 
			      (float)(h->total-1));
  
  /* Calculate the expected values for the histogram.
   * Note that the magic 0.5 correction appears again.
   * Calculating difference between distribution functions for Gaussian 
   * would be correct but hard.
   */
  hsize     = h->max - h->min + 1;
  h->expect = (float *) MallocOrDie(sizeof(float) * hsize);
  for (idx = 0; idx < hsize; idx++)
    h->expect[idx] = 0.;

  for (sc = h->min; sc <= h->max; sc++)
    {
      delta = (float) sc + 0.5 - h->param[GAUSS_MEAN];
      h->expect[sc - h->min] =
	(float) h->total * ((1. / (h->param[GAUSS_SD] * sqrt(2.*3.14159))) * 
        (exp(-1.* delta*delta / (2. * h->param[GAUSS_SD] * h->param[GAUSS_SD]))));
    }

  /* Calculate the goodness-of-fit (within region that was fitted)
   */
  h->chisq = 0.;
  nbins    = 0;
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
      {
	delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
	h->chisq += delta * delta / h->expect[sc-h->min];
	nbins++;
      }
	/* -1 d.f. for normalization; -2 d.f. for two free parameters */
  if (nbins > 3)
    h->chip = (float) IncompleteGamma((double)(nbins-3)/2., 
				      (double) h->chisq/2.);
  else
    h->chip = 0.;		

  return 1;
}


/* Function: GaussianSetHistogram()
 * 
 * Purpose:  Instead of fitting the histogram to a Gaussian,
 *           simply set the Gaussian parameters from an external source.
 */
void
GaussianSetHistogram(struct histogram_s *h, float mean, float sd)
{
  int   sc;
  int   hsize, idx;
  int   nbins;
  float delta;

  UnfitHistogram(h);
  h->fit_type          = HISTFIT_GAUSSIAN;
  h->param[GAUSS_MEAN] = mean;
  h->param[GAUSS_SD]   = sd;

  /* Calculate the expected values for the histogram.
   */
  hsize     = h->max - h->min + 1;
  h->expect = (float *) MallocOrDie(sizeof(float) * hsize);
  for (idx = 0; idx < hsize; idx++)
    h->expect[idx] = 0.;

  /* Note: ideally we'd use the Gaussian distribution function
   * to find the histogram occupancy in the window sc..sc+1. 
   * However, the distribution function is hard to calculate.
   * Instead, estimate the histogram by taking the density at sc+0.5.
   */
  for (sc = h->min; sc <= h->max; sc++)
    { 
      delta = ((float)sc + 0.5) - h->param[GAUSS_MEAN];
      h->expect[sc - h->min] =
	(float) h->total * ((1. / (h->param[GAUSS_SD] * sqrt(2.*3.14159))) * 
	    (exp(-1.*delta*delta / (2. * h->param[GAUSS_SD] * h->param[GAUSS_SD]))));
    }

  /* Calculate the goodness-of-fit (within whole region)
   */
  h->chisq = 0.;
  nbins    = 0;
  for (sc = h->lowscore; sc <= h->highscore; sc++)
    if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
      {
	delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
	h->chisq += delta * delta / h->expect[sc-h->min];
	nbins++;
      }
	/* -1 d.f. for normalization */
  if (nbins > 1)
    h->chip = (float) IncompleteGamma((double)(nbins-1)/2., 
				      (double) h->chisq/2.);
  else
    h->chip = 0.;		
}



/* Function: ExtremeValueP()
 * 
 * Purpose:  Calculate P(S>x) according to an extreme
 *           value distribution, given x and the parameters
 *           of the distribution (characteristic
 *           value mu, decay constant lambda).
 *           
 * Args:     x      = score
 *           mu     = characteristic value of extreme value distribution
 *           lambda = decay constant of extreme value distribution
 *           
 * Return:   P(S>x)
 */                   
double
ExtremeValueP(float x, float mu, float lambda)
{
 return (1.0 - exp(-1.0* exp(-1.0 * lambda * (x - mu))));
}


/* Function: ExtremeValueP2()
 * 
 * Purpose:  Calculate P(S>x) in a database of size N,
 *           using P(S>x) for a single sequence, according
 *           to a Poisson distribution.
 *
 * Args:     x      = score
 *           mu     = characteristic value of extreme value distribution
 *           lambda = decay constant of extreme value distribution
 *           N      = number of trials (number of sequences)
 *
 * Return:   P(S>x) for database of size N
 */
double
ExtremeValueP2(float x, float mu, float lambda, int N)
{
  return (1.0 - exp(-1.0 * N * ExtremeValueP(x,mu,lambda)));
}

/* Function: ExtremeValueE()
 * 
 * Purpose:  Calculate E(S>x) in a database of size N,
 *           using P(S>x) for a single sequence: simply np.
 *
 * Args:     x      = score
 *           mu     = characteristic value of extreme value distribution
 *           lambda = decay constant of extreme value distribution
 *           N      = number of trials (number of sequences)
 *
 * Return:   E(S>x) for database of size N
 */
double
ExtremeValueE(float x, float mu, float lambda, int N)
{
  return (double)N * ExtremeValueP(x,mu,lambda);
}


/* Function: EVDrandom()
 * 
 * Purpose:  Randomly sample an x from an EVD.
 *           Trivially done by the transformation method, since
 *           the distribution is analytical:
 *              x = \mu - \frac{\log \left[ -\log P(S<x) \right]}{\lambda}
 *           where P(S<x) is sampled uniformly on 0 < P(S<x) < 1.
 */
float
EVDrandom(float mu, float lambda)
{
  float p = 0.0;

  /* Very unlikely, but possible,
   * that sre_random() would give us exactly 0 or 1 
   */
  while (p == 0. || p == 1.) p = sre_random(); 
  return mu - log(-1. * log(p)) / lambda;
} 
 



/***********************************************************************
 *
 * Additional optimization of EVD fits -- not used in production code
 * 
 *********************************************************************** 
 */
void
evd_bruteforce(struct histogram_s *h, int lowbound, int highbound)
{
  float bestmu, bestlambda, bestchisq;
  float mu, lambda, chisq;

  bestmu     = h->param[EVD_MU];
  bestlambda = h->param[EVD_LAMBDA];
  bestchisq  = evd_goodness(h, lowbound, highbound, bestmu, bestlambda);

  for (mu = h->param[EVD_MU] - 0.5; mu <= h->param[EVD_MU] + 0.5; mu += 0.01)
    for (lambda = h->param[EVD_LAMBDA] - 0.05; lambda <= h->param[EVD_LAMBDA] + 0.05; lambda += 0.001)
      {
	chisq = evd_goodness(h, lowbound, highbound, mu, lambda);
	if (chisq < bestchisq) 
	  {
	    bestmu = mu;
	    bestlambda = lambda;
	    bestchisq = chisq;
	  }
      }

  h->param[EVD_MU]     = bestmu;
  h->param[EVD_LAMBDA] = bestlambda;
}


/* Function: evd_goodness()
 * 
 * Purpose:  Calculate the chi-square for a given histogram for
 *           a given mu and lambda. Returns the chi-square value.
 *           
 *           Warning: uses histogram->expect and ->chisq as temp space;
 *           histogram->expect must already be properly allocated.
 */
float
evd_goodness(struct histogram_s *h, int lowbound, int highbound,
	     float mu, float lambda)
{
  int   sc;			/* loop index for score (histogram bin)  */
  int   nbins;			/* number of bins included in chi-square */
  float delta;

  /* Calculate the expected values for the histogram
   * within fitted region.
   */
  for (sc = lowbound; sc <= highbound; sc++)
    {
      h->expect[sc - h->min] =
	ExtremeValueE((float)(sc), mu, lambda, h->total) -
	ExtremeValueE((float)(sc+1), mu, lambda, h->total);
      assert(h->expect[sc - h->min] > 0.);
    }

  /* Calculate the chi-square.
   */
  h->chisq = 0.;
  nbins = 0;
  for (sc = lowbound; sc <= highbound; sc++)
    if (h->expect[sc-h->min] >= 5. && h->histogram[sc-h->min] >= 5)
    {
      delta = (float) h->histogram[sc-h->min] - h->expect[sc-h->min];
      h->chisq += delta * delta / h->expect[sc-h->min];
      nbins++;
    }
  if (nbins < 3) h->chisq = 9999999.;
  return h->chisq;
}

