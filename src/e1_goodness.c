/* e1_goodness - funtions for goodness of fit test of 
 *               a model to a geometric distribution of insert lengths
 * 
 * Contents:
 *
 * ER, Mon Sep  2 09:05:36 EDT 2013 [Janelia] 
 * SVN $Id:$
 */

#include <string.h>
#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_random.h"
#include "esl_ratematrix.h"
#include "esl_rootfinder.h"
#include "esl_stats.h"
#include "esl_vectorops.h"

#include "e2.h"
#include "e1_goodness.h"
#include "e1_simulate.h"
#include "ratematrix.h"

static double geometric_cdf(double x, double p);
static double geometric_generic_cdf(double x, void *params);

int
GoodnessGeometric(FILE *fp, ESL_HISTOGRAM *h, double *eta, char *errbuf, int verbose)
{
  double  b;
  double  bfit;
  double  params[1];
  double  G, Gp, X2, X2p;
  double *x;
  int     n;
  
  if (fp == NULL) return eslOK;
  
  esl_histogram_GetData(h, &x, &n);
  geometric_Fit(x, n, &bfit);
  if (eta) b = *eta;
  else     b =  bfit;
 
  params[0] = b;
  esl_histogram_SetExpect(h, geometric_generic_cdf, &params);
  
  printf("h->n %d\n", (int)h->n);
  if (h->n > 0) {
    esl_histogram_Write(stdout, h);
    esl_histogram_Plot(fp, h);
    esl_histogram_Goodness(h, 0, NULL, &G, &Gp, &X2, &X2p);
    if (eta) printf("bernoulli: given* = %f fit %f\n", b, bfit);
    else     printf("bernoulli: fit* %f\n", b);
    printf("G   = %f  p = %f\n", G, Gp);
    printf("X^2 = %f  p = %f\n", X2, X2p);
  }
  
  return eslOK;
}

int 
geometric_Fit(double *x, int n, double *ret_b)
{
  double b;
  double meanL;
 
  esl_stats_DMean(x, n, &meanL, NULL);

  b = meanL/(1.0+meanL);
  *ret_b = b;

  return eslOK;
}


static double 
geometric_cdf(double x, double p)
{
  double y;
  y = log(p);
  return 1.0 - exp( (x+1.0)*y );
}

/* Function:  geometric_generic_cdf()
 *
 * Purpose:   Generic-API version of CDF function.
 */
static double
geometric_generic_cdf(double x, void *params) 
{
  double *p = (double *) params;
  return geometric_cdf(x, p[0]);
}



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
