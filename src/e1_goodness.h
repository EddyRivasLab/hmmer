/* goodness - funtions for goodness of fit test of 
 *            a model to a geometric distribution of insert lengths
 */
#ifndef E1GOODNESS_INCLUDED
#define E1GOODNESS_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_histogram.h"

extern int GoodnessGeometric(FILE *fp, ESL_HISTOGRAM *h, double *eta, char *errbuf, int verbose);
extern int geometric_Fit(double *xv, int n, double *ret_b);

#endif /*E1GOODNESS_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
