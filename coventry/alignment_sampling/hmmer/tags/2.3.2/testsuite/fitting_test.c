/* fitting_test.c
 * 17 June 1997 (see notebook)
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "structs.h"
#include "funcs.h"
#include "squid.h"

#include "globals.h"

int
main(int argc, char **argv)
{
  int   n;			/* number of EVD samples */
  float p1, p2;
  struct histogram_s *histog;   
  int   i,j;
  float x;
  int   seed;
  int   do_evd, set, fit_evd, show_hist;

  p1 = atof(argv[1]);		/* mu or mean */
  p2 = atof(argv[2]);		/* lambda or sd */
  n  = atoi(argv[3]);		/* # of histograms */
  do_evd = atoi(argv[4]);	/* 1 to sample EVD; 0 to sample Gaussian */
  set = atoi(argv[5]);		/* 1 to set instead of fit the dist  */
  fit_evd = atoi(argv[6]);	/* 1 to fit EVD; 0 to fit Gaussian  */
  show_hist = atoi(argv[7]);	/* 1 to show histogram */

  seed = (int) time ((time_t *) NULL); 
  sre_srandom(seed);

  for (j = 0; j < n; j++)
    {
      histog = AllocHistogram(-200, 200, 100);
      for (i = 0; i < 2500; i++)
	{
	  if (do_evd) x = EVDrandom(p1, p2);
	  else        x = Gaussrandom(p1, p2);

	  assert(x > -100.);
	  assert(x < 100.);
	  AddToHistogram(histog, x);
	}

      if (set && fit_evd)
	ExtremeValueSetHistogram(histog, p1, p2);
      else if (set && !fit_evd)
	GaussianSetHistogram(histog, p1, p2);
      else if (!set && fit_evd)
	ExtremeValueFitHistogram(histog, 9999.);
      else 
	GaussianFitHistogram(histog, 9999.);

      printf("%f\n", histog->chip);

      if (show_hist)
	PrintASCIIHistogram(stdout, histog); 
     
      FreeHistogram(histog);
    }

  return 0;
}
