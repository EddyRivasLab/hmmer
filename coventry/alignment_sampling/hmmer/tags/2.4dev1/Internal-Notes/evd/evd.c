#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int
main(int argc, char **argv)
{
  double mu;
  double lambda;
  double x;
  double px;
  double psx;

  mu     = atof(argv[1]);
  lambda = atof(argv[2]);

  for (x = -5.; x <= 20.; x += 0.1)
    {
				/* density function P(S=x) */
      px = lambda * exp(-1.*lambda * (x - mu) - exp(-1* lambda * (x - mu)));
				/* distribution function P(S < x) */
      psx= exp(-1* exp(-1*lambda * (x - mu)));

      printf("%f %f\n", x, px);

    }
  printf("&\n");
}
