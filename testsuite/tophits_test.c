/* tophits_test.c
 * SRE, Tue Oct 28 08:03:10 1997 [Newton Institute, Cambridge UK]
 * 
 * Test driver for tophits.c. Returns 0 if everything is OK.
 * 
 * Options:
 *    -v   Verbose; print stuff.
 */

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "structs.h"
#include "funcs.h"
#include "globals.h"
#include "squid.h"

static char banner[] = "\
tophits_test : internal verification of tophits.c";

static char usage[] = "\
Usage: tophits_test [-options]\n\
  Available options are:\n\
    -h     : help; display this usage info\n\
    -s <x> : set random seed to <x>\n\
    -v     : be verbose (default is to simply exit with status 1 or 0)\n\
"; 

static char experts[] = "\
\n";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE },
  { "-s", TRUE, sqdARG_INT  }, 
  { "-v", TRUE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  struct tophit_s *hit;         /* hit list                     */
  int i,j;			/* counters                     */
  int nsamples;			/* option: # of random "scores" */
  int be_verbose;		/* option: TRUE to show output  */
  int seed;			/* option: random number seed   */
  int paramH;			/* option: H parameter          */
  int paramA;			/* option: A parameter          */
  double *list;                 /* list of "scores"             */
  double tmp;			/* used for swapping            */
  float score, score2;

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */

  /*********************************************** 
   * Parse command line
   ***********************************************/
  be_verbose = FALSE;
  seed       = (int) time ((time_t *) NULL);
  paramH     = 100;
  paramA     = 10;
  nsamples   = 1000;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-s") == 0) { seed       = atoi(optarg); }
    else if (strcmp(optname, "-v") == 0) { be_verbose = TRUE;         }
    else if (strcmp(optname, "-h") == 0) {
      HMMERBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }
  if (argc - optind != 0)
    Die("Incorrect number of arguments.\n%s\n", usage);

  sre_srandom(seed);
  if (be_verbose)
    printf("%d\tSEED\n", seed);

  /*********************************************** 
   * Generate three tiers of numbers:
   *    paramA - really good scores, 1000-2000
   *    paramH - good scores, 100-200
   *    nsamples - paramH - paramA: bad scores, 10-20         
   * then shuffle.
   ***********************************************/

  list = MallocOrDie (sizeof(double) * nsamples);
  for (i = 0; i < paramA; i++)
    list[i] = 1000. + 1000. * sre_random();
  for (; i < paramA + paramH; i++)
    list[i] = 100. + 100. * sre_random(); 
  for (; i < nsamples; i++)
    list[i] = 10. + 10. * sre_random();

  for (i = 0; i < nsamples; i++)
    {
      j = CHOOSE(nsamples);
      tmp = list[j];
      list[j] = list[i];
      list[i] = tmp;
    }

  if (be_verbose)
    for (i = 0; i < nsamples; i++)
      printf("%8.2f\tTest set\n", list[i]);

  /*********************************************** 
   * Test of FullSortTophits().
   *     Fill up a hit list with random numbers;
   *     FullSort it;
   *     check that all top H are >= 100 and sorted.
   ***********************************************/

  hit = AllocTophits(100);
  for (i = 0; i < nsamples; i++)
    RegisterHit(hit, list[i], 0., (float) list[i], 0., 0., 
		NULL, NULL, NULL, /* name, acc, desc */
		0,0,0,
		0,0,0,
		0,0,
		NULL);
  FullSortTophits(hit);

  if (be_verbose)
    {
      for (i = 0; i < hit->num; i++)
	{
	  GetRankedHit(hit, i,  NULL, &score, NULL, NULL, 
		       NULL, NULL, NULL, /* name, acc, desc */
		       NULL, NULL, NULL, 
		       NULL, NULL, NULL, 
		       NULL, NULL,
		       NULL);
	  printf("%8.2f  FullSort()\n", score);
	}
    }

  for (i = 0; i < hit->num-1; i++)
    {
      GetRankedHit(hit, i,  NULL, &score, NULL, NULL, 
		   NULL, NULL, NULL, /* name, acc, desc */
		   NULL, NULL, NULL, 
		   NULL, NULL, NULL,
		   NULL, NULL,
		   NULL);
      GetRankedHit(hit, i+1,NULL, &score2,NULL, NULL, 
		   NULL, NULL, NULL, /* name, acc, desc */
		   NULL, NULL, NULL, 
		   NULL, NULL, NULL, 
		   NULL, NULL,
		   NULL);
      if (score < score2)
	Die("FullSortTophits() fails test: order wrong");
      if (i < paramA && score < 1000.)
	Die("FullSortTophits() fails test: lost a number");
      if (i < paramA + paramH && score < 100.)
	Die("FullSortTophits() fails test: lost a number");
    }

  FreeTophits(hit);
  free(list);
  SqdClean();
  
  if (be_verbose) printf("tophits_test is OK\n");
  return 0;
}
