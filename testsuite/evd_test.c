/* evd_test.c
 * SRE, Wed Nov 12 11:17:27 1997 [St. Louis]
 * cp evd_test.c ../src/testdriver.c; cd ../src; make testdriver
 * 
 * Test driver for EVD distribution support in histogram.c
 * Generates random EVD samples; fits them; checks fitted mu, lambda
 * against parametric mu, lambda. If they differ badly, calls Die(). 
 * If OK, returns EXIT_SUCCESS.
 * 
 * RCS $Id$
 */


#include <stdio.h>
#include <time.h>
#include <math.h>

#include "structs.h"
#include "funcs.h"
#include "globals.h"
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static char banner[] = "\
evd_test : testing of EVD code in histogram.c";

static char usage[] = "\
Usage: testdriver [-options]\n\
  Available options are:\n\
  -h              : help; display this usage info\n\
  -c <x>          : censor data below <x>\n\
  -e <n>          : sample <n> times from EVD\n\
  -g <n>          : add <n> Gaussian samples of \"noise\"\n\
  -n <n>          : set number of trials to <n>\n\
  -s <n>          : set random seed to <n>\n\
  -v              : be verbose (default is to simply exit with status 1 or 0)\n\
  --xmgr <file>   : save graphical data to <file>\n\
  --hist          : fit to histogram instead of raw samples\n\
  --loglog <file> : save log log regression line to <file>\n\
  --regress       : do old-style linear regression fit, not ML\n\
  --mu <x>        : set EVD mu to <x>\n\
  --lambda <x>    : set EVD lambda to <x>\n\
  --mean <x>      : set Gaussian mean to <x>\n\
  --sd   <x>      : set Gaussian std. dev. to <x>\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-h",       TRUE,  sqdARG_NONE  },
  { "-c",       TRUE,  sqdARG_FLOAT },
  { "-e",       TRUE,  sqdARG_INT },
  { "-g",       TRUE,  sqdARG_INT },
  { "-n",       TRUE,  sqdARG_INT },
  { "-s",       TRUE,  sqdARG_INT   }, 
  { "-v",       TRUE,  sqdARG_NONE  },
  { "--xmgr",   FALSE, sqdARG_STRING},
  { "--hist",   FALSE, sqdARG_NONE},
  { "--loglog", FALSE, sqdARG_STRING},
  { "--regress",FALSE, sqdARG_NONE},
  { "--mu",     FALSE, sqdARG_FLOAT},
  { "--lambda", FALSE, sqdARG_FLOAT},
  { "--mean",   FALSE, sqdARG_FLOAT},
  { "--sd",     FALSE, sqdARG_FLOAT},
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  struct histogram_s *h;        /* histogram structure          */    
  int ntrials;			/* number of different fits     */
  int be_verbose;               /* option: TRUE to show output  */
  int seed;                     /* option: random number seed   */
  int   nevd;                   /* # of samples from EVD        */
  float mu;			/* EVD mu parameter             */
  float lambda;			/* EVD lambda parameter         */
  int   ngauss;			/* # of samples from Gaussian   */
  float mean;			/* Gaussian "noise" mean        */
  float sd;			/* Gaussian "noise" std. dev.   */
  float x;			/* a random sample              */
  int   i, idx;
  float *val;			/* array of samples             */
  float mlmu;			/* estimate of mu               */
  float mllambda;		/* estimate of lambda           */

  char *xmgrfile;               /* output file for XMGR graph data */
  char *logfile;                /* output file for regression line */
  FILE *xmgrfp;                 /* open output file                */
  FILE *logfp;                  /* open log log file               */
  int   do_ml;			/* TRUE to do a max likelihood fit */
  int   fit_hist;		/* TRUE to fit histogram instead of samples */
  int   censoring;		/* TRUE to left-censor the data    */
  float censorlevel;		/* value to censor at              */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_inuse(&histid1);
  fprintf(stderr, "[... memory debugging is ON ...]\n");
#endif
  
  /*********************************************** 
   * Parse command line
   ***********************************************/
  be_verbose = FALSE;
  seed       = (int) time ((time_t *) NULL);
  ntrials    = 1;
  nevd       = 1000;
  mu         = -20.0;
  lambda     = 0.4;
  ngauss     = 0;
  mean       = 20.;
  sd         = 20.;
  xmgrfile   = NULL;
  logfile    = NULL;
  xmgrfp     = NULL;
  logfp      = NULL;
  do_ml      = TRUE;
  censoring  = FALSE;
  censorlevel= 0.;
  fit_hist   = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-e")       == 0) { nevd       = atoi(optarg); } 
    else if (strcmp(optname, "-c")       == 0) { censoring  = TRUE;
                                                 censorlevel= atof(optarg); } 
    else if (strcmp(optname, "-g")       == 0) { ngauss     = atoi(optarg); } 
    else if (strcmp(optname, "-n")       == 0) { ntrials    = atoi(optarg); }
    else if (strcmp(optname, "-s")       == 0) { seed       = atoi(optarg); }
    else if (strcmp(optname, "-v")       == 0) { be_verbose = TRUE;         }
    else if (strcmp(optname, "--xmgr")   == 0) { xmgrfile   = optarg; }
    else if (strcmp(optname, "--hist")   == 0) { fit_hist   = TRUE; }
    else if (strcmp(optname, "--loglog") == 0) { logfile    = optarg; }
    else if (strcmp(optname, "--regress")== 0) { do_ml      = FALSE; }
    else if (strcmp(optname, "--mu")     == 0) { mu         = atof(optarg); } 
    else if (strcmp(optname, "--lambda") == 0) { lambda     = atof(optarg); } 
    else if (strcmp(optname, "--mean")   == 0) { mean       = atof(optarg); } 
    else if (strcmp(optname, "--sd")     == 0) { sd         = atof(optarg); } 
    else if (strcmp(optname, "-h")       == 0) {
      Banner(stdout, banner);
      puts(usage);
      exit(0);
    }
  }
  if (argc - optind != 0)
    Die("Incorrect number of arguments.\n%s\n", usage);

  sre_srandom(seed);

  /****************************************************************
   * Print options
   ****************************************************************/

  if (be_verbose)  
    {
      puts("--------------------------------------------------------");
      printf("EVD samples    = %d\n", nevd);
      printf("mu, lambda     = %f, %f\n", mu, lambda);
      if (ngauss > 0) {
	printf("Gaussian noise = %d\n", ngauss);
	printf("mean, sd       = %f, %f\n", mean, sd); 
      }
      if (censoring) printf("pre-censoring  = ON, at %f\n", censorlevel);
      printf("total trials   = %d\n", ntrials);
      printf("random seed    = %d\n", seed);
      printf("fit method     = %s\n", do_ml ? "ML" : "linear regression");
      printf("fit is to      = %s\n", fit_hist ? "histogram" : "list");
      puts("--------------------------------------------------------");
    }
  
  if (xmgrfile != NULL) 
    if ((xmgrfp = fopen(xmgrfile, "w")) == NULL)
      Die("Failed to open output file %s", xmgrfile);
  if (logfile != NULL) 
    if ((logfp = fopen(logfile, "w")) == NULL)
      Die("Failed to open output file %s", logfile);

  /* Generate random EVD "signal" (and Gaussian "noise") 
   * samples and put them in the histogram
   */
  while (ntrials--)
    {
      val = MallocOrDie(sizeof(double) * (nevd+ngauss));
      h   = AllocHistogram(-20, 20, 10);

				/* EVD signal */
      idx = 0;
      for (i = 0; i < nevd; i++)
	{
	  x = EVDrandom(mu, lambda);
	  if (! censoring || x > censorlevel)
	    {
	      AddToHistogram(h, x);
	      val[idx] = x;
	      idx++;
	    }
	}
				/* Gaussian noise */
      for (; i < nevd + ngauss; i++)
	{
	  x = Gaussrandom(mean, sd);
	  if (! censoring || x > censorlevel)
	    {
	      AddToHistogram(h, x);
	      val[idx] = x;
	      idx++;
	    }
	}

      if (do_ml)
	{

	  if (censoring)
	    {
	      if (be_verbose)
		printf("I have censored the data at %f: %d observed, %d censored\n", censorlevel, idx, (nevd+ngauss)-idx);

	      EVDCensoredFit(val, NULL, idx, 
			     (nevd+ngauss)-idx, censorlevel,
			     &mlmu, &mllambda); 
	      ExtremeValueSetHistogram(h, (float) mlmu, (float) mllambda, 2); 
	    }
	  else
	    {
	      if (fit_hist)
		{
		  ExtremeValueFitHistogram(h, TRUE, 20.);  
		}
	      else
		{
		  EVDMaxLikelyFit(val, NULL, idx, &mlmu, &mllambda); 
		  ExtremeValueSetHistogram(h, (float) mlmu, (float) mllambda, 2); 
		}
	    }
	}
      else
	EVDBasicFit(h);

      if (be_verbose) {
	printf("%f\tmu\n",     h->param[EVD_MU]);
	printf("%f\tlambda\n", h->param[EVD_LAMBDA]);
	printf("%f\t%% error on mu\n", 
	       fabs(100. * (h->param[EVD_MU] - mu) / mu));
	printf("%f\t%% error on lambda\n", 
	       fabs(100. * (h->param[EVD_LAMBDA] - lambda) / lambda));
	printf("%f\tchi-squared P value\n", h->chip);
      }
      if (xmgrfp != NULL) PrintXMGRHistogram(xmgrfp, h);
      /*      if (xmgrfp != NULL) PrintXMGRDistribution(xmgrfp, h); */
      if (logfp  != NULL) PrintXMGRRegressionLine(logfp, h);

      /* Generate the expected lines: sets 5,7 of xmgrfile (manually delete 4,6)
       *                              set 3 of loglogfile  (manually delete 2)
       */
      ExtremeValueSetHistogram(h, mu, lambda, 0);
      if (xmgrfp != NULL) PrintXMGRHistogram(xmgrfp, h);
      /*      if (xmgrfp != NULL) PrintXMGRDistribution(xmgrfp, h); */
      if (logfp  != NULL) PrintXMGRRegressionLine(logfp, h);

      /* Do the internal test.
       * Criterion: on a 1000 sample EVD of u = -40 and lambda = 0.4,
       * estimate u to within +/- 2 and lambda to within +/- 0.05. 
       */
      if (fabs(h->param[EVD_MU] - mu) > 2.)
	Die("evd_test: tolerance to mu exceeded (%f)", 
	    fabs(h->param[EVD_MU] - mu));
      if (fabs(h->param[EVD_LAMBDA] - lambda) > 0.05)
	Die("evd_test: tolerance to lambda exceeded (%f)",
	    fabs(h->param[EVD_LAMBDA] - lambda));

      FreeHistogram(h);
      free(val);
    }

#ifdef MEMDEBUG
  current_size = malloc_inuse(&histid2);
  if (current_size != orig_size) Die("evd_test failed memory test");
  else fprintf(stderr, "[No memory leaks.]\n");
#endif
  
  if (xmgrfp != NULL) fclose(xmgrfp);
  if (logfp != NULL)  fclose(logfp);
  return EXIT_SUCCESS;
}
