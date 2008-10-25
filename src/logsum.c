/* p7_FLogsum() function used in the Forward() algorithm.
 * 
 * Contents:
 *    1. Floating point log sum.
 *    2. Benchmark driver.
 *    3. Unit tests.
 *    4. Test driver.
 *    5. Example.
 *    6. Copyright and license information.
 *
 * Exegesis:
 * 
 * Internally, HMMER3 profile scores are in nats: floating point
 * log-odds probabilities, with the log odds taken relative to
 * background residue frequencies, and the log to the base e.
 * 
 * The Forward algorithm needs to calculate sums of probabilities.
 * Given two log probabilities s1 and s2, where s1 = \log
 * \frac{p_1}{f_1}, and s2 = \log \frac{p_2}{f_2}, we need to
 * calculate s3 = \log \frac{p_1 + p_2}{f_3}.
 * 
 * The Forward algorithm guarantees that f_1 = f_2 = f_3, because it
 * is always concerned with summing terms that describe different
 * parses of the same target sequence prefix, and the product of the
 * background frequencies for the same sequence prefix is a constant.
 * 
 * The naive solution is s3 = log(e^{s1} + e^{s2}). This requires
 * expensive calls to log() and exp().
 * 
 * A better solution is s3 = s1 + log(1 + e^{s2-s1}). s1 should be the
 * greater, so s2-s1 is negative. For sufficiently small s2 << s1,
 * e^{s2-s1} becomes less than the machine's FLT_EPSILON, and s3 ~=
 * s1. (This is at about s2-s1 < -15.9, for the typical FLT_EPSILON of
 * 1.2e-7.)
 * 
 * With some loss of accuracy, we can precalculate log(1 + e^{s2-s1})
 * for a discretized range of differences (s2-s1), and compute s3 = s1
 * + table_lookup(s2-s1). This is what HMMER's p7_FLogsum() function
 * does.
 * 
 * SRE, Wed Jul 11 11:00:57 2007 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"
#include <math.h>
#include "hmmer.h"

static float flogsum_lookup[p7_LOGSUM_TBL];

/*****************************************************************
 *= 1. floating point log sum
 *****************************************************************/

/* Function:  p7_FLogsumInit()
 * Synopsis:  Initialize the p7_Logsum() function.
 * Incept:    SRE, Thu Apr 10 08:46:23 2008 [Janelia]
 *
 * Purpose:   Initialize the lookup table for <p7_FLogsum()>. 
 *            This function must be called once before any
 *            call to <p7_FLogsum()>.
 *            
 *            The precision of the lookup table is determined
 *            by the compile-time <p7_LOGSUM_TBL> constant.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_FLogsumInit(void)
{
  static int firsttime = TRUE;
  if (!firsttime) return eslOK;
  firsttime = FALSE;

  int i;
  for (i = 0; i < p7_LOGSUM_TBL; i++) 
    flogsum_lookup[i] = log(1. + exp((double) -i / p7_INTSCALE));
  return eslOK;
}

/* Function:  p7_FLogsum()
 * Synopsis:  Approximate $\log(e^a + e^b)$.
 * Incept:    SRE, Fri Jul 13 15:30:39 2007 [Janelia]
 *
 * Purpose:   Returns a fast table-driven approximation to
 *            $\log(e^a + e^b)$.
 *            
 *            Either <a> or <b> (or both) may be $-\infty$,
 *            but neither may be $+\infty$ or <NaN>.
 *
 * Note:      This function is a critical optimization target, because
 *            it's in the inner loop of generic Forward() algorithms.
 */
float
p7_FLogsum(float a, float b)
{
  const float max = ESL_MAX(a, b);
  const float min = ESL_MIN(a, b);
#if 0
  return (min == -eslINFINITY || (max-min) >= 15.7f) ? max : max + log(1.0 + exp(min-max));  /* SRE: While debugging SSE impl. Remember to remove! */
#endif
  return (min == -eslINFINITY || (max-min) >= 15.7f) ? max : max + flogsum_lookup[(int)((max-min)*p7_INTSCALE)];
} 

/* Function:  p7_FLogsumError()
 * Synopsis:  Compute absolute error in probability from Logsum.
 * Incept:    SRE, Sun Aug  3 10:22:18 2008 [Janelia]
 *
 * Purpose:   Compute the absolute error in probability space
 *            resulting from <p7_FLogsum()>'s table lookup 
 *            approximation: approximation result - exact result.
 *                                                  
 *            This is of course computable analytically for
 *            any <a,b> given <p7_LOGSUM_TBL>; but the function
 *            is useful for some routines that want to determine
 *            if <p7_FLogsum()> has been compiled in its
 *            exact slow mode for debugging purposes. Testing
 *            <p7_FLogsumError(-0.4, -0.5) > 0.0001>
 *            for example, suffices to detect that the function
 *            is compiled in its fast approximation mode given
 *            the defaults. 
 */
float
p7_FLogsumError(float a, float b)
{
  float approx = p7_FLogsum(a,b);
  float exact  = log(exp(a) + exp(b));
  return (exp(approx) - exp(exact));
}


/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7LOGSUM_BENCHMARK
/* gcc -o benchmark -g -O2 -I. -L. -I../easel -L../easel -Dp7LOGSUM_BENCHMARK logsum.c -leasel -lm
 * ./benchmark
 */
/* All times in units of nanoseconds/iteration: cpu time * 10.
 * All times derived from 1e8 iterations (-N 100000000) unless stated.
 * All runs on my workstation, a 3.2GHz Xeon.
 * Times in brackets are difference from baseline.  
 * To get baselines, comment out the appropriate Logsum() call and recompile.
 * 
 * Floating point:   gcc -g -O2
 *                   ---------      
 *   baseline:        274.5
 *   p7_FLogsum()     293.2  [18.7]
 *  
 * Integer version:             
 *   baseline:        269.9                                       
 *   p7_ILogsum()     271.8   [1.9]
 */
#include "p7_config.h"

#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",    0 },
  { "-i",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "run the integer version",                 0 },
  { "-v",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual results",     0 },
  { "-N",        eslARG_INT,"100000000",NULL,"n>0", NULL,  NULL, NULL, "number of trials",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "benchmark driver for logsum functions()";

static float 
naive1(float s1, float s2)
{
  return log(exp(s1) + exp(s2));
}

static float 
naive2(float s1, float s2)
{
  if (s1 > s2) return s1 + log(1 + exp(s2-s1));
  else         return s2 + log(1 + exp(s1-s2));
}

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r       = esl_randomness_Create(42);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  int             N       = esl_opt_GetInteger(go, "-N");
  int             i;

  if (esl_opt_GetBoolean(go, "-i"))
    {
      int  x, z;

      p7_ILogsumInit();
      esl_stopwatch_Start(w);
      for (z = 0, i = 0; i < N; i++)
	{
	  x = z - esl_random(r) * 7000;

	  if (esl_opt_GetBoolean(go, "-v"))  
	    printf("%d %d %d \n", z, x, p7_ILogsum(x, z));

	  z = p7_ILogsum(x,z);  
	  z -= 119;
	}
      esl_stopwatch_Stop(w);
    }
  else
    {
      float  x, z;

      p7_FLogsumInit();
      esl_stopwatch_Start(w);
      for (z = 0., i = 0; i < N; i++)
	{
	  x = z - esl_random(r) * 7.;

	  if (esl_opt_GetBoolean(go, "-v"))  
	    printf("%g %g %g %g %g\n", z, x, p7_FLogsum(x, z), naive1(x,z), fabs(p7_FLogsum(x, z) - naive1(x,z)));

	  z  = p7_FLogsum(x, z);       
	  /* z = naive2(x,y); */
	  z -= 0.1187;		/* empirically balancing z near 0 */
	}
      esl_stopwatch_Stop(w);
  
    }
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7LOGSUM_BENCHMARK*/
/*-------------------- end, benchmark ---------------------------*/


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7LOGSUM_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

static void
utest_FLogsumError(ESL_GETOPTS *go, ESL_RANDOMNESS *r)
{
  int     N          = esl_opt_GetInteger(go, "-N");
  float   maxval     = esl_opt_GetReal(go, "-S");
  int     be_verbose = esl_opt_GetBoolean(go, "-v");
  float   maxerr = 0.0;
  float   avgerr = 0.0;
  int     i;
  float   a,b,result,exact,err;

  for (i = 0; i < N; i++)
    {
      a = (esl_random(r) - 0.5) * maxval * 2.; /* uniform draws on -maxval..maxval */
      b = (esl_random(r) - 0.5) * maxval * 2.; 

      exact  = log(exp(a) + exp(b));
      result = p7_FLogsum(a,b);
      err    = fabs(exact-result) / maxval;

      avgerr += err;
      maxerr = ESL_MAX(maxerr, err);

      if (be_verbose)
	printf("%8.4f %8.4f %8.4f %8.4f %8.4f\n", a, b, exact, result, err);
    }
  avgerr /= (float) N;

  if (be_verbose) {
    printf("average error = %f\n", avgerr);
    printf("max error     = %f\n", maxerr);
  }

  if (maxerr > 0.0001) esl_fatal("maximum error of %f is too high: logsum unit test fails", maxerr);
  if (avgerr > 0.0001) esl_fatal("average error of %f is too high: logsum unit test fails", avgerr);
}

static void
utest_FLogsumSpecials(void)
{
  char *msg = "logsum specials unit test failed";

  if (p7_FLogsum(0.0,          -eslINFINITY) !=          0.0) esl_fatal(msg);
  if (p7_FLogsum(-eslINFINITY,          0.0) !=          0.0) esl_fatal(msg);
  if (p7_FLogsum(-eslINFINITY, -eslINFINITY) != -eslINFINITY) esl_fatal(msg);
}
#endif /*p7LOGSUM_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/

/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
#ifdef p7LOGSUM_TESTDRIVE
/*
  gcc -o logsum_utest -msse2 -g -Wall -I. -L. -I../easel -L../easel -Dp7LOGSUM_TESTDRIVE logsum.c -leasel -lm
  ./logsum_utest
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",               0},
  {"-N",  eslARG_INT,    "1000", NULL, "n>0",NULL, NULL, NULL, "number of samples",                 0},
  {"-S",  eslARG_REAL,   "20.0", NULL, "x>0",NULL, NULL, NULL, "maximum operand value",             0},
  {"-r",  eslARG_NONE,     NULL, NULL, NULL, NULL, NULL, NULL, "use arbitrary random number seed",  0},
  {"-s",  eslARG_INT,      "42", NULL, "n>0",NULL, NULL, NULL, "random number seed",                0},
  {"-v",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show verbose output",               0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for logsum.c";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = NULL;

  if (esl_opt_GetBoolean(go, "-r")) r = esl_randomness_CreateTimeseeded();
  else                              r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  p7_FLogsumInit();

  utest_FLogsumError(go, r);
  utest_FLogsumSpecials();

  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7LOGSUM_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/


/*****************************************************************
 * 5. Example.
 *****************************************************************/
#ifdef p7LOGSUM_EXAMPLE
/* gcc -o example -g -O2 -I. -L. -I../easel -L../easel -Dp7LOGSUM_EXAMPLE logsum.c -leasel -lm
 * ./example -0.5 -0.5
 */
#include "p7_config.h"
#include "easel.h"
#include "hmmer.h"

int
main(int argc, char **argv)
{
  float a = atof(argv[1]);
  float b = atof(argv[2]);
  float result;

  p7_FLogsumInit();
  result = p7_FLogsum(a, b);
  printf("p7_FLogsum(%f,%f) = %f\n", a, b, result);

  result = log(exp(a) + exp(b));
  printf("log(e^%f + e^%f) = %f\n", a, b, result);

  printf("Absolute error in probability: %f\n", p7_FLogsumError(a,b));
  return eslOK;
}
#endif /*p7LOGSUM_EXAMPLE*/
/*--------------------- end, example ----------------------------*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/


