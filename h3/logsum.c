/* The p7_FLogsum() function used in the Forward() algorithm.
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
 * Contents:      
 * 
 * SRE, Wed Jul 11 11:00:57 2007 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <math.h>

#include "hmmer.h"

static float flogsum_lookup[p7_LOGSUM_TBL];
static int   ilogsum_lookup[p7_LOGSUM_TBL];

void
p7_FLogsumInit(void)
{
  static int firsttime = TRUE;
  if (!firsttime) return;
  firsttime = FALSE;

  int i;
  for (i = 0; i < p7_LOGSUM_TBL; i++) 
    flogsum_lookup[i] = log(1. + exp((double) -i / p7_INTSCALE));
  return;
}

/* Function:  
 * Synopsis:  
 * Incept:    SRE, Fri Jul 13 15:30:39 2007 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 *
 * Note:      This function is a critical optimization target, because
 *            it's in the inner loop of Forward() algorithms.
 */
float
p7_FLogsum(float s1, float s2)
{
#if 0
  return (log(exp(s1) + exp(s2))); /* SRE: While debugging SSE impl. Remember to remove! */
#endif
  const float max = ESL_MAX(s1, s2);
  const float min = ESL_MIN(s1, s2);
  return  (min == -eslINFINITY || (max-min) >= 15.7f) ? max : max + flogsum_lookup[(int)((max-min)*p7_INTSCALE)];

} 


void 
p7_ILogsumInit(void)
{
  static int firsttime = TRUE;
  if (!firsttime)  return;
  firsttime = FALSE;
    
  int i;
  for (i = 0; i < p7_LOGSUM_TBL; i++) 
    ilogsum_lookup[i] = rint(p7_INTSCALE * (log(1.+exp((double) -i/p7_INTSCALE))));
}


int 
p7_ILogsum(int s1, int s2)
{
  const int max = ESL_MAX(p7_IMPOSSIBLE, ESL_MAX(s1, s2));
  const int min = ESL_MIN(s1, s2);
  return  (min <= p7_IMPOSSIBLE || (max-min) >= p7_LOGSUM_TBL) ? max : max + ilogsum_lookup[max-min];
} 



/*****************************************************************
 * Benchmark driver.
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



