/* h4_logsum() function used in log-space Forward/Backward algorithms.
 * 
 * Contents:
 *    1. Floating point log-sum-exp-2arg: LSE2
 *    2. Debugging/development tools
 *    3. Benchmark driver
 *    4. Unit tests
 *    5. Test driver
 *    6. Example
 */
#include "h4_config.h"

#include <math.h>

#include "easel.h"

/* h4LOGSUM_SCALE defines the precision of the calculation; the
 * default of 500.0 means discretizing in bins of 0.002 bits.
 * h4LOGSUM_TBL defines the size of the lookup table; the
 * default of 11500 means entries are calculated for differences of 0
 * to 23.000 bits (when h4LOGSUM_SCALE is 500.0).  2^{-h4LOGSUM_TBL /
 * h4LOGSUM_SCALE} should be on the order of the machine FLT_EPSILON,
 * 1.2e-7 for IEEE754 floating point.
 */
#define h4LOGSUM_SCALE 500.0f
#define h4LOGSUM_TBL   11500

static float h4_flogsum_lookup[h4LOGSUM_TBL];  // h4LOGSUM_TBL=11500: (A-B) = 0..23 bits, steps of 0.002
static int   h4_logsum_initialized = FALSE;    // A flag to allow us to crash out of logsum() if lookup table wasn't initialized
static int   h4_logsum_max         = FALSE;    // Some debugging tests force logsum() to do max(), and we need a flag to get the slow/exact mode to do it

/*****************************************************************
 *# 1. floating point log-sum-exp-2arg: LSE2
 *****************************************************************/

/* Function:  h4_logsum_Init()
 * Synopsis:  Initialize the h4_logsum() lookup table
 *
 * Purpose:   Initialize the lookup table for <h4_logsum()>. 
 *            This function must be called once before any
 *            call to <h4_logsum()>.
 *            
 *            The precision of the lookup table is determined
 *            by the compile-time <h4LOGSUM_TBL> constant.
 *
 * Returns:   <eslOK> on success.
 */
int
h4_logsum_Init(void)
{
  int i;

  if (h4_logsum_initialized) return eslOK;

  // I used to add +0.5 to i to pre-round the table and reduce
  // discretization error by 2x; but then A+A gives maximum error,
  // violating principle of least surprise.
  for (i = 0; i < h4LOGSUM_TBL; i++) 
    h4_flogsum_lookup[i] = (float) log2(1. + exp2(- ((double) i) / h4LOGSUM_SCALE)); 
  h4_logsum_initialized = TRUE;
  h4_logsum_max         = FALSE;
  return eslOK;
}

/* Function:  h4_logsum_Reinit()
 * Synopsis:  Reinitialize the h4_logsum() lookup table
 *
 * Purpose:   Force a reinitialization of the lookup table.  Calling
 *            <h4_logsum_Init()> a second time has no effect, because
 *            of an optimization that checks the
 *            <h4_logsum_initialized> flag to avoid wasting time
 *            recalculating the table. When a caller is running a unit
 *            test that uses <h4_logsum_InitMax()> to hack the table
 *            to convert F/B sums to Viterbi max's, it wants to reset
 *            the table when it's done, using this function.
 */
int
h4_logsum_Reinit(void)
{
  h4_logsum_initialized = FALSE;
  return h4_logsum_Init();
}


/* Function:  h4_logsum()
 * Synopsis:  Approximate $\log_2(2^a + 2^b)$.
 *
 * Purpose:   Returns a fast table-driven approximation to
 *            $\log_2(2^a + 2^b)$.
 *            
 *            Either <a> or <b> (or both) may be $-\infty$,
 *            but neither may be $+\infty$ or <NaN>.
 *
 * Note:      This function is a critical optimization target, because it's
 *            in the inner loop of sparse and reference Forward/Backward
 *            algorithms.
 *            
 *            Compiling with the <h4LOGSUM_SLOWEXACT> flag bypasses
 *            the table-driven approximation and uses the exact 
 *            calculation instead; useful for debugging.
 *            
 *            Because this is in critical path, we don't test
 *            h4_logsum_max, like we do in slow/exact mode. Instead,
 *            debugging tools that want this function to yield max()
 *            instead of logsum() use h4_logsum_InitMax() to
 *            initialize the lookup table to all 0. This way, we don't
 *            put any other instructions in the critical path just to
 *            get that rarely-used debugging functionality.
 */
float
h4_logsum(float a, float b)
{
  const float max = ESL_MAX(a, b);
  const float min = ESL_MIN(a, b);      // needs to work when a or b is -inf, so fabs(a-b) or abs((int) a-b) don't work.

  ESL_DASSERT1(( h4_logsum_initialized ));

#ifdef h4LOGSUM_SLOWEXACT
  return (h4_logsum_max || min == -eslINFINITY || (max-min) >= 23.0f) ? max : max + log2f(1.0 + exp2f(min-max));  
#else
  return                  (min == -eslINFINITY || (max-min) >= 23.0f) ? max : max + h4_flogsum_lookup[(int)((max-min)*h4LOGSUM_SCALE)]; // roundf() costs too much time.
#endif
} 





/*****************************************************************
 * 2. Debugging/development tools
 *****************************************************************/

/* Function:  h4_logsum_IsSlowExact()
 * Synopsis:  Return TRUE if compiled for slow but exact calculation.
 *
 * Purpose:   When debugging, especially when comparing scores or DP
 *            cell values from a prob-space implementation against a
 *            log-space implementation, it can be useful to compile
 *            the code with <h4LOGSUM_SLOWEXACT>. This replaces the
 *            table-driven approximation with the slow but exact max +
 *            log2(1.0 + exp2(min-max)) calculation, to avoid score
 *            differences arising from the logsum approximation.  When
 *            unit tests compare scores, they need to know how
 *            stringently to do the comparison: they can call this
 *            function to decide whether to compare stringently.
 *            
 * Returns:   <TRUE> or <FALSE>.
 */
int
h4_logsum_IsSlowExact(void)
{
#ifdef h4LOGSUM_SLOWEXACT
  return TRUE;
#else
  return FALSE;
#endif
}


/* Function:  h4_logsum_InitMax()
 * Synopsis:  Initialize the lookup table for max rather than logsum
 *
 * Purpose:   Make h4_logsum(a,b) calls return max(a,b), by setting all
 *            lookup table values to 0. This converts any log-space
 *            Forward/Inside (or Backward/Outside) algorithm to a
 *            Viterbi/CYK algorithm, at no cost in speed. Useful for debugging and unit
 *            tests, where we can make sure that a complex
 *            path-summing algorithm is correctly scoring at least the
 *            single optimal path, using the native code for that
 *            algorithm.
 *            
 *            Caller must call <h4_logsum_Reinit()> when it's
 *            finished, to rebuild the lookup table correctly. For
 *            example, a unit test will call <h4_logsum_InitMax()> on
 *            entry, and <h4_logsum_Reinit()> on exit. Calling
 *            <h4_logsum_Init()> doesn't work, because of an internal
 *            flag (<logsum_initialized>) that prevents <_Init()> from
 *            wasting time rebuilding a table that's already valid.
 */
int
h4_logsum_InitMax(void)
{
  int i;
  
  for (i = 0; i < h4LOGSUM_TBL; i++)
    h4_flogsum_lookup[i] = 0.;
  h4_logsum_initialized = TRUE;
  h4_logsum_max         = TRUE;   // so we get max() even if in slow/exact mode
  return eslOK;
}
  

/* logsum_exact()
 *
 * Purpose:   Calculates $\log_2(2^a + 2^b)$ exactly, without the fast
 *            table driven approximation. Rarely used - just for some
 *            debugging and unit testing situations. To get
 *            h4_logsum() itself to use the exact calculation,
 *            compile with the <h4LOGSUM_SLOWEXACT> flag.
 */
#if defined h4LOGSUM_TESTDRIVE || defined h4LOGSUM_BENCHMARK || defined h4LOGSUM_EXAMPLE
static float 
logsum_exact(float a, float b)
{
  const float max = ESL_MAX(a, b);
  const float min = ESL_MIN(a, b);
  return (min == -eslINFINITY || (max-min) >= 23.f) ? max : max + log2f(1.0 + exp2f(min-max));  
}
#endif /* TESTDRIVE || BENCHMARK || EXAMPLE */




/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef h4LOGSUM_BENCHMARK

#include "h4_config.h"

#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "general.h"
#include "logsum.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",    0 },
  { "-n",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "naive time: log2f(exp2f(A)+exp2f(B))",    0 },
  { "-x",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "exact time: A + log2f(1+exp2f(-(A-B)))",  0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",           0 },
  { "-v",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual results",     0 },
  { "-N",        eslARG_INT,"100000000",NULL,"n>0", NULL,  NULL, NULL, "number of trials",                        0 },
  { "--version", eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 0, argc, argv, "logsum benchmark driver", "[-options]");
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  int             N       = esl_opt_GetInteger(go, "-N");
  int             i;
  float          *A, *B, *C;

  h4_logsum_Init();

  /* Create the problem: sample N values A,B on interval -10,10 */
  A = malloc(sizeof(float) * N);
  B = malloc(sizeof(float) * N);
  C = malloc(sizeof(float) * N);
  for (i = 0; i < N; i++)
    {
      A[i] = esl_random(rng) * 20. - 10.;
      B[i] = esl_random(rng) * 20. - 10.;
    }
  
  /* Run */
  esl_stopwatch_Start(w);
  if (esl_opt_GetBoolean(go, "-n"))
    {
      for (i = 0; i < N; i++)
	C[i] = log2f(exp2f(A[i]) + exp2f(B[i]));
    }
  else if (esl_opt_GetBoolean(go, "-x"))
    {
      for (i = 0; i < N; i++)
	C[i] = logsum_exact(A[i], B[i]);
    }
  else
    {
      for (i = 0; i < N; i++)
	C[i] = h4_logsum(A[i], B[i]);       
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("max sum = %.2f\n", esl_vec_FMax(C, N));

  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  free(A);
  free(B);
  free(C);
  return 0;
}
#endif /*h4LOGSUM_BENCHMARK*/
/*-------------------- end, benchmark ---------------------------*/




/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef h4LOGSUM_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

/* utest_error()
 * 
 * Analysis shows that maximum absolute error in the table-driven approximation
 * is ~ w/2, for bin width w = 1/h4LOGSUM_SCALE bits [xref logsum.md]. Test that
 * this is so, by comparing a large number of exact vs h4_logsum() calculations.
 * 
 * Also, I've designed it so that logsum(A+A) = A+1, so test that too.
 * See logsum.md for notes, and for why I don't try to reduce
 * discretization error on delta from w to w/2 by rounding delta to
 * nearest table bin, nor by prerounding the lookup table entries to
 * logsum(i+w/2). I bet some clever joker in the future will try to
 * halve the error this way; catch them.
 */
static void
utest_error(void)
{
  char            msg[] = "logsum numerical error unit test failed";
  ESL_RANDOMNESS *rng   = esl_randomness_Create(0); 
  int     N      = 1000000;
  float   maxval = 20.0;   // test operands uniformly distributed on -20..20
  float   maxerr = 0.;
  float   errlim = 1.01 * 1.0 / (2. * (float) h4LOGSUM_SCALE);   // see notes in logsum.md for approx analysis. The 1.01 allows a little slop for the approx + floating point roundoff.
  int     i;
  float   a,b,result,exact,abserr;

  for (i = 0; i < N; i++)
    {
      a = (esl_random(rng) - 0.5) * maxval * 2.; /* uniform draws on -maxval..maxval */
      b = (esl_random(rng) - 0.5) * maxval * 2.; 

      // Verify that logsum(A,B) is within expected absolute error.
      exact  = logsum_exact(a,b);
      result = h4_logsum(a,b);
      abserr = fabs(exact-result);

      if (abserr > errlim) esl_fatal(msg);
      if (abserr > maxerr) maxerr = abserr;

      // Verify that logsum(A,A) = A+1 exactly.
      result = h4_logsum(a,a);
      exact  = a + 1.;
      abserr = fabs(exact-result);
      if (abserr != 0.0) esl_fatal(msg);  // exact comparison is ok; I designed this to be true.
    }
  esl_randomness_Destroy(rng);
}


/* utest_bounds()
 * 
 * h4_logsum() needs to work when either of its arguments is -\infty.
 */
static void
utest_bounds(void)
{
  char *msg = "logsum bounds unit test failed";

  if (h4_logsum(0.0,          -eslINFINITY) !=          0.0) esl_fatal(msg);
  if (h4_logsum(-eslINFINITY,          0.0) !=          0.0) esl_fatal(msg);
  if (h4_logsum(-eslINFINITY, -eslINFINITY) != -eslINFINITY) esl_fatal(msg);
}
#endif /*p7LOGSUM_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/



/*****************************************************************
 * 6. Test driver.
 *****************************************************************/
#ifdef h4LOGSUM_TESTDRIVE
#include "h4_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "general.h"
#include "logsum.h"

static ESL_OPTIONS options[] = {
  /* name          type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",         eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",               0},
  { "--version", eslARG_NONE,    NULL, NULL, NULL,  NULL, NULL, NULL, "show HMMER version number",         0 },  
  { 0,0,0,0,0,0,0,0,0,0},
};

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go         = h4_CreateDefaultApp(options, 0, argc, argv, "logsum.c test driver", "[-options]");
  
  fprintf(stderr, "## %s\n", argv[0]);

  h4_logsum_Init();

  utest_error();
  utest_bounds();

  esl_getopts_Destroy(go);

  fprintf(stderr, "#  status = ok\n");
  return eslOK;
}
#endif /*h4LOGSUM_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/




/*****************************************************************
 * 6. Example.
 *****************************************************************/
#ifdef h4LOGSUM_EXAMPLE

#include "h4_config.h"

#include "easel.h"

#include "logsum.h"

int
main(int argc, char **argv)
{
  float a = atof(argv[1]);
  float b = atof(argv[2]);
  float r_approx, r_exact;

  h4_logsum_Init();

  r_approx = h4_logsum(a, b);
  printf("h4_logsum(%f,%f) = %f\n", a, b, r_approx);

  r_exact = logsum_exact(a,b);
  printf("log2(2^%f + 2^%f) = %f\n", a, b, r_exact);

  printf("Absolute error: %g\n", fabs(r_approx - r_exact));
  return eslOK;
}
#endif /*h4LOGSUM_EXAMPLE*/
/*--------------------- end, example ----------------------------*/



