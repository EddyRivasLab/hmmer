/* Accuracy metrics for state paths (alignments)
 * 
 * Contents: 
 *    1. The <P7_TRACE_METRICS> object
 *    2. p7_trace_metrics() calculation
 *    3. Unit tests
 *    4. Test driver
 *    5. Example
 *    6. Copyright and license information
 */

#include "p7_config.h"

#include "easel.h"

#include "base/p7_trace.h"

#include "misc/p7_trace_metrics.h"

/*****************************************************************
 * 1. The P7_TRACE_METRICS object
 *****************************************************************/

/* Function:  p7_trace_metrics_Create()
 * Synopsis:  Create a new trace metrics object.
 *
 * Purpose:   Create a new trace metrics object, with
 *            its counts all initialized to zero; return
 *            a pointer to it. 
 *
 * Args:      (void)
 *
 * Returns:   ptr to new <P7_TRACE_METRICS> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_TRACE_METRICS *
p7_trace_metrics_Create(void)
{
  P7_TRACE_METRICS *tm = NULL;
  int               status;

  ESL_ALLOC(tm, sizeof(P7_TRACE_METRICS));
  p7_trace_metrics_Zero(tm);
  return tm;
  
 ERROR:
  return NULL;
}

/* Function:  p7_trace_metrics_Zero()
 * Synopsis:  Zero all the counts in a trace metrics object.
 *
 * Purpose:   Zero the counts in a trace metrics object <tm>.
 */
int
p7_trace_metrics_Zero(P7_TRACE_METRICS *tm)
{
  tm->state_tp  = 0;
  tm->state_fp  = 0;
  tm->state_fn  = 0;
  tm->align_tp  = 0;
  tm->align_fp  = 0;
  tm->align_fn  = 0;
  tm->region_tp = 0;
  tm->region_fp = 0;
  tm->region_fn = 0;
  tm->edge_tp   = 0;
  tm->edge_fp   = 0;
  tm->edge_fn   = 0;
  
  return eslOK;
}

/* Function:  p7_trace_metrics_Destroy()
 * Synopsis:  Frees a trace metrics object.
 */
void
p7_trace_metrics_Destroy(P7_TRACE_METRICS *tm)
{
  if (tm) free(tm);
}
/*----------- end, P7_TRACE_METRICS object ----------------------*/


/*****************************************************************
 * 2. p7_trace_metrics() calculation
 *****************************************************************/

/* Function:  p7_trace_metrics()
 * Synopsis:  Calculate and accumulate tp/fp/fn state path accuracy counts. 
 *
 * Purpose:   Given a reference trace <reftr> that we assume is a gold
 *            standard, and a test trace <testtr> to compare to that
 *            reference; calculate accuracy metrics and increment 
 *            corresponding counters in the metrics structure <tm>.
 *            
 *            This is an increment (accumulation of counts) so that
 *            caller may sum up average accuracy statistics over many
 *            state path comparisons.
 *            
 *            Four measures are calculated: state, alignment, region,
 *            and edge accuracies. 
 *            
 *            For state accuracy, every state in the state paths are
 *            compared; matches (true positives) agree completely, in
 *            state type, position i in the target sequence, and (for
 *            main model states) position k in the query profile.
 *            Three traceback positions always match, by construction:
 *            <S->N> and the terminal <T>. False positives are in the
 *            test path but not the reference; false negatives are in
 *            the reference but not the test. The sum <TP+FN> is the
 *            number of states in the reference path <reftr->N>;
 *            <TP+FP> is <testtr->N>.
 *            
 *            For alignment accuracy, only main model residue-emitting
 *            states are counted, and we evaluate the exact
 *            correspondence of Mk/Ik states to residues i. To count
 *            as a match (true positive), i, k, and statetype must all
 *            agree (except that Glocal/local distinction is ignored;
 *            alignment to ML versus MG is counted as identical, for
 *            example). The sum <TP+FN> is the number of M/I states in
 *            the reference path; <TP+FP> is the number in the test
 *            path.
 *            
 *            For region accuracy, again only main model
 *            residue-emitting states are counted, but we only
 *            evaluate whether the state path has correctly assigned
 *            residues i as homologous (in the model) or not; the
 *            exact correspondence to an Mk/Ik state is not evaluated.
 *            To count as a match, both the reference and test trace
 *            must have assigned residue i to a main model state
 *            (M/I). False positives are residues assigned as
 *            homologous only in the test path; false negatives are
 *            only in the reference path. The sum <TP+FN> is the
 *            number of M/I states in the reference path, and <TP+FP>
 *            is the number in the test path.
 *            
 *            For edge accuracy, only B and E states are counted.  To
 *            match, there must be a corresponding B/E state at the
 *            same target sequence position <i> in both reference and
 *            test state path. False positives are edges called only
 *            in the test path; false negatives are edges only in the
 *            reference path. The sum <TP+FN> is the number of B/E
 *            states in the reference path; <TP+FP> is the number in
 *            the test path.
 *
 * Args:      reftr  - reference state path
 *            testtr - test state path
 *            tm     - RESULT: accumulated accuracy metrics
 *
 * Returns:   <eslOK> on success, and <tm> has been updated.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_trace_metrics(const P7_TRACE *reftr, const P7_TRACE *testtr, P7_TRACE_METRICS *tm)
{
  int zr, zt;			/* position in reference, test trace */

  for (zr = zt = 0; zr < reftr->N; zr++)
    {
      /* Increment zt, trying to align test to the current position of the reference    */
      if (reftr->i[zr]) {  /* if reference emits i: move test to where it emitted i too */
	while (testtr->i[zt] < reftr->i[zr])
	  { /* everything we must skip over in the test is a) nonemitting and b) a false positive */
	    if (testtr->st[zt] == p7T_B || testtr->st[zt] == p7T_E) tm->edge_fp++;
	    tm->state_fp++;
	    zt++;
	  }
      } else { /* if reference is mute: either find a corresponding mute state in test, or halt and wait on T or next emitter */		
	while (testtr->i[zt] == 0 && ! (testtr->st[zt] == reftr->st[zr] && testtr->k[zt] == reftr->k[zr]) && zt < testtr->N-1)
	  {  /* again, anything extra we see in the test trace is a false positive */
	    if (testtr->st[zt] == p7T_B || testtr->st[zt] == p7T_E) tm->edge_fp++;
	    tm->state_fp++;
	    zt++;
	  }
      }

      /* After the zt update, we can be in one of the following states:
       *   if reference trace emits: then test trace also emits
       *     same position i. (though not necessarily w/ same state or k)
       *   if reference trace is on a mute state: test trace may emit or not.
       *      if test trace emits: then its position i is the next one that
       *          the reference trace will reach next; it is paused,
       *          waiting for reference to get thru mute states and catch up.
       *      if not, then either: reference and test are on same mute state;
       *          or the test trace is already done, on its terminal T state.
       *
       *   corollary: if test trace is on a mute state, it must either match 
       *      the reference mute state, or it must be done (on T).
       */
      /* All-states comparison */
      if (reftr->st[zr] == testtr->st[zt] && reftr->i[zr]  == testtr->i[zt]  && reftr->k[zr]  == testtr->k[zt]) 
	  tm->state_tp++;
      else {
	tm->state_fn++;
	if (reftr->i[zr])  tm->state_fp++; /* checking i[z] makes sure we only count each i once */
      }

      /* Alignment comparison */
      if (reftr->i[zr])	/* this makes sure we only evaluate alignment of each i once. we can go through zr loop more than once w/ testtr->i[zt] on the same i */
	{
	  if ( p7_trace_IsMain(reftr->st[zr]))
	    {
	      if ( reftr->k[zr] == testtr->k[zt] && 
		   ( (p7_trace_IsM(reftr->st[zr]) && p7_trace_IsM(testtr->st[zt])) ||
		     (p7_trace_IsI(reftr->st[zr]) && p7_trace_IsI(testtr->st[zt]))) )
		tm->align_tp++;
	      else {
		tm->align_fn++;	
		if (p7_trace_IsMain(testtr->st[zt])) tm->align_fp++; /* not only did testtr fail to correctly align i, it aligned it somewhere else that's wrong (a false pos, in addition to being a false neg) */
	      }
	    }
	  else { /* reference emits it w/ NJC; if test puts it in model, it's a fp */
	    if (p7_trace_IsMain(testtr->st[zt])) tm->align_fp++;
	  }
	}
      
      /* Region comparison */
      if (reftr->i[zr]) {
	if      (p7_trace_IsMain(reftr->st[zr]) && p7_trace_IsMain(testtr->st[zt])) tm->region_tp++;
	else if (p7_trace_IsMain(reftr->st[zr]))                    	            tm->region_fn++;
	else if (p7_trace_IsMain(testtr->st[zt]))	                            tm->region_fp++;
      }

      /* Edge comparison */
      if (reftr->st[zr] == p7T_B) {
	if (testtr->st[zt] == p7T_B) tm->edge_tp++;
	else                         tm->edge_fn++;
      }
      if (reftr->st[zr] == p7T_E) {
	if (testtr->st[zt] == p7T_E) tm->edge_tp++;
	else                         tm->edge_fn++;
      }
      
      /* Advance zt unless it's supposed to wait */
      if (reftr->i[zr] || ( zt < testtr->N-1 && testtr->i[zt] == 0)) zt++;
    }
  return eslOK;
}
/*------------ end, p7_trace_metrics() calculation --------------*/


/*****************************************************************
 * 2. Debugging/development tools
 *****************************************************************/

int
p7_trace_metrics_Dump(FILE *ofp, P7_TRACE_METRICS *tm)
{
  fprintf(ofp, "metric      TP      FP      FN    sens     ppv\n");
  fprintf(ofp, "-------  ------- ------- ------- ------- -------\n");

  fprintf(ofp, "state    %7d %7d %7d %7.4f %7.4f\n", 
	  tm->state_tp, tm->state_fp, tm->state_fn,
	  (double) tm->state_tp / ((double) tm->state_tp + (double) tm->state_fn),
	  (double) tm->state_tp / ((double) tm->state_tp + (double) tm->state_fp));
  fprintf(ofp, "align    %7d %7d %7d %7.4f %7.4f\n", 
	  tm->align_tp, tm->align_fp, tm->align_fn,
	  (double) tm->align_tp / ((double) tm->align_tp + (double) tm->align_fn),
	  (double) tm->align_tp / ((double) tm->align_tp + (double) tm->align_fp));
  fprintf(ofp, "region   %7d %7d %7d %7.4f %7.4f\n", 
	  tm->region_tp, tm->region_fp, tm->region_fn,
	  (double) tm->region_tp / ((double) tm->region_tp + (double) tm->region_fn),
	  (double) tm->region_tp / ((double) tm->region_tp + (double) tm->region_fp));
  fprintf(ofp, "edge     %7d %7d %7d %7.4f %7.4f\n", 
	  tm->edge_tp, tm->edge_fp, tm->edge_fn,
	  (double) tm->edge_tp / ((double) tm->edge_tp + (double) tm->edge_fn),
	  (double) tm->edge_tp / ((double) tm->edge_tp + (double) tm->edge_fp));
  return eslOK;
}


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7TRACE_METRICS_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "base/p7_bg.h"
#include "build/modelsample.h"
#include "search/modelconfig.h"
#include "misc/emit.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_viterbi.h"

/* The "consistency" test samples <N> random profiles; from each it
 * emits a sequence in multihit glocal/local <L> mode, and keeps the
 * emitted trace as a reference; it then realigns the emitted sequence
 * by ReferenceViterbi and takes the Viterbi trace as a test path.
 * 
 * For each reference/test path pair, it runs p7_trace_metrics() and
 * validates that the expected TP+FP, TP+FN sums add up to what they
 * should (see header documentation on p7_trace_metrics()).
 */
static void
utest_consistency(ESL_RANDOMNESS *rng, int alphatype, int M, int L, int N)
{
  char              msg[]  = "trace_metrics consistency unit test failed";
  ESL_ALPHABET     *abc    = esl_alphabet_Create(alphatype);
  ESL_SQ           *sq     = esl_sq_CreateDigital(abc);
  P7_BG            *bg     = p7_bg_Create(abc);
  P7_HMM           *hmm    = NULL;
  P7_PROFILE       *gm     = p7_profile_Create(M, abc);
  P7_TRACE         *reftr  = p7_trace_Create();
  P7_TRACE         *testtr = p7_trace_Create();
  P7_REFMX         *vmx    = p7_refmx_Create(100, 100); /* will grow as needed */
  P7_TRACE_METRICS *tm     = p7_trace_metrics_Create();
  int               refct[p7T_NSTATETYPES];
  int               testct[p7T_NSTATETYPES];
  int               idx;
  int               test_nMI, ref_nMI;
  int               test_nBE, ref_nBE;

  for (idx = 0; idx < N; idx++)
    {
      /* Zero the metrics: we're going to test one at a time. */
      if ( p7_trace_metrics_Zero(tm) != eslOK) esl_fatal(msg);

      /* Sample a random HMM, and config a profile (multihit local/glocal) from it */
      if ( p7_hmm_Sample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
      if ( p7_profile_Config(gm, hmm, bg)   != eslOK) esl_fatal(msg);
      if ( p7_profile_SetLength(gm, L)      != eslOK) esl_fatal(msg);

      /* Emit a sequence and keep its trace as the reference */
      if ( p7_profile_SetLength(gm, L)                 != eslOK) esl_fatal(msg);
      if ( p7_ProfileEmit(rng, hmm, gm, bg, sq, reftr) != eslOK) esl_fatal(msg);

      /* Viterbi alignment of that sequence back to the model gives us a test state path */
      if ( p7_profile_SetLength(gm, sq->n)                            != eslOK) esl_fatal(msg);
      if ( p7_bg_SetLength(bg, sq->n)                                 != eslOK) esl_fatal(msg); 
      if ( p7_ReferenceViterbi(sq->dsq, sq->n, gm, vmx, testtr, NULL) != eslOK) esl_fatal(msg);

      /* calculate the metrics... */
      if ( p7_trace_metrics(reftr, testtr, tm) != eslOK) esl_fatal(msg);

      //p7_trace_DumpAnnotated(stdout, reftr, gm, sq->dsq);
      //p7_trace_DumpAnnotated(stdout, testtr, gm, sq->dsq);
      //p7_trace_metrics_Dump(stdout, tm);

      /* independently, tote up state use counts, used for consistency checks */
      if ( p7_trace_GetStateUseCounts(reftr, refct)   != eslOK) esl_fatal(msg);
      if ( p7_trace_GetStateUseCounts(testtr, testct) != eslOK) esl_fatal(msg);
      test_nMI = testct[p7T_MG] + testct[p7T_ML] + testct[p7T_IG] + testct[p7T_IL];
      ref_nMI  =  refct[p7T_MG] +  refct[p7T_ML] +  refct[p7T_IG] +  refct[p7T_IL];
      test_nBE = testct[p7T_B]  + testct[p7T_E];
      ref_nBE  =  refct[p7T_B]  +  refct[p7T_E];

      /* consistency checks as documented in p7_trace_metrics() header */
      if (tm->state_tp  + tm->state_fn  !=  reftr->N) esl_fatal(msg); 
      if (tm->state_tp  + tm->state_fp  != testtr->N) esl_fatal(msg); 
      if (tm->align_tp  + tm->align_fn  !=  ref_nMI)  esl_fatal(msg);
      if (tm->align_tp  + tm->align_fp  != test_nMI)  esl_fatal(msg);
      if (tm->region_tp + tm->region_fn !=  ref_nMI)  esl_fatal(msg);
      if (tm->region_tp + tm->region_fp != test_nMI)  esl_fatal(msg);
      if (tm->edge_tp   + tm->edge_fn   !=  ref_nBE)  esl_fatal(msg);
      if (tm->edge_tp   + tm->edge_fp   != test_nBE)  esl_fatal(msg);

      /* additionally, both traces should Validate()... of course. */
      if (p7_trace_Validate( reftr, abc, sq->dsq, NULL) != eslOK) esl_fatal(msg);
      if (p7_trace_Validate(testtr, abc, sq->dsq, NULL) != eslOK) esl_fatal(msg);

      p7_refmx_Reuse(vmx);
      p7_trace_Reuse(reftr);
      p7_trace_Reuse(testtr);
      esl_sq_Reuse(sq);
    }  

  p7_trace_Destroy(reftr);
  p7_trace_Destroy(testtr);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
}

#endif /*p7TRACE_METRICS_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/


/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
/* gcc -o p7_trace_metrics_utest -Dp7TRACE_METRICS_TESTDRIVE -I . -I ../easel -L. -L../easel p7_trace_metrics.c -lhmmer -leasel -lm 
 */
#ifdef p7TRACE_METRICS_TESTDRIVE


#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 }, /* this test should always succeed; is safe to use any rng seed */
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for trace metrics";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go        = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r         = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             alphatype = eslAMINO;
  int             M         = 20;
  int             L         = 20;
  int             N         = 10;

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  utest_consistency(r, alphatype, M, L, N);

  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  fprintf(stderr, "#  status = ok\n");
  return 0;
}
#endif /*p7TRACE_METRICS_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/



/*****************************************************************
 * 5. Example.
 *****************************************************************/
#ifdef p7TRACE_METRICS_EXAMPLE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "configured mean seq length for profile",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "example of calculating accuracy metrics for comparison of two traces";


/* Test: emit two traces from the same profile. 
 *       compare them.
 */
int
main(int argc, char **argv)
{
  ESL_GETOPTS      *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_RANDOMNESS   *rng     = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  char             *hmmfile = esl_opt_GetArg(go, 1);
  int               L       = esl_opt_GetInteger(go, "-L");
  ESL_ALPHABET     *abc     = NULL;
  P7_HMMFILE       *hfp     = NULL;
  P7_HMM           *hmm     = NULL;
  P7_BG            *bg      = NULL;
  P7_PROFILE       *gm      = NULL;
  P7_TRACE         *reftr   = p7_trace_Create();
  P7_REFMX         *vmx     = NULL;
  P7_TRACE         *testtr  = p7_trace_Create();
  ESL_SQ           *sq      = NULL;
  P7_TRACE_METRICS *tm      = p7_trace_metrics_Create();
  char              errbuf[eslERRBUFSIZE];
  int               status;

  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",                       status, hmmfile, errbuf);  

  status = p7_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT)   p7_Fail("Bad file format in HMM file %s:\n%s\n",          hfp->fname, hfp->errbuf);
  else if (status == eslEINCOMPAT) p7_Fail("HMM in %s is not in the expected %s alphabet\n", hfp->fname, esl_abc_DecodeType(abc->type));
  else if (status == eslEOF)       p7_Fail("Empty HMM file %s? No HMM data found.\n",        hfp->fname);
  else if (status != eslOK)        p7_Fail("Unexpected error in reading HMMs from %s\n",     hfp->fname);
  p7_hmmfile_Close(hfp);

  bg = p7_bg_Create(abc);                
  gm = p7_profile_Create(hmm->M, abc);   

  p7_profile_ConfigUniglocal(gm, hmm, bg, L); 
  sq = esl_sq_CreateDigital(abc);

  /* Emit a sequence, and keep its trace as the "reference" */
  p7_ProfileEmit(rng, hmm, gm, bg, sq, reftr);

  /* Viterbi alignment of that sequence back to the model: Viterbi trace is the "test" */
  vmx = p7_refmx_Create(gm->M, sq->n);

  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, sq->n);
  p7_bg_SetLength(bg, sq->n);
  p7_ReferenceViterbi(sq->dsq, sq->n, gm, vmx, testtr, NULL);

  p7_trace_DumpAnnotated(stdout, reftr,  gm, sq->dsq);
  p7_trace_DumpAnnotated(stdout, testtr, gm, sq->dsq);

  /* Calculate trace accuracy metrics. */
  p7_trace_metrics(reftr, testtr, tm);

  p7_trace_metrics_Dump(stdout, tm);

  p7_trace_metrics_Destroy(tm);
  p7_refmx_Destroy(vmx);
  p7_trace_Destroy(reftr);
  p7_trace_Destroy(testtr);
  esl_sq_Destroy(sq);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7TRACE_METRICS_EXAMPLE*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
