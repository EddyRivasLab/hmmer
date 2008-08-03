/* The Forward filter implementation; SSE version.
 * 
 * SIMD vectorized, striped, interleaved, one-row O(M),
 * probability-space single precision implementation of
 * the Forward algorithm.
 * 
 * Calculates the Forward score in limited range (-127..128 bits).
 * Will not underflow for local alignment, but may overflow
 * on high-scoring target sequences.
 * 
 * Contents:
 *   1. Forward filter implementation.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 *   
 * SRE, Sun Aug  3 09:24:12 2008 [waiting out a thunderstorm, Einstein's, St. Louis]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_sse.h"


/*****************************************************************
 * 1. The p7_ForwardFilter() DP implementation.
 *****************************************************************/

/* Function:  p7_ForwardFilter()
 * Synopsis:  Calculates Forward score, vewy vewy fast, with limited upper range.
 * Incept:    SRE, Thu Dec 13 08:54:07 2007 [Janelia]
 *
 * Purpose:   Calculates the Forward score for sequence <dsq> of length <L> 
 *            residues, using optimized profile <om>, and a preallocated
 *            one-row DP matrix <ox>. Return the Forward score (in nats)
 *            in <ret_sc>.
 *            
 *            The Forward score may overflow, and will, on
 *            high-scoring sequences. Range is limited to -88 to +88
 *            nats (-127 to 127 bits). Scores will not underflow, for
 *            models configured in local mode, within HMMER's design
 *            limits ($L \leq 10^5$; $M \leq 10^4$).
 *            
 *            The model must be in a local mode; other modes cannot
 *            guarantee that we will not underflow.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: Forward score (in nats)          
 *
 * Returns:   <eslOK> on success;
 *            <eslERANGE> if the score overflows, in which case
 *            <*ret_sc> is also set <eslINFINITY>, and the sequence
 *            may be treated as a high-scoring hit.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small, or if the profile
 *            isn't in local alignment mode.
 */
int
p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m128 mpv, dpv, ipv;   /* previous row values                                       */
  register __m128 sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128 dcv;		   /* delayed storage of D(i,q+1)                               */
  register __m128 xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128 xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  __m128   zerov;		   /* splatted 0.0's in a vector                                */
  float    xN, xE, xB, xC, xJ;	   /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over quads 0..nq-1                                */
  int j;			   /* counter over DD iterations (4 is full serialization)      */
  int Q       = p7O_NQF(om->M);	   /* segment length: # of vectors                              */
  __m128 *dp  = ox->dpf[0];        /* using {MDI}MX(q) macro requires initialization of <dp>    */
  __m128 *rp;			   /* will point at om->rf[x] for residue x[i]                  */
  __m128 *tp;			   /* will point into (and step thru) om->tf                    */


  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ4)                                 ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  if (om->mode != p7_LOCAL && om->mode != p7_UNILOCAL) ESL_EXCEPTION(eslEINVAL, "Fast filter only works for local alignment");
  ox->M  = om->M;

  /* Initialization.
   */
  zerov = _mm_setzero_ps();
  for (q = 0; q < Q; q++)
    MMXo(q) = IMXo(q) = DMXo(q) = zerov;
  xE    = 0.;
  xN    = 1.;
  xJ    = 0.;
  xB    = om->xf[p7O_N][p7O_MOVE];
  xC    = 0.;

#if p7_DEBUGGING
  if (ox->debugging) p7_omx_DumpFloatRow(ox, TRUE, 0, 9, 5, xE, xN, xJ, xB, xC);	/* logify=TRUE, <rowi>=0, width=8, precision=5*/
#endif

  for (i = 1; i <= L; i++)
    {
      rp    = om->rf[dsq[i]];
      tp    = om->tf;
      dcv   = _mm_setzero_ps();
      xEv   = _mm_setzero_ps();
      xBv   = _mm_set1_ps(xB);

      /* Right shifts by 4 bytes. 4,8,12,x becomes x,4,8,12.  Shift zeros on.
       */
      mpv = MMXo(Q-1);  mpv = _mm_shuffle_ps(mpv, mpv, _MM_SHUFFLE(2, 1, 0, 0));   mpv = _mm_move_ss(mpv, zerov);
      dpv = DMXo(Q-1);  dpv = _mm_shuffle_ps(dpv, dpv, _MM_SHUFFLE(2, 1, 0, 0));   dpv = _mm_move_ss(dpv, zerov);
      ipv = IMXo(Q-1);  ipv = _mm_shuffle_ps(ipv, ipv, _MM_SHUFFLE(2, 1, 0, 0));   ipv = _mm_move_ss(ipv, zerov);
      
      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
	  sv   =                _mm_mul_ps(xBv, *tp);  tp++;
	  sv   = _mm_add_ps(sv, _mm_mul_ps(mpv, *tp)); tp++;
	  sv   = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++;
	  sv   = _mm_add_ps(sv, _mm_mul_ps(dpv, *tp)); tp++;
	  sv   = _mm_mul_ps(sv, *rp);                  rp++;
	  xEv  = _mm_add_ps(xEv, sv);
	  
	  /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
	   * {MDI}MX(q) is then the current, not the prev row
	   */
	  mpv = MMXo(q);
	  dpv = DMXo(q);
	  ipv = IMXo(q);

	  /* Do the delayed stores of {MD}(i,q) now that memory is usable */
	  MMXo(q) = sv;
	  DMXo(q) = dcv;

	  /* Calculate the next D(i,q+1) partially: M->D only;
           * delay storage, holding it in dcv
	   */
	  dcv   = _mm_mul_ps(sv, *tp); tp++;

	  /* Calculate and store I(i,q) */
	  sv     =                _mm_mul_ps(mpv, *tp);  tp++;
	  sv     = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++;
	  IMXo(q) = _mm_mul_ps(sv, *rp);                  rp++;
	}	  

      /* Now the DD paths. We would rather not serialize them but 
       * in an accurate Forward calculation, we have few options.
       */
      /* dcv has carried through from end of q loop above; store it 
       * in first pass, we add M->D and D->D path into DMX
       */
      /* We're almost certainly're obligated to do at least one complete 
       * DD path to be sure: 
       */
      dcv    = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
      dcv    = _mm_move_ss(dcv, zerov);
      DMXo(0) = zerov;
      tp     = om->tf + 7*Q;	/* set tp to start of the DD's */
      for (q = 0; q < Q; q++) 
	{
	  DMXo(q) = _mm_add_ps(dcv, DMXo(q));	
	  dcv    = _mm_mul_ps(DMXo(q), *tp); tp++; /* extend DMXo(q), so we include M->D and D->D paths */
	}

      /* now. on small models, it seems best (empirically) to just go
       * ahead and serialize. on large models, we can do a bit better,
       * by testing for when dcv (DD path) accrued to DMXo(q) is below
       * machine epsilon for all q, in which case we know DMXo(q) are all
       * at their final values. The tradeoff point is (empirically) somewhere around M=100,
       * at least on my desktop. We don't worry about the conditional here;
       * it's outside any inner loops.
       */
      if (om->M < 100)
	{			/* Fully serialized version */
	  for (j = 1; j < 4; j++)
	    {
	      dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
	      dcv = _mm_move_ss(dcv, zerov);
	      tp = om->tf + 7*Q;	/* set tp to start of the DD's */
	      for (q = 0; q < Q; q++) 
		{
		  DMXo(q) = _mm_add_ps(dcv, DMXo(q));	
		  dcv    = _mm_mul_ps(dcv, *tp);   tp++; /* note, extend dcv, not DMXo(q); only adding DD paths now */
		}	    
	    }
	} 
      else
	{			/* Slightly parallelized version, but which incurs some overhead */
	  for (j = 1; j < 4; j++)
	    {
	      register __m128 cv;	/* keeps track of whether any DD's change DMXo(q) */

	      dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(2, 1, 0, 0));
	      dcv = _mm_move_ss(dcv, zerov);
	      tp  = om->tf + 7*Q;	/* set tp to start of the DD's */
	      cv  = zerov;
	      for (q = 0; q < Q; q++) 
		{
		  sv     = _mm_add_ps(dcv, DMXo(q));	
		  cv     = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMXo(q))); /* remember if DD paths changed any DMXo(q): *without* conditional branch */
		  DMXo(q) = sv;	                                    /* store new DMXo(q) */
		  dcv    = _mm_mul_ps(dcv, *tp);   tp++;            /* note, extend dcv, not DMXo(q); only adding DD paths now */
		}	    
	      if (! _mm_movemask_ps(cv)) break; /* DD's didn't change any DMXo(q)? Then we're done, break out. */
	    }
	}

      /* Add D's to xEv */
      for (q = 0; q < Q; q++) xEv = _mm_add_ps(DMXo(q), xEv);

      /* Finally the "special" states, which start from Mk->E (->C, ->J->B) */
      /* The following incantation is a horizontal sum of xEv's elements  */
      /* These must follow DD calculations, because D's contribute to E in Forward
       * (as opposed to Viterbi)
       */
      xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(0, 3, 2, 1)));
      xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(1, 0, 3, 2)));
      _mm_store_ss(&xE, xEv);

      xN =  xN * om->xf[p7O_N][p7O_LOOP];
      xC = (xC * om->xf[p7O_C][p7O_LOOP]) +  (xE * om->xf[p7O_E][p7O_MOVE]);
      xJ = (xJ * om->xf[p7O_J][p7O_LOOP]) +  (xE * om->xf[p7O_E][p7O_LOOP]);
      xB = (xJ * om->xf[p7O_J][p7O_MOVE]) +  (xN * om->xf[p7O_N][p7O_MOVE]);
      /* and now xB will carry over into next i, and xC carries over after i=L */

#if p7_DEBUGGING
      if (ox->debugging) p7_omx_DumpFloatRow(ox, TRUE, i, 9, 5, xE, xN, xJ, xB, xC);	/* logify=TRUE, <rowi>=i, width=8, precision=5*/
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T, and flip back to log space (nats) */
  /* On overflow, xC is inf or nan (nan arises because inf*0 = nan). */
  /* On an underflow (which shouldn't happen), we counterintuitively return infinity:
   * the effect of this is to force the caller to rescore us with full range.
   */
  if       (isnan(xC))      { *ret_sc = eslINFINITY; return eslERANGE; }
  else if  (xC == 0.0)      { *ret_sc = eslINFINITY; return eslERANGE; } /* on underflow, force caller to rescore us! */
  else if  (isinf(xC) == 1) { *ret_sc = eslINFINITY; return eslERANGE; }
  else                      *ret_sc = log(xC * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}
/*------------------ end, p7_ForwardFilter() --------------------*/


/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7FWDFILTER_BENCHMARK
/* -c and -x options are useful for debugging and characterization
 * of precision and range. For all the optimized implementations,
 * -c compares to the "generic" implementation of the same scoring
 * algorithm (MSV, Viterbi, Forward). -x compares to a well-trusted
 * implementation that achieves exactly the same scores, possibly
 * by some sort of emulation mode.
 * 
 * Here, for ForwardFilter, -c compares to GForward()
 * scores. p7_ForwardFilter() calculates "exact" scores in probability
 * space, whereas p7_GForward() uses a table-driven log-sum
 * approximation; comparison of the two scores reveals the imprecision
 * in p7_GForward()'s approximation.
 * 
 * -x also compares to GForward() scores, but requires that you first
 * go over to logsum.c::p7_FLogsum(), uncomment to have it calculate
 * log(exp(s1) + exp(s2) (slow but exact), and recompile the driver.
 * (If you don't do this, the driver will exit with a warning.)
 */

/* 
   gcc -o benchmark-fwdfilter -std=gnu99 -g -Wall -msse2 -I.. -L.. -I../../easel -L../../easel -Dp7FWDFILTER_BENCHMARK fwdfilter.c -lhmmer -leasel -lm 
   icc -o benchmark-fwdfilter -O3 -static -I.. -L.. -I../../easel -L../../easel -Dp7FWDFILTER_BENCHMARK fwdfilter.c -lhmmer -leasel -lm 

   ./benchmark-fwdfilter <hmmfile>             runs benchmark
   ./benchmark-fwdfilter -b <hmmfile>          gets baseline time to subtract: just random seq generation
   ./benchmark-fwdfilter -c -N100 <hmmfile>    compare scores of SSE to generic impl
   ./benchmark-fwdfilter -x -N100 <hmmfile>    test that scores match trusted implementation.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "baseline timing: don't run DP at all",             0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-x", "compare scores to generic implementation (debug)", 0 }, 
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                  0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-c", "equate scores to trusted implementation (debug)",  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for SSE ForwardFilter implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox      = NULL;
  P7_GMX         *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc1, sc2;

  p7_FLogsumInit();

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  if (esl_opt_GetBoolean(go, "-x") && p7_FLogsumError(-0.4, -0.5) > 0.0001)
    p7_Fail("-x here requires p7_Logsum() recompiled in slow exact mode");

  ox = p7_omx_Create(gm->M, 0, 0); /* a one-row matrix */
  gx = p7_gmx_Create(gm->M, L);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      if (! esl_opt_GetBoolean(go, "-b")) 
	{
	  p7_ForwardFilter(dsq, L, om, ox, &sc1);   

	  if (esl_opt_GetBoolean(go, "-c") || esl_opt_GetBoolean(go, "-x"))
	  {
	    p7_GForward(dsq, L, gm, gx, &sc2); 
	    printf("%.4f %.4f\n", sc1, sc2);  
	  }
	}
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);

  free(dsq);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7FWDFILTER_BENCHMARK*/
/*--------------------- end, benchmark driver -------------------*/




/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7FWDFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* ForwardFilter() unit test:
 * compare to GForward() scores and ForwardParser() scores.
 */
static void
utest_forward_filter(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M, 0, L);
  P7_GMX      *gx  = p7_gmx_Create(M, L);
  float tolerance;
  float sc1, sc2, sc3;

  if (p7_FLogsumError(-0.4, -0.5) > 0.0001) tolerance = 1.0;  /* weaker test.   */
  else tolerance = 0.0001;   /* stronger test: FLogsum() is in slow exact mode. */

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_ForwardFilter(dsq, L, om, ox, &sc1);
      p7_GForward     (dsq, L, gm, gx, &sc2);
      p7_ForwardParser(dsq, L, om, ox, &sc3);

      if (fabs(sc1-sc2) > tolerance)
	esl_fatal("ForwardFilter unit test failed: GForward differs (%.2f, %.2f)", sc1, sc2);
      if (fabs(sc1-sc3) > 0.0001)
	esl_fatal("ForwardFilter unit test failed: ForwardParser differs (%.2f, %.2f)", sc1, sc3);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7FWDFILTER_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/

/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7FWDFILTER_TESTDRIVE
/* 
   gcc -g -Wall -msse2 -std=gnu99 -I.. -L.. -I../../easel -L../../easel -o fwdfilter_utest -Dp7FWDFILTER_TESTDRIVE fwdfilter.c -lhmmer -leasel -lm
   ./fwdfilter_utest
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the SSE implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = NULL;
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  p7_FLogsumInit();

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("ForwardFilter() tests, DNA\n");
  utest_forward_filter(r, abc, bg, M, L, N);   
  utest_forward_filter(r, abc, bg, 1, L, 10);  
  utest_forward_filter(r, abc, bg, M, 1, 10);  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("ForwardFilter() tests, protein\n");
  utest_forward_filter(r, abc, bg, M, L, N);   
  utest_forward_filter(r, abc, bg, 1, L, 10);  
  utest_forward_filter(r, abc, bg, M, 1, 10);  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7FWDFILTER_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/

/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7FWDFILTER_EXAMPLE
/* A minimal example.
   Also useful for debugging on small HMMs and sequences.

   gcc -g -Wall -msse2 -std=gnu99 -I.. -L.. -I../../easel -L../../easel -o example -Dp7FWDFILTER_EXAMPLE fwdfilter.c -lhmmer -leasel -lm
   ./example <hmmfile> <seqfile>
 */ 
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "impl_sse.h"

int 
main(int argc, char **argv)
{
  char           *hmmfile = argv[1];
  char           *seqfile = argv[2];
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox      = NULL;
  P7_GMX         *gx      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           sc;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
  if  (esl_sqio_Read(sqfp, sq) != eslOK) p7_Fail("Failed to read sequence");

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  /* allocate DP matrices, both a generic and an optimized one */
  ox = p7_omx_Create(gm->M, 0, sq->n);
  gx = p7_gmx_Create(gm->M, sq->n);

  /* Useful to place and compile in for debugging: 
     p7_oprofile_Dump(stdout, om);      dumps the optimized profile
     p7_omx_SetDumpMode(ox, TRUE);      makes the fast DP algorithms dump their matrices
     p7_gmx_Dump(stdout, gx);           dumps a generic DP matrix
  */

  p7_ForwardFilter(sq->dsq, sq->n, om, ox, &sc);  printf("forward filter score: %.2f nats\n", sc);
  p7_GForward     (sq->dsq, sq->n, gm, gx, &sc);  printf("forward (generic):    %.2f nats\n", sc);
  p7_ForwardParser(sq->dsq, sq->n, om, ox, &sc);  printf("forward parser:       %.2f nats\n", sc);

  /* now in a real app, you'd need to convert raw nat scores to final bit
   * scores, by subtracting the null model score and rescaling.
   */

  /* cleanup */
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  return 0;
}
#endif /*p7FWDFILTER_EXAMPLE*/
/*-------------------- end, example -----------------------------*/




/*****************************************************************
 * @LICENSE@
 *****************************************************************/
