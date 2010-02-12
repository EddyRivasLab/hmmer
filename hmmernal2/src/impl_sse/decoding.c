/* Posterior decoding algorithms; SSE versions.
 * 
 * Contents:
 *   1. Posterior decoding algorithms.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 *   
 * SRE, Mon Aug 18 08:15:50 2008 [Janelia]
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
 * 1. Posterior decoding algorithms.
 *****************************************************************/

/* Function:  p7_Decoding()
 * Synopsis:  Posterior decoding of residue assignment.
 * Incept:    SRE, Fri Aug  8 14:29:42 2008 [UA217 to SFO]
 *
 * Purpose:   Identical to <p7_GDecoding()> except that <om>, <oxf>,
 *            <oxb> are SSE optimized versions. See <p7_GDecoding()>
 *            documentation for more info.
 *
 * Args:      om   - profile (must be the same that was used to fill <oxf>, <oxb>).
 *            oxf  - filled Forward matrix 
 *            oxb  - filled Backward matrix
 *            pp   - RESULT: posterior decoding matrix.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_Decoding(const P7_OPROFILE *om, const P7_OMX *oxf, P7_OMX *oxb, P7_OMX *pp)
{
  register __m128 *ppv;
  register __m128 *fv;
  register __m128 *bv;
  register __m128  totrv;
  register __m128  scv;
  int    L  = oxf->L;
  int    M  = om->M;
  int    Q  = p7O_NQF(M);	
  int    i,q;
  float  totp_reciprocal = 1.0 / (oxf->xmx[p7X_C][L] * om->xf[p7O_C][p7O_MOVE]);

  ppv = pp->dpf[0];
  for (q = 0; q < Q; q++) {
    *ppv = _mm_setzero_ps(); ppv++;
    *ppv = _mm_setzero_ps(); ppv++;
    *ppv = _mm_setzero_ps(); ppv++;
  }
  pp->xmx[p7X_E][0] = 0.0;
  pp->xmx[p7X_N][0] = 0.0;
  pp->xmx[p7X_J][0] = 0.0;
  pp->xmx[p7X_C][0] = 0.0;
  pp->xmx[p7X_B][0] = 0.0;

  totrv = _mm_set1_ps(totp_reciprocal);
  for (i = 1; i <= L; i++)
    {
      ppv =  pp->dpf[i];
      fv  = oxf->dpf[i];
      bv  = oxb->dpf[i];
      scv = _mm_set1_ps(oxf->xmx[p7X_SCALE][i]);

      for (q = 0; q < Q; q++)
	{
	  /* M */
	  *ppv = _mm_mul_ps(*fv,  *bv);
	  *ppv = _mm_mul_ps(*ppv,  scv);
	  *ppv = _mm_mul_ps(*ppv,  totrv);
	  ppv++;  fv++;  bv++;

	  /* D */
	  *ppv = _mm_setzero_ps();
	  ppv++;  fv++;  bv++;

	  /* I */
	  *ppv = _mm_mul_ps(*fv,  *bv);
	  *ppv = _mm_mul_ps(*ppv,  scv);
	  *ppv = _mm_mul_ps(*ppv,  totrv);
	  ppv++;  fv++;  bv++;
	}
      pp->xmx[p7X_E][i] = 0.0;
      pp->xmx[p7X_N][i] = oxf->xmx[p7X_N][i-1] * oxb->xmx[p7X_N][i] * om->xf[p7O_N][p7O_LOOP] * totp_reciprocal;
      pp->xmx[p7X_J][i] = oxf->xmx[p7X_J][i-1] * oxb->xmx[p7X_J][i] * om->xf[p7O_J][p7O_LOOP] * totp_reciprocal;
      pp->xmx[p7X_C][i] = oxf->xmx[p7X_C][i-1] * oxb->xmx[p7X_C][i] * om->xf[p7O_C][p7O_LOOP] * totp_reciprocal;
      pp->xmx[p7X_B][i] = 0.0;
    }
  return eslOK;
}

/* Function:  p7_DomainDecoding()
 * Synopsis:  Posterior decoding of domain location.
 * Incept:    SRE, Tue Aug  5 08:39:07 2008 [Janelia]
 *
 * Purpose:   Identical to <p7_GDomainDecoding()> except that <om>, <oxf>,
 *            <oxb> are SSE optimized versions. See <p7_GDomainDecoding()>
 *            documentation for more info.
 *
 * Args:      gm   - profile
 *            oxf  - filled Forward matrix
 *            oxb  - filled Backward matrix
 *            ddef - container for the results.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_DomainDecoding(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef)
{
  int   L               = oxf->L;
  float totp_reciprocal = 1.0 / (oxf->xmx[p7X_C][L] * om->xf[p7O_C][p7O_MOVE]);
  float njcp;
  int   i;

  ddef->btot[0] = 0.0;
  ddef->etot[0] = 0.0;
  ddef->mocc[0] = 0.0;
  for (i = 1; i <= L; i++)
    {
      ddef->btot[i] = ddef->btot[i-1] +
	(oxf->xmx[p7X_B][i-1] * oxb->xmx[p7X_B][i-1] * oxf->xmx[p7X_SCALE][i-1] * totp_reciprocal);
      ddef->etot[i] = ddef->etot[i-1] +
	(oxf->xmx[p7X_E][i]   * oxb->xmx[p7X_E][i]   * oxf->xmx[p7X_SCALE][i]   * totp_reciprocal);

      njcp  = oxf->xmx[p7X_N][i-1] * oxb->xmx[p7X_N][i] * om->xf[p7O_N][p7O_LOOP] * totp_reciprocal;
      njcp += oxf->xmx[p7X_J][i-1] * oxb->xmx[p7X_J][i] * om->xf[p7O_J][p7O_LOOP] * totp_reciprocal;
      njcp += oxf->xmx[p7X_C][i-1] * oxb->xmx[p7X_C][i] * om->xf[p7O_C][p7O_LOOP] * totp_reciprocal;
      ddef->mocc[i] = 1. - njcp;
    }
  ddef->L = oxf->L;
  return eslOK;
}
/*------------------ end, posterior decoding --------------------*/

/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/
#ifdef p7DECODING_BENCHMARK
/*
   icc  -O3 -static -o decoding_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_BENCHMARK decoding.c -lhmmer -leasel -lm 
   ./decoding_benchmark <hmmfile>     
                    RRM_1 (M=72)       Caudal_act (M=136)     SMC_N (M=1151)
                 -----------------    ------------------     ---------------
   21 Aug 08      3.52u (409 Mc/s)     15.36u (177 Mc/s)     318.78u (72.2 Mc/s)

   The length dependency probably indicates L1 cache missing; because we're 
   manipulating 3 matrices at the same time, we can't fit the calculation 
   in cache.
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
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                  0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for posterior residue decoding, SSE version";

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
  P7_OMX         *fwd     = NULL;
  P7_OMX         *bck     = NULL;
  P7_OMX         *pp      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc;
  double          Mcs;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);                 p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);    p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);    p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  fwd = p7_omx_Create(gm->M, L, L);
  bck = p7_omx_Create(gm->M, L, L);
  pp  = p7_omx_Create(gm->M, L, L);

  esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  p7_Forward (dsq, L, om, fwd,      &fsc);
  p7_Backward(dsq, L, om, fwd, bck, &bsc);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    p7_Decoding(om, fwd, bck, pp);              
  esl_stopwatch_Stop(w);

  Mcs = (double) N * (double) L * (double) gm->M * 1e-6 / (double) w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(pp);
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
#endif /*p7DECODING_BENCHMARK*/
/*------------------ end, benchmark driver ----------------------*/




/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7DECODING_TESTDRIVE

#endif /*p7DECODING_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7DECODING_TESTDRIVE

#endif /*p7DECODING_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/




/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7DECODING_EXAMPLE

#endif /*p7DECODING_EXAMPLE*/
/*------------------------ example ------------------------------*/



/*****************************************************************
 * @LICENSE@
 *****************************************************************/

