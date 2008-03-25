/* "null2" model: biased composition correction
 * 
 * See p7_domaindef.c -- null2 correction of per-seq and per-domain
 * scores is embedded within p7_domaindef's logic; we split it out
 * to a separate file because it's so important.
 * 
 * SRE, Thu Feb 28 09:51:27 2008 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"
#include "easel.h"
#include "esl_vectorops.h"
#include "hmmer.h"

#define MMX(i,k)      (dp[(i)][(k) * p7G_NSCELLS + p7G_M])
#define IMX(i,k)      (dp[(i)][(k) * p7G_NSCELLS + p7G_I])
#define DMX(i,k)      (dp[(i)][(k) * p7G_NSCELLS + p7G_D])
#define XMX(i,s)      (xmx[(i) * p7G_NXCELLS + (s)])


/* Function:  p7_Null2Corrections()
 * Synopsis:  Calculate biased composition corrections
 * Incept:    SRE, Thu Feb 28 09:52:28 2008 [Janelia]
 *
 * Purpose:   Calculate biased composition corrections for per-domain
 *            and per-sequence scores, for scoring an envelope $i..j$
 *            with model <gm>. The envelope $i..j$ is defined by
 *            passing a <dsq> pointer to residue $i-1$, and an
 *            envelope length <Ld>, where <Ld> $ = j-i+1$.
 *            
 *            (Why not pass the original full-length <dsq> and $i,j$
 *            envelope coords, instead of the offset <dsq> pointer?
 *            Because our caller has already calculated small Forward
 *            and Backward matrices and posterior probabilities for
 *            the offset <dsq>, so we want to work within this offset
 *            coordinate frame.)
 *            
 *            Caller provides a calculated matrix of posterior
 *            probabilities for each residue in <pp>, for the chunk of
 *            digital sequence in the envelope: e.g. the same offset
 *            <dsq+i-1> pointer, over length <Ld>. The null2
 *            correction uses these to calculate a log odds emission
 *            score vector as the posterior weighted average over all
 *            emission vectors used in alignments (inclusive of N,C,J)
 *            in the envelope <i>..<j>.
 *            
 *            The caller provides <noverlap>, the number of residues
 *            that overlap with a previous envelope. This is used to
 *            avoid overcounting null2 corrections against the
 *            per-sequence score. If no preceding envelope overlaps
 *            <i>..<j>, <noverlap> is 0.
 *            
 *            Caller also provides a workspace DP matrix <wrk> which
 *            contains at least two rows.
 *
 *            The calculation optionally returns the null2 log odds
 *            scores, and two corrections. The null2 score vector
 *            <opt_null2> is a <0..K-1> array containing $\log
 *            \frac{f_d{x}}{f{x}}$ scores; caller provides allocated
 *            space for this score vector, if it wants to see it. The
 *            per-domain correction <*opt_domcorrection> is the score
 *            (in nats) that the caller should subtract from the
 *            domain envelope score for <i>..<j>). The per-sequence
 *            correction <*opt_seqcorrection> is the score (in nats)
 *            that the caller should subtract from the overall
 *            per-sequence score, in addition to the null1 score.
 *            
 *            If the envelope <i..j> is independent (nonoverlapping)
 *            with any previous domain, these two corrections are the
 *            same number; when the envelope overlaps with a previous
 *            envelope, the per-seq correction only accounts for the
 *            last <Ld - noverlap> residues of the envelope, even
 *            though the null model is calculated over the entire
 *            envelope. This isn't well principled; it's just a way of
 *            making sure each residue only gets corrected once in the
 *            per-sequence score.
 *
 * Args:      gm                - profile, in any mode, target length model set to <L>
 *            dsq               - offset digital seq being scored; <dsq+i-1>; <1..Ld>
 *            Ld                - length of domain envelope
 *            noverlap          - number of residues that overlap with a previous envelope
 *            pp                - posterior prob matrix, for <gm> against domain envelope in <dsq>
 *            wrk               - workspace, a DP matrix with at least two rows.
 *            opt_null2         - optRETURN: null2 log odds scores; <0..K-1>; caller allocated space
 *            opt_domcorrection - optRETURN: per-domain score correction (in nats) to subtract
 *            opt_seqcorrection - optRETURN: per-seq score correction (in nats) to subtract
 *
 * Returns:   <eslOK> on success; <opt_null2> contains the null2 scores;
 *            <*opt_domcorrection>, <*opt_seqcorrection> contain
 *            the corrections.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
p7_Null2Corrections(const P7_PROFILE *gm, const ESL_DSQ *dsq, int Ld, int noverlap,
		    const P7_GMX *pp, P7_GMX *wrk,
		    float *opt_null2, float *opt_domcorrection, float *opt_seqcorrection)
{
  int      M   = gm->M;
  float  **dp  = wrk->dp;
  float   *xmx = wrk->xmx;
  float    my_null2[p7_MAXABET];
  float   *null2  = NULL;
  float    domsc, seqsc;
  float    omega;
  float    xfactor;
  int      x;			/* over symbols 0..K-1                       */
  int      i;			/* over offset envelope dsq positions 1..Ld  */
  int      k;			/* over model M states 1..M, I states 1..M-1 */

  /* Where are we storing the null2 log odds scores? Either caller provided
   * space, or we use our own. */
  null2 = (opt_null2 != NULL ? opt_null2 : my_null2);

  /* Calculate expected # of times that each emitting state was used
   * in generating the Ld residues in this domain.
   * The 0 row in <wrk> is used to hold these numbers.
   */
  esl_vec_FCopy(pp->dp[1],            (M+1)*p7G_NSCELLS, wrk->dp[0]); 
  esl_vec_FCopy(pp->xmx+p7G_NXCELLS,  p7G_NXCELLS,       wrk->xmx);  
  for (i = 2; i <= Ld; i++)
    {
      esl_vec_FAdd(wrk->dp[0], pp->dp[i],             (M+1)*p7G_NSCELLS);
      esl_vec_FAdd(wrk->xmx,   pp->xmx+i*p7G_NXCELLS, p7G_NXCELLS);
    }
  
  /* Convert those expected #'s to log frequencies; these we'll use as
   * the log posterior weights.
   */
  esl_vec_FLog(wrk->dp[0], (M+1)*p7G_NSCELLS);
  esl_vec_FLog(wrk->xmx,   p7G_NXCELLS);

  esl_vec_FIncrement(wrk->dp[0], (M+1)*p7G_NSCELLS, -log((float)Ld));
  esl_vec_FIncrement(wrk->xmx,   p7G_NXCELLS,       -log((float)Ld));

  /* Calculate null2's log odds emission probabilities, by taking
   * posterior weighted sum over all emission vectors used in paths
   * explaining the domain.  This is dog-slow; a point for future
   * optimization
   */
  xfactor = XMX(0,p7G_N);
  xfactor = p7_FLogsum(xfactor, XMX(0,p7G_C));
  xfactor = p7_FLogsum(xfactor, XMX(0,p7G_J));

  for (x = 0; x < gm->abc->K; x++)
    { 
      esl_vec_FCopy(wrk->dp[0], (M+1)*p7G_NSCELLS, wrk->dp[1]);
      for (k = 1; k < M; k++)
	{
	  MMX(1,k) += p7P_MSC(gm, k, x);
	  IMX(1,k) += p7P_ISC(gm, k, x);
	}
      MMX(1,M) += p7P_MSC(gm, k, x);

      null2[x] = esl_vec_FLogSum(wrk->dp[1], (M+1)*p7G_NSCELLS);
      null2[x] = p7_FLogsum(null2[x], xfactor);
    }
  /* now null2[x] = \log \frac{f_d(x)}{f_0(x)} for all x in alphabet,
   * 0..K-1, where f_d(x) are the ad hoc "null2" residue frequencies
   * for this envelope.
   */

  /* Tote up the null2 path scores for the envelope */
  noverlap = (noverlap > Ld ? Ld : noverlap); /* just checking. */
  domsc = seqsc = 0.0;
  for (i = 1; i <= noverlap; i++) domsc += null2[dsq[i]];
  for ( ;     i <= Ld;       i++) seqsc += null2[dsq[i]];
  domsc += seqsc;

  /* Sum over the null1, null2 paths for this domain */
  omega = 1.0f / 256.0f;
  if (opt_domcorrection != NULL) *opt_domcorrection = p7_FLogsum( log (1. - omega), log(omega) + domsc);
  if (opt_seqcorrection != NULL) *opt_seqcorrection = p7_FLogsum( log (1. - omega), log(omega) + seqsc);
  return eslOK;
}





/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/

/* A speed benchmark for the null2 calculation. */
/* Code is related to the domaindef benchmark.  */

/* With gcc -g code, RRM_1 model: 
 *    8.74u - 2.50u = 6.24u     original version
 *    5.15u - 2.51u = 2.64u     rearranged to precalculate expected # of state usage in envelope
 *    5.12u - 5.02u = 0.10u     rearranged to take post probs as argument, and to use esl_vec_LogSum
 */


#ifdef p7NULL2_BENCHMARK
/* gcc -o null2_benchmark -g -Wall -I../easel -L../easel -I. -L. -Dp7NULL2_BENCHMARK null2.c -lhmmer -leasel -lm
 * ./null2_benchmark <hmmfile> 
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_alphabet.h"
#include "esl_stopwatch.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-N",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "number of sampled sequences",                      0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "length of sampled sequences",                      0 },
  { "-b",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "baseline timing - no null2 processing",            0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                  0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "timing benchmark for the null2 correction calculation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go          = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_RANDOMNESS *r           = NULL;
  ESL_STOPWATCH  *w           = esl_stopwatch_Create();
  char           *hmmfile     = esl_opt_GetArg(go, 1);
  P7_HMMFILE     *hfp         = NULL;
  P7_HMM         *hmm         = NULL;
  ESL_ALPHABET   *abc         = NULL;
  P7_BG          *bg          = NULL;
  P7_PROFILE     *gm          = NULL;
  P7_GMX         *fwd         = NULL;
  P7_GMX         *bck         = NULL;
  ESL_DSQ        *dsq         = NULL;
  int             N           = esl_opt_GetInteger(go, "-N");
  int             L           = esl_opt_GetInteger(go, "-L");
  int             do_baseline = esl_opt_GetBoolean(go, "-b");
  int             idx;

  /* Set up the RNG */
  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  /* Read in one HMM */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);

  /* Other initial allocations */
  dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  fwd  = p7_gmx_Create(gm->M, L);
  bck  = p7_gmx_Create(gm->M, L);
  p7_FLogsumInit();

  esl_stopwatch_Start(w);
  for (idx = 0; idx < N; idx++)
    {
      esl_rnd_xfIID(r, bg->f, abc->K, L, dsq); /* sample a random digital seq of length L */

      p7_GForward (dsq, L, gm, fwd, NULL); 
      p7_GBackward(dsq, L, gm, bck, NULL);       

      p7_PosteriorDecoding(L, gm, fwd, bck, bck); /* <bck> now contains posterior probs */
      
      if (! do_baseline) 
	p7_Null2Corrections(gm, dsq, L, 0, bck, fwd, NULL, NULL, NULL);	/* use <fwd> as workspace */
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "CPU time:   ");

  free(dsq);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7NULL2_BENCHMARK*/


/*****************************************************************
 * Unit tests
 *****************************************************************/
#ifdef p7NULL2_TESTDRIVE

static void
utest_correct_normalization(ESL_RANDOMNESS *r, P7_PROFILE *gm, P7_BG *bg, ESL_DSQ *dsq, int L, P7_GMX *fwd, P7_GMX *bck)
{
  char *msg = "normalization unit test failed";
  float null2[p7_MAXABET];
  float sum;
  int   x;

  esl_rnd_xfIID(r, bg->f, gm->abc->K, L, dsq); /* sample a random digital seq of length L */

  p7_GForward (dsq, L, gm, fwd, NULL); 
  p7_GBackward(dsq, L, gm, bck, NULL);       
  p7_PosteriorDecoding(L, gm, fwd, bck, bck); /* <bck> now contains posterior probs */
  p7_Null2Corrections(gm, dsq, L, 0, bck, fwd, null2, NULL, NULL);	/* use <fwd> as workspace */

  /* Convert null2 from lod score to frequencies f_d  */
  for (x = 0; x < gm->abc->K; x++)
    null2[x] = exp(null2[x]) * bg->f[x];

  sum = esl_vec_FSum(null2, gm->abc->K);
  if (sum < 0.99 || sum > 1.01) esl_fatal(msg);
}  
  

#endif /*p7NULL2_TESTDRIVE*/


/*****************************************************************
 * Test driver
 *****************************************************************/

#ifdef p7NULL2_TESTDRIVE
/* gcc -o null2_utest -g -Wall -I../easel -L../easel -I. -L. -Dp7NULL2_TESTDRIVE null2.c -lhmmer -leasel -lm
 * ./null2_utest
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_alphabet.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "length of sampled sequences",                      0 },
  { "-M",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "length of sampled HMM",                            0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                  0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options]";
static char banner[] = "unit test driver for the null2 correction calculation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go          = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r           = NULL;
  ESL_ALPHABET   *abc         = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm         = NULL;
  P7_BG          *bg          = NULL;
  P7_PROFILE     *gm          = NULL;
  P7_GMX         *fwd         = NULL;
  P7_GMX         *bck         = NULL;
  ESL_DSQ        *dsq         = NULL;
  int             M           = esl_opt_GetInteger(go, "-M");
  int             L           = esl_opt_GetInteger(go, "-L");

  /* Set up the RNG */
  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  /* Sample a random HMM */
  p7_hmm_Sample(r, M, abc, &hmm);

  /* Configure a profile from the sampled HMM */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);

  /* Other initial allocations */
  dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  fwd  = p7_gmx_Create(gm->M, L);
  bck  = p7_gmx_Create(gm->M, L);
  p7_FLogsumInit();

  utest_correct_normalization(r, gm, bg, dsq, L, fwd, bck);

  free(dsq);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7NULL2_TESTDRIVE*/


  

  
   
