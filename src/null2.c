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


/* Function:  p7_null2_ByExpectation()
 * Synopsis:  Calculate biased composition corrections
 * Incept:    SRE, Thu Feb 28 09:52:28 2008 [Janelia]
 *
 * Purpose:   Calculate biased composition corrections for per-domain
 *            and per-sequence scores, for scoring an envelope $i..j$
 *            with model <gm>. The envelope $i..j$ is defined by
 *            passing a <dsq> pointer to residue $i-1$, and an
 *            envelope length <Ld>, where <Ld> $ = j-i+1$.
 *            
 *            The expectation method is applied to envelopes in
 *            simple, well resolved regions (regions containing just a
 *            single envelope, where no stochastic traceback
 *            clustering was required).
 *            
 *            Watch out for some offset coord issues in passing
 *            arguments. <dsq> is the original full length sequence,
 *            <1..L>, and <ienv,jenv> are the $i..j$ boundaries of the
 *            envelope. However, the posterior probability matrix <pp>
 *            has been calculated by the caller for only the envelope,
 *            so its rows are numbered <1..Ld>, for envelope
 *            <ienv..jenv> of length <Ld=jenv-ienv+1>.
 *            
 *            Caller provides a calculated matrix of posterior
 *            probabilities for each residue in <pp>, for the chunk of
 *            digital sequence in the envelope: e.g. for an offset
 *            <dsq+i-1> over length <Ld>. The null2 correction uses
 *            these to calculate a log odds emission score vector as
 *            the posterior weighted average over all emission vectors
 *            used in alignments (inclusive of N,C,J) in the envelope
 *            <i>..<j>, and this null2 is applied to the entire
 *            envelope.
 *            
 *            Caller also provides a workspace DP matrix <wrk> which
 *            contains at least two rows.
 *
 *            The calculation can optionally return the null2 log odds
 *            scores in <opt_null2> as a <0..Kp-1> array containing
 *            $\log \frac{f'{x}}{f{x}}$ scores; caller provides
 *            allocated space for this score vector, if it wants to
 *            see it.
 *
 * Args:      ddef              - where the n2sc null2 score vector is
 *            gm                - profile, in any mode, target length model set to <L>
 *            dsq               - digital seq being scored; <1..L>
 *            ienv              - start of envelope being scored
 *            jenv              - end of envelope being scored                    
 *            pp                - posterior prob matrix, for <gm> against domain envelope <dsq+i-1> (offset)
 *            wrk               - workspace, a DP matrix with at least two rows.
 *            opt_null2         - optRETURN: null2 log odds scores; <0..Kp-1>; caller allocated space
 *
 * Returns:   <eslOK> on success; <opt_null2> contains the null2 scores.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_null2_ByExpectation(P7_DOMAINDEF *ddef, const P7_PROFILE *gm, const ESL_DSQ *dsq, int ienv, int jenv,
		       const P7_GMX *pp, P7_GMX *wrk, float *opt_null2)
{
  int      M      = gm->M;
  int      Ld     = jenv-ienv+1;
  float   *null2  = NULL;
  float  **dp     = wrk->dp;
  float   *xmx    = wrk->xmx;
  float    my_null2[p7_MAXCODE];
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
   * explaining the domain.
   * This is dog-slow; a point for future optimization.
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

 /* make sure it has valid scores for all degeneracies;
  * this is averaging the log odds scores; maybe it should be averaging the probs;
  * but it won't make much difference, degeneracies are rare, we're
  * really just avoiding complete cockup on degenerate residues
  */
  esl_abc_FAvgScVec(gm->abc, null2);

  /* Now assign the null2 scores to residues i..j in the envelope */
  for (i = ienv; i <= jenv; i++)
    ddef->n2sc[i] = null2[dsq[i]];

  /* ta-da */
  return eslOK;
}


/* Function:  p7_null2_BySampling()
 * Synopsis:  Assign null2 scores to an envelope by the sampling method.
 * Incept:    SRE, Thu May  1 10:00:43 2008 [Janelia]
 *
 * Purpose:   Assign null2 scores to the region <ireg..jreg> in the
 *            domain definition object <ddef> by the stochastic traceback
 *            sampling method: construct position-specific null2
 *            scores for each residue in the region by integrating 
 *            over many multidomain null2 models, sampled by
 *            stochastic traceback; each stochastic traceback defines
 *            one or more domains in the region, and a null2 model is
 *            assigned to each individual domains by averaging the
 *            emission distributions of states used in that domain.
 *            
 *            <gm> is the profile; <dsq> is the complete sequence
 *            <1..L>; <gx> is a matrix allocated for enough size that
 *            we can use it for a Forward calculation of <gm> against
 *            a region of length <Lr=jreg-ireg+1>.
 *
 * Args:      
 *
 * Returns:   <eslOK> on success, and the <ddef->n2sc> scores are set
 *            for region <i..j>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_null2_BySampling(P7_DOMAINDEF *ddef, const P7_PROFILE *gm, const ESL_DSQ *dsq, int ireg, int jreg, P7_GMX *gx)
{
  float *suse     = NULL;	/* state usage. suse[0] = any I state; suse[1..M] = indexed M states  */
  float *nsc      = NULL;	/* null scores assigned to each residue i..j; indexed 0..Lr-1         */
  int    Lr       = jreg-ireg+1;/* length of the whole region in residues */
  int    Ld;			/* length of a particular domain, as defined in a sampled trace       */
  float  sc;
  int    t, d;
  float  null2[p7_MAXCODE];	/* a null2 model frequency vector, 0..Kp-1   */
  int    x, k, z, pos;
  int    status;
 
  ESL_ALLOC(suse,     sizeof(float) * (gm->M+1));
  ESL_ALLOC(nsc,      sizeof(float) * Lr);

  p7_GForward(dsq+ireg-1, Lr, gm, gx, &sc); /* this call will disappear, merged into domain definition */

  for (t = 0; t < ddef->nsamples; t++)
    {
      esl_vec_FSet(nsc, Lr, 0.);

      p7_GStochasticTrace(ddef->r, dsq+ireg-1, Lr, gm, gx, ddef->tr);
      p7_trace_Index(ddef->tr);

      pos = 1;
      for (d = 0; d < ddef->tr->ndom; d++) /* For each segment (domain) in this sampled trace: */
	{
	  Ld = ddef->tr->sqto[d] - ddef->tr->sqfrom[d] + 1;
	  
	  /* Calculate emitting state usage in this particular segment: */
	  esl_vec_FSet(suse, gm->M+1, 0.);
	  for (z = ddef->tr->tfrom[d]; z <= ddef->tr->tto[d]; z++) 
	    {
	      if      (ddef->tr->st[z] == p7T_M) suse[ddef->tr->k[z]] += 1.0;
	      else if (ddef->tr->st[z] == p7T_I) suse[0] += 1.0;
	    }
	  esl_vec_FScale(suse, gm->M+1, 1.0 / (float) Ld); /* now suse[] are usage frequencies */
	  
	  /* from those state usages, calculate null2 emission probs
	   * as posterior weighted sum over all emission vectors used
	   * in path segment 
	   */
	  esl_vec_FSet(null2, gm->abc->K, 0.);
	  for (x = 0; x < gm->abc->K; x++)
	    for (k = 1; k <= gm->M; k++)
	      if (suse[k] > 0.) null2[x] +=  exp(p7P_MSC(gm, k, x)) * suse[k];
	  if (suse[0] > 0.)	/* nasty gotcha potential here: watch out for M=1 models with no insert states at all */
	    for (x = 0; x < gm->abc->K; x++)
	      null2[x] += exp(p7P_ISC(gm, 1, x)) * suse[0];    /* uses I1 only: assumes all I emissions are identical */
	  /* null2[x] is now f'(x) / f(x) ratio (not log) */

	  esl_abc_FAvgScVec(gm->abc, null2); /* make sure it has valid scores for all degeneracies: averaged probs */
	  /* note an inconsistency: in ByExpectation, we averaged the log scores; here we avg the ratios; it doesn't
           * seem worth reconciling this small issue, degenerate residues are rare */

	  /* residues outside domains get bumped +1: because f'(x) = f(x), so f'(x)/f(x) = 1 in these segments */
	  for (; pos < ddef->tr->sqfrom[d]; pos++) nsc[pos-1] += 1.0;
	  
	  /* Residues inside domains get bumped by their null2 ratio */
	  for (; pos <= ddef->tr->sqto[d]; pos++)  nsc[pos-1] += null2[dsq[ireg+pos-1]];
	}
      /* the remaining residues in the region outside any domains get +1 */
      for (; pos <= Lr; pos++) nsc[pos-1] += 1.0;	  


      /* increment the position-specific null2 scores in the region by this nsc sample; 
       * note offset -- nsc is 0..Lr-1 for ireg..jreg 
       */
      esl_vec_FAdd(ddef->n2sc+ireg, nsc, Lr);

      /* get ready to take another sample. */
      p7_trace_Reuse(ddef->tr);
    }	  

  /* Convert the tot_nsc[0..Lr-1] ratios to log odds scores for the region */
  for (pos = ireg; pos <= jreg; pos++) 
    ddef->n2sc[pos] = log(ddef->n2sc[pos] / (float) ddef->nsamples);

  free(suse);
  free(nsc);
  return eslOK;

 ERROR:
  if (suse != NULL) free(suse);
  if (nsc  != NULL) free(nsc);
  return status;
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
#include "esl_randomseq.h"
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
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq); /* sample a random digital seq of length L */

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

  esl_rsq_xfIID(r, bg->f, gm->abc->K, L, dsq); /* sample a random digital seq of length L */

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


  

  
   
