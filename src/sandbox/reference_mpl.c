#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "hmmer.h"

#include "sandbox/reference_mpl_fwd.h"

int
p7_ReferenceMPLSearch(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,
		      P7_COORDS2 *dom)
{
  P7_REFMX        *fwd     = p7_refmx_Create(gm->M, L);
  P7_REFMX        *mpl     = p7_refmx_Create(gm->M, L);
  P7_TRACE        *tr      = p7_trace_Create();
  P7_COORDS2_HASH *hashtbl = p7_coords2_hash_Create(0,0,0);
  double          remprob = 1.0;
  float          *wrk     = NULL; /* reusable workspace needed for stochastic trace */
  int      max_iterations = 100000;
  int      iteration;
  float    fwdsc, sc;
  double   mplprob, bestprob;
  int      keyidx, bestidx;
  int      status;

  double   lossthresh = log(0.001);
  int      d;

  /* Calculate a Forward matrix; and the Forward score,
   * the normalization factor for MPL probabilities.
   */
  p7_ReferenceForward(dsq, L, gm, fwd,     &fwdsc);


  /* Labelling implied by the Viterbi path is a good
   * initial guess. Initialize best so far: (bestprob, dom).
   */
  p7_ReferenceViterbi(dsq, L, gm, mpl, tr, &sc);
  p7_trace_Index(tr);
  p7_coords2_SetFromTrace(dom, tr);                 /* sets <dom>      */
  p7_coords2_hash_Store(hashtbl, dom, &bestidx);    /* stores <dom> in the hash */
  p7_refmx_Reuse(mpl);

  p7_ReferenceMPLForward(dsq, L, gm, dom->arr, dom->n, mpl, &sc);   
  bestprob = exp(sc - fwdsc);	                                         /* sets <bestprob> */
  remprob  = 1.0 - bestprob;

  printf("Vit %6.2f %8.2g %8.2g ", sc, bestprob, remprob);
  printf("%2d ", dom->n);
  for (d = 0; d < dom->n; d++) printf("%4d %4d ", dom->arr[d].n1, dom->arr[d].n2);
  printf("\n");

  if (bestprob > remprob) goto DONE; /* we're already done: guaranteed we have the MPL */

  p7_trace_Reuse(tr);
  p7_refmx_Reuse(mpl);

  /* Sample paths from the posterior, to sample labellings.
   */
  for (iteration = 1; iteration <= max_iterations; iteration++)
    {
      p7_reference_trace_Stochastic(rng, &wrk, gm, fwd, tr);
      p7_trace_Index(tr);
      p7_coords2_SetFromTrace(dom, tr);
      status = p7_coords2_hash_Store(hashtbl, dom, &keyidx);
      /* That status is either eslOK or eslEDUP. */

      if (status == eslOK)	/* eslOK = the trace implies a new labeling we haven't seen yet */
	{
	  p7_ReferenceMPLForward(dsq, L, gm, dom->arr, dom->n, mpl, &sc);
	  mplprob  = exp(sc - fwdsc);
	  remprob -= mplprob;

	  printf("NEW %6.2f %8.2g %8.2g ", sc, mplprob, remprob);
	  printf("%2d ", dom->n);
	  for (d = 0; d < dom->n; d++) printf("%4d %4d ", dom->arr[d].n1, dom->arr[d].n2);
	  printf("\n");

	  if (mplprob > bestprob) 
	    {
	      bestprob = mplprob;
	      bestidx  = keyidx;
	    }
	}
      else 
	{
	  printf("dup %6s %8s %8s ", "-", "-", "-");
	  printf("%2d ", dom->n);
	  for (d = 0; d < dom->n; d++) printf("%4d %4d ", dom->arr[d].n1, dom->arr[d].n2);
	  printf("\n");
	}

      if (bestprob > remprob) goto DONE;
      if (iteration * log(1.0 - bestprob) < lossthresh) goto DONE;

      p7_trace_Reuse(tr);
      p7_refmx_Reuse(mpl);
    }
  
 DONE:
  p7_coords2_hash_Get(hashtbl, bestidx, dom);

  if (wrk) free(wrk);
  p7_coords2_hash_Destroy(hashtbl);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(mpl);
  p7_refmx_Destroy(fwd);
  return eslOK;
}


/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7REFERENCE_MPL_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "development of MPL";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(0);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_COORDS2     *dom     = p7_coords2_Create(0,0);
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* Read in one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);
  esl_sqfile_Close(sqfp);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);

  /* Set the profile and null model's target length models */
  p7_bg_SetLength     (bg, sq->n);
  p7_profile_SetLength(gm, sq->n);

  /* Do it. */
  p7_ReferenceMPLSearch(rng, sq->dsq, sq->n, gm, dom);

  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_MPL_EXAMPLE*/




