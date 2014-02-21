#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "hmmer.h"
#include "sandbox/reference_asc_forward.h"

int
p7_coords2_SetAnchorsFromTrace(const P7_REFMX *pp, const P7_TRACE *tr, P7_COORDS2 *anch)
{
  float *dpc;
  float  ppv;
  int z;
  int i,k,s;
  int d = 0;
  float best_pp = -1.;
  int status;

  if (( status = p7_coords2_GrowTo(anch, tr->ndom) ) != eslOK) goto ERROR;

  for (z = 0; z < tr->N; z++)
    {
      if (tr->i[z]) i = tr->i[z]; /* keeps track of last i emitted, for when a D state is best */

      if (p7_trace_IsMain(tr->st[z]))
	{
	  k   = tr->k[z];
	  dpc = pp->dp[i] + k * p7R_NSCELLS;
	  ppv = 0.0;
	  for (s = 0; s < p7R_NSCELLS; s++) ppv += dpc[s];
	  if (ppv > best_pp)
	    {
	      anch->arr[d].n1 = i;
	      anch->arr[d].n2 = k;
	      best_pp         = ppv;
	    }
	}
      else if (tr->st[z] == p7T_E)
	{
	  d++;
	  best_pp = -1.;
	}
    }

  anch->n    = d;
  anch->dim1 = pp->L;
  anch->dim2 = pp->M;
  return eslOK;
  
 ERROR:
  return status;
}


int
p7_ReferenceASCSearch(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_COORDS2 *anch)
{
  P7_REFMX       *rxf     = p7_refmx_Create(gm->M, L);
  P7_REFMX       *rxd     = p7_refmx_Create(gm->M, L);
  P7_REFMX       *rxau    = p7_refmx_Create(gm->M, L);
  P7_REFMX       *rxad    = p7_refmx_Create(gm->M, L);
  P7_TRACE       *tr      = p7_trace_Create();
  P7_COORD2_HASH *hashtbl = p7_coord2_hash_Create(0,0,0);
  float          *wrk     = NULL;                        /* tmp work space needed by stochastic traceback */
  float           fwdsc;
  float           asc, vsc;
  float           ascprob;
  float           best_ascprob;
  int32_t         keyidx, best_keyidx;
  int             iteration;
  int             max_iterations = 1000;
  int             d;
  float           lossthresh = log(0.001);
  int             am_done    = FALSE;
  int             status;
  
  /* Calculate a Forward matrix, for sampling traces;
   *   the Forward score, for normalization;
   *   and a Decoding matrix, for posterior (i,k) probabilities.
   */
  p7_ReferenceForward (dsq, L, gm, rxf, &fwdsc);
  p7_ReferenceBackward(dsq, L, gm, rxd, NULL);   
  p7_ReferenceDecoding(dsq, L, gm, rxf, rxd, rxd);   

  /* Labelling implied by the Viterbi path is a good initial guess. */
  p7_ReferenceViterbi(dsq, L, gm, rxau, tr, &vsc);
  p7_trace_Index(tr);
  p7_coords2_SetAnchorsFromTrace(rxd, tr, anch);
  p7_coord2_hash_Store(hashtbl, anch->arr, anch->n, &best_keyidx);

  p7_refmx_Reuse(rxau);
  p7_trace_Reuse(tr);

  p7_ReferenceASCForward(dsq, L, gm, anch->arr, anch->n, rxau, rxad, &asc);
  best_ascprob = exp(asc - fwdsc);

  printf("# Forward score: %6.2f nats\n", fwdsc);
  printf("# Viterbi score: %6.2f nats\n", vsc);

  printf("VIT %6.2f %8.4g ", asc, best_ascprob);
  printf("%2d ", anch->n);
  for (d = 0; d < anch->n; d++) printf("%4d %4d ", anch->arr[d].n1, anch->arr[d].n2);
  printf("\n");

  /* Sample paths from the posterior, to sample anchor sets
   */
  for (iteration = 1; iteration <= max_iterations; iteration++)
    {
      p7_reference_trace_Stochastic(rng, &wrk, gm, rxf, tr);
      p7_trace_Index(tr);
      p7_coords2_SetAnchorsFromTrace(rxd, tr, anch);
      
      status = p7_coord2_hash_Store(hashtbl, anch->arr, anch->n, &keyidx);
      /* status = eslOK means it's a new anchor set;
       *          eslEDUP means we've seen this anchor set before
       */

      if (status == eslOK)
	{
	  p7_ReferenceASCForward(dsq, L, gm, anch->arr, anch->n, rxau, rxad, &asc);
	  ascprob  = exp(asc - fwdsc);

	  printf("NEW %6.2f %8.4g ", asc, ascprob);
	  printf("%2d ", anch->n);
	  for (d = 0; d < anch->n; d++) printf("%4d %4d ", anch->arr[d].n1, anch->arr[d].n2);
	  printf("\n");

	  if (ascprob > best_ascprob)
	    {
	      best_ascprob = ascprob;
	      best_keyidx  = keyidx;
	    }

	  p7_refmx_Reuse(rxau);
	  p7_refmx_Reuse(rxad);
	}
      else 
	{
	  printf("dup %6s %8s %8s ", "-", "-", "-");
	  printf("%2d ", anch->n);
	  for (d = 0; d < anch->n; d++) printf("%4d %4d ", anch->arr[d].n1, anch->arr[d].n2);
	  printf("\n");
	}
      
      if (! am_done && iteration * log(1.0 - best_ascprob) < lossthresh) 
	{
	  am_done = TRUE;
	  printf("### I think I'm done.\n");
	}

      p7_coords2_Reuse(anch);
      p7_trace_Reuse(tr);
    }

 DONE:
  p7_coord2_hash_Get(hashtbl, L, best_keyidx, anch);

  free(wrk);
  p7_coord2_hash_Destroy(hashtbl);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(rxad);
  p7_refmx_Destroy(rxau);
  p7_refmx_Destroy(rxd);
  p7_refmx_Destroy(rxf);
  return eslOK;
}




/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7REFERENCE_ASC_EXAMPLE
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
static char banner[] = "example of ASC search mplementation";

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
  P7_COORDS2     *anch    = p7_coords2_Create(0,0);
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
 
  /* Read one sequence */
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
  p7_ReferenceASCSearch(rng, sq->dsq, sq->n, gm, anch);


  p7_coords2_Destroy(anch);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_ASC_FWD_EXAMPLE*/
