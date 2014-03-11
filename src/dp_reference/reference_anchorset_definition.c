#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "hmmer.h"
#include "dp_reference/reference_asc_forward.h"
#include "dp_reference/reference_anchorset_definition.h"

static int dump_xstats_asc(FILE *ofp, P7_XSTATS_ASC *stats);
static int compare_anchorset_to_trace(P7_TRACE *tr, P7_COORDS2 *anch, P7_XSTATS_ASC *stats);


int
p7_coords2_SetAnchorsFromTrace(const P7_REFMX *pp, const P7_TRACE *tr, P7_COORDS2 *anch)
{
  float *dpc;
  float  ppv;
  int z;
  int k,s;
  int d = 0;
  int   i       = 0;
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
p7_ReferenceASCSearch(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_COORDS2 *anch, P7_XSTATS_ASC *stats)
{
  P7_REFMX       *rxf     = p7_refmx_Create(gm->M, L);
  P7_REFMX       *rxd     = p7_refmx_Create(gm->M, L);
  P7_REFMX       *rxau    = p7_refmx_Create(gm->M, L);
  P7_REFMX       *rxad    = p7_refmx_Create(gm->M, L);
  P7_TRACE       *vtr     = p7_trace_Create();
  P7_TRACE       *tr      = p7_trace_Create();
  P7_COORD2_HASH *hashtbl = p7_coord2_hash_Create(0,0,0);
  float          *wrk     = NULL;                        /* tmp work space needed by stochastic traceback */
  float           fwdsc;
  float           asc, vsc, best_asc;
  float           ascprob,  best_ascprob;
  int32_t         keyidx, best_keyidx;
  int             iteration;
  int             max_iterations = 1000;
  int             d;
  float           lossthresh = log(0.001);
  int             am_done    = FALSE;
  int             be_verbose = TRUE;
  int             status;
  
  stats->tot_iterations       = max_iterations;
  stats->tot_asc_calculations = 0;
  stats->n_to_find            = -1;
  stats->nsamples_in_best     = 0;	/* Viterbi path doesn't count as a sample */
  stats->best_is_viterbi      = TRUE;	/* until it's not */  
  stats->best_found_late      = FALSE;

  /* Calculate a Forward matrix, for sampling traces;
   *   the Forward score, for normalization;
   *   and a Decoding matrix, for posterior (i,k) probabilities.
   */
  p7_ReferenceForward (dsq, L, gm, rxf, &fwdsc);
  p7_ReferenceBackward(dsq, L, gm, rxd, NULL);   
  p7_ReferenceDecoding(dsq, L, gm, rxf, rxd, rxd);   

  /* Labelling implied by the Viterbi path is a good initial guess. */
  p7_ReferenceViterbi(dsq, L, gm, rxau, vtr, &vsc);
  p7_trace_Index(vtr);
  p7_coords2_SetAnchorsFromTrace(rxd, vtr, anch);
  p7_coord2_hash_Store(hashtbl, anch->arr, anch->n, &best_keyidx);

  p7_refmx_Reuse(rxau);
  /* don't reuse Viterbi trace! we need it again at end - just keep it */

  p7_ReferenceASCForward(dsq, L, gm, anch->arr, anch->n, rxau, rxad, &best_asc);
  best_ascprob            = exp(best_asc - fwdsc);

  stats->fsc         = fwdsc;
  stats->vsc         = vsc;
  stats->vit_asc     = best_asc;
  stats->vit_ascprob = best_ascprob;

  if (be_verbose) 
    {
      printf("# Forward score: %6.2f nats\n", fwdsc);
      printf("# Viterbi score: %6.2f nats\n", vsc);

      printf("VIT %6.2f %8.4g ", asc, best_ascprob);
      printf("%2d ", anch->n);
      for (d = 0; d < anch->n; d++) printf("%4d %4d ", anch->arr[d].n1, anch->arr[d].n2);
      printf("\n");
    }

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
	  
	  if (! am_done) stats->tot_asc_calculations++;

	  if (be_verbose) 
	    {
	      printf("NEW %6.2f %8.4g ", asc, ascprob);
	      printf("%2d ", anch->n);
	      for (d = 0; d < anch->n; d++) printf("%4d %4d ", anch->arr[d].n1, anch->arr[d].n2);
	      printf("\n");
	    }

	  if (ascprob > best_ascprob)
	    {
	      best_asc         = asc;
	      best_ascprob     = ascprob;
	      best_keyidx      = keyidx;

	      stats->nsamples_in_best = 1;
	      stats->best_is_viterbi  = FALSE;
	      if (am_done) stats->best_found_late = TRUE;
	    }

	  p7_refmx_Reuse(rxau);
	  p7_refmx_Reuse(rxad);
	}
      else 
	{
	  if (keyidx == best_keyidx) stats->nsamples_in_best++;

	  if (be_verbose)
	    {
	      printf("dup %6s %8s %8s ", "-", "-", "-");
	      printf("%2d ", anch->n);
	      for (d = 0; d < anch->n; d++) printf("%4d %4d ", anch->arr[d].n1, anch->arr[d].n2);
	      printf("\n");
	    }
	}
      
      if (! am_done && iteration * log(1.0 - best_ascprob) < lossthresh) 
	{
	  am_done          = TRUE;
	  stats->n_to_find = iteration;
	  if (be_verbose) printf("### I think I'm done.\n");
	}

      p7_coords2_Reuse(anch);
      p7_trace_Reuse(tr);
    }

  // DONE:
  p7_coord2_hash_Get(hashtbl, L, best_keyidx, anch);
  stats->solution_not_found = ( am_done ? FALSE : TRUE );
  stats->best_asc     = best_asc;
  stats->best_ascprob = best_ascprob;
  compare_anchorset_to_trace(vtr, anch, stats);

  if (be_verbose)
    {
      printf("BEST %6.2f %8.4g ", stats->best_asc, stats->best_ascprob);
      printf("%2d ", anch->n);
      for (d = 0; d < anch->n; d++) printf("%4d %4d ", anch->arr[d].n1, anch->arr[d].n2);
      printf("\n");
    }

  free(wrk);
  p7_coord2_hash_Destroy(hashtbl);
  p7_trace_Destroy(vtr);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(rxad);
  p7_refmx_Destroy(rxau);
  p7_refmx_Destroy(rxd);
  p7_refmx_Destroy(rxf);
  return eslOK;
}


/* trace must be indexed */
static int
compare_anchorset_to_trace(P7_TRACE *tr, P7_COORDS2 *anch, P7_XSTATS_ASC *stats)
{
  int ad;
  int td              = 0;
  int anch_in_this_td = 0;

  stats->anch_outside    = 0;
  stats->anch_unique     = 0;
  stats->anch_multiple   = 0;

  stats->dom_zero        = 0;
  stats->dom_one         = 0;
  stats->dom_multiple    = 0;

  /* For n domains in tr:
   *   they can either be hit 0 times, 1 time, or 2+ times by anchors.
   * For m anchors in anch:
   *   they can either fall outside any domain, uniquely in a domain, or multiply in a domain.
   */
  for (ad = 0; ad < anch->n; ad++)
    {
      if   (anch->arr[ad].n1 < tr->sqfrom[td] || td == tr->ndom) 
	stats->anch_outside++;
      else if (anch->arr[ad].n1 >= tr->sqfrom[td] && anch->arr[ad].n1 <= tr->sqto[td])
	anch_in_this_td++;
      else 
	{
	  /* we have to advance <td>, and try again */
	  if      (anch_in_this_td == 0) { stats->dom_zero++; }
	  else if (anch_in_this_td == 1) { stats->anch_unique++; stats->dom_one++; }
	  else if (anch_in_this_td > 1)  { stats->anch_multiple += anch_in_this_td; stats->dom_multiple++; }
	  anch_in_this_td = 0;
	  td++;
	  ad--;			/* forces reevaluation of <ad> when we go back around; a bit hacky! */
	}
    }
  
  /* we're out of anchors. If td == tr->ndom, we also know we
   * handled what happened with anchors in the last <td>. But if
   * td == tr->ndom-1, we haven't yet resolved what happened with final <td> yet,
   * and if td is even smaller, we have some dom_zero's to count.
   */
  for (; td < tr->ndom; td++)
    {
      if      (anch_in_this_td == 0) { stats->dom_zero++; }
      else if (anch_in_this_td == 1) { stats->anch_unique++; stats->dom_one++; }
      else if (anch_in_this_td > 1)  { stats->anch_multiple += anch_in_this_td; stats->dom_multiple++; }
      anch_in_this_td = 0;
    }

  ESL_DASSERT1(( stats->dom_zero     + stats->dom_one     + stats->dom_multiple  == tr->ndom ));
  ESL_DASSERT1(( stats->anch_outside + stats->anch_unique + stats->anch_multiple == anch->n  ));

  return eslOK;
}

static int
dump_xstats_asc(FILE *ofp, P7_XSTATS_ASC *stats)
{
  fprintf(ofp, "# Stats on the ASC solution:\n");
  fprintf(ofp, "# tot_iterations       = %d\n",    stats->tot_iterations);
  fprintf(ofp, "# tot_asc_calculations = %d\n",    stats->tot_asc_calculations);
  fprintf(ofp, "# n_to_find            = %d\n",    stats->n_to_find);
  fprintf(ofp, "# Viterbi score (nats) = %.2f\n",  stats->vsc);
  fprintf(ofp, "# Forward score (nats) = %.2f\n",  stats->fsc);
  fprintf(ofp, "# best_asc             = %.2f\n",  stats->best_asc);
  fprintf(ofp, "# vit_asc              = %.2f\n",  stats->vit_asc);
  fprintf(ofp, "# vit_ascprob          = %6.4f\n", stats->vit_ascprob);
  fprintf(ofp, "# best_ascprob         = %6.4f\n", stats->best_ascprob);
  fprintf(ofp, "# nsamples_in_best     = %d\n",    stats->nsamples_in_best);
  fprintf(ofp, "# best_is_viterbi      = %s\n",    (stats->best_is_viterbi    ? "TRUE" : "FALSE"));
  fprintf(ofp, "# best_found_late      = %s\n",    (stats->best_found_late    ? "TRUE" : "FALSE"));
  fprintf(ofp, "# solution_not_found   = %s\n",    (stats->solution_not_found ? "TRUE" : "FALSE"));
  fprintf(ofp, "# anch_outside         = %d\n",    stats->anch_outside);
  fprintf(ofp, "# anch_unique          = %d\n",    stats->anch_unique);
  fprintf(ofp, "# anch_multiple        = %d\n",    stats->anch_multiple);
  fprintf(ofp, "# dom_zero             = %d\n",    stats->dom_zero);
  fprintf(ofp, "# dom_one              = %d\n",    stats->dom_one);
  fprintf(ofp, "# dom_multiple         = %d\n",    stats->dom_multiple);
  return eslOK;
}

/*****************************************************************
 * x. Statistics collection driver.
 *****************************************************************/
#ifdef p7REFERENCE_ASC_STATS
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type           default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,   FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",                   0 },
  { "-Z",          eslARG_INT,      "1",  NULL, NULL,   NULL,  NULL, NULL, "set sequence # to <n>, for E-value calculations",        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "statistics collection on anchor set constraint calculations";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go         = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng        = esl_randomness_Create(0);
  char           *hmmfile    = esl_opt_GetArg(go, 1);
  P7_HMMFILE     *hfp        = NULL;
  ESL_ALPHABET   *abc        = NULL;
  char           *seqfile    = esl_opt_GetArg(go, 2);
  ESL_SQ         *sq         = NULL;
  int             format     = eslSQFILE_UNKNOWN;
  ESL_SQFILE     *sqfp       = NULL;
  P7_BG          *bg         = NULL;
  P7_HMM         *hmm        = NULL;
  P7_PROFILE     *gm         = NULL;           /* profile in H4's standard dual-mode local/glocal */
  P7_OPROFILE    *om         = NULL;
  P7_FILTERMX    *fx         = p7_filtermx_Create(100);
  P7_CHECKPTMX   *cx         = p7_checkptmx_Create(100, 100, ESL_MBYTES(p7_RAMLIMIT));
  P7_SPARSEMASK  *sm         = p7_sparsemask_Create(100, 100);
  P7_COORDS2     *anch       = p7_coords2_Create(0,0);
  P7_XSTATS_ASC   stats;
  int             Z          = esl_opt_GetInteger(go, "-Z");
  float           nullsc;
  int             status;

  /* Read in one HMM. Set alphabet to whatever the HMM's alphabet is. */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Configure vector, dual-mode, and local-only profiles from HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);

  om = p7_oprofile_Create(hmm->M, abc);
  p7_oprofile_Convert(gm, om);
  
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* For each target sequence... */
  while (( status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_bg_SetLength           (bg,  sq->n);
      p7_profile_SetLength      (gm,  sq->n);
      p7_oprofile_ReconfigLength(om,  sq->n);

      if (( status = p7_pipeline_AccelerationFilter(sq->dsq, sq->n, om, bg, fx, cx, sm)) == eslOK)
	{
	  printf("%-15s  %-30s  ", gm->name, sq->name);

	  p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

	  p7_ReferenceASCSearch(rng, sq->dsq, sq->n, gm, anch, &stats);
	     
	  printf("%8.2f %10.4g %8.2f %10.4g ",
		 stats.fsc, 
		 (double) Z * esl_exp_surv   ( (stats.fsc - nullsc) / eslCONST_LOG2, gm->evparam[p7_FTAU], gm->evparam[p7_FLAMBDA]),
		 stats.vsc, 
		 (double) Z * esl_gumbel_surv( (stats.vsc - nullsc) / eslCONST_LOG2, gm->evparam[p7_VMU],  gm->evparam[p7_VLAMBDA]));
		   
		 
	  printf("%8.2f %6.4f %3s ", 
		 stats.vit_asc,
		 stats.vit_ascprob,
		 (stats.best_is_viterbi    ? "YES" : "no"));

	  printf("%8.2f %6.4f %6.4f ", 
		 stats.best_asc,
		 stats.best_ascprob,
		 (float) stats.nsamples_in_best / (float) stats.tot_iterations);
	  
	  printf("%4d %4d %3s %3s ", 
		 stats.n_to_find,
		 stats.tot_asc_calculations,
		 (stats.best_found_late    ? "YES" : "no"),
		 (stats.solution_not_found ? "YES" : "no"));
	    
	  printf("%4d %4d %4d ", 
		 stats.anch_outside,
		 stats.anch_unique,
		 stats.anch_multiple);

	  printf("%4d %4d %4d ",
		 stats.dom_zero,
		 stats.dom_one,
		 stats.dom_multiple);

	  printf("\n");
	}

      p7_filtermx_Reuse(fx);
      p7_checkptmx_Reuse(cx);
      p7_sparsemask_Reuse(sm);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));
  else if (status != eslEOF)     p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(cx);
  p7_filtermx_Destroy(fx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7REFERENCE_ASC_STATS*/




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
  P7_XSTATS_ASC   stats;
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
  p7_ReferenceASCSearch(rng, sq->dsq, sq->n, gm, anch, &stats);

  dump_xstats_asc(stdout, &stats);

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
