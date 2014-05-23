/* Most probable anchor set (MPAS) algorithm
 * Reference implementation: test/development.
 * 
 * Contents:
 *    1. MPAS algorithm 
 *    2. Making trial anchor sets (from traces, for example)
 *    3. Internal (static) functions used by MPAS
 *    4. Statistics collection driver
 *    5. Example
 *    6. License and copyright information.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "base/p7_anchorhash.h"
#include "base/p7_anchors.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_asc_fwdback.h"
#include "dp_reference/reference_trace.h"

#include "dp_reference/reference_anchors.h"


#include "search/p7_mpas.h"

static int dump_current_anchorset(FILE *ofp, const P7_ANCHORS *anch);


/*****************************************************************
 * 1. MPAS algorithm
 *****************************************************************/

/* Function:  p7_reference_Anchors()
 * Synopsis:  The most probable anchor set (MPAS) algorithm.
 *
 * Purpose:   Find the most probable anchor set for comparison of 
 *            query profile <gm> to target digital sequence <dsq>
 *            of length <L>.
 *            
 *            An "anchor set" defines the number and locations of
 *            domains in the target sequence (i.e. homologous
 *            subsequences that each match the model). Each "anchor" is
 *            a correspondence (i0,k0) between a residue i0 in the
 *            sequence and a match state k0 in the profile.  The
 *            probability mass of a domain is defined as the sum of
 *            all paths that pass through the domain's anchor (either
 *            using the glocal MG or the local ML match state). The
 *            probability mass of a domain structure of D domains is
 *            defined as the sum of all paths that pass through each
 *            of the D anchors in separate domains (i.e., paths that
 *            look like SN..NB..(Mk0)..EJ..JB..(Mk0)..EC..CT.)
 *            
 *            The MPAS algorithm depends on stochastic path sampling,
 *            so caller provides a random number generator <rng>.
 *
 *            Caller has already done what we call "First Analysis" on
 *            the sequence: standard Viterbi, Forward/Backward, and
 *            posterior Decoding. From this, caller provides the
 *            Forward matrix <rxf>, the Decoding matrix <rxd>, and the
 *            Viterbi trace <tr>. The matrices will remain unchanged,
 *            returned just as they came in. <tr> will be overwritten
 *            by newly sampled paths.
 *            
 *            Caller provides two workspaces. <byp_wrk> is a "bypass"
 *            workspace that the stochastic trace algorithm knows how
 *            to deal with. A caller will generally set <float *wrk =
 *            NULL>, pass <&wrk> for the argument, and <free(wrk)>
 *            when done. The other workspace is an empty hash table <ah>
 *            that we need for fast checking of which anchor sets
 *            we've already calculated scores for.
 *
 *            For retrieving the result, caller provides two empty
 *            matrices <afu>, <afd>; an empty anchor set <anch>; and
 *            <asc_sc>. Upon return, <afu> and <afd> are the ASC UP
 *            and ASC DOWN matrices for the most probable anchor set;
 *            <anch> is that anchor set; and <asc_sc> is the raw
 *            anchor-set-constrained score in nats.
 *            
 *            Optionally, caller may override default parameters by
 *            passing in a <prm> structure. Pass <NULL> to use
 *            defaults.
 *
 *            Also optionally, caller may provide a pointer to a
 *            <stats> structure, to collect and retrieve statistics
 *            about the MPAS algorithm's performance. Pass <NULL> if
 *            this information isn't needed.
 *            
 *            The statistics in <stats> are collected in two phases.
 *            The first phase is collected here. A second phase
 *            requires having both the Viterbi trace and the anchor
 *            set, in order to compare the two. That requires having
 *            two trace structures around, one for sampling in the
 *            MPAS algorithm and one for remembering the Viterbi
 *            trace. We split the collection into these two phases
 *            because we don't want to force the caller to always have
 *            two trace structures, even if the caller isn't
 *            collecting the extra statistics; this way, the caller
 *            that is interested, can make its own copy of the Viterbi
 *            path. The <stats> structure has flags for whether it has
 *            the data for phase 1, phase 2, or both.
 *            
 *            Collecting <stats> does come with some computational
 *            overhead; avoid it in production code (of course, the
 *            reference implementation isn't supposed to be in
 *            production code in the first place).
 *
 * Args:      rng       : random number generator
 *            dsq       : target sequence, digital (1..L) 
 *            L         : length of <dsq>
 *            gm        : query profile, dual-mode, (1..M)
 *            rxf       : Forward matrix for <gm> x <dsq> comparison
 *            rxd       : Decoding matrix for <gm> x <dsq> comparison
 *            tr        : Viterbi trace for <gm> x <dsq>, then used as space for sampled paths  
 *            byp_wrk   : BYPASS: workspace array for stochastic tracebacks
 *            ah        : empty hash table for alternative <anch> solutions
 *            afu       : empty matrix to be used for ASC UP calculations
 *            afd       : empty matrix to be used for ASC DOWN calculations
 *            anch      : empty anchor data structure to hold result
 *            ret_asc   : score of the most probable anchor set (raw, nats)
 *            prm       : OPTIONAL : non-default control parameters
 *            stats     : OPTIONAL : detailed data collection for development, debugging (or NULL)
 *
 * Returns:   <eslOK> on success, and:
 *            anch    : contains most probable anchor set
 *            afu     : contains ASC Forward UP matrix for MPAS solution
 *            afd     : contains ASC Forward DOWN matrix for it
 *            ret_asc : is the ASC Forward score of it, raw (nats)
 *            stats   : if provided, contains statistics on the optimization that was done
 *            
 *            ah      : undefined (contains hashed data on all anchorsets that were tried)
 *            tr      : undefined (contains last path that that algorithm happened to sample)
 *            rng     : internal state changed (because we sampled numbers from it)
 */
int
p7_reference_Anchors(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
		     const P7_REFMX *rxf, const P7_REFMX *rxd,
		     P7_TRACE *tr,  float **byp_wrk,  P7_ANCHORHASH *ah,
		     P7_REFMX *afu, P7_REFMX *afd, P7_ANCHORS *anch,  float *ret_asc,
		     P7_MPAS_PARAMS *prm, P7_MPAS_STATS *stats)
{
  float           fwdsc          = P7R_XMX(rxf, L, p7R_C) + gm->xsc[p7P_C][p7P_MOVE];  // if this lookup seems iffy, efficiency-wise, make fwdsc an argument; caller has it
  int             max_iterations = (prm ? prm->max_iterations      : 1000);
  float           loglossthresh  = (prm ? log(prm->loss_threshold) : log(0.001));
  int             be_verbose     = (prm ? prm->be_verbose          : FALSE);
  float           best_asc       = -eslINFINITY;
  int             iteration      = 0;   // Counter for stochastic samples. Is 0 for the Viterbi trace, 1st path thru main loop
  float           asc,    best_ascprob;
  int32_t         keyidx, best_keyidx;
  int32_t         asc_keyidx;	        // keyidx that our current <asc>, <afu>, <afd> are for
  int             status;
  
  /* Contract checks / argument validation */
  ESL_DASSERT1(( rxf->type == p7R_FORWARD  ));
  ESL_DASSERT1(( rxd->type == p7R_DECODING ));
  ESL_DASSERT1(( rxf->L == L && rxf->M == gm->M ));
  ESL_DASSERT1(( rxd->L == L && rxd->M == gm->M ));
  ESL_DASSERT1(( tr->L  == L && tr->M  == gm->M ));
  /* afu, afd will be checked, reallocated, initialized in ASCForward() calls */

  /* Initialize optional <stats> */
  if (stats) {
    p7_mpas_stats_Init(stats);
    p7_trace_Score(tr, dsq, gm, &(stats->vsc));
    stats->fsc = fwdsc;
  }

  if (be_verbose) printf("# Forward score: %6.2f nats\n", fwdsc);

  /* MPAS algorithm main loop */
  while (1) // ends on convergence tests at end of loop
    {
      /* First time thru, <tr> is the Viterbi path, from caller; after that, it's a stochastic sampled path */
      p7_reference_anchors_SetFromTrace(rxd, tr, anch);
      status = p7_anchorhash_Store(ah, anch, &keyidx);

      /* <status> is either eslOK or eslEDUP.
       *    <eslOK>   = <anch> is new and needs to be scored;
       *    <eslEDUP> = <anch> is a duplicate, doesn't need scoring, only counting
       */
      if (status == eslOK)
	{
	  /* Get the ASC score for this new anchor set */
	  p7_refmx_Reuse(afu);
	  p7_refmx_Reuse(afd);
	  p7_ReferenceASCForward(dsq, L, gm, anch->a, anch->D, afu, afd, &asc);
	  asc_keyidx = keyidx;
	  
	  if (stats) 
	    {
	      stats->tot_asc_calculations++;
	      if (iteration == 0) {
		stats->vit_asc     = asc;
		stats->vit_ascprob = exp(asc - fwdsc);
	      } 	      
	    }

	  if (be_verbose)
	    {
	      if      (iteration == 0) printf("VIT      %5d %6.2f %10.4g ", keyidx, asc, exp(asc-fwdsc));
	      else if (asc > best_asc) printf("NEW/BEST %5d %6.2f %10.4g ", keyidx, asc, exp(asc-fwdsc));
	      else                     printf("NEW      %5d %6.2f %10.4g ", keyidx, asc, exp(asc-fwdsc));
	      dump_current_anchorset(stdout, anch);
	    }

	  /* If it's better than our best solution so far, store it. */
	  if (asc > best_asc)
	    {
	      best_asc     = asc;
	      best_ascprob = exp(asc - fwdsc);
	      best_keyidx  = keyidx;
	      
	      if (stats) 
		{
		  if (iteration > 0) stats->best_is_viterbi = FALSE;
		  stats->nsamples_in_best = 1;
		}
	    }
	}	  
      else if (status == eslEDUP)
	{
	  if (stats && keyidx == best_keyidx) stats->nsamples_in_best++;

	  if (be_verbose)
	    {
	      printf("dup      %5d %6s %10s ", keyidx, "-", "-");
	      dump_current_anchorset(stdout, anch);
	    }
	}

      /* Convergence / termination tests. See note [1] at end of file for details.  */
      if (best_ascprob >= 0.5 || iteration * log(1.0 - best_ascprob) < loglossthresh)
	break;
      if (iteration == max_iterations) {
	if (stats) stats->solution_not_found = TRUE;
	break; 
      }

      /* Otherwise, sample a new path and loop around again */
      p7_trace_Reuse(tr);
      p7_reference_trace_Stochastic(rng, byp_wrk, gm, rxf, tr);

      iteration++;
    }

  /* <anch> contains the solution corresponding to <keyidx>;
   * <afu>, <afd> contain sol'n corresponding to <asc_keyidx>.
   * If either doesn't match best_keyidx we need to update it.
   */
  if (keyidx     != best_keyidx) p7_anchorhash_Get(ah, best_keyidx, anch);
  if (asc_keyidx != best_keyidx) p7_ReferenceASCForward(dsq, L, gm, anch->a, anch->D, afu, afd, &asc);

  if (stats) 
    {
      stats->tot_iterations = iteration;
      stats->best_asc       = asc;
      stats->best_ascprob   = exp(asc - fwdsc);
      stats->has_part1      = TRUE;
    }

  if (be_verbose)
    {
      printf("WINNER   %5d %6.2f %10.4g ", best_keyidx, best_asc, exp(best_asc-fwdsc));
      dump_current_anchorset(stdout, anch);
    }

  *ret_asc = asc;
  return eslOK;
}

/*****************************************************************
 * Footnotes 
 * Additional note [1]. On the convergence tests.
 *
 * The three convergence tests for the MPAS algorithm work
 * as follows.
 * 
 * 1. i * log(1 - p) < log(t)
 *    Assume that by sampling paths from posterior probability
 *    P(\pi), we're sampling anchor sets A from P(A). (This is
 *    almost, but not quite true, so the test is a heuristic,
 *    not an actual proof.) Suppose there's an AS that's better
 *    than our current best, with probability p >
 *    best_ascprob. The probability that we haven't sampled such
 *    an AS yet, in i iterations, is < (1-best_ascprob)^i.  If
 *    this probability falls below some low/negligible
 *    threshold, we can declare that we've probably found the
 *    solution already: hence, i log (1-best_ascprob) < log t.
 * 
 *    If this test succeeds, we don't know that afu, afd, asc
 *    are on the best solution, so we may need to recalculate
 *    them one last time before being done; hence the 
 *    keyidx != best_keyidx test and the recalculation. (The
 *    other way to do it is to keep the best matrices somewhere,
 *    but we would rather burn a small amount of time than a lot
 *    of memory.
 *
 * 2. best_ascprob >= 0.5 
 *    If this AS happens to dominate, we know we're done (probably
 *    quickly); moreover, this will be triggered immediately after we
 *    calculated a new ASC score, so we know that what's in afu, afd,
 *    asc is from that solution, and keyidx == best_keyidx.
 *    
 * 3. iteration == max_iterations
 *    I believe the MPAS problem is NP-hard. The algorithm usually
 *    works in reasonable time by testing the probabilities above (and
 *    accepting the error probability in test [1]). But we may need to
 *    stop the search and end with the best solution we've found so
 *    far.
 *****************************************************************/

/*-------------------- end, MPAS algorithm ----------------------*/


/*****************************************************************
 * 2. Making trial anchor sets
 *****************************************************************/

/* Function:  p7_reference_anchors_SetFromTrace()
 * Synopsis:  Make an anchor set from a traceback path.
 *
 * Purpose:   Given a path <tr>, and posterior decoding matrix <pp>, for
 *            every domain in <tr> choose the best anchor i,k by
 *            choosing the match state (ML+MG, marginalized) with
 *            highest posterior probability. Put the anchor
 *            coordinates and count into <anch>, an allocated, empty
 *            structure provided by the caller.
 *            
 *            <anch> may be reallocated here, if needed.
 *            
 *            <tr> does not need to be indexed.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure. 
 */
int
p7_reference_anchors_SetFromTrace(const P7_REFMX *pp, const P7_TRACE *tr, P7_ANCHORS *anch)
{
  const float *dpc;
  int          z;		/* index in trace position */
  float        ppv;
  float        best_ppv = -1.;
  int          status;

  p7_anchors_Reuse(anch);
  anch->D = 1;  // we'll increment this every time we find an E state

  for (z = 0; z < tr->N; z++)
    {
      if (p7_trace_IsM(tr->st[z]))
	{
	  dpc = pp->dp[tr->i[z]] + tr->k[z] * p7R_NSCELLS;
	  ppv = dpc[p7R_ML] + dpc[p7R_MG];
	  if (ppv > best_ppv)
	    {
	      anch->a[anch->D].i0 = tr->i[z];
	      anch->a[anch->D].k0 = tr->k[z];
	      best_ppv   = ppv;
	    }
	}
      else if (tr->st[z] == p7T_E)
	{
	  anch->D++;
	  best_ppv = -1.;
	  if ((status = p7_anchors_Grow(anch)) != eslOK) goto ERROR; /* Make sure we have room for another domain */
	}
    }

  anch->D--;  			// because it's D+1 at end of the loop above
  p7_anchor_SetSentinels(anch->a, anch->D, tr->L, tr->M);
  return eslOK;
  
 ERROR:
  return status;
}
/*-------------- end, trial anchor set making -------------------*/



/*****************************************************************
 * 3. Internal (static) functions used by MPAS
 *****************************************************************/

static int
dump_current_anchorset(FILE *ofp, const P7_ANCHORS *anch)
{
  int d;

  fprintf(ofp, "%2d ", anch->D);
  for (d = 1; d <= anch->D; d++)
    fprintf(ofp, "%4d %4d ", anch->a[d].i0, anch->a[d].k0);
  fprintf(ofp, "\n");
  return eslOK;
}
/*--------------------- end, static functions -------------------*/




/*****************************************************************
 * 4. Statistics collection driver.
 *****************************************************************/
#ifdef p7REFERENCE_ANCHORS_STATS
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
  { "-s",          eslARG_INT,      "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                          0 },
  { "-Z",          eslARG_INT,      "1",  NULL, NULL,   NULL,  NULL, NULL, "set sequence # to <n>, for E-value calculations",        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "statistics collection on anchor set constraint calculations";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create( esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  P7_HMMFILE     *hfp     = NULL;
  ESL_ALPHABET   *abc     = NULL;
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_SQ         *sq      = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  ESL_SQFILE     *sqfp    = NULL;
  P7_BG          *bg      = NULL;
  P7_HMM         *hmm     = NULL;
  P7_PROFILE     *gm      = NULL;           /* profile in H4's standard dual-mode local/glocal */
  P7_OPROFILE    *om      = NULL;
  P7_FILTERMX    *fx      = p7_filtermx_Create(100);
  P7_CHECKPTMX   *cx      = p7_checkptmx_Create(100, 100, ESL_MBYTES(p7_RAMLIMIT));
  P7_SPARSEMASK  *sm      = p7_sparsemask_Create(100, 100);
  P7_ANCHORS     *anch    = p7_anchors_Create();
  P7_REFMX       *rxf     = p7_refmx_Create(100,100);
  P7_REFMX       *rxd     = p7_refmx_Create(100,100);
  P7_REFMX       *afu     = p7_refmx_Create(100,100);
  P7_REFMX       *afd     = p7_refmx_Create(100,100);
  P7_TRACE       *tr      = p7_trace_Create();
  P7_TRACE       *vtr     = p7_trace_Create();
  P7_ANCHORHASH  *ah      = p7_anchorhash_Create();
  float          *wrk     = NULL;
  int             Z       = esl_opt_GetInteger(go, "-Z");
  P7_MPAS_STATS   stats;
  float           nullsc, vsc, fsc, asc;
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

	  p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf, vtr, &vsc);
	  p7_reference_trace_Viterbi(gm, rxf, tr);         // a second copy; we have no p7_trace_Copy() function yet
	  p7_ReferenceForward (sq->dsq, sq->n, gm, rxf, &fsc);
	  p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd, NULL);   
	  p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd);   

	  p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, tr, &wrk, ah,
			       afu, afd, anch, &asc, NULL, &stats);
	     
	  p7_trace_Index(vtr);
	  p7_mpas_stats_CompareAS2Trace(&stats, anch, vtr);

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
	  
	  printf("%5d %4d %3s ", 
		 stats.tot_iterations,
		 stats.tot_asc_calculations,
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
	  
	  p7_refmx_Reuse(rxf);
	  p7_refmx_Reuse(rxd);
	  p7_refmx_Reuse(afd);
	  p7_refmx_Reuse(afu);
	  p7_trace_Reuse(tr);
	  p7_trace_Reuse(vtr);
	  p7_anchorhash_Reuse(ah);
	  p7_anchors_Reuse(anch);
	}

      p7_filtermx_Reuse(fx);
      p7_checkptmx_Reuse(cx);
      p7_sparsemask_Reuse(sm);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));
  else if (status != eslEOF)     p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  if (wrk) free(wrk);
  p7_anchors_Destroy(anch);
  p7_anchorhash_Destroy(ah);
  p7_trace_Destroy(tr);
  p7_trace_Destroy(vtr);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(rxf);
  p7_refmx_Destroy(rxd);
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
#endif /*p7REFERENCE_ANCHORS_STATS*/

/*----------------- end, statistics driver ----------------------*/





/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7REFERENCE_ANCHORS_EXAMPLE
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
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of running MPAS (most probable anchor set) algorithm";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
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
  P7_ANCHORS     *anch    = p7_anchors_Create();
  P7_ANCHORHASH  *hashtbl = p7_anchorhash_Create();
  P7_REFMX       *rxf     = NULL;
  P7_REFMX       *rxd     = NULL;
  P7_REFMX       *afu     = NULL;
  P7_REFMX       *afd     = NULL;
  P7_TRACE       *tr      = NULL;
  P7_TRACE       *vtr     = NULL;
  float          *wrk     = NULL;
  P7_MPAS_PARAMS  prm;
  P7_MPAS_STATS   stats;
  float           fsc, vsc, asc;
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

  /* Allocate DP matrices and tracebacks */
  rxf = p7_refmx_Create(gm->M, sq->n);
  rxd = p7_refmx_Create(gm->M, sq->n);
  tr  = p7_trace_Create();
  vtr = p7_trace_Create();
  afu = p7_refmx_Create(gm->M, sq->n);
  afd = p7_refmx_Create(gm->M, sq->n);

  /* First pass analysis */
  p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf, vtr, &vsc);
  p7_reference_trace_Viterbi(gm, rxf, tr);         // a second copy; we have no p7_trace_Copy() function yet
  p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc);
  p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd, NULL);   
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd);   

  /* Customize parameters */
  prm.max_iterations = 1000;
  prm.loss_threshold = 0.001;
  prm.be_verbose     = TRUE;

  /* Do it. */
  p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, tr, &wrk, ah,
		       afu, afd, anch, &asc, &prm, &stats);

  p7_trace_Index(vtr);
  p7_mpas_stats_CompareAS2Trace(&stats, anch, vtr);
  p7_mpas_stats_Dump(stdout, &stats);


  p7_anchors_Destroy(anch);
  p7_anchorhash_Destroy(ah);
  if (wrk) free(wrk);
  p7_trace_Destroy(tr);
  p7_trace_Destroy(vtr);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(rxd);
  p7_refmx_Destroy(rxf);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_ANCHORS_EXAMPLE*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
