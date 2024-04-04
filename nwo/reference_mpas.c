/* Reference implementation of the most probable anchor set (MPAS) algorithm
 *
 * The MPAS algorithm identifies an optimal domain annotation of the sequence;
 * i.e. how many regions in the sequence are homologous to the query profile, and
 * where. It is an ensemble-based algorithm, marginalizing over all alignments
 * consistent with a given annotation, rather than depending on any one particular
 * alignment. A domain annotation for D domains is defined by an "anchor set" of D
 * anchors i0,k0. A path is consistent with an annotation if each and every
 * homologous region in the path aligns an Mk0 state (ML or MG) to residue x_i0.
 *
 * Contents:
 *   1. MPAS algorithm
 *   2. Internal (static) functions
 *   3. Statistics collection driver
 *   4. Example
 * 
 * See also:
 *     h4_mpas.*        Data structures shared between reference, sparse MPAS
 *     reference_asc.*  Reference ASC dynamic programming
 */
#include <h4_config.h>

#include "easel.h"
#include "esl_random.h"

#include "h4_anchorhash.h"
#include "h4_anchorset.h"
#include "h4_mpas.h"
#include "h4_path.h"
#include "h4_refmx.h"

#include "reference_asc.h"
#include "reference_dp.h"
#include "reference_mpas.h"

/*****************************************************************
 * 1. MPAS algorithm
 *****************************************************************/

/* Function:  h4_reference_MPAS()
 * Synopsis:  The most probable anchor set (MPAS) algorithm.
 *
 * Purpose:   Find the most probable anchor set for comparison of query profile <hmm>
 *            in mode <mo> to target digital sequence <dsq> of length <L>.
 *            
 *            An "anchor set" defines the number and locations of domains in the
 *            target sequence (i.e. homologous subsequences that each match the
 *            model). Each "anchor" is a correspondence (i0,k0) between a residue i0
 *            in the sequence and a match state k0 in the profile.  The probability
 *            mass of a domain is defined as the sum of all paths that pass through
 *            the domain's anchor (either using the glocal MG or the local ML match
 *            state). The probability mass of a domain structure of D domains is
 *            defined as the sum of all paths that pass through each of the D anchors
 *            in separate domains (i.e., paths that look like
 *            SN..NB..(Mk0)..EJ..JB..(Mk0)..EC..CT.)
 *            
 *            The MPAS algorithm depends on stochastic path sampling, so caller
 *            provides a random number generator <rng>.
 *
 *            Caller has already done what we call "First Analysis" on the sequence:
 *            standard Viterbi, Forward/Backward, and posterior Decoding. From this,
 *            caller provides the Forward matrix <rxf>, the Decoding matrix <rxd>,
 *            and the Viterbi path <vpi>. The matrices will remain unchanged,
 *            returned just as they came in. <vpi> will be overwritten by newly
 *            sampled paths.
 *            
 *            Caller provides two workspaces. <byp_wrk> is a "bypass" workspace that
 *            the stochastic trace algorithm knows how to deal with. A caller will
 *            generally set <float *wrk = NULL>, pass <&wrk> for the argument, and
 *            <free(wrk)> when done. The other workspace is an empty hash table <ah>
 *            that we need for fast checking of which anchor sets we've already
 *            calculated scores for. Upon return, <ah> contains hashed data on all
 *            anchorsets that were sampled and scored; these data are used in some
 *            testing experiments, but not in production code.
 *
 *            For retrieving the result, caller provides two empty matrices <afu>,
 *            <afd>; an empty anchor set <anch>; and <asc_sc>. Upon return, <afu> and
 *            <afd> are the ASC UP and ASC DOWN matrices for the most probable anchor
 *            set; <anch> is that anchor set; and <asc_sc> is the raw
 *            anchor-set-constrained score in bits.
 *            
 *            Optionally, caller may override default parameters by passing in a
 *            <prm> structure. Pass <NULL> to use defaults.
 *
 *            Also optionally, caller may provide a pointer to a <stats> structure,
 *            to collect and retrieve statistics about the MPAS algorithm's
 *            performance. Pass <NULL> if this information isn't needed.
 *            
 *            The statistics in <stats> are collected in two phases.  The first phase
 *            is collected here. A second phase requires having both the Viterbi
 *            trace and the anchor set, in order to compare the two. That requires
 *            having two trace structures around, one for sampling in the MPAS
 *            algorithm and one for remembering the Viterbi trace. We split the
 *            collection into these two phases because we don't want to force the
 *            caller to always have two trace structures, even if the caller isn't
 *            collecting the extra statistics; this way, the caller that is
 *            interested, can make its own copy of the Viterbi path. The <stats>
 *            structure has flags for whether it has the data for phase 1, phase 2,
 *            or both.
 *            
 *            Collecting <stats> does come with some computational overhead; avoid it
 *            in production code (of course, the reference implementation isn't
 *            supposed to be in production code in the first place).
 *
 * Args:      rng       : random number generator                [internal state changed]
 *            dsq       : target sequence, digital (1..L) 
 *            L         : length of <dsq>
 *            hmm       : query profile, dual-mode, (1..M)
 *            mo        : comparison mode, with length set
 *            rxf       : Forward matrix for <hmm> x <dsq> comparison
 *            rxd       : Decoding matrix for <hmm> x <dsq> comparison
 *            pi        : Viterbi trace for <hmm> x <dsq>, then used as space for sampled paths   [caller-provided space; reused/resized as needed here]
 *            byp_wrk   : BYPASS: workspace array for stochastic tracebacks
 *            ah        : empty hash table for alternative <anch> solutions   [caller-provided space; reused/resized as needed here]
 *            afu       : empty matrix to be used for ASC UP calculations     [caller-provided space; reused/resized as needed here]
 *            afd       : empty matrix to be used for ASC DOWN calculations   [caller-provided space; reused/resized as needed here]
 *            anch      : empty anchor data structure to hold result          [caller-provided space; reused/resized as needed here]
 *            ret_asc   : score of the most probable anchor set (bits)
 *            prm       : OPTIONAL : non-default control parameters (or NULL)
 *            stats     : OPTIONAL : detailed data collection for development, debugging (or NULL)
 *
 * Returns:   <eslOK> on success, and:
 *            anch    : contains most probable anchor set
 *            afu     : contains ASC Forward UP matrix for MPAS solution
 *            afd     : contains ASC Forward DOWN matrix for it
 *            ret_asc : is the ASC Forward score of it, raw (bits)
 *            stats   : if provided, contains statistics on the optimization that was done
 *            ah      : contains hashed data on all anchorsets that were sampled 
 *
 *            pi      : undefined (contains last path that algorithm happened to sample)
 *            rng     : internal state changed (because we sampled numbers from it)
 */
int
h4_reference_MPAS(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo,
                  const H4_REFMX *rxf, const H4_REFMX *rxd, H4_PATH *pi,
                  float **byp_wrk,  H4_ANCHORHASH *ah,
                  H4_REFMX *afu, H4_REFMX *afd, H4_ANCHORSET *anch, float *ret_asc,
                  const H4_MPAS_PARAMS *prm, H4_MPAS_STATS *stats)
{
  float           fwdsc          = H4R_XMX(rxf, L, h4R_C) + mo->xsc[h4_C][h4_MOVE] - mo->nullsc;  
  int             max_iterations = (prm ? prm->max_iterations      : h4_MPAS_MAX_ITERATIONS);
  float           loglossthresh  = (prm ? log(prm->loss_threshold) : log(h4_MPAS_LOSS_THRESHOLD));
  int             nmax_sampling  = (prm ? prm->nmax_sampling       : h4_MPAS_NMAX_SAMPLING);
  int             be_verbose     = (prm ? prm->be_verbose          : h4_MPAS_BE_VERBOSE);
  float           best_asc       = -eslINFINITY;
  int             iteration      = 0;      // Counter for stochastic samples. Is 0 for the Viterbi trace, 1st path thru main loop
  int             satisfied      = FALSE;  // Flips to TRUE when termination conditions satisfied. If <nmax_sampling> is TRUE, we continue sampling anyway.
  float           asc,    best_ascprob;
  int32_t         keyidx, best_keyidx;
  int32_t         asc_keyidx;	        // keyidx that our current <asc>, <afu>, <afd> are for
  int             status;
  
  /* Contract checks / argument validation */
  ESL_DASSERT1(( rxf->type == h4R_FORWARD  ));
  ESL_DASSERT1(( rxd->type == h4R_DECODING ));
  ESL_DASSERT1(( rxf->L == L && rxf->M == hmm->M ));
  ESL_DASSERT1(( rxd->L == L && rxd->M == hmm->M ));
  /* afu, afd will be checked, reallocated, initialized in ASCForward() calls */

  /* Initialize optional <stats> */
  if (stats) {
    h4_mpas_stats_Reuse(stats);
    h4_path_Score(pi, dsq, hmm, mo, &(stats->vsc));
    stats->fsc = fwdsc;
  }

  /* Prepare caller-provided empty/reused allocations 
   *   <anch> Reuse()'d in anchorset_from_path()
   *   <afu>,<afd> Reuse()'d in DP algorithms
   */
  h4_anchorhash_Reuse(ah);
  
  if (be_verbose) printf("# Forward score: %6.2f bits\n", fwdsc);

  /* MPAS algorithm main loop */
  while (1) // ends on convergence tests at end of loop
    {
      /* First time thru, <pi> is the Viterbi path, from caller; after that, it's a stochastic sampled path */
      h4_reference_mpas_path2anchors(pi, rxd, anch);
      status = h4_anchorhash_Store(ah, anch, 0, &keyidx);

      /* now <status> is either eslOK or eslEDUP.
       *    <eslOK>   = <anch> is new and needs to be scored;
       *    <eslEDUP> = <anch> is a duplicate, doesn't need scoring, only counting
       */
      if (status == eslOK)
	{
	  /* Get the ASC score for this new anchor set */
	  h4_reference_asc_Forward(dsq, L, hmm, mo, anch, afu, afd, &asc);
	  asc_keyidx = keyidx;
	  
	  if (stats) 
	    {
	      stats->tot_prob += exp2(asc-fwdsc);
	      stats->tot_asc_calculations++;
	      if (iteration == 0) {  
		stats->vit_asc     = asc;
		stats->vit_ascprob = exp2(asc - fwdsc);
	      } 	      
	    }

	  if (be_verbose)
	    {
	      if      (iteration == 0) printf("VIT      %5d %6.2f %10.4g ", keyidx, asc, exp2(asc-fwdsc));
	      else if (asc > best_asc) printf("NEW/BEST %5d %6.2f %10.4g ", keyidx, asc, exp2(asc-fwdsc));
	      else                     printf("NEW      %5d %6.2f %10.4g ", keyidx, asc, exp2(asc-fwdsc));
	      h4_anchorset_DumpOneLine(stdout, anch);
	    }

	  /* If it's better than our best solution so far, store it. */
	  if (asc > best_asc)
	    {
	      best_asc     = asc;
	      best_ascprob = exp2(asc - fwdsc);
	      best_keyidx  = keyidx;
	      
	      if (stats) 
		{
		  if (iteration > 0) stats->best_is_viterbi = FALSE;
		  if (satisfied)     stats->late_solution   = TRUE;  // Oops. We already thought we found the "best" anchorset yet here's a better one.
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
	      h4_anchorset_DumpOneLine(stdout, anch);
	    }
	}

      /* Convergence / termination tests. See note [1] at end of file for details.  
       * <iteration> ranges from 0..max_iterations; iteration 0 is the Viterbi path
       * and not a stochastic sample; the <iteration> in the termination condition
       * is the # of stochastic samples taken.
       */
      if (best_ascprob >= 0.5 || iteration * log(1.0 - best_ascprob) < loglossthresh)
	{ 
	  if (!satisfied && be_verbose) printf(" --- termination condition satisfied --- \n");
	  satisfied = TRUE;           // Usually we stop, having satisfied the termination tests.
	  if (! nmax_sampling) break; // But in testing, we can set <nmax_sampling> flag to see if 
	}                             // termination tests work, and no better solution is found later.
      if (iteration == max_iterations) {
	if (!satisfied && stats) stats->solution_not_found = TRUE;
	break; 
      }

      /* Otherwise, sample a new path and loop around again */
      h4_anchorset_Reuse(anch);
      h4_path_Reuse(pi);
      h4_reference_StochasticTrace(rng, byp_wrk, hmm, mo, rxf, pi);

      iteration++;
    }

  /* Now <anch> contains the solution corresponding to <keyidx>;
   * <afu>, <afd> contain sol'n corresponding to <asc_keyidx>.
   * If either doesn't match best_keyidx we need to update it.
   */
  if (keyidx     != best_keyidx) h4_anchorhash_Get(ah, best_keyidx, 0, anch);
  if (asc_keyidx != best_keyidx) h4_reference_asc_Forward(dsq, L, hmm, mo, anch, afu, afd, &asc);

  if (stats) 
    {
      stats->tot_iterations = iteration;
      stats->best_asc       = asc;
      stats->best_ascprob   = exp2(asc - fwdsc);
      stats->has_part1      = TRUE;
    }

  if (be_verbose)
    {
      printf("WINNER   %5d %6.2f %10.4g ", best_keyidx, best_asc, exp2(best_asc-fwdsc));
      h4_anchorset_DumpOneLine(stdout, anch);
    }

  *ret_asc = asc;
  return eslOK;
}



/*****************************************************************
 * 2. Internal (static) functions
 *****************************************************************/

/* Function:  h4_reference_mpas_path2anchors()
 * Synopsis:  Heuristically choose a good anchor set from a path.
 * Incept:    SRE, Mon 11 Mar 2024
 *
 * Purpose:   Given a path <pi>, and posterior decoding matrix <rxd>,
 *            for every domain in <pi> choose an anchor i,k by
 *            choosing the match state (ML+MG, marginalized) with
 *            highest posterior probability. Put the anchor
 *            coordinates and count into <anch>, an allocated,
 *            possible empty structure provided by the caller.
 *
 *            It's possible for a glocal domain to have no MG state at
 *            all, and only use {DI} states. Such a domain is
 *            "unanchorable" and skipped. If all domains in the path
 *            are unanchorable, the resulting anchorset is empty
 *            (D=0).
 *            
 *            <anch> will be reinitialized here and may be reallocated
 *            if needed.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure. 
 *
 * Note:      Needs to be here in reference_mpas, not h4_path or
 *            h4_anchorset, because it depends on the reference matrix
 *            structure H4_REFMX.  We'll have a separate function for
 *            the sparse DP version.
 */
int
h4_reference_mpas_path2anchors(const H4_PATH *pi, const H4_REFMX *rxd, H4_ANCHORSET *anch)
{
  const float *dpc;
  int          i,k,z,r;	
  int          i0,k0;
  float        ppv;
  float        best_ppv = -1.;
  int          status;

  if ((status = h4_anchorset_Reuse(anch))                        != eslOK) return status;
  if ((status = h4_anchorset_SetSentinels(anch, rxd->L, rxd->M)) != eslOK) return status;

  i = 1;
  for (z = 0; z < pi->Z; z++)
    {
      if      (pi->st[z] == h4P_N) i += pi->rle[z]-1;
      else if (pi->st[z] == h4P_G) k = 1;          
      else if (pi->st[z] == h4P_L) k = pi->rle[z]; 
      else if (h4_path_IsM(pi->st[z]))
        {
          for (r = 0; r < pi->rle[z]; r++)
            {
              dpc = rxd->dp[i] + k*h4R_NSCELLS;
              ppv = dpc[h4R_ML] + dpc[h4R_MG];
              if (ppv > best_ppv) { i0 = i; k0 = k; best_ppv = ppv; }
              i++; k++;
            }
        }
      else if (h4_path_IsI(pi->st[z])) { i += pi->rle[z]; }
      else if (h4_path_IsD(pi->st[z])) { k += pi->rle[z]; }
      else if (pi->st[z] == h4P_J || pi->st[z] == h4P_C)
        {
          if (best_ppv != -1.) {  // if it's still -1, this domain is unanchorable
            if ((status = h4_anchorset_Add(anch, i0, k0)) != eslOK) return status;
          }
          best_ppv = -1.;
          i += pi->rle[z]-1;
        }

    }
  return eslOK;
}


/*****************************************************************
 * 3. Statistics driver
 *****************************************************************/
#ifdef h4REFERENCE_MPAS_STATS

// TK

#endif // h4REFERENCE_MPAS_STATS


/*****************************************************************
 * 4. Example
 *****************************************************************/
#ifdef h4REFERENCE_MPAS_EXAMPLE
#include <h4_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "general.h"
#include "h4_hmmfile.h"
#include "h4_profile.h"
#include "h4_mode.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                                  docgroup*/
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump report on all sampled anchorsets",                  0 },
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",                   0 },
  { "-k",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "keep sampling until n_max: test termination conditions", 0 },
  { "-n",        eslARG_INT,   "1000", NULL, NULL,   NULL,  NULL, NULL, "maximum number of samples (n_max)",                      0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                          0 },
  { "-t",        eslARG_REAL, "0.001", NULL, NULL,   NULL,  NULL, NULL, "loss threshold",                                         0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show HMMER version info",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of running MPAS (most probable anchor set) algorithm";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  ESL_SQ         *sq      = NULL;
  H4_ANCHORSET   *anch    = h4_anchorset_Create(0,0,0);
  H4_ANCHORHASH  *ah      = h4_anchorhash_Create();
  H4_REFMX       *rxf     = NULL;
  H4_REFMX       *rxd     = NULL;
  H4_REFMX       *afu     = NULL;
  H4_REFMX       *afd     = NULL;
  H4_PATH        *vpi     = h4_path_Create();
  H4_PATH        *pi      = NULL;
  float          *wrk     = NULL;
  H4_MPAS_PARAMS *prm     = h4_mpas_params_Create();
  H4_MPAS_STATS  *stats   = h4_mpas_stats_Create();
  float           fsc, vsc, asc;
  int             status;

  /* Customize MPAS algorithm parameters */
  prm->max_iterations = esl_opt_GetInteger(go, "-n");
  prm->loss_threshold = esl_opt_GetReal   (go, "-t");
  prm->nmax_sampling  = esl_opt_GetBoolean(go, "-k");
  prm->be_verbose     = TRUE;

  /* Read in one HMM */
  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);
 
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such sequence file %s", seqfile);
  else if (status == eslEFORMAT)   esl_fatal("Format of sequence file %s unrecognized.", seqfile);
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);
 
  /* Read one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      esl_fatal("Unexpected error %d reading sequence file %s", status, sqfp->filename);
  esl_sqfile_Close(sqfp);

  /* Set target length model */
  h4_mode_SetLength(mo, sq->n);

  /* Allocate DP matrices */
  rxf = h4_refmx_Create(hmm->M, sq->n);
  rxd = h4_refmx_Create(hmm->M, sq->n);
  afu = h4_refmx_Create(hmm->M, sq->n);
  afd = h4_refmx_Create(hmm->M, sq->n);

  /* First pass analysis */
  h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxf, vpi, &vsc);
  pi = h4_path_Clone(vpi);      // we need to keep the Viterbi path for this example, not just overwrite it in MPAS
  h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc);
  h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxd,      NULL);   
  h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxd, rxd);   

  /* Do it. */
  h4_reference_MPAS(rng, sq->dsq, sq->n, hmm, mo, rxf, rxd, pi, &wrk, ah,
                    afu, afd, anch, &asc, prm, stats);

  h4_mpas_stats_Dump(stdout, stats);

  if (esl_opt_GetBoolean(go, "-a"))
    {
      int key;

      h4_anchorset_Reuse(anch);
      h4_refmx_Reuse(afu);
      h4_refmx_Reuse(afd);

      for (key = 0; key < ah->nkeys; key++)
	{
	  h4_anchorhash_Get(ah, key, 0, anch);
	  h4_reference_asc_Forward(sq->dsq, sq->n, hmm, mo, anch, afu, afd, &asc);

	  printf("anchorset %8d %6.4f %6.4f ",
		 key+1, 
		 exp2(asc - fsc), 
		 (float) ah->key_count[key] / (float) (stats->tot_iterations + 1));
	  h4_anchorset_DumpOneLine(stdout, anch);

	  h4_anchorset_Reuse(anch);
	  h4_refmx_Reuse(afu);
	  h4_refmx_Reuse(afd);
	}
    }


  free(wrk);
  h4_mpas_stats_Destroy(stats);
  h4_mpas_params_Destroy(prm);
  h4_anchorset_Destroy(anch);
  h4_anchorhash_Destroy(ah);
  h4_path_Destroy(pi);    h4_path_Destroy(vpi);
  h4_refmx_Destroy(afd);  h4_refmx_Destroy(afu);
  h4_refmx_Destroy(rxd);  h4_refmx_Destroy(rxf);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif //h4REFERENCE_MPAS_EXAMPLE
