/* Collect data on the effectiveness of using ensemble calculations
 * as opposed to optimal alignment; search a sequence database with
 * one profile, collect various statistics.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_histogram.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "sandbox/reference_mpl_fwd.h"

static ESL_OPTIONS options[] = {
  /* name           type           default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,   FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",                   0 },
  { "-Z",          eslARG_INT,      "1",  NULL, NULL,   NULL,  NULL, NULL, "set sequence # to <n>, for E-value calculations",        0 },
  { "--histplot",  eslARG_OUTFILE, NULL,  NULL, NULL,   NULL,  NULL, NULL, "output histograms to <f> (xmgrace format)",              0 },
  { "--lhistplot", eslARG_OUTFILE, NULL,  NULL, NULL,   NULL,  NULL, NULL, "output local-only histograms to <f> (xmgrace format)",   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "testing importance of ensemble calculations";

static int acceleration_filter(ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_BG *bg,
			       P7_FILTERMX *fx, P7_CHECKPTMX *cx, P7_SPARSEMASK *sm);
static int update_histograms(P7_REFMX *pp, P7_COORD2 *dom, int ndom, ESL_HISTOGRAM *indom, ESL_HISTOGRAM *outdom);
static int count_nin_nout_above(P7_REFMX *pp, P7_COORD2 *dom, int ndom, float thresh, int *opt_nin, int *opt_nout);
static int count_nin_nout_below(P7_REFMX *pp, P7_COORD2 *dom, int ndom, float thresh, int *opt_nin, int *opt_nout);


int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go         = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
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
  P7_PROFILE     *lgm        = NULL;           /* profile in local-only mode, emulating H3        */
  P7_OPROFILE    *om         = NULL;
  P7_FILTERMX    *fx         = p7_filtermx_Create(100);
  P7_CHECKPTMX   *cx         = p7_checkptmx_Create(100, 100, ESL_MBYTES(p7_SPARSIFY_RAMLIMIT));
  P7_SPARSEMASK  *sm         = p7_sparsemask_Create(100, 100);
  P7_REFMX       *vit        = p7_refmx_Create(100, 100);
  P7_REFMX       *fwd        = p7_refmx_Create(100, 100);
  P7_REFMX       *bck        = p7_refmx_Create(100, 100);
  P7_REFMX       *pp         = p7_refmx_Create(100, 100);
  P7_REFMX       *mpl        = p7_refmx_Create(100, 100);
  P7_TRACE       *tr         = p7_trace_Create();
  P7_COORD2      *dom        = NULL;
  ESL_HISTOGRAM  *invit      = NULL;
  ESL_HISTOGRAM  *outvit     = NULL;
  char           *histfile   = esl_opt_GetString(go, "--histplot");
  FILE           *histfp     = NULL;
  P7_REFMX       *vit_l      = p7_refmx_Create(100, 100);
  P7_REFMX       *fwd_l      = p7_refmx_Create(100, 100);
  P7_REFMX       *bck_l      = p7_refmx_Create(100, 100);
  P7_REFMX       *pp_l       = p7_refmx_Create(100, 100);
  P7_REFMX       *mpl_l      = p7_refmx_Create(100, 100);
  P7_TRACE       *tr_l       = p7_trace_Create();
  P7_COORD2      *dom_l      = NULL;
  ESL_HISTOGRAM  *invit_l    = NULL;
  ESL_HISTOGRAM  *outvit_l   = NULL;
  char           *histfile_l = esl_opt_GetString(go, "--lhistplot");
  FILE           *histfp_l   = NULL;
  int             Z          = esl_opt_GetInteger(go, "-Z");
  float           fsc,   vsc,   bsc,   mplsc;
  float           fsc_l, vsc_l, bsc_l, mplsc_l;
  int             nintotal,   nin90,   nout90,   nin50,   nout50,   nin0;
  int             nintotal_l, nin90_l, nout90_l, nin50_l, nout50_l, nin0_l;
  int             d;
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

  lgm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(lgm, hmm, bg, 100);

  om = p7_oprofile_Create(hmm->M, abc);
  p7_oprofile_Convert(gm, om);
  
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* Histograms for results */
  invit  = esl_histogram_Create(0.0, 1.0, 0.1);
  outvit = esl_histogram_Create(0.0, 1.0, 0.1);
  if (histfile && (histfp = fopen(histfile, "w")) == NULL) p7_Fail("Failed to open histogram file %s for writing.", histfile);

  invit_l  = esl_histogram_Create(0.0, 1.0, 0.1);
  outvit_l = esl_histogram_Create(0.0, 1.0, 0.1);
  if (histfile_l && (histfp_l = fopen(histfile_l, "w")) == NULL) p7_Fail("Failed to open local histogram file %s for writing.", histfile_l);

  /* For each target sequence... */
  while (( status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_bg_SetLength           (bg,  sq->n);
      p7_profile_SetLength      (gm,  sq->n);
      p7_profile_SetLength      (lgm, sq->n);
      p7_oprofile_ReconfigLength(om,  sq->n);

      if (( status = acceleration_filter(sq->dsq, sq->n, om, bg, fx, cx, sm)) == eslOK)
	{
	  printf("%15s  %30s  ", gm->name, sq->name);

	  p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

	  /* glocal/local dual-mode */
	  p7_ReferenceViterbi (sq->dsq, sq->n, gm, vit, tr,  &vsc);
	  p7_ReferenceForward (sq->dsq, sq->n, gm, fwd,      &fsc);   
	  p7_ReferenceBackward(sq->dsq, sq->n, gm, bck,      &bsc);
	  p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, bck, pp);

	  /* Create a coord2 list, crudely. We can do better later. */
	  p7_trace_Index(tr);
	  ESL_ALLOC(dom, sizeof(P7_COORD2) * tr->ndom);
	  nintotal = 0;
	  for (d = 0; d < tr->ndom; d++)
	    {
	      //printf("domain %2d:  seq pos %5d..%-5d  model pos %5d..%5d\n",
	      //d+1, tr->sqfrom[d], tr->sqto[d], tr->hmmfrom[d], tr->hmmto[d]);

	      dom[d].n1 = tr->sqfrom[d];
	      dom[d].n2   = tr->sqto[d];
	      nintotal    += tr->sqto[d] - tr->sqfrom[d] + 1;
	    }
	     
	  p7_ReferenceMPLForward(sq->dsq, sq->n, gm, dom, tr->ndom, mpl, &mplsc);
	     
	  update_histograms   (pp, dom, tr->ndom, invit, outvit);
	  count_nin_nout_above(pp, dom, tr->ndom, 0.9,   &nin90, &nout90); /* # of res labeled as in domain w/ pp >= 0.9; and outside domain yet pp >= 0.9 */
	  count_nin_nout_above(pp, dom, tr->ndom, 0.5,   &nin50, &nout50); /* ditto, for 0.5 threshold */
	  count_nin_nout_below(pp, dom, tr->ndom, 0.1,   &nin0,  NULL);    /* # of residues labeled by Viterbi as in a domain, yet pp < 0.1 */
	  
	  printf("%10.4g %10.4g %10.4f %10.4f %10.4f %10.4g ",
		 (double) Z * esl_gumbel_surv( (vsc - nullsc) / eslCONST_LOG2, gm->evparam[p7_VMU],  gm->evparam[p7_VLAMBDA]),
		 (double) Z * esl_exp_surv   ( (fsc - nullsc) / eslCONST_LOG2, gm->evparam[p7_FTAU], gm->evparam[p7_FLAMBDA]),
		 vsc, 
		 mplsc,
		 fsc, 
		 exp(mplsc - fsc));

	  printf("%5d %5d %5d %5d ",  nintotal, nin90, nin50, nin0);
	  printf("%5d %5d %5d ",      (int) sq->n - nintotal, nout90, nout50);


	  /* now repeat, w/ local-only model */
	  p7_ReferenceViterbi (sq->dsq, sq->n, lgm, vit_l, tr_l,  &vsc_l);
	  p7_ReferenceForward (sq->dsq, sq->n, lgm, fwd_l,        &fsc_l);   
	  p7_ReferenceBackward(sq->dsq, sq->n, lgm, bck_l,        &bsc_l);
	  p7_ReferenceDecoding(sq->dsq, sq->n, lgm, fwd_l, bck_l, pp_l);

	  p7_trace_Index(tr_l);
	  ESL_ALLOC(dom_l, sizeof(P7_COORD2) * tr_l->ndom);
	  nintotal_l = 0;
	  for (d = 0; d < tr_l->ndom; d++)
	    {
	      //printf("local %2d:  seq pos %5d..%-5d  model pos %5d..%5d\n",
	      //d+1, tr_l->sqfrom[d], tr_l->sqto[d], tr_l->hmmfrom[d], tr_l->hmmto[d]);
	      dom_l[d].n1 = tr_l->sqfrom[d];
	      dom_l[d].n2   = tr_l->sqto[d];
	      nintotal_l    += tr_l->sqto[d] - tr_l->sqfrom[d] + 1;
	    }

	  p7_ReferenceMPLForward(sq->dsq, sq->n, lgm, dom_l, tr_l->ndom, mpl_l, &mplsc_l);
	     
	  update_histograms   (pp_l, dom_l, tr_l->ndom, invit_l, outvit_l);
	  count_nin_nout_above(pp_l, dom_l, tr_l->ndom, 0.9,   &nin90_l, &nout90_l); /* # of res labeled as in domain w/ pp >= 0.9; and outside domain yet pp >= 0.9 */
	  count_nin_nout_above(pp_l, dom_l, tr_l->ndom, 0.5,   &nin50_l, &nout50_l); /* ditto, for 0.5 threshold */
	  count_nin_nout_below(pp_l, dom_l, tr_l->ndom, 0.1,   &nin0_l,  NULL);    /* # of residues labeled by Viterbi as in a domain, yet pp < 0.1 */
	  
	  printf("%10.4g %10.4g %10.4f %10.4f %10.4f %10.4g ",
		 (double) Z * esl_gumbel_surv( (vsc_l - nullsc) / eslCONST_LOG2, gm->evparam[p7_VMU],  gm->evparam[p7_VLAMBDA]),
		 (double) Z * esl_exp_surv   ( (fsc_l - nullsc) / eslCONST_LOG2, gm->evparam[p7_FTAU], gm->evparam[p7_FLAMBDA]),
		 vsc_l, 
		 mplsc_l,
		 fsc_l, 
		 exp(mplsc_l - fsc_l));

	  printf("%5d %5d %5d %5d ",  nintotal_l, nin90_l, nin50_l, nin0_l);
	  printf("%5d %5d %5d ",      (int) sq->n - nintotal_l, nout90_l, nout50_l);
	  printf("\n");


	  p7_trace_Reuse(tr);
	  p7_refmx_Reuse(vit);
	  p7_refmx_Reuse(fwd);
	  p7_refmx_Reuse(bck);
	  p7_refmx_Reuse(pp);
	  p7_refmx_Reuse(mpl);

	  p7_trace_Reuse(tr_l);
	  p7_refmx_Reuse(vit_l);
	  p7_refmx_Reuse(fwd_l);
	  p7_refmx_Reuse(bck_l);
	  p7_refmx_Reuse(pp_l);
	  p7_refmx_Reuse(mpl_l);
	}

      p7_filtermx_Reuse(fx);
      p7_checkptmx_Reuse(cx);
      p7_sparsemask_Reuse(sm);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));
  else if (status != eslEOF)     p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);


  /* Results output
   */
  if (histfp) 
    {
      esl_histogram_Plot(histfp, invit);
      esl_histogram_Plot(histfp, outvit);
    }
  if (histfp_l) 
    {
      esl_histogram_Plot(histfp_l, invit_l);
      esl_histogram_Plot(histfp_l, outvit_l);
    }

  if (histfp) fclose(histfp);
  esl_histogram_Destroy(outvit);
  esl_histogram_Destroy(invit);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(mpl);
  p7_refmx_Destroy(pp);
  p7_refmx_Destroy(vit);
  p7_refmx_Destroy(bck);
  p7_refmx_Destroy(fwd);
  p7_sparsemask_Destroy(sm);

  if (histfp_l) fclose(histfp_l);
  esl_histogram_Destroy(outvit_l);
  esl_histogram_Destroy(invit_l);
  p7_trace_Destroy(tr_l);
  p7_refmx_Destroy(mpl_l);
  p7_refmx_Destroy(pp_l);
  p7_refmx_Destroy(vit_l);
  p7_refmx_Destroy(bck_l);
  p7_refmx_Destroy(fwd_l);

  p7_checkptmx_Destroy(cx);
  p7_filtermx_Destroy(fx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(lgm);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;

 ERROR:
  return status;
}


/* acceleration_filter()
 * Simplified version of the standard HMMER acceleration filter of [Eddy11].
 *
 * Purpose:  Implements a simpified version of the HMMER acceleration
 *           pipeline in <p7_Pipeline()>: only the vector acceleration
 *           steps, and without the bias filter. Compares vector
 *           profile <om> against digital sequence <dsq> of length
 *           <L>, calculating bit scores against null model
 *           <bg>. Returns <eslOK> if the comparison passes the
 *           acceleration filters, and <eslFAIL> if it does not.
 * 
 *           Caller provides vector DP matrices: <fx> for SSV/MSV and
 *           <cx> for checkpointed Forward/Backward. These DP matrices
 *           can have any previous allocation size because they will
 *           be reallocated as needed by the vector DP routines.  Upon
 *           return, <fx> contains the Viterbi filter data from
 *           <p7_ViterbiFilter()>, and <cx> contains the checkpointed
 *           DP Forward/Backward matrix from <p7_ForwardFilter()> and
 *           <p7_BackwardFilter()>.
 *           
 *           Caller also provides a sparse mask object <sm>, which
 *           upon successful return will contain the sparse mask. It
 *           too can be provided in any valid allocation size, because
 *           it is reallocated as needed.
 *           
 *           Acceleration filter thresholds are set to defaults: P
 *           $\leq$ 0.02 for MSV, P $\leq$ 1e-3 for Viterbi, and P
 *           $\leq$ 1e-5 for Forward filters.
 *           
 * Args:     dsq  : digital sequence, 1..L
 *           L    : length of <dsq>
 *           om   : vector profile, with length parameterization set
 *           bg   : null model, with length parameterization set
 *           fx   : vector MSV/Vit DP matrix, at any valid alloc size
 *           cx   : vector checkpointed Fwd/Bck DP mx, any valid size
 *           sm   : RETURN: sparse mask
 *
 * Returns:  <eslOK> if <dsq> is a high-scoring target that passes the
 *           acceleration filter thresholds. <sm> contains the sparse
 *           DP mask for the <dsq>/<om> comparison. <fx> and <cx> DP
 *           matrices are valid too, if caller wants to use them for 
 *           something.
 *           
 *           <eslFAIL> if <dsq> is a low-scoring target that fails the
 *           filter thresholds. Now the data in <fx>, <cx>, <sm> are
 *           undefined though they remain validly allocated and can be
 *           reused.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      J12/120.
 * 
 * To do:     We could put this in p7_Pipeline itself and streamline code
 *            there. We could add pipeline stats as an optional
 *            argument, add the bias filter, add an optional input arg
 *            for the F1/F2/F3 parameters, and do something about
 *            SCAN_MODELS mode -- perhaps just making two versions of
 *            this call, one for search mode and one for scan mode.
 */
static int
acceleration_filter(ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_BG *bg,
		    P7_FILTERMX *fx, P7_CHECKPTMX *cx, P7_SPARSEMASK *sm)
{
  float  usc, vitsc, fwdsc;
  float  nullsc;
  float  seq_score;
  double P;
  float  F1 = 0.02;
  float  F2 = 1e-3;
  float  F3 = 1e-5;

  if (L == 0) return eslFAIL;

  p7_bg_NullOne(bg, dsq, L, &nullsc);

  p7_MSVFilter(dsq, L, om, fx, &usc);
  seq_score = (usc - nullsc) / eslCONST_LOG2;
  P = esl_gumbel_surv(seq_score, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);
  if (P > F1) return eslFAIL;

  p7_ViterbiFilter(dsq, L, om, fx, &vitsc);  
  seq_score = (vitsc - nullsc) / eslCONST_LOG2;
  P  = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
  if (P > F2) return eslFAIL;
  
  p7_ForwardFilter(dsq, L, om, cx, &fwdsc);
  seq_score = (fwdsc - nullsc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  if (P > F3) return eslFAIL;
 
  p7_BackwardFilter(dsq, L, om, cx, sm, p7_SPARSIFY_THRESH);

  return eslOK;
}
  


static int
update_histograms(P7_REFMX *pp, P7_COORD2 *dom, int ndom, ESL_HISTOGRAM *indom, ESL_HISTOGRAM *outdom)
{
  float phomology;
  int   i;
  int   d = 0;			/* index of current or next domain in <dom> */

  for (i = 1; i <= pp->L; i++)
    {
      phomology = 1.0 - (P7R_XMX(pp, i, p7R_N) + P7R_XMX(pp, i, p7R_JJ) + P7R_XMX(pp, i, p7R_CC)); /* JJ,CC, not J,C because we only want emitting J,C */
      
      /* Test for whether i is in the current or next domain:
       *   if there's a domain d to come, but i isn't in it: i < dom[d].n1
       *   if i is beyond all domain segments:               d >= ndom
       *   and the d >= ndom test has to be first, to avoid evaluating dom[d].n1 for invalid d
       */
      if (d >= ndom || i < dom[d].n1) esl_histogram_Add(outdom, phomology); /* i is outside annotated domains */
      else                               esl_histogram_Add(indom,  phomology); /* i is within annotated domain   */
      
      if (d < ndom && i == dom[d].n2) d++;
    }
  return eslOK;
}  

static int 
count_nin_nout_above(P7_REFMX *pp, P7_COORD2 *dom, int ndom, float thresh, int *opt_nin, int *opt_nout)
{
  float phomology;
  int   i;
  int   nin  = 0;
  int   nout = 0;
  int   d    = 0;			/* index of current or next domain in <dom> */


  for (i = 1; i <= pp->L; i++)
    {
      phomology = 1.0 - (P7R_XMX(pp, i, p7R_N) + P7R_XMX(pp, i, p7R_JJ) + P7R_XMX(pp, i, p7R_CC)); /* JJ,CC, not J,C because we only want emitting J,C */
      
      if (phomology >= thresh)
	{
	  if (d >= ndom || i < dom[d].n1) nout++;
	  else                               nin++;
	}

      if (d < ndom && i == dom[d].n2) d++;
    }

  if (opt_nin)  *opt_nin  = nin;
  if (opt_nout) *opt_nout = nout;
  return eslOK;
}

static int 
count_nin_nout_below(P7_REFMX *pp, P7_COORD2 *dom, int ndom, float thresh, int *opt_nin, int *opt_nout)
{
  float phomology;
  int   i;
  int   nin  = 0;
  int   nout = 0;
  int   d    = 0;			/* index of current or next domain in <dom> */

  for (i = 1; i <= pp->L; i++)
    {
      phomology = 1.0 - (P7R_XMX(pp, i, p7R_N) + P7R_XMX(pp, i, p7R_JJ) + P7R_XMX(pp, i, p7R_CC)); /* JJ,CC, not J,C because we only want emitting J,C */
      
      if (phomology < thresh)
	{
	  if (d >= ndom || i < dom[d].n1) nout++;
	  else                               nin++;
	}

      if (d < ndom && i == dom[d].n2) d++;
    }

  if (opt_nin)  *opt_nin  = nin;
  if (opt_nout) *opt_nout = nout;
  return eslOK;
}



/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 * 
 * Incept: SRE, Tue Jan  7 09:40:47 2014 [Janelia Farm] 
 *****************************************************************/



