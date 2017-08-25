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
//#include "sandbox/reference_mpl_fwd.h"

static ESL_OPTIONS options[] = {
  /* name           type           default  env  range  toggles reqs incomp  help                                       docgroup*/
  { (char *) "-h",          eslARG_NONE,   FALSE,  NULL, NULL,   NULL,  NULL, NULL,(char *)  "show brief help on version and usage",                   0 },
  { (char *) "-Z",          eslARG_INT,     (char *)  "1",  NULL, NULL,   NULL,  NULL, NULL,(char *)  "set sequence # to <n>, for E-value calculations",        0 },
  { (char *) "--histplot",  eslARG_OUTFILE, NULL,  NULL, NULL,   NULL,  NULL, NULL,(char *)  "output histograms to <f> (xmgrace format)",              0 },
  {(char *)  "--lhistplot", eslARG_OUTFILE, NULL,  NULL, NULL,   NULL,  NULL, NULL, (char *) "output local-only histograms to <f> (xmgrace format)",   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "testing importance of ensemble calculations";

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
  
  P7_ENGINE      *eng        = NULL;           // Overthruster has <fx>, <cx>, <sm>.

  P7_REFMX       *vit        = p7_refmx_Create(100, 100);
  P7_REFMX       *fwd        = p7_refmx_Create(100, 100);
  P7_REFMX       *bck        = p7_refmx_Create(100, 100);
  P7_REFMX       *pp         = p7_refmx_Create(100, 100);
  //  P7_REFMX       *mpl        = p7_refmx_Create(100, 100);
  P7_TRACE       *tr         = p7_trace_Create();
  P7_COORD2      *dom        = NULL;
  ESL_HISTOGRAM  *invit      = NULL;
  ESL_HISTOGRAM  *outvit     = NULL;
  char           *histfile   = esl_opt_GetString(go, (char *) "--histplot");
  FILE           *histfp     = NULL;
  P7_REFMX       *vit_l      = p7_refmx_Create(100, 100);
  P7_REFMX       *fwd_l      = p7_refmx_Create(100, 100);
  P7_REFMX       *bck_l      = p7_refmx_Create(100, 100);
  P7_REFMX       *pp_l       = p7_refmx_Create(100, 100);
  //  P7_REFMX       *mpl_l      = p7_refmx_Create(100, 100);
  P7_TRACE       *tr_l       = p7_trace_Create();
  P7_COORD2      *dom_l      = NULL;
  ESL_HISTOGRAM  *invit_l    = NULL;
  ESL_HISTOGRAM  *outvit_l   = NULL;
  char           *histfile_l = esl_opt_GetString(go, (char *) "--lhistplot");
  FILE           *histfp_l   = NULL;
  int             Z          = esl_opt_GetInteger(go, (char *) "-Z");
  float           fsc,   vsc,   bsc;
  float           fsc_l, vsc_l, bsc_l;
  int             nintotal,   nin90,   nout90,   nin50,   nout50,   nin0;
  int             nintotal_l, nin90_l, nout90_l, nin50_l, nout50_l, nin0_l;
  int             d;
  float           nullsc;
  int             status;
 
  /* Read in one HMM. Set alphabet to whatever the HMM's alphabet is. */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail((char *) "Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail((char *) "Failed to read HMM");
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
  if      (status == eslENOTFOUND) p7_Fail((char *) "No such file.");
  else if (status == eslEFORMAT)   p7_Fail((char *) "Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail((char *) "Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail((char *) "Open failed, code %d.", status);

  /* Histograms for results */
  invit  = esl_histogram_Create(0.0, 1.0, 0.1);
  outvit = esl_histogram_Create(0.0, 1.0, 0.1);
  if (histfile && (histfp = fopen(histfile, "w")) == NULL) p7_Fail((char *) "Failed to open histogram file %s for writing.", histfile);

  invit_l  = esl_histogram_Create(0.0, 1.0, 0.1);
  outvit_l = esl_histogram_Create(0.0, 1.0, 0.1);
  if (histfile_l && (histfp_l = fopen(histfile_l, "w")) == NULL) p7_Fail((char *) "Failed to open local histogram file %s for writing.", histfile_l);

  eng  = p7_engine_Create(abc, /*params=*/NULL, /*stats=*/NULL, gm->M, /*L_hint=*/400);

  /* For each target sequence... */
  while (( status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_bg_SetLength           (bg,  sq->n);
      p7_profile_SetLength      (gm,  sq->n);
      p7_profile_SetLength      (lgm, sq->n);
      p7_oprofile_ReconfigLength(om,  sq->n);

      if (( status = p7_engine_Overthruster(eng, sq->dsq, sq->n, om, bg)) == eslOK)
	{
	  printf("%-15s  %-30s  ", gm->name, sq->name);

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

	      dom[d].n1  = tr->sqfrom[d];
	      dom[d].n2  = tr->sqto[d];
	      nintotal   += tr->sqto[d] - tr->sqfrom[d] + 1;
	    }
	     
	  //p7_ReferenceMPLForward(sq->dsq, sq->n, gm, dom, tr->ndom, mpl, &mplsc);
	     
	  update_histograms   (pp, dom, tr->ndom, invit, outvit);
	  count_nin_nout_above(pp, dom, tr->ndom, 0.9,   &nin90, &nout90); /* # of res labeled as in domain w/ pp >= 0.9; and outside domain yet pp >= 0.9 */
	  count_nin_nout_above(pp, dom, tr->ndom, 0.5,   &nin50, &nout50); /* ditto, for 0.5 threshold */
	  count_nin_nout_below(pp, dom, tr->ndom, 0.1,   &nin0,  NULL);    /* # of residues labeled by Viterbi as in a domain, yet pp < 0.1 */
	  
	  printf("%10.4g %10.4g %10.4f %10.4f %10.4f %10.4g ",
		 (double) Z * esl_gumbel_surv( (vsc - nullsc) / eslCONST_LOG2, gm->evparam[p7_VMU],  gm->evparam[p7_VLAMBDA]),
		 (double) Z * esl_exp_surv   ( (fsc - nullsc) / eslCONST_LOG2, gm->evparam[p7_FTAU], gm->evparam[p7_FLAMBDA]),
		 vsc, 
		 0.0, // was mplsc
 		 fsc, 
		 0.0);  // was exp(mplsc - fsc));

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

	  //p7_ReferenceMPLForward(sq->dsq, sq->n, lgm, dom_l, tr_l->ndom, mpl_l, &mplsc_l);
	     
	  update_histograms   (pp_l, dom_l, tr_l->ndom, invit_l, outvit_l);
	  count_nin_nout_above(pp_l, dom_l, tr_l->ndom, 0.9,   &nin90_l, &nout90_l); /* # of res labeled as in domain w/ pp >= 0.9; and outside domain yet pp >= 0.9 */
	  count_nin_nout_above(pp_l, dom_l, tr_l->ndom, 0.5,   &nin50_l, &nout50_l); /* ditto, for 0.5 threshold */
	  count_nin_nout_below(pp_l, dom_l, tr_l->ndom, 0.1,   &nin0_l,  NULL);      /* # of residues labeled by Viterbi as in a domain, yet pp < 0.1 */
	  
	  printf("%10.4g %10.4g %10.4f %10.4f %10.4f %10.4g ",
		 (double) Z * esl_gumbel_surv( (vsc_l - nullsc) / eslCONST_LOG2, gm->evparam[p7_VMU],  gm->evparam[p7_VLAMBDA]),
		 (double) Z * esl_exp_surv   ( (fsc_l - nullsc) / eslCONST_LOG2, gm->evparam[p7_FTAU], gm->evparam[p7_FLAMBDA]),
		 vsc_l, 
		 0.0,    // was mplsc_l
		 fsc_l, 
                 0.0);  // was exp(mplsc_l - fsc_l)

	  printf("%5d %5d %5d %5d ",  nintotal_l, nin90_l, nin50_l, nin0_l);
	  printf("%5d %5d %5d ",      (int) sq->n - nintotal_l, nout90_l, nout50_l);
	  printf("\n");


	  p7_trace_Reuse(tr);
	  p7_refmx_Reuse(vit);
	  p7_refmx_Reuse(fwd);
	  p7_refmx_Reuse(bck);
	  p7_refmx_Reuse(pp);
	  //p7_refmx_Reuse(mpl);

	  p7_trace_Reuse(tr_l);
	  p7_refmx_Reuse(vit_l);
	  p7_refmx_Reuse(fwd_l);
	  p7_refmx_Reuse(bck_l);
	  p7_refmx_Reuse(pp_l);
	  //p7_refmx_Reuse(mpl_l);
	}

      p7_engine_Reuse(eng);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) p7_Fail((char *) "Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));
  else if (status != eslEOF)     p7_Fail((char *) "Unexpected error %d reading sequence file %s", status, sqfp->filename);


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
  //p7_refmx_Destroy(mpl);
  p7_refmx_Destroy(pp);
  p7_refmx_Destroy(vit);
  p7_refmx_Destroy(bck);
  p7_refmx_Destroy(fwd);

  if (histfp_l) fclose(histfp_l);
  esl_histogram_Destroy(outvit_l);
  esl_histogram_Destroy(invit_l);
  p7_trace_Destroy(tr_l);
  //p7_refmx_Destroy(mpl_l);
  p7_refmx_Destroy(pp_l);
  p7_refmx_Destroy(vit_l);
  p7_refmx_Destroy(bck_l);
  p7_refmx_Destroy(fwd_l);

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



