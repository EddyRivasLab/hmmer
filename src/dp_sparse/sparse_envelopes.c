#include "p7_config.h"

#include "easel.h"

#include "base/p7_envelopes.h" 
#include "base/p7_anchors.h"
#include "base/p7_profile.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_asc_fwdback.h"


int
p7_sparse_Envelopes(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,
		    const P7_ANCHOR *anch, int D, 
		    P7_SPARSEMX *asf, const P7_SPARSEMX *asd,
		    P7_ENVELOPES *env)
{
  const P7_SPARSEMASK *sm = asd->sm;
  const float         *xD = asd->xmx;
  const float         *xF = asf->xmx;
  float theta   = 0.5;
  float epsilon = 0.005;
  P7_ANCHOR reanch[3];  
  float pB, pE, pL, pG;
  int   d,g,i,s;
  int   status;

  if (( status = p7_envelopes_Reinit(env, D)) != eslOK) goto ERROR;
  for (d = 1; d <= D; d++)
    {
      env->arr[d].i0     = anch[d].i0;
      env->arr[d].k0     = anch[d].k0;
      env->arr[d].alia   = 0;
      env->arr[d].alib   = 0;
      env->arr[d].ka     = 0;
      env->arr[d].kb     = 0;
      env->arr[d].ia     = 0;
      env->arr[d].ib     = 0;
      env->arr[d].oea    = 0;
      env->arr[d].oeb    = 0;
      env->arr[d].env_sc = 0.;
      env->arr[d].flags  = 0;
    }
  env->D = D;
  env->L = L;
  env->M = gm->M;
  p7_envelope_SetSentinels(env->arr, D, L, gm->M);  


  d = 1;
  for (g = 1; g <= sm->S; g++)
    {
      if (env->arr[d].i0 > sm->seg[g].ib) continue;

      pB = pL = pG = pE = 0.;
      for (i = sm->seg[g].ia-1; i <= sm->seg[g].ib; i++, xD += p7S_NXCELLS, xF += p7S_NXCELLS)
	{
	  if (i == env->arr[d].i0) 
	    { 
	      if (pG >= pL) env->arr[d].flags |= p7E_IS_GLOCAL;
	      d++; 
	      pB = pL = pG = pE = 0.; 
	    }

	  if (env->arr[d-1].i0 >= sm->seg[g].ia)  // in DOWN sector
	    {
	      pE += xD[p7S_E];
	      if (1.0 - pE < theta && ! env->arr[d-1].ib)  
		env->arr[d-1].ib = i;
	      if (1.0 - pE < epsilon   && ! env->arr[d-1].oeb) 
		{
		  env->arr[d-1].oeb = i;
		  s = (d == D+1 ? p7S_C : p7S_J);
		  if (xD[s] < 1.0 - (2.*epsilon)) env->arr[d-1].flags &= (~p7E_ENVSC_APPROX);
		  
		  env->arr[d-1].env_sc += xF[s];
		  env->arr[d-1].env_sc += gm->xsc[p7P_C][p7P_LOOP] * (L-env->arr[d].oeb) + gm->xsc[p7P_C][p7P_MOVE];
		}
	    }

	  /* In UP sector.
	   *  Once we've identified ia, we know we've identified oea too,
           *  because oea <= ia (because epsilon <= theta).
	   */
	  if (env->arr[d].i0 <= sm->seg[g].ib && ! env->arr[d].ia)
	    {
	      pB += xD[p7S_B];
	      pL += xD[p7S_L];
	      pG += xD[p7S_G];

	      if (pB >= theta)
		env->arr[d].ia = i+1;

	      if (pB >= epsilon   && ! env->arr[d].oea)
		{
		  env->arr[d].oea = i+1;
		  s = (d == 1 ? p7S_N : p7S_J);
		  if (xD[s] >= 1.0 - 2.*epsilon) env->arr[d].flags |= p7E_ENVSC_APPROX;
		  
		  env->arr[d].env_sc -= xF[s];
		  env->arr[d].env_sc += gm->xsc[p7P_N][p7P_LOOP] * (env->arr[d].oea-1);
		}

	    }
	}
    }
  

  /* For any domains <d> for which the env score approximation is invalid,
   * use ASC ForwardSeg() to recalculate.
   */
  for (d = 1; d <= D; d++)
    if (! (env->arr[d].flags & p7E_ENVSC_APPROX))
      {
	reanch[1].i0 = anch[d].i0;
	reanch[1].k0 = anch[d].k0;
	p7_anchor_SetSentinels(reanch, 1, L, gm->M);

	/* TODO: is <asf> allocation ok? do we need to reuse() it? */
	p7_sparse_asc_ForwardSeg(dsq, L, gm, reanch, 1, 1, sm,
				 env->arr[d].oea-1, env->arr[d].oea, env->arr[d].oeb,
				 asf, 0., -eslINFINITY, -eslINFINITY, asf->dp, asf->xmx,
				 NULL, NULL, &(env->arr[d].env_sc), NULL, NULL, NULL, NULL);

	env->arr[d].env_sc += gm->xsc[p7P_C][p7P_LOOP] * (L-env->arr[d].oeb) + gm->xsc[p7P_C][p7P_MOVE];
	// TODO: that could be -inf * 0; FIXME; ALSO ABOVE, same problem
      }


  return eslOK;

 ERROR:
  return status;
}





/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7SPARSE_ENVELOPES_EXAMPLE

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
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "include all cells in sparse mx",                   0 },
  { "-s",         eslARG_INT,    "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
 
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of sparse envelope definition";

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
  P7_OPROFILE    *om      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_CHECKPTMX   *cx      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  P7_SPARSEMX    *sxf     = NULL;
  P7_SPARSEMX    *sxd     = NULL;
  P7_TRACE       *tr      = NULL;
  P7_SPARSEMX    *asf     = NULL;
  P7_SPARSEMX    *asb     = NULL;
  P7_SPARSEMX    *asd     = NULL;
  P7_ANCHORS     *anch    = p7_anchors_Create();
  P7_ANCHORS     *vanch   = p7_anchors_Create();
  P7_ANCHORHASH  *ah      = p7_anchorhash_Create();
  P7_ENVELOPES   *env     = p7_envelopes_Create();
  float          *wrk     = NULL;
  float           fsc, vsc, asc, asc_f, asc_b;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
 
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Get a sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create (hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config   (gm, hmm, bg);
  p7_oprofile_Convert (gm, om);

  p7_bg_SetLength     (bg, sq->n);
  p7_profile_SetLength(gm, sq->n);
  p7_oprofile_ReconfigLength(om, sq->n);

  /* We need a sparse mask <sm>.
   * To get it, run checkpointed Fwd/Bck/Decoding
   */
  cx = p7_checkptmx_Create(hmm->M, sq->n, ESL_MBYTES(32));
  sm = p7_sparsemask_Create(gm->M, sq->n);
  if (esl_opt_GetBoolean(go, "-a")) 
    p7_sparsemask_AddAll(sm);
  else {
    p7_ForwardFilter (sq->dsq, sq->n, om, cx, /*fsc=*/NULL);
    p7_BackwardFilter(sq->dsq, sq->n, om, cx, sm, p7_SPARSEMASK_THRESH_DEFAULT);
  }

  /* Allocate DP matrices, traceback */
  sxf  = p7_sparsemx_Create(sm);
  sxd  = p7_sparsemx_Create(sm);
  asf  = p7_sparsemx_Create(sm);
  asb  = p7_sparsemx_Create(sm);
  asd  = p7_sparsemx_Create(sm);
  tr   = p7_trace_Create();

  /* First pass analysis */
  p7_SparseViterbi (sq->dsq, sq->n, gm, sm,  sxf, tr, &vsc);
  p7_SparseForward (sq->dsq, sq->n, gm, sm,  sxf,     &fsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, sm,  sxd,     NULL);
  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxd,     sxd);


  /* MPAS */
  p7_sparse_anchors_SetFromTrace(sxd, tr, vanch);
  p7_trace_Reuse(tr);

  p7_sparse_Anchors(rng, sq->dsq, sq->n, gm,
		    vsc, fsc, sxf, sxd, vanch,
		    tr, &wrk, ah,
		    asf, anch, &asc, 
		    NULL);

  /* Remaining ASC calculations; MPAS did <asf> */
  p7_sparse_asc_Backward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asb, &asc_b);
  p7_sparse_asc_Decoding(sq->dsq, sq->n, gm, anch->a, anch->D, asc_f, asf, asb, asd);
  
  p7_spascmx_DumpSpecials(stdout, asd, anch->a, anch->D);

  /* Envelope determination. */
  p7_sparse_Envelopes(sq->dsq, sq->n, gm, anch->a, anch->D, asf, asd, env);

  p7_envelopes_Dump(stdout, env);

  if (wrk) free(wrk);
  p7_envelopes_Destroy(env);
  p7_anchorhash_Destroy(ah);
  p7_anchors_Destroy(anch);
  p7_anchors_Destroy(vanch);
  p7_sparsemx_Destroy(asd);
  p7_sparsemx_Destroy(asb);
  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemx_Destroy(sxf);
  p7_trace_Destroy(tr);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(cx);
  p7_profile_Destroy(gm); 
  p7_oprofile_Destroy(om);
  esl_sqfile_Close(sqfp); 
  esl_sq_Destroy(sq);
  p7_hmmfile_Close(hfp);  
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_ENVELOPES_EXAMPLE*/
/*---------------------- end, example ---------------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
