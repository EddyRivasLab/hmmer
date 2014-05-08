/* Reference implementation of anchor set constrained (ASC) alignment;
 * by maximum expected gain (MEG).
 * 
 * All reference implementation code is for development and
 * testing. It is not used in HMMER's main executables. Production
 * code uses sparse dynamic programming.
 * 
 * Contents:
 *    1. ASC MEG alignment, fill.
 *    2. ASC MEG alignment, trace.
 *    x. Copyright and license information
 */
#include "p7_config.h"

#include "easel.h"

#include "base/p7_envelopes.h"
#include "base/p7_profile.h"

#include "dp_reference/p7_refmx.h"

/*****************************************************************
 * 1. ASC MEG alignment, fill.
 *****************************************************************/

/* If we were using Plan9, we could do this as a "1-state" model,
 * 1 DP cell per i,k. But since we need to prohibit D->I,I->D, we
 * can't.
 * 
 * Because we fit the calculation into one matrix, we can't 
 * initialize an ia-1 top row for the up matrices, because 
 * we may be adjacent to an ib row of an immediately preceding
 * domain.
 */

/* Function:  p7_reference_ASC_Align()
 * Synopsis:  Maximum expected gain, anchor-set-constrained alignment.
 *
 * Purpose:   Calculates an anchor-set-constrained (ASC) maximum expected
 *            gain (MEG) alignment, given profile <gm>, envelope data <env>,
 *            and ASC posterior decoding matrix UP/DOWN pair <apu>/<apd>.
 *
 *            Caller provides <mx> as space for the DP calculation; this
 *            can be of any allocated size, and will be reallocated if needed.
 *
 * Args:      gm   : profile
 *            env  : envelopes (env->n of them)
 *            apu  : ASC Decoding UP matrix
 *            apd  : ASC Decoding DOWN matrix
 *            mx   : MEG alignment matrix space
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_reference_ASC_Align(const P7_PROFILE *gm, const P7_ENVELOPES *env, const P7_REFMX *apu, const P7_REFMX *apd, P7_REFMX *mx)
{
  int    M = env->M;
  int    L = env->L;
  float *tsc;
  float *dpp;
  float *dpc;
  float *ppp;
  float *xc, *xp;
  float  mlv, dlv, xB, xE;
  int    d,i,k,s;
  int    status;

  /* Contract checks / argument validation */
  ESL_DASSERT1(( gm->M == M ));
  ESL_DASSERT1(( apu->M == M && apu->L == L && apu->type == p7R_ASC_DECODE_UP));
  ESL_DASSERT1(( apd->M == M && apd->L == L && apd->type == p7R_ASC_DECODE_DOWN));
  
  /* Reallocation if needed */
  if ( (status = p7_refmx_GrowTo(mx, M, L)) != eslOK) return status;
  mx->M    = M;
  mx->L    = L;
  mx->type = p7R_ASC_ALIGN;
#if eslDEBUGLEVEL >= 1
  p7_refmx_Zero(mx, M, L);
#endif


  for (d = 0; d < env->n; d++)
    {
      /*****************************************************************
       * UP sector for domain d: 
       *   i= ia..i0-1; k = 1..k0-1.
       *   if ia==i0, or k0==1, there are 0 cells in UP sector.
       */

      /* ia row is boundary case; and it may not even exist */
      xB  = 0.;                         // as if it were ia-1
      i   = env->arr[d].ia;             // for notation 'clarity'
      dpc = mx->dp[i];                  //  we need dpc outside block below, so transition to DOWN works
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; // now dpc is on k=1.
      if (i < env->arr[d].i0)
	{
	  /* more setup */
	  tsc = gm->tsc;                    // on k=0 (w/ LM, GM entries off-by-one: for k=1)
	  ppp = apu->dp[i] + p7R_NSCELLS;   // on k=1 of UP decoding matrix, current row
	  dlv = -eslINFINITY;

	  for (k = 1; k < env->arr[d].k0; k++)
	    {                                         
	      dpc[p7R_ML] = mlv = ppp[p7R_ML] + ppp[p7R_MG]
		+ P7_DELTAT(xB, ESL_MAX(tsc[p7P_LM], tsc[p7P_GM]));  // only way into M on 1st row is an entry

	      tsc += p7P_NTRANS;  

	      dpc[p7R_IL]       = -eslINFINITY;                      // no way into I
	      dpc[p7R_DL]       = dlv;
	      dlv               = ESL_MAX( P7_DELTAT(mlv, tsc[p7P_MD]), P7_DELTAT(dlv, tsc[p7P_DD]));

	      dpc += p7R_NSCELLS; 
	      ppp += p7R_NSCELLS; 
	    }

	  /* update xB */
	  ppp = apd->dp[i] + (M+1) * p7R_NSCELLS;  
	  xc  = mx->dp[i]  + (M+1) * p7R_NSCELLS;  
	  xc[p7R_B] = xB = (d == 0 ? xB + ppp[p7R_N] : xB + ppp[p7R_JJ]);	  
	}

      /* remaining rows of UP sector (if any) */
      for (i = env->arr[d].ia + 1; i < env->arr[d].i0; i++)
	{
	  /* setup, as above */
	  tsc = gm->tsc;
	  ppp = apu->dp[i] + p7R_NSCELLS;
	  dpp = mx->dp[i-1];               // k=0. No valid values here. We will ++ before access.
	  dpc = mx->dp[i];
	  for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; // now dpc is on k=1.
	  dlv = -eslINFINITY;

	  for (k = 1; k < env->arr[d].k0; k++)
	    { 
	      mlv = dpc[p7R_ML] = ppp[p7R_ML] + ppp[p7R_MG] + 
		ESL_MAX( ESL_MAX ( P7_DELTAT(dpp[p7R_ML], tsc[p7P_MM]),
				   P7_DELTAT(dpp[p7R_IL], tsc[p7P_IM])),
			 ESL_MAX ( P7_DELTAT(dpp[p7R_DL], tsc[p7P_DM]), 
				   P7_DELTAT( xB, ESL_MAX(tsc[p7P_LM], tsc[p7P_GM]))));

	      tsc += p7P_NTRANS;
	      dpp += p7R_NSCELLS;
	      
	      dpc[p7R_IL] = ppp[p7R_IL] + ppp[p7R_IG] +
		ESL_MAX( P7_DELTAT(dpp[p7R_ML], tsc[p7P_MI]),
			 P7_DELTAT(dpp[p7R_IL], tsc[p7P_II]));
		
	      dpc[p7R_DL] = dlv;
	      dlv = ESL_MAX( P7_DELTAT(mlv, tsc[p7P_MD]), P7_DELTAT(dlv, tsc[p7P_DD]));

	      dpc += p7R_NSCELLS;
	      ppp += p7R_NSCELLS;
	    }

	  ppp = apd->dp[i] + (M+1) * p7R_NSCELLS;  
	  xc  = mx->dp[i]  + (M+1) * p7R_NSCELLS;  
	  xc[p7R_B] = xB = (d == 0 ? xB + ppp[p7R_N] : xB + ppp[p7R_JJ]);	  
	} /* ends loop over rows i in UP sector d */


      /* The very last cell we calculated was the cell diagonal
       * from the anchor. We need that cell for initializing DOWN.
       * <dpc> just stepped past it; step it back. Because we
       * were careful to always initialize dpc even when UP sector
       * has no valid cells, this always works (so long as we don't
       * dereference it, because it still might point to
       * invalid data).
       */
      dpp = dpc - p7R_NSCELLS;

      /* now setup on the anchor cell i0,k0 */
      i   = env->arr[d].i0;
      tsc = gm->tsc     + (env->arr[d].k0-1) * p7R_NSCELLS;   
      dpc = mx->dp[i]   + (env->arr[d].k0-1) * p7R_NSCELLS; 
      ppp = apd->dp[i]  +  env->arr[d].k0    * p7R_NSCELLS;   // ppp[] on k0
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; // now dpc is on k=k0.

      if (env->arr[d].ia < env->arr[d].i0)
	dpc[p7R_ML] = mlv = ppp[p7R_ML] + ppp[p7R_MG] + 
	  ESL_MAX( ESL_MAX ( P7_DELTAT(dpp[p7R_ML], tsc[p7P_MM]),
			     P7_DELTAT(dpp[p7R_IL], tsc[p7P_IM])),
		   ESL_MAX ( P7_DELTAT(dpp[p7R_DL], tsc[p7P_DM]), 
			     P7_DELTAT( xB, ESL_MAX(tsc[p7P_LM], tsc[p7P_GM]))));
      else // special case, no cells in UP sector: we must enter anchor on LMk/GMk entry
	dpc[p7R_ML] = mlv = ppp[p7R_ML] + ppp[p7R_MG] + P7_DELTAT( xB, ESL_MAX(tsc[p7P_LM], tsc[p7P_GM]));

      dpc[p7R_IL] = -eslINFINITY;
      dpc[p7R_DL] = -eslINFINITY;
      tsc        += p7P_NTRANS;
      xE          = mlv;
      dlv         = P7_DELTAT(mlv, tsc[p7P_MD]);
      dpc        += p7R_NSCELLS;

      /* Remainder of anchor row k0 in DOWN 
       * is only reachable on deletion paths from anchor.
       */
      for (k = env->arr[d].k0+1; k <= M; k++)
	{
	  dpc[p7R_ML] = dpc[p7R_IL] = -eslINFINITY;
	  tsc        += p7P_NTRANS;
	  dpc[p7R_DL] = dlv;
	  dlv         = P7_DELTAT(mlv, tsc[p7P_MD]);
	  dpc        += p7R_NSCELLS;
	}

      /* dpc now sits on start of specials in mx, which is where
       * we want to store xE; we can get away with dpc[E] here
       */
      dpc[p7R_E] = xE;

      /* rest of the DOWN sector for domain d:
       *  i = i0..ib; k=k0..M.
       */
      for (i = env->arr[d].i0+1; i <= env->arr[d].ib; i++)
	{
	  tsc = gm->tsc      + (env->arr[d].k0-1) * p7P_NTRANS;
	  dpp = mx->dp[i-1]  + (env->arr[d].k0-1) * p7R_NSCELLS;
	  dpc = mx->dp[i]    + (env->arr[d].k0-1) * p7R_NSCELLS;
	  ppp = apd->dp[i]   +  env->arr[d].k0    * p7R_NSCELLS;
	  for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;
	  dlv = xE = -eslINFINITY;

	  for (k = env->arr[d].k0; k <= M; k++)
	    {
	      dpc[p7R_ML] = mlv = ppp[p7R_ML] + ppp[p7R_MG] + 
		ESL_MAX( ESL_MAX ( P7_DELTAT(dpp[p7R_ML], tsc[p7P_MM]),
				   P7_DELTAT(dpp[p7R_IL], tsc[p7P_IM])),
		                   P7_DELTAT(dpp[p7R_DL], tsc[p7P_DM])); 
	      tsc += p7P_NTRANS;
	      dpp += p7R_NSCELLS;

	      dpc[p7R_IL] = ppp[p7R_IL] + ppp[p7R_IG] +
		ESL_MAX( P7_DELTAT(dpp[p7R_ML], tsc[p7P_MI]),
			 P7_DELTAT(dpp[p7R_IL], tsc[p7P_II]));
		
	      xE = ESL_MAX(xE, ESL_MAX(mlv, dlv));

	      dpc[p7R_DL] = dlv;
	      dlv         = ESL_MAX( P7_DELTAT(mlv, tsc[p7P_MD]), 
				     P7_DELTAT(dlv, tsc[p7P_DD]));

	      dpc += p7R_NSCELLS;
	      ppp += p7R_NSCELLS;
	    }

	  xp   = dpp + p7R_NSCELLS;
	  xc   = dpc;
	  ppp  = apd->dp[i] + (M+1) * p7R_NSCELLS;  
	  xc[p7R_E] = (d == env->n-1 ? 
		       ESL_MAX(xE, xp[p7R_E] + ppp[p7R_CC]) :
		       ESL_MAX(xE, xp[p7R_E] + ppp[p7R_JJ]));
	}
    } /* end loop over domains d */

  return eslOK;
}


/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7REFERENCE_ASC_ALIGN_EXAMPLE
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
static char banner[] = "example of ASC MEG alignment step";

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
  P7_COORDS2     *anch    = p7_coords2_Create(0,0);
  P7_REFMX       *rxf     = NULL;
  P7_REFMX       *rxd     = NULL;
  P7_REFMX       *afu     = NULL;
  P7_REFMX       *afd     = NULL;
  P7_REFMX       *apu     = NULL;
  P7_REFMX       *apd     = NULL;
  P7_TRACE       *tr      = NULL;
  float          *wrk     = NULL;
  P7_ENVELOPES   *env     = p7_envelopes_Create(0,0);
  P7_COORDS2_HASH *hashtbl = p7_coords2_hash_Create(0,0,0);
  float           fsc, vsc, asc, asc_b;
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
  afu = p7_refmx_Create(gm->M, sq->n);
  afd = p7_refmx_Create(gm->M, sq->n);

  /* First pass analysis */
  p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf,  tr, &vsc);
  p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc);
  p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd, NULL);   
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd);   

  /* MPAS algorithm gets us an anchor set */
  p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, tr, &wrk, hashtbl,
		       afu, afd, anch, &asc, NULL, NULL);

  
  /* We no longer need rxf and rxd. 
   * Use their space for apu/apd pair, which will briefly
   * hold ASC Backward matrices, then get used for ASC Decoding.
   */
  apu = rxf; p7_refmx_Reuse(apu);
  apd = rxd; p7_refmx_Reuse(apd);

  p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch->arr, anch->n, apu, apd, &asc_b);
  
  /* ASC Decoding takes afu/afd and abu/abd as input;
   * overwrites abu/abd with decoding matrices
   */
  p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->arr, anch->n, afu, afd, apu, apd, apu, apd);

  /* Envelope calculation needs to get four matrices:
   * ASC Decoding pair, apu/apd, and it will leave these constant;
   * ASC Forward pair,  afu/afd, and it will overwrite these.
   */
  p7_reference_Envelopes(sq->dsq, sq->n, gm, anch->arr, anch->n, apu, apd, afu, afd, env);

  p7_envelopes_Dump(stdout, env);

  /* MEG alignment step uses afu as its matrix; apu/apd decoding matrices as its input */
  p7_refmx_Reuse(afu);
  p7_reference_ASC_Align(gm, env, apu, apd, afu);

  p7_refmx_Dump(stdout, afu);

  p7_envelopes_Destroy(env);
  p7_coords2_hash_Destroy(hashtbl);
  if (wrk) free(wrk);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(rxd);
  p7_refmx_Destroy(rxf);
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
#endif /*p7REFERENCE_ASC_ALIGN_EXAMPLE*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
