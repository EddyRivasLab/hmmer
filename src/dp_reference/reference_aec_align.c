/* Reference implementation of anchor/envelope constrained (AEC) alignment;
 * by maximum expected gain (MEG).
 * 
 * All reference implementation code is for development and
 * testing. It is not used in HMMER's main executables. Production
 * code uses sparse dynamic programming.
 * 
 * Contents:
 *    1. AEC MEG alignment, fill.
 *    2. Unit tests.
 *    3. Test driver.
 *    4. Example.
 *    5. Copyright and license information
 */
#include "p7_config.h"

#include "easel.h"

#include "base/p7_envelopes.h"
#include "base/p7_profile.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_aec_align.h"
#include "dp_reference/reference_aec_trace.h"

/*****************************************************************
 * 1. AEC MEG alignment, fill.
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

/* Function:  p7_reference_AEC_Align()
 * Synopsis:  Maximum expected gain, anchor/envelope-constrained alignment.
 *
 * Purpose:   Calculates an anchor/envelope-constrained (AEC) maximum expected
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
 *            tr   : RETURN: MEG path
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_reference_AEC_Align(const P7_PROFILE *gm, P7_ENVELOPES *env, const P7_REFMX *apu, const P7_REFMX *apd, P7_REFMX *mx, P7_TRACE *tr)
{
  int    M = env->M;
  int    L = env->L;
  float *tsc;
  float *dpp;
  float *dpc;
  float *ppp;
  float *xc, *xp;
  float  mlv, dlv, xB, xE, xN;
  float  mgv, dgv;
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
  mx->type = p7R_AEC_MEG_ALIGN;

  /* Initialize rows 0..ia[0]-1.
   * Note that we don't really have to do this. We know,
   * from the envelope constraint, that the path is SNNNNN...
   * up to ia(0); only at ia(0) do we start worrying about
   * other transitions. We work it all out without shortcuts
   * because we're the reference implementation.
   */
  xN = 0.;
  for (i = 0; i < env->arr[0].ia; i++)
    {
      dpc = mx->dp[i];
      for (s = 0; s < (M+1)*p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;

      xc  = dpc;
      ppp = apd->dp[i] + (M+1)*p7R_NSCELLS;        
      xc[p7R_E]      = -eslINFINITY;
      xc[p7R_N] = xN = (i == 0 ? 0. : xN + ppp[p7R_N]);
      xc[p7R_J]      = -eslINFINITY;
      xc[p7R_B]      = -eslINFINITY;
      xc[p7R_L]      = -eslINFINITY;
      xc[p7R_G]      = -eslINFINITY;
      xc[p7R_C]      = -eslINFINITY;
      xc[p7R_JJ]     = -eslINFINITY;
      xc[p7R_CC]     = -eslINFINITY;
    }
  /* As we leave: xN retains the value for N(ia-1), which we need below. */
  xc[p7R_B] = xB = xN;
  xc[p7R_L]      = ( env->arr[0].flags & p7E_IS_GLOCAL ? -eslINFINITY : xB );
  xc[p7R_G]      = ( env->arr[0].flags & p7E_IS_GLOCAL ? xB : -eslINFINITY );


  for (d = 0; d < env->n; d++)
    {
      /*****************************************************************
       * UP sector for domain d: rows i= ia[d]..i0[d]-1
       * It's possible to have zero rows in an UP sector, when ia==i0
       */

      /* ia row is boundary case; and it may not even exist.
       * xB is already set to B(ia-1)
       */
      i   = env->arr[d].ia;             // for notation 'clarity'
      dpc = mx->dp[i];                  //  we need dpc outside block below, so transition to DOWN works
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; // now dpc is on k=1.
      if (i < env->arr[d].i0)
	{
	  tsc = gm->tsc;                    // on k=0 (w/ LM, GM entries off-by-one: for k=1)
	  ppp = apu->dp[i] + p7R_NSCELLS;   // on k=1 of UP decoding matrix, current row
	  dlv = dgv = -eslINFINITY;

	  if (env->arr[d].flags & p7E_IS_GLOCAL)
	    {
	      for (k = 1; k < env->arr[d].k0; k++)
		{
		  dpc[p7R_ML]       = -eslINFINITY;
		  dpc[p7R_MG] = mgv = ppp[p7R_ML] + ppp[p7R_MG] + P7_DELTAT(xB, tsc[p7P_GM]);  // only way into M on 1st row is an entry
		  tsc += p7P_NTRANS;  
		  dpc[p7R_IL] = dpc[p7R_IG] = -eslINFINITY;                      // no way into I
		  dpc[p7R_DL] = -eslINFINITY;
		  dpc[p7R_DG] = dgv;
		  dgv         = ESL_MAX( P7_DELTAT(mgv, tsc[p7P_MD]), P7_DELTAT(dgv, tsc[p7P_DD]));

		  dpc += p7R_NSCELLS; 
		  ppp += p7R_NSCELLS; 
		}
	    }
	  else
	    {
	      for (k = 1; k < env->arr[d].k0; k++)
		{                                         
		  dpc[p7R_ML] = mlv = ppp[p7R_ML] + ppp[p7R_MG] + P7_DELTAT(xB, tsc[p7P_LM]);
		  dpc[p7R_MG]       = -eslINFINITY;

		  tsc += p7P_NTRANS;  
		  
		  dpc[p7R_IL] = dpc[p7R_IG] = -eslINFINITY;                      // no way into I
		  dpc[p7R_DL] = dlv;
		  dpc[p7R_DG] = -eslINFINITY;
		  dlv         = ESL_MAX( P7_DELTAT(mlv, tsc[p7P_MD]), P7_DELTAT(dlv, tsc[p7P_DD]));

		  dpc += p7R_NSCELLS; 
		  ppp += p7R_NSCELLS; 
		}
	    }

	  for (k = env->arr[d].k0; k <= M; k++)   // Initialize unused remainder of UP row in standard REFMX
	    for (s = 0; s < p7R_NSCELLS; s++) 
	      *dpc++ = -eslINFINITY;

	  /* Specials, for first row ia.
	   */
	  ppp = apd->dp[i] + (M+1) * p7R_NSCELLS;  
	  xc  = dpc;
	  xc[p7R_E]      = -eslINFINITY;
	  xc[p7R_N]      = (d == 0 ? xB + ppp[p7R_N] : -eslINFINITY);
	  xc[p7R_J]      = (d == 0 ? -eslINFINITY    : xB + ppp[p7R_JJ]);
	  xc[p7R_B] = xB = (d == 0 ? xc[p7R_N] : xc[p7R_J]);
	  xc[p7R_L]      = (env->arr[d].flags & p7E_IS_GLOCAL ? -eslINFINITY : xB);
	  xc[p7R_G]      = (env->arr[d].flags & p7E_IS_GLOCAL ? xB : -eslINFINITY);
	  xc[p7R_C]      = -eslINFINITY;
	  xc[p7R_JJ]     = -eslINFINITY;
	  xc[p7R_CC]     = -eslINFINITY;
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
	  dlv = dgv = -eslINFINITY;
	  
	  if (env->arr[d].flags & p7E_IS_GLOCAL)
	    {
	      for (k = 1; k < env->arr[d].k0; k++)
		{ 
		  dpc[p7R_ML] = -eslINFINITY;
		  dpc[p7R_MG] = mgv = ppp[p7R_ML] + ppp[p7R_MG] + 
		    ESL_MAX( ESL_MAX ( P7_DELTAT(dpp[p7R_MG], tsc[p7P_MM]),
				       P7_DELTAT(dpp[p7R_IG], tsc[p7P_IM])),
			     ESL_MAX ( P7_DELTAT(dpp[p7R_DG], tsc[p7P_DM]), 
				       P7_DELTAT(         xB, tsc[p7P_GM])));

		  tsc += p7P_NTRANS;
		  dpp += p7R_NSCELLS;
	      
		  dpc[p7R_IL] = -eslINFINITY;
		  dpc[p7R_IG] = ppp[p7R_IL] + ppp[p7R_IG] +
		    ESL_MAX( P7_DELTAT(dpp[p7R_MG], tsc[p7P_MI]),
			     P7_DELTAT(dpp[p7R_IG], tsc[p7P_II]));
		
		  dpc[p7R_DL] = -eslINFINITY;
		  dpc[p7R_DG] = dgv;
		  dgv         = ESL_MAX( P7_DELTAT(mgv, tsc[p7P_MD]), P7_DELTAT(dgv, tsc[p7P_DD]));

		  dpc += p7R_NSCELLS;
		  ppp += p7R_NSCELLS;
		}
	    }
	  else 
	    {
	      for (k = 1; k < env->arr[d].k0; k++)
		{ 
		  dpc[p7R_ML] = mlv = ppp[p7R_ML] + ppp[p7R_MG] + 
		    ESL_MAX( ESL_MAX ( P7_DELTAT(dpp[p7R_ML], tsc[p7P_MM]),
				       P7_DELTAT(dpp[p7R_IL], tsc[p7P_IM])),
			     ESL_MAX ( P7_DELTAT(dpp[p7R_DL], tsc[p7P_DM]), 
				       P7_DELTAT(         xB, tsc[p7P_LM])));
		  dpc[p7R_MG] = -eslINFINITY;

		  tsc += p7P_NTRANS;
		  dpp += p7R_NSCELLS;
	      
		  dpc[p7R_IL] = ppp[p7R_IL] + ppp[p7R_IG] +
		    ESL_MAX( P7_DELTAT(dpp[p7R_ML], tsc[p7P_MI]),
			     P7_DELTAT(dpp[p7R_IL], tsc[p7P_II]));
		  dpc[p7R_IG] = -eslINFINITY;
		
		  dpc[p7R_DL] = dlv;
		  dpc[p7R_DG] = -eslINFINITY;
		  dlv = ESL_MAX( P7_DELTAT(mlv, tsc[p7P_MD]), P7_DELTAT(dlv, tsc[p7P_DD]));

		  dpc += p7R_NSCELLS;
		  ppp += p7R_NSCELLS;
		}
	    }

	  for (k = env->arr[d].k0; k <= M; k++)   // Initialize unused remainder of UP row in standard REFMX
	    for (s = 0; s < p7R_NSCELLS; s++) 
	      *dpc++ = -eslINFINITY;

	  ppp = apd->dp[i] + (M+1) * p7R_NSCELLS;  
	  xc  = dpc;
	  xc[p7R_E]      = -eslINFINITY;
	  xc[p7R_N]      = (d == 0 ? xB + ppp[p7R_N] : -eslINFINITY);
	  xc[p7R_J]      = (d == 0 ?    -eslINFINITY : xB + ppp[p7R_JJ]);
	  xc[p7R_B] = xB = (d == 0 ? xc[p7R_N] : xc[p7R_J]);
	  xc[p7R_L]      = (env->arr[d].flags & p7E_IS_GLOCAL ? -eslINFINITY : xB);
	  xc[p7R_G]      = (env->arr[d].flags & p7E_IS_GLOCAL ? xB : -eslINFINITY);
	  xc[p7R_C]      = -eslINFINITY;
	  xc[p7R_JJ]     = -eslINFINITY;
	  xc[p7R_CC]     = -eslINFINITY;
	} /* ends loop over rows i in UP sector d */
      /* as we leave: xB = L|G for row i0-1 */


      /*****************************************************************
       * Now the DOWN sector: rows i = i0[d]..ib[d]
       */
      /* now setup on the anchor cell i0,k0 */
      i   = env->arr[d].i0;
      dpp = mx->dp[i-1] + (env->arr[d].k0-1) * p7R_NSCELLS;   
      tsc = gm->tsc     + (env->arr[d].k0-1) * p7P_NTRANS;
      ppp = apd->dp[i]  +  env->arr[d].k0    * p7R_NSCELLS;   // ppp[] on k0
      dpc = mx->dp[i];
      for (s = 0; s < p7R_NSCELLS * env->arr[d].k0; s++) *dpc++ = -eslINFINITY; // dpc[] now on k0

      if (env->arr[d].flags & p7E_IS_GLOCAL)
	{
	  dpc[p7R_ML] = -eslINFINITY;
	  if (env->arr[d].ia < env->arr[d].i0)
	    dpc[p7R_MG] = mgv = ppp[p7R_ML] + ppp[p7R_MG] + 
	      ESL_MAX( ESL_MAX ( P7_DELTAT(dpp[p7R_MG], tsc[p7P_MM]),
				 P7_DELTAT(dpp[p7R_IG], tsc[p7P_IM])),
		       ESL_MAX ( P7_DELTAT(dpp[p7R_DG], tsc[p7P_DM]), 
				 P7_DELTAT(         xB, tsc[p7P_GM])));
	  else // special case, no cells in UP sector: we must enter anchor on LMk/GMk entry
	    dpc[p7R_MG] = mgv = ppp[p7R_ML] + ppp[p7R_MG] + P7_DELTAT( xB, tsc[p7P_LM]);


	  tsc        += p7P_NTRANS;
	  dpc[p7R_IL] = dpc[p7R_IG] = -eslINFINITY;
	  dpc[p7R_DL] = dpc[p7R_DG] = -eslINFINITY;

	  xE          = mgv;
	  dgv         = P7_DELTAT(mgv, tsc[p7P_MD]);
	  dpc        += p7R_NSCELLS;
	}
      else
	{
	  if (env->arr[d].ia < env->arr[d].i0)
	    dpc[p7R_ML] = mlv = ppp[p7R_ML] + ppp[p7R_MG] + 
	      ESL_MAX( ESL_MAX ( P7_DELTAT(dpp[p7R_ML], tsc[p7P_MM]),
				 P7_DELTAT(dpp[p7R_IL], tsc[p7P_IM])),
		       ESL_MAX ( P7_DELTAT(dpp[p7R_DL], tsc[p7P_DM]), 
				 P7_DELTAT(         xB, tsc[p7P_LM])));
	  else // special case, no cells in UP sector: we must enter anchor on LMk/GMk entry
	    dpc[p7R_ML] = mlv = ppp[p7R_ML] + ppp[p7R_MG] + P7_DELTAT( xB, tsc[p7P_LM]);

	  dpc[p7R_MG] = -eslINFINITY;

	  tsc        += p7P_NTRANS;

	  dpc[p7R_IL] = dpc[p7R_IG] = -eslINFINITY;
	  dpc[p7R_DL] = dpc[p7R_DG] = -eslINFINITY;

	  xE          = mlv;
	  dlv         = P7_DELTAT(mlv, tsc[p7P_MD]);
	  dpc        += p7R_NSCELLS;
	}
      /* now dpc[] is on k0+1 */

      /* Remainder of anchor row i0 in DOWN 
       * is only reachable on deletion paths from anchor.
       */
      if (env->arr[d].flags & p7E_IS_GLOCAL)
	for (k = env->arr[d].k0+1; k <= M; k++)
	  {
	    dpc[p7R_ML] = dpc[p7R_MG] = -eslINFINITY;
	    dpc[p7R_IL] = dpc[p7R_IG] = -eslINFINITY;
	    tsc        += p7P_NTRANS;
	    dpc[p7R_DL] = -eslINFINITY;
	    dpc[p7R_DG] = dgv;
	    dgv         = P7_DELTAT(dgv, tsc[p7P_DD]);
	    dpc        += p7R_NSCELLS;
	  }
      else 
	for (k = env->arr[d].k0+1; k <= M; k++)
	  {
	    dpc[p7R_ML] = dpc[p7R_MG] = -eslINFINITY;
	    dpc[p7R_IL] = dpc[p7R_IG] = -eslINFINITY;
	    tsc        += p7P_NTRANS;
	    dpc[p7R_DL] = dlv;
	    dpc[p7R_DG] = -eslINFINITY;
	    dlv         = P7_DELTAT(dlv, tsc[p7P_DD]);
	    dpc        += p7R_NSCELLS;
	  }

      /* dpc now sits on start of specials in mx */
      xc = dpc;
      xc[p7R_E]  = xE;
      xc[p7R_N]  = -eslINFINITY;
      xc[p7R_J]  = (d == env->n-1 ? -eslINFINITY : xE);  // Only reachable by E->J, on this first row of DOWN sector.
      xc[p7R_B]  = -eslINFINITY;
      xc[p7R_L]  = -eslINFINITY;
      xc[p7R_G]  = -eslINFINITY;
      xc[p7R_C]  = (d == env->n-1 ? xE : -eslINFINITY);
      xc[p7R_JJ] = -eslINFINITY;
      xc[p7R_CC] = -eslINFINITY;

      /* rest of the DOWN sector rows for domain d:
       *  i = i0..ib; k=k0..M.
       */
      for (i = env->arr[d].i0+1; i <= env->arr[d].ib; i++)
	{
	  tsc = gm->tsc      + (env->arr[d].k0-1) * p7P_NTRANS;
	  dpp = mx->dp[i-1]  + (env->arr[d].k0-1) * p7R_NSCELLS;
	  ppp = apd->dp[i]   +  env->arr[d].k0    * p7R_NSCELLS;

	  dpc = mx->dp[i];
	  for (s = 0; s < p7R_NSCELLS * env->arr[d].k0; s++) *dpc++ = -eslINFINITY; // dpc[] now on k0
	  dlv = dgv = xE = -eslINFINITY;

	  if (env->arr[d].flags & p7E_IS_GLOCAL)
	    {
	      for (k = env->arr[d].k0; k <= M; k++)
		{
		  dpc[p7R_ML] = -eslINFINITY;
		  dpc[p7R_MG] = mgv = ppp[p7R_ML] + ppp[p7R_MG] + 
		    ESL_MAX( ESL_MAX ( P7_DELTAT(dpp[p7R_MG], tsc[p7P_MM]),
				       P7_DELTAT(dpp[p7R_IG], tsc[p7P_IM])),
		                       P7_DELTAT(dpp[p7R_DG], tsc[p7P_DM])); 
		  tsc += p7P_NTRANS;
		  dpp += p7R_NSCELLS;

		  dpc[p7R_IL] = -eslINFINITY;
		  dpc[p7R_IG] = ppp[p7R_IL] + ppp[p7R_IG] +
		    ESL_MAX( P7_DELTAT(dpp[p7R_MG], tsc[p7P_MI]),
			     P7_DELTAT(dpp[p7R_IG], tsc[p7P_II]));

		  dpc[p7R_DL] = -eslINFINITY;
		  dpc[p7R_DG] = dgv;
		  dgv         = ESL_MAX( P7_DELTAT(mgv, tsc[p7P_MD]), 
					 P7_DELTAT(dgv, tsc[p7P_DD]));

		  dpc += p7R_NSCELLS;
		  ppp += p7R_NSCELLS;

		}
	      /* dpc is on specials, +1 past [M] cells, hence the hacky -p7R_NSCELLS here: */
	      xE = ESL_MAX( dpc[p7R_MG-p7R_NSCELLS], dpc[p7R_DG-p7R_NSCELLS]);
	    }
	  else
	    {
	      for (k = env->arr[d].k0; k <= M; k++)
		{
		  dpc[p7R_ML] = mlv = ppp[p7R_ML] + ppp[p7R_MG] + 
		    ESL_MAX( ESL_MAX ( P7_DELTAT(dpp[p7R_ML], tsc[p7P_MM]),
				       P7_DELTAT(dpp[p7R_IL], tsc[p7P_IM])),
			               P7_DELTAT(dpp[p7R_DL], tsc[p7P_DM])); 
		  dpc[p7R_MG] = -eslINFINITY;
		  tsc += p7P_NTRANS;
		  dpp += p7R_NSCELLS;

		  dpc[p7R_IL] = ppp[p7R_IL] + ppp[p7R_IG] +
		    ESL_MAX( P7_DELTAT(dpp[p7R_ML], tsc[p7P_MI]),
			     P7_DELTAT(dpp[p7R_IL], tsc[p7P_II]));
		  dpc[p7R_IG] = -eslINFINITY;
		
		  xE = ESL_MAX(xE, ESL_MAX(mlv, dlv));

		  dpc[p7R_DL] = dlv;
		  dpc[p7R_DG] = -eslINFINITY;
		  dlv         = ESL_MAX( P7_DELTAT(mlv, tsc[p7P_MD]), 
					 P7_DELTAT(dlv, tsc[p7P_DD]));

		  dpc += p7R_NSCELLS;
		  ppp += p7R_NSCELLS;
		}
	    }

	  /* Specials, at the end of each row. */
	  xp   = dpp + p7R_NSCELLS;
	  xc   = dpc;
	  ppp  = apd->dp[i] + (M+1) * p7R_NSCELLS;  

	  xc[p7R_E]  = xE;
	  xc[p7R_N]  = -eslINFINITY;
	  xc[p7R_J]  = (d == env->n-1 ? -eslINFINITY : ESL_MAX(xE, xp[p7R_J] + ppp[p7R_JJ]));
	  xc[p7R_B]  = -eslINFINITY;
	  xc[p7R_L]  = -eslINFINITY;
	  xc[p7R_G]  = -eslINFINITY;
	  xc[p7R_C]  = (d == env->n-1 ? ESL_MAX(xE, xp[p7R_J] + ppp[p7R_JJ]) : -eslINFINITY);
	  xc[p7R_JJ] = -eslINFINITY;
	  xc[p7R_CC] = -eslINFINITY;
	} /* end loop over i up to ib(d) */

      /* Now, interdomain rows ib[d]+1 to ia[d]-1. 
       * We know we're going to stay in either J|C on these rows.
       */
      if (d < env->n-1)
	{ // Here we know we're in J, between two domains d and d+1:
	  /* a little overwriting of row ib(d), which might transition directly to an ia(d+1) */
	  xB = xc[p7R_J];

	  /* Now, from ib(d)+1 to ia(d+1)-1 (or L), we know we stay
	   * in J|C. It's possible that ib(d)+1 == ia, in which case
	   * we won't do anything in this loop.
	   */
	  for (i = env->arr[d].ib+1; i < env->arr[d+1].ia; i++)
	    {
	      dpc = mx->dp[i];
	      for (s = 0; s < (M+1)*p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;

	      xc = dpc;
	      ppp = apd->dp[i] + (M+1)*p7R_NSCELLS;        
	      
	      xc[p7R_E]      = -eslINFINITY;
	      xc[p7R_N]      = -eslINFINITY;
	      xc[p7R_J] = xB = (i == 0 ? 0. : xB + ppp[p7R_JJ]);
	      xc[p7R_B]      = -eslINFINITY;
	      xc[p7R_L]      = -eslINFINITY;
	      xc[p7R_G]      = -eslINFINITY;
	      xc[p7R_C]      = -eslINFINITY;
	      xc[p7R_JJ]     = -eslINFINITY;
	      xc[p7R_CC]     = -eslINFINITY;
	    }
	  xc[p7R_B] = xB;
	  xc[p7R_L] = ( env->arr[d].flags & p7E_IS_GLOCAL ? -eslINFINITY : xB );
	  xc[p7R_G] = ( env->arr[d].flags & p7E_IS_GLOCAL ? xB : -eslINFINITY );
	}
      else
	{ // Here we know there's no more domains. We're in C until the end. No B.
	  xB = xc[p7R_C];  // we're going to misuse "xB" instead of declaring a new variable for this
	  for (i = env->arr[d].ib+1; i <= L; i++)
	    {
	      dpc = mx->dp[i];
	      for (s = 0; s < (M+1)*p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;

	      xc = dpc;
	      ppp = apd->dp[i] + (M+1)*p7R_NSCELLS;        
	      
	      xc[p7R_E]      = -eslINFINITY;
	      xc[p7R_N]      = -eslINFINITY;
	      xc[p7R_J]      = -eslINFINITY;
	      xc[p7R_B]      = -eslINFINITY;
	      xc[p7R_L]      = -eslINFINITY;
	      xc[p7R_G]      = -eslINFINITY;
	      xc[p7R_C] = xB = xB + ppp[p7R_CC];   
	      xc[p7R_JJ]     = -eslINFINITY;
	      xc[p7R_CC]     = -eslINFINITY;
	    }
	}

    } /* end loop over domains d */

  return p7_reference_aec_trace_MEG(gm, env, apd, mx, tr);
}
/*---------------- end, AEC/MEG alignment fill ------------------*/




/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef p7REFERENCE_AEC_ALIGN_TESTDRIVE
#include "hmmer.h"


/* "crashtestdummy" test
 * Compare randomly selected profile to sequences sampled
 * from that profile. 
 * 
 * The test is not very stringent, because we don't know the
 * true alignment. This is just a test that nothing obviously
 * horrible happens, like a crash, or obviously invalid data
 * structures.
 * 
 * Extended from reference_envelopes.c unit tests.
 *
 * We test:
 *    1. Coordinates of each envelope/alignment are coherent:
 *       1 <= oea <= ia <= alia <= i0 <= alib <= ib <= oeb <= L
 *       1 <= ka <= k0 <= kb <= M
 *       
 *    2. Envelopes do not overlap (assuming default threshold of
 *       0.5 when defining them):
 *         ia(d) > ib(d-1)  for d = 1..D-1
 *       (Outer envelopes, in contrast, can overlap.)
 *       
 *    3. envsc(d) <= asc_sc <= fwdsc.
 *    
 *    4. If D=1 (single domain) in both the generated trace
 *       and the inferred envelopes, and the domain coords in 
 *       the trace are encompassed by the outer envelope,
 *       then envsc(d) >= generated trace score.
 *
 *    5. MEG trace passes validation.
 *
 *    6. Indexed trace domains agree with envelope data.
 */
static void
utest_crashtestdummy(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N)
{
  char             msg[] = "reference_aec_align:: crash test dummy failed";
  ESL_SQ          *sq    = esl_sq_CreateDigital(abc);
  P7_BG           *bg    = p7_bg_Create(abc);
  P7_HMM          *hmm   = NULL;
  P7_PROFILE      *gm    = p7_profile_Create(M, abc);
  P7_TRACE        *gtr   = p7_trace_Create();            // generated trace
  P7_TRACE        *vtr   = p7_trace_Create();            // Viterbi trace
  P7_TRACE        *tr    = p7_trace_Create();            // MEG trace
  P7_REFMX        *rxf   = p7_refmx_Create(M, 20);       // Fwd, Vit ~~> ASC Decode UP
  P7_REFMX        *rxd   = p7_refmx_Create(M, 20);       // Bck, Decode ~~> ASC Decode DOWN
  P7_REFMX        *afu   = p7_refmx_Create(M, 20);       // ASC Fwd UP
  P7_REFMX        *afd   = p7_refmx_Create(M, 20);       // ASC Fwd DOWN
  P7_REFMX        *apu   = rxf;                          // for 'clarity' we use two names for this mx
  P7_REFMX        *apd   = rxd;                          //   ... and this one too.
  float           *wrk   = NULL;
  P7_COORDS2      *anch  = p7_coords2_Create(0,0);
  P7_COORDS2_HASH *htbl  = p7_coords2_hash_Create(0,0,0);
  P7_ENVELOPES    *env   = p7_envelopes_Create(0,0);
  float            tol   = 0.001;
  char             errbuf[eslERRBUFSIZE];
  float  gsc, fsc, asc;
  int    idx;
  int    d;
  
  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      /* Emit sequence from model, using an arbitrary length model of <M>;
       * restrict the emitted sequence length to 6M, arbitrarily, to 
       * keep it down to something reasonable.
       */
      if ( p7_profile_SetLength(gm, M) != eslOK) esl_fatal(msg);
      do {
	esl_sq_Reuse(sq);
	if (p7_ProfileEmit(rng, hmm, gm, bg, sq, gtr) != eslOK) esl_fatal(msg);
      } while (sq->n > M * 6); 
      if (p7_trace_Index   (gtr)                      != eslOK) esl_fatal(msg);
      if (p7_trace_Score   (gtr, sq->dsq, gm, &gsc)   != eslOK) esl_fatal(msg);

      /* Reset the length model to the actual length sq->n, then
       * put it through the domain postprocessing analysis pipeline
       */
      if ( p7_profile_SetLength(gm, sq->n)                          != eslOK) esl_fatal(msg);
     
      /* First pass analysis */
      if ( p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf, vtr, NULL) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd,      NULL) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd)  != eslOK) esl_fatal(msg);

      /* Anchor determination (MPAS algorithm) */
      if ( p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, vtr, &wrk, htbl,
				afu, afd, anch, &asc, NULL, NULL)  != eslOK) esl_fatal(msg);

      /* Reuse rxf,rxd as apu, apd; finish ASC analysis with Backward, Decoding */
      p7_refmx_Reuse(apu);  p7_refmx_Reuse(apd);
      if ( p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch->arr, anch->n, apu, apd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->arr, anch->n, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(msg);

      /* Envelope calculation */
      if ( p7_reference_Envelopes(sq->dsq, sq->n, gm, anch->arr, anch->n, apu, apd, afu, afd, env) != eslOK) esl_fatal(msg);

      //p7_envelopes_Dump(stdout, env);

      p7_refmx_Reuse(afu);
      if ( p7_reference_AEC_Align(gm, env, apu, apd, afu, tr) != eslOK) esl_fatal(msg);

      //p7_refmx_Dump(stdout, afu);
      //p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);

      /* Test 1. Coords of each domain are coherent */
      if (anch->n != env->n) esl_fatal(msg);
      for (d = 0; d < anch->n; d++)
	{
	  if (! (1 <= env->arr[d].oea &&
		 env->arr[d].oea  <= env->arr[d].ia    &&
		 env->arr[d].ia   <= env->arr[d].alia  &&
		 env->arr[d].alia <= env->arr[d].i0    &&
		 env->arr[d].i0   <= env->arr[d].alib  &&
		 env->arr[d].alib <= env->arr[d].ib    &&
		 env->arr[d].ib   <= env->arr[d].oeb   &&
		 env->arr[d].oeb  <= sq->n)) esl_fatal(msg);
	  if (! (1 <= env->arr[d].ka &&
		 env->arr[d].ka <= env->arr[d].k0 &&
		 env->arr[d].k0 <= env->arr[d].kb &&
		 env->arr[d].kb <= gm->M)) esl_fatal(msg);
	}

      /* Test 2. Envelopes do not overlap. */
      for (d = 1; d < anch->n; d++)
	if (! (env->arr[d].ia > env->arr[d-1].ib)) esl_fatal(msg);

      /* Test 3. envsc(d) <= asc_sc <= fwdsc */
      for (d = 0; d < anch->n; d++)
	if (! (env->arr[d].env_sc <= asc+tol && asc <= fsc+tol)) esl_fatal(msg);

      /* Test 4, only on D=1 case with generated trace's domain 
       * encompassed by the outer envelope 
       */
      if (gtr->ndom == 1 &&  anch->n   == 1 && 
	  gtr->sqfrom[0] >= env->arr[0].oea &&
	  gtr->sqto[0]   <= env->arr[0].oeb)
	if (! ( env->arr[0].env_sc >= gsc)) esl_fatal(msg);

      /* Test 5. MEG trace passes validation */
      if (p7_trace_Validate(tr,  abc, sq->dsq, errbuf)  != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
      if (p7_trace_Validate(gtr, abc, sq->dsq, errbuf)  != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
      if (p7_trace_Validate(vtr, abc, sq->dsq, errbuf)  != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);

      /* Test 6. Indexed domains agree with envelope data */
      p7_trace_Index(tr);
      if (tr->ndom != env->n) esl_fatal(msg);
      for (d = 0; d < env->n; d++)
	if (! ( tr->sqfrom[d]  == env->arr[d].alia &&
		tr->sqto[d]    == env->arr[d].alib &&
		tr->hmmfrom[d] == env->arr[d].ka   &&
		tr->hmmto[d]   == env->arr[d].kb)) esl_fatal(msg);


      p7_envelopes_Reuse(env);
      p7_coords2_Reuse(anch);
      p7_coords2_hash_Reuse(htbl);
      p7_refmx_Reuse(rxf); p7_refmx_Reuse(rxd);
      p7_refmx_Reuse(afu); p7_refmx_Reuse(afd);
      p7_trace_Reuse(gtr); p7_trace_Reuse(vtr); p7_trace_Reuse(tr);
      esl_sq_Reuse(sq);
    }
      
  if (wrk) free(wrk);
  p7_envelopes_Destroy(env);
  p7_coords2_Destroy(anch);
  p7_coords2_hash_Destroy(htbl);
  p7_refmx_Destroy(afu); p7_refmx_Destroy(afd);
  p7_refmx_Destroy(rxf); p7_refmx_Destroy(rxd);
  p7_trace_Destroy(vtr); p7_trace_Destroy(gtr); p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_sq_Destroy(sq);
}

/* "singlemulti" test.
 * 
 * Use p7_modelsample_SinglePathedASC() to create a
 * profile/sequence/anchorset comparison that has only a single
 * possible path when anchor set constrained. Now the expected true
 * path and envelope(s) are known, from that single path, so we can compare.
 * 
 * In order to guarantee only one possible path, while allowing
 * multiple domains, the profile is limited to glocal-only.
 *
 * This is an extension of a very similar unit test from
 * reference_envelopes.c.
 * 
 * We test:
 *     1. The MEG trace is identical to the generated path.
 *     1. Trace and envelopes agree on number of domains.
 *     2. For each domain, oea==ia, oeb==ib, and these coords
 *        agree with the trace.
 *     3. In the case of a single domain (D=1), the envelope
 *        score == the trace score.
 */
static void
utest_singlemulti(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N)
{
  char          msg[] = "reference_aec_align:: singlemulti unit test failed";
  P7_BG        *bg    = p7_bg_Create(abc);
  P7_HMM       *hmm   = NULL;
  P7_PROFILE   *gm    = NULL;
  ESL_DSQ      *dsq   = NULL;
  int           L;
  P7_TRACE     *gtr   = NULL;
  P7_TRACE     *tr    = p7_trace_Create();
  P7_COORD2    *anch  = NULL;
  int           D;
  P7_REFMX     *afu   = p7_refmx_Create(M, 20);
  P7_REFMX     *afd   = p7_refmx_Create(M, 20);
  P7_REFMX     *apu   = p7_refmx_Create(M, 20);
  P7_REFMX     *apd   = p7_refmx_Create(M, 20);
  P7_ENVELOPES *env   = p7_envelopes_Create(0,0);
  float         gsc;
  float         tol   = 0.001;
  int           idx;
  int           d;
  
  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_SinglePathedASC(rng, M, bg, &hmm, &gm, &dsq, &L, &gtr, &anch, &D, &gsc) != eslOK) esl_fatal(msg);

      if ( p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCBackward(dsq, L, gm, anch, D, apu, apd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(msg);

      if ( p7_reference_Envelopes (dsq, L, gm, anch, D, apu, apd, afu, afd, env)      != eslOK) esl_fatal(msg);

      p7_refmx_Reuse(afu);
      if ( p7_reference_AEC_Align( gm, env, apu, apd, afu, tr) != eslOK) esl_fatal(msg);
      if ( p7_trace_Index(tr)                                  != eslOK) esl_fatal(msg);

      //p7_refmx_Dump(stdout, afu);
      //p7_trace_DumpAnnotated(stdout, gtr, gm, dsq);
      //p7_trace_DumpAnnotated(stdout,  tr, gm, dsq);

      /* Test 1. The traces are identical. */
      if ( p7_trace_Compare(gtr, tr, 0.0) != eslOK) esl_fatal(msg);

      /* Test 2. Domain #'s agree */
      if (! (gtr->ndom == D && env->n == D)) esl_fatal(msg);
      if (! (tr->ndom  == D))                esl_fatal(msg);

      /* Test 3. Envelope coords (and outer env coords) match trace */
      for (d = 0; d < D; d++)
	{
	  if (! (env->arr[d].alia == gtr->sqfrom[d] &&
		 env->arr[d].alia ==  tr->sqfrom[d] &&
		 env->arr[d].alia == env->arr[d].ia &&
		 env->arr[d].alia == env->arr[d].oea)) esl_fatal(msg);

	  if (! (env->arr[d].alib == gtr->sqto[d]   &&
		 env->arr[d].alib ==  tr->sqto[d]   &&
		 env->arr[d].alib == env->arr[d].ib &&
		 env->arr[d].alib == env->arr[d].oeb)) esl_fatal(msg);

	  if (! (env->arr[d].ka == gtr->hmmfrom[d] &&
		 env->arr[d].ka ==  tr->hmmfrom[d])) esl_fatal(msg);

	  if (! (env->arr[d].kb == gtr->hmmto[d] &&
		 env->arr[d].kb ==  tr->hmmto[d])) esl_fatal(msg);
	}

      /* Test 4. If D == 1, envelope score == trace score. */
      if (D == 1 &&  esl_FCompare(env->arr[0].env_sc, gsc, tol) != eslOK) esl_fatal(msg);

      p7_envelopes_Reuse(env);
      p7_refmx_Reuse(afu); p7_refmx_Reuse(afd);
      p7_refmx_Reuse(apu); p7_refmx_Reuse(apd);
      p7_trace_Reuse(tr);

      free(dsq);
      free(anch);
      p7_trace_Destroy(gtr);
      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
    }
     
  p7_trace_Destroy(tr);
  p7_envelopes_Destroy(env);
  p7_refmx_Destroy(afu); p7_refmx_Destroy(afd);
  p7_refmx_Destroy(apu); p7_refmx_Destroy(apd);
  p7_bg_Destroy(bg);
}


#endif /*p7REFERENCE_AEC_ALIGN_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/


/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7REFERENCE_AEC_ALIGN_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for AEC/MEG alignment inference";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  int             M    = 20;

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_crashtestdummy(rng, M, abc, 10);  
  utest_singlemulti   (rng, M, abc, 10);

  fprintf(stderr, "#  status = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}

#endif /*p7REFERENCE_AEC_ALIGN_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/



/*****************************************************************
 * 4. Example
 *****************************************************************/
#ifdef p7REFERENCE_AEC_ALIGN_EXAMPLE
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
  p7_trace_Reuse(tr);
  p7_reference_AEC_Align(gm, env, apu, apd, afu, tr);

  //p7_refmx_Dump(stdout, afu);
  //p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);
  p7_envelopes_Dump(stdout, env);

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
#endif /*p7REFERENCE_AEC_ALIGN_EXAMPLE*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
