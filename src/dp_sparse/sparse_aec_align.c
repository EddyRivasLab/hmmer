/* Anchor/envelope constrained (AEC) alignment.
 *
 * Given envelopes and anchors for our domains in the target sequence,
 * now we determine an optimal alignment for each domain. Each of
 * these alignments is global with respect to the envelope-defined
 * subsequence ia(d)..ib(d); constrained to use only i,k supercells in
 * the sparse mask; and constrained to pass through the domain anchor
 * i0(d),k0(d).
 * 
 * The optimal alignment is
 *     \hat{\pi} = \argmax{\pi} \sum_{i=ia(d)}^{ib(d)} P(\pi_i),
 * where P(\pi_i) is the posterior probability of generating x_i with
 * state \pi_i, as obtained by sparse ASC posterior decoding.  This
 * optimization criterion has been called a "label gain function"
 * [HamadaAsai2012]. It was first introduced by [Kall2005].
 * 
 * Contents:
 *   1. 
 *   x. Footnotes.
 */
#if 0

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_anchors.h"
#include "base/p7_envelopes.h"
#include "base/p7_trace.h"

#include "dp_sparse/p7_sparsemx.h"


static int aec_fill(const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_ENVELOPES *env, int d,
		    const float **ret_ppp, float **ret_dpc, int *ret_g, float *ret_xX);

int
p7_sparse_aec_Align(const P7_PROFILE *gm, const P7_SPARSEMX *asd, 
		    P7_ENVELOPES *env, P7_SPARSEMX *aec, P7_TRACE *tr)
{
  const P7_SPARSEMASK *sm = asd->sm;
  const float *ppp = asd->dp;
  float       *dpc = NULL;
  int          M = env->M;
  int          L = env->L;
  int          g;
  int          d;
  float        xX;

  /* Contract checks / argument validation */
  ESL_DASSERT1(( gm->M == M ));
  ESL_DASSERT1(( sm->M == M && sm->L == L ));
  ESL_DASSERT1(( asd->type == p7S_ASC_DECODE ));
  ESL_DASSERT1(( aec->type == p7S_UNSET ));
  ESL_DASSERT1(( tr->N == 0 ));

  /* Reallocation if needed. */
  // TODO. 
  aec->type = p7S_AEC_ALIGN;
  dpc       = aec->dp;

  xX = 0.;
  for (g = 1, d = 1; d <= env->D; d++)
    aec_fill(gm, sm, env, d, &ppp, &dpc, &g, &xX);

  return eslOK;
}

static int
aec_fill(const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_ENVELOPES *env, int d,
	 const float **ret_ppp, float **ret_dpc, int *ret_g, float *ret_xX)
{
  const float *tsc       = gm->tsc;  // sets up TSC() macro, access to profile's transitions      
  const float *ppp       = *ret_ppp;
  float       *dpc       = *ret_dpc;
  const float *dpp       = NULL;
  float       *last_dpc  = NULL;
  int          is_glocal = (env->arr[d].flags & p7E_IS_GLOCAL);
  int          g         = *ret_g;
  float        xX        = *ret_xX;
  float dc;
  int   i,k,y,z;

  /* <ppp> maintenance (1): skip leading rows in segment, catch up to row ia(d) */
  while ( env->arr[d].i0 > sm->seg[g].ib) g++;                                           // Find seg g that anchor d is in. Must succeed for some g <= sm->S.
  for (i = sm->seg[g].ia; i < env->arr[d].ia; i++)                                       // Skip ASC rows ia(g)..ia(d)-1, which must be in UP.
    for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d].k0; z++)                       // UP row includes any supercells 1..k0(d)-1
      ppp += p7S_NSCELLS;      

  /* UP sector */
  for (i = env->arr[d].ia; i < env->arr[d].i0; i++)
    {
      /* <ppp> maintenance (2): skip leading DOWN supercells when i is ASC DOWN+UP row */
      if (env->arr[d-1].i0 >= sm->seg[g].ia) {                               // If there's an anchor above us in this seg, ASC matrix row includes DOWN.
	z = 0; while (z < sm->n[i] && sm->k[i][z] < env->arr[d-1].k0) z++;   // DOWN is k0(d-1)..M, so skip 1..k0(d-1)-1 in the sparse list
	for (; z < sm->n[i]; z++) ppp += p7S_NSCELLS;
      }

      /* UP top row ia(d) is initialization. 
       * We must enter G->Mk on this row, and only on this row.
       * There is no prev row (not directly, anyway), so no entries from M/I.
       */
      if (i == env->arr[d].ia) 
	{
	  dpp      = NULL;
	  last_dpc = dpc;
	  dc       = -eslINFINITY;

	  for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d].k0; z++, dpc += p7S_NSCELLS, ppp += p7S_NSCELLS)
	    {
	      dpc[p7S_ML] = ppp[p7S_ML] + ppp[p7S_MG] + 
	 	            (is_glocal ? ESL_MAX( P7_DELTAT(dc, TSC(p7P_DM, k-1)),
						  P7_DELTAT(xX, TSC(p7P_GM, k-1))) :
		                         ESL_MAX( P7_DELTAT(dc, TSC(p7P_DM, k-1)),
						  P7_DELTAT(xX, TSC(p7P_LM, k-1))));
	      dpc[p7S_IL] = -eslINFINITY;                                
	      dpc[p7S_DL] = dc;

	      /* Advance calculation of next Dk+1. No Mk->E transition because we're in UP. */
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1)
		dc = ESL_MAX( P7_DELTAT( dpc[p7S_ML], TSC(p7P_MD, k)),
			      P7_DELTAT( dpc[p7S_DL], TSC(p7P_DD, k)));
	      else 
		dc = -eslINFINITY;
	    }
	}
      

      /* UP's remaining rows, ia+1..i0-1 */
      else
	{
	  dpp      = last_dpc;
	  last_dpc = dpc;
	  dc       = -eslINFINITY;

	  for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d].k0; z++, dpc += p7S_NSCELLS, ppp += p7S_NSCELLS)
	    {
	      k = sm->k[i][z]; 

	      /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	      while (y < sm->n[i-1] && sm->k[i-1][y]  < k-1) { y++; dpp += p7S_NSCELLS; }
	      if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) 
		dpc[p7S_ML] = ppp[p7S_ML] + ppp[p7S_MG] + 
		              ESL_MAX ( ESL_MAX( P7_DELTAT(dpp[p7S_MG], TSC(p7P_MM, k-1)),
						 P7_DELTAT(dpp[p7S_IG], TSC(p7P_IM, k-1))),
					         P7_DELTAT(dpp[p7S_DG], TSC(p7P_DM, k-1)));
	      else
		dpc[p7S_ML] = -eslINFINITY;

	      /* Try to find cell i-1,k; then compute I(i,k) from it */
	      if (y < sm->n[i-1] && sm->k[i-1][y] < k) { y++; dpp += p7S_NSCELLS; }
	      if (y < sm->n[i-1] && sm->k[i-1][y] == k) 
		dpc[p7S_IL] = ppp[p7S_IL] + ppp[p7S_IL] +
		              ESL_MAX( P7_DELTAT(dpp[p7S_MG], TSC(p7P_MI, k)),
				       P7_DELTAT(dpp[p7S_IG], TSC(p7P_II, k)));
	      else
		dpc[p7S_IL] = -eslINFINITY;

	      /* Delayed store of Dk */
	      dpc[p7S_DL] = dc;
	      
	      /* Advance calculation of next Dk+1 */
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) 
		dc = ESL_MAX( P7_DELTAT( dpc[p7S_ML], TSC(p7P_MD, k)), P7_DELTAT( dpc[p7S_DL], TSC(p7P_DD, k)));
	      else 
		dc = -eslINFINITY;
	    }
	}
    } 
  /* end of UP sector calculation */

  /* Normal handoff to the anchor cell is relatively easy:
   * we know bottom UP corner is -1 supercell from where <dpc> is now,
   * because we just computed it.
   */

  /* Now the DOWN sector */
  for (i = env->arr[d].i0; i <= env->arr[d].ib; i++)
    {
      if (i == env->arr[d].i0) 
	{
	  dpp       = dpc - p7S_NSCELLS;
	  last_dpc  = dpc;

	  z = 0; while (sm->k[i][z] < env->arr[d-1].k0) z++;  // here we know we must reach the anchor cell k0; sm->n[i] boundary check unnecessary
	  k = env->arr[d-1].k0;

	  if (env->arr[d].i0 == env->arr[d].ia)
	    dpc[p7S_ML] = ppp[p7S_ML] + ppp[p7S_MG] + 
                          (is_glocal ? P7_DELTAT(xX, TSC(p7P_GM, k-1)) : 
                                       P7_DELTAT(xX, TSC(p7P_LM, k-1)));
	  else
	    dpc[p7S_ML] = ppp[p7S_ML] + ppp[p7S_MG] + 
	                  ESL_MAX( ESL_MAX( P7_DELTAT( dpp[p7S_MG], TSC(p7P_MM, k-1)),   
					    P7_DELTAT( dpp[p7S_IG], TSC(p7P_IM, k-1))),  
					    P7_DELTAT( dpp[p7S_DG], TSC(p7P_DM, k-1)));
	  dpc[p7S_IL] = -eslINFINITY;
	  dpc[p7S_DL] = -eslINFINITY;

	  xX = (is_glocal ? P7_DELTAT( P7_DELTAT(dpc[p7S_ML], TSC(p7P_MD, k)), TSC(p7P_DGE, k)) :
		            dpc[p7S_ML]);

	  /* Advance calculation of next Dk+1.
	   * We*/
	  if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) 
	    dc = P7_DELTAT( dpc[p7S_ML], TSC(p7P_MD, k));
	  else 
	    dc = -eslINFINITY;

	  /* Remainder of initial row i0 can only be reached by delete paths,
	   * and we've already dealt with Mk->E exit 
	   */
	  for (z = z+1; z < sm->n[i]; z++, dpc += p7S_NSCELLS, ppp += p7S_NSCELLS )
	    {
	      dpc[p7S_ML] = -eslINFINITY;
	      dpc[p7S_IL] = -eslINFINITY;
	      dpc[p7S_DL] = dc;

	      /* Advance calculation of next Dk+1 */
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) 
		dc = P7_DELTAT( dc, TSC(p7P_DD, k));
	      else 
		dc = -eslINFINITY;
	    } // end loop over sparse supercells z on initial row i0, glocal path
	} // end initialization of row i0


      else // for all other DOWN rows i, other than initialization on i0:
	{
	  dpp      = last_dpc;
	  last_dpc = dpc;
	  dc       = -eslINFINITY;

	  z = 0; while (z < sm->n[i]   && sm->k[i][z]   < env->arr[d-1].k0) z++;
	  y = 0; while (y < sm->n[i-1] && sm->k[i-1][z] < env->arr[d-1].k0) y++;
	  for (; z < sm->n[i]; z++, dpc += p7S_NSCELLS, ppp += p7S_NSCELLS )
	    {
	      k = sm->k[i][z]; 

	      /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	      while (y < sm->n[i-1] && sm->k[i-1][y]  < k-1) { y++; dpp += p7S_NSCELLS; }
	      if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) 
		dpc[p7S_ML] = ppp[p7S_ML] + ppp[p7S_MG] + 
		              ESL_MAX ( ESL_MAX( P7_DELTAT(dpp[p7S_MG], TSC(p7P_MM, k-1)),
						 P7_DELTAT(dpp[p7S_IG], TSC(p7P_IM, k-1))),
			                         P7_DELTAT(dpp[p7S_DG], TSC(p7P_DM, k-1)));
	      else
		dpc[p7S_ML] = -eslINFINITY;

	      /* Try to find cell i-1,k; then compute I(i,k) from it */
	      if (y < sm->n[i-1] && sm->k[i-1][y] < k) { y++; dpp += p7S_NSCELLS; }
	      if (y < sm->n[i-1] && sm->k[i-1][y] == k) 
		dpc[p7S_IL] = ppp[p7S_IL] + ppp[p7S_IL] +
		              ESL_MAX( P7_DELTAT(dpp[p7S_MG], TSC(p7P_MI, k)),
				       P7_DELTAT(dpp[p7S_IG], TSC(p7P_II, k)));
	      else
		dpc[p7S_IL] = -eslINFINITY;

	      /* Delayed store of Dk */
	      dpc[p7S_DL] = dc;
		  
	      xX =  (is_glocal ? P7_DELTAT( P7_DELTAT(dpc[p7S_ML], TSC(p7P_MD, k)), TSC(p7P_DGE, k)) :
		                 dpc[p7S_ML]);

	      /* Advance calculation of next Dk+1 */
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) 
		dc = ESL_MAX( P7_DELTAT( dpc[p7S_ML], TSC(p7P_MD, k)),
			      P7_DELTAT( dpc[p7S_DL], TSC(p7P_DD, k)));
	      else 
		dc = -eslINFINITY;
	    } // end recursion over supercells k in DOWN row i


	  /* <ppp> maintenance (3): skip trailing UP supercells when i is ASC DOWN/UP row */
	  if (env->arr[d].i0 <= sm->seg[g].ib)                                  // If there's an anchor below us, ASC row includes UP
	    for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d].k0; z++)      // UP is 1..k0(d)-1 
	      ppp += p7S_NSCELLS;
	}
    } // end of DOWN sector calculation


  /* <ppp> maintenance (4): If <d> is the last anchor in segment, skip trailing rows ib(d)+1..ib(g), which must be DOWN only rows */
  if (env->arr[d].i0 > sm->seg[g].ib)
    for (i = env->arr[d].ib+1; i <= sm->seg[g].ib; i++)
      {
	z = 0; while (z < sm->n[i] && sm->k[i][z] < env->arr[d-1].i0) z++;   // DOWN is k0(d-1)..M, so skip 1..k0(d-1)-1
	for (; z < sm->n[i]; z++) ppp += p7S_NSCELLS;
      }

  *ret_ppp = ppp;
  *ret_dpc = dpc;
  *ret_g   = g;
  *ret_xX  = xX;
  return eslOK;
}


/*****************************************************************
 * x. Choice selection functions for traceback
 *****************************************************************/


/* <dpc> and <ppc> point at DP supercell i,k in AEC mx, ASC PP mx, respectively.
 * z is index of that supercell in sparse mask (i.e. k = k[i][z])
 */

static inline int
select_down_m(const P7_SPARSEMASK *sm, const P7_ENVELOPES *env, int d, int i, const float **mod_dpc, int *mod_z)
{
  const float *dpp = *mod_dpc;
  int          z   = *mod_z;
  int          k   = sm->k[i][z];
  int          gstate[3] = { p7T_MG, p7T_IG, p7T_DG };  
  int          lstate[3] = { p7T_ML, p7T_IL, p7T_DL };
  float        path[3];

  while (sm->k[i][z]   >= env->arr[d].k0) { z--; dpp -= p7S_NSCELLS; }  // skip back over rest of DOWN(i) row. z>=0 check not needed, because k0>=1.
  z = sm->n[i-1];
  while (sm->k[i-1][z] >  k-1)            { z--; dpp -= p7S_NSCELLS; }  // skip back to i-1,k-1 cell on DOWN(i-1) row. no z>=0 check, because i-1,k-1 supercell must exist.
  ESL_DASSERT1(( z >= 0 && sm->k[i-1][z] == k-1 ));                     // connected supercell i-1,k-1 must exist, by construction.
  /* Now dpp is the i-1, k-1 supercell, which is k[i-1][y] in the sparse mask for row i-1  */

  path[0] = dpp[p7S_ML];
  path[1] = dpp[p7S_IL];
  path[2] = dpp[p7S_DL];
  
  *mod_dpc = dpp;
  *mod_z   = z;
  return (env->arr[d].flags & p7E_IS_GLOCAL ? gstate[esl_vec_FArgMax(path, 3)] : lstate[esl_vec_FArgMax(path, 3)]);
}

static inline int
select_up_m(const P7_SPARSEMASK *sm, const P7_ENVELOPES *env, int d, int i, const float **mod_dpc, int *mod_z)
{
  const float *dpp       = *mod_dpc;
  int          z         = *mod_z;
  int          k         =  sm->k[i][z];
  int          gstate[3] = { p7T_MG, p7T_IG, p7T_DG };  
  int          lstate[3] = { p7T_ML, p7T_IL, p7T_DL };
  float        path[3];
  
  dpp -= z * p7S_NSCELLS;                                                // skip remainder of UP(i)
  z    = sm->n[i-1]; 
  while (sm->k[i-1][z] >= env->arr[d].k0)   z--;                         // skip z back until we're in the UP sector. z>=0 check not needed, we know there's at least 1 supercell 
  while (sm->k[i-1][z] >  k-1)            { z--; dpp -= p7S_NSCELLS; }   // skip back to i-1,k-1 cell on UP(i-1) row; again z>=0 check not needed
  ESL_DASSERT1(( z >= 0 && sm->k[i-1][z] == k-1 ));

  path[0] = dpp[p7S_ML];
  path[1] = dpp[p7S_IL];
  path[2] = dpp[p7S_DL];
  
  *mod_dpc = dpp;
  *mod_z   = z;
  return (env->arr[d].flags & p7E_IS_GLOCAL ? gstate[esl_vec_FArgMax(path, 3)] : lstate[esl_vec_FArgMax(path, 3)]);
}	       
  

static inline int
select_down_i(const P7_SPARSEMASK *sm, const P7_ENVELOPES *env, int d, int i, const float **mod_dpc, int *mod_z)
{
  const float *dpp = *mod_dpc;
  int          z   = *mod_z;
  int          k   = sm->k[i][z];

  while (z >= 0 && sm->k[i][z]  >= env->arr[d].k0) { z--; dpp--; }  // skip back over rest of DOWN(i).
  z = sm->n[i-1]; while (sm->k[i-1][z] >  k)       { z--; dpp--; }  // skip back to i-1,k cell on DOWN(i-1)
  ESL_DASSERT1(( z >= 0 && sm->k[i-1][z] == k ));                   // connected supercell i-1,k must exist, and must have >=1 value >= -inf, by construction.

  *mod_dpc = dpp;
  *mod_z   = z;
  if (dpp[p7S_ML] >= dpp[p7S_IL]) { return (env->arr[d].flags & p7E_IS_GLOCAL ? p7T_MG : p7T_ML); }
  else                            { return (env->arr[d].flags & p7E_IS_GLOCAL ? p7T_IG : p7T_IL); }
}

static inline int
select_up_i(const P7_SPARSEMASK *sm, const P7_ENVELOPES *env, int d, int i, const float **mod_dpc, int *mod_z)
{
  const float *dpp = *mod_dpc;
  int          z   = *mod_z;
  int          k   = sm->k[i][z];
  
  dpp -= z * p7S_NSCELLS;                                                // skip remainder of UP(i)
  z    = sm->n[i-1]; 
  while (sm->k[i-1][z] >= env->arr[d].k0)   z--;                         // skip z back until we're in the UP sector. z>=0 check not needed, we know there's at least 1 supercell 
  while (sm->k[i-1][z] >  k)              { z--; dpp -= p7S_NSCELLS; }   // skip back to i-1,k cell on UP(i-1) row; again z>=0 check not needed
  ESL_DASSERT1(( z >= 0 && sm->k[i-1][z] == k ));

  *mod_dpc = dpp;
  *mod_z   = z;
  if (dpp[p7S_ML] >= dpp[p7S_IL]) { return (env->arr[d].flags & p7E_IS_GLOCAL ? p7T_MG : p7T_ML); }
  else                            { return (env->arr[d].flags & p7E_IS_GLOCAL ? p7T_IG : p7T_IL); }
}

static inline int
select_d(const P7_ENVELOPES *env, int d, const float **mod_dpc);
{
  const float *dpp;

  *mod_dpc = dpp = *mod_dpc - p7S_NSCELLS;
  
  if (dpp[p7S_ML] >= dpp[p7S_DL]) { return (env->arr[d].flags & p7E_IS_GLOCAL ? p7T_MG : p7T_ML); }
  else                            { return (env->arr[d].flags & p7E_IS_GLOCAL ? p7T_DG : p7T_DL); }
}



/* Mk->E only occurs on final row ib(d), in a DOWN sector.
 *    kp is the sparse cell list k[ib(d)]. (passing kp, nk terser than passing i and sm)
 *    nk is length of kp[] list
 *    k0 is k0(d) anchor for this domain.
 *   *mod_dpc points at last sparse cell on row ib(d);
 *            upon return, it's been updated to the supercell of the connected M.
 *    upon return, *ret_z is the index of the connected M in kp[].
 */
static inline int
select_e(const P7_SPARSEMASK *sm, const P7_ENVELOPES *env, int d, int i, const float **mod_dpc, int *ret_z)
{
  const float *dpp      = *mod_dpc;
  const float *best_dpp = dpp;
  int          best_z   = sm->n[i]-1;
  int          z;

  for (z = sm->n[i]-1; sm->k[i][z] >= env->arr[d].k0; z--, dpp -= p7S_NSCELLS)
    if (dpp[p7S_ML] > best_dpp[p7S_ML]) { best_dpp = dpp; best_z = z; }

  *mod_dpc = best_dpp;
  *ret_z   = best_z;
  return ( env->arr[d].flags & p7E_IS_GLOCAL ? p7T_MG : p7T_ML );
}



int
aec_trace(const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_ENVELOPES *env, P7_TRACE *tr)
{
  int sprv;
  int d;

  if ((status = p7_trace_Append(tr, p7T_T, 0, 0)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, 0, 0)) != eslOK) return status;

  for (d = env->D; d >= 1; d--)  
    {
      k = gm->M+1;

      /* Residues ib(d)+1 .. ia(d+1)-1 are assigned to C | J.
       * Sentinel at ia(D+1) = L+1, so we're clean, don't need a special case there.
       */
      for (i = env->arr[d+1].ia-1; i > env->arr[d].ib; i--)
	if ((status = p7_trace_Append(tr, (d == env->D ? p7T_C : p7T_J), 0, i)) != eslOK) return status;
      if ((status = p7_trace_Append(tr, p7T_E, 0, i)) != eslOK) return status;
      scur = p7T_E;
      // i now ib(d); scur now E; k is now 

      /* DOWN matrix: trace back until we reach the anchor, which always exists */
      while (i > env->arr[d].i0 || k > env->arr[d].k0)
	{
	  switch (scur) {
	  case p7T_ML: case p7T_MG:   sprv = select_down_m(sm, env, d, i, &dpc, &z);  i--; k--;              break;
	  case p7T_IL: case p7T_IG:   sprv = select_down_i(sm, env, d, i, &dpc, &z);  i--;                   break;
	  case p7T_DL: case p7T_DG:   sprv = select_d     (    env, d,    &dpc);           k--; z--;         break;   
	  case p7T_E:                 sprv = select_e     (sm, env, d, i, &dpc, &z);       k = sm->k[i][z];  break;  
	  default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in sparse AEC traceback");
	  }
	  
	  /* Right wing unfolding, on glocal Mk->E exit. Can only happen from DOWN. */
	  if ( scur == p7T_E && (env->arr[d].flags & p7E_IS_GLOCAL))
	    for (k2 = gm->M; k2 > k; k2--)
	      if ( (status = p7_trace_Append(tr, p7T_DG, k2, i)) != eslOK) return status;
      
	  /* Append selected state - then move there, and iterate. */
	  if ( (status = p7_trace_Append(tr, sprv, k, i)) != eslOK) return status;
	  scur = sprv;
	}

      /* Now we're on the anchor. i == i0; k == k0; scur == M. */

      /* UP matrix. Trace back until we reach Mk, ia(d).
       *   If ia(d) == i0, this block is simply skipped; there's no UP sector, B->{LG}->Mk0,i0 entry into anchor.
       */
      while (i > env->arr[d].ia  ||  (scur == p7T_DG || scur == p7T_DL))  // in other words: stop when we get to M state of row ia(d).
	{
	  switch (scur) {
	  case p7T_ML: case p7T_MG:  sprv = select_up_m(sm, env, d, i, &dpc, &z);  i--; k--;      break;
	  case p7T_IL: case p7T_IG:  sprv = select_up_i(sm, env, d, i, &dpc, &z);  i--;           break;
	  case p7T_DL: case p7T_DG:  sprv = select_d   (    env, d,    &dpc);           k--; z--; break;
	  default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in sparse AEC traceback");
	  }
	  
	  /* Append selected state - then move there, and iterate. */
	  if ( (status = p7_trace_Append(tr, sprv, k, i)) != eslOK) return status;
	  scur = sprv;
	}


      sprv = (env->arr[d].flags & p7E_IS_GLOCAL ? p7T_G : p7T_L );
      if ( sprv == p7T_G )  // Left wing unfolding, on glocal G->Mk entry.
	for (; k >= 1; k--)
	  if ( (status = p7_trace_Append(tr, p7T_DG, k, i)) != eslOK) return status;
      if ( (status = p7_trace_Append(tr, sprv,                     0, i)) != eslOK) return status;
      if ( (status = p7_trace_Append(tr, p7T_B,                    0, i)) != eslOK) return status;
      if ( (status = p7_trace_Append(tr, (d == 1 ? p7T_N : p7T_J), 0, i)) != eslOK) return status;
      // i is now ia(d), and we're on B.
    }

  for (i = env->arr[1].ia-1; i >= 1; i--)
    if ((status = p7_trace_Append(tr, p7T_N, 0, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_S, 0, 0))   != eslOK) return status;
  
  tr->M = sm->M;
  tr->L = sm->L;
  return p7_trace_Reverse(tr);
}



/*****************************************************************
 * x. Exegesis.
 ***************************************************************** 
 * 

 */


/*****************************************************************
 * x. Exegesis, footnotes
 ***************************************************************** 
 * 
 * [1] HACKY MISUSE OF P7_SPARSEMX
 * 
 * Should go back through and design data structure that works for all
 * three styles of sparse calculation: standard, ASC, and AEC. Using
 * current sparsemx for AEC is wasteful. In AEC we only need three
 * states M,I,D and no specials. (Perhaps less someday; see [2].)  We
 * only use the L cells (even for glocal paths) and don't even waste
 * time initializing the G cells.
 * 
 * 
 * 
 * [2] ONE-STATE DP SOMEDAY
 * 
 * In principle, AEC only needs one DP cell per i,k. The only reason
 * we need 3 states is that we're Plan7, thus prohibiting D->I/I->D
 * transition by construction. If we someday move back to Plan9,
 * revisit.
 * 
 * 
 * 
 * [3] DOMAINS CAN BE ADJACENT
 * 
 * We know there must be at least one unused row between any two
 * sparse segments (ia(g) - ib(g-1) > 1). Sparse DP calculations
 * depend on this when they initialize an ia(g)-1 row in an UP sector.
 * In AEC, however, two domains can be adjacent: ia(d) - ib(d-1) >= 1.
 * UP sectors have to be initialized on their top row ia(d).
 * 
 * 
 * 
 * [4] CUMULATIVE OR INDIVIDUAL DOMAIN MEG SCORES?
 * 
 * Right now, xX is accumulating the total MEG score for target
 * sequence.  It would be easy to get the MEG score for each domain,
 * if we want, especially since we have the DP calculation broken into
 * individual domains.
 * 
 * 
 * 
 * [5] BEAUTIFICATION CAMPAIGN SOMEDAY
 * 
 * All this stuff with pointers, <dpc> and such, is ugly. Stepping by
 * non-unit values, p7S_NSCELLS and such, is ugly. Boundary
 * checking on UP and DOWN sparse row traversal: ugly.
 * 
 * May be cleaner and more efficient if DP matrices were arrays of
 * structures, dp[c].s, where <c> counts over sparse supercells and a
 * supercell is a fixed-size structure w/ elements <s>.
 * 
 * Might also create sparse ASC traversal aid w/ sm->nu[i] (number of
 * cells on UP row) and sm->ku[i] (pointing to the first cell), and
 * analogous for sm->nd[i], sm->kd[i]. 
 * 
 * SPARSEMASK could contain ptrs to anchors and envelopes, normally
 * NULL. If anchors is non-NULL, it's a sparse ASC matrix, w/
 * traversal aids set for ASC DP. If envelopes non-NULL, it's a sparse
 * AEC matrix.
 *
 */

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

#endif /*0*/
