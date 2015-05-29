/* Anchor/envelope constrained (AEC) alignment.
 *
 * Given envelopes and anchors for our domains in the target sequence,
 * now we determine an optimal alignment for each domain.
 *
 * Each alignment is global with respect to the envelope-defined
 * subsequence ia(d)..ib(d); constrained to use only i,k supercells in
 * the sparse mask; and constrained to pass through the domain anchor
 * i0(d),k0(d).
 * 
 * The optimality criterion is "maximum expected accuracy".  We
 * maximize the expected accuracy of the assignment of profile states
 * to each aligned residue; i.e.
 *     \hat{\pi} = \argmax{\pi} \sum_{i=ia(d)}^{ib(d)} P(\pi_i),
 * where P(\pi_i) is the posterior probability of generating x_i with
 * state \pi_i, as obtained by sparse ASC posterior decoding.  This
 * criterion has been called a "label gain function"
 * [HamadaAsai2012]. It was first introduced by [Kall2005].
 * 
 * Contents:
 *   1. API wrapper, p7_sparse_aec_Align().
 *   2. DP recursion, aec_fill().
 *   3. Choice selection functions for the traceback.
 *   4. Traceback, aec_trace().
 *   5. Footnotes.
 *   6. Unit tests.
 *   7. Test driver.
 *   8. Example.
 *   9. Copyright and licence information.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"


#include "base/p7_profile.h"
#include "base/p7_envelopes.h"
#include "base/p7_trace.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spaecmx.h"

static int aec_fill (const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_ENVELOPES *env, int d,   const float **mod_ppp, float **mod_dpc, int *mod_g, float *mod_xX);
static int aec_trace(const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_ENVELOPES *env,  const float *dpc, P7_TRACE *tr);


/*****************************************************************
 * 1. API wrapper, p7_sparse_aec_Align()
 *****************************************************************/

/* Function:  p7_sparse_aec_Align()
 * Synopsis:  Anchor/envelope constrained optimal alignment of profile to sequence.
 *
 * Purpose:   Given a profile <gm>, the sparse posterior decoding matrix
 *            <asd> that we got from aligning <gm> to a target
 *            sequence, and the envelope set <env> that defines
 *            domains and anchors that we've called on the target
 *            sequence; we perform a maximum expected accuracy
 *            anchor/envelope-constrained alignment to each domain.
 *            
 *            This requires a sparse AEC DP matrix, <aec>. Caller
 *            provides us with this empty object. It will be
 *            reallocated if needed. Upon return, the <aec> object
 *            contains the filled sparse AEC DP matrix.
 *            
 *            The optimal alignment trace for the complete sequence is
 *            stored in <tr>. Caller provides us an empty <tr> object, 
 *            and it will be reallocated as needed too.
 *            
 * Args:      gm  : profile
 *            asd : sparse posterior decoding matrix for gm/seq comparison
 *            env : envelope set determined for seq
 *            aec : caller-provided sparse DP matrix for AEC alignment fill
 *            tr  : RETURN: optimal alignment
 *
 * Returns:   <eslOK> on success. <tr> contains the optimal alignment path
 *            for the complete sequence. <aec> contains the sparse
 *            AEC DP matrix.
 *
 * Throws:    <eslEMEM> on allocation failure. State of <tr>, <aec> 
 *            undefined.
 */
int
p7_sparse_aec_Align(const P7_PROFILE *gm, const P7_SPARSEMX *asd, 
		    P7_ENVELOPES *env, P7_SPARSEMX *aec, P7_TRACE *tr)
{
  const P7_SPARSEMASK *sm  = asd->sm;   
  const float         *ppp = asd->dp;
  float               *dpc = NULL;
  int    g;
  int    d;
  float  xX;                            

  /* Contract checks / argument validation */
  ESL_DASSERT1(( gm->M == env->M ));
  ESL_DASSERT1(( sm->M == env->M && sm->L == env->L ));
  ESL_DASSERT1(( asd->type == p7S_ASC_DECODE ));
  ESL_DASSERT1(( aec->type == p7S_UNSET ));
  ESL_DASSERT1(( tr->N == 0 ));

  /* Reallocation if needed. TODO */
  aec->type = p7S_AEC_ALIGN;
  aec->sm   = sm;               // This should be in p7_spaecmx_Resize(), when we have it.
  dpc       = aec->dp;
  

  /* AEC matrix only needs main state supercells (no specials), and it
   * only needs them for ia(d)..i0/k0(d)..ib(d) for each domain d.
   */
  xX = 0.;
  g  = 0;
  for (d = 1; d <= env->D; d++)
    aec_fill(gm, sm, env, d, &ppp, &dpc, &g, &xX);

  //p7_spaecmx_Dump(stdout, aec, env);  

  aec_trace(gm, sm, env, dpc, tr);     // dpc is +1 from last supercell. aec_trace knows to back it up -1 supercell, because it does that for each successive domain anyway
  return eslOK;
}




/*****************************************************************
 * 2. DP recursion, aec_fill()
 *****************************************************************/

/* aec_fill()
 * Fill stage of AEC DP, one domain <d> at a time.
 * 
 * Maintaining state in posterior decoding matrix <asd> and the AEC DP
 * matrix <aec> is detail-oriented, so listen up. We only need main
 * states, no specials, so this routine only deals with ptrs into
 * <asd> and <aec> main state memory: <ppp> and <dpc>, respectively.
 * 
 * <d>   : Caller makes D calls to aec_fill(), d=1..D.
 * <ppp> : If <d> is in the same segment as the previous one, then <ppp>
 *         is on next ASC row, ib(d-1)+1. If <d> is in a new segment,
 *         <ppp> is at start of first row of that segment g, ia(g).
 *         Starts at asd->dp.
 * <dpc> : Start of DP supercells for domain <d> in AEC matrix. Starts at
 *         aec->dp.
 * <xX>  : Running total posterior prob for all homologous rows in alignment;
 *         (nonhomologous rows, outside domains, are not counted).
 *         Starts at 0.
 *
 * At end, when all domains have been looped over, <dpc> is +1 from
 * the final supercell of the AEC matrix.  <dpc> - <aec->dp> is then
 * the number of supercells.
 */
static int
aec_fill(const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_ENVELOPES *env, int d, 
	 const float **mod_ppp, float **mod_dpc, int *mod_g, float *mod_xX)
{
  const float *tsc       = gm->tsc;    // sets up TSC() macro, access to profile's transitions      
  const float *ppp       = *mod_ppp;   
  float       *dpc       = *mod_dpc;
  const float *dpp       = NULL;
  float       *last_dpc  = NULL;
  int          is_glocal = (env->arr[d].flags & p7E_IS_GLOCAL);
  int          i         = env->arr[d-1].ib+1;                  // for case where <d> is in current <g>; if not, initialization below will reset i. For d=1, d-1 sentinel gives i=1.
  int          g         = *mod_g;
  float        xX        = *mod_xX;
  float        dc;
  int          k,y,z;

  /* <ppp> maintenance (1): skip leading rows <i>, catch <i> and <ppp> up to row ia(d), where <dpc> already is.
   *
   * Two cases: 
   *   1. <d> is first domain in a new segment <g>. 
   *        In this case, we advance <g> to the correct segment, and set i=ia(g).
   *        In the ASC matrix <ppp>, we skip rows ia(g)..ia(d)-1; we know these are only UP sector rows.
   *   2. <d> isn't the first domain in <g>, we're still working on a <g>.
   *        Now we leave <g> as it is and we use i=ib(d-1)+1, which was already set above
   *        In the ASC matrix <ppp>, we skip rows ib(d-1)+1..ia(d)-1.
   *        We know caller has given us i == ib(d-1)+1.
   *        We know each row has both a DOWN and an UP sector that <ppp> must be skipped past.
   * 
   * When this maintenance is complete:
   *    i = ia(d)
   *    <dpc> is on the AEC supercell for first sparse supercell in UP(i) 
   *    <ppp> is on first ASC supercell on row i (which still might be in a DOWN sector and need more skipping)
   */
  while ( env->arr[d].i0 > sm->seg[g].ib)                                           // Find seg that anchor d is in. Must succeed for some g <= sm->S. 
    { g++; i = sm->seg[g].ia; }                                                     // <ppp> has no storage in anchorless segments.        

  if (env->arr[d-1].i0 < sm->seg[g].ia)                                             // If <d> is 1st domain in seg <g>: 
    {                                                                               // then for ASC rows ia(g)..ia(d)-1, ASC is only in UP sector; <ppp> advances over UP(i) rows.
      for ( ; i < env->arr[d].ia; i++)      
	for (z = 0;  z < sm->n[i] && sm->k[i][z] < env->arr[d].k0;  z++)            // UP row includes any supercells 1..k0(d)-1
	  ppp += p7S_NSCELLS;  
    }
  else                                                                              // If <d> isn't 1st domain in segment <g>:
    {                                                                               // then for ASC rows ib(d-1)+1..ia(d)-1, must advance ppp past DOWN(i) and UP(i)
      for ( ; i < env->arr[d].ia; i++) {
	for (z = sm->n[i]-1; z >= 0       && sm->k[i][z] >= env->arr[d-1].k0; z--)  // Order of access doesn't matter here; only # of sparse supercells in DOWN(i)
	  ppp += p7S_NSCELLS;  
	for (z = 0;          z < sm->n[i] && sm->k[i][z] <  env->arr[d].k0;   z++)  // UP row includes any supercells 1..k0(d)-1
	  ppp += p7S_NSCELLS; 
      }
    }



  /* Recursion starts with the UP sector for domain <d>.
   *    This runs from ia(d)..i0(d)-1,  1..k0(d)-1.
   */
  for (i = env->arr[d].ia; i < env->arr[d].i0; i++)
    {
      /* <ppp> maintenance (2): skip leading DOWN supercells when i is ASC DOWN+UP row */
      if (env->arr[d-1].i0 >= sm->seg[g].ia)                                     // If there's an anchor above us in this seg, ASC matrix row includes DOWN.
	for (z = sm->n[i]-1; z >= 0 && sm->k[i][z] >= env->arr[d-1].k0; z--)     // DOWN is k0(d-1)..M, and the order we skip them doesn't matter, only the order.
	  ppp += p7S_NSCELLS;    

      /* UP top row ia(d) is initialization. 
       * We must enter G->Mk on this row, and *only* on this row.
       * There are no entries from M/I on any prev row, but we can do M->D along the row.
       */
      if (i == env->arr[d].ia) 
	{
	  dpp      = NULL;
	  last_dpc = dpc;
	  dc       = -eslINFINITY;

	  for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d].k0; z++, dpc += p7S_NSCELLS, ppp += p7S_NSCELLS)
	    {
	      k = sm->k[i][z];
	      dpc[p7S_ML] = ppp[p7S_ML] + ppp[p7S_MG] + ( is_glocal ?  P7_DELTAT(xX, TSC(p7P_GM, k-1)) : P7_DELTAT(xX, TSC(p7P_LM, k-1)));
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

	  y = 0;
	  for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d].k0; z++, dpc += p7S_NSCELLS, ppp += p7S_NSCELLS)
	    {
	      k = sm->k[i][z]; 

	      /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	      while (y < sm->n[i-1] && sm->k[i-1][y]  < k-1) { y++; dpp += p7S_NSCELLS; }
	      if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) 
		dpc[p7S_ML] = ppp[p7S_ML] + ppp[p7S_MG] + 
		              ESL_MAX ( ESL_MAX( P7_DELTAT(dpp[p7S_ML], TSC(p7P_MM, k-1)),
						 P7_DELTAT(dpp[p7S_IL], TSC(p7P_IM, k-1))),
					         P7_DELTAT(dpp[p7S_DL], TSC(p7P_DM, k-1)));
	      else
		dpc[p7S_ML] = -eslINFINITY;

	      /* Try to find cell i-1,k; then compute I(i,k) from it */
	      if (y < sm->n[i-1] && sm->k[i-1][y] < k) { y++; dpp += p7S_NSCELLS; }
	      if (y < sm->n[i-1] && sm->k[i-1][y] == k) 
		dpc[p7S_IL] = ppp[p7S_IL] + ppp[p7S_IG] +
		              ESL_MAX( P7_DELTAT(dpp[p7S_ML], TSC(p7P_MI, k)),
				       P7_DELTAT(dpp[p7S_IL], TSC(p7P_II, k)));
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

	  z = 0; while (sm->k[i][z] < env->arr[d].k0) z++;  // here we know we must reach the anchor cell k0; sm->n[i] boundary check unnecessary
	  k = env->arr[d].k0;

	  if (env->arr[d].i0 == env->arr[d].ia)                                          // This handles the special case of no UP sector, in which case we MUST enter at the anchor.
	    dpc[p7S_ML] = ppp[p7S_ML] + ppp[p7S_MG] + 
                          (is_glocal ? P7_DELTAT(xX, TSC(p7P_GM, k-1)) : 
                                       P7_DELTAT(xX, TSC(p7P_LM, k-1)));
	  else
	    dpc[p7S_ML] = ppp[p7S_ML] + ppp[p7S_MG] +                                    // Otherwise, you can't enter at the anchor at all, because ali is global on ia..ib. 
	                  ESL_MAX( ESL_MAX( P7_DELTAT( dpp[p7S_ML], TSC(p7P_MM, k-1)),   
			        	    P7_DELTAT( dpp[p7S_IL], TSC(p7P_IM, k-1))),  
		                            P7_DELTAT( dpp[p7S_DL], TSC(p7P_DM, k-1)));
	  dpc[p7S_IL] = -eslINFINITY;
	  dpc[p7S_DL] = -eslINFINITY;

	  xX = (is_glocal ? P7_DELTAT( P7_DELTAT(dpc[p7S_ML], TSC(p7P_MD, k)), TSC(p7P_DGE, k)) : dpc[p7S_ML]);
	  
	  /* Advance calculation of next Dk+1. */
	  if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) dc = P7_DELTAT( dpc[p7S_ML], TSC(p7P_MD, k));
	  else                                        dc = -eslINFINITY;

	  dpc += p7S_NSCELLS;
	  ppp += p7S_NSCELLS;

	  /* Remainder of initial row i0 can only be reached by delete paths,
	   * and we've already dealt with Mk->E exit 
	   */
	  for (z = z+1; z < sm->n[i]; z++, dpc += p7S_NSCELLS, ppp += p7S_NSCELLS )
	    {
	      dpc[p7S_ML] = -eslINFINITY;
	      dpc[p7S_IL] = -eslINFINITY;
	      dpc[p7S_DL] = dc;

	      /* Advance calculation of next Dk+1 */
	      k = sm->k[i][z];
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) dc = P7_DELTAT( dc, TSC(p7P_DD, k));
	      else                                  	  dc = -eslINFINITY;
	    } // end loop over sparse supercells z on initial row i0, glocal path
	} // end initialization of row i0


      else // for all other DOWN rows i, other than initialization on i0:
	{
	  dpp      = last_dpc;
	  last_dpc = dpc;
	  dc       = -eslINFINITY;

	  z = 0; while (z < sm->n[i]   && sm->k[i][z]   < env->arr[d].k0) z++;
	  y = 0; while (y < sm->n[i-1] && sm->k[i-1][y] < env->arr[d].k0) y++;
	  for (; z < sm->n[i]; z++, dpc += p7S_NSCELLS, ppp += p7S_NSCELLS )
	    {
	      k = sm->k[i][z]; 

	      /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	      while (y < sm->n[i-1] && sm->k[i-1][y]  < k-1) { y++; dpp += p7S_NSCELLS; }
	      if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) 
		dpc[p7S_ML] = ppp[p7S_ML] + ppp[p7S_MG] + 
		              ESL_MAX ( ESL_MAX( P7_DELTAT(dpp[p7S_ML], TSC(p7P_MM, k-1)),
						 P7_DELTAT(dpp[p7S_IL], TSC(p7P_IM, k-1))),
			                         P7_DELTAT(dpp[p7S_DL], TSC(p7P_DM, k-1)));
	      else
		dpc[p7S_ML] = -eslINFINITY;

	      /* Try to find cell i-1,k; then compute I(i,k) from it */
	      if (y < sm->n[i-1] && sm->k[i-1][y] < k) { y++; dpp += p7S_NSCELLS; }
	      if (y < sm->n[i-1] && sm->k[i-1][y] == k) 
		dpc[p7S_IL] = ppp[p7S_IL] + ppp[p7S_IG] +
		              ESL_MAX( P7_DELTAT(dpp[p7S_ML], TSC(p7P_MI, k)),
				       P7_DELTAT(dpp[p7S_IL], TSC(p7P_II, k)));
	      else
		dpc[p7S_IL] = -eslINFINITY;

	      /* Delayed store of Dk */
	      dpc[p7S_DL] = dc;
		  
	      xX =  ESL_MAX( xX,
			     (is_glocal ? P7_DELTAT( P7_DELTAT(dpc[p7S_ML], TSC(p7P_MD, k)), TSC(p7P_DGE, k)) : dpc[p7S_ML]));

	      /* Advance calculation of next Dk+1 */
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) 
		dc = ESL_MAX( P7_DELTAT( dpc[p7S_ML], TSC(p7P_MD, k)),
			      P7_DELTAT( dpc[p7S_DL], TSC(p7P_DD, k)));
	      else 
		dc = -eslINFINITY;
	    } // end recursion over supercells k in DOWN row i
	}


      /* <ppp> maintenance (3): skip trailing UP supercells when i is ASC DOWN/UP row */
      if (env->arr[d+1].i0 <= sm->seg[g].ib)                                // If there's an anchor below us, ASC row includes UP
	for (z = 0; z < sm->n[i] && sm->k[i][z] < env->arr[d+1].k0; z++)    // UP is 1..k0(d)-1 
	  ppp += p7S_NSCELLS;
    } // end of DOWN sector calculation


  /* <ppp> maintenance (4): If <d> is the last anchor in segment, skip trailing rows ib(d)+1..ib(g), which must be DOWN only rows */
  if (env->arr[d+1].i0 > sm->seg[g].ib)
    for (i = env->arr[d].ib+1; i <= sm->seg[g].ib; i++)
      for (z = sm->n[i]-1; z >= 0 && sm->k[i][z] >= env->arr[d].k0; z--)  // DOWN is k0(d)..M; order doesn't matter, so we can skip it backwards.
	ppp += p7S_NSCELLS;

  *mod_ppp = ppp;  // <ppp> either sits on start of next row (if there are more domains in this segment), or start of next segment. 
  *mod_dpc = dpc;  // <dpc> sits on ia(d) for next domain in the AEC DP
  *mod_g   = g;    // <g> is the segment we're in, or just finished.
  *mod_xX  = xX;   // <xX> is the total labeling posterior probability so far, for rows included in domains (nonhomologous rows are excluded from xX).
  return eslOK;
}


/*****************************************************************
 * 3. Choice selection functions for traceback
 *****************************************************************/

/* On entry to a select* function, <dpc> points at the current DP supercell i,k in AEC mx.
 * The select* function updates it to point to the new supercell that we trace back to.
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

  while (z >= 0 && sm->k[i][z]   >= env->arr[d].k0) { z--; dpp -= p7S_NSCELLS; }  // skip back over rest of DOWN(i) row. 
  z = sm->n[i-1]-1;
  while (sm->k[i-1][z] >  k-1)                      { z--; dpp -= p7S_NSCELLS; }  // skip back to i-1,k-1 cell on DOWN(i-1) row. no z>=0 check, because i-1,k-1 supercell must exist.
  ESL_DASSERT1(( z >= 0 && sm->k[i-1][z] == k-1 ));                               // connected supercell i-1,k-1 must exist, by construction.
  /* Now dpp is the i-1, k-1 supercell, which is k[i-1][y] in the sparse mask for row i-1  */
  path[0] = dpp[p7S_ML];
  path[1] = dpp[p7S_IL];
  path[2] = dpp[p7S_DL];
  
  //printf("down_m choosing from M=%.4f  I=%.4f  D=%.4f\n", path[0], path[1], path[2]);

  *mod_dpc = dpp;
  *mod_z   = z;
  return (env->arr[d].flags & p7E_IS_GLOCAL ? gstate[esl_vec_FArgMax(path, 3)] : lstate[esl_vec_FArgMax(path, 3)]);
}

/* up_m also gets used from the anchor cell, which is actually in the DOWN sector,
 * but traceback works almost the same as if it were in UP.
 */
static inline int
select_up_m(const P7_SPARSEMASK *sm, const P7_ENVELOPES *env, int d, int i, const float **mod_dpc, int *mod_z)
{
  const float *dpp       = *mod_dpc;
  int          z         = *mod_z;
  int          k         =  sm->k[i][z];
  int          gstate[3] = { p7T_MG, p7T_IG, p7T_DG };  
  int          lstate[3] = { p7T_ML, p7T_IL, p7T_DL };
  float        path[3];
  
  if (i < env->arr[d].i0) dpp -= (z+1) * p7S_NSCELLS;                    // skip remainder of UP(i), putting <dpp> on last supercell of prev row...
  else                    dpp -= p7S_NSCELLS;                            //  ... unless we're on anchor cell, in which case we know we just need to move -1 supercell.
  z    = sm->n[i-1]-1;                                                   // now starting from last sparse supercell on prev row i-1...
  while (sm->k[i-1][z] >= env->arr[d].k0)   z--;                         // skip z back until we're in the UP sector. z>=0 check not needed, we know there's at least 1 supercell 
  while (sm->k[i-1][z] >  k-1)            { z--; dpp -= p7S_NSCELLS; }   // skip back to i-1,k-1 cell on UP(i-1) row; again z>=0 check not needed
  ESL_DASSERT1(( z >= 0 && sm->k[i-1][z] == k-1 ));

  path[0] = dpp[p7S_ML];
  path[1] = dpp[p7S_IL];
  path[2] = dpp[p7S_DL];
  
  //printf("up_m choosing from M=%.4f  I=%.4f  D=%.4f\n", path[0], path[1], path[2]);

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

  while (z >= 0 && sm->k[i][z]  >= env->arr[d].k0) { z--; dpp -= p7S_NSCELLS; }  // skip back over rest of DOWN(i).
  z = sm->n[i-1]-1; while (sm->k[i-1][z] >  k)     { z--; dpp -= p7S_NSCELLS; }  // skip back to i-1,k cell on DOWN(i-1)
  ESL_DASSERT1(( z >= 0 && sm->k[i-1][z] == k ));                                // connected supercell i-1,k must exist, and must have >=1 value >= -inf, by construction.

  //printf("down_i choosing from M=%.4f  I=%.4f\n", dpp[p7S_ML], dpp[p7S_IL]);

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
  
  dpp -= (z+1) * p7S_NSCELLS;                                            // skip remainder of UP(i); we can't be the anchor cell (unlike select_up_m)
  z    = sm->n[i-1]-1; 
  while (sm->k[i-1][z] >= env->arr[d].k0)   z--;                         // skip z back until we're in the UP sector. z>=0 check not needed, we know there's at least 1 supercell 
  while (sm->k[i-1][z] >  k)              { z--; dpp -= p7S_NSCELLS; }   // skip back to i-1,k cell on UP(i-1) row; again z>=0 check not needed
  ESL_DASSERT1(( z >= 0 && sm->k[i-1][z] == k ));

  //printf("up_i choosing from M=%.4f  I=%.4f\n", dpp[p7S_ML], dpp[p7S_IL]);

  *mod_dpc = dpp;
  *mod_z   = z;
  if (dpp[p7S_ML] >= dpp[p7S_IL]) { return (env->arr[d].flags & p7E_IS_GLOCAL ? p7T_MG : p7T_ML); }
  else                            { return (env->arr[d].flags & p7E_IS_GLOCAL ? p7T_IG : p7T_IL); }
}

static inline int
select_d(const P7_ENVELOPES *env, int d, const float **mod_dpc)
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
 * Convention still holds here, even though E isn't explicitly stored
 * in DP matrix: *mod_dpc points at last sparse cell we calculated, which
 * is either the M that started domain d+1, or +1 off edge of matrix for d=D+1.
 *    Upon return, it's been updated to the supercell of the connected M.
 *    upon return, *ret_z is the index of the connected M in kp[].
 */
static inline int
select_e(const P7_SPARSEMASK *sm, const P7_ENVELOPES *env, int d, int i, const float **mod_dpc, int *ret_z)
{
  const float *dpp      = *mod_dpc - p7S_NSCELLS;   // now *dpp is on last sparse supercell of row ib(d).
  const float *best_dpp = dpp;
  int          best_z   = sm->n[i]-1;
  int          z;

  for (z = sm->n[i]-1; z >= 0 && sm->k[i][z] >= env->arr[d].k0; z--, dpp -= p7S_NSCELLS)
    if (dpp[p7S_ML] > best_dpp[p7S_ML]) { best_dpp = dpp; best_z = z; }

  *mod_dpc = best_dpp;
  *ret_z   = best_z;
  return ( env->arr[d].flags & p7E_IS_GLOCAL ? p7T_MG : p7T_ML );
}



/*****************************************************************
 * 4. Traceback, aec_trace()
 *****************************************************************/

static int
aec_trace(const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_ENVELOPES *env, const float *dpc, P7_TRACE *tr)
{
  int scur,sprv;
  int d;
  int i,k,k2;
  int z = 0;     // initialization is solely to quiet static code analyzers
  int status;

  if ((status = p7_trace_Append(tr, p7T_T, 0, 0)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, 0, 0)) != eslOK) return status;

  for (d = env->D; d >= 1; d--)  
    {
      k    = gm->M+1;

      /* Residues ib(d)+1 .. ia(d+1)-1 are assigned to C | J.
       * Sentinel at ia(D+1) = L+1, so we don't need a special case there.
       */
      for (i = env->arr[d+1].ia-1; i > env->arr[d].ib; i--)
	if ((status = p7_trace_Append(tr, (d == env->D ? p7T_C : p7T_J), 0, i)) != eslOK) return status;
      if ((status = p7_trace_Append(tr, p7T_E, 0, i)) != eslOK) return status;
      scur = p7T_E;
      // i now ib(d); scur now E; k is M+1 awaiting traceback from E

      /* DOWN matrix: trace back until we reach the anchor, which always exists */
      while (i > env->arr[d].i0 || k > env->arr[d].k0)
	{
	  //printf("I am now on i=%d, k=%d and scur=%d(%s)\n", i, k, scur, p7_trace_DecodeStatetype(scur));

	  switch (scur) {
	  case p7T_ML: case p7T_MG:   sprv = select_down_m(sm, env, d, i, &dpc, &z);  i--; k--;              break;
	  case p7T_IL: case p7T_IG:   sprv = select_down_i(sm, env, d, i, &dpc, &z);  i--;                   break;
	  case p7T_DL: case p7T_DG:   sprv = select_d     (    env, d,    &dpc);           k--; z--;         break;   
	  case p7T_E:                 sprv = select_e     (sm, env, d, i, &dpc, &z);       k = sm->k[i][z];  break;  
	  default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in sparse AEC traceback");
	  }
	  
	  /* When we trace E back, deducing Mk->E exit: record kb boundary, and 
	   * do right wing unfolding.
	   */
	  if ( scur == p7T_E )
	    {
	      if (env->arr[d].flags & p7E_IS_GLOCAL) 
		{
		  env->arr[d].kb = gm->M;
		  for (k2 = gm->M; k2 > k; k2--)
		    if ( (status = p7_trace_Append(tr, p7T_DG, k2, i)) != eslOK) return status;		  
		}
	      else env->arr[d].kb = k;
	    }
      
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

      /* We're on ia(d), in the optimal Mk. There may be still supercells on this dpc row < k,
       * but we're obligated to leave dpc on first supercell, i.e. +1 from ib(d-1), so select_e
       * will be able to work.
       * (We can only have supercells < k if ia(d) is in UP sector, which it usually is.
       *  Watch out for rare case where ia==i0, in which case first supercell is the anchor itself.)
       */
      if (env->arr[d].ia < env->arr[d].i0) dpc -= z * p7S_NSCELLS; 

      sprv = (env->arr[d].flags & p7E_IS_GLOCAL ? p7T_G : p7T_L );
      if ( sprv == p7T_G )  // Left wing unfolding, on glocal G->Mk entry.
	{
	  for (; k > 1; k--)
	    if ( (status = p7_trace_Append(tr, p7T_DG, k-1, i))             != eslOK) return status;
	  env->arr[d].ka = 1;
	}
      else env->arr[d].ka = k;
	    
      if ( (status = p7_trace_Append(tr, sprv,                     0, i)) != eslOK) return status;
      if ( (status = p7_trace_Append(tr, p7T_B,                    0, i)) != eslOK) return status;
      if ( (status = p7_trace_Append(tr, (d == 1 ? p7T_N : p7T_J), 0, i)) != eslOK) return status;
      // i is now ia(d), and we're on B, and dpc is +1 from ib(d-1)
    }

  for (i = env->arr[1].ia-1; i >= 1; i--)
    if ((status = p7_trace_Append(tr, p7T_N, 0, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_S, 0, 0))   != eslOK) return status;
  
  tr->M = sm->M;
  tr->L = sm->L;
  return p7_trace_Reverse(tr);
}





/*****************************************************************
 * 5. Footnotes
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
 */

/*****************************************************************
 * 6. Unit tests
 *****************************************************************/
#ifdef p7SPARSE_AEC_ALIGN_TESTDRIVE

#include "hmmer.h"

/* "compare_reference" test
 *
 */
static void
utest_compare_reference(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N, int be_verbose)
{
  char            msg[]  = "sparse_aec_align:: compare_reference unit test failed";
  ESL_SQ          *sq    = esl_sq_CreateDigital(abc);
  P7_BG           *bg    = p7_bg_Create(abc);
  P7_HMM          *hmm   = NULL;
  P7_PROFILE      *gm    = p7_profile_Create(M, abc);
  P7_ANCHORS      *anch  = p7_anchors_Create();
  P7_TRACE        *gtr   = p7_trace_Create();            // generated trace
  P7_TRACE        *rtr   = p7_trace_Create();            // trace from reference AEC
  P7_TRACE        *str   = p7_trace_Create();            // trace from sparse AEC
  P7_REFMX        *rfu   = p7_refmx_Create(M, 20);       // ASC Forward UP;  and AEC alignment.
  P7_REFMX        *rfd   = p7_refmx_Create(M, 20);       // ASC Forward DOWN
  P7_REFMX        *rpu   = p7_refmx_Create(M, 20);       // ASC posterior decoding UP
  P7_REFMX        *rpd   = p7_refmx_Create(M, 20);       // ASC posterior decoding DOWN
  P7_ENVELOPES    *renv  = p7_envelopes_Create();        // Envelopes calculated by reference impl
  P7_ENVELOPES    *senv  = p7_envelopes_Create();        // Envelopes calculated by sparse impl
  P7_SPARSEMASK   *sm    = p7_sparsemask_Create(M, M);
  P7_SPARSEMX     *asf   = p7_sparsemx_Create(NULL);     // sparse ASC Forward
  P7_SPARSEMX     *asx   = p7_sparsemx_Create(NULL);     // sparse ASC Backward, and AEC alignment
  P7_SPARSEMX     *asd   = p7_sparsemx_Create(NULL);     // sparse ASC Decoding
  float fsc;
  int   idx;

  /* Sample a random standard profile, no special constraints. */
  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);

  /* For N different test sequences... */
  for (idx = 0; idx < N; idx++)
    {
      /* ... emit sequence from model, using an arbitrary length model of <M>;
       * restrict the emitted sequence length to 6M, arbitrarily, to 
       * keep it down to something reasonable.
       */
      if ( p7_profile_SetLength(gm, M) != eslOK) esl_fatal(msg);
      do {
	esl_sq_Reuse(sq);
	if (p7_ProfileEmit(rng, hmm, gm, bg, sq, gtr) != eslOK) esl_fatal(msg);
      } while (sq->n > M * 6); 
      if ( p7_profile_SetLength(gm, sq->n) != eslOK) esl_fatal(msg);

      /* Use generated trace to create a plausible anchor set; doesn't need to be MPAS */
      if ( p7_trace_Index(gtr)                        != eslOK) esl_fatal(msg);
      if ( p7_anchors_SampleFromTrace(anch, rng, gtr) != eslOK) esl_fatal(msg);

      /* Reference calculations */
      if ( p7_ReferenceASCForward (sq->dsq, sq->n, gm, anch->a, anch->D, rfu, rfd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch->a, anch->D, rpu, rpd, NULL)               != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->a, anch->D, rfu, rfd, rpu, rpd, rpu, rpd) != eslOK) esl_fatal(msg);  // Reference ASC decoding can overwrite backwards matrix
      if ( p7_reference_Envelopes (sq->dsq, sq->n, gm, anch->a, anch->D, rpu, rpd, rfu, rfd, renv)     != eslOK) esl_fatal(msg);

      if ( p7_refmx_Reuse(rfu)                                   != eslOK) esl_fatal(msg);
      if ( p7_reference_AEC_Align( gm, renv, rpu, rpd, rfu, rtr) != eslOK) esl_fatal(msg);
      if ( p7_trace_Index(rtr)                                   != eslOK) esl_fatal(msg);

      /* Mark all cells in sparse mask */
      if ( p7_sparsemask_Reinit(sm, M, sq->n) != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_AddAll(sm)           != eslOK) esl_fatal(msg);

      /* Sparse calculations */
      if ( p7_sparse_asc_Forward (sq->dsq, sq->n, gm, anch->a, anch->D, sm, asf, &fsc)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Backward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asx, NULL)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Decoding(sq->dsq, sq->n, gm, anch->a, anch->D, fsc, asf, asx, asd) != eslOK) esl_fatal(msg);
      if ( p7_sparse_Envelopes   (sq->dsq, sq->n, gm, anch->a, anch->D, asf, asd, senv)     != eslOK) esl_fatal(msg);

      if ( p7_sparsemx_Reuse(asx)                       != eslOK) esl_fatal(msg);
      if ( p7_sparse_aec_Align(gm, asd, senv, asx, str) != eslOK) esl_fatal(msg);
      if ( p7_trace_Index(str)                          != eslOK) esl_fatal(msg);
   
      if (be_verbose) p7_trace_DumpAnnotated(stdout, rtr, gm, sq->dsq);
      if (be_verbose) p7_trace_DumpAnnotated(stdout, str, gm, sq->dsq);

      /* Test 1. Traces are identical */
      if ( p7_trace_Compare(rtr, str, 0.0) != eslOK) esl_fatal(msg);

      
      p7_anchors_Reuse(anch);
      p7_sparsemask_Reuse(sm);
      p7_envelopes_Reuse(senv);  p7_envelopes_Reuse(renv);
      p7_trace_Reuse(str);       p7_trace_Reuse(rtr);        p7_trace_Reuse(gtr);
      p7_refmx_Reuse(rfu);       p7_refmx_Reuse(rfd);
      p7_refmx_Reuse(rpu);       p7_refmx_Reuse(rpd);
      p7_sparsemx_Reuse(asf);    p7_sparsemx_Reuse(asx);     p7_sparsemx_Reuse(asd);
      esl_sq_Reuse(sq);
    }

  p7_anchors_Destroy(anch);
  p7_sparsemask_Destroy(sm);
  p7_envelopes_Destroy(senv); p7_envelopes_Destroy(renv);
  p7_trace_Destroy(str);      p7_trace_Destroy(rtr);         p7_trace_Destroy(gtr);
  p7_sparsemx_Destroy(asf);   p7_sparsemx_Destroy(asx);      p7_sparsemx_Destroy(asd);
  p7_refmx_Destroy(rfu);      p7_refmx_Destroy(rfd);
  p7_refmx_Destroy(rpu);      p7_refmx_Destroy(rpd);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_sq_Destroy(sq);
}


/* "singlemulti" test
 * 
 * Create <N> profile/seq/anchorset comparisons that have only one
 * possible path once it's anchor set constrained, using
 * p7_modelsample_SinglePathedASC(). For each, create a sparse mask that
 * includes that path (50% of the time by dirtying around the path;
 * 50% of the time by adding all cells). AEC alignment must recover
 * exactly that path.
 * 
 * To achieve this, the model is in multiglocal mode (no local
 * alignments). Anchor set constraint is what allows it to be
 * multidomain, because the anchor set defines a unique number of
 * domains. The key idea is that the profile generates a fixed number
 * of match states (possibly involving D state usage, using transition
 * probs of 0 to force using either M or D deterministically at each
 * k), and match states cannot generate some residue X. We construct a
 * target sequence that uses X for all N/J/C/I-aligned residues, and
 * non-X for match states. 
 * 
 * We test:
 *   1. The MEG trace is identical to the generated path.
 *   2. For each domain envelope, oa=ia, ob=ib, and these
 *      coords agree with the trace.
 *   3. In the case of a single domain (D=1), envelope score == 
 *      trace score.    
 */
static void
utest_singlemulti(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc, int N, int be_verbose)
{
  char           msg[] = "sparse_aec_align :: singlemulti unit test failed";
  P7_BG         *bg    = p7_bg_Create(abc);
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = NULL;
  ESL_DSQ       *dsq   = NULL;
  int            L;
  P7_TRACE      *gtr   = NULL;
  P7_ANCHOR     *anch  = NULL;
  int            D;
  float          gsc, fsc, bsc;
  P7_SPARSEMASK *sm    = NULL;
  P7_SPARSEMX   *asf   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asb   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asd   = p7_sparsemx_Create(NULL);
  P7_ENVELOPES  *env   = p7_envelopes_Create();
  P7_SPARSEMX   *aec   = NULL;
  P7_TRACE      *tr    = p7_trace_Create();
  float          tol   = 0.001;
  int            idx;
  int            d;

  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_SinglePathedASC(rng, M, bg, &hmm, &gm, &dsq, &L, &gtr, &anch, &D, &gsc) != eslOK) esl_fatal(msg);
      
      if (be_verbose) p7_trace_DumpAnnotated(stdout, gtr, gm, dsq);

      if (( sm = p7_sparsemask_Create(M, L)) == NULL) esl_fatal(msg);
      if (idx%2) { if ( p7_sparsemask_AddAll(sm)                 != eslOK) esl_fatal(msg); }
      else       { if ( p7_sparsemask_SetFromTrace(sm, rng, gtr) != eslOK) esl_fatal(msg); }

      if ( p7_sparse_asc_Forward (dsq, L, gm, anch, D, sm, asf, &fsc)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Backward(dsq, L, gm, anch, D, sm, asb, &bsc)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Decoding(dsq, L, gm, anch, D, fsc, asf, asb, asd) != eslOK) esl_fatal(msg);

      if ( p7_sparse_Envelopes(dsq, L, gm, anch, D, asf, asd, env) != eslOK) esl_fatal(msg);

      if (be_verbose) p7_spascmx_Dump(stdout, asd, anch, D);

      aec = asf;
      p7_sparsemx_Reuse(aec);
      if ( p7_sparse_aec_Align(gm, asd, env, aec, tr) != eslOK) esl_fatal(msg);
      if ( p7_trace_Index(tr)                         != eslOK) esl_fatal(msg);
      
      if (be_verbose) p7_spaecmx_Dump(stdout, aec, env);
      if (be_verbose) p7_trace_DumpAnnotated(stdout, tr, gm, dsq);

      /* Data objects should be valid, of course */
      if ( p7_spaecmx_Validate(aec, env, NULL)       != eslOK) esl_fatal(msg);
      if ( p7_trace_Validate(tr, gm->abc, dsq, NULL) != eslOK) esl_fatal(msg);

      /* Test 1. Traces are identical. */
      if ( p7_trace_Compare(gtr, tr, 0.0) != eslOK) esl_fatal(msg);

      /* Test 2. Domain #'s agree */
      if (! (gtr->ndom == D && env->D == D)) esl_fatal(msg);
      if (! (tr->ndom  == D))                esl_fatal(msg);

      /* Test 3. Envelope coords (and outer env coords) match trace */
      for (d = 1; d <= D; d++)
	{
	  if (! (env->arr[d].ia == gtr->sqfrom[d-1] &&   // beware: domains in trace are 0..ndom-1, off by one from env's 1..D
		 env->arr[d].ia ==  tr->sqfrom[d-1] &&
		 env->arr[d].ia == env->arr[d].oa)) esl_fatal(msg);

	  if (! (env->arr[d].ib == gtr->sqto[d-1]   &&
		 env->arr[d].ib ==  tr->sqto[d-1]   &&
		 env->arr[d].ib == env->arr[d].ob)) esl_fatal(msg);

	  if (! (env->arr[d].ka == gtr->hmmfrom[d-1] &&
		 env->arr[d].ka ==  tr->hmmfrom[d-1])) esl_fatal(msg);

	  if (! (env->arr[d].kb == gtr->hmmto[d-1] &&
		 env->arr[d].kb ==  tr->hmmto[d-1])) esl_fatal(msg);
	}
      
      /* Test 4. If D == 1, envelope score == trace score. */
      if (D == 1 &&  esl_FCompare(env->arr[1].env_sc, gsc, tol) != eslOK) esl_fatal(msg);

      p7_envelopes_Reuse(env);
      p7_sparsemx_Reuse(asf);
      p7_sparsemx_Reuse(asb);
      p7_sparsemx_Reuse(asd);
      p7_trace_Reuse(tr);

      p7_sparsemask_Destroy(sm);
      free(dsq);
      free(anch);
      p7_trace_Destroy(gtr);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);
    }


  p7_trace_Destroy(tr);
  p7_envelopes_Destroy(env);
  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(asb);
  p7_sparsemx_Destroy(asd);
  p7_bg_Destroy(bg);
}


#endif /*p7SPARSE_AEC_ALIGN_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/


/*****************************************************************
 * 7. Test driver
 *****************************************************************/
#ifdef p7SPARSE_AEC_ALIGN_TESTDRIVE
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
static char banner[] = "unit test driver for sparse AEC/MEG alignment inference";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  int             M    = 20; 

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_compare_reference(rng, M, abc, 10, /*be_verbose=*/FALSE);
  utest_singlemulti      (rng, M, abc, 10, /*be_verbose=*/FALSE); 

  fprintf(stderr, "#  status = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_AEC_ALIGN_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/



/*****************************************************************
 * 8. Example
 *****************************************************************/
#ifdef p7SPARSE_AEC_ALIGN_EXAMPLE

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
  { "-s",        eslARG_INT,      "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
 
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of 'back end' domain processing, sparse (production) version";

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
  float           fsc, vsc, asc_f, asc_b;
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
   * We get it from checkpointed Fwd/Bck/Decoding, the last step in the acceleration pipeline. 
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

  /* First pass analysis 
   * Uses two sparse matrices: <sxf>, <sxd>,
   * and also gets the unconstrained Viterbi trace, <tr>
   */
  p7_SparseViterbi (sq->dsq, sq->n, gm, sm,  sxf, tr, &vsc);
  p7_SparseForward (sq->dsq, sq->n, gm, sm,  sxf,     &fsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, sm,  sxd,     NULL);
  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxd,     sxd);

  /* MPAS algorithm 
   * <sxf>, <sxd> stay unchanged; after MPAS they can be reused (?)
   * <asf>   = ASC Forward matrix
   * <anch>  = optimal anchor set
   * <asc_f> = ASC Forward score.
   */
  p7_sparse_anchors_SetFromTrace(sxd, tr, vanch);
  p7_trace_Reuse(tr);
  p7_sparse_Anchors(rng, sq->dsq, sq->n, gm,
		    vsc, fsc, sxf, sxd, vanch,
		    tr, &wrk, ah,
		    asf, anch, &asc_f, 
		    NULL);


  /* Remaining ASC calculations; MPAS did <asf>. 
   * We don't need <asb>, except that current implementation doesn't
   * allow decoding to overwrite the backwards mx.
   */
  p7_sparse_asc_Backward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asb, &asc_b);
  p7_sparse_asc_Decoding(sq->dsq, sq->n, gm, anch->a, anch->D, asc_f, asf, asb, asd);
  
  //  p7_spascmx_DumpSpecials(stdout, asd, anch->a, anch->D);
  //p7_spascmx_Dump(stdout, asd, anch->a, anch->D); 

  /* Envelope determination. */
  p7_sparse_Envelopes(sq->dsq, sq->n, gm, anch->a, anch->D, asf, asd, env);

  p7_envelopes_Dump(stdout, env);

  /* AEC alignment. */
  p7_trace_Reuse(tr);
  p7_sparsemx_Reuse(sxf);
  p7_sparse_aec_Align(gm, asd, env, sxf, tr);

  if ( p7_spaecmx_Validate(sxf, env, NULL) != eslOK) esl_fatal("AEC matrix failed validation");

  p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);

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
#endif /*p7SPARSE_AEC_ALIGN_EXAMPLE */


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
