/* Sparse envelope scoring.
 * A modified version of Forward.
 */

#include "p7_config.h"

#include "easel.h"

#include "hmmer.h"
#include "p7_sparsemx.h"

/* Function:  p7_SparseEnvScore()
 * Synopsis:  Calculate an envelope Forward score.
 *
 * Purpose:   Calculate the score of an envelope <iae..ibe>, <kae..kbe>
 *            in the comparison of profile <gm> to sequence <dsq> of
 *            length <L>. That is, calculate the Forward score of this
 *            sequence as if it had only a single domain, constrained
 *            to lie within residues <iae..ibe> on the sequence and
 *            positions <kae..kbe> on the profile.
 *            
 *            Caller provides a sparse matrix <sx> for the
 *            calculation, initializaed the same way as a full sparse
 *            Forward calculation. We currently abuse this space for
 *            efficiency of the envelope calculation. Upon return, the
 *            contents of <sx> are undefined, except that it still
 *            retains its reference to a sparse mask <sx->sm>.
 *            
 *            The raw envelope score, in nats, is returned in
 *            <*ret_envsc>.
 *
 * Args:      dsq       - digital sequence, 1..L
 *            L         - length of dsq
 *            gm        - profile
 *            iae,ibe   - envelope coord bounds on sequence, 1<=iae<=ibe<=L
 *            kae,kbe   - envelope coord bounds on model, 1<=kae<=kbe<=M
 *            sx        - caller-allocated matrix for dsq x gm sparse comparison
 *            ret_envsc - RETURN: envelope score in nats
 *
 * Returns:   <eslOK> on success, and <*ret_envsc> is the envelope score.
 *
 * Throws:    (no abnormal error conditions)
 * 
 * Notes:  
 *   1. This is a modified version of sparse_fwdback.c::p7_SparseForward().
 *      The changes:
 *         - J state is ignored (treated as if tE->J=0). Even if the model
 *           <gm> is in multihit mode, envelope calculation is restricted
 *           to scoring a single domain.
 *         - 1..iae-1 are forced to be emitted by N; ibe+1..L are forced
 *           to be emitted by C.
 *         - Only sparse cells within <kae..kbe> envelope on each row are
 *           considered.
 *
 *   2. DP matrix <sx> has no structure other than being big enough
 *      for the entire sparse matrix for 1..L; therefore we know it is
 *      big enough for iae..ibe; so we just start immediately,
 *      rather than trying to position ourself at iae and leaving
 *      the correct # of initialized DP cells for 1..iae-1.
 *      We only store cells in the intersection of the sparse mask and 
 *      the envelope. 
 *      
 *   3. By construction of envelopes, we know iae..ibe must lie 
 *      within a single sparse segment. We do not have to check
 *      for rows with no included cells.  
 */

int
p7_SparseEnvScore(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
		  int iae, int ibe, int kae, int kbe,
		  P7_SPARSEMX *sx, float *ret_envsc)
{
  P7_SPARSEMASK *sm  = sx->sm;
  float         *xpc = sx->xmx;
  float         *dpc = sx->dp;
  float         *dpp;
  float         *last_dpc;
  const float const *tsc = gm->tsc;
  const float const *rsc;
  float          xE, xB, xL, xG, xN, xC;
  float          mlc, mgc;
  float          dlc, dgc;
  int            i,k;
  int            y,z;

  /* Initialize on <iae-1>. We know iae is in a sparse segment. */
  xpc[p7S_E]  =      -eslINFINITY;
  xpc[p7S_N]  = xN = (iae > 1 ? (iae-1) * gm->xsc[p7P_N][p7P_LOOP] : 0.0);
  xpc[p7S_J]  =      -eslINFINITY;   // J prohibited in envelope scoring
  xpc[p7S_B]  = xB = xN + gm->xsc[p7P_N][p7P_MOVE];
  xpc[p7S_L]  = xL = xB + gm->xsc[p7P_B][0];
  xpc[p7S_G]  = xG = xB + gm->xsc[p7P_B][1];
  xpc[p7S_C]  = xC = -eslINFINITY;
  xpc[p7S_JJ] =      -eslINFINITY;
  xpc[p7S_CC] =      -eslINFINITY;
  xpc += p7S_NXCELLS;

  /* If iae-1 is a stored sparse row, we need to initialize main states.
   * Remember, we only store the intersection of mask and envelope */
  if (sm->n[iae-1])   
    {
      dpp = dpc;
      y = 0; while (y < sm->n[iae-1] && sm->k[iae-1][y] < kae) y++; /* skip k<kae cells */
      for (; y < sm->n[iae-1] && sm->k[iae-1][y] <= kbe; y++)     /* i.e., for the intersection of mask and kae..kbe: */
	{
	  dpc[p7S_ML] = dpc[p7S_MG] = -eslINFINITY;
	  dpc[p7S_IL] = dpc[p7S_IG] = -eslINFINITY;
	  dpc[p7S_DL] = dpc[p7S_DG] = -eslINFINITY;
	  dpc += p7S_NSCELLS;
	}
    }
  /* Now xpc,dpc are on row iae. 
   * If there are sparse cells on iae-1, dpp is on the start of the row;
   * if not, dpp is invalid, but neither will it be accessed in the code below.
   */

  /* main recursion */
  for (i = iae; i <= ibe; i++)
    {
      rsc = gm->rsc[dsq[i]];
      last_dpc = dpc;
      dlc = dgc = xE = -eslINFINITY;

      y = 0; while (y < sm->n[i-1] && sm->k[i-1][y] < kae) y++; // enforce k>=kae envelope on prev row: skip to first valid sparse cell
      z = 0; while (z < sm->n[i]   && sm->k[i]  [z] < kae) z++; // ditto for k>=kae envelope on current row
      for ( ; z < sm->n[i] && sm->k[i][z] <= kbe; z++)          // i.e. for cells in the intersection of sparse mask and kae..kbe:
	{
	  k = sm->k[i][z];	/* next sparse cell to calculate: (i,k), kae<=k<=kbe */

	  /* Try to find i-1,k-1 in envelope, then compute M(i,k) from it */
	  mlc = xL + TSC(p7P_LM, k-1);
	  mgc = xG + TSC(p7P_GM, k-1);
	  while (y < sm->n[i-1] && sm->k[i-1][y]  < k-1) { y++; dpp += p7S_NSCELLS; } // tricky, the dpp increment: works because we can prove that y is either already at end (n[i-1]), or always in the mask/env intersection
	  if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) {                            // also tricky; note that for k=kae, k-1 will not have been stored on prev row (out of env), cannot access it; but we already initialized y to have k[i-1][y] >= kae
	    mlc = p7_FLogsum( p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MM, k-1),
					  dpp[p7R_IL] + TSC(p7P_IM, k-1)),
			      p7_FLogsum( dpp[p7R_DL] + TSC(p7P_DM, k-1),
					  mlc));        
	    mgc = p7_FLogsum( p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MM, k-1),
					  dpp[p7R_IG] + TSC(p7P_IM, k-1)),
			      p7_FLogsum( dpp[p7R_DG] + TSC(p7P_DM, k-1),
					  mgc));
	  }
	  dpc[p7S_ML] = mlc = MSC(k) + mlc;
	  dpc[p7S_MG] = mgc = MSC(k) + mgc;

	  /* Try to find cell i-1,k; then compute I(i,k) from it */
	  while (y < sm->n[i-1] && sm->k[i-1][y]  < k) { y++; dpp += p7S_NSCELLS; }
	  if    (y < sm->n[i-1] && sm->k[i-1][y] == k) {
	    dpc[p7S_IL] = p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MI,k),  dpp[p7R_IL] + TSC(p7P_II, k)); // + ISC(k) removed here; put back in if inserts have nonzero score
	    dpc[p7S_IG] = p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MI,k),  dpp[p7R_IG] + TSC(p7P_II, k)); // ditto
	  } else { dpc[p7S_IL] = dpc[p7S_IG] = -eslINFINITY; }

	  /* Local exit paths */
	  xE = p7_FLogsum(xE, p7_FLogsum(mlc, dlc));

	  /* delayed store of Dk; advance calculation of next D_k+1 */
	  dpc[p7S_DL] = dlc;
	  dpc[p7S_DG] = dgc;
	  if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) { /* is there a (i,k+1) cell to our right? */
	    dlc = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
	    dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	  } else { 		/* if not, we must add the {MD}k->Dk+1..E glocal exit path! (even from internal sparse cells, not just the last sparse cell)  */
	    xE  = p7_FLogsum( xE, TSC(p7P_DGE, k) + p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k)));
	    dlc = dgc = -eslINFINITY; // D paths do not propagate thru unmarked cells
	  } 
	  dpc += p7S_NSCELLS;
	} // end loop over sparse cells on current row

      xpc[p7S_E]      = xE;                             // we accumulated all MLk->...E exits, and all MGk->..E exits from MGk that had no DGk+1 to transit to
      xpc[p7S_N] = xN = xN + gm->xsc[p7P_N][p7P_LOOP];
      xpc[p7S_J]      = -eslINFINITY;                   // J prohibited in envelope calculation
      xpc[p7S_B] = xB = xN + gm->xsc[p7P_N][p7P_MOVE];
      xpc[p7S_L] = xL = xB + gm->xsc[p7P_B][0];
      xpc[p7S_G] = xG = xB + gm->xsc[p7P_B][1];
      xpc[p7S_C] = xC = p7_FLogsum( xE + gm->xsc[p7P_E][p7P_MOVE], xC + gm->xsc[p7P_C][p7P_LOOP]);
      xpc[p7S_JJ]     = -eslINFINITY;
      xpc[p7S_CC]     = -eslINFINITY;
      xpc += p7S_NXCELLS;
      dpp = last_dpc;		/* now dpp is valid on the current row, for sure */
    } // end loop over rows i=iae..ibe

  sx->type = p7S_ENVSCORE;
  *ret_envsc = xC 
               + gm->xsc[p7P_C][p7P_MOVE]  // C->T  
               + (sm->L == ibe ? 0.0f : (sm->L - ibe) * gm->xsc[p7P_C][p7P_LOOP]); // ibe+1..L are tCC's. Guard against 0*-inf=NaN.
  return eslOK;
}
