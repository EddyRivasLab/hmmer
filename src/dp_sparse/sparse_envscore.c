/* Sparse envelope scoring; a modified version of Forward.
 * 
 * The "envelope score" is the score of a single domain as if it were
 * the only domain in the target sequence.
 * 
 * Specifically, the envelope score is a restricted Forward score,
 * disallowing any core model paths outside the envelope
 * iae..ibe/kae..kbe, and disallowing multiple passes within the
 * envelope by disallowing (not evaluating) J state uses.  Thus
 * 1..iae-1 are forced to be accounted for by N state, and ibe+1..L by
 * C state; M1..Mkae-1 and Mkbe+1..m are unused.  (Dk's can be used on
 * glocal retracted entry/exit.)
 * 
 * The calculation is done by implicitly treating the intersection of
 * the envelope and the given sparse mask as the effective sparse
 * mask.  Sparse DP matrices from an envelope calculation have only
 * the cells in that intersection; if you were to navigate around in a
 * sparse DP envelope matrix, you'd have to make the same implicit
 * intersection as you traverse it.
 * 
 * Contents:
 *   1. Exact calculation, by sparse DP
 *   2. Approximate calculation, from a sparse Forward matrix already in hand
 *   3. Debugging tools
 *   4. Unit tests
 *   5. Test driver
 *   6. Example
 */
#include <p7_config.h>

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dsq.h"

#include "base/p7_profile.h"

#include "misc/logsum.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_envscore.h"

/*****************************************************************
 * 1. Exact calculation, by sparse DP
 *****************************************************************/

/* Function:  p7_SparseEnvscore()
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
 *            calculation, initialized the same way as a full sparse
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
p7_SparseEnvscore(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
		  int iae, int ibe, int kae, int kbe,
		  const P7_SPARSEMASK *sm, P7_SPARSEMX *sx, float *ret_envsc)
{
  const float * const tsc = gm->tsc;
  const float * rsc;
  float         *xpc;
  float         *dpc;
  float         *dpp;
  float         *last_dpc;
  float          xE, xB, xL, xG, xN, xC;
  float          mlc, mgc;
  float          dlc, dgc;
  int            i,k;
  int            y,z;
  int            status;

  /* Contract checks / arg validation */
  ESL_DASSERT1( (sm->L == L) );
  ESL_DASSERT1( (sm->M == gm->M) );
  ESL_DASSERT1( (iae >= 1 && iae <= sm->L) );
  ESL_DASSERT1( (ibe >= 1 && ibe <= sm->L && ibe >= iae) );
  ESL_DASSERT1( (kae >= 1 && kae <= sm->M) );
  ESL_DASSERT1( (kbe >= 1 && kbe <= sm->M && kbe >= kae) );

  /* Assure that <sx> is allocated large enough; we might be reusing
   * it.  (We could be even more efficient, since we only store the
   * intersection of the mask and the envelope, but here we're
   * allocating for the entire mask.)
   */
  if ( (status = p7_sparsemx_Reinit(sx, sm)) != eslOK) return status;
  sx->type = p7S_ENVSCORE;

  /* Setting these ptrs has to come after realloc we may have just done: */
  xpc = sx->xmx;
  dpc = sx->dp;

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
   * if not, dpp is invalid, but neither will it be accessed in the 
   * first pass through the loop below; then it'll get set to last_dpc.
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
	  if (k < kbe && z < sm->n[i]-1 && sm->k[i][z+1] == k+1) { /* is there a (i,k+1) cell to our right? */
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

  *ret_envsc = xC 
               + gm->xsc[p7P_C][p7P_MOVE]  // C->T  
               + (sm->L == ibe ? 0.0f : (sm->L - ibe) * gm->xsc[p7P_C][p7P_LOOP]); // ibe+1..L are tCC's. Guard against 0*-inf=NaN.
  return eslOK;
}
/*--------------- end, exact DP calculation ---------------------*/



/*****************************************************************
 * 2. Approximate calculation, from sparse Fwd mx already in hand
 *****************************************************************/

/* Function:  p7_sparsemx_SparseEnvscoreApprox()
 * Synopsis:  Envelope score by fast Forward endpoint approximation method.
 *
 * Purpose:   Implements envelope scoring by the fast Forward endpoint
 *            approximation method. Returns the envelope score of
 *            envelope <iae..ibe> in the target sequence for which
 *            we've calculated a Forward matrix <sxf>, comparing it to
 *            profile <gm>. The "envelope score" is defined as the
 *            Forward score for the entire sequence <1..L>, subject to
 *            the condition that there is only a single domain that lies
 *            somewhere within <iae..ibe> on the target sequence.
 *
 *            <*ret_envsc> is a raw Forward score in nats.  The caller
 *            still needs to include the null model(s) contribution to
 *            this score and convert it to bits, before calculating
 *            P-values or reporting it in output.
 *            
 *            If the envelope contains negligible probability mass in
 *            the core model, the calculation here is subject to
 *            numerical error; the function returns <eslEINACCURATE>
 *            in this case, but still sets <*ret_envsc> to what it
 *            estimated. A caller that wants the right answer must check
 *            for the <eslINACCURATE>/<eslOK> status return.
 *            
 *            Caller has determined that the envelope <iae..ibe>
 *            contains only a single domain (negligible probability
 *            mass in any more than that); otherwise, this technique
 *            does not work. 
 *            
 *            The model must be in either standard multihit
 *            (<tCC=tJJ=tNN>, <tEJ=tEC>) mode or standard unihit mode
 *            (<tCC=tNN>, <tEC=0>, <tEJ=-inf>) for the approximation
 *            to be valid. (Even if valid, it may be inaccurate.)
 *            
 *            For a single-domain envelope <iae..ibe>, multihit with
 *            equal CC,JJ,TT transition scores, the score <delta>
 *            contributed by the <iae..ibe> region is
 *            
 *            $\delta = C_j 
 *                      + \log \[ 1 - e^{C_{i-1} + L_e \tau_{CC} - C_j} \]
 *                      - \log \[ e^{N_{i-1}} + e^{C_{i-1}} \]$
 *            
 *            where $L_e$ is the envelope length <ibe-iae+1>. We then
 *            add the <S->N...> contribution for <iae-1> flanking
 *            residues and <...C->T> contribution for <L-ibe> flanking
 *            residues to convert <delta> to a raw Forward score for
 *            the whole sequence, constrained to a single domain in
 *            the envelope.
 *
 * Args:      gm        - profile used to create the Forward matrix
 *            sxf       - Forward matrix calculated for gm x seq comparison
 *            iae       - left endpoint of envelope, 1..L on target seq
 *            ibe       - right endpoint of envelope, 1..L on target seq
 *            ret_envsc - RETURN: raw envelope Forward score, in nats
 *
 * Returns:   <eslOK> on success, and <*ret_envsc> contains the raw 
 *            envelope Forward score in nats.
 *            
 *            <eslEINACCURATE> if the value in <*ret_envsc> may be inaccurate
 *            because of numerical error. Caller may want to call the exact
 *            envelope score calculation in this case.
 *
 * Xref:      SRE:J10/43-45.
 */
int
p7_SparseEnvscoreApprox(P7_PROFILE *gm, P7_SPARSEMX *sxf, int iae, int ibe, float *ret_envsc)
{
  const float *xc = sxf->xmx;
  float Ci1, Ni1, Cj;
  float deltaC;
  float Le;
  int   g, ia, ib, last_ib;
  float envsc;
  int   is_unihit = (gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY ? TRUE : FALSE );

  /* Contract check / arg validation */
  ESL_DASSERT1( (sxf->type == p7S_FORWARD) );
  ESL_DASSERT1( ( gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY ||  gm->xsc[p7P_E][p7P_LOOP] == gm->xsc[p7P_E][p7P_MOVE]) ); /* unihit vs. multihit standard parameterizations */

  /* Accessing the {CN}(iae-1) and C(ibe) entries is complicated because of sparse storage; please forgive. */
  /* Find first stored special row i>=iae; get C(i-1), N(i-1) from there; reset iae=i */
  for (g = 1; g <= sxf->sm->S; g++)
    {
      ia = sxf->sm->seg[g].ia;
      ib = sxf->sm->seg[g].ib;
      
      if      (iae < ia)   {                             iae = ia;  Ci1 = xc[p7S_C]; Ni1 = xc[p7S_N]; break; }
      else if (iae <= ib)  { xc += (iae-ia)*p7S_NXCELLS; ia  = iae; Ci1 = xc[p7S_C]; Ni1 = xc[p7S_N]; break; }
      else    xc += (ib-ia+2) * p7S_NXCELLS;
    }
  /* now xc is on row ia-1, in segment g, which ends at ib. check if ib is in that segment. */
  
  /* Find last stored row j<=ibe, get C(j) from there; reset ibe=j */
  if (ib < iae)  { *ret_envsc = -eslINFINITY; return eslOK; }
  if (ibe <= ib) {
    xc += (ibe - ia + 1) * p7S_NXCELLS;
    Cj = xc[p7S_C];
  } else {
    xc     += (ib-ia+1) * p7S_NXCELLS; // xc was on ia-1. Move it to (ib) in that previous segment
    last_ib = ib;	 	       /* Now we're looking for ibe' = max_i i <= ibe such that xc[i] is stored  */
    for (g = g+1; g <= sxf->sm->S; g++)
      {
	ia = sxf->sm->seg[g].ia;
	ib = sxf->sm->seg[g].ib;
	
	if      (ibe < ia-1) { ibe = last_ib;                 Cj = xc[p7S_C]; break; }
	else if (ibe <= ib)  { xc += (ibe-ia+2)*p7S_NXCELLS;  Cj = xc[p7S_C]; break; }
	else  {
	  xc += (ib-ia+2)*p7S_NXCELLS;
	  last_ib = ib;
	}
      }
    if (g == sxf->sm->S+1) { ibe = last_ib; Cj = xc[p7S_C]; }
  }
  /* now iae,ibe may have been moved, to correspond to outermost stored rows in the envelope */


  /* First part of envsc is Cj + log(1-exp(-deltaC)), and the log term needs
   * special numerical handling; using lim x->0 1-e^-x = x for small deltaC,
   * lim x->0 log (1-x) = -x for large deltaC
   */
  envsc  = Cj;			/* first term. */
  Le     = ibe - iae + 1;
  deltaC = Cj - Ci1 - Le * gm->xsc[p7P_C][p7P_LOOP];

  if      (deltaC  < 0) ESL_EXCEPTION(eslEINCONCEIVABLE, "no, something's wrong, this intermediate term is >= 0 by construction");
  if      (deltaC == 0)                  envsc = -eslINFINITY;
  else if (deltaC < eslSMALLX1)          envsc += logf(deltaC);        // logf(deltaC) -> -eslINFINITY, for small deltaC ->0
  else if (exp(-1.*deltaC) < eslSMALLX1) envsc -= expf(-1.*deltaC);    // expf(-deltaC) -> 0, for large deltaC -> inf
  else                                   envsc += logf(1.-exp(-1.*deltaC));

  /* third term is a standard log-sum-exp-2, may as well use our fast approx */
  envsc -= (is_unihit ? Ni1 : p7_FLogsum(Ci1, Ni1) );
  
  /* left flank: envsc already includes ...N->B->...; add S->N..N flank. S->N is 1.0.  */
  envsc += gm->xsc[p7P_N][p7P_LOOP] * (iae-1);
  
  /* right flank: envsc includes E->C; add C..C->T flank */
  envsc += gm->xsc[p7P_C][p7P_LOOP] * (sxf->sm->L - ibe);
  envsc += gm->xsc[p7P_C][p7P_MOVE];

  *ret_envsc = envsc;

  /* if Cj and Ci1 path term are very close, deltaC ~ 0, and the accuracy of
   * our calculation is compromised by a catatrophic cancellation that I
   * don't see how to wriggle out of. Warn that the calculation is inaccurate
   * by setting eslEINACCURATE flag.
   */
  return (deltaC >= 0.1f ? eslOK : eslEINACCURATE);
}
/*-------------- end, approximate method ------------------------*/

/*****************************************************************
 * x. Debugging tools
 *****************************************************************/

/* Function:  p7_sparse_envscore_Dump()
 * Synopsis:  Dump a P7_SPARSEMX from an envscore calculation.
 *
 * Purpose:   To stream <ofp>, dump a P7_SPARSEMX <sxe> that was calculated by
 *            <p7_SparseEnvscore()> for envelope <iae..ibe,kae..kbe>. 
 * 
 *            This is an appropriately specialized version of
 *            <p7_sparsemx_Dump()>; the envscore calculation only
 *            calculates and stores the intersection of the sparse
 *            mask and the envelope, so the resulting <P7_SPARSEMX> is
 *            not suitable for normal sparsemx functions.
 */
int
p7_sparse_envscore_Dump(FILE *ofp, P7_SPARSEMX *sxe, int iae, int ibe, int kae, int kbe)
{
  int   storder[p7S_NSCELLS] = { p7S_ML, p7S_IL, p7S_DL, p7S_MG, p7S_IG, p7S_DG };
  const P7_SPARSEMASK *sm  = sxe->sm;
  float               *dpc = sxe->dp;
  float               *xc  = sxe->xmx;
  int width     = 9;
  int precision = 4;
  int i,k,x,y,z;
  int zb;
  int nin;

  fprintf(ofp, "       ");
  for (k = kae; k <= kbe;         k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0;   x < p7S_NXCELLS;  x++) fprintf(ofp, "%*s ", width, p7_sparsemx_DecodeSpecial(x));
  fprintf(ofp, "\n");

  fprintf(ofp, "       ");
  for (k = kae; k <= kbe;        k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0;   x < p7S_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");

  for (i = iae-1; i <= ibe; i++)
    {
      zb  = 0;  while (zb < sm->n[i] && sm->k[i][zb] < kae) zb++;

      for (x = 0; x < p7S_NSCELLS; x++)
	{
	  fprintf(ofp, "%3d %s ", i, p7_sparsemx_DecodeState(storder[x]));
	  for (z = zb, k = kae; k <= kbe; k++)
	    {
	      while (z < sm->n[i] && sm->k[i][z]  < k) z++; 
	      if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + (z-zb)*p7S_NSCELLS + storder[x]));
	      else                                     fprintf(ofp, "%*s ",   width, ".....");
	    }
	  if (x==0) for (y = 0; y < p7S_NXCELLS; y++) fprintf(ofp, "%*.*f ", width, precision, xc[y]);
	  fputc('\n', ofp);
	}
      fputc('\n', ofp);

      nin = 0; while (zb < sm->n[i] && sm->k[i][zb] <= kbe) { nin++; zb++; }

      dpc += nin * p7R_NSCELLS;
      xc  += p7R_NXCELLS;
    }
  return eslOK;
}

/* Function:  p7_sparse_envscore_IntersectedMask()
 * Synopsis:  Create an explicit P7_SPARSEMASK as a mask/envelope intersection.
 *
 * Purpose:   The envelope score calculation implicitly intersects the
 *            sparse mask and the envelope. As an independent test of
 *            the <p7_SparseEnvscore()> calculation, you can
 *            explicitly intersect the mask and the envelope, and run
 *            a standard <p7_SparseForward()> calculation with that
 *            intersection; you get the same result. The
 *            intersected-mask unit test uses this. We expose this
 *            function outside the unit tests because it's useful for
 *            other debugging (including it in the example main(), for
 *            example).
 *
 * Args:      osm - the original/old sparse mask
 *            iae,ibe,kae,kbe - envelope coords
 *            ret_sm - RETURN: newly created mask that's the intersection.
 *
 * Returns:   <eslOK> on success. <*ret_sm> is allocated here, and should
 *            be free'd by the caller using <p7_sparsemask_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation failure. Now <*ret_sm> is <NULL>.
 */
int
p7_sparse_envscore_IntersectedMask(P7_SPARSEMASK *osm, int iae, int ibe, int kae, int kbe, P7_SPARSEMASK **ret_sm)
{
  P7_SPARSEMASK *sm;
  int            i,z;

  ESL_DASSERT1( (iae >= 1 && iae <= osm->L) );
  ESL_DASSERT1( (ibe >= 1 && ibe <= osm->L && ibe >= iae) );
  ESL_DASSERT1( (kae >= 1 && kae <= osm->M) );
  ESL_DASSERT1( (kbe >= 1 && kbe <= osm->M && kbe >= kae) );

  if ( (sm = p7_sparsemask_Create(osm->M, osm->L)) == NULL) { *ret_sm = NULL; return eslEMEM; }

  /* Remember, API for Add() is designed for production use, in a vector Backwards;
   * so we go back through i's and z's, and pass (q,r) vector coords
   */
  for (i = ibe; i >= iae; i--)
    if (osm->n[i])
      {
	p7_sparsemask_StartRow(sm, i);
	z = osm->n[i]-1; 	
	while (z >= 0 && osm->k[i][z] >  kbe)   z--;
	while (z >= 0 && osm->k[i][z] >= kae) { p7_sparsemask_Add(sm, (osm->k[i][z]-1)%osm->Q, (osm->k[i][z]-1)/osm->Q); z--; }
	p7_sparsemask_FinishRow(sm);
      }
  p7_sparsemask_Finish(sm);
  
  *ret_sm = sm; 
  return eslOK;
}

/*--------------- end, debugging tools --------------------------*/


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7SPARSE_ENVSCORE_TESTDRIVE
#include "esl_vectorops.h"

#include "base/p7_bg.h"
#include "base/p7_masstrace.h"

#include "build/modelsample.h"
#include "search/modelconfig.h"
#include "misc/emit.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"

#include "dp_sparse/sparse_viterbi.h"
#include "dp_sparse/sparse_fwdback.h"
#include "dp_sparse/sparse_decoding.h"
#include "dp_sparse/sparse_trace.h"
#include "dp_sparse/sparse_masstrace.h"


/* The unit tests use two different alternative techniques to
 * calculate an envelope-restricted Forward score.
 * 
 * The intersected_mask test constructs a new sparse mask by
 * explicitly intersecting the envelope and the original mask.
 * 
 * The zeroed_emission test constructs a new target sequence and new
 * HMM, where x_1..iae-1,iae+1..L are a residue that the HMM can't
 * emit, and M_1..kae-1,kbe+1..M have e(x)=1.0 for a residue that the
 * target sequence doesn't contain.
 * 
 * Either way, a recalculation of sparse Forward with the (possibly
 * modified) mask on the (possibly modified) target sequence results
 * in a envelope-restricted Forward score. This new Forward score may
 * still contain J uses; this can be prevented by testing against
 * unihit profiles that prohibit J's. 
 * 
 * Let:
 *    e1 = an exact env score w/ original target, original mask,
 *         on some envelope iae..ibe/kae..kbe
 *    a1 = approx env score for same, calculated from the original
 *         Forward matrix
 *    f2 = new Fwd score, recalculated by one of the two methods
 *         of restricting paths to the envelope
 *    e2 = new exact env score, recalculated, using the modified mask or 
 *         modified target seq
 *    a2 = new approx env score, obtained from the recalculated Fwd 
 *         matrix that f2 came from
 *         
 * To see how these scores relate, here's a table of what additional
 * parts of the path ensemble they include:
 * 
 *                                     e1    a1     f2     e2    a2
 * can use J (multihit within env)     no  YES[1] YES[1]   no  YES[1]
 * can use M outside kae..kbe          no  YES[2]   no     no    no
 * can use core outside iae..ibe       no  YES[3]   no     no    no
 * 
 *  [1] If model is unihit, these become "no".
 *  [2] If kae..kbe is 1..M (glocal), this becomes "no".
 *  [3] The approx score may include path mass for an overlapping
 *      ali that starts before i (so is not in Ci-1 yet) and ends before j 
 *      (so is in Cj); Cj-C{i-1} is an overestimate of domain
 *      mass inside the envelope for this reason.
 *  
 * So:
 *   in the general case:  (e1 == e2) <= (f2 == a2) <= a1
 *   w/ unihit profile:    (e1 == f2 == e2 == a2)   <= a1
 *   
 * [xref SRE:J10/68]
 * 
 * The zeroed-emission test has a special wrinkle, arising from a
 * convention we use in sparse DP. The convention is that exit paths
 * through Dk's are counted for glocal (Mk->DGk+1...E) even for
 * non-marked DG's, whereas paths through DLk's are only counted when
 * they're marked. The zeroed-emission test disallows paths outside
 * the envelope by using strategically placed prob 0 emissions; this
 * strategy does not prohibit probability mass that flows through
 * unmarked DLk..->E paths, as in MLkbe-> DLkbe+1...E. Therefore to
 * make our tests work we have to prohibit local paths through
 * unmarked DLk's: for example, either use glocal-only, or set
 * t(Mk->Dk+1)=0 for all k. [xref J10/70-71]
 */

/* First step of the zeroed-emissions utest is to sample an HMM in
 * which all A/C {MI}k emissions are 0; this works for both DNA and AMINO
 * models. Assume that in the digital alphabet A=0,C=1 (this is true
 * for both DNA and AMINO).
 */
static void
prohibit_all_ac_in_hmm(P7_HMM *hmm)
{
  int k;
  for (k = 1; k <= hmm->M; k++)
    {
      hmm->mat[k][0] = hmm->mat[k][1] = 0.0f;
      esl_vec_FNorm(hmm->mat[k], hmm->abc->K);

      hmm->ins[k][0] = hmm->ins[k][1] = 0.0f;
      esl_vec_FNorm(hmm->ins[k], hmm->abc->K);
    }
}

/* and we also prohibit any C's in the sampled target sequence,
 * so the post-doctoring HMM (with ek(C)=1 for Mk states outside
 * the model envelope) can't generate any residue in the target.
 * (C's can't be generated by the HMM we use, but they could be 
 * generated by the background distribution in NN/CC/JJ emissions).
 */
static void
prohibit_all_c_in_seq(ESL_SQ *sq)
{
  int i;
  for (i = 1; i <= sq->n; i++)
    if (sq->dsq[i] == 1) sq->dsq[i] = 3; /* D for AMINO; G for DNA. Doesn't matter what we turn it to, really, as long as it's not C */
}


/* and because of the detail of nonzero exit mass through unmarked
 * DLks (see above; and J10/70-71), when the zeroed-emission test
 * uses local paths, it prohibits all prob mass in deletes, by 
 * setting t(Mk->Dk+1)=0 for all k.
 */
static void 
prohibit_delete_transitions(P7_HMM *hmm)
{
  int k;
  for (k = 1; k < hmm->M; k++)
    {
      hmm->t[k][p7H_MD] = 0.0;
      esl_vec_FNorm(P7H_TMAT(hmm, k), p7H_NTMAT);
    }
}
    

/* After we've sampled a target seq from the A/C-prohibited profile,
 * and done sparse calculations against it to determine envelopes, now
 * for each envelope we construct a modified dsq that has A's 
 * for all i=1..iae-1,ibe+1..L outside the sequence envelope.
 * We aren't too concerned with runtime, in the unit tests, to
 * we allocate that dsq here and don't worry about alloc/free cycles.
 */
static void
create_envdoctored_dsq(const ESL_SQ *sq, int iae, int ibe, ESL_DSQ **ret_dsq)
{
  ESL_DSQ *dsq = NULL;
  int      i;

  esl_dsq_Clone(sq->dsq, sq->n, &dsq);
  for (i = 1;     i < iae;    i++) dsq[i] = 0;
  for (i = ibe+1; i <= sq->n; i++) dsq[i] = 0;
  *ret_dsq = dsq;
}

/* and simiilarly, create a new profile, in which we've prohibited
 * emissions in Mk=1..kae-1,kbe..M by setting ek(C) = 1.0, and ek(all
 * other residues) = -inf. Because the target seq contains no C's,
 * this prohibits paths from using any of the model outside the
 * kae..kbe envelope. Configure the new profile the same way the old
 * <ogm> was.
 */
static void
create_envdoctored_profile(const P7_HMM *ohmm, const P7_PROFILE *ogm, const P7_BG *bg, int kae, int kbe, P7_PROFILE **ret_gm)
{
  P7_HMM     *new = p7_hmm_Clone(ohmm);
  P7_PROFILE *gm  = p7_profile_Create(ogm->M, ogm->abc);
  int     k;

  for (k = 1;     k < kae;     k++) { esl_vec_FSet(new->mat[k], new->abc->K, 0.0); new->mat[k][1] = 1.0; }
  for (k = kbe+1; k <= new->M; k++) { esl_vec_FSet(new->mat[k], new->abc->K, 0.0); new->mat[k][1] = 1.0; }
  
  p7_profile_ConfigCustom(gm, new, bg, ogm->L, ogm->nj, ogm->pglocal);
  p7_hmm_Destroy(new);
  *ret_gm = gm;
}


/* the two utests, intersected-mask and zeroed-emissions, have very
 * similar structures, so we have one engine that runs both of them,
 * in different combinations of unihit/multihit, glocal/dual. The
 * first arg to the engine is TRUE for the intersected-mask test,
 * FALSE for the zeroed-emissions test.
 */
static void
utest_engine(int do_intersected_mask, char *msg, ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int do_unihit, int do_glocal)
{
  int            do_zeroed_emissions = (do_intersected_mask ? FALSE : TRUE ); /* exclusive, one or the other */
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = p7_profile_Create(M, abc);
  P7_OPROFILE   *om    = p7_oprofile_Create(M, abc);
  ESL_SQ        *sq    = esl_sq_CreateDigital(abc);
  P7_CHECKPTMX  *ox    = p7_checkptmx_Create(M, L, ESL_MBYTES(32));
  P7_SPARSEMASK *sm    = p7_sparsemask_Create(M, L);
  P7_SPARSEMX   *sx    = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxf   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxb   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxd   = p7_sparsemx_Create(NULL);
  P7_TRACE      *vtr   = p7_trace_CreateWithPP();
  P7_MASSTRACE  *mt    = p7_masstrace_CreateSlim(M, L);
  P7_SPARSEMASK *sm_restricted  = NULL;
  P7_PROFILE    *gm_restricted  = NULL;
  ESL_DSQ       *dsq_restricted = NULL;
  float          fsc;
  float          e1_sc, a1_sc, f2_sc, e2_sc, a2_sc;
  int            iae,ibe,kae,kbe;
  int            idx;
  int            d;
  char           errbuf[eslERRBUFSIZE];
  float          tol  = 0.001;
  float          tol2 = 1.0;	/* relaxed tolerance for Approx() calculations that report they're inaccurate */
  int            a1_status, a2_status;

  /* Sample a profile. Config as requested for the utest version. */
  if (p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if (do_zeroed_emissions) {
    prohibit_all_ac_in_hmm(hmm);
    if (! do_glocal) prohibit_delete_transitions(hmm);
  }
  if      (do_unihit && do_glocal) { if ( p7_profile_ConfigUniglocal(gm, hmm, bg, L)           != eslOK) esl_fatal(msg); }
  else if (do_unihit)              { if ( p7_profile_ConfigCustom   (gm, hmm, bg, L, 0.0, 0.5) != eslOK) esl_fatal(msg); }
  else if (do_glocal)              { if ( p7_profile_ConfigGlocal   (gm, hmm, bg, L)           != eslOK) esl_fatal(msg); }
  else                             { if ( p7_profile_Config         (gm, hmm, bg)              != eslOK) esl_fatal(msg); }
  if (p7_oprofile_Convert(gm, om) != eslOK) esl_fatal(msg);
  
  for (idx = 0; idx < N; idx++)
    {
      /* Sample a sequence from the profile */
      if ( p7_profile_SetLength(gm, L)  != eslOK) esl_fatal(msg);   /* config to generate mean length of L (length was probably reset by last emitted seq) */
      do {
	esl_sq_Reuse(sq);
	p7_ProfileEmit(rng, hmm, gm, bg, sq, NULL);
      } while (sq->n < 1 || sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
      if (do_zeroed_emissions) prohibit_all_c_in_seq(sq);
      if ( p7_profile_SetLength      (gm, sq->n) != eslOK) esl_fatal(msg);
      if ( p7_oprofile_ReconfigLength(om, sq->n) != eslOK) esl_fatal(msg);

      /* Fwd/Bck local filter to calculate the sparse mask */
      if ( p7_checkptmx_Reinit(ox, M, sq->n)                                 != eslOK) esl_fatal(msg);
      if ( p7_ForwardFilter (sq->dsq, sq->n, om, ox, /*fsc=*/NULL)           != eslOK) esl_fatal(msg);
      if ( p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH) != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_Validate(sm, errbuf)                                != eslOK) esl_fatal("%s\n  %s", msg, errbuf);

      /* Sparse DP calculations */
      if ( p7_SparseViterbi   (sq->dsq, sq->n, gm, sm, sx,  vtr, NULL) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward   (sq->dsq, sq->n, gm, sm, sxf,     &fsc)  != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward  (sq->dsq, sq->n, gm, sm, sxb,      NULL) != eslOK) esl_fatal(msg);
      if ( p7_SparseDecoding  (sq->dsq, sq->n, gm, sxf, sxb, sxd)      != eslOK) esl_fatal(msg);
      p7_sparsemx_TracePostprobs(sxd, vtr); /* selecting domain anchors requires pp annotation of the trace */
      p7_trace_Index(vtr);
      p7_sparsemx_Reuse(sx);	/* we'll reuse it for mass trace dp below, like we do in production pipeline */

      //p7_sparsemx_Dump(stdout, sxf);

      /* Test each domain */
      for (d = 0; d < vtr->ndom; d++)
	{
	  if (p7_SparseMasstrace(sq->dsq, sq->n, gm, sxf, sxb, vtr, vtr->anch[d], 0.1, sx, mt, &iae, &ibe, &kae, &kbe) != eslOK) esl_fatal(msg);
	  p7_sparsemx_Reuse(sx);

	  a1_status = p7_SparseEnvscoreApprox(gm, sxf, iae, ibe, &a1_sc); /* EnvscoreApprox can return eslOK or eslEINACCURATE */
	  if (a1_status != eslOK && a1_status != eslEINACCURATE) esl_fatal(msg);

	  if (p7_SparseEnvscore(sq->dsq, sq->n, gm, iae, ibe, kae, kbe, sm, sx, &e1_sc) != eslOK) esl_fatal(msg);
	  //printf("# Exact: for domain %d, iae..ibe=%d..%d, kae..kbe=%d..%d\n", d, iae, ibe, kae, kbe);
	  //p7_sparse_envscore_Dump(stdout, sx, iae, ibe, kae, kbe);
	  p7_sparsemx_Reuse(sx);

	  /* Impose the restrictions, appropriate to our test */
	  if (do_intersected_mask)
	    {			                                                          /* In the intersected-mask utest: */
	      p7_sparse_envscore_IntersectedMask(sm, iae, ibe, kae, kbe, &sm_restricted); /* we intersect the sparse mask with the envelope... */
	      dsq_restricted = sq->dsq;                                                   /* but no change to <dsq> target... */
	      gm_restricted  = gm;	                                                  /* nor to the profile */

	      if ( p7_sparsemask_Validate(sm_restricted, errbuf) != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
	    }
	  else if (do_zeroed_emissions)
	    {                                                                    /* In the zeroed-emission utest: */
	      sm_restricted = sm;                                                /* no change to the sparse mask... */
	      create_envdoctored_dsq(sq, iae, ibe, &dsq_restricted);             /* but seq has all A outside, no C inside envelope */
	      create_envdoctored_profile(hmm, gm, bg, kae, kbe, &gm_restricted); /* and profile has ek(C)=1 outside, ek(AC)=0 inside envelope */
	    }

	  /* Now recalculate on the utest's envelope-restricted versions: */
	  if (p7_SparseForward(dsq_restricted, sq->n, gm_restricted, sm_restricted, sx, &f2_sc) != eslOK) esl_fatal(msg);
	  //p7_sparsemx_Dump(stdout, sx);
	  a2_status = p7_SparseEnvscoreApprox(gm_restricted, sx, iae, ibe, &a2_sc);
	  if (a2_status != eslOK && a2_status != eslEINACCURATE) esl_fatal(msg);
	  p7_sparsemx_Reuse(sx);
	  if (p7_SparseEnvscore(dsq_restricted, sq->n, gm_restricted, iae, ibe, kae, kbe, sm_restricted, sx, &e2_sc) != eslOK) esl_fatal(msg);
	  //p7_sparse_envscore_Dump(stdout, sx, iae, ibe, kae, kbe);
	  
	  if (do_unihit) 
	    {
	      if (esl_FCompare(e1_sc, f2_sc, /*rtol=*/0.0, tol)                              != eslOK) esl_fatal("%s\n e1 (%f) should equal f2 (%f)",             msg, e1_sc, f2_sc);
	      if (esl_FCompare(e1_sc, e2_sc, /*rtol=*/0.0, tol)                              != eslOK) esl_fatal("%s\n e1 (%f) should equal e2 (%f)",             msg, e1_sc, e2_sc);
	      if (esl_FCompare(e1_sc, a2_sc, /*rtol=*/0.0, (a2_status == eslOK? tol : tol2)) != eslOK) esl_fatal("%s\n e1 (%f) should equal a2 (%f)",             msg, e1_sc, a2_sc);
	      if ( e1_sc > a1_sc + (a1_status == eslOK? tol: tol2))                                    esl_fatal("%s\n  (e1,f2,e2,a2) (%f) should be <= a1 (%f)", msg, e1_sc, a1_sc);
	    }
	  else           
	    {
	      if ( esl_FCompare(e1_sc, e2_sc, /*rtol=*/0.0, tol)                              != eslOK) esl_fatal("%s\n  e1 (%f) should equal e2 (%f)\n",         msg, e1_sc, e2_sc);
	      if ( esl_FCompare(f2_sc, a2_sc, /*rtol=*/0.0, (a2_status == eslOK? tol : tol2)) != eslOK) esl_fatal("%s\n  f2 (%f) should equal a2 (%f)\n",         msg, f2_sc, a2_sc);
	      if ( e2_sc > f2_sc + tol)                                                                 esl_fatal("%s\n  (e1,e2) (%f) should be <= (f2,a2) (%f)", msg, e2_sc, f2_sc);
	      if ( f2_sc > a1_sc + (a1_status == eslOK? tol : tol2))                                    esl_fatal("%s\n  (f2,a2) (%f) should be <= a1 (%f)",      msg, f2_sc, a1_sc);
	    }

	  if (do_intersected_mask) 
	    {
	      p7_sparsemask_Destroy(sm_restricted);
	    }
	  else
	    {
	      free(dsq_restricted);
	      p7_profile_Destroy(gm_restricted);
	    }
	  p7_sparsemx_Reuse(sx);
	  p7_masstrace_Reuse(mt);
	}
      
      p7_trace_Reuse(vtr);
      p7_sparsemx_Reuse(sxd);
      p7_sparsemx_Reuse(sxb);
      p7_sparsemx_Reuse(sxf);
      p7_sparsemask_Reuse(sm);
    }

  p7_masstrace_Destroy(mt);
  p7_trace_Destroy(vtr);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sx);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(ox);
  esl_sq_Destroy(sq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}

/* Finally, the unit tests.
 */
static void 
utest_intersected_mask(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int do_unihit, int do_glocal)
{ 
  char msg[] = "sparse envscore, intersected-mask test failed";
  utest_engine(TRUE, msg, rng, abc, bg, M, L, N, do_unihit, do_glocal);
}
static void 
utest_zeroed_emissions(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int do_unihit, int do_glocal)
{
  char msg[] = "sparse envscore, zeroed-emissions test failed";
  utest_engine(FALSE, msg, rng, abc, bg, M, L, N, do_unihit, do_glocal);
}
#endif /*p7SPARSE_ENVSCORE_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/


/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7SPARSE_ENVSCORE_TESTDRIVE
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,     "70", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,     "35", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,      "5", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for sparse envelope scores";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg   = p7_bg_Create(abc);
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  utest_zeroed_emissions(r, abc, bg, M, L, N, /*do_unihit=*/FALSE, /*do_glocal=*/FALSE); 
  utest_zeroed_emissions(r, abc, bg, M, L, N, /*do_unihit=*/TRUE,  /*do_glocal=*/FALSE); 
  utest_zeroed_emissions(r, abc, bg, M, L, N, /*do_unihit=*/FALSE, /*do_glocal=*/TRUE); 
  utest_zeroed_emissions(r, abc, bg, M, L, N, /*do_unihit=*/TRUE,  /*do_glocal=*/TRUE); 

  utest_intersected_mask(r, abc, bg, M, L, N, /*do_unihit=*/FALSE, /*do_glocal=*/FALSE); 
  utest_intersected_mask(r, abc, bg, M, L, N, /*do_unihit=*/TRUE,  /*do_glocal=*/FALSE); 
  utest_intersected_mask(r, abc, bg, M, L, N, /*do_unihit=*/FALSE, /*do_glocal=*/TRUE); 
  utest_intersected_mask(r, abc, bg, M, L, N, /*do_unihit=*/TRUE,  /*do_glocal=*/TRUE); 

  fprintf(stderr, "#  status = ok\n");

  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7SPARSE_ENVSCORE_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/


/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef p7SPARSE_ENVSCORE_EXAMPLE

#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",              0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of envelope scoring, sparse dual implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_HMM         *hmm     = NULL;
  ESL_SQ         *sq      = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_CHECKPTMX   *ox      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  P7_SPARSEMX    *sx      = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxf     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxb     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxd     = p7_sparsemx_Create(NULL);
  P7_TRACE       *vtr     = p7_trace_CreateWithPP();
  P7_MASSTRACE   *mt      = NULL;
  P7_SPARSEMASK  *sm_intersected = NULL;
  float           fsc, vsc, bsc;
  int             iae, ibe, kae, kbe;
  float           e1_sc, a1_sc, f2_sc, e2_sc, a2_sc;
  int             d;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Open sequence database */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Read in one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);
  esl_sqfile_Close(sqfp);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_oprofile_Convert(gm, om);

  /* Set the profile and null model's target length models */
  p7_bg_SetLength           (bg, sq->n);
  p7_profile_SetLength      (gm, sq->n);
  p7_oprofile_ReconfigLength(om, sq->n);

  /* Use f/b filter to create sparse mask */
  ox = p7_checkptmx_Create(hmm->M, sq->n, ESL_MBYTES(32));
  sm = p7_sparsemask_Create(gm->M, sq->n);
  p7_ForwardFilter (sq->dsq, sq->n, om, ox, /*fsc=*/NULL);
  p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH);
  
  /* Sparse DP calculations */
  p7_SparseViterbi (sq->dsq, sq->n, gm, sm, sx, vtr, &vsc);
  p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf,     &fsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxb,     &bsc);
  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxb, sxd);

  p7_sparsemx_TracePostprobs(sxd, vtr); /* selecting domain anchors requires pp annotation of the trace */
  p7_trace_Index(vtr);
  p7_sparsemx_Reuse(sx);

  /* Look at each domain: */
  mt = p7_masstrace_CreateSlim(gm->M, sq->n);
  for (d = 0; d < vtr->ndom; d++)
    {
      p7_SparseMasstrace(sq->dsq, sq->n, gm, sxf, sxb, vtr, vtr->anch[d], 0.1, sx, mt, &iae, &ibe, &kae, &kbe);
      p7_sparsemx_Reuse(sx);

      p7_SparseEnvscoreApprox(gm, sxf, iae, ibe, &a1_sc);
      p7_SparseEnvscore(sq->dsq, sq->n, gm, iae, ibe, kae, kbe, sm, sx, &e1_sc);
      p7_sparsemx_Reuse(sx);

      p7_sparse_envscore_IntersectedMask(sm, iae, ibe, kae, kbe, &sm_intersected);

      p7_SparseForward(sq->dsq, sq->n, gm, sm_intersected, sx, &f2_sc);
      p7_SparseEnvscoreApprox(gm, sx, iae, ibe, &a2_sc);
      p7_sparsemx_Reuse(sx);

      p7_SparseEnvscore(sq->dsq, sq->n, gm, iae, ibe, kae, kbe, sm_intersected, sx, &e2_sc);
      p7_sparsemx_Reuse(sx);

      printf("Domain %d:\n", d);
      printf("   ali ia..ib            = %d..%d\n", vtr->sqfrom[d],  vtr->sqto[d]);
      printf("   ali ka..kb            = %d..%d\n", vtr->hmmfrom[d], vtr->hmmto[d]);
      printf("   anchor i0,k0,st0      = %d,%d,%s; %.4f\n", vtr->i[vtr->anch[d]], vtr->k[vtr->anch[d]], p7_trace_DecodeStatetype(vtr->st[vtr->anch[d]]), vtr->pp[vtr->anch[d]]);
      printf("   env iae..ibe          = %d..%d\n", iae, ibe);
      printf("   env kae..kbe          = %d..%d\n", kae, kbe);
      printf("   envelope score (e1)   = %.6f\n", e1_sc);
      printf("   approx env score (a1) = %.6f\n", a1_sc);
      printf("   intersected, f2       = %.6f\n", f2_sc);
      printf("   intersected, e2       = %.6f\n", e2_sc);
      printf("   intersected, a2       = %.6f\n", a2_sc);

      p7_masstrace_Reuse(mt);
      p7_sparsemask_Destroy(sm_intersected);
    }


  /* Cleanup */
  esl_sq_Destroy(sq);
  p7_masstrace_Destroy(mt);
  p7_trace_Destroy(vtr);
  p7_sparsemx_Destroy(sx);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(ox);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}

#endif /*p7SPARSE_ENVSCORE_EXAMPLE*/
/*---------------- end, example main()  ------------------------*/
