/* Production implementation of anchor set constrained (ASC) sparse
 * dynamic programming.
 *
 * Contents:
 *   1. Sparse ASC Forward
 *   2. Sparse ASC Forward, segmental version
 */


/*****************************************************************
 * 1. Sparse ASC Forward 
 *****************************************************************/

/* Function:  h4_sparse_asc_Forward()
 * Synopsis:  Sparse anchor set constrained Forward algorithm.
 *
 * Purpose:   Run the sparse anchor set constrained (ASC) Forward
 *            algorithm, comparing profile <hmm> in comparison mode
 *            <mo> to target sequence <dsq> of length <L>, dually
 *            constrained by anchor set <anch> and sparse mask
 *            <sm>. Fills in the sparse ASC Forward matrix <asf> and
 *            optionally returns the sparse raw ASC Forward bitscore
 *            in <opt_sc>.
 *
 * Args:      dsq    : target digital sequence, 1..L
 *            L      : length of <dsq>
 *            hmm    : query profile
 *            mo     : comparison mode
 *            anch   : anchor set 
 *            sm     : sparse mask that constrains the DP
 *            asf    : RESULT : sparse ASC Forward DP matrix (reused/realloc'ed here)
 *            opt_sc : optRETURN : sparse ASC Forward score
 *
 * Returns:   <eslOK> on success.
 *  
 * Throws:    <eslEMEM> on allocation failure.
 */
int
h4_sparse_asc_Forward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo,
                      const H4_ANCHORSET *anch,  const H4_SPARSEMASK *sm, H4_SPARSEMX *asf, float *opt_sc)
{
  float *dpc;
  float *xc;
  float xN  = 0.;
  float xJ  = -eslINFINITY;
  float xC  = -eslINFINITY;
  int   d   = 1;
  int   g;
  int   ngap;
  int   status;

  /* Argument validation */
  ESL_DASSERT1(( sm->L == L     ));
  ESL_DASSERT1(( sm->M == hmm->M ));

  /* Reallocation, if needed */
  if (( status = h4_spascmx_Reinit(asf, anch, sm)) != eslOK) return status;
  asf->type = h4S_ASC_FWD;

  /* Iterate the ForwardSeg algorithm over each segment in turn. 
   */
  dpc = asf->dp;
  xc  = asf->xmx;
  for (g = 1; g <= sm->S; g++)
    {
      ngap = (g == 1 ? sm->seg[g].ia-1 : sm->seg[g].ia - sm->seg[g-1].ib - 1);
      h4_sparse_asc_ForwardSeg(dsq, L, gm, anch, D, d, sm, ngap, sm->seg[g].ia, sm->seg[g].ib, 
			       asf, xN, xJ, xC, dpc, xc, 
			       &xN, &xJ, &xC, &d, &dpc, &xc, NULL);
    }
  
  /* Account for the last ib(S)+1..L residues in ...CCCC path.
   */
  if (L > sm->seg[sm->S].ib) 
    xC +=  (float) (L - sm->seg[sm->S].ib) * gm->xsc[p7P_C][p7P_LOOP];
  if (opt_sc) *opt_sc = xC + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}
/*--------------- end, sparse ASC Forward -----------------------*/

