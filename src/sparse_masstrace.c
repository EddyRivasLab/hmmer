

/* 
 *  The probability mass that flows through state X at position i,k is calculated recursively:
 *  
 *  \rho^X(i,k) = \sum_Y \rho^Y(i_Y,k_Y) * 
 *                          exp \left(  \alpha^X(i,k) + \tau_k(XY) + \epsilon(Y_k_Y, x_i_Y) - \alpha^Y(i_Y, k_Y)} \right)
 *  where the notation is as follows:
 *    X = current state
 *    Y = state that X can transition to (MID only; exclude E)
 *    i_Y,k_Y    = condensed notation for next i,k. i_M,k_M = i+1,k+1; i_I,k_I = i+1,k; i_D,k_D = i,k+1.
 *    \alpha()   = Forward matrix, log space in nats
 *    \tau_k(XY) = log transition prob, for example log t(Mk->Mk+1) is \tau_k(MM)
 *    \epsilon() = emission scores, for example e(Mk+1, x_i+1)
 *    
 *  A little useful guidance:
 *    
 *    this thing: \alpha^X(i,k) + \tau_k(XY) + \epsilon(Y_k_Y, x_i_Y)
 *    is the XY component (edge) of the Forward recursion into state Y, accounting for the X->Y state transition.
 *    
 *    this thing: \alpha^Y(i_Y, k_Y)} 
 *    is the sum of all such edges into Y: i.e., the Forward recursion.
 *    
 *    if these were in probability space, the ratio of one edge into Y
 *    over all edges into Y gives us the _edge weight_: the fraction
 *    of posterior probability mass that flows back that edge. The sum
 *    of edge weights into Y is 1. Since everything's in null-scaled
 *    log space, the ratio shows up as a subtraction and an exp().
 *    
 *    For a Y=Mk match state, one such edge weight is the B->Mk edge;
 *    mass that flows back that edge disappears from our recursion.
 *    Essentially, it's this lossage that we're tracking back: how far
 *    we have to trace back in i,k before enough mass is lost from the
 *    main state alignment ensemble that the mass we're still tracing
 *    back drops below the <massthresh> threshold.
 *    
 *    In a stochastic traceback, we would trace back from Y
 *    stochastically, according to its edge weight distribution
 *    (inclusive of the B->Mk edges). Notationally, we implement a
 *    stochastic traceback as a _push_ from a current state Y to a
 *    next state X, whereas the mass traceback is implemented as a
 *    _pull_ to a current state X from connected states Y.
 */


int
p7_sparse_masstrace_Backward(const P7_SPARSEMX *fwd, P7_SPARSEMX *mass, P7_TRACE *tr, int z, float massthresh, int *ret_iae, int *ret_kae)
{
  const P7_SPARSEMASK *sm = fwd->sm;
  int   st0 = tr->st[z];
  int   k0  = tr->k[z];
  int   i0  = tr->i[z];
  float pp0 = tr->pp[z];
  int   iae = i0;
  int   kae = k0;
  float *rhoc;			/* ptr into current row i rho (trace mass) cells */
  float *rhon;			/* ptr into next row i+1 rho (trace mass) cells */
  const float *dpc;		/* ptr into current row i Forward DP cells */
  const float *dpn;		/* ptr into next row i+1 Forward DP cells */
  float *last_rhoc;
  const float *last_dpc;
  const float *rsn;		/* residue score vector for next row i+1. Enables MSN(k) notation macro. */
  const float *tsc = gm->tsc;   /* transition score vector. Enables TSC(k) notation macro. */
  int   i, k, v,w,s;
  float  rowmass;
  float *kmass           = NULL;
  int    iae_proven_done = FALSE;
  int    kae_proven_done = FALSE;

  /* we should avoid this tmp alloc. find space for kmass in sparsemx. */
  ESL_ALLOC(kmass, sizeof(float) * (sm->M));
  for (k = 0; k<= sm->M; k++) 
    kmass[k] = 0.;

  dpc  = last_dpc  = fwd->dp  + p7S_NSCELLS*((sm->k[i] + sm->n[i] - 1) - sm->kmem);      // unwarranted pointer-arithmetic chumminess with sparsemask, sparsemx: assumes that dp[] and kmem[] are EXACTLY 1:1 with each other
  rhoc = last_rhoc = mass->dp + p7S_NSCELLS*((sm->k[i] + sm->n[i] - 1) - sm->kmem);      // sm->k[i+1]-1 doesn't work because of i=L edge case (no k[L+1]) 

  /* ok, so, the following is insane.
   * note that local, glocal paths don't cross; 
   * if st0 is {MDI}L, only {MDI}L can be reached; analogous for {MDI}G.
   * so we only need to compute the L or G interleaved half of the sparse matrix, depending on st0
   * so...
   * the supercells are [ ML MG IL IG DL DG ]
   * dpc, rhoc currently point at the final supercell on row i: specifically at its ML cell
   * so... 
   * if st0 is {MDI}G glocal, offset ptrs by +1 cell; now our ptrs are on MG 
   * because of relative indexing in dp (everything done by stepping exactly in supercell units),
   * this single offset suffices to keep dpc,dpn,rhoc,rhon on L vs. G in all the rest of the code here
   * we access everything with L indices, but with the offset, these will exactly be the G cells
   * seriously, trust me.
   */
  if (st0 == p7T_MG || st0 == p7T_IG || st0 == p7T_DG)
    { dpc += 1; rhoc += 1; }

  /* special case the first row (i0): initialize on i0,k0, then pull delete path to the left. */
  // first, find k0, bravely assuming that k0 MUST be in the sparse list for i0
  for (v = sm->n[i0]-1; sm->k[i][v] != k0; v--)
    {
      rhoc[p7S_ML] = rhoc[p7S_IL] = rhoc[p7S_DL] = 0.;
      dpc  -= p7S_NSCELLS; 
      rhoc -= p7S_NSCELLS; 
    }

  /* now v is the index of k0 in row i0's sparse cell list;
   *     *dpc is the sparse DP supercell \alpha(i0,k0) for X = {ML MG IL IG DL DG}
   *     *rho is \rho(i0,k0) supercell
   */
  switch (st0) {
  case p7T_ML: case p7T_MG: rhoc[p7S_ML] = pp0; rhoc[p7S_IL] = 0.;  rhoc[p7S_DL] = 0.;  break;
  case p7T_IL: case p7T_IG: rhoc[p7S_ML] = 0.;  rhoc[p7S_IL] = pp0; rhoc[p7S_DL] = 0.;  break;
  case p7T_DL: case p7T_DL: rhoc[p7S_ML] = 0.;  rhoc[p7S_IL] = 0.;  rhoc[p7S_DL] = pp0; break;
  default:     ESL_EXCEPTION(eslEINCONCEIVABLE, "you know it is");
  }
  kmass[k0] += pp0;
  dpc       -= p7S_NSCELLS;
  rhoc      -= p7S_NSCELLS;
  
  /* now pull to the left. If we didn't start on a D, or
   * if we don't have contiguous supercells on the row, this
   * is all an expensive way of zeroing the row: but it's
   * clearer this way than a lot of special case branching,
   * especially since it's obvious where the same ->D case is
   * in the main recursion later.
   */ 
  for (v = v-1; v >= 0; v--) 
    {
      k = sm->k[i][v];
      if (sm->k[i][v+1] == k+1) {
	rhoc[p7S_ML] = rhoc[p7S_ML+p7S_NSCELLS] * exp( dpc[p7S_ML] + TSC(p7P_MD, k) - dpc[p7S_DL+p7S_NSCELLS]);
	rhoc[p7S_IL] = 0.0f;
	rhoc[p7S_DL] = rhoc[p7S_DL+p7S_NSCELLS] * exp( dpc[p7S_DL] + TSC(p7P_DD, k) - dpc[p7S_DL+p7S_NSCELLS]);
      } else { rhoc[p7S_ML] = rhoc[p7S_IL] = rhoc[p7S_DL] = 0.0f; }

      kmass[k] += rhoc[p7S_ML] + rhoc[p7S_DL];
      if (kmass[k] >= massthresh) kae = k; 
      dpc      -= p7S_NSCELLS;
      rhoc     -= p7S_NSCELLS;
    }

  
  /* The main recursion.  */
  for (i = i0-1; i >= 1; i--)
    {
      dpn     = last_dpc;
      rhon    = last_rhoc;
      rsn     = gm->rsc[dsq[i+1]];  // MSN() notation macro now valid
      rowmass = 0.;

      last_dpc  = dpc;
      last_rhoc = rhoc;

      w = sm->n[i+1]-1;
      v = sm->n[i] - 1; while (sm->k[i][v] > k0) { v--; dpc -= p7S_NSCELLS; rhoc -= p7S_NSCELLS; }
      if (v < 0) break;		/* no cells on row at all? trace mass can't flow back any further then; we're done for sure. */
      for (; v >= 0; v--)	/* for all sparse k on row, such that k <= k0. if no such k exist, this code doesn't execute, and mass flow is done */
	{
	  k = sm->k[i][v];

	  while (w >= 0 && sm->k[i+1][w]  > k+1) { w--; dpn -= p7S_NSCELLS; rhon -= p7S_NSCELLS; }
	  if    (w >= 0 && sm->k[i+1][w] == k+1) {
	    rhoc[p7S_ML]  = rhon[p7S_ML] * exp( dpc[p7S_ML] + TSC(p7P_MM, k) + MSN(k+1) - dpn[p7S_ML]);
	    rhoc[p7S_IL]  = rhon[p7S_ML] * exp( dpc[p7S_IL] + TSC(p7P_IM, k) + MSN(k+1) - dpn[p7S_ML]);
	    rhoc[p7S_DL]  = rhon[p7S_ML] * exp( dpc[p7S_DL] + TSC(p7P_DM, k) + MSN(k+1) - dpn[p7S_ML]);
	  } else { rhoc[p7S_ML] = rhoc[p7S_IL] = rhoc[p7S_DL] = 0.0f; }

	  while (w >= 0 && sm->k[i+1][w]  > k) { w--; dpn -= p7S_NSCELLS; rhon -= p7S_NSCELLS; }
	  if    (w >= 0 && sm->k[i+1][w] == k) { 
	    rhoc[p7S_ML] += rhon[p7S_IL] * exp( dpc[p7S_ML] + TSC(p7P_MI, k) - dpn[p7S_IL]); // insert scores ISN(k) assumed to be zero
	    rhoc[p7S_IL] += rhon[p7S_IL] * exp( dpc[p7S_IL] + TSC(p7P_II, k) - dpn[p7S_IL]);
	  }

	  if (v < sm->n[i]-1 && sm->k[i][v+1] == k+1) {
	    rhoc[p7S_ML] += rhoc[p7S_ML+p7S_NSCELLS] * exp( dpc[p7S_ML] + TSC(p7P_MD, k) - dpc[p7S_DL+p7S_NSCELLS]);
	    rhoc[p7S_DL] += rhoc[p7S_DL+p7S_NSCELLS] * exp( dpc[p7S_DL] + TSC(p7P_DD, k) - dpc[p7S_DL+p7S_NSCELLS]);
	  }

	  kmass[k] += rhoc[p7S_ML] + rhoc[p7S_DL]; /* kmass[k] is a lower bound on how much probability mass is flowing leftwards thru this column  */
	  if (k < kae && kmass[k] >= massthresh) kae = k; 
	  if (kmass[k] + rowmass < massthresh) kae_proven_done = TRUE; /* kmass[k] + rowmass is the upper bound on what can flow leftward thru k */
	  rowmass  += rhoc[p7S_ML] + rhoc[p7S_IL]; /* how much total probability mass is flowing upwards through this row  */

	  rhoc -= p7S_NSCELLS;
	  dpc  -= p7S_NSCELLS;
	}

      if (rowmass < massthresh)  iae_proven_done = TRUE;
      if (iae_proven_done && kae_proven_done) break;
      iae = i;
    }
  
  free(kmass);
  *ret_iae = iae;
  *ret_kae = kae;
  return eslOK;

 ERROR:
  if (kmass) free(kmass);
  return status;
}

