/* Mass trace algorithm, sparse DP version:
 * finding bounds of a domain envelope.
 * 
 * Contents:
 *   1. Upwards recursion (iae/kae starts)
 *   2. Downwards recursion (ibe/kbe ends)
 *   3. API, wrapper for up/down recursions
 *   4. Unit tests
 *   5. Test driver
 *   6. Example
 *   7. Copyright, license information.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "base/p7_profile.h"
#include "base/p7_masstrace.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_masstrace.h"

/*****************************************************************
 * 1. The upwards recursion.
 *****************************************************************/

static int
masstrace_up(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *fwd,
	     int i0, int k0, int st0, float massthresh, 
	     P7_SPARSEMX *mass, P7_MASSTRACE *mt, int *ret_iae, int *ret_kae)
{
  const P7_SPARSEMASK *sm = fwd->sm;
  int   iae = i0;		/* RETURN: envelope start coord in sequence (1..i0) */
  int   kae = k0;		/* RETURN: envelope start coord in profile (1..k0) */
  float *rhoc;			/* ptr into current row i rho (trace mass) cells */
  float *rhon;			/* ptr into next row i+1 rho (trace mass) cells */
  const float *dpc;		/* ptr into current row i Forward DP cells */
  const float *dpn;		/* ptr into next row i+1 Forward DP cells */
  float *last_rhoc;
  const float *last_dpc;
  const float *rsn;		/* residue score vector for next row i+1. Enables MSN(k) notation macro. */
  const float *tsc = gm->tsc;   /* transition score vector. Enables TSC(k) notation macro. */
  int   i, k, v,w;
  float  rowmass;
  int    iae_proven_done = FALSE;
  int    kae_proven_done = FALSE;

  /* dpc, rhoc initialized to point at last sparse supercell on row i0 */
  dpc  = fwd->dp  + p7S_NSCELLS*((sm->k[i0] + sm->n[i0] - 1) - sm->kmem);      // unwarranted pointer-arithmetic chumminess with sparsemask, sparsemx: assumes that dp[] and kmem[] are EXACTLY 1:1 with each other
  rhoc = mass->dp + p7S_NSCELLS*((sm->k[i0] + sm->n[i0] - 1) - sm->kmem);      // sm->k[i+1]-1 doesn't work because of i=L edge case (no k[L+1]) 

  /* ok, so, the following is not pretty.
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
   * we access everything with L indices, but with the offset, these will be exactly the G cells!
   * seriously, trust me.
   */
  if (st0 == p7T_MG || st0 == p7T_IG || st0 == p7T_DG) { dpc += 1; rhoc += 1; }
  last_dpc  = dpc;
  last_rhoc = rhoc;

  /* Another wrinkle: on glocal paths, we know kae=1 and kbe=M.
   * The code below is going to neglect the prob mass that flows
   * through wing-retracted DGks; so for glocal mass trace, it
   * will calculate kmass[] distribution exclusive of the wing
   * retractions. On a glocal path, we must be sure not to modify kae
   * after this initialization.
   */
  if (p7_trace_IsGlocal(st0)) kae = 1;

  /* Initialization of first row, i0:
   */
  /* first, skip backwards to find k0, bravely assuming that k0 MUST be in the sparse list for i0 */
  for (v = sm->n[i0]-1; sm->k[i0][v] != k0; v--)
    { dpc  -= p7S_NSCELLS; rhoc -= p7S_NSCELLS;  }

  /* now v is the index of k0 in row i0's sparse cell list;
   *     *dpc is the sparse DP supercell \alpha(i0,k0) for X = {ML MG IL IG DL DG}
   *     *rho is \rho(i0,k0) supercell
   */
  switch (st0) {
  case p7T_ML: case p7T_MG: rhoc[p7S_ML] = 1.; rhoc[p7S_IL] = 0.;  rhoc[p7S_DL] = 0.;  break;
  case p7T_IL: case p7T_IG: rhoc[p7S_ML] = 0.; rhoc[p7S_IL] = 1.;  rhoc[p7S_DL] = 0.;  break;
  case p7T_DL: case p7T_DG: rhoc[p7S_ML] = 0.; rhoc[p7S_IL] = 0.;  rhoc[p7S_DL] = 1.;  break;
  default:     ESL_EXCEPTION(eslEINCONCEIVABLE, "you know it is");
  }
  mt->kmass[k0]  = 1.0;
  if (mt->imass) mt->imass[i0] = 1.0;
  dpc  -= p7S_NSCELLS;
  rhoc -= p7S_NSCELLS;
  
  /* now pull to the left to finish the initialization of row i0.  If
   * we didn't start on a D, or if we don't have contiguous supercells
   * on the row, this amounts to an expensive way of zeroing the row:
   * but it's clearer this way than a lot of special case branching,
   * especially since it's obvious where the same ->D case is in the
   * main recursion later.
   */ 
  for (v = v-1; v >= 0; v--) 
    {
      k = sm->k[i0][v];
      if (sm->k[i0][v+1] == k+1) { /* if v,v+1 sparse cells are contiguous k,k+1, we can propagate D path leftwards */
	rhoc[p7S_ML] = P7_MASSTRACE_INCREMENT(rhoc[p7S_DL+p7S_NSCELLS],  dpc[p7S_ML] + TSC(p7P_MD, k), dpc[p7S_DL+p7S_NSCELLS]);
	rhoc[p7S_IL] = 0.0f;
	rhoc[p7S_DL] = P7_MASSTRACE_INCREMENT(rhoc[p7S_DL+p7S_NSCELLS],  dpc[p7S_DL] + TSC(p7P_DD, k), dpc[p7S_DL+p7S_NSCELLS]);
      } else { rhoc[p7S_ML] = rhoc[p7S_IL] = rhoc[p7S_DL] = 0.0f; } /* else we can't, and no probability mass reaches these cells (nor any others on the row) */

      mt->kmass[k] += rhoc[p7S_ML] + rhoc[p7S_DL];
      if (k < kae && mt->kmass[k] >= massthresh) kae = k; /* k<kae test suffices to prevent changing kae on glocal path */
      dpc  -= p7S_NSCELLS;
      rhoc -= p7S_NSCELLS;
    }
  /* Now dpc, rhoc are on last supercell of row i0-1. */
  
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
      v = sm->n[i] - 1;
      while (v >= 0 && sm->k[i][v] > k0) { v--; dpc -= p7S_NSCELLS; rhoc -= p7S_NSCELLS; }
      if    (v < 0) break;	/* no cells on row at all? trace mass can't flow back any further then; we're done for sure. */
      for (; v >= 0; v--)	/* for all sparse k on row, such that k <= k0. if no such k exist, this code doesn't execute, and mass flow is done */
	{
	  k = sm->k[i][v];

	  /* Try to find the M(i+1,k+1) cell on row i+1. If it exists: apportion its mass to our current cells {MID}ik */
	  while (w >= 0 && sm->k[i+1][w]  > k+1) { w--; dpn -= p7S_NSCELLS; rhon -= p7S_NSCELLS; }
	  if    (w >= 0 && sm->k[i+1][w] == k+1 && k < k0) {  // note k<k0 test; for k=k0, Mk+1 was never initialized in rhon[]
	    rhoc[p7S_ML]  = P7_MASSTRACE_INCREMENT(rhon[p7S_ML], dpc[p7S_ML] + TSC(p7P_MM, k) + MSN(k+1), dpn[p7S_ML]);
	    rhoc[p7S_IL]  = P7_MASSTRACE_INCREMENT(rhon[p7S_ML], dpc[p7S_IL] + TSC(p7P_IM, k) + MSN(k+1), dpn[p7S_ML]);
	    rhoc[p7S_DL]  = P7_MASSTRACE_INCREMENT(rhon[p7S_ML], dpc[p7S_DL] + TSC(p7P_DM, k) + MSN(k+1), dpn[p7S_ML]);
	  } else { rhoc[p7S_ML] = rhoc[p7S_IL] = rhoc[p7S_DL] = 0.0f; }

	  /* Try to find I(i+1,k) cell on row i+1; if exists, apportion its mass */
	  while (w >= 0 && sm->k[i+1][w]  > k) { w--; dpn -= p7S_NSCELLS; rhon -= p7S_NSCELLS; }
	  if    (w >= 0 && sm->k[i+1][w] == k) { 
	    rhoc[p7S_ML] += P7_MASSTRACE_INCREMENT(rhon[p7S_IL], dpc[p7S_ML] + TSC(p7P_MI, k), dpn[p7S_IL]); // insert scores ISN(k) assumed to be zero
	    rhoc[p7S_IL] += P7_MASSTRACE_INCREMENT(rhon[p7S_IL], dpc[p7S_IL] + TSC(p7P_II, k), dpn[p7S_IL]);
	  }

	  /* If v+1 is k+1 ... then v,v+1 are contiguous supercells ... and D(i,k+1) exists, so apportion its mass. */
	  if (v < sm->n[i]-1 && sm->k[i][v+1] == k+1 && k < k0) { // k<k0 test: don't look for Dk0+1
	    rhoc[p7S_ML] += P7_MASSTRACE_INCREMENT(rhoc[p7S_DL+p7S_NSCELLS], dpc[p7S_ML] + TSC(p7P_MD, k), dpc[p7S_DL+p7S_NSCELLS]);
	    rhoc[p7S_DL] += P7_MASSTRACE_INCREMENT(rhoc[p7S_DL+p7S_NSCELLS], dpc[p7S_DL] + TSC(p7P_DD, k), dpc[p7S_DL+p7S_NSCELLS]);
	  }

	  mt->kmass[k] += rhoc[p7S_ML] + rhoc[p7S_DL]; /* kmass[k] is a lower bound on how much probability mass is flowing leftwards thru this column  */
	  if (k < kae && mt->kmass[k] >= massthresh) kae = k; /* k<kae condition suffices to prevent changing preset kae=1 on glocal path */
	  if (kae == 1 || mt->kmass[k] + rowmass < massthresh) kae_proven_done = TRUE; /* kmass[k] + rowmass is the upper bound on what can flow leftward thru k */
	  rowmass  += rhoc[p7S_ML] + rhoc[p7S_IL]; /* accumulate total probability mass that's still flowing upwards through this row  */

	  rhoc -= p7S_NSCELLS;
	  dpc  -= p7S_NSCELLS;
	}
      if (rowmass < massthresh)  iae_proven_done = TRUE; else iae = i;     /* we keep decrementing iae until we drop below our mass threshold */

      if      (mt->imass)                          mt->imass[i] = rowmass; /* if we're storing <imass>, we do all rows. */
      else if (iae_proven_done && kae_proven_done) break;                  /* if not, we stop recursing up in i when we've provably obtained iae,kae */
    }

  *ret_iae   = iae;
  *ret_kae   = kae;
  return eslOK;
}
/*--------------- end, upwards recursion ------------------------*/


/*****************************************************************
 * 2. The downwards recursion.
 *****************************************************************/

static int
masstrace_down(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *bck,
	       int i0, int k0, int st0, float massthresh,
	       P7_SPARSEMX *mass, P7_MASSTRACE *mt, int *ret_ibe, int *ret_kbe)
{
  const P7_SPARSEMASK *sm = bck->sm; 
  int          ibe = i0;	/* RETURN: envelope end coord in sequence (i0..L) */
  int          kbe = k0;	/* RETURN: envelope end coord in profile (10..M) */
  float       *rhoc;  		/* ptr that steps through sparse mass supercells on current row i */
  float       *rhop;		/* ptr that steps through sparse mass supercells on previous row i-1 */
  float       *last_rhoc;	/* used to set the next rhop */
  const float *dpc;		/* ptr that steps through sparse Forward supercells on current row i */
  const float *dpp;		/* ptr that steps through sparse Forward supercells on previous row i-1 */
  const float *last_dpc;	/* used to set the next dpp */
  const float *rsc;		/* residue score vector on current row. Enables MSC(k) notation macro */
  const float *tsc   = gm->tsc; /* transition score vector. Enables TSC(k) notation macro */
  float        rowmass;
  int ibe_proven_done = FALSE;
  int kbe_proven_done = FALSE;
  int v;			/* v index steps through list of sparse supercells k on current row i    */
  int w;			/* w index steps through list of sparse supercells k on previous row i-1 */
  int k;			/* index in profile positions 1..M */
  int i;			/* index in sequence positions 1..L */

  /* dpc, rhoc set to first sparse supercell on row i0. */
  dpc  = bck->dp  + p7S_NSCELLS*(sm->k[i0] - sm->kmem);  // unwarranted ptr-arithmetic chumminess with sparse matrix layout
  rhoc = mass->dp + p7S_NSCELLS*(sm->k[i0] - sm->kmem);
  
  /* See comment in _up() about the following magic: */
  if (st0 == p7T_MG || st0 == p7T_IG || st0 == p7T_DG)
    { dpc += 1; rhoc += 1; }
  last_dpc  = dpc;  /* these two lines must be placed after the magic above, of course */
  last_rhoc = rhoc;

  /* See comment in _up(): for glocal paths we know kbe=M; we're going to neglect wing-retracted exits thru DG in code below */
  if (p7_trace_IsGlocal(st0)) kbe = gm->M;

  /* Start initializing the first row, i0. 
   * Up() is responsible for k<k0. We are (Down() is) responsible for k>k0.
   * In both Up() and Down(), (s0,i0,k0) anchor itself is set to 1.0, and
   * it doesn't hurt to redo that here.
   */
  /* first, skip to k0. We trust that it exists in sparse cell list: caller said so. */
  for (v = 0; sm->k[i0][v] != k0; v++) { dpc += p7S_NSCELLS; rhoc += p7S_NSCELLS; }

  /* now v is the index of k0 in row i0's sparse cell list, and rhoc/dpc are on the
   * anchor cell
   */
  switch (st0) {
  case p7T_ML: case p7T_MG: rhoc[p7S_ML] = 1.; rhoc[p7S_IL] = 0.;  rhoc[p7S_DL] = 0.;  break;
  case p7T_IL: case p7T_IG: rhoc[p7S_ML] = 0.; rhoc[p7S_IL] = 1.;  rhoc[p7S_DL] = 0.;  break;
  case p7T_DL: case p7T_DG: rhoc[p7S_ML] = 0.; rhoc[p7S_IL] = 0.;  rhoc[p7S_DL] = 1.;  break;
  default:     ESL_EXCEPTION(eslEINCONCEIVABLE, "you know it is");
  }
  mt->kmass[k0]  = 1.0;
  if (mt->imass) mt->imass[i0] = 1.0;
  dpc       += p7S_NSCELLS;
  rhoc      += p7S_NSCELLS;

  /* now pull to the right, on the delete path.
   * if we started on an I, this is just an expensive way to zero remaining sparse cells k>k0 
   */
  for (v = v+1; v < sm->n[i0]; v++) 
    {
      k = sm->k[i0][v];
      rhoc[p7S_ML] = rhoc[p7S_IL] = 0.;
      if (sm->k[i0][v-1] == k-1) 
	{
	  rhoc[p7S_DL] = 
	    P7_MASSTRACE_INCREMENT(rhoc[p7S_ML-p7S_NSCELLS], dpc[p7S_DL] + TSC(p7P_MD, k-1), dpc[p7S_ML-p7S_NSCELLS]) + // yes, those array indices are negative;
	    P7_MASSTRACE_INCREMENT(rhoc[p7S_DL-p7S_NSCELLS], dpc[p7S_DL] + TSC(p7P_DD, k-1), dpc[p7S_DL-p7S_NSCELLS]);  // but that's valid C,  seriously.
	}
      else rhoc[p7S_DL] = 0.;

      mt->kmass[k] += rhoc[p7S_ML] + rhoc[p7S_DL];
      if (k > kbe && mt->kmass[k] >= massthresh) kbe = k; /* k>kbe test suffices to avoid changing preset kbe=M on glocal path */
      dpc  += p7S_NSCELLS;
      rhoc += p7S_NSCELLS;
    }

  /* The main recursion */
  for (i = i0 + 1; i <= L; i++)
    {
      dpp     = last_dpc;
      rhop    = last_rhoc;
      rsc     = gm->rsc[dsq[i]];  // MSC(k) notation macro now valid
      rowmass = 0.;
      
      last_dpc  = dpc;
      last_rhoc = rhoc;

      v = w = 0;
      while (v <  sm->n[i] && sm->k[i][v] < k0) { v++; dpc += p7S_NSCELLS; rhoc += p7S_NSCELLS; }  // for Down pass, k >= k0
      if    (v == sm->n[i]) break;  // no cells on row i at all? trace mass can't flow any more. break completely out (of i loop), we're done.
      for (; v <  sm->n[i]; v++)
	{
	  k = sm->k[i][v];

	  /* Try to find supercell {MID}(i-1,k-1). If it exists, each cell apportions some of its mass to M(i,k) */
	  while (w < sm->n[i-1] && sm->k[i-1][w]  < k-1) { w++; dpp += p7S_NSCELLS; rhop += p7S_NSCELLS; }
	  if    (w < sm->n[i-1] && sm->k[i-1][w] == k-1 && k > k0) // k > k0 test is there to prevent looking at k0-1 supercell
	    {
	      rhoc[p7S_ML] = 
		P7_MASSTRACE_INCREMENT(rhop[p7S_ML], dpc[p7S_ML] + TSC(p7P_MM, k-1) + MSC(k), dpp[p7S_ML]) +
		P7_MASSTRACE_INCREMENT(rhop[p7S_IL], dpc[p7S_ML] + TSC(p7P_IM, k-1) + MSC(k), dpp[p7S_IL]) +
		P7_MASSTRACE_INCREMENT(rhop[p7S_DL], dpc[p7S_ML] + TSC(p7P_DM, k-1) + MSC(k), dpp[p7S_DL]);
	    }
	  else rhoc[p7S_ML] = 0.;

	  /* Try to find supercell {MID}(i-1,k), which is either w or w+1 if it exists; if so, its MI cells apportion mass to I(i,k) */
	  while (w < sm->n[i-1] && sm->k[i-1][w]  < k) { w++; dpp += p7S_NSCELLS; rhop += p7S_NSCELLS; }
	  if    (w < sm->n[i-1] && sm->k[i-1][w] == k && k < gm->M)  // k=M check because Im doesn't exist, and we need to avoid a -inf - -inf = nan
	    {
	      rhoc[p7S_IL] = 
		P7_MASSTRACE_INCREMENT(rhop[p7S_ML], dpc[p7S_IL] + TSC(p7P_MI, k), dpp[p7S_ML]) +   // here we're assuming ISC(k)=0
		P7_MASSTRACE_INCREMENT(rhop[p7S_IL], dpc[p7S_IL] + TSC(p7P_II, k), dpp[p7S_IL]);    // ditto

	    }
	  else rhoc[p7S_IL] = 0.;

	  /* If v-1 is k-1, then v-1,v are contiguous, and supercell i,k-1) exists; if so, its MD cells appportion mass to D(i,k)  */
	  if (v > 0 && sm->k[i][v-1] == k-1 && k > k0) 
	    {
	      rhoc[p7S_DL] = 
		P7_MASSTRACE_INCREMENT(rhoc[p7S_ML-p7S_NSCELLS], dpc[p7S_DL] + TSC(p7P_MD, k-1), dpc[p7S_ML-p7S_NSCELLS]) +  // yes, those array indices are negative;
		P7_MASSTRACE_INCREMENT(rhoc[p7S_DL-p7S_NSCELLS], dpc[p7S_DL] + TSC(p7P_DD, k-1), dpc[p7S_DL-p7S_NSCELLS]);   // but that's valid C,  seriously.
	    }
	  else rhoc[p7S_DL] = 0.;

	  mt->kmass[k] += rhoc[p7S_ML] + rhoc[p7S_DL];       // lower bound on how much mass is flowing right, through column k (partial sum, rows i0..i). Don't count I.
	  if (k > kbe && mt->kmass[k] > massthresh) kbe = k; // update kbe envelope end bound: rightmost k that satisfies threshold; k>kbe test avoids changing preset kbe=M on glocal
	  if (kbe == gm->M || mt->kmass[k] + rowmass < massthresh) kbe_proven_done = TRUE;  // *upper* bound on mass flowing right, by adding total mass that's still to the left of k 
	  rowmass  += rhoc[p7S_ML] + rhoc[p7S_IL]; // how much mass is still flowing down, through this row i. Don't count D.

	  rhoc += p7S_NSCELLS;  // advance rhoc, dpc ptrs by one supercell;
	  dpc  += p7S_NSCELLS;  // we will figure out its k index when the v loop rolls around now...
	}

      if (rowmass < massthresh) ibe_proven_done = TRUE; else ibe = i;

      if      (mt->imass) mt->imass[i] = rowmass;
      else if (ibe_proven_done && kbe_proven_done) break;
    } // end loop over i=i0..L

  *ret_ibe   = ibe;
  *ret_kbe   = kbe;
  return eslOK;
}
/*--------------- end, downwards recursion ----------------------*/


/*****************************************************************
 * 3. API, wrapping the Up and Down recursions
 *****************************************************************/



/* Function:  p7_SparseMasstrace()
 * Synopsis:  Calculate envelope coords, given a domain's anchor point.
 *
 * Purpose:   Determine envelope coords for an individual domain by the
 *            "mass trace" algorithm.
 *
 *            The caller has compared profile <gm> to digital sequence
 *            <dsq> of length <L>, obtaining Forward matrix <fwd> and
 *            Backward matrix <bck>, for a sparse mask that both of
 *            those matrices have copies of (<fwd->sm == bck->sm>);
 *            and we have obtained an optimal trace <tr> for that
 *            comparison (probably a Viterbi trace), which we're going
 *            to use to define the individual domains we find in
 *            <dsq>. 
 *
 *            Consider the individual domain defined by an anchor
 *            defined by position <z> in trace <tr> (see below). For
 *            each possible trace that contains the anchor, define
 *            (ia,ib,ka,kb) as the bounds of the domain on the
 *            sequence and the profile. For a glocal path, ka=1,kb=M (i.e.
 *            wing-retracted paths through terminal DG's count).
 *            
 *            (We uniquely define a domain by its "anchor": a triplet
 *            (i0,k0,st0={MDI}{LG}), such that any subpath
 *            ...B->...->anchor->...E ... is considered to be the
 *            "same" domain, when we're looking at alternative paths
 *            of the same comparison in an ensemble.)
 * 
 *            Summed over the partial posterior trace ensemble of
 *            traces containing this domain anchor, calculate
 *            cumulative distributions away from the anchor:
 *               P(ia <= i), the prob that the domain starts at least as far out as i, i<=i0
 *               P(ib >= i), prob that domain ends at least as far out as i, i>=i0
 *               and similarly P(ka <= k), P(kb >= k).
 *               
 *            For us to run this calculation, caller provides us with
 *            a sparse DP matrix <mass>, and a <P7_MASSTRACE> object
 *            <mt> which will store the necessary histograms. These objects
 *            can be reused from a previous calculation, even a smaller one;
 *            they will be reallocated if needed, and reinitialized here.
 *               
 *            Define envelope coords as:
 *               iae = \min_{i \leq i0} P(ia \leq i) \geq massthresh
 *               ibe = \max_{i \geq i0} P(ib \geq i) \geq massthresh
 *               kae = \min_{k \leq k0} P(ka \leq k) \geq massthresh
 *               kbe = \max_{k \geq k0} P(kb \geq k) \geq massthresh
 *            
 *            i.e. the widest start/end points that contain at least a
 *            posterior probability of <massthresh>, (Put another way:
 *            the probability that the domain starts at some i < iae
 *            is less than <massthresh>.) These probabilities are
 *            conditional on having the domain at all (e.g., the
 *            envelope *given* that the domain is present; our sums
 *            are over the partial ensemble of traces containing this
 *            domain's anchor triplet).
 *            
 *            Return the envelope coords in <*ret_iae>, <*ret_ibe>,
 *            <*ret_kae>, <*ret_kbe>.
 *
 * Args:      dsq     - digital sequence 1..L
 *            L       - length of digital sequence
 *            gm      - profile
 *            fwd     - sparse Forward matrix that caller calculated for gm x dsq comparison
 *            bck     - sparse Backward matrix caller calc'ed 
 *            tr      - optimal path caller wants to use to define domain structure of <dsq>
 *            z       - index in <tr>'s path of an anchor triplet (i0,k0,st0) = (<tr->i[z]>, <tr->k[z]>, <tr->st[z]>);
 *                      the anchor defines an individual domain, for purposes of identifying the "same" domain in 
 *                      different paths
 *            mass    - sparse DP matrix for the mass trace calculation; can be reused, will be reallocated if needed.
 *            mt      - cumulative histogram space needed for the calculation; can be reused, will be reallocated if needed.
 *            ret_iae - RETURN: iae envelope coord
 *            ret_ibe - RETURN: ibe envelope coord
 *            ret_kae - RETURN: kae envelope coord
 *            ret_kbe - RETURN: kbe envelope coord
 *
 * Returns:   <eslOK> on success. <mass> may have been reallocated, and
 *            it now contains the $\rho$ values from the sparse mass
 *            trace recursions, both up and down. <mt> may have been
 *            reallocated, and it now contains the cumulative $P(ia
 *            \leq i)$, etc. distributions away from the anchor.
 *            The envelope coords are in <*ret_iae,*ret_ibe,*ret_kae,*ret_kbe>.
 *
 * Throws:    <eslEMEM> if a reallocation (of <mass> or <mt>) fails.
 */
int
p7_SparseMasstrace(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *fwd, const P7_SPARSEMX *bck, const P7_TRACE *tr, int z, float massthresh, 
		   P7_SPARSEMX *mass, P7_MASSTRACE *mt,
		   int *ret_iae, int *ret_ibe, int *ret_kae, int *ret_kbe)
{
  int   st0 = tr->st[z];	/* anchor point's state (p7T_{MDI}{LG}) */
  int   k0  = tr->k[z];		/* anchor point's k position in profile (1..M) */
  int   i0;
  int   status;

  /* contract check / arg validation */
  ESL_DASSERT1 ( (fwd->type == p7S_FORWARD) );
  ESL_DASSERT1 ( (bck->type == p7S_BACKWARD) );
  ESL_DASSERT1 ( (fwd->sm == bck->sm) );
  ESL_DASSERT1 ( (fwd->sm->L == L) );
  ESL_DASSERT1 ( (fwd->sm->M == gm->M) );
  ESL_DASSERT1 ( (tr->L == L) );
  ESL_DASSERT1 ( (tr->M == gm->M) );

  /* Find i0, the last emitted i (i.e. watch out for st0=D; D states have no tr->i[] */
  while (z && ! tr->i[z]) z--;
  i0 = tr->i[z];

  /* We might be reusing <mass> dp matrix; make sure it's big enough.
   * Set its type now, so it can be Dumped or otherwise analyzed, if we need to.
   * No need to zero it out.
   */
  if ( (status = p7_sparsemx_Reinit(mass, fwd->sm)) != eslOK) return status;
  mass->type = p7S_MASSTRACE;

  /* We might be reusing <mt> workspace. Make sure it's big enough. Then zero it, which sets its L,M. */
  if ( (status = p7_masstrace_GrowTo(mt, gm->M, L)) != eslOK) return status;
  p7_masstrace_Zero  (mt, gm->M, L); 

  /* Up() recursion finds <iae>,<kae>; fills upper left quadrant of <mass> and <mt> (i<=i0, k<=k0); 
   * Down() recurstion finds <ibe>,<kbe>; fills lower right quadrant of <mass> and <mt> (i>=i0, k>=k0) 
   */
  if ( (status = masstrace_up  (dsq, L, gm, fwd, i0, k0, st0, massthresh, mass, mt, ret_iae, ret_kae)) != eslOK) return status;
  if ( (status = masstrace_down(dsq, L, gm, bck, i0, k0, st0, massthresh, mass, mt, ret_ibe, ret_kbe)) != eslOK) return status;

  /* Algorithm neglects glocal wing retraction kmass thru DGks; but on
   * glocal paths we know kmass[k]=1. When the <mt> is being used for
   * its distributions (i.e. when it's not Slim; when mt->imass is
   * non-NULL), then clean up this detail, so unit tests, etc don't
   * mess up on it. (Which explains the weird-looking test on
   * mt->imass before setting kmass; this is not a typo.)
   */
  if (mt->imass && (st0 == p7T_MG || st0 == p7T_IG || st0 == p7T_DG))
    { esl_vec_FSet(mt->kmass+1, mt->M, 1.0); }

  mt->i0  = i0;
  mt->k0  = k0;
  mt->st0 = st0;
  return eslOK;
}
/*----------------- end, API wrapper ----------------------------*/


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7SPARSE_MASSTRACE_TESTDRIVE

#include "base/p7_bg.h"
#include "base/p7_prior.h"

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

static void
utest_approx_masstrace(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L)
{
  char           msg[]  = "sparse masstrace, approx-masstrace unit test failed";
  P7_PRIOR      *pri    = NULL;
  P7_HMM        *hmm    = NULL;
  P7_PROFILE    *gm     = p7_profile_Create(M, abc);
  ESL_SQ        *sq     = esl_sq_CreateDigital(abc);       /* space for generated (homologous) target seqs              */
  P7_OPROFILE   *om     = p7_oprofile_Create(M, abc);
  P7_CHECKPTMX  *ox     = NULL;
  P7_TRACE      *vtr    = p7_trace_CreateWithPP(); /* domain anchor selection in trace_Index() requires a pp-annotated trace */
  P7_TRACE      *str    = p7_trace_Create();
  P7_SPARSEMASK *sm     = NULL;
  P7_SPARSEMX   *sxv    = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxf    = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxb    = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxd    = p7_sparsemx_Create(NULL);
  float         *wrk    = NULL;	/* reusable scratch workspace needed by stochastic trace */
  P7_MASSTRACE  *mte    = NULL;
  P7_MASSTRACE  *mta    = NULL;
  int            d;
  int            ntr    = 10000;
  int            i, ntrc;
  float          tol    = 0.02;
  int            i0,k0,st0;	/* anchor triplet */
  int            iae,ibe,kae,kbe;
  int            z;
  float          fsc;
  char           errbuf[eslERRBUFSIZE];

  /* Sample a profile from prior. Config as usual, multihit dial-mode. 
   * We sample from the prior, as opposed to uniform p7_modelsample(), 
   * in order to get a more realistic, info-rich profile. This utest
   * needs high posterior prob anchors in its domain alignment paths,
   * to be able to sample domain endpoints effectively.
   */
  if       (abc->type == eslAMINO) pri = p7_prior_CreateAmino();
  else if  (abc->type == eslDNA)   pri = p7_prior_CreateNucleic();
  else if  (abc->type == eslRNA)   pri = p7_prior_CreateNucleic();
  else                             pri = p7_prior_CreateLaplace(abc);

  if ( p7_modelsample_Prior(rng, M, abc, pri, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)               != eslOK) esl_fatal(msg);
  if ( p7_oprofile_Convert(gm, om)                  != eslOK) esl_fatal(msg);

  /* Generate (sample) a sequence from the profile */
  if ( p7_profile_SetLength(gm, L)                 != eslOK) esl_fatal(msg);
  do {
    esl_sq_Reuse(sq);
    p7_ProfileEmit(rng, hmm, gm, bg, sq, NULL);
  } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
  if ( p7_profile_SetLength(gm, sq->n)       != eslOK) esl_fatal(msg);
  if ( p7_oprofile_ReconfigLength(om, sq->n) != eslOK) esl_fatal(msg);

  /* Fwd/Bck local filter to calculate the sparse mask */
  if (  (ox = p7_checkptmx_Create(M, sq->n, ESL_MBYTES(32)))    == NULL) esl_fatal(msg);
  if (  (sm = p7_sparsemask_Create(M, sq->n))                   == NULL) esl_fatal(msg);

  if ( p7_ForwardFilter (sq->dsq, sq->n, om, ox, /*fsc=*/NULL) != eslOK) esl_fatal(msg);
  if ( p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSEMASK_THRESH_DEFAULT)           != eslOK) esl_fatal(msg);

   /* Sparse DP calculations, and exact posterior decoding */
  if ( p7_SparseViterbi   (sq->dsq, sq->n, gm, sm, sxv, vtr, NULL) != eslOK) esl_fatal(msg);
  if ( p7_SparseForward   (sq->dsq, sq->n, gm, sm, sxf,     &fsc)  != eslOK) esl_fatal(msg);
  if ( p7_SparseBackward  (sq->dsq, sq->n, gm, sm, sxb,      NULL) != eslOK) esl_fatal(msg);
  if ( p7_SparseDecoding  (sq->dsq, sq->n, gm, sxf, sxb, sxd)      != eslOK) esl_fatal(msg);
  p7_sparsemx_TracePostprobs(sxd, vtr); /* selecting domain anchors requires pp annotation of the trace */
  p7_trace_Index(vtr);
  p7_sparsemx_Reuse(sxv);	/* we'll reuse it for mass trace dp below, like we do in production pipeline */

  /* Choose one domain at random (we know there's at least one) */
  if (!vtr->ndom)   esl_fatal(msg);
  d = esl_rnd_Roll(rng, vtr->ndom);
  if (!vtr->anch[d]) esl_fatal(msg);

  /* Determine the anchor triplet (CountTrace needs it, whereas SparseMasstrace works it out for itself */
  for (i0=0, z = vtr->anch[d]; z && !i0; z--) i0 = vtr->i[z];
  k0  = vtr->k[vtr->anch[d]];
  st0 = vtr->st[vtr->anch[d]];

  //p7_trace_DumpAnnotated(stdout, vtr, gm, sq->dsq);

  /* Stochastic trace ensemble approximation to mass trace */
  ntrc = 0;
  if ( ( mta = p7_masstrace_Create(gm->M, sq->n))   == NULL)  esl_fatal(msg);
  if (         p7_masstrace_Zero(mta, gm->M, sq->n) != eslOK) esl_fatal(msg);
  for (i = 0; i < ntr; i++)
    {
      p7_sparse_trace_Stochastic(rng, &wrk, gm, sxf, str);
      //p7_trace_DumpAnnotated(stdout, str, gm, sq->dsq);
      p7_masstrace_CountTrace(str, i0, k0, st0, mta, &ntrc);
      p7_trace_Reuse(str);
    }
  p7_masstrace_FinishCount(mta, ntrc);

  /* Exact mass trace algorithm */
  if ( ( mte = p7_masstrace_Create(gm->M, sq->n)) == NULL)  esl_fatal(msg);
  if ( p7_SparseMasstrace(sq->dsq, sq->n, gm, sxf, sxb, vtr, vtr->anch[d], 0.1, sxv, mte, &iae, &ibe, &kae, &kbe) != eslOK) esl_fatal(msg);
  
  //printf("max = %f\n", p7_masstrace_GetMaxAbsDiff(mte, mta));
  //p7_masstrace_Dump(stdout, mte);
  //p7_masstrace_Dump(stdout, mta);

  /* Tests */
  if ( p7_masstrace_Validate(mte, errbuf)  != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
  if ( p7_masstrace_Validate(mta, errbuf)  != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
  if ( p7_masstrace_Compare(mte, mta, tol) != eslOK) esl_fatal(msg);

  /* Clean up */
  p7_masstrace_Destroy(mta);
  p7_masstrace_Destroy(mte);
  if (wrk) free(wrk);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxv);
  p7_sparsemask_Destroy(sm);
  p7_trace_Destroy(str);
  p7_trace_Destroy(vtr);
  p7_checkptmx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_prior_Destroy(pri);
  esl_sq_Destroy(sq);
}
#endif /*p7SPARSE_MASSTRACE_TESTDRIVE*/
/*------------------ end, unit tests ----------------------------*/



/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7SPARSE_MASSTRACE_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

/* This is a stochastic unit test, with some probabilistic variability expected, so we use a fixed RNG seed 
 * to avoid rare expected excursions causing rare failures.
 */
static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for sparse mass trace algorithms";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg   = p7_bg_Create(abc);
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  utest_approx_masstrace(r, abc, bg, M, L);

  fprintf(stderr, "#  status = ok\n");

  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}

#endif /*p7SPARSE_MASSTRACE_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/


/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef p7SPARSE_MASSTRACE_EXAMPLE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",              0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of Forward/Backward, sparse dual implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
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
  P7_SPARSEMASK  *sm      = NULL;
  P7_SPARSEMX    *sxv     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxf     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxb     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxd     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxm     = p7_sparsemx_Create(NULL);
  P7_MASSTRACE   *mt      = p7_masstrace_Create(100,100); /* M,L hints. Will be grown. */
  P7_TRACE       *tr      = p7_trace_CreateWithPP();
  float           vsc, fsc;
  int             iae,ibe;
  int             kae,kbe;
  int             d,z;
  int             status;

  ESL_RANDOMNESS *rng     = esl_randomness_Create(42);
  float          *wrk     = NULL;
  P7_TRACE       *str     = p7_trace_Create();
  int             ntr     = 10000;
  int             i, ntrx;
  int             i0,k0,st0;

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
  p7_profile_Config(gm, hmm, bg);

  /* Allocate bands, matrices */
  sm  = p7_sparsemask_Create(gm->M, sq->n);
  p7_sparsemask_AddAll(sm);

  /* Set the profile and null model's target length models */
  p7_bg_SetLength           (bg, sq->n);
  p7_profile_SetLength      (gm, sq->n);

  /* Sparse DP calculations */
  p7_SparseViterbi (sq->dsq, sq->n, gm, sm, sxv, tr, &vsc);
  p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf,     &fsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxb,     NULL);
  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxb, sxd);
  p7_sparsemx_TracePostprobs(sxd, tr);
  p7_trace_Index(tr);

  //p7_sparsemx_Dump(stdout, sxf);
  //p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);
  
  for (d = 0; d < tr->ndom; d++)
    {
      /* Determine the anchor triplet (CountTrace needs it, whereas SparseMasstrace works it out for itself */
      for (i0=0, z = tr->anch[d]; z && !i0; z--) i0 = tr->i[z];
      k0  = tr->k[tr->anch[d]];
      st0 = tr->st[tr->anch[d]];

      /* stochastic ensemble approximation */
      ntrx = 0;
      p7_masstrace_GrowTo(mt, gm->M, sq->n);
      p7_masstrace_Zero(mt, gm->M, sq->n);
      for (i = 0; i < ntr; i++)
	{
	  p7_sparse_trace_Stochastic(rng, &wrk, gm, sxf, str);
	  p7_masstrace_CountTrace(str, i0, k0, st0, mt, &ntrx);
	  p7_trace_Reuse(str);
	}
      p7_masstrace_FinishCount(mt, ntrx);
      p7_masstrace_PlotImass(stdout, mt);
      p7_masstrace_Reuse(mt);

      /* Mass trace calculation */
      p7_SparseMasstrace(sq->dsq, sq->n, gm, sxf, sxb, tr, tr->anch[d], 0.1, sxm, mt, &iae, &ibe, &kae, &kbe);

      p7_masstrace_PlotImass(stdout, mt);
      //p7_masstrace_PlotKmass(stdout, mt);

      //      printf("# domain %3d  iali: %d..%d [%daa]  ienv: %d..%d [%daa]  kali: %d..%d [%daa]  kenv: %d..%d [%daa]\n",
      //	     d,
      //	     tr->sqfrom[d],  tr->sqto[d],  tr->sqto[d]-tr->sqfrom[d]+1, iae, ibe, ibe-iae+1,
      //	     tr->hmmfrom[d], tr->hmmto[d], tr->hmmto[d]-tr->hmmfrom[d]+1, kae, kbe, kbe-kae+1);
      p7_sparsemx_Reuse(sxm);
      p7_masstrace_Reuse(mt);
    }
  
  /* Cleanup */
  p7_trace_Destroy(str);
  esl_randomness_Destroy(rng);
  free(wrk);

  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_sparsemx_Destroy(sxv);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemx_Destroy(sxm);
  p7_sparsemask_Destroy(sm);
  p7_masstrace_Destroy(mt);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
/*------------------ end, example driver ------------------------*/
#endif /*p7SPARSE_MASSTRACE_EXAMPLE*/



/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

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
 *    over all edges into Y gives us the fraction of posterior
 *    probability mass that flows back that edge: call this the _edge
 *    weight_. The sum of all edge weights into Y is 1. Since
 *    everything's in null-scaled log space, the ratio shows up as a
 *    subtraction and an exp().
 *    
 *    For a Y=Mk match state, one such incoming edge weight is the
 *    B->Mk edge; mass that flows back that edge disappears from our
 *    recursion.  Essentially, it's this lossage that we're tracking
 *    back: how far do we have to trace back in i,k before enough mass
 *    is lost from the main state alignment ensemble that the mass
 *    we're still tracing back drops below the <massthresh> threshold?
 *    
 *    In a stochastic traceback, we would trace back from Y
 *    stochastically, according to its edge weight distribution
 *    (inclusive of the B->Mk edges). Notationally, we implement a
 *    stochastic traceback as a _push_ from a current state Y to a
 *    next state X, whereas the mass traceback is implemented as a
 *    _pull_ to a current state X from connected states Y.
 */


