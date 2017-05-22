/* Production implementation of anchor set constrained (ASC) Forward
 * and Backward, using sparse dynamic programming.
 * 
 * Contents:
 *    1. Sparse ASC Forward
 *    2. Sparse ASC Forward, segmental version
 *    3. Sparse ASC Backward
 *    4. Sparse ASC Decoding
 *    5. Footnotes
 *    6. Unit tests
 *    7. Test driver
 *    8. Example
 */

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_anchors.h"

#include "misc/logsum.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spascmx.h"

#include "dp_sparse/sparse_asc_fwdback.h"

/*****************************************************************
 * 1. Sparse ASC Forward 
 *****************************************************************/

/* Function:  p7_sparse_asc_Forward()
 * Synopsis:  Sparse anchor set constrained Forward algorithm.
 *
 * Purpose:   Run the sparse anchor set constrained (ASC) Forward
 *            algorithm, comparing profile <gm> to target sequence
 *            <dsq> of length <L>, constrained both by the anchor set
 *            array <anch> of <1..D> anchors, and by the sparse mask
 *            <sm>. Fills in the sparse ASC Forward matrix <asf> and
 *            optionally returns the sparse ASC Forward score (raw, in
 *            nats) in <opt_sc>.
 *            
 *            Caller provides <asf> space, but it can be allocated to
 *            anything; it is reallocated here as needed.
 *
 * Args:      dsq    : target digital sequence, 1..L
 *            L      : length of <dsq>
 *            gm     : query profile
 *            anch   : anchor set array, 1..D
 *            D      : number of anchors in <anch>
 *            sm     : sparse mask that constrains the DP
 *            asf    : RESULT : sparse ASC Forward DP matrix
 *            opt_sc : optRETURN : sparse ASC Forward score
 *
 * Returns:   <eslOK> on success.
 *  
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_sparse_asc_Forward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D, 
			  const P7_SPARSEMASK *sm, P7_SPARSEMX *asf, float *opt_sc)
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

  /* Contract checks, argument validation */
  ESL_DASSERT1(( sm->L == L     ));
  ESL_DASSERT1(( sm->M == gm->M ));

  /* Assure that <asf> is allocated large enough. 
   */
  if (( status = p7_spascmx_Resize(asf, sm, anch, D)) != eslOK) return status;
  asf->type = p7S_ASC_FWD;

  /* Iterate the ForwardSeg algorithm over each segment in turn. 
   */
  dpc = asf->dp;
  xc  = asf->xmx;
  for (g = 1; g <= sm->S; g++)
    {
      ngap = (g == 1 ? sm->seg[g].ia-1 : sm->seg[g].ia - sm->seg[g-1].ib - 1);
      p7_sparse_asc_ForwardSeg(dsq, L, gm, anch, D, d, sm, ngap, sm->seg[g].ia, sm->seg[g].ib, 
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





/*****************************************************************
 * 2. Sparse ASC Forward Segment
 *****************************************************************/


/* Function:  p7_sparse_asc_ForwardSeg()
 * Synopsis:  Computes sparse ASC Forward on a single sparse segment
 *
 * Purpose:   This routine is used in three ways:
 *            
 *            First, we use it in segmental divide and conquer in the
 *            fastest MPAS algorithm. There we need to calculate (and
 *            re-calculate) the ASC Forward matrix for a sparse
 *            segment, using different guesses for anchor sets within
 *            it.
 *            
 *            Second, we compose the overall ASC Forward algorithm
 *            (for one known anchor set) by calling ForwardSeg on each
 *            segment in turn.
 *            
 *            Third, to calculate envelope scores, we may need to do
 *            an ASC Forward calculation for a single anchor between
 *            two specific coordinates <ia..ib> that lie within some
 *            sparse segment.
 *            
 *            The first and second uses require precise handling of
 *            information about the starting and ending values at
 *            N/J/C states and position of ptrs into the ASC Forward
 *            DP matrix because we're building up a complete matrix,
 *            one segment at a time (and in the MPAS case, rewriting
 *            the currently last segment many times).
 *            
 *            The third use (envelope scores) doesn't need to build up
 *            a complete matrix, but does require using ia/ib coords
 *            that are within a segment g, rather than being exactly
 *            that sparse segment's endpoints.
 *            
 *            So here's how we do all that:
 *            
 *            We're comparing profile <gm> to digital target sequence
 *            <dsq> of length <L>, constrained both by sparse mask
 *            <sm> and by anchor set <anch> containing <D> domain
 *            anchors. 
 *            
 *            In D&C MPAS, <anch> is a prefix up to and including the
 *            current segment. In full ASC Forward <anch> is a
 *            complete set. In envelope scores, <anch> is a single
 *            anchor (with sentinels) and <D=1>.
 *            
 *            Caller provides sparse ASC matrix <asf> and (unlike our
 *            other DP implementations), we do no reallocation
 *            here. Caller is responsible for making sure that <asf>
 *            has sufficient memory to hold the calculation. This is
 *            because the demands of complete-sequence ASC Forward
 *            vs. segmental divide and conquer MPAS vs. envelope
 *            scoring are different.
 *            
 *            Caller must also provide a bunch of state information
 *            about where we are in an ongoing DP calculation, from
 *            previous segments we've calculated:
 *            
 *            <d> is the index (1..D,D+1) of the next domain anchor in
 *            <anch> that we can reach in <dsq> (even if it is
 *            downstream of the current segment, and thus not reached
 *            in this segment); if it is D+1, that means we're already
 *            past all anchors (and no anchors occur in this segment).
 *            <anch> can either be the anchor set for the entire
 *            sequence (in which case <d> may be somewhere inside that
 *            array; if there is no anchor in the segment, <i0(d) >
 *            ib>), or it may be an anchor set only for this
 *            segment (in which case <d==1>; if there is no anchor in
 *            the segment <d==1> and <D==0>). 
 *            
 *            <ngap> is the number of residues since the last segment
 *            we calculated. This is either <ia-1> (if this is first
 *            or only segment), or <ia(g) - ib(g-1) - 1> (if there was
 *            a previous segment <g-1>, as in MPAS or Forward).

 *            <ia>..<ib> are the coords of the segment we're to calculate.
 *            For MPAS and full Forward, they are <seg[g].ia..seg[g].ib>
 *            for the current segment <g>. For envelope scores, they are
 *            <seg[g].ia <= ia <= ib <= seg[g].ib> coords within whatever
 *            segment the envelope is in (and we don't need to know <g>).
 *            
 *            <xN>, <xJ>, <xC> are the ASC Forward values at states
 *            N,J,C on the last row of the previous segment <ib(g-1)>.
 *            For envelope scores, or for the first segment, <xN=0>
 *            and <xJ=xC=-inf>.
 *            
 *            <dpc>, <xc> are pointers into the sparse ASC DP matrix
 *            to the first main supercell and first special supercell
 *            for this segment, respectively (i.e. the next supercells
 *            after the previous segment's calculation). (We need this
 *            state information because addressing in sparse matrices
 *            is implicit; we only know where we are by sequentially
 *            traversing the matrix.) For envelope scores, or for the
 *            first segment, <dpc=asf->dp> and <xc=asf->xmx>.
 *            
 *            After the DP calculation for this segment is complete,
 *            we optionally get state information back, for use at the
 *            next call (if any).  <opt_xN>, <opt_xJ>, <opt_xC>
 *            contain the N/J/C Forward values at <ib> (env scores
 *            don't need these). <opt_d> optionally contains the index
 *            of the next anchor we can reach (neither MPAS nor env
 *            scores need this). <opt_dpn> and <opt_xn> are pointers
 *            to the next supercells we can store in <asf> (env scores
 *            don't need these).  It is safe for these to be ptrs to
 *            the same variables as the input state variables,
 *            allowing input to be overwritten with the output.
 *            
 *            For initialization conditions on the first segment in a
 *            sequence, we would pass <d=1>, <ngap = seg[1].ia-1>, <ia
 *            = seg[1].ia>, <ib = seg[1].ib>, <xN=0>, <xJ=xC=-inf>,
 *            <dpc=asf->dp>, and <xc=asf->xmx>.
 *            
 *            For envelope score of <oa..ob>, we create an anchor
 *            set containing the one anchor for our envelope (<D=1>),
 *            and pass <d=1>, <ngap = oa-1>, <ia=oa>, <ib=ob>,
 *            <xN=0>, <xJ=xC=-inf>, <dpc=asf->dp>, and <xc=asf->xmx>.
 *            The only state information that an envelope score needs
 *            to get back is <opt_xC>.
 *            
 *            In a complete sequence Forward, for subsequent segments,
 *            we pass the state information we got from the previous
 *            segment's return. In a segmental divide-and-conquer
 *            MPAS, to recalculate the same segment with a different
 *            anchor set, we can call as many times as we like using
 *            the same state information, which results in overwriting
 *            previous calculation of this segment.
 *            
 *            When we calculate the last segment <g=sm->S>, <opt_xC>
 *            is the C value at <ib(S)>. To obtain a complete-sequence
 *            ASC Forward or envelope score, caller then must add CC
 *            terms for residues <ib+1..L> (there are L-ib of
 *            them), and the <t_CT> term.
 *            
 *            Upon return, <asf> has been extended to include the DP
 *            matrix storage for this segment, and <ret_asc> contains
 *            the sparse ASC Forward segment score. This segment score
 *            is the log sum of all paths from the starting N or J
 *            state at <ia(g)-1> to the ending N or J state at <ib(g)>.
 *            
 *            The segment score makes the assumption that
 *            <t_NN>=<t_JJ>= <t_CC>. However, nothing else in this
 *            routine makes that assumption. Segmental divide and
 *            conquer MPAS requires the assumption (and uses segment
 *            scores). Complete-sequence ASC Forward does not require
 *            the assumption (and does not use segment scores).
 *
 * Args:      dsq     : digital sequence
 *            L       : length of <dsq>
 *            gm      : profile
 *            anch    : array of anchors 1..D, with sentinels at 0,D+1
 *            D       : number of anchors in <anch>
 *            d       : next anchor that we can reach in this or subsequent segments; 1..D+1 (D+1 means no anchors in this or subseq segments)
 *            sm      : sparse mask
 *            ngap    : number of previous residues since last segment: ia-1 or ia-ib(g-1)-1
 *            ia      : start coord of segment: ia(g), or oa for an envelope score
 *            ib      : end coord of segment:   ib(g), or ob for an envelope score
 *            asf     : sparse ASC matrix that already contains calculation for segments 1..g-1;
 *                      caller guarantees that <asf> is allocated for the DP calculation of segment g
 *            xN      : last N value at ib(g-1), or 0 for g=1
 *            xJ      : last J value at ib(g-1), or -eslINFINITY for g=1
 *            xC      : last C value at ib(g-1), or -eslINFINITY for g=1
 *            dpc     : ptr to next main supercell in <asf> where this segment g is to start storage
 *            xc      : ptr to next special supercell in <asf> where segment g is to start storage
 *            opt_xN  : RETURN: N value at ib(g)
 *            opt_xJ  : RETURN: J value at ib(g)
 *            opt_xC  : RETURN: C value at ib(g)
 *            opt_d   : optRETURN: next anchor we can see. ib(g) < ia(g+1) <= i0(*ret_d)
 *            opt_dpn : RETURN: ptr to next main supercell in <asf>, where next segment starts storage.
 *                              For last segment <g==sm->S>, may point outside valid memory.    
 *            opt_xn  : RETURN: ptr to next special supercell in <asf>, where next segment will start storage. 
 *                              For last segment <g==sm->S>, may point outside valid memory.
 *            opt_asc : optRETURN: sparse ASC Forward segment score. Invalid unless t_NN=t_JJ=t_CC.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_sparse_asc_ForwardSeg(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
			 const P7_ANCHOR     *anch, int D, int d,
			 const P7_SPARSEMASK *sm,   int ngap, int ia, int ib,
			 P7_SPARSEMX *asf, float xN, float xJ, float xC, float *dpc, float *xc,
			 float *opt_xN, float *opt_xJ, float *opt_xC, int *opt_d, float **opt_dpn, float **opt_xn, float *opt_asc)
{
  float const *tsc = gm->tsc;	 // sets up TSC() macro, access to profile's transitions      
  float const *rsc = NULL;	 // will be set up for MSC(), ISC() macros for residue scores for each row i
  int    ndown     = 0;          //   ... this bumps to 1 when i reaches an anchor; row # of DOWN sector; back to 0 at seg end
  float *dpp       = NULL;       // ptr for stepping thru DP cells on a previous row. Gets initialized by...
  float *last_down = NULL;       //   ... remembering the <dpc> start of each DOWN...
  float *last_up   = NULL;       //   ... or UP sector as we compute them.
  int    Ds        = 0;          // # of anchors we see in this segment; if we see 0 we need to know that for some special casing
  float  mlc, mgc;               // tmp variables for computing MLk, MGk score
  float  dlc, dgc;               //   ... and DLk, DGk.
  float  xE,xB,xL,xG;            // current special state values on this row; specials are not always stored
  int    i;                      // index for dsq[i] seq position/rows 0..L
  int    k;                      // index of model positions 1..M
  int    z;                      // index for sparse cell indices k[i][z] along a row
  int    y;			 //   ... and for previous row, k[i-1][y].
  int    s;                      // index of individual states
  float  xNJbase;                // Save xN | xJ at ia(g)-1, whichever one wasn't -inf

  
  ESL_DASSERT1(( asf->type == p7S_ASC_FWD ));  // Caller already allocated <asf> and set its type


  /* Initialize on row ia(g)-1
   * Follows "always calculate, sometimes store" rule for syncing numerical error w/ SparseForward(); see footnote [1].
   */
  xNJbase = ESL_MAX(xN, xJ);                                                  // Either xN == -inf || xJ == -inf, depending on whether we've seen an anchor yet in dsq
  xE = -eslINFINITY;
  xN = (ngap == 0 ? xN : xN + (float) ngap * gm->xsc[p7P_N][p7P_LOOP]);       // ngap=0 only occurs for g=1, ia(1)=1 (1st residue is in 1st seg)
  xJ = (ngap == 0 ? xJ : xJ + (float) ngap * gm->xsc[p7P_J][p7P_LOOP]);       // ngap==0 test is there because t_NN may be -inf; must avoid 0 * -inf = NaN
  xB = ESL_MAX(xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]); // Either xN == -inf || xJ == -inf, depending on whether we've seen an anchor yet in dsq
  xL = xB + gm->xsc[p7P_B][0];                                                // [0] is B->L 
  xG = xB + gm->xsc[p7P_B][1];                                                // [1] is B->G
  xC = (ngap == 0 ? xC : xC + (float) ngap * gm->xsc[p7P_C][p7P_LOOP]);       // Up 'til 1st anchor is reached, xJ=xC=-inf; afterwards, xN=-inf. 

  if ( anch[d].i0 <= ib)  // suffices to test if there's an anchor in this segment or not
    {
      *xc++ = xE;  
      *xc++ = xN; 
      *xc++ = xJ;  
      *xc++ = xB; 
      *xc++ = xL; 
      *xc++ = xG;  
      *xc++ = xC;
      *xc++ = -eslINFINITY; // JJ; only valid in a decoding mx
      *xc++ = -eslINFINITY; // CC; ditto
    }

  for (i = ia; i <= ib; i++)
    {
      if      (i == anch[d].i0) { ndown = 1; d++; Ds++; }  // When i reaches next anchor, bump d to next domain index, and DOWN sector becomes active.
      else if (ndown)           { ndown++;              }  // Counting ndown lets us determine when we're in 1st vs. subsequent rows; useful for boundary conditions on DP 

      xE = -eslINFINITY;

      /*****************************************************************
       * Computing row i when it is in a DOWN sector
       *****************************************************************/
      if (ndown == 1)                                // topmost row of a DOWN sector is a special case for initialization.
	{
	  rsc       = gm->rsc[dsq[i]];	             // now MSC(k), ISC(k) residue score macros work 
	  dpp       = last_up;                       // yes, UP: this initialization of i0,k0 in DOWN is the connecting thread from UP to DOWN.
	  last_down = dpc;
	  y = z     = 0;

	  while ( sm->k[i][z] < anch[d-1].k0) z++;   // skip z to the anchor cell for row i. (Which MUST exist; we don't even need to check z<n[i])
	  k = sm->k[i][z];

	  mlc = xL + TSC(p7P_LM, k-1);
	  mgc = xG + TSC(p7P_GM, k-1);
	  while (dpp && y < sm->n[i-1] && sm->k[i-1][y] <  k-1) { y++; dpp += p7S_NSCELLS; }  // dpp is an UP row. i-1,k-1 may not exist (UP sector can be empty), hence the check for non-NULL dpp
	  if    (dpp && y < sm->n[i-1] && sm->k[i-1][y] == k-1) {                             //   ... but usually, there's an UP sector with a i-1,k-1 cell 
	    mlc = p7_FLogsum( p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MM, k-1),
					  dpp[p7R_IL] + TSC(p7P_IM, k-1)),
			      p7_FLogsum( dpp[p7R_DL] + TSC(p7P_DM, k-1),
					  mlc));
	    mgc = p7_FLogsum( p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MM, k-1),
					  dpp[p7R_IG] + TSC(p7P_IM, k-1)),
			      p7_FLogsum( dpp[p7R_DG] + TSC(p7P_DM, k-1),
					  mgc));
	  }
	  *dpc++ = mlc = MSC(k) + mlc;
	  *dpc++ = mgc = MSC(k) + mgc;
	  *dpc++ = -eslINFINITY;                      // IL(i0,k0)
	  *dpc++ = -eslINFINITY;                      // IG(i0,k0)
	  xE     = p7_FLogsum(xE, mlc);               // dlc is -inf, not included in sum
	  *dpc++ = -eslINFINITY;                      // DL(i0,k0)
	  *dpc++ = -eslINFINITY;                      // DG(i0,k0)
	  if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {                     // is there a (i,k+1) cell to our right? 
	    dlc =  mlc + TSC(p7P_MD, k);                                    //  ... then ->Dk+1 is precalculated as usual.
	    dgc =  mgc + TSC(p7P_MD, k);                                    // Since we know we're the first cell (the anchor) we initialize dlc, dgc here.
	  } else {                                                          // If not, we MUST add {MD}l->Dk+1..E glocal exit path, even from internal sparse cells - not just last cell! 
	    xE  = p7_FLogsum( xE, mgc + TSC(p7P_DGE, k) + TSC(p7P_MD, k));
	    dlc = dgc = -eslINFINITY;
	  }
	  
	  for (z = z+1; z < sm->n[i]; z++)                                  // Now the rest of the anchor row for DOWN. Can only be a ->DDD->E path.
	    {
	      k = sm->k[i][z];
	      
	      *dpc++ = -eslINFINITY;
	      *dpc++ = -eslINFINITY;
	      *dpc++ = -eslINFINITY;
	      *dpc++ = -eslINFINITY;
	      xE     = p7_FLogsum(xE, dlc);
	      *dpc++ = dlc;
	      *dpc++ = dgc;
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {               // see comment above for ->Dk+1 precalculation.
		dlc = dlc + TSC(p7P_DD, k);
		dgc = dgc + TSC(p7P_DD, k);
	      } else {   
		xE  = p7_FLogsum( xE, dgc + TSC(p7P_DGE, k) + TSC(p7P_DD, k)); 
		dlc = dgc = -eslINFINITY;
	      }
	    }
	}

      else if (ndown)                                             // Main recursion for DOWN rows. 
	{
	  rsc = gm->rsc[dsq[i]];	                          // now MSC(k), ISC(k) residue score macros work 
	  y = z = 0; 
	  dlc = dgc = -eslINFINITY;

	  dpp        = last_down;
	  last_down  = dpc;

	  while (z < sm->n[i]   && sm->k[i][z]   < anch[d-1].k0) z++;   // skip sparse cells that aren't in DOWN sector (which is k0..M) 
	  while (y < sm->n[i-1] && sm->k[i-1][y] < anch[d-1].k0) y++;   //   .. and on prev row too.
	  for (; z < sm->n[i]; z++)                                 // then calculate the rest of the row.
	    {                                  
	      k = sm->k[i][z];  // for notational convenience

	      /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	      mlc = mgc = -eslINFINITY;
	      while (y < sm->n[i-1] && sm->k[i-1][y] < k-1)        { y++; dpp += p7S_NSCELLS; }     // skip cells that exist in sparse ASC matrix, but aren't (i-1,k-1)
	      if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) {
		mlc = p7_FLogsum( p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MM, k-1),
					      dpp[p7R_IL] + TSC(p7P_IM, k-1)),
				              dpp[p7R_DL] + TSC(p7P_DM, k-1));
		mgc = p7_FLogsum( p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MM, k-1),
					      dpp[p7R_IG] + TSC(p7P_IM, k-1)),
				              dpp[p7R_DG] + TSC(p7P_DM, k-1));
	      }
	      *dpc++ = mlc = MSC(k) + mlc;
	      *dpc++ = mgc = MSC(k) + mgc;

	      /* Try to find cell i-1,k; then compute I(i,k) from it */
	      if (y < sm->n[i-1] && sm->k[i-1][y] < k) { y++; dpp += p7S_NSCELLS; }
	      if (y < sm->n[i-1] && sm->k[i-1][y] == k) {
		*dpc++ = p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MI,k),  dpp[p7R_IL] + TSC(p7P_II, k)); // +ISC(k) if we weren't enforcing it to zero
		*dpc++ = p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MI,k),  dpp[p7R_IG] + TSC(p7P_II, k)); // ditto
	      } else {
		*dpc++ = -eslINFINITY;
		*dpc++ = -eslINFINITY;
	      }

	      /* Local exit paths */
	      xE = p7_FLogsum(xE, p7_FLogsum(mlc, dlc));

	      /* Delayed store of Dk. */
	      *dpc++ = dlc;
	      *dpc++ = dgc;

	      /* Advance calculation of next D_k+1 */
	      if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) { 
		dlc = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
		dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	      } else {
		xE  = p7_FLogsum( xE, TSC(p7P_DGE, k) + p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k))); 
		dlc = dgc = -eslINFINITY;
	      }
	    } // end loop over z: we've finished sparse row i for DOWN sector.
	} // end of <ndown> block, calculation of DOWN row i
	  

      /*****************************************************************
       * Computing row i when it's in an UP sector 
       *****************************************************************/

      if (anch[d].i0 <= ib)              // d may be D+1 here: if so, sentinel makes the the comparison to seg[g].ib fail
	{
	  if (ndown == 1)                // The first row of an UP matrix for subsequent domains in a segment is the anchor row i0.
	    {                            // All sparse cells in this row are unreachable, initialized to -inf. They only get used for G->DDDD->Mk,i+1 decoding, wing unfolding.
	      last_up = dpc;
	      for (z = 0; z < sm->n[i] && sm->k[i][z] < anch[d].k0; z++)
		for (s = 0; s < p7S_NSCELLS; s++)
		  *dpc++ = -eslINFINITY;
	    }

	  else if (i == ia)                    // The first row of UP(d) when d is the first domain in a segment
	    {                                  // is the first row of the segment. Cells on this row can be reached
	      rsc       = gm->rsc[dsq[i]];     // by {GL}->Mk entries, followed by {Dk,Mk}->Dk+1 delete transitions.
	      dpp       = NULL;                // The (ndown==1) code must precede this, to deal with a case of an anchor on first row of a segment, which means a nonexistent UP sector and immediate init in DOWN
	      last_up   = dpc;
	      dlc = dgc = -eslINFINITY;

	      for (z=0, y=0; z < sm->n[i] && sm->k[i][z] < anch[d].k0; z++)  // for all sparse cells in UP sector on this row
		{
		  k = sm->k[i][z];

		  *dpc++ = mlc = xL  + TSC(p7P_LM, k-1) + MSC(k);
		  *dpc++ = mgc = xG  + TSC(p7P_GM, k-1) + MSC(k);
		  *dpc++       = -eslINFINITY;
		  *dpc++       = -eslINFINITY;
		  *dpc++       = dlc;
		  *dpc++       = dgc;
		  // Advance calculation of next D_k+1, if it's there
		  if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {  // is there an (i,k+1) cell to our right?
		    dlc = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
		    dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
		  } else 
		    dlc = dgc = -eslINFINITY;
		}
	    }

	  else                      // all other cases: we're within UP(d), not on its first row; 
	    {                       // so i-1 row exists in UP(d) and sparse mask; we may look at dp[i-1] and sm->k[i-1]
	      rsc = gm->rsc[dsq[i]];	                   // now MSC(k), ISC(k) residue score macros work 
	      dpp       = last_up;
	      last_up   = dpc;
	      dlc = dgc = -eslINFINITY;

	      for (z=0, y=0; z < sm->n[i] && sm->k[i][z] < anch[d].k0; z++)  // for all sparse cells in UP sector on this row
		{
		  k = sm->k[i][z];

		  // Try to find cell i-1,k-1. Then compute M(i,k) from it.
		  mlc = xL  + TSC(p7P_LM, k-1);
		  mgc = xG  + TSC(p7P_GM, k-1);
		  while (y < sm->n[i-1] && sm->k[i-1][y] < k-1)  { y++; dpp+= p7S_NSCELLS; }
		  if    (y < sm->n[i-1] && sm->k[i-1][y] == k-1) {
		    mlc = p7_FLogsum( p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MM, k-1),
						  dpp[p7R_IL] + TSC(p7P_IM, k-1)),
				      p7_FLogsum( dpp[p7R_DL] + TSC(p7P_DM, k-1),
						  mlc));
		    mgc = p7_FLogsum( p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MM, k-1),
						  dpp[p7R_IG] + TSC(p7P_IM, k-1)),
				      p7_FLogsum( dpp[p7R_DG] + TSC(p7P_DM, k-1),
						  mgc));
		  }
		  *dpc++ = mlc = MSC(k) + mlc;
		  *dpc++ = mgc = MSC(k) + mgc;

		  // Try to find cell i-1, k. Then compute I(i,k) from it.
		  if (y < sm->n[i-1] && sm->k[i-1][y] < k) { y++; dpp += p7S_NSCELLS; }
		  if (y < sm->n[i-1] && sm->k[i-1][y] == k) {
		    *dpc++ = p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MI,k),  dpp[p7R_IL] + TSC(p7P_II, k)); // +ISC(k) if we weren't enforcing it to zero
		    *dpc++ = p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MI,k),  dpp[p7R_IG] + TSC(p7P_II, k)); // ditto
		  } else {
		    *dpc++ = -eslINFINITY;
		    *dpc++ = -eslINFINITY;
		  }

		  // Delayed store of Dk
		  *dpc++ = dlc;
		  *dpc++ = dgc;

		  // Advance calculation of next D_k+1, if it's there
		  if (z < sm->n[i]-1 && sm->k[i][z+1] == k+1) {  // is there an (i,k+1) cell to our right?
		    dlc = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
		    dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
		  } else 
		    dlc = dgc = -eslINFINITY;
	      
		} // End of loop over z, sparse cells on row i
	    } // End of dealing with all but the first row in an UP sector
	} // End of UP sector.


      /*****************************************************************
       * Specials. Always calculated; sometimes stored.
       *****************************************************************/

      /* To prevent discrepancies caused by roundoff error
       * accumulation, calculation of specials must mirror what
       * SparseForward() does, reinitializing every ia(g)-1 even if
       * that segment g isn't stored by ASC. Thus we separate
       * *calculation* of specials (by SparseForward() rules) from
       * *storage* of specials (by ASC rules).  For full explanation
       * of why this is so, see footnote [1].
       */

      /* Again, this follows "always calculate, sometimes store": see footnote [1] */
      /* xE was already collected above */
      if (Ds == 0)    xN = xN + gm->xsc[p7P_N][p7P_LOOP];                             // which might still be -inf, if xN started as -inf in this seg
      else            xN = -eslINFINITY;
      if (ndown != 1) xJ = p7_FLogsum(xE + gm->xsc[p7P_E][p7P_LOOP], xJ + gm->xsc[p7P_J][p7P_LOOP]);
      else            xJ =            xE + gm->xsc[p7P_E][p7P_LOOP];                  // block J->J propagation across anchor i0
      if (ndown != 1) xC = p7_FLogsum(xE + gm->xsc[p7P_E][p7P_MOVE], xC + gm->xsc[p7P_C][p7P_LOOP]);  
      else            xC =            xE + gm->xsc[p7P_E][p7P_MOVE];
      xB = p7_FLogsum(xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]);  // only one term is finite, so the FLogsum() is overkill
      xL = xB + gm->xsc[p7P_B][0];                                                
      xG = xB + gm->xsc[p7P_B][1];                                                


      if (anch[d].i0 <= ib || ndown)
	{
	  *xc++ = xE;
	  *xc++ = xN;
	  *xc++ = xJ;
	  *xc++ = xB;
	  *xc++ = xL;
	  *xc++ = xG;
	  *xc++ = xC;
	  *xc++ = -eslINFINITY; // JJ; only valid in a decoding mx
	  *xc++ = -eslINFINITY; // CC; ditto
	}
    }

  if (opt_xN)  *opt_xN  = xN;
  if (opt_xJ)  *opt_xJ  = xJ;
  if (opt_xC)  *opt_xC  = xC;
  if (opt_d)   *opt_d   = d;
  if (opt_dpn) *opt_dpn = dpc;
  if (opt_xn)  *opt_xn  = xc;
  if (opt_asc) *opt_asc = (xN != -eslINFINITY && Ds==0) ?  xN - xNJbase : xJ - xNJbase;  // Usually J(ib) - (J(ia-1) || N(ia-1)) suffices, but watch out for N..N empty path ending in N(ib)
  return eslOK;
}
/*------------ end, sparse ASC Forward segment ------------------*/





/*****************************************************************
 * 3. Sparse ASC Backward
 *****************************************************************/

/* Numerical error matching w/ sparse backward: necessary? */

int
p7_sparse_asc_Backward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
		       const P7_ANCHOR *anch, int D, const P7_SPARSEMASK *sm,
		       P7_SPARSEMX *asb, float *opt_sc)
{
  const float *tsc = gm->tsc;        // ptr into model transition score vector; enables TSC(k,s) macro shorthand 
  const float *rsc = NULL;	     // ptr into emission scores on current row i; enables MSC(k),ISC(k) macro shorthand 
  const float *rsn;		     // ptr into emission scores on next row i+1; enables MSN(k),ISN(k) 
  const float *last_rsc = NULL;      
  float       *dpc;		     // ptr stepping through current DP row i sparse cells 
  float       *dpn = NULL;	     // ptr stepping through next DP row (i+1) sparse cells; NULL when no i+1 row stored 
  float       *last_dpc;
  float       *xc;
  int          i,k,y,z;
  int          g;                    // index over sparse segments, 1..S
  int          d;                    // index of the next anchor we will see as we decrement i. d=D..1,(0); i0(d) <= i
  int          s;
  float        xC,xG,xL,xB,xJ,xN,xE; // tmp vars for special cells on this row i
  float        dlc,dgc;		     // tmp vars for D vals on current row i  
  float        mln,mgn,iln,ign;	     // tmp vars for M+e, I+e on next row i+1 
  float        mlc = 0.;             // initialization is solely to quiet static code analyzers
  float        mgc = 0.;             //  ... ditto.
  int          nup = 0;              // counter for how many rows into an UP sector we are
  int          ngap;                 // number of residues in an intersegment gap
  int64_t      asc_ncells;
  int          asc_nrows;
  int          status;

  /* Assure that <asb> is allocated large enough.
   */
  if ((status = p7_spascmx_Resize(asb, sm, anch, D)) != eslOK) return status;
  asb->type = p7S_ASC_BCK;

  /* TODO:
   * We need a more efficient way to find the end of the dp and xmx
   * sparse arrays in <asb>. We end up calling p7_spascmx_MinSizeof()
   * three times (each _Resize()) calls it), but we really only need
   * to do it once. But for now, this is fine; right now we're more
   * concerned with correctness of p7_sparse_asc_Backward() than
   * efficiency.
   */
  p7_spascmx_MinSizeof(sm, anch, D, &asc_ncells, &asc_nrows);

  dpc = asb->dp  + (asc_ncells - 1) * p7S_NSCELLS;
  xc  = asb->xmx + (asc_nrows  - 1) * p7S_NXCELLS;

  d   = D;
  xC  = gm->xsc[p7P_C][p7P_MOVE];   // initialize with C->T on i=L. Because this is C(L), this dictates that xC is a partial "precalculation" of value C(i), including tCC
  xJ  = -eslINFINITY;
  xN  = -eslINFINITY;
  for (g = sm->S; g >= 1; g--)
    {

      /* Starting a new segment, we're on ib(g) row.
       * Our convention in starting the i loop below is that:
       *     xC/J/N already include contribution of tCC/tJJ/tNN to reach row i from i+1
       *     xL/G have been accumulated on row i+1
       * And our convention for finishing the previous segment is that 
       *     xC/J/N include terms to reach ia(g+1)-1 (or L)
       *     
       * So in initializing for a new row ib, we must account for
       * transitions from ib(g) to ia(g+1)-1 (or L). That means
       * assessing a number of tXX transitions equal to the number of 
       * residues in the intersegment gap.
       */
      ngap = (g == sm->S ?  L - sm->seg[g].ib : sm->seg[g+1].ia - sm->seg[g].ib - 1);   // xC/J/N must account for ib->..->ia-1 transitions; ib+1..ia-1 residues
      xC +=  (ngap == 0  ?  0. : (float) ngap * gm->xsc[p7P_C][p7P_LOOP]);              // beware 0 * -inf = NaN; hence test on ngap == 0
      xJ +=  (ngap == 0  ?  0. : (float) ngap * gm->xsc[p7P_J][p7P_LOOP]);
      xN +=  (ngap == 0  ?  0. : (float) ngap * gm->xsc[p7P_N][p7P_LOOP]);     
      xL =   -eslINFINITY;                                                              // xL/xG remain -inf while in the last DOWN sector in the segment,
      xG =   -eslINFINITY;                                                              // until we reach that sector's anchor. After that, xL/xG will accumulate
                                                                                        // {LG}->Mi+1,k paths back from Mk's in UP row, or from anchor cell in DOWN.
      nup = 0;

      for (i = sm->seg[g].ib; i >= sm->seg[g].ia; i--)
	{
	  rsn = last_rsc;                     // MSN(),ISN() macros will now work for row i+1 
	  rsc = last_rsc = gm->rsc[dsq[i]];   // MSC(),ISC() macros will now work for row i

	  /*****************************************************************
	   * Specials first.
	   * We have to "always calculate, sometimes store"... but here, we've
	   * already partially calculated xC/J/N by including the tXX contribution
	   * that gets us back to this row i; and if we're in a segment with no
	   * anchors, that's the only contribution. So we can go directly to storage code:
	   */
	  if (nup || anch[d].i0 >= sm->seg[g].ia) // if there's at least one anchor in this segment...
	    {
	      ESL_DASSERT1(( xc >= asb->xmx ));

	      xc[p7S_CC] = -eslINFINITY;
	      xc[p7S_JJ] = -eslINFINITY;
	      xc[p7S_C]  = xC;                // C(i) can only be reached from C(i+1) tC->C (or tC->T for i=L), and we already included that.
	      xc[p7S_G]  = xG;
	      xc[p7S_L]  = xL;
	      xc[p7S_B]  = xB = p7_FLogsum( xL + gm->xsc[p7P_B][0],        xG + gm->xsc[p7P_B][1]);         // on i=ib segment start, evaluates to -inf, because xL=xG=-inf there
	      xc[p7S_J]  = xJ = p7_FLogsum( xJ, 		           xB + gm->xsc[p7P_J][p7P_MOVE]);  // on i=ib, evaluates to xJ (only path is J->J(i+1))
	      xc[p7S_N]  = xN = p7_FLogsum( xN, 	                   xB + gm->xsc[p7P_N][p7P_MOVE]);  // ditto xN
	      xc[p7S_E]  = xE = p7_FLogsum( xJ + gm->xsc[p7P_E][p7P_LOOP], xC + gm->xsc[p7P_E][p7P_MOVE]);  // on i=ib, i<L, E->{JC} are both possible; on i=L, evaluates to xC + tEC.
	      xc -= p7S_NXCELLS;
	    }

	  /*****************************************************************
	   * UP sector.
	   *****************************************************************/

	  /* Three cases to deal with:
           *   1. Initialization case, when i is the final (i.e. bottommost) row of an UP sector. 
           *      Here we have no next UP row; instead we connect bottom corner to
           *      DOWN anchor cell ML/MG, via mlc,mgc
           *
           *   2. Termination case, when i is an anchor row. 
           *      No explicit path can use anchor row in UP sector.
           *      Only the D's of a G->D1..Dk-1->Mk entry path can be used,
           *      but they aren't explicitly stored; they're unfolded in decoding, later.
           *      So: all values are -inf.
	   *
           *   3. Main recursion.
           *   
           * When anchors occur on adjacent rows, the UP matrix for
           * the first one consists of just a single row, and that
           * must be treated as the termination case; only an unfolded
           * G->DDDk-1 entry path can pass along such a row.
           *
	   * As we enter: 
	   *    dpc      is on last supercell of row i.
	   *    last_dpc happens to be on last supercell of previous row, but
	   *             we don't use it.
	   *    mlc,mgc  are values from the anchor cell diagonal i+1,k+1 from dpc
	   * 
	   * As we finish the row:
	   *   dpc        is at the last supercell of the next row we need to calculate.
	   *              That's either a DOWN sector on this same row i, or (if we're the first,
	   *              topmost UP sector in the segment) the previous row i-1 of this same 
	   *              UP sector.
	   *   last_dpc   has been set to the last supercell of row i, where <dpc> started.     
	   *   dpn        is at the last supercell of row i+1 of next UP or DOWN sector
	   *   xL, xG     have accumulated L->MLk, G->MGk entry probability sum
	   */

	  /* Termination case: on anchor row i0, all UP sector values are -inf.
	   * Must test for this first, because a row can be both first and terminal.
	   */
	  if (nup && i == anch[d].i0)
	    {
	      dpn      = last_dpc;
	      last_dpc = dpc;

	      xL  = xG = -eslINFINITY;
	      z   = sm->n[i]-1;   while (z >= 0 && sm->k[i][z]   >= anch[d+1].k0) z--;
	      y   = sm->n[i+1]-1; while (y >= 0 && sm->k[i+1][y] >= anch[d+1].k0) y--;
	      for (; z >= 0; z--) 
		{
		  for (s = 0; s < p7S_NSCELLS; s++)
		    dpc[s] = -eslINFINITY;
		  dpc -= p7S_NSCELLS;
		}
	      for (; y >= 0; y--)
		dpn -= p7S_NSCELLS;
	    }
	  /* Initialization case */
	  else if (nup == 1)  
	    {
	      last_dpc = dpc;
	      dpn      = NULL;

	      dlc = dgc = -eslINFINITY;
	      xL  = xG  = -eslINFINITY;
	      z = sm->n[i]-1;
	      while (z >= 0 && sm->k[i][z] >= anch[d+1].k0) z--;
	      if    (z >= 0 && sm->k[i][z] == anch[d+1].k0-1)             
		{
		  k = sm->k[i][z];
		  dpc[p7S_ML] =        mlc + TSC(p7P_MM, k);   // mlc/mgc were precalculated and held over from 
		  dpc[p7S_MG] =        mgc + TSC(p7P_MM, k);   //   when we calculated the anchor cell in DOWN;
		  dpc[p7S_IL] =        mlc + TSC(p7P_IM, k);   //   we already included the emission score there.
		  dpc[p7S_IG] =        mgc + TSC(p7P_IM, k);
		  dpc[p7S_DL] = dlc =  mlc + TSC(p7P_DM, k);
		  dpc[p7S_DG] = dgc =  mgc + TSC(p7P_DM, k);
		  xL          = dpc[p7S_ML] + MSC(k) + TSC(p7P_LM, k-1);
		  xG          = dpc[p7S_MG] + MSC(k) + TSC(p7P_GM, k-1);
		  z--;
		  dpc -= p7S_NSCELLS;
		}
	      
	      for ( ; z >= 0; z--)
		{
		  k = sm->k[i][z];
		  if (sm->k[i][z+1] != k+1) { dlc = dgc = -eslINFINITY; }

		  dpc[p7S_ML] = dlc + TSC(p7P_MD, k);
		  dpc[p7S_MG] = dgc + TSC(p7P_MD, k);
		  dpc[p7S_IL] = -eslINFINITY;
		  dpc[p7S_IG] = -eslINFINITY;
		  dpc[p7S_DL] = dlc = dlc + TSC(p7P_DD, k);
		  dpc[p7S_DG] = dgc = dgc + TSC(p7P_DD, k);
		  xL          = p7_FLogsum(xL, dpc[p7S_ML] + MSC(k) + TSC(p7P_LM, k-1));
		  xG          = p7_FLogsum(xG, dpc[p7S_MG] + MSC(k) + TSC(p7P_GM, k-1));
		  dpc -= p7S_NSCELLS;
		}
	    }



	  /* main recursion eqns for an UP sector row i */
	  else if (nup)
	    {
	      dpn      = last_dpc;
	      last_dpc = dpc;

	      dlc = dgc = -eslINFINITY;
	      xL  = xG  = -eslINFINITY;
	      z   = sm->n[i]-1;   while (z >= 0 && sm->k[i][z]   >= anch[d+1].k0) z--;
	      y   = sm->n[i+1]-1; while (y >= 0 && sm->k[i+1][y] >= anch[d+1].k0) y--;

	      for ( ; z >= 0; z--)
		{
		  k = sm->k[i][z];

		  /* Try to find and pick up mln, mgn from i+1,k+1 */
		  while (y >= 0 && sm->k[i+1][y]  > k+1) { y--; dpn -= p7S_NSCELLS; }  // if y is already on i+1,k+1, it doesn't move
		  if    (y >= 0 && sm->k[i+1][y] == k+1) {
		    mln = dpn[p7S_ML] + MSN(k+1);
		    mgn = dpn[p7S_MG] + MSN(k+1);
		  } else { mln = mgn = -eslINFINITY; }
		      
		  /* Try to pick up iln,ign from i+1,k */
		  while (y >= 0 && sm->k[i+1][y]  > k)  { y--; dpn -= p7S_NSCELLS; }  // if y is already on i+1,k, it doesn't move
		  if    (y >= 0 && sm->k[i+1][y] == k) {
		    iln = dpn[p7S_IL]; // + ISN(k), if it weren't fixed to 0
		    ign = dpn[p7S_IG]; // + ISN(k), ditto
		  } else { iln = ign = -eslINFINITY; }

		  /* Check if dlc,dgc have become invalid because there is no k+1 sparse cell */
		  if (z < sm->n[i]-1 && sm->k[i][z+1] != k+1) 
		    dlc = dgc = -eslINFINITY;

		  /* M(i,k) calculations need to use dgc,dlc before we
		   * change them, while they're still D(i,k+1). 
		   */
		  dpc[p7S_ML] = mlc = p7_FLogsum( p7_FLogsum(mln + TSC(p7P_MM, k),      // ML(i,k) =   ML(i+1,k+1) * t(k,MM) * e_M(k+1, x_i+1)      | mln = log[ ML(i+1,k+1) * e_M(k+1,x_i+1)]
							     iln + TSC(p7P_MI, k)),     //           + IL(i+1,k)   * t(k,MI) * e_I(k, x_i+1)        | iln = log[ IL(i+1,k)   * e_I(k,  x_i+1)]
						             dlc + TSC(p7P_MD, k));     //           + DL(i,  k+1) * t(k,DD)                        | dlc = DL(i,k+1), wrapped around from prev loop iteration. This has tricky boundary conditions at k=kbc, and k=kbc=M
		  dpc[p7S_MG] = mgc = p7_FLogsum( p7_FLogsum(mgn + TSC(p7P_MM, k),      // MG(i,k) is essentially the same recursion.
							     ign + TSC(p7P_MI, k)),     
						             dgc + TSC(p7P_MD, k));     
		  
		  /* Accumulate xG, xL as we sweep over the current row; 
		   * they will get used to initialize the next row we calculate (i-1). Note
		   * how we need emission scores on current row for this,
		   * and on next row for the mgn/mln calculations above.
		   */
		  xG = p7_FLogsum( xG, dpc[p7S_MG] + MSC(k) + TSC(p7P_GM, k-1));   // t(k-1,p7P_GM) = left wing retraction, entry at Mk, stored off-by-one at k-1 
		  xL = p7_FLogsum( xL, dpc[p7S_ML] + MSC(k) + TSC(p7P_LM, k-1));   // t(k-1,p7P_LM) = uniform local entry at Mk, off-by-one storage at k-1 

		  /* I(i,k) calculations */
		  dpc[p7S_IL] = p7_FLogsum( mln + TSC(p7P_IM, k), iln + TSC(p7P_II, k));
		  dpc[p7S_IG] = p7_FLogsum( mgn + TSC(p7P_IM, k), ign + TSC(p7P_II, k));
	  
		  /* D(i,k) calculations; 
		   * we can start storing results in <dpc> and decrementing;
		   * and we're done with dgc,dlc, so can update them for
		   * time thru the k loop.
		   * Dg has no exit to E; ...DGm->E paths have been wrapped into the right wing retraction in t(k,MGE).
		   */
		  dpc[p7S_DL] = dlc = p7_FLogsum(mln + TSC(p7P_DM, k), dlc  + TSC(p7P_DD, k));
		  dpc[p7S_DG] = dgc = p7_FLogsum(mgn + TSC(p7P_DM, k), dgc  + TSC(p7P_DD, k)); 

		  dpc -= p7S_NSCELLS;
		}
	      
	      // dpn may be on first supercell on row i+1, or -1 from that;
	      // we need it to be -1 for sure. If we go to a row i in a DOWN sector
	      // next, the DOWN sector code depends on dpn moving smoothly from this
	      // UP row to the end of the DOWN row.
	      while (y >= 0) { y--; dpn -= p7S_NSCELLS; }
	    }

	  /* DOWN sector: 
           *    Compute sparse supercells k0(d)..M 
	   *    Can only exit from DOWN, not enter; paths back from E are computed,
	   *    but not paths back to G,L.
	   *    Exception: can enter the anchor cell i0,k0.
	   */

	  // We depend on the following state of things:
	  //   dpc : currently points at the last sparse cell on DOWN row i


	  if ( anch[d].i0 >= sm->seg[g].ia)  // d may be 0 here, and if so, sentinel i0(0) makes the test fail
	    {

	      /* Last row of DOWN sector: initialization case. 
	       * No next row i+1. Cells only reachable from E, D 
	       */
	      if (nup == 1 || i == sm->seg[g].ib)  
		{                                  
		  if (! nup) last_dpc = dpc;
		  dpn = NULL;
		  dlc = -eslINFINITY;
		  dgc = xE + TSC(p7P_DGE, sm->k[i][sm->n[i]-1]);       // tDGE wing-retracted exit term needed on any Dk with no adjacent Dk+1.
		                                                       // tDGEk is Dk+1..Dm->E path sum; DGE[M] and [M-1] are 0
		  for (z = sm->n[i]-1; z >= 0 && sm->k[i][z] >= anch[d].k0; z--)
		    {
		      k = sm->k[i][z];

		      if (z < sm->n[i]-1 && sm->k[i][z+1] != k+1) {    // If supercell k has no adjacent k+1, reinitialize D paths
			dlc = -eslINFINITY;
			dgc = xE + TSC(p7P_DGE, k);
		      }

		      dpc[p7S_ML] = mlc = p7_FLogsum(dlc + TSC(p7P_MD, k), xE);
		      dpc[p7S_MG] = mgc =            dgc + TSC(p7P_MD, k);
		      dpc[p7S_IL] =       -eslINFINITY;
		      dpc[p7S_IG] =       -eslINFINITY;
		      dpc[p7S_DL] = dlc = p7_FLogsum(dlc + TSC(p7P_DD, k), xE);
		      dpc[p7S_DG] = dgc =            dgc + TSC(p7P_DD, k);        // Yes it's correct. Suppose k=M. DD term is DM->M+1 edge = 0. Initialization w/ DGE was also 0. Total = 0.
                                                                                  // Now suppose we start on an internal k. DD term is Dk->Dk+1; DGE is Dk+1..Dm->E; that's the complete Dk->E exit path

		      dpc -= p7S_NSCELLS;
		    } // end loop over z. Calculation of all sparse cells k on current DOWN row i completed.
                      // dpc now -1 from first DOWN sparse supercell on row i: i.e. last supercell of prev row (in up or down sector, depending)
		} // end code block for initialization of the last row of a DOWN sector

	      /* Any other row i of DOWN sector.
	       * We have a next row i+1. 
	       */
	      else  // We're interior in a DOWN sector, with at least one DOWN row below us at i+1, so all i+1 indices (in sm->n[] for example) are safe
                {   // Cells k = k0(d)..M must be computed.
		    // From a DOWN sector, paths can end from the model, but not start in it,
		    // so we pull {MD}->E transitions backward from xE, but we don't evaluate B->{LG}->Mk
		  if (! nup) {
		    dpn      = last_dpc;  // if we're the last DOWN sector in segment, we're the first access to dpn; set it.
		    last_dpc = dpc;      
		  }

		  dlc = -eslINFINITY;
		  dgc = xE + TSC(p7P_DGE, sm->k[i][sm->n[i]-1]);  // ?? 
		  y   = sm->n[i+1] - 1;   
		  for (z = sm->n[i]-1; z >= 0 && sm->k[i][z] >= anch[d].k0; z--)
		    {
		      k = sm->k[i][z];

		      /* Try to find and pick up mln, mgn from i+1,k+1 */
		      while (y >= 0 && sm->k[i+1][y]  > k+1) { y--; dpn -= p7S_NSCELLS; }  // if y is already on i+1,k+1, it doesn't move
		      if    (y >= 0 && sm->k[i+1][y] == k+1) {
			mln = dpn[p7S_ML] + MSN(k+1);
			mgn = dpn[p7S_MG] + MSN(k+1);
		      } else { mln = mgn = -eslINFINITY; }
		      
		      /* Try to pick up iln,ign from i+1,k */
		      while (y >= 0 && sm->k[i+1][y]  > k)  { y--; dpn -= p7S_NSCELLS; }  // if y is already on i+1,k, it doesn't move
		      if    (y >= 0 && sm->k[i+1][y] == k) {
			iln = dpn[p7S_IL]; // + ISN(k), if it weren't fixed to 0
			ign = dpn[p7S_IG]; // + ISN(k), ditto
		      } else { iln = ign = -eslINFINITY; }

		      /* Check if dlc,dgc are invalid and need to be reinitialized */
		      if (z < sm->n[i]-1 && sm->k[i][z+1] != k+1) {
			dlc = -eslINFINITY;
			dgc = xE + TSC(p7P_DGE, k);
		      }

		      /* M(i,k) calculations need to use dgc,dlc before we
		       * change them, while they're still D(i,k+1). 
		       */
		      dpc[p7S_ML] = mlc = p7_FLogsum( p7_FLogsum(mln + TSC(p7P_MM, k),      // ML(i,k) =   ML(i+1,k+1) * t(k,MM) * e_M(k+1, x_i+1)      | mln = log[ ML(i+1,k+1) * e_M(k+1,x_i+1)]
								 iln + TSC(p7P_MI, k)),     //           + IL(i+1,k)   * t(k,MI) * e_I(k, x_i+1)        | iln = log[ IL(i+1,k)   * e_I(k,  x_i+1)]
						      p7_FLogsum(dlc + TSC(p7P_MD, k),      //           + DL(i,  k+1) * t(k,DD)                        | dlc = DL(i,k+1), wrapped around from prev loop iteration. This has tricky boundary conditions at k=kbc, and k=kbc=M
								 xE));                      //           +  E(i)       * t(MkE)=1.0
		      dpc[p7S_MG] = mgc = p7_FLogsum( p7_FLogsum(mgn + TSC(p7P_MM, k),      // MG(i,k) is essentially the same recursion, without a transition to E.
								 ign + TSC(p7P_MI, k)),     
						                 dgc + TSC(p7P_MD, k));     

		      /* I(i,k) calculations */
		      dpc[p7S_IL] = p7_FLogsum( mln + TSC(p7P_IM, k), iln + TSC(p7P_II, k));
		      dpc[p7S_IG] = p7_FLogsum( mgn + TSC(p7P_IM, k), ign + TSC(p7P_II, k));

		      /* D(i,k) calculations; 
		       * we can start storing results in <dpc> and decrementing;
		       * and we're done with dgc,dlc, so can update them for
		       * time thru the k loop.
		       * Dg has no exit to E; ...DGm->E paths have been wrapped into the right wing retraction in t(k,MGE).
		       */
		      dpc[p7S_DL] = dlc = p7_FLogsum( p7_FLogsum(mln  + TSC(p7P_DM, k), 
								 dlc  + TSC(p7P_DD, k)),
						      xE);                   
		      dpc[p7S_DG] = dgc = p7_FLogsum( mgn  + TSC(p7P_DM, k),  
						      dgc  + TSC(p7P_DD, k)); 

		      dpc -= p7S_NSCELLS;

		    } // end loop over sparse cells k in DOWN sector.
		      // dpc is now one supercell beyond: i.e. on last supercell of UP
		} // end of code block when i is an interior DOWN row, where we have a valid i+1 row to look at
	      


	      /* The handoff: from anchor cell (top left corner of
               * DOWN sector) to bottom right corner of UP sector. All
               * paths must pass thru anchor cell; they either start
               * exactly there (via {LG}->M) or they continue back to
               * the lower right corner of UP. We precalculate
               * mlc/mgc, which will be used to initialize that UP
               * supercell when we get to the prev row i-1; we also
               * precalculate xL/xG for i-1 which can only come from
               * the anchor cell when we're in DOWN (as opposed to UP
               * sector supercells, all of which can go back to xL/xG.
               * Meanwhile (because paths have to go thru anchor cell)
               * we set C/J/N to -inf. J,N will be restarted on row
               * i-1 when we pull from L,G back to B back to J,N.
	       */
	      if (i == anch[d].i0) 
		{
		  mlc += MSC(anch[d].k0);               // now, don't touch mlc/mgc again until UP initialization of prev row
		  mgc += MSC(anch[d].k0);
		  xL = mlc + TSC(p7P_LM, anch[d].k0-1);
		  xG = mgc + TSC(p7P_GM, anch[d].k0-1);
		  xJ = -eslINFINITY;
		  xC = -eslINFINITY;
		  xN = -eslINFINITY;
		}
	    } // end of all DOWN sector calculations for row i
	  
	  /* Partial precalculation of xC/J/N for row i-1 we'll do next.
	   * When loop starts new row i, xC/xJ/xN must already include 
	   * tCC/JJ/NN transition that reaches that row
	   */
	  xC += gm->xsc[p7P_C][p7P_LOOP];
	  xJ += gm->xsc[p7P_J][p7P_LOOP];
	  xN += gm->xsc[p7P_N][p7P_LOOP];

	  if      (i == anch[d].i0)   { nup = 1; d--; }
	  else if (nup)               { nup++;        }
	} // end loop over rows i of a segment

      /* Finally the ia(g)-1 row specials */
      if (nup)
	{
	  xc[p7S_CC] = -eslINFINITY;
	  xc[p7S_JJ] = -eslINFINITY;
	  xc[p7S_C]  = xC;   
	  xc[p7S_G]  = xG;
	  xc[p7S_L]  = xL;
	  xc[p7S_B]  = xB = p7_FLogsum( xL + gm->xsc[p7P_B][0],        xG + gm->xsc[p7P_B][1]);         // on i=ib segment start, evaluates to -inf, because xL=xG=-inf there
	  xc[p7S_J]  = xJ = p7_FLogsum( xJ, 		           xB + gm->xsc[p7P_J][p7P_MOVE]);  // on i=ib, evaluates to xJ (only path is J->J(i+1))
	  xc[p7S_N]  = xN = p7_FLogsum( xN, 	                   xB + gm->xsc[p7P_N][p7P_MOVE]);  // ditto xN
	  xc[p7S_E]  = xE = p7_FLogsum( xJ + gm->xsc[p7P_E][p7P_LOOP], xC + gm->xsc[p7P_E][p7P_MOVE]);  // on i=ib, i<L, E->{JC} are both possible; on i=L, evaluates to xC + tEC.
	  xc -= p7S_NXCELLS;
	}

      last_rsc = NULL;
    } // end loop over segments g

  /* Account for 1..ia(1)-1 residues in NNNN... path */
  if ( sm->seg[1].ia > 1 )
    xN += (float) (sm->seg[1].ia - 1) * gm->xsc[p7P_N][p7P_LOOP];
  if (opt_sc) *opt_sc = xN;	 // S->N = 1.0, so no transition score needs to be added.
  return eslOK;
}
/*-------------- end, sparse ASC Backward -----------------------*/




/*****************************************************************
 * 4. Sparse ASC Decoding
 *****************************************************************/

// Access pattern deliberately mirrors that of Backward.
// Future optimization: merge w/ Backward
//  takes Fwd as arg - helps in knowing matrix size

// May need renormalization.

// Merge w/ fwdback. Make one file, asc.c, w/ all the ASC unit testing.
// sparse_fwdback,decoding ditto: call that one firstpass.c?

// Unlike other Decoding, this does NOT allow overwrite of asb.
//   don't worry - will not need it, after optimization by merging w/ Backwards

int
p7_sparse_asc_Decoding(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
		       const P7_ANCHOR *anch, int D,
		       float totsc, const P7_SPARSEMX *asf, P7_SPARSEMX *asb, P7_SPARSEMX *asd)
{
  const P7_SPARSEMASK *sm        = asf->sm;   // sparse mask, defining what cells are in the DP matrix
  const float         *dpcF;                  // ptr into... main supercells of Forward sparse ASC matrix on current row i
  const float         *xcF;                   // ... specials of Forward mx, current row i
  const float         *dpcB;                  // ... main supercells of Backward mx, current row i
  const float         *dpnB;                  // ... ... and next row, i+1; needed for G->Mk wing unfolding
  const float         *last_dpcB;             // ... ... and the start of each row; used for reinitializing <dpnB>
  const float         *xcB;                   // ... specials of Backward mx, curr row i
  float               *dpcD;                  // ... main supercells of Decoding matrix
  float               *xcD;                   // ... specials of Decoding matrix
  const float         *tsc       = gm->tsc;   // ptr into model transition score vector; enables TSC(k,s) macro shorthand 
  const float         *rsn;                   // ptr into emission scores for next row i+1; enables MSN() macro; needed for G->Mk wing unfolding
  float                xFG;                   // to remember Fwd G(i) value, for G->D...D->Mk wing unfolding
  float                delta;                 // to accumulate probability of G->D...D->Mk unfolded paths
  int                  nup;		      // counts UP sector rows; 1=bottommost
  int                  i;                     // index for dp matrix rows/sequence positions in dsq; 1..L
  int                  d;                     // index for anchors; 1..D
  int                  g;                     // index for sparse segments; 1..sm->S
  int                  z;                     // index for sparse supercells on a current row; 0..sm->n[i]-1
  int                  y;                     //  ... and for a next row i+1; 0..sm->n[i+1]-1
  int                  s;                     // index for states in a supercell
  int64_t asc_ncells;
  int     asc_nrows;
  int     status;


  /* Assure that <asd> is allocated large enough.
   */
  if ((status = p7_spascmx_Resize(asd, sm, anch, D)) != eslOK) return status;
  asd->type = p7S_ASC_DECODE;

  /* TODO:
   * Same as ASC Backward, 
   * we need a more efficient way to find the end of the dp and xmx
   * sparse arrays in <asd>. We end up calling p7_spascmx_MinSizeof()
   * three times (each _Resize()) calls it), but we really only need
   * to do it once. But for now, this is fine; right now we're more
   * concerned with correctness than efficiency.
   */
  p7_spascmx_MinSizeof(sm, anch, D, &asc_ncells, &asc_nrows);

  dpcF = asf->dp  + (asc_ncells - 1) * p7S_NSCELLS;
  dpcB = asb->dp  + (asc_ncells - 1) * p7S_NSCELLS;
  dpcD = asd->dp  + (asc_ncells - 1) * p7S_NSCELLS;

  xcF  = asf->xmx + (asc_nrows  - 1) * p7S_NXCELLS;
  xcB  = asb->xmx + (asc_nrows  - 1) * p7S_NXCELLS;
  xcD  = asd->xmx + (asc_nrows  - 1) * p7S_NXCELLS;

  d  = D;
  for (g = sm->S; g >= 1; g--)
    {
      if (anch[d].i0 < sm->seg[g].ia) continue;         // no anchor in seg? Then no sparse storage for this seg, skip it.
      
      nup = 0;
      for (i = sm->seg[g].ib; i >= sm->seg[g].ia; i--)
	{
	  ESL_DASSERT1(( dpcF >= asf->dp-p7S_NSCELLS ));   // we can run out of supercells before running out of rows. When we run out of supercells, dpc{FBD} are oob at -6.
	  ESL_DASSERT1(( dpcB >= asb->dp-p7S_NSCELLS ));  
	  ESL_DASSERT1(( dpcD >= asd->dp-p7S_NSCELLS ));  
	  ESL_DASSERT1(( xcF  >= asf->xmx ));
	  ESL_DASSERT1(( xcB  >= asb->xmx ));
	  ESL_DASSERT1(( xcD  >= asd->xmx ));
	  

	  /* Specials */
	  xcD[p7S_C]  = exp(xcF[p7S_C]  + xcB[p7S_C] - totsc);
	  xcD[p7S_G]  = exp(xcF[p7S_G]  + xcB[p7S_G] - totsc);
	  xcD[p7S_L]  = exp(xcF[p7S_L]  + xcB[p7S_L] - totsc);
	  xcD[p7S_B]  = exp(xcF[p7S_B]  + xcB[p7S_B] - totsc);
	  xcD[p7S_J]  = exp(xcF[p7S_J]  + xcB[p7S_J] - totsc);
	  xcD[p7S_N]  = exp(xcF[p7S_N]  + xcB[p7S_N] - totsc);
	  xcD[p7S_E]  = exp(xcF[p7S_E]  + xcB[p7S_E] - totsc);

	  xFG  = xcF[p7S_G];     // Remember Fwd G(i) score. We'll need it for G->Mk wing unfolding in UP sector rows i
	  xcF -= p7S_NXCELLS;    // Back xcF up early. Now we can decode J->J(i) emission, C->C(i) emission 
	  xcD[p7S_CC] = (i == anch[d].i0 ? 0. : exp(xcF[p7S_C] + gm->xsc[p7P_C][p7P_LOOP] + xcB[p7S_C] - totsc));  // i0 test because anchor i0 cannot be emitted by anything but Mk0
	  xcD[p7S_JJ] = (i == anch[d].i0 ? 0. : exp(xcF[p7S_J] + gm->xsc[p7P_J][p7P_LOOP] + xcB[p7S_J] - totsc));

	  xcB -= p7S_NXCELLS;
	  xcD -= p7S_NXCELLS;


	  /* UP sector */
	  if (nup == 1)
	    {
	      last_dpcB = dpcB;
	      rsn       = gm->rsc[dsq[i+1]];                                                                                // i+1 must exist in seq, because UP sector can't include i+1
	      delta     = exp(xFG + TSC(p7P_GM, anch[d+1].k0-1) + MSN(anch[d+1].k0) + dpcB[p7S_MG + p7S_NSCELLS] - totsc);  // tricksy: dpcB+p7S_NSCELLS happens to be the anchor cell i0,k0!

	      z = sm->n[i]-1;  while (z >= 0 && sm->k[i][z] >= anch[d+1].k0) z--;   // skip to last sparse supercell in UP row i
	      for (; z >= 0; z--)
		{
		  for (s = 0; s < p7S_NSCELLS; s++)
		    dpcD[s] = exp(dpcF[s] + dpcB[s] - totsc);
		  dpcD[p7S_DG] += delta;   // G->Mk0 wing unfolding: each sparse cell on last row of UP has a single special-case connection to MGk0(i0) anchor

		  dpcD -= p7S_NSCELLS;
		  dpcF -= p7S_NSCELLS;
		  dpcB -= p7S_NSCELLS;
		}
	    }
	  else if (nup)  // for all UP rows except the bottommost one; here we know row i+1 exists in UP
	    {
	      dpnB      = last_dpcB;
	      last_dpcB = dpcB;
	      rsn       = gm->rsc[dsq[i+1]]; 
	      delta     = 0.;

	      z = sm->n[i]-1;   while (z >= 0 && sm->k[i][z]   >= anch[d+1].k0) z--;
	      y = sm->n[i+1]-1; while (y >= 0 && sm->k[i+1][y] >= anch[d+1].k0) y--;  // n[i+1] exists, because we dealt with bottommost row (nup=1) as special case above
	      for ( ; z >= 0; z--)
		{
		  /* Update delta, backing up over all supercells on NEXT row i+1 such that j > k: i.e. such that G->DDDDk->Mj path crosses/includes Dk */
		  while (y >= 0 && sm->k[i+1][y] > sm->k[i][z]) {
		    delta += exp(xFG + TSC(p7P_GM, sm->k[i+1][y]-1) + MSN(sm->k[i+1][y]) + dpnB[p7S_MG] - totsc);
		    y--;
		    dpnB -= p7S_NSCELLS;
		  }
		  
		  for (s = 0; s < p7S_NSCELLS; s++)
		    dpcD[s] = exp(dpcF[s] + dpcB[s] - totsc);
		  dpcD[p7S_DG] += delta;  // 

		  dpcD  -= p7S_NSCELLS;
		  dpcF  -= p7S_NSCELLS;
		  dpcB  -= p7S_NSCELLS;
		}
	    }

	  
	  /* Is row i in a DOWN sector? */
	  if (anch[d].i0 >= sm->seg[g].ia) 
	    {
	      for (z = sm->n[i]-1; z >= 0 && sm->k[i][z] >= anch[d].k0; z--)
		{
		  for (s = 0; s < p7S_NSCELLS; s++)
		    dpcD[s] = exp(dpcF[s] + dpcB[s] - totsc);
		  dpcD  -= p7S_NSCELLS;
		  dpcF  -= p7S_NSCELLS;
		  dpcB  -= p7S_NSCELLS;
		}
	    }
	  
	  if      (i == anch[d].i0) { nup = 1; d--; }
	  else if (nup)             { nup++;        }
	} // end loop over rows i for one segment g
      
      /* Now we're on ia-1 specials. */
      xcD[p7S_CC] = 0.;                                    // There must be an anchor below us; C is now impossible
      xcD[p7S_JJ] = exp(xcF[p7S_J]  + xcB[p7S_J] - totsc); // E is impossible, so no E->J mass; all prob mass in J(ia-1) must come from ia-2, w/ emission at ia-1
      xcD[p7S_C]  = 0.;                                    // Don't bother calculating it, we know it's zero by construction.
      xcD[p7S_G]  = exp(xcF[p7S_G]  + xcB[p7S_G] - totsc); 
      xcD[p7S_L]  = exp(xcF[p7S_L]  + xcB[p7S_L] - totsc);
      xcD[p7S_B]  = exp(xcF[p7S_B]  + xcB[p7S_B] - totsc);
      xcD[p7S_J]  = xcD[p7S_JJ];
      xcD[p7S_N]  = exp(xcF[p7S_N]  + xcB[p7S_N] - totsc);
      xcD[p7S_E]  = 0.;                                    // Known by construction.
      xcF -= p7S_NXCELLS;
      xcB -= p7S_NXCELLS;
      xcD -= p7S_NXCELLS;
    } // end loop over segments g

  ESL_DASSERT1(( dpcF == asf->dp-p7S_NSCELLS  ));  
  ESL_DASSERT1(( dpcB == asb->dp-p7S_NSCELLS  ));  
  ESL_DASSERT1(( dpcD == asd->dp-p7S_NSCELLS  ));  
  ESL_DASSERT1(( xcF  == asf->xmx-p7S_NXCELLS ));  
  ESL_DASSERT1(( xcB  == asb->xmx-p7S_NXCELLS ));  
  ESL_DASSERT1(( xcD  == asd->xmx-p7S_NXCELLS ));  

  p7_spascmx_Renormalize(asd, anch, D);
  return eslOK;
}




/*****************************************************************
 * 5. Footnotes
 *****************************************************************
 *
 * [1] NUMERICAL ERROR MATCHED TO SPARSE FORWARD/BACKWARD
 * 
 * The sparse Forward and Backward algorithms (sparse_fwdback.c) skip
 * across intersegment intervals, where the path is forced to be
 * NN..NN/JJ..JJ/CC..CC, and reinitialize at every ia(g)-1 row by
 * adding ng*t_XX for <ng> residues in the interval. This is a small
 * speed optimization.
 * 
 * It also reduces numerical roundoff error. Note that np != \sum_n p
 * in floating point arithmetic. Multiplication is nigh-exact 
 * but summation accumulates error.
 * 
 * It's important for sparse ASC Forward/Backward to do exactly the
 * same calculation, because we will compare ASC Forward scores to
 * overall Forward scores to calculate the posterior probability of an
 * anchor set, and we can't tolerate numerical error giving us asc_f >
 * fsc, which results in posterior probs > 1. 
 * 
 * Hence we decouple *calculation* of the special states from the
 * *storage* of their values. The calculation mirrors (unconstrained)
 * sparse Forward/Backward, whereas storage follows the rules of ASC.
 * This works fine because ASC stores specials for a strict subset of
 * the segments that the unconstrained sparse algorithms store
 * (specifically, segments with >= 1 anchor in them).
 * 
 * Both p7_SparseForward() and p7_sparse_asc_Forward() are still
 * subject to numerical roundoff error accumulation, but they are
 * subject to *exactly the same* error accumulation, so we don't
 * get asc_f > fsc.
 * 
 * xref SRE:J13/150.
 *
 *------------------ end, footnotes -----------------------------*/



/*****************************************************************
 * 6. Unit tests
 *****************************************************************/
#ifdef p7SPARSE_ASC_FWDBACK_TESTDRIVE
#include "hmmer.h"

/* "generation" test.
 * 
 * Sample a profile, generate some synthetic homologous sequences from
 * it. For each seq, run the standard analysis pipeline: vector local checkpointed
 * decoding to get a sparse mask, and sparse segmental MPAS to get 
 * the anchor set. Finish the sparse ASC analysis w/ Backward, Decoding.
 *
 * Not a very stringent test, since we have no guarantees or constraints on paths or
 * scores: basically just a test for obvious crashes and invalid DP cells.
 * Advantage, though, is that we're generating a thorough and unconstrained sample
 * of profiles and sequences.
 *
 * Tests:
 *   - ASC fwd and bck scores equal (within tolerance)
 *   - F/B/D sparse ASC matrices pass Validation
 */
static void
utest_generation(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int L, int N)
{
  char           msg[] = "sparse_asc_fwdback :: generation unit test failed";
  P7_BG         *bg    = p7_bg_Create(abc);
  P7_PRIOR      *pri   = NULL;
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = p7_profile_Create(M, abc);
  P7_OPROFILE   *om    = p7_oprofile_Create(M, abc);
  ESL_SQ        *sq    = esl_sq_CreateDigital(abc);
  P7_CHECKPTMX  *cx    = p7_checkptmx_Create(M, 100, ESL_MBYTES(32));
  P7_SPARSEMASK *sm    = p7_sparsemask_Create(M, 100);
  P7_TRACE      *tr    = p7_trace_Create();
  P7_SPARSEMX   *sxf   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxd   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asf   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asb   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asd   = p7_sparsemx_Create(NULL);
  P7_ANCHORS    *anch  = p7_anchors_Create();
  P7_ANCHORS    *vanch = p7_anchors_Create();
  P7_ANCHORHASH *ah    = p7_anchorhash_Create();
  float         *wrk   = NULL;
  float          vsc, fsc, asc_f, asc_b;
  int            idx;
  char           errbuf[eslERRBUFSIZE];
  float          tol   = 0.01;

  
  if      (abc->type == eslAMINO) pri = p7_prior_CreateAmino();
  else if (abc->type == eslDNA)   pri = p7_prior_CreateNucleic();
  else if (abc->type == eslRNA)   pri = p7_prior_CreateNucleic();
  else                            pri = p7_prior_CreateLaplace(abc);

  /* Sample a random profile. Do it from prior, so it has good/realistic scores. */
  if ( p7_modelsample_Prior(rng, M, abc, pri, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)               != eslOK) esl_fatal(msg);
  if ( p7_oprofile_Convert(gm, om)                  != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      /* Sample a synthetic homologous sequence from profile, of mean length <L> */
      if ( p7_profile_SetLength(gm, L)  != eslOK) esl_fatal(msg);   /* config to generate mean length of L (length was probably reset by last emitted seq) */
      do {
	esl_sq_Reuse(sq);
	p7_ProfileEmit(rng, hmm, gm, bg, sq, NULL);
      } while (sq->n > L * 10); /* allow many domains, but keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */

      /* Set profile length models appropriately for new seq length  */
      if ( p7_bg_SetLength           (bg, sq->n) != eslOK) esl_fatal(msg);
      if ( p7_profile_SetLength      (gm, sq->n) != eslOK) esl_fatal(msg);
      if ( p7_oprofile_ReconfigLength(om, sq->n) != eslOK) esl_fatal(msg);

      /* Use vector local checkpointed F/B/D to set sparse mask */
      if ( p7_checkptmx_Reinit(cx, gm->M, sq->n)                             != eslOK) esl_fatal(msg);
      if ( p7_ForwardFilter (sq->dsq, sq->n, om, cx, /*fsc=*/NULL)           != eslOK) esl_fatal(msg);
      if ( p7_BackwardFilter(sq->dsq, sq->n, om, cx, sm, p7_SPARSIFY_THRESH) != eslOK) esl_fatal(msg);

      /* First pass sparse analysis */
      if ( p7_SparseViterbi (sq->dsq, sq->n, gm, sm,  sxf, tr, &vsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward (sq->dsq, sq->n, gm, sm,  sxf,     &fsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward(sq->dsq, sq->n, gm, sm,  sxd,     NULL) != eslOK) esl_fatal(msg);
      if ( p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxd,     sxd)  != eslOK) esl_fatal(msg);

      /* MPAS to get the anchor set and sparse ASC Forward */
      if ( p7_sparse_anchors_SetFromTrace(sxd, tr, vanch) != eslOK) esl_fatal(msg);
      if ( p7_trace_Reuse(tr)                             != eslOK) esl_fatal(msg);
      if ( p7_sparse_Anchors(rng, sq->dsq, sq->n, gm,
			     vsc, fsc, sxf, sxd, vanch,
			     tr, &wrk, ah, 
			     asf, anch, &asc_f, 
			     NULL)     != eslOK) esl_fatal(msg);

      if ( p7_anchors_Validate(anch, sq->n, gm->M, errbuf)    != eslOK) esl_fatal("%s\n   %s\n", msg, errbuf);

      /* Remainder of sparse ASC analysis */
      if ( p7_sparse_asc_Backward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asb, &asc_b)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Decoding(sq->dsq, sq->n, gm, anch->a, anch->D, asc_f, asf, asb, asd) != eslOK) esl_fatal(msg);

      /* Tests */
      if (! diagfp)
	{
	  if ( p7_spascmx_Validate(asf, anch->a, anch->D, errbuf) != eslOK) esl_fatal("%s\n   %s\n", msg, errbuf);
	  if ( p7_spascmx_Validate(asb, anch->a, anch->D, errbuf) != eslOK) esl_fatal("%s\n   %s\n", msg, errbuf);
	  if ( p7_spascmx_Validate(asd, anch->a, anch->D, errbuf) != eslOK) esl_fatal("%s\n   %s\n", msg, errbuf);

	  if ( esl_FCompareAbs(asc_f, asc_b, tol) != eslOK) esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%20g\n", asc_f - asc_b);

      p7_anchors_Reuse(vanch);
      p7_anchors_Reuse(anch);
      p7_anchorhash_Reuse(ah);
      p7_sparsemx_Reuse(asf);
      p7_sparsemx_Reuse(asb);
      p7_sparsemx_Reuse(asd);
      p7_sparsemx_Reuse(sxf);
      p7_sparsemx_Reuse(sxd);
      p7_trace_Reuse(tr);
      p7_sparsemask_Reuse(sm);
      esl_sq_Reuse(sq);
    }


  if (wrk) free(wrk);
  p7_anchors_Destroy(vanch);
  p7_anchors_Destroy(anch);
  p7_anchorhash_Destroy(ah);
  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(asb);
  p7_sparsemx_Destroy(asd);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxd);
  p7_trace_Destroy(tr);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(cx);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_prior_Destroy(pri);
  p7_bg_Destroy(bg);
}

/* "compare_reference" unit test.
 * 
 * When we include all i,k supercells in the sparse mask, then sparse
 * ASC DP calculations give the same results as reference ASC
 * calculations, even at the level of individual DP cells.
 * 
 * Sample a random profile of length <M>.  Generate <N> sequences from
 * that profile, using a length model of <L> during sampling.  For
 * each sampled sequence, make a sparse mask that contains all i,k
 * cells. Make an anchor set from the generating trace (the anchor set
 * just needs to be reasonable, not optimal). Do both reference ASC
 * and sparse ASC Fwd, Bck, Decoding.
 * 
 * Then:
 *   1. Reference and sparse Fwd, Bck scores must be identical (within 
 *      numerical tolerance).
 *   2. Cells of reference and sparse DP matrices (F, B, D) have identical values
 *      (within numerical tolerance).
 *   3. Sparse ASC matrix structures (F, B, D) all pass Validate().
 */
static void
utest_compare_reference(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int L, int N)
{
  char           msg[] = "sparse_asc_fwdback :: compare_reference unit test failed";
  P7_BG         *bg    = p7_bg_Create(abc);
  P7_HMM        *hmm   = NULL;
  ESL_SQ        *sq    = esl_sq_CreateDigital(abc);
  P7_PROFILE    *gm    = p7_profile_Create(M, abc);
  P7_TRACE      *gtr   = p7_trace_Create();
  P7_SPARSEMASK *sm    = p7_sparsemask_Create(M, L);
  P7_ANCHORS    *anch  = p7_anchors_Create();
  P7_REFMX      *afu   = p7_refmx_Create(100,100);
  P7_REFMX      *afd   = p7_refmx_Create(100,100);
  P7_REFMX      *abu   = p7_refmx_Create(100,100);
  P7_REFMX      *abd   = p7_refmx_Create(100,100);
  P7_REFMX      *apu   = p7_refmx_Create(100,100);
  P7_REFMX      *apd   = p7_refmx_Create(100,100);
  P7_SPARSEMX   *asf   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asb   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asd   = p7_sparsemx_Create(NULL);
  float          fsc_r, bsc_r;
  float          fsc, bsc;
  int            idx;
  char           errbuf[eslERRBUFSIZE];
  float          tol   = 0.01;

  /* Sample a profile. Config as usual, multihit dual-mode local/glocal. */
  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      /* Generate (sample) a sequence from the profile */
      if ( p7_profile_SetLength(gm, L)  != eslOK) esl_fatal(msg);   /* config to generate mean length of L (length was probably reset by last emitted seq) */
      do {
	esl_sq_Reuse(sq);
	p7_ProfileEmit(rng, hmm, gm, bg, sq, gtr);
      } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
      if ( p7_profile_SetLength(gm, sq->n) != eslOK) esl_fatal(msg);

      /* Mark all cells in sparse mask */
      if ( p7_sparsemask_Reinit(sm, M, sq->n) != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_AddAll(sm)           != eslOK) esl_fatal(msg);

      //p7_trace_DumpAnnotated(stdout, gtr, gm, sq->dsq);
      
      /* Use generating trace to create a plausible anchor set */
      if ( p7_trace_Index(gtr)                        != eslOK) esl_fatal(msg);
      if ( p7_anchors_SampleFromTrace(anch, rng, gtr) != eslOK) esl_fatal(msg);

      //p7_anchors_Dump(stdout, anch);

      /* reference ASC calculations */
      if ( p7_ReferenceASCForward (sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, &fsc_r)             != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch->a, anch->D, abu, abd, &bsc_r)             != eslOK) esl_fatal(msg);
      if ( p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, abu, abd, apu, apd) != eslOK) esl_fatal(msg);
      
      /* sparse ASC calculations */
      if ( p7_sparse_asc_Forward (sq->dsq, sq->n, gm, anch->a, anch->D, sm, asf, &fsc)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Backward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asb, &bsc)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Decoding(sq->dsq, sq->n, gm, anch->a, anch->D, fsc, asf, asb, asd) != eslOK) esl_fatal(msg);

      //p7_refmx_Dump(stdout, abu);
      //p7_refmx_Dump(stdout, abd);
      //p7_spascmx_Dump(stdout, asb, anch->a, anch->D);

      /* comparisons */
      if (! diagfp) 
	{
	  if ( p7_spascmx_CompareReference(asf, anch->a, anch->D, afu, afd, tol) != eslOK) esl_fatal(msg); // Do Compare()'s  first. Makes debugging easier.
	  if ( p7_spascmx_CompareReference(asb, anch->a, anch->D, abu, abd, tol) != eslOK) esl_fatal(msg); // Better to fail on bad cell in DP matrix, then just see bad score diff.
	  if ( p7_spascmx_CompareReference(asd, anch->a, anch->D, apu, apd, tol) != eslOK) esl_fatal(msg); 

	  if ( p7_spascmx_Validate(asf, anch->a, anch->D, errbuf) != eslOK) esl_fatal("%s\n   %s\n", msg, errbuf);
	  if ( p7_spascmx_Validate(asb, anch->a, anch->D, errbuf) != eslOK) esl_fatal("%s\n   %s\n", msg, errbuf);
	  if ( p7_spascmx_Validate(asd, anch->a, anch->D, errbuf) != eslOK) esl_fatal("%s\n   %s\n", msg, errbuf);

	  if ( esl_FCompareAbs(fsc, fsc_r, tol) != eslOK) esl_fatal(msg);
	  if ( esl_FCompareAbs(bsc, bsc_r, tol) != eslOK) esl_fatal(msg);
	  if ( esl_FCompareAbs(fsc, bsc,   tol) != eslOK) esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%20g %20g %20g\n", fsc_r-fsc, bsc_r-bsc, fsc-bsc);

      p7_sparsemask_Reuse(sm);
      p7_anchors_Reuse(anch);
      p7_trace_Reuse(gtr);
      p7_refmx_Reuse(afu);    p7_refmx_Reuse(afd);   
      p7_refmx_Reuse(abu);    p7_refmx_Reuse(abd);
      p7_refmx_Reuse(apu);    p7_refmx_Reuse(apd);
      p7_sparsemx_Reuse(asf); 
      p7_sparsemx_Reuse(asb);
      p7_sparsemx_Reuse(asd);
      esl_sq_Reuse(sq);
    }

  p7_sparsemx_Destroy(asf); 
  p7_sparsemx_Destroy(asb);
  p7_sparsemx_Destroy(asd);
  p7_refmx_Destroy(afu);    p7_refmx_Destroy(afd);
  p7_refmx_Destroy(abu);    p7_refmx_Destroy(abd);
  p7_refmx_Destroy(apu);    p7_refmx_Destroy(apd);
  p7_anchors_Destroy(anch);
  p7_trace_Destroy(gtr);
  p7_sparsemask_Destroy(sm);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_sq_Destroy(sq);
}



/* "singlesingle"
 * 
 * The "singlesingle" path uses p7_modelsample_SinglePathedSeq()
 * to sample a contrived profile/seq pair that has only a one possible
 * path with probability 1: that is, P(\pi | x, H) = 1.
 * 
 * This is restricted to a single domain and glocal alignment.  Local
 * alignment and multihit paths are not tested.
 * 
 * To get an anchor set, we randomly choose an anchor MG state in the
 * trace's single domain.
 * 
 * To get a sparse mask, we use p7_sparsemask_SetFromTrace() to
 * include all of the trace's cells, plus a randomized scattering of
 * others.
 * 
 * Because only a single path is possible, it doesn't matter whether
 * we are sparse or reference, ASC or unconstrained, Viterbi or
 * Forward: all scores are the same. We don't have to test all
 * possible combinations, because we use the "singlesingle" style test
 * elsewhere. Here it suffices to compare to the trace score.
 *
 * Thus:
 *    1. Sparse ASC fwd score = trace score
 *    2. Sparse ASC bck score, ditto.
 *    3. Sparse ASC decoding matrix has 1.0 where trace passes thru.
 *    
 *****************************************************************
 * Analysis of stochastic error:
 *
 * tsc-fsc difference is exactly zero. Order of summation is exactly
 * the same.
 */
static void
utest_singlesingle(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int N)
{
  char           msg[] = "sparse_asc_fwdback :: singlesingle unit test failed";
  P7_BG         *bg    = p7_bg_Create(abc);
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = NULL;
  ESL_DSQ       *dsq   = NULL;
  P7_TRACE      *tr    = NULL;
  P7_ANCHOR     *anch  = NULL;
  P7_SPARSEMASK *sm    = NULL;
  P7_SPARSEMX   *asf   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asb   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asd   = p7_sparsemx_Create(NULL);
  P7_REFMX      *rxd   = p7_refmx_Create(100, 100);
  int           D,L;
  int           idx;
  float         tsc, fsc, bsc;
  float         tol = 0.0001;
  char          errbuf[eslERRBUFSIZE];

  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_SinglePathedSeq(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &tsc) != eslOK) esl_fatal(msg);

      if ((sm = p7_sparsemask_Create(M, L))         == NULL)  esl_fatal(msg);
      if ( p7_sparsemask_SetFromTrace(sm, rng, tr)  != eslOK) esl_fatal(msg);

      if ( p7_sparse_asc_Forward (dsq, L, gm, anch, D, sm, asf, &fsc)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Backward(dsq, L, gm, anch, D, sm, asb, &bsc)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Decoding(dsq, L, gm, anch, D, fsc, asf, asb, asd) != eslOK) esl_fatal(msg);

      /* To check that ASC decoding left us with 1.0's in cells along
       * <tr> and 0.0 elsewhere, we use some trickery. We count the
       * trace into a reference matrix.  Because we only count 1 trace
       * in, now <rxd> elements are only 0.0, 1.0 as we need them; no
       * normalization is needed. Then we call the
       * spascmx_CompareReference() routine. Even though that routine
       * is designed for comparison to reference ASC (with two
       * matrices UP and DOWN), it will work in this case when we pass
       * the same matrix <rxd> for both.
       */
      if ( p7_refmx_GrowTo   (rxd, gm->M, L)               != eslOK) esl_fatal(msg);      
      if ( p7_refmx_SetType  (rxd, gm->M, L, p7R_DECODING) != eslOK) esl_fatal(msg);
      if ( p7_refmx_SetValues(rxd, 0.0)                    != eslOK) esl_fatal(msg);
      if ( p7_refmx_CountTrace(tr, rxd)                    != eslOK) esl_fatal(msg);  

      //p7_trace_DumpAnnotated(stdout, tr, gm, dsq);
      //p7_spascmx_Dump(stdout, asf, anch, D);

      if (!diagfp)
	{
	  if ( p7_spascmx_CompareReference(asd, anch, D, rxd, rxd, tol) != eslOK) esl_fatal(msg); // TODO: should tol be 0.0?

	  if ( p7_spascmx_Validate(asf, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf);
	  if ( p7_spascmx_Validate(asb, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf); 
	  if ( p7_spascmx_Validate(asd, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf); 

	  if (esl_FCompareAbs(tsc, fsc, tol) != eslOK) esl_fatal(msg);
	  if (esl_FCompareAbs(tsc, bsc, tol) != eslOK) esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%20g %20g\n", tsc-fsc, tsc-bsc);

      p7_sparsemx_Reuse(asf);
      p7_sparsemx_Reuse(asb);
      p7_sparsemx_Reuse(asd);
      p7_refmx_Reuse(rxd);
      free(anch);
      p7_trace_Destroy(tr);
      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
      p7_sparsemask_Destroy(sm);
      free(dsq);
    }

  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(asb);
  p7_sparsemx_Destroy(asd);
  p7_refmx_Destroy(rxd);
  p7_bg_Destroy(bg);
}


/* "multisingle" unit test
 * 
 * Sample a contrived profile/seq/anchor triplet such that all
 * possible paths must pass thru the anchor. Thus, standard and ASC
 * F/B give the same answer. 
 * 
 * In the "multisingle" version of this idea, we have to disallow
 * local paths and multiple domains, but we allow N/J/C emission and
 * insert states. The contrived model is uniglocal, with one chosen
 * match state k0 that emits anchor residue X with probability 1 and
 * cannot emit any other residue. The contrived sequence has only a
 * single X in it. Because the model is forced to visit Mk0 exactly
 * once (because it's uniglocal), that single X must be aligned to
 * Mk0.
 * 
 * The sparse mask is constructed in two different ways. 
 *  (1) By using p7_sparsemask_SetFromTrace() to include the
 *      generating trace's cells, plus a randomized scattering
 *      of others.
 *  (2) By using local checkpointed F/B/Decoding as in normal
 *      HMMER pipeline.
 *      
 * Runs <N> different <gm>,<dsq> sampled comparisons, alternating
 * which of the two sparse mask construction alternatives it uses,
 * starting with SetFromTrace().
 * 
 * Then:
 *     1. Sparse Fwd score = Sparse ASC Fwd score.
 *     2. ditto for Bck scores
 *     3. Viterbi score >= sampled trace score
 *     4. Fwd/Bck scores >= Viterbi score
 */
static void
utest_multisingle(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int N)
{
  char           msg[] = "sparse_asc_fwdback :: multisingle unit test failed";
  P7_BG         *bg    = p7_bg_Create(abc);
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = NULL;
  ESL_DSQ       *dsq   = NULL;
  int            L;
  P7_TRACE      *gtr   = NULL;
  P7_ANCHOR     *anch  = NULL;
  int            D;
  float          gsc;
  P7_OPROFILE   *om  = NULL;
  P7_CHECKPTMX  *cx  = NULL;
  P7_SPARSEMASK *sm  = NULL;
  P7_SPARSEMX   *sxf = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxb = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asf = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asb = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asd = p7_sparsemx_Create(NULL);
  int            idx;
  float          vsc, fsc, bsc, asc_f, asc_b;
  char           errbuf[eslERRBUFSIZE];
  float          tol = 0.002;

  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_AnchoredUni(rng, M, bg, &hmm, &gm, &dsq, &L, &gtr, &anch, &D, &gsc) != eslOK) esl_fatal(msg);

      if ((sm = p7_sparsemask_Create(gm->M, L)) == NULL) esl_fatal(msg);
      
      if (idx%2) {
	if ((om = p7_oprofile_Create(hmm->M, abc)) == NULL)  esl_fatal(msg);
	if ( p7_oprofile_Convert(gm, om)           != eslOK) esl_fatal(msg);
	if ( p7_oprofile_ReconfigMultihit(om, L)   != eslOK) esl_fatal(msg);
	om->mode = p7_LOCAL;

	if ((cx = p7_checkptmx_Create(hmm->M, L, ESL_MBYTES(32)))      == NULL)  esl_fatal(msg);
	if ( p7_ForwardFilter (dsq, L, om, cx, /*fsc=*/NULL)           != eslOK) esl_fatal(msg);
	if ( p7_BackwardFilter(dsq, L, om, cx, sm, p7_SPARSIFY_THRESH) != eslOK) esl_fatal(msg);

	p7_oprofile_Destroy(om);
	p7_checkptmx_Destroy(cx);
      } else
	if ( p7_sparsemask_SetFromTrace(sm, rng, gtr) != eslOK) esl_fatal(msg);


      if ( p7_SparseViterbi      (dsq, L, gm,          sm, sxf, NULL, &vsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward      (dsq, L, gm,          sm, sxf,       &fsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward     (dsq, L, gm,          sm, sxb,       &bsc) != eslOK) esl_fatal(msg);

      if ( p7_sparse_asc_Forward (dsq, L, gm, anch, D, sm, asf,     &asc_f)   != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Backward(dsq, L, gm, anch, D, sm, asb,     &asc_b)   != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Decoding(dsq, L, gm, anch, D, asc_f, asf, asb, asd)  != eslOK) esl_fatal(msg);
      
      if (! diagfp)
	{
	  if ( p7_spascmx_Validate(asf, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf);
	  if ( p7_spascmx_Validate(asb, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf); 
	  if ( p7_spascmx_Validate(asd, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf); 

	  if (esl_FCompareAbs( fsc, asc_f, tol) != eslOK) esl_fatal(msg);
	  if (esl_FCompareAbs( fsc, asc_b, tol) != eslOK) esl_fatal(msg);
	  if (gsc > vsc+tol)                              esl_fatal(msg);
	  if (vsc > fsc)                                  esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%20g %20g %20g %20g %20g\n", fsc-asc_f, fsc-asc_b, bsc-asc_b, vsc-gsc, fsc-vsc);
  
      p7_sparsemx_Reuse(asf);
      p7_sparsemx_Reuse(asb);
      p7_sparsemx_Reuse(asd);
      p7_sparsemx_Reuse(sxf);
      p7_sparsemx_Reuse(sxb);

      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
      free(dsq);
      p7_trace_Destroy(gtr);
      free(anch);
      p7_sparsemask_Destroy(sm);
    }

  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(asb);
  p7_sparsemx_Destroy(asd);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_bg_Destroy(bg);
}


/* "multipath_local" test
 * 
 * Like multisingle, we contrive a profile/seq/anchorset triplet 
 * for which all paths must pass through the anchor set, so standard
 * and ASC F/B give the same score.
 * 
 * This variant of the multi* tests local alignment transitions.  In
 * order to do that, it disallows N/C/J emissions, insert states, and
 * multihit.
 * 
 * The contrivance is that only M states can generate residues (N/C/J
 * and I are disallowed). We choose a special Mk0 state in the model
 * as the anchor, and only this state can generate a special residue X
 * (chosen randomly) with nonzero probability. Meanwhile the target
 * sequence has exactly 1 X residue at position i0. Mk0/x_i0 are
 * therefore forced to align.
 * 
 * As in multisingle, the sparse mask is constructed either by using
 * p7_sparsemask_SetFromTrace(), or by using local checkpointed
 * F/B/Decoding as in normal HMMER pipeline.
 *      
 * Runs <N> different <gm>,<dsq> sampled comparisons, alternating
 * which of the two sparse mask construction alternatives it uses,
 * starting with SetFromTrace().
 * 
 * Then:
 *   1. Fwd score = ASC Fwd score
 *   2. Bck score = ASC Bck score
 *   3. Fwd/Bck scores >= Viterbi score
 */
static void
utest_multipath_local(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int N)
{
  char           msg[] = "sparse_asc_fwdback :: multipath_local unit test failed";
  P7_BG         *bg    = p7_bg_Create(abc);
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = NULL;
  ESL_DSQ       *dsq   = NULL;
  int            L;
  P7_TRACE      *gtr   = NULL;
  P7_ANCHOR     *anch  = NULL;
  int            D;
  float          gsc;
  P7_OPROFILE   *om  = NULL;
  P7_CHECKPTMX  *cx  = NULL;
  P7_SPARSEMASK *sm  = NULL;
  P7_SPARSEMX   *sxf = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxb = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asf = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asb = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asd = p7_sparsemx_Create(NULL);
  int            idx;
  float          vsc, fsc, bsc, asc_f, asc_b;
  char           errbuf[eslERRBUFSIZE];
  float          tol = 0.0001;

  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_AnchoredLocal(rng, M, bg, &hmm, &gm, &dsq, &L, &gtr, &anch, &D, &gsc) != eslOK) esl_fatal(msg);

      if ((sm = p7_sparsemask_Create(gm->M, L)) == NULL) esl_fatal(msg);
      
      if (idx%2) {
	if ((om = p7_oprofile_Create(hmm->M, abc)) == NULL)  esl_fatal(msg);
	if ( p7_oprofile_Convert(gm, om)           != eslOK) esl_fatal(msg);
	if ( p7_oprofile_ReconfigMultihit(om, 0)   != eslOK) esl_fatal(msg);
	om->mode = p7_LOCAL;

	if ((cx = p7_checkptmx_Create(hmm->M, L, ESL_MBYTES(32)))      == NULL)  esl_fatal(msg);
	if ( p7_ForwardFilter (dsq, L, om, cx, /*fsc=*/NULL)           != eslOK) esl_fatal(msg);
	if ( p7_BackwardFilter(dsq, L, om, cx, sm, p7_SPARSIFY_THRESH) != eslOK) esl_fatal(msg);

	p7_oprofile_Destroy(om);
	p7_checkptmx_Destroy(cx);
      } else 
	if ( p7_sparsemask_SetFromTrace(sm, rng, gtr) != eslOK) esl_fatal(msg);
  
      if ( p7_SparseViterbi      (dsq, L, gm,          sm, sxf, NULL, &vsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward      (dsq, L, gm,          sm, sxf,       &fsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward     (dsq, L, gm,          sm, sxb,       &bsc) != eslOK) esl_fatal(msg);

      if ( p7_sparse_asc_Forward (dsq, L, gm, anch, D, sm, asf,     &asc_f)  != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Backward(dsq, L, gm, anch, D, sm, asb,     &asc_b)  != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Decoding(dsq, L, gm, anch, D, asc_f, asf, asb, asd) != eslOK) esl_fatal(msg);
   
      if (! diagfp)
	{
	  if ( p7_spascmx_Validate(asf, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf);
	  if ( p7_spascmx_Validate(asb, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf); 
	  if ( p7_spascmx_Validate(asd, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf); 

	  if (esl_FCompareAbs( fsc, asc_f, tol) != eslOK) esl_fatal(msg);
	  if (esl_FCompareAbs( bsc, asc_b, tol) != eslOK) esl_fatal(msg);
	  if (gsc > vsc+tol)                              esl_fatal(msg);
	  if (vsc > fsc)                                  esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%20g %20g %20g %20g %20g\n", fsc-asc_f, bsc-asc_b, asc_f-asc_b, vsc-gsc, fsc-vsc);
  
      p7_sparsemx_Reuse(asf);
      p7_sparsemx_Reuse(asb);
      p7_sparsemx_Reuse(asd);
      p7_sparsemx_Reuse(sxf);
      p7_sparsemx_Reuse(sxb);

      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
      free(dsq);
      p7_trace_Destroy(gtr);
      free(anch);
      p7_sparsemask_Destroy(sm);
    }

  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(asb);
  p7_sparsemx_Destroy(asd);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_bg_Destroy(bg);
}

/* "multimulti" test
 * 
 * All three "multi*" tests contrive a profile/seq/anchorset triplet
 * for which all paths must pass through the anchor set. Thus, standard
 * and ASC F/B give the same scores.
 * 
 * In the "multimulti" variant, multihit is allowed; the price is 
 * disallowing N/C/J emissions, I states, and local alignment paths.
 * The model is L=0 multiglocal.
 * 
 * The trick is as follows. Choose a special residue X. Only the
 * anchor state Mk0 can emit this residue with nonzero probability.
 * The target sequence has exactly D X residues, one per domain.
 * Since all D of them must be accounted for, and the only way to
 * align to them is with Mk0, all paths must visit the D anchors.
 * 
 * As in multisingle, the sparse mask is constructed either by using
 * p7_sparsemask_SetFromTrace(), or by using local checkpointed
 * F/B/Decoding as in normal HMMER pipeline.
 *      
 * Runs <N> different <gm>,<dsq> sampled comparisons, alternating
 * which of the two sparse mask construction alternatives it uses,
 * starting with SetFromTrace().
 * 
 * Then:
 *   1. Fwd score = ASC Fwd score
 *   2. Bck score = ASC Bck score
 *   3. Fwd/Bck scores >= Viterbi score
 */
static void
utest_multimulti(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int N)
{
  char           msg[] = "sparse_asc_fwdback :: multimulti unit test failed";
  P7_BG         *bg    = p7_bg_Create(abc);
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = NULL;
  ESL_DSQ       *dsq   = NULL;
  int            L;
  P7_TRACE      *gtr   = NULL;
  P7_ANCHOR     *anch  = NULL;
  int            D;
  float          gsc;
  P7_OPROFILE   *om  = NULL;
  P7_CHECKPTMX  *cx  = NULL;
  P7_SPARSEMASK *sm  = NULL;
  P7_SPARSEMX   *sxf = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxb = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asf = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asb = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asd = p7_sparsemx_Create(NULL);
  int            idx;
  float          vsc, fsc, bsc, asc_f, asc_b;
  char           errbuf[eslERRBUFSIZE];
  float          tol = 0.01;

  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_AnchoredMulti(rng, M, bg, &hmm, &gm, &dsq, &L, &gtr, &anch, &D, &gsc) != eslOK) esl_fatal(msg);

      if ((sm = p7_sparsemask_Create(gm->M, L)) == NULL) esl_fatal(msg);
      
      if (idx%2) {
	if ((om = p7_oprofile_Create(hmm->M, abc)) == NULL)  esl_fatal(msg);
	if ( p7_oprofile_Convert(gm, om)           != eslOK) esl_fatal(msg);
	if ( p7_oprofile_ReconfigMultihit(om, 0)   != eslOK) esl_fatal(msg);
	om->mode = p7_LOCAL;

	if ((cx = p7_checkptmx_Create(hmm->M, L, ESL_MBYTES(32)))      == NULL)  esl_fatal(msg);
	if ( p7_ForwardFilter (dsq, L, om, cx, /*fsc=*/NULL)           != eslOK) esl_fatal(msg);
	if ( p7_BackwardFilter(dsq, L, om, cx, sm, p7_SPARSIFY_THRESH) != eslOK) esl_fatal(msg);

	p7_oprofile_Destroy(om);
	p7_checkptmx_Destroy(cx);
      } else 
	if ( p7_sparsemask_SetFromTrace(sm, rng, gtr) != eslOK) esl_fatal(msg);
  
      if ( p7_SparseViterbi      (dsq, L, gm,          sm, sxf, NULL, &vsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward      (dsq, L, gm,          sm, sxf,       &fsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward     (dsq, L, gm,          sm, sxb,       &bsc) != eslOK) esl_fatal(msg);

      if ( p7_sparse_asc_Forward (dsq, L, gm, anch, D, sm, asf,     &asc_f)  != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Backward(dsq, L, gm, anch, D, sm, asb,     &asc_b)  != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Decoding(dsq, L, gm, anch, D, asc_f, asf, asb, asd) != eslOK) esl_fatal(msg);
   
      if (! diagfp)
	{
	  if ( p7_spascmx_Validate(asf, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf);
	  if ( p7_spascmx_Validate(asb, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf); 
	  if ( p7_spascmx_Validate(asd, anch, D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf); 

	  if (esl_FCompareAbs( fsc,   asc_f, tol) != eslOK) esl_fatal(msg);
	  if (esl_FCompareAbs( bsc,   asc_b, tol) != eslOK) esl_fatal(msg);
	  if (esl_FCompareAbs( asc_f, asc_b, tol) != eslOK) esl_fatal(msg);
	  if (gsc > vsc+tol)                                esl_fatal(msg);
	  if (vsc > fsc)                                    esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%20g %20g %20g %20g %20g\n", fsc-asc_f, bsc-asc_b, asc_f-asc_b, vsc-gsc, fsc-vsc);
  
      p7_sparsemx_Reuse(asf);
      p7_sparsemx_Reuse(asb);
      p7_sparsemx_Reuse(asd);
      p7_sparsemx_Reuse(sxf);
      p7_sparsemx_Reuse(sxb);

      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
      free(dsq);
      p7_trace_Destroy(gtr);
      free(anch);
      p7_sparsemask_Destroy(sm);
    }

  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(asb);
  p7_sparsemx_Destroy(asd);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_bg_Destroy(bg);
}

/* "emulated_viterbi" utest
 * 
 * If we choose an anchor set implied by the unconstrained sparse
 * Viterbi path, then sparse ASC Viterbi path is identical to that
 * path.
 * 
 * We can use this fact to test Forward and Backward, even though we
 * don't have (or need) a sparse ASC Viterbi implementation.  Any F/B
 * implementation can be converted to Viterbi scoring by making the
 * p7_FLogsum() function do a max instead of a log-sum. Our logsum
 * implementation includes the necessary machinery for a temporary
 * switch in its logic, by using p7_logsum_InitMax() and
 * p7_logsum_Reinit().
 *
 * Holds true for any model in any configuration, and any sequence,
 * and any sparse mask.
 */
static void
utest_emulated_viterbi(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int L, int N)
{
  char           msg[] = "sparse_asc_fwdback :: emulated_viterbi unit test failed";
  P7_BG         *bg    = p7_bg_Create(abc);
  P7_PRIOR      *pri   = NULL;
  P7_HMM        *hmm;
  P7_PROFILE    *gm    = p7_profile_Create(M, abc);
  ESL_SQ        *sq    = esl_sq_CreateDigital(abc);
  int            idx;
  P7_OPROFILE   *om;
  P7_CHECKPTMX  *cx;
  P7_SPARSEMASK *sm    = p7_sparsemask_Create(M,L);
  P7_SPARSEMX   *sxf   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxd   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asf   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *asb   = p7_sparsemx_Create(NULL);
  P7_TRACE      *vtr   = p7_trace_Create();
  P7_ANCHORS    *anch  = p7_anchors_Create();
  float          vsc, fsc, asc, asc_b;
  char           errbuf[eslERRBUFSIZE];
  float          tol   = 0.0001;

  if      (abc->type == eslAMINO) pri = p7_prior_CreateAmino();
  else if (abc->type == eslDNA)   pri = p7_prior_CreateNucleic();
  else if (abc->type == eslRNA)   pri = p7_prior_CreateNucleic();
  else                            pri = p7_prior_CreateLaplace(abc);

  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_Prior(rng, M, abc, pri, &hmm) != eslOK) esl_fatal(msg);
      if ( p7_profile_Config(gm, hmm, bg)               != eslOK) esl_fatal(msg);

      if ( p7_profile_SetLength(gm, L)  != eslOK) esl_fatal(msg);   /* config to generate mean length of L (length was probably reset by last emitted seq) */
      do {
	esl_sq_Reuse(sq);
	p7_ProfileEmit(rng, hmm, gm, bg, sq, NULL);
      } while (sq->n > L * 10); /* allow many domains, but keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
      if ( p7_profile_SetLength(gm, sq->n) != eslOK) esl_fatal(msg);

      if ((om = p7_oprofile_Create(hmm->M, abc))                             == NULL)  esl_fatal(msg);
      if ( p7_oprofile_Convert(gm, om)                                       != eslOK) esl_fatal(msg);
      if ((cx = p7_checkptmx_Create(hmm->M, sq->n, ESL_MBYTES(32)))          == NULL)  esl_fatal(msg);
      if ( p7_ForwardFilter (sq->dsq, sq->n, om, cx, /*fsc=*/NULL)           != eslOK) esl_fatal(msg);
      if ( p7_BackwardFilter(sq->dsq, sq->n, om, cx, sm, p7_SPARSIFY_THRESH) != eslOK) esl_fatal(msg);

      if ( p7_SparseViterbi (sq->dsq, sq->n, gm, sm, sxf, vtr, &vsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf,      NULL) != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxd,      NULL) != eslOK) esl_fatal(msg);
      if ( p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxd, sxd)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_anchors_SetFromTrace(sxd, vtr, anch)            != eslOK) esl_fatal(msg);

      p7_logsum_InitMax();    // Zero the logsum lookup table; now FLogsum() calls will compute max()
      if ( p7_SparseForward      (sq->dsq, sq->n, gm,                   sm, sxf, &fsc)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Forward (sq->dsq, sq->n, gm, anch->a, anch->D, sm, asf, &asc)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Backward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asb, &asc_b)    != eslOK) esl_fatal(msg);
      p7_logsum_Reinit();     // Reset the logsum lookup table
      
      if ( ! diagfp)
	{
	  if ( p7_spascmx_Validate(asf, anch->a, anch->D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf);
	  if ( p7_spascmx_Validate(asb, anch->a, anch->D, errbuf) != eslOK) esl_fatal("%s\n  %s\n", msg, errbuf); 

	  if (esl_FCompareAbs( vsc, fsc,   tol) != eslOK) esl_fatal(msg);
	  if (esl_FCompareAbs( fsc, asc,   tol) != eslOK) esl_fatal(msg);
	  if (esl_FCompareAbs( fsc, asc_b, tol) != eslOK) esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%4d %4d %20g %20g %20g %20g\n", (int) sq->n, anch->D, vsc, vsc-fsc, fsc-asc, fsc-asc_b);

      
      p7_anchors_Reuse(anch);
      p7_trace_Reuse(vtr);
      p7_sparsemx_Reuse(sxf);
      p7_sparsemx_Reuse(sxd);	
      p7_sparsemx_Reuse(asf);	
      p7_sparsemx_Reuse(asb);	
      p7_sparsemask_Reuse(sm);
      esl_sq_Reuse(sq);
      p7_profile_Reuse(gm);

      p7_oprofile_Destroy(om);
      p7_checkptmx_Destroy(cx);
      p7_hmm_Destroy(hmm);
    }



  p7_anchors_Destroy(anch);
  p7_trace_Destroy(vtr);
  p7_sparsemx_Destroy(sxf);  
  p7_sparsemx_Destroy(sxd);  
  p7_sparsemx_Destroy(asf);  
  p7_sparsemx_Destroy(asb);  
  p7_sparsemask_Destroy(sm);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_prior_Destroy(pri);
  p7_bg_Destroy(bg);
}

#endif /*p7SPARSE_ASC_FWDBACK_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/




/*****************************************************************
 * 7. Test driver
 *****************************************************************/
#ifdef p7SPARSE_ASC_FWDBACK_TESTDRIVE

#include "p7_config.h"

#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

#if p7_DEVELOPMENT
#define p7_RNGSEED "0"
#else
#define p7_RNGSEED "42"
#endif

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",     eslARG_NONE,        FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",     eslARG_INT,    p7_RNGSEED, NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",     eslARG_INT,          "10", NULL, NULL,  NULL,  NULL, NULL, "mean sequence length to generate",               0 },
  { "-M",     eslARG_INT,          "10", NULL, NULL,  NULL,  NULL, NULL, "length of sampled profile",                      0 },
  { "-N",     eslARG_INT,          "10", NULL, NULL,  NULL,  NULL, NULL, "number of samples to attempt per unit test",     0 },
  { "--diag", eslARG_STRING,       NULL, NULL, NULL,  NULL,  NULL, NULL, "dump data on a utest's chance failure rate",     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for sparse ASC Forward/Backward dynamic programming";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  int             M    = 10;
  int             L    = 10;
  int             N    = esl_opt_GetInteger(go, "-N");

  if (esl_opt_IsOn(go, "--diag"))
    { 
      char *which = esl_opt_GetString(go, "--diag");

      if      (strcmp(which, "generation")        == 0) utest_generation       (stdout, rng, abc, M, L, N);
      else if (strcmp(which, "compare_reference") == 0) utest_compare_reference(stdout, rng, abc, M, L, N);
      else if (strcmp(which, "singlesingle")      == 0) utest_singlesingle     (stdout, rng, abc, M,    N);
      else if (strcmp(which, "multisingle")       == 0) utest_multisingle      (stdout, rng, abc, M,    N);
      else if (strcmp(which, "multipath_local")   == 0) utest_multipath_local  (stdout, rng, abc, M,    N);
      else if (strcmp(which, "multimulti")        == 0) utest_multimulti       (stdout, rng, abc, M,    N);
      else if (strcmp(which, "viterbi")           == 0) utest_emulated_viterbi (stdout, rng, abc, M, L, N);
      else esl_fatal("--diag takes: compare_reference, singlesingle, multisingle");
    }
  else
    {
      fprintf(stderr, "## %s\n", argv[0]);
      fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

      utest_generation       (NULL, rng, abc, M, L, N);
      utest_compare_reference(NULL, rng, abc, M, L, N);
      utest_singlesingle     (NULL, rng, abc, M,    N);
      utest_multisingle      (NULL, rng, abc, M,    N);
      utest_multipath_local  (NULL, rng, abc, M,    N);
      utest_multimulti       (NULL, rng, abc, M,    N);
      utest_emulated_viterbi (NULL, rng, abc, M, L, N);

      fprintf(stderr, "#  status = ok\n");
    }

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_ASC_FWDBACK_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/


/*****************************************************************
 * 8. Example
 *****************************************************************/
#ifdef p7SPARSE_ASC_FWDBACK_EXAMPLE

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
static char banner[] = "example of sparse ASC Forward";

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
  P7_TRACE       *vtr     = NULL;
  P7_REFMX       *rxf     = NULL;
  P7_REFMX       *rxd     = NULL;
  P7_REFMX       *afu     = NULL;
  P7_REFMX       *afd     = NULL;
  P7_SPARSEMX    *asf     = NULL;
  P7_SPARSEMX    *asb     = NULL;
  P7_SPARSEMX    *asd     = NULL;
  P7_ANCHORS     *anch    = p7_anchors_Create();
  float          *wrk     = NULL;
  P7_ANCHORHASH  *ah      = p7_anchorhash_Create();
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
    p7_BackwardFilter(sq->dsq, sq->n, om, cx, sm, p7_SPARSIFY_THRESH);
  }

  /* We need an anchor set <anch>.
   * To get it, run the reference prototype code;
   * (we don't have the MPAS algorithm in its sparse production form yet)
   * TODO: now we do. Convert these to sparse versions.
   */
  vtr = p7_trace_Create();
  rxf = p7_refmx_Create(gm->M, sq->n);
  rxd = p7_refmx_Create(gm->M, sq->n);
  afu = p7_refmx_Create(gm->M, sq->n);
  afd = p7_refmx_Create(gm->M, sq->n);


  p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxf, vtr, &vsc);
  p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc);   
  p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd,      NULL);   
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd, rxd);   

  p7_reference_Anchors(rng, sq->dsq, sq->n, gm, rxf, rxd, vtr, &wrk, ah,
		       afu, afd, anch, &asc, NULL, NULL);


  /* Finally...
   * Run sparse ASC Forward.
   */
  asf = p7_sparsemx_Create(sm);
  asb = p7_sparsemx_Create(sm);
  asd = p7_sparsemx_Create(sm);

  p7_sparse_asc_Forward (sq->dsq, sq->n, gm, anch->a, anch->D, sm, asf, &asc_f);
  p7_sparse_asc_Backward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asb, &asc_b);
  p7_sparse_asc_Decoding(sq->dsq, sq->n, gm, anch->a, anch->D, asc_f, asf, asb, asd);

  //p7_spascmx_Dump(stdout, asb, anch->a, anch->D);

  p7_spascmx_Validate(asf, anch->a, anch->D, NULL);
  p7_spascmx_Validate(asb, anch->a, anch->D, NULL);
  p7_spascmx_Validate(asd, anch->a, anch->D, NULL);

  printf("Reference ASC fwd score = %.2f nats\n", asc);
  printf("Sparse ASC fwd score    = %.2f nats\n", asc_f);
  printf("Sparse ASC bck score    = %.2f nats\n", asc_b);

  p7_anchorhash_Destroy(ah);
  if (wrk) free(wrk);
  p7_sparsemx_Destroy(asb);
  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(asd);
  p7_anchors_Destroy(anch);
  p7_trace_Destroy(vtr);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(cx);
  p7_refmx_Destroy(afu);   p7_refmx_Destroy(afd);
  p7_refmx_Destroy(rxf);   p7_refmx_Destroy(rxd);
  p7_profile_Destroy(gm);  p7_oprofile_Destroy(om);
  esl_sqfile_Close(sqfp);  esl_sq_Destroy(sq);
  p7_hmmfile_Close(hfp);   p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_ASC_FWDBACK_EXAMPLE*/
/*---------------------- end, example ---------------------------*/


