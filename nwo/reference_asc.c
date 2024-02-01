/* Reference implementation of anchor-set-constrained (ASC) dynamic programming.
 *
 * The ASC technique sums over all paths consistent with a particular domain
 * annotation. A domain annotation for D domains is defined by an "anchor set" of D
 * anchors i0,k0. Each domain in a valid path must align an Mk0 state (ML or MG) to
 * residue x_i0.
 *
 * This is the reference implementation, used for testing. It is not used in HMMER's
 * main programs. The production code uses sparse ASC DP.
 *
 * These routines use the table-driven log-sum-exp approximation; caller must have
 * initialized with <h4_logsum_init()> before they will work.
 * 
 * Contents:
 *   1. ASC Forward
 *   2. ASC Backward
 *   3. ASC Decoding
 *   4. Unit tests
 *   5. Test driver
 *   6. Example
 */
#include <h4_config.h>

#include "easel.h"

#include "logsum.h"
#include "h4_profile.h"
#include "h4_mode.h"
#include "h4_anchorset.h"
#include "h4_refmx.h"



/*****************************************************************
 * 1. ASC Forward
 *****************************************************************/

/* Function:  h4_reference_asc_Forward()
 * Synopsis:  Anchor-set-constrained (ASC) Forward algorithm.
 * Incept:    SRE, Fri 01 Jan 2021
 *
 * Purpose:   The anchor-set-constrained (ASC) Forward
 *            algorithm. Compare digital sequence <dsq> of length <L>
 *            to profile <hmm> in mode <mo>, constrained by anchor set
 *            <anch>. Fill in ASC UP/DOWN matrices <mxu/mxd>, and
 *            optionally return the ASC Forward score, in bits, in
 *            <*opt_sc>.
 *
 *            DP matrices <mxu> and <mxd> can be provided as any
 *            allocated size; reused and reallocated here as needed.
 *
 * Args:      dsq     - digital sequence, 1..L
 *            L       - length of dsq
 *            hmm     - profile HMM
 *            mo      - comparison mode, with length set
 *            anch    - anchorset
 *            mxu     - RESULT: ASC UP matrix (reused/reallocated as needed)
 *            mxd     - RESULT: ASC DOWN matrix (ditto)
 *            opt_sc  - optRETURN: ASC fwd bit score
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on re-allocation failure.
 */
int
h4_reference_asc_Forward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo,
                         const H4_ANCHORSET *anch, H4_REFMX *mxu, H4_REFMX *mxd, float *opt_sc)
{
  const float *tsc;		// ptr for stepping thru profile's transition parameters 
  const float *rsc;		// ptr for stepping thru profile's emission parameters   
  float *dpp;                   // ptr into main states ({MID}{LG}) for prev row...      
  float *dpc;			//   ... and a current row.                              
  float *xp;			// ptr into specials (ENJBLGC) for a previous row...     
  float *xc;			//   ... and a current row.                              
  int    d;                     // counter over domains 1..D
  int    i;			// counter over sequence positions (rows): 0,1..L
  int    k;			// counter over model positions (columns): 0,1..M
  int    s;			// counter over states                           
  float  dlv, dgv;		//   ... and for DL, DG scores                   
  float  xE;			// tmp var for accumulating E score for a row    
  int    M = hmm->M;		// for a bit more clarity, less dereference clutter      
  int    status;

  /* Don't try to contract check the <mo> length model like we do elsewhere.
   * Envelope determination calls ASC Forward on arbitrary subsequences to
   * recalculate envelope scores.
   */

  /* reallocation, if needed */
  if ( (status = h4_refmx_GrowTo(mxu, M, L)) != eslOK) return status;
  if ( (status = h4_refmx_GrowTo(mxd, M, L)) != eslOK) return status;
  mxu->M    = mxd->M    = M;
  mxu->L    = mxd->L    = L;
  mxu->type = h4R_ASC_FWD_UP;
  mxd->type = h4R_ASC_FWD_DOWN;

#if eslDEBUGLEVEL > 0
  h4_refmx_SetValues(mxu, -eslINFINITY);
  h4_refmx_SetValues(mxd, -eslINFINITY);
#endif

  /* Initialize i=0..i0(1)-1 specials. 
   * All specials are kept in the <mxd> DOWN matrix. 
   * UP specials are unused.
   */
  xp = NULL;
  for (i = 0; i < anch->a[1].i0; i++)  // works even for D=0 case, because then a[1] is the sentinel, with a[1].i0 = L+1
    {
      xc = mxd->dp[i] + (M+1) * h4R_NSCELLS;
      xc[h4R_E]  = -eslINFINITY;
      xc[h4R_N]  = (i == 0 ? 0. : xp[h4R_N] + mo->xsc[h4_N][h4_LOOP]);
      xc[h4R_J]  = -eslINFINITY;
      xc[h4R_B]  = xc[h4R_N] + mo->xsc[h4_N][h4_MOVE];
      xc[h4R_L]  = xc[h4R_B] + mo->xsc[h4_B][0]; 
      xc[h4R_G]  = xc[h4R_B] + mo->xsc[h4_B][1]; 
      xc[h4R_C]  = -eslINFINITY;
      xc[h4R_JJ] = -eslINFINITY;
      xc[h4R_CC] = -eslINFINITY;
      xp = xc;
    }


  /* Iterate over domains d=1..D: */
  for (d = 1; d <= anch->D; d++)
    {
      /*****************************************************************
       * Part 1. UP matrix sector for domain d
       *    In an UP sector, we can enter the model, but not exit;
       *    so we make use of L,G state values from a previous row.
       *
       *    The UP sector includes:
       *       i = i0(d-1) .. i0(d)-1   (for d=1, sentinel i0(0)=0 makes first UP sector 0..i0(1)-1)
       *       k = 1 to k0(d)-1       
       *
       *    Supercells in the top UP row, i0(d-1), are only used in
       *    decoding, for delete states in G->DD->{MI}k wing unfolding.
       *    So here, in fwd, we'll initialize row i0(d-1) to all -inf,
       *    which simplifies the boundary condition for the next row.
       *    
       *    There's a tricky case when i0(d) = i0(d-1)+1, with anchors
       *    on adjacent rows. We need to connect the lower right UP supercell
       *    to the anchor supercell in DOWN; we do this by making sure <dpc>
       *    has a well-defined final state, at i0(d-1),k0. When we go to 
       *    connect to DOWN, we'll back it up one supercell, setting dpp = dpc - h4R_NSCELLS.
       *                                                                         
       *    Because the H4_REFMX has a column k=0 allocated for us, we may 
       *    as well use it; so we also initialize supercells in k=0 at
       *    every UP row to -inf, again to simplify boundary conditions.
       *    
       *    It's possible to have no cells in the UP matrix at all, if
       *    k0(d) = 1. In that case, the way we handle initialization,
       *    the combination of the above two points is important: 
       *    we initialize k=0 on the previous row, and dpc is left hovering
       *    just after it.
       *****************************************************************/

      /* Initialization of top row, i0(d-1), to -inf's */
      dpc = mxu->dp[ anch->a[d-1].i0 ];      // for d=1 that's dp[0]
      for (s = 0; s < h4R_NSCELLS * anch->a[d].k0; s++) *dpc++ = -eslINFINITY;   
      /* in the special case of i0(d)=i0(d-1)+1,
       * dpc is now sitting on i0(d)-1, k0(d)): supercell above the anchor, 
       * where it needs to be to make the connection to DOWN.
       */

      
      /* Now we recurse for remaining rows, i0(d-1)+1..i0[d].i-1 row.
       * (It's possible that no such rows exist, depending where that next anchor is.)
       */
      for (i = anch->a[d-1].i0+1 ; i < anch->a[d].i0; i++)
	{
	  rsc = hmm->rsc[dsq[i]] + 1;                               // Start <rsc> at k=1 on row i (+1 skips k=0)
	  tsc = hmm->tsc[0];                                        // Start <tsc> at k=0, where off-by-one {LG}->M transition scores are
	  xp  = mxd->dp[i-1] + (M+1) * h4R_NSCELLS;                 // <xp> set to the specials of the previous row, i-1
	  dpp = mxu->dp[i-1];                                       // Start <dpp> on k=0, which we know has been init'ed to -inf above.
	  dpc = mxu->dp[i];                                         // Start <dpc> on k=0 of row i...
	  for (s = 0; s < h4R_NSCELLS; s++) *dpc++ = -eslINFINITY;  //   ... initialize that supercell, and now <dpc> is on k=1.
	  dlv = dgv = -eslINFINITY;

	  for (k = 1; k < anch->a[d].k0; k++)
	    {
              dpc[h4R_ML] = *rsc + h4_logsum( h4_logsum(dpp[h4R_ML] + tsc[h4_MM],
                                                        dpp[h4R_IL] + tsc[h4_IM]),
                                              h4_logsum(dpp[h4R_DL] + tsc[h4_DM],
                                                        xp[h4R_L]   + tsc[h4_LM]));
              dpc[h4R_MG] = *rsc + h4_logsum( h4_logsum(dpp[h4R_MG] + tsc[h4_MM],
                                                        dpp[h4R_IG] + tsc[h4_IM]),
                                              h4_logsum(dpp[h4R_DG] + tsc[h4_DM],
                                                        xp[h4R_G]   + tsc[h4_GM]));

	      rsc++;                    // rsc advances to next match emission score, at k+1
	      tsc += h4_NTSC;           // tsc advances to transitions for states k
	      dpp += h4R_NSCELLS;	// dpp advances to states k

              dpc[h4R_IL] = h4_logsum( h4_logsum( dpp[h4R_ML] + tsc[h4_MI],
                                                  dpp[h4R_IL] + tsc[h4_II]),
                                                  dpp[h4R_DL] + tsc[h4_DI]);
              dpc[h4R_IG] = h4_logsum( h4_logsum( dpp[h4R_MG] + tsc[h4_MI],
                                                  dpp[h4R_IG] + tsc[h4_II]),
                                       h4_logsum( dpp[h4R_DG] + tsc[h4_DI],
                                                  xp[h4R_G]   + tsc[h4_GI]));  // funky wing-retracted G->D1..Dk->Ik entry path

              /* no E state update, if you're comparing to unconstrained Forward code; can't reach E from UP upstream side */

              /* Delete state, deferred storage trick */
	      dpc[h4R_DL] = dlv;
	      dpc[h4R_DG] = dgv;
              dlv = h4_logsum( h4_logsum(dpc[h4R_ML] + tsc[h4_MD],
                                         dpc[h4R_IL] + tsc[h4_ID]),
                                         dpc[h4R_DL] + tsc[h4_DD]);
              dgv = h4_logsum( h4_logsum(dpc[h4R_MG] + tsc[h4_MD],
                                         dpc[h4R_IG] + tsc[h4_ID]),
                                         dpc[h4R_DG] + tsc[h4_DD]);

              dpc += h4R_NSCELLS;    // dpc advances to what will be states k, when we roll around the for loop
	    }
	}

      /* The very last cell we calculated was the cell diagonal from
       * the anchor (i.e. where the DOWN matrix is going to start), and
                 * we want that diagonal cell for initializing DOWN.
       * <dpc> just stepped past our desired cell; step back, to set <dpp>.
       * Critically, this works even if the UP matrix was empty.
       */
      dpp = dpc - h4R_NSCELLS;

      /*****************************************************************
       * Part 2. DOWN matrix sector for domain d
       *    
       *    In a DOWN sector, we can exit the model, but we can only
       *    enter at the anchor state itself; so we collect xE on each row,
       *    and we handle the special case of anchor state entry when we
       *    initialize that supercell. After each main row is calculated,
       *    we use xE to set the specials on that same row.
       *
       *    DOWN sector d includes:
       *       i = i0(d) to i0(d+1)-1  (for d=D, sentinel i0(D+1) makes last DOWN sector i0(D)..L)
       *       k = k0(d) to M
       *       
       *    We'll initialize the top row with a partial DP calc, then
       *    do the remaining rows with the full calculation.  (With
       *    the UP matrix, we could initialize a prev row to -inf, but
       *    with the DOWN matrices, they are exactly abutting when we
       *    squeeze them down into two-matrix form.) 
       *    
       *    We do, however, take advantage of being able to initialize 
       *    the supercell at k0(d)-1 to -inf.
       *****************************************************************/

      /* Start with k0[d]-1 on first row, and set all cells to -inf*/
      i   = anch->a[d].i0;
      tsc = hmm->tsc[0] + (anch->a[d].k0-1) * h4_NTSC;         // Start <tsc> on k0-1, i.e. k-1 relative to start of calculation
      rsc = hmm->rsc[dsq[i]] + anch->a[d].k0;                  // <rsc> starts on scores for k0
      xp  = mxd->dp[i-1] + (M+1) * h4R_NSCELLS;                // <xp> on specials for i0(d)-1
      dpc = mxd->dp[i] + (anch->a[d].k0-1) * h4R_NSCELLS;      // We will initialize k0(d)-1 supercell to -inf's...
      for (s = 0; s < h4R_NSCELLS; s++) *dpc++ = -eslINFINITY; //   ... and now <dpc> now sits on k0(d)

      /* Then calculate the anchor cell (i0,k0) as an UP calc, using
       * <dpp> which was already set, above. 
       */
      dpc[h4R_ML] = *rsc + h4_logsum( h4_logsum(dpp[h4R_ML] + tsc[h4_MM],
                                                dpp[h4R_IL] + tsc[h4_IM]),
                                      h4_logsum(dpp[h4R_DL] + tsc[h4_DM],
                                                xp[h4R_L]   + tsc[h4_LM]));
      dpc[h4R_MG] = *rsc + h4_logsum( h4_logsum(dpp[h4R_MG] + tsc[h4_MM],
                                                dpp[h4R_IG] + tsc[h4_IM]),
                                      h4_logsum(dpp[h4R_DG] + tsc[h4_DM],
                                                xp[h4R_G]   + tsc[h4_GM]));

      tsc   += h4_NTSC;
      dpc[h4R_IL] = -eslINFINITY;
      dpc[h4R_IG] = -eslINFINITY;
      dpc[h4R_DL] = -eslINFINITY;
      dpc[h4R_DG] = -eslINFINITY;
      dlv = dpc[h4R_ML] + tsc[h4_MD];
      dgv = dpc[h4R_MG] + tsc[h4_MD];

      /* xE initialization counts exits from anchor cell.
       * Unlike the rest of the top row, MG/ML exits from the anchor cell
       * need to be calculated. Also, it has to watch out for the
       * glocal exit case when the anchor cell (unusually) sits on k=M.
       */
      xE  = (anch->a[d].k0 == M ? h4_logsum(dpc[h4R_ML], dpc[h4R_MG]) : dpc[h4R_ML]);
      dpc += h4R_NSCELLS;

      /* Now we initialize the rest of the top row i0(d) from k=k0(d)+1 to M,
       * which is only reachable on deletion paths from the anchor.
       */
      for (k = anch->a[d].k0+1; k <= M; k++)
	{
	  dpc[h4R_ML] = -eslINFINITY; // ML. No entry, and unreachable from other cells too. 
	  dpc[h4R_MG] = -eslINFINITY; // MG. Ditto.
	  dpc[h4R_IL] = -eslINFINITY; // IL. Not reachable on top row. 
	  dpc[h4R_IG] = -eslINFINITY; // IG. Ditto.
	  dpc[h4R_DL] = dlv;          // DL. Customary delayed store of prev calculation.
	  dpc[h4R_DG] = dgv;          // DG. Ditto.

	  tsc   += h4_NTSC;
	  
	  xE  = (k == M ?                                  // Glocal exit included if k==M.
		 h4_logsum( xE, h4_logsum( dlv, dgv)) :    // We know all non-anchor-cell M's are -inf on top row, so 
		 h4_logsum( xE, dlv));			   //   we don't include MG/ML in these sums.
	  
	  dlv    = dlv + tsc[h4_DD];
	  dgv    = dgv + tsc[h4_DD];

          dpc += h4R_NSCELLS;    // dpc advances to what will be states k, when we roll around the for loop
        }

      /* dpc now sits on the start of the specials, in mxd */
      xc = dpc;
      xc[h4R_E]  = xE;
      xc[h4R_N]  = -eslINFINITY; 
      xc[h4R_J]  = xc[h4R_E] + mo->xsc[h4_E][h4_LOOP];
      xc[h4R_B]  = xc[h4R_J] + mo->xsc[h4_J][h4_MOVE];
      xc[h4R_L]  = xc[h4R_B] + mo->xsc[h4_B][0]; 
      xc[h4R_G]  = xc[h4R_B] + mo->xsc[h4_B][1]; 
      xc[h4R_C]  = xc[h4R_E] + mo->xsc[h4_E][h4_MOVE];
      xc[h4R_JJ] = -eslINFINITY;
      xc[h4R_CC] = -eslINFINITY;

      /* Now we can do the remaining rows in the Down sector of domain d. */
      for (i = anch->a[d].i0+1 ; i < anch->a[d+1].i0; i++)         // sentinel i0(D+1) = L+1 
	{
	  rsc = hmm->rsc[dsq[i]] + anch->a[d].k0;                   // Start <rsc> on (x_i, k0) 
	  tsc = hmm->tsc[0]      + (anch->a[d].k0-1) * h4_NTSC;     // Start <tsc> on (k0-1)
	  dpp = mxd->dp[i-1]     + (anch->a[d].k0-1) * h4R_NSCELLS; // Start <dpp> on (i-1, k0-1) 
	  dpc = mxd->dp[i]       + (anch->a[d].k0-1) * h4R_NSCELLS; // Start <dpc> on (i,   k0-1)... 
	  for (s = 0; s < h4R_NSCELLS; s++) *dpc++ = -eslINFINITY;  //  ... and initialize the k0-1 cells to -inf... 
                                                           	    //  ... so, now dpc is on k0
	  dlv = dgv = xE = -eslINFINITY;

  	  for (k = anch->a[d].k0; k <= M; k++) 
	    {				  
	      dpc[h4R_ML] = *rsc + h4_logsum( h4_logsum(dpp[h4R_ML] + tsc[h4_MM],     // at first column (k0_d), only IL/IG have finite scores,
                                                        dpp[h4R_IL] + tsc[h4_IM]),    //    because we initialized k0-1 dpp's to -inf
                                                        dpp[h4R_DL] + tsc[h4_DM]);
              dpc[h4R_MG] = *rsc + h4_logsum( h4_logsum(dpp[h4R_MG] + tsc[h4_MM],
                                                        dpp[h4R_IG] + tsc[h4_IM]),
                                                        dpp[h4R_DG] + tsc[h4_DM]);
              
	      rsc++;		    // rsc advances to next match emission score (k+1)
	      tsc += h4_NTSC;       // tsc advances to transitions in states k     
	      dpp += h4R_NSCELLS;   // dpp advances to cells for states k          

              dpc[h4R_IL] = h4_logsum( h4_logsum( dpp[h4R_ML] + tsc[h4_MI],
                                                  dpp[h4R_IL] + tsc[h4_II]),
                                                  dpp[h4R_DL] + tsc[h4_DI]);
              dpc[h4R_IG] = h4_logsum( h4_logsum( dpp[h4R_MG] + tsc[h4_MI],
                                                  dpp[h4R_IG] + tsc[h4_II]),
                                                  dpp[h4R_DG] + tsc[h4_DI]);

	      xE  = (k == M ?
		     h4_logsum( xE, h4_logsum( h4_logsum(dpc[h4R_ML], dlv), h4_logsum(dpc[h4R_MG], dgv))) : // k=M includes glocal exits  
		     h4_logsum( xE, h4_logsum(dpc[h4R_ML], dlv)));                                          // k<M allows local exit only 

              dpc[h4R_DL] = dlv;                                               // DL. Customary delayed store.
	      dpc[h4R_DG] = dgv;                                               //   ... ditto for DG store.
              dlv = h4_logsum( h4_logsum(dpc[h4R_ML] + tsc[h4_MD],             // Precalculation of DL for next k.
                                         dpc[h4R_IL] + tsc[h4_ID]),
                                         dpc[h4R_DL] + tsc[h4_DD]);
              dgv = h4_logsum( h4_logsum(dpc[h4R_MG] + tsc[h4_MD],             //   ... ditto for DG calculation.
                                         dpc[h4R_IG] + tsc[h4_ID]),
                                         dpc[h4R_DG] + tsc[h4_DD]);
              dpc += h4R_NSCELLS;
	    }

          /*****************************************************************
	   *  Having finished and stored the DOWN calculation on row i, with value xE,
           *  we can calculate and store the specials - also in the DOWN matrix.
           *  dpc ptr is already on the special state storage.
	   *****************************************************************/

	  xc = dpc;
	  xp = dpp + h4R_NSCELLS;	
	  xc[h4R_E]  = xE;		
	  xc[h4R_N]  = -eslINFINITY; 
	  xc[h4R_J]  = h4_logsum( xp[h4R_J] + mo->xsc[h4_J][h4_LOOP], xc[h4R_E] + mo->xsc[h4_E][h4_LOOP] );
	  xc[h4R_B]  = xc[h4R_J] + mo->xsc[h4_J][h4_MOVE]; 
	  xc[h4R_L]  = xc[h4R_B] + mo->xsc[h4_B][0]; 
	  xc[h4R_G]  = xc[h4R_B] + mo->xsc[h4_B][1]; 
	  xc[h4R_C]  = h4_logsum( xp[h4R_C] + mo->xsc[h4_C][h4_LOOP], xc[h4R_E] + mo->xsc[h4_E][h4_MOVE] );
	  xc[h4R_JJ] = -eslINFINITY;                                                                           
	  xc[h4R_CC] = -eslINFINITY;       
	} /* end loop over rows i of DOWN sector for domain d */

       } /* end loop over domains d=0..D-1; DP calculation complete. */

  /* As we leave the DP recursion, <xc> is still sitting on the
   * special states for the last row L... even for the edge case
   * of D=0 (and the edge case L=0 which must also have D=0).
   */
  if (opt_sc) *opt_sc = xc[h4R_C] + mo->xsc[h4_C][h4_MOVE] - mo->nullsc;
  return eslOK;
}
/*-------------------- end, ASC Forward -------------------------*/



/*****************************************************************
 * 2. ASC Backward
 *****************************************************************/

/* Function:  h4_reference_asc_Backward()
 * Synopsis:  Anchor-set-constrained (ASC) Forward algorithm.
 * Incept:    SRE, Sun 07 Feb 2021
 *
 * Purpose:   The anchor-set-constrained (ASC) Backward
 *            algorithm. Compare digital sequence <dsq> of length <L>
 *            to profile <hmm> in mode <mo>, constrained to sum only
 *            over those paths that have <D> domains that use the
 *            anchor set <anch>. Return the ASC Backward score, in
 *            bits, in <*opt_sc>.
 *
 *            Caller provides two reference DP matrices <mxu> and
 *            <mxd>. They can be of any allocated size; they will be
 *            reallocated here as needed. Upon return, <mxu> and <mxd>
 *            contain the ASC Backward UP and DOWN matrices,
 *            respectively. In domain analysis, they will be needed
 *            later for posterior decoding.
 *
 * Args:      input:
 *            dsq     - digital sequence, 1..L
 *            L       - length of dsq
 *            hmm     - profile HMM
 *            mo      - comparison mode, with length set
 *            anch    - anchorset
 *
 *            output (allocated object provided):
 *            mxu     - ASC UP matrix
 *            mxd     - ASC DOWN matrix
 *
 *            output:
 *            opt_sc  - ASC bck bit score
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on re-allocation failure.
 */
int
h4_reference_asc_Backward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo,
                          const H4_ANCHORSET *anch, H4_REFMX *mxu, H4_REFMX *mxd, float *opt_sc)
{
 const float *tsc;		// ptr for stepping thru profile's transition parameters 
 const float *rsc;		// ptr into emission scores for residue dsq[i] on current row i
 const float *rsn;		// ptr into emission scores for residue dsq[i+1] on next row i+1 
 float *dpc;                    // ptr into main states ({MID}{LG}) for current row i
 float *dpn;			//   ... and next row i+1
 float *xc;			// ptr into specials (ENJBLGC) for current row
 float *xn = NULL;              //   ... and next row i+1
 int    M  = hmm->M;	        // for a bit more clarity, less dereference clutter      
 int    D  = anch->D;	        //   ... ditto
 int    d;                      // counter over domains 1..D
 int    i;			// counter over sequence positions (rows): 0,1..L
 int    k;			// counter over model positions (columns): 0,1..M
 int    s;			// counter over states                           
 int    iend;                   // same as anch->a[d+1].i0-1, for clarity: last row in current sector
 float  mgc, mlc;
 float  mgn, mln;
 float  dgn, dln;
 float  ign, iln;
 float  xE;			// tmp var for accumulating E score for a row    
 float  xG, xL;
 float  xJ, xN;
 int    status;

 /* contract checks / arg validation */
 ESL_DASSERT1( ( mo->L == L || mo->L == 0) );   // length model in comparison mode is either L (usually) or 0 (some unit tests)

 /* reallocation, if needed */
 if ( (status = h4_refmx_GrowTo(mxu, M, L)) != eslOK) return status;
 if ( (status = h4_refmx_GrowTo(mxd, M, L)) != eslOK) return status;
 mxu->M    = mxd->M    = M;
 mxu->L    = mxd->L    = L;
 mxu->type = h4R_ASC_BCK_UP;
 mxd->type = h4R_ASC_BCK_DOWN;

#if eslDEBUGLEVEL > 0
  h4_refmx_SetValues(mxu, -eslINFINITY);
  h4_refmx_SetValues(mxd, -eslINFINITY);
#endif

 /* initialize i=L..i0(D) specials, backwards from L */
 for (i = L; i >= anch->a[D].i0; i--)     // L=0,D=0 case: this would be 0..0 (thanks to i0(0)=0 sentinel)
    {
      xc = mxd->dp[i] + (M+1) * h4R_NSCELLS;
      xc[h4R_CC] = -eslINFINITY;
      xc[h4R_JJ] = -eslINFINITY;
      xc[h4R_C]  = (i == L ? mo->xsc[h4_C][h4_MOVE] : xn[h4R_C] + mo->xsc[h4_C][h4_LOOP]);
      xc[h4R_G]  = -eslINFINITY;
      xc[h4R_L]  = -eslINFINITY;
      xc[h4R_B]  = -eslINFINITY;
      xc[h4R_J]  = -eslINFINITY;
      xc[h4R_N]  = -eslINFINITY;
      xc[h4R_E]  = xc[h4R_C] + mo->xsc[h4_E][h4_MOVE];
      xn = xc;
    }

  /* The code below is designed to be easily convertible to one-row memory efficient DP, if needed */
  for (d = D; d >= 1; d--)
    {
      /* DOWN matrix (d)
       *   i = i0(d)..i0(d+1)-1
       *   k = k0(d)..M
       * calculated Backward. 
       * In the DOWN matrix, paths can end from the model, but not start in it,
       * so we evaluate {MD}->E transitions backward, but we don't evaluate 
       * B->{LG}->Mk
       */
      iend = anch->a[d+1].i0-1;                                         // <i==iend> tests used below for boundary conditions on last DOWN row
      for (i = iend; i >= anch->a[d].i0; i--)                           // at d=D, this is (L down to i0(D)), thanks to i0(D+1)=L+1 sentinel
	{
	  rsn  = (i == iend ? NULL : hmm->rsc[dsq[i+1]] + M);           // residue scores on next row; start at k=M
	  tsc  = hmm->tsc[M];                                           // transition scores: start at M
	  dpc  = mxd->dp[i] + M * h4R_NSCELLS;                          // current row of DP matrix: start at M
	  xE   = dpc[h4R_NSCELLS+h4R_E];                                // pick up the xE score; specials start at M+1, hence the h4R_NSCELLS bump here
	  dpn  = (i == iend ? NULL : mxd->dp[i+1] + M * h4R_NSCELLS);   // next row of DP matrix: start at M

	  mgn = dgn = -eslINFINITY;
	  mln = dln = -eslINFINITY;
	  ign = iln = -eslINFINITY;
	  xL  = -eslINFINITY;

          for (k = M; k >= anch->a[d].k0; k--)
	    {
	      if (i != iend) {           // in one-row memory-efficient dp, dpc could be same as dpn, so:
		ign = dpn[h4R_IG];       // pick up I scores, before storing anything in these cells
		iln = dpn[h4R_IL];       // if insert scores were non-zero, we would add rsn[I] here
	      }

	      /* M calculations. Storage deferred for one-row reasons. */
	      mgc = (k == M ? xE :                             // at k=M, MG->E is possible, and it happens to be the only transition that's possible
		     h4_logsum( h4_logsum(mgn + tsc[h4_MM],    // mgn (+ *rsn) was picked up in last k loop, so now it's i+1,k+1
                                          ign + tsc[h4_MI]),   // ign was just picked up, so it's i+1,k
                                          dgn + tsc[h4_MD]));  // dgn is remembered from prev loop, so now it's i,k+1
	      mlc =  h4_logsum( h4_logsum(mln + tsc[h4_MM],   
                                          iln + tsc[h4_MI]),
                                h4_logsum(dln + tsc[h4_MD],
					  xE));

              /* Must do IG|IL before DG|DL, to use dgn/dln values before we overwrite them */
	      dpc[h4R_IG] = h4_logsum( h4_logsum(mgn + tsc[h4_IM],
                                                 ign + tsc[h4_II]),
                                                 dgn + tsc[h4_ID]);
	      dpc[h4R_IL] = h4_logsum( h4_logsum(mln + tsc[h4_IM],
                                                 iln + tsc[h4_II]),
                                                 dln + tsc[h4_ID]);

              dpc[h4R_DG] = dgn = (k == M ?  xE :
                                  h4_logsum( h4_logsum(mgn + tsc[h4_DM],
                                                       ign + tsc[h4_DI]),
                                                       dgn + tsc[h4_DD]));
	      dpc[h4R_DL] = dln = h4_logsum( h4_logsum(mln + tsc[h4_DM],
                                                       iln + tsc[h4_DI]),
                                             h4_logsum(dln + tsc[h4_DD],
                                                       xE));
	      
	      if (i != iend) {              // pick up M[i+1][k] values, add residue emission to them;
		mgn =  dpn[h4R_MG] + *rsn;  // when we loop around, these become M[i+1][k+1] values we need for DP
		mln =  dpn[h4R_ML] + *rsn;
		rsn--;
		dpn -= h4R_NSCELLS;
	      } 

	      dpc[h4R_MG] = mgc;           // now that we've picked up mgn/mln, safe to store MG,ML
	      dpc[h4R_ML] = mlc;
	      
	      tsc -= h4_NTSC;           
	      dpc -= h4R_NSCELLS;
	    }
	} 
      /* end of the DOWN sector.
       * now mgc/mlc are the scores in the a[d].(i0,k0) anchor cell; tsc is on k0-1
       */
      rsc = hmm->rsc[ dsq[anch->a[d].i0]] + anch->a[d].k0;
      mgn = mgc + *rsc;  
      mln = mlc + *rsc;
      xG  = mgn + tsc[h4_GM];
      xL  = mln + tsc[h4_LM];
      xJ = xN = -eslINFINITY;

      /* UP sector */
      iend = anch->a[d].i0-1;                           // <i == iend> tests for boundary condition on last UP row, shorthand for <i == anch[d].i0-1>
      for (i = iend; i > anch->a[d-1].i0; i--)          // at d=1, this is i0(1)-1 down to 1, thanks to i0(0)=0 sentinel
	{
	  xc = mxd->dp[i] + (M+1) * h4R_NSCELLS;        // on specials, which are in DOWN matrix
	  xc[h4R_CC] = -eslINFINITY;                    // CC,JJ are only used in decoding matrices
	  xc[h4R_JJ] = -eslINFINITY;
	  xc[h4R_C]  = -eslINFINITY;                    // C is now unreachable, when anchor set constrained.
	  xc[h4R_G]  = xG;                              // xG was accumulated during prev row; G->Mk wing unfolded
	  xc[h4R_L]  = xL;                              // xL accumulated on prev row
	  xc[h4R_B]  = h4_logsum(xG + mo->xsc[h4_B][1],  xL + mo->xsc[h4_B][0]); 
	  xc[h4R_J]  = xJ = h4_logsum(xJ + mo->xsc[h4_J][h4_LOOP], xc[h4R_B] + mo->xsc[h4_J][h4_MOVE]);
	  xc[h4R_N]  = xN = h4_logsum(xN + mo->xsc[h4_N][h4_LOOP], xc[h4R_B] + mo->xsc[h4_N][h4_MOVE]);
	  xc[h4R_E]  = xc[h4R_J] + mo->xsc[h4_E][h4_LOOP];  

	  tsc = hmm->tsc[anch->a[d].k0-1];                                            // transition scores: start at anch[d].k0-1
	  dpc = mxu->dp[i] + (anch->a[d].k0-1) * h4R_NSCELLS;                         // on anch[d].k0-1
	  dpn = (i == iend ? NULL : mxu->dp[i+1] + (anch->a[d].k0-1) * h4R_NSCELLS);  // on anch[d].k0-1
	  rsc = hmm->rsc[dsq[i]] + (anch->a[d].k0-1);
	  rsn = (i == iend ? NULL : hmm->rsc[dsq[i+1]] + (anch->a[d].k0-1));

	  xG  = xL  = -eslINFINITY; 
	  dgn = dln = -eslINFINITY;
	  ign = iln = -eslINFINITY;
	  if (i < iend) mgn = mln = -eslINFINITY;       // not on iend, because there we allow mgn/mln to carry over from anchor cell 

	  /* The recursion is the same as for the DOWN sector, so only differences are commented on: */
	  for (k = anch->a[d].k0-1; k >= 1; k--)
	    {
	      if (i < iend) {           
		ign = dpn[h4R_IG];       
		iln = dpn[h4R_IL];       
	      }

	      /* M calculations include no E contributions: can't do M->E in UP matrix */
	      mgc =  h4_logsum( h4_logsum(mgn + tsc[h4_MM],    
                                          ign + tsc[h4_MI]),   
				          dgn + tsc[h4_MD]);
	      mlc =  h4_logsum( h4_logsum(mln + tsc[h4_MM],   
                                          iln + tsc[h4_MI]),
				          dln + tsc[h4_MD]);
		     
              /* must calc IG|IL before DG|DL to use dgn|dln before overwriting them */
	      dpc[h4R_IG] = h4_logsum( h4_logsum(mgn + tsc[h4_IM],
                                                 ign + tsc[h4_II]),
                                                 dgn + tsc[h4_ID]);
	      dpc[h4R_IL] = h4_logsum( h4_logsum(mln + tsc[h4_IM],
                                                 iln + tsc[h4_II]),
                                                 dln + tsc[h4_ID]);
	      
              /* pull back to xG, xL: UP sector specific  */
	      xG = h4_logsum(xG, mgc + *rsc + tsc[h4_GM - h4_NTSC]);
              xG = h4_logsum(xG, dpc[h4R_IG] + tsc[h4_GI]);          // this is the G->D1..Dk->Ik path in Plan9
	      xL = h4_logsum(xL, mlc + *rsc + tsc[h4_LM - h4_NTSC]);
	      rsc--;

	      /* no D->E contributions from xE in UP matrix */
	      dpc[h4R_DG] = dgn = h4_logsum( h4_logsum(mgn + tsc[h4_DM],
                                                       ign + tsc[h4_DI]),
                                                       dgn + tsc[h4_DD]);
              dpc[h4R_DL] = dln = h4_logsum( h4_logsum(mln + tsc[h4_DM],
                                                       iln + tsc[h4_DI]),
                                                       dln + tsc[h4_DD]);
	      
	      if (i < iend) {       
		mgn =  dpn[h4R_MG] + *rsn;  // when we loop around, these become M[i+1][k+1] values we need for DP
		mln =  dpn[h4R_ML] + *rsn;
		rsn--;
		dpn -= h4R_NSCELLS;
	      } else mgn = mln = -eslINFINITY;

	      dpc[h4R_MG] = mgc;           // now that we've picked up mgn/mln, safe to store MG,ML
	      dpc[h4R_ML] = mlc;
	      
	      tsc -= h4_NTSC;           
	      dpc -= h4R_NSCELLS;
	    }
	} /* end backwards loop over i for UP matrix d */
      /* i is now on anch[d-1].i, or 0 */

      dpc = mxu->dp[i];
      for (s = 0; s < h4R_NSCELLS * anch->a[d].k0; s++) *dpc++ = -eslINFINITY;   

      xc = mxd->dp[i] + (M+1) * h4R_NSCELLS;   // on specials, which are in DOWN matrix
      xc[h4R_CC] = -eslINFINITY;               // CC,JJ are only used in decoding matrices
      xc[h4R_JJ] = -eslINFINITY;
      xc[h4R_C]  = -eslINFINITY;               // C is now unreachable, when anchor set constrained.
      xc[h4R_G]  = xG;                         // xG was accumulated during prev row; G->Mk wing unfolded
      xc[h4R_L]  = xL;                         // xL accumulated on prev row
      xc[h4R_B]  = h4_logsum(xG + mo->xsc[h4_B][1],  xL + mo->xsc[h4_B][0]); 
      xc[h4R_J]  = xJ = h4_logsum(xJ + mo->xsc[h4_J][h4_LOOP], xc[h4R_B] + mo->xsc[h4_J][h4_MOVE]);
      xc[h4R_N]  = xN = h4_logsum(xN + mo->xsc[h4_N][h4_LOOP], xc[h4R_B] + mo->xsc[h4_N][h4_MOVE]);
      xc[h4R_E]  = xc[h4R_J] + mo->xsc[h4_E][h4_LOOP];  
    } /* end loop over domains d */

  if (opt_sc) *opt_sc = xN - mo->nullsc;
  return eslOK;
}


/*****************************************************************
 * 3. ASC Decoding
 *****************************************************************/

/* Function:  h4_reference_asc_Decoding()
 * Synopsis:  Anchor-set-constrained (ASC) Decoding algorithm.
 * Incept:    SRE, Fri 12 Feb 2021
 *
 * Purpose:   The anchor-set-constrained ASC Decoding algorithm.
 *            Given calculated Forward and Backward up/down matrices <afu>/<afd>
 *            and <abu>/<abd>, for a comparison of profile <hmm> to sequence <dsq> of length <L>,
 *            in comparison mode <mo>, constrained by anchor set <anch>: compute the
 *            posterior decoding matrices <apu> and <apd>.
 *
 *            Caller provides preallocated matrices <apu>,<apd>. They
 *            can be of any contents and any allocated size; they will
 *            be reinitialized and reallocated here.
 *
 *            Alternatively, the caller can overwrite the ASC Backward matrices
 *            by passing <abu>/<abd> for <apu>/apd>. That is, the call can look like
 *            <(... afu, afd, abu, abd, abu, abd)>. 
 *
 * Args:      input:
 *            dsq     - digital sequence, 1..L
 *            L       - length of dsq
 *            hmm     - profile HMM
 *            mo      - comparison mode, with length set
 *            anch    - anchorset
 *            afu     - ASC Forward UP matrix
 *            afd     - ASC Forward DOWN matrix
 *            abu     - ASC Backward UP matrix
 *            abd     - ASC Backward DOWN matrix
 *
 *            output (allocated object provided):
 *            apu     - ASC UP matrix
 *            apd     - ASC DOWN matrix
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on re-allocation failure.
 *
 * Xref:      H10/28: on <const H4_REFMX *abu, ... H4_REFMX *abu> being valid in C.
 */
int
h4_reference_asc_Decoding(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_ANCHORSET *anch,
                          const H4_REFMX *afu, const H4_REFMX *afd, const H4_REFMX *abu, const H4_REFMX *abd,
                          H4_REFMX *apu, H4_REFMX *apd)
{
  const float *rsc;		// ptr to current row's residue emission scores in <hmm>
  float       *fwdp;		// ptr into row of Forward DP matrix, UP or DOWN (<afu> or <afd>)
  float       *bckp;		// ptr into row of Backward DP matrix (<abu> or <abd>)
  float       *ppp;		// ptr into row of Decoding DP matrix (<apu> or <apd>)
  const int    M = hmm->M;	// for clarity, pull out model's size
  const int    D = anch->D;     //   ... similarly, # of domains
  int          d, i, k, s;	// indices for D domains, L residues, M nodes, and 6 states (or a whole row of them)
  int          iend;		// tmp var for start or end of a chunk in UP or DOWN matrices 
  float        totsc;		// overall Backward (=Forward) score, our normalization factor  
  float        denom;		// sum of pp for all emitting states for this residue, for renorm
  float        delta;		// piece of probability alloted to Dj in G->D1..Dk-1->Mk and G->D1..Dk->Ik wings
  float        xJ, xC;		// used to keep J,C fwd scores from prev row, for decoding JJ/CC
  float        xG;		// for clarity, tmp var, G(i-1) pulled out for wing unfolding
  int          status;

  /* Contract checks / argument validation */
  ESL_DASSERT1( (afu->type == h4R_ASC_FWD_UP)   );
  ESL_DASSERT1( (afd->type == h4R_ASC_FWD_DOWN) );
  ESL_DASSERT1( (abu->type == h4R_ASC_BCK_UP)   );
  ESL_DASSERT1( (abd->type == h4R_ASC_BCK_DOWN) );
  ESL_DASSERT1( (afu->L == L && afd->L == L && abu->L == L && abd->L == L) );
  ESL_DASSERT1( (afu->M == M && afd->M == M && abu->M == M && abd->M == M) );
 
  /* Reallocation, if needed. 
   * Caller is allowed to overwrite abu -> apu, abd -> apd
   */
  if ( apu != abu && ( status = h4_refmx_GrowTo(apu, M, L)) != eslOK) return status;
  if ( apd != abd && ( status = h4_refmx_GrowTo(apd, M, L)) != eslOK) return status;
  apu->L = apd->L = L;
  apu->M = apd->M = M;
  apu->type = h4R_ASC_DECODE_UP;
  apd->type = h4R_ASC_DECODE_DOWN;

#if eslDEBUGLEVEL > 0
  if (abu != apu && abd != apd) {           // easier to debug if non-ASC cells are -inf
    h4_refmx_SetValues(apu, -eslINFINITY);  // but we can only do this if we're not overwriting the bck mx
    h4_refmx_SetValues(apd, -eslINFINITY);
  }
#endif

  /* Initialize specials on rows 1..i0[1]-1 
   * We've above the first anchor, so only S->N->B->{LG} is possible in specials. 
   * We pick up totsc from row 0 of backwards.
   */
  for (i = 0; i < anch->a[1].i0; i++)
    {
      fwdp  = afd->dp[i] + (M+1) * h4R_NSCELLS;
      bckp  = abd->dp[i] + (M+1) * h4R_NSCELLS;
      ppp   = apd->dp[i] + (M+1) * h4R_NSCELLS;

      if (i == 0) totsc = bckp[h4R_N]; 

      ppp[h4R_JJ] = 0.0;	
      ppp[h4R_CC] = 0.0; 
      ppp[h4R_E]  = 0.0; 
      ppp[h4R_N]  = (i == 0 ? 1.0 : exp2f(fwdp[h4R_N] + bckp[h4R_N] - totsc));
      ppp[h4R_J]  = 0.0;  
      ppp[h4R_B]  = exp2f(fwdp[h4R_B] + bckp[h4R_B] - totsc);
      ppp[h4R_L]  = exp2f(fwdp[h4R_L] + bckp[h4R_L] - totsc);
      ppp[h4R_G]  = exp2f(fwdp[h4R_G] + bckp[h4R_G] - totsc); 
      ppp[h4R_C]  = 0.0;

      *(apu->dp[i]) = ppp[h4R_N];  // that's a hack. We'll need a rowwise sum over emitters for renormalizaion.
                                   // stash partial sum in k=0,ML slot of UP, which is otherwise unused
    }


  for (d = 1; d <= D; d++)
    {
      /* UP sector */
      /* Row iend-1 (0, for example, but also the top edge of each UP matrix) is a boundary case:
       * only reachable by wing retraction. Initialize it all to zero.
       */
      iend = anch->a[d-1].i0+1;
      for (s = 0; s < anch->a[d].k0*h4R_NSCELLS; s++) apu->dp[iend-1][s] = 0.;

      for (i = iend; i < anch->a[d].i0; i++)
	{
	  /* Wing unfolding of G->D1..Dk-1->MGk and G->D1..Dk->IGk paths:
           * increment PREVIOUS row i-1 with these.
	   * 
	   * In Forward/Backward we used G->{M,I}Gk entry directly, but in
           * decoding we need the posterior probs of the D's. Each Dj in the
           * PREVIOUS row gets an added correction <delta>, which is
           * the sum of all G->Mk|Ik paths that run through it, j<k.
           * This step is the only reason we need <rsc>, hence <dsq>,
           * in posterior decoding.
	   */
	  bckp  = abu->dp[i]       + (anch->a[d].k0-1) * h4R_NSCELLS;  // bckp on i, k0-1
	  ppp   = apu->dp[i-1]     + (anch->a[d].k0-1) * h4R_NSCELLS;  // ppp starts on i+1, k0-1 (PREVIOUS row)
	  rsc   = hmm->rsc[dsq[i]] + (anch->a[d].k0-1);
    	  xG    = *(afd->dp[i-1] + (M+1) * h4R_NSCELLS + h4R_G);       // ugh. I don't see any good way of avoiding this reach out into memory
	  delta = 0.0;
	  for (k = anch->a[d].k0-1; k >= 1; k--)
	    {
              delta       += exp2f(xG + hmm->tsc[k][h4_GI] + bckp[h4R_IG] - totsc);          // add G->D1..Dk->Ik
	      ppp[h4R_DG] += delta;                                                         // store at k
	      delta       += exp2f(xG + hmm->tsc[k-1][h4_GM] + *rsc + bckp[h4R_MG] - totsc); // add G->D1..Dk-1->Mk, delay store
	      ppp  -= h4R_NSCELLS;
	      bckp -= h4R_NSCELLS;
              rsc--;
	    }

          fwdp  = afu->dp[i] + h4R_NSCELLS;  // e.g. [i][1]
	  bckp  = abu->dp[i] + h4R_NSCELLS;  //  ... ditto
	  ppp   = apu->dp[i];
	  denom = *ppp;	 // pick up what we stashed earlier; we're now going to finish row i and renormalize if needed 
	  for (s = 0; s < h4R_NSCELLS; s++) *ppp++ = 0.0;

	  /* Main decoding recursion UP cells in row i
	   */
	  for (k = 1; k < anch->a[d].k0; k++)
	    {
	      ppp[h4R_ML] = exp2f(fwdp[h4R_ML] + bckp[h4R_ML] - totsc); denom += ppp[h4R_ML];
	      ppp[h4R_MG] = exp2f(fwdp[h4R_MG] + bckp[h4R_MG] - totsc); denom += ppp[h4R_MG];
	      ppp[h4R_IL] = exp2f(fwdp[h4R_IL] + bckp[h4R_IL] - totsc); denom += ppp[h4R_IL];
	      ppp[h4R_IG] = exp2f(fwdp[h4R_IG] + bckp[h4R_IG] - totsc); denom += ppp[h4R_IG];
	      ppp[h4R_DL] = exp2f(fwdp[h4R_DL] + bckp[h4R_DL] - totsc);
	      ppp[h4R_DG] = exp2f(fwdp[h4R_DG] + bckp[h4R_DG] - totsc);

	      fwdp += h4R_NSCELLS;
	      bckp += h4R_NSCELLS;
	      ppp  += h4R_NSCELLS;
	    }

          /* Because of numerical error, we force rowwise renormalization.
           * The tricky bit here is that i could be emitted either in UP or DOWN sector.
           * So to collect the sum for the denominator, we need the partial sum from DOWN (and its JJ/CC specials),
           * which is why we have that denom-stashing hack;
           * and when we renormalize row i here, we do UP, DOWN, and specials (in DOWN).
           */
#ifndef h4REFERENCE_ASC_TESTDRIVE  // don't renormalize when unit testing, in case it hides real problems
	  denom = 1.0 / denom;		                                                // multiplication may be faster than division
	  ppp   = apu->dp[i] + h4R_NSCELLS;                                             // that's k=1 in UP
	  for (s = 0; s < (anch->a[d].k0-1) * h4R_NSCELLS; s++) *ppp++ *= denom;        // UP matrix row i renormalized
	  ppp   = apd->dp[i] + (anch->a[d].k0) * h4R_NSCELLS;                           // that's k0 in DOWN
	  for (s = 0; s < (M - anch->a[d].k0 + 1) * h4R_NSCELLS; s++) *ppp++ *= denom;  // DOWN matrix row i renormalized
	  for (s = 0; s < h4R_NXCELLS; s++) ppp[s] *= denom;
#endif
        } // end loop over i's in UP sector for domain d

      /* One last wing-retracted row:
       * On last row of each UP chunk, the i0(d)-1 row, those DG's
       * are reachable on a G(i)->D1...Dk0-1->MG(i0,k0) path 
       * to the anchor cell itself that we need to unfold;
       * and the bck value at the anchor is in the DOWN matrix.
       * The G->D1..Dk->Ik path isn't relevant here because the anchor is ML|MG.
       */
      i     = anch->a[d].i0;	
      k     = anch->a[d].k0;
      xG    = *(afd->dp[i-1] + (M+1) * h4R_NSCELLS + h4R_G); // ugly, sorry. Reach out and get the xG value on row i0-1.
      bckp  = abd->dp[i]      + k * h4R_NSCELLS;             // *bckp is the anchor cell i0,k0 (in DOWN)
      rsc   = hmm->rsc[dsq[i]] + k;
      delta = exp2f(xG + hmm->tsc[k-1][h4_GM] + *rsc + bckp[h4R_MG] - totsc);  // k-1 because tGM is stored off-by-one
      ppp   = apu->dp[i-1] + h4R_NSCELLS;                    // ppp starts on i-1,k=1 cells
      for (k = 1; k < anch->a[d].k0; k++) 
	{
	  ppp[h4R_DG] += delta;   // unlike the main wing retraction above, there's only one path going through these 
	  ppp += h4R_NSCELLS;     //   DGk's (the G->Mk0 entry into the anchor cell), so there's no retro-accumulation of delta
	}

      /* main recursion for the DOWN sector */
      xJ = xC = -eslINFINITY;
      for (i = anch->a[d].i0; i < anch->a[d+1].i0; i++)
	{
	  fwdp  = afd->dp[i] + anch->a[d].k0 * h4R_NSCELLS;
	  bckp  = abd->dp[i] + anch->a[d].k0 * h4R_NSCELLS;
	  ppp   = apd->dp[i] + anch->a[d].k0 * h4R_NSCELLS;
	  denom = 0.0;
	  
	  for (k = anch->a[d].k0; k <= M; k++)
	    {
	      ppp[h4R_ML] = exp2f(fwdp[h4R_ML] + bckp[h4R_ML] - totsc); denom += ppp[h4R_ML];
	      ppp[h4R_MG] = exp2f(fwdp[h4R_MG] + bckp[h4R_MG] - totsc); denom += ppp[h4R_MG];
	      ppp[h4R_IL] = exp2f(fwdp[h4R_IL] + bckp[h4R_IL] - totsc); denom += ppp[h4R_IL];
	      ppp[h4R_IG] = exp2f(fwdp[h4R_IG] + bckp[h4R_IG] - totsc); denom += ppp[h4R_IG];
	      ppp[h4R_DL] = exp2f(fwdp[h4R_DL] + bckp[h4R_DL] - totsc);
	      ppp[h4R_DG] = exp2f(fwdp[h4R_DG] + bckp[h4R_DG] - totsc);

	      fwdp += h4R_NSCELLS;
	      bckp += h4R_NSCELLS;
	      ppp  += h4R_NSCELLS;
	    }
	  /* fwdp, bckp, ppp now all sit at M+1, start of specials */

  	  ppp[h4R_JJ] = (d == D ? 0.0 : exp2f(xJ + mo->xsc[h4_J][h4_LOOP] + bckp[h4R_J] - totsc)); xJ = fwdp[h4R_J]; denom += ppp[h4R_JJ];
	  ppp[h4R_CC] = (d  < D ? 0.0 : exp2f(xC + mo->xsc[h4_C][h4_LOOP] + bckp[h4R_C] - totsc)); xC = fwdp[h4R_C]; denom += ppp[h4R_CC];
	  ppp[h4R_E]  = exp2f(fwdp[h4R_E] + bckp[h4R_E] - totsc); 
	  ppp[h4R_N]  = 0.0;
	  ppp[h4R_J]  = (d == D ? 0.0 : exp2f(fwdp[h4R_J] + bckp[h4R_J] - totsc));
	  ppp[h4R_B]  = exp2f(fwdp[h4R_B] + bckp[h4R_B] - totsc); 
	  ppp[h4R_L]  = exp2f(fwdp[h4R_L] + bckp[h4R_L] - totsc);
	  ppp[h4R_G]  = exp2f(fwdp[h4R_G] + bckp[h4R_G] - totsc);  
	  ppp[h4R_C]  = (d  < D ? 0.0 : exp2f(fwdp[h4R_C] + bckp[h4R_C] - totsc));

          /* Even with forced renormalization, I don't think you can
	   * do much better than the error inherent in the default
	   * h4_logsum() calculation. For example, the calculation of
	   * ppp[C] above can give a pp of > 1.0, even if the sum of
	   * all emitters on the row is denom=1.0, because of logsum
	   * table approximation error in fwdp[C].
	   */
	  if (d < D) *(apu->dp[i]) = denom;  // hack: stash denom, which we'll pick up again when we do the next UP matrix.
	  //#ifndef p7REFERENCE_ASC_DECODING_TESTDRIVE
	  if (d == D) 
	    { // UP matrices only go through anch[D].i-1, so for last DOWN chunk, we need to renormalize.
	      denom = 1.0 / denom;
	      ppp = apd->dp[i] + (anch->a[d].k0) * h4R_NSCELLS;                           // that's k0 in DOWN
	      for (s = 0; s < (M - anch->a[d].k0 + 1) * h4R_NSCELLS; s++) *ppp++ *= denom; 
	      for (s = 0; s < h4R_NXCELLS; s++) ppp[s] *= denom;
	    }
	  //#endif
	} /* end loop over i's in DOWN chunk for domain d*/
    } /* end loop over domains d=0..D-1 */
  return eslOK;
}


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef h4REFERENCE_ASC_TESTDRIVE

#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "ascmx.h"
#include "emit.h"
#include "general.h"
#include "modelsample.h"
#include "reference_dp.h"

/* The unit tests create a variety of contrived comparisons where we can compare
 * standard DP matrices and results to ASC counterparts. Routines for creating these
 * contrived comparisons are in modelsample.c.
 * 
 * Most of the tests (all but "singlemulti") create comparisons where the standard
 * and the ASC calculations have identical path ensembles, by forcing the standard DP
 * calculation to use the anchor cells even without constraint by the ASC
 * algorithm. Forcing this requires unusual models, but there are various ways to do
 * it, and different ways allows us to test different aspects of the model, as
 * summarized below (also see SRE:J13/59):
 *                                                                     same   
 *                      multiple   N/C/J     I             multiple   std,ASC
 *  Unit test            paths?    emits?  emits?  Local?  domains?   ensemble?
 *  ---------------     --------  -------  ------  ------  --------   ---------
 *  singlepath              -        -      1/0       -        -        YES
 *  singlesingle            -       YES     YES       -        -        YES
 *  singlemulti             -       YES     YES       -       YES         -
 *  ensemble uni          YES       YES     YES       -        -        YES
 *  ensemble local        YES        -       -       YES       -          -
 *  ensemble multi        YES        -       -        -       YES         -  
 * 
 * The three single pathed tests verify that the decoding matrix has pp=1.0 along the
 * path, 0.0 elsewhere. (The singlemulti test is special: the ASC Decoding matrix
 * has only one nonzero path, but standard DP has an ensemble.) Ensemble
 * tests compare standard to ASC matrices.
 * 
 * When comparing ASC matrices to standard Decoding matrices cell-by-cell, if there
 * is an ensemble, we have to marginalize UP+DOWN before doing the comparison.
 */


/* "singlepath" test
 *
 * Create a profile that has only one possible path with $P(\pi |
 * \theta) = 1$, when using a uniglocal L=0 alignment mode.
 *
 * The idea is to set one transition prob to 1.0 for each state.
 * t_II have to be 0.0, so insert lengths are either 0 or 1.
 *
 * Now:
 *   1. Viterbi = Fwd = Bck = ASC Fwd = ASC Bck scores.
 *   2. Viterbi path = emitted path
 *   3. Viterbi DP matrix cells = Forward cells
 *   4. ASC Forward U/D matrix cells = Fwd cells, within ASC regions
 *   5. (ditto for ASC Backward, Bck)
 *   6. ASC Decoding matrix has 1.0 for all cells in path
 *   7. ASC Decoding matrix = std Decoding, all cells in ASC regions.
 */
static void
utest_singlepath(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M)
{
  char          failmsg[] = "reference_asc singlepath unit test failed";
  H4_PROFILE   *hmm       = NULL;
  H4_MODE      *mo        = h4_mode_Create();
  ESL_SQ       *sq        = esl_sq_CreateDigital(abc);
  H4_PATH      *pi        = h4_path_Create();           // emitted path
  H4_PATH      *vpi       = h4_path_Create();           // Viterbi path
  H4_ANCHORSET *anch      = NULL;
  H4_REFMX     *rxv       = h4_refmx_Create(100,100);   // reference Viterbi DP matrix
  H4_REFMX     *rxf       = h4_refmx_Create(100,100);   // reference Forward DP matrix
  H4_REFMX     *rxb       = h4_refmx_Create(100,100);   // reference Backward DP matrix
  H4_REFMX     *rxd       = h4_refmx_Create(100,100);   // reference Decoding DP matrix
  H4_REFMX     *afu       = h4_refmx_Create(100,100);   // ASC Forward UP matrix
  H4_REFMX     *afd       = h4_refmx_Create(100,100);   // ASC Forward DOWN matrix
  H4_REFMX     *abu       = h4_refmx_Create(100,100);   // ASC Backward UP matrix
  H4_REFMX     *abd       = h4_refmx_Create(100,100);   // ASC Backward DOWN matrix
  H4_REFMX     *apu       = h4_refmx_Create(100,100);   // ASC Decoding UP matrix
  H4_REFMX     *apd       = h4_refmx_Create(100,100);   // ASC Decoding DOWN matrix
  float         tsc, vsc, fsc, bsc, asc_f, asc_b;
  float         abstol    = 0.0001;

  if ( h4_modelsample_SinglePath(rng, abc, M, &hmm)   != eslOK) esl_fatal(failmsg);
  if ( h4_mode_SetUniglocal(mo)                       != eslOK) esl_fatal(failmsg);
  if ( h4_mode_SetLength(mo, 0)                       != eslOK) esl_fatal(failmsg);
  
  if ( h4_emit(rng, hmm, mo, sq, pi)                  != eslOK) esl_fatal(failmsg);
  if ( h4_path_Score(pi, sq->dsq, hmm, mo, &tsc)      != eslOK) esl_fatal(failmsg);

  if (( anch = h4_anchorset_Create(0, sq->n, hmm->M)) == NULL)  esl_fatal(failmsg);
  if ( h4_anchorset_SampleFromPath(rng, pi, anch)        != eslOK) esl_fatal(failmsg);

  if ( h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxv, vpi, &vsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxb,      &bsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxb, rxd)  != eslOK) esl_fatal(failmsg);

  if ( h4_reference_asc_Forward (sq->dsq, sq->n, hmm, mo, anch, afu, afd, &asc_f) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch, abu, abd, &asc_b) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, abu, abd, apu, apd) != eslOK) esl_fatal(failmsg);

#if 0
  esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, FALSE);
  h4_profile_Dump(stdout, hmm);
  h4_mode_Dump(stdout, mo);
  //h4_refmx_Dump(stdout, rxf);
  //h4_refmx_Dump(stdout, rxb);
  h4_refmx_Dump(stdout, rxd);
  //h4_refmx_Dump(stdout, afu);
  //h4_refmx_Dump(stdout, afd);
  //h4_refmx_Dump(stdout, abu);
  //h4_refmx_Dump(stdout, abd);
  h4_refmx_Dump(stdout, apu);
  h4_refmx_Dump(stdout, apd);

  h4_path_Dump(stdout, pi);
  //h4_path_Dump(stdout, vpi);
  h4_anchorset_Dump(stdout, anch);

  printf("tsc   = %.3f\n", tsc);
  printf("fsc   = %.3f\n", fsc);
  printf("bsc   = %.3f\n", bsc);
  printf("asc_f = %.3f\n", asc_f);
  printf("asc_b = %.3f\n", asc_b);
#endif

  if ( h4_path_Compare(pi, vpi)                 != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, vsc,   0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, fsc,   0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, bsc,   0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, asc_f, 0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, asc_f, 0.0, abstol) != eslOK) esl_fatal(failmsg);

  /* to compare Viterbi to Fwd matrix, we have to hack around a safety check in the structures */
  rxv->type = h4R_FORWARD;
  if (h4_refmx_Compare(rxf, rxv, abstol) != eslOK) esl_fatal(failmsg);
  rxv->type = h4R_VITERBI;

  if ( h4_ascmx_fb_compare_std(rxf, afu, afd, anch, abstol) != eslOK) esl_fatal(failmsg);
  if ( h4_ascmx_fb_compare_std(rxb, abu, abd, anch, abstol) != eslOK) esl_fatal(failmsg);

  if ( h4_ascmx_pp_compare_path(pi,  apu, apd, anch, abstol)  != eslOK) esl_fatal(failmsg);
  if ( h4_ascmx_pp_compare_std (rxd, apu, apd, anch, abstol)  != eslOK) esl_fatal(failmsg);
  if ( h4_ascmx_pp_Validate    (apu, apd, anch, abstol, NULL) != eslOK) esl_fatal(failmsg);

  h4_refmx_Destroy(afu); h4_refmx_Destroy(afd); 
  h4_refmx_Destroy(abu); h4_refmx_Destroy(abd);
  h4_refmx_Destroy(apu); h4_refmx_Destroy(apd);
  h4_refmx_Destroy(rxf); h4_refmx_Destroy(rxb);
  h4_refmx_Destroy(rxd); h4_refmx_Destroy(rxv);
  h4_path_Destroy(pi);   h4_path_Destroy(vpi);
  h4_anchorset_Destroy(anch);
  esl_sq_Destroy(sq); 
  h4_profile_Destroy(hmm); h4_mode_Destroy(mo);
}


/* "singlesingle" test
 *
 * Create a profile/sequence pair that has only one possible path with $P(\pi | x,
 * \theta) = 1$, when using uniglocal mode.
 *
 * Unlike the "singlepath" test, now we can use L>0 (N,C states) and tII > 0 (inserts
 * of length > 1).
 *
 * See modelsample.c::single_path_seq_engine() for explanation of how the case is
 * constructed.
 * 
 * Now:
 *   1. path score = Viterbi = Fwd = Bck = ASC Fwd = ASC Bck scores.
 *   2. Viterbi path = emitted path
 *   3. ASC Decoding matrix has 1.0 for all cells in path
 *   4. ASC Decoding matrix = standard Decoding mx, cell by cell
 */
static void
utest_singlesingle(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M)
{
  char          failmsg[] = "reference_asc singlesingle unit test failed";
  H4_PROFILE   *hmm       = NULL;
  H4_MODE      *mo        = NULL;
  ESL_SQ       *sq        = NULL;
  H4_PATH      *pi        = NULL;
  H4_PATH      *vpi       = h4_path_Create();          
  H4_ANCHORSET *anch      = NULL;
  H4_REFMX     *rxv       = h4_refmx_Create(100,100);  
  H4_REFMX     *rxf       = h4_refmx_Create(100,100);  
  H4_REFMX     *rxb       = h4_refmx_Create(100,100);  
  H4_REFMX     *rxd       = h4_refmx_Create(100,100);  
  H4_REFMX     *afu       = h4_refmx_Create(100,100);  
  H4_REFMX     *afd       = h4_refmx_Create(100,100);  
  H4_REFMX     *abu       = h4_refmx_Create(100,100);   
  H4_REFMX     *abd       = h4_refmx_Create(100,100);   
  H4_REFMX     *apu       = h4_refmx_Create(100,100);   
  H4_REFMX     *apd       = h4_refmx_Create(100,100);   
  float         tsc, vsc, fsc, bsc, asc_f, asc_b;
  float         abstol    = 0.0001;

  if ( h4_modelsample_SinglePathSeq(rng, abc, M, &hmm, &sq, &mo, &pi, &anch, &tsc)  != eslOK) esl_fatal(failmsg);

  if ( h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxv, vpi, &vsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxb,      &bsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxb, rxd)  != eslOK) esl_fatal(failmsg);

  if ( h4_reference_asc_Forward (sq->dsq, sq->n, hmm, mo, anch, afu, afd, &asc_f) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch, abu, abd, &asc_b) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, abu, abd, apu, apd) != eslOK) esl_fatal(failmsg);

#if 0
  esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, FALSE);
  h4_profile_Dump(stdout, hmm);
  h4_mode_Dump(stdout, mo);
 
  h4_refmx_Dump(stdout, rxf);
  //  h4_refmx_Dump(stdout, rxb);
  //  h4_refmx_Dump(stdout, rxd);
  h4_refmx_Dump(stdout, afu);
  h4_refmx_Dump(stdout, afd);
  //  h4_refmx_Dump(stdout, abu);
  //  h4_refmx_Dump(stdout, abd);
  //  h4_refmx_Dump(stdout, apu);
  //  h4_refmx_Dump(stdout, apd);

  h4_path_Dump(stdout, pi);
  //h4_path_Dump(stdout, vpi);
  h4_anchorset_Dump(stdout, anch);

  printf("tsc   = %.3f\n", tsc);
  printf("fsc   = %.3f\n", fsc);
  printf("bsc   = %.3f\n", bsc);
  printf("asc_f = %.3f\n", asc_f);
  printf("asc_b = %.3f\n", asc_b);
#endif

  if ( h4_path_Compare(pi, vpi)                 != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, vsc,   0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, fsc,   0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, bsc,   0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, asc_f, 0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, asc_b, 0.0, abstol) != eslOK) esl_fatal(failmsg);

  if ( h4_ascmx_pp_Validate(apu, apd, anch, abstol, NULL)    != eslOK) esl_fatal(failmsg);
  if ( h4_ascmx_pp_compare_path(pi,  apu, apd, anch, abstol) != eslOK) esl_fatal(failmsg);
  if ( h4_ascmx_pp_compare_std (rxd, apu, apd, anch, abstol) != eslOK) esl_fatal(failmsg);

  h4_refmx_Destroy(afu);   h4_refmx_Destroy(afd);
  h4_refmx_Destroy(abu);   h4_refmx_Destroy(abd);
  h4_refmx_Destroy(apu);   h4_refmx_Destroy(apd);
  h4_refmx_Destroy(rxf);   h4_refmx_Destroy(rxb);
  h4_refmx_Destroy(rxd);   h4_refmx_Destroy(rxv);
  h4_path_Destroy(pi);     h4_path_Destroy(vpi);
  h4_profile_Destroy(hmm); h4_mode_Destroy(mo);
  h4_anchorset_Destroy(anch);
  esl_sq_Destroy(sq);
}


/* "singlemulti" test
 *
 * Create a profile/sequence/anchorset triplet that has only one possible
 * path with $P(\pi | x, A, \theta) = 1$, when using a glocal mode.
 *
 * Unlike the "singlesingle" test, now we can use multihit mode.
 *
 * With standard (not ASC) algorithms, more than one path may be
 * possible.  Therefore reference Forward, Backward, and even Viterbi
 * scores can exceed ASC F/B scores and the path score.
 *
 * See <h4_modelsample_SinglePathedASC()> and its
 * <single_path_seq_engine()> for explanation of how the case is
 * constructed.
 * 
 * Now:
 *   1. path score = ASC Fwd = ASC Bck score.
 *   2. ASC Decoding matrix has 1.0 for all cells in path
 */
static void
utest_singlemulti(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M)
{
  char          failmsg[] = "reference_asc singlemulti unit test failed";
  H4_PROFILE   *hmm       = NULL;
  H4_MODE      *mo        = NULL;
  ESL_SQ       *sq        = NULL;
  H4_PATH      *pi        = NULL;
  H4_PATH      *vpi       = h4_path_Create();          
  H4_ANCHORSET *anch      = NULL;
  H4_REFMX     *rxv       = h4_refmx_Create(100,100);  
  H4_REFMX     *rxf       = h4_refmx_Create(100,100);  
  H4_REFMX     *rxb       = h4_refmx_Create(100,100);  
  H4_REFMX     *rxd       = h4_refmx_Create(100,100);  
  H4_REFMX     *afu       = h4_refmx_Create(100,100);  
  H4_REFMX     *afd       = h4_refmx_Create(100,100);  
  H4_REFMX     *abu       = h4_refmx_Create(100,100);  
  H4_REFMX     *abd       = h4_refmx_Create(100,100);  
  H4_REFMX     *apu       = h4_refmx_Create(100,100);  
  H4_REFMX     *apd       = h4_refmx_Create(100,100);  
  float         tsc, vsc, fsc, bsc, asc_f, asc_b;
  float         abstol    = 0.001;

  if ( h4_modelsample_SinglePathASC(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &tsc)  != eslOK) esl_fatal(failmsg);

  if ( h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxv, vpi, &vsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxb,      &bsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxb, rxd)  != eslOK) esl_fatal(failmsg);

  if ( h4_reference_asc_Forward (sq->dsq, sq->n, hmm, mo, anch, afu, afd, &asc_f) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch, abu, abd, &asc_b) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, abu, abd, apu, apd) != eslOK) esl_fatal(failmsg);

#if 0
  esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, FALSE);
  h4_profile_Dump(stdout, hmm);
  h4_mode_Dump(stdout, mo);

  h4_refmx_Dump(stdout, rxf);
  h4_refmx_Dump(stdout, rxb);
  h4_refmx_Dump(stdout, rxd);
  h4_refmx_Dump(stdout, afu);
  h4_refmx_Dump(stdout, afd);
  h4_refmx_Dump(stdout, abu);
  h4_refmx_Dump(stdout, abd);
  h4_refmx_Dump(stdout, apu);
  h4_refmx_Dump(stdout, apd);

  h4_path_Dump(stdout, pi);
  h4_path_Dump(stdout, vpi);
  h4_anchorset_Dump(stdout, anch);

  printf("tsc   = %.3f\n", tsc);
  printf("vsc   = %.3f\n", vsc);
  printf("fsc   = %.3f\n", fsc);
  printf("bsc   = %.3f\n", bsc);
  printf("asc_f = %.3f\n", asc_f);
  printf("asc_b = %.3f\n", asc_b);
#endif

  if ( esl_FCompare(tsc, asc_f, 0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompare(tsc, asc_b, 0.0, abstol) != eslOK) esl_fatal(failmsg);
 
  if ( h4_ascmx_pp_Validate(apu, apd, anch, abstol, NULL)    != eslOK) esl_fatal(failmsg);
  if ( h4_ascmx_pp_compare_path(pi,  apu, apd, anch, abstol) != eslOK) esl_fatal(failmsg);

  h4_refmx_Destroy(afu);   h4_refmx_Destroy(afd);
  h4_refmx_Destroy(abu);   h4_refmx_Destroy(abd);
  h4_refmx_Destroy(apu);   h4_refmx_Destroy(apd);
  h4_refmx_Destroy(rxf);   h4_refmx_Destroy(rxb);
  h4_refmx_Destroy(rxd);   h4_refmx_Destroy(rxv);
  h4_path_Destroy(pi);     h4_path_Destroy(vpi);
  h4_profile_Destroy(hmm); h4_mode_Destroy(mo);
  h4_anchorset_Destroy(anch);
  esl_sq_Destroy(sq);
}


/* "ensemble" tests
 *
 * Contrive a profile/sequence/anchorset such that the ensemble of all
 * valid paths must pass through the anchor set, so standard F/B score
 * = ASC F/B score.
 *
 * The general idea is to put an X in the sequence at each anchor
 * point i0 (where X is a randomly chosen residue). Match state
 * emissions are set so only Mk0 can generate X, and all other Mk have
 * e_k(X) = 0 and can't.
 *
 * Three different versions of the test:
 *
 * "uni" test (which=0), with h4_modelsample_AnchoredUni():
 *     Allows N/C/J and I, unlike the other tests. 
 *     Profile in unihit glocal mode. Transitions t_{k0-1} to M are
 *     set to 1.0 and e_k0(X) = 1.0, so Mk0 must be used once and it
 *     can't align to anything but the one X.
 *     
 * "local" test (which=1), with h4_modelsample_AnchoredLocal():
 *     Allows local paths, unlike the other tests.
 *     Profile in unihit local/glocal L=0 mode. Only M states can
 *     emit; N/C/J disallowed (by L=0) and I disallowed (by setting
 *     tXI transitions to 0). Only Mk0 can emit X with nonzero
 *     probability, so even local paths must align Mk0 to the
 *     one X.  (Here e_k0(X)>0 is sufficient, unlike the other two
 *     versions; and t_{k0-1} transitions don't need to force to Mk0.)
 *
 * "multi" test (which=2), with h4_modelsample_AnchoredMulti():
 *     Allows multihit, unlike the other tests.
 *     Profile in multihit glocal L=0 mode. Only M states can emit;
 *     N/C/J and I disallowed. Transitions t_{k0-1} to M are set to
 *     1.0 for force Mk0 to be used, e_k0(X) = 1.0, and e_k(X) = 0 for
 *     all other k != k0. Now Mk0 must be used in each domain, Mk0
 *     cannot align to anything but X, and no other state besides Mk0
 *     can generate X: all paths must use D Mk0,X anchors for D domains.
 *
 * Now:
 *     1. Fwd/Bck score = ASC Fwd/Bck score
 *     2. ASC Decoding matrices = std Decoding matrix, cell by cell
 * 
 * And also, obviously:
 *     2. Viterbi score >= score of sampled path that created sq
 *     3. Fwd/Bck score >= Viterbi score
 */
static void
utest_ensemble(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int which)
{
  char          failmsg0[] = "reference_asc ensemble uni unit test failed";
  char          failmsg1[] = "reference_asc ensemble local unit test failed";
  char          failmsg2[] = "reference_asc ensemble multi test failed";
  char         *failmsg;
  H4_PROFILE   *hmm        = NULL;
  H4_MODE      *mo         = NULL;
  ESL_SQ       *sq         = NULL;
  H4_PATH      *pi         = NULL;
  H4_PATH      *vpi        = h4_path_Create();          
  H4_ANCHORSET *anch       = NULL;
  H4_REFMX     *rxv        = h4_refmx_Create(100,100);  
  H4_REFMX     *rxf        = h4_refmx_Create(100,100);  
  H4_REFMX     *rxb        = h4_refmx_Create(100,100);  
  H4_REFMX     *rxd        = h4_refmx_Create(100,100);  
  H4_REFMX     *afu        = h4_refmx_Create(100,100);  
  H4_REFMX     *afd        = h4_refmx_Create(100,100);  
  H4_REFMX     *abu        = h4_refmx_Create(100,100);  
  H4_REFMX     *abd        = h4_refmx_Create(100,100);  
  H4_REFMX     *apu        = h4_refmx_Create(100,100);  
  H4_REFMX     *apd        = h4_refmx_Create(100,100);  
  float         tsc, vsc, fsc, bsc, asc_f, asc_b;
  float         fwd_atol;
  float         vit_atol;
  float         dec_atol;

  ESL_DASSERT1(( which >= 0 && which <= 2 ));

  switch (which) {
  case 0:  failmsg = failmsg0;  fwd_atol = (h4_logsum_IsSlowExact() ? 0.001  : 0.002);  vit_atol = 1e-5;  dec_atol = 0.001; break;
  case 1:  failmsg = failmsg1;  fwd_atol = (h4_logsum_IsSlowExact() ? 0.001  : 0.002);  vit_atol = 1e-5;  dec_atol = 0.001;  break;
  case 2:  failmsg = failmsg2;  fwd_atol = (h4_logsum_IsSlowExact() ? 0.0001 : 0.01);   vit_atol = 1e-5;  dec_atol = (h4_logsum_IsSlowExact() ? 0.001 : 0.01); break;
  }
  
  switch (which) {
  case 0: if ( h4_modelsample_AnchoredUni  (rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &tsc)  != eslOK) esl_fatal(failmsg); break;
  case 1: if ( h4_modelsample_AnchoredLocal(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &tsc)  != eslOK) esl_fatal(failmsg); break;
  case 2: if ( h4_modelsample_AnchoredMulti(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &tsc)  != eslOK) esl_fatal(failmsg); break;
  }

  if ( h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxv, vpi, &vsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxb,      &bsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxb, rxd)  != eslOK) esl_fatal(failmsg);

  if ( h4_reference_asc_Forward (sq->dsq, sq->n, hmm, mo, anch, afu, afd, &asc_f) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_asc_Backward(sq->dsq, sq->n, hmm, mo, anch, abu, abd, &asc_b) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_asc_Decoding(sq->dsq, sq->n, hmm, mo, anch, afu, afd, abu, abd, apu, apd) != eslOK) esl_fatal(failmsg);
  
#if 0
  esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, FALSE);
  h4_profile_Dump(stdout, hmm);
  h4_mode_Dump(stdout, mo);

  h4_refmx_Dump(stdout, rxf);
  h4_refmx_Dump(stdout, rxb);
  h4_refmx_Dump(stdout, rxd);
  h4_refmx_Dump(stdout, afu);
  h4_refmx_Dump(stdout, afd);
  h4_refmx_Dump(stdout, abu);
  h4_refmx_Dump(stdout, abd);
  h4_refmx_Dump(stdout, apu);
  h4_refmx_Dump(stdout, apd);

  h4_path_Dump(stdout, pi);
  h4_path_Dump(stdout, vpi);
  h4_anchorset_Dump(stdout, anch);

  printf("tsc   = %.3f\n", tsc);
  printf("vsc   = %.3f\n", vsc);
  printf("fsc   = %.3f\n", fsc);
  printf("bsc   = %.3f\n", bsc);
  printf("asc_f = %.3f\n", asc_f);
  printf("asc_b = %.3f\n", asc_b);
#endif

  if (esl_FCompare(fsc, bsc,   0.0, fwd_atol) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(fsc, asc_f, 0.0, fwd_atol) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(fsc, asc_b, 0.0, fwd_atol) != eslOK) esl_fatal(failmsg);
  if (tsc > vsc+vit_atol)                                  esl_fatal(failmsg);
  if (vsc > fsc)                                           esl_fatal(failmsg);

  if ( h4_ascmx_pp_Validate(apu, apd, anch, dec_atol, NULL)     != eslOK) esl_fatal(failmsg);
  if ( h4_ascmx_pp_compare_std (rxd, apu, apd, anch, dec_atol)  != eslOK) esl_fatal(failmsg);

  h4_refmx_Destroy(afu);    h4_refmx_Destroy(afd);
  h4_refmx_Destroy(abu);    h4_refmx_Destroy(abd);
  h4_refmx_Destroy(apu);    h4_refmx_Destroy(apd);
  h4_refmx_Destroy(rxf);    h4_refmx_Destroy(rxb);
  h4_refmx_Destroy(rxd);    h4_refmx_Destroy(rxv);
  h4_path_Destroy(vpi);     h4_path_Destroy(pi);
  h4_profile_Destroy(hmm);  h4_mode_Destroy(mo);
  h4_anchorset_Destroy(anch);
  esl_sq_Destroy(sq);
}
#endif // h4REFERENCE_ASC_TESTDRIVE

/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef h4REFERENCE_ASC_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                   docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show brief help summary",                 0 },
  { "-s",         eslARG_INT,     "0", NULL, NULL,  NULL,   NULL,  NULL, "set random number generator seed",        0 },
  { "-M",         eslARG_INT,    "10", NULL, NULL,  NULL,   NULL,  NULL, "set test profile length",                 0 },
  { "-N",         eslARG_INT,   "100", NULL, NULL,  NULL,"--diag", NULL, "number of times to run utest for --diag", 0 },
  { "--diag",     eslARG_STRING, NULL, NULL, NULL,  NULL,   NULL,  NULL, "dump data on utest <s> failure rate",     0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show HMMER version number",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv, "test driver for reference_asc", "[-options]");
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);
  int             M   = esl_opt_GetInteger(go, "-M");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng)); 

  utest_singlepath  (/*diagfp=*/NULL, rng, abc, M);
  utest_singlesingle(           NULL, rng, abc, M);
  utest_singlemulti (           NULL, rng, abc, M);

  utest_ensemble(NULL, rng, abc, M, 0); // 0 = ensemble uni test
  utest_ensemble(NULL, rng, abc, M, 1); // 1 = ensemble local test
  utest_ensemble(NULL, rng, abc, M, 2); // 2 = ensemble multi test

  fprintf(stderr, "#  status   = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}

#endif // h4REFERENCE_ASC_TESTDRIVE

/***************************************************************** 
 * 6. Example
 *****************************************************************/
#ifdef h4REFERENCE_ASC_EXAMPLE

#include <h4_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_refmx.h"

#include "general.h"
#include "reference_asc.h"
#include "reference_dp.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show HMMER version info",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile> <ndom> [<i0> <k0>]...";
static char banner[] = "example of ASC reference implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, -1, argc, argv, banner, usage);  // -1 means allow any number of cmdline args
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  H4_ANCHORSET   *anch    = NULL;
  H4_REFMX       *mxu     = NULL;
  H4_REFMX       *mxd     = NULL;
  H4_REFMX       *fwd     = NULL;
  int             D, d, i, i0, k0;
  float           fsc, asc_fsc;
  int             status;

  /* Read one profile */
  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);

  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format unrecognized.");
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);
 
  /* Read one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      esl_fatal("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  /* Read anchorset from command line */
  D = strtol( esl_opt_GetArg(go, 3), NULL, 10);
  anch = h4_anchorset_Create(0, sq->n, hmm->M);
  for (i = 4, d = 1; d <= D; d++)
    {
      i0 = strtol( esl_opt_GetArg(go, i), NULL, 10); i++;
      k0 = strtol( esl_opt_GetArg(go, i), NULL, 10); i++;
      h4_anchorset_Add(anch, i0, k0);
    }

  /* Allocate UP and DOWN matrices */
  mxu = h4_refmx_Create(hmm->M, sq->n);
  mxd = h4_refmx_Create(hmm->M, sq->n);
  fwd = h4_refmx_Create(hmm->M, sq->n);

  h4_mode_SetLength(mo, sq->n);

  /* Run ASC Forward */
  h4_reference_asc_Forward(sq->dsq, sq->n, hmm, mo, anch, mxu, mxd, &asc_fsc);
  printf("ASC Fwd score: %.3f bits\n", asc_fsc);

  h4_refmx_Dump(stdout, mxu);
  h4_refmx_Dump(stdout, mxd);

  /* Run standard (unconstrained) Forward for comparison */
  h4_reference_Forward(sq->dsq, sq->n, hmm, mo, fwd, &fsc);
  h4_refmx_Dump(stdout, fwd);
  printf("Fwd score: %.3f bits\n", fsc);

  /* Cleanup */
  h4_hmmfile_Close(hfp);
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  h4_refmx_Destroy(mxu);
  h4_refmx_Destroy(mxd);
  h4_refmx_Destroy(fwd);
  h4_anchorset_Destroy(anch);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  exit(0);
}
#endif // h4REFERENCE_ASC_EXAMPLE
