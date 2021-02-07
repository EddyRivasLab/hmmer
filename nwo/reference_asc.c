/* Reference implementations of anchor-set-constrained (ASC) dynamic
 * programming.
 *
 * All reference implementation code is for testing. It is not used in
 * HMMER's main executables. The production code uses sparse ASC DP.
 *
 * Contents:
 *   1. ASC Forward
 *   x. Unit tests
 *   x. Test driver
 *   x. Example
 */
#include "h4_config.h"

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
 *            to profile <hmm> in mode <mo>, constrained to sum only
 *            over those paths that have <D> domains that use the
 *            anchor set <anch>. Return the ASC Forward raw score, in
 *            bits, in <*opt_sc>.
 *
 *            Caller provides two reference DP matrices <mxu> and
 *            <mxd>. They can be of any allocated size; they will be
 *            reallocated here as needed. Upon return, <mxu> and <mxd>
 *            contain the ASC Forward UP and DOWN matrices,
 *            respectively. In domain analysis, they will be needed
 *            later for posterior decoding.
 *
 *            
 *
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
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
  if (opt_sc) *opt_sc = xc[h4R_C] + mo->xsc[h4_C][h4_MOVE]; /* C->T */
  return eslOK;
}
/*-------------------- end, ASC Forward -------------------------*/


/*****************************************************************
 * x. Unit tests
 *****************************************************************/
#ifdef h4REFERENCE_ASC_TESTDRIVE

#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "emit.h"
#include "general.h"
#include "modelsample.h"
#include "reference_dp.h"


/* ascmatrix_compare()
 * 
 * Compare a standard DP matrix <std> to ASC UP and DOWN matrices
 * <ascu> and <ascd>, using anchor set <anch>.  Compare all valid
 * values in the ASC matrices to their counterparts in the standard
 * matrix. If any absolute difference exceeds <tolerance>, return
 * <eslFAIL>. Else return <eslOK>.
 * 
 * Similar to h4_refmx_Compare(). See notes there for why we prefer
 * absolute, not relative difference (short version: error accumulates
 * more as an absolute # per term added; DP cells with values close to
 * zero may have large relative errors).
 */
static int
ascmatrix_compare(H4_REFMX *std, H4_REFMX *ascu, H4_REFMX *ascd, H4_ANCHORSET *anch, float abstol)
{
  char  msg[] = "comparison of ASC DP F|B matrix to standard DP F|B matrix failed";
  int   M     = std->M;
  int   L     = std->L;
  float val;
  int   d,i,k,s;

  /* Contract checks */
  ESL_DASSERT1( (ascu->M == M && ascu->L == L));
  ESL_DASSERT1( (ascd->M == M && ascd->L == L));
  ESL_DASSERT1( ((std->type == h4R_FORWARD  && ascu->type == h4R_ASC_FWD_UP && ascd->type == h4R_ASC_FWD_DOWN) ||
		 (std->type == h4R_BACKWARD && ascu->type == h4R_ASC_BCK_UP && ascd->type == h4R_ASC_BCK_DOWN)));

  for (d = 1, i = 0; i <= L; i++)
    {
      if (i == anch->a[d].i0) d++;  // d = idx of current UP matrix; d-1 is current DOWN
      for (k = 1; k <= M; k++)      //   ... alternatively, you can think d = idx of next anchor we'll reach (1..D; D+1 when we're past final anchor)
	for (s = 0; s < h4R_NSCELLS; s++)
	  {
            // It's possible for an i,k to be valid in both DOWN and UP matrices, in which case
            // some paths are in each; we have to sum to get what's in the unconstrained matrix.
            // Thus we traverse the afu/afd matrices in an unusual pattern, sweeping out all i,k
            // and testing which ones are valid in afu and afd.
	    val = -eslINFINITY;
	    if (d > 1        && k >= anch->a[d-1].k0) val = h4_logsum(val, H4R_MX(ascd,i,k,s)); // DOWN (valid rows start at a[1].i0. Start looking at DOWN at first anchor.)
	    if (d <= anch->D && k <  anch->a[d].k0)   val = h4_logsum(val, H4R_MX(ascu,i,k,s)); // UP   (valid rows end at a[D].i0. Stop looking at UP when we're at last anchor.)
	    
	    if (val != -eslINFINITY && esl_FCompareNew(H4R_MX(std,i,k,s), val, 0.0, abstol) != eslOK)
	      ESL_FAIL(eslFAIL, NULL, msg);
	  }

      for (s = 0; s < h4R_NXCELLS; s++)
	{
	  /* special exception: ignore L state in this test, for Backwards.
	   * reference backward can have valid prob for {MI}(k<k0) on
	   * the anchor row, and via M->I, reference can reach an Mk
	   * on any row i that ASC can't reach, in addition to Mk that
	   * ASC can reach; thus L can legitimately differ between ASC
	   * and reference Backward. [H5/26]
	   */
	  if (std->type == h4R_BACKWARD && s == h4R_L) continue;

	  val = H4R_XMX(ascd,i,s);
	  if (val != -eslINFINITY && esl_FCompareNew(H4R_XMX(std,i,s), val, 0.0, abstol) != eslOK)
	    ESL_FAIL(eslFAIL, NULL, msg);
	}
    }
  return eslOK;
}




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
 *   2. Viterbi trace = emitted trace
 *   3. Viterbi DP matrix cells = Forward cells
 *   4. ASC Forward U/D matrix cells = Fwd cells, within ASC regions
 *   5. (ditto for ASC Backward, Bck)
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
  H4_REFMX     *afu       = h4_refmx_Create(100,100);   // ASC Forward UP matrix
  H4_REFMX     *afd       = h4_refmx_Create(100,100);   // ASC Forward DOWN matrix
  float         tsc, vsc, fsc, asc_f;
  float         abstol    = 0.0001;

  if ( h4_modelsample_SinglePath(rng, abc, M, &hmm)   != eslOK) esl_fatal(failmsg);
  if ( h4_mode_SetUniglocal(mo)                       != eslOK) esl_fatal(failmsg);
  if ( h4_mode_SetLength(mo, 0)                       != eslOK) esl_fatal(failmsg);
  
  if ( h4_emit(rng, hmm, mo, sq, pi)                  != eslOK) esl_fatal(failmsg);
  if ( h4_path_Score(pi, sq->dsq, hmm, mo, &tsc)      != eslOK) esl_fatal(failmsg);

  if (( anch = h4_anchorset_Create(0, sq->n, hmm->M)) == NULL)  esl_fatal(failmsg);
  if ( h4_anchorset_SetFromPath(rng, pi, anch)        != eslOK) esl_fatal(failmsg);

  if ( h4_reference_Viterbi(sq->dsq, sq->n, hmm, mo, rxv, vpi, &vsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Forward(sq->dsq, sq->n, hmm, mo, rxf,      &fsc) != eslOK) esl_fatal(failmsg);

  if ( h4_reference_asc_Forward(sq->dsq, sq->n, hmm, mo, anch, afu, afd, &asc_f) != eslOK) esl_fatal(failmsg);

#if 0
  h4_refmx_Dump(stdout, rxf);
  h4_refmx_Dump(stdout, afu);
  h4_refmx_Dump(stdout, afd);

  h4_path_Dump(stdout, pi);
  h4_path_Dump(stdout, vpi);
  h4_anchorset_Dump(stdout, anch);

  printf("tsc   = %.3f\n", tsc);
  printf("asc_f = %.3f\n", asc_f);
  printf("diff  = %.3f\n", tsc - asc_f);
#endif

  if ( h4_path_Compare(pi, vpi)                 != eslOK) esl_fatal(failmsg);
  if ( esl_FCompareNew(tsc, vsc,   0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompareNew(tsc, fsc,   0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompareNew(tsc, asc_f, 0.0, abstol) != eslOK) esl_fatal(failmsg);

  /* to compare Viterbi to Fwd matrix, we have to hack around a safety check in the structures */
  rxv->type = h4R_FORWARD;
  if (h4_refmx_Compare(rxf, rxv, abstol) != eslOK) esl_fatal(failmsg);
  rxv->type = h4R_VITERBI;

  if ( ascmatrix_compare(rxf, afu, afd, anch, abstol) != eslOK) esl_fatal(failmsg);

  h4_refmx_Destroy(afd);
  h4_refmx_Destroy(afu);
  h4_refmx_Destroy(rxf);
  h4_refmx_Destroy(rxv);
  h4_anchorset_Destroy(anch);
  h4_path_Destroy(pi);
  h4_path_Destroy(vpi);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
}


/* "singlesingle" test
 *
 * Create a profile/sequence pair that has only one possible path
 * with $P(\pi | x, \theta) = 1$, when using uniglocal mode.
 *
 * Unlike the "singlepath" test, now we can use L>0 (N,C states)
 * and tII > 0 (inserts of length > 1).
 *
 * See h4_modelsample_SinglePathSeq() for explanation of how
 * the case is constructed.
 * 
 * Now:
 *   1. trace score = Viterbi = Fwd = Bck = ASC Fwd = ASC Bck scores.
 *   2. Viterbi trace = emitted trace
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
  H4_REFMX     *afu       = h4_refmx_Create(100,100);  
  H4_REFMX     *afd       = h4_refmx_Create(100,100);  
  float         tsc, vsc, fsc, asc_f;
  float         abstol    = 0.0001;

  if ( h4_modelsample_SinglePathSeq(rng, abc, M, &hmm, &sq, &mo, &pi, &anch, &tsc)  != eslOK) esl_fatal(failmsg);

  if ( h4_reference_Viterbi(sq->dsq, sq->n, hmm, mo, rxv, vpi, &vsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Forward(sq->dsq, sq->n, hmm, mo, rxf,      &fsc) != eslOK) esl_fatal(failmsg);

  if ( h4_reference_asc_Forward(sq->dsq, sq->n, hmm, mo, anch, afu, afd, &asc_f) != eslOK) esl_fatal(failmsg);

#if 0
  h4_refmx_Dump(stdout, rxf);
  h4_refmx_Dump(stdout, afu);
  h4_refmx_Dump(stdout, afd);

  h4_path_Dump(stdout, pi);
  h4_path_Dump(stdout, vpi);
  h4_anchorset_Dump(stdout, anch);

  printf("tsc   = %.3f\n", tsc);
  printf("asc_f = %.3f\n", asc_f);
  printf("diff  = %.3f\n", tsc - asc_f);
#endif

  if ( h4_path_Compare(pi, vpi)                 != eslOK) esl_fatal(failmsg);
  if ( esl_FCompareNew(tsc, vsc,   0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompareNew(tsc, fsc,   0.0, abstol) != eslOK) esl_fatal(failmsg);
  if ( esl_FCompareNew(tsc, asc_f, 0.0, abstol) != eslOK) esl_fatal(failmsg);

  h4_refmx_Destroy(afd);
  h4_refmx_Destroy(afu);
  h4_refmx_Destroy(rxf);
  h4_refmx_Destroy(rxv);
  h4_anchorset_Destroy(anch);
  h4_path_Destroy(pi);
  h4_path_Destroy(vpi);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
}


/* "singlemulti" test
 *
 * Create a profile/sequence/anchorset triplet that has only one possible
 * path with $P(\pi | x, A, \theta) = 1$, when using a glocal mode.
 *
 * Unlike the "singlesingle" test, now we can use multihit mode.
 *
 * See <h4_modelsample_SinglePathedASC()> for explanation of how the
 * case is constructed.
 * 
 * Now trace score = ASC Fwd = ASC Bck score.
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
  H4_REFMX     *afu       = h4_refmx_Create(100,100);  
  H4_REFMX     *afd       = h4_refmx_Create(100,100);  
  float         tsc, vsc, fsc, asc_f;
  float         abstol    = 0.0001;

  if ( h4_modelsample_SinglePathASC(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &tsc)  != eslOK) esl_fatal(failmsg);

  if ( h4_reference_Viterbi(sq->dsq, sq->n, hmm, mo, rxv, vpi, &vsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Forward(sq->dsq, sq->n, hmm, mo, rxf,      &fsc) != eslOK) esl_fatal(failmsg);

  if ( h4_reference_asc_Forward(sq->dsq, sq->n, hmm, mo, anch, afu, afd, &asc_f) != eslOK) esl_fatal(failmsg);

#if 0
  h4_refmx_Dump(stdout, rxf);
  h4_refmx_Dump(stdout, afu);
  h4_refmx_Dump(stdout, afd);

  h4_path_Dump(stdout, pi);
  h4_path_Dump(stdout, vpi);
  h4_anchorset_Dump(stdout, anch);

  printf("tsc   = %.3f\n", tsc);
  printf("asc_f = %.3f\n", asc_f);
  printf("diff  = %.3f\n", tsc - asc_f);
#endif

  if ( esl_FCompareNew(tsc, asc_f, 0.0, abstol) != eslOK) esl_fatal(failmsg);

  h4_refmx_Destroy(afd);
  h4_refmx_Destroy(afu);
  h4_refmx_Destroy(rxf);
  h4_refmx_Destroy(rxv);
  h4_anchorset_Destroy(anch);
  h4_path_Destroy(pi);
  h4_path_Destroy(vpi);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
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
 * 
 * And also, obviously:
 *     2. Viterbi score >= score of sampled trace that created sq
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
  H4_REFMX     *afu        = h4_refmx_Create(100,100);  
  H4_REFMX     *afd        = h4_refmx_Create(100,100);  
  float         tsc, vsc, fsc, asc_f;
  float         fwd_atol;
  float         vit_atol;

  ESL_DASSERT1(( which >= 0 && which <= 2 ));

  switch (which) {
  case 0:  failmsg = failmsg0;  fwd_atol = (h4_logsum_IsSlowExact() ? 0.001  : 0.002);  vit_atol = 1e-5;  break;
  case 1:  failmsg = failmsg1;  fwd_atol = (h4_logsum_IsSlowExact() ? 0.001  : 0.002);  vit_atol = 1e-5;  break;
  case 2:  failmsg = failmsg2;  fwd_atol = (h4_logsum_IsSlowExact() ? 0.0001 : 0.01);   vit_atol = 1e-5;  break;
  }
  
  switch (which) {
  case 0: if ( h4_modelsample_AnchoredUni  (rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &tsc)  != eslOK) esl_fatal(failmsg); break;
  case 1: if ( h4_modelsample_AnchoredLocal(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &tsc)  != eslOK) esl_fatal(failmsg); break;
  case 2: if ( h4_modelsample_AnchoredMulti(rng, abc, M, &hmm, &sq, &anch, &mo, &pi, &tsc)  != eslOK) esl_fatal(failmsg); break;
  }

  if ( h4_reference_Viterbi(sq->dsq, sq->n, hmm, mo, rxv, vpi, &vsc) != eslOK) esl_fatal(failmsg);
  if ( h4_reference_Forward(sq->dsq, sq->n, hmm, mo, rxf,      &fsc) != eslOK) esl_fatal(failmsg);

  if ( h4_reference_asc_Forward(sq->dsq, sq->n, hmm, mo, anch, afu, afd, &asc_f) != eslOK) esl_fatal(failmsg);
  
#if 0
  h4_profile_Dump(stdout, hmm);
  h4_mode_Dump(stdout, mo);

  esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, FALSE);

  h4_refmx_Dump(stdout, rxf);
  h4_refmx_Dump(stdout, afu);
  h4_refmx_Dump(stdout, afd);

  h4_path_Dump(stdout, pi);
  h4_path_Dump(stdout, vpi);
  h4_anchorset_Dump(stdout, anch);

  printf("tsc   = %.3f\n", tsc);
  printf("vsc   = %.3f\n", vsc);
  printf("fsc   = %.3f\n", fsc);
  printf("asc_f = %.3f\n", asc_f);
  printf("diff  = %.3f\n", fsc - asc_f);
#endif

  if (esl_FCompareNew(fsc, asc_f, 0.0, fwd_atol) != eslOK) esl_fatal(failmsg);
  if (tsc > vsc+vit_atol)                                  esl_fatal(failmsg);
  if (vsc > fsc)                                           esl_fatal(failmsg);

  h4_refmx_Destroy(afu);
  h4_refmx_Destroy(afd);
  h4_refmx_Destroy(rxf);
  h4_refmx_Destroy(rxv);
  h4_anchorset_Destroy(anch);
  h4_path_Destroy(vpi);
  h4_path_Destroy(pi);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
}
#endif // h4REFERENCE_ASC_TESTDRIVE

/*****************************************************************
 * x. Test driver
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
 * x. Example
 *****************************************************************/

#ifdef h4REFERENCE_ASC_EXAMPLE

#include "h4_config.h"

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
static char banner[] = "example of ASC Forward reference implementation";

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
