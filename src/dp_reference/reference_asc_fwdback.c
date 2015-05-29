/* Reference implementation of anchor set constrained (ASC) Forward
 * and Backward.
 * 
 * All reference implementation code is for development and
 * testing. It is not used in HMMER's main executables. Production
 * code uses sparse dynamic programming.
 * 
 * Contents:
 *    1. ASC Forward
 *    2. ASC Backward
 *    3. Footnotes
 *    4. Unit tests
 *    5. Test driver
 *    6. Example
 *    7. Copyright and license information
 */


#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_anchors.h"

#include "misc/logsum.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_asc_fwdback.h"

/*****************************************************************
 * 1. ASC Forward
 *****************************************************************/


/* Function:  p7_ReferenceASCForward()
 * Synopsis:  Calculate an anchor-set-constrained (ASC) Forward score.
 *
 * Purpose:   The anchor set constrained (ASC) Forward algorithm. Given
 *            digital sequence <dsq> of length <L>, profile <gm> to
 *            compare it to, and an anchor set <anch> for <D> domains;
 *            calculate ASC Forward score, and return it in
 *            <*opt_sc>.
 *            
 *            Caller provides two DP matrices <mxu> and <mxd>. They
 *            can be of any allocated size, and they will be
 *            reallocated if needed. Upon return, <mxu> and <mxd>
 *            contain the ASC Forward DP matrices for the UP and DOWN
 *            sectors of the calculation, respectively. In domain
 *            analysis, they will be needed later for posterior
 *            decoding.
 *
 *            Caller must have initialized at least once (per program
 *            invocation) with a <p7_FLogsumInit()> call, because this
 *            function uses <p7_FLogsum()>.
 *            
 * Args:      dsq    : digital target sequence 1..L
 *            L      : length of <dsq>
 *            gm     : profile
 *            anch   : array of (i,k) anchors defining <dsq>'s domain structure (0) 1..D (D+1) with sentinels
 *            D      : number of anchors in <anch> array -- # of domains
 *            mxu    : UP matrix
 *            mxd    : DOWN matrix
 *            opt_sc : optRETURN - ASC Forward lod score, in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_ReferenceASCForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D,
		       P7_REFMX *mxu, P7_REFMX *mxd, float *opt_sc)
{
  const float *tsc;		/* ptr for stepping thru profile's transition parameters */
  const float *rsc;		/* ptr for stepping thru profile's emission parameters   */
  float *dpp;                   /* ptr into main states ({MID}{LG}) for prev row...      */
  float *dpc;			/*   ... and a current row.                              */
  float *xp;			/* ptr into specials (ENJBLGC) for a previous row...     */
  float *xc;			/*   ... and a current row.                              */
  int   d;			/* counter over domains: 0..D-1                          */
  int   i;			/* counter over sequence positions (rows): 0,1..L        */
  int   k;			/* counter over model positions (columns): 0,1..M        */
  int   s;			/* counter over states                                   */
  float mlv, mgv;		/* tmp variables for ML, MG scores...                    */
  float dlv, dgv;		/*   ... and for DL, DG scores                           */
  float xE;			/* tmp var for accumulating E score for a row            */
  int   M = gm->M;		/* for a bit more clarity, less dereference clutter      */
  int   status;

  /* Don't try to contract check the <gm> length model like we do elsewhere.
   * Envelope determination calls ASC Forward on arbitrary subsequences to
   * recalculate envelope scores.
   */

  /* reallocation, if needed */
  if ( (status = p7_refmx_GrowTo(mxu, M, L)) != eslOK) return status;
  if ( (status = p7_refmx_GrowTo(mxd, M, L)) != eslOK) return status;
  mxu->M    = mxd->M    = M;
  mxu->L    = mxd->L    = L;
  mxu->type = p7R_ASC_FWD_UP;
  mxd->type = p7R_ASC_FWD_DOWN;

  /* Initialize i=0..anch[0].i-1 specials. 
   * All specials are stored in DOWN matrix.
   */
  xp   = NULL;
  for (i = 0; i < anch[1].i0; i++)   // Note this even works for D=0/L=0 case, because then anch[D+1] = L+1 sentinel, and we'll do row 0 (only) here
    {
      xc = mxd->dp[i] + (M+1) * p7R_NSCELLS; 
      xc[p7R_E]  = -eslINFINITY;
      xc[p7R_N]  = (i == 0 ? 0. : xp[p7R_N] + gm->xsc[p7P_N][p7P_LOOP]);
      xc[p7R_J]  = -eslINFINITY;
      xc[p7R_B]  = xc[p7R_N] + gm->xsc[p7P_N][p7P_MOVE];
      xc[p7R_L]  = xc[p7R_B] + gm->xsc[p7P_B][0]; 
      xc[p7R_G]  = xc[p7R_B] + gm->xsc[p7P_B][1]; 
      xc[p7R_C]  = -eslINFINITY;
      xc[p7R_JJ] = -eslINFINITY;
      xc[p7R_CC] = -eslINFINITY;
      xp = xc;
    }

  /* Iterate over domains d=1..D: */
  for (d = 1; d <= D; d++)
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
       *    decoding, for delete states in G->D1DDDk-1->Mk wing unfolding.
       *    So here, in fwd, we'll initialize row i0(d-1) to all -inf,
       *    which somewhat simplifies the boundary condition for the next row.
       *    
       *    There's a potentially tricky case when i0(d) = i0(d-1)+1, anchors
       *    on adjacent rows. We need to connect the lower right UP supercell
       *    to the anchor supercell in DOWN; we do this by making sure <dpc>
       *    has a well-defined final state, at i0(d-1),k0. When we go to 
       *    connect to DOWN, we'll back it up one supercell, setting dpp = dpc - p7R_NSCELLS.
       *                                                                         
       *    Because the REFMX has a column k=0 allocated for us, we may 
       *    as well use it; so we also initialize supercells in k=0 at
       *    every UP row to -inf, again to somewhat simplify boundary conditions.
       *    
       *    It's possible to have no cells in the UP matrix at all, if
       *    k0(d) = 1. In that case, the way we handle initialization,
       *    the combination of the above two points is important: 
       *    we initialize k=0 on the previous row, and dpc is left hovering
       *    just after it.
       *****************************************************************/

      /* Initialization of top row, i0(d-1), to -inf's */
      dpc = mxu->dp[ anch[d-1].i0 ];                                          // that's dp[0] for d=1 
      for (s = 0; s < p7R_NSCELLS * anch[d].k0; s++) *dpc++ = -eslINFINITY;   
      /* in the special case of i0(d)=i0(d-1)+1,
       * dpc is now sitting on i0(d)-1, k0(d).k): supercell above the anchor, 
       * where it needs to be to make the connection to DOWN.
       */

      /* Now we recurse for remaining rows, i0(d-1)+1..i0[d].i-1 row.
       * (It's possible that no such rows exist, depending where that next anchor is.)
       */
      for (i = anch[d-1].i0+1 ; i < anch[d].i0; i++)
	{
	  rsc = gm->rsc[dsq[i]] + p7P_NR;                           // Start <rsc> at k=1 on row i 
	  tsc = gm->tsc;                                            // Start <tsc> at k=0, where off-by-one {LG}->M transition scores are
	  xp  = mxd->dp[i-1] + (M+1) * p7R_NSCELLS;                 // <xp> set to the specials of the previous row, i-1
	  dpp = mxu->dp[i-1];                                       // Start <dpp> on k=0, which we know has been init'ed to -inf above.
	  dpc = mxu->dp[i];                                         // Start <dpc> on k=0 of row i...
	  for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;  //   ... initialize that cell, and now <dpc> is on k=1.
	  dlv = dgv = -eslINFINITY;

	  for (k = 1; k < anch[d].k0; k++)
	    {
	      mlv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_ML) + *(tsc + p7P_MM),
							   *(dpp+p7R_IL) + *(tsc + p7P_IM)),
						p7_FLogsum(*(dpp+p7R_DL) + *(tsc + p7P_DM),
							      xp[p7R_L]  + *(tsc + p7P_LM)));
	      mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_MG) + *(tsc + p7P_MM),
							   *(dpp+p7R_IG) + *(tsc + p7P_IM)),
						p7_FLogsum(*(dpp+p7R_DG) + *(tsc + p7P_DM),
							      xp[p7R_G]  + *(tsc + p7P_GM)));

	      rsc++;              
	      tsc += p7P_NTRANS;  
	      dpp += p7R_NSCELLS;	

	      *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_ML) + *(tsc + p7P_MI), *(dpp + p7R_IL) + *(tsc + p7P_II));
	      *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_MG) + *(tsc + p7P_MI), *(dpp + p7R_IG) + *(tsc + p7P_II));
	      rsc++;	

	      *dpc++ = dlv;
	      *dpc++ = dgv;
	      dlv = p7_FLogsum( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD));
	      dgv = p7_FLogsum( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD));
	    }
	}

      /* The very last cell we calculated was the cell diagonal from
       * the anchor (i.e. where the DOWN matrix is going to start), and
       * we want that diagonal cell for initializing DOWN.
       * <dpc> just stepped past our desired cell; step back.
       * This works even if the UP matrix was empty.
       */
      dpp = dpc - p7R_NSCELLS;

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

      /* Start with anch[d].k-1 on first row, and set all cells to -inf*/
      i   = anch[d].i0;
      tsc = gm->tsc +         (anch[d].k0-1) * p7P_NTRANS;    // Start <tsc> on k0-1, i.e. k-1 relative to start of calculation
      rsc = gm->rsc[dsq[i]] + (anch[d].k0    * p7P_NR);       // <rsc> is on scores for k0
      xp  = mxd->dp[i-1] + (M+1) * p7R_NSCELLS;               // <xp> on specials for anch.i-1
      dpc = mxd->dp[i] + (anch[d].k0-1) * p7R_NSCELLS;        // We will initialize k0(d)-1 supercell...
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;//   ... and now <dpc> now sits on k0(d)

      /* Then calculate the anchor cell (i0,k0) as an UP calc, using
       * <dpp> which was already set, above. 
       */
      mlv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_ML) + *(tsc + p7P_MM),
						   *(dpp+p7R_IL) + *(tsc + p7P_IM)),
					p7_FLogsum(*(dpp+p7R_DL) + *(tsc + p7P_DM),
						      xp[p7R_L]  + *(tsc + p7P_LM)));
      mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_MG) + *(tsc + p7P_MM),
						   *(dpp+p7R_IG) + *(tsc + p7P_IM)),
					p7_FLogsum(*(dpp+p7R_DG) + *(tsc + p7P_DM),
						      xp[p7R_G]  + *(tsc + p7P_GM)));

      tsc   += p7P_NTRANS;
      *dpc++ = -eslINFINITY;	/* IL */
      *dpc++ = -eslINFINITY;	/* IG */
      *dpc++ = -eslINFINITY;	/* DL */
      *dpc++ = -eslINFINITY;	/* DG */
      dlv = mlv + *(tsc + p7P_MD);
      dgv = mgv + *(tsc + p7P_MD);

      /* xE initialization counts exits from the anchor cell we calculated.
       * Unlike the rest of the top row, MG/ML exits from the anchor cell
       * need to be calculated. Also, it has to watch out for the
       * glocal exit case when the anchor cell (unusually) sits on M
       * itself.
       */
      xE  = (anch[d].k0 == M ? p7_FLogsum(mlv, mgv) : mlv);

      /* Now we initialize the rest of the top row i0(d) from k=k0(d)+1 to M,
       * which is only reachable on deletion paths from the anchor.
       */
      for (k = anch[d].k0+1; k <= M; k++)
	{
	  *dpc++ = -eslINFINITY; // ML. No entry, and unreachable from other cells too. 
	  *dpc++ = -eslINFINITY; // MG. Ditto.
	  *dpc++ = -eslINFINITY; // IL. Not reachable on top row. 
	  *dpc++ = -eslINFINITY; // IG. Ditto.
	  *dpc++ = dlv;          // DL. Customary delayed store of prev calculation.
	  *dpc++ = dgv;          // DG. Ditto.

	  tsc   += p7P_NTRANS;
	  
	  xE  = (k == M ?                                  // Glocal exit included if k==M.
		 p7_FLogsum( xE, p7_FLogsum( dlv, dgv)) :  // We know all non-anchor-cell M's are -inf on top row, so 
		 p7_FLogsum( xE, dlv));			   //  we don't include M in these sums.
	  
	  dlv    = dlv + *(tsc + p7P_DD);
	  dgv    = dgv + *(tsc + p7P_DD);
	}

      /* dpc now sits on the start of the specials, in mxd */
      xc = dpc;
      xc[p7R_E]  = xE;
      xc[p7R_N]  = -eslINFINITY; 
      xc[p7R_J]  = xc[p7R_E] + gm->xsc[p7P_E][p7P_LOOP];
      xc[p7R_B]  = xc[p7R_J] + gm->xsc[p7P_J][p7P_MOVE];
      xc[p7R_L]  = xc[p7R_B] + gm->xsc[p7P_B][0]; 
      xc[p7R_G]  = xc[p7R_B] + gm->xsc[p7P_B][1]; 
      xc[p7R_C]  = xc[p7R_E] + gm->xsc[p7P_E][p7P_MOVE];
      xc[p7R_JJ] = -eslINFINITY;
      xc[p7R_CC] = -eslINFINITY;

      /* Now we can do the remaining rows in the Down sector of domain d. */
      for (i = anch[d].i0+1 ; i < anch[d+1].i0; i++)               // sentinel i0(D+1) = L+1 
	{
	  rsc = gm->rsc[dsq[i]] + anch[d].k0     * p7P_NR;         // Start <rsc> on (x_i, anchor_k, MAT) */
	  tsc = gm->tsc         + (anch[d].k0-1) * p7P_NTRANS;	   // Start <tsc> on (anchor_k-1), to pick up LMk,GMk entries 
	  dpp = mxd->dp[i-1]    + (anch[d].k0-1) * p7R_NSCELLS;	   // Start <dpp> on (i-1, anchor_k-1) 
	  dpc = mxd->dp[i]      + (anch[d].k0-1) * p7R_NSCELLS;	   // Start <dpc> on (i, anchor_k-1)... 
	  for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; //  ... and initialize the k-1 cells to -inf... 
                                                           	   //  ... so, now dpc is on anchor_k.
	  dlv = dgv = xE = -eslINFINITY;

  	  for (k = anch[d].k0; k <= M; k++) 
	    {				  
	      mlv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_ML) + *(tsc + p7P_MM),
							   *(dpp+p7R_IL) + *(tsc + p7P_IM)),
						           *(dpp+p7R_DL) + *(tsc + p7P_DM));
	      mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_MG) + *(tsc + p7P_MM),
							   *(dpp+p7R_IG) + *(tsc + p7P_IM)),
					 	           *(dpp+p7R_DG) + *(tsc + p7P_DM));

	      rsc++;                // rsc advances to insert score for position k 
	      tsc += p7P_NTRANS;    // tsc advances to transitions in states k     
	      dpp += p7R_NSCELLS;   // dpp advances to cells for states k          

	      *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_ML) + *(tsc + p7P_MI), *(dpp + p7R_IL) + *(tsc + p7P_II));
	      *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_MG) + *(tsc + p7P_MI), *(dpp + p7R_IG) + *(tsc + p7P_II));
	      rsc++;		    // rsc advances to next match state emission  

	      xE  = (k == M ?
		     p7_FLogsum( xE, p7_FLogsum( p7_FLogsum(mlv, dlv), p7_FLogsum(mgv, dgv))) : // k=M includes glocal exits  
		     p7_FLogsum( xE, p7_FLogsum(mlv, dlv)));                                    // k<M allows local exit only 

	      *dpc++ = dlv;                                                    // DL. Customary delayed store.
	      *dpc++ = dgv;                                                    //   ... ditto for DG store.
	      dlv = p7_FLogsum( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD)); // Precalculation of DL for next k.
	      dgv = p7_FLogsum( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD)); //   ... ditto for DG calculation.
	    }

	  /*****************************************************************
	   *  Having finished and stored the DOWN calculation on row i, with value xE,
           *  we can calculate and store the specials - also in the DOWN matrix.
           *    dpc[] is already on the special state storage.
	   *****************************************************************/

	  xc = dpc;
	  xp = dpp + p7R_NSCELLS;	
	  xc[p7R_E]  = xE;		
	  xc[p7R_N]  = -eslINFINITY; 
	  xc[p7R_J]  = p7_FLogsum( xp[p7R_J] + gm->xsc[p7P_J][p7P_LOOP], xc[p7R_E] + gm->xsc[p7P_E][p7P_LOOP] );
	  xc[p7R_B]  = xc[p7R_J] + gm->xsc[p7P_J][p7P_MOVE]; 
	  xc[p7R_L]  = xc[p7R_B] + gm->xsc[p7P_B][0]; 
	  xc[p7R_G]  = xc[p7R_B] + gm->xsc[p7P_B][1]; 
	  xc[p7R_C]  = p7_FLogsum( xp[p7R_C] + gm->xsc[p7P_C][p7P_LOOP], xc[p7R_E] + gm->xsc[p7P_E][p7P_MOVE] );
	  xc[p7R_JJ] = -eslINFINITY;                                                                           
	  xc[p7R_CC] = -eslINFINITY;       
	} /* end loop over rows i of DOWN sector for domain d */

    } /* end loop over domains d=0..D-1; DP calculation complete. */

  /* As we leave the DP recursion, <xc> is still sitting on the
   * special states for the last row L... even for the edge case
   * of D==0 (and the edge case L=0 which must also have D==0).
   */
  if (opt_sc) *opt_sc = xc[p7R_C] + gm->xsc[p7P_C][p7P_MOVE]; /* C->T */
  return eslOK;
}
/*-------------------- end, ASC Forward -------------------------*/

/*****************************************************************
 * 2. ASC Backward
 *****************************************************************/

/* Function:  p7_ReferenceASCBackward()
 * Synopsis:  Calculate an anchor-set-constrained (ASC) Backward score.
 *
 * Purpose:   The anchor set constrained (ASC) Backward algorithm.
 *            Given digital sequence <dsq> of length <L>, profile <gm> to
 *            compare it to, and an anchor set <anch> for <D> domains;
 *            calculate ASC Backward score, and return it in
 *            <*opt_sc>.
 *            
 *            Caller provides two DP matrices <abu> and <abd>. They
 *            can be of any allocated size, and they will be
 *            reallocated if needed. Upon return, <abu> and <abd>
 *            contain the ASC Backward DP matrices for the UP and DOWN
 *            sectors of the calculation, respectively. In domain
 *            analysis, they will be needed later for posterior
 *            decoding.
 * 
 *            Caller must have initialized at least once (per program
 *            invocation) with a <p7_FLogsumInit()> call, because this
 *            function uses <p7_FLogsum()>.
 *            
 *            <anch> and <D> might be data in a <P7_ANCHORS> list
 *            management container: for example, for <P7_ANCHORS *dom>,
 *            you would pass <dom->a> and <dom->D>.
 *
 * Args:      dsq    : digital target sequence 1..L
 *            L      : length of <dsq>
 *            gm     : profile
 *            anch   : array of (i0,k0) anchors defining <dsq>'s domain structure; (0) 1..D (D+1) with sentinels
 *            D      : number of anchors in <anch> array -- # of domains
 *            abu    : UP matrix
 *            abd    : DOWN matrix
 *            opt_sc : optRETURN - ASC Backward lod score, in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_ReferenceASCBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D,
			P7_REFMX *abu, P7_REFMX *abd, float *opt_sc)
{
  int    M = gm->M;
  const float *tsc;		/* ptr into transition scores of <gm> */
  const float *rsc;		/* ptr into emission scores of <gm> for residue dsq[i] on current row i  */
  const float *rsn;		/* ptr into emission scores of <gm> for residue dsq[i+1] on next row i+1 */
  float *dpc, *dpn;		/* ptrs into DP matrix for current row i, next row i+1  */
  float *xc;			/* ptr to specials on current row; specials are stored in DOWN, <abd> */
  int    d;                   	/* counter over domains 0..D-1 */
  int    i;			/* counter over sequence positions 0.1..L (DP rows) */
  int    k;			/* counter over model positions 0.1..M (DP columns) */
  int    s;
  int    iend;
  float  mgc, mlc;
  float  mgn, mln;
  float  dgn, dln;
  float  ign, iln;
  float  xE;
  float  xG,  xL;
  float  xC, xJ, xN;
  int    status;

  /* contract checks / arg validation */
  ESL_DASSERT1( ( gm->L == L || gm->L == 0) ); /* length model in profile is either L (usually) or 0 (some unit tests) */

  /* reallocation, if needed */
  if ( (status = p7_refmx_GrowTo(abu, gm->M, L)) != eslOK) return status;
  if ( (status = p7_refmx_GrowTo(abd, gm->M, L)) != eslOK) return status;
  abu->M    = abd->M    = M;
  abu->L    = abd->L    = L;
  abu->type = p7R_ASC_BCK_UP;
  abd->type = p7R_ASC_BCK_DOWN;

  for (i = L; i >= anch[D].i0; i--)     // L=0,D=0 case: this would be 0..0 (thanks to i0(0)=0 sentinel)
    {
      xc = abd->dp[i] + (M+1) * p7R_NSCELLS;
      xc[p7R_CC] = -eslINFINITY;
      xc[p7R_JJ] = -eslINFINITY;
      xc[p7R_C]  = xC = (i == L ? gm->xsc[p7P_C][p7P_MOVE] : xC + gm->xsc[p7P_C][p7P_LOOP]);
      xc[p7R_G]  = -eslINFINITY;
      xc[p7R_L]  = -eslINFINITY;
      xc[p7R_B]  = -eslINFINITY;
      xc[p7R_J]  = -eslINFINITY;
      xc[p7R_N]  = -eslINFINITY;
      xc[p7R_E]  = xC + gm->xsc[p7P_E][p7P_MOVE];  // cppcheck thinks xC can be uninitialized here, but it can't; xC is initialized at i==L, 1st pass thru loop
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
      iend = anch[d+1].i0-1;                                            // <i==iend> tests used below for boundary conditions on last DOWN row
      for (i = iend; i >= anch[d].i0; i--)                              // at d=D, this is (L down to i0(D)), thanks to i0(D+1)=L+1 sentinel
	{
	  rsn  = (i == iend ? NULL : gm->rsc[dsq[i+1]] + M * p7P_NR);   // residue scores on next row; start at M
	  tsc  = gm->tsc + M * p7P_NTRANS;                              // transition scores: start at M
	  dpc  = abd->dp[i]   + M * p7R_NSCELLS;                        // current row of DP matrix: start at M
	  xE   = dpc[p7R_NSCELLS+p7R_E];                                // pick up the xE score; specials start at M+1, hence the p7R_NSCELLS bump here
	  dpn  = (i == iend ? NULL : abd->dp[i+1] + M * p7R_NSCELLS);   // next row of DP matrix: start at M

	  mgn = dgn = -eslINFINITY;
	  mln = dln = -eslINFINITY;
	  ign = iln = -eslINFINITY;
	  xL  = -eslINFINITY;

	  for (k = M; k >= anch[d].k0; k--)
	    {
	      if (i != iend) {           // in one-row memory-efficient dp, dpc could be same as dpn, so:
		ign = dpn[p7R_IG];       // pick up I scores, before storing anything in these cells
		iln = dpn[p7R_IL];       // if insert scores were non-zero, we would add rsn[p7R_I] here
	      }

	      /* M calculations. Storage deferred for one-row reasons. */
	      mgc = (k == M ? xE :                                // at k=M, MG->E is possible, and it happens to be the only transition that's possible
		     p7_FLogsum( p7_FLogsum(mgn + tsc[p7P_MM],    // mgn (+ *rsn) was picked up in last k loop, so now it's i+1,k+1
					    ign + tsc[p7P_MI]),   // ign was just picked up, so it's i+1,k
        				    dgn + tsc[p7P_MD]));  // dgn is remembered from prev loop, so now it's i,k+1
	      mlc =  p7_FLogsum( p7_FLogsum(mln + tsc[p7P_MM],   
					    iln + tsc[p7P_MI]),
				 p7_FLogsum(dln + tsc[p7P_MD],
					    xE));
		     
	      dpc[p7R_DG] = dgn = (k == M ?  xE :  p7_FLogsum(mgn + tsc[p7P_DM], dgn + tsc[p7P_DD]));
	      dpc[p7R_DL] = dln =  p7_FLogsum( xE, p7_FLogsum(mln + tsc[p7P_DM], dln + tsc[p7P_DD]));
	      
	      dpc[p7R_IG] = p7_FLogsum(mgn + tsc[p7P_IM], ign + tsc[p7P_II]);
	      dpc[p7R_IL] = p7_FLogsum(mln + tsc[p7P_IM], iln + tsc[p7P_II]);
	      
	      if (i != iend) {              // pick up M[i+1][k] values, add residue emission to them;
		mgn =  dpn[p7R_MG] + *rsn;  // when we loop around, these become M[i+1][k+1] values we need for DP
		mln =  dpn[p7R_ML] + *rsn;
		rsn -= p7P_NR;
		dpn -= p7R_NSCELLS;
	      } 

	      dpc[p7R_MG] = mgc;           // now that we've picked up mgn/mln, safe to store MG,ML
	      dpc[p7R_ML] = mlc;
	      
	      tsc -= p7P_NTRANS;           
	      dpc -= p7R_NSCELLS;
	    }
	}
      /* mgc/mlc are the scores in the anchor cell anch[d].i,k */
      /* tsc is on anch[d].k-1 */
      rsc = gm->rsc[ dsq[anch[d].i0]] + anch[d].k0 * p7P_NR;
      mgn = mgc + *rsc;  
      mln = mlc + *rsc;
      xG  = mgn + tsc[p7P_GM];
      xL  = mln + tsc[p7P_LM];

      xJ = xN = -eslINFINITY;

      /* UP matrix */
      iend = anch[d].i0-1;                              // <i == iend> tests for boundary condition on last UP row, shorthand for <i == anch[d].i0-1>
      for (i = iend; i > anch[d-1].i0; i--)             // at d=1, this is i0(1)-1 down to 1, thanks to i0(0)=0 sentinel
	{
	  xc = abd->dp[i] + (M+1) * p7R_NSCELLS;        // on specials, which are in DOWN matrix
	  xc[p7R_CC] = -eslINFINITY;                    // CC,JJ are only used in decoding matrices
	  xc[p7R_JJ] = -eslINFINITY;
	  xc[p7R_C]  = -eslINFINITY;                    // C is now unreachable, when anchor set constrained.
	  xc[p7R_G]  = xG;                              // xG was accumulated during prev row; G->Mk wing unfolded
	  xc[p7R_L]  = xL;                              // xL accumulated on prev row
	  xc[p7R_B]  = p7_FLogsum(xG + gm->xsc[p7P_B][1],  xL + gm->xsc[p7P_B][0]); 
	  xc[p7R_J]  = xJ = p7_FLogsum(xJ + gm->xsc[p7P_J][p7P_LOOP], xc[p7R_B] + gm->xsc[p7P_J][p7P_MOVE]);
	  xc[p7R_N]  = xN = p7_FLogsum(xN + gm->xsc[p7P_N][p7P_LOOP], xc[p7R_B] + gm->xsc[p7P_N][p7P_MOVE]);
	  xc[p7R_E]  = xc[p7R_J] + gm->xsc[p7P_E][p7P_LOOP];  

	  tsc = gm->tsc    + (anch[d].k0-1) * p7P_NTRANS;                                    // transition scores: start at anch[d].k-1
	  dpc = abu->dp[i] + (anch[d].k0-1) * p7R_NSCELLS;                                   // on anch[d].k-1
	  dpn = (i == anch[d].i0-1 ? NULL : abu->dp[i+1] + (anch[d].k0 - 1) * p7R_NSCELLS);  // on anch[d].k-1
	  rsc = gm->rsc[dsq[i]] + (anch[d].k0-1) * p7P_NR;
	  rsn = (i == anch[d].i0-1 ? NULL : gm->rsc[dsq[i+1]] + (anch[d].k0-1) * p7P_NR);

	  xG  = xL  = -eslINFINITY; 
	  dgn = dln = -eslINFINITY;
	  ign = iln = -eslINFINITY;
	  if (i < iend) mgn = mln = -eslINFINITY;       // not on iend, because there we allow mgn/mln to carry over from anchor cell 

	  /* The recursion is the same as for the DOWN matrix, so only differences are commented on: */
	  for (k = anch[d].k0-1; k >= 1; k--)
	    {
	      if (i < iend) {           
		ign = dpn[p7R_IG];       
		iln = dpn[p7R_IL];       
	      }

	      /* M calculations include no E contributions: can't do M->E in UP matrix */
	      mgc =  p7_FLogsum( p7_FLogsum(mgn + tsc[p7P_MM],    
					    ign + tsc[p7P_MI]),   
				            dgn + tsc[p7P_MD]);
	      mlc =  p7_FLogsum( p7_FLogsum(mln + tsc[p7P_MM],   
					    iln + tsc[p7P_MI]),
				            dln + tsc[p7P_MD]);
		     
	      xG = p7_FLogsum(xG, mgc + *rsc + tsc[p7P_GM - p7P_NTRANS]);
	      xL = p7_FLogsum(xL, mlc + *rsc + tsc[p7P_LM - p7P_NTRANS]);
	      rsc -= p7P_NR;

	      /* same for no D->E contributions in UP matrix */
	      dpc[p7R_DG] = dgn =  p7_FLogsum(mgn + tsc[p7P_DM], dgn + tsc[p7P_DD]);
	      dpc[p7R_DL] = dln =  p7_FLogsum(mln + tsc[p7P_DM], dln + tsc[p7P_DD]);
	      
	      dpc[p7R_IG] = p7_FLogsum(mgn + tsc[p7P_IM], ign + tsc[p7P_II]);
	      dpc[p7R_IL] = p7_FLogsum(mln + tsc[p7P_IM], iln + tsc[p7P_II]);
	      
	      if (i < iend) {       
		mgn =  dpn[p7R_MG] + *rsn;  // when we loop around, these become M[i+1][k+1] values we need for DP
		mln =  dpn[p7R_ML] + *rsn;
		rsn -= p7P_NR;
		dpn -= p7R_NSCELLS;
	      } else mgn = mln = -eslINFINITY;

	      dpc[p7R_MG] = mgc;           // now that we've picked up mgn/mln, safe to store MG,ML
	      dpc[p7R_ML] = mlc;
	      
	      tsc -= p7P_NTRANS;           
	      dpc -= p7R_NSCELLS;
	    }


	} /* end backwards loop over i for UP matrix d */
      /* i is now on anch[d-1].i, or 0 */

      dpc = abu->dp[i];
      for (s = 0; s < p7R_NSCELLS * anch[d].k0; s++) *dpc++ = -eslINFINITY;   

      xc = abd->dp[i] + (M+1) * p7R_NSCELLS;   // on specials, which are in DOWN matrix
      xc[p7R_CC] = -eslINFINITY;  // CC,JJ are only used in decoding matrices
      xc[p7R_JJ] = -eslINFINITY;
      xc[p7R_C]  = -eslINFINITY;  // C is now unreachable, when anchor set constrained.
      xc[p7R_G]  = xG;            // xG was accumulated during prev row; G->Mk wing unfolded
      xc[p7R_L]  = xL;            // xL accumulated on prev row
      xc[p7R_B]  = p7_FLogsum(xG + gm->xsc[p7P_B][1],  xL + gm->xsc[p7P_B][0]); 
      xc[p7R_J]  = xJ = p7_FLogsum(xJ + gm->xsc[p7P_J][p7P_LOOP], xc[p7R_B] + gm->xsc[p7P_J][p7P_MOVE]);
      xc[p7R_N]  = xN = p7_FLogsum(xN + gm->xsc[p7P_N][p7P_LOOP], xc[p7R_B] + gm->xsc[p7P_N][p7P_MOVE]);
      xc[p7R_E]  = xc[p7R_J] + gm->xsc[p7P_E][p7P_LOOP];  
    } /* end loop over domains d */

  if (opt_sc) *opt_sc = xN;
  return eslOK;
}


/*****************************************************************
 * 3. Footnotes
 ***************************************************************** 
 *
 * [1] xC VS xJ IN FORWARD; xN VS xJ IN BACKWARD
 * 
 * In a standard Forward algorithm, we can't look ahead to know
 * whether we're going to have an additional domain later in the
 * sequence, so both xC and xJ can have valid scores. In an ASC
 * algorithm, we do know, because we have the anchors -- so in
 * principle, we could calculate ASC Forward so that xC == -inf until
 * we hit the last anchor, then xJ == -inf after that. We don't do
 * that, because we can't always do it, and we don't want to create
 * unnecessary discrepancies with other Forward algorithms. Same goes
 * for xN/xJ in the backward direction.
 */


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7REFERENCE_ASC_FWDBACK_TESTDRIVE
#include "hmmer.h"

/* ascmatrix_compare()
 * 
 * Compare a standard DP matrix <std> to ASC UP and DOWN matrices
 * <ascu> and <ascd>, using anchor set <anch> for <D> domains.
 * Compare all valid values in the ASC matrices to their counterparts
 * in the standard matrix. If any absolute difference exceeds <tolerance>,
 * return <eslFAIL>. Else return <eslOK>.
 * 
 * Similar to p7_refmx_Compare(). See notes there for why we prefer
 * absolute, not relative difference (short version: error accumulates
 * more as an absolute # per term added; DP cells with values close to
 * zero may have large relative errors).
 */
static int
ascmatrix_compare(P7_REFMX *std, P7_REFMX *ascu, P7_REFMX *ascd, P7_ANCHOR *anch, int D, float tolerance)
{
  char  msg[] = "comparison of ASC DP F|B matrix to standard DP F|B matrix failed";
  int   M     = std->M;
  int   L     = std->L;
  float val;
  int   d,i,k,s;

  /* Contract checks */
  ESL_DASSERT1( (ascu->M == M && ascu->L == L));
  ESL_DASSERT1( (ascd->M == M && ascd->L == L));
  ESL_DASSERT1( ((std->type == p7R_FORWARD  && ascu->type == p7R_ASC_FWD_UP && ascd->type == p7R_ASC_FWD_DOWN) ||
		 (std->type == p7R_BACKWARD && ascu->type == p7R_ASC_BCK_UP && ascd->type == p7R_ASC_BCK_DOWN)));


  for (d = 1, i = 0; i <= L; i++)
    {
      if (i == anch[d].i0) d++;  // d = idx of current UP matrix; d-1 is current DOWN
      for (k = 1; k <= M; k++)   //   ... alternatively, you can think d = idx of next anchor we'll reach (1..D)
	for (s = 0; s < p7R_NSCELLS; s++)
	  {
	    val = -eslINFINITY;
	    if (k >= anch[d-1].k0) val = p7_FLogsum(val, P7R_MX(ascd,i,k,s)); // DOWN
	    if (k <  anch[d].k0)   val = p7_FLogsum(val, P7R_MX(ascu,i,k,s)); // UP
	    
	    if (val != -eslINFINITY && esl_FCompareAbs(val, P7R_MX(std,i,k,s), tolerance) != eslOK)
	      ESL_FAIL(eslFAIL, NULL, msg);
	  }

      for (s = 0; s < p7R_NXCELLS; s++)
	{
	  val = P7R_XMX(ascd,i,s);
	  if (val != -eslINFINITY && esl_FCompareAbs(val, P7R_XMX(std,i,s), tolerance) != eslOK)
	    ESL_FAIL(eslFAIL, NULL, msg);
	}
    }
  return eslOK;
}



/* The "singlepath" test uses p7_modelsample_SinglePathed() to sample
 * a contrived HMM that has only one possible path with probability 1;
 * that is, P(\pi | H) = 1. A uniglocal L=0 profile of that HMM also
 * has only one possible path. Any match state in that path is suitable
 * as an anchor, for comparison/testing of ASC algorithms. 
 *
 * Because there is literally only a single path, every cell in the
 * valid area of the ASC DP UP/DOWN matrices is identical (within 
 * floating point epsilon) to the standard DP matrix cells, for both
 * Forward/Backward. The same is not true of other unit tests, where
 * some prefixes/suffixes that ASC disallows have finite mass in standard Fwd/Bck but 
 * cannot reach finished paths.
 * 
 * For this situation:
 *   1. Viterbi = Fwd = Bck = ASC Fwd = ASC Bck scores.
 *   2. Viterbi trace == trace that generated the test sequence.
 *   3. Viterbi DP matrix cells == Forward DP matrix cells
 *   4. ASC Forward U/D matrix cells == std F matrix cells, within ASC regions.
 *   5. Ditto for ASC Backward cells compared to standard Backward.
 *   
 *****************************************************************
 * May fail stochastically. Analysis: SRE:J13/126
 *  ./reference_fwdback_utest --diag singlepath -N 10000
 * Field 1 = trace sc - vsc
 * Field 2 = trace sc - fsc
 * Field 3 = trace sc - bsc
 * Field 4 = trace sc - asc_f
 * Field 5 = trace_sc - asc_b                     
 * 
 * Because the trace score is summed in the same order as a Viterbi or
 * Forwards DP calculation, we expect exactly 0 difference for 1,2,4.
 * Backwards DP summation order, however, gives rise to roundoff
 * error.  Because this is a single-pathed test, we expect no
 * significant difference between default and exact logsum calculations.
 * 
 * Default:
 *    Zero difference for vsc, fsc, asc_f.
 *    trace_sc - bsc   : mean -2e-8 sd 1e-6 => 10 sigma and some => 1e-4
 *    trace_sc - asc_b : mean -2e-8 sd 1e-6 
 * 
 * Exact logsum:
 *    Zero difference for vsc, fsc, asc_f.
 *    trace_sc - bsc   : mean  2e-8 sd 1e-6 => 10 sigma and some => 1e-4
 *    trace_sc - asc_b : mean -2e-8 sd 1e-6
 *   
 */
static void
utest_singlepath(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback singlepath unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = p7_profile_Create(M, bg->abc);
  ESL_SQ     *sq        = esl_sq_CreateDigital(bg->abc);
  P7_TRACE   *tr        = p7_trace_Create();
  P7_TRACE   *vtr       = p7_trace_Create();
  P7_ANCHOR  *anch      = malloc(sizeof(P7_ANCHOR) * 3);  // *3, because of sentinels
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  int         nm, z, which;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       tolerance = 0.0001;
  int         status;
  
  if (gm  == NULL || sq  == NULL || tr  == NULL || bg  == NULL || vtr == NULL || rxv == NULL ||
      rxf == NULL || rxb == NULL || afu == NULL || afd == NULL || abu == NULL || abd == NULL || anch == NULL) 
    esl_fatal(failmsg);

  /* Sample single-pathed HMM, create uniglocal L=0 profile */
  if ((status = p7_modelsample_SinglePathed(rng, M, bg->abc, &hmm)) != eslOK) esl_fatal(failmsg);
  if ( p7_profile_ConfigUniglocal(gm, hmm, bg, 0)                   != eslOK) esl_fatal(failmsg); 

  /* Emit sequence and trace. Get the trace score. */
  if (p7_ProfileEmit(rng, hmm, gm, bg, sq, tr) != eslOK) esl_fatal(failmsg);
  if (p7_trace_Index(tr)                       != eslOK) esl_fatal(failmsg);
  if (p7_trace_Score(tr, sq->dsq, gm, &sc)     != eslOK) esl_fatal(failmsg);

  /* Create anchor set. We know it's uniglocal, so D=1. */
  D  = 1;
  nm = 0;
  for (z = 1; z < tr->N; z++) if (p7_trace_IsM(tr->st[z])) nm++;
  which = 1 + esl_rnd_Roll(rng, nm);
  for (z = 1; z < tr->N; z++) {
    if (p7_trace_IsM(tr->st[z])) {
      which--;
      if (which == 0) {
	anch[1].i0 = tr->i[z];
	anch[1].k0 = tr->k[z];
	break;
      }
    }
  }
  p7_anchor_SetSentinels(anch, D, tr->L, tr->M);

  if ((status = p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxv, vtr, &vsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceBackward(sq->dsq, sq->n, gm, rxb,      &bsc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceASCForward (sq->dsq, sq->n, gm, anch, D, afu, afd, &asc_f)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch, D, abu, abd, &asc_b)) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Vit:\n"); p7_refmx_Dump(stdout, rxv);
  //printf("### Reference Fwd:\n"); p7_refmx_Dump(stdout, rxf);
  //printf("### Reference Bck:\n"); p7_refmx_Dump(stdout, rxb);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);
  //printf("### ASC Bck UP:\n");    p7_refmx_Dump(stdout, abu);
  //printf("### ASC Bck DOWN:\n");  p7_refmx_Dump(stdout, abd);
  //p7_trace_DumpAnnotated(stdout, tr,  gm, sq->dsq);
  //p7_trace_DumpAnnotated(stdout, vtr, gm, sq->dsq);

  if (p7_trace_Compare(tr, vtr, 0)        != eslOK) esl_fatal(failmsg);  // generated trace = Viterbi trace

  if (! diagfp) 
    {
      if (esl_FCompareAbs(sc, vsc,   tolerance) != eslOK) esl_fatal(failmsg);  // generated trace score = Viterbi score
      if (esl_FCompareAbs(sc, fsc,   tolerance) != eslOK) esl_fatal(failmsg);  //  ... = Forward score
      if (esl_FCompareAbs(sc, bsc,   tolerance) != eslOK) esl_fatal(failmsg);  //  ... = Backward score
      if (esl_FCompareAbs(sc, asc_f, tolerance) != eslOK) esl_fatal(failmsg);  //  ... = ASC Forward score
      if (esl_FCompareAbs(sc, asc_b, tolerance) != eslOK) esl_fatal(failmsg);  //  ... = ASC Backward score

      /* to compare Viterbi to Fwd matrix, we have to hack around a safety check in the structures */
      rxv->type = p7R_FORWARD;
      if (p7_refmx_Compare(rxf, rxv, tolerance) != eslOK) esl_fatal(failmsg);
      rxv->type = p7R_VITERBI;

      if (ascmatrix_compare(rxf, afu, afd, anch, D, tolerance) != eslOK) esl_fatal(failmsg);
      if (ascmatrix_compare(rxb, abu, abd, anch, D, tolerance) != eslOK) esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g %20g %20g %20g\n", sc-vsc, sc-fsc, sc-bsc, sc-asc_f, sc-asc_b);
  
  free(anch);
  esl_sq_Destroy(sq);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(abu);
  p7_refmx_Destroy(abd);
  p7_refmx_Destroy(rxb);
  p7_refmx_Destroy(rxf);
  p7_refmx_Destroy(rxv);
  p7_trace_Destroy(vtr);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}



/* The "singlesingle" test (for single path, single domain) uses
 * p7_modelsample_SinglePathedSeq() to sample a contrived profile/seq
 * pair that has only one possible path with probability 1 (that is,
 * P(\pi | x, H) = 1). This must be a glocal single domain
 * path. Choose any anchor A in that domain. 
 * 
 * Then:
 *   1. Viterbi = Fwd = Bck = ASC Fwd = ASC Bck scores.
 *   2. Viterbi trace = trace that generated the test sequence. 
 * 
 *****************************************************************
 * May fail stochastically. Analysis: SRE:J13/126
 *  ./reference_fwdback_utest --diag singlesingle -N 10000
 * Field 1 = trace sc - vsc
 * Field 2 = trace sc - fsc
 * Field 3 = trace sc - bsc
 * Field 4 = trace sc - asc_f
 * Field 5 = trace_sc - asc_b                     
 * 
 * Because the trace score is summed in the same order as a Viterbi or
 * Forwards DP calculation, we expect exactly 0 difference for 1,2,4,
 * but we use a tolerance anyway because we haven't proven this
 * exact order is architecture independent. Backwards DP summation order, 
 * however, gives rise to roundoff error.  Because this is a single-pathed 
 * test, we expect no significant difference between default and exact
 * logsum calculations.
 * 
 * Default:
 *    Zero difference for vsc, fsc, asc_f.
 *    trace_sc - bsc   : mean -4e-8 sd 6e-6  range -0.0002 .. 0.0001 => 10x range => 0.001
 *    trace_sc - asc_b : mean -4e-8 sd 6e-6 
 * 
 * Exact logsum:
 *    Zero difference for vsc, fsc, asc_f.
 *    trace_sc - bsc   : mean -2e-8 sd 4e-6 range -0.0002 .. 0.0001  => 10x range => 0.001
 *    trace_sc - asc_b : mean -2e-8 sd 4e-6
 */
static void
utest_singlesingle(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback singlesingle unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_TRACE   *vtr       = p7_trace_Create();
  P7_ANCHOR  *anch      = NULL;
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       tolerance = 0.001;
  int         status;
  
  if (bg  == NULL || vtr == NULL || rxv == NULL ||
      rxf == NULL || rxb == NULL || afu == NULL ||
      afd == NULL || abu == NULL || abd == NULL)   esl_fatal(failmsg);

  if ((status = p7_modelsample_SinglePathedSeq(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceViterbi (dsq, L, gm, rxv, vtr, &vsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceForward (dsq, L, gm, rxf,      &fsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceBackward(dsq, L, gm, rxb,      &bsc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd, &asc_f)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd, &asc_b)) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Vit:\n"); p7_refmx_Dump(stdout, rxv);
  //printf("### Reference Fwd:\n"); p7_refmx_Dump(stdout, rxf);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);
  //p7_trace_DumpAnnotated(stdout, tr,  gm, dsq);
  //p7_trace_DumpAnnotated(stdout, vtr, gm, dsq);

  if (p7_trace_Compare(tr, vtr, 0)        != eslOK) esl_fatal(failmsg);

  if (! diagfp)
    {
      if (esl_FCompareAbs(sc, vsc,   tolerance) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(sc, fsc,   tolerance) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(sc, bsc,   tolerance) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(sc, asc_f, tolerance) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(sc, asc_b, tolerance) != eslOK) esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g %20g %20g %20g\n", sc-vsc, sc-fsc, sc-bsc, sc-asc_f, sc-asc_b);

  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(abu);
  p7_refmx_Destroy(abd);
  p7_refmx_Destroy(rxb);
  p7_refmx_Destroy(rxf);
  p7_refmx_Destroy(rxv);
  p7_trace_Destroy(vtr);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}



/* The "singlemulti" test (single path, multiple domains) samples
 * a contrived profile/seq/anchorset triplet that has only one 
 * path when anchor set constrained; that is, P(\pi | x, A, H) = 1.
 * This must be a glocal path, and may contain multiple domains.
 * 
 * Then:
 *   1. Trace score = ASC Fwd score = ASC Bck score.
 *   
 *****************************************************************
 * May fail stochastically. Analysis: SRE:J13/126.
 *   ./reference_fwdback_utest --diag singlemulti -N 10000
 * Field 1 = trace_sc - asc_f
 * Field 2 = trace_sc - asc_b
 * 
 * Trace score is summed in same order as ASC Forward, so we expect
 * exactly zero difference for (tsc - asc_f); however, because we haven't
 * proven this is architecture independent, we still use a tolerance.
 * Backwards is subject to roundoff error accumulation. Because this
 * is a singlepath test, we expect no significant difference between 
 * default and logsum calculations.
 * 
 * Default: 
 *   Zero difference for asc_f.
 *   tsc - asc_b:  mean 6e-8, sd 1.2e-5, range -0.0004 .. 0.0004 => 50x range => 0.002
 *   
 * Exact:
 *   Zero difference for asc_f
 *   tsc - asc_b: mean 3e-7, sd 1.3e-5, range -0.0004 .. 0.0005 => 50x range => 0.002
 */
static void
utest_singlemulti(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback singlemulti unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_ANCHOR  *anch      = NULL;
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, asc_f, asc_b;
  float       tolerance = 0.002;
  int         status;
  
  if (bg  == NULL || afu == NULL || afd == NULL || abu == NULL || abd == NULL) 
    esl_fatal(failmsg);

  if ((status = p7_modelsample_SinglePathedASC(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd, &asc_f)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd, &asc_b)) != eslOK) esl_fatal(failmsg);

  //p7_trace_DumpAnnotated(stdout, tr, gm, dsq);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);

  if (!diagfp)
    {
      if (esl_FCompareAbs(sc, asc_f, tolerance) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(sc, asc_b, tolerance) != eslOK) esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g\n", sc-asc_f, sc-asc_b);

  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(abu);
  p7_refmx_Destroy(abd);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}



/* The "multisingle" test (multiple path, single domain) samples
 * a contrived profile/seq/anchor triplet for which all paths must
 * pass through the anchor. Thus, standard and ASC F/B give the same
 * results. 
 *
 * To make this work, the contrived model is uniglocal, with a chosen
 * match state k0 that emits an anchor residue X with probability 1,
 * and the contrived sequence has a single X at the anchor position
 * (thus forcing all paths to pass through it, even without an ASC
 * algorithm).
 *
 * Then: 
 *    1. Fwd score = ASC Fwd score = Bck score = ASC Bck score.
 *    2. Viterbi score >= sampled trace score.
 *    3. Fwd/Bck scores >= Viterbi score.
 *
 *****************************************************************
 * May fail stochastically. Analysis: SRE:J13/126
 *  ./reference_fwdback_utest --diag multisingle -N 10000
 * Field 1 = fsc - bsc; same as reference_fwdback_utest
 * Field 2 = fsc - asc_f
 * Field 3 = bsc - asc_b
 * Field 4 = vsc - sampled trace score. 
 * Field 5 = fsc - vsc.                 
 * 
 * Because ASC and reference Fwd, Bck are computed in identical orders,
 * we expect 2,3 to be exactly zero with no roundoff error. We still use
 * a tolerance though, because I haven't proven this must be
 * true across all platforms. Field1 is the F-B roundoff error.
 * Field 4 is usually large, because the sampled trace is suboptimal,
 * but sometimes the trace is identical to the optimal one (in which 
 * case we get zero with no roundoff error -- or, yes it's true,
 * a trace that's "better" than Viterbi, because of numerical roundoff
 * error -- so we need to use a tolerance.
 * Field 5 should always be >=0 because this is a multipath test.
 * Also, since this is a multipath test, default logsums do have 
 * higher error than exact logsums.
 *
 * Default:
 *   Zero difference for 2,3.
 *   fsc - bsc: mean -2e-6, s.d. 0.0002 range -0.0008 .. 0.0009 =>  0.002
 *   vsc - tsc: minimum -1e-6                                   =>  1e-5
 *   fsc - vsc: minimum 0.04
 *
 * Exact:
 *   Zero difference for 2,3
 *   fsc - bsc: mean -5e-8 sd 2.4e-6 range -6e-5 .. 6e-5 => 0.001
 *   vsc - tsc: min -1e-6                                => 1e-5
 *   fsc - vsc: min 0
 * 
 */
static void
utest_multisingle(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback multisingle unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_ANCHOR  *anch      = NULL;
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       ftol      = ( p7_logsum_IsSlowExact() ? 0.001 : 0.002 );
  float       vtol      = 1e-5;
  int         status;
  
  if (bg  == NULL || rxf == NULL || rxb == NULL || afu == NULL ||
      rxv == NULL || afd == NULL || abu == NULL || abd == NULL)   esl_fatal(failmsg);

  if ((status = p7_modelsample_AnchoredUni(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceViterbi    (dsq, L, gm,          rxv, NULL, &vsc))   != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceForward    (dsq, L, gm,          rxf ,      &fsc))   != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceBackward   (dsq, L, gm,          rxb,       &bsc))   != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd,  &asc_f)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd,  &asc_b)) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Fwd:\n"); p7_refmx_Dump(stdout, rxf);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);
  //printf("FWD = BCK = ASC_FWD = ASC_BCK = %.2f\n", fsc);

  if (!diagfp) 
    {
      if (esl_FCompareAbs(fsc, bsc,   ftol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(fsc, asc_f, ftol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(bsc, asc_b, ftol) != eslOK) esl_fatal(failmsg);
      if (sc  > vsc+vtol)                             esl_fatal(failmsg);
      if (vsc > fsc)                                  esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g %20g %20g %20g\n", fsc-bsc, fsc-asc_f, bsc-asc_b, vsc-sc, fsc-vsc);
  
  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(abu);
  p7_refmx_Destroy(abd);
  p7_refmx_Destroy(rxb);
  p7_refmx_Destroy(rxf);
  p7_refmx_Destroy(rxv);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}

/* The "multipath_local" test (multiple path, local alignment allowed), like
 * "multisingle", samples a contrived profile/seq/anchorset triplet
 * for which all paths must pass through the anchor set. Thus,
 * standard and ASC F/B give the same results.
 * 
 * This version of the various multi* tests is unique in that it 
 * tests local alignment transitions.
 * 
 * The contrivance is similar for both the multilocal and multimulti
 * tests. Only M states can generate residues. Prevent any
 * N/C/J/I emission (because these can emit any residue) by setting
 * L=0 length model and all tMI = 0.  Model is unihit L=0 dual-mode:
 * local/glocal transitions are allowed, but only for one domain.
 * Choose a special residue X; only the anchor state Mk0 can
 * generate this residue with nonzero probability. The target sequence
 * has exactly 1 X residue, at the anchor.
 * 
 * Then:
 *    1. Fwd score = ASC Fwd score = Bck score = ASC Bck score.
 *    2. Viterbi score >= sampled trace score.
 *    3. Fwd/Bck scores >= Viterbi score.
 *    
 *****************************************************************
 * May fail stochastically. Analysis: SRE:J13/126
 * ./reference_fwdback_utest --diag multipath_local -N 10000
 * Field 1 = fsc - bsc; same as reference_fwdback_utest
 * Field 2 = fsc - asc_f
 * Field 3 = bsc - asc_b
 * Field 4 = vsc - sampled trace score. 
 * Field 5 = fsc - vsc.                 
 *                 
 * Same issues as the multisingle test; see notes there, above.
 * 
 * Default:
 *  Zero difference for 2,3.
 *  fsc-bsc : mean 8e-7 sd 0.0002 range -0.0006 .. 0.0007 => 0.002
 *  vsc-tsc : min -1e-6                                   => 1e-5
 *  fsc-vsc : min 5e-6
 *  
 * Exact:
 *  Zero difference for 2,3
 *  fsc-bsc : mean 3e-9 sd 4e-7 range -2e-6 .. 2e-6  => 0.0001
 *  vsc-tsc : min -1e-6                              => 1e-5
 *  fsc-vsc : min 3e-6                 
 */
static void
utest_multipath_local(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback multipath_local unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_ANCHOR  *anch      = NULL;
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       ftol      = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.002 );
  float       vtol      = 1e-5;
  int         status;
  
  if (bg  == NULL || rxf == NULL || rxb == NULL || afu == NULL ||
      rxv == NULL || afd == NULL || abu == NULL || abd == NULL)   esl_fatal(failmsg);

  if ((status = p7_modelsample_AnchoredLocal(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceViterbi    (dsq, L, gm,          rxv, NULL, &vsc))   != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceForward    (dsq, L, gm,          rxf,       &fsc))   != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceBackward   (dsq, L, gm,          rxb,       &bsc))   != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd,  &asc_f)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd,  &asc_b)) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Fwd:\n"); p7_refmx_Dump(stdout, rxf);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);

  if (!diagfp) 
    {
      if (esl_FCompareAbs(fsc, bsc,   ftol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(fsc, asc_f, ftol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(bsc, asc_b, ftol) != eslOK) esl_fatal(failmsg);
      if (sc  > vsc+vtol)                             esl_fatal(failmsg);
      if (vsc > fsc)                                  esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g %20g %20g %20g\n", fsc-bsc, fsc-asc_f, bsc-asc_b, vsc-sc, fsc-vsc);
  
  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(abu);
  p7_refmx_Destroy(abd);
  p7_refmx_Destroy(rxb);
  p7_refmx_Destroy(rxf);
  p7_refmx_Destroy(rxv);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}


/* The "multimulti" test (multiple path, multiple domain), like
 * "multisingle", samples a contrived profile/seq/anchorset triplet
 * for which all paths must pass through the anchor set. Thus,
 * standard and ASC F/B give the same results.
 * 
 * To get multiple domains, this uses a different contrivance than
 * the one used in the "multisingle" test.
 * 
 * The contrivance: Only M states can generate residues. Prevent any
 * N/C/J/I emission (because these can emit any residue) by setting
 * L=0 length model and all tMI = 0.  Model is multihit glocal;
 * multiple domains are allowed, but only in glocal mode.
 * Choose a special residue X; only the anchor state Mk0 can
 * generate this residue with nonzero probability. The target sequence
 * has exactly D X residues, one per domain.
 * 
 * Now:
 *    1. Fwd score = ASC Fwd score = Bck score = ASC Bck score.
 *    2. Viterbi score >= sampled trace score.
 *    3. Fwd/Bck scores >= Viterbi score.
 *    
 *****************************************************************
 * May fail stochastically. Analysis: SRE:J13/126.
 * ./reference_fwdback_utest --diag multimulti -N 10000
 * Field 1:  fsc-bsc   Same as reference_fwdback_utest
 * Field 2 = fsc - asc_f
 * Field 3 = bsc - asc_b
 * Field 4 = vsc - sampled trace score. 
 * Field 5 = fsc - vsc.            
 * 
 * Same numerical issues as in multisingle; see comments there.
 * 
 * Default:
 *    Zero difference for 2,3
 *    fwd-bck: mean -3e-6 sd 0.0002 range -0.001 .. 0.001  => 0.01
 *    vsc-tsc: min -1e-6                                   => 1e-5
 *    fsc-vsc: min 0
 *    
 * Exact:
 *    Zero difference for 2,3
 *    fwd-bck: mean -2e-8 sd 1e-6 range -2e-5 .. 1e-5      => 0.0001
 *    vsc-tsc: min -2e-6                                   => 1e-5
 *    fsc-vsc: min 0
 */
static void
utest_multimulti(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback multimulti unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_ANCHOR  *anch      = NULL;
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       ftol      = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  float       vtol      = 1e-5;
  int         status;
  
  if (bg  == NULL || rxf == NULL || rxb == NULL || afu == NULL ||
      rxv == NULL || afd == NULL || abu == NULL || abd == NULL)   esl_fatal(failmsg);

  if ((status = p7_modelsample_AnchoredMulti(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceViterbi    (dsq, L, gm,          rxv, NULL, &vsc))   != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceForward    (dsq, L, gm,          rxf,       &fsc))   != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceBackward   (dsq, L, gm,          rxb,       &bsc))   != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd,  &asc_f)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd,  &asc_b)) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Fwd:\n"); p7_refmx_Dump(stdout, rxf);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);

  if (!diagfp) 
    {
      if (esl_FCompareAbs(fsc, bsc,   ftol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(fsc, asc_f, ftol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(bsc, asc_b, ftol) != eslOK) esl_fatal(failmsg);
      if (sc  > vsc+vtol)                             esl_fatal(failmsg);
      if (vsc > fsc)                                  esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g %20g %20g %20g\n", fsc-bsc, fsc-asc_f, bsc-asc_b, vsc-sc, fsc-vsc);

  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(abu);
  p7_refmx_Destroy(abd);
  p7_refmx_Destroy(rxb);
  p7_refmx_Destroy(rxf);
  p7_refmx_Destroy(rxv);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}
#endif /*p7REFERENCE_ASC_FWDBACK_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/


/*****************************************************************
 * 5. Test driver.
 *****************************************************************/
#ifdef p7REFERENCE_ASC_FWDBACK_TESTDRIVE

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
  { "-N",     eslARG_INT,         "100", NULL, NULL,  NULL,  NULL, NULL, "number of times to run utest for --diag",        0 },
  { "--diag", eslARG_STRING,       NULL, NULL, NULL,  NULL,  NULL, NULL, "dump data on a utest's chance failure rate",     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for reference ASC Forward/Backward dynamic programming";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  int             M    = 10;
  int             N    = esl_opt_GetInteger(go, "-N");

  if (esl_opt_IsOn(go, "--diag"))
    {  // --diag is for studying chance failure rates, error distributions
      char *which = esl_opt_GetString(go, "--diag");

      if      (strcmp(which, "singlepath")      == 0) while (N--) utest_singlepath     (stdout, rng, M, abc);
      else if (strcmp(which, "singlesingle")    == 0) while (N--) utest_singlesingle   (stdout, rng, M, abc);
      else if (strcmp(which, "singlemulti")     == 0) while (N--) utest_singlemulti    (stdout, rng, M, abc);
      else if (strcmp(which, "multisingle")     == 0) while (N--) utest_multisingle    (stdout, rng, M, abc);
      else if (strcmp(which, "multipath_local") == 0) while (N--) utest_multipath_local(stdout, rng, M, abc);
      else if (strcmp(which, "multimulti")      == 0) while (N--) utest_multimulti     (stdout, rng, M, abc);
      else esl_fatal("--diag takes: singlepath, singlesingle, singlemulti, multisingle, multipath_local, multimulti");
    }
  else  // Running the unit tests is what we usually do:
    {
      fprintf(stderr, "## %s\n", argv[0]);
      fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

      utest_singlepath     (NULL, rng, M, abc);
      utest_singlesingle   (NULL, rng, M, abc);
      utest_singlemulti    (NULL, rng, M, abc);
      utest_multisingle    (NULL, rng, M, abc);
      utest_multipath_local(NULL, rng, M, abc);
      utest_multimulti     (NULL, rng, M, abc);

      fprintf(stderr, "#  status = ok\n");
    }

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_ASC_FWDBACK_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/



/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef p7REFERENCE_ASC_FWDBACK_EXAMPLE
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
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "stop after doing the Viterbi anchors, no manual",  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile> <ndom> [<i0> <k0>]...";
static char banner[] = "example of ASC Forward reference implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, -1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_REFMX       *vit     = p7_refmx_Create(100, 100);
  P7_REFMX       *fwd     = p7_refmx_Create(100, 100);
  P7_REFMX       *bck     = p7_refmx_Create(100, 100);
  P7_REFMX       *pp      = p7_refmx_Create(100, 100);
  P7_REFMX       *mxu     = p7_refmx_Create(100, 100);
  P7_REFMX       *mxd     = p7_refmx_Create(100, 100);
  P7_TRACE       *tr      = p7_trace_Create();
  P7_ANCHORS     *anchv   = p7_anchors_Create();
  P7_ANCHORS     *anchm   = p7_anchors_Create();
  int             D;
  int             d,i;
  float           fsc, vsc, asc, asc_b;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
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
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config   (gm, hmm, bg);
  p7_bg_SetLength     (bg, sq->n);
  p7_profile_SetLength(gm, sq->n);

  /* Read anchor coords from command line */
  D = strtol( esl_opt_GetArg(go, 3), NULL, 10);
  p7_anchors_Resize(anchm, D);
  for (i = 4, d = 1; d <= D; d++)
    {
      anchm->a[d].i0 = strtol( esl_opt_GetArg(go, i), NULL, 10); i++;
      anchm->a[d].i0 = strtol( esl_opt_GetArg(go, i), NULL, 10); i++;
    }
  anchm->D = D;
  p7_anchor_SetSentinels(anchm->a, D, sq->n, gm->M);

  p7_ReferenceViterbi (sq->dsq, sq->n, gm, vit, tr, &vsc);
  p7_ReferenceForward (sq->dsq, sq->n, gm, fwd, &fsc);   
  p7_ReferenceBackward(sq->dsq, sq->n, gm, bck, NULL);   
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, bck, pp);   

  p7_refmx_DumpBestDecoding(stdout, sq->dsq, sq->n, gm, pp);
  //p7_trace_Dump(stdout, tr);

  p7_reference_anchors_SetFromTrace(pp, tr, anchv);
  p7_ReferenceASCForward(sq->dsq, sq->n, gm, anchv->a, anchv->D, mxu, mxd, &asc);

  //p7_refmx_Dump(stdout, mxu);
  //p7_refmx_Dump(stdout, mxd);

  p7_refmx_Reuse(mxu);
  p7_refmx_Reuse(mxd);
  p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anchv->a, anchv->D, mxu, mxd, &asc_b);

  p7_refmx_Dump(stdout, mxu);
  p7_refmx_Dump(stdout, mxd);

  printf("%-20s VIT   %6.2f %6.2f %6.2f %6.2f %8.4g ", sq->name, vsc, asc, asc_b, fsc, exp(asc-fsc));
  printf("%2d ", anchv->D);
  for (d = 1; d <= anchv->D; d++) printf("%4d %4d ", anchv->a[d].i0, anchv->a[d].k0);
  printf("\n");



  if (! esl_opt_GetBoolean(go, "-v")) 
    {
      p7_refmx_Reuse(mxu);
      p7_refmx_Reuse(mxd);
      p7_ReferenceASCForward(sq->dsq, sq->n, gm, anchm->a, anchm->D, mxu, mxd, &asc);

      printf("%-20s YOURS %6s %6.2f %6.2f %8.4g ", sq->name, "", asc, fsc, exp(asc-fsc));
      printf("%2d ", anchm->D);
      for (d = 1; d <= anchm->D; d++) printf("%4d %4d ", anchm->a[d].i0, anchm->a[d].k0);
      printf("\n");
    }

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_anchors_Destroy(anchm);
  p7_anchors_Destroy(anchv);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(mxu);
  p7_refmx_Destroy(mxd);
  p7_refmx_Destroy(pp);
  p7_refmx_Destroy(bck);
  p7_refmx_Destroy(fwd);
  p7_refmx_Destroy(vit);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_ASC_FWDBACK_EXAMPLE*/




/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
