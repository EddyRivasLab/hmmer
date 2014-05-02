/* Reference implementation of anchor set constrained (ASC) Forward
 * and Backward.
 * 
 * All reference implementation code is for development and
 * testing. It is not used in HMMER's main executables. Production
 * code uses sparse dynamic programming.
 * 
 * Contents:
 *    1. ASC Forward.
 *    2. ASC Backward.
 *    3. Unit tests.
 *    4. Test driver.
 *    5. Example.
 *    6. Copyright and license information
 */


#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_coords2.h"

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
 *            The two coords in <anch>, <anch[].n1> and <anch[].n2>,
 *            are assigned to (i,k) pairs (in that order). The anchors
 *            in <anch> must be sorted in order of increasing sequence
 *            position <i>.
 *            
 *            <anch> and <D> might be data in a <P7_COORDS2> list
 *            management container: for example, for <P7_COORDS2 *dom>,
 *            you would pass <dom->arr> and <dom->n>.
 *
 * Args:      dsq    : digital target sequence 1..L
 *            L      : length of <dsq>
 *            gm     : profile
 *            anch   : array of (i,k) anchors defining <dsq>'s domain structure
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
p7_ReferenceASCForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_COORD2 *anch, int D,
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
  int   iend;			/* tmp var for row index for end of a DOWN sector        */
  int   M = gm->M;		/* for a bit more clarity, less dereference clutter      */
  int   status;

  /* Don't try to contract check the <gm> length model like we do elsewhere.
   * Envelope determination calls ASC Forward on arbitrary subsequences to
   * recalculate envelope scores.
   */

  /* reallocation, if needed */
  if ( (status = p7_refmx_GrowTo(mxu, M, L)) != eslOK) return status;
  if ( (status = p7_refmx_GrowTo(mxd, M, L)) != eslOK) return status;
#if eslDEBUGLEVEL >= 1
  p7_refmx_Zero(mxu, M, L);
  p7_refmx_Zero(mxd, M, L);
#endif
  mxu->M    = mxd->M    = M;
  mxu->L    = mxd->L    = L;
  mxu->type = p7R_ASC_FWD_UP;
  mxd->type = p7R_ASC_FWD_DOWN;

  /* Initialize i=0..anch[0].i-1 specials. 
   * All specials are stored in DOWN matrix.
   * You can think of this as a DOWN matrix for
   * a boundary condition anchor at 0,M+1: 
   * i.e. i = 0..anch[0].n1-1, k=M+1..M (nothing).
   */
  iend = (D == 0) ? 1 : anch[0].n1;
  xp   = NULL;
  for (i = 0; i < iend; i++)
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

  /* Iterate over domains d=0..D-1: */
  for (d = 0; d < D; d++)
    {
      /*****************************************************************
       * Part 1. UP matrix sector for domain d
       *    In an UP sector, we can enter the model, but not exit;
       *    so we make use of L,G state values from a previous row.
       *
       *    The UP sector includes:
       *       i = anch[d-1].i+1 to anch[d].i-1 (1..anch[d].i-1 for d==0)
       *       k = 1 to anch[d].k-1
       *
       *    We'll initialize row anch[d-1].i; then do the remaining rows.
       *    
       *    It's possible for the UP matrix to be empty (no cells), when
       *    anch[d].i == anch[d-1].i+1, including the case of anch[0].i = 1.
       *    It's also possible for the UP matrix to only consist of
       *    the initialization column k=0, when anch[d].k == 1.
       *****************************************************************/

      /* Initialization of previous row, anch[d-1].i (or, 0 for d=0) */
      i   = (d == 0 ? 0 : anch[d-1].n1);                                      // i is always our current row index, 0.1..L 
      dpc = mxu->dp[i];                                        
      for (s = 0; s < p7R_NSCELLS * anch[d].n2; s++) *dpc++ = -eslINFINITY;   // initialize previous row above UP matrix (0..anch[d].k-1) to -inf
      /* dpc is now sitting on (anch[d].i-1, anch[d].k): supercell above the anchor. 
       * We rely on that position, if UP matrix has no cells.
       */

      /* Now we recurse for remaining rows, down to anch[d].i-1 row.
       * (It's possible that no such rows exist, depending where that next anchor is.)
       */
      for (i = i+1 ; i < anch[d].n1; i++)
	{
	  rsc = gm->rsc[dsq[i]] + p7P_NR;                           // Start <rsc> at k=1 on row i 
	  tsc = gm->tsc;                                            // Start <tsc> at k=0, where off-by-one {LG}->M transition scores are
	  xp  = mxd->dp[i-1] + (M+1) * p7R_NSCELLS;                 // <xp> set to the specials of the previous row, i-1
	  dpp = mxu->dp[i-1];                                       // Start <dpp> on k=0, which we know has been init'ed to -inf above.
	  dpc = mxu->dp[i];                                         // Start <dpc> on k=0 of row i...
	  for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;  //   ... initialize that cell, and now <dpc> is on k=1.
	  dlv = dgv = -eslINFINITY;

	  for (k = 1; k < anch[d].n2; k++)
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
       *    In a DOWN sector, we can exit the model, but not enter,
       *    so we collect xE on each row,
       *    and use it to set the specials for that row.
       *    
       *    The DOWN sector includes:
       *       i = anch[d].i to anch[d+1].i-1  (anch[d].i to L, for last domain d=D-1)
       *       k = anch[d].k to M
       *       
       *    We'll initialize the top row with a partial DP calc, then
       *    do the remaining rows with the full calculation.  (With
       *    the UP matrix, we could initialize a prev row to -inf, but
       *    with the DOWN matrices, they are exactly abutting when we
       *    squeeze them down into two-matrix form.) 
       *    
       *    The top row starts with the anchor cell, which is
       *    initialized with one final UP calculation that allows
       *    entry exactly on the anchor.
       *****************************************************************/

      /* Start with anch[d].k-1 on first row, and set all cells to -inf*/
      i   = anch[d].n1;
      tsc = gm->tsc +         (anch[d].n2-1) * p7P_NTRANS;    // Start <tsc> on anch.k, i.e. k-1 relative to start of calculation
      rsc = gm->rsc[dsq[i]] + (anch[d].n2    * p7P_NR);       // <rsc> is on scores for anch.k
      xp  = mxd->dp[i-1] + (M+1) * p7R_NSCELLS;               // <xp> on specials for anch.i-1
      dpc = mxd->dp[i] + (anch[d].n2-1) * p7R_NSCELLS;
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;// <dpc> now sits on anch.k

      /* Then calculate the anchor cell (anch.i, anch.k) as an UP calc, using
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
      xE  = (anch[d].n2 == M ? p7_FLogsum(mlv, mgv) : mlv);

      /* Initialization of the rest of the top row from k=anch.k+1 to M,
       * which is only reachable on deletion paths from the anchor.
       */
      for (k = anch[d].n2+1; k <= M; k++)
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
      xc[p7R_J]  = (d == D-1 ? -eslINFINITY : xc[p7R_E] + gm->xsc[p7P_E][p7P_LOOP]);
      xc[p7R_B]  = xc[p7R_J] + gm->xsc[p7P_J][p7P_MOVE];
      xc[p7R_L]  = xc[p7R_B] + gm->xsc[p7P_B][0]; 
      xc[p7R_G]  = xc[p7R_B] + gm->xsc[p7P_B][1]; 
      xc[p7R_C]  = (d == D-1 ? xc[p7R_E] + gm->xsc[p7P_E][p7P_MOVE] : -eslINFINITY);
      xc[p7R_JJ] = -eslINFINITY;
      xc[p7R_CC] = -eslINFINITY;

      /* Now we can do the remaining rows in the Down sector of domain d. */
      iend = (d < D-1 ? anch[d+1].n1 : L+1);
      for (i = i+1 ; i < iend; i++)
	{
	  rsc = gm->rsc[dsq[i]] + anch[d].n2     * p7P_NR;         // Start <rsc> on (x_i, anchor_k, MAT) */
	  tsc = gm->tsc         + (anch[d].n2-1) * p7P_NTRANS;	   // Start <tsc> on (anchor_k-1), to pick up LMk,GMk entries 
	  dpp = mxd->dp[i-1]    + (anch[d].n2-1) * p7R_NSCELLS;	   // Start <dpp> on (i-1, anchor_k-1) 
	  dpc = mxd->dp[i]      + (anch[d].n2-1) * p7R_NSCELLS;	   // Start <dpc> on (i, anchor_k-1)... 
	  for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; //  ... and initialize the k-1 cells to -inf... 
                                                           	   //  ... so, now dpc is on anchor_k.
	  dlv = dgv = xE = -eslINFINITY;

  	  for (k = anch[d].n2; k <= M; k++) 
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
	  xc[p7R_J]  = (d == D-1 ? -eslINFINITY : p7_FLogsum( xp[p7R_J] + gm->xsc[p7P_J][p7P_LOOP], xc[p7R_E] + gm->xsc[p7P_E][p7P_LOOP]));
	  xc[p7R_B]  = xc[p7R_J] + gm->xsc[p7P_J][p7P_MOVE]; 
	  xc[p7R_L]  = xc[p7R_B] + gm->xsc[p7P_B][0]; 
	  xc[p7R_G]  = xc[p7R_B] + gm->xsc[p7P_B][1]; 
	  xc[p7R_C]  = (d == D-1 ? p7_FLogsum( xp[p7R_C] + gm->xsc[p7P_C][p7P_LOOP], xc[p7R_E] + gm->xsc[p7P_E][p7P_MOVE]) : -eslINFINITY);
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
 *            The two coords in <anch>, <anch[].n1> and <anch[].n2>,
 *            are assigned to (i,k) pairs (in that order). The anchors
 *            in <anch> must be sorted in order of increasing sequence
 *            position <i>.
 *            
 *            <anch> and <D> might be data in a <P7_COORDS2> list
 *            management container: for example, for <P7_COORDS2 *dom>,
 *            you would pass <dom->arr> and <dom->n>.
 *
 * Args:      dsq    : digital target sequence 1..L
 *            L      : length of <dsq>
 *            gm     : profile
 *            anch   : array of (i,k) anchors defining <dsq>'s domain structure
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
p7_ReferenceASCBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_COORD2 *anch, int D,
			P7_REFMX *abu, P7_REFMX *abd, float *opt_sc)
{
  const float *tsc;		/* ptr into transition scores of <gm> */
  const float *rsc;		/* ptr into emission scores of <gm> for residue dsq[i] on current row i  */
  const float *rsn;		/* ptr into emission scores of <gm> for residue dsq[i+1] on next row i+1 */
  float *dpc, *dpn;		/* ptrs into DP matrix for current row i, next row i+1  */
  float *xc;			/* ptr to specials on current row; specials are stored in DOWN, <abd> */
  int    d;                   	/* counter over domains 0..D-1 */
  int    i;			/* counter over sequence positions 0.1..L (DP rows) */
  int    k;			/* counter over model positions 0.1..M (DP columns) */
  int    iend;
  float  mgc, mlc;
  float  mgn, mln;
  float  dgn, dln;
  float  ign, iln;
  float  xE;
  float  xG,  xL;
  float  xC, xJ, xN;
  int    M = gm->M;
  int    status;

  /* contract checks / arg validation */
  ESL_DASSERT1( ( gm->L == L || gm->L == 0) ); /* length model in profile is either L (usually) or 0 (some unit tests) */

  /* reallocation, if needed */
  if ( (status = p7_refmx_GrowTo(abu, gm->M, L)) != eslOK) return status;
  if ( (status = p7_refmx_GrowTo(abd, gm->M, L)) != eslOK) return status;
#if eslDEBUGLEVEL >= 1
  p7_refmx_Zero(abu, M, L);
  p7_refmx_Zero(abd, M, L);
#endif
  abu->M    = abd->M    = M;
  abu->L    = abd->L    = L;
  abu->type = p7R_ASC_BCK_UP;
  abd->type = p7R_ASC_BCK_DOWN;

  iend = (D == 0 ? 0 : anch[D-1].n1);
  for (i = L; i >= iend; i--)
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
      xc[p7R_E]  = xC + gm->xsc[p7P_E][p7P_MOVE];
    }

  /* The code below is designed to be easily convertible to one-row memory efficient DP, if needed */
  for (d = D-1; d >= 0; d--)
    {
      /* DOWN matrix.
       *   i = anch[d].i .. anch[d+1].i-1
       *   k = anch[d].k .. M
       *   calculated Backward. 
       * In the DOWN matrix, paths can end from the model, but not start in it,
       * so we evaluate {MD}->E transitions backward, but we don't evaluate 
       * B->{LG}->Mk
       */
      iend = (d == D-1 ? L : anch[d+1].n1-1);
      for (i = iend; i >= anch[d].n1; i--)
	{
	  rsn  = (i == iend ? NULL : gm->rsc[dsq[i+1]] + M * p7P_NR);   // residue scores on next row; start at M
	  tsc  = gm->tsc + M * p7P_NTRANS;                              // transition scores: start at M
	  dpc  = abd->dp[i]   + M * p7R_NSCELLS;                        // current row of DP matrix: start at M
	  xE   = dpc[p7R_NSCELLS+p7R_E];                                // pick up the xE score; specials start at M+1, hence the p7R_NSCELLS bump here
	  dpn  = (i == iend ? NULL : abd->dp[i+1] + M * p7R_NSCELLS);   // next row of DP matrix: start at M

	  mgn = dgn = -eslINFINITY;
	  mln = dln = -eslINFINITY;
	  ign = iln = -eslINFINITY;
	  xG  = xL  = -eslINFINITY;

	  for (k = M; k >= anch[d].n2; k--)
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
      rsc = gm->rsc[ dsq[anch[d].n1]] + anch[d].n2 * p7P_NR;
      mgn = mgc + *rsc;  
      mln = mlc + *rsc;
      xG  = mgn + tsc[p7P_GM];
      xL  = mln + tsc[p7P_LM];

      xJ = xN = -eslINFINITY;

      /* UP matrix */
      iend = (d == 0 ? 1 : anch[d-1].n1+1);
      for (i = anch[d].n1-1; i >= iend; i--)
	{
	  xc = abd->dp[i] + (M+1) * p7R_NSCELLS;   // on specials, which are in DOWN matrix
	  xc[p7R_CC] = -eslINFINITY;  // CC,JJ are only used in decoding matrices
	  xc[p7R_JJ] = -eslINFINITY;
	  xc[p7R_C]  = -eslINFINITY;  // C is now unreachable, when anchor set constrained.
	  xc[p7R_G]  = xG;            // xG was accumulated during prev row; G->Mk wing unfolded
	  xc[p7R_L]  = xL;            // xL accumulated on prev row
	  xc[p7R_B]  = p7_FLogsum(xG + gm->xsc[p7P_B][1],  xL + gm->xsc[p7P_B][0]); 
	  xc[p7R_J]  = xJ = (d == 0 ? -eslINFINITY : p7_FLogsum(xJ + gm->xsc[p7P_J][p7P_LOOP], xc[p7R_B] + gm->xsc[p7P_J][p7P_MOVE]));
	  xc[p7R_N]  = xN = (d  > 0 ? -eslINFINITY : p7_FLogsum(xN + gm->xsc[p7P_N][p7P_LOOP], xc[p7R_B] + gm->xsc[p7P_N][p7P_MOVE]));
	  xc[p7R_E]  = xc[p7R_J] + gm->xsc[p7P_E][p7P_LOOP];  

	  tsc = gm->tsc    + (anch[d].n2-1) * p7P_NTRANS;                                    // transition scores: start at anch[d].k-1
	  dpc = abu->dp[i] + (anch[d].n2-1) * p7R_NSCELLS;                                   // on anch[d].k-1
	  dpn = (i == anch[d].n1-1 ? NULL : abu->dp[i+1] + (anch[d].n2 - 1) * p7R_NSCELLS);  // on anch[d].k-1
	  rsc = gm->rsc[dsq[i]] + (anch[d].n2-1) * p7P_NR;
	  rsn = (i == anch[d].n1-1 ? NULL : gm->rsc[dsq[i+1]] + (anch[d].n2-1) * p7P_NR);

	  xG  = xL  = -eslINFINITY; 
	  dgn = dln = -eslINFINITY;
	  ign = iln = -eslINFINITY;
	  if (i < anch[d].n1-1) mgn = mln = -eslINFINITY; /* allow mgn/mln to carry over from anchor cell */

	  /* The recursion is the same as for the DOWN matrix, so only differences are commented on: */
	  for (k = anch[d].n2-1; k >= 1; k--)
	    {
	      if (i < anch[d].n1-1) {           
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
	      
	      if (i < anch[d].n1-1) {       
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

      xc = abd->dp[i] + (M+1) * p7R_NSCELLS;   // on specials, which are in DOWN matrix
      xc[p7R_CC] = -eslINFINITY;  // CC,JJ are only used in decoding matrices
      xc[p7R_JJ] = -eslINFINITY;
      xc[p7R_C]  = -eslINFINITY;  // C is now unreachable, when anchor set constrained.
      xc[p7R_G]  = xG;            // xG was accumulated during prev row; G->Mk wing unfolded
      xc[p7R_L]  = xL;            // xL accumulated on prev row
      xc[p7R_B]  = p7_FLogsum(xG + gm->xsc[p7P_B][1],  xL + gm->xsc[p7P_B][0]); 
      xc[p7R_J]  = xJ = (d == 0 ? -eslINFINITY : p7_FLogsum(xJ + gm->xsc[p7P_J][p7P_LOOP], xc[p7R_B] + gm->xsc[p7P_J][p7P_MOVE]));
      xc[p7R_N]  = xN = (d  > 0 ? -eslINFINITY : p7_FLogsum(xN + gm->xsc[p7P_N][p7P_LOOP], xc[p7R_B] + gm->xsc[p7P_N][p7P_MOVE]));
      xc[p7R_E]  = xc[p7R_J] + gm->xsc[p7P_E][p7P_LOOP];  

    } /* end loop over domains d */

  if (opt_sc) *opt_sc = xN;
  return eslOK;
}

/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7REFERENCE_ASC_FWDBACK_TESTDRIVE
#include "hmmer.h"



/* Compare a standard DP matrix <std> to ASC UP and DOWN matrices
 * <ascu> and <ascd>, using anchor set <anch> for <D> domains.
 * Compare all valid values in the ASC matrices to their counterparts
 * in the standard matrix for equality within absolute (not relative)
 * tolerance <epsilon>. (Ignore cells that are valid in the standard
 * matrix, but not in the ASC matrices.) Return <eslOK> if all
 * cell comparisons succeed; <eslFAIL> if not. 
 * 
 * Similar to p7_refmx_Compare(). See notes there for why we prefer
 * absolute, not relative epsilon (short version: error accumulates
 * more as an absolute # per term added; DP cells with values close to
 * zero may have large relative errors).
 */
static int
ascmatrix_compare(P7_REFMX *std, P7_REFMX *ascu, P7_REFMX *ascd, P7_COORD2 *anch, int D, float epsilon)
{
  int M         = std->M;
  int L         = std->L;
  int killmenow = FALSE;
  int d,i,k,s,iz;
#ifdef p7_DEBUGGING
  killmenow = TRUE;
#endif

  ESL_DASSERT1( (ascu->M == M && ascu->L == L));
  ESL_DASSERT1( (ascd->M == M && ascd->L == L));
  ESL_DASSERT1( ((std->type == p7R_FORWARD  && ascu->type == p7R_ASC_FWD_UP && ascd->type == p7R_ASC_FWD_DOWN) ||
		 (std->type == p7R_BACKWARD && ascu->type == p7R_ASC_BCK_UP && ascd->type == p7R_ASC_BCK_DOWN)));


  /* i=0..anch[0].i-1 specials are a boundary case... */
  iz = (D == 0 ? 1 : anch[0].n1);
  for (i = 0; i < iz; i++)
    {
      /* In Backwards, don't look at J or L cells. */
      if (std->type == p7R_FORWARD) 
	{
	  for (s = 0; s < p7R_NXCELLS; s++)
	    if (P7R_XMX(ascd,i,s) != -eslINFINITY && esl_FCompareAbs( P7R_XMX(std,i,s), P7R_XMX(ascd,i,s), epsilon) == eslFAIL) 
	      { if (killmenow) abort(); return eslFAIL; }
	}
      else
	{
	  if (esl_FCompareAbs( P7R_XMX(std,i,p7R_N), P7R_XMX(ascd,i,p7R_N), epsilon) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	  if (esl_FCompareAbs( P7R_XMX(std,i,p7R_B), P7R_XMX(ascd,i,p7R_B), epsilon) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	  if (esl_FCompareAbs( P7R_XMX(std,i,p7R_G), P7R_XMX(ascd,i,p7R_G), epsilon) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	}
    }

  for (d = 0; d < D; d++)
    {
      /* UP sector for domain d */
      iz = (d == 0 ? 1  : anch[d-1].n1+1);
      for (i = iz; i < anch[d].n1; i++)
	{
	  for (k = 1; k < anch[d].n2; k++)
	    {
	      for (s = 0; s < p7R_NSCELLS; s++)
		if (P7R_MX(ascu,i,k,s) != -eslINFINITY && esl_FCompareAbs( P7R_MX(std,i,k,s), P7R_MX(ascu,i,k,s), epsilon) == eslFAIL)
		  { if (killmenow) abort(); return eslFAIL; }
	    }
	}

      /* DOWN sector for domain d */
      iz = ( d < D-1 ? anch[d+1].n1 : L+1);
      for (i = anch[d].n1; i < iz; i++)
	{
	  for (k = anch[d].n2; k <= M; k++)
	    {
	      for (s = 0; s < p7R_NSCELLS; s++)
		if (P7R_MX(ascd,i,k,s) != -eslINFINITY && esl_FCompareAbs( P7R_MX(std,i,k,s), P7R_MX(ascd,i,k,s), epsilon) == eslFAIL) 
		  { if (killmenow) abort(); return eslFAIL; }
	    }

	  /* In specials, ASC has some known discrepancies from standard implementation.
	   * Standard implementation fills in C and J identically because it doesn't know where last domain is,
	   * and continues into J->B->{LG}. In ASC, J/B/L/G are -inf after final domain d. 
	   */
	  if (d < D-1) 
	    {
	      for (s = 0; s < p7R_NXCELLS; s++)
		if (P7R_XMX(ascd,i,s) != -eslINFINITY && s != p7R_C && esl_FCompareAbs( P7R_XMX(std,i,s), P7R_XMX(ascd,i,s), epsilon) == eslFAIL) 
		  { if (killmenow) abort(); return eslFAIL; }
	    }
	  else 
	    {
	      if (esl_FCompareAbs( P7R_XMX(std,i,p7R_E), P7R_XMX(ascd,i,p7R_E), epsilon) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	      if (esl_FCompareAbs( P7R_XMX(std,i,p7R_C), P7R_XMX(ascd,i,p7R_C), epsilon) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	    }
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
 */
static void
utest_singlepath(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback singlepath unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = p7_profile_Create(M, bg->abc);
  ESL_SQ     *sq        = esl_sq_CreateDigital(bg->abc);
  P7_TRACE   *tr        = p7_trace_Create();
  P7_TRACE   *vtr       = p7_trace_Create();
  P7_COORD2  *anch      = malloc(sizeof(P7_COORD2));
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
  float       epsilon   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  int         status;
  
  if (gm  == NULL || sq  == NULL || tr  == NULL ||
      bg  == NULL || vtr == NULL || rxv == NULL ||
      rxf == NULL || rxb == NULL || afu == NULL ||
      afd == NULL || abu == NULL || abd == NULL || anch == NULL)   esl_fatal(failmsg);

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
	anch[0].n1 = tr->i[z];
	anch[0].n2 = tr->k[z];
	break;
      }
    }
  }

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

  if (esl_FCompare(sc, vsc,   epsilon) != eslOK) esl_fatal(failmsg);  // generated trace score = Viterbi score
  if (esl_FCompare(sc, fsc,   epsilon) != eslOK) esl_fatal(failmsg);  //  ... = Forward score
  if (esl_FCompare(sc, bsc,   epsilon) != eslOK) esl_fatal(failmsg);  //  ... = Backward score
  if (esl_FCompare(sc, asc_f, epsilon) != eslOK) esl_fatal(failmsg);  //  ... = ASC Forward score
  if (esl_FCompare(sc, asc_b, epsilon) != eslOK) esl_fatal(failmsg);  //  ... = ASC Backward score
  if (p7_trace_Compare(tr, vtr, 0) != eslOK) esl_fatal(failmsg);      // generated trace = Viterbi trace

  /* to compare Viterbi to Fwd matrix, we have to hack around a safety check in the structures */
  rxv->type = p7R_FORWARD;
  if (p7_refmx_Compare(rxf, rxv, epsilon) != eslOK) esl_fatal(failmsg);
  rxv->type = p7R_VITERBI;

  if (ascmatrix_compare(rxf, afu, afd, anch, D, epsilon) != eslOK) esl_fatal(failmsg);
  if (ascmatrix_compare(rxb, abu, abd, anch, D, epsilon) != eslOK) esl_fatal(failmsg);
  
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
 */
static void
utest_singlesingle(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback singlesingle unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_TRACE   *vtr       = p7_trace_Create();
  P7_COORD2  *anch      = NULL;
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       epsilon   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
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

  if (esl_FCompare(sc, vsc, epsilon)   != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(sc, fsc, epsilon)   != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(sc, bsc, epsilon)   != eslOK) esl_fatal(failmsg);
  if (p7_trace_Compare(tr, vtr, 0)     != eslOK) esl_fatal(failmsg);

  if (esl_FCompare(sc, asc_f, epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(sc, asc_b, epsilon) != eslOK) esl_fatal(failmsg);

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
 */
static void
utest_singlemulti(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback singlemulti unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_COORD2  *anch      = NULL;
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, asc_f, asc_b;
  float       epsilon   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  int         status;
  
  if (bg  == NULL || afu == NULL || afd == NULL || abu == NULL || abd == NULL) 
    esl_fatal(failmsg);

  if ((status = p7_modelsample_SinglePathedASC(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd, &asc_f)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd, &asc_b)) != eslOK) esl_fatal(failmsg);

  //p7_trace_DumpAnnotated(stdout, tr, gm, dsq);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);

  if (esl_FCompare(sc, asc_f, epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(sc, asc_b, epsilon) != eslOK) esl_fatal(failmsg);

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
 */
static void
utest_multisingle(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback multisingle unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_COORD2  *anch      = NULL;
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       epsilon   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
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

  if (esl_FCompare(fsc, bsc,   epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(fsc, asc_f, epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(bsc, asc_b, epsilon) != eslOK) esl_fatal(failmsg);
  if (sc  > vsc+epsilon)                          esl_fatal(failmsg);
  if (vsc > fsc+epsilon)                          esl_fatal(failmsg);
  
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
 */
static void
utest_multipath_local(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback multipath_local unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_COORD2  *anch      = NULL;
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       epsilon   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
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

  if (esl_FCompare(fsc, bsc,   epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(fsc, asc_f, epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(bsc, asc_b, epsilon) != eslOK) esl_fatal(failmsg);
  if (sc  > vsc+epsilon)                          esl_fatal(failmsg);
  if (vsc > fsc+epsilon)                          esl_fatal(failmsg);
  
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
 */
static void
utest_multimulti(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_fwdback multimulti unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_COORD2  *anch      = NULL;
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       epsilon   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
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

  if (esl_FCompare(fsc, bsc,   epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(fsc, asc_f, epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(bsc, asc_b, epsilon) != eslOK) esl_fatal(failmsg);
  if (sc  > vsc+epsilon)                          esl_fatal(failmsg);
  if (vsc > fsc+epsilon)                          esl_fatal(failmsg);

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
 * 4. Test driver.
 *****************************************************************/
#ifdef p7REFERENCE_ASC_FWDBACK_TESTDRIVE

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
static char banner[] = "unit test driver for reference ASC Forward/Backward dynamic programming";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  int             M    = 10;

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_singlepath     (rng, M, abc);
  utest_singlesingle   (rng, M, abc);
  utest_singlemulti    (rng, M, abc);
  utest_multisingle    (rng, M, abc);
  utest_multipath_local(rng, M, abc);
  utest_multimulti     (rng, M, abc);

  fprintf(stderr, "#  status = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_ASC_FWDBACK_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/



/*****************************************************************
 * 5. Example
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
  P7_COORDS2     *anchv   = p7_coords2_Create(0,0);
  P7_COORDS2     *anchm   = p7_coords2_Create(0,0);
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
  p7_coords2_GrowTo(anchm, D);
  for (i = 4, d = 0; d < D; d++)
    {
      anchm->arr[d].n1 = strtol( esl_opt_GetArg(go, i), NULL, 10); i++;
      anchm->arr[d].n2 = strtol( esl_opt_GetArg(go, i), NULL, 10); i++;
    }
  anchm->n    = D;

  p7_ReferenceViterbi (sq->dsq, sq->n, gm, vit, tr, &vsc);
  p7_ReferenceForward (sq->dsq, sq->n, gm, fwd, &fsc);   
  p7_ReferenceBackward(sq->dsq, sq->n, gm, bck, NULL);   
  p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, bck, pp);   

  p7_refmx_DumpBestDecoding(stdout, sq->dsq, sq->n, gm, pp);
  //p7_trace_Dump(stdout, tr);

  p7_reference_anchors_SetFromTrace(pp, tr, anchv);
  p7_ReferenceASCForward(sq->dsq, sq->n, gm, anchv->arr, anchv->n, mxu, mxd, &asc);

  //p7_refmx_Dump(stdout, mxu);
  //p7_refmx_Dump(stdout, mxd);

  p7_refmx_Reuse(mxu);
  p7_refmx_Reuse(mxd);
  p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anchv->arr, anchv->n, mxu, mxd, &asc_b);

  p7_refmx_Dump(stdout, mxu);
  p7_refmx_Dump(stdout, mxd);

  printf("%-20s VIT   %6.2f %6.2f %6.2f %6.2f %8.4g ", sq->name, vsc, asc, asc_b, fsc, exp(asc-fsc));
  printf("%2d ", anchv->n);
  for (d = 0; d < anchv->n; d++) printf("%4d %4d ", anchv->arr[d].n1, anchv->arr[d].n2);
  printf("\n");



  if (! esl_opt_GetBoolean(go, "-v")) 
    {
      p7_refmx_Reuse(mxu);
      p7_refmx_Reuse(mxd);
      p7_ReferenceASCForward(sq->dsq, sq->n, gm, anchm->arr, anchm->n, mxu, mxd, &asc);

      printf("%-20s YOURS %6s %6.2f %6.2f %8.4g ", sq->name, "", asc, fsc, exp(asc-fsc));
      printf("%2d ", anchm->n);
      for (d = 0; d < anchm->n; d++) printf("%4d %4d ", anchm->arr[d].n1, anchm->arr[d].n2);
      printf("\n");
    }

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_coords2_Destroy(anchm);
  p7_coords2_Destroy(anchv);
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
