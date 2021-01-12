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
      dpc = mxu->dp[ anch->a[d-1].i0 ];      // that's dp[0] for d=1 
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
      dpc += h4R_NSCELLS;

      /* xE initialization counts exits from anchor cell.
       * Unlike the rest of the top row, MG/ML exits from the anchor cell
       * need to be calculated. Also, it has to watch out for the
       * glocal exit case when the anchor cell (unusually) sits on k=M.
       */
      xE  = (anch->a[d].k0 == M ? h4_logsum(dpc[h4R_ML], dpc[h4R_MG]) : dpc[h4R_ML]);

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
	  xc[h4R_J]  = h4_logsum( xp[h4_J] + mo->xsc[h4_J][h4_LOOP], xc[h4_E] + mo->xsc[h4_E][h4_LOOP] );
	  xc[h4R_B]  = xc[h4_J] + mo->xsc[h4_J][h4_MOVE]; 
	  xc[h4R_L]  = xc[h4_B] + mo->xsc[h4_B][0]; 
	  xc[h4R_G]  = xc[h4_B] + mo->xsc[h4_B][1]; 
	  xc[h4R_C]  = h4_logsum( xp[h4_C] + mo->xsc[h4_C][h4_LOOP], xc[h4_E] + mo->xsc[h4_E][h4_MOVE] );
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



/*****************************************************************
 * x. Test driver
 *****************************************************************/


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
