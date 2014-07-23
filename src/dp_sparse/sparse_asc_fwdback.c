/* Production implementation of anchor set constrained (ASC) Forward
 * and Backward, using sparse dynamic programming.
 * 
 * Contents:
 *    1. Sparse ASC Forward
 *    2. Sparse ASC Backward
 *    3. Footnotes
 *    4. Unit tests
 *    5. Test driver
 *    6. Example
 *    7. Copyright and license information
 */

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_coords2.h"

#include "misc/logsum.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spascmx.h"


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
  float const *tsc = gm->tsc;	 // sets up TSC() macro, access to profile's transitions      
  float const *rsc = NULL;	 // will be set up for MSC(), ISC() macros for residue scores for each row i
  int    g         = 1;	         // index of next or current segment. When we enter it, and while we're in it, in_seg = TRUE. s=1..S, with sentinels.
  int    in_seg    = FALSE;      //   ... this bumps to TRUE when we see ia(g), first row of segment; to FALSE on ib(g), last row.
  int    d         = 1;          // idx of next domain anchor we reach. Row i may be in sector UP(d), DOWN(d-1).
  int    ndown     = 0;          //   ... this bumps to 1 when i reaches an anchor; row # of DOWN sector; back to 0 at seg end
  float *dpc       = NULL;       // all indexing is relative! This ptr walks through the sparse ASC matrix main cells.
  float *xc        = NULL;       //   ... and this one, through the specials.
  float *dpp       = NULL;       // ptr for stepping thru DP cells on a previous row. Gets initialized by...
  float *last_down = NULL;       //   ... remembering the <dpc> start of each DOWN...
  float *last_up   = NULL;       //   ... or UP sector as we compute them.
  float  mlc, mgc;               // tmp variables for computing MLk, MGk score
  float  dlc, dgc;               //   ... and DLk, DGk.
  float  xE,xN,xJ,xB,xL,xG,xC;   // current special state values on this row; specials are not always stored
  int i;                         // index for dsq[i] seq position/rows 0..L
  int k;                         // index of model positions 1..M
  int z;                         // index for sparse cell indices k[i][z] along a row
  int y;			 //   ... and for previous row, k[i-1][y].
  int s;                         // index of individual states
  int ngap;                      // number of residues in an intersegment interval, all must be accounted for by N, J, or C, depending on what anchor we're on
  int status;

  /* Contract checks, argument validation */
  ESL_DASSERT1(( sm->L == L     ));
  ESL_DASSERT1(( sm->M == gm->M ));

  /* Assure that <asf> is allocated large enough. 
   */
  if (( status = p7_spascmx_Reinit(asf, sm, anch, D)) != eslOK) return status;
  asf->type = p7S_ASC_FWD;
  
  dpc = asf->dp;   // you can't initialize <dpc> until AFTER the _Reinit() call above, because Reinit might move the memory.
  xc  = asf->xmx;  //   ... ditto for <xc>, as well as other tmp ptrs into the DP matrix memory.
  xN  = 0.;
  for (i = 0; i <= L; i++)
    {
      /* Mechanics of traversing a sparse ASC matrix, row by row */
      if      (i == anch[d].i0) { ndown = 1; d++; }            // when i reaches next anchor; bump d to next domain index, and DOWN sector is active...
      else if (sm->n[i] == 0)   { ndown = 0;      }            //  ... until when we reach end of segment, when DOWN becomes inactive again.
      else if (ndown)           { ndown++;        }            // counting ndown lets us determine when we're in 1st, 2nd row; useful for boundary conditions on UP.

      if      (i >  sm->seg[g].ib)  { in_seg  = FALSE; g++; }  // g bumps, to start expecting to see start of segment <g> next.
      else if (i == sm->seg[g].ia)  { in_seg  = TRUE;       }  // g might be S+1, but this is safe because of sentinel seg[S+1].ia=ib=L+2      

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
	  *dpc++ = -eslINFINITY;       // IL(i0,k0)
	  *dpc++ = -eslINFINITY;       // IG(i0,k0)
	  xE     = p7_FLogsum(xE, mlc);               // dlc is -inf, not included in sum
	  *dpc++ = -eslINFINITY;       // DL(i0,k0)
	  *dpc++ = -eslINFINITY;       // DG(i0,k0)
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

	  while (z < sm->n[i] && sm->k[i][z] < anch[d-1].k0) z++;   // skip sparse cells that aren't in DOWN sector (which is k0..M) 
	  for (; z < sm->n[i]; z++)                                 // then calculate the rest of the row.
	    {                                  
	      k = sm->k[i][z];  // for notational convenience

	      /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	      mlc = mgc = -eslINFINITY;
	      while (y < sm->n[i-1] && sm->k[i-1][y] < anch[d-1].k0) y++;                           // skip cells on prev row that aren't in DOWN at all
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

      if (in_seg && anch[d].i0 <= sm->seg[g].ib)                             // d may be D+1 here: if so, sentinel makes the the comparison to seg[g].ib fail
	{
	  if (ndown == 1)       // The first row of an UP matrix for subsequent domains in a segment is the anchor row i0.
	    {                   // All sparse cells in this row are unreachable, initialized to -inf. They only get used for G->DDDD->Mk,i+1 decoding, wing unfolding.
	      last_up = dpc;
	      for (z = 0; z < sm->n[i] && sm->k[i][z] < anch[d].k0; z++)
		for (s = 0; s < p7S_NSCELLS; s++)
		  *dpc++ = -eslINFINITY;
	    }

	  else if (i == sm->seg[g].ia)         // The first row of UP(d) when d is the first domain in a segment
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
       * *calculation* of specials (by SparseForward()) rules from
       * *storage* of specials (by ASC rules).  For full explanation
       * of why this is so, see footnote [1].
       */
      if (i == sm->seg[g].ia-1)
	{
	  ngap =  (g == 1 ? i : sm->seg[g].ia - sm->seg[g-1].ib - 1);

	  if (d == 1)                  xN = xN + (ngap == 0 ? 0. : (float) ngap * gm->xsc[p7P_N][p7P_LOOP]);
	  else                         xN = -eslINFINITY;

	  if      (d == 1 || d == D+1) xJ = -eslINFINITY;
	  else                         xJ = p7_FLogsum(xJ + (float) ngap * gm->xsc[p7P_J][p7P_LOOP], 
						       xE +                gm->xsc[p7P_E][p7P_LOOP]);

	  if      (d == D+1)           xC = p7_FLogsum(xC + (float) ngap * gm->xsc[p7P_C][p7P_LOOP], 
						       xE +                gm->xsc[p7P_E][p7P_MOVE]);
	  else                         xC = -eslINFINITY;

	  xB = (d == 1 ? xN + gm->xsc[p7P_N][p7P_MOVE] : xJ + gm->xsc[p7P_J][p7P_MOVE]);
	  xL = xB + gm->xsc[p7P_B][0]; /* B->L */
	  xG = xB + gm->xsc[p7P_B][1]; /* B->G */
	}
      else if (i >= sm->seg[g].ia && i <= sm->seg[g].ib)
	{
	  xN = (d == 1 ? xN + gm->xsc[p7P_N][p7P_LOOP] : -eslINFINITY);

	  if      (d == 1 || d == D+1) xJ = -eslINFINITY;
	  else if (ndown == 1)         xJ = xE + gm->xsc[p7P_E][p7P_LOOP];  // don't propagate J->J across DOWN sector boundaries
	  else                         xJ = p7_FLogsum(xJ + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_LOOP]);

	  xB = (d == 1 ? xN + gm->xsc[p7P_N][p7P_MOVE] : xJ + gm->xsc[p7P_J][p7P_MOVE]);
	  xL = xB + gm->xsc[p7P_B][0]; /* B->L */
	  xG = xB + gm->xsc[p7P_B][1]; /* B->G */
	  xC = (d == D+1 ? p7_FLogsum(xE + gm->xsc[p7P_E][p7P_MOVE],  xC + gm->xsc[p7P_C][p7P_LOOP]) : -eslINFINITY);
	}
      
      /* ASC only stores specials on ia(g)-1..ib(g) for segments g that contain >= 1 anchor. simple to say, a little awkward to test for: */
      if ( (i >= sm->seg[g].ia-1 && anch[d].i0 <= sm->seg[g].ib) || ndown)   // d is next anchor we'll see. i0(d) <= ib until we're in last DOWN sector in seg.
	{                                                                    // in last DOWN sector in seg, i0(d) is in a later seg... hence we need the <ndown> test too here
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
    } // end of i=0..L. Sparse ASC DP matrix complete. Mission accomplished.

  ngap = L - sm->seg[sm->S].ib;
  xC   = (ngap == 0 ? xC : xC + (float) ngap * gm->xsc[p7P_C][p7P_LOOP]);
  if (opt_sc) *opt_sc = xC + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}
/*--------------- end, sparse ASC Forward -----------------------*/



/*****************************************************************
 * 2. Sparse ASC Backward
 *****************************************************************/

int
p7_sparse_asc_Backward(void)
{

  return eslOK;
}


/*-------------- end, sparse ASC Backward -----------------------*/


/*****************************************************************
 * 3. Footnotes
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
 * 3. Unit tests
 *****************************************************************/
#ifdef p7SPARSE_ASC_FWDBACK_TESTDRIVE
#include "hmmer.h"

/* "compare_reference" unit test.
 * 
 * When we include all i,k supercells in the sparse mask, then sparse
 * DP calculations give the same results as reference calculations,
 * even at the level of individual DP cells.
 * 
 * Sample a random profile of length <M>.  Generate <N> sequences from
 * that profile, using a length model of <L> during sampling.  For
 * each sampled sequence, Make a sparse mask that contains all i,k
 * cells. Make an anchor set from the generating trace (the anchor set
 * just needs to be reasonable, not optimal).
 * 
 * Then:
 *   1. Reference and sparse Forward scores must be identical (within 
 *      numerical tolerance).
 *   2. Cells of reference and sparse DP matrix have identical values
 *      (within numerical tolerance).
 *   3. Sparse Fwd matrix structure passes Validate().
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
  P7_SPARSEMX   *asf   = p7_sparsemx_Create(NULL);
  float          sc1, sc2;
  int            idx;
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

      /* reference ASC forward calculation */
      if ( p7_ReferenceASCForward(sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, &sc1) != eslOK) esl_fatal(msg);
      
      /* sparse ASC forward calculation */
      if ( p7_sparse_asc_Forward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asf, &sc2)   != eslOK) esl_fatal(msg);

      //p7_refmx_Dump(stdout, afu);
      //p7_refmx_Dump(stdout, afd);
      //p7_spascmx_Dump(stdout, asf, anch->a, anch->D);

      /* comparisons */
      if (! diagfp) 
	{
	  if ( p7_spascmx_CompareReference(asf, anch->a, anch->D, afu, afd, tol) != eslOK) esl_fatal(msg); // test this first; makes debugging easier, if there's a bad difference, make it fail on the bad cell in the DP matrix
	  if ( esl_FCompareAbs(sc1, sc2, tol) != eslOK) esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%20g\n", sc1-sc2);

      p7_sparsemask_Reuse(sm);
      p7_anchors_Reuse(anch);
      p7_trace_Reuse(gtr);
      p7_refmx_Reuse(afu);
      p7_refmx_Reuse(afd);
      p7_sparsemx_Reuse(asf);
      esl_sq_Reuse(sq);
    }

  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(afd);
  p7_anchors_Destroy(anch);
  p7_trace_Destroy(gtr);
  p7_sparsemask_Destroy(sm);
  p7_sparsemx_Destroy(asf);
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
  char           failmsg[] = "sparse_asc_fwdback :: singlesingle unit test failed";
  P7_BG         *bg        = p7_bg_Create(abc);
  P7_HMM        *hmm       = NULL;
  P7_PROFILE    *gm        = NULL;
  ESL_DSQ       *dsq       = NULL;
  P7_TRACE      *tr        = NULL;
  P7_ANCHOR     *anch      = NULL;
  P7_SPARSEMASK *sm        = NULL;
  P7_SPARSEMX   *asf       = p7_sparsemx_Create(NULL);
  int   D,L;
  int   idx;
  float tsc, fsc;
  float tol = 0.0001;

  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_SinglePathedSeq(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &tsc) != eslOK) esl_fatal(failmsg);

      if ((sm = p7_sparsemask_Create(M, L))            == NULL) esl_fatal(failmsg);
      if ( p7_sparsemask_SetFromTrace(sm, rng, tr)    != eslOK) esl_fatal(failmsg);

      if ( p7_sparse_asc_Forward(dsq, L, gm, anch, D, sm, asf, &fsc) != eslOK) esl_fatal(failmsg);

      //p7_trace_DumpAnnotated(stdout, tr, gm, dsq);
      //p7_spascmx_Dump(stdout, asf, anch, D);
      
      if (!diagfp)
	{
	  if (esl_FCompareAbs(tsc, fsc, tol) != eslOK) esl_fatal(failmsg);
	}
      else
	fprintf(diagfp, "%20g\n", tsc-fsc);

      p7_sparsemx_Reuse(asf);

      free(anch);
      p7_trace_Destroy(tr);
      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
      p7_sparsemask_Destroy(sm);
      free(dsq);
    }

  p7_sparsemx_Destroy(asf);
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
  P7_SPARSEMX   *asf = p7_sparsemx_Create(NULL);
  int            idx;
  float          vsc, fsc, asc_f;
  float          tol = 0.0001;

  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_AnchoredUni(rng, M, bg, &hmm, &gm, &dsq, &L, &gtr, &anch, &D, &gsc) != eslOK) esl_fatal(msg);

      if ((sm = p7_sparsemask_Create(gm->M, L)) == NULL) esl_fatal(msg);
      
      if (idx%2) {
	if ((om = p7_oprofile_Create(hmm->M, abc)) == NULL)  esl_fatal(msg);
	if ( p7_oprofile_Convert(gm, om)           != eslOK) esl_fatal(msg);
	if ( p7_oprofile_ReconfigMultihit(om, L)   != eslOK) esl_fatal(msg);
	om->mode = p7_LOCAL;

	if ((cx = p7_checkptmx_Create(hmm->M, L, ESL_MBYTES(32))) == NULL) esl_fatal(msg);
	if ( p7_ForwardFilter (dsq, L, om, cx, /*fsc=*/NULL) != eslOK) esl_fatal(msg);
	if ( p7_BackwardFilter(dsq, L, om, cx, sm, p7_SPARSEMASK_THRESH_DEFAULT) != eslOK) esl_fatal(msg);

	p7_oprofile_Destroy(om);
	p7_checkptmx_Destroy(cx);
      } else
	if ( p7_sparsemask_SetFromTrace(sm, rng, gtr) != eslOK) esl_fatal(msg);


      if ( p7_SparseViterbi     (dsq, L, gm,          sm, sxf, NULL, &vsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward     (dsq, L, gm,          sm, sxf,       &fsc) != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Forward(dsq, L, gm, anch, D, sm, asf,     &asc_f) != eslOK) esl_fatal(msg);
      
      if (! diagfp)
	{
	  if (esl_FCompareAbs( fsc, asc_f, tol) != eslOK) esl_fatal(msg);
	  // bck score test here
	  if (gsc > vsc+tol)                              esl_fatal(msg);
	  if (vsc > fsc)                                  esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%20g %20g %20g\n", fsc-asc_f, vsc-gsc, fsc-vsc);
  
      p7_sparsemx_Reuse(asf);
      p7_sparsemx_Reuse(sxf);

      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
      free(dsq);
      p7_trace_Destroy(gtr);
      free(anch);
      p7_sparsemask_Destroy(sm);
    }

  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(sxf);
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
  P7_SPARSEMX   *asf = p7_sparsemx_Create(NULL);
  int            idx;
  float          vsc, fsc, asc_f;
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

	if ((cx = p7_checkptmx_Create(hmm->M, L, ESL_MBYTES(32))) == NULL) esl_fatal(msg);
	if ( p7_ForwardFilter (dsq, L, om, cx, /*fsc=*/NULL) != eslOK) esl_fatal(msg);
	if ( p7_BackwardFilter(dsq, L, om, cx, sm, p7_SPARSEMASK_THRESH_DEFAULT) != eslOK) esl_fatal(msg);

	p7_oprofile_Destroy(om);
	p7_checkptmx_Destroy(cx);
      } else 
	if ( p7_sparsemask_SetFromTrace(sm, rng, gtr) != eslOK) esl_fatal(msg);
  
      if ( p7_SparseViterbi     (dsq, L, gm,          sm, sxf, NULL, &vsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward     (dsq, L, gm,          sm, sxf,       &fsc) != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Forward(dsq, L, gm, anch, D, sm, asf,     &asc_f) != eslOK) esl_fatal(msg);
   
      if (! diagfp)
	{
	  if (esl_FCompareAbs( fsc, asc_f, tol) != eslOK) esl_fatal(msg);
	  // bck score test here
	  if (gsc > vsc+tol)                              esl_fatal(msg);
	  if (vsc > fsc)                                  esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%20g %20g %20g\n", fsc-asc_f, vsc-gsc, fsc-vsc);
  
      p7_sparsemx_Reuse(asf);
      p7_sparsemx_Reuse(sxf);

      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
      free(dsq);
      p7_trace_Destroy(gtr);
      free(anch);
      p7_sparsemask_Destroy(sm);
    }

  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(sxf);
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
  P7_SPARSEMX   *asf = p7_sparsemx_Create(NULL);
  int            idx;
  float          vsc, fsc, asc_f;
  float          tol = 0.0001;

  for (idx = 0; idx < N; idx++)
    {
      if ( p7_modelsample_AnchoredMulti(rng, M, bg, &hmm, &gm, &dsq, &L, &gtr, &anch, &D, &gsc) != eslOK) esl_fatal(msg);

      if ((sm = p7_sparsemask_Create(gm->M, L)) == NULL) esl_fatal(msg);
      
      if (idx%2) {
	if ((om = p7_oprofile_Create(hmm->M, abc)) == NULL)  esl_fatal(msg);
	if ( p7_oprofile_Convert(gm, om)           != eslOK) esl_fatal(msg);
	if ( p7_oprofile_ReconfigMultihit(om, 0)   != eslOK) esl_fatal(msg);
	om->mode = p7_LOCAL;

	if ((cx = p7_checkptmx_Create(hmm->M, L, ESL_MBYTES(32))) == NULL) esl_fatal(msg);
	if ( p7_ForwardFilter (dsq, L, om, cx, /*fsc=*/NULL) != eslOK) esl_fatal(msg);
	if ( p7_BackwardFilter(dsq, L, om, cx, sm, p7_SPARSEMASK_THRESH_DEFAULT) != eslOK) esl_fatal(msg);

	p7_oprofile_Destroy(om);
	p7_checkptmx_Destroy(cx);
      } else 
	if ( p7_sparsemask_SetFromTrace(sm, rng, gtr) != eslOK) esl_fatal(msg);
  
      if ( p7_SparseViterbi     (dsq, L, gm,          sm, sxf, NULL, &vsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward     (dsq, L, gm,          sm, sxf,       &fsc) != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Forward(dsq, L, gm, anch, D, sm, asf,     &asc_f) != eslOK) esl_fatal(msg);
   
      if (! diagfp)
	{
	  if (esl_FCompareAbs( fsc, asc_f, tol) != eslOK) esl_fatal(msg);
	  // bck score test here
	  if (gsc > vsc+tol)                              esl_fatal(msg);
	  if (vsc > fsc)                                  esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%20g %20g %20g\n", fsc-asc_f, vsc-gsc, fsc-vsc);
  
      p7_sparsemx_Reuse(asf);
      p7_sparsemx_Reuse(sxf);

      p7_hmm_Destroy(hmm);
      p7_profile_Destroy(gm);
      free(dsq);
      p7_trace_Destroy(gtr);
      free(anch);
      p7_sparsemask_Destroy(sm);
    }

  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(sxf);
  p7_bg_Destroy(bg);
}

static void
utest_viterbi(FILE *diagfp, ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, int M, int L, int N)
{
  char           msg[] = "sparse_asc_fwdback :: viterbi unit test failed";
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
  P7_TRACE      *vtr   = p7_trace_Create();
  P7_ANCHORS    *anch  = p7_anchors_Create();
  float          vsc, fsc, asc;
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

      if ((om = p7_oprofile_Create(hmm->M, abc))                                       == NULL)  esl_fatal(msg);
      if ( p7_oprofile_Convert(gm, om)                                                 != eslOK) esl_fatal(msg);
      if ((cx = p7_checkptmx_Create(hmm->M, sq->n, ESL_MBYTES(32)))                    == NULL)  esl_fatal(msg);
      if ( p7_ForwardFilter (sq->dsq, sq->n, om, cx, /*fsc=*/NULL)                     != eslOK) esl_fatal(msg);
      if ( p7_BackwardFilter(sq->dsq, sq->n, om, cx, sm, p7_SPARSEMASK_THRESH_DEFAULT) != eslOK) esl_fatal(msg);

      if ( p7_SparseViterbi (sq->dsq, sq->n, gm, sm, sxf, vtr, &vsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf,      NULL) != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxd,      NULL) != eslOK) esl_fatal(msg);
      if ( p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxd, sxd)      != eslOK) esl_fatal(msg);
      if ( p7_sparse_anchors_SetFromTrace(sxd, vtr, anch)            != eslOK) esl_fatal(msg);

      p7_logsum_InitMax();    // Zero the logsum lookup table; now FLogsum() calls will compute max()
      if ( p7_SparseForward     (sq->dsq, sq->n, gm,                   sm, sxf, &fsc) != eslOK) esl_fatal(msg);
      if ( p7_sparse_asc_Forward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asf, &asc) != eslOK) esl_fatal(msg);
      p7_logsum_Reinit();     // Reset the logsum lookup table
      
      if ( ! diagfp)
	{
	  if (esl_FCompareAbs( vsc, fsc, tol) != eslOK) esl_fatal(msg);
	  if (esl_FCompareAbs( fsc, asc, tol) != eslOK) esl_fatal(msg);
	}
      else
	fprintf(diagfp, "%4d %4d %20g %20g %20g\n", (int) sq->n, anch->D, vsc, vsc-fsc, fsc-asc);

      
      p7_anchors_Reuse(anch);
      p7_trace_Reuse(vtr);
      p7_sparsemx_Reuse(sxf);
      p7_sparsemx_Reuse(sxd);	
      p7_sparsemx_Reuse(asf);	
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
  p7_sparsemask_Destroy(sm);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
}

#endif /*p7SPARSE_ASC_FWDBACK_TESTDRIVE*/

/*------------------- end, unit tests ---------------------------*/



/*****************************************************************
 * 4. Test driver
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

      if      (strcmp(which, "compare_reference") == 0) utest_compare_reference(stdout, rng, abc, M, L, N);
      else if (strcmp(which, "singlesingle")      == 0) utest_singlesingle     (stdout, rng, abc, M,    N);
      else if (strcmp(which, "multisingle")       == 0) utest_multisingle      (stdout, rng, abc, M,    N);
      else if (strcmp(which, "multipath_local")   == 0) utest_multipath_local  (stdout, rng, abc, M,    N);
      else if (strcmp(which, "multimulti")        == 0) utest_multimulti       (stdout, rng, abc, M,    N);
      else if (strcmp(which, "viterbi")           == 0) utest_viterbi          (stdout, rng, abc, M, L, N);
      else esl_fatal("--diag takes: compare_reference, singlesingle, multisingle");
    }
  else
    {
      fprintf(stderr, "## %s\n", argv[0]);
      fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

      utest_viterbi          (NULL, rng, abc, M, L, N);

      utest_compare_reference(NULL, rng, abc, M, L, N);
      utest_singlesingle     (NULL, rng, abc, M,    N);
      utest_multisingle      (NULL, rng, abc, M,    N);
      utest_multipath_local  (NULL, rng, abc, M,    N);
      utest_multimulti       (NULL, rng, abc, M,    N);

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
 * 5. Example
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
  P7_ANCHORS     *anch    = p7_anchors_Create();
  float          *wrk     = NULL;
  P7_ANCHORHASH  *ah      = p7_anchorhash_Create();
  float           fsc, vsc, asc, asc_sparse;
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
    p7_BackwardFilter(sq->dsq, sq->n, om, cx, sm, p7_SPARSEMASK_THRESH_DEFAULT);
  }

  /* We need an anchor set <anch>.
   * To get it, run the reference prototype code;
   * (we don't have the MPAS algorithm in its sparse production form yet)
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


  //p7_refmx_Dump(stdout, afu);
  //p7_refmx_Dump(stdout, afd);

  /* Finally...
   * Run sparse ASC Forward.
   */
  asf = p7_sparsemx_Create(sm);
  p7_sparse_asc_Forward(sq->dsq, sq->n, gm, anch->a, anch->D, sm, asf, &asc_sparse);

  

  //p7_spascmx_Dump(stdout, asf, anch->arr, anch->n);

  printf("Reference ASC fwd score = %.2f nats\n", asc);
  printf("Sparse ASC fwd score    = %.2f nats\n", asc_sparse);

  p7_anchorhash_Destroy(ah);
  if (wrk) free(wrk);
  p7_sparsemx_Destroy(asf);
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




/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
