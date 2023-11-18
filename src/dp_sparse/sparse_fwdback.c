/* Sparse dynamic programming.
 * 
 * Contents:
 *   1. Sparse Forward
 *   2. Sparse Backward
 *   3. Benchmark driver
 *   4. Unit tests
 *   5. Test driver
 *   6. Example 
 */
#include <p7_config.h>

#include <string.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "base/p7_profile.h"

#include "misc/logsum.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_fwdback.h"

/*****************************************************************
 * 1. Sparse Forward
 *****************************************************************/

/* Function:  p7_SparseForward()
 * Synopsis:  Forward algorithm, in sparse DP.
 *
 * Purpose:   Compare profile <gm> to digital sequence <dsq> (of length
 *            <L>), by the Forward algorithm, using sparse dynamic
 *            programming, restricted to the mask <sm>.  Fill in the
 *            sparse DP matrix <sx>, and the Forward score is also
 *            optionally returned in <opt_sc>.
 *            
 *            <sx> can be reused from previous calculations, even
 *            smaller ones; see <p7_sparsemx_Reuse()>. If necessary,
 *            it will be reallocated here, to be large enough for the
 *            <gm->M> by <L> calculation restricted to masked cells
 *            <sm>.
 *
 * Args:      dsq     -  target sequence 1..L
 *            L       -  length of <dsq>
 *            gm      -  profile
 *            sm      -  sparse mask
 *            sx      -  Forward matrix to fill (may be reallocated here)
 *            opt_sc  -  optRETURN: raw Forward score, in nats
 *
 * Returns:   <eslOK> on success. Matrix <sx> may be reallocated, and is filled.
 *            The raw Forward score in nats is optionally in <*opt_sc>.
 *            
 * Throws:    <eslEMEM> if we fail to reallocate a large enough <sx>.           
 */
int
p7_SparseForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_SPARSEMX *sx, float *opt_sc)
{
  float const *tsc    = gm->tsc;	 /* sets up TSC() macro, access to profile's transitions      */
  float       *xpc;	                 /* ptr that steps through current special cells              */
  float       *dpc;	                 /* ptr to step thru current row i main DP cells              */
  float       *dpp;			 /* ptr to step thru previous row i-1 main DP cells           */
  float       *last_dpc;		 /* used to reinit dpp after each sparse row computation      */
  float const *rsc;			 /* will be set up for MSC(), ISC() macros for residue scores */
  int          ng;
  float        xE, xN, xJ, xB, xL, xG, xC;  /* tmp scores on special states. only stored when in row bands, and on ia-1 before a seg */
  float        mlc, mgc;		 /* temporary score calculations M(i,k)                       */
  float        dlc, dgc;		 /* precalculated D(i,k+1) value on current row               */
  int         *kc = sm->k[0];		 /* <kc> points to list of sparse indices k for curr row i    */
  int         *kp;			 /* <kp> points to the previous row's sparse cell index list  */
  int          i,k;	      	         /* i,k row,col (seq position, profile position) cell coords  */
  int          y,z;			 /* indices in lists of k coords on prev, current row         */
  int          status;

  /* Contract checks on arguments */
  ESL_DASSERT1( (sm->L == L) );
  ESL_DASSERT1( (sm->M == gm->M) );

  /* Assure that <sx> is allocated large enough (we might be reusing it).
   * Set its type now, so we can Dump/Validate/etc. during debugging this routine, if needed.
   */
  if ( (status = p7_sparsemx_Reinit(sx, sm)) != eslOK) return status;
  sx->type = p7S_FORWARD;

  xN  = 0.0f;
  xJ  = -eslINFINITY;
  xC  = -eslINFINITY;
  ng  = 0;
  xpc = sx->xmx;
  dpc = sx->dp;
  for (i = 1; i <= L; i++)
    {
      if (! sm->n[i]) { ng++; continue; }   /* skip rows that have no included cells */

      /* Reinitialize and store specials for row ia-1 just outside sparsified segment */
      if (i == 1 || ng) {   // need for i=1 here may be confusing. here's why: We only reach i == 1 test when sm->n[1] > 0, i.e. seg[1].ia = 1, ia-1 = 0, ng = 0
	*xpc++ = xE = -eslINFINITY;
	*xpc++ = xN  = xN + ( ng ? ng * gm->xsc[p7P_N][p7P_LOOP] : 0.0); /* test ng, because we must watch out for 0*-inf special case */
	*xpc++ = xJ  = xJ + ( ng ? ng * gm->xsc[p7P_J][p7P_LOOP] : 0.0);
	*xpc++ = xB  = p7_FLogsum( xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]);
	*xpc++ = xL  = xB + gm->xsc[p7P_B][0]; /* B->L */
	*xpc++ = xG  = xB + gm->xsc[p7P_B][1]; /* B->G */
	*xpc++ = xC  = xC + ( ng ? ng * gm->xsc[p7P_C][p7P_LOOP] : 0.0);
	*xpc++       = -eslINFINITY; /* JJ: this space only used in a Decoding matrix. */
	*xpc++       = -eslINFINITY; /* CC: this space only used in a Decoding matrix. */
	ng = 0;
      }

      rsc = gm->rsc[dsq[i]];	/* now MSC(k), ISC(k) residue score macros work */
      last_dpc = dpc;		/* remember where dpc started; dpp will be set here after we finish each row calculation */

      kp = kc;                  /* last row we did becomes prev row now; ready to step through k indices of previous row's sparse cells */
      kc = sm->k[i];		/* ditto for current row i */

      dlc = dgc = xE = -eslINFINITY;
      for (z=0, y=0; z < sm->n[i]; z++) /* Iterate over the one or more sparse cells (i,k) that we calculate on this row. */
	{
	  k = kc[z]; /* next sparse cell to calculate: (i,k) */
	  
	  /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	  mlc = xL  + TSC(p7P_LM, k-1);
	  mgc = xG  + TSC(p7P_GM, k-1);
	  while (y < sm->n[i-1] && kp[y] < k-1)  { y++; dpp+=p7S_NSCELLS; }
	  if    (y < sm->n[i-1] && kp[y] == k-1) {
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

	  /* Try to find cell i-1,k; then compute I(i,k) from it */
	  while (y < sm->n[i-1] && kp[y] < k)  { y++; dpp+=p7S_NSCELLS; }
	  if    (y < sm->n[i-1] && kp[y] == k) {
	    *dpc++ = p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MI,k),  dpp[p7R_IL] + TSC(p7P_II, k)); // +ISC(k) if we weren't enforcing it to zero
	    *dpc++ = p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MI,k),  dpp[p7R_IG] + TSC(p7P_II, k)); // ditto
	  } else {
	    *dpc++ = -eslINFINITY;
	    *dpc++ = -eslINFINITY;
	  }
	    
	  /* local exit paths */
	  xE = p7_FLogsum(xE, p7_FLogsum(mlc, dlc));

	  /* delayed store of Dk; advance calculation of next D_k+1 */
	  *dpc++ = dlc;
	  *dpc++ = dgc;
	  if (z < sm->n[i]-1 && kc[z+1] == k+1) { /* is there a (i,k+1) cell to our right? */
	    dlc = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
	    dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	  } else {   /* if not, we MUST add {MD}l->Dk+1..E glocal exit path, even from internal sparse cells - not just last cell! */
	    xE  = p7_FLogsum( xE, TSC(p7P_DGE, k) + p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k))); // this line comes first BEFORE you set dgc to -inf, fool
	    dlc = dgc = -eslINFINITY;
	  }
	}
      *xpc++ = xE;  /* we already accumulated all Mk->E exits, both local and glocal */
      *xpc++ = xN = xN + gm->xsc[p7P_N][p7P_LOOP];
      *xpc++ = xJ = p7_FLogsum( xJ + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_LOOP]);
      *xpc++ = xB = p7_FLogsum( xJ + gm->xsc[p7P_J][p7P_MOVE],  xN + gm->xsc[p7P_N][p7P_MOVE]);
      *xpc++ = xL = xB + gm->xsc[p7P_B][0]; /* B->L */
      *xpc++ = xG = xB + gm->xsc[p7P_B][1]; /* B->G */
      *xpc++ = xC = p7_FLogsum( xE + gm->xsc[p7P_E][p7P_MOVE],  xC + gm->xsc[p7P_C][p7P_LOOP]);
      *xpc++      = -eslINFINITY; /* JJ: this space only used in a Decoding matrix. */
      *xpc++      = -eslINFINITY; /* CC: this space only used in a Decoding matrix. */

      /* now dpc is on the start of the next sparsified row */
      dpp = last_dpc;
    }

  if (opt_sc != NULL) *opt_sc = xC + ( ng ? ng *  gm->xsc[p7P_C][p7P_LOOP] : 0.0f) + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}
/*--------------- end, sparse Forward  --------------------------*/




/*****************************************************************
 * 2. Sparse Backward
 *****************************************************************/

/* Function:  p7_SparseBackward()
 * Synopsis:  Backward algorithm, in sparse DP.
 *
 * Purpose:   Compare profile <gm> to digital sequence <dsq> (of length
 *            <L>) by the Backward algorithm, using sparse dynamic
 *            programming, as constrained by the sparse mask <sm>.
 *            Fill in the sparse DP Backward matrix <sx> and
 *            (optionally) return the overall raw Backward score in
 *            <*opt_sc>.
 *
 * Args:      dsq     - digital target sequence 1..L
 *            L       - length of <dsq>
 *            gm      - profile
 *            sm      - sparse mask
 *            sx      - Backward matrix to fill; may get reallocated here
 *            opt_sc  - optRETURN: raw Backwards score, in nats
 *                      
 * Returns:   <eslOK> on success. Matrix <sx> is filled (and may be
 *            reallocated), and the raw Backward score is in <*opt_sc>.
 *            
 * Throws:    <eslEMEM> if we fail to reallocated a large enough <sx>.           
 */
int
p7_SparseBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_SPARSEMX *sx, float *opt_sc)
{
  const float    *tsc = gm->tsc; // ptr into model transition score vector; enables TSC(k,s) macro shorthand 
  const float    *rsc = NULL;	 // ptr into emission scores on current row i; enables MSC(k),ISC(k) macro shorthand 
  const float    *rsn;		 // ptr into emission scores on next row i+1; enables MSN(k),ISN(k) 
  float          *xp;            // ptr that steps through current cells 
  float          *dpc;		 // ptr stepping through current DP row i sparse cells 
  float          *dpn = NULL;	 // ptr stepping through next DP row (i+1) sparse cells; NULL when no i+1 row stored 
  float          *last_dpc;	 // remembers where dpc started, so dpn can be set for next i 
  int    i,k,y,z;
  float  xC,xG,xL,xB,xJ,xN,xE;	 // tmp vars for special cells on this row i
  float  dlc,dgc;		 // tmp vars for D vals on current row i  
  float  mln,mgn,iln,ign;	 // tmp vars for M+e, I+e on next row i+1 
  int    status;

  /* Contract checks on arguments */
  ESL_DASSERT1( (sm->L == L) );
  ESL_DASSERT1( (sm->M == gm->M) );

  /* Assure that <sx> is allocated large enough (we might be reusing it).
   * Set its type now, so we can Dump/Validate/etc. during debugging this routine, if needed.
   */
  if ( (status = p7_sparsemx_Reinit(sx, sm)) != eslOK) return status;
  sx->type = p7S_BACKWARD;

  /* In Backwards, we traverse DP matrix backwards; init ptrs to last elements */
  xp  = sx->xmx + (sm->nrow + sm->S - 1)*p7S_NXCELLS; // last supercell in xmx 
  dpc = sx->dp  + (sm->ncells-1)*p7S_NSCELLS;	      // last supercell in dp  

  /* xC/N/J follow a convention that as we start a row i, they already include
   * all contribution from row i+1; i.e. they are calculated as a lookahead, 
   * then immediately stored.
   */
  xC = gm->xsc[p7P_C][p7P_MOVE];  // xC[L] can only be reached via C->T transition 
  xJ = xN = -eslINFINITY;	  // J,N unreached on L, and indeed until at least first ib (last stored row)
  xG = xL = -eslINFINITY;	  // L,G can only be reached from a stored next row i+1.

  for (i = sm->L; i >= 1; i--)
    {
      /* if this row isn't stored, then it can only extend CC/JJ/NN:
       * extend them to prev row (remember precalculation convention)
       */
      if (! sm->n[i]) 
	{ 
	  xC += gm->xsc[p7P_C][p7P_LOOP];
	  xJ += gm->xsc[p7P_J][p7P_LOOP];
	  xN += gm->xsc[p7P_N][p7P_LOOP];
	  continue;
	}

      /* Else we continue: this row i is stored, both as a main row (in dpc) and specials (in xmx) */
      xp[p7S_CC] = -eslINFINITY;     // CC only stored in a decoding matrix 
      xp[p7S_JJ] = -eslINFINITY;     // ditto JJ 
      xp[p7S_C]  = xC;               // deferred store of xC. 
      xp[p7S_G]  = xG;               // on i=ib segment start, xG cannot be reached; it will be -inf
      xp[p7S_L]  = xL;               // ditto xL
      xp[p7S_B]  = xB = p7_FLogsum( xL + gm->xsc[p7P_B][0],        xG + gm->xsc[p7P_B][1]);         // on i=ib segment start, evaluates to -inf
      xp[p7S_J]  = xJ = p7_FLogsum( xJ, 		           xB + gm->xsc[p7P_J][p7P_MOVE]);  // on i=ib, evaluates to xJ
      xp[p7S_N]  = xN = p7_FLogsum( xN, 	                   xB + gm->xsc[p7P_N][p7P_MOVE]);  // ditto xN
      xp[p7S_E]  = xE = p7_FLogsum( xJ + gm->xsc[p7P_E][p7P_LOOP], xC + gm->xsc[p7P_E][p7P_MOVE]);  // on i=ib, i<L, E->{JC} are both possible; i=L, evaluates to xC + tEC.
      xp -= p7S_NXCELLS;

      last_dpc = dpc;               // last_dpc will remember where dpc was (at end of row i); before we loop back to new i, dpn will be set to last_dpc.
      rsn      = rsc;               // from previous row. invalid if there is no i+1 row (indeed, will be NULL for i=L) but doesn't matter; won't be accessed in that case
      rsc      = gm->rsc[dsq[i]];   // MSC(k),ISC(k) macros now work

      dlc = -eslINFINITY; 
      dgc = xE + TSC(p7P_DGE, sm->k[i][sm->n[i]-1]);
      xG  = xL = -eslINFINITY;               // prepare to accumulate xG, xL logsums
      y = (i == sm->L ? -1 : sm->n[i+1]-1);  // if n[i+1] were allocated and set to a 0 sentinel, we could avoid this branch
      for (z = sm->n[i]-1; z >= 0; z--)
	{
	  k = sm->k[i][z];

	  /* try to pick up mln, mgn from i+1,k+1 (and add their emission score) */
	  while (y >= 0 && sm->k[i+1][y]  > k+1) { y--; dpn -= p7S_NSCELLS; } // note if y were already on i+1,k+1, it doesn't move
	  if    (y >= 0 && sm->k[i+1][y] == k+1) {
	    mln = dpn[p7S_ML] + MSN(k+1);
	    mgn = dpn[p7S_MG] + MSN(k+1);
	  } else { mln = mgn = -eslINFINITY; }
	  
	  /* try to pick up iln,ign from i+1,k */
	  while (y >= 0 && sm->k[i+1][y]  > k) { y--; dpn -= p7S_NSCELLS; } // note if y were already on i+1,k+1, it doesn't move
	  if    (y >= 0 && sm->k[i+1][y] == k) {
	    iln = dpn[p7S_IL]; // + ISN(k), if it weren't fixed to 0
	    ign = dpn[p7S_IG]; // + ISN(k), ditto
	  } else { iln = ign = -eslINFINITY; }

	  /* see if dlc, dgc are valid from prev calculation on this row. if not reinit. note glocal path, wing-retracted. */
	  if (z < sm->n[i]-1 && sm->k[i][z+1] != k+1) { 
	    dlc = -eslINFINITY; 
	    dgc = xE + TSC(p7P_DGE, k);
	  }
	  
	  /* M(i,k) calculations need to use dgc,dlc before we
	   * change them, while they're still D(i,k+1). 
	   */
	  dpc[p7S_ML] = p7_FLogsum( p7_FLogsum(mln + TSC(p7P_MM, k),      // ML(i,k) =   ML(i+1,k+1) * t(k,MM) * e_M(k+1, x_i+1)      | mln = log[ ML(i+1,k+1) * e_M(k+1,x_i+1)]
					       iln + TSC(p7P_MI, k)),     //           + IL(i+1,k)   * t(k,MI) * e_I(k, x_i+1)        | iln = log[ IL(i+1,k)   * e_I(k,  x_i+1)]
				    p7_FLogsum(dlc + TSC(p7P_MD, k),      //           + DL(i,  k+1) * t(k,DD)                        | dlc = DL(i,k+1), wrapped around from prev loop iteration. This has tricky boundary conditions at k=kbc, and k=kbc=M
					       xE));                      //           +  E(i)       * t(MkE)=1.0
	  dpc[p7S_MG] = p7_FLogsum( p7_FLogsum(mgn + TSC(p7P_MM, k),      // MG(i,k) is essentially the same recursion, without a transition to E.
					       ign + TSC(p7P_MI, k)),     
					       dgc + TSC(p7P_MD, k));     
	      
	  /* Accumulate xG, xL as we sweep over the current row; 
	   * they get used to initialize the next row (i-1). Note
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
	  dpc[p7S_DL] = dlc = p7_FLogsum( p7_FLogsum(mln  + TSC(p7P_DM, k), 
						     dlc  + TSC(p7P_DD, k)),
					             xE);                   
	  dpc[p7S_DG] = dgc = p7_FLogsum( mgn  + TSC(p7P_DM, k),  
					  dgc  + TSC(p7P_DD, k)); 

	  dpc -= p7S_NSCELLS;
	} // end loop over sparse cells k on this row.
      dpn = last_dpc;

      /* precalculate what xC/xJ/xN will be on previous row i-1... these values get stored as we roll around i loop */
      xC += gm->xsc[p7P_C][p7P_LOOP];
      xJ += gm->xsc[p7P_J][p7P_LOOP];
      xN += gm->xsc[p7P_N][p7P_LOOP];

      /* If we are terminating a segment (if this i is an ia), we need to 
       * store specials for row ia-1. This works for i=1 and storing the
       * final i=0 special row too.
       */
      if (sm->n[i-1] == 0) 
	{
	  xp[p7S_CC] = -eslINFINITY;	/* CC only stored in a Decoding matrix. */
	  xp[p7S_JJ] = -eslINFINITY;    /* JJ only stored in a Decoding matrix. */
	  xp[p7S_C]  = xC;
	  xp[p7S_G]  = xG;
	  xp[p7S_L]  = xL;
	  xp[p7S_B]  = xB = p7_FLogsum( xL + gm->xsc[p7P_B][0],        xG + gm->xsc[p7P_B][1]);
	  xp[p7S_J]  = xJ = p7_FLogsum( xJ,     	               xB + gm->xsc[p7P_J][p7P_MOVE]);
	  xp[p7S_N]  = xN = p7_FLogsum( xN,		               xB + gm->xsc[p7P_N][p7P_MOVE]);
	  xp[p7S_E]  = xE = p7_FLogsum( xJ + gm->xsc[p7P_E][p7P_LOOP], xC + gm->xsc[p7P_E][p7P_MOVE]);
	  xp -= p7S_NXCELLS;

	  xG = xL = -eslINFINITY; /* when we start the next segment at an ib, these will get stored */
	}
    } // end loop over i

  if (opt_sc) *opt_sc = xN;	/* tS->N is 1.0, no cost. */
  return eslOK;
}
/*----------------- end, sparse Backwards ----------------------*/



/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/
#ifdef p7SPARSE_FWDBACK_BENCHMARK

/* As of 13 Aug 2012, on cluster node:
 * # Caudal_act (M= 143)
 * # Sparse Forward:      13.8 Mc/s
 * # Sparse Backward:     13.1 Mc/s
 * # Sparse Decoding:     11.3 Mc/s
 * # Sparse Viterbi:      64.3 Mc/s
 */
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,   "2000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for sparse dual local/glocal Forward/Backward implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  P7_SPARSEMX    *sxv     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxf     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxb     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxd     = p7_sparsemx_Create(NULL);
  P7_TRACE       *tr      = p7_trace_Create();
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, F_time, B_time, V_time, D_time;
  double          F_speed, B_speed, V_speed, D_speed;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, L);

  sm  = p7_sparsemask_Create(gm->M, L);
  p7_sparsemask_AddAll(sm);

  /* Baseline time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Forward benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_SparseForward  (dsq, L, gm, sm, sxf, &sc);
      p7_sparsemx_Reuse (sxf);
    }
  esl_stopwatch_Stop(w);
  F_time  = w->user - base_time;
  F_speed = (double) N * (double) L * (double) gm->M * 1e-6 / F_time;

  /* Backward benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_SparseBackward (dsq, L, gm, sm, sxb, &sc);
      p7_sparsemx_Reuse (sxb);
    }
  esl_stopwatch_Stop(w);
  B_time  = w->user - base_time;
  B_speed = (double) N * (double) L * (double) gm->M * 1e-6 / B_time;

  /* Viterbi benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_SparseViterbi  (dsq, L, gm, sm, sxv, tr, &sc);
      p7_sparsemx_Reuse (sxv);
      p7_trace_Reuse    (tr);
    }
  esl_stopwatch_Stop(w);
  V_time  = w->user - base_time;
  V_speed = (double) N * (double) L * (double) gm->M * 1e-6 / V_time;

  /* Decoding benchmark time (by differencing out F/B part) */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_SparseForward  (dsq, L, gm, sm, sxf, &sc);
      p7_SparseBackward (dsq, L, gm, sm, sxb, &sc);
      p7_SparseDecoding (dsq, L, gm, sxf, sxb, sxd);

      p7_sparsemx_Reuse (sxf);
      p7_sparsemx_Reuse (sxb);
      p7_sparsemx_Reuse (sxd);
      p7_trace_Reuse    (tr);
    }
  esl_stopwatch_Stop(w);
  D_time  = (w->user - base_time) - F_time - B_time;
  D_speed = (double) N * (double) L * (double) gm->M * 1e-6 / D_time;

  printf("# %s (M= %d)\n", gm->name, gm->M);
  printf("# Sparse Forward:  %8.1f Mc/s\n", F_speed);
  printf("# Sparse Backward: %8.1f Mc/s\n", B_speed);
  printf("# Sparse Decoding: %8.1f Mc/s\n", D_speed);
  printf("# Sparse Viterbi:  %8.1f Mc/s\n", V_speed);

  free(dsq);
  p7_sparsemx_Destroy(sxv);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemask_Destroy(sm);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_FWDBACK_BENCHMARK*/
/*------------------ end, benchmark -----------------------------*/



/*****************************************************************
 * 4. Unit tests
 *****************************************************************/

#ifdef p7SPARSE_FWDBACK_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_dsq.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "base/p7_bg.h"

#include "build/modelsample.h"
#include "build/p7_builder.h"
#include "build/seqmodel.h"
#include "search/modelconfig.h"
#include "misc/emit.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"

#include "dp_sparse/sparse_viterbi.h"
#include "dp_sparse/sparse_trace.h"
#include "dp_sparse/sparse_decoding.h"

#include "dp_reference/reference_fwdback.h"
#include "dp_reference/reference_viterbi.h"
#include "dp_reference/reference_trace.h"
#include "dp_reference/reference_decoding.h"

/*  The "randomseq" utest compares a randomly sampled profile to
 *  random iid sequences by sparse DP (including fwd/bck filter step
 *  to generate the sparse mask), and tests:
 *     
 *  1. Forward score is always >= Viterbi
 *  2. Forward score = Backwards score
 *  3. F,B,V matrices pass Validate()
 *  
 *  It also runs full reference DP, and tests:
 *  
 *  4. Reference F,B score >= sparse F,B score 
 *     (but not necessarily very close, esp on iid seqs)
 *  5. Reference V score >= sparse V score.
 *  6. All values v in reference, sparse F,B,V matrices 
 *     satisfy v_ref >= v_sparse
 *     
 *  These are relatively weak constraints, but one advantage of the
 *  unit test is that it tests the entire sparse calculation sequence,
 *  including fwd/bck filtering to generate the sparse mask, which can
 *  then include the most tricky indexing situations (multisegment
 *  sparse masks). Contrast to the "compare-reference" utest, which
 *  guarantees strong constraints (exact match to the reference
 *  implementation, within roundoff error) but only by using "full"
 *  sparse masks with all cells marked; or to the "singlepath" utest,
 *  which also guarantees strong constraints, but only by restricting
 *  to a special singlepath profile with uniglocal,L=0, which cannot
 *  generate a multisegment sparse path.
 */
static void
utest_randomseq(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char           msg[]  = "sparse fwdback 'randomseq' unit test failed";
  P7_HMM        *hmm    = NULL;
  P7_PROFILE    *gm     = p7_profile_Create(M, abc);
  P7_OPROFILE   *om     = p7_oprofile_Create(M, abc);
  ESL_DSQ       *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_CHECKPTMX  *ox     = p7_checkptmx_Create(M, L, ESL_MBYTES(32));
  P7_SPARSEMASK *sm     = p7_sparsemask_Create(M, L);
  P7_SPARSEMX   *sxv    = p7_sparsemx_Create(sm);
  P7_SPARSEMX   *sxf    = p7_sparsemx_Create(sm);
  P7_SPARSEMX   *sxb    = p7_sparsemx_Create(sm);
  P7_REFMX      *rxv    = p7_refmx_Create(M, L);
  P7_REFMX      *rxf    = p7_refmx_Create(M, L);
  P7_REFMX      *rxb    = p7_refmx_Create(M, L);
  float          fsc_s, bsc_s, vsc_s; /* sparse DP scores for fwd, bck, vit    */
  float          fsc_r, bsc_r, vsc_r; /* reference DP scores for fwd, bck, vit */
  int            idx;
  float          tol  = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
   
  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);
  if ( p7_oprofile_Convert(gm, om)       != eslOK) esl_fatal(msg);
  if ( p7_profile_SetLength(gm, L)       != eslOK) esl_fatal(msg);
  if ( p7_oprofile_ReconfigLength(om, L) != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      /* Sample a random iid sequence */
      if (esl_rsq_xfIID(rng, bg->f, abc->K, L, dsq)      != eslOK) esl_fatal(msg);

      /* Fwd/Bck local filter to calculate the sparse mask */
      if ( p7_checkptmx_Reinit(ox, M, L)                             != eslOK) esl_fatal(msg);
      if ( p7_ForwardFilter (dsq, L, om, ox, /*fsc=*/NULL)           != eslOK) esl_fatal(msg);
      if ( p7_BackwardFilter(dsq, L, om, ox, sm, p7_SPARSIFY_THRESH) != eslOK) esl_fatal(msg);

      /* Sparse DP calculations */
      if ( p7_SparseViterbi   (dsq, L, gm, sm, sxv, NULL, &vsc_s) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward   (dsq, L, gm, sm, sxf, &fsc_s)       != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward  (dsq, L, gm, sm, sxb, &bsc_s)       != eslOK) esl_fatal(msg);
  
      /* Tests (1)-(3) */
      if (fsc_s < vsc_s+tol)                         esl_fatal(msg); /* 1. Fwd score >= Vit score */
      if (fabs(fsc_s - bsc_s) > tol)                 esl_fatal(msg); /* 2. Fwd score == Bck score */
      if ( p7_sparsemx_Validate(sxv, NULL) != eslOK) esl_fatal(msg); /* 3. F,B,V matrices must Validate() */
      if ( p7_sparsemx_Validate(sxf, NULL) != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Validate(sxb, NULL) != eslOK) esl_fatal(msg);
      
      /* Reference DP calculations */
      if (p7_ReferenceForward (dsq, L, gm, rxf,       &fsc_r)  != eslOK) esl_fatal(msg);
      if (p7_ReferenceBackward(dsq, L, gm, rxb,       &bsc_r)  != eslOK) esl_fatal(msg);
      if (p7_ReferenceViterbi (dsq, L, gm, rxv, NULL, &vsc_r)  != eslOK) esl_fatal(msg);

      /* Tests (4)-(6) */
      if (fsc_s > fsc_r+tol) esl_fatal(msg); /* 4. Reference F,B score >= sparse F,B score */
      if (vsc_s > vsc_r+tol) esl_fatal(msg); /* 5. Reference V score >= sparse V score     */
      if ( p7_sparsemx_CompareReferenceAsBound(sxv, rxv, tol) != eslOK) esl_fatal(msg); /* 6. All matrix values satisfy v_ref >= v_sparse */
      if ( p7_sparsemx_CompareReferenceAsBound(sxf, rxf, tol) != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_CompareReferenceAsBound(sxb, rxb, tol) != eslOK) esl_fatal(msg);
      
      /* reuse DP matrices and mask */
      if ( p7_refmx_Reuse(rxv)     != eslOK) esl_fatal(msg);
      if ( p7_refmx_Reuse(rxf)     != eslOK) esl_fatal(msg);
      if ( p7_refmx_Reuse(rxb)     != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxv)  != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxf)  != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxb)  != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_Reuse(sm) != eslOK) esl_fatal(msg); 
    }

  free(dsq);
  p7_sparsemask_Destroy(sm);
  p7_sparsemx_Destroy(sxv);  p7_sparsemx_Destroy(sxf);  p7_sparsemx_Destroy(sxb);
  p7_refmx_Destroy(rxv);     p7_refmx_Destroy(rxf);     p7_refmx_Destroy(rxb);
  p7_checkptmx_Destroy(ox);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
  p7_hmm_Destroy(hmm);
}



/* "compare_reference" unit test.
 * 
 *  Tests that the sparse implementations of Viterbi, Forward,
 *  Backward, and Decoding give the same results as
 *  reference implementation, when all cells are included
 *  in the sparse mask.  
 *  
 *  Samples a random profile and compare it to 'homologous'
 *  sequences, generated from the profile, with the sparse
 *  mask set to mark all cells, and tests:
 *  
 * 1. Reference and sparse scores must be identical (V,F,B), 
 *    within a tolerance;
 * 2. Reference and sparse viterbi traces must be identical;
 * 3. All sparse V,F,B,D DP matrix structures must Validate();
 * 4. Sparse Viterbi trace structure must Validate();
 * 5. Reference and sparse matrices are identical within tolerance. 
 */
static void
utest_compare_reference(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char           msg[]  = "sparse fwdback, compare-reference unit test failed";
  P7_HMM        *hmm    = NULL;
  P7_PROFILE    *gm     = p7_profile_Create(M, abc);
  ESL_SQ        *sq     = esl_sq_CreateDigital(abc);       /* space for generated (homologous) target seqs              */
  P7_SPARSEMASK *sm     = p7_sparsemask_Create(M, L);
  P7_SPARSEMX   *sxv    = p7_sparsemx_Create(sm);
  P7_SPARSEMX   *sxf    = p7_sparsemx_Create(sm);
  P7_SPARSEMX   *sxb    = p7_sparsemx_Create(sm);
  P7_SPARSEMX   *sxd    = p7_sparsemx_Create(sm);
  P7_REFMX      *rxv    = p7_refmx_Create(M, L);
  P7_REFMX      *rxf    = p7_refmx_Create(M, L);
  P7_REFMX      *rxb    = p7_refmx_Create(M, L);
  P7_REFMX      *rxd    = p7_refmx_Create(M, L);
  P7_TRACE      *rtr    = p7_trace_Create();
  P7_TRACE      *str    = p7_trace_Create();
  int            idx;
  float          vsc_s, fsc_s, bsc_s;
  float          vsc_r, fsc_r, bsc_r;
  float          tol = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  char           errbuf[eslERRBUFSIZE];
  
  /* Sample a profile. Config as usual, multihit dual-mode local/glocal. */
  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      /* Generate (sample) a sequence from the profile */
      if ( p7_profile_SetLength(gm, L)  != eslOK) esl_fatal(msg);   /* config to generate mean length of L (length was probably reset by last emitted seq) */
      do {
	esl_sq_Reuse(sq);
	p7_ProfileEmit(rng, hmm, gm, bg, sq, NULL);
      } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
      if ( p7_profile_SetLength(gm, sq->n) != eslOK) esl_fatal(msg);

      /* Mark all cells in sparse mask */
      if ( p7_sparsemask_Reinit(sm, M, sq->n) != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_AddAll(sm)           != eslOK) esl_fatal(msg);

      /* Sparse DP calculations  */
      if ( p7_SparseViterbi (sq->dsq, sq->n, gm, sm, sxv, str, &vsc_s) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf,      &fsc_s) != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxb,      &bsc_s) != eslOK) esl_fatal(msg);
      if ( p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxb, sxd)        != eslOK) esl_fatal(msg);

      /* Reference DP calculations */
      if ( p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxv, rtr, &vsc_r) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc_r) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceBackward(sq->dsq, sq->n, gm, rxb,      &bsc_r) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxb, rxd)    != eslOK) esl_fatal(msg);

      //p7_sparsemx_Dump(stdout, sxd);
      //p7_refmx_Dump(stdout, rxd);

      /* Tests */
      if ( esl_FCompare(vsc_s, vsc_r, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg); /* (1): reference and sparse scores (V,F,B) identical within tolerance */
      if ( esl_FCompare(fsc_s, fsc_r, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg);
      if ( esl_FCompare(bsc_s, bsc_r, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg);
      if ( p7_trace_CompareLoosely(rtr, str, sq->dsq)    != eslOK) esl_fatal(msg); /* (2): reference, sparse Viterbi tracebacks identical; see notes on p7_trace_CompareLoosely */
      if ( p7_sparsemx_Validate(sxv, NULL)               != eslOK) esl_fatal(msg); /* (3): All sparse DP matrices Validate() (V,F,B,D) */
      if ( p7_sparsemx_Validate(sxf, NULL)               != eslOK) esl_fatal(msg); /*      (the <NULL> arg is an optional <errbuf>)    */
      if ( p7_sparsemx_Validate(sxb, NULL)               != eslOK) esl_fatal(msg);       
      if ( p7_sparsemx_Validate(sxd, errbuf)             != eslOK) esl_fatal(errbuf);
      if ( p7_trace_Validate(str, abc, sq->dsq, NULL)    != eslOK) esl_fatal(msg); /* (4): Sparse DP Viterbi trace must Validate() */
      if ( p7_sparsemx_CompareReference(sxv, rxv, tol)   != eslOK) esl_fatal(msg); /* (5): Sparse and reference DP matrices all identical within tolerance */
      if ( p7_sparsemx_CompareReference(sxf, rxf, tol)   != eslOK) esl_fatal(msg); /* (5): Sparse and reference DP matrices all identical within tolerance */
      if ( p7_sparsemx_CompareReference(sxb, rxb, tol)   != eslOK) esl_fatal(msg); /* (5): Sparse and reference DP matrices all identical within tolerance */
      if ( p7_sparsemx_CompareReference(sxd, rxd, tol)   != eslOK) esl_fatal(msg); /* (5): Sparse and reference DP matrices all identical within tolerance */
      
      if ( p7_trace_Reuse(str)     != eslOK) esl_fatal(msg);
      if ( p7_trace_Reuse(rtr)     != eslOK) esl_fatal(msg);
      if ( p7_refmx_Reuse(rxv)     != eslOK) esl_fatal(msg);
      if ( p7_refmx_Reuse(rxf)     != eslOK) esl_fatal(msg);
      if ( p7_refmx_Reuse(rxb)     != eslOK) esl_fatal(msg);
      if ( p7_refmx_Reuse(rxd)     != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxv)  != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxf)  != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxb)  != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxd)  != eslOK) esl_fatal(msg);
      if ( esl_sq_Reuse(sq)        != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_Reuse(sm) != eslOK) esl_fatal(msg); 
    }

  p7_trace_Destroy(rtr);
  p7_trace_Destroy(str);
  p7_sparsemask_Destroy(sm);
  p7_sparsemx_Destroy(sxv);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxd);
  p7_refmx_Destroy(rxv);
  p7_refmx_Destroy(rxf);
  p7_refmx_Destroy(rxb);
  p7_refmx_Destroy(rxd);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}



/* The "reference-constrained" utest compares a randomly sampled
 * profile against homologous sequences (i.e. sampled from the same
 * profile) by sparse DP, using a sparse mask that includes all the
 * cells in the reference Viterbi trace plus a sprinkling of
 * additional random sparse supercells. Thus by construction, we know
 * that the sparse Viterbi score must equal the reference, and the
 * sparse Forward and Backward scores must be >= the reference Viterbi
 * score. 
 *
 * This test is important because it bounds the sparse DP scores from
 * *below*, which is nontrivial. If sparse DP missed a high scoring
 * path, it could elude other unit tests like 'randomseq' where
 * we only have the sparse DP score bounded from above by the 
 * reference DP.
 * 
 * Tests:
 * 1. Reference and sparse Viterbi scores are equal (within tolerance)
 * 2. Reference and sparse Viterbi traces are identical 
 * 3. Sparse F,B scores are >= reference Viterbi score, because we
 *    know the reference V path is included in the sparse mask.
 *
 * Also, more general test criteria apply, shared with 'randomseq' for 
 * example:
 * 
 * 4. Sparse F score == B
 * 5. All values v in reference, sparse V,F,B matrices satisfy v_r >= v_s
 * 6. Sparse V,F,B,D matrices Validate()
 */
static void
utest_reference_constrained(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char           msg[]  = "sparse fwdback, reference-constrained unit test failed";
  P7_HMM        *hmm    = NULL;
  P7_PROFILE    *gm     = p7_profile_Create(M, abc);
  ESL_SQ        *sq     = esl_sq_CreateDigital(abc);       /* space for generated (homologous) target seqs              */
  P7_REFMX      *rxv    = p7_refmx_Create(M, L);
  P7_REFMX      *rxf    = p7_refmx_Create(M, L);
  P7_REFMX      *rxb    = p7_refmx_Create(M, L);
  P7_TRACE      *rtr    = p7_trace_Create();
  P7_TRACE      *str    = p7_trace_Create();
  P7_SPARSEMASK *sm     = p7_sparsemask_Create(M, L);
  P7_SPARSEMX   *sxv    = p7_sparsemx_Create(sm);
  P7_SPARSEMX   *sxf    = p7_sparsemx_Create(sm);
  P7_SPARSEMX   *sxb    = p7_sparsemx_Create(sm);
  P7_SPARSEMX   *sxd    = p7_sparsemx_Create(sm);
  int            idx;
  float          vsc_s, fsc_s, bsc_s;
  float          vsc_r, fsc_r, bsc_r;
  float          tol = ( p7_logsum_IsSlowExact() ? 0.001 : 0.01 ); /* numerical error on V scores is a little higher in this test */

  /* Sample a profile. 
   * Config as usual: multihit dual-mode local/glocal, so all paths in it are valid.
   */
  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      /* Generate (sample) a sequence from the profile */
      if ( p7_profile_SetLength(gm, L)  != eslOK) esl_fatal(msg);   /* config to generate mean length of L (length was probably reset by last emitted seq) */
      do {
	esl_sq_Reuse(sq);
	p7_ProfileEmit(rng, hmm, gm, bg, sq, NULL);
      } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
      if ( p7_profile_SetLength(gm, sq->n) != eslOK) esl_fatal(msg);

      /* Reference DP calculations, including a reference Viterbi traceback */
      if ( p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,       &fsc_r) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceBackward(sq->dsq, sq->n, gm, rxb,       &bsc_r) != eslOK) esl_fatal(msg);
      if ( p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxv, rtr,  &vsc_r) != eslOK) esl_fatal(msg);
      
      /* Use the reference Viterbi trace to create a sparse mask */
      if ( p7_sparsemask_Reinit(sm, M, sq->n)  != eslOK) esl_fatal(msg);
      p7_sparsemask_SetFromTrace(sm, rng, rtr);

      /* Sparse DP calculations, in which we know the reference Viterbi trace will be scored */
      if ( p7_SparseViterbi (sq->dsq, sq->n, gm, sm, sxv, str, &vsc_s) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf,      &fsc_s) != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxb,      &bsc_s) != eslOK) esl_fatal(msg);
      if ( p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxb, sxd)        != eslOK) esl_fatal(msg);

      //p7_trace_DumpAnnotated(stdout, str, gm, sq->dsq);
      //p7_trace_DumpAnnotated(stdout, rtr, gm, sq->dsq);
      //p7_sparsemask_Dump(stdout, sm);
      //p7_sparsemx_Dump(stdout, sxv);
      //p7_refmx_Dump(stdout, rxv);

      /* Tests */
      if ( esl_FCompare_old(vsc_s, vsc_r, tol)                    != eslOK) esl_fatal(msg); /* (1) sparse, reference V scores equal */
      if ( p7_trace_CompareLoosely(str, rtr, sq->dsq)         != eslOK) esl_fatal(msg); /* (2) sparse, reference V traces identical; see notes on p7_trace_CompareLoosely */
      if ( ! (fsc_s - vsc_r + tol > 0.0f))                              esl_fatal(msg); /* (3) sparse F,B scores >= reference V score. */
      if ( esl_FCompare_old(fsc_s, bsc_s, tol)                    != eslOK) esl_fatal(msg); /* (4) sparse F score = B score */
      if ( p7_sparsemx_CompareReferenceAsBound(sxv, rxv, tol) != eslOK) esl_fatal(msg); /* (5) All V,F,B matrix values satisfy v_ref >= v_sparse */
      if ( p7_sparsemx_CompareReferenceAsBound(sxf, rxf, tol) != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_CompareReferenceAsBound(sxb, rxb, tol) != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Validate(sxv, NULL)                    != eslOK) esl_fatal(msg); /* (6): All sparse DP matrices Validate() (V,F,B,D) */
      if ( p7_sparsemx_Validate(sxf, NULL)                    != eslOK) esl_fatal(msg); /*      (the <NULL> arg is an optional <errbuf>)    */
      if ( p7_sparsemx_Validate(sxb, NULL)                    != eslOK) esl_fatal(msg);       
      if ( p7_sparsemx_Validate(sxd, NULL)                    != eslOK) esl_fatal(msg);

      /* reuse DP matrices, mask, and trace */
      if ( p7_sparsemask_Reuse(sm) != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxv)  != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxf)  != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxb)  != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reuse(sxd)  != eslOK) esl_fatal(msg);
      if ( p7_refmx_Reuse(rxv)     != eslOK) esl_fatal(msg);
      if ( p7_refmx_Reuse(rxf)     != eslOK) esl_fatal(msg);
      if ( p7_refmx_Reuse(rxb)     != eslOK) esl_fatal(msg);
      if ( p7_trace_Reuse(rtr)     != eslOK) esl_fatal(msg);
      if ( p7_trace_Reuse(str)     != eslOK) esl_fatal(msg);
    }

  p7_trace_Destroy(rtr);
  p7_trace_Destroy(str);
  p7_sparsemask_Destroy(sm);
  p7_sparsemx_Destroy(sxv);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxd);
  p7_refmx_Destroy(rxv);
  p7_refmx_Destroy(rxf);
  p7_refmx_Destroy(rxb);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}


/* The 'singlepath' utest samples a specially configured profile that
 * only allows one possible statepath. Generate a sequence from that
 * path, and set the sparse mask to include the path, plus marking a
 * dusting of additional randomly chosen sparse cells; then use sparse
 * DP routines to compare the profile to the seq. Test:
 * 
 * 1. Sparse V score = sparse F score = sparse B score
 * 2. Sparse V score = trace score of the generated seq
 * 3. Sparse Viterbi trace = generated trace
 * 
 * The key observation in this utest is that V=F=B scores, because we
 * know there is only a single possible path.
 * 
 * In order to guarantee a single path, the sampled profile must be
 * uniglocal L=0 (local mode would allow different local subseq
 * entry/exit; multihit or L>0 would allow different alignments of the
 * model to a target seq. Thus many paths are not exercised by the
 * test. Most importantly, because the profile aligns globally to the
 * sequence it generates, the sparse mask can't be multisegment; every
 * row i has at least one sparse cell marked.
 */
static void
utest_singlepath(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int N)
{
  char           msg[] = "sparse fwdback singlepath unit test failed";
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = p7_profile_Create(M, abc);
  ESL_SQ        *sq    = esl_sq_CreateDigital(abc);
  P7_TRACE      *gtr   = p7_trace_Create();           /* generated trace */
  P7_TRACE      *vtr   = p7_trace_Create();	      /* viterbi trace */
  P7_SPARSEMASK *sm    = p7_sparsemask_Create(M, M);  /* exact initial alloc size doesn't matter */
  P7_SPARSEMX   *sxv   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxf   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxb   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxd   = p7_sparsemx_Create(NULL);
  float          tsc, fsc, vsc, bsc;
  float          tol   = 1e-3;
  int            idx;

  for (idx = 0; idx < N; idx++)
    {
      /* Create a profile that has only a single possible path (including
       * emissions) thru it; requires configuring in uniglocal mode w/ L=0
       */
      if ( p7_modelsample_SinglePathed(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
      if ( p7_profile_ConfigUniglocal(gm, hmm, bg, 0)     != eslOK) esl_fatal(msg);

      /* Sample that sequence and path; get its trace score */
      if ( p7_ProfileEmit(rng, hmm, gm, bg, sq, gtr)    != eslOK) esl_fatal(msg);
      if ( p7_trace_Score(gtr, sq->dsq, gm, &tsc)       != eslOK) esl_fatal(msg);
  
      //p7_trace_DumpAnnotated(stdout, gtr, gm, sq->dsq);

      /* Build a randomized sparse mask around that trace */
      if ( p7_sparsemask_Reinit(sm, M, sq->n)  != eslOK) esl_fatal(msg); 
      p7_sparsemask_SetFromTrace(sm, rng, gtr);

      //p7_sparsemask_Dump(stdout, sm);

      /* Run DP routines, collect scores that should all match trace score */
      if ( p7_SparseViterbi  (sq->dsq, sq->n, gm, sm, sxv, vtr, &vsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseForward  (sq->dsq, sq->n, gm, sm, sxf,      &fsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseBackward (sq->dsq, sq->n, gm, sm, sxb,      &bsc) != eslOK) esl_fatal(msg);
      if ( p7_SparseDecoding (sq->dsq, sq->n, gm, sxf, sxb, sxd)      != eslOK) esl_fatal(msg);
  
      //p7_sparsemx_Dump(stdout, sxv);
      //p7_trace_DumpAnnotated(stdout, gtr, gm, sq->dsq);
      //p7_trace_DumpAnnotated(stdout, vtr, gm, sq->dsq);

      /* Since only a single path is possible, trace score and Fwd score match */
      if ( esl_FCompare(tsc, fsc, /*rtol=*/0.0, tol)   != eslOK) esl_fatal(msg); 
      if ( esl_FCompare(tsc, bsc, /*rtol=*/0.0, tol)   != eslOK) esl_fatal(msg);
      if ( esl_FCompare(tsc, vsc, /*rtol=*/0.0, tol)   != eslOK) esl_fatal(msg);
      if ( p7_trace_Compare(gtr, vtr, 0.0f)            != eslOK) esl_fatal(msg); // 0.0 is <pptol> arg, unused, because neither trace has PP annotation
  
      esl_sq_Reuse(sq);
      p7_trace_Reuse(vtr);
      p7_trace_Reuse(gtr);
      p7_sparsemask_Reuse(sm);
      p7_sparsemx_Reuse(sxf);
      p7_sparsemx_Reuse(sxv);
      p7_profile_Reuse(gm);
      p7_hmm_Destroy(hmm);
    }
  
  p7_sparsemx_Destroy(sxv);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemask_Destroy(sm);
  p7_trace_Destroy(vtr);
  p7_trace_Destroy(gtr);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
}


/* The 'internal-glocal-exit' utest tests for a peculiar quirk of
 * sparse DP: that MGk->...->E glocal exits may cross unmarked Dk+1.Dm
 * supercells, so any time we have an {DM}Gk sparse cell that does not
 * have DGk+1 marked, DP routines must allow an internal
 * wing-retracted {DM}Gk->E exit.
 * 
 * The utest creates a handcrafted simple example of comparing a 40aa
 * model to a 38aa truncated target, where the final transition needs
 * to cross an unmarked cell. It then tests the usual sort of stuff
 * against the reference:
 * 1. Viterbi score for sparse, reference are identical within tolerance
 * 2. Sparse F score = B
 * 3. Reference F score = B 
 * 4. Viterbi traces for sparse, reference are identical
 */
static void
utest_internal_glocal_exit(void)
{
  char          msg[]   = "sparse fwdback, internal-glocal-exit utest failed";
  ESL_ALPHABET *abc     = esl_alphabet_Create(eslAMINO);
  char          qseq[]  = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY";
  char          tseq[]  = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTV";
  int           M       = strlen(qseq);
  int           L       = strlen(tseq);
  ESL_DSQ      *qsq     = NULL;
  ESL_DSQ      *tsq     = NULL;
  float         popen   = 0.1;
  float         pextend = 0.4;
  P7_BUILDER   *bld     = p7_builder_Create(NULL, abc);
  P7_BG        *bg      = p7_bg_Create(abc);
  P7_HMM       *hmm     = NULL;
  P7_PROFILE   *gm      = p7_profile_Create(M, abc);
  P7_SPARSEMASK *sm     = p7_sparsemask_Create(M,L);
  P7_SPARSEMX   *sxv    = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxf    = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxb    = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxd    = p7_sparsemx_Create(NULL);
  P7_REFMX      *rxv    = p7_refmx_Create(M, L);
  P7_REFMX      *rxf    = p7_refmx_Create(M, L);
  P7_REFMX      *rxb    = p7_refmx_Create(M, L);
  P7_REFMX      *rxd    = p7_refmx_Create(M, L);
  P7_TRACE      *rtr    = p7_trace_Create();
  P7_TRACE      *str    = p7_trace_Create();
  float          vsc_s, fsc_s, bsc_s, tsc_s;
  float          vsc_r, fsc_r, bsc_r, tsc_r;
  int            i;
  float          tol = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  int            Q   = sm->Q; 

  /* Create the 40aa A-YA-Y test model */
  if ( esl_dsq_Create(abc, qseq, &qsq)                                                   != eslOK) esl_fatal(msg);
  if ( p7_builder_LoadScoreSystem(bld, "BLOSUM62", popen, pextend, bg)                   != eslOK) esl_fatal(msg); 
  if ( p7_Seqmodel(abc, qsq, M, "aatest", bld->Q, bg->f, bld->popen, bld->pextend, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_hmm_SetComposition(hmm)                                                        != eslOK) esl_fatal(msg);
  if ( p7_hmm_SetConsensus(hmm, NULL)                                                    != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)                                                    != eslOK) esl_fatal(msg);

  /* Create the 38aa truncated test target seq */
  if ( esl_dsq_Create(abc, tseq, &tsq)                                                   != eslOK) esl_fatal(msg);
  
  /* Create a sparse mask that includes the main diagonal and the last column. 
   * This mask exercises the potential bug.
   */
  for (i = L; i >= 1; i--)  // remember, sparsemask construction is constrained by backwards pass of the backwards filter
    {
      /* Add cells k=M (last column), k=i (diagonal) to each row.
       * API requires that they're added in backwards order (M, then i).
       * API also requires that they're added not as <k> coord but as <q,r> slot coords (see sparsemask.c docs)
       */
      if ( p7_sparsemask_StartRow(sm, i)                   != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_Add(sm, (M-1)%Q, (M-1)/Q) != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_Add(sm, (i-1)%Q, (i-1)/Q) != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_FinishRow(sm)                     != eslOK) esl_fatal(msg);
    }
  if ( p7_sparsemask_Finish(sm)          != eslOK) esl_fatal(msg);
  
  /* Reference DP calculations */
  if ( p7_ReferenceViterbi (tsq, L, gm, rxv, rtr, &vsc_r) != eslOK) esl_fatal(msg);
  if ( p7_ReferenceForward (tsq, L, gm, rxf,      &fsc_r) != eslOK) esl_fatal(msg);
  if ( p7_ReferenceBackward(tsq, L, gm, rxb,      &bsc_r) != eslOK) esl_fatal(msg);
  if ( p7_ReferenceDecoding(tsq, L, gm, rxf, rxb, rxd)    != eslOK) esl_fatal(msg);

  /* Sparse DP calculations */
  if ( p7_SparseViterbi (tsq, L, gm, sm, sxv, str, &vsc_s) != eslOK) esl_fatal(msg);
  if ( p7_SparseForward (tsq, L, gm, sm, sxf,      &fsc_s) != eslOK) esl_fatal(msg);
  if ( p7_SparseBackward(tsq, L, gm, sm, sxb,      &bsc_s) != eslOK) esl_fatal(msg);
  if ( p7_SparseDecoding(tsq, L, gm, sxf, sxb, sxd)        != eslOK) esl_fatal(msg);

  if ( p7_trace_Score(str, tsq, gm, &tsc_s) != eslOK) esl_fatal(msg);
  if ( p7_trace_Score(rtr, tsq, gm, &tsc_r) != eslOK) esl_fatal(msg);

  /* Tests */
  if ( esl_FCompare(vsc_s, vsc_r, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg); 
  if ( esl_FCompare(fsc_s, bsc_s, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg);
  if ( esl_FCompare(fsc_r, bsc_r, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg);
  if ( p7_trace_Compare(str, rtr, 0.0f)              != eslOK) esl_fatal(msg); // 0.0 is <pptol> arg, unused, because neither trace has PP annotation

  p7_trace_Destroy(rtr);
  p7_trace_Destroy(str);
  p7_sparsemask_Destroy(sm);
  p7_sparsemx_Destroy(sxv);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxd);
  p7_refmx_Destroy(rxv);
  p7_refmx_Destroy(rxf);
  p7_refmx_Destroy(rxb);
  p7_refmx_Destroy(rxd);
  p7_builder_Destroy(bld);
  p7_bg_Destroy(bg);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);  
  free(qsq);
  free(tsq);
}
  
#endif /*p7SPARSE_FWDBACK_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/


/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7SPARSE_FWDBACK_TESTDRIVE

#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,     "20", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for sparse DP of dual-mode profile";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg   = p7_bg_Create(abc);
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  utest_randomseq            (r, abc, bg, M, L, N);
  utest_compare_reference    (r, abc, bg, M, L, N);
  utest_reference_constrained(r, abc, bg, M, L, N);
  utest_singlepath           (r, abc, bg, M,    N);
  utest_internal_glocal_exit();

  fprintf(stderr, "#  status = ok\n");

  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}

#endif /*p7SPARSE_FWDBACK_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/



/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef p7SPARSE_FWDBACK_EXAMPLE

#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",              0 },
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "include all cells in sparse mx",                    0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                     0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Backward DP matrix for examination",           0 },
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump posterior decoding matrix for examination",    0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Forward DP matrix for examination",            0 },
  { "-M",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump sparse mask for examination",                  0 },
  { "-S",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump a stochastic trace for examination",           0 },
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Viterbi trace for examination",                0 },
  { "-V",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Viterbi DP matrix for examination",            0 },
  { "--diplot",  eslARG_OUTFILE,FALSE, NULL, NULL,   NULL,  NULL, NULL, "save domain inference plot to <f>",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of Forward/Backward, sparse dual implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
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
  P7_OPROFILE    *om      = NULL;
  P7_CHECKPTMX   *ox      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  P7_SPARSEMX    *sxv     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxf     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxb     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxd     = p7_sparsemx_Create(NULL);
  P7_TRACE       *tr      = p7_trace_CreateWithPP();
  char           *difile  = NULL;
  float           fsc, vsc, bsc;
  float           nullsc;
  char            errbuf[eslERRBUFSIZE];
  int             status;

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
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_oprofile_Convert(gm, om);

  /* Set the profile and null model's target length models */
  p7_bg_SetLength           (bg, sq->n);
  p7_profile_SetLength      (gm, sq->n);
  p7_oprofile_ReconfigLength(om, sq->n);

  /* Use f/b filter to create sparse mask */
  ox = p7_checkptmx_Create(hmm->M, sq->n, ESL_MBYTES(32));
  sm  = p7_sparsemask_Create(gm->M, sq->n);
  if (esl_opt_GetBoolean(go, "-a"))  
    p7_sparsemask_AddAll(sm);
  else {
    p7_ForwardFilter (sq->dsq, sq->n, om, ox, /*fsc=*/NULL);
    p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH);
  }
  
  /* Sparse DP calculations */
  p7_SparseViterbi (sq->dsq, sq->n, gm, sm, sxv, tr, &vsc);
  p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf,     &fsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxb,     &bsc);
  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxb, sxd);
  p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

  p7_sparsemx_TracePostprobs(sxd, tr);
  p7_trace_Index(tr);
  
  if (esl_opt_GetBoolean(go, "-B")) p7_sparsemx_Dump(stdout, sxb);
  if (esl_opt_GetBoolean(go, "-D")) p7_sparsemx_Dump(stdout, sxd);
  if (esl_opt_GetBoolean(go, "-F")) p7_sparsemx_Dump(stdout, sxf);
  if (esl_opt_GetBoolean(go, "-M")) p7_sparsemask_Dump(stdout, sm);
  if (esl_opt_GetBoolean(go, "-T")) p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);
  if (esl_opt_GetBoolean(go, "-V")) p7_sparsemx_Dump(stdout, sxv);

  if (( difile = esl_opt_GetString(go, "--diplot")) != NULL)
    {
      FILE *difp;
      if ((difp = fopen(difile, "w")) == NULL) esl_fatal("failed to open file %s\n", difile);
      p7_sparsemx_PlotDomainInference(difp, sxd, 1, sq->n, tr);
      fclose(difp);
    }

  if (esl_opt_GetBoolean(go, "-S")) {
    p7_trace_Reuse(tr);
    p7_sparse_trace_Stochastic(rng, NULL, gm, sxf, tr);
    if ( p7_trace_Validate(tr, abc, sq->dsq, errbuf) != eslOK) esl_fatal("stochastic trace validation failed: %s", errbuf);
    p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);
  }

  if ( p7_sparsemx_Validate( sxv, errbuf)  != eslOK) esl_fatal("sxv validation failed: %s", errbuf);
  if ( p7_sparsemx_Validate( sxf, errbuf)  != eslOK) esl_fatal("sxf validation failed: %s", errbuf);
  if ( p7_sparsemx_Validate( sxb, errbuf)  != eslOK) esl_fatal("sxb validation failed: %s", errbuf);
  if ( p7_sparsemx_Validate( sxd, errbuf)  != eslOK) esl_fatal("sxd validation failed: %s", errbuf);

  printf("target sequence:      %s\n",         sq->name);
  printf("vit raw score:        %.4f nats\n",  vsc);
  printf("fwd raw score:        %.4f nats\n",  fsc);
  printf("bck raw score:        %.4f nats\n",  bsc);
  printf("null score:           %.2f nats\n",  nullsc);
  printf("per-seq score:        %.2f bits\n",  (fsc - nullsc) / eslCONST_LOG2);

  /* Cleanup */
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_sparsemx_Destroy(sxv);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(ox);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_FWDBACK_EXAMPLE*/

