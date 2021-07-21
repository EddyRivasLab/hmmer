/* Sparse DP routines.
 *
 * Contents:
 *   1. Viterbi
 *   2. Forward
 *   3. Backward
 *   4. Posterior decoding
 *   5. Selection functions for Viterbi traceback
 *   6. Selection functions for stochastic traceback
 *   7. Traceback engine and API
 *   8. Unit tests
 *   9. Test driver
 *  10. Example
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_sparsemask.h"
#include "h4_sparsemx.h"

#include "logsum.h"
#include "sparse_dp.h"

/*****************************************************************
 * 1. Sparse Viterbi
 *****************************************************************/

/* Function:  h4_sparse_Viterbi()
 * Synopsis:  Viterbi algorithm, sparse DP.
 *
 * Purpose:   Compare profile <hmm> in mode <mo> to digital sequence
 *            <dsq> of length <L>, by the Viterbi algorithm, using
 *            sparse dynamic programming constrained by the sparse
 *            mask <sm>.  Fill in the sparse Viterbi matrix <sx>;
 *            (optionally) trace back the optimal path and return it
 *            in the trace structure <opt_pi> if the caller provides
 *            one; and (optionally) return the Viterbi raw score in
 *            bits in <*opt_vsc>.
 *            
 *            <sx> and <opt_pi> are provided as existing space;
 *            reused/reallocated as needed.
 *            
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of <dsq>
 *            hmm     - profile
 *            sm      - sparse mask
 *            sx      - Viterbi matrix to fill (reallocated/reused)
 *            opt_pi  - optRESULT: optimal path (reallocated/reused); or NULL if unwanted
 *            opt_vsc - optRETURN: raw Viterbi score in bits; or NULL
 *
 * Returns:   <eslOK> on success; <sx> is the Viterbi DP matrix;
 *            <opt_pi>, if non-NULL, contains Viterbi path; and
 *            <*opt_vsc> optionally contains the raw Viterbi score.
 *
 * Throws:    <eslEMEM> if reallocation of <sx> or <opt_pi> fails.
 */
int
h4_sparse_Viterbi(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMASK *sm, H4_SPARSEMX *sx, H4_PATH *opt_pi, float *opt_vsc)
{
  float         *xpc; 	                 // ptr that steps through current special cells   
  float         *dpc;	                 // ptr to step thru current row i main DP cells 
  float         *dpp;			 // ptr to step thru previous row i-1 main DP cells
  float         *last_dpc;		 // used to reinit dpp after each sparse row computation
  int            ng;                     // counts residues outside sparsified segments 
  float          xE, xN, xJ, xB, xL, xG, xC;  // tmp scores on special states. only stored when in row bands, and on ia-1 before a seg
  float          mlc, mgc;		 // temporary score calculations M(i,k)
  float          ilc, igc;               //   ... and I(i,k)
  float          dlc, dgc;		 // precalculated D(i,k+1) value on current row
  int           *kc = sm->k[0];		 // <kc> points to the list of sparse cell indices k for current row i
  int           *kp;			 // <kp> points to the previous row's sparse cell index list
  int            i,k;	      	         // i,k row,col (seq position, profile position) cell coords
  int            y,z;			 // indices in lists of k coords on prev, current row
  int            status;

  /* Argument validation */
  ESL_DASSERT1( (sm->L == L) );
  ESL_DASSERT1( (sm->M == hmm->M) );

  /* Assure that <sx> is allocated large enough (we might be reusing it).
   * Set its type now, so we can Dump/Validate/etc. during debugging this routine, if needed.
   */
  if ( (status = h4_sparsemx_Reinit(sx, sm)    ) != eslOK) return status;
  sx->type = h4S_VITERBI;

  xN  = 0.0f;
  xJ  = -eslINFINITY;
  xC  = -eslINFINITY;
  ng  = 0;
  xpc = sx->xmx;
  dpc = sx->dp;
  for (i = 1; i <= L; i++)
    {
      if (! sm->n[i]) { ng++; continue; }   // skip rows that have no included cells

      /* Reinitialize and store specials for row ia-1 just outside a sparsified segment */
      if (i == 1 || ng) {
	*xpc++ = xE = -eslINFINITY;
	*xpc++ = xN  = xN + ( ng ? ng * mo->xsc[h4_N][h4_LOOP] : 0.0); // test ng, because we must watch out for 0*-inf special case 
	*xpc++ = xJ  = xJ + ( ng ? ng * mo->xsc[h4_J][h4_LOOP] : 0.0);
	*xpc++ = xB  = ESL_MAX( xN + mo->xsc[h4_N][h4_MOVE], xJ + mo->xsc[h4_J][h4_MOVE]);
	*xpc++ = xL  = xB + mo->xsc[h4_B][0]; // B->L 
	*xpc++ = xG  = xB + mo->xsc[h4_B][1]; // B->G 
	*xpc++ = xC  = xC + ( ng ? ng * mo->xsc[h4_C][h4_LOOP] : 0.0);
	*xpc++       = -eslINFINITY; // JJ: this space only used in a Decoding matrix.
	*xpc++       = -eslINFINITY; // CC: this space only used in a Decoding matrix.
	ng = 0;
      }

      last_dpc = dpc;		// remember where dpc started; dpp will be set here after we finish each row calculation 
      kp = kc;                  // last row we did becomes prev row now; ready to step through k indices of previous row's sparse cells 
      kc = sm->k[i];		// ditto for current row i 
      dlc = dgc = xE = -eslINFINITY;
      for (z=0, y=0; z < sm->n[i]; z++) // Iterate over the one or more sparse cells (i,k) that we calculate on this row. 
	{
	  k = kc[z]; // index of next sparse cell to calculate: (i,k) 
	  
	  /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	  mlc = xL  + hmm->tsc[k-1][h4_LM];
	  mgc = xG  + hmm->tsc[k-1][h4_GM];
	  while (y < sm->n[i-1] && kp[y]  < k-1) { y++; dpp+=h4S_NSCELLS; }
	  if    (y < sm->n[i-1] && kp[y] == k-1) {
	    mlc = ESL_MAX( ESL_MAX( dpp[h4S_ML] + hmm->tsc[k-1][h4_MM],
				    dpp[h4S_IL] + hmm->tsc[k-1][h4_IM]),
			   ESL_MAX( dpp[h4S_DL] + hmm->tsc[k-1][h4_DM],
				    mlc));        
	    mgc = ESL_MAX( ESL_MAX( dpp[h4S_MG] + hmm->tsc[k-1][h4_MM],
				    dpp[h4S_IG] + hmm->tsc[k-1][h4_IM]),
			   ESL_MAX( dpp[h4S_DG] + hmm->tsc[k-1][h4_DM],
				    mgc));
	  }
	  *dpc++ = mlc = hmm->rsc[dsq[i]][k] + mlc;
	  *dpc++ = mgc = hmm->rsc[dsq[i]][k] + mgc;

	  /* Try to find cell i-1,k; then compute I(i,k) from it */
	  while (y < sm->n[i-1] && kp[y] < k)  { y++; dpp+=h4S_NSCELLS; }
	  if    (y < sm->n[i-1] && kp[y] == k) {
	    *dpc++ = ilc = ESL_MAX( ESL_MAX( dpp[h4S_ML] + hmm->tsc[k][h4_MI],
                                             dpp[h4S_IL] + hmm->tsc[k][h4_II]),
                                             dpp[h4S_DL] + hmm->tsc[k][h4_DI]);  // +ISC(k), if we weren't enforcing it to zero
	    *dpc++ = igc = ESL_MAX( ESL_MAX( dpp[h4S_MG] + hmm->tsc[k][h4_MI],
                                             dpp[h4S_IG] + hmm->tsc[k][h4_II]),
                                    ESL_MAX( dpp[h4S_DG] + hmm->tsc[k][h4_DI],
                                             xG          + hmm->tsc[k][h4_GI])); // ditto
	  } else {
	    *dpc++ = ilc = -eslINFINITY;
	    *dpc++ = igc = xG + hmm->tsc[k][h4_GI];
	  }
	    
	  /* local exit paths, with implicit transition prob 1.0.
           * In Plan9 models, there are edge cases where Ik->Dk+1->E is an optimal path.
           */
	  xE = ESL_MAX(xE, ESL_MAX(mlc, dlc));

	  /* delayed store of Dk; advance calculation of next D_k+1 */
	  *dpc++ = dlc;
	  *dpc++ = dgc;
          if (z < sm->n[i]-1 && kc[z+1] == k+1) { // is there a (i,k+1) cell to our RIGHT? (not a typo) 
	    dlc = ESL_MAX( ESL_MAX( mlc + hmm->tsc[k][h4_MD],
                                    ilc + hmm->tsc[k][h4_ID]),
                                    dlc + hmm->tsc[k][h4_DD]);
	    dgc = ESL_MAX( ESL_MAX( mgc + hmm->tsc[k][h4_MD],
                                    igc + hmm->tsc[k][h4_ID]),
                                    dgc + hmm->tsc[k][h4_DD]);
	  } else if (k==hmm->M) { // no? is it because we're k==M?
            xE = ESL_MAX(xE, ESL_MAX(mgc, dgc));
            dlc = dgc = -eslINFINITY;
          } else { // if not, then evaluate {MD}Gk->Dk+1..E glocal exit path from internal sparse cells with no k+1 neighbor
	    xE  = ESL_MAX( xE,
                           hmm->tsc[k][h4_DGE] + ESL_MAX( ESL_MAX( mgc + hmm->tsc[k][h4_MD],
                                                                   igc + hmm->tsc[k][h4_ID]),
                                                                   dgc + hmm->tsc[k][h4_DD]));  
	    dlc = dgc = -eslINFINITY;
	  }
	}

      *xpc++ = xE;  // we already max'ed over all Mk->E exits, both local and glocal
      *xpc++ = xN = xN + mo->xsc[h4_N][h4_LOOP];
      *xpc++ = xJ = ESL_MAX( xJ + mo->xsc[h4_J][h4_LOOP],  xE + mo->xsc[h4_E][h4_LOOP]);
      *xpc++ = xB = ESL_MAX( xJ + mo->xsc[h4_J][h4_MOVE],  xN + mo->xsc[h4_N][h4_MOVE]);
      *xpc++ = xL = xB + mo->xsc[h4_B][0]; // B->L 
      *xpc++ = xG = xB + mo->xsc[h4_B][1]; // B->G 
      *xpc++ = xC = ESL_MAX( xE + mo->xsc[h4_E][h4_MOVE],  xC + mo->xsc[h4_C][h4_LOOP]);
      *xpc++      = -eslINFINITY; // JJ: this space only used in a Decoding matrix 
      *xpc++      = -eslINFINITY; // CC: this space only used in a Decoding matrix 

      /* now dpc is on the start of the next sparsified row */
      dpp = last_dpc;
    }

  if (opt_vsc) *opt_vsc = xC + ( ng ? ng * mo->xsc[h4_C][h4_LOOP] : 0.0f) + mo->xsc[h4_C][h4_MOVE];
  if (opt_pi && xC != -eslINFINITY) return h4_sparse_ViterbiTrace(hmm, mo, sx, opt_pi);
  else                              return eslOK;
}
/*-------------------- end, Viterbi -----------------------------*/



/*****************************************************************
 * 2. Sparse Forward
 *****************************************************************/

/* Function:  h4_sparse_Forward()
 * Synopsis:  Forward algorithm, sparse DP
 *
 * Purpose:   Compare profile <hmm> in mode <mo> to digital sequence
 *            <dsq> of length <L>, by the sparse Forward DP algorithm,
 *            constrained by the sparse mask <sm>.  Fill in the sparse
 *            Forward matrix <sx>, and (optionally) return the Forward
 *            raw score in bits in <*opt_fsc>.
 *            
 *            <sx> is provided as existing space; reused/reallocated
 *            as needed.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of <dsq>
 *            hmm     - profile
 *            mo      - comparison mode
 *            sm      - sparse mask
 *            sx      - Forward matrix to fill (reallocated/reused)
 *            opt_fsc - optRETURN: raw Forward score in bits; or NULL
 *
 * Returns:   <eslOK> on success. 
 *
 * Throws:    <eslEMEM> if reallocation of <sx> fails.
 */
int
h4_sparse_Forward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMASK *sm, H4_SPARSEMX *sx, float *opt_fsc)
{
  float         *xpc; 	                 // ptr that steps through current special cells   
  float         *dpc;	                 // ptr to step thru current row i main DP cells 
  float         *dpp;			 // ptr to step thru previous row i-1 main DP cells
  float         *last_dpc;		 // used to reinit dpp after each sparse row computation
  int            ng;                     // counts residues outside sparsified segments 
  float          xE, xN, xJ, xB, xL, xG, xC;  // tmp scores on special states. only stored when in row bands, and on ia-1 before a seg
  float          mlc, mgc;		 // temporary score calculations M(i,k)
  float          ilc, igc;               //   ... and I(i,k)
  float          dlc, dgc;		 // precalculated D(i,k+1) value on current row
  int           *kc = sm->k[0];		 // <kc> points to the list of sparse cell indices k for current row i
  int           *kp;			 // <kp> points to the previous row's sparse cell index list
  int            i,k;	      	         // i,k row,col (seq position, profile position) cell coords
  int            y,z;			 // indices in lists of k coords on prev, current row
  int            status;

  /* Argument validation */
  ESL_DASSERT1( (sm->L == L) );
  ESL_DASSERT1( (sm->M == hmm->M) );

  /* Assure that <sx> is allocated large enough (we might be reusing it).
   * Set its type now, so we can Dump/Validate/etc. during debugging this routine, if needed.
   */
  if ( (status = h4_sparsemx_Reinit(sx, sm)) != eslOK) return status;
  sx->type = h4S_FORWARD;

  xN  = 0.0f;
  xJ  = -eslINFINITY;
  xC  = -eslINFINITY;
  ng  = 0;
  xpc = sx->xmx;
  dpc = sx->dp;
  for (i = 1; i <= L; i++)
    {
      if (! sm->n[i]) { ng++; continue; }   // skip rows that have no included cells

      /* Reinitialize and store specials for row ia-1 just outside a sparsified segment */
      if (i == 1 || ng) {
	*xpc++ = xE = -eslINFINITY;
	*xpc++ = xN  = xN + ( ng ? ng * mo->xsc[h4_N][h4_LOOP] : 0.0); // test ng, because we must watch out for 0*-inf special case 
	*xpc++ = xJ  = xJ + ( ng ? ng * mo->xsc[h4_J][h4_LOOP] : 0.0);
	*xpc++ = xB  = h4_logsum( xN + mo->xsc[h4_N][h4_MOVE], xJ + mo->xsc[h4_J][h4_MOVE]);
	*xpc++ = xL  = xB + mo->xsc[h4_B][0]; // B->L 
	*xpc++ = xG  = xB + mo->xsc[h4_B][1]; // B->G 
	*xpc++ = xC  = xC + ( ng ? ng * mo->xsc[h4_C][h4_LOOP] : 0.0);
	*xpc++       = -eslINFINITY; // JJ: this space only used in a Decoding matrix.
	*xpc++       = -eslINFINITY; // CC: this space only used in a Decoding matrix.
	ng = 0;
      }

      last_dpc = dpc;		// remember where dpc started; dpp will be set here after we finish each row calculation 
      kp = kc;                  // last row we did becomes prev row now; ready to step through k indices of previous row's sparse cells 
      kc = sm->k[i];		// ditto for current row i 
      dlc = dgc = xE = -eslINFINITY;
      for (z=0, y=0; z < sm->n[i]; z++) // Iterate over the one or more sparse cells (i,k) that we calculate on this row. 
	{
	  k = kc[z]; // index of next sparse cell to calculate: (i,k) 
	  
	  /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	  mlc = xL  + hmm->tsc[k-1][h4_LM];
	  mgc = xG  + hmm->tsc[k-1][h4_GM];
	  while (y < sm->n[i-1] && kp[y]  < k-1) { y++; dpp+=h4S_NSCELLS; }
	  if    (y < sm->n[i-1] && kp[y] == k-1) {
	    mlc = h4_logsum( h4_logsum( dpp[h4S_ML] + hmm->tsc[k-1][h4_MM],
				        dpp[h4S_IL] + hmm->tsc[k-1][h4_IM]),
			     h4_logsum( dpp[h4S_DL] + hmm->tsc[k-1][h4_DM],
				        mlc));        
	    mgc = h4_logsum( h4_logsum( dpp[h4S_MG] + hmm->tsc[k-1][h4_MM],
				        dpp[h4S_IG] + hmm->tsc[k-1][h4_IM]),
			     h4_logsum( dpp[h4S_DG] + hmm->tsc[k-1][h4_DM],
				        mgc));
	  }
	  *dpc++ = mlc = hmm->rsc[dsq[i]][k] + mlc;
	  *dpc++ = mgc = hmm->rsc[dsq[i]][k] + mgc;

	  /* Try to find cell i-1,k; then compute I(i,k) from it */
	  while (y < sm->n[i-1] && kp[y] < k)  { y++; dpp+=h4S_NSCELLS; }
	  if    (y < sm->n[i-1] && kp[y] == k) {
	    *dpc++ = ilc = h4_logsum( h4_logsum( dpp[h4S_ML] + hmm->tsc[k][h4_MI],
                                                 dpp[h4S_IL] + hmm->tsc[k][h4_II]),
                                                 dpp[h4S_DL] + hmm->tsc[k][h4_DI]);  
	    *dpc++ = igc = h4_logsum( h4_logsum( dpp[h4S_MG] + hmm->tsc[k][h4_MI],
                                                 dpp[h4S_IG] + hmm->tsc[k][h4_II]),
                                      h4_logsum( dpp[h4S_DG] + hmm->tsc[k][h4_DI],
                                                 xG          + hmm->tsc[k][h4_GI])); 
	  } else {
	    *dpc++ = ilc = -eslINFINITY;
	    *dpc++ = igc = xG + hmm->tsc[k][h4_GI];
	  }
	    
	  /* local exit paths, with implicit transition prob 1.0.
           * In Plan9 models, there are edge cases where Ik->Dk+1->E is an optimal path.
           */
	  xE = h4_logsum(xE, h4_logsum(mlc, dlc));

	  /* delayed store of Dk; advance calculation of next D_k+1 */
	  *dpc++ = dlc;
	  *dpc++ = dgc;
          if (z < sm->n[i]-1 && kc[z+1] == k+1) { // is there a (i,k+1) cell to our RIGHT? (not a typo) 
	    dlc = h4_logsum( h4_logsum( mlc + hmm->tsc[k][h4_MD],
                                        ilc + hmm->tsc[k][h4_ID]),
                                        dlc + hmm->tsc[k][h4_DD]);
	    dgc = h4_logsum( h4_logsum( mgc + hmm->tsc[k][h4_MD],
                                        igc + hmm->tsc[k][h4_ID]),
                                        dgc + hmm->tsc[k][h4_DD]);
	  } else if (k==hmm->M) { // no? is it because we're k==M?
            xE = h4_logsum(xE, h4_logsum(mgc, dgc));
            dlc = dgc = -eslINFINITY;
          } else { // if not, then evaluate {MD}Gk->Dk+1..E glocal exit path from internal sparse cells
	    xE  = h4_logsum(xE,
                            hmm->tsc[k][h4_DGE] + h4_logsum( h4_logsum( mgc + hmm->tsc[k][h4_MD],
                                                                        igc + hmm->tsc[k][h4_ID]),
                                                                        dgc + hmm->tsc[k][h4_DD])); 
	    dlc = dgc = -eslINFINITY;
	  }
	}

      *xpc++ = xE;  // we already max'ed over all Mk->E exits, both local and glocal
      *xpc++ = xN = xN + mo->xsc[h4_N][h4_LOOP];
      *xpc++ = xJ = h4_logsum( xJ + mo->xsc[h4_J][h4_LOOP],  xE + mo->xsc[h4_E][h4_LOOP]);
      *xpc++ = xB = h4_logsum( xJ + mo->xsc[h4_J][h4_MOVE],  xN + mo->xsc[h4_N][h4_MOVE]);
      *xpc++ = xL = xB + mo->xsc[h4_B][0]; // B->L 
      *xpc++ = xG = xB + mo->xsc[h4_B][1]; // B->G 
      *xpc++ = xC = h4_logsum( xE + mo->xsc[h4_E][h4_MOVE],  xC + mo->xsc[h4_C][h4_LOOP]);
      *xpc++      = -eslINFINITY; // JJ: this space only used in a Decoding matrix 
      *xpc++      = -eslINFINITY; // CC: this space only used in a Decoding matrix 

      /* now dpc is on the start of the next sparsified row */
      dpp = last_dpc;
    }

  if (opt_fsc) *opt_fsc = xC + ( ng ? ng * mo->xsc[h4_C][h4_LOOP] : 0.0f) + mo->xsc[h4_C][h4_MOVE];
  return eslOK;
}
/*-------------------- end, Forward -----------------------------*/

/*****************************************************************
 * 3. Sparse Backward
 *****************************************************************/

/* Function:  h4_sparse_Backward()
 * Synopsis:  Backward algorithm, sparse DP.
 *
 * Purpose:   Compare profile <hmm> in mode <mo> to digital sequence
 *            <dsq> of length <L> by the sparse Backward DP algorithm,
 *            constrained by the sparse mask <sm>.  Fill in the sparse
 *            DP Backward matrix <sx> and (optionally) return the
 *            overall raw Backward score in <*opt_sc>.
 *
 * Args:      dsq     - digital target sequence 1..L
 *            L       - length of <dsq>
 *            hmm     - profile
 *            mo      - comparison mode
 *            sm      - sparse mask
 *            sx      - Backward matrix to fill (reallocated/reused)
 *            opt_sc  - optRETURN: raw Backwards score in bits
 *                      
 * Returns:   <eslOK> on success. 
 *            
 * Throws:    <eslEMEM> if reallocation of <sx> fails
 */
int
h4_sparse_Backward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMASK *sm, H4_SPARSEMX *sx, float *opt_sc)
{
  float          *xc;            // ptr that steps through current special cells 
  float          *dpc;		 // ptr stepping through current DP row i sparse cells 
  float          *dpn = NULL;	 // ptr stepping through next DP row (i+1) sparse cells; NULL when no i+1 row stored 
  float          *last_dpc;	 // remembers where dpc started, so dpn can be set for next i 
  int    i,k,y,z;
  float  xC,xG,xL,xB,xJ,xN,xE;	 // tmp vars for special cells on this row i
  float  dlc,dgc;		 // tmp vars for D vals on current row i  
  float  mln,mgn,iln,ign;	 // tmp vars for M+e, I+e on next row i+1 
  int    status;

  /* Argument validation */
  ESL_DASSERT1( (sm->L == L) );
  ESL_DASSERT1( (sm->M == hmm->M) );

  /* Assure <sx> is allocated large enough (we might be reusing it).
   * Set its type now, so we can Dump/Validate/etc. during debugging this routine, if needed.
   */
  if ( (status = h4_sparsemx_Reinit(sx, sm)) != eslOK) return status;
  sx->type = h4S_BACKWARD;

  /* In Backwards, we traverse DP matrix backwards; init ptrs to last elements */
  xc  = sx->xmx + (sm->nrow + sm->S - 1)*h4S_NXCELLS; // last supercell in xmx 
  dpc = sx->dp  + (sm->ncells-1)*h4S_NSCELLS;	      // last supercell in dp  

  /* xC/N/J follow a convention that as we start a row i, they already include
   * all contribution from row i+1; i.e. they are calculated as a lookahead, 
   * then immediately stored.
   */
  xC = mo->xsc[h4_C][h4_MOVE];    // xC[L] can only be reached via C->T transition 
  xJ = xN = -eslINFINITY;	  // J,N unreached on L, and indeed until at least first ib (last stored row)
  xG = xL = -eslINFINITY;	  // L,G can only be reached from a stored next row i+1.

  for (i = sm->L; i >= 1; i--)
    {
      /* if this row isn't stored, then it can only extend CC/JJ/NN:
       * extend them to prev row (remember precalculation convention)
       */
      if (! sm->n[i]) 
	{ 
	  xC += mo->xsc[h4_C][h4_LOOP];
	  xJ += mo->xsc[h4_J][h4_LOOP];
	  xN += mo->xsc[h4_N][h4_LOOP];
	  continue;
	}

      /* Else we continue: this row i is stored, both as a main row (in dpc) and specials (in xmx) */
      xc[h4S_CC] = -eslINFINITY;     // CC only stored in a decoding matrix 
      xc[h4S_JJ] = -eslINFINITY;     // ditto JJ 
      xc[h4S_C]  = xC;               // deferred store of xC. 
      xc[h4S_G]  = xG;               // on i=ib segment start, xG cannot be reached; it will be -inf
      xc[h4S_L]  = xL;               // ditto xL
      xc[h4S_B]  = xB = h4_logsum( xL + mo->xsc[h4_B][0],       xG + mo->xsc[h4_B][1]);        // on i=ib segment start, evaluates to -inf
      xc[h4S_J]  = xJ = h4_logsum( xJ, 		                xB + mo->xsc[h4_J][h4_MOVE]);  // on i=ib, evaluates to xJ
      xc[h4S_N]  = xN = h4_logsum( xN, 	                        xB + mo->xsc[h4_N][h4_MOVE]);  // ditto xN
      xc[h4S_E]  = xE = h4_logsum( xJ + mo->xsc[h4_E][h4_LOOP], xC + mo->xsc[h4_E][h4_MOVE]);  // on i=ib, i<L, E->{JC} are both possible; i=L, evaluates to xC + tEC.
      xc -= h4S_NXCELLS;

      last_dpc = dpc;                        // last_dpc will remember where dpc was (at end of row i); before we loop back to new i, dpn will be set to last_dpc.
      xG  = xL = -eslINFINITY;               // prepare to accumulate xG, xL logsums
      y = (i == sm->L ? -1 : sm->n[i+1]-1);  // if n[i+1] were allocated and set to a 0 sentinel, we could avoid this branch
      for (z = sm->n[i]-1; z >= 0; z--)
	{
	  k = sm->k[i][z];
          if (k == hmm->M)  // k=M is special
            {
              dpc[h4S_ML] = dpc[h4S_MG] = xE;
              dpc[h4S_IL] = dpc[h4S_IG] = -eslINFINITY;
              dpc[h4S_DL] = dlc = dpc[h4S_DG] = dgc = xE;
              xG          = dpc[h4S_MG] + hmm->rsc[dsq[i]][k] + hmm->tsc[k-1][h4_GM];
              xL          = dpc[h4S_ML] + hmm->rsc[dsq[i]][k] + hmm->tsc[k-1][h4_LM];
            }
          else
            {
              /* try to pick up mln, mgn from i+1,k+1 (and add their emission score) */
              while (y >= 0 && sm->k[i+1][y]  > k+1) { y--; dpn -= h4S_NSCELLS; } // note if y were already on i+1,k+1, it doesn't move
              if    (y >= 0 && sm->k[i+1][y] == k+1) {
                mln = dpn[h4S_ML] + hmm->rsc[dsq[i+1]][k+1];
                mgn = dpn[h4S_MG] + hmm->rsc[dsq[i+1]][k+1];
              } else { mln = mgn = -eslINFINITY; }
	  
              /* try to pick up iln,ign from i+1,k */
              while (y >= 0 && sm->k[i+1][y]  > k) { y--; dpn -= h4S_NSCELLS; } // note if y were already on i+1,k+1, it doesn't move
              if    (y >= 0 && sm->k[i+1][y] == k) {
                iln = dpn[h4S_IL]; 
                ign = dpn[h4S_IG]; 
              } else { iln = ign = -eslINFINITY; }

              /* see if dlc, dgc are valid from prev calculation on this row. if not reinit. note glocal path, wing-retracted. */
              if (z == sm->n[i]-1 || sm->k[i][z+1] != k+1) { 
                dlc = -eslINFINITY; 
                dgc = xE + hmm->tsc[k][h4_DGE];
              }
	  
              /* M(i,k) calculations need to use dgc,dlc before we
               * change them, while they're still D(i,k+1). 
               */
              dpc[h4S_ML] = h4_logsum( h4_logsum(mln + hmm->tsc[k][h4_MM],      // ML(i,k) =   ML(i+1,k+1) * t(k,MM) * e_M(k+1, x_i+1)      | mln = log[ ML(i+1,k+1) * e_M(k+1,x_i+1)]
                                                 iln + hmm->tsc[k][h4_MI]),     //           + IL(i+1,k)   * t(k,MI)                        | iln = log[ IL(i+1,k)   * e_I(k,  x_i+1)]
                                       h4_logsum(dlc + hmm->tsc[k][h4_MD],      //           + DL(i,  k+1) * t(k,MD)                        | dlc = DL(i,k+1), wrapped around from prev loop iteration. This has tricky boundary conditions at k=kbc, and k=kbc=M
                                                 xE));                          //           +  E(i)       * t(MkE)=1.0
              dpc[h4S_MG] = h4_logsum( h4_logsum(mgn + hmm->tsc[k][h4_MM],      // MG(i,k) is essentially the same recursion, without a transition to E.
                                                 ign + hmm->tsc[k][h4_MI]),     
                                                 dgc + hmm->tsc[k][h4_MD]);     
	      
              /* I(i,k) calculations */
              dpc[h4S_IL] = h4_logsum( h4_logsum( mln + hmm->tsc[k][h4_IM],
                                                  iln + hmm->tsc[k][h4_II]),
                                                  dlc + hmm->tsc[k][h4_ID]);
              dpc[h4S_IG] = h4_logsum( h4_logsum( mgn + hmm->tsc[k][h4_IM],
                                                  ign + hmm->tsc[k][h4_II]),
                                                  dgc + hmm->tsc[k][h4_ID]);

              /* Accumulate xG, xL as we sweep over the current row; 
               * they get used to initialize the next row (i-1). Note
               * how we need emission scores on current row for this,
               * versus next row for the mgn/mln calculations above.
               */
              xG = h4_logsum( h4_logsum( xG,
                                         dpc[h4S_MG] + hmm->rsc[dsq[i]][k] + hmm->tsc[k-1][h4_GM]), 
                                         dpc[h4S_IG] + hmm->tsc[k][h4_GI]);
              xL = h4_logsum( xL, dpc[h4S_ML] + hmm->rsc[dsq[i]][k] + hmm->tsc[k-1][h4_LM]);   
	  
              /* D(i,k) calculations; 
               * we can start storing results in <dpc> and decrementing;
               * and we're done with dgc,dlc, so can update them for
               * time thru the k loop.
               */
              dpc[h4S_DL] = dlc = h4_logsum( h4_logsum( mln  + hmm->tsc[k][h4_DM],
                                                        iln  + hmm->tsc[k][h4_DI]),
                                             h4_logsum( dlc  + hmm->tsc[k][h4_DD],
                                                        xE));                   
              dpc[h4S_DG] = dgc = h4_logsum( h4_logsum( mgn  + hmm->tsc[k][h4_DM],
                                                        ign  + hmm->tsc[k][h4_DI]),
                                                        dgc  + hmm->tsc[k][h4_DD]);
            } // end k<M cases

	  dpc -= h4S_NSCELLS;
	} // end loop over sparse cells k on this row.
      dpn = last_dpc;

      /* precalculate what xC/xJ/xN will be on previous row i-1... these values get stored as we roll around i loop */
      xC += mo->xsc[h4_C][h4_LOOP];
      xJ += mo->xsc[h4_J][h4_LOOP];
      xN += mo->xsc[h4_N][h4_LOOP];

      /* If we are terminating a segment (if this i is an ia), we need to 
       * store specials for row ia-1. This works for i=1 and storing the
       * final i=0 special row too.
       */
      if (sm->n[i-1] == 0) 
	{
	  xc[h4S_CC] = -eslINFINITY;	/* CC only stored in a Decoding matrix. */
	  xc[h4S_JJ] = -eslINFINITY;    /* JJ only stored in a Decoding matrix. */
	  xc[h4S_C]  = xC;
	  xc[h4S_G]  = xG;
	  xc[h4S_L]  = xL;
	  xc[h4S_B]  = xB = h4_logsum( xL + mo->xsc[h4_B][0],       xG + mo->xsc[h4_B][1]);
	  xc[h4S_J]  = xJ = h4_logsum( xJ,     	                    xB + mo->xsc[h4_J][h4_MOVE]);
	  xc[h4S_N]  = xN = h4_logsum( xN,		            xB + mo->xsc[h4_N][h4_MOVE]);
	  xc[h4S_E]  = xE = h4_logsum( xJ + mo->xsc[h4_E][h4_LOOP], xC + mo->xsc[h4_E][h4_MOVE]);
	  xc -= h4S_NXCELLS;

	  xG = xL = -eslINFINITY; /* when we start the next segment at an ib, these will get stored */
	}
    } // end loop over i

  if (opt_sc) *opt_sc = xN;	/* tS->N is 1.0, no cost. */
  return eslOK;
}
/*-------------------- end, Backward ----------------------------*/


/*****************************************************************
 * 4. Sparse posterior decoding
 *****************************************************************/

/* Function:  h4_sparse_Decoding()
 * Synopsis:  Posterior decoding algorithm, in sparse DP.
 *
 * Purpose:   Given Forward matrix <sxf> and Backward matrix <sxb> that
 *            have been computed for a comparison of profile <hmm> in 
 *            mode <mo> against digital sequence <dsq> of length <L>;
 *            fill in the posterior decoding matrix <sxd>.
 *            
 *            <sxb> and <sxd> can point to the same structure, in
 *            which case posterior decoding overwrites the Backward
 *            matrix; a trick that the caller might use to save
 *            allocating a new matrix. (Can't do the same with
 *            Forward, because of a detail involving CC/JJ transition
 *            decoding.)
 * 
 *            If <sxd> is an independent matrix (i.e. not overwriting
 *            <sxb>), it will be reinitialized (and possibly reallocated)
 *            to use the same sparse mask as <sxf> and <sxb>. 
 *            
 *            {M/I}(i,k) is the prob that state {MI}k generated
 *            residue <x_i>.  D(i,k) is the prob that a state path
 *            used state Dk after already generating residue <x_i>
 *            with a previous M or I state. {BLG}(i) is the
 *            probability of being in a {BLG} state just as we start a
 *            new domain at <x_i+1>. E(i) is the probability of being
 *            in the end state on row <i>, having ended the domain on
 *            this row.  These all arise from standard posterior
 *            decoding eqn's.
 *            
 *            Watch out for {NJC}, though, the states that emit on
 *            transition; these are decoded in a slightly nonstandard
 *            way. {NN/CC/JJ}(i) is the probability that we emitted
 *            residue <i> on an {NN/CC/JJ} transition. {NCJ}(i) is the
 *            probability that the state path uses {NCJ} on (i), which
 *            is a sum of emitting <i> and an initial mute transition
 *            {S->N,E->C,E->J}.
 *
 * Args:      dsq  - digital sequence, 1..L
 *            L    - length of <dsq>
 *            hmm  - profile
 *            mo   - comparison mode
 *            sxf  - Forward matrix, computed by caller
 *            sxb  - Backward matrix, computed by caller
 *            sxd  - Decoding matrix to fill (realloc/reused as needed)
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if a reallocation of <sxd> fails.
 */
int
h4_sparse_Decoding(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMX *sxf, const H4_SPARSEMX *sxb, H4_SPARSEMX *sxd)
{
  const H4_SPARSEMASK *sm = sxf->sm;
  const float   *dpf      = sxf->dp;
  float         *dpb      = sxb->dp;
  float         *dpd;
  const float   *xf       = sxf->xmx;
  float         *xb       = sxb->xmx;
  float         *xd;
  float         *dpd2;
  const float   *rsc;
  float          totsc  = xf[h4S_N] + xb[h4S_N];  // first ia-1 is not necessarily row 0, so xb[N] alone does not suffice!
  float          norm;                // except for numerical roundoff error accumulation, so we also explicitly renormalize each row.
  float          xC,xJ,xG;            // Forward values remembered from prev row i-1.
  int            i,k,y,z;
  float          delta;		      // additions to DGk's, resulting from unfolding the wing-retracted entry/exit paths */
#ifndef h4SPARSE_DECODING_TESTDRIVE   // see note on renormalization further below. 
  int            x;
#endif
  int            status;

  /* Argument validation */
  ESL_DASSERT1 ( (sxf->type == h4S_FORWARD) );
  ESL_DASSERT1 ( (sxb->type == h4S_BACKWARD) );
  ESL_DASSERT1 ( (sxf->sm == sxb->sm) );
  ESL_DASSERT1 ( (sm->L   == L) );
  ESL_DASSERT1 ( (sm->M   == hmm->M) );

  /* Assure that <sxd> is allocated large enough (we might be reusing it),
   * either because we're overwriting <sxb> (which we know is large enough),
   * or because we reinitialize it (if it isn't already using <sm>). 
   */
  if (sxd != sxb && (status = h4_sparsemx_Reinit(sxd, sm)) != eslOK) return status;
  sxd->type = h4S_DECODING;

  xJ = xC = -eslINFINITY;
  dpd = sxd->dp;
  xd  = sxd->xmx;
  for (i = 1; i <= L; i++)
    {
      rsc = hmm->rsc[dsq[i]];	

      /* i=ia-1, initialization of a sparse segment; special storage */
      if (sm->n[i] && !sm->n[i-1])
	{
	  norm       = 0.0f;
	  xd[h4S_E]  = 0.0f;                                
	  xd[h4S_N]  = exp2f(xf[h4S_N] + xb[h4S_N] - totsc);                  norm += xd[h4S_N];
	  xd[h4S_J]  = exp2f(xf[h4S_J] + xb[h4S_J] - totsc);  
	  xd[h4S_B]  = exp2f(xf[h4S_B] + xb[h4S_B] - totsc);
	  xd[h4S_L]  = exp2f(xf[h4S_L] + xb[h4S_L] - totsc);
	  xd[h4S_G]  = exp2f(xf[h4S_G] + xb[h4S_G] - totsc);  xG = xf[h4S_G];
	  xd[h4S_C]  = exp2f(xf[h4S_C] + xb[h4S_C] - totsc);
	  xd[h4S_JJ] = xd[h4S_J];                            xJ = xf[h4S_J]; norm += xd[h4S_JJ];
	  xd[h4S_CC] = xd[h4S_C];                            xC = xf[h4S_C]; norm += xd[h4S_CC];
	  
#ifndef h4SPARSE_DECODING_TESTDRIVE /* see note on renormalization further below. */
	  norm = 1.0f/norm;
	  for (x = 0; x < h4S_NXCELLS; x++) xd[x] *= norm;
#endif
	  xf += h4S_NXCELLS;
	  xb += h4S_NXCELLS;
	  xd += h4S_NXCELLS;
	}

      /* For each sparse cell z (k=k[i][z]) on row: 
       * DG is special, because we have to unfold the wing-retracted entries and exits. 
       *
       * First, exits: unfolding right wing retractions 
       *
       * Caution: Although we're unfolding wings correctly, by design,
       * we're only storing that probability mass in i,k cells in the
       * sparse mask. The sparse mask was created by local-only
       * decoding in ForwardFilter(). Suppose there's a high
       * probability G->DDDD->Mk entry, for example. The D's won't
       * necessarily be in the mask! Thus, sparse decoding does NOT
       * give you a complete decoding of D's. I don't think we care;
       * but I'm leaving a bread crumb here, just in case.
       */
      dpd2 = dpd;
      for (delta = 0.0f, z = 0; z < sm->n[i]; z++)
	{
	  k = sm->k[i][z];
	  dpd[h4S_DG] = delta + exp2f(dpf[h4S_DG] + dpb[h4S_DG] - totsc);   // because we might be overwriting backwards mx with decoding, this is the last time we can access dpb[h4S_DG] values on the row! 
	  if (z == sm->n[i]-1 || sm->k[i][z+1] != k+1) // No cell to our right? then {MD}Gk->Dk+1..Dm->E path delta added to all subsequent stored Dk+1..m */
	    delta += exp2f( xb[h4S_E] + hmm->tsc[k][h4_DGE] - totsc + h4_logsum( h4_logsum( dpf[h4S_MG] + hmm->tsc[k][h4_MD],
                                                                                            dpf[h4S_IG] + hmm->tsc[k][h4_ID]),
                                                                                            dpf[h4S_DG] + hmm->tsc[k][h4_DD]));
          dpd += h4S_NSCELLS; dpf += h4S_NSCELLS; dpb += h4S_NSCELLS;
	}
      /* Second, entries: unfolding left wing retractions. The DG's for a G->D..D->Mk path are on the PREVIOUS ROW from the Mk! */
      /* All Mk on current row contributes a delta; and each delta applies to all 1..k-1 on prev row. Sparsity on both rows makes this tricky */
      /* and remember, we need the residue score for e(x_i,MGk) as well as the backwards score MGk,i */
      for (delta = 0.0f, y = sm->n[i-1]-1, z = sm->n[i]-1; z >= 0; z--)
	{
	  k      = sm->k[i][z];
	  dpb   -= h4S_NSCELLS;
          delta += exp2f(xG + hmm->tsc[k][h4_GI]           + dpb[h4S_IG] - totsc); // G->D1..Dk->Ik entry path
	  while (y >= 0 && sm->k[i-1][y] >= k) { dpd2 -= h4S_NSCELLS; dpd2[h4S_DG] += delta; y--; }
	  delta += exp2f(xG + hmm->tsc[k-1][h4_GM] + rsc[k] + dpb[h4S_MG] - totsc); // G->D1..Dk-1->Mk entry path, added to all stored D1..Dk-1 
	}
      while (y >= 0) { dpd2 -= h4S_NSCELLS; dpd2[h4S_DG] += delta; y--; } /* dpd2 now sits on first sparse cell of prev row i-1. */
      /* note that dpb ran up and back down; dpf and dpd only ran up, and need to be run back down */
      dpf -= sm->n[i]*h4S_NSCELLS;
      dpd -= sm->n[i]*h4S_NSCELLS;

      norm = 0.0;
      for (z = 0; z < sm->n[i]; z++, dpd += h4S_NSCELLS, dpf += h4S_NSCELLS, dpb += h4S_NSCELLS)
	{
	  dpd[h4S_ML] = exp2f(dpf[h4S_ML] + dpb[h4S_ML] - totsc); norm += dpd[h4S_ML];
	  dpd[h4S_MG] = exp2f(dpf[h4S_MG] + dpb[h4S_MG] - totsc); norm += dpd[h4S_MG];
	  dpd[h4S_IL] = exp2f(dpf[h4S_IL] + dpb[h4S_IL] - totsc); norm += dpd[h4S_IL];
	  dpd[h4S_IG] = exp2f(dpf[h4S_IG] + dpb[h4S_IG] - totsc); norm += dpd[h4S_IG];
	  dpd[h4S_DL] = exp2f(dpf[h4S_DL] + dpb[h4S_DL] - totsc);  // nonemitters don't count toward normalization
	  // DG was already done above.
	}
      
      /* specials on each stored row */
      if (sm->n[i]) {
	xd[h4S_JJ] = exp2f(  xJ      + xb[h4S_J] + mo->xsc[h4_J][h4_LOOP] - totsc); xJ = xf[h4S_J]; norm += xd[h4S_JJ];  // JJ,CC calculations must come before J,C; they depend on xb[J,C], which we may overwrite when we calc J,C
	xd[h4S_CC] = exp2f(  xC      + xb[h4S_C] + mo->xsc[h4_C][h4_LOOP] - totsc); xC = xf[h4S_C]; norm += xd[h4S_CC];
	xd[h4S_E]  = exp2f(xf[h4S_E] + xb[h4S_E] - totsc);                                  
	xd[h4S_N]  = exp2f(xf[h4S_N] + xb[h4S_N] - totsc);                                          norm += xd[h4S_N];
	xd[h4S_J]  = exp2f(xf[h4S_J] + xb[h4S_J] - totsc);  
	xd[h4S_B]  = exp2f(xf[h4S_B] + xb[h4S_B] - totsc);                                  
	xd[h4S_L]  = exp2f(xf[h4S_L] + xb[h4S_L] - totsc);                                  
	xd[h4S_G]  = exp2f(xf[h4S_G] + xb[h4S_G] - totsc);                          xG = xf[h4S_G];                              
	xd[h4S_C]  = exp2f(xf[h4S_C] + xb[h4S_C] - totsc);                   

	/* Renormalization.
	 * 
	 * Roundoff error accumulation in F/B is significant. For large
	 * target seqs, it isn't unusual to have a whole nat of
	 * difference in overall fwd vs. bck raw score; for example, fn3
	 * vs. TITIN_HUMAN. Since we only use the bck score (at i=0) for
	 * normalization above, pp's would have very large systematic
	 * error. 
	 * 
	 * To squash this error in production code, we renormalize rows.
	 * Because renormalization can hide real errors, we don't
	 * renormalize when we've compiled the code for unit testing.
	 * Default unit tests don't run large enough M/L to create a
	 * lot of error accumulation.
	 *
	 * (A similar note appears in reference_decoding.)
	 */
#ifndef h4SPARSE_DECODING_TESTDRIVE	
	norm = 1.0f/norm;
	dpd -= sm->n[i]*h4S_NSCELLS;  /* back up to start of row again */
	for (x = 0; x < sm->n[i]*h4S_NSCELLS; x++) *dpd++ *= norm;  
	for (x = 0; x < h4S_NXCELLS;          x++) xd[x]  *= norm;
#endif
	  
	xd += h4S_NXCELLS;
	xf += h4S_NXCELLS;
	xb += h4S_NXCELLS;
      }
    }
  return eslOK;
}
/*--------------- end, posterior decoding -----------------------*/



/*****************************************************************
 *  5. Selection functions for Viterbi traceback
 *****************************************************************/

static inline int
v_select_ml(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np, const float *xp, int *ret_z)
{
  ESL_UNUSED(r);
  int   y = 0;
  while (y < np && kp[y]  < k-1) { y++; dpp += h4S_NSCELLS; } 
  if    (y < np && kp[y] == k-1) 
    {
      float path[4];

      path[0] =   xp[h4S_L] + hmm->tsc[k-1][h4_LM];
      path[1] = dpp[h4S_ML] + hmm->tsc[k-1][h4_MM];
      path[2] = dpp[h4S_IL] + hmm->tsc[k-1][h4_IM];
      path[3] = dpp[h4S_DL] + hmm->tsc[k-1][h4_DM];

      *ret_z = y;
      return (h4P_L + esl_vec_FArgMax(path, 4));  // assumes L-ML-IL-DL order in h4P_*
    }
  else { *ret_z = 0; return h4P_L; }
}
static inline int
v_select_mg(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np, const float *xp, int *ret_z)
{
  ESL_UNUSED(r);
  int   y = 0;
  while (y < np && kp[y]  < k-1) { y++; dpp += h4S_NSCELLS; } 
  if    (y < np && kp[y] == k-1) 
    {
      float path[4];

      path[0] =   xp[h4S_G] + hmm->tsc[k-1][h4_GM];
      path[1] = dpp[h4S_MG] + hmm->tsc[k-1][h4_MM];
      path[2] = dpp[h4S_IG] + hmm->tsc[k-1][h4_IM];
      path[3] = dpp[h4S_DG] + hmm->tsc[k-1][h4_DM];

      *ret_z = y;
      return (h4P_G + esl_vec_FArgMax(path, 4));  // assumes G-MG-IG-DG order in h4P_*
    }
  else { *ret_z = 0; return h4P_G; }
}
static inline int
v_select_il(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np, int *ret_z)
{
  ESL_UNUSED(r);
  float path[3];
  int   y = 0;
  while (kp[y] != k) { y++; dpp += h4S_NSCELLS; } // a little brave; we know an appropriate sparse cell exists on prv row, else we couldn't reach I on cur 
  path[0] = dpp[h4S_ML] + hmm->tsc[k][h4_MI];
  path[1] = dpp[h4S_IL] + hmm->tsc[k][h4_II];
  path[2] = dpp[h4S_DL] + hmm->tsc[k][h4_DI];
  *ret_z = y;
  return (h4P_ML + esl_vec_FArgMax(path, 3));    // assumes ML-IL-DL order in h4P_*
}
static inline int
v_select_ig(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np, const float *xp, int *ret_z)
{
  ESL_UNUSED(r);
  int   y = 0;
  while (y < np && kp[y]  < k) { y++; dpp += h4S_NSCELLS; } 
  if    (y < np && kp[y] == k)  // because of G->Ik wing-retracted entry, we need to check whether sparse cell k exists above us, unlike plan7
    {
      float path[4];

      path[0] =  xp[h4S_G]  + hmm->tsc[k][h4_GI];
      path[1] = dpp[h4S_MG] + hmm->tsc[k][h4_MI];
      path[2] = dpp[h4S_IG] + hmm->tsc[k][h4_II];
      path[3] = dpp[h4S_DG] + hmm->tsc[k][h4_DI];

      *ret_z = y;
      return (h4P_G + esl_vec_FArgMax(path, 4));  // assumes G-MG-IG-DG order in h4P_*
    }
  else { *ret_z = 0; return h4P_G; }
}
static inline int
v_select_dl(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, int k, const float *dpp)
{
  ESL_UNUSED(r);
  float path[3];
  path[0] = dpp[h4S_ML] + hmm->tsc[k-1][h4_MD]; // more bravery. fact that we're tracing back from DL means that sparse cell k-1 must exist in z-1 
  path[1] = dpp[h4S_IL] + hmm->tsc[k-1][h4_ID];
  path[2] = dpp[h4S_DL] + hmm->tsc[k-1][h4_DD];
  return (h4P_ML + esl_vec_FArgMax(path, 3));   // assumes ML-IL-DL order in h4P_*
}
static inline int
v_select_dg(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, int k, const float *dpp)
{
  ESL_UNUSED(r);
  float path[3];
  path[0] = dpp[h4S_MG] + hmm->tsc[k-1][h4_MD]; // more bravery. fact that we're tracing back from DL means that sparse cell k-1 must exist in z-1 
  path[1] = dpp[h4S_IG] + hmm->tsc[k-1][h4_ID];
  path[2] = dpp[h4S_DG] + hmm->tsc[k-1][h4_DD];
  return (h4P_MG + esl_vec_FArgMax(path, 3));   // assumes MG-IG-DG order in h4P_*
}
static inline int
v_select_j(ESL_RANDOMNESS *r, const H4_MODE *mo, const float *xc)
{
  ESL_UNUSED(r);
  return ( (*(xc-h4S_NXCELLS+h4S_J) + mo->xsc[h4_J][h4_LOOP] >  // i.e. xp[J] on prv row i-1. 
                          xc[h4S_E] + mo->xsc[h4_E][h4_LOOP]) ? h4P_J : h4P_E);
}
static inline int
v_select_c(ESL_RANDOMNESS *r, const H4_MODE *mo, const float *xc)
{
  ESL_UNUSED(r);
  return ( (*(xc-h4S_NXCELLS+h4S_C) + mo->xsc[h4_C][h4_LOOP] >  // i.e. xp[C] on prv row i-1. 
                          xc[h4S_E] + mo->xsc[h4_E][h4_MOVE]) ? h4P_C : h4P_E);
}
static inline int
v_select_e(ESL_RANDOMNESS *r, float *wrk, const H4_PROFILE *hmm, const float *dpp, const int *kp, int np, int *ret_z)
{
  ESL_UNUSED(r);
  ESL_UNUSED(wrk);
  float max  = -eslINFINITY;
  float pathsc;
  int   smax = -1;
  int   zmax = -1;
  int   k,z;

  for (z = 0; z < np; z++)
    { 
      k = kp[z];
      if (dpp[h4S_ML] >= max) { max = dpp[h4S_ML]; smax = h4P_ML; zmax = z; }   
      if (dpp[h4S_DL] >= max) { max = dpp[h4S_DL]; smax = h4P_DL; zmax = z; }  // In Plan9, ML->IL->DL->E can be a Viterbi path. In Plan7, DL never ends a path.
      
      if (k == hmm->M) // final M cell doesn't use DGE wing folding, and it has no I.
        {
          if (dpp[h4S_MG] > max) { max = dpp[h4S_MG]; smax = h4P_MG; zmax = z; }
          if (dpp[h4S_DG] > max) { max = dpp[h4S_DG]; smax = h4P_DG; zmax = z; }
        }
      else if (z == np-1 || kp[z+1] != k+1) // sparse cell k with no following cell k+1: check glocal exit Mk->...->E across unmarked cells
	{ 
	  pathsc = dpp[h4S_MG] + hmm->tsc[k][h4_MD] + hmm->tsc[k][h4_DGE]; if (pathsc > max) { max = pathsc; smax = h4P_MG; zmax = z; }
          pathsc = dpp[h4S_IG] + hmm->tsc[k][h4_ID] + hmm->tsc[k][h4_DGE]; if (pathsc > max) { max = pathsc; smax = h4P_IG; zmax = z; }
          pathsc = dpp[h4S_DG] + hmm->tsc[k][h4_DD] + hmm->tsc[k][h4_DGE]; if (pathsc > max) { max = pathsc; smax = h4P_DG; zmax = z; }
	}
      dpp += h4S_NSCELLS;
    }
  *ret_z = zmax;
  return smax;
}
static inline int
v_select_b(ESL_RANDOMNESS *r, const H4_MODE *mo, const float *xc)
{
  ESL_UNUSED(r);
  return ( (xc[h4S_J] + mo->xsc[h4_J][h4_MOVE] > 
            xc[h4S_N] + mo->xsc[h4_N][h4_MOVE]) ? h4P_J : h4P_N);
}
/*----------- end, Viterbi traceback selections -----------------*/


/***************************************************************** 
 * 6. Selection functions for stochastic traceback
 *****************************************************************/

/* TK TK */
/*---------- end, stochastic traceback selections ---------------*/



/***************************************************************** 
 * 7. Traceback engine and API
 *****************************************************************/

static int
sparse_traceback_engine(ESL_RANDOMNESS *rng, float *wrk, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMX *sx, H4_PATH *pi)
{
  const H4_SPARSEMASK *sm = sx->sm;
  const float *dpc;             // points at start of row i
  const float *dpp;             // points at start of row i-1
  const float *xc;              // points to current specials, row i
  const float *xp;              // points to prv specials, row i-1
  int            k  = 0;	// current coord in profile consensus
  int            i;	        // current coord in sequence (that dpc and xc are on)
  int            bump_i;        // TRUE to do i-- 
  int            scur, sprv;
  int            z;		// current position in n[i] sparse k array entries on current row dp[]
  int            g;             // what segment we're in, 1..sm->S
  int            status;

  if ((status = h4_path_Reuse(pi)) != eslOK) return status;

  int (*select_ml)(ESL_RANDOMNESS *rng,             const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np, const float *xp, int *ret_z);
  int (*select_mg)(ESL_RANDOMNESS *rng,             const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np, const float *xp, int *ret_z);
  int (*select_il)(ESL_RANDOMNESS *rng,             const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np,                  int *ret_z);
  int (*select_ig)(ESL_RANDOMNESS *rng,             const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np, const float *xp, int *ret_z);
  int (*select_dl)(ESL_RANDOMNESS *rng,             const H4_PROFILE *hmm, int k, const float *dpp);
  int (*select_dg)(ESL_RANDOMNESS *rng,             const H4_PROFILE *hmm, int k, const float *dpp);
  int (*select_e) (ESL_RANDOMNESS *rng, float *wrk, const H4_PROFILE *hmm,        const float *dpp, const int *kp, int np,                  int *ret_z);
  int (*select_j) (ESL_RANDOMNESS *rng,             const H4_MODE *mo, const float *xc);
  int (*select_c) (ESL_RANDOMNESS *rng,             const H4_MODE *mo, const float *xc);
  int (*select_b) (ESL_RANDOMNESS *rng,             const H4_MODE *mo, const float *xc);
  
  if (sx->type == h4S_VITERBI)	// configure for Viterbi optimal traceback 
    {
      select_ml = &v_select_ml;      select_mg = &v_select_mg;
      select_il = &v_select_il;      select_ig = &v_select_ig;
      select_dl = &v_select_dl;      select_dg = &v_select_dg;
      select_j  = &v_select_j;       select_c  = &v_select_c;
      select_e  = &v_select_e;       select_b  = &v_select_b;
    }
#if 0 // TK TK
  else if (sx->type == h4S_FORWARD) // configure for stochastic traceback 
    {
      select_ml = &sto_select_ml;      select_mg = &sto_select_mg;
      select_il = &sto_select_il;      select_ig = &sto_select_ig;
      select_dl = &sto_select_dl;      select_dg = &sto_select_dg;
      select_j  = &sto_select_j;       select_c  = &sto_select_c;
      select_e  = &sto_select_e;       select_b  = &sto_select_b;
    }
#endif
  else ESL_EXCEPTION(eslEINCONCEIVABLE, "neither Forward or Viterbi?");


  /* Initialization to ib[S] row; ib[S]+1..L are assigned to C */
  i    = sm->seg[sm->S].ib;  
  scur = h4P_C;
  dpc  = sx->dp  + (sm->ncells - sm->n[i]) * h4S_NSCELLS;   // to *start* of row ib[S], so we can index with z
  xc   = sx->xmx + (sm->nrow+sm->S-1) * h4S_NXCELLS;        // to last xc supercell, which is ib[S]
  if ((status = h4_path_AppendElement(pi, h4P_C, sm->L-i+1)) != eslOK) return status;  // L-i residues, +1 for emit-on-transition

  for (g = sm->S; g >= 1; g--)
    {
      while ( i >= sm->seg[g].ia || ! h4_path_IsX(scur))
        {
          bump_i = FALSE;
          if (i > 0) {
            dpp    = dpc - (sm->n[i-1] * h4S_NSCELLS);
            xp     = xc - h4S_NXCELLS;
          }

          switch (scur) {
          case h4P_ML: sprv = (*select_ml)(rng, hmm, k, dpp, sm->k[i-1], sm->n[i-1], xp, &z); bump_i = TRUE; k--;      break;
          case h4P_MG: sprv = (*select_mg)(rng, hmm, k, dpp, sm->k[i-1], sm->n[i-1], xp, &z); bump_i = TRUE; k--;      break;
          case h4P_IL: sprv = (*select_il)(rng, hmm, k, dpp, sm->k[i-1], sm->n[i-1],     &z); bump_i = TRUE;           break;
          case h4P_IG: sprv = (*select_ig)(rng, hmm, k, dpp, sm->k[i-1], sm->n[i-1], xp, &z); bump_i = TRUE;           break;
          case h4P_DL: sprv = (*select_dl)(rng, hmm, k, dpc + (z-1)*h4S_NSCELLS);                            k--; z--; break; 
          case h4P_DG: sprv = (*select_dg)(rng, hmm, k, dpc + (z-1)*h4S_NSCELLS);                            k--; z--; break;
          case h4P_N:  sprv = h4P_N;                                                                                   break;
          case h4P_J:  sprv = (*select_j)(rng, mo, xc);                                                                break; 
          case h4P_C:  sprv = (*select_c)(rng, mo, xc);                                                                break; 
          case h4P_E:  sprv = (*select_e)(rng, wrk, hmm, dpc, sm->k[i], sm->n[i], &z);                  k=sm->k[i][z]; break;
          case h4P_B:  sprv = (*select_b)(rng, mo, xc);                                                                break; 
          case h4P_L:  sprv = h4P_B;                                                                                   break;
          case h4P_G:  sprv = h4P_B;                                                                                   break;
          default:     ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
          }
          if (h4_path_IsX(sprv) && scur == sprv) bump_i = TRUE;
      
          /* Glocal B->G->{MI}k left wing retraction entry: unfold it */
          if (sprv == h4P_G && k > 0) {
            if ((status = h4_path_AppendElement(pi, h4P_DG, k)) != eslOK) return status;
            k=0;
          }
          /* Glocal {MID}Gk->E right wing retraction: off last sparse cell k, {MID}k->Dk+1->E */
          if (scur == h4P_E && k < hmm->M && (sprv == h4P_MG || sprv == h4P_IG || sprv == h4P_DG)) {
            if ((status = h4_path_AppendElement(pi, h4P_DG, hmm->M-k)) != eslOK) return status;
          }
      
          if (sprv == h4P_L) status = h4_path_AppendElement(pi, h4P_L, k+1); // if we just backed into an L state, record the k of L->Mk in its rle[].
          else               status = h4_path_Append       (pi, sprv);       // note that this appends G if that's our scur
          if (status != eslOK) return status;

          if (bump_i)
            {
              i--;
              dpc = dpp;
              xc  = xp;
            } 
          scur = sprv;
        }

      // Now we're on ia[g]-1, in a N|J|C state, and we haven't emitted ia[g]-1 itself yet.
      // Residues ib[g-1]+1..ia[g]-1 are between segments, emitted by the same N|J|C state.
      if (sprv != h4P_N)
        {
          if ((status = h4_path_AppendSeveral(pi, scur, sm->seg[g].ia - sm->seg[g-1].ib - 1)) != eslOK) return status; 
          i     = sm->seg[g-1].ib;
          dpc  -= sm->n[i] * h4S_NSCELLS;   
          dpp  = dpc - (sm->n[i-1] * h4S_NSCELLS);
          xc  -= h4S_NXCELLS;              
          xp   = xc - h4S_NXCELLS;
        }
    }
  // now i=0 if we've accounted for all seq; else, on last residue that N needs to account for.
  if ((status = h4_path_AppendSeveral(pi, scur, i)) != eslOK) return status;  
  return h4_path_Reverse(pi);
}

/* Function:  h4_sparse_ViterbiTrace()
 * Synopsis:  Traceback of a sparse Viterbi DP matrix.
 *
 * Purpose:   Caller has filled sparse Viterbi matrix <sx> with
 *            <h4_sparse_Viterbi()>, in a comparison involving profile
 *            <hmm> in comparison mode <mo>. Perform a Viterbi
 *            traceback, storing the result in <pi>, which the caller
 *            provides.
 *
 * Args:      hmm - profile
 *            mo  - comparison mode
 *            sx  - filled sparse Viterbi DP matrix
 *            pi  - path structure to store the result in (reused/reallocated here as needed)
 *
 * Returns:   <eslOK> on success, and <pi> contains the optimal path.
 *
 * Throws:    <eslEMEM> on allocation failure; <pi> may need to grow to
 *            accommodate the trace.
 * 
 *            <eslEINVAL> on various coding, validation failures. 
 */
int
h4_sparse_ViterbiTrace(const H4_PROFILE *hmm, const H4_MODE *mo, const H4_SPARSEMX *sx, H4_PATH *pi)
{
  return sparse_traceback_engine(NULL, NULL, hmm, mo, sx, pi);
}

/*------------- end, traceback engine and API -------------------*/


/***************************************************************** 
 * 8. Unit tests
 *****************************************************************/
#ifdef h4SPARSE_DP_TESTDRIVE

#include "esl_sq.h"

#include "h4_checkptmx.h"
#include "h4_refmx.h"

#include "emit.h"
#include "fbfilter.h"
#include "modelsample.h"
#include "reference_dp.h"


/* utest_generation()
 *
 * The "generation" unit test emits a sequence from a randomly sampled
 * profile HMM, then compares them by sparse DP, including a fwd/bck
 * filter step to generate the sparse mask. It tests that:
 *     
 *  1. Fwd, Bck score >= Vit
 *  2. Fwd = Bck score
 *  3. F,B,V,D matrices pass Validate()
 *  4. Score of Vit path = Vit score
 *
 *  It also runs reference DP for comparison, and tests:
 *  
 *  5. Reference F,B,V scores >= sparse F,B,V score 
 *  6. All values v in reference, sparse F,B,V matrices 
 *     satisfy v_ref >= v_sparse
 */
static void
utest_generation(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, int M, int L, int N)
{
  char           msg[]  = "sparse_dp:: generation unit test failed";
  H4_PROFILE    *hmm    = NULL;
  H4_MODE       *mo     = h4_mode_Create();
  ESL_SQ        *sq     = esl_sq_CreateDigital(abc);
  H4_CHECKPTMX  *cpx    = h4_checkptmx_Create(M, L, ESL_MBYTES(32));
  H4_SPARSEMASK *sm     = h4_sparsemask_Create(M, L);
  H4_SPARSEMX   *sxv    = h4_sparsemx_Create(NULL);
  H4_SPARSEMX   *sxf    = h4_sparsemx_Create(NULL);
  H4_SPARSEMX   *sxb    = h4_sparsemx_Create(NULL);
  H4_SPARSEMX   *sxd    = h4_sparsemx_Create(NULL);
  H4_REFMX      *rxv    = h4_refmx_Create(M, L);
  H4_REFMX      *rxf    = h4_refmx_Create(M, L);
  H4_REFMX      *rxb    = h4_refmx_Create(M, L);
  H4_REFMX      *rxd    = h4_refmx_Create(M, L);
  H4_PATH       *epi    = h4_path_Create();
  H4_PATH       *vpi    = h4_path_Create();
  float          esc;                 // score of emitted seq path
  float          fsc_s, bsc_s, vsc_s; // sparse DP scores for fwd, bck, vit
  float          fsc_r, bsc_r, vsc_r; // reference DP scores for fwd, bck, vit
  float          psc_s;               // score of viterbi path
  int            idx;
  float          tol  = ( h4_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  char           errbuf[eslERRBUFSIZE];
   
  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  /* choose a random alignment mode: */
  switch (esl_rnd_Roll(rng, 6)) {
  case 0: h4_mode_SetDefault(mo);   break;
  case 1: h4_mode_SetLocal(mo);     break;
  case 2: h4_mode_SetGlocal(mo);    break;
  case 3: h4_mode_SetUnihit(mo);    break;
  case 4: h4_mode_SetUnilocal(mo);  break;
  case 5: h4_mode_SetUniglocal(mo); break;
  }

  for (idx = 0; idx < N; idx++)
    {
      /* Generate a sequence */
      if ( h4_mode_SetLength(mo, L)                   != eslOK) esl_fatal(msg);                              // set length model to L for emission
      do { if ( h4_emit(rng, hmm, mo, sq, epi)        != eslOK) esl_fatal(msg); } while (sq->n > (M+L)*3);   // keep seq length reasonable
      if ( h4_mode_SetLength(mo, sq->n)               != eslOK) esl_fatal(msg);   // reset length model to this sampled seq's length before scoring or DP 
      if ( h4_path_Score(epi, sq->dsq, hmm, mo, &esc) != eslOK) esl_fatal(msg);   // emitted path score. Viterbi can improve on this.

      /* Fwd/Bck local filter to calculate the sparse mask */
      if ( h4_fwdfilter(sq->dsq, sq->n, hmm, mo, cpx, NULL)                  != eslOK) esl_fatal(msg);
      if ( h4_bckfilter(sq->dsq, sq->n, hmm, mo, cpx, sm, h4SPARSIFY_THRESH) != eslOK) esl_fatal(msg);
  
      /* Sparse DP calculations */
      if ( h4_sparse_Viterbi (sq->dsq, sq->n, hmm, mo, sm, sxv, vpi, &vsc_s) != eslOK) esl_fatal(msg);
      if ( h4_sparse_Forward (sq->dsq, sq->n, hmm, mo, sm, sxf,      &fsc_s) != eslOK) esl_fatal(msg);
      if ( h4_sparse_Backward(sq->dsq, sq->n, hmm, mo, sm, sxb,      &bsc_s) != eslOK) esl_fatal(msg);
      if ( h4_sparse_Decoding(sq->dsq, sq->n, hmm, mo, sxf, sxb, sxd)        != eslOK) esl_fatal(msg);
      if ( h4_path_Score(vpi, sq->dsq, hmm, mo, &psc_s)                      != eslOK) esl_fatal(msg);

      /* Reference DP calculations */
      if ( h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxv, NULL, &vsc_r)  != eslOK) esl_fatal(msg);
      if ( h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,       &fsc_r)  != eslOK) esl_fatal(msg);
      if ( h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxb,       &bsc_r)  != eslOK) esl_fatal(msg);
      if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxb, rxd)      != eslOK) esl_fatal(msg);
      
      /* Tests */  
      if ( fsc_s < vsc_s+tol)                             esl_fatal(msg);                      // 1. Fwd score >= Vit score 
      if ( esl_FCompare(fsc_s, bsc_s, 0.0, tol) != eslOK) esl_fatal(msg);                      // 2. Fwd score == Bck score 
      if ( h4_sparsemx_Validate(sxv, errbuf)    != eslOK) esl_fatal("%s:\n %s", msg, errbuf);  // 3. F,B,V,D matrices Validate() 
      if ( h4_sparsemx_Validate(sxf, errbuf)    != eslOK) esl_fatal("%s:\n %s", msg, errbuf);
      if ( h4_sparsemx_Validate(sxb, errbuf)    != eslOK) esl_fatal("%s:\n %s", msg, errbuf);
      if ( h4_sparsemx_Validate(sxd, errbuf)    != eslOK) esl_fatal("%s:\n %s", msg, errbuf);
      if ( esl_FCompare(vsc_s, psc_s, 0.0, tol) != eslOK) esl_fatal(msg);                      // 4. V score = Vit path score
      if ( fsc_s > fsc_r+tol)                             esl_fatal(msg);                      // 5. Reference F,B,V score >= sparse F,B,V score 
      if ( bsc_s > bsc_r+tol)                             esl_fatal(msg);                      
      if ( vsc_s > vsc_r+tol)                             esl_fatal(msg);                      
      if ( h4_sparsemx_CompareReferenceAsBound(sxv, rxv, tol) != eslOK) esl_fatal(msg);        // 6. All matrix values satisfy v_ref >= v_sparse 
      if ( h4_sparsemx_CompareReferenceAsBound(sxf, rxf, tol) != eslOK) esl_fatal(msg);
      if ( h4_sparsemx_CompareReferenceAsBound(sxb, rxb, tol) != eslOK) esl_fatal(msg);
      if ( h4_path_Validate(vpi, hmm->M, sq->n, errbuf) != eslOK) esl_fatal("%s:\n %s", msg, errbuf);

      if ( esl_sq_Reuse(sq) != eslOK) esl_fatal(msg);
    }


  h4_path_Destroy(epi);      h4_path_Destroy(vpi);
  h4_sparsemx_Destroy(sxv);  h4_sparsemx_Destroy(sxf);  h4_sparsemx_Destroy(sxb);  h4_sparsemx_Destroy(sxd);
  h4_refmx_Destroy(rxv);     h4_refmx_Destroy(rxf);     h4_refmx_Destroy(rxb);     h4_refmx_Destroy(rxd);
  h4_sparsemask_Destroy(sm);
  h4_checkptmx_Destroy(cpx);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
}


/* utest_compare_reference()
 * 
 * The "compare_reference" unit test tests that the sparse
 * implementations of Viterbi, Forward, Backward, and Decoding give
 * the same results as reference implementation, when all cells are
 * included in the sparse mask.
 *  
 * Samples a random profile and compare it to 'homologous' sequences,
 * generated from the profile, with the sparse mask set to mark all
 * cells, and tests:
 *  
 * 1. Reference and sparse scores must be identical (V,F,B), 
 *    within a tolerance;
 * 2. Reference and sparse viterbi traces must be identical;
 * 3. All sparse V,F,B,D DP matrix structures must Validate();
 * 4. Sparse Viterbi trace structure must Validate();
 * 5. Reference and sparse matrices are identical within tolerance. 
 */
static void
utest_compare_reference(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, int M, int L, int N)
{
  char           msg[]  = "sparse_dp:: compare_reference unit test failed";
  H4_PROFILE    *hmm    = NULL;
  H4_MODE       *mo     = h4_mode_Create();
  ESL_SQ        *sq     = esl_sq_CreateDigital(abc);       /* space for generated (homologous) target seqs              */
  H4_SPARSEMASK *sm     = h4_sparsemask_Create(M, L);
  H4_SPARSEMX   *sxv    = h4_sparsemx_Create(sm);
  H4_SPARSEMX   *sxf    = h4_sparsemx_Create(sm);
  H4_SPARSEMX   *sxb    = h4_sparsemx_Create(sm);
  H4_SPARSEMX   *sxd    = h4_sparsemx_Create(sm);
  H4_REFMX      *rxv    = h4_refmx_Create(M, L);
  H4_REFMX      *rxf    = h4_refmx_Create(M, L);
  H4_REFMX      *rxb    = h4_refmx_Create(M, L);
  H4_REFMX      *rxd    = h4_refmx_Create(M, L);
  H4_PATH       *rpi    = h4_path_Create();
  H4_PATH       *spi    = h4_path_Create();
  int            idx;
  float          vsc_s, fsc_s, bsc_s;
  float          vsc_r, fsc_r, bsc_r;
  float          tol = ( h4_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  char           errbuf[eslERRBUFSIZE];
  
  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  /* choose a random alignment mode: */
  switch (esl_rnd_Roll(rng, 6)) {
  case 0: h4_mode_SetDefault(mo);   break;
  case 1: h4_mode_SetLocal(mo);     break;
  case 2: h4_mode_SetGlocal(mo);    break;
  case 3: h4_mode_SetUnihit(mo);    break;
  case 4: h4_mode_SetUnilocal(mo);  break;
  case 5: h4_mode_SetUniglocal(mo); break;
  }

  for (idx = 0; idx < N; idx++)
    {
      /* Generate a sequence */
      if ( h4_mode_SetLength(mo, L)                   != eslOK) esl_fatal(msg);                              // set length model to L for emission
      do { if ( h4_emit(rng, hmm, mo, sq, NULL)       != eslOK) esl_fatal(msg); } while (sq->n > (M+L)*3);   // keep seq length reasonable
      if ( h4_mode_SetLength(mo, sq->n)               != eslOK) esl_fatal(msg);   // reset length model to this sampled seq's length before scoring or DP 

      /* Mark all cells in sparse mask */
      if ( h4_sparsemask_Reinit(sm, hmm->M, sq->n) != eslOK) esl_fatal(msg);
      if ( h4_sparsemask_AddAll(sm)                != eslOK) esl_fatal(msg);

      /* Sparse DP calculations  */
      if ( h4_sparse_Viterbi (sq->dsq, sq->n, hmm, mo, sm, sxv, spi, &vsc_s) != eslOK) esl_fatal(msg);
      if ( h4_sparse_Forward (sq->dsq, sq->n, hmm, mo, sm, sxf,      &fsc_s) != eslOK) esl_fatal(msg);
      if ( h4_sparse_Backward(sq->dsq, sq->n, hmm, mo, sm, sxb,      &bsc_s) != eslOK) esl_fatal(msg);
      if ( h4_sparse_Decoding(sq->dsq, sq->n, hmm, mo, sxf, sxb, sxd)        != eslOK) esl_fatal(msg);

      /* Reference DP calculations */
      if ( h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxv, rpi, &vsc_r) != eslOK) esl_fatal(msg);
      if ( h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc_r) != eslOK) esl_fatal(msg);
      if ( h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxb,      &bsc_r) != eslOK) esl_fatal(msg);
      if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxb, rxd)    != eslOK) esl_fatal(msg);

      /* Tests */
      if ( esl_FCompare(vsc_s, vsc_r, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg);                      // (1): reference and sparse scores (V,F,B) identical within tolerance
      if ( esl_FCompare(fsc_s, fsc_r, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg);
      if ( esl_FCompare(bsc_s, bsc_r, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg);
      if ( h4_path_Compare(rpi, spi)                     != eslOK) esl_fatal(msg);                      // (2): reference, sparse Viterbi tracebacks identical; see notes on p7_trace_CompareLoosely
      if ( h4_sparsemx_Validate(sxv, errbuf)             != eslOK) esl_fatal("%s:\n  %s", msg, errbuf); // (3): All sparse DP matrices Validate() (V,F,B,D)
      if ( h4_sparsemx_Validate(sxf, errbuf)             != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
      if ( h4_sparsemx_Validate(sxb, errbuf)             != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
      if ( h4_sparsemx_Validate(sxd, errbuf)             != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);
      if ( h4_path_Validate(spi, hmm->M, sq->n, errbuf)  != eslOK) esl_fatal("%s:\n  %s", msg, errbuf); // (4): Sparse DP Viterbi trace must Validate() 
      if ( h4_path_Validate(rpi, hmm->M, sq->n, errbuf)  != eslOK) esl_fatal("%s:\n  %s", msg, errbuf); 
      if ( h4_sparsemx_CompareReference(sxv, rxv, tol)   != eslOK) esl_fatal(msg);                      // (5): Sparse and reference DP matrices all identical within tolerance
      if ( h4_sparsemx_CompareReference(sxf, rxf, tol)   != eslOK) esl_fatal(msg); 
      if ( h4_sparsemx_CompareReference(sxb, rxb, tol)   != eslOK) esl_fatal(msg);
      if ( h4_sparsemx_CompareReference(sxd, rxd, tol)   != eslOK) esl_fatal(msg);
      
      if ( esl_sq_Reuse(sq) != eslOK) esl_fatal(msg);
    }

  h4_path_Destroy(rpi);      h4_path_Destroy(spi);
  h4_refmx_Destroy(rxv);     h4_refmx_Destroy(rxf);     h4_refmx_Destroy(rxb);     h4_refmx_Destroy(rxd);
  h4_sparsemx_Destroy(sxv);  h4_sparsemx_Destroy(sxf);  h4_sparsemx_Destroy(sxb);  h4_sparsemx_Destroy(sxd);
  h4_sparsemask_Destroy(sm);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
}
#endif //h4SPARSE_DP_TESTDRIVE
/*------------------- end, unit tests ---------------------------*/



/***************************************************************** 
 * 9. Test driver
 *****************************************************************/
#ifdef h4SPARSE_DP_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                   docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show brief help summary",                   0 },
  { "-s",         eslARG_INT,     "0", NULL, NULL,  NULL,   NULL,  NULL, "set random number generator seed",          0 },
  { "-L",         eslARG_INT,    "10", NULL, NULL,  NULL,   NULL,  NULL, "set length model for sampled seqs",         0 },
  { "-M",         eslARG_INT,    "10", NULL, NULL,  NULL,   NULL,  NULL, "set test profile length",                   0 },
  { "-N",         eslARG_INT,   "100", NULL, NULL,  NULL,   NULL,  NULL, "number of profile/seq comparisons to test", 0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,   NULL,  NULL, "show HMMER version number",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv, "test driver for reference_aec", "[-options]");
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);
  int             L   = esl_opt_GetInteger(go, "-L");
  int             M   = esl_opt_GetInteger(go, "-M");
  int             N   = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng)); 

  utest_generation       (rng, abc, M, L, N);
  utest_compare_reference(rng, abc, M, L, N);

  fprintf(stderr, "#  status   = ok\n");

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif //h4SPARSE_DP_TESTDRIVE
/*------------------- end, test driver --------------------------*/




/***************************************************************** 
 * 10. Example
 *****************************************************************/
#ifdef h4SPARSE_DP_EXAMPLE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_checkptmx.h"
#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_sparsemask.h"
#include "h4_sparsemx.h"

#include "fbfilter.h"
#include "general.h"
#include "sparse_dp.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                  docgroup*/
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "include all cells in sparse mx",          0 },
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help",                         0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Backward DP matrix for examination", 0 },
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Decoding DP matrix for examination", 0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Forward DP matrix for examination",  0 },
  { "-M",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump sparse mask for examination",        0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Viterbi path for examination",       0 }, 
  { "-V",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Viterbi DP matrix for examination",  0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show HMMER version info",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of using the sparse Vit/Fwd/Bck/Decoding DP implementations";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  H4_CHECKPTMX   *cpx     = h4_checkptmx_Create(100, 100, ESL_MBYTES(32));
  H4_SPARSEMASK  *sm      = h4_sparsemask_Create(100, 100);
  H4_SPARSEMX    *sxv     = h4_sparsemx_Create(NULL);
  H4_SPARSEMX    *sxf     = h4_sparsemx_Create(NULL);
  H4_SPARSEMX    *sxb     = h4_sparsemx_Create(NULL);
  H4_SPARSEMX    *sxd     = h4_sparsemx_Create(NULL);
  H4_PATH        *vpi     = h4_path_Create();
  float           ffsc, vsc, fsc, bsc, pisc;
  int             status;

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);

  if ( esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, /*env=*/NULL, &sqfp) != eslOK) esl_fatal("couldn't open sequence file %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      h4_mode_SetLength(mo, sq->n);

      if (esl_opt_GetBoolean(go, "-a"))
        {
          if ( h4_sparsemask_Reinit(sm, hmm->M, sq->n) != eslOK) esl_fatal("sparsemask reinit failed");
          if ( h4_sparsemask_AddAll(sm)                != eslOK) esl_fatal("sparsemask add all failed");
          ffsc = -eslINFINITY;
        }
      else 
        {
          if ( h4_fwdfilter(sq->dsq, sq->n, hmm, mo, cpx, &ffsc)                 != eslOK) esl_fatal("fwdfilter failed");
          if ( h4_bckfilter(sq->dsq, sq->n, hmm, mo, cpx, sm, h4SPARSIFY_THRESH) != eslOK) esl_fatal("bckfilter failed");
        }
          
      if (esl_opt_GetBoolean(go, "-M")) h4_sparsemask_Dump(stdout, sm);

      if ( h4_sparse_Viterbi (sq->dsq, sq->n, hmm, mo, sm, sxv, vpi, &vsc) != eslOK) esl_fatal("sparse Viterbi failed");
      if ( h4_sparse_Forward (sq->dsq, sq->n, hmm, mo, sm, sxf,      &fsc) != eslOK) esl_fatal("sparse Forward failed");
      if ( h4_sparse_Backward(sq->dsq, sq->n, hmm, mo, sm, sxb,      &bsc) != eslOK) esl_fatal("sparse Backward failed");
      if ( h4_sparse_Decoding(sq->dsq, sq->n, hmm, mo, sxf, sxb, sxd)      != eslOK) esl_fatal("sparse Decoding failed");
      if ( h4_path_Score(vpi, sq->dsq, hmm, mo, &pisc)                     != eslOK) esl_fatal("path score failed");

      if (esl_opt_GetBoolean(go, "-V")) h4_sparsemx_Dump(stdout, sxv);
      if (esl_opt_GetBoolean(go, "-F")) h4_sparsemx_Dump(stdout, sxf);
      if (esl_opt_GetBoolean(go, "-B")) h4_sparsemx_Dump(stdout, sxb);
      if (esl_opt_GetBoolean(go, "-D")) h4_sparsemx_Dump(stdout, sxd);
      if (esl_opt_GetBoolean(go, "-P")) h4_path_Dump(stdout, vpi);

      printf("%12s  fwd filter raw score: %.4f bits\n", sq->name, ffsc);
      printf("%12s  sparse Vit raw score: %.4f bits\n", sq->name, vsc);
      printf("%12s  sparse Fwd raw score: %.4f bits\n", sq->name, fsc);
      printf("%12s  sparse Bck raw score: %.4f bits\n", sq->name, bsc);
      printf("%12s  Vit path raw score:   %.4f bits\n", sq->name, pisc);
      
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);

  h4_sparsemask_Destroy(sm);
  h4_sparsemx_Destroy(sxv);
  h4_sparsemx_Destroy(sxf);
  h4_sparsemx_Destroy(sxb);
  h4_sparsemx_Destroy(sxd);
  h4_path_Destroy(vpi);
  h4_checkptmx_Destroy(cpx);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif //h4SPARSE_DP_EXAMPLE
/*-------------------- end, example ---------------------------- */
