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

#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_sparsemask.h"
#include "h4_sparsemx.h"

#include "logsum.h"

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
  if (opt_pi && (status = h4_path_Reuse(opt_pi)) != eslOK) return status;
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
  if (opt_pi && xC != -eslINFINITY) return eslOK; //h4_sparse_ViterbiTrace(hmm, sx, opt_pi);
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
 *            <dsq> of length <L>, by the Forward algorithm, using
 *            sparse dynamic programming constrained by the sparse
 *            mask <sm>.  Fill in the sparse Forward matrix <sx>,
 *            and (optionally) return the Forward raw score in
 *            bits in <*opt_fsc>.
 *            
 *            <sx> is provided as existing space; reused/reallocated
 *            as needed.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of <dsq>
 *            hmm     - profile
 *            sm      - sparse mask
 *            sx      - Forward matrix to fill (reallocated/reused)
 *            opt_fsc - optRETURN: raw Forward score in bits; or NULL
 *
 * Returns:   <eslOK> on success. <sx> is the Forward DP matrix;
 *            raw Forward bitscore is optionally in <*opt_fsc>.
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
                                                 dpp[h4S_DL] + hmm->tsc[k][h4_DI]);  // +ISC(k), if we weren't enforcing it to zero
	    *dpc++ = igc = h4_logsum( h4_logsum( dpp[h4S_MG] + hmm->tsc[k][h4_MI],
                                                 dpp[h4S_IG] + hmm->tsc[k][h4_II]),
                                      h4_logsum( dpp[h4S_DG] + hmm->tsc[k][h4_DI],
                                                 xG          + hmm->tsc[k][h4_GI])); // ditto
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
v_select_il(ESL_RANDOMNESS *r, const H4_PROFILE *gm, int k, const float *dpp, const int *kp, int np, int *ret_z)
{
  ESL_UNUSED(r);
  int   state[3] = { h4P_MG, h4P_IG, h4P_DG };
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
      int   state[4] = { h4P_MG, h4P_IG, h4P_DG, h4P_G };
      float path[4];

      path[0] = dpp[h4S_ML] + hmm->tsc[k][h4_MI];
      path[1] = dpp[h4S_IL] + hmm->tsc[k][h4_II];
      path[2] = dpp[h4S_DL] + hmm->tsc[k][h4_DI];
      path[3] =  xp[h4S_G]  + hmm->tsc[k][h4_GI];
      *ret_z = y;
      return state[esl_vec_FArgMax(path, 4)]
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
  path[1] = dpp[h4S_DG] + hmm->tsc[k-1][h4_DD];
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
  int   z;

  for (z = 0; z < np; z++)
    { 
      if (dpp[h4S_ML] >= max) { max = dpp[h4S_ML]; smax = h4P_ML; zmax = z; }   // don't need to check DLk->E path; these can't occur in a Viterbi path 
      
      if (kp[z] == hmm->M) // final M cell doesn't use DGE wing folding, and it has no I.
        {
          if (dpp[h4S_MG] > max) { max = dpp[h4_MG]; smax = h4P_MG; zmax = z; }
          if (dpp[h4S_DG] > max) { max = dpp[h4_DG]; smax = h4P_DG; zmax = z; }
        }
      else if (z == np-1 || kp[z+1] != kp[z]+1) // sparse cell k with no following cell k+1: check glocal exit Mk->...->E across unmarked cells
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
sparse_traceback_engine(ESL_RANDOMNESS *rng, float *wrk, const H4_PROFILE *hmm, H4_MODE *mo, const H4_SPARSEMX *sx, H4_PATH *pi)
{
  const H4_SPARSEMASK *sm = sx->sm;
  int            k  = 0;	// current coord in profile consensus
  int            i;	        // current coord in sequence (that snxt and dp are on)
  int            scur, snxt;
  int            z;		// current position in n[i] sparse k array entries on current row dp[]
  int            g;             // what segment we're in, 1..sm->S
  int            status;

  ESL_DASSERT1(( pi->N == 0 ));   // path must be empty - _Reuse() was called on it by our parent

  int (*select_ml)(ESL_RANDOMNESS *rng,             const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np, const float *xp, int *ret_z);
  int (*select_mg)(ESL_RANDOMNESS *rng,             const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np, const float *xp, int *ret_z);
  int (*select_il)(ESL_RANDOMNESS *rng,             const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np,                  int *ret_z);
  int (*select_ig)(ESL_RANDOMNESS *rng,             const H4_PROFILE *hmm, int k, const float *dpp, const int *kp, int np,                  int *ret_z);
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
  i    = sm->ib[sm->S];  
  scur = h4P_C;
  xc   = sx->xmx + (sm->nrow+sm->S-1) * h4S_NXCELLS;        // to last xc supercell, which is ib[S]
  dp   = sx->dp  + (sm->ncells - sm->n[i]) * h4S_NSCELLS;   // to *start* of row ib[S], so we can index with z
  if ((status = h4_path_AppendElement(pi, h4P_C, L-i+1)) != eslOK) return status;  // L-i residues, +1 for emit-on-transition

  for (g = sm->S; g >= 1; g--)
    {
      while ( i >= sm->ia[g] || (scur != h4P_N && scur != h4P_J))
        {
          switch (scur) {
          case h4P_ML: snxt = (*select_ml)(rng, hmm, k, dp, sm->k[i-1], sm->n[i-1], xc-h4S_NXCELLS, &z); i--; k--;      xc -= h4S_NXCELLS; break;
          case h4P_MG: snxt = (*select_mg)(rng, hmm, k, dp, sm->k[i-1], sm->n[i-1], xc-h4S_NXCELLS, &z); i--; k--;      xc -= h4S_NXCELLS; break;
          case h4P_IL: snxt = (*select_il)(rng, hmm, k, dp, sm->k[i-1], sm->n[i-1],                 &z); i--;           xc -= h4S_NXCELLS; break;
          case h4P_IG: snxt = (*select_ig)(rng, hmm, k, dp, sm->k[i-1], sm->n[i-1],                 &z); i--;           xc -= h4S_NXCELLS; break;
          case h4P_DL: snxt = (*select_dl)(rng, hmm, k, dp + (z-1)*h4S_NSCELLS);                              k--; z--;                    break; 
          case h4P_DG: snxt = (*select_dg)(rng, hmm, k, dp + (z-1)*h4S_NSCELLS);                              k--; z--;                    break;
          case h4P_N:  snxt = h4P_N;                                                                                                       break;
          case h4P_J:  snxt = (*select_j)(rng, mo, xc);                                                                                    break; 
          case h4P_C:  snxt = (*select_c)(rng, mo, xc);                                                                                    break; // connect to E(i), C(i-1). E(i) valid if n[i]>0; and if E(i) valid, C(i-1) must be stored too. if E(i) invalid, connect to C, and it doesn't matter if xc is valid or not
          case h4P_E:  snxt = (*select_e)(rng, wrk, hmm, dp, sm->k[i], sm->n[i], &z);                         k=sm->k[i][z];               break;
          case h4P_B:  snxt = (*select_b)(rng, mo, xc);                                                                                    break; // {NJ}(i) -> B(i). If we reached B(i), xc[i] valid, so NJ must also be valid.
          case h4P_L:  snxt = h4P_B;                                                                                                       break;
          case h4P_G:  snxt = h4P_B;                                                                                                       break;
          default:     ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
          }
      
          /* Glocal B->G->{MI}k left wing retraction entry: unfold it */
          if (snxt == h4P_G && k > 0) {
            if (status = h4_path_AppendElement(pi, h4P_DG, k) != eslOK) return status;
            k=0;
          }
          /* Glocal Mk->E right wing retraction: off last sparse cell k, {MID}k->Dk+1->E */
          if (scur == h4P_E && k < hmm->M) {
            if (status = h4_path_AppendElement(pi, h4P_DG, hmm->M-k) != eslOK) return status;
          }
      
          if (snxt == h4P_L) status = h4_path_AppendElement(pi, h4P_L, k+1); // if we just backed into an L state, record the k of L->Mk in its rle[].
          else               status = h4_path_Append       (pi, snxt);       // note that this appends G if that's our scur
          if (status != eslOK) return status;

          if (snxt == h4P_ML || snxt == h4P_MG || snxt == h4P_IL || snxt == h4P_IG) 
            dp -= sm->n[i] * h4S_NSCELLS; 
          else if ( h4_path_IsX(snxt) && scur == snxt) // deferred i decrement for NN/JJ, because we have to check for NN/JJ emit as opposed to {NJ}->B init
            {
              i--;
              dp -= sm->n[i] * h4S_NSCELLS;  // if we've reached i=ia[g]-1 and are about to finish a segment, this does nothing (no dp row there)
              xc -= h4S_NXCELLS;             //   ... whereas this puts xc on ia[g]-1 (xc does have a supercell there)
            }
          scur = snxt;
        }

      // Now we're on ia[g]-1, in an N|J state. Residues ib[g-1]+1..ia[g]-1 are between segments, emitted by the same N|J state.
      i = sm->ib[g-1];
      dp -= sm->n[i] * h4S_NSCELLS;   // put dp on ib[g-1]
      xc -= h4S_NXCELLS;              // and xc
      if ((status = h4_path_AppendSeveral(pi, scur, sm->ia[g] - sm->ib[g-1] - 1) != eslOK) return status;
    }

  return h4_path_Reverse(pi);
}
/*------------- end, traceback engine and API -------------------*/


/***************************************************************** 
 * 8. Unit tests
 *****************************************************************/
/* TK TK */
/*------------------- end, unit tests ---------------------------*/



/***************************************************************** 
 * 9. Test driver
 *****************************************************************/
/* TK TK */
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
  /* name           type      default  env  range  toggles reqs incomp  help                        docgroup*/
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "include all cells in sparse mx",     0 },
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help",                    0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Vit DP matrix for examination", 0 },
  { "-V",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Fwd DP matrix for examination", 0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show HMMER version info",            0 },
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
  H4_PATH        *vpi     = h4_path_Create();
  float           ffsc, vsc, fsc;
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

      if ( h4_sparse_Viterbi(sq->dsq, sq->n, hmm, mo, sm, sxv, vpi, &vsc)    != eslOK) esl_fatal("sparse Viterbi failed");
      if ( h4_sparse_Forward(sq->dsq, sq->n, hmm, mo, sm, sxf,      &fsc)    != eslOK) esl_fatal("sparse Forward failed");
          
      if (esl_opt_GetBoolean(go, "-F")) h4_sparsemx_Dump(stdout, sxf);
      if (esl_opt_GetBoolean(go, "-V")) h4_sparsemx_Dump(stdout, sxv);

      printf("%12s  fwd filter raw score: %.2f bits\n", sq->name, ffsc);
      printf("%12s  sparse Vit raw score: %.2f bits\n", sq->name, vsc);
      printf("%12s  sparse Fwd raw score: %.2f bits\n", sq->name, fsc);

      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);

  h4_sparsemask_Destroy(sm);
  h4_sparsemx_Destroy(sxv);
  h4_sparsemx_Destroy(sxf);
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
