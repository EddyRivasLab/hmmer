/* Sparse DP routines.
 *
 * Contents:
 *   1. Viterbi
 * 
 */
#include "h4_config.h"

#include "easel.h"

#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_sparsemask.h"
#include "h4_sparsemx.h"


/*****************************************************************
 * 1. Sparse Viterbi
 *****************************************************************/

/* Function:  h4_sparse_Viterbi()
 * Synopsis:  Viterbi optimal path algorithm, in sparse DP.
 *
 * Purpose:   Compare profile <hmm> in mode <mo> to digital sequence
 *            <dsq> of length <L>, by the Viterbi algorithm, using
 *            sparse dynamic programming, as constrained by the sparse
 *            mask <sm>.  Fill in the sparse Viterbi matrix <sx>;
 *            (optionally) trace back the optimal path and return it
 *            in the trace structure <opt_pi> if the caller provides
 *            one; and (optionally) return the Viterbi raw score in
 *            bits in <*opt_vsc>.
 *            
 *            <sx> and <opt_pi> are provided as existing space; 
 *            they will be reused/reallocated as needed.
 *            
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of <dsq>
 *            gm      - profile
 *            sm      - sparse mask
 *            sx      - Viterbi matrix to fill (reallocated/reused)
 *            opt_pi  - optRESULT: optimal path (reallocated/reused); or NULL if unwanted
 *            opt_vsc - optRETURN: raw Viterbi score in bits; or NULL
 *
 * Returns:   <eslOK> on success; <opt_pi>, if non-NULL, contains Viterbi path;
 *            and <*opt_vsc> optionally contains the raw Viterbi score.
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
  if ( (status = h4_path_Reuse(opt_pi))      != eslOK) return status;
  if ( (status = h4_sparsemx_Reinit(sx, sm)) != eslOK) return status;
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
                                             dpp[h4S_DG] + hmm->tsc[k][h4_DI]);  // ditto
	  } else {
	    *dpc++ = ilc = -eslINFINITY;
	    *dpc++ = igc = -eslINFINITY;
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
	  } else {  		// if not, we MUST consider {MD}Gk->Dk+1..E glocal exit path, even from internal sparse cells - not just last cell! 
	    xE  = ESL_MAX( xE,
                           hmm->tsc[k][h4_DGE] + ESL_MAX( ESL_MAX( mgc + hmm->tsc[k][h4_MD],
                                                                   igc + hmm->tsc[k][h4_ID]),
                                                                   dgc + hmm->tsc[k][h4_DD]));  // yes, the D path can contribute; we only use wing-retraction on sparse cells k where k+1 is unmarked; if k=M, for example, we must check D->E
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

  xC += ( ng ? ng * mo->xsc[h4_C][h4_LOOP] : 0.0f) + mo->xsc[h4_C][h4_MOVE];

  if (opt_vsc) *opt_vsc = xC;
  if (opt_pi && xC != -eslINFINITY) return eslFAIL; //h4_sparse_ViterbiTrace(hmm, sx, opt_pi);
  else                              return eslOK;
}
/*-------------------- end, Viterbi -----------------------------*/
