/* Sparse dynamic programming: Viterbi implementation.
 * 
 * Contents:
 *   1. Sparse Viterbi
 *   2. Copyright and license information
 */
#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_trace.h"
#include "dp_sparse/sparse_viterbi.h"

/*****************************************************************
 * 1. Sparse Viterbi
 *****************************************************************/

/* Function:  p7_SparseViterbi()
 * Synopsis:  Viterbi optimal path algorithm, in sparse DP.
 *
 * Purpose:   Compare profile <gm> to digital sequence <dsq> of length <L>,
 *            by the Viterbi algorithm, using sparse dynamic programming,
 *            as constrained by the sparse mask <sm>.
 *            Fill in the sparse Viterbi matrix <sx>; (optionally) trace
 *            back the optimal path and return it in the trace structure <opt_tr>
 *            if the caller provides one; and (optionally) return the 
 *            Viterbi raw score in nats in <*opt_sc>.
 *            
 *            <sx> can be reused from previous calculations, even
 *            smaller ones; see <p7_sparsemx_Reuse()>. If necessary,
 *            it will be reallocated here, to be large enough for the
 *            <gm->M> by <L> calculation restricted to masked cells
 *            <sm>.
 *            
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of <dsq>
 *            gm      - profile
 *            sm      - sparse mask
 *            sx      - Viterbi matrix to fill; (may be reallocated here)
 *            opt_tr  - optRESULT: trace structure with optimal traceback; or NULL if caller doesn't want it
 *            opt_sc  - optRETURN: raw Viterbi score in nats; or NULL if result unwanted
 *
 * Returns:   <eslOK> on success; <opt_tr>, if non-NULL, contains the optimal traceback;
 *            and <*opt_sc> optionally contains the raw Viterbi score.
 */
int
p7_SparseViterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMASK *sm, P7_SPARSEMX *sx,  P7_TRACE *opt_tr, float *opt_sc)
{
  float const   *tsc    = gm->tsc;	 /* sets up TSC() macro, access to profile's transitions */
  float const   *rsc;			 /* will be set up for MSC(), ISC() macros for residue scores */
  float         *xpc; 	                 /* ptr that steps through current special cells   */
  float         *dpc;	                 /* ptr to step thru current row i main DP cells */
  float         *dpp;			 /* ptr to step thru previous row i-1 main DP cells */
  float         *last_dpc;		 /* used to reinit dpp after each sparse row computation */
  int            ng;
  float          xE, xN, xJ, xB, xL, xG, xC;  /* tmp scores on special states. only stored when in row bands, and on ia-1 before a seg */
  float          mlc, mgc;		 /* temporary score calculations M(i,k)         */
  float          dlc, dgc;		 /* precalculated D(i,k+1) value on current row */
  int           *kc = sm->k[0];		 /* <kc> points to the list of sparse cell indices k for current row i */
  int           *kp;			 /* <kp> points to the previous row's sparse cell index list */
  int            i,k;	      	         /* i,k row,col (seq position, profile position) cell coords */
  int            y,z;			 /* indices in lists of k coords on prev, current row */
  int            status;

  /* Contract checks on arguments */
  ESL_DASSERT1( (sm->L == L) );
  ESL_DASSERT1( (sm->M == gm->M) );

  /* Assure that <sx> is allocated large enough (we might be reusing it).
   * Set its type now, so we can Dump/Validate/etc. during debugging this routine, if needed.
   */
  if ( (status = p7_sparsemx_Reinit(sx, sm)) != eslOK) return status;
  sx->type = p7S_VITERBI;

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
      if (i == 1 || ng) {
	*xpc++ = xE = -eslINFINITY;
	*xpc++ = xN  = xN + ( ng ? ng * gm->xsc[p7P_N][p7P_LOOP] : 0.0); /* test ng, because we must watch out for 0*-inf special case */
	*xpc++ = xJ  = xJ + ( ng ? ng * gm->xsc[p7P_J][p7P_LOOP] : 0.0);
	*xpc++ = xB  = ESL_MAX( xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]);
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
	  while (y < sm->n[i-1] && kp[y]  < k-1) { y++; dpp+=p7S_NSCELLS; }
	  if    (y < sm->n[i-1] && kp[y] == k-1) {
	    mlc = ESL_MAX( ESL_MAX( dpp[p7R_ML] + TSC(p7P_MM, k-1),
				    dpp[p7R_IL] + TSC(p7P_IM, k-1)),
			   ESL_MAX( dpp[p7R_DL] + TSC(p7P_DM, k-1),
				    mlc));        
	    mgc = ESL_MAX( ESL_MAX( dpp[p7R_MG] + TSC(p7P_MM, k-1),
				    dpp[p7R_IG] + TSC(p7P_IM, k-1)),
			   ESL_MAX( dpp[p7R_DG] + TSC(p7P_DM, k-1),
				    mgc));
	  }
	  *dpc++ = mlc = MSC(k) + mlc;
	  *dpc++ = mgc = MSC(k) + mgc;

	  /* Try to find cell i-1,k; then compute I(i,k) from it */
	  while (y < sm->n[i-1] && kp[y] < k)  { y++; dpp+=p7S_NSCELLS; }
	  if    (y < sm->n[i-1] && kp[y] == k) {
	    *dpc++ = ISC(k) + ESL_MAX( dpp[p7R_ML] + TSC(p7P_MI,k),  dpp[p7R_IL] + TSC(p7P_II, k));
	    *dpc++ = ISC(k) + ESL_MAX( dpp[p7R_MG] + TSC(p7P_MI,k),  dpp[p7R_IG] + TSC(p7P_II, k));
	  } else {
	    *dpc++ = -eslINFINITY;
	    *dpc++ = -eslINFINITY;
	  }
	    
	  /* local exit paths (a F/V difference here: in V, no Dk->E path can win */
	  xE = ESL_MAX(xE, mlc);

	  /* delayed store of Dk; advance calculation of next D_k+1 */
	  *dpc++ = dlc;
	  *dpc++ = dgc;
	  if (z < sm->n[i]-1 && kc[z+1] == k+1) { /* is there a (i,k+1) cell to our right? */
	    dlc = ESL_MAX( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
	    dgc = ESL_MAX( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	  } else {  		/* if not, we MUST consider {MD}Gk->Dk+1..E glocal exit path, even from internal sparse cells - not just last cell! */
	    xE  = ESL_MAX( xE,  TSC(p7P_DGE, k) + ESL_MAX( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k)));  // yes, the D path can contribute; we only use wing-retraction on sparse cells k where k+1 is unmarked; if k=M, for example, we must check D->E
	    dlc = dgc = -eslINFINITY;
	  }
	}

      *xpc++ = xE;  // we already max'ed over all Mk->E exits, both local and glocal
      *xpc++ = xN = xN + gm->xsc[p7P_N][p7P_LOOP];
      *xpc++ = xJ = ESL_MAX( xJ + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_LOOP]);
      *xpc++ = xB = ESL_MAX( xJ + gm->xsc[p7P_J][p7P_MOVE],  xN + gm->xsc[p7P_N][p7P_MOVE]);
      *xpc++ = xL = xB + gm->xsc[p7P_B][0]; /* B->L */
      *xpc++ = xG = xB + gm->xsc[p7P_B][1]; /* B->G */
      *xpc++ = xC = ESL_MAX( xE + gm->xsc[p7P_E][p7P_MOVE],  xC + gm->xsc[p7P_C][p7P_LOOP]);
      *xpc++      = -eslINFINITY; /* JJ: this space only used in a Decoding matrix. */
      *xpc++      = -eslINFINITY; /* CC: this space only used in a Decoding matrix. */

      /* now dpc is on the start of the next sparsified row */
      dpp = last_dpc;
    }

  xC += ( ng ? ng *  gm->xsc[p7P_C][p7P_LOOP] : 0.0f) + gm->xsc[p7P_C][p7P_MOVE];

  if (opt_sc) *opt_sc = xC;
  if (opt_tr && xC != -eslINFINITY) return p7_sparse_trace_Viterbi(gm, sx, opt_tr);
  else                              return eslOK;
}
/*-------------------- end, Viterbi -----------------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
