/* Reference implementation of Viterbi scoring and alignment.
 *   dual-mode (local/glocal)
 *   quadratic memory (not banded, not checkpointed)
 *   standard C code (not striped, not vectorized)
 *   
 * The reference implementation is for testing and debugging, and to
 * provide an example simpler than our production DP code, which
 * layers on some more complicated techniques (banding, vectorization,
 * checkpointing).
 * 
 * Contents:
 *   1. Viterbi DP fill.
 *   2. Viterbi optimal traceback.
 *   3. Copyright and license information.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_refmx.h"

static int reference_viterbi_traceback(const P7_PROFILE *gm, const P7_REFMX *rmx, P7_TRACE *tr);

/*****************************************************************
 * 1. Viterbi DP fill.
 *****************************************************************/

/* Function:  p7_ReferenceViterbi()
 * Synopsis:  Reference implementation of the Viterbi algorithm.
 *
 * Purpose:   Given a target sequence <dsq> of length <L> residues, a
 *            query profile <gm>, and a DP matrix <rmx> that the
 *            caller has allocated for a <gm->M> by <L> comparison,
 *            do the Viterbi optimal alignment algorithm.
 *            Return the Viterbi score in nats in <*opt_sc>, if
 *            caller provides it. 
 *            Return the Viterbi optimal trace in <*opt_tr>, if
 *            caller provides an allocated trace structure.
 *            
 * Args:      dsq     - digital target sequence 1..L
 *            L       - length of <dsq> in residues
 *            gm      - query profile
 *            rmx     - DP matrix, allocated gm->M by L
 *            opt_tr  - optRETURN: trace structure for Viterbi alignment, or NULL
 *            opt_sc  - optRETURN: Viterbi raw score in nats, or NULL
 *
 * Returns:   <eslOK> on success. <rmx> contains the Viterbi matrix.
 */
int
p7_ReferenceViterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, P7_TRACE *opt_tr, float *opt_sc)
{
  int    M          = gm->M;
  float *dpc, *dpp;
  const float *tsc;		/* ptr for stepping thru profile's transition parameters */
  const float *rsc;		/* ptr for stepping thru profile's emission parameters   */
  int    i, k, s;
  float  mlv, mgv;	      /* ML,MG cell values on current row   */
  float  dlv, dgv; 	      /* pushed-ahead DL,DG cell k+1 values */
  float  xE, xN, xJ, xB, xL, xG;
  
  /* Initialization of the zero row. */
  dpc = rmx->dp[0];
  for (s = 0; s < (M+1) * p7R_NSCELLS; s++)
    *dpc++ = -eslINFINITY; 	                               // all M,I,D; k=0..M
  *dpc++ = -eslINFINITY;	                               // E
  *dpc++ = 0.0;			                               // N
  *dpc++ = -eslINFINITY;                                       // J
  *dpc++ = gm->xsc[p7P_N][p7P_MOVE];                           // B
  *dpc++ = xL = gm->xsc[p7P_N][p7P_MOVE] + gm->xsc[p7P_B][0];  // L
  *dpc++ = xG = gm->xsc[p7P_N][p7P_MOVE] + gm->xsc[p7P_B][1];  // G
  *dpc        = -eslINFINITY;                                  // C
  /* *dpc must end and stay on C state, to handle L=0 case where the recursion body below doesn't run */

  /* Main DP recursion */
  for (i = 1; i <= L; i++)
    {
      /* Initialization for a new row */
      rsc = gm->rsc[dsq[i]] + p7P_NR;	/* this ptr steps through the row's emission scores 1..M. skip k=0 */
      tsc = gm->tsc;			/* this ptr steps through profile's transition scores 0..M         */

      dpp = rmx->dp[i-1];               /* previous row dpp is already set, and at k=0 */
      dpc = rmx->dp[i];                 /* current DP row, skip k=0, start at k=1.  */
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;

      dlv = dgv = -eslINFINITY;
      xE  =       -eslINFINITY;

      /* Main inner loop of the recursion */
      for (k = 1; k < M; k++)
	{
	  /* match states MLk, MGk */
	  mlv = *dpc++ = *rsc + ESL_MAX( ESL_MAX(*(dpp+p7R_ML) + *(tsc + p7P_MM),
						 *(dpp+p7R_IL) + *(tsc + p7P_IM)),
					 ESL_MAX(*(dpp+p7R_DL) + *(tsc + p7P_DM),
				 	         xL            + *(tsc + p7P_LM)));

	  mgv = *dpc++ = *rsc + ESL_MAX( ESL_MAX(*(dpp+p7R_MG) + *(tsc + p7P_MM),
						 *(dpp+p7R_IG) + *(tsc + p7P_IM)),
 					 ESL_MAX(*(dpp+p7R_DG) + *(tsc + p7P_DM),
						 xG            + *(tsc + p7P_GM)));

	  rsc++;                /* rsc advances to insert score for position k */
	  tsc += p7P_NTRANS;    /* tsc advances to transitions in states k     */
	  dpp += p7R_NSCELLS;	/* dpp advances to cells for states k          */

	  /* Insert state calculations ILk, IGk. */
	  *dpc++ = *rsc + ESL_MAX( *(dpp + p7R_ML) + *(tsc + p7P_MI), *(dpp + p7R_IL) + *(tsc + p7P_II));
	  *dpc++ = *rsc + ESL_MAX( *(dpp + p7R_MG) + *(tsc + p7P_MI), *(dpp + p7R_IG) + *(tsc + p7P_II));
	  rsc++;		/* rsc advances to next match state emission   */

	  /* E state update; local paths only, transition prob 1.0 in implicit probability model, Dk->E can never be optimal */
	  xE  = ESL_MAX( mlv, xE);

	  /* Delete state, deferred storage trick */
	  *dpc++ = dlv;
	  *dpc++ = dgv;
	  dlv = ESL_MAX( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD));
	  dgv = ESL_MAX( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD));
	}

      /* k=M node is unrolled and handled separately. No I state, and glocal exits. */
      mlv = *dpc++ = *rsc + ESL_MAX( ESL_MAX(*(dpp+p7R_ML) + *(tsc + p7P_MM),
					     *(dpp+p7R_IL) + *(tsc + p7P_IM)),
				     ESL_MAX(*(dpp+p7R_DL) + *(tsc + p7P_DM),
					     xL            + *(tsc + p7P_LM)));

      mgv = *dpc++ = *rsc + ESL_MAX( ESL_MAX(*(dpp+p7R_MG) + *(tsc + p7P_MM),
					     *(dpp+p7R_IG) + *(tsc + p7P_IM)),
				     ESL_MAX(*(dpp+p7R_DG) + *(tsc + p7P_DM),
					     xG            + *(tsc + p7P_GM)));
      dpp  += p7R_NSCELLS; 

      /* I_M state doesn't exist      */
      *dpc++ = -eslINFINITY;	/* IL */
      *dpc++ = -eslINFINITY;	/* IG */

      /* E state update now includes glocal exits: transition prob 1.0 from MG_m; DLk->E still can't be optimal, but DGk->E can */
      xE  = ESL_MAX( ESL_MAX( mgv, dgv),
		     ESL_MAX( xE,  mlv));
      
      /* D_M state: deferred storage only */
      *dpc++ = dlv;
      *dpc++ = dgv;
    
      /* row i is now finished, and dpc[] is positioned exactly on first special state, E */
      dpp += p7R_NSCELLS;    /* now dpp[] is also positioned exactly on first special, E */
      
      *dpc++ = xE;		/* E */
      *dpc++ = xN = *(dpp + p7R_N) + gm->xsc[p7P_N][p7P_LOOP]; /* N */
      *dpc++ = xJ = ESL_MAX( *(dpp + p7R_J) + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_LOOP]); /* J */
      *dpc++ = xB = ESL_MAX(            xN  + gm->xsc[p7P_N][p7P_MOVE],  xJ + gm->xsc[p7P_J][p7P_MOVE]); /* B */
      *dpc++ = xL = xB  + gm->xsc[p7P_B][0]; /* L */
      *dpc++ = xG = xB  + gm->xsc[p7P_B][1]; /* G */
      *dpc        = ESL_MAX( *(dpp + p7R_C) + gm->xsc[p7P_C][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_MOVE]); /* C */
    }
  /* Done with all rows i. As we leave, dpc is still sitting on the xC value for i=L ... including even the L=0 case */
  
  if (opt_sc) *opt_sc = *dpc + gm->xsc[p7P_C][p7P_MOVE];
  rmx->M = M;
  rmx->L = L;
  if (opt_tr) return reference_viterbi_traceback(gm, rmx, opt_tr);
  else        return eslOK;
}
/*------------------ end, viterbi DP fill -----------------------*/


/*****************************************************************
 * 2. Viterbi traceback
 *****************************************************************/

static inline int
select_ml(const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  int   state[4] = { p7T_ML, p7T_IL, p7T_DL, p7T_L };
  float path[4];

  path[0] = P7R_MX(rmx, i-1, k-1, p7R_ML) + P7P_TSC(gm, k-1, p7P_MM);
  path[1] = P7R_MX(rmx, i-1, k-1, p7R_IL) + P7P_TSC(gm, k-1, p7P_IM);
  path[2] = P7R_MX(rmx, i-1, k-1, p7R_DL) + P7P_TSC(gm, k-1, p7P_DM);
  path[3] = P7R_XMX(rmx, i-1, p7R_L)     + P7P_TSC(gm, k-1, p7P_LM);
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
select_mg(const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  int   state[4] = { p7T_MG, p7T_IG, p7T_DG, p7T_G };
  float path[4];

  path[0] = P7R_MX(rmx, i-1, k-1, p7R_MG) + P7P_TSC(gm, k-1, p7P_MM);
  path[1] = P7R_MX(rmx, i-1, k-1, p7R_IG) + P7P_TSC(gm, k-1, p7P_IM);
  path[2] = P7R_MX(rmx, i-1, k-1, p7R_DG) + P7P_TSC(gm, k-1, p7P_DM);
  path[3] = P7R_XMX(rmx, i-1, p7R_G)      + P7P_TSC(gm, k-1, p7P_GM);
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
select_il(const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  float path[2];

  path[0] = P7R_MX(rmx, i-1, k, p7R_ML) + P7P_TSC(gm, k, p7P_MI);
  path[1] = P7R_MX(rmx, i-1, k, p7R_IL) + P7P_TSC(gm, k, p7P_II);
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_IL);
}

static inline int
select_ig(const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  float path[2];

  path[0] = P7R_MX(rmx, i-1, k, p7R_MG) + P7P_TSC(gm, k, p7P_MI);
  path[1] = P7R_MX(rmx, i-1, k, p7R_IG) + P7P_TSC(gm, k, p7P_II);
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_IG);
}


static inline int
select_dl(const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  float path[2];

  path[0] = P7R_MX(rmx, i, k-1, p7R_ML) + P7P_TSC(gm, k-1, p7P_MD);
  path[1] = P7R_MX(rmx, i, k-1, p7R_DL) + P7P_TSC(gm, k-1, p7P_DD);
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_DL);
}

static inline int
select_dg(const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  float path[2];

  path[0] = P7R_MX(rmx, i, k-1, p7R_MG) + P7P_TSC(gm, k-1, p7P_MD);
  path[1] = P7R_MX(rmx, i, k-1, p7R_DG) + P7P_TSC(gm, k-1, p7P_DD);
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_DG);
}

static inline int
select_e(const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int *ret_k)
{
  float max  = -eslINFINITY;
  int   smax = -1;
  int   kmax = -1;
  int   k;

  for (k = 1; k <= gm->M; k++)
    if (P7R_MX(rmx, i, k, p7R_ML) > max) { max = P7R_MX(rmx, i, k, p7R_ML); smax = p7T_ML; kmax = k; }

  if (P7R_MX(rmx, i, gm->M, p7R_MG) > max) { max = P7R_MX(rmx, i, gm->M, p7R_MG); smax = p7T_MG; kmax = k; }
  if (P7R_MX(rmx, i, gm->M, p7R_DG) > max) { max = P7R_MX(rmx, i, gm->M, p7R_DG); smax = p7T_DG; kmax = k; }

  *ret_k = kmax;
  return smax;
}

static inline int
select_n(int i)
{
  return ((i==0) ? p7T_S : p7T_N);
}

static inline int
select_j(const P7_PROFILE *gm, const P7_REFMX *rmx, int i)
{
  float path[2];
  path[0] = P7R_XMX(rmx, i-1, p7R_J) + gm->xsc[p7P_J][p7P_LOOP];
  path[1] = P7R_XMX(rmx, i,   p7R_E) + gm->xsc[p7P_E][p7P_LOOP];
  return ( (path[0] > path[1]) ? p7T_J : p7T_E);
}

static inline int
select_b( const P7_PROFILE *gm, const P7_REFMX *rmx, int i)
{
  float path[2];
  path[0] = P7R_XMX(rmx, i, p7R_N) + gm->xsc[p7P_N][p7P_MOVE];
  path[1] = P7R_XMX(rmx, i, p7R_J) + gm->xsc[p7P_J][p7P_MOVE];
  return ( (path[0] > path[1]) ? p7T_N : p7T_J);
}

static inline int
select_c(const P7_PROFILE *gm, const P7_REFMX *rmx, int i)
{
  float path[2];
  path[0] = P7R_XMX(rmx, i-1, p7R_C) + gm->xsc[p7P_C][p7P_LOOP];
  path[1] = P7R_XMX(rmx, i,   p7R_E) + gm->xsc[p7P_E][p7P_MOVE];
  return ( (path[0] > path[1]) ? p7T_C : p7T_E);
}

static int
reference_viterbi_traceback(const P7_PROFILE *gm, const P7_REFMX *rmx, P7_TRACE *tr)
{
  int   i    = rmx->L;
  int   k    = 0;
  int   sprv = p7T_C;
  int   scur;
  int   status;

  if ((status = p7_trace_Append(tr, p7T_T, k, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, k, i)) != eslOK) return status;

  while (sprv != p7T_S)
    {
      switch (sprv) {
      case p7T_ML: scur = select_ml(gm, rmx, i, k); k--; i--; break;
      case p7T_MG: scur = select_mg(gm, rmx, i, k); k--; i--; break;
      case p7T_IL: scur = select_il(gm, rmx, i, k);      i--; break;
      case p7T_IG: scur = select_ig(gm, rmx, i, k);      i--; break;
      case p7T_DL: scur = select_dl(gm, rmx, i, k); k--;      break;
      case p7T_DG: scur = select_dg(gm, rmx, i, k); k--;      break;
      case p7T_E:  scur = select_e (gm, rmx, i, &k);          break;
      case p7T_N:  scur = select_n (         i    );          break;
      case p7T_J:  scur = select_j (gm, rmx, i);              break;
      case p7T_B:  scur = select_b (gm, rmx, i);              break;
      case p7T_L:  scur = p7T_B;                              break;
      case p7T_G:  scur = p7T_B;                              break;
      case p7T_C:  scur = select_c (gm, rmx, i);              break;
      default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
      }

      /* A glocal B->G->Mk wing-retraction entry: unfold it */
      if (scur == p7T_G) {
	while (k > 1) {
	  if ( (status = p7_trace_Append(tr, p7T_DG, k-1, i)) != eslOK) return status;
	  k--;
	}
      }

      if ( (status = p7_trace_Append(tr, scur, k, i)) != eslOK) return status;

      /* For NCJ, we had to defer i decrement. */
      if ( (scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;
      sprv = scur;
    }
  
  tr->M = rmx->M;
  tr->L = rmx->L;
  return p7_trace_Reverse(tr);
}
/*-------------- end, viterbi traceback -------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/


