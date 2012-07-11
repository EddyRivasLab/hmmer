/* 
 * 
 * Contents:
 *   1. Sparse Forward
 *   2. Sparse Backward
 *   3. Sparse Viterbi
 *   4. Sparse Viterbi trace
 *   5. Benchmark driver
 *   6. Unit tests
 *   7. Test driver
 *   8. Example 
 *   9. Notes
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_sparsemx.h"


/*****************************************************************
 * 1. Sparse Forward
 *****************************************************************/

int
p7_SparseForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_SPARSEMX *sx, float *opt_sc)
{
  P7_SPARSEMASK *sm   = sx->sm;
  float         *xpc  = sx->xmx;	 /* ptr that steps through current special cells   */
  float         *dpc  = sx->dp;	         /* ptr to step thru current row i main DP cells */
  float         *dpp;			 /* ptr to step thru previous row i-1 main DP cells */
  float         *last_dpc;		 /* used to reinit dpp after each sparse row computation */
  float const *tsc    = gm->tsc;	 /* sets up TSC() macro, access to profile's transitions */
  float const *rsc;			 /* will be set up for MSC(), ISC() macros for residue scores */
  int          ng;
  float        xE, xN, xJ, xB, xL, xG, xC;  /* tmp scores on special states. only stored when in row bands, and on ia-1 before a seg */
  float        mlc, mgc;		 /* temporary score calculations M(i,k)         */
  float        dlc, dgc;		 /* precalculated D(i,k+1) value on current row */
  int         *kc = sm->k[0];		 /* <kc> points to the list of sparse cell indices k for current row i */
  int         *kp;			 /* <kp> points to the previous row's sparse cell index list */
  int          i,k;	      	         /* i,k row,col (seq position, profile position) cell coords */
  int          y,z;			 /* indices in lists of k coords on prev, current row */


#ifdef P7_DEBUGGING  
  if (L != sx->L) ESL_EXCEPTION("L, sx->L disagree: sparse matrix wasn't allocated or reinitialized for this sequence");
#endif

  xN = 0.0f;
  xJ = -eslINFINITY;
  xC = -eslINFINITY;
  ng = 0;
  for (i = 1; i <= L; i++)
    {
      if (! sm->n[i]) { ng++; continue; }   /* skip rows that have no included cells */

      /* Reinitialize and store specials for row ia-1 just outside sparsified segment */
      if (i == 1 || ng) {
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

      kp = kc;                /* last row we did becomes prev row now; ready to step through k indices of previous row's sparse cells */
      kc = sm->k[i];		/* ditto for current row i */

      dlc = dgc = xE = -eslINFINITY;
      for (z = sm->n[i]-1, y = sm->n[i-1]-1; z >= 0; z--) /* Iterate over the one or more sparse cells (i,k) that we calculate on this row. Remember, sparsemask is stored in reverse order! */
	{
	  k = kc[z]; /* next sparse cell to calculate: (i,k) */
	  
	  /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	  while (y >= 0 && kp[y] < k-1) { y--; dpp+=p7S_NSCELLS; }
	  mlc = xL  + TSC(p7P_LM, k-1);
	  mgc = xG  + TSC(p7P_GM, k-1);
	  if (y >= 0 && kp[y] == k-1) {
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
	  while (y >= 0 && kp[y] < k) { y--; dpp+=p7S_NSCELLS; }
	  if (y >= 0 && kp[y] == k) {
	    *dpc++ = ISC(k) + p7_FLogsum( dpp[p7R_ML] + TSC(p7P_MI,k),  dpp[p7R_IL] + TSC(p7P_II, k));
	    *dpc++ = ISC(k) + p7_FLogsum( dpp[p7R_MG] + TSC(p7P_MI,k),  dpp[p7R_IG] + TSC(p7P_II, k));
	  } else {
	    *dpc++ = -eslINFINITY;
	    *dpc++ = -eslINFINITY;
	  }
	    
	  /* local exit paths */
	  xE = p7_FLogsum(xE, p7_FLogsum(mlc, dlc));

	  /* delayed store of Dk; advance calculation of next D_k+1 */
	  *dpc++ = dlc;
	  *dpc++ = dgc;
	  if (z >= 1 && kc[z-1] == k+1) { /* is there a (i,k+1) cell to our right? */
	    dlc = p7_FLogsum( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
	    dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	  } else if (z == 0) {             /* last sparse cell on row? we need dgc to complete MGk->E path, below */
	    dlc = -eslINFINITY;
	    dgc = p7_FLogsum( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	  } else {
	    dlc = dgc = -eslINFINITY;
	  }
	}

      *xpc++ = xE = p7_FLogsum( xE, dgc + TSC(p7P_DGE, k));   // glocal exit path(s) added on, from last k cell D(k)->E
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

  sx->type = p7S_FORWARD;
  if (opt_sc != NULL) *opt_sc = xC + ( ng ? ng *  gm->xsc[p7P_C][p7P_LOOP] : 0.0f) + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}
/*--------------- end, sparse Forward  --------------------------*/


/*****************************************************************
 * 3. Sparse Viterbi
 *****************************************************************/
static int sparse_viterbi_traceback(const P7_PROFILE *gm, const P7_SPARSEMX *sx, P7_TRACE *tr);

int
p7_SparseViterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_SPARSEMX *sx, float *opt_sc, P7_TRACE *opt_tr)
{
  P7_SPARSEMASK *sm   = sx->sm;
  float         *xpc  = sx->xmx;	 /* ptr that steps through current special cells   */
  float         *dpc  = sx->dp;	         /* ptr to step thru current row i main DP cells */
  float         *dpp;			 /* ptr to step thru previous row i-1 main DP cells */
  float         *last_dpc;		 /* used to reinit dpp after each sparse row computation */
  float const   *tsc    = gm->tsc;	 /* sets up TSC() macro, access to profile's transitions */
  float const   *rsc;			 /* will be set up for MSC(), ISC() macros for residue scores */
  int            ng;
  float          xE, xN, xJ, xB, xL, xG, xC;  /* tmp scores on special states. only stored when in row bands, and on ia-1 before a seg */
  float          mlc, mgc;		 /* temporary score calculations M(i,k)         */
  float          dlc, dgc;		 /* precalculated D(i,k+1) value on current row */
  int           *kc = sm->k[0];		 /* <kc> points to the list of sparse cell indices k for current row i */
  int           *kp;			 /* <kp> points to the previous row's sparse cell index list */
  int            i,k;	      	         /* i,k row,col (seq position, profile position) cell coords */
  int            y,z;			 /* indices in lists of k coords on prev, current row */


#ifdef P7_DEBUGGING  
  if (L != sx->L) ESL_EXCEPTION("L, sx->L disagree: sparse matrix wasn't allocated or reinitialized for this sequence");
#endif

  xN = 0.0f;
  xJ = -eslINFINITY;
  xC = -eslINFINITY;
  ng = 0;
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
      for (z = sm->n[i]-1, y = sm->n[i-1]-1; z >= 0; z--) /* Iterate over the one or more sparse cells (i,k) that we calculate on this row. Remember, sparsemask is stored in reverse order! */
	{
	  k = kc[z]; /* next sparse cell to calculate: (i,k) */
	  
	  /* Try to find cell i-1,k-1; then compute M(i,k) from it */
	  while (y >= 0 && kp[y] < k-1) { y--; dpp+=p7S_NSCELLS; }
	  mlc = xL  + TSC(p7P_LM, k-1);
	  mgc = xG  + TSC(p7P_GM, k-1);
	  if (y >= 0 && kp[y] == k-1) {
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
	  while (y >= 0 && kp[y] < k) { y--; dpp+=p7S_NSCELLS; }
	  if (y >= 0 && kp[y] == k) {
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
	  if (z >= 1 && kc[z-1] == k+1) { /* is there a (i,k+1) cell to our right? */
	    dlc = ESL_MAX( mlc + TSC(p7P_MD, k), dlc + TSC(p7P_DD, k));
	    dgc = ESL_MAX( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	  } else if (z == 0) {             /* last sparse cell on row? we need dgc to complete MGk->E path, below */
	    dlc = -eslINFINITY;
	    dgc = ESL_MAX( mgc + TSC(p7P_MD, k), dgc + TSC(p7P_DD, k));
	  } else {
	    dlc = dgc = -eslINFINITY;
	  }
	}

      *xpc++ = xE = ESL_MAX( xE, dgc + TSC(p7P_DGE, k));   // glocal exit path(s) added on, from last k cell D(k)->E; dgc includes Mk exit as Mk->Dk+1->E
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

  sx->type = p7S_VITERBI;
  xC += ( ng ? ng *  gm->xsc[p7P_C][p7P_LOOP] : 0.0f) + gm->xsc[p7P_C][p7P_MOVE];

  if (opt_sc) *opt_sc = xC;
  if (opt_tr && xC != -eslINFINITY) return sparse_viterbi_traceback(gm, sx, opt_tr);
  else                              return eslOK;
}




/*-------------------- end, Viterbi -----------------------------*/

static inline int
select_ml(const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, float *xp, int *ret_z)
{
  int   y = np-1;
  while (y >= 0 && kp[y] != k-1) { y--; dpp += p7S_NSCELLS; } /* dpp[] is in normal order; kp is reversed */
  if (y >= 0 && kp[y] == k-1) 
    {
      int   state[4] = { p7T_ML, p7T_IL, p7T_DL, p7T_L };
      float path[4];

      path[0] = dpp[p7S_ML] + P7P_TSC(gm, k-1, p7P_MM);
      path[1] = dpp[p7S_IL] + P7P_TSC(gm, k-1, p7P_IM);
      path[2] = dpp[p7S_DL] + P7P_TSC(gm, k-1, p7P_DM);
      path[3] =   xp[p7B_L] + P7P_TSC(gm, k-1, p7P_LM);
      *ret_z = y;
      return state[esl_vec_FArgMax(path, 4)];
    }
  else { *ret_z = 0; return p7T_L; }
}
static inline int
select_mg(const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, float *xp, int *ret_z)
{
  int   y = np-1;
  while (y >= 0 && kp[y] != k-1) { y--; dpp += p7S_NSCELLS; } /* dpp[] is in normal order; kp is reversed */
  if (y >= 0 && kp[y] == k-1) 
    {
      int   state[4] = { p7T_MG, p7T_IG, p7T_DG, p7T_G };
      float path[4];

      path[0] = dpp[p7S_MG] + P7P_TSC(gm, k-1, p7P_MM);
      path[1] = dpp[p7S_IG] + P7P_TSC(gm, k-1, p7P_IM);
      path[2] = dpp[p7S_DG] + P7P_TSC(gm, k-1, p7P_DM);
      path[3] =   xp[p7B_G] + P7P_TSC(gm, k-1, p7P_LM);
      *ret_z = y;
      return state[esl_vec_FArgMax(path, 4)];
    }
  else { *ret_z = 0; return p7T_G; }
}
static inline int
select_il(const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, int *ret_z)
{
  float path[2];
  int y = np-1;

  while (kp[y] != k) { y--; dpp += p7S_NSCELLS; } /* a little brave; we know an appropriate sparse cell exists on prv row, else we couldn't reach I on cur */
  path[0] = dpp[p7S_ML] + P7P_TSC(gm, k, p7P_MI);
  path[1] = dpp[p7S_IL] + P7P_TSC(gm, k, p7P_II);
  *ret_z = y;
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_IL);
}
static inline int
select_ig(const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, int *ret_z)
{
  float path[2];
  int y = np-1;

  while (kp[y] != k) { y--; dpp += p7S_NSCELLS; } /* a little brave; we know an appropriate sparse cell exists on prv row, else we couldn't reach I on cur */
  path[0] = dpp[p7S_MG] + P7P_TSC(gm, k, p7P_MI);
  path[1] = dpp[p7S_IG] + P7P_TSC(gm, k, p7P_II);
  *ret_z = y;
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_IG);
}
static inline int
select_dl(const P7_PROFILE *gm, int k, const float *dpc)
{
  float path[2];
  path[0] = dpc[p7S_ML] + P7P_TSC(gm, k-1, p7P_MD); /* more bravery. fact that we're tracing back from DL means that sparse cell k-1 must exist in z+1 */
  path[1] = dpc[p7S_DL] + P7P_TSC(gm, k-1, p7P_DD);
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_DL);
}
static inline int
select_dg(const P7_PROFILE *gm, int k, const float *dpc)
{
  float path[2];
  path[0] = dpc[p7S_MG] + P7P_TSC(gm, k-1, p7P_MD); /* more bravery. fact that we're tracing back from DL means that sparse cell k-1 must exist in z+1 */
  path[1] = dpc[p7S_DG] + P7P_TSC(gm, k-1, p7P_DD);
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_DG);
}
static inline int
select_j(const P7_PROFILE *gm, const float *xc)
{
  float path[2];
  path[0] = *(xc-p7S_NXCELLS+p7S_J) + gm->xsc[p7P_J][p7P_LOOP]; /* i.e. xp[p7S_J] on prv row i-1. */
  path[1] = xc[p7S_E]               + gm->xsc[p7P_E][p7P_LOOP];
  return ( (path[0] > path[1]) ? p7T_J : p7T_E);
}
static inline int
select_c(const P7_PROFILE *gm, const float *xc)
{
  float path[2];
  path[0] = *(xc-p7S_NXCELLS+p7S_C) + gm->xsc[p7P_C][p7P_LOOP]; /* i.e. xp[p7S_C] on prv row i-1. */
  path[1] =   xc[p7S_E]             + gm->xsc[p7P_E][p7P_MOVE];
  return ( (path[0] > path[1]) ? p7T_C : p7T_E);
}
static inline int
select_e(const P7_PROFILE *gm, const float *dpp, int *kp, int np, int *ret_z)
{
  float max  = -eslINFINITY;
  int   smax = -1;
  int   zmax = -1;
  int   z;

  /* oi vey; suppose previous row stored cells k= 2, 7, 8, 11
   * dpp[] points to these in forward order:  [ 2 7 8 11 ]
   * but kp[0..np-1] is in REVERSE order [11 8 7 2] because we store SPARSEMASK reversed
   */
  for (z = np-1; z > 0; z--)
    {       /* don't need to check DL->E path; these can't occur in a Viterbi path */
      if (dpp[p7S_ML] >= max) { max = dpp[p7S_ML]; smax = p7T_ML; zmax = z; }
      dpp += p7S_NSCELLS;
    }
  /* now dpp[] is on the final sparse supercell on this row, and z=0; k=kp[0] */
  if (dpp[p7S_ML] >= max) { max = dpp[p7S_ML]; smax = p7T_ML; zmax = 0; }
  
  /* Glocal Mk,Dk->E: Mkb + t(MD,kb->kb+1) + t(Dkb+1->E) wing retraction 
   * remember DGE is stored off by one; TSC(gm,kb,DGE) is t(Dkb+1->E) wing retraction 
   * for this to work on boundary condition kb=M, requires TSC(gm,M,DGE) = TSC(gm,M,MD) = TSC(gm,M,DD) = 0.0 
   * for this to work on boundary condition kb=M-1, requires TSC(gm,M-1,DGE) = 0.0 
   * and those boundary conditions are enforced: see modelconfig.c 
   */
  if ( dpp[p7S_MG] + P7P_TSC(gm, kp[0], p7P_MD) + P7P_TSC(gm, kp[0], p7P_DGE) >= max) { max = dpp[p7S_MG]; smax = p7T_MG; zmax = 0; }
  if ( dpp[p7S_DG] + P7P_TSC(gm, kp[0], p7P_DD) + P7P_TSC(gm, kp[0], p7P_DGE) >= max) { max = dpp[p7S_DG]; smax = p7T_DG; zmax = 0; }
  *ret_z = zmax;
  return smax;
}
static inline int
select_b(const P7_PROFILE *gm, const float *xc)
{
  float path[2];
  path[0] = xc[p7S_J] + gm->xsc[p7P_J][p7P_MOVE];
  path[1] = xc[p7S_N] + gm->xsc[p7P_N][p7P_MOVE];
  return ( (path[0] > path[1]) ? p7T_J : p7T_N);
}


/* Function:  
 * Synopsis:  
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEINVAL> if trace object isn't suitable, like if it
 *            was already used and not reinitialized (with _Reuse()).
 */
static int
sparse_viterbi_traceback(const P7_PROFILE *gm, const P7_SPARSEMX *sx, P7_TRACE *tr)
{
  P7_SPARSEMASK *sm = sx->sm;
  int            k  = 0;	/* current coord in profile consensus */
  int            i  = sm->L;	/* current coord in sequence (that snxt is on) */
  float         *dp;		/* points to main model DP cells for next valid row i */
  int            ip = sm->L;    /* current coord that <dp> is on (cur or prev row) */
  float         *xc;		/* points to xmx[i] special cells for current row (if it was stored; or prev stored row <i) */
  int            xc_on_i;       /* TRUE if <xc> is on row i; FALSE if xc points at a special row < i */
  int            scur, snxt;
  int            k2;		/* extra k counter while extending wings */
  int            z;		/* current position in n[i] sparse k array entries on current row dp[] */
  int            status;

#ifdef p7_DEBUGGING
  if (tr->N) ESL_EXCEPTION(eslEINVAL, "trace isn't empty - forgot to Reuse()?");
#endif

  /* <dp> points to the main sparse row we're tracing to, when we can
   * trace to an <snxt> in the model, and <ip> is the index of that
   * row; except for initiation conditions, ip is either i or i-1.
   * <dp> decrements by a row (i.e. by n[ip-1] supercells) when <ip>
   * decrements, which is when <snxt> accounts for x_i (i.e. if snxt =
   * M,I or an NN/CC/JJ emit.  Main model cells are stored for any row
   * i with n[i]>0.
   * 
   * <xc> points to the current row i, or an earlier row <i. xc_on_i
   * is TRUE when <xc> is on a stored special row i. Specials are
   * stored not only on rows with n[i]>0, but also on the row ia-1
   * immediately preceding a segment of rows i=ia..ib all with n[i]>0.
   * <xc> decrements whenever i decrements, which is when <scur> was
   * an M or I, or on an NN/CC/JJ emit.
   * 
   * When we emit NN/CC/JJ, we decrement both i and ip, but the
   * i decrement has to be deferred until after we attach the 
   * <snxt> state to the trace, because it's explaining i, not i-1.
   *
   * We build the trace backwards, and reverse it when we're done.
   * Remember, trace_Append() will filter k,i coords appropriately, only storing
   * them where it makes sense for the particular state, so it's harmless to
   * always pass current k,i even for nonemitting states.
   */

  xc      = sx->xmx + (sm->nrow+sm->nseg-1)*p7S_NXCELLS; /* initialized to last stored row, an end-of-seg ib; may be <L */
  xc_on_i = (sm->n[sm->L] ? TRUE : FALSE);		 /* if last row is in segment, stored, ib==L, then xc points to stored special row */

  dp      = sx->dp  + (sm->ncells - sm->n[sm->L]) * p7S_NSCELLS; /* <dp> is initialized on row ip=L, which might be empty */
      
  if ((status = p7_trace_Append(tr, p7T_T, k, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, k, i)) != eslOK) return status;
  scur = p7T_C;
  while (scur != p7T_S)
    {
      switch (scur) {
      case p7T_ML: snxt = select_ml(gm, k, dp, sm->k[i-1], sm->n[i-1], xc-p7S_NXCELLS, &z); i--; k--;      xc -= p7S_NXCELLS; break;
      case p7T_MG: snxt = select_mg(gm, k, dp, sm->k[i-1], sm->n[i-1], xc-p7S_NXCELLS, &z); i--; k--;      xc -= p7S_NXCELLS; break;
      case p7T_IL: snxt = select_il(gm, k, dp, sm->k[i-1], sm->n[i-1],                 &z); i--;           xc -= p7S_NXCELLS; break;
      case p7T_IG: snxt = select_ig(gm, k, dp, sm->k[i-1], sm->n[i-1],                 &z); i--;           xc -= p7S_NXCELLS; break;
      case p7T_DL: snxt = select_dl(gm, k, dp + (sm->n[i]-z-2)*p7S_NSCELLS);                     k--; z++;                    break; // something is wrong here.
      case p7T_DG: snxt = select_dg(gm, k, dp + (sm->n[i]-z-2)*p7S_NSCELLS);                     k--; z++;                    break;
      case p7T_N:  snxt =    (i == 0 ? p7T_S : p7T_N);                                                                        break;
      case p7T_J:  snxt = (sm->n[i] ? select_j(gm, xc) : p7T_J);                                                              break; 
      case p7T_C:  snxt = (sm->n[i] ? select_c(gm, xc) : p7T_C);                                                              break; // connect to E(i), C(i-1). E(i) valid if n[i]>0; and if E(i) valid, C(i-1) must be stored too. if E(i) invalid, connect to C, and it doesn't matter if xc is valid or not
      case p7T_E:  snxt = select_e(gm, dp, sm->k[i], sm->n[i], &z);                              k=sm->k[i][z];               break;
      case p7T_B:  snxt = select_b(gm, xc);                                                                                   break; // {NJ}(i) -> B(i). If we reached B(i), xc[i] valid, so NJ must also be valid.
      case p7T_L:  snxt = p7T_B;                                                                                              break;
      case p7T_G:  snxt = p7T_B;                                                                                              break;
      default:     ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
      }

      if (snxt == p7T_ML || snxt == p7T_MG || snxt == p7T_IL || snxt == p7T_IG) 
	{ ip--; dp -= sm->n[ip] * p7S_NSCELLS; }
	
      /* Glocal B->G->Mk left wing retraction entry: unfold it */
      if (snxt == p7T_G) 
	for (; k >= 1; k--) 
	  if ( (status = p7_trace_Append(tr, p7T_DG, k, i)) != eslOK) return status;
      /* Glocal Mk->E right wing retraction: off last sparse cell k, Mk->Dk+1->E or Dk->Dk+1->E */
      if (scur == p7T_E && (snxt == p7T_MG || snxt == p7T_DG))
	for (k2 = gm->M; k2 > k; k2--) 
	  if ( (status = p7_trace_Append(tr, p7T_DG, k2, i)) != eslOK) return status;
      
      /* Append the <snxt> state */
      if ( (status = p7_trace_Append(tr, snxt, k, i)) != eslOK) return status;

      /* NN,CC,JJ have a deferred i decrement, because of emit-on-transition 
       * note that the only way to leave a segment is via a JJ or NN, so
       * the only place we need to set xc_on_i FALSE is here.
       */
      if ( (snxt == p7T_N || snxt == p7T_J || snxt == p7T_C) && scur == snxt) 
	{
	  i--;
	  if (xc_on_i) xc -= p7S_NXCELLS;
	  xc_on_i = (sm->n[i] || sm->n[i+1]) ? TRUE : FALSE; 

	  ip--;
	  dp -= sm->n[ip] * p7S_NSCELLS; 
	}

      scur = snxt;
    }

  tr->M = sm->M;
  tr->L = sm->L;
  return p7_trace_Reverse(tr);
}




/*****************************************************************
 * 5. Benchmark driver
 *****************************************************************/

#ifdef p7SPARSE_FWDBACK_BENCHMARK
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "p7_sparsemx.h"

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
  P7_SPARSEMX    *sxf     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  /* Initialize log-sum calculator */
  impl_Init();
  p7_FLogsumInit();

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, L);

  sm  = p7_sparsemask_Create(gm->M, L, 0.0);
  p7_sparsemask_AddAll(sm);
  sxf = p7_sparsemx_Create(sm);

  /* Baseline time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_SparseForward  (dsq, L, gm, sxf, &sc);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_sparsemx_Destroy(sxf);
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
 * 6. Unit tests
 *****************************************************************/
#ifdef p7SPARSE_FWDBACK_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_randomseq.h"

static void
sparsemask_set_from_trace(ESL_RANDOMNESS *rng, P7_SPARSEMASK *sm, P7_TRACE *tr)
{
  char  msg[] = "sparse fwdback, sparse mask creation from trace failed";
  float cellprob = 0.5;
  float rowprob  = 0.2;
  int   i,k,z;
  int   status;
  
  z = tr->N-1;
  for (i = sm->L; i >= 1; i--) /* sparsemask api requires building it backwards */
    {
      while (tr->i[z] != i) z--; /* find trace position that generated this residue. */
    
      if ( (status = p7_sparsemask_StartRow(sm, i)) != eslOK) esl_fatal(msg);

      /* If this residue was emitted by the model, at least that cell
       * must be present; thus the row must be present. 
       * Tricky: in addition to the actual emitting cell i,k, we may
       * also need to add one or more delete cells i,k-1... 
       */
      if (p7_trace_IsM(tr->st[z]) || p7_trace_IsI(tr->st[z])) 
	{
	  while (p7_trace_IsD(tr->st[z+1])) z++;
	  
	  for (k = sm->M; k > tr->k[z]; k--) 
	    if (esl_random(rng) < cellprob)
	      if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) esl_fatal(msg); 

	  while (p7_trace_IsD(tr->st[z])) {
	    k = tr->k[z]; 
	    if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) esl_fatal(msg); 
	    z--;
	  }

	  k = tr->k[z];
	  if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) esl_fatal(msg); 
	  
	  for (k = k-1; k >= 1; k--)
	    if (esl_random(rng) < cellprob)
	      if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) esl_fatal(msg); 
	}
      else
	{
	  if (esl_random(rng) < rowprob)
	    for (k = sm->M; k >= 1; k--)
	      if (esl_random(rng) < cellprob)
		if ((status = p7_sparsemask_Add(sm, (k-1)%sm->Q, (k-1)/sm->Q)) != eslOK) esl_fatal(msg); /* append to k[i] list, increment n[i] count, reallocating as needed; doesn't deal w/ segments (nrow,nseg,i[]) */
	}

      if ((status = p7_sparsemask_FinishRow(sm)) != eslOK) esl_fatal(msg);
    }
  if ( (status = p7_sparsemask_Finish(sm)) != eslOK) esl_fatal(msg);
  return;
}




static void
utest_compare_reference(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char           msg[]  = "sparse fwdback, reference comparison unit test failed";
  P7_HMM        *hmm    = NULL;
  P7_PROFILE    *gm     = p7_profile_Create(M, abc);
  ESL_SQ        *sq     = esl_sq_CreateDigital(abc);       /* space for generated (homologous) target seqs              */
  ESL_DSQ       *dsqmem = malloc(sizeof(ESL_DSQ) * (L+2)); /* space for randomly generated target sequences of length L */
  ESL_DSQ       *dsq    = NULL;				 /* ptr into either dsqmem (nonhomologous targets) or sq->dsq (homologous targets) */
  P7_SPARSEMASK *sm     = p7_sparsemask_Create(M, L, 0.0);
  P7_SPARSEMX   *sx     = p7_sparsemx_Create(sm);
  P7_REFMX      *rmx    = p7_refmx_Create(M, L);
  P7_REFMX      *sxcopy = p7_refmx_Create(M, L);
  int            tL;
  int            idx;
  float          fsc1, fsc2;
  float          tol   =  (p7_logsum_IsSlowExact() ? 1e-5 : 0.01);
  
  if ( p7_hmm_Sample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)   != eslOK) esl_fatal(msg);

  for (idx = 0; idx < N; idx++)
    {
      if ( p7_profile_SetLength(gm, L)  != eslOK) esl_fatal(msg);   /* config to generate mean length of L (length may have been reset to tL by last emitted seq) */

      /* A 50:50 mix of generated (homologous) and random (nonhomologous) sequences */
      if (esl_rnd_Roll(rng, 2)) 
	{ /* fixed-length random emission */
	  esl_rsq_xfIID(rng, bg->f, abc->K, L, dsqmem);  
	  dsq = dsqmem;  
	  tL = L;     
	}
      else  
	{ /* variable-length seq emission: length config on <gm>,<om> will change to <tL> before DP */
	  do {
	    esl_sq_Reuse(sq);
	    p7_ProfileEmit(rng, hmm, gm, bg, sq, NULL);
	  } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
	  dsq = sq->dsq; 
	  tL = sq->n; 
	}

      if ( p7_profile_SetLength(gm, tL)         != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_Reinit(sm, M, tL, 0.0) != eslOK) esl_fatal(msg);
      if ( p7_sparsemask_AddAll(sm)             != eslOK) esl_fatal(msg);

      if ( p7_refmx_GrowTo(rmx,    M, tL)  != eslOK) esl_fatal(msg);
      if ( p7_refmx_GrowTo(sxcopy, M, tL)  != eslOK) esl_fatal(msg);
      if ( p7_sparsemx_Reinit(sx, sm)      != eslOK) esl_fatal(msg);

      p7_SparseForward   (dsq, tL, gm, sx,  &fsc1);
      p7_ReferenceForward(dsq, tL, gm, rmx, &fsc2);

      p7_sparsemx_Copy2Reference(sx, sxcopy);

      //printf("%3d %.4f %.4f\n", idx, fsc1, fsc2);
      //p7_sparsemx_Dump(stdout, sx);
      //p7_refmx_Dump(stdout, sxcopy);
      //p7_refmx_Dump(stdout, rmx);

      if ( fabs(fsc1-fsc2)                    > tol)    esl_fatal(msg);
      if ( p7_refmx_Validate(sxcopy, NULL)    != eslOK) esl_fatal(msg);
      if ( p7_refmx_Compare(rmx, sxcopy, tol) != eslOK) esl_fatal(msg);
      
      esl_sq_Reuse(sq);
      p7_refmx_Reuse(rmx);
      p7_refmx_Reuse(sxcopy);
      p7_sparsemask_Reuse(sm);
      p7_sparsemx_Reuse(sx);
    }

  free(dsqmem);
  p7_sparsemask_Destroy(sm);
  p7_sparsemx_Destroy(sx);
  p7_refmx_Destroy(rmx);
  p7_refmx_Destroy(sxcopy);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}

static void
utest_singlepath(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int N)
{
  char           msg[] = "sparse fwdback singlepath unit test failed";
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = p7_profile_Create(M, abc);
  ESL_SQ        *sq    = esl_sq_CreateDigital(abc);
  P7_TRACE      *gtr   = p7_trace_Create();           /* generated trace */
  P7_TRACE      *vtr   = p7_trace_Create();	      /* viterbi trace */
  P7_SPARSEMASK *sm    = NULL;
  P7_SPARSEMX   *sxf   = NULL;
  P7_SPARSEMX   *sxv   = NULL;
  float          tsc, fsc, vsc;
  float          tol   = 1e-4;
  int            idx;

  for (idx = 0; idx < N; idx++)
    {
      /* Create a profile that has only a single possible path (including
       * emissions) thru it; requires configuring in uniglocal mode w/ L=0
       */
      if ( p7_hmm_SampleSinglePathed(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
      if ( p7_profile_ConfigUniglocal(gm, hmm, bg, 0)   != eslOK) esl_fatal(msg);

      /* Sample that sequence and path; get its trace score */
      if ( p7_ProfileEmit(rng, hmm, gm, bg, sq, gtr)    != eslOK) esl_fatal(msg);
      if ( p7_trace_Score(gtr, sq->dsq, gm, &tsc)       != eslOK) esl_fatal(msg);
  
      //p7_trace_DumpAnnotated(stdout, gtr, gm, sq->dsq);

      /* Build a randomized sparse mask around that trace */
      if  (sm) { if (   p7_sparsemask_Reinit(sm, M, sq->n, 0.0)  != eslOK) esl_fatal(msg); }
      else     { if ( (sm = p7_sparsemask_Create(M, sq->n, 0.0)) == NULL) esl_fatal(msg); }
      sparsemask_set_from_trace(rng, sm, gtr);

      if  (sxf) { if (   p7_sparsemx_Reinit(sxf, sm) != eslOK) esl_fatal(msg); }
      else      { if ( (sxf = p7_sparsemx_Create(sm)) == NULL) esl_fatal(msg); }

      if  (sxv) { if (   p7_sparsemx_Reinit(sxv, sm) != eslOK) esl_fatal(msg); }
      else      { if ( (sxv = p7_sparsemx_Create(sm)) == NULL) esl_fatal(msg); }

      //p7_sparsemask_Dump(stdout, sm);

      /* Run DP routines, collect scores that should all match trace score */
      if ( p7_SparseForward(sq->dsq, sq->n, gm, sxf, &fsc)      != eslOK) esl_fatal(msg);
      if ( p7_SparseViterbi(sq->dsq, sq->n, gm, sxv, &vsc, vtr) != eslOK) esl_fatal(msg);
  
      //p7_sparsemx_Dump(stdout, sxv);
      p7_trace_DumpAnnotated(stdout, gtr, gm, sq->dsq);
      p7_trace_DumpAnnotated(stdout, vtr, gm, sq->dsq);

      /* Since only a single path is possible, trace score and Fwd score match */
      if ( esl_FCompareAbs(tsc, vsc, tol)   != eslOK) esl_fatal(msg);
      if ( esl_FCompareAbs(tsc, fsc, tol)   != eslOK) esl_fatal(msg);
      if ( p7_trace_Compare(gtr, vtr, 0.0f) != eslOK) esl_fatal(msg); // 0.0 is <pptol> arg, unused, because neither trace has PP annotation
  
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
  p7_sparsemask_Destroy(sm);
  p7_trace_Destroy(vtr);
  p7_trace_Destroy(gtr);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
}


#endif /*p7SPARSE_FWDBACK_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/


/*****************************************************************
 * 7. Test driver
 *****************************************************************/
#ifdef p7SPARSE_FWDBACK_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "p7_sparsemx.h"

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

  p7_FLogsumInit();
  impl_Init();

  utest_compare_reference(r, abc, bg, M, L, N);
  utest_singlepath       (r, abc, bg, M,    N);

  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}

#endif /*p7SPARSE_FWDBACK_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/



/*****************************************************************
 * 8. Example
 *****************************************************************/
#ifdef p7SPARSE_FWDBACK_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_regexp.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_sparsemx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",              0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Forward DP matrix for examination",            0 },
  { "-V",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Viterbi DP matrix for examination",            0 },
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Viterbi trace for examination",                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of Forward/Backward, sparse dual implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
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
  P7_SPARSEMASK  *sm      = NULL;
  P7_SPARSEMX    *sxf     = NULL;
  P7_SPARSEMX    *sxv     = NULL;
  P7_TRACE       *tr      = p7_trace_Create();
  float           fsc, vsc;
  float           nullsc;
  int             status;

  /* Initialize log-sum calculator */
  impl_Init();
  p7_FLogsumInit();

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
  p7_profile_Config(gm, hmm, bg);

  /* Allocate bands, matrices */
  sm  = p7_sparsemask_Create(gm->M, sq->n, 0.0);
  p7_sparsemask_AddAll(sm);
  sxf = p7_sparsemx_Create(sm);
  sxv = p7_sparsemx_Create(sm);

  /* Set the profile and null model's target length models */
  p7_bg_SetLength           (bg, sq->n);
  p7_profile_SetLength      (gm, sq->n);

  /* Sparse forward calculation */
  p7_SparseForward (sq->dsq, sq->n, gm, sxf, &fsc);
  p7_SparseViterbi (sq->dsq, sq->n, gm, sxf, &vsc, tr);
  p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

 if (esl_opt_GetBoolean(go, "-F")) p7_sparsemx_Dump(stdout, sxf);
 if (esl_opt_GetBoolean(go, "-V")) p7_sparsemx_Dump(stdout, sxv);
 if (esl_opt_GetBoolean(go, "-T")) p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);


  printf("target sequence:      %s\n",         sq->name);
  printf("vit raw score:        %.4f nats\n",  vsc);
  printf("fwd raw score:        %.4f nats\n",  fsc);
  printf("null score:           %.2f nats\n",  nullsc);
  printf("per-seq score:        %.2f bits\n",  (fsc - nullsc) / eslCONST_LOG2);

  /* Cleanup */
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxv);
  p7_sparsemask_Destroy(sm);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_FWDBACK_EXAMPLE*/
/*------------------ end, example driver ------------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
