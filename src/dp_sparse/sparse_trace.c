/* Traceback routines, both Viterbi and stochastic.
 * The traceback machinery for V and Sto traces is the same; only the
 * select*() choice functions differ. 
 * 
 * Contents:
 *   1. Choice selection functions for Viterbi optimal traces
 *   2. Choice selection functions for stochastic traces of a Forward matrix
 *   3. Traceback engine, shared by Viterbi and stochastic
 *   4. Exposed API, wrappers around the engine
 *   5. Example
 *   6. Copyright and license information.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "base/p7_profile.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_trace.h"

/*****************************************************************
 * 1. Choice selection functions for Viterbi traces
 *****************************************************************/

static inline int
v_select_ml(ESL_RANDOMNESS *r, const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, float *xp, int *ret_z)
{
  ESL_UNUSED(r);
  int   y = 0;
  while (y < np && kp[y]  < k-1) { y++; dpp += p7S_NSCELLS; } 
  if    (y < np && kp[y] == k-1) 
    {
      int   state[4] = { p7T_ML, p7T_IL, p7T_DL, p7T_L };
      float path[4];

      path[0] = dpp[p7S_ML] + P7P_TSC(gm, k-1, p7P_MM);
      path[1] = dpp[p7S_IL] + P7P_TSC(gm, k-1, p7P_IM);
      path[2] = dpp[p7S_DL] + P7P_TSC(gm, k-1, p7P_DM);
      path[3] =   xp[p7S_L] + P7P_TSC(gm, k-1, p7P_LM);
      *ret_z = y;
      return state[esl_vec_FArgMax(path, 4)];
    }
  else { *ret_z = 0; return p7T_L; }
}
static inline int
v_select_mg(ESL_RANDOMNESS *r, const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, float *xp, int *ret_z)
{
  ESL_UNUSED(r);
  int   y = 0;
  while (y < np && kp[y]  < k-1) { y++; dpp += p7S_NSCELLS; } 
  if    (y < np && kp[y] == k-1) 
    {
      int   state[4] = { p7T_MG, p7T_IG, p7T_DG, p7T_G };
      float path[4];

      path[0] = dpp[p7S_MG] + P7P_TSC(gm, k-1, p7P_MM);
      path[1] = dpp[p7S_IG] + P7P_TSC(gm, k-1, p7P_IM);
      path[2] = dpp[p7S_DG] + P7P_TSC(gm, k-1, p7P_DM);
      path[3] =   xp[p7S_G] + P7P_TSC(gm, k-1, p7P_GM);
      *ret_z = y;
      return state[esl_vec_FArgMax(path, 4)];
    }
  else { *ret_z = 0; return p7T_G; }
}
static inline int
v_select_il(ESL_RANDOMNESS *r, const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, int *ret_z)
{
  ESL_UNUSED(r);
  float path[2];
  int   y = 0;
  while (kp[y] != k) { y++; dpp += p7S_NSCELLS; } /* a little brave; we know an appropriate sparse cell exists on prv row, else we couldn't reach I on cur */
  path[0] = dpp[p7S_ML] + P7P_TSC(gm, k, p7P_MI);
  path[1] = dpp[p7S_IL] + P7P_TSC(gm, k, p7P_II);
  *ret_z = y;
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_IL);
}
static inline int
v_select_ig(ESL_RANDOMNESS *r, const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, int *ret_z)
{
  ESL_UNUSED(r);
  float path[2];
  int y = 0;
  while (kp[y] != k) { y++; dpp += p7S_NSCELLS; } /* a little brave; we know an appropriate sparse cell exists on prv row, else we couldn't reach I on cur */
  path[0] = dpp[p7S_MG] + P7P_TSC(gm, k, p7P_MI);
  path[1] = dpp[p7S_IG] + P7P_TSC(gm, k, p7P_II);
  *ret_z = y;
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_IG);
}
static inline int
v_select_dl(ESL_RANDOMNESS *r, const P7_PROFILE *gm, int k, const float *dpp)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = dpp[p7S_ML] + P7P_TSC(gm, k-1, p7P_MD); /* more bravery. fact that we're tracing back from DL means that sparse cell k-1 must exist in z+1 */
  path[1] = dpp[p7S_DL] + P7P_TSC(gm, k-1, p7P_DD);
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_DL);
}
static inline int
v_select_dg(ESL_RANDOMNESS *r, const P7_PROFILE *gm, int k, const float *dpp)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = dpp[p7S_MG] + P7P_TSC(gm, k-1, p7P_MD); /* more bravery. fact that we're tracing back from DL means that sparse cell k-1 must exist in z+1 */
  path[1] = dpp[p7S_DG] + P7P_TSC(gm, k-1, p7P_DD);
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_DG);
}
static inline int
v_select_j(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const float *xc)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = *(xc-p7S_NXCELLS+p7S_J) + gm->xsc[p7P_J][p7P_LOOP]; /* i.e. xp[p7S_J] on prv row i-1. */
  path[1] = xc[p7S_E]               + gm->xsc[p7P_E][p7P_LOOP];
  return ( (path[0] > path[1]) ? p7T_J : p7T_E);
}
static inline int
v_select_c(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const float *xc)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = *(xc-p7S_NXCELLS+p7S_C) + gm->xsc[p7P_C][p7P_LOOP]; /* i.e. xp[p7S_C] on prv row i-1. */
  path[1] =   xc[p7S_E]             + gm->xsc[p7P_E][p7P_MOVE];
  return ( (path[0] > path[1]) ? p7T_C : p7T_E);
}
static inline int
v_select_e(ESL_RANDOMNESS *r, float *wrk, const P7_PROFILE *gm, const float *dpp, int *kp, int np, int *ret_z)
{
  ESL_UNUSED(r);
  ESL_UNUSED(wrk);
  float max  = -eslINFINITY;
  float pathsc;
  int   smax = -1;
  int   zmax = -1;
  int   z;

  /* Glocal Mk,Dk->E: Mkb + t(MD,kb->kb+1) + t(Dkb+1->E) wing retraction 
   * remember DGE is stored off by one; TSC(gm,kb,DGE) is t(Dkb+1->E) wing retraction 
   * for this to work on boundary condition kb=M, requires TSC(gm,M,DGE) = TSC(gm,M,MD) = TSC(gm,M,DD) = 0.0 
   * for this to work on boundary condition kb=M-1, requires TSC(gm,M-1,DGE) = 0.0 
   * and those boundary conditions are enforced: see modelconfig.c 
   * 
   * Dk->E paths don't need to be checked in Viterbi;
   *   in local, any Mk->Dk+1->...E path must be inferior to Mk->E
   *   in glocal, an Mk->Dk+1->...E path can win, but all such paths are already checked by wing-retracted Mk->...E exits
   */
  for (z = 0; z < np-1; z++)
    {       /* don't need to check DL->E path; these can't occur in a Viterbi path */
      if (dpp[p7S_ML] >= max) { max = dpp[p7S_ML]; smax = p7T_ML; zmax = z; }
      
      if (kp[z+1] != kp[z]+1) 
	{ 	/* sparse cell k with no following cell k+1: then a glocal exit Mk->...->E across unmarked cells is possible */
	  pathsc = dpp[p7S_MG] + P7P_TSC(gm, kp[z], p7P_MD) + P7P_TSC(gm, kp[z], p7P_DGE);
	  if (pathsc > max) { max = pathsc; smax = p7T_MG; zmax = z; }
	}
      dpp += p7S_NSCELLS;
    }
  /* last cell np-1 is out of loop because we don't want to bump dpp after it */
  if (dpp[p7S_ML] > max) { max = dpp[p7S_ML]; smax = p7T_ML; zmax = z; }
  pathsc = dpp[p7S_MG] + P7P_TSC(gm, kp[z], p7P_MD) + P7P_TSC(gm, kp[z], p7P_DGE);  if (pathsc > max) { max = pathsc; smax = p7T_MG; zmax = z; }  
  pathsc = dpp[p7S_DG] + P7P_TSC(gm, kp[z], p7P_DD) + P7P_TSC(gm, kp[z], p7P_DGE);  if (pathsc > max) { max = pathsc; smax = p7T_DG; zmax = z; }  
  *ret_z = zmax;
  return smax;
}
static inline int
v_select_b(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const float *xc)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = xc[p7S_J] + gm->xsc[p7P_J][p7P_MOVE];
  path[1] = xc[p7S_N] + gm->xsc[p7P_N][p7P_MOVE];
  return ( (path[0] > path[1]) ? p7T_J : p7T_N);
}
/*-------------- end, Viterbi choice functions ------------------*/



/*****************************************************************
 * 2. Choice selection functions for stochastic trace of Fwd matrix
 *****************************************************************/
/* These are essentially copies of the Viterbi choice functions,
 * except that instead of simply choosing the maximum path, we
 * renormalize the log probs in path[] (with <esl_vec_FLogNorm()>)
 * to create a probability vector, and sample from that (with
 * <esl_rnd_FChoose()>.
 */
static inline int
sto_select_ml(ESL_RANDOMNESS *rng, const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, float *xp, int *ret_z)
{
  int   y = 0;
  while (y < np && kp[y]  < k-1) { y++; dpp += p7S_NSCELLS; } 
  if    (y < np && kp[y] == k-1) 
    {
      int   state[4] = { p7T_ML, p7T_IL, p7T_DL, p7T_L };
      float path[4];

      path[0] = dpp[p7S_ML] + P7P_TSC(gm, k-1, p7P_MM);
      path[1] = dpp[p7S_IL] + P7P_TSC(gm, k-1, p7P_IM);
      path[2] = dpp[p7S_DL] + P7P_TSC(gm, k-1, p7P_DM);
      path[3] =   xp[p7S_L] + P7P_TSC(gm, k-1, p7P_LM);
      esl_vec_FLogNorm(path, 4);
      *ret_z = y;
      return state[esl_rnd_FChoose(rng, path, 4)];
    }
  else { *ret_z = 0; return p7T_L; }
}
static inline int
sto_select_mg(ESL_RANDOMNESS *rng, const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, float *xp, int *ret_z)
{
  int   y = 0;
  while (y < np && kp[y]  < k-1) { y++; dpp += p7S_NSCELLS; } 
  if    (y < np && kp[y] == k-1) 
    {
      int   state[4] = { p7T_MG, p7T_IG, p7T_DG, p7T_G };
      float path[4];

      path[0] = dpp[p7S_MG] + P7P_TSC(gm, k-1, p7P_MM);
      path[1] = dpp[p7S_IG] + P7P_TSC(gm, k-1, p7P_IM);
      path[2] = dpp[p7S_DG] + P7P_TSC(gm, k-1, p7P_DM);
      path[3] =   xp[p7S_G] + P7P_TSC(gm, k-1, p7P_GM);
      esl_vec_FLogNorm(path, 4);
      *ret_z = y;
      return state[esl_rnd_FChoose(rng,path,4)];
    }
  else { *ret_z = 0; return p7T_G; }
}
static inline int
sto_select_il(ESL_RANDOMNESS *rng, const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, int *ret_z)
{
  float path[2];
  int   y = 0;
  while (kp[y] != k) { y++; dpp += p7S_NSCELLS; } /* a little brave; we know an appropriate sparse cell exists on prv row, else we couldn't reach I on cur */
  path[0] = dpp[p7S_ML] + P7P_TSC(gm, k, p7P_MI);
  path[1] = dpp[p7S_IL] + P7P_TSC(gm, k, p7P_II);
  esl_vec_FLogNorm(path, 2);
  *ret_z = y;
  return ( (esl_rnd_FChoose(rng,path,2) == 0) ? p7T_ML : p7T_IL);
}
static inline int
sto_select_ig(ESL_RANDOMNESS *rng, const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, int *ret_z)
{
  float path[2];
  int y = 0;
  while (kp[y] != k) { y++; dpp += p7S_NSCELLS; } /* a little brave; we know an appropriate sparse cell exists on prv row, else we couldn't reach I on cur */
  path[0] = dpp[p7S_MG] + P7P_TSC(gm, k, p7P_MI);
  path[1] = dpp[p7S_IG] + P7P_TSC(gm, k, p7P_II);
  esl_vec_FLogNorm(path, 2);
  *ret_z = y;
  return ( (esl_rnd_FChoose(rng,path,2) == 0) ? p7T_MG : p7T_IG);
}
static inline int
sto_select_dl(ESL_RANDOMNESS *rng, const P7_PROFILE *gm, int k, const float *dpp)
{
  float path[2];
  path[0] = dpp[p7S_ML] + P7P_TSC(gm, k-1, p7P_MD); /* more bravery. fact that we're tracing back from DL means that sparse cell k-1 must exist in z+1 */
  path[1] = dpp[p7S_DL] + P7P_TSC(gm, k-1, p7P_DD);
  esl_vec_FLogNorm(path, 2);
  return ( (esl_rnd_FChoose(rng,path,2) == 0) ? p7T_ML : p7T_DL);
}
static inline int
sto_select_dg(ESL_RANDOMNESS *rng, const P7_PROFILE *gm, int k, const float *dpp)
{
  float path[2];
  path[0] = dpp[p7S_MG] + P7P_TSC(gm, k-1, p7P_MD); /* more bravery. fact that we're tracing back from DL means that sparse cell k-1 must exist in z+1 */
  path[1] = dpp[p7S_DG] + P7P_TSC(gm, k-1, p7P_DD);
  esl_vec_FLogNorm(path, 2);
  return ( (esl_rnd_FChoose(rng,path,2) == 0) ? p7T_MG : p7T_DG);
}
static inline int
sto_select_j(ESL_RANDOMNESS *rng, const P7_PROFILE *gm, const float *xc)
{
  float path[2];
  path[0] = *(xc-p7S_NXCELLS+p7S_J) + gm->xsc[p7P_J][p7P_LOOP]; /* i.e. xp[p7S_J] on prv row i-1. */
  path[1] = xc[p7S_E]               + gm->xsc[p7P_E][p7P_LOOP];
  esl_vec_FLogNorm(path, 2);
  return ( (esl_rnd_FChoose(rng,path,2) == 0) ? p7T_J : p7T_E);
}
static inline int
sto_select_c(ESL_RANDOMNESS *rng, const P7_PROFILE *gm, const float *xc)
{
  float path[2];
  path[0] = *(xc-p7S_NXCELLS+p7S_C) + gm->xsc[p7P_C][p7P_LOOP]; /* i.e. xp[p7S_C] on prv row i-1. */
  path[1] =   xc[p7S_E]             + gm->xsc[p7P_E][p7P_MOVE];
  esl_vec_FLogNorm(path, 2);
  return ( (esl_rnd_FChoose(rng,path,2) == 0) ? p7T_C : p7T_E);
}
static inline int
sto_select_e(ESL_RANDOMNESS *rng, float *wrk, const P7_PROFILE *gm, const float *dpp, int *kp, int np, int *ret_z)
{
  /* The <wrk> vector is a scratch space with room for at least p7S_NSCELLS*np floats */
  int   w,z;
  int   state[6] = { p7T_ML, p7T_MG, p7T_BOGUS, p7T_BOGUS, p7T_DL, p7T_DG }; /* w%6 in profile states mapped to trace statetypes */

  for (z = 0; z < np; z++, dpp += p7S_NSCELLS)
    {
      wrk[z*p7S_NSCELLS + p7S_ML] = dpp[p7S_ML];
      wrk[z*p7S_NSCELLS + p7S_DL] = dpp[p7S_DL];
      
      if (z == np-1 || kp[z+1] != kp[z]+1) 
	{ /* sparse cell k with no following cell k+1: then a glocal exit {MD}k->...->E across unmarked cells is possible */      
	  wrk[z*p7S_NSCELLS + p7S_MG] = dpp[p7S_MG] + P7P_TSC(gm, kp[z], p7P_MD) + P7P_TSC(gm, kp[z], p7P_DGE);
	  wrk[z*p7S_NSCELLS + p7S_DG] = dpp[p7S_DG] + P7P_TSC(gm, kp[z], p7P_DD) + P7P_TSC(gm, kp[z], p7P_DGE);
	}
      else
	{ /* sparse cell with a following cell k+1; no glocal exit, instead prob flows along {MD}k->Dk+1 path  */
	  wrk[z*p7S_NSCELLS + p7S_MG] = -eslINFINITY;
	  wrk[z*p7S_NSCELLS + p7S_DG] = -eslINFINITY;
	}

      wrk[z*p7S_NSCELLS + p7S_IL] = -eslINFINITY;
      wrk[z*p7S_NSCELLS + p7S_IG] = -eslINFINITY;
    }

  /* wrk has p7S_NSCELLS*np values, some of which are -inf; normalize and choose one; then sort out z,s components of its index */
  esl_vec_FLogNorm(wrk, np*p7S_NSCELLS);
  w = esl_rnd_FChoose(rng,wrk,np*p7S_NSCELLS);
  *ret_z = w / p7S_NSCELLS;
  return state[w%p7S_NSCELLS];
}
static inline int
sto_select_b(ESL_RANDOMNESS *rng, const P7_PROFILE *gm, const float *xc)
{
  float path[2];
  path[0] = xc[p7S_J] + gm->xsc[p7P_J][p7P_MOVE];
  path[1] = xc[p7S_N] + gm->xsc[p7P_N][p7P_MOVE];
  esl_vec_FLogNorm(path, 2);
  return ( (esl_rnd_FChoose(rng,path,2) == 0) ? p7T_J : p7T_N);
}
/*------------- end, stochastic choice functions ----------------*/



/*****************************************************************
 * 3. Traceback engine, shared by V and Sto tracing
 *****************************************************************/

/* 
 * Can work either as a complete-sequence traceback of a forward
 * matrix, or a "segment trace" of a single segment. (A segment trace
 * is a hack used in the sparse MPAS algorithm; see sparse_anchors.c,
 * including footnote [4] in that .c file.) 
 * 
 * For a normal (complete-sequence) traceback: pass <g=0>. (And then
 * it doesn't matter what you pass for <dp>, <xc>, so you may as well
 * pass NULL for them.)
 * 
 * For a segment traceback, pass <g=1..G>, the segment number you want
 * to trace. Now <dp>, <xc> can be used as hints to speed up access to
 * the sparse matrix: if the caller knows a ptr to the start of the
 * main matrix row for the last row <ib> of the segment, pass that as
 * <dp>; if it knows a ptr to the start of the specials for the last
 * row of the segment g, pass that as <xc>. If these aren't known, no
 * problem, pass <dp> and <xc> as <NULL>, and they'll be determined by
 * a slower method.
 * 
 * Optionally, can return the number of domains in the trace in <*opt_D>.
 * 
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, 
 *            in <p7_trace_Append()> for example.
 */
static int
sparse_traceback_engine(ESL_RANDOMNESS *rng, float *wrk, const P7_PROFILE *gm, const P7_SPARSEMX *sx, P7_TRACE *tr,
			int g, float *dp, float *xc, int *opt_D)
{
  const P7_SPARSEMASK *sm = sx->sm;
  int            k  = 0;	/* current coord in profile consensus */
  int            i;	        /* current coord in sequence (that snxt is on) */
  int            ip;            /* current coord that <dp> is on (cur or prev row) */
  int            xc_on_i;       /* TRUE if <xc> is on row i; FALSE if xc points at a special row < i */
  int            scur, snxt;
  int            k2;		/* extra k counter while extending wings */
  int            z;		/* current position in n[i] sparse k array entries on current row dp[] */
  int            istart;        /* row we're tracing back to: 0 for complete, ia-1 for segment trace */
  int            D  = 0;	/* counts number of domains in the trace, by counting E states. */
  int            status;

  ESL_DASSERT1(( tr->N == 0 ));   // trace must be empty - _Reuse() was called on it

  int (*select_ml)(ESL_RANDOMNESS *rng,             const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, float *xp, int *ret_z);
  int (*select_mg)(ESL_RANDOMNESS *rng,             const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np, float *xp, int *ret_z);
  int (*select_il)(ESL_RANDOMNESS *rng,             const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np,            int *ret_z);
  int (*select_ig)(ESL_RANDOMNESS *rng,             const P7_PROFILE *gm, int k, const float *dpp, int *kp, int np,            int *ret_z);
  int (*select_dl)(ESL_RANDOMNESS *rng,             const P7_PROFILE *gm, int k, const float *dpp);
  int (*select_dg)(ESL_RANDOMNESS *rng,             const P7_PROFILE *gm, int k, const float *dpp);
  int (*select_j) (ESL_RANDOMNESS *rng,             const P7_PROFILE *gm, const float *xc);
  int (*select_c) (ESL_RANDOMNESS *rng,             const P7_PROFILE *gm, const float *xc);
  int (*select_e) (ESL_RANDOMNESS *rng, float *wrk, const P7_PROFILE *gm,        const float *dpp, int *kp, int np,            int *ret_z);
  int (*select_b) (ESL_RANDOMNESS *rng,             const P7_PROFILE *gm, const float *xc);
  
  if (sx->type == p7S_VITERBI)	/* configure for Viterbi optimal traceback */
    {
      select_ml = &v_select_ml;      select_mg = &v_select_mg;
      select_il = &v_select_il;      select_ig = &v_select_ig;
      select_dl = &v_select_dl;      select_dg = &v_select_dg;
      select_j  = &v_select_j;       select_c  = &v_select_c;
      select_e  = &v_select_e;       select_b  = &v_select_b;
    }
  else if (sx->type == p7S_FORWARD) /* configure for stochastic traceback */
    {
      select_ml = &sto_select_ml;      select_mg = &sto_select_mg;
      select_il = &sto_select_il;      select_ig = &sto_select_ig;
      select_dl = &sto_select_dl;      select_dg = &sto_select_dg;
      select_j  = &sto_select_j;       select_c  = &sto_select_c;
      select_e  = &sto_select_e;       select_b  = &sto_select_b;
    }



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
  if (g)  // segment trace version.
    {
      if (dp == NULL)
	{
	  dp = sx->dp;
	  for (i =  1; i < sm->seg[g].ia; i++)
	    dp += p7S_NSCELLS * sm->n[i];
	}
      if (xc == NULL)
	{
	  xc = sx->xmx;
	  int g2;
	  for (g2 = 1; g2 < g; g++)
	    xc += p7S_NXCELLS * (sm->seg[g].ib - sm->seg[g].ia+2);  // ia-1..ib = ib-ia+2 to skip segment
	}

      i = ip = sm->seg[g].ib;   
      xc_on_i = TRUE;   // we know ib is stored, it's in a segment.		
      if ((status = p7_trace_Append(tr, p7T_J, k, i)) != eslOK) return status;   // _Append() will actually set i[0] = 0, not i; there's always one mute J state.
      scur = p7T_J;                                                              // If we do use a J to generate ib, that gets handled and added by the main recursion.
      istart = sm->seg[g].ia;
    }
  else  // complete-sequence trace version
    {
      i = ip  = sm->L;
      xc      = sx->xmx + (sm->nrow+sm->S-1)*p7S_NXCELLS; /* initialized to last stored row, an end-of-seg ib; may be <L */
      xc_on_i = (sm->n[sm->L] ? TRUE : FALSE);	          /* if last row is in segment, stored, ib==L, then xc points to stored special row */

      dp      = sx->dp  + (sm->ncells - sm->n[sm->L]) * p7S_NSCELLS; /* <dp> is initialized on row ip=L, which might be empty */

      if ((status = p7_trace_Append(tr, p7T_T, k, i)) != eslOK) return status;
      if ((status = p7_trace_Append(tr, p7T_C, k, i)) != eslOK) return status;
      scur   = p7T_C;
      istart = 0;
    }

  /* Main recursion. */
  while ( i > istart || (scur != p7T_N && scur != p7T_J))    // Normally we trace to N(0). Segment trace goes to {N|J}(ia-1). This while() logic handles both.
    {
      switch (scur) {
      case p7T_ML: snxt = (*select_ml)(rng, gm, k, dp, sm->k[i-1], sm->n[i-1], xc-p7S_NXCELLS, &z); i--; k--;      xc -= p7S_NXCELLS; break;
      case p7T_MG: snxt = (*select_mg)(rng, gm, k, dp, sm->k[i-1], sm->n[i-1], xc-p7S_NXCELLS, &z); i--; k--;      xc -= p7S_NXCELLS; break;
      case p7T_IL: snxt = (*select_il)(rng, gm, k, dp, sm->k[i-1], sm->n[i-1],                 &z); i--;           xc -= p7S_NXCELLS; break;
      case p7T_IG: snxt = (*select_ig)(rng, gm, k, dp, sm->k[i-1], sm->n[i-1],                 &z); i--;           xc -= p7S_NXCELLS; break;
      case p7T_DL: snxt = (*select_dl)(rng, gm, k, dp + (z-1)*p7S_NSCELLS);                              k--; z--;                    break; 
      case p7T_DG: snxt = (*select_dg)(rng, gm, k, dp + (z-1)*p7S_NSCELLS);                              k--; z--;                    break;
      case p7T_N:  snxt = p7T_N;                                                                                                      break;
      case p7T_J:  snxt = (sm->n[i] ? (*select_j)(rng, gm, xc) : p7T_J);                                                              break; 
      case p7T_C:  snxt = (sm->n[i] ? (*select_c)(rng, gm, xc) : p7T_C);                                                              break; // connect to E(i), C(i-1). E(i) valid if n[i]>0; and if E(i) valid, C(i-1) must be stored too. if E(i) invalid, connect to C, and it doesn't matter if xc is valid or not
      case p7T_E:  snxt = (*select_e)(rng, wrk, gm, dp, sm->k[i], sm->n[i], &z);                         k=sm->k[i][z];   D++;        break;
      case p7T_B:  snxt = (*select_b)(rng, gm, xc);                                                                                   break; // {NJ}(i) -> B(i). If we reached B(i), xc[i] valid, so NJ must also be valid.
      case p7T_L:  snxt = p7T_B;                                                                                                      break;
      case p7T_G:  snxt = p7T_B;                                                                                                      break;
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

  if (g == 0)  // only for a complete trace, we add the S state. Segment trace ends with N or J.
    if ((status = p7_trace_Append(tr, p7T_S, 0, 0)) != eslOK) return status;

  tr->M = sm->M;
  tr->L = sm->L;
  if (opt_D) *opt_D = D;
  return p7_trace_Reverse(tr);  // trace_Reverse doesn't care if it's a complete trace or a segment trace.
}
/*--------------- end, traceback engine -------------------------*/


/*****************************************************************
 * 4. Exposed API, wrappers around the trace engine.
 *****************************************************************/

/* Function:  p7_sparse_trace_Viterbi()
 * Synopsis:  Traceback of a sparse Viterbi DP matrix.
 *
 * Purpose:   Caller has filled sparse Viterbi matrix <sx> with
 *            <p7_SparseViterbi()>, in a comparison involving profile
 *            <gm>. Perform a Viterbi traceback, storing the result in
 *            <tr>, which the caller provides. 
 *
 * Args:      gm - profile
 *            sx - filled sparse Viterbi DP matrix
 *            tr - trace structure to store the result in
 *
 * Returns:   <eslOK> on success, and <tr> contains the optimal trace.
 *
 * Throws:    <eslEMEM> on allocation failure; <tr> may need to grow to
 *            accommodate the trace.
 * 
 *            <eslEINVAL> on various coding, contract check failures. 
 */
int
p7_sparse_trace_Viterbi(const P7_PROFILE *gm, const P7_SPARSEMX *sx, P7_TRACE *tr)
{
  return (sparse_traceback_engine(NULL, NULL, gm, sx, tr, 0, NULL, NULL, NULL));
}

/* Function:  p7_sparse_trace_Stochastic()
 * Synopsis:  Stochastic traceback of a sparse Forward DP matrix.
 *
 * Purpose:   Caller has filled sparse Forward matrix <sx> with
 *            <p7_SparseForward()>, in a comparison involving profile
 *            <gm>. Using random number generator <rng>, perform a
 *            stochastic traceback, storing the sampled trace in <tr>,
 *            which the caller provides.
 *            
 *            The calculation requires a temporary scratch space of up
 *            to <M*p7S_NSCELLS> floats. In order to minimize
 *            alloc/free cycles, caller can provide an existing
 *            scratch space of any size, to be reallocated if needed,
 *            using Easel's "bypass" idiom. That is: if you pass
 *            <NULL> for <wrk_byp>, then workspace is allocated and
 *            free'd internally; this is a little slower, but saves
 *            you having to declare (and free) a variable. If you
 *            already have an allocated space <wrk>, you can pass it
 *            as <&wrk>. If <wrk==NULL> and you pass <&wrk>, then it
 *            will be allocated here; this provides a clean
 *            initialization mechanism for loops (i.e. init <wrk=NULL>
 *            and pass <&wrk> to the routine inside a loop, and it
 *            will be allocated by the first call, and reused by
 *            subsequent ones.).
 *
 * Args:      rng     - random number generator
 *            wrk_byp - ptr to workspace to reuse, or <NULL> to alloc internally instead.
 *                      <M*p7S_NSCELLS> floats are needed in this workspace; if <*wrk_byp>
 *                      is smaller, it will be reallocated; if <*wrk_byp> is <NULL> it
 *                      will be allocated.
 *            gm      - profile
 *            sxf     - filled sparse Forward DP matrix
 *            tr      - trace structure to store the result in
 *
 * Returns:   <eslOK> on success.
 *            <tr> contains the sampled trace, and may have been internally reallocated.
 *            If <wrk_byp> was non-<NULL>, <*wrk_byp> workspace may have been reallocated.
 *
 * Throws:    <eslEMEM> on allocation failure; <tr> may need to grow to
 *            accommodate the trace.
 * 
 *            <eslEINVAL> on various coding, contract check failures. 
 */
int
p7_sparse_trace_Stochastic(ESL_RANDOMNESS *rng, float **wrk_byp, const P7_PROFILE *gm, const P7_SPARSEMX *sxf, P7_TRACE *tr)
{
  float *wrk = NULL;
  int    status;

  ESL_DASSERT1( (sxf->type == p7S_FORWARD)  );

  if      (esl_byp_IsInternal(wrk_byp)) { ESL_ALLOC  (wrk,      sxf->sm->M * p7S_NSCELLS * sizeof(float));                 }
  else if (esl_byp_IsReturned(wrk_byp)) { ESL_ALLOC  (*wrk_byp, sxf->sm->M * p7S_NSCELLS * sizeof(float)); wrk = *wrk_byp; }
  else if (esl_byp_IsProvided(wrk_byp)) { ESL_REALLOC(*wrk_byp, sxf->sm->M * p7S_NSCELLS * sizeof(float)); wrk = *wrk_byp; }

  status = sparse_traceback_engine(rng, wrk, gm, sxf, tr, 0, NULL, NULL, NULL);

  if  (esl_byp_IsInternal(wrk_byp)) { free(wrk); }
  return status;

 ERROR:
  return status;
}



/* Function:  p7_sparse_trace_SegmentStochastic()
 * Synopsis:  Stochastic traceback of a single sparse segment. Used in MPAS.
 *
 * Purpose:   Caller wants a stochastic traceback for segment <g> of a
 *            sparse Forward matrix <sxf>. Caller provides random number
 *            generator <rng>, and the profile <gm> that was used in
 *            calculating <sxf>. 
 *            
 *            It's possible (indeed likely) that the segment contains
 *            no domains. This case gets handled specially, and the
 *            caller must provide the segment forward score <sF> for
 *            this handling (see footnote [6] in sparse_anchors.c for
 *            full explanation). Rather than wasting time constructing
 *            the empty path in <tr> explicitly, we simply return an
 *            empty <tr> (i.e. one with <tr->N = 0>). Caller must
 *            check for this case.
 *            
 *            Stochastic tracebacks require some auxiliary working
 *            memory, of at least <M*p7S_NSCELLS> floats. Caller may
 *            provide an existing scratch space via <wrk_byp>,
 *            following Easel's bypass idiom. You can pass <NULL>, and
 *            the necessary work space will be allocated and
 *            deallocated internally. See
 *            <p7_sparse_trace_Stochastic()> documentation for more
 *            information.
 *            
 *            If caller happens to know where row <ib> (i.e. the start
 *            point of the traceback) starts in the sparse Forward
 *            matrix storage in <sxf>, it should provide a ptr to the
 *            main matrix row <ib> via <dpc>, and to the specials via
 *            <xc>. If these are unknown, pass <NULL> for both, and
 *            they will be determined (at some cost of time).
 *            
 *            The resulting trace, stored in <tr>, is a hacked
 *            "segment trace". It starts with N|J and ends with J, and
 *            accounts for residues ia..ib. See footnote [4] in
 *            sparse_anchors.c. Although footnote [4] says that
 *            segment traces allow empty paths (all JJ...JJ), because
 *            the caller is responsible for predetermining whether the
 *            empty path is sampled, this function only returns
 *            segment traces with D>0.
 *
 * Returns:   <eslOK> on success. <tr> contains a segment trace
 *            sampled from the deltaS + p0 ensemble. If the empty 
 *            path p0 was sampled, <tr> is empty, with <tr->N=0>.
 *
 * Throws:    <eslEMEM> on allocation failure (in <wrk_byp>, or 
 *            in p7_trace_Append(), for example).
 * 
 *            <eslNOHALT> if rejection sampling in the deltaS ensemble 
 *            appears to be caught in an infinite loop. 
 */
int
p7_sparse_trace_SegmentStochastic(ESL_RANDOMNESS *rng, float **wrk_byp, const P7_PROFILE *gm, const P7_SPARSEMX *sxf, float sF, int g, float *dpc, float *xc, P7_TRACE *tr)
{
  float *wrk      = NULL;
  float  p0       = gm->xsc[p7P_J][p7P_J] * (sxf->sm->seg[g].ib - sxf->sm->seg[g].ia + 1);
  float  pempty   = exp(p0 - sF);
  int    failsafe = 0;
  int    D;
  int    status;

  /* Contract checks. */
  ESL_DASSERT1(( gm->xsc[p7P_N][p7P_N] == gm->xsc[p7P_J][p7P_J] ));   // Use of segmented MPAS requires t_NN = t_JJ = tCC. See footnote [1] in sparse_anchors.c.
  ESL_DASSERT1(( gm->xsc[p7P_N][p7P_N] == gm->xsc[p7P_C][p7P_C] ));
  ESL_DASSERT1(( pempty <= 0.5 ));                                    // See footnote [6] in sparse_anchors.c for detailed explanation.
  ESL_DASSERT1(( tr->N  == 0 ));                                      // Caller has done a _Reuse(), which we depend on.
  ESL_DASSERT1(( sxf->type == p7S_FORWARD ));
  
  /* Idiomatic handling of a "bypass" variable. wrk_byp=NULL: internal. *wrk_byp=NULL: allocate one for me and pass it back. Else: use the wrk buffer I'm giving you. */
  if      (esl_byp_IsInternal(wrk_byp)) { ESL_ALLOC  (wrk,      sxf->sm->M * p7S_NSCELLS * sizeof(float));                 }
  else if (esl_byp_IsReturned(wrk_byp)) { ESL_ALLOC  (*wrk_byp, sxf->sm->M * p7S_NSCELLS * sizeof(float)); wrk = *wrk_byp; }
  else if (esl_byp_IsProvided(wrk_byp)) { ESL_REALLOC(*wrk_byp, sxf->sm->M * p7S_NSCELLS * sizeof(float)); wrk = *wrk_byp; }

  if (esl_random(rng) >= pempty)           // else, we sample the p0 empty path, leave tr->N=0, and return success.
    do {
      status = sparse_traceback_engine(rng, wrk, gm, sxf, tr, g, dpc, xc, &D);
      if (failsafe++ >= 100000) ESL_XEXCEPTION(eslENOHALT, "so much for footnote [6]");
    } while (status == eslOK && D == 0);   // reject empty JJ...JJ traces. 
  
  if  (esl_byp_IsInternal(wrk_byp)) { free(wrk); }
  return status;

 ERROR:
  if  (esl_byp_IsInternal(wrk_byp)) { free(wrk); }
  return status;
}



/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7SPARSE_TRACE_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",              0 },
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "include all cells in sparse mx",                    0 },
  { "-n",        eslARG_INT,   "1000", NULL, NULL,   NULL,  NULL, NULL, "number of stochastic traces to sample",             0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                     0 },
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump posterior decoding matrix for examination",    0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Forward DP matrix for examination",            0 },
  { "-M",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump sparse mask for examination",                  0 },
  { "-S",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump stochastic traces for examination",            0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of using stochastic tracebacks to approximate posterior decoding";

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
  P7_SPARSEMX    *sxf     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxb     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxd     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxd2    = NULL;
  P7_TRACE       *tr      = p7_trace_Create();
  float          *wrk     = NULL;
  int             N       = esl_opt_GetInteger(go, "-n");
  int             idx;
  float           fsc;
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
    p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSEMASK_THRESH_DEFAULT);
  }
  
  /* Sparse DP calculations */
  p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf,  &fsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxb,  NULL);
  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxb, sxd);

  if (esl_opt_GetBoolean(go, "-M")) p7_sparsemask_Dump(stdout, sm);
  if (esl_opt_GetBoolean(go, "-F")) p7_sparsemx_Dump(stdout, sxf);

  /* Collect N stochastic traces, count them into <sxd> */
  sxd2 = p7_sparsemx_Create(sm);
  p7_sparsemx_Zero(sxd2);
  for (idx = 0; idx < N; idx++)
    {
      p7_sparse_trace_Stochastic(rng, &wrk, gm, sxf, tr);
      if (esl_opt_GetBoolean(go, "-S")) p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);

      p7_sparsemx_CountTrace(tr, sxd2);

      p7_trace_Reuse(tr);
    }

  /* Renormalize <sxd2> */
  esl_vec_FScale(sxd2->dp,   sxd2->sm->ncells*p7S_NSCELLS,            1./(float)N);
  esl_vec_FScale(sxd2->xmx, (sxd2->sm->nrow+sxd2->sm->S)*p7S_NXCELLS, 1./(float)N);
  
  if (esl_opt_GetBoolean(go, "-D")) p7_sparsemx_Dump(stdout, sxd2);

  p7_sparsemx_CompareDecoding(sxd, sxd2, 0.01);


  /* Cleanup */
  if (wrk) free(wrk);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_sparsemx_Destroy(sxf);
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
#endif /*p7SPARSE_TRACE_EXAMPLE*/
/*------------------ end, example driver ------------------------*/




/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/ 
