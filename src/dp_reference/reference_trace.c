/* Traceback routines for Viterbi, MEG alignment, and stochastic
 * sampling, for the reference implementation. All tracebacks use the
 * same machinery (the reference_traceback_engine()), using different
 * select*() functions.
 * 
 * Contents:
 *   1. Choice selection for Viterbi traces
 *   2. Choice selection for stochastic traces
 *   3. Traceback engine
 *   4. Exposed API, wrappers around the engine
 *   5. Copyright and license information
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_trace.h"



/*****************************************************************
 * 1. Choice selection functions for Viterbi traces
 *****************************************************************/

static inline int
v_select_ml(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  int   state[4] = { p7T_ML, p7T_IL, p7T_DL, p7T_L };
  float path[4];

  path[0] = P7R_MX(rmx, i-1, k-1, p7R_ML) + P7P_TSC(gm, k-1, p7P_MM);
  path[1] = P7R_MX(rmx, i-1, k-1, p7R_IL) + P7P_TSC(gm, k-1, p7P_IM);
  path[2] = P7R_MX(rmx, i-1, k-1, p7R_DL) + P7P_TSC(gm, k-1, p7P_DM);
  path[3] = P7R_XMX(rmx, i-1, p7R_L)     + P7P_TSC(gm, k-1, p7P_LM);
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
v_select_mg(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  int   state[4] = { p7T_MG, p7T_IG, p7T_DG, p7T_G };
  float path[4];

  path[0] = P7R_MX(rmx, i-1, k-1, p7R_MG) + P7P_TSC(gm, k-1, p7P_MM);
  path[1] = P7R_MX(rmx, i-1, k-1, p7R_IG) + P7P_TSC(gm, k-1, p7P_IM);
  path[2] = P7R_MX(rmx, i-1, k-1, p7R_DG) + P7P_TSC(gm, k-1, p7P_DM);
  path[3] = P7R_XMX(rmx, i-1, p7R_G)      + P7P_TSC(gm, k-1, p7P_GM);
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
v_select_il(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  float path[2];

  path[0] = P7R_MX(rmx, i-1, k, p7R_ML) + P7P_TSC(gm, k, p7P_MI);
  path[1] = P7R_MX(rmx, i-1, k, p7R_IL) + P7P_TSC(gm, k, p7P_II);
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_IL);
}

static inline int
v_select_ig(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  float path[2];

  path[0] = P7R_MX(rmx, i-1, k, p7R_MG) + P7P_TSC(gm, k, p7P_MI);
  path[1] = P7R_MX(rmx, i-1, k, p7R_IG) + P7P_TSC(gm, k, p7P_II);
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_IG);
}


static inline int
v_select_dl(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  float path[2];

  path[0] = P7R_MX(rmx, i, k-1, p7R_ML) + P7P_TSC(gm, k-1, p7P_MD);
  path[1] = P7R_MX(rmx, i, k-1, p7R_DL) + P7P_TSC(gm, k-1, p7P_DD);
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_DL);
}

static inline int
v_select_dg(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  float path[2];

  path[0] = P7R_MX(rmx, i, k-1, p7R_MG) + P7P_TSC(gm, k-1, p7P_MD);
  path[1] = P7R_MX(rmx, i, k-1, p7R_DG) + P7P_TSC(gm, k-1, p7P_DD);
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_DG);
}

static inline int
v_select_e(ESL_RANDOMNESS *r, float *wrk, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int *ret_k)
{
  ESL_UNUSED(r);
  ESL_UNUSED(wrk);
  float max  = -eslINFINITY;
  int   smax = -1;
  int   kmax = -1;
  int   k;

  for (k = 1; k <= gm->M; k++)
    if (P7R_MX(rmx, i, k, p7R_ML) > max) { max = P7R_MX(rmx, i, k, p7R_ML); smax = p7T_ML; kmax = k; }

  if (P7R_MX(rmx, i, gm->M, p7R_MG) > max) { max = P7R_MX(rmx, i, gm->M, p7R_MG); smax = p7T_MG; kmax = gm->M; }
  if (P7R_MX(rmx, i, gm->M, p7R_DG) > max) {                                      smax = p7T_DG; kmax = gm->M; }

  *ret_k = kmax;
  return smax;
}

static inline int
v_select_j(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = P7R_XMX(rmx, i-1, p7R_J) + gm->xsc[p7P_J][p7P_LOOP];
  path[1] = P7R_XMX(rmx, i,   p7R_E) + gm->xsc[p7P_E][p7P_LOOP];
  return ( (path[0] > path[1]) ? p7T_J : p7T_E);
}

static inline int
v_select_b( ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = P7R_XMX(rmx, i, p7R_N) + gm->xsc[p7P_N][p7P_MOVE];
  path[1] = P7R_XMX(rmx, i, p7R_J) + gm->xsc[p7P_J][p7P_MOVE];
  return ( (path[0] > path[1]) ? p7T_N : p7T_J);
}

static inline int
v_select_c(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = P7R_XMX(rmx, i-1, p7R_C) + gm->xsc[p7P_C][p7P_LOOP];
  path[1] = P7R_XMX(rmx, i,   p7R_E) + gm->xsc[p7P_E][p7P_MOVE];
  return ( (path[0] > path[1]) ? p7T_C : p7T_E);
}
/*----------- end, viterbi selection functions ------------------*/



/*****************************************************************
 * 2. Selection functions for stochastic trace of Fwd matrix
 *****************************************************************/


static inline int
sto_select_ml(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  int   state[4] = { p7T_ML, p7T_IL, p7T_DL, p7T_L };
  float path[4];

  path[0] = P7R_MX(rmx, i-1, k-1, p7R_ML) + P7P_TSC(gm, k-1, p7P_MM);
  path[1] = P7R_MX(rmx, i-1, k-1, p7R_IL) + P7P_TSC(gm, k-1, p7P_IM);
  path[2] = P7R_MX(rmx, i-1, k-1, p7R_DL) + P7P_TSC(gm, k-1, p7P_DM);
  path[3] = P7R_XMX(rmx, i-1, p7R_L)     + P7P_TSC(gm, k-1, p7P_LM);
  esl_vec_FLogNorm(path,4);
  return state[esl_rnd_FChoose(r, path, 4)];
}

static inline int
sto_select_mg(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  int   state[4] = { p7T_MG, p7T_IG, p7T_DG, p7T_G };
  float path[4];

  path[0] = P7R_MX(rmx, i-1, k-1, p7R_MG) + P7P_TSC(gm, k-1, p7P_MM);
  path[1] = P7R_MX(rmx, i-1, k-1, p7R_IG) + P7P_TSC(gm, k-1, p7P_IM);
  path[2] = P7R_MX(rmx, i-1, k-1, p7R_DG) + P7P_TSC(gm, k-1, p7P_DM);
  path[3] = P7R_XMX(rmx, i-1, p7R_G)      + P7P_TSC(gm, k-1, p7P_GM);
  esl_vec_FLogNorm(path,4);
  return state[esl_rnd_FChoose(r, path, 4)];
}

static inline int
sto_select_il(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  float path[2];

  path[0] = P7R_MX(rmx, i-1, k, p7R_ML) + P7P_TSC(gm, k, p7P_MI);
  path[1] = P7R_MX(rmx, i-1, k, p7R_IL) + P7P_TSC(gm, k, p7P_II);
  esl_vec_FLogNorm(path, 2);
  return ( (esl_rnd_FChoose(r,path,2) == 0) ? p7T_ML : p7T_IL);
}

static inline int
sto_select_ig(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  float path[2];

  path[0] = P7R_MX(rmx, i-1, k, p7R_MG) + P7P_TSC(gm, k, p7P_MI);
  path[1] = P7R_MX(rmx, i-1, k, p7R_IG) + P7P_TSC(gm, k, p7P_II);
  esl_vec_FLogNorm(path, 2);
  return ( (esl_rnd_FChoose(r,path,2) == 0) ? p7T_MG : p7T_IG);
}


static inline int
sto_select_dl(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  float path[2];

  path[0] = P7R_MX(rmx, i, k-1, p7R_ML) + P7P_TSC(gm, k-1, p7P_MD);
  path[1] = P7R_MX(rmx, i, k-1, p7R_DL) + P7P_TSC(gm, k-1, p7P_DD);
  esl_vec_FLogNorm(path, 2);
  return ( (esl_rnd_FChoose(r,path,2) == 0) ? p7T_ML : p7T_DL);
}

static inline int
sto_select_dg(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  float path[2];

  path[0] = P7R_MX(rmx, i, k-1, p7R_MG) + P7P_TSC(gm, k-1, p7P_MD);
  path[1] = P7R_MX(rmx, i, k-1, p7R_DG) + P7P_TSC(gm, k-1, p7P_DD);
  esl_vec_FLogNorm(path, 2);
  return ( (esl_rnd_FChoose(r,path,2) == 0) ? p7T_MG : p7T_DG);
}

static inline int
sto_select_e(ESL_RANDOMNESS *r, float *wrk, const P7_PROFILE *gm, const P7_REFMX *rmx, int i, int *ret_k)
{
  ESL_UNUSED(r);
  ESL_UNUSED(wrk);
  int   k;

  wrk[0]                                     = P7R_MX(rmx, i, gm->M, p7R_DG); /* 0   : DGm->E exit           */
  for (k = 1; k <= gm->M; k++)  wrk[k]       = P7R_MX(rmx, i, k,     p7R_DL); /* 1..M: DLk->E exits for k=1..M */
  for (k = 1; k <= gm->M; k++)  wrk[k+gm->M] = P7R_MX(rmx, i, k,     p7R_ML); /* M+1..2M : MLk->E exits for k=1..M */
  wrk[2*gm->M+1]                             = P7R_MX(rmx, i, gm->M, p7R_MG); /* 2M+1    : MGk->E exit */

  esl_vec_FLogNorm(wrk, 2*gm->M+2);
  k = esl_rnd_FChoose(r, wrk, 2*gm->M+2);
  if      (k == 2*gm->M+1) { *ret_k = gm->M;   return p7T_MG; }
  else if (k == 0)         { *ret_k = gm->M;   return p7T_DG; }
  else if (k <= gm->M)     { *ret_k = k;       return p7T_DL; }
  else                     { *ret_k = k-gm->M; return p7T_ML; }
}

static inline int
sto_select_j(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = P7R_XMX(rmx, i-1, p7R_J) + gm->xsc[p7P_J][p7P_LOOP];
  path[1] = P7R_XMX(rmx, i,   p7R_E) + gm->xsc[p7P_E][p7P_LOOP];
  esl_vec_FLogNorm(path, 2);
  return (  (esl_rnd_FChoose(r,path,2) == 0) ? p7T_J : p7T_E);
}

static inline int
sto_select_b( ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = P7R_XMX(rmx, i, p7R_N) + gm->xsc[p7P_N][p7P_MOVE];
  path[1] = P7R_XMX(rmx, i, p7R_J) + gm->xsc[p7P_J][p7P_MOVE];
  esl_vec_FLogNorm(path, 2);
  return (  (esl_rnd_FChoose(r,path,2) == 0) ? p7T_N : p7T_J);
}

static inline int
sto_select_c(ESL_RANDOMNESS *r, const P7_PROFILE *gm, const P7_REFMX *rmx, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = P7R_XMX(rmx, i-1, p7R_C) + gm->xsc[p7P_C][p7P_LOOP];
  path[1] = P7R_XMX(rmx, i,   p7R_E) + gm->xsc[p7P_E][p7P_MOVE];
  esl_vec_FLogNorm(path, 2);
  return (  (esl_rnd_FChoose(r,path,2) == 0) ? p7T_C : p7T_E);
}
/*--------- end, stochastic selection functions -----------------*/




/*****************************************************************
 * 3. Traceback engine
 *****************************************************************/

static int 
reference_trace_engine(ESL_RANDOMNESS *rng, float *wrk, const P7_PROFILE *gm, const P7_REFMX *rmx, P7_TRACE *tr)
{
  int   i    = rmx->L;
  int   k    = 0;
  int   sprv = p7T_C;
  int   scur;
  int   status;

  int (*select_ml)(ESL_RANDOMNESS *,          const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_mg)(ESL_RANDOMNESS *,          const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_il)(ESL_RANDOMNESS *,          const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_ig)(ESL_RANDOMNESS *,          const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_dl)(ESL_RANDOMNESS *,          const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_dg)(ESL_RANDOMNESS *,          const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_e) (ESL_RANDOMNESS *, float *, const P7_PROFILE *, const P7_REFMX *, int, int *);
  int (*select_j) (ESL_RANDOMNESS *,          const P7_PROFILE *, const P7_REFMX *, int);
  int (*select_b) (ESL_RANDOMNESS *,          const P7_PROFILE *, const P7_REFMX *, int);
  int (*select_c) (ESL_RANDOMNESS *,          const P7_PROFILE *, const P7_REFMX *, int);
  
  /* if the target sequence is impossible: leave trace alone, with tr->N=0, our convention for an impossible trace */
  if (P7R_XMX(rmx, rmx->L, p7R_C) == -eslINFINITY) { tr->M = rmx->M; tr->L = rmx->L; return eslOK; }

  if (rmx->type == p7R_VITERBI)	/* configure for Viterbi optimal traceback */
    {
      select_ml = &v_select_ml;      select_mg = &v_select_mg;
      select_il = &v_select_il;      select_ig = &v_select_ig;
      select_dl = &v_select_dl;      select_dg = &v_select_dg;
      select_j  = &v_select_j;       select_c  = &v_select_c;
      select_e  = &v_select_e;       select_b  = &v_select_b;
    }
  else if (rmx->type == p7R_FORWARD) /* configure for stochastic traceback */
    {
      select_ml = &sto_select_ml;      select_mg = &sto_select_mg;
      select_il = &sto_select_il;      select_ig = &sto_select_ig;
      select_dl = &sto_select_dl;      select_dg = &sto_select_dg;
      select_j  = &sto_select_j;       select_c  = &sto_select_c;
      select_e  = &sto_select_e;       select_b  = &sto_select_b;
    }

  if ((status = p7_trace_Append(tr, p7T_T, k, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, k, i)) != eslOK) return status;
  while (sprv != p7T_S)
    {
      switch (sprv) {
      case p7T_ML: scur = (*select_ml)(rng,      gm, rmx, i, k); k--; i--; break;
      case p7T_MG: scur = (*select_mg)(rng,      gm, rmx, i, k); k--; i--; break;
      case p7T_IL: scur = (*select_il)(rng,      gm, rmx, i, k);      i--; break;
      case p7T_IG: scur = (*select_ig)(rng,      gm, rmx, i, k);      i--; break;
      case p7T_DL: scur = (*select_dl)(rng,      gm, rmx, i, k); k--;      break;
      case p7T_DG: scur = (*select_dg)(rng,      gm, rmx, i, k); k--;      break;
      case p7T_E:  scur = (*select_e) (rng, wrk, gm, rmx, i, &k);          break;
      case p7T_N:  scur = (i==0 ? p7T_S : p7T_N);                          break;
      case p7T_J:  scur = (*select_j) (rng, gm, rmx, i);                   break;
      case p7T_B:  scur = (*select_b) (rng, gm, rmx, i);                   break;
      case p7T_L:  scur = p7T_B;                                           break;
      case p7T_G:  scur = p7T_B;                                           break;
      case p7T_C:  scur = (*select_c) (rng, gm, rmx, i);                   break;
      default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
      }

      /* A glocal B->G->Mk wing-retraction entry: unfold it */
      if (scur == p7T_G) {
	while (k) {
	  if ( (status = p7_trace_Append(tr, p7T_DG, k, i)) != eslOK) return status;
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
/*---------------- end, traceback engine ------------------------*/

/*****************************************************************
 * 4. Exposed API, wrappers around the trace engine 
 *****************************************************************/ 

int
p7_reference_trace_Viterbi(const P7_PROFILE *gm, const P7_REFMX *rmx, P7_TRACE *tr)
{
  return (reference_trace_engine(NULL, NULL, gm, rmx, tr));
}

int
p7_reference_trace_Stochastic(ESL_RANDOMNESS *rng, float **wrk_byp, const P7_PROFILE *gm,  const P7_REFMX *rmx, P7_TRACE *tr)
{
  float *wrk = NULL;
  int    status;

  if      (esl_byp_IsInternal(wrk_byp)) { ESL_ALLOC  (wrk,      (2*gm->M+2) * sizeof(float));                 }
  else if (esl_byp_IsReturned(wrk_byp)) { ESL_ALLOC  (*wrk_byp, (2*gm->M+2) * sizeof(float)); wrk = *wrk_byp; }
  else if (esl_byp_IsProvided(wrk_byp)) { ESL_REALLOC(*wrk_byp, (2*gm->M+2) * sizeof(float)); wrk = *wrk_byp; }

  status = reference_trace_engine(rng, wrk, gm, rmx, tr);

  if  (esl_byp_IsInternal(wrk_byp)) { free(wrk); }
  return status;
 ERROR:
  return status;
}

/*----------------- end, API wrappers ---------------------------*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/

