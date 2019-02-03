/* Traceback routines for Viterbi, MEG alignment, and stochastic
 * sampling, for the reference implementation. All tracebacks use the
 * same machinery (the reference_traceback_engine()), using different
 * select*() functions.
 * 
 * Contents: 
 *    1. Selection functions for Viterbi tracebacks
 *    2. Selection functions for stochastic tracebacks
 *    3. Traceback engine
 *    4. Exposed API, wrappers around the engine
 */

#include "h4_config.h"

#include "easel.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_path.h"
#include "h4_refmx.h"


/*****************************************************************
 * 1. Selection functions for Viterbi tracebacks
 *****************************************************************/

static inline int8_t
v_select_ml(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[4] = { h4P_ML, h4P_IL, h4P_DL, h4P_L };
  float  path[4];

  path[0] = H4R_MX(rmx, i-1, k-1, h4R_ML) + hmm->tsc[k-1][h4_MM];
  path[1] = H4R_MX(rmx, i-1, k-1, h4R_IL) + hmm->tsc[k-1][h4_IM];
  path[2] = H4R_MX(rmx, i-1, k-1, h4R_DL) + hmm->tsc[k-1][h4_DM];
  path[3] = H4R_XMX(rmx, i-1, h4R_L)      + hmm->tsc[k-1][h4_LM];
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int8_t
v_select_mg(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[4] = { h4P_MG, h4P_IG, h4P_DG, h4P_G };
  float  path[4];

  path[0] = H4R_MX(rmx, i-1, k-1, h4R_MG) + hmm->tsc[k-1][h4_MM];
  path[1] = H4R_MX(rmx, i-1, k-1, h4R_IG) + hmm->tsc[k-1][h4_IM];
  path[2] = H4R_MX(rmx, i-1, k-1, h4R_DG) + hmm->tsc[k-1][h4_DM];
  path[3] = H4R_XMX(rmx, i-1, h4R_G)      + hmm->tsc[k-1][h4_GM];
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int8_t
v_select_il(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[3] = { h4P_ML, h4P_IL, h4P_DL };
  float  path[3];

  path[0] = H4R_MX(rmx, i-1, k, h4R_ML) + hmm->tsc[k][h4_MI];
  path[1] = H4R_MX(rmx, i-1, k, h4R_IL) + hmm->tsc[k][h4_II];
  path[2] = H4R_MX(rmx, i-1, k, h4R_DL) + hmm->tsc[k][h4_DI];
  return state[esl_vec_FArgMax(path, 3)];
}

static inline int8_t
v_select_ig(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[3] = { h4P_MG, h4P_IG, h4P_DG };
  float  path[3];

  path[0] = H4R_MX(rmx, i-1, k, h4R_MG) + hmm->tsc[k][h4_MI];
  path[1] = H4R_MX(rmx, i-1, k, h4R_IG) + hmm->tsc[k][h4_II];
  path[2] = H4R_MX(rmx, i-1, k, h4R_DG) + hmm->tsc[k][h4_DI];
  return state[esl_vec_FArgMax(path, 3)];
}

static inline int8_t
v_select_dl(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[3] = { h4P_ML, h4P_IL, h4P_DL };
  float  path[3];

  path[0] = H4R_MX(rmx, i, k-1, h4R_ML) + hmm->tsc[k-1][h4_MD];
  path[1] = H4R_MX(rmx, i, k-1, h4R_IL) + hmm->tsc[k-1][h4_ID];
  path[2] = H4R_MX(rmx, i, k-1, h4R_DL) + hmm->tsc[k-1][h4_DD];
  return state[esl_vec_FArgMax(path, 3)];
}

static inline int8_t
v_select_dg(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rmx, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[3] = { h4P_MG, h4P_IG, h4P_DG };
  float  path[3];

  path[0] = H4R_MX(rmx, i, k-1, h4R_MG) + hmm->tsc[k-1][h4_MD];
  path[1] = H4R_MX(rmx, i, k-1, h4R_IG) + hmm->tsc[k-1][h4_ID];
  path[2] = H4R_MX(rmx, i, k-1, h4R_DG) + hmm->tsc[k-1][h4_DD];
  return state[esl_vec_FArgMax(path, 3)];
}

static inline int8_t
v_select_e(ESL_RANDOMNESS *r, float *wrk, const H4_PROFILE *hmm, const H4_REFMX *rmx, int i, int *ret_k)
{
  ESL_UNUSED(r);
  ESL_UNUSED(wrk);
  float  max  = -eslINFINITY;
  int8_t smax = -1;
  int    kmax = -1;
  int    k;

  for (k = 1; k <= hmm->M; k++)
    if (H4R_MX(rmx, i, k, h4R_ML) > max) { max = H4R_MX(rmx, i, k, h4R_ML); smax = h4P_ML; kmax = k; }

  if (H4R_MX(rmx, i, hmm->M, h4R_MG) > max) { max = H4R_MX(rmx, i, hmm->M, h4R_MG); smax = h4P_MG; kmax = hmm->M; }
  if (H4R_MX(rmx, i, hmm->M, h4R_DG) > max) {                                       smax = h4P_DG; kmax = hmm->M; }

  *ret_k = kmax;
  return smax;
}

static inline int8_t
v_select_j(ESL_RANDOMNESS *r, const H4_MODE *mo, const H4_REFMX *rmx, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = H4R_XMX(rmx, i-1, h4R_J) + mo->xsc[h4_J][h4_LOOP];
  path[1] = H4R_XMX(rmx, i,   h4R_E) + mo->xsc[h4_E][h4_LOOP];
  return ( (path[0] > path[1]) ? h4P_J : h4P_E);
}

static inline int8_t
v_select_b( ESL_RANDOMNESS *r, const H4_MODE *mo, const H4_REFMX *rmx, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = H4R_XMX(rmx, i, h4R_N) + mo->xsc[h4_N][h4_MOVE];
  path[1] = H4R_XMX(rmx, i, h4R_J) + mo->xsc[h4_J][h4_MOVE];
  return ( (path[0] > path[1]) ? h4P_N : h4P_J);
}

static inline int8_t
v_select_c(ESL_RANDOMNESS *r, const H4_MODE *mo, const H4_REFMX *rmx, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = H4R_XMX(rmx, i-1, h4R_C) + mo->xsc[h4_C][h4_LOOP];
  path[1] = H4R_XMX(rmx, i,   h4R_E) + mo->xsc[h4_E][h4_MOVE];
  return ( (path[0] > path[1]) ? h4P_C : h4P_E);
}





/*****************************************************************
 * 2. Selection functions for stochastic tracebacks
 *****************************************************************/ 

/*--- TK TK TK  src/dp_reference/reference_trace.c ---*/




/*****************************************************************
 * 3. Traceback engine
 *****************************************************************/ 

static int 
reference_trace_engine(ESL_RANDOMNESS *rng, float *wrk, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_REFMX *rmx, H4_PATH *pi)
{
  int    i    = rmx->L;
  int    k    = 0;
  int8_t sprv = h4P_C;
  int8_t scur;
  int    status;

  int8_t (*select_ml)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_mg)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_il)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_ig)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_dl)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_dg)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_e) (ESL_RANDOMNESS *, float *, const H4_PROFILE *, const H4_REFMX *, int, int *);
  int8_t (*select_j) (ESL_RANDOMNESS *,          const H4_MODE *,    const H4_REFMX *, int);
  int8_t (*select_b) (ESL_RANDOMNESS *,          const H4_MODE *,    const H4_REFMX *, int);
  int8_t (*select_c) (ESL_RANDOMNESS *,          const H4_MODE *,    const H4_REFMX *, int);
  
  /* if the target sequence is impossible: leave trace alone, with tr->Z=0, our convention for an impossible trace */
  if (H4R_XMX(rmx, rmx->L, h4R_C) == -eslINFINITY)  return eslOK; 

  if (rmx->type == h4R_VITERBI)	/* configure for Viterbi optimal traceback */
    {
      select_ml = &v_select_ml;      select_mg = &v_select_mg;
      select_il = &v_select_il;      select_ig = &v_select_ig;
      select_dl = &v_select_dl;      select_dg = &v_select_dg;
      select_j  = &v_select_j;       select_c  = &v_select_c;
      select_e  = &v_select_e;       select_b  = &v_select_b;
    }
#if 0 // NWO TK
  else if (rmx->type == h4R_FORWARD) /* configure for stochastic traceback */
    {
      select_ml = &sto_select_ml;      select_mg = &sto_select_mg;
      select_il = &sto_select_il;      select_ig = &sto_select_ig;
      select_dl = &sto_select_dl;      select_dg = &sto_select_dg;
      select_j  = &sto_select_j;       select_c  = &sto_select_c;
      select_e  = &sto_select_e;       select_b  = &sto_select_b;
    }
#endif

  if ((status = h4_path_Append(pi, h4P_C)) != eslOK) return status;
  while (sprv != h4P_S)
    {
      switch (sprv) {
      case h4P_ML: scur = (*select_ml)(rng,      hmm, rmx, i, k); k--; i--; break;
      case h4P_MG: scur = (*select_mg)(rng,      hmm, rmx, i, k); k--; i--; break;
      case h4P_IL: scur = (*select_il)(rng,      hmm, rmx, i, k);      i--; break;
      case h4P_IG: scur = (*select_ig)(rng,      hmm, rmx, i, k);      i--; break;
      case h4P_DL: scur = (*select_dl)(rng,      hmm, rmx, i, k); k--;      break;
      case h4P_DG: scur = (*select_dg)(rng,      hmm, rmx, i, k); k--;      break;
      case h4P_E:  scur = (*select_e) (rng, wrk, hmm, rmx, i, &k);          break;
      case h4P_N:  scur = (i==0 ? h4P_S : h4P_N);                           break;
      case h4P_J:  scur = (*select_j) (rng, mo, rmx, i);                    break;
      case h4P_B:  scur = (*select_b) (rng, mo, rmx, i);                    break;
      case h4P_L:  scur = h4P_B;                                            break;
      case h4P_G:  scur = h4P_B;                                            break;
      case h4P_C:  scur = (*select_c) (rng, mo, rmx, i);                    break;
      default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
      }

      /* A glocal B->G->Mk wing-retraction entry: unfold it */
      if (scur == h4P_G) {
	while (k) {
	  if ( (status = h4_path_Append(pi, h4P_DG)) != eslOK) return status;
	  k--;
	}
      }

      if ( (status = h4_path_Append(pi, scur)) != eslOK) return status;

      if ( (scur == h4P_N || scur == h4P_J || scur == h4P_C) && scur == sprv) i--;
      sprv = scur;
    }
  
  return h4_path_Reverse(pi);
}



int
h4_reference_trace_Viterbi(const H4_PROFILE *hmm, const H4_MODE *mo, const H4_REFMX *rmx, H4_PATH *pi)
{
  return reference_trace_engine(NULL, NULL, hmm, mo, rmx, pi);
}

