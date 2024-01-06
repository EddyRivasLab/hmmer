/* P7_ENVELOPES manages an array of "envelopes", information about
 * each detected domain in a target sequence.
 * 
 * Contents:
 *    1. P7_ENVELOPE structure (one envelope)
 *    2. P7_ENVELOPES object (array of envelopes)
 */
#include <p7_config.h>

#include <stdlib.h>

#include "easel.h"

#include "base/p7_envelopes.h"


/*****************************************************************
 * 1. The P7_ENVELOPE data structure
 *****************************************************************/

int
p7_envelope_SetSentinels(P7_ENVELOPE *env, int D, int L, int M)
{
  /* Anchor i0,k0 sentinels are set exactly as
   * p7_anchor_SetSentinels() does; see explanation over there.
   *
   * Likewise, in envelope determination we want maximum envelope
   * bounds to be i0(d-1)+1..i0(d)..i0(d+1)-1; i.e. 
   *    i0(d-1) + 1 <= oa(d) <= ia(d) <= i0(d), and
   *    i0(d)       <= ib(d) <= ob(d) <= i0(d+1)-1
   * so for that too, we need i0(0) = 0, i0(D+1) = L+1.
   */
  env[0].i0   = 0;
  env[D+1].i0 = L+1;
  env[0].k0   = M+1;
  env[D+1].k0 = 0;

  /* I don't think other sentinel values matter, but here's a
   * guess what they would be if they did matter; make i 
   * coords = i0, and k coords = k0.
   */
  env[0].oa    = env[0].ob    = 0;
  env[D+1].oa  = env[D+1].ob  = L+1;

  env[0].ia    = env[0].ib    = 0;
  env[D+1].ia  = env[D+1].ib  = L+1;

  env[0].ka   = env[0].kb   = M+1;
  env[D+1].ka = env[D+1].kb = 0;

  env[0].env_sc = env[D+1].env_sc = 0.0;
  env[0].flags  = env[D+1].flags  = 0;
  return eslOK;
}


/*****************************************************************
 * 2. The P7_ENVELOPES object.
 *****************************************************************/

P7_ENVELOPES *
p7_envelopes_Create(void)
{
  P7_ENVELOPES *envs             = NULL;
  int           default_nalloc   = 8;    // i.e. for up to D=6, because of sentinels
  int           default_nredline = 64;   //          ...to D=62
  int           status;

  ESL_ALLOC(envs, sizeof(P7_ENVELOPES));
  envs->arr      = NULL;
  envs->D        = 0;
  envs->L        = 0;
  envs->M        = 0;
  envs->nalloc   = default_nalloc;
  envs->nredline = default_nredline;
  
  ESL_ALLOC(envs->arr, sizeof(P7_ENVELOPE) * envs->nalloc);
  return envs;

 ERROR:
  p7_envelopes_Destroy(envs);
  return NULL;
}

int
p7_envelopes_Reinit(P7_ENVELOPES *envs, int D)
{
  int status;

  if ( D+2 > envs->nalloc) {   // grow?
    ESL_REALLOC(envs->arr, sizeof(P7_ENVELOPE) * (D+2));
    envs->nalloc = D+2;
  }                           
  else if (envs->nalloc > envs->nredline && D+2 <= envs->nredline) { // shrink?
    ESL_REALLOC(envs->arr, sizeof(P7_ENVELOPE) * envs->nredline);
    envs->nalloc = envs->nredline;
  }

  envs->D = 0;
  envs->L = 0;
  envs->M = 0;
  return eslOK;

 ERROR:
  return status;
}

int
p7_envelopes_Reuse(P7_ENVELOPES *envs)
{
  int status;

  if (envs->nalloc > envs->nredline)
    {
      ESL_REALLOC(envs->arr, sizeof(P7_ENVELOPE) * envs->nredline);
      envs->nalloc = envs->nredline;
    }
  
  envs->D = 0;
  envs->L = 0;
  envs->M = 0;
  return eslOK;

 ERROR:
  return status;
}

void
p7_envelopes_Destroy(P7_ENVELOPES *envs)
{
  if (envs) 
    {
      if (envs->arr) free(envs->arr);
      free(envs);
    }
  return;
}


/*****************************************************************
 * 2. Debugging and development tools
 *****************************************************************/

int
p7_envelopes_Dump(FILE *ofp, P7_ENVELOPES *env)
{
  int   d;
  float bias;
  float Sdd;
  float tCC   = (float) (env->L) / (float) (env->L+3);
  float omega = log(1./256.);  // hardcoded - FIXME

  esl_dataheader(ofp, 5, "dom", 5, "ia", 5, "ib", 
		 5, "i0", 5, "k0", 5, "ka", 5, "kb",
		 5, "oa", 5, "ob",
		 6, "env_sc", 6, "null2", 6, "bias", 6, "Sdd",
		 3, "app", 3, "glo",
		 0);

  for (d = 1; d <= env->D; d++)
    {
      bias = log(1 + exp(env->arr[d].null2_sc + omega));
      Sdd  = eslCONST_LOG2R * 
	(env->arr[d].env_sc - ( (float) env->L * log(tCC) + log(1-tCC)) - bias);

      fprintf(ofp, "%-5d %5d %5d %5d %5d %5d %5d %5d %5d %6.2f %6.2f %6.2f %6.2f %3s %3s\n",
	      d,
	      env->arr[d].ia,  env->arr[d].ib,
	      env->arr[d].i0,  env->arr[d].k0,
	      env->arr[d].ka,  env->arr[d].kb,
	      env->arr[d].oa,  env->arr[d].ob,
	      env->arr[d].env_sc,
	      env->arr[d].null2_sc,
	      bias, 
	      Sdd,
	      (env->arr[d].flags & p7E_ENVSC_APPROX ? "YES" : "n"),
	      (env->arr[d].flags & p7E_IS_GLOCAL    ? "YES" : "n"));
    }

  fprintf(ofp, "# Total domains = %d\n",    env->D);
  fprintf(ofp, "# M,L           = %d,%d\n", env->M, env->L);
  fprintf(ofp, "# nalloc        = %d\n",    env->nalloc);
  fprintf(ofp, "# nredline      = %d\n",    env->nredline);
  return eslOK;
}
