/* P7_ENVELOPES manages an array of "envelopes", information about
 * each detected domain in a target sequence.
 * 
 * Contents:
 *    1. P7_ENVELOPE structure (one envelope)
 *    2. P7_ENVELOPES object (array of envelopes)
 *    x. Copyright and license information
 */
#include "p7_config.h"

#include <stdlib.h>

#include "easel.h"

#include "base/p7_envelopes.h"


/*****************************************************************
 * 1. The P7_ENVELOPE data structure
 *****************************************************************/

int
p7_envelope_SetSentinels(P7_ENVELOPE *env, int D, int L, int M)
{
  /* Anchor i0,k0 sentinels are set exactly as p7_anchor_SetSentinels() does; see explanation over there */
  env[0].i0   = 0;
  env[D+1].i0 = L+1;
  env[0].k0   = M+1;
  env[D+1].k0 = 0;

  /* I don't think other sentinel values matter, but here's a
   * guess what they would be if they mattered; make i 
   * coords = i0, and k coords = k0.
   */
  env[0].oea  = env[0].oeb  = 0;
  env[0].ia   = env[0].ib   = 0;

  env[D+1].oea  = env[D+1].oeb  = L+1;
  env[D+1].ia   = env[D+1].ib   = L+1;

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
  int e;

  fprintf(ofp, "#%3s %5s %5s %5s %5s %5s %5s %5s %5s %6s %3s %3s\n",
	  "dom", "ia",  "ib", "i0",  "k0", "ka",  "kb", 
	  "oea",  "oeb", "env_sc", "app", "glo");
  fprintf(ofp, "#%3s %5s %5s %5s %5s %5s %5s %5s %5s %6s %3s %3s\n",
	  "---", "-----",  "-----", "-----",  "-----", "-----",  "-----", 
	  "-----",  "-----", "------", "---", "---");
  for (e = 1; e <= env->D; e++)
    fprintf(ofp, "%-4d %5d %5d %5d %5d %5d %5d %5d %5d %6.2f %3s %3s\n",
	    e,
	    env->arr[e].ia,   env->arr[e].ib,
	    env->arr[e].i0,   env->arr[e].k0,
	    env->arr[e].ka,   env->arr[e].kb,
	    env->arr[e].oea,  env->arr[e].oeb,
	    env->arr[e].env_sc,
	    (env->arr[e].flags & p7E_ENVSC_APPROX ? "YES" : "n"),
	    (env->arr[e].flags & p7E_IS_GLOCAL    ? "YES" : "n"));

  fprintf(ofp, "# Total domains = %d\n",    env->D);
  fprintf(ofp, "# M,L           = %d,%d\n", env->M, env->L);
  fprintf(ofp, "# nalloc        = %d\n",    env->nalloc);
  fprintf(ofp, "# nredline      = %d\n",    env->nredline);
  return eslOK;
}
