/* P7_ENVELOPES
 */
#include "p7_config.h"

#include <stdlib.h>

#include "easel.h"

#include "base/p7_envelopes.h"


/*****************************************************************
 * 1. The P7_ENVELOPES object.
 *****************************************************************/

P7_ENVELOPES *
p7_envelopes_Create(int32_t nalloc, int32_t nredline)
{
  P7_ENVELOPES *envs = NULL;
  int           status;

  ESL_ALLOC(envs, sizeof(P7_ENVELOPES));
  envs->arr      = NULL;
  envs->n        = 0;
  envs->L        = 0;
  envs->M        = 0;
  envs->nalloc   = (nalloc   > 0 ? nalloc   : 8);
  envs->nredline = (nredline > 0 ? nredline : 64);
  
  ESL_ALLOC(envs->arr, sizeof(P7_ENVELOPE) * envs->nalloc);
  return envs;

 ERROR:

  p7_envelopes_Destroy(envs);
  return NULL;
}

int
p7_envelopes_GrowTo(P7_ENVELOPES *envs, int32_t nalloc)
{
  int status;

  if (envs->nalloc >= nalloc) return eslOK;
  
  ESL_REALLOC(envs->arr, sizeof(P7_ENVELOPE) * nalloc);
  envs->nalloc = nalloc;
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
  
  envs->n = 0;
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

  fprintf(ofp, "#%3s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %6s %3s %3s\n",
	  "dom", "ia",  "ib", "i0",  "k0", "alia", "alib", "ka",  "kb", 
	  "oea",  "oeb", "env_sc", "app", "glo");
  fprintf(ofp, "#%3s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %6s %3s %3s\n",
	  "---", "-----",  "-----", "-----",  "-----", "-----", "-----", "-----",  "-----", 
	  "-----",  "-----", "------", "---", "---");
  for (e = 0; e < env->n; e++)
    fprintf(ofp, "%-4d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %6.2f %3s %3s\n",
	    e+1,
	    env->arr[e].ia,   env->arr[e].ib,
	    env->arr[e].i0,   env->arr[e].k0,
	    env->arr[e].alia, env->arr[e].alib,
	    env->arr[e].ka,   env->arr[e].kb,
	    env->arr[e].oea,  env->arr[e].oeb,
	    env->arr[e].env_sc,
	    (env->arr[e].flags & p7E_ENVSC_APPROX ? "YES" : "n"),
	    (env->arr[e].flags & p7E_IS_GLOCAL    ? "YES" : "n"));

  fprintf(ofp, "# Total domains = %d\n", env->n);
  fprintf(ofp, "# M,L           = %d,%d\n", env->M, env->L);
  fprintf(ofp, "# nalloc        = %d\n", env->nalloc);
  fprintf(ofp, "# nredline      = %d\n", env->nredline);
  return eslOK;
}
