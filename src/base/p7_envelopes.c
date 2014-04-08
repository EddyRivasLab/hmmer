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
