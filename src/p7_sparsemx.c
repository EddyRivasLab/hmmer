#include "p7_config.h"

#include "easel.h"

P7_SPARSEMASK *
p7_sparsemask_CreateFull(int M, int L)
{
  P7_SPARSEMASK *sm            = NULL;
  int            i;

  ESL_ALLOC(sm, sizeof(P7_SPARSEMASK));
  sm->L      = L;
  sm->M      = M;

  sm->kin    = NULL;
  sm->nin    = NULL;
  sm->kmem   = NULL;

  sm->ralloc   = L+1;
  sm->kalloc   = (int64_t) L * (int64_t) M;

  ESL_ALLOC(sm->kin,  sizeof(int *) * sm->ralloc);
  ESL_ALLOC(sm->nin,  sizeof(int)   * sm->ralloc);
  ESL_ALLOC(sm->kmem, sizeof(int)   * sm->kalloc);

  sm->kin[0] = NULL;
  sm->nin[0] = 0;
  sm->ncells = 0;

  for (i = 1; i <= L; i++) 
    {
      sm->kin[i] = sm->kmem + sm->ncells;
      for (k = 1; k <= sm->M; k++) sm->kin[i][k-1] = k; 
      sm->nin[i]  = sm->M;
      sm->ncells += sm->M;
    }
  return sm;

 ERROR:
  p7_sparsemask_Destroy(sm);
  return NULL;
}

void
p7_sparsemask_Destroy(P7_SPARSEMX *sm)
{
  if (sm) {
    if (sm->kin)  free(sm->kin);
    if (sm->nin)  free(sm->nin);
    if (sm->kmem) free(sm->kmem);
    free(sm);
  }
}


    

