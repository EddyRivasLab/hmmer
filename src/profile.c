

P7_PROFILE *
p7_profile_Create(int M)
{
  P7_PROFILE *gm;
  int         x;

  gm          = MallocOrDie(sizeof(P7_PROFILE));
  gm->tsc     = MallocOrDie (7     *           sizeof(int *));
  gm->msc     = MallocOrDie (MAXCODE   *       sizeof(int *));
  gm->isc     = MallocOrDie (MAXCODE   *       sizeof(int *)); 
  gm->tsc[0]  = MallocOrDie ((7*M)     *       sizeof(int));
  gm->msc[0]  = MallocOrDie ((MAXCODE*(M+1)) * sizeof(int));
  gm->isc[0]  = MallocOrDie ((MAXCODE*M) *     sizeof(int));
  gm->bsc     = MallocOrDie ((M+1) *           sizeof(int));
  gm->esc     = MallocOrDie ((M+1) *           sizeof(int));
  
  for (x = 1; x < MAXCODE; x++) {
    gm->msc[x] = gm->msc[0] + x * (M+1);
    gm->isc[x] = gm->isc[0] + x * M;
  }
  for (x = 0; x < 7; x++)
    gm->tsc[x] = gm->tsc[0] + x * M;

  gm->M    = M;
  return gm;
}


void
p7_profile_Crutch(P7_HMM *hmm, P7_PROFILE *gm)
{
  int k, x;

  for (x = 0; x < 7; x++)
    for (k = 0; k < hmm->M; k++)
      gm->tsc[x][k] = hmm->tsc[x][k];
	
  for (x = 0; x < MAXCODE; x++)
    for (k = 0; k <= hmm->M; k++)
      gm->msc[x][k] = hmm->msc[x][k];

  for (x = 0; x < MAXCODE; x++)
    for (k = 0; k < hmm->M; k++)      
      gm->isc[x][k] = hmm->isc[x][k];

  for (k = 0; k < 4; k++)
    for (x = 0; x < 2; x++)
      gm->xsc[k][x] = hmm->xsc[k][x];

  for (k = 0; k <= hmm->M; k++) {
    gm->bsc[k] = hmm->bsc[k];
    gm->esc[k] = hmm->esc[k];
  }

  gm->M  = hmm->M;
  return;
}

void
p7_profile_Destroy(P7_PROFILE *gm)
{
  free(gm->esc);
  free(gm->bsc);
  free(gm->isc[0]);
  free(gm->msc[0]);
  free(gm->tsc[0]);
  free(gm->isc);
  free(gm->msc);
  free(gm->tsc);
  free(gm);
  return;
}
  
