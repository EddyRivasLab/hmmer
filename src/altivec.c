#ifdef ALTIVEC
void
AllocPlan7Body(struct plan7_s *hmm, int M) 
{
  int k, x;

  hmm->M = M;

  hmm->rf     = MallocOrDie ((M+2) * sizeof(char));
  hmm->cs     = MallocOrDie ((M+2) * sizeof(char));
  hmm->ca     = MallocOrDie ((M+2) * sizeof(char));
  hmm->map    = MallocOrDie ((M+1) * sizeof(int));

  hmm->t      = MallocOrDie (M     *           sizeof(float *));
  hmm->tsc    = MallocOrDie (7     *           sizeof(int *));
  hmm->mat    = MallocOrDie ((M+1) *           sizeof(float *));
  hmm->ins    = MallocOrDie (M     *           sizeof(float *));
  hmm->msc    = MallocOrDie (MAXCODE   *       sizeof(int *));
  hmm->isc    = MallocOrDie (MAXCODE   *       sizeof(int *)); 
  
  hmm->t[0]   = MallocOrDie ((7*M)     *       sizeof(float));
  /* Allocate extra memory so tsc[TMM,TIM,TDM,TMD,TDD] start on the 
   * 16-byte cache boundary, and tsc[TMI,TII] start
   * 12 bytes offset from the boundary. 
   */
  hmm->tsc_mem = MallocOrDie (((7*(M+16)))  *   sizeof(int));
  hmm->mat[0] = MallocOrDie ((MAXABET*(M+1)) * sizeof(float));
  hmm->ins[0] = MallocOrDie ((MAXABET*M) *     sizeof(float));
  /* Allocate extra mem. to make sure all members of msc,isc start
   * on 12-byte offsets from cache boundary.
   */
  hmm->msc_mem = MallocOrDie ((MAXCODE*(M+1+16)) * sizeof(int));
  hmm->isc_mem = MallocOrDie ((MAXCODE*(M+16)) *   sizeof(int));

  /* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   */
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * MAXABET;
    if (k < M) {
      hmm->ins[k] = hmm->ins[0] + k * MAXABET;
      hmm->t[k]   = hmm->t[0]   + k * 7;
    }
  }
  
  /* align tsc pointers */
  hmm->tsc[TMM] = (int *) (((((unsigned long int) hmm->tsc_mem) + 15) & (~0xf)));
  hmm->tsc[TMI] = (int *) (((((unsigned long int) hmm->tsc_mem) + (M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  hmm->tsc[TMD] = (int *) (((((unsigned long int) hmm->tsc_mem) + 2*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TIM] = (int *) (((((unsigned long int) hmm->tsc_mem) + 3*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TII] = (int *) (((((unsigned long int) hmm->tsc_mem) + 4*(M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  hmm->tsc[TDM] = (int *) (((((unsigned long int) hmm->tsc_mem) + 5*(M+12)*sizeof(int) + 15) & (~0xf)));
  hmm->tsc[TDD] = (int *) (((((unsigned long int) hmm->tsc_mem) + 6*(M+12)*sizeof(int) + 15) & (~0xf)));

  for (x = 0; x < MAXCODE; x++) {
    hmm->msc[x] = (int *) (((((unsigned long int)hmm->msc_mem) + x*(M+1+12)*sizeof(int) + 15) & (~0xf)) + 12);
    hmm->isc[x] = (int *) (((((unsigned long int)hmm->isc_mem) + x*(M+12)*sizeof(int) + 15) & (~0xf)) + 12);
  }
  /* tsc[0] is used as a boundary condition sometimes [Viterbi()],
   * so set to -inf always.
   */
  for (x = 0; x < 7; x++)
    hmm->tsc[x][0] = -INFTY;

  hmm->begin  = MallocOrDie  ((M+1) * sizeof(float));
  hmm->bsc_mem= MallocOrDie  ((M+1+12) * sizeof(int));
  hmm->end    = MallocOrDie  ((M+1) * sizeof(float));
  hmm->esc_mem= MallocOrDie  ((M+1+12) * sizeof(int));

  hmm->bsc = (int *) (((((unsigned long int) hmm->bsc_mem) + 15) & (~0xf)) + 12);
  hmm->esc = (int *) (((((unsigned long int) hmm->esc_mem) + 15) & (~0xf)) + 12);
  
  return;
}  


#endif /*ALTIVEC*/
