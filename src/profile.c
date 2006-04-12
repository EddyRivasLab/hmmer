

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
  
/* Function:  p7_profile_Dump()
 * Incept:    SRE, Thu Apr  6 13:46:25 2006 [AA890 enroute to Boston]
 *
 * Purpose:   Debugging: print log-odds scores of a configured plan7
 *            profile <gm> to a stream <fp>, in roughly the same format
 *            as a save file.  
 */
void
p7_profile_Dump(FILE *fp, P7_PROFILE *gm)
{
  char buf[p7_MAX_SC_TXTLEN];
  int k;			/* counter for nodes */
  int x;			/* counter for symbols */
  int ts;			/* counter for state transitions */
  
  fprintf(fp, "N: %6s ", p7_score2txt(gm->xsc[XTN][MOVE]), buf);
  fprintf(fp, "%6s\n",   p7_score2txt(gm->xsc[XTN][LOOP]), buf);

  for (k = 1; k <= gm->M; k++)
    {
				/* Line 1: k, match emissions */
      fprintf(fp, " %5d ", k);
      for (x = 0; x < gm->abc->K; x++) 
        fprintf(fp, "%6s ", p7_score2txt(gm->msc[x][k]), buf);
      fputs("\n", fp);
				/* Line 2: insert emissions */
      fprintf(fp, "       ");
      for (x = 0; x < gm->abc->K; x++) 
	fprintf(fp, "%6s ", (k < gm->M) ? p7_score2txt(gm->isc[x][k], buf) : "*");
      fputs("\n", fp);
				/* Line 3: transition probs; begin, end */
      fprintf(fp, "       ");
      for (ts = 0; ts < 7; ts++)
	fprintf(fp, "%6s ", (k < gm->M) ? p7_score2txt(gm->tsc[ts][k], buf) : "*"); 
      fprintf(fp, "%6s ", p7_score2txt(gm->bsc[k]), buf);
      fprintf(fp, "%6s ", p7_score2txt(gm->esc[k]), buf);
      fputs("\n", fp);
    }
  fprintf(fp, "E: %6s ", p7_score2txt(gm->xsc[XTE][MOVE]));
  fprintf(fp, "%6s\n",   p7_score2txt(gm->xsc[XTE][LOOP])); 

  fprintf(fp, "J: %6s ", p7_score2txt(gm->xsc[XTJ][MOVE]));
  fprintf(fp, "%6s\n",   p7_score2txt(gm->xsc[XTJ][LOOP])); 

  fprintf(fp, "C: %6s ", p7_score2txt(gm->xsc[XTC][MOVE]));
  fprintf(fp, "%6s\n",   p7_score2txt(gm->xsc[XTC][LOOP])); 

  fputs("//\n", fp);
}



/* Function:  p7_profile_GetTransition()
 * Incept:    SRE, Wed Apr 12 14:20:18 2006 [St. Louis]
 *
 * Purpose:   Convenience function that looks up a transition score in
 *            profile <gm> for a transition from state type <st1> in
 *            node <k1> to state type <st2> in node <k2>. For unique
 *            state types that aren't in nodes (<p7_STS>, for example), the
 *            <k> value is ignored; would be customarily passed as 0.
 *            Return the transition score in <ret_tsc>.
 *            
 * Returns:   <eslOK> on success.            
 * 
 * Throws:    <eslEINVAL> if a nonexistent transition is requested.
 */
int
p7_profile_GetTransition(P7_PROFILE *gm, char st1, int k1, char st2, int k2,
			 int *ret_tsc)
{
  int status = eslOK;
  int tsc    = 0;
  *ret_tsc   = 0;

  switch (st1) {
    case p7_STS: /* S->N is p=1 */
    break;	

  case p7_STN:
    switch (st2) {
    case p7_STB: tsc =  gm->xsc[p7_XTN][p7_MOVE]; 
    case p7_STN: tsc =  gm->xsc[p7_XTN][p7_LOOP]; 
    default:     status = eslEINVAL;
    }
    break;

  case p7_STB:
    switch (st2) {
    case p7_STM: tsc = gm->bsc[k2]; 
    default:     status = eslEINVAL;
    }
    break;

  case p7_STM:
    switch (st2) {
    case p7_STM: tsc = gm->tsc[p7_TMM][k1];
    case p7_STI: tsc = gm->tsc[p7_TMI][k1];
    case p7_STD: tsc = gm->tsc[p7_TMD][k1];
    case p7_STE: tsc = gm->esc[k1];
    default:     status = eslEINVAL;
    }
    break;

  case p7_STD:
    switch (st2) {
    case p7_STM: tsc = gm->tsc[p7_TDM][k1]; 
    case p7_STD: tsc = gm->tsc[p7_TDD][k1];
    default:     status = eslEINVAL;
    }
    break;

  case p7_STI:
    switch (st2) {
    case p7_STM: tsc = gm->tsc[p7_TIM][k1];
    case p7_STI: tsc = gm->tsc[p7_TII][k1];
    default:     status = eslEINVAL;
    }
    break;

  case p7_STE:
    switch (st2) {
    case p7_STC: tsc = gm->xsc[p7_XTE][p7_MOVE]; 
    case p7_STJ: tsc = gm->xsc[p7_XTE][p7_LOOP]; 
    default:     status = eslEINVAL;
    }
    break;

  case p7_STJ:
    switch (st2) {
    case p7_STB: tsc = gm->xsc[p7_XTJ][p7_MOVE]; 
    case p7_STJ: tsc = gm->xsc[p7_XTJ][p7_LOOP]; 
    default:     status = eslEINVAL;
    }
    break;

  case p7_STC:
    switch (st2) {
    case p7_STT:  tsc = gm->xsc[p7_XTC][p7_MOVE]; 
    case p7_STC:  tsc = gm->xsc[p7_XTC][p7_LOOP]; 
    default:      status = eslEINVAL;
    }
    break;

  case p7_STT:   
    break;

  default:
    esl_error(eslEINVAL, "illegal state typde %d in traceback", st1);
    return eslEINVAL;
  }

  if (status != eslOK)
    {
      esl_error(status, "illegal %s->%s transition", 
		p7_hmm_Statetype(st1), p7_hmm_Statetype(st2));
      return status;
    }
  return eslOK;
}
