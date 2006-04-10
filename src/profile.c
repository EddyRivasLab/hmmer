

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

