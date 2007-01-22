/* Routines for the P7_PROFILE structure - a Plan 7 search profile
 *                                         
 * SRE, Thu Jan 11 15:16:47 2007 [Janelia] [Sufjan Stevens, Illinois]
 * SVN $Id$
 */


/* Function:  p7_profile_Create()
 * Incept:    SRE, Thu Jan 11 15:53:28 2007 [Janelia]
 *
 * Purpose:   Creates a profile of <M> nodes, for digital alphabet <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    NULL on allocation error.
 *
 * Xref:      STL11/125.
 */
P7_PROFILE *
p7_profile_Create(int M, ESL_ALPHABET *abc)
{
  P7_PROFILE *gm = NULL;
  int         x;
  int         status;

  /* level 0 */
  ESL_ALLOC(gm, sizeof(P7_PROFILE));
  gm->tsc   = gm->msc = gm->isc = NULL;
  gm->bsc   = gm->esc = NULL;
  gm->begin = gm->end = NULL;
  
  /* level 1 */
  ESL_ALLOC(gm->tsc, sizeof(int *) * 7);           gm->tsc[0] = NULL;
  ESL_ALLOC(gm->msc, sizeof(int *) * (abc->Kp-1)); gm->msc[0] = NULL;
  ESL_ALLOC(gm->isc, sizeof(int *) * (abc->Kp-1)); gm->isc[0] = NULL;
  ESL_ALLOC(gm->bsc, sizeof(int) * (M+1));      
  ESL_ALLOC(gm->esc, sizeof(int) * (M+1));      

  /* Begin, end may eventually disappear in production
   * code, but we need them in research code for now to
   * be able to emulate & test HMMER2 configurations.
   */
  ESL_ALLOC(gm->begin, sizeof(float) * (M+1));
  ESL_ALLOC(gm->end,   sizeof(float) * (M+1));
  
  /* level 2 */
  ESL_ALLOC(gm->tsc[0], sizeof(int) * 7*M);    
  ESL_ALLOC(gm->msc[0], sizeof(int) * (abc->Kp-1) * (M+1));
  ESL_ALLOC(gm->isc[0], sizeof(int) * (abc->Kp-1) * M);
  for (x = 1; x < p7_MAXCODE; x++) {
    gm->msc[x] = gm->msc[0] + x * (M+1);
    gm->isc[x] = gm->isc[0] + x * M;
  }
  for (x = 0; x < 7; x++)
    gm->tsc[x] = gm->tsc[0] + x * M;

  gm->M    = M;
  gm->abc  = abc;
  gm->hmm  = NULL;
  gm->bg   = NULL;
  return gm;

 ERROR:
  p7_profile_Destroy(gm);
  return NULL;
}

/* Function:  p7_profile_Destroy()
 * Incept:    SRE, Thu Jan 11 15:54:17 2007 [Janelia]
 *
 * Purpose:   Frees a profile <gm>.
 *
 * Returns:   (void).
 *
 * Xref:      STL11/125.
 */
void
p7_profile_Destroy(P7_PROFILE *gm)
{
  if (gm != NULL) {
    if (gm->tsc   != NULL && gm->tsc[0] != NULL) free(gm->tsc[0]);
    if (gm->msc   != NULL && gm->msc[0] != NULL) free(gm->msc[0]);
    if (gm->isc   != NULL && gm->isc[0] != NULL) free(gm->isc[0]);
    if (gm->tsc   != NULL) free(gm->tsc);
    if (gm->msc   != NULL) free(gm->msc);
    if (gm->isc   != NULL) free(gm->isc);
    if (gm->bsc   != NULL) free(gm->bsc);
    if (gm->esc   != NULL) free(gm->esc);
    if (gm->begin != NULL) free(gm->begin);
    if (gm->end   != NULL) free(gm->end);
  }
  free(gm);
  return;
}
