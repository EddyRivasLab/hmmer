/* plan7.c
 * The Plan 7 HMM data structure, P7_HMM
 * 
 * SRE, Sat Nov 16 14:19:56 1996
 * SVN $Id: plan7.c 1487 2005-11-12 20:40:57Z eddy $
 */

#include "p7_config.h"		/* must be included first */

#include <stdio.h>
#include <string.h>		/* strcpy(), strlen()             */
#include <time.h>		/* p7_hmm_SetCtime() calls time() */

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_vectorops.h>

#include "plan7.h"		


/* Function:  p7_hmm_Create()
 * Incept:    SRE, Fri Mar 31 14:07:43 2006 [St. Louis]
 *
 * Purpose:   Allocate a <P7_HMM> of <M> nodes, for symbol
 *            alphabet <abc>, and return a pointer to it.
 *            
 *            The HMM only keeps a copy of the <abc> alphabet
 *            pointer. The caller is responsible for providing the
 *            alphabet, keeping it around while the HMM is in use,
 *            and (eventually) free'ing the alphabet when it's
 *            not needed any more. (Basically, just a step removed
 *            from keeping the alphabet as a global.)
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_HMM *
p7_hmm_Create(int M, ESL_ALPHABET *abc) 
{
  P7_HMM *hmm = NULL;

  if (hmm = p7_hmm_CreateShell();
  p7_hmm_CreateBody(hmm, M, abc);
  return hmm;
}  

/* Function:  p7_hmm_CreateShell()
 * Incept:    SRE, Fri Mar 31 14:09:45 2006 [St. Louis]
 *
 * Purpose:   Allocate the shell of a <P7_HMM>: everything that
 *            doesn't depend on knowing the number of nodes M. 
 *            
 *            HMM input (<hmmio.c>) uses two-step shell/body allocation
 *            because it has to read data before it reads the model
 *            size M.
 *
 * Returns:   a pointer to the new <P7_HMM> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_HMM *
p7_hmm_CreateShell(void) 
{
  P7_HMM *hmm = NULL;

  ESL_ALLOC(hmm, sizeof(P7_HMM));
  hmm->abc      = NULL;
  hmm->M        = 0;
  hmm->t        = NULL;
  hmm->mat      = NULL;
  hmm->ins      = NULL;

  hmm->null     = NULL;
  hmm->p1       = 0.;

  /* xt[][] left un-init */
  hmm->begin    = NULL;
  hmm->end      = NULL;

  hmm->gm       = NULL;
  hmm->om       = NULL;
  hmm->lscore   = 0.;

  hmm->name     = NULL;
  hmm->acc      = NULL;
  hmm->desc     = NULL;
  hmm->rf       = NULL;
  hmm->cs       = NULL;
  hmm->ca       = NULL;
  hmm->comlog   = NULL; 
  hmm->ctime    = NULL;
  hmm->map      = NULL;
  hmm->checksum = 0;

  hmm->tpri     = NULL;
  hmm->mpri     = NULL;
  hmm->ipri     = NULL;

  hmm->ga1 = hmm->ga2 = 0.;
  hmm->tc1 = hmm->tc2 = 0.;
  hmm->nc1 = hmm->nc2 = 0.;

  hmm->lvstats  = NULL;
  hmm->lfstats  = NULL;
  hmm->gvstats  = NULL;
  hmm->gfstats  = NULL;
  
  hmm->flags    = 0;

 CLEANEXIT:
  return hmm;
}  

/* Function:  p7_hmm_CreateBody()
 * Incept:    SRE, Fri Mar 31 14:24:44 2006 [St. Louis]
 *
 * Purpose:   Given an allocated shell <hmm>, and a now-known number
 *            of nodes <M> and alphabet <abc>, allocate
 *            the remainder of it for that many nodes.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, the entire
 *            HMM is free'd (including the shell).
 */
int
p7_hmm_CreateBody(P7_HMM *hmm, int M, ESL_ALPHABET *abc) 
{
  int k, x;

  hmm->abc    = abc;
  hmm->M      = M;

  ESL_ALLOC(hmm->t,      M                   * sizeof(float *));
  ESL_ALLOC(hmm->mat,    (M+1)               * sizeof(float *));
  ESL_ALLOC(hmm->ins,    M                   * sizeof(float *));
  ESL_ALLOC(hmm->t[0],   (7*M)               * sizeof(float));
  ESL_ALLOC(hmm->mat[0], (hmm->abc->K*(M+1)) * sizeof(float));
  ESL_ALLOC(hmm->ins[0], (hmm->abc->K*M)     * sizeof(float));

  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * hmm->abc->K;
    if (k < M) {
      hmm->ins[k] = hmm->ins[0] + k * hmm->abc->K;
      hmm->t[k]   = hmm->t[0]   + k * 7;
    }
  }

  ESL_ALLOC(hmm->null,  hmm->abc->K * sizeof(float));
  esl_vec_FSet(hmm->null, hmm->abc->K, 1./(float)hmm->abc->K); /* uniform */

  ESL_ALLOC(hmm->begin, (M+1) * sizeof(float));
  ESL_ALLOC(hmm->end,   (M+1) * sizeof(float));

  ESL_ALLOC(hmm->rf,  (M+2) * sizeof(char));
  ESL_ALLOC(hmm->cs,  (M+2) * sizeof(char));
  ESL_ALLOC(hmm->ca,  (M+2) * sizeof(char));
  ESL_ALLOC(hmm->map, (M+1) * sizeof(int));
  
  p7_hmm_ZeroCounts(hmm);

 CLEANEXIT:
  if (status != eslOK) p7_hmm_Destroy(hmm);
  return status;
}  


/* Function:  p7_hmm_Destroy()
 * Incept:    SRE, Fri Mar 31 15:13:25 2006 [St. Louis]
 *
 * Purpose:   Frees both the shell and body of an <hmm>.
 *            The <hmm> may be damaged (incompletely allocated),
 *            or even <NULL>.
 *
 * Note:      remember, leave abc alone. We just have a copy of its
 *            pointer, that application gave us.
 *
 * Returns:   (void).
 */
void
p7_hmm_Destroy(P7_HMM *hmm)
{
  if (hmm == NULL) return;

  if (hmm->mat     != NULL) {
    if (hmm->mat[0] != NULL) free(hmm->mat[0]);
    free(hmm->mat);
  }
  if (hmm->ins     != NULL) {
    if (hmm->ins[0] != NULL) free(hmm->ins[0]);
    free(hmm->ins);
  }
  if (hmm->t != NULL) {
    if (hmm->t[0] != NULL) free(hmm->t[0]);
    free(hmm->t);
  }
  if (hmm->null    != NULL) free(hmm->null);
  if (hmm->begin   != NULL) free(hmm->begin);
  if (hmm->end     != NULL) free(hmm->end);

  if (hmm->name    != NULL) free(hmm->name);
  if (hmm->acc     != NULL) free(hmm->acc);
  if (hmm->desc    != NULL) free(hmm->desc);
  if (hmm->rf      != NULL) free(hmm->rf);
  if (hmm->cs      != NULL) free(hmm->cs);
  if (hmm->ca      != NULL) free(hmm->ca);
  if (hmm->comlog  != NULL) free(hmm->comlog);
  if (hmm->ctime   != NULL) free(hmm->ctime);
  if (hmm->map     != NULL) free(hmm->map);

  if (hmm->tpri    != NULL) free(hmm->tpri);
  if (hmm->mpri    != NULL) free(hmm->mpri);
  if (hmm->ipri    != NULL) free(hmm->ipri);

  if (hmm->lvstats != NULL) p7_evinfo_Destroy(hmm->lvstats);
  if (hmm->lfstats != NULL) p7_evinfo_Destroy(hmm->lfstats);
  if (hmm->gvstats != NULL) p7_evinfo_Destroy(hmm->gvstats);
  if (hmm->gfstats != NULL) p7_evinfo_Destroy(hmm->gfstats);

  free(hmm);
  return;
}

/* Function: p7_hmm_ZeroCounts()
 * 
 * Purpose:  Zeros the counts/probabilities fields in a model (both
 *           core and configured form).  
 *           Leaves null model alone;
 *           invalidates any profile or optimized profile by dropping
 *           HASBITS flag, but leaves the memory allocated;
 *           invalidates any statistical fits by dropping STATS flags,
 *           but leaves memory allocated.
 */
void
p7_hmm_ZeroCounts(P7_HMM *hmm)
{
  int k;
  for (k = 0; k < hmm->M; k++)
    {
      esl_vec_FSet(hmm->t[k],   7,           0.);  /* t[0] uses only TMM,TMD */
      esl_vec_FSet(hmm->mat[k], hmm->abc->K, 0.);  /* mat[0] unused */
      esl_vec_FSet(hmm->ins[k], hmm->abc->K, 0.);  /* ins[0] unused */
    }
  esl_vec_FSet(hmm->mat[hmm->M], hmm->abc->K, 0.);

  esl_vec_FSet(hmm->begin, hmm->M+1, 0.); /* include unused begin[0] */
  esl_vec_FSet(hmm->end,   hmm->M+1, 0.); /* include unused begin[0] */
  for (k = 0; k < 4; k++)
    esl_vec_FSet(hmm->xt[k], 2, 0.);
  
  hmm->mode   =   P7_NO_MODE;
  hmm->flags  &= ~PLAN7_HASPROB;	/* invalidates probabilities        */
  hmm->flags  &= ~PLAN7_HASBITS;	/* invalidates scores               */
  hmm->flags  &= ~PLAN7_STATS_LV;	/* invalidates local Viterbi stats  */
  hmm->flags  &= ~PLAN7_STATS_LF;	/* invalidates local Forward stats  */
  hmm->flags  &= ~PLAN7_STATS_GV;	/* invalidates glocal Viterbi stats */
  hmm->flags  &= ~PLAN7_STATS_GF;	/* invalidates glocal Forward stats */
  return;
}


/*****************************************************************
 * 2. Convenience functions for setting annotation in the HMM.
 *****************************************************************/ 

/* Function: p7_hmm_SetName()
 * 
 * Purpose:  Set or change the name of a Plan7 HMM to <name>.
 *           Any trailing whitespace (including newline) is chopped off.     
 *      
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslEMEM> on allocation error, and original name (if any) 
 *           remains.
 */
int
p7_hmm_SetName(P7_HMM *hmm, char *name)
{
  int   status;
  void *tmp;
  int   n;

  if (name == NULL) 
    {
      if (hmm->name != NULL) free(hmm->name); 
      hmm->name = NULL;
    }
  else
    {
      n = strlen(name);

      if      (hmm->name == NULL)   
	ESL_ALLOC (hmm->name,      sizeof(char)*(n+1));
      else if (n > strlen(hmm->name))
	ESL_RALLOC(hmm->name, tmp, sizeof(char)*(n+1));

      strcpy(hmm->name, name);
      esl_strchop(hmm->name);
    }
 CLEANEXIT:
  return status;
}


/* Function: p7_hmm_SetAccession()
 * 
 * Purpose:  Set of change the accession number of a Plan7 HMM to <acc>. 
 *           Trailing whitespace (including newline) is chopped.     

 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslEMEM> on allocation error, and original name (if any) 
 *           remains.
 */
int
p7_hmm_SetAccession(P7_HMM *hmm, char *acc)
{
  int   status;
  void *tmp;
  int   n;

  if (acc == NULL) 
    {
      if (hmm->acc != NULL) free(hmm->acc); 
      hmm->acc = NULL;
    }
  else
    {
      n = strlen(acc);

      if      (hmm->acc == NULL)     
	ESL_ALLOC (hmm->acc,      sizeof(char)*(n+1));
      else if (n > strlen(hmm->acc)) 
	ESL_RALLOC(hmm->acc, tmp, sizeof(char)*(n+1));

      strcpy(hmm->acc, acc);
      esl_strchop(hmm->acc);
    }
  status = eslOK;
 CLEANEXIT:
  return status;
}

/* Function: p7_hmm_SetDescription()
 * 
 * Purpose:  Change the description line of a Plan7 HMM. 
 *           Trailing whitespace (including newline) is chopped.
 */
void
p7_hmm_SetDescription(P7_HMM *hmm, char *desc)
{
  int   status;
  void *tmp;
  int   n;

  if (desc == NULL) 
    {
      if (hmm->desc != NULL) free(hmm->desc); 
      hmm->desc = NULL;
      hmm->flags &= ~PLAN7_DESC;
    }
  else
    {
      n = strlen(desc);

      if      (hmm->desc == NULL)     
	ESL_ALLOC (hmm->desc,      sizeof(char)*(n+1));
      else if (n > strlen(hmm->desc)) 
	ESL_RALLOC(hmm->desc, tmp, sizeof(char)*(n+1));

      strcpy(hmm->desc, desc);
      esl_strchop(hmm->desc);
      hmm->flags |= PLAN7_DESC;
    }
  status = eslOK;
 CLEANEXIT:
  return status;
}

/* Function: p7_hmm_Comlog()
 * Date:     SRE, Wed Oct 29 09:57:30 1997 [TWA 721 over Greenland] 
 * 
 * Purpose:  Concatenate command line options and append as a line in the
 *           command line log. Command line log is multiline, with each line
 *           ending in newline char, except for last line.
 *           
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.          
 */
int
p7_hmm_Comlog(P7_HMM *hmm, int argc, char **argv)
{
  int   status;
  void *tmp;
  int   len;
  int   i;

  /* figure out length of added command line, and (re)allocate comlog */
  len = argc;	/* account for 1 space or \n per arg */
  for (i = 0; i < argc; i++)
    len += strlen(argv[i]);

  if (hmm->comlog != NULL)
    {
      len += strlen(hmm->comlog); /* last comlog already ends w/ \n. */
      ESL_RALLOC(hmm->comlog, tmp, sizeof(char)* (len+1));
    }
  else
    {
      ESL_ALLOC(hmm->comlog, sizeof(char)* (len+1));
      *(hmm->comlog) = '\0'; /* need this to make strcat work */
    }

  strcat(hmm->comlog, "\n");
  for (i = 0; i < argc; i++)
    {
      strcat(hmm->comlog, argv[i]);
      if (i < argc-1) strcat(hmm->comlog, " ");
    }

  status = eslOK;
 CLEANEXIT:
  return status;
}

/* Function: p7_hmm_SetCtime()
 * Date:     SRE, Wed Oct 29 11:53:19 1997 [TWA 721 over the Atlantic]
 * 
 * Purpose:  Set the <ctime> field in a new HMM to the current time.
 *           This function is not reentrant and not threadsafe, because
 *           it calls the nonreentrant ANSI C ctime() function.
 * 
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.
 */
int
p7_hmm_SetCtime(P7_HMM *hmm)
{
  time_t date = time(NULL);
  if (hmm->ctime != NULL) free(hmm->ctime);
  if (esl_strdup(ctime(&date), -1, &(hmm->ctime)) != eslOK)
    { hmm->ctime = NULL; return eslEMEM; }
  esl_strchop(hmm->ctime, -1);
  return eslOK;
}


/* Function: p7_hmm_SetNull()
 * 
 * Purpose:  Set the null model emission probabilities of an HMM.
 *           Null model <p1> is set by length modeling; see P7ReconfigLength().
 */
int
p7_hmm_SetNull(P7_HMM *hmm, float *null, int K)
{
  if (K != hmm->abc->K) 
    ESL_ERROR(eslEINVAL, "null model, hmm alphabet sizes differ");
  esl_vec_FCopy(hmm->null, null, hmm->abc->K);
  return eslOK;
}



/* Function:  p7_hmm_Rescale()
 * Incept:    Steve Johnson, 3 May 2004
 *            eweights code incorp: SRE, Thu May 20 10:34:03 2004 [St. Louis]
 *
 * Purpose:   Scale a counts-based HMM by some factor, for
 *            adjusting counts to a new effective sequence number.
 *            Only affects the core probability model (<t>, <ins>, and <mat>).
 *
 * Args:      hmm        - counts based HMM.
 *            scale      - scaling factor (e.g. eff_nseq/nseq); 1.0=no scaling.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_hmm_Rescale(P7_HMM *hmm, float scale)
{
  int k;

  for (k = 1; k <= hmm->M; k++) 
    esl_vec_FScale(hmm->mat[k], hmm->abc->K, scale);
  for (k = 1; k <  hmm->M; k++) 
    esl_vec_FScale(hmm->ins[k], hmm->abc->K, scale);
  for (k = 0; k <  hmm->M; k++) /* including begin->M1 and begin->D1 */
    esl_vec_FScale(hmm->t[k],   7,      scale);
  return eslOK;
}


/* Function: p7_hmm_Renormalize()
 * 
 * Purpose:  Take an HMM in counts form, and renormalize
 *           all probability vectors in the core probability model. Enforces
 *           Plan7 restrictions on nonexistent transitions. Sets the
 *           <PLAN7_HASPROB> flag. 
 *
 *           Leaves other flags (stats and profile) alone, so caller
 *           needs to be wary. Renormalizing a probability model that
 *           has stats and profile scores wouldn't usually invalidate
 *           those data; and if we're renormalizing a counts model, we
 *           shouldn't have stats or profile scores yet anyway.
 *           
 * Args:     hmm - the model to renormalize.
 *                 
 * Return:   <eslOK> on success.
 */                          
int
p7_hmm_Renormalize(P7_HMM *hmm)
{
  int   k;			/* counter for model position */
  float d;			/* denominator */

  for (k = 1; k <= hmm->M; k++)  /* match emissions: 1..M */
    esl_vec_FNorm(hmm->mat[k], hmm->abc->K);
  for (k = 1; k < hmm->M; k++)	/* insert emissions: 1..M-1 */
    esl_vec_FNorm(hmm->ins[k], hmm->abc->K);
  for (k = 1; k < hmm->M; k++)	/* transitions: 1..M-1 */
    {
      esl_vec_FNorm(hmm->t[k],   3);	/* match  */
      esl_vec_FNorm(hmm->t[k]+3, 2);	/* insert */
      esl_vec_FNorm(hmm->t[k]+5, 2);	/* delete */
    }

  hmm->t[0][TMI] = 0.;          /* make sure... */
  esl_vec_FNorm(hmm->t[0], 3);  /* begin transitions; TMM, TMD are valid */

  /* Enforce nonexistent but allocated transitions: */
  esl_vec_FSet(hmm->mat[0], hmm->abc->K, 0.); /* mat[0] */
  esl_vec_FSet(hmm->ins[0], hmm->abc->K, 0.); /* ins[0] */
  esl_vec_FSet(hmm->t[0]+3, 2, 0.); /* t[0] delete */
  esl_vec_FSet(hmm->t[0]+5, 2, 0.); /* t[0] insert */

  hmm->flags |= PLAN7_HASPROB;	/* set the probabilities OK flag */
  return eslOK;
}
  

    
/* Function:  p7_hmm_Dump()
 * Incept:    SRE, Sun May  8 08:34:55 2005 [St. Louis]
 *
 * Purpose:   Debugging: dump the probabilities (or counts) from an HMM.
 */
void
p7_hmm_Dump(FILE *fp, P7_HMM *hmm)
{
  int k;			/* counter for nodes */
  int x;			/* counter for symbols */
  int ts;			/* counter for state transitions */
  
  fprintf(fp, "B->M1,B->D1: %5.1f %5.1f\n", hmm->t[0][TMM], hmm->t[0][TMD]);
  for (k = 1; k <= hmm->M; k++)
    {
				/* Line 1: k, match emissions */
      fprintf(fp, " %5d ", k);
      for (x = 0; x < hmm->abc->K; x++) 
        fprintf(fp, "%5.1f ", hmm->mat[k][x]);
      fputs("\n", fp);
				/* Line 2: insert emissions */
      fprintf(fp, "       ");
      for (x = 0; x < hmm->abc->K; x++) 
	fprintf(fp, "%5.1f ", (k < hmm->M) ? hmm->ins[k][x] : 0);
      fputs("\n", fp);
				/* Line 3: transition probs */
      fprintf(fp, "       ");
      for (ts = 0; ts < 7; ts++)
	fprintf(fp, "%5.1f ", (k < hmm->M) ? hmm->t[k][ts] : 0); 
      fputs("\n", fp);
    }
  fputs("//\n", fp);
}



/************************************************************
 * @LICENSE@
 ************************************************************/


