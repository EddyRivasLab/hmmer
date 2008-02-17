/* The Plan7 core HMM data structure.
 * 
 * Contents:
 *   1. The P7_HMM object: allocation, initialization, destruction.
 *   2. Convenience routines for setting fields in an HMM.
 *   3. Renormalization and rescaling counts in core HMMs.
 *   4. Debugging and development code.
 *   5. Other routines in the API.
 *   6. Unit tests.
 *   7. Test driver. 
 *   8. Copyright and license.
 * 
 * SRE, Mon Jan  1 16:20:29 2007 [Casa de Gatos] [Verdi, La Traviata]
 * SVN $Id$
 */

#include "p7_config.h"		/* must be included first */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"
#include "esl_random.h"
#include "esl_dirichlet.h"

#include "hmmer.h"




/*****************************************************************
 * 1. The P7_HMM object: allocation, initialization, destruction.
 *****************************************************************/

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
p7_hmm_Create(int M, const ESL_ALPHABET *abc) 
{
  P7_HMM *hmm = NULL;

  if ((hmm = p7_hmm_CreateShell()) == NULL) return NULL;
  p7_hmm_CreateBody(hmm, M, abc);
  return hmm;
}  


/* Function:  p7_hmm_CreateShell()
 * Incept:    SRE, Fri Mar 31 14:09:45 2006 [St. Louis]
 *
 * Purpose:   Allocate the shell of a <P7_HMM>: everything that
 *            doesn't depend on knowing the number of nodes M. 
 *            
 *            HMM input (<hmmio.c>) uses two-step shell/body
 *            allocation because it has to read for a ways from the
 *            HMM file before it reads the model size M or the
 *            alphabet type.
 *
 * Returns:   a pointer to the new <P7_HMM> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_HMM *
p7_hmm_CreateShell(void) 
{
  P7_HMM *hmm = NULL;
  int     z;
  int     status;

  ESL_ALLOC(hmm, sizeof(P7_HMM));
  hmm->M        = 0;
  hmm->t        = NULL;
  hmm->mat      = NULL;
  hmm->ins      = NULL;

  hmm->name     = NULL;
  hmm->acc      = NULL;
  hmm->desc     = NULL;
  hmm->rf       = NULL;
  hmm->cs       = NULL;
  hmm->ca       = NULL;
  hmm->comlog   = NULL; 
  hmm->nseq     = 0;
  hmm->eff_nseq = 0;
  hmm->ctime    = NULL;
  hmm->map      = NULL;
  hmm->checksum = 0;

  for (z = 0; z < p7_NCUTOFFS; z++) hmm->cutoff[z]  = 0.0f;
  for (z = 0; z < p7_NEVPARAM; z++) hmm->evparam[z] = 0.0f;

  hmm->offset   = 0;
  hmm->flags    = 0;
  hmm->abc      = NULL;
  return hmm;

 ERROR:
  return NULL;
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
 * Throws:    <eslEMEM> on allocation failure; in this case, the HMM
 *            is likely corrupted, and the caller should destroy it.
 */
int
p7_hmm_CreateBody(P7_HMM *hmm, int M, const ESL_ALPHABET *abc) 
{
  int k;
  int status;

  hmm->abc = abc;
  hmm->M   = M;

  /* level 1 */
  ESL_ALLOC(hmm->t,    (M+1) * sizeof(float *));
  ESL_ALLOC(hmm->mat,  (M+1) * sizeof(float *));
  ESL_ALLOC(hmm->ins,  (M+1) * sizeof(float *));
  hmm->t[0]   = NULL;
  hmm->mat[0] = NULL;
  hmm->ins[0] = NULL;

  /* level 2 */
  ESL_ALLOC(hmm->t[0],   (p7H_NTRANSITIONS*(M+1)) * sizeof(float));
  ESL_ALLOC(hmm->mat[0], (abc->K*(M+1))           * sizeof(float));
  ESL_ALLOC(hmm->ins[0], (abc->K*(M+1))           * sizeof(float));
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * hmm->abc->K;
    hmm->ins[k] = hmm->ins[0] + k * hmm->abc->K;
    hmm->t[k]   = hmm->t[0]   + k * p7H_NTRANSITIONS;
  }

  /* Enforce conventions on unused but allocated distributions, so
   * Compare() tests succeed unless memory was corrupted.
   */
  if ((status = p7_hmm_Zero(hmm)) != eslOK) goto ERROR;
  hmm->mat[0][0]    = 1.0;
  hmm->t[0][p7H_DM] = 1.0; 

  /* Optional allocation, status flag dependent */
  if (hmm->flags & p7H_RF)  ESL_ALLOC(hmm->rf,  (M+2) * sizeof(char));
  if (hmm->flags & p7H_CS)  ESL_ALLOC(hmm->cs,  (M+2) * sizeof(char));
  if (hmm->flags & p7H_CA)  ESL_ALLOC(hmm->ca,  (M+2) * sizeof(char));
  if (hmm->flags & p7H_MAP) ESL_ALLOC(hmm->map, (M+1) * sizeof(int));
  
  return eslOK;

 ERROR:
  return status;
}  


/* Function:  p7_hmm_Destroy()
 * Incept:    SRE, Fri Mar 31 15:13:25 2006 [St. Louis]
 *
 * Purpose:   Frees both the shell and body of an <hmm>.
 *            Works even if the <hmm> is damaged (incompletely allocated)
 *            or even <NULL>.
 *
 * Note:      Remember, leave reference pointers like abc, gm, and
 *            bg alone. These are under the application's control not ours.
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

  if (hmm->name    != NULL) free(hmm->name);
  if (hmm->acc     != NULL) free(hmm->acc);
  if (hmm->desc    != NULL) free(hmm->desc);
  if (hmm->rf      != NULL) free(hmm->rf);
  if (hmm->cs      != NULL) free(hmm->cs);
  if (hmm->ca      != NULL) free(hmm->ca);
  if (hmm->comlog  != NULL) free(hmm->comlog);
  if (hmm->ctime   != NULL) free(hmm->ctime);
  if (hmm->map     != NULL) free(hmm->map);

  free(hmm);
  return;
}

/* Function:  p7_hmm_CopyParameters()
 * Incept:    SRE, Fri May  4 14:10:17 2007 [Janelia]
 *
 * Purpose:   Copy parameters of <src> to <dest>. The HMM <dest> must
 *            be allocated by the caller for the same 
 *            alphabet and M as <src>. 
 *            
 *            Both core and search model parameters are copied.
 *
 *            No annotation is copied.  This is because several
 *            annotation fields are variable-length strings that
 *            require individual allocations.  The
 *            <p7_hmm_CopyParameters()> function is for cases where we
 *            have to repeatedly reset the parameters of a model - for
 *            example, in entropy weighting.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_hmm_CopyParameters(const P7_HMM *src, P7_HMM *dest)
{
  int k;
  for (k = 0; k <= src->M; k++) {
    esl_vec_FCopy(src->t[k],   p7H_NTRANSITIONS, dest->t[k]);
    esl_vec_FCopy(src->mat[k], src->abc->K,      dest->mat[k]);
    esl_vec_FCopy(src->ins[k], src->abc->K,      dest->ins[k]);
  }
  return eslOK;
}

/* Function:  p7_hmm_Clone()
 * Incept:    SRE, Fri Jan 26 15:34:42 2007 [Janelia]
 *
 * Purpose:   Duplicates an hmm.
 * 
 *            Note: does not duplicate the objects the HMM refers to,
 *            if any (profile, null model, or alphabet); only copies
 *            the reference pointers.
 * 
 * Returns:   a pointer to the duplicate.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_HMM *
p7_hmm_Clone(const P7_HMM *hmm)
{
  int     status;
  P7_HMM *new = NULL;
  int     z;

  if ((new = p7_hmm_Create(hmm->M, hmm->abc)) == NULL) goto ERROR;
  p7_hmm_CopyParameters(hmm, new);
  
  if (hmm->name != NULL    && (status = esl_strdup(hmm->name,   -1, &(new->name)))   != eslOK) goto ERROR;
  if (hmm->acc  != NULL    && (status = esl_strdup(hmm->acc,    -1, &(new->acc)))    != eslOK) goto ERROR;
  if (hmm->desc != NULL    && (status = esl_strdup(hmm->desc,   -1, &(new->desc)))   != eslOK) goto ERROR;
  if (hmm->flags & p7H_RF  && (status = esl_strdup(hmm->rf,     -1, &(new->rf)))     != eslOK) goto ERROR;
  if (hmm->flags & p7H_CS  && (status = esl_strdup(hmm->cs,     -1, &(new->cs)))     != eslOK) goto ERROR;
  if (hmm->flags & p7H_CA  && (status = esl_strdup(hmm->ca,     -1, &(new->ca)))     != eslOK) goto ERROR;
  if (hmm->comlog != NULL  && (status = esl_strdup(hmm->comlog, -1, &(new->comlog))) != eslOK) goto ERROR;
  if (hmm->ctime  != NULL  && (status = esl_strdup(hmm->ctime,  -1, &(new->ctime)))  != eslOK) goto ERROR;
  if (hmm->flags & p7H_MAP) {
    ESL_ALLOC(new->map, sizeof(int) * (hmm->M+1));
    esl_vec_ICopy(hmm->map, hmm->M+1, new->map);
  }
  new->nseq     = hmm->nseq;
  new->eff_nseq = hmm->eff_nseq;
  new->checksum = hmm->checksum;

  for (z = 0; z < p7_NEVPARAM; z++) new->evparam[z] = hmm->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) new->cutoff[z]  = hmm->cutoff[z];

  new->offset   = hmm->offset;
  new->flags    = hmm->flags;
  new->abc      = hmm->abc;
  return new;

 ERROR:
  if (new != NULL) p7_hmm_Destroy(new);
  return NULL;
}

/* Function:  p7_hmm_Scale()
 * Incept:    SRE, Fri May  4 14:19:33 2007 [Janelia]
 *
 * Purpose:   Given a counts-based model <hmm>, scale core
 *            by a multiplicative factor of <scale>. Used in
 *            absolute sequence weighting.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_hmm_Scale(P7_HMM *hmm, double scale)
{
  int k;

  for (k = 0; k <= hmm->M; k++) {
    esl_vec_FScale(hmm->t[k],   p7H_NTRANSITIONS, scale);  
    esl_vec_FScale(hmm->mat[k], hmm->abc->K,      scale);  
    esl_vec_FScale(hmm->ins[k], hmm->abc->K,      scale);  
  }
  return eslOK;
}


/* Function:  p7_hmm_Zero()
 * Incept:    SRE, Mon Jan  1 16:32:59 2007 [Casa de Gatos]
 *
 * Purpose:   Zeroes the counts/probabilities fields in core model.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_hmm_Zero(P7_HMM *hmm)
{
  int k;

  for (k = 0; k <= hmm->M; k++) {
    esl_vec_FSet(hmm->t[k],   p7H_NTRANSITIONS, 0.);  
    esl_vec_FSet(hmm->mat[k], hmm->abc->K, 0.);  
    esl_vec_FSet(hmm->ins[k], hmm->abc->K, 0.);  
  }
  return eslOK;
}



/* Function:  p7_hmm_DescribeStatetype()
 * Incept:    SRE, Mon Jan  1 18:47:34 2007 [Casa de Gatos]
 *
 * Purpose:   Returns the state type in text, as a string of length 1 
 *            (2 if you count NUL). For example, <p7_Statetype(p7T_S)>
 *            returns "S".
 */
char *
p7_hmm_DescribeStatetype(char st)
{
  switch (st) {
  case p7T_M: return "M";
  case p7T_D: return "D";
  case p7T_I: return "I";
  case p7T_S: return "S";
  case p7T_N: return "N";
  case p7T_B: return "B";
  case p7T_E: return "E";
  case p7T_C: return "C";
  case p7T_T: return "T";
  case p7T_J: return "J";
  default:     return "?";
  }
}




/*****************************************************************
 * 2. Convenience routines for setting fields in an HMM.
 *****************************************************************/ 

/* Function: p7_hmm_SetName()
 * Incept:   SRE, Mon Jan  1 16:53:23 2007 [Casa de Gatos]
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

  if (name == NULL) {
    if (hmm->name != NULL) free(hmm->name); 
    hmm->name = NULL;
  } else {
    n = strlen(name);
    ESL_RALLOC(hmm->name, tmp, sizeof(char)*(n+1));
    strcpy(hmm->name, name);
    if ((status = esl_strchop(hmm->name, n)) != eslOK) goto ERROR;
  }
  return eslOK;

 ERROR:
  return status;
}

/* Function: p7_hmm_SetAccession()
 * Incept:   SRE, Mon Jan  1 16:53:53 2007 [Casa de Gatos]
 * 
 * Purpose:  Set or change the accession number of a Plan7 HMM to <acc>,
 *           and raise the <P7_ACC> flag. Trailing whitespace (including newline) 
 *           is chopped.  
 *           
 *           If <acc> is <NULL>, unset the HMM's accession (if any) and drop 
 *           the <P7_ACC> flag.
 *
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

  if (acc == NULL) {
    if (hmm->acc != NULL) free(hmm->acc); 
    hmm->acc = NULL;
    hmm->flags &= ~p7H_ACC;
  } else {
    n = strlen(acc);
    ESL_RALLOC(hmm->acc, tmp, sizeof(char)*(n+1));
    strcpy(hmm->acc, acc);
    if ((status = esl_strchop(hmm->acc, n)) != eslOK) goto ERROR;
    hmm->flags |= p7H_ACC;
  }
  return eslOK;

 ERROR:
  return status;
}

/* Function: p7_hmm_SetDescription()
 * Incept:   SRE, Mon Jan  1 16:59:28 2007 [Casa de Gatos]
 * 
 * Purpose:  Set or change the description line of a Plan7 HMM. 
 *           Trailing whitespace (including newline) is chopped.
 */
int
p7_hmm_SetDescription(P7_HMM *hmm, char *desc)
{
  int   status;
  void *tmp;
  int   n;

  if (desc == NULL) 
    {
      if (hmm->desc != NULL) free(hmm->desc); 
      hmm->desc   = NULL;
      hmm->flags &= ~p7H_DESC;
    }
  else
    {
      n = strlen(desc);
      ESL_RALLOC(hmm->desc, tmp, sizeof(char)*(n+1));
      strcpy(hmm->desc, desc);
      if ((status = esl_strchop(hmm->desc, n)) != eslOK) goto ERROR;
      hmm->flags |= p7H_DESC;
    }
  return eslOK;

 ERROR:
  return status;
}

/* Function: p7_hmm_AppendComlog()
 * Incept:   SRE, Mon Jan  1 18:23:42 2007 [Casa de Gatos]
 * 
 * Purpose:  Concatenate command line options and append as a new line in the
 *           command line log. Command line log is multiline, with each line
 *           ending in newline char, except for last line.
 *           
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.          
 */
int
p7_hmm_AppendComlog(P7_HMM *hmm, int argc, char **argv)
{
  int   status;
  void *tmp;
  int   n;
  int   i;

  /* figure out length of added command line, and (re)allocate comlog */
  n = argc-1;	/* account for 1 space per arg, except last one */
  for (i = 0; i < argc; i++)
    n += strlen(argv[i]);

  if (hmm->comlog != NULL) {
    n += strlen(hmm->comlog) + 1; /* +1 for the \n we're going to add to the old comlog */
    ESL_RALLOC(hmm->comlog, tmp, sizeof(char)* (n+1));
    strcat(hmm->comlog, "\n");
  } else {
    ESL_ALLOC(hmm->comlog, sizeof(char)* (n+1));
    *(hmm->comlog) = '\0'; /* need this to make strcat work */
  }

  for (i = 0; i < argc-1; i++)
    {
      strcat(hmm->comlog, argv[i]);
      strcat(hmm->comlog, " ");
    }
  strcat(hmm->comlog, argv[argc-1]);
  return eslOK;

 ERROR:
  return status;
}

/* Function: p7_hmm_SetCtime()
 * Date:     SRE, Wed Oct 29 11:53:19 1997 [TWA 721 over the Atlantic]
 * 
 * Purpose:  Set the <ctime> field in a new HMM to the current time.
 *
 *           This function is not reentrant and not threadsafe, because
 *           it calls the nonreentrant ANSI C ctime() function.
 * 
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure. <eslESYS> if the <time()>
 *           system call fails to obtain the calendar time.
 */
int
p7_hmm_SetCtime(P7_HMM *hmm)
{
  int    status;
  char  *s = NULL;
  time_t date;

  if ((date   = time(NULL))                       == -1) { status = eslESYS; goto ERROR; }
  if ((status = esl_strdup(ctime(&date), -1, &s)) != eslOK) goto ERROR;
  if ((status = esl_strchop(s, -1))               != eslOK) goto ERROR;
  
  if (hmm->ctime != NULL) free(hmm->ctime);
  hmm->ctime = s;
  return eslOK;

 ERROR:
  if (s != NULL) free(s);
  return status;
}
/*---------------- end, internal-setting routines ---------------*/




/*****************************************************************
 * 3. Renormalization and rescaling counts in core HMMs.
 *****************************************************************/ 

/* Function:  p7_hmm_Rescale()
 * Incept:    Steve Johnson, 3 May 2004
 *            eweights code incorp: SRE, Thu May 20 10:34:03 2004 [St. Louis]
 *
 * Purpose:   Scale a counts-based core HMM by some factor, for
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

  for (k = 0; k <= hmm->M; k++) {
    esl_vec_FScale(hmm->mat[k], hmm->abc->K, scale);
    esl_vec_FScale(hmm->ins[k], hmm->abc->K, scale);
    esl_vec_FScale(hmm->t[k],   7,           scale);
  }
  return eslOK;
}

/* Function: p7_hmm_Renormalize()
 * Incept:   SRE, Mon Jan  1 18:39:42 2007 [Casa de Gatos]
 * 
 * Purpose:  Take a core HMM in counts form, and renormalize
 *           all probability vectors in the core probability model. Enforces
 *           Plan7 restrictions on nonexistent transitions.
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

  for (k = 0; k <= hmm->M; k++) {
    esl_vec_FNorm(hmm->mat[k], hmm->abc->K);
    esl_vec_FNorm(hmm->ins[k], hmm->abc->K);
    esl_vec_FNorm(P7H_TMAT(hmm, k), p7H_NTMAT);	/* TMX */
    esl_vec_FNorm(P7H_TDEL(hmm, k), p7H_NTDEL);	/* TIX */
    esl_vec_FNorm(P7H_TINS(hmm, k), p7H_NTINS);	/* TDX */
  }
  /* If t[M][TD*] distribution was all zeros, we just made TDD nonzero. Oops.
   * Re-enforce t's on that final delete state. */
  hmm->t[hmm->M][p7H_DM] = 1.0;
  hmm->t[hmm->M][p7H_DD] = 0.0;

  /* Rare: if t[M][TM*] distribution was all zeros (all final transitions
   * were D_M -> E) then we just made nonexistent M_M->D_M+1 transition nonzero.
   * Fix that too.
   */
  if (hmm->t[hmm->M][p7H_MD] > 0.) {
    hmm->t[hmm->M][p7H_MD] = 0.;
    hmm->t[hmm->M][p7H_MM] = 0.5;
    hmm->t[hmm->M][p7H_MI] = 0.5;
  }

  return eslOK;
}
  
/*****************************************************************
 * 4. Debugging and development code
 *****************************************************************/

/* Function:  p7_hmm_Dump()
 * Incept:    SRE, Mon Jan  1 18:44:15 2007 [Casa de Gatos]
 *
 * Purpose:   Debugging: dump the probabilities (or counts) from a core HMM.
 * 
 * Returns:   <eslOK> on success.
 */
int
p7_hmm_Dump(FILE *fp, P7_HMM *hmm)
{
  int k;			/* counter for nodes */
  int x;			/* counter for symbols */
  int ts;			/* counter for state transitions */
  
  for (k = 0; k <= hmm->M; k++)
    {				/* Line 1: k, match emissions */
      fprintf(fp, " %5d ", k);
      for (x = 0; x < hmm->abc->K; x++) 
        fprintf(fp, "%9.4f ", hmm->mat[k][x]);
      fputs("\n", fp);
				/* Line 2: insert emissions */
      fprintf(fp, "       ");
      for (x = 0; x < hmm->abc->K; x++) 
	fprintf(fp, "%9.4f ", hmm->ins[k][x]);
      fputs("\n", fp);
				/* Line 3: transition probs */
      fprintf(fp, "       ");
      for (ts = 0; ts < 7; ts++)
	fprintf(fp, "%9.4f ", hmm->t[k][ts]); 
      fputs("\n", fp);
    }
  fputs("//\n", fp);
  return eslOK;
}

/* Function:  p7_hmm_Sample()
 * Incept:    SRE, Sat Jan  6 13:43:03 2007 [Casa de Gatos]
 *
 * Purpose:   Creates a random HMM of length <M> nodes,
 *            for alphabet <abc>, obtaining randomness from
 *            <r>.
 * 
 *            Probably only useful for debugging.
 *            
 * Note:      Compare p7_hmm_Renormalize(), which has a similar
 *            structure, except it normalizes instead of
 *            sampling each probability vector.           
 *
 * Returns:   <eslOK> on success, and the new hmm is returned
 *            through <ret_hmm); caller is responsible for 
 *            freeing this object with <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_hmm_Sample(ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[random HMM created by sampling]";
  int     k;
  int     status;

  hmm = p7_hmm_Create(M, abc);
  if (hmm == NULL) { status = eslEMEM; goto ERROR; }
  
  for (k = 0; k <= M; k++)
    {
      if (k > 0) esl_dirichlet_FSampleUniform(r, abc->K, hmm->mat[k]);
      esl_dirichlet_FSampleUniform(r, abc->K, hmm->ins[k]);
      esl_dirichlet_FSampleUniform(r, 3,      hmm->t[k]);
      esl_dirichlet_FSampleUniform(r, 2,      hmm->t[k]+3);
      if (k > 0) esl_dirichlet_FSampleUniform(r, 2,      hmm->t[k]+5);
    }
  /* Node M is special: no transitions to D, transitions to M
   * are interpreted as transitions to E. Overwrite a little of
   * what we did in node M.
   */
  esl_dirichlet_FSampleUniform(r, 2, hmm->t[M]);    /* TMM,TMI only */
  hmm->t[M][p7H_MD] = 0.;	
  hmm->t[M][p7H_DM] = 1.0;
  hmm->t[M][p7H_DD] = 0.0;
  
  /* Add mandatory annotation
   */
  p7_hmm_SetName(hmm, "sampled-hmm");
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  hmm->nseq     = 0;
  hmm->eff_nseq = 0;
  p7_hmm_SetCtime(hmm);
  hmm->checksum = 0;

  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;

}

/* Function:  p7_hmm_SampleUngapped()
 * Incept:    SRE, Thu Jan 25 09:38:30 2007 [Janelia]
 *
 * Purpose:   Same as <p7_hmm_Sample()>, except all 
 *            M $\rightarrow$ M transitions are 1.0:
 *            an ungapped model. Useful for testing 
 *            as a limit case.
 *            
 * Returns:   <eslOK> on success, and the new hmm is returned
 *            through <ret_hmm); caller is responsible for 
 *            freeing this object with <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      STL11/140
 */
int
p7_hmm_SampleUngapped(ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  int     k;
  int     status;

  if ((status = p7_hmm_Sample(r, M, abc, &hmm)) != eslOK) goto ERROR;
  for (k = 0; k <= M; k++) {
    hmm->t[k][p7H_MM] = 1.0;
    hmm->t[k][p7H_MD] = 0.0;
    hmm->t[k][p7H_MI] = 0.0;
  }
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}

/* Function:  esl_hmm_SampleEnumerable()
 * Incept:    SRE, Wed Apr 18 09:38:09 2007 [Janelia]
 *
 * Purpose:   Sample a random HMM with random emission and 
 *            transition probabilities with the exception that
 *            all transitions to insert are zero. This makes
 *            it possible to create a model with a finite,
 *            easily enumerable sequence space (all seqs of
 *            length $\leq M).
 *            
 *            To achieve this in the profile as well as the core HMM,
 *            the caller must configure a unihit mode
 *            (<p7_ProfileConfig(hmm, bg, gm, p7_UNILOCAL)> or
 *            <p7_UNIGLOCAL>), and a target length of zero
 *            (<p7_ReconfigLength(gm, 0)>).
 *            
 *            Useful for debugging and validating Forward/Viterbi
 *            algorithms.
 *            
 * Returns:   <eslOK> on success. The newly allocated hmm is returned through
 *            <ret_hmm>. The caller is responsible for freeing this object
 *            with <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_hmm_SampleEnumerable(ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, P7_HMM **ret_hmm)
{
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[random enumerable HMM created by sampling]";
  int     k;
  float   tmp[2];
  int     status;
  
  hmm = p7_hmm_Create(M, abc);
  if (hmm == NULL) { status = eslEMEM; goto ERROR; }

  for (k = 0; k <= M; k++)
    {
      if (k > 0) esl_dirichlet_FSampleUniform(r, abc->K, hmm->mat[k]); /* match emission probs  */
      esl_dirichlet_FSampleUniform(r, abc->K, hmm->ins[k]);            /* insert emission probs */
      esl_dirichlet_FSampleUniform(r, 2,      tmp);       
      hmm->t[k][p7H_MM] = tmp[0];
      hmm->t[k][p7H_MI] = 0.;
      hmm->t[k][p7H_MD] = tmp[1];
      hmm->t[k][p7H_IM] = 1.;                                          /* I transitions irrelevant since I's are unreached. */
      hmm->t[k][p7H_II] = 0.;
      if (k > 0) esl_dirichlet_FSampleUniform(r, 2,      hmm->t[k]+5); /* delete transitions to M,D */
    }

  /* Node M is special: no transitions to D, transitions to M
   * are interpreted as transitions to E. Overwrite a little of
   * what we did in node M.
   */
  hmm->t[M][p7H_MM] = 1.;
  hmm->t[M][p7H_MD] = 0.;	
  hmm->t[M][p7H_DM] = 1.;
  hmm->t[M][p7H_DD] = 0.;
  
  /* Add mandatory annotation
   */
  p7_hmm_SetName(hmm, "sampled-hmm");
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  hmm->nseq     = 0;
  hmm->eff_nseq = 0;
  p7_hmm_SetCtime(hmm);
  hmm->checksum = 0;

  /* SRE DEBUGGING */
  p7_hmm_Validate(hmm, NULL, 0.0001);

  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}

/* Function:  p7_hmm_SampleUniform()
 * Incept:    SRE, Thu Feb 22 10:04:19 2007 [Janelia]
 *
 * Purpose:   Sample a model that uses uniform transition probabilities,
 *            determined by <tmi>, <tii>, <tmd>, and <tdd>,
 *            the probabilistic equivalent of gap-open/gap-extend for
 *            inserts, deletes.
 *            
 *            Useful for testing expected behavior on single-sequence
 *            models, where transitions are position-independent.
 *
 * Returns:   <eslOK> on success, and the new hmm is returned
 *            through <ret_hmm); caller is responsible for 
 *            freeing this object with <p7_hmm_Destroy()>.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      J1/5.
 */
int
p7_hmm_SampleUniform(ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, 
		     float tmi, float tii, float tmd, float tdd,
		     P7_HMM **ret_hmm)
{
  int     status;
  P7_HMM *hmm    = NULL;
  char   *logmsg = "[HMM with uniform transitions, random emissions]";
  int     k;

  hmm = p7_hmm_Create(M, abc);
  if (hmm == NULL) { status = eslEMEM; goto ERROR; }
  
  for (k = 0; k <= M; k++)
    {
      if (k > 0) esl_dirichlet_FSampleUniform(r, abc->K, hmm->mat[k]);
      esl_dirichlet_FSampleUniform(r, abc->K, hmm->ins[k]);
      hmm->t[k][p7H_MM] = 1.0 - tmi - tmd;
      hmm->t[k][p7H_MI] = tmi;
      hmm->t[k][p7H_MD] = tmd;
      hmm->t[k][p7H_IM] = 1.0 - tii;
      hmm->t[k][p7H_II] = tii;
      hmm->t[k][p7H_DM] = 1.0 - tdd;
      hmm->t[k][p7H_DD] = tdd;
    }

  /* Deal w/ special stuff at node 0, M, overwriting some of what we
   * just did. 
   */
  hmm->t[M][p7H_MM] = 1.0 - tmi;
  hmm->t[M][p7H_MD] = 0.;
  hmm->t[M][p7H_DM] = 1.0;
  hmm->t[M][p7H_DD] = 0.;
  
  /* Add mandatory annotation
   */
  p7_hmm_SetName(hmm, "sampled-hmm");
  p7_hmm_AppendComlog(hmm, 1, &logmsg);
  hmm->nseq     = 0;
  hmm->eff_nseq = 0;
  p7_hmm_SetCtime(hmm);
  hmm->checksum = 0;

  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}



/* Function:  p7_hmm_Compare()
 * Incept:    SRE, Sat Jan  6 14:14:58 2007 [Casa de Gatos]
 *
 * Purpose:   Compare two HMMs <h1> and <h2> to each other;
 *            return <eslOK> if they're identical, and <eslFAIL>
 *            if they differ. Floating-point probabilities are 
 *            compared for equality within a fractional tolerance
 *            <tol>. 
 */
int
p7_hmm_Compare(P7_HMM *h1, P7_HMM *h2, float tol)
{
  int k, z;
  
  if (h1->abc->type != h2->abc->type) return eslFAIL;
  if (h1->M         != h2->M)         return eslFAIL;
  if (h1->flags     != h2->flags)     return eslFAIL;
  
  for (k = 0; k <= h1->M; k++)	/* (it's safe to include 0 here.) */
    {
      if (esl_vec_FCompare(h1->mat[k], h2->mat[k], h1->abc->K, tol) != eslOK) return eslFAIL;
      if (esl_vec_FCompare(h1->ins[k], h2->ins[k], h1->abc->K, tol) != eslOK) return eslFAIL;
      if (esl_vec_FCompare(h1->t[k],   h2->t[k],   7,          tol) != eslOK) return eslFAIL;
    }

  if (strcmp(h1->name,   h2->name)   != 0) return eslFAIL;
  if (strcmp(h1->comlog, h2->comlog) != 0) return eslFAIL;
  if (strcmp(h1->ctime,  h2->ctime)  != 0) return eslFAIL;
  if (h1->nseq     != h2->nseq)            return eslFAIL;
  if (h1->eff_nseq != h2->eff_nseq)        return eslFAIL;
  if (h1->checksum != h2->checksum)        return eslFAIL;

  if ((h1->flags & p7H_ACC)  && strcmp(h1->acc,  h2->acc)  != 0) return eslFAIL;
  if ((h1->flags & p7H_DESC) && strcmp(h1->desc, h2->desc) != 0) return eslFAIL;
  if ((h1->flags & p7H_RF)   && strcmp(h1->rf,   h2->rf)   != 0) return eslFAIL;
  if ((h1->flags & p7H_CS)   && strcmp(h1->cs,   h2->cs)   != 0) return eslFAIL;
  if ((h1->flags & p7H_CA)   && strcmp(h1->ca,   h2->ca)   != 0) return eslFAIL;
  if ((h1->flags & p7H_MAP)  && esl_vec_ICompare(h1->map, h2->map, h1->M+1) != 0) return eslFAIL;

  if (h1->flags & p7H_GA) {
    if (esl_FCompare(h1->cutoff[p7_GA1], h2->cutoff[p7_GA1], tol) != eslOK) return eslFAIL;
    if (esl_FCompare(h1->cutoff[p7_GA2], h2->cutoff[p7_GA2], tol) != eslOK) return eslFAIL;
  }
  if (h1->flags & p7H_TC) {
    if (esl_FCompare(h1->cutoff[p7_TC1], h2->cutoff[p7_TC1], tol) != eslOK) return eslFAIL;
    if (esl_FCompare(h1->cutoff[p7_TC2], h2->cutoff[p7_TC2], tol) != eslOK) return eslFAIL;
  }
  if (h1->flags & p7H_NC) {
    if (esl_FCompare(h1->cutoff[p7_NC1], h2->cutoff[p7_NC1], tol) != eslOK) return eslFAIL;
    if (esl_FCompare(h1->cutoff[p7_NC2], h2->cutoff[p7_NC2], tol) != eslOK) return eslFAIL;
  }

  if (h1->flags & p7H_STATS) {
    for (z = 0; z < p7_NEVPARAM; z++)
      if (esl_FCompare(h1->evparam[z], h2->evparam[z], tol) != eslOK) return eslFAIL;
  }

  return eslOK;
}

/* Function:  p7_hmm_Validate()
 * Incept:    SRE, Sat Jan  6 14:43:00 2007 [Casa de Gatos]
 *
 * Purpose:   Validates the internals of the HMM structure <hmm>.
 * 
 *            Probability vectors are validated to sum up to
 *            within a fractional tolerance <tol> of 1.0.
 *
 *            Probably only useful for debugging and development,
 *            not production code.
 *
 * Returns:   <eslOK> if <hmm> internals look fine.
 *            Returns <eslFAIL> if something is wrong.
 */
int
p7_hmm_Validate(P7_HMM *hmm, char *errbuf, float tol)
{
  int status;
  int k;

  if (hmm            == NULL)       ESL_XFAIL(eslFAIL, errbuf, "HMM is a null pointer");
  if (hmm->M         <  1)          ESL_XFAIL(eslFAIL, errbuf, "HMM has M < 1");
  if (hmm->abc       == NULL)       ESL_XFAIL(eslFAIL, errbuf, "HMM has no alphabet reference");
  if (hmm->abc->type == eslUNKNOWN) ESL_XFAIL(eslFAIL, errbuf, "HMM's alphabet is set to unknown");
  
  for (k = 0; k <= hmm->M; k++)
    {
      if (esl_vec_FValidate(hmm->mat[k], hmm->abc->K, tol, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "mat[%d] fails pvector validation", k);
      if (esl_vec_FValidate(hmm->ins[k], hmm->abc->K, tol, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "ins[%d] fails pvector validation", k);
      if (esl_vec_FValidate(hmm->t[k],   3,           tol, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "t_M[%d] fails pvector validation", k);
      if (esl_vec_FValidate(hmm->t[k]+3, 2,           tol, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "t_I[%d] fails pvector validation", k);
      if (esl_vec_FValidate(hmm->t[k]+5, 2,           tol, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "t_D[%d] fails pvector validation", k);
    }
  if (hmm->t[hmm->M][p7H_MD] != 0.0) ESL_XFAIL(eslFAIL, errbuf, "TMD should be 0 for last node");
  if (hmm->t[hmm->M][p7H_DM] != 1.0) ESL_XFAIL(eslFAIL, errbuf, "TDM should be 1 for last node");
  if (hmm->t[hmm->M][p7H_DD] != 0.0) ESL_XFAIL(eslFAIL, errbuf, "TDD should be 0 for last node");

  /* Don't be strict about mandatory name, comlog, ctime for now in development */
  /*  if (hmm->name     == NULL) return eslFAIL; */
  /*  if (hmm->comlog   == NULL) return eslFAIL; */
  /*  if (hmm->ctime    == NULL) return eslFAIL;  */
  if (hmm->nseq     <  0 )   ESL_XFAIL(eslFAIL, errbuf, "invalid nseq");
  if (hmm->eff_nseq <  0 )   ESL_XFAIL(eslFAIL, errbuf, "invalid eff_nseq");
  if (hmm->checksum <  0 )   ESL_XFAIL(eslFAIL, errbuf, "invalid checksum");

  if (  (hmm->flags & p7H_ACC)  && hmm->acc  == NULL) ESL_XFAIL(eslFAIL, errbuf, "accession null but p7H_ACC flag is up");
  if (! (hmm->flags & p7H_ACC)  && hmm->acc  != NULL) ESL_XFAIL(eslFAIL, errbuf, "accession present but p7H_ACC flag is down");
  if (  (hmm->flags & p7H_DESC) && hmm->desc == NULL) ESL_XFAIL(eslFAIL, errbuf, "description null but p7H_DESC flag is up");
  if (! (hmm->flags & p7H_DESC) && hmm->desc != NULL) ESL_XFAIL(eslFAIL, errbuf, "description present but p7H_DESC flag is down");
  if (hmm->flags & p7H_RF) {
    if (hmm->rf == NULL || strlen(hmm->rf) != hmm->M+1) ESL_XFAIL(eslFAIL, errbuf, "p7H_RF flag up, but rf string is invalid");
  } else 
    if (hmm->rf != NULL)                                ESL_XFAIL(eslFAIL, errbuf, "p7H_RF flag down, but rf string is present");
  if (hmm->flags & p7H_CS) {
    if (hmm->cs == NULL || strlen(hmm->cs) != hmm->M+1) ESL_XFAIL(eslFAIL, errbuf, "p7H_CS flag up, but cs string is invalid");
  } else 
    if (hmm->cs != NULL)                                ESL_XFAIL(eslFAIL, errbuf, "p7H_CS flag down, but cs string is present");
  if (hmm->flags & p7H_CA) {
    if (hmm->ca == NULL || strlen(hmm->ca) != hmm->M+1) ESL_XFAIL(eslFAIL, errbuf, "p7H_CA flag up, but ca string is invalid");
  } else 
    if (hmm->ca != NULL)                                ESL_XFAIL(eslFAIL, errbuf, "p7H_CA flag down, but ca string is present");
  if (  (hmm->flags & p7H_MAP) && hmm->map == NULL)      ESL_XFAIL(eslFAIL, errbuf, "p7H_MAP flag up, but map string is null");
  if (! (hmm->flags & p7H_MAP) && hmm->map != NULL)      ESL_XFAIL(eslFAIL, errbuf, "p7H_MAP flag down, but map string is present");

  if (hmm->flags & p7H_STATS) {
    if (hmm->evparam[p7_LAMBDA] <= 0.) ESL_XFAIL(eslFAIL, errbuf, "lambda parameter can't be negative");
  }

  return eslOK;

 ERROR:
  return status;
}
/*------------- end of debugging/development code ----------------*/




/*****************************************************************
 * 5. Other routines in the API.
 *****************************************************************/

/* Function:  p7_hmm_CalculateOccupancy()
 * Incept:    SRE, Mon Jan 22 08:10:05 2007 [Janelia]
 *
 * Purpose:   Calculate a vector <occ[1..M]> containing probability
 *            that each match state is used in a sampled path through
 *            the model. Caller provides allocated space (<M+1> floats)
 *            for <occ>.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_hmm_CalculateOccupancy(const P7_HMM *hmm, float *occ)
{
  int k;

  occ[0] = 0.;			                    /* no M_0 state */
  occ[1] = hmm->t[0][p7H_MI] + hmm->t[0][p7H_MM];   /* initialize w/ 1 - B->D_1 */
  for (k = 2; k <= hmm->M; k++)
    occ[k] = occ[k-1] * (hmm->t[k-1][p7H_MM] + hmm->t[k-1][p7H_MI]) +
      (1.0-occ[k-1]) * hmm->t[k-1][p7H_DM];
  return eslOK;
}

/*---------------- end of the rest of the API -------------------*/




/*****************************************************************
 * 6. Unit tests.
 *****************************************************************/
#ifdef p7HMM_TESTDRIVE

/* The occupancy test is based on the principle that
 * the stationary match occupancy probability in a random HMM 
 * converges to 0.6, for long enough M (STL11/138)
 */
static void
utest_occupancy(P7_HMM *hmm)
{
  char  *msg = "modelconfig.c::calculate_occupancy() unit test failed";
  float *occ;
  float  x;

  occ = malloc(sizeof(float) * (hmm->M+1));
  p7_hmm_CalculateOccupancy(hmm, occ);
  x = esl_vec_FSum(occ+1, hmm->M) / (float) hmm->M;
  if (esl_FCompare(x, 0.6, 0.1) != eslOK)           esl_fatal(msg);
  free(occ);
  return;
}

#endif /*p7HMM_TESTDRIVE*/
/*---------------------- end of unit tests -----------------------*/


/*****************************************************************
 * 7. Test driver.
 *****************************************************************/

#ifdef p7HMM_TESTDRIVE

#include "p7_config.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "hmmer.h"

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r    = esl_randomness_Create(42);
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm  = NULL;
  int             M    = 100;

  p7_hmm_Sample(r, M, abc, &hmm);

  utest_occupancy(hmm);

  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  exit(0); /* success */
}

#endif /*p7HMM_TESTDRIVE*/
/*-------------------- end of test driver ---------------------*/


/************************************************************
 * @LICENSE@
 ************************************************************/

