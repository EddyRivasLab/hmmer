/* The Plan7 core HMM data structure.
 * 
 * Contents:
 *   1. The P7_HMM object: allocation, initialization, destruction.
 *   2. Convenience routines for setting fields in an HMM.
 *   3. Renormalization and rescaling counts in core HMMs.
 *   4. Unit tests.
 *   5. Test driver. 
 *   6. Copyright and license.
 * 
 * SRE, Mon Jan  1 16:20:29 2007 [Casa de Gatos] [Verdi, La Traviata; Maria Callas]
 * SVN $Id$
 */

#include "p7_config.h"		/* must be included first */

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_vectorops.h>


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
p7_hmm_Create(int M, ESL_ALPHABET *abc) 
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
  int     status;

  ESL_ALLOC(hmm, sizeof(P7_HMM));
  hmm->abc      = NULL;
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
  hmm->ctime    = NULL;
  hmm->map      = NULL;
  hmm->checksum = 0;

  hmm->ga1 = hmm->ga2 = 0.;
  hmm->tc1 = hmm->tc2 = 0.;
  hmm->nc1 = hmm->nc2 = 0.;

  hmm->flags    = 0;
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
 * Throws:    <eslEMEM> on allocation failure; in this case, the entire
 *            HMM is free'd (including the shell).
 */
int
p7_hmm_CreateBody(P7_HMM *hmm, int M, ESL_ALPHABET *abc) 
{
  int k, x;
  int status;

  hmm->abc    = abc;
  hmm->M      = M;

  ESL_ALLOC(hmm->t,      M              * sizeof(float *));
  ESL_ALLOC(hmm->mat,    (M+1)          * sizeof(float *));
  ESL_ALLOC(hmm->ins,    M              * sizeof(float *));
  ESL_ALLOC(hmm->t[0],   (7*M)          * sizeof(float));
  ESL_ALLOC(hmm->mat[0], (abc->K*(M+1)) * sizeof(float));
  ESL_ALLOC(hmm->ins[0], (abc->K*M)     * sizeof(float));

  for (k = 1; k < M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * hmm->abc->K;
    hmm->ins[k] = hmm->ins[0] + k * hmm->abc->K;
    hmm->t[k]   = hmm->t[0]   + k * 7;
  }
  hmm->mat[M] = hmm->mat[0] + M * hmm->abc->K;

  if (hmm->flags & p7_RF)  ESL_ALLOC(hmm->rf,  (M+2) * sizeof(char));
  if (hmm->flags & p7_CS)  ESL_ALLOC(hmm->cs,  (M+2) * sizeof(char));
  if (hmm->flags & p7_CA)  ESL_ALLOC(hmm->ca,  (M+2) * sizeof(char));
  if (hmm->flags & p7_MAP) ESL_ALLOC(hmm->map, (M+1) * sizeof(int));
  
  if ((status = p7_hmm_Zero(hmm)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  return status;
}  


/* Function:  p7_hmm_Destroy()
 * Incept:    SRE, Fri Mar 31 15:13:25 2006 [St. Louis]
 *
 * Purpose:   Frees both the shell and body of an <hmm>.
 *            Works even if the <hmm> is damaged (incompletely allocated)
 *            or even <NULL>.
 *
 * Note:      Remember, leave abc alone. It's just a reference ptr
 *            that the application gave us.
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

/* Function:  p7_hmm_Zero()
 * Incept:    SRE, Mon Jan  1 16:32:59 2007 [Casa de Gatos]
 *
 * Purpose:   Zeroes the counts/probabilities fields in a model;
 *            drop the p7_HASPROBS flag (probs no longer valid).
 *
 * Returns:   <eslOK> on success.
 */
int
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

  hmm->flags  &= ~p7_HASPROB;	/* invalidates probabilities        */
  hmm->flags  &= ~p7_HASBITS;	/* invalidates scores               */
  return eslOK;
}

/* Function:  p7_hmm_Dump()
 * Incept:    SRE, Mon Jan  1 18:44:15 2007 [Casa de Gatos]
 *
 * Purpose:   Debugging: dump the probabilities (or counts) from an HMM.
 * 
 * Returns:   <eslOK> on success.
 */
int
p7_hmm_Dump(FILE *fp, P7_HMM *hmm)
{
  int k;			/* counter for nodes */
  int x;			/* counter for symbols */
  int ts;			/* counter for state transitions */
  
  fprintf(fp, "B->M1,B->D1: %5.1f %5.1f\n", hmm->t[0][p7_TMM], hmm->t[0][p7_TMD]);
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
	fprintf(fp, "%5.1f ", (k < hmm->M) ? hmm->ins[k][x] : 0.);
      fputs("\n", fp);
				/* Line 3: transition probs */
      fprintf(fp, "       ");
      for (ts = 0; ts < 7; ts++)
	fprintf(fp, "%5.1f ", (k < hmm->M) ? hmm->t[k][ts] : 0.); 
      fputs("\n", fp);
    }
  fputs("//\n", fp);
  return eslOK;
}


/* Function:  p7_hmm_DescribeStatetype()
 * Incept:    SRE, Mon Jan  1 18:47:34 2007 [Casa de Gatos]
 *
 * Purpose:   Returns the state type in text, as a string of length 1 
 *            (2 if you count NUL). For example, <p7_Statetype(p7_STS)>
 *            returns "S".
 */
char *
p7_hmm_DescribeStatetype(char st)
{
  switch (st) {
  case p7_STM: return "M";
  case p7_STD: return "D";
  case p7_STI: return "I";
  case p7_STS: return "S";
  case p7_STN: return "N";
  case p7_STB: return "B";
  case p7_STE: return "E";
  case p7_STC: return "C";
  case p7_STT: return "T";
  case p7_STJ: return "J";
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
    hmm->flags &= ~p7_ACC;
  } else {
    n = strlen(acc);
    ESL_RALLOC(hmm->acc, tmp, sizeof(char)*(n+1));
    strcpy(hmm->acc, acc);
    if ((status = esl_strchop(hmm->acc, n)) != eslOK) goto ERROR;
    hmm->flags |= p7_ACC;
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
      hmm->flags &= ~p7_DESC;
    }
  else
    {
      n = strlen(desc);
      ESL_RALLOC(hmm->desc, tmp, sizeof(char)*(n+1));
      strcpy(hmm->desc, desc);
      if ((status = esl_strchop(hmm->desc, n)) != eslOK) goto ERROR;
      hmm->flags |= p7_DESC;
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


/*****************************************************************
 * 3. Renormalization and rescaling counts in core HMMs.
 *****************************************************************/ 

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

  for (k = 1; k <= hmm->M; k++)  esl_vec_FScale(hmm->mat[k], hmm->abc->K, scale);
  for (k = 1; k <  hmm->M; k++)  esl_vec_FScale(hmm->ins[k], hmm->abc->K, scale);
  for (k = 0; k <  hmm->M; k++)  esl_vec_FScale(hmm->t[k],   7,           scale);
  return eslOK;
}

/* Function: p7_hmm_Renormalize()
 * Incept:   SRE, Mon Jan  1 18:39:42 2007 [Casa de Gatos]
 * 
 * Purpose:  Take an HMM in counts form, and renormalize
 *           all probability vectors in the core probability model. Enforces
 *           Plan7 restrictions on nonexistent transitions. Raises the
 *           <p7_HASPROB> flag. 
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

  for (k = 1; k <= hmm->M; k++) 
    esl_vec_FNorm(hmm->mat[k], hmm->abc->K);
  for (k = 1; k < hmm->M; k++)
    esl_vec_FNorm(hmm->ins[k], hmm->abc->K);
  for (k = 1; k < hmm->M; k++) {
      esl_vec_FNorm(hmm->t[k],   3);	/* match  */
      esl_vec_FNorm(hmm->t[k]+3, 2);	/* insert */
      esl_vec_FNorm(hmm->t[k]+5, 2);	/* delete */
  }

  /* The t[0] is special: TM* are the begin transitions
   */
  hmm->t[0][p7_TMI] = 0.;       /* make sure TMI = 0. */
  esl_vec_FNorm(hmm->t[0], 3);  /* begin transitions; TMM, TMD are valid */
  esl_vec_FSet(hmm->t[0]+3, 2, 0.); /* t[0] delete */
  esl_vec_FSet(hmm->t[0]+5, 2, 0.); /* t[0] insert */
  esl_vec_FSet(hmm->mat[0], hmm->abc->K, 0.); /* mat[0] */
  esl_vec_FSet(hmm->ins[0], hmm->abc->K, 0.); /* ins[0] */

  hmm->flags |= p7_HASPROB;	/* set the probabilities OK flag */
  return eslOK;
}
  


/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7HMM_TESTDRIVE

static void
utest_foo(void)
{

  return;
}


#endif /*p7HMM_TESTDRIVE*/
/*---------------------- end of unit tests -----------------------*/




/*****************************************************************
 * 5. Test driver.
 *****************************************************************/

#ifdef p7HMM_TESTDRIVE

#include <p7_config.h>
#include <p7_hmm.h>

int
main(int argc, char **argv)
{
  
  exit(0); /* success */
}

#endif /*p7HMM_TESTDRIVE*/
/*-------------------- end of test driver ---------------------*/


/************************************************************
 * @LICENSE@
 ************************************************************/

