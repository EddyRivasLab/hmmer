/* Routines for the P7_PROFILE structure - Plan 7's search profile
 *                                         
 *    1. The P7_PROFILE object: allocation, initialization, destruction.
 *    2. Access methods.
 *    3. Debugging and development tools.
 *    4. Unit tests.
 *    5. Test driver.
 *    6. Example.
 *
 * See also: 
 *   modelconfig.c : routines that configure a profile given an HMM
 */
#include "p7_config.h"

#include <stdio.h>
#include <string.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "base/p7_profile.h"
#include "base/p7_trace.h"

/*****************************************************************
 * 1. The P7_PROFILE object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  p7_profile_Create()
 * Synopsis:  Allocates a profile.
 *
 * Purpose:   Allocates for a profile of up to <M> nodes, for digital
 *            alphabet <abc>.
 *            
 *            Because this function might be in the critical path (in
 *            hmmscan, for example), we leave much of the model
 *            uninitialized, including scores and length model
 *            probabilities. The <p7_profile_Config()> call is what
 *            sets these.
 *            
 *            The reference pointer <gm->abc> is set to <abc>.
 *
 * Returns:   a pointer to the newly allocated profile.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_PROFILE *
p7_profile_Create(int allocM, const ESL_ALPHABET *abc)
{
  P7_PROFILE *gm = NULL;
  int         x;
  int         status;

  /* level 0 */
  ESL_ALLOC(gm, sizeof(P7_PROFILE));
  gm->tsc       = NULL;
  gm->rsc       = NULL;
  gm->name      = NULL;
  gm->acc       = NULL;
  gm->desc      = NULL;
  gm->rf        = NULL;
  gm->mm        = NULL;
  gm->cs        = NULL;
  gm->consensus = NULL;

  /* level 1 */
  ESL_ALLOC(gm->tsc,       sizeof(float)   * (allocM+1) * p7P_NTRANS); /* 0..M */
  ESL_ALLOC(gm->rsc,       sizeof(float *) * abc->Kp);
  ESL_ALLOC(gm->rf,        sizeof(char)    * (allocM+2)); /* yes, +2: each is (0)1..M, +trailing \0  */
  ESL_ALLOC(gm->mm,        sizeof(char)    * (allocM+2));
  ESL_ALLOC(gm->cs,        sizeof(char)    * (allocM+2));
  ESL_ALLOC(gm->consensus, sizeof(char)    * (allocM+2));
  gm->rsc[0] = NULL;
  
  /* level 2 */
  ESL_ALLOC(gm->rsc[0], sizeof(float) * abc->Kp * (allocM+1) * p7P_NR);
  for (x = 1; x < abc->Kp; x++) 
    gm->rsc[x] = gm->rsc[0] + x * (allocM+1) * p7P_NR;

  /* Initialization of tsc[0], including removal of I0.  tsc[k-1,LM],tsc[k-1,GM] will be configured + overwritten later */
  esl_vec_FSet(gm->tsc, p7P_NTRANS, -eslINFINITY);  
  /* tsc[M] initialized and Im removed when we know actual M : see modelconfig.c */

  for (x = 0; x < abc->Kp; x++) {        
    P7P_MSC(gm, 0, x) = -eslINFINITY;                  /* no emissions from nonexistent M_0... */
    P7P_ISC(gm, 0, x) = -eslINFINITY;                  /* nor I_0...                           */
    /* I_M is initialized in profile config, when we know actual M, not just allocated max M   */
  }
  x = esl_abc_XGetGap(abc);	                       /* no emission can emit/score gap characters */
  esl_vec_FSet(gm->rsc[x], (allocM+1)*p7P_NR, -eslINFINITY);
  x = esl_abc_XGetMissing(abc);	                       /* no emission can emit/score missing data characters */
  esl_vec_FSet(gm->rsc[x], (allocM+1)*p7P_NR, -eslINFINITY);


  /* Set remaining info  */
  gm->M                = 0;
  gm->allocM           = allocM;
  gm->L                = -1;  	   /* "unset" flag */
  gm->nj               = -1.0f;    /* "unset" flag */
  gm->pglocal          = -1.0f;    /* "unset" flag */

  gm->roff             = -1;
  gm->eoff             = -1;
  gm->offs[p7_MOFFSET] = -1;
  gm->offs[p7_FOFFSET] = -1;
  gm->offs[p7_POFFSET] = -1;

  gm->name             = NULL;
  gm->acc              = NULL;
  gm->desc             = NULL;
  gm->rf[0]            = 0;     /* RF line is optional annotation; this flags that it's not set yet */
  gm->mm[0]            = 0;     /* likewise for MM annotation line */
  gm->cs[0]            = 0;     /* likewise for CS annotation line */
  gm->consensus[0]     = 0;
  
  for (x = 0; x < p7_NEVPARAM; x++) gm->evparam[x] = p7_EVPARAM_UNSET;
  for (x = 0; x < p7_NCUTOFFS; x++) gm->cutoff[x]  = p7_CUTOFF_UNSET;
  for (x = 0; x < p7_MAXABET;  x++) gm->compo[x]   = p7_COMPO_UNSET;

  gm->max_length  = -1;		/* "unset" */
  gm->abc         = abc;

  return gm;

 ERROR:
  p7_profile_Destroy(gm);
  return NULL;
}


/* Function:  p7_profile_Copy()
 * Synopsis:  Copy a profile.
 *
 * Purpose:   Copies profile <src> to profile <dst>, where <dst>
 *            has already been allocated to be of sufficient size,
 *            and has the same alphabet.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <dst> is too small 
 *            to fit <src> or is for a different alphabet.
 */
int
p7_profile_Copy(const P7_PROFILE *src, P7_PROFILE *dst)
{
  int x,z;
  int status;

  if (src->M         >   dst->allocM)   ESL_EXCEPTION(eslEINVAL, "destination profile is too small to hold a copy of source profile");
  if (src->abc->type != dst->abc->type) ESL_EXCEPTION(eslEINVAL, "destination profile has different alphabet than source");

  dst->M = src->M;
  esl_vec_FCopy(src->tsc, (src->M+1)*p7P_NTRANS, dst->tsc);
  for (x = 0; x < src->abc->Kp;   x++) esl_vec_FCopy(src->rsc[x], (src->M+1)*p7P_NR, dst->rsc[x]);
  for (x = 0; x < p7P_NXSTATES;   x++) esl_vec_FCopy(src->xsc[x], p7P_NXTRANS,       dst->xsc[x]);

  dst->L           = src->L;
  dst->nj          = src->nj;
  dst->pglocal     = src->pglocal;

  if (dst->name) free(dst->name);   
  if (dst->acc)  free(dst->acc);    
  if (dst->desc) free(dst->desc);   
  if ((status = esl_strdup(src->name, -1, &(dst->name)))      != eslOK) return status; 
  if ((status = esl_strdup(src->acc,  -1, &(dst->acc)))       != eslOK) return status; 
  if ((status = esl_strdup(src->desc, -1, &(dst->desc)))      != eslOK) return status; 

  strcpy(dst->rf,        src->rf);         /* RF is optional: if it's not set, *rf=0, and strcpy still works fine */
  strcpy(dst->mm,        src->mm);         /* MM is also optional annotation */
  strcpy(dst->cs,        src->cs);         /* CS is also optional annotation */
  strcpy(dst->consensus, src->consensus);  /* consensus though is always present on a valid profile */

  for (z = 0; z < p7_NEVPARAM; z++) dst->evparam[z] = src->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) dst->cutoff[z]  = src->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) dst->compo[z]   = src->compo[z];

  for (x = 0; x < p7_NOFFSETS; ++x) dst->offs[x] = src->offs[x];
  dst->roff        = src->roff;
  dst->eoff        = src->eoff;

  dst->max_length  = src->max_length;
  return eslOK;
}


/* Function:  p7_profile_Clone()
 * Synopsis:  Duplicates a profile.
 *
 * Purpose:   Duplicate profile <gm>; return a pointer
 *            to the newly allocated copy.
 *            
 * Returns:   ptr to new clone.
 * 
 * Throws:    <NULL> on allocation failure.
 */
P7_PROFILE *
p7_profile_Clone(const P7_PROFILE *gm)
{
  P7_PROFILE *g2 = NULL;
  int         status;

  if ((g2 = p7_profile_Create(gm->allocM, gm->abc)) == NULL) return NULL;
  if ((status = p7_profile_Copy(gm, g2)) != eslOK) goto ERROR;
  return g2;
  
 ERROR:
  p7_profile_Destroy(g2);
  return NULL;
}



/* Function:  p7_profile_Reuse()
 * Synopsis:  Prepare profile to be re-used for a new HMM.
 *
 * Purpose:   Prepare profile <gm>'s memory to be re-used
 *            for a new HMM.
 */
int
p7_profile_Reuse(P7_PROFILE *gm)
{
  /* name, acc, desc annotation is dynamically allocated for each HMM */
  if (gm->name != NULL) { free(gm->name); gm->name = NULL; }
  if (gm->acc  != NULL) { free(gm->acc);  gm->acc  = NULL; }
  if (gm->desc != NULL) { free(gm->desc); gm->desc = NULL; }

  /* set annotations to empty strings */
  gm->rf[0]        = 0;
  gm->mm[0]        = 0;
  gm->cs[0]        = 0;
  gm->consensus[0] = 0;
      
  /* reset some other things, but leave the rest alone. */
  gm->M       = 0;
  gm->L       = -1;
  gm->nj      = -1.0f;
  gm->pglocal = -1.0f;

  gm->offs[p7_MOFFSET] = -1;
  gm->offs[p7_FOFFSET] = -1;
  gm->offs[p7_POFFSET] = -1;
  gm->roff             = -1;
  gm->eoff             = -1;

  gm->max_length       = -1;
  return eslOK;
}


/* Function:  p7_profile_Sizeof()
 * Synopsis:  Return the allocated size of a P7_PROFILE.
 *
 * Purpose:   Return the allocated size of a <P7_PROFILE>, in bytes.
 * 
 *            Because we don't know the allocated size of some
 *            annotation fields (name, acc, desc), we actually
 *            calculate the minimum allocation size, counting these
 *            strings as strlen(s)+1. This suffices because Sizeof()'s
 *            most critical use is in determining the minimum
 *            necessary size of a serialized structure. When we
 *            instead use Sizeof() to report memory consumption, a
 *            rough number is adequate.
 *            
 *            Very roughly, 276*M bytes; 50KB per typical M=200 model;
 *            30MB for a design limit M=100K model. Dominated by
 *            transition/emission parameters.
 */
size_t
p7_profile_Sizeof(P7_PROFILE *gm)
{
  size_t n = 0;

  /* these mirror malloc()'s in p7_profile_Create(); maintain one:one correspondence for maintainability */
  n += sizeof(P7_PROFILE);
  n += sizeof(float)   * (gm->allocM+1) * p7P_NTRANS;           /* gm->tsc       */
  n += sizeof(float *) * gm->abc->Kp;	                        /* gm->rsc       */
  n += sizeof(float)   * gm->abc->Kp * (gm->allocM+1) * p7P_NR; /* gm->rsc[0]    */
  n += sizeof(char)    * (gm->allocM+2);	                /* gm->rf        */
  n += sizeof(char)    * (gm->allocM+2);                        /* gm->mm        */
  n += sizeof(char)    * (gm->allocM+2);	                /* gm->cs        */
  n += sizeof(char)    * (gm->allocM+2);	                /* gm->consensus */

  if (gm->name) n += sizeof(char) * (strlen(gm->name) + 1);
  if (gm->acc)  n += sizeof(char) * (strlen(gm->acc)  + 1);
  if (gm->desc) n += sizeof(char) * (strlen(gm->desc) + 1);

  return n;
}


/* Function:  p7_profile_Destroy()
 * Synopsis:  Frees a profile.
 *
 * Purpose:   Frees a profile <gm>.
 *
 * Returns:   (void).
 */
void
p7_profile_Destroy(P7_PROFILE *gm)
{
  if (gm) {
    if (gm->rsc && gm->rsc[0]) free(gm->rsc[0]);
    if (gm->tsc)       free(gm->tsc);
    if (gm->rsc)       free(gm->rsc);
    if (gm->name)      free(gm->name);
    if (gm->acc)       free(gm->acc);
    if (gm->desc)      free(gm->desc);
    if (gm->rf)        free(gm->rf);
    if (gm->mm)        free(gm->mm);
    if (gm->cs)        free(gm->cs);
    if (gm->consensus) free(gm->consensus);
    free(gm);
  }
  return;
}


/*****************************************************************
 * 2. Access methods.
 *****************************************************************/

/* Function:  p7_profile_IsLocal()
 * Synopsis:  Return TRUE if profile is in a local alignment mode
 *
 * Purpose:   Return <TRUE> if profile is in a pure local alignment mode (only),
 *            not dual-mode or glocal-only mode.
 */
int
p7_profile_IsLocal(const P7_PROFILE *gm)
{
  return (gm->pglocal == 0.0f ? TRUE : FALSE);
}

int
p7_profile_IsGlocal(const P7_PROFILE *gm)
{
  return (gm->pglocal == 1.0f ? TRUE : FALSE);
}

/* Function:  p7_profile_IsMultihit()
 * Synopsis:  Return TRUE if profile is in a multihit alignment mode.
 *
 * Purpose:   Return <TRUE> if profile is in a multihit alignment mode.
 */
int
p7_profile_IsMultihit(const P7_PROFILE *gm)
{
  if (gm->nj > 0.0) return TRUE;
  return FALSE;
}





/* Function:  p7_profile_GetT()
 *
 * Purpose:   Convenience function that looks up and returns a
 *            transition score in profile <gm> for a transition from
 *            state type <st1> in node <k1> to state type <st2> in
 *            node <k2>. For unique state types that aren't in nodes
 *            (<p7T_S>, for example), the <k> value is ignored, and
 *            it is customarily passed as 0.
 *
 *            Assumes that the transition is valid. If caller has any
 *            question about this, it needs to validate the transition
 *            itself. 
 *            
 *            Shouldn't be used in time-critical code.
 */
float
p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1, char st2, int k2)
{
  switch (st1) {
  case p7T_S: return 0.0;
  case p7T_N: return (st2 == p7T_B  ? gm->xsc[p7P_N][p7P_MOVE] : gm->xsc[p7P_N][p7P_LOOP]);
  case p7T_B: return (st2 == p7T_L  ? gm->xsc[p7P_B][0]        : gm->xsc[p7P_B][1]);
  case p7T_L: return P7P_TSC(gm, k2-1, p7P_LM);
  case p7T_G: return (st2 == p7T_MG ? gm->xsc[p7P_G][0]        : gm->xsc[p7P_G][1]);
  case p7T_E: return (st2 == p7T_C  ? gm->xsc[p7P_E][p7P_MOVE] : gm->xsc[p7P_E][p7P_LOOP]); 
  case p7T_C: return (st2 == p7T_T  ? gm->xsc[p7P_C][p7P_MOVE] : gm->xsc[p7P_C][p7P_LOOP]); 
  case p7T_J: return (st2 == p7T_B  ? gm->xsc[p7P_J][p7P_MOVE] : gm->xsc[p7P_J][p7P_LOOP]); 

  case p7T_ML: 
  case p7T_MG:
    switch (st2) {
    case p7T_ML: case p7T_MG: return P7P_TSC(gm, k1, p7P_MM);
    case p7T_IL: case p7T_IG: return P7P_TSC(gm, k1, p7P_MI);
    case p7T_DL: case p7T_DG: return P7P_TSC(gm, k1, p7P_MD);
    case p7T_E:               return 0.0;
    }

  case p7T_DL:
  case p7T_DG:
    switch (st2) {
    case p7T_ML: case p7T_MG:  return P7P_TSC(gm, k1, p7P_DM);
    case p7T_DL: case p7T_DG:  return P7P_TSC(gm, k1, p7P_DD);
    case p7T_E:                return 0.0;
    }

  case p7T_IL: 
  case p7T_IG:  
    switch (st2) {
    case p7T_ML: case p7T_MG: return P7P_TSC(gm, k1, p7P_IM);
    case p7T_IL: case p7T_IG: return P7P_TSC(gm, k1, p7P_II);
    }
  }    
  /*NOTREACHED*/
  return 0.0;
}


/*****************************************************************
 * 3. Debugging and development tools.
 *****************************************************************/

/* Function:  p7_profile_Dump()
 * Synopsis:  Dump a P7_PROFILE structure to stream, for examination.
 *
 * Purpose:   Write the contents of <gm> in terse
 *            debugging/examination tabular format, to the open output
 *            stream <fp>.
 */
int
p7_profile_Dump(FILE *fp, P7_PROFILE *gm)
{
  int width     = 9;
  int precision = 4;
  int k;
  int x;

  fputs("# profile object dump\n", fp);
  fprintf(fp, "   name      = %s\n", gm->name);
  fprintf(fp, "   length M  = %d\n", gm->M);
  fprintf(fp, "   alloc M   = %d\n", gm->allocM);
  fprintf(fp, "   accession = %s\n", gm->acc  ? gm->acc  : "[null accession]");
  fprintf(fp, "   desc      = %s\n", gm->desc ? gm->desc : "[null description]");
  fprintf(fp, "   alphabet  = %s\n", esl_abc_DecodeType(gm->abc->type));
  fputs("\n", fp);

  fputs("# configuration parameters\n", fp);
  fprintf(fp, "   length model L  = %d\n",   gm->L);
  fprintf(fp, "   multihit     nj = %.1f\n", gm->nj);
  fprintf(fp, "   glocal       pg = %.2f\n", gm->pglocal);
  fputs("\n", fp);

  fputs("# special profile states\n", fp);
  fprintf(fp, "state  %*s  %*s\n", width,  "LOOP/0", width,  "MOVE/1");
  fprintf(fp, "-----  %*s  %*s\n", width, "-------", width, "-------");
  fprintf(fp, "E      %*.*f  %*.*f\n",                 width, precision, gm->xsc[p7P_E][0], width, precision, gm->xsc[p7P_E][1]);
  fprintf(fp, "N      %*.*f  %*.*f\n",                 width, precision, gm->xsc[p7P_N][0], width, precision, gm->xsc[p7P_N][1]);
  fprintf(fp, "J      %*.*f  %*.*f\n",                 width, precision, gm->xsc[p7P_J][0], width, precision, gm->xsc[p7P_J][1]);
  fprintf(fp, "C      %*.*f  %*.*f\n",                 width, precision, gm->xsc[p7P_C][0], width, precision, gm->xsc[p7P_C][1]);
  fprintf(fp, "B      %*.*f  %*.*f  // B->L  B->G\n",  width, precision, gm->xsc[p7P_B][0], width, precision, gm->xsc[p7P_B][1]);
  fprintf(fp, "G      %*.*f  %*.*f  // G->M1 G->D1\n", width, precision, gm->xsc[p7P_G][0], width, precision, gm->xsc[p7P_G][1]);
  fputs("\n", fp);

  fputs("# state transition parameters\n", fp);
  fputs("pos   cons ", fp);
  for (x = 0; x < p7P_NTRANS; x++) fprintf(fp, "%*s ", width, p7_profile_DecodeT(x));
  fputs("\n", fp);

  fputs("----- ---- ", fp);
  for (x = 0; x < p7P_NTRANS; x++) fprintf(fp, "%*s ", width, "-------");
  fputs("\n", fp);

  for (k = 0; k <= gm->M; k++) 	/* transitions are valid for k=1..M-1, but stored for 0..M */
    {
      fprintf(fp, "%5d    %c ", k, gm->consensus[k]);
      for (x = 0; x < p7P_NTRANS; x++) fprintf(fp, "%*.*f ", width, precision, P7P_TSC(gm, k, x));
      fputs("\n", fp);
    }
  fputs("\n", fp);

  fputs("# match state emission parameters\n", fp);
  fputs("pos   cons   rf   cs ", fp);
  for (x = 0; x < gm->abc->K; x++) fprintf(fp, "%*s%c ", width-1, " ", gm->abc->sym[x]);
  fputs("\n", fp);

  fputs("----- ---- ---- ---- ", fp);
  for (x = 0; x < gm->abc->K; x++) fprintf(fp, "%*s ", width, "-------");
  fputs("\n", fp);

  for (k = 0; k <= gm->M; k++) 
    {
      fprintf(fp, "%5d    %c    %c    %c ", k, gm->consensus[k], gm->rf[0] ? gm->rf[k] : '-', gm->cs[0] ? gm->cs[k] : '-');
      for (x = 0; x < gm->abc->K; x++) fprintf(fp, "%*.*f ", width, precision, P7P_MSC(gm, k, x));
      fputs("\n", fp);
    }
  fputs("\n", fp);

  fputs("# insert state emission parameters\n", fp);
  fputs("# pos cons ", fp);
  for (x = 0; x < gm->abc->K; x++) fprintf(fp, "%*s%c ", width-1, " ", gm->abc->sym[x]);
  fputs("\n", fp);

  fputs("----- ---- ", fp);
  for (x = 0; x < gm->abc->K; x++) fprintf(fp, "%*s ", width, "-------");
  fputs("\n", fp);

  for (k = 0; k <= gm->M; k++) 	/* no insert state I_M, so k=1..M-1 */
    {
      fprintf(fp, "%5d    %c ", k, gm->consensus[k]);
      for (x = 0; x < gm->abc->K; x++) fprintf(fp, "%*.*f ", width, precision, P7P_ISC(gm, k, x));
      fputs("\n", fp);
    }
  fputs("\n", fp);


  fputs("# model composition\n", fp);
  for (x = 0; x < gm->abc->K; x++) fprintf(fp, "%*s%c ", width-1, " ", gm->abc->sym[x]);
  fputs("\n", fp);
  for (x = 0; x < gm->abc->K; x++) fprintf(fp, "%*s ", width, "-------");
  fputs("\n", fp);
  for (x = 0; x < gm->abc->K; x++) fprintf(fp, "%*.*f ", width, precision, gm->compo[x]);
  fputs("\n\n", fp);

  if (gm->evparam[0] == p7_EVPARAM_UNSET) 
    fputs("# statistical parameters: UNSET\n", fp);
  else 
    {
      fputs("# statistical parameters\n", fp);
      fprintf(fp, "   msv mu     = %.4f\n", gm->evparam[p7_MMU]);
      fprintf(fp, "   msv lambda = %.4f\n", gm->evparam[p7_MLAMBDA]);
      fprintf(fp, "   vit mu     = %.4f\n", gm->evparam[p7_VMU]);
      fprintf(fp, "   vit lambda = %.4f\n", gm->evparam[p7_VLAMBDA]);
      fprintf(fp, "   fwd tau    = %.4f\n", gm->evparam[p7_FTAU]);
      fprintf(fp, "   fwd lambda = %.4f\n", gm->evparam[p7_FLAMBDA]);
    }
  fputs("\n", fp);
  
  fputs("# bit score cutoffs\n", fp);
  if (gm->cutoff[p7_GA1] == p7_CUTOFF_UNSET) fprintf(fp, " GA1 = %6s   GA2 = %6s\n",   "unset", "unset");
  else                                       fprintf(fp, " GA1 = %6.1f GA2 = %6.1f\n", gm->cutoff[p7_GA1], gm->cutoff[p7_GA2]);
  if (gm->cutoff[p7_TC1] == p7_CUTOFF_UNSET) fprintf(fp, " TC1 = %6s   TC2 = %6s\n",   "unset", "unset");
  else                                       fprintf(fp, " TC1 = %6.1f TC2 = %6.1f\n", gm->cutoff[p7_TC1], gm->cutoff[p7_TC2]);
  if (gm->cutoff[p7_NC1] == p7_CUTOFF_UNSET) fprintf(fp, " NC1 = %6s   NC2 = %6s\n",   "unset", "unset");
  else                                       fprintf(fp, " NC1 = %6.1f NC2 = %6.1f\n", gm->cutoff[p7_NC1], gm->cutoff[p7_NC2]);
  fputs("\n", fp);
  
  fputs("# derived calculations\n", fp);
  fprintf(fp, "   max_length (nhmmer) = %d\n", gm->max_length);
  fputs("\n", fp);
  
  fputs("#### END\n", fp);
  return eslOK;
}


/* Function:  p7_profile_DecodeT()
 * Synopsis:  Returns printable string representation of transition parameter index code.
 */
char *
p7_profile_DecodeT(int tidx)
{
  switch (tidx) {
  case p7P_MM:  return "tMM";
  case p7P_IM:  return "tIM";
  case p7P_DM:  return "tDM";
  case p7P_LM:  return "tLMk";
  case p7P_GM:  return "tGMk";
  case p7P_MD:  return "tMD";
  case p7P_DD:  return "tDD";
  case p7P_MI:  return "tMI";
  case p7P_II:  return "tII";   
  case p7P_DGE: return "tDgE";
  }
  return NULL;
}


  


/* Function:  p7_profile_Validate()
 *
 * Purpose:   Validates the internals of the generic profile structure
 *            <gm>.
 *            
 * Returns:   <eslOK> if <gm> internals look fine. Returns <eslFAIL>
 *            if something is wrong, and leaves an error message in
 *            <errbuf> if caller passed it non-<NULL>.
 */
int
p7_profile_Validate(const P7_PROFILE *gm, char *errbuf, float tol)
{
  int     k;
  int     M      = gm->M;
  double *pstart = NULL;
  int     status;

  ESL_ALLOC(pstart, sizeof(double) * (gm->M+1));

  /* Validate tsc[0] boundary condition: 
   * tLM, tGM, tDGE are valid transitions at nonexistent node k=0, 
   * because of their off-by-one storage (i.e. tsc[k-1,LM] = tLk->M)
   */
  if (P7P_TSC(gm, 0, p7P_MM)  != -eslINFINITY ||
      P7P_TSC(gm, 0, p7P_IM)  != -eslINFINITY ||
      P7P_TSC(gm, 0, p7P_DM)  != -eslINFINITY ||
      // LM, GM skipped past...
      P7P_TSC(gm, 0, p7P_MD)  != -eslINFINITY ||
      P7P_TSC(gm, 0, p7P_DD)  != -eslINFINITY ||
      P7P_TSC(gm, 0, p7P_MI)  != -eslINFINITY ||
      P7P_TSC(gm, 0, p7P_II)  != -eslINFINITY) ESL_XFAIL(eslFAIL, errbuf, "transition probs at 0 not set properly");
      // DGE skipped.
      
  /* Validate tsc[M] boundary conditions.
   *  t(Mm->D) = 0   as an initialization condition to make Backward work
   *  t(Dm->D) = 0   ditto
   *  t(DGE,k) = t(Dk+1->..->Dm->E) = 0 at k=M and k=M-1
   */  
  if (P7P_TSC(gm, M, p7P_MM)  != -eslINFINITY ||
      P7P_TSC(gm, M, p7P_IM)  != -eslINFINITY ||
      P7P_TSC(gm, M, p7P_DM)  != -eslINFINITY ||
      P7P_TSC(gm, M, p7P_LM)  != -eslINFINITY ||
      P7P_TSC(gm, M, p7P_GM)  != -eslINFINITY ||
      /*... MD DD ... */
      P7P_TSC(gm, M, p7P_MI)  != -eslINFINITY ||
      P7P_TSC(gm, M, p7P_II)  != -eslINFINITY) ESL_XFAIL(eslFAIL, errbuf, "transition probs at M not set properly");

  if (P7P_TSC(gm, M, p7P_MD)  != 0.0f ||  
      P7P_TSC(gm, M, p7P_DD)  != 0.0f ||  
      P7P_TSC(gm, M, p7P_DGE) != 0.0f)  ESL_XFAIL(eslFAIL, errbuf, "transition probs at M not set properly");

  if (M>1 && P7P_TSC(gm, M-1, p7P_DGE) != 0.0f)  ESL_XFAIL(eslFAIL, errbuf, "transition probs at M not set properly");

  /* Validate local entry distribution.
   * this is an implicit probability distribution,
   * corresponding to the implicit local alignment model, and we have
   * to calculate the M(M+1)/2 fragment probabilities accordingly.
   */
  pstart[0] = 0.0;
  for (k = 1; k <= gm->M; k++) {
    pstart[k] = exp(P7P_TSC(gm, k-1, p7P_LM)) * (gm->M - k + 1);     /* multiply p_ij by the number of exits j; note off-by-one storage, L->Mk in tsc[k-1] */
    if (pstart[k] > 1.0 && pstart[k] <= 1.0 + tol) pstart[k] = 1.0;  // There's a special case, where only one M state is occupiable, with entry prob 1.0/(M-k+1). (M-k+1)*exp(log(1.0/M-k+1)) can give 1+epsilon. 
  }
  if (esl_vec_DValidate(pstart, gm->M+1, tol, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "local entry distribution is not normalized properly");

  /* Validate the glocal entry distribution.  This is an explicit
   * probability distribution, corresponding to left wing retraction.
   * To make it sum to one, we also have to add in the probability
   * of the mute cycle G->D1...Dm->E; we put this in pstart[0].
   * [SRE:J9/98].
   */
  for (k = 1; k <= gm->M; k++)
    pstart[k] = exp(P7P_TSC(gm, k-1, p7P_GM));
  /* mute path probability */
  p7_profile_GetMutePathLogProb(gm, &(pstart[0]));
  pstart[0] = exp(pstart[0]);	/* this will often underflow, but when it does, it means that path's negligible in p, and the validation of the sum will still succeed */

  if (esl_vec_DValidate(pstart, gm->M+1, tol, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "glocal entry distribution is not normalized properly");

  free(pstart);
  return eslOK;

 ERROR:
  if (pstart != NULL) free(pstart);
  return eslFAIL;
}


/* Function:  p7_profile_GetMutePathLogProb()
 * Synopsis:  Get the ln prob of G->D1...Dm->E path.
 *
 * Purpose:   Calculate the log probability (in nats) of the mute glocal
 *            path G->D1...Dm->E, a probability mass that profile
 *            configuration removes from the model (by wing retraction
 *            and removal of the D1 state).  We sometimes need to add
 *            this term back, especially in some debugging routines
 *            that check for perfect summation over paths.
 *            
 *            For example, in <p7_profile_Validate()>, we check that
 *            the wing-retracted glocal entry probability distribution
 *            sums to one, and in some "enumeration" unit tests, we
 *            check that the total sum $\sum_x P(x \mid M) = 1.0$ over
 *            sequences $x$.
 * 
 * Args:      gm           - configured profile
 *            ret_mute_lnp - RETURN: log_e (nat) score of G->D1..Dm->E path.    
 *
 * Returns:   <eslOK> on success, and <*ret_mute_lnp> is the answer.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      [SRE:J9/98,100].
 */
int
p7_profile_GetMutePathLogProb(const P7_PROFILE *gm, double *ret_mute_lnp)
{
  int    k;
  double mute_lnp = 0.0;

  mute_lnp = gm->xsc[p7P_G][1];                                   /* G->D1 */
  for (k = 1; k < gm->M; k++) mute_lnp += P7P_TSC(gm, k, p7P_DD); /* D1->D2,...Dm-1->Dm; Dm->E is 1.0 */
  *ret_mute_lnp = mute_lnp;
  return eslOK;
}


/* Function:  p7_profile_Compare()
 * Synopsis:  Compare two profiles for equality.
 *
 * Purpose:   Compare two profiles <gm1> and <gm2> to each other.
 *            Return <eslOK> if they're identical, and <eslFAIL> if
 *            they differ. Floating-point probabilities are 
 *            compared for equality within a fractional tolerance
 *            <tol>.  Only compares the scores, not any annotation
 *            on the profiles.
 */
int
p7_profile_Compare(P7_PROFILE *gm1, P7_PROFILE *gm2, float tol)
{
  int x;

  if (gm1->M       != gm2->M)       return eslFAIL;
  if (gm1->L       != gm2->L)       return eslFAIL;
  if (gm1->nj      != gm2->nj)      return eslFAIL;
  if (gm1->pglocal != gm2->pglocal) return eslFAIL;

  if (esl_vec_FCompare(gm1->tsc, gm2->tsc, (gm1->M+1)*p7P_NTRANS, tol)     != eslOK) return eslFAIL;
  for (x = 0; x < gm1->abc->Kp; x++) 
    if (esl_vec_FCompare(gm1->rsc[x], gm2->rsc[x], (gm1->M+1)*p7P_NR, tol) != eslOK) return eslFAIL;

  for (x = 0; x < p7P_NXSTATES; x++)
    if (esl_vec_FCompare(gm1->xsc[x], gm2->xsc[x], p7P_NXTRANS, tol)       != eslOK) return eslFAIL;

  return eslOK;
}



/* Function:  p7_profile_SameAsSSV()
 * Synopsis:  Round a generic profile's scores to give SSV scores.
 *
 * Purpose:   Set a generic profile's scores so that the reference
 *            Viterbi implementation will give the same score as
 *            <p7_SSVFilter()>, by rescaling and rounding using a
 *            factor of <scale_b> (that is, pass <om->scale_b> from
 *            the vector profile you're matching to).
 *
 *            Caller must first set <gm> to unihit local mode. Then
 *            this routine sets: all t_MM scores = 0; all other core
 *            transitions = -inf; all <t_BMk> entries uniformly <log
 *            2/(M(M+1))>; <tCC, tNN, tJJ> scores 0 (caller will apply
 *            the 2 nat approximation); rounded in the same way as the
 *            8-bit limited precision.
 *            
 *            To convert a generic Viterbi score <gsc> calculated with
 *            this profile to a nat score that should match
 *            SSVFilter() exactly, do <(gsc / om->scale_b) - 2.0>.
 *            
 *            <gm> is irrevocably altered. You can only call this once
 *            on a given <gm>.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_profile_SameAsSSV(P7_PROFILE *gm, float scale_b)
{
  int    k,x;
  float  tbm = roundf(scale_b * (log(2.0f / ((float) gm->M * (float) (gm->M+1)))));

  /* Transitions */
  esl_vec_FSet(gm->tsc, p7P_NTRANS * gm->M, -eslINFINITY);
  for (k = 1; k <  gm->M; k++) P7P_TSC(gm, k, p7P_MM) = 0.0f;
  for (k = 0; k <  gm->M; k++) P7P_TSC(gm, k, p7P_LM) = tbm;
  
  /* Emissions */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 0; k <= gm->M; k++)
      {
	gm->rsc[x][k*2]   = (gm->rsc[x][k*2] <= -eslINFINITY) ? -eslINFINITY : roundf(scale_b * gm->rsc[x][k*2]);
	gm->rsc[x][k*2+1] = 0;	/* insert scores are implicitly zero in the filters. */
      }	

   /* Specials */
  for (k = 0; k < p7P_NXSTATES; k++)
    for (x = 0; x < p7P_NXTRANS; x++)
      gm->xsc[k][x] = (gm->xsc[k][x] <= -eslINFINITY) ? -eslINFINITY : roundf(scale_b * gm->xsc[k][x]);

  /* NN, CC, JJ hardcoded 0 in limited precision */
  gm->xsc[p7P_N][p7P_LOOP] =  gm->xsc[p7P_J][p7P_LOOP] =  gm->xsc[p7P_C][p7P_LOOP] = 0;

  return eslOK;
}


/* Function:  p7_profile_SameAsVF()
 * Synopsis:  Round a generic profile to match ViterbiFilter scores.
 *
 * Purpose:   Round all the scores in a generic (lspace) <P7_PROFILE> <gm> in
 *            exactly the same way that the scores in a
 *            <P7_OPROFILE>, by rescaling and rounding using a factor
 *            of <scale_w> (that is, pass <om->scale_w> from the vector 
 *            profile you're matching). Then we can test that two profiles
 *            give identical internal scores in testing, say,
 *            <p7_ViterbiFilter()> against <p7_GViterbi()>. 
 *            
 *            The 3nat approximation is used; NN=CC=JJ=0, and 3 nats are
 *            subtracted at the end to account for their contribution.
 *            
 *            To convert a generic Viterbi score <gsc> calculated with this profile
 *            to a nat score that should match ViterbiFilter() exactly,
 *            do <(gsc / om->scale_w) - 3.0>.
 *
 *            <gm> must be the same profile that <om> was constructed from.
 * 
 *            <gm> is irrevocably altered by this call. 
 *            
 *            Do not call this more than once on any given <gm>! 
 *
 * Args:      <gm>      - generic profile that <om> was built from.          
 *            <scale_w> - i.e. <om->scale_w> from a vector profile we're matching to
 *
 * Returns:   <eslOK> on success.
 */
int
p7_profile_SameAsVF(P7_PROFILE *gm, float scale_w)
{
  int k;
  int x;

  /* Transitions */
  /* <= -eslINFINITY test is used solely to silence compiler about floating point comparison.
   * really testing == -eslINFINITY 
   */
  for (x = 0; x < gm->M*p7P_NTRANS; x++)
    gm->tsc[x] = (gm->tsc[x] <= -eslINFINITY) ? -eslINFINITY : roundf(scale_w * gm->tsc[x]);
  
  /* Enforce the rule that no II can be 0; max of -1 */
  for (x = p7P_II; x < gm->M*p7P_NTRANS; x += p7P_NTRANS) 
    if (gm->tsc[x] == 0.0) gm->tsc[x] = -1.0;

  /* Emissions */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 0; k <= gm->M; k++)
      {
	gm->rsc[x][k*2]   = (gm->rsc[x][k*2]   <= -eslINFINITY) ? -eslINFINITY : roundf(scale_w * gm->rsc[x][k*2]);
	gm->rsc[x][k*2+1] = 0.0;	/* insert score: VF makes it zero no matter what. */
      }	

  /* Specials */
  for (k = 0; k < p7P_NXSTATES; k++)
    for (x = 0; x < p7P_NXTRANS; x++)
      gm->xsc[k][x] = (gm->xsc[k][x] <= -eslINFINITY) ? -eslINFINITY : roundf(scale_w * gm->xsc[k][x]);

  /* 3nat approximation: NN, CC, JJ hardcoded 0 in limited precision */
  gm->xsc[p7P_N][p7P_LOOP] =  gm->xsc[p7P_J][p7P_LOOP] =  gm->xsc[p7P_C][p7P_LOOP] = 0.0;
  return eslOK;
}




/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7PROFILE_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_random.h"

#include "base/p7_bg.h"
#include "build/modelsample.h"
#include "search/modelconfig.h"

static void
utest_Compare(void)
{
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(42);
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm  = NULL;
  P7_BG          *bg   = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_PROFILE     *gm2  = NULL;
  int             M    = 200;
  int             L    = 400;

  p7_modelsample(r, M, abc, &hmm); /* master and worker's sampled profiles are identical */
  bg  = p7_bg_Create(abc);

  gm  = p7_profile_Create(hmm->M, abc);
  gm2 = p7_profile_Create(hmm->M, abc);

  p7_profile_Config(gm,  hmm, bg);
  p7_profile_Config(gm2, hmm, bg);

  p7_profile_SetLength(gm,  L);
  p7_profile_SetLength(gm2, L);

  if (p7_profile_Compare(gm, gm2, 0.001) != eslOK) p7_Die("identical profile comparison failed");
  
  p7_profile_Destroy(gm);
  p7_profile_Destroy(gm2);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return;
}


#endif /*p7PROFILE_TESTDRIVE*/

/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7PROFILE_TESTDRIVE
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",              0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for p7_profile.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  fprintf(stderr, "## %s\n", argv[0]);

  utest_Compare();

  fprintf(stderr, "#  status = ok\n");

  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7PROFILE_TESTDRIVE*/


/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef p7PROFILE_EXAMPLE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "very verbose: dump profile object",                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "example for profile object code";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  float           ftol    = 1e-4;        /* floating-point tolerance for checking parameters against expected probs or summing to 1 */
  char            errbuf[eslERRBUFSIZE];

  /* Read in one HMM; sets alphabet to the HMM's alphabet */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Set up a null model */
  bg = p7_bg_Create(abc);

  /* Allocate and configure a profile from HMM and null model */
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, 400);    /* 400 is arbitrary here; this is whatever your target seq length L is */
  
  printf("profile memory consumed: %" PRId64 " bytes\n", (int64_t) p7_profile_Sizeof(gm));

  /* Debugging tools allow dumping, validating the object */
  if (p7_profile_Validate(gm, errbuf, ftol) != eslOK) p7_Fail("profile validation failed\n  %s\n", errbuf);
  
  if (esl_opt_GetBoolean(go, "--vv"))
    p7_profile_Dump(stdout, gm);

  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /* p7PROFILE_EXAMPLE */

