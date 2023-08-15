/* E2_GMX implementation: a generic dynamic programming matrix
 *
 * Contents:
 *   1. The <E2_GMX> object
 *   2. Debugging aids
 *   3. Unit tests
 *   4. Test driver
 *   5. Copyright and license information
 */

#include <stdio.h>		/* FILE */
#include <string.h>

#include "p7_config.h"
#include "hmmer.h"

#include "e2.h"
#include "e2_gmx.h"

static int print_transition(FILE *ofp, E2_GMX *gx, enum e2g_cells_e e2cell, int i, int jstart, int jend, int kstart, int kend, int flags, int width, int precision);

/*****************************************************************
 *= 1. The <E2_GMX> object.
 *****************************************************************/

/* Function:  e2_gmx_Create()
 * Synopsis:  Allocate a new <E2_GMX>.
 *
 * Purpose:   Allocate a reusable, resizeable <E2_GMX> for 
 *            sequences up to length <allocL+2>.
 *            
 *            We've set this up so it should be easy to allocate
 *            aligned memory, though we're not doing this yet.
 *
 * Returns:   a pointer to the new <E2_GMX>.
 *
 * Throws:    <NULL> on allocation error.
 */
E2_GMX *
e2_gmx_Create(int M, int allocL1, int allocL2)
{
  E2_GMX *gx = NULL;
  int     r,x;
  int     status;

  ESL_ALLOC(gx, sizeof(E2_GMX));
  gx->dp_mem = NULL;
  gx->dp     = NULL;
 
  gx->allocL = (allocL1+1) * (allocL2+1);
  gx->allocW = (M+1) * e2G_NSCELLS + e2G_NXCELLS;
  gx->allocN = (int64_t) gx->allocL * (int64_t) gx->allocW;

  ESL_ALLOC(gx->dp_mem,  sizeof(float)   * gx->allocN);
  ESL_ALLOC(gx->dp,      sizeof(float *) * gx->allocN);
  for (r = 0; r < gx->allocL; r++)
    gx->dp[r] = gx->dp_mem + (r * gx->allocW);
  
  /* Initialize memory that's allocated but unused, only to keep
   * valgrind and friends happy.
   */
  for (r = 0; r < gx->allocL; r++) 
    for (x = 0; x < e2G_NSCELLS; x++)
      gx->dp[r][x] = -eslINFINITY;

  gx->validL = gx->allocL;
  gx->M      = M;
  gx->Lrow   = allocL1;
  gx->Lcol   = allocL2;
  gx->L      = (allocL1+1)*(allocL2+1);
  gx->type   = E2_UNSET;
  gx->ncells = (uint64_t) (gx->allocL);
  return gx;

 ERROR:
  if (gx != NULL) e2_gmx_Destroy(gx);
  return NULL;
}

/* Function:  e2_gmx_GrowTo()
 * Synopsis:  Assure that DP matrix is big enough.
 *
 * Purpose:   Assures that a DP matrix <gx> is allocated
 *            for  a sequence of length up to <L>; 
 *            reallocates if necessary.
 *            
 *            This function does not respect the configured
 *            <RAMLIMIT>; it will allocate what it's told to
 *            allocate. 
 *
 * Returns:   <eslOK> on success, and <gx> may be reallocated upon
 *            return; any data that may have been in <gx> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <gx> must be assumed to be invalidated.
 */
int
e2_gmx_GrowTo(E2_GMX *gx, int M, int L1, int L2)
{
  int      W        = (M+1) * e2G_NSCELLS + e2G_NXCELLS;
  int      L        = (L1+1)*(L2+1);
  uint64_t N        = (int64_t) L * (int64_t) W;
  int      do_reset = FALSE;
  int      r,x;
  int      status;

  /* are we already big enough? */
  if (W <= gx->allocW && L <= gx->validL) return eslOK;
  
  /* must we reallocate the matrix cells? */
  if (N > gx->allocN)
    {
      ESL_REALLOC(gx->dp_mem, sizeof(float) * N);
      gx->allocN = N;
      do_reset   = TRUE;
    }
  
  /* must we reallocate the row pointers? */
  if (L > gx->allocL)
    {
      ESL_REALLOC(gx->dp, sizeof(float *) * L);
      gx->allocL = L;
      do_reset   = TRUE;
    }

  /* must we widen the rows? */
  if (W > gx->allocW) do_reset = TRUE;

  /* must we set some more valid row pointers? */
  if (L > gx->validL) do_reset = TRUE;

  /* resize rows, reset valid row pointers */
  if (do_reset)
    {
      gx->allocW = W;
      gx->validL = ESL_MIN(gx->allocL, (int) ( gx->allocN / (uint64_t) gx->allocW));
      for (r = 0; r < gx->validL; r++)
	{
	  gx->dp[r] = gx->dp_mem + (r * gx->allocW);
	  for (x = 0; x < e2G_NSCELLS; x++)  
	    gx->dp[r][x] = -eslINFINITY;
	}
    }
  gx->M = M;
  gx->L = (L1+1)*(L2+1);
  return eslOK;

 ERROR:
  return status;
}

/* Function:  e2_gmx_Sizeof()
 * Synopsis:  Returns the allocation size of DP matrix, in bytes.
 */
size_t 
e2_gmx_Sizeof(E2_GMX *gx)
{
  size_t n = 0;
  
  n += sizeof(E2_GMX);
  n += gx->allocN * sizeof(float);     /* dp_mem */
  n += gx->allocL * sizeof(float *);	/* dp[]   */
  return n;
}



/* Function:  e2_gmx_Reuse()
 * Synopsis:  Recycle a generic DP matrix.
 *
 * Purpose:   Recycles <gx> for reuse.
 *
 * Returns:   <eslOK> on success.
 */
int
e2_gmx_Reuse(E2_GMX *gx)
{
  /* not much to do here. The memory rearrangement for a new seq is all in GrowTo(). */
  gx->Lcol = 0;
  gx->L    = 0;
  gx->type = E2_UNSET;
   return eslOK;
}


/* Function:  e2_gmx_Destroy()
 * Synopsis:  Frees a DP matrix.
 *
 * Purpose:   Frees a <E2_GMX>.
 *
 * Returns:   (void)
 */
void
e2_gmx_Destroy(E2_GMX *gx)
{
  if (gx == NULL) return;
  if (gx->dp      != NULL)  free(gx->dp);
  if (gx->dp_mem  != NULL)  free(gx->dp_mem);
  free(gx);
  return;
}

/*****************************************************************
 * 2. Debugging aids
 *****************************************************************/
/* Function:  e2_gmx_Compare()
 * Synopsis:  Compare two DP matrices for equality within given tolerance.
 *
 * Purpose:   Compare all the values in DP matrices <gx1> and <gx2> using
 *            <esl_FCompare()> and relative epsilon <tolerance>. If any
 *            value pairs differ by more than the acceptable <tolerance>
 *            return <eslFAIL>.  If all value pairs are identical within
 *            tolerance, return <eslOK>. 
 */
int
e2_gmx_Compare(E2_GMX *gx1, E2_GMX *gx2, float tolerance)
{
  int i,j,x;
  int k, s;

  if (gx1->Lrow != gx2->Lrow) return eslFAIL;
  if (gx1->Lcol != gx2->Lcol) return eslFAIL;
  if (gx1->type != gx2->type) return eslFAIL;
  
  for (i = 0; i <= gx1->Lrow; i++)
  {
      for (j = 1; j <= gx1->Lcol; j++) /* k=0 is a boundary; doesn't need to be checked */
      {
	x = ID(i,j,gx1->Lcol);
 
	for (k = 0; k <= gx1->M; k++)   
	  for (s = 0; s < e2G_NSCELLS; s++)
	    if ( esl_FCompareAbs(E2G_MX(gx1,x,k,s), E2G_MX(gx2,x,k,s), tolerance) == eslFAIL) return eslFAIL;
	for (s = 0; s < e2G_NXCELLS; s++)
	  if ( esl_FCompareAbs(E2G_XMX(gx1,i,s), E2G_XMX(gx2,i,s), tolerance) == eslFAIL)   return eslFAIL; 
	
      }
  }
  return eslOK;	
}



/* Function:  e2_gmx_Dump()
 * Synopsis:  Dump a DP matrix to a stream, for diagnostics.
 *
 * Purpose:   Dump matrix <gx> to stream <fp> for diagnostics.
 *
 *            <flags> control some optional output behaviors, as follows:
 *              | <e2_HIDE_SPECIALS> | don't show scores for <ENJBC> states  |
 *              | <e2_SHOW_LOG>      | <gx> is in probs; show as log probs   |
 */
int
e2_gmx_Dump(FILE *ofp, E2_GMX *gx, int flags)
{
  return e2_gmx_DumpWindow(ofp, gx, 0, gx->Lrow, 0, gx->Lcol, 0, gx->M, flags);
}

/* Function:  e2_gmx_DumpWindow()
 * Synopsis:  Dump a window of a DP matrix to a stream for diagnostics.
 *
 * Purpose:   Dump a window of matrix <gx> to stream <fp> for diagnostics,
 *            from row <istart> to <iend>, from column <kstart> to <kend>.
 *            
 *            Asking for <0..L,0..M> with <flags=e2_SHOW_SPECIALS> is the
 *            same as <e2_gmx_Dump()>.
 *            
 *            <flags> control some optional output behaviors, as follows:
 *              | <e2_HIDE_SPECIALS> | don't show scores for <ENJBC> states  |
 *              | <e2_SHOW_LOG>      | <gx> is in probs; show as log probs   |
 *  
 * Returns:   <eslOK> on success.
 */
int
e2_gmx_DumpWindow(FILE *ofp, E2_GMX *gx, int istart, int iend, int jstart, int jend, int kstart, int kend, int flags)
{
  int   width     = 9;
  int   precision = 4;
  int   i,j;
  int   k,s,x;

  /* Header */
  fprintf(ofp, "     ");
  for (j = jstart; j <= jend; j++) fprintf(ofp, "%*d ", width, j);
  for (j = jstart; j <= jend; j++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");
  
  fprintf(ofp, "       ");
  for (k = kstart; k <= kend;  k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < e2G_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");
  
  /* DP matrix data */
  for (i = istart; i <= iend; i++)
    for (s = 0; s < e2G_NSCELLS; s++)
      print_transition(ofp, gx, s, i, jstart, jend, kstart, kend, flags, width, precision);
   
  return eslOK;
}

int
e2_gmx_NTFromTag(enum e2g_cells_e *ret_e2cell, char *tag)
{
 
  if (ret_e2cell == NULL) return eslFAIL;

  if      (strcmp(tag, "BB")==0)  *ret_e2cell = e2G_BB;
  else if (strcmp(tag, "IB")==0)  *ret_e2cell = e2G_IB;
  else if (strcmp(tag, "SS")==0)  *ret_e2cell = e2G_SS;
  else if (strcmp(tag, "DS")==0)  *ret_e2cell = e2G_DS;
  else if (strcmp(tag, "IS")==0)  *ret_e2cell = e2G_IS;
  else if (strcmp(tag, "SD")==0)  *ret_e2cell = e2G_SD;
  else if (strcmp(tag, "DD")==0)  *ret_e2cell = e2G_DD;
  else if (strcmp(tag, "ID")==0)  *ret_e2cell = e2G_ID;
  else if (strcmp(tag, "BI")==0)  *ret_e2cell = e2G_BI;
  else if (strcmp(tag, "SI")==0)  *ret_e2cell = e2G_SI;
  else if (strcmp(tag, "DI")==0)  *ret_e2cell = e2G_DI;
  else if (strcmp(tag, "II")==0)  *ret_e2cell = e2G_II;

  else return eslFAIL;

  return eslOK;
}
 
int
e2_gmx_NTtag(enum e2g_cells_e e2cell, char **ret_tag)
{
  char tag[STRSIZE];

  if      (e2cell == e2G_BB) strcpy(tag, "BB");
  else if (e2cell == e2G_IB) strcpy(tag, "IB");
  else if (e2cell == e2G_SS) strcpy(tag, "SS");
  else if (e2cell == e2G_DS) strcpy(tag, "DS");
  else if (e2cell == e2G_IS) strcpy(tag, "IS");
  else if (e2cell == e2G_SD) strcpy(tag, "SD");
  else if (e2cell == e2G_DD) strcpy(tag, "DD");
  else if (e2cell == e2G_ID) strcpy(tag, "ID");
  else if (e2cell == e2G_BI) strcpy(tag, "BI");
  else if (e2cell == e2G_SI) strcpy(tag, "SI");
  else if (e2cell == e2G_DI) strcpy(tag, "DI");
  else if (e2cell == e2G_II) strcpy(tag, "II");

  else return eslFAIL;

  if (ret_tag != NULL) *ret_tag = tag;
  return eslOK;
}

/* Function:  p7_refmx_DecodeSpecial()
 * Synopsis:  Convert special state code to string for debugging output.
 *
 * Purpose:   Given a special state code (such as <p7R_E>), return a
 *            string label, suitable for state label in debugging output.
 *
 * Args:      type  - special state code, such as <p7R_E>
 *
 * Returns:   ptr to a static string representation
 *
 * Throws:    (no abnormal error conditions)
 */
char *
e2_gmx_DecodeSpecial(int type)
{
  switch (type) {
  case e2G_EE:  return "EE";
  case e2G_N1:  return "N1";
  case e2G_N2:  return "N2";
  case e2G_J1:  return "J1";
  case e2G_J2:  return "J2";
  case e2G_C1:  return "C1";
  case e2G_C2:  return "C2";
  case e2G_JJ1: return "JJ1";
  case e2G_JJ2: return "JJ2";
  case e2G_CC1: return "CC1";
  case e2G_CC2: return "CC2";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_REFMX special state code %d\n", type);
  return NULL;
}
char *
e2_gmx_DecodeState(enum e2g_cells_e e2cell)
{
  switch (e2cell) {
  case e2G_BB: return "BB";
  case e2G_IB: return "IB";
  case e2G_SS: return "SS";
  case e2G_DS: return "DS";
  case e2G_IS: return "IS";
  case e2G_SD: return "SD";
  case e2G_DD: return "DD";
  case e2G_ID: return "ID";
  case e2G_BI: return "BI";
  case e2G_SI: return "SI";
  case e2G_DI: return "DI";
  case e2G_II: return "II";
  case e2G_ii: return "ii";
 }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such E2_GMX state code %d\n", e2cell);
  return NULL;

}

int
print_transition(FILE *ofp, E2_GMX *gx, enum e2g_cells_e s, int i, int jstart, int jend, int kstart, int kend, int flags, int width, int precision)
{
  char *tag = NULL;
  int   j,x,k;
  float val;

  e2_gmx_NTtag(s, &tag);
  
  fprintf(ofp, "%3d %s ", i, tag);
  for (j = jstart; j <= jend; j++)  
    {
      x = ID(i,j,gx->Lcol);

      for (k = kstart; k <= kend; k++) {
	val = gx->dp[x][k*e2G_NSCELLS + s];
	if (flags & e2_SHOW_LOG) val = log(val);
	fprintf(ofp, "%*.*f ", width, precision, val);
      }
    }
  fprintf(ofp, "\n");

  return eslOK;
}


/*------------- end, (most) debugging tools ---------------------*/


/*****************************************************************
 * 3. Validation of matrices.
 *****************************************************************/
/* Gets its own section because it's a little fiddly and detailed.
 * See e2_gmx.h for notes on the patterns that get tested.
 */
static inline int
validate_dimensions(E2_GMX *gmx, char *errbuf)
{
  if ( gmx->M <= 0)               ESL_FAIL(eslFAIL, errbuf, "nonpositive M");
  if ( gmx->L <= 0)               ESL_FAIL(eslFAIL, errbuf, "nonpositive L");
  if ( gmx->L+1    > gmx->validL) ESL_FAIL(eslFAIL, errbuf, "L+1 is larger than validL");
  if ( gmx->validL > gmx->allocL) ESL_FAIL(eslFAIL, errbuf, "validL larger than allocL");
  if ( (gmx->M+1)*e2G_NSCELLS+e2G_NXCELLS < gmx->allocW) ESL_FAIL(eslFAIL, errbuf, "M is too large for allocW");
  return eslOK;
}

static inline int 
validate_mainstate(E2_GMX *gmx, int i, int k, int s, float val, char *errbuf)
{
  if (E2G_MX(gmx, i, k, s) != val) ESL_FAIL(eslFAIL, errbuf, "expected %f at i=%d, k=%d, %s", val, i, k, e2_gmx_DecodeState(s));
  return eslOK;
}

static inline int
validate_special(E2_GMX *gmx, int i, int s, float val, char *errbuf)
{
  if (E2G_XMX(gmx, i, s) != val) ESL_FAIL(eslFAIL, errbuf, "expected %f at i=%d, %s", val, i, e2_gmx_DecodeSpecial(s));
  return eslOK;
}

static inline int
validate_no_nan(E2_GMX *gmx, char *errbuf)
{
  int i,k,s;
  for (i = 0; i <= gmx->L; i++)
    for (k = 0; k <= gmx->M; k++)
      {
	for (s = 0; s < e2G_NSCELLS; s++)
	  if (isnan(E2G_MX(gmx, i, k, s))) ESL_FAIL(eslFAIL, errbuf, "found NaN at i=%d, k=%d, %s", i, k, e2_gmx_DecodeState(s));      
	for (s = 0; s < e2G_NXCELLS; s++)
	  if (isnan(E2G_XMX(gmx, i, s)))   ESL_FAIL(eslFAIL, errbuf, "found NaN at i=%d, %s", i, e2_gmx_DecodeSpecial(s));      
     }
  return eslOK;
}


static inline int
validate_column_zero(E2_GMX *gmx, float val, char *errbuf)
{
  int i, s, status;
  for (i = 0; i <= gmx->L; i++)
    for (s = 0; s < e2G_NSCELLS; s++)
      if (( status = validate_mainstate(gmx, i, 0, s, val, errbuf)) != eslOK) return status;
  return eslOK;
}

static inline int
validate_forward(E2_GMX *gmx, char *errbuf)
{
  int k,s,i;
  int status;

 if (( status = validate_column_zero(gmx, -eslINFINITY, errbuf)) != eslOK) return status;

   /* Row i=0 */
  for (k = 1; k <= gmx->M; k++)
    for (s = 0; s <= e2G_NSCELLS; s++)
      if ( (status = validate_mainstate(gmx, 0, k, s, -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_EE,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_N1,     0.0f,         errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_N2,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_J1,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_J2,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_C1,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_C2,     -eslINFINITY, errbuf)) != eslOK) return status;

  /* CC/JJ are only for decoding; Forward matrix has them as -inf */
  for (i = 0; i <= gmx->L; i++)
    {
      if ( ( status = validate_special(gmx, 0, e2G_NN2,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_JJ1,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_JJ2,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_CC1,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_CC2,    -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}

static inline int
validate_backward(E2_GMX *gmx, char *errbuf)
{
  int k,s,i;
  int status;

 if (( status = validate_column_zero(gmx, -eslINFINITY, errbuf)) != eslOK) return status;

  /* Row i=0 */
  for (k = 1; k <= gmx->M; k++)
    for (s = 0; s < e2G_NSCELLS; s++)
      if ( (status = validate_mainstate(gmx, 0, k, s, -eslINFINITY, errbuf)) != eslOK) return status;

  if ( (status = validate_special  (gmx, gmx->L,    e2G_N1,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (gmx, gmx->L,    e2G_N2,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (gmx, gmx->L,    e2G_J1,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (gmx, gmx->L,    e2G_J2,  -eslINFINITY, errbuf)) != eslOK) return status;

  /* CC/JJ are only for decoding; Backward matrix has them as -inf */
  for (i = 0; i <= gmx->L; i++)
    {
      if ( ( status = validate_special(gmx, 0, e2G_NN2,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_JJ1,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_JJ2,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_CC1,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_CC2,    -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}

static inline int
validate_decoding(E2_GMX *gmx, char *errbuf)
{
  int k,i,s;
  int status;

  /* All of k=0 column is 0.0 */
  if (( status = validate_column_zero(gmx, 0.0, errbuf)) != eslOK) return status;  

 /* Row i=0 */
  for (k = 1; k <= gmx->M; k++)
    for (s = 0; s <= e2G_NSCELLS; s++)
      if ( (status = validate_mainstate(gmx, 0, k, s, -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_EE,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_N2,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_J1,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_J2,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_C1,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(gmx, 0, e2G_C2,     -eslINFINITY, errbuf)) != eslOK) return status;

  if ( (status = validate_special  (gmx, gmx->L,    e2G_N1,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (gmx, gmx->L,    e2G_N2,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (gmx, gmx->L,    e2G_J1,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (gmx, gmx->L,    e2G_J2,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (gmx, gmx->L,    e2G_C1,  -eslINFINITY, errbuf)) != eslOK) return status;
 
  /* and throughout, CC/JJ are only for decoding; Alignment matrix has them as -inf */
  for (i = 0; i <= gmx->L; i++)
    {
      if ( ( status = validate_special(gmx, 0, e2G_JJ1,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_JJ2,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_CC1,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_CC2,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(gmx, 0, e2G_NN2,    -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}

static inline int
validate_alignment(E2_GMX *gmx, char *errbuf)
{
  int k,s;
  int status;

  if (( status = validate_column_zero(gmx, -eslINFINITY, errbuf)) != eslOK) return status;  

  /* Row i=0 */
  for (k = 1; k <= gmx->M; k++)
    for (s = 0; s <= e2G_NSCELLS; s++)
      if ( (status = validate_mainstate(gmx, 0, k, s, -eslINFINITY, errbuf)) != eslOK) return status;

  return eslOK;
}


/* Function:  p7_refmx_Validate()
 * Synopsis:  Validates a reference DP matrix.
 *
 * Purpose:   Validates the internals of the
 *            DP matrix <gmx>. Returns <eslOK> if
 *            it passes. Returns <eslFAIL> if it fails,
 *            and sets <errbuf> to contain an explanation,
 *            if caller provided an <errbuf>.
 *
 * Args:      gmx    - reference DP matrix to validate.
 *            errbuf - allocated space for err msg, or NULL
 *
 * Returns:   <eslOK> on success; <eslFAIL> on failure, and
 *            <errbuf> (if provided) contains explanation.
 *
 * Xref:      See notes in <p7_refmx.h> for an explanation
 *            of the particular patterns we're checking.
 */
int
e2_gmx_Validate(E2_GMX *gmx, char *errbuf)
{
  int status;

  if ( (status = validate_dimensions(gmx, errbuf)) != eslOK) return status;
  if ( (status = validate_no_nan    (gmx, errbuf)) != eslOK) return status;
  
  switch (gmx->type) {
  case E2_FORWARD:   if ( (status = validate_forward    (gmx,               errbuf)) != eslOK) return status;  break;
  case E2_VITERBI:   if ( (status = validate_forward    (gmx,               errbuf)) != eslOK) return status;  break; // Viterbi has same pattern as Forward
  case E2_BACKWARD:  if ( (status = validate_backward   (gmx,               errbuf)) != eslOK) return status;  break;
  case E2_DECODING:  if ( (status = validate_decoding   (gmx,               errbuf)) != eslOK) return status;  break;     
  case E2_ALIGNMENT: if ( (status = validate_alignment  (gmx,               errbuf)) != eslOK) return status;  break; 
  case E2_UNSET:     if ( (status = validate_column_zero(gmx, -eslINFINITY, errbuf)) != eslOK) return status;  break;
  default:            ESL_FAIL(eslFAIL, errbuf, "no such reference DP algorithm type %d", gmx->type);
  }
  return eslOK;
}
/*------------------ end, validation ----------------------------*/
