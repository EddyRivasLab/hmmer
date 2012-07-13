/* P7_REFMX implementation: dynamic programming matrix for reference 
 * implementation of dual-mode (local/glocal) alignment.
 *
 * See also: reference_fwdback.c
 * 
 * Contents:
 *   1. The <P7_REFMX> object
 *   2. Debugging aids
 *   3. Validation 
 *   4. Copyright and license information.
 */

#include "p7_config.h"

#include "hmmer.h"
#include "p7_refmx.h"


/*****************************************************************
 * 1. The <P7_REFMX> object
 *****************************************************************/

/* Function:  p7_refmx_Create()
 * Synopsis:  Create a new <P7_REFMX> DP matrix.
 *
 * Purpose:   Create a new <P7_REFMX> matrix for a model of
 *            length <M> and a target sequence of length
 *            <L>.
 *
 * Args:      M - model length
 *            L - target sequence length
 * 
 * Returns:   ptr to the new <P7_REFMX>.
 *
 * Throws:    <NULL> on any allocation failure.
 */
P7_REFMX *
p7_refmx_Create(int M, int L)
{
  P7_REFMX *rmx = NULL;
  int      r,x;
  int      status;

  ESL_ALLOC(rmx, sizeof(P7_REFMX));
  rmx->dp_mem = NULL;
  rmx->dp     = NULL;

  rmx->allocR = L+1;
  rmx->allocW = (M+1) * p7R_NSCELLS + p7R_NXCELLS;
  rmx->allocN = (int64_t) rmx->allocR * (int64_t) rmx->allocW;

  ESL_ALLOC(rmx->dp_mem, sizeof(float  ) * rmx->allocN);
  ESL_ALLOC(rmx->dp,     sizeof(float *) * rmx->allocR);
  for (r = 0; r < rmx->allocR; r++)
    rmx->dp[r] = rmx->dp_mem + (r * rmx->allocW);

  /* Initialize all k=0 cells to -inf. */
  for (r = 0; r < rmx->allocR; r++)
    for (x = 0; x < p7R_NSCELLS; x++)
      rmx->dp[r][x] = -eslINFINITY;

  rmx->validR = rmx->allocR;
  rmx->M      = 0;
  rmx->L      = 0;
  rmx->type   = p7R_UNSET;
  return rmx;

 ERROR:
  if (rmx) p7_refmx_Destroy(rmx);
  return NULL;
}

/* Function:  p7_refmx_GrowTo()
 * Synopsis:  Efficiently reallocate a <P7_REFMX>.
 *
 * Purpose:   Efficiently reallocate the matrix <rmx> to a new
 *            DP problem size, for model length <M> and target 
 *            sequence length <L>. Reuse the existing allocation
 *            as much as possible, to minimize reallocation calls.
 *
 * Args:      rmx - existing DP matrix
 *            M   - new model length
 *            L   - new target sequence length
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on memory allocation failure.
 */
int
p7_refmx_GrowTo(P7_REFMX *rmx, int M, int L)
{
  int      W        = (M+1) * p7R_NSCELLS + p7R_NXCELLS;
  int      R        = L+1;
  uint64_t N        = (int64_t) R * (int64_t) W;
  int      do_reset = FALSE;
  int      r,x;
  int      status;

  /* are we already big enough? */
  if (W <= rmx->allocW && R <= rmx->validR) return eslOK;

  /* must we reallocate the matrix cells? */
  if (N > rmx->allocN)
    {
      ESL_REALLOC(rmx->dp_mem, sizeof(float) * N);
      rmx->allocN = N;
      do_reset    = TRUE;
    }
  
  /* must we reallocate the row pointers? */
  if (R > rmx->allocR)
    {
      ESL_REALLOC(rmx->dp, sizeof(float *) * R);
      rmx->allocR = R;
      do_reset    = TRUE;
    }

  /* must we widen the rows? */
  if (W > rmx->allocW) do_reset = TRUE;

  /* must we set some more valid row pointers? */
  if (R > rmx->validR) do_reset = TRUE;

  /* resize rows, reset valid row pointers */
  if (do_reset)
    {
      rmx->allocW = W;
      rmx->validR = ESL_MIN(rmx->allocR, (int) ( rmx->allocN / (uint64_t) rmx->allocW));
      for (r = 0; r < rmx->validR; r++)
	{
	  rmx->dp[r] = rmx->dp_mem + (r * rmx->allocW);
	  for (x = 0; x < p7R_NSCELLS; x++)  
	    rmx->dp[r][x] = -eslINFINITY;
	}
    }
  rmx->M = 0;
  rmx->L = 0;
  return eslOK;

 ERROR:
  return status;
}
	

/* Function:  p7_refmx_Sizeof()
 * Synopsis:  Returns current allocated size of reference DP matrix, in bytes.
 */
size_t
p7_refmx_Sizeof(const P7_REFMX *rmx)
{
  size_t n = sizeof(P7_REFMX);
  n += rmx->allocN * sizeof(float);     /* dp_mem */
  n += rmx->allocR * sizeof(float *);	/* dp[]   */
  return n;
}

/* Function:  p7_refmx_MinSizeof()
 * Synopsis:  Returns minimum required size of reference DP matrix, in bytes.
 * 
 * Purpose:   Calculate and return the minimal required size of 
 *            a reference DP matrix for a comparison of a 
 *            profile of length <M> to a sequence of length <L>.
 *            
 *            This could overflow, on machines with 32b size_t. 
 */
size_t
p7_refmx_MinSizeof(int M, int L)
{
  size_t n = sizeof(P7_REFMX);
  n += sizeof(float)   * (size_t) (L+1) * (size_t) ( (M+1)*p7R_NSCELLS + p7R_NXCELLS); // dp_mem
  n += sizeof(float *) * (size_t) (L+1);                                               // dp[]
  return n;
}

/* Function:  p7_refmx_Reuse()
 * Synopsis:  Finish using a <P7_REFMX> without deallocation.
 *
 * Purpose:   Caller says it is done using <rmx> for now, but is
 *            soon going to use it again for a new problem;
 *            so reinitialize, saving deallocation/reallocation.
 *            Equiv to <p7_refmx_Destroy(); p7_refmx_Create()> in 
 *            effect, but faster.
 *
 * Args:      rmx - matrix that will be reused.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_refmx_Reuse(P7_REFMX *rmx)
{
  rmx->M    = 0;
  rmx->L    = 0;
  rmx->type = p7R_UNSET;
  return eslOK;
}


/* Function:  p7_refmx_Destroy()
 * Synopsis:  Free a <P7_REFMX>.
 *
 * Purpose:   Free the matrix <rmx>.
 */
void
p7_refmx_Destroy(P7_REFMX *rmx)
{
  if (!rmx) return;
  if (rmx->dp_mem) free(rmx->dp_mem);
  if (rmx->dp)     free(rmx->dp);
  free(rmx);
}


/*****************************************************************
 * 2. Debugging aids
 *****************************************************************/

/* Function:  p7_refmx_Compare()
 * Synopsis:  Compare two DP matrices for equality (within given tolerance)
 *
 * Purpose:   Compare all the values in DP matrices <rx1>, <rx2> for
 *            equality within absolute epsilon <tolerance>, using
 *            <esl_FCompareAbs()> calls. Return <eslOK> if all cell
 *            comparisons succeed; <eslFAIL> if not.
 *            
 *            Absolute difference comparison is preferred over
 *            relative differences. Numerical error accumulation in DP
 *            scales more with the number of terms than their
 *            magnitude. DP cells with values close to zero (and hence
 *            small absolute differences) may reasonably have large
 *            relative differences.
 */
int
p7_refmx_Compare(const P7_REFMX *rx1, const P7_REFMX *rx2, float tolerance)
{
  int i,k,s;
  int killmenow = FALSE;
#ifdef p7_DEBUGGING
  killmenow = TRUE;	/* Setting <killmenow> during debugging helps find where a comparison fails */
#endif /*p7_DEBUGGING*/
  
  if (rx1->M    != rx2->M)    { if (killmenow) abort(); return eslFAIL; }
  if (rx1->L    != rx2->L)    { if (killmenow) abort(); return eslFAIL; }
  if (rx1->type != rx2->type) { if (killmenow) abort(); return eslFAIL; }
  
  for (i = 0; i <= rx1->L; i++)
    {
      for (k = 0; k <= rx1->M; k++)   
	for (s = 0; s < p7R_NSCELLS; s++)
	  if ( esl_FCompareAbs(P7R_MX(rx1,i,k,s), P7R_MX(rx2,i,k,s), tolerance) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
      for (s = 0; s < p7R_NXCELLS; s++)
	if ( esl_FCompareAbs(P7R_XMX(rx1,i,s), P7R_XMX(rx2,i,s), tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
    }
  return eslOK;	
}

/* Function:  p7_refmx_CompareLocal()
 * Synopsis:  Compare two DP matrices (local paths only) for equality 
 *
 * Purpose:   A variant of <p7_refmx_Compare()> that compares only 
 *            cells on local paths. <MG,IG,DG,G> states are excluded.
 *            
 *            This gets used in unit tests that compare Backwards
 *            matrices computed by the ForwardFilter() with the
 *            reference implementation. You can't expect these Bck
 *            matrices to compare completely equal.  In the fwdfilter,
 *            glocal paths are all -inf by construction (they're not
 *            evaluated at all).  In the reference impl, in local
 *            mode, backwards values are finite along glocal paths
 *            (G/MG/IG/DG) because the zero transition that prohibits
 *            the path is the B->G transition, which isn't evaluated
 *            in the recursion until the *end* of glocal paths. 
 */
int
p7_refmx_CompareLocal(const P7_REFMX *rx1, const P7_REFMX *rx2, float tolerance)
{
  int i,k;
  int killmenow = FALSE;
#ifdef p7_DEBUGGING
  killmenow = TRUE;	/* Setting <killmenow> during debugging helps find where a comparison fails */
#endif /*p7_DEBUGGING*/
  
  if (rx1->M    != rx2->M)    { if (killmenow) abort(); return eslFAIL; }
  if (rx1->L    != rx2->L)    { if (killmenow) abort(); return eslFAIL; }
  if (rx1->type != rx2->type) { if (killmenow) abort(); return eslFAIL; }
  
  for (i = 0; i <= rx1->L; i++)
    {
      for (k = 0; k <= rx1->M; k++)   
	{
	  if ( esl_FCompareAbs(P7R_MX(rx1,i,k,p7R_ML), P7R_MX(rx2,i,k,p7R_ML), tolerance) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	  if ( esl_FCompareAbs(P7R_MX(rx1,i,k,p7R_IL), P7R_MX(rx2,i,k,p7R_IL), tolerance) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	  if ( esl_FCompareAbs(P7R_MX(rx1,i,k,p7R_DL), P7R_MX(rx2,i,k,p7R_DL), tolerance) == eslFAIL) { if (killmenow) abort(); return eslFAIL; }
	}
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_E),  P7R_XMX(rx2,i,p7R_E),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_N),  P7R_XMX(rx2,i,p7R_N),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_J),  P7R_XMX(rx2,i,p7R_J),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_B),  P7R_XMX(rx2,i,p7R_B),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_L),  P7R_XMX(rx2,i,p7R_L),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_C),  P7R_XMX(rx2,i,p7R_C),  tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_JJ), P7R_XMX(rx2,i,p7R_JJ), tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
      if ( esl_FCompareAbs(P7R_XMX(rx1,i,p7R_CC), P7R_XMX(rx2,i,p7R_CC), tolerance) == eslFAIL)   { if (killmenow) abort(); return eslFAIL; }
    }
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
p7_refmx_DecodeSpecial(int type)
{
  switch (type) {
  case p7R_E:  return "E";
  case p7R_N:  return "N";
  case p7R_J:  return "J";
  case p7R_B:  return "B";
  case p7R_L:  return "L";
  case p7R_G:  return "G";
  case p7R_C:  return "C";
  case p7R_JJ: return "JJ";
  case p7R_CC: return "CC";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_REFMX special state code %d\n", type);
  return NULL;
}

/* Function:  p7_refmx_DecodeState()
 * Synopsis:  Convert main state code to string for debugging output.
 *
 * Purpose:   Given a main state code (such as <p7R_MG>), return a
 *            string label, suitable as a state label in debugging
 *            output.
 *
 * Args:      type - main state code, such as <p7R_MG>
 *
 * Returns:   ptr to static string representation.
 */
char *
p7_refmx_DecodeState(int type)
{
  switch (type) {
  case p7R_ML: return "ML";
  case p7R_MG: return "MG";
  case p7R_IL: return "IL";
  case p7R_IG: return "IG";
  case p7R_DL: return "DL";
  case p7R_DG: return "DG";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such P7_REFMX main state code %d\n", type);
  return NULL;
}
    

/* Function:  p7_refmx_Dump()
 * Synopsis:  Dump a <P7_REFMX> for examination.
 *
 * Purpose:   Print the contents of <rmx> to stream <ofp>, for
 *            examination/debugging.
 */
int
p7_refmx_Dump(FILE *ofp, P7_REFMX *rmx)
{
  return p7_refmx_DumpWindow(ofp, rmx, 0, rmx->L, 0, rmx->M);
}

int
p7_refmx_DumpWindow(FILE *ofp, P7_REFMX *rmx, int istart, int iend, int kstart, int kend)
{
  int   width     = 9;
  int   precision = 4;
  int   i,k,x;

  /* Header */
  fprintf(ofp, "       ");
  for (k = kstart; k <= kend;       k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0;      x < p7R_NXCELLS; x++) fprintf(ofp, "%*s ", width, p7_refmx_DecodeSpecial(x));
  fprintf(ofp, "\n");

  fprintf(ofp, "       ");
  for (k = kstart; k <= kend;  k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < p7R_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");

   /* DP matrix data */
  for (i = istart; i <= iend; i++)
    {
      fprintf(ofp, "%3d ML ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_ML]);
      for (x = 0;  x <  p7R_NXCELLS; x++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][ (rmx->M+1) * p7R_NSCELLS + x]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IL ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_IL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DL ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_DL]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d MG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_MG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d IG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_IG]);
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d DG ", i);
      for (k = kstart; k <= kend;    k++)  fprintf(ofp, "%*.*f ", width, precision, rmx->dp[i][k * p7R_NSCELLS + p7R_DG]);
      fprintf(ofp, "\n\n");
  }
  return eslOK;
}


/* a private hack for making heatmaps */
int
p7_refmx_DumpCSV(FILE *fp, P7_REFMX *pp, int istart, int iend, int kstart, int kend)
{
  int   width     = 7;
  int   precision = 5;
  int   i, k;
  float val;

  fprintf(fp, "i,");
  for (k = kend; k >= kstart; k--)
    fprintf(fp, "%-d%s", k, k==kstart ? "\n" : ",");

  for (i = istart; i <= iend; i++)
    {
      fprintf(fp, "%-d,", i);
      for (k = kend; k >= kstart; k--)
	{
	  val = 
	    pp->dp[i][k * p7R_NSCELLS + p7R_ML] + pp->dp[i][k * p7R_NSCELLS + p7R_MG] + 
 	    pp->dp[i][k * p7R_NSCELLS + p7R_IL] + pp->dp[i][k * p7R_NSCELLS + p7R_IG] + 
	    pp->dp[i][k * p7R_NSCELLS + p7R_DL] + pp->dp[i][k * p7R_NSCELLS + p7R_DG];

	  fprintf(fp, "%*.*f%s", width, precision, val, k==kstart ? "\n" : ", ");
	}
    }
  return eslOK;
}
/*------------- end, (most) debugging tools ---------------------*/


/*****************************************************************
 * 3. Validation of matrices.
 *****************************************************************/
/* Gets its own section because it's a little fiddly and detailed.
 * See p7_refmx.h for notes on the patterns that get tested.
 */

static inline int
validate_dimensions(P7_REFMX *rmx, char *errbuf)
{
  if ( rmx->M <= 0)               ESL_FAIL(eslFAIL, errbuf, "nonpositive M");
  if ( rmx->L <= 0)               ESL_FAIL(eslFAIL, errbuf, "nonpositive L");
  if ( rmx->L+1    > rmx->validR) ESL_FAIL(eslFAIL, errbuf, "L+1 is larger than validR");
  if ( rmx->validR > rmx->allocR) ESL_FAIL(eslFAIL, errbuf, "validR larger than allocR");
  if ( (rmx->M+1)*p7R_NSCELLS+p7R_NXCELLS < rmx->allocW) ESL_FAIL(eslFAIL, errbuf, "M is too large for allocW");
  return eslOK;
}

static inline int 
validate_mainstate(P7_REFMX *rmx, int i, int k, int s, float val, char *errbuf)
{
  if (P7R_MX(rmx, i, k, s) != val) ESL_FAIL(eslFAIL, errbuf, "expected %f at i=%d, k=%d, %s", val, i, k, p7_refmx_DecodeState(s));
  return eslOK;
}

static inline int
validate_special(P7_REFMX *rmx, int i, int s, float val, char *errbuf)
{
  if (P7R_XMX(rmx, i, s) != val) ESL_FAIL(eslFAIL, errbuf, "expected %f at i=%d, %s", val, i, p7_refmx_DecodeSpecial(s));
  return eslOK;
}

static inline int
validate_no_nan(P7_REFMX *rmx, char *errbuf)
{
  int i,k,s;
  for (i = 0; i <= rmx->L; i++)
    for (k = 0; k <= rmx->M; k++)
      {
	for (s = 0; s < p7R_NSCELLS; s++)
	  if (isnan(P7R_MX(rmx, i, k, s))) ESL_FAIL(eslFAIL, errbuf, "found NaN at i=%d, k=%d, %s", i, k, p7_refmx_DecodeState(s));      
	for (s = 0; s < p7R_NXCELLS; s++)
	  if (isnan(P7R_XMX(rmx, i, s)))   ESL_FAIL(eslFAIL, errbuf, "found NaN at i=%d, %s", i, p7_refmx_DecodeSpecial(s));      
      }
  return eslOK;
}

static inline int
validate_column_zero(P7_REFMX *rmx, float val, char *errbuf)
{
  int i, s, status;
  for (i = 0; i <= rmx->L; i++)
    for (s = 0; s < p7R_NSCELLS; s++)
      if (( status = validate_mainstate(rmx, i, 0, s, val, errbuf)) != eslOK) return status;
  return eslOK;
}


static inline int
validate_forward(P7_REFMX *rmx, char *errbuf)
{
  int i,k,s;
  int status;

  if (( status = validate_column_zero(rmx, -eslINFINITY, errbuf)) != eslOK) return status;

  /* Row i=0 */
  for (k = 1; k <= rmx->M; k++)
    for (s = 0; s <= p7R_NSCELLS; s++)
      if ( (status = validate_mainstate(rmx, 0, k, s, -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_E,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_N,     0.0f,         errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_J,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_C,     -eslINFINITY, errbuf)) != eslOK) return status;

  /* Row i=1 has some additional expected -inf's, different from the remaining rows 2..L */
  for (k = 1; k <= rmx->M; k++) {
    if ( (status = validate_mainstate(rmx, 1, k, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
    if ( (status = validate_mainstate(rmx, 1, k, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
  }

  /* Rows 2..L, plus the row i=1 cells we didn't just check */
  for (i = 1; i <= rmx->L; i++)
    {
      if ((status =    validate_mainstate(rmx, i,      1, p7R_DL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i,      1, p7R_DG, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i, rmx->M, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i, rmx->M, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
    }

  /* CC/JJ are only for decoding; Forward matrix has them as -inf */
  for (i = 0; i <= rmx->L; i++)
    {
      if ( ( status = validate_special(rmx, 0, p7R_JJ,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(rmx, 0, p7R_CC,    -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}

static inline int
validate_backward(P7_REFMX *rmx, char *errbuf)
{
  int i,k,s;
  int status;

  if (( status = validate_column_zero(rmx, -eslINFINITY, errbuf)) != eslOK) return status;  

  /* Row i=0 */
  for (k = 1; k <= rmx->M; k++)
    for (s = 0; s < p7R_NSCELLS; s++)
      if ( (status = validate_mainstate(rmx, 0, k, s, -eslINFINITY, errbuf)) != eslOK) return status;

  /* Rows 1..L */
  for (i = 1; i <= rmx->L; i++)
    {
      if ((status = validate_mainstate(rmx, i, rmx->M, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status = validate_mainstate(rmx, i, rmx->M, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
    }

  /* Row i=L has some additional expected -inf's, different from the remaining rows */
  for (k = 1; k <= rmx->M; k++) {
    if ( (status = validate_mainstate(rmx, rmx->L, k, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
    if ( (status = validate_mainstate(rmx, rmx->L, k, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
  }
  if ( (status = validate_special  (rmx, rmx->L,    p7R_N,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_J,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_B,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_L,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_G,  -eslINFINITY, errbuf)) != eslOK) return status;

  /* CC/JJ are only for decoding; Backward matrix has them as -inf */
  for (i = 0; i <= rmx->L; i++)
    {
      if ( ( status = validate_special(rmx, 0, p7R_JJ,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(rmx, 0, p7R_CC,    -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}

static inline int
validate_decoding(P7_REFMX *rmx, char *errbuf)
{
  int i,k,s;
  int status;

  /* All of k=0 column is 0.0 */
  if (( status = validate_column_zero(rmx, 0.0, errbuf)) != eslOK) return status;  

  /* Row i=0: main states all 0.0  */
  for (k = 1; k <= rmx->M; k++)
    for (s = 0; s <= p7R_NSCELLS; s++)
      if ( (status = validate_mainstate(rmx, 0, k, s, 0.0f, errbuf)) != eslOK) return status;
  /*          and E,J,C are unreached */
  if ( ( status = validate_special(rmx, 0, p7R_E, 0.0f, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_J, 0.0f, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_C, 0.0f, errbuf)) != eslOK) return status;

  /* Rows 1..L */
  for (i = 1; i <= rmx->L; i++)
    {
      if ((status = validate_mainstate(rmx, i,      1, p7R_DL, 0.0f, errbuf)) != eslOK) return status;
      if ((status = validate_mainstate(rmx, i,      1, p7R_DG, 0.0f, errbuf)) != eslOK) return status;
      if ((status = validate_mainstate(rmx, i, rmx->M, p7R_IL, 0.0f, errbuf)) != eslOK) return status;
      if ((status = validate_mainstate(rmx, i, rmx->M, p7R_IG, 0.0f, errbuf)) != eslOK) return status;
    }

  /* Rows i=1,L have some additional zeros */
  for (k = 1; k < rmx->M; k++) {
    if ((status = validate_mainstate(rmx,      1, k, p7R_IL, 0.0f, errbuf)) != eslOK) return status;
    if ((status = validate_mainstate(rmx,      1, k, p7R_IG, 0.0f, errbuf)) != eslOK) return status;
    if ((status = validate_mainstate(rmx, rmx->L, k, p7R_IL, 0.0f, errbuf)) != eslOK) return status;
    if ((status = validate_mainstate(rmx, rmx->L, k, p7R_IG, 0.0f, errbuf)) != eslOK) return status;
  }
  /* on row 1, JJ and CC cannot be reached */
  if ((status = validate_special  (rmx,      1, p7R_JJ, 0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx,      1, p7R_CC, 0.0f, errbuf)) != eslOK) return status;
  /* on row L, only ...E->C(CC)->T specials can be reached. */
  if ((status = validate_special  (rmx, rmx->L, p7R_N,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx, rmx->L, p7R_J,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx, rmx->L, p7R_B,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx, rmx->L, p7R_L,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx, rmx->L, p7R_G,  0.0f, errbuf)) != eslOK) return status;
  if ((status = validate_special  (rmx, rmx->L, p7R_JJ, 0.0f, errbuf)) != eslOK) return status;
  return eslOK;
}

static inline int
validate_alignment(P7_REFMX *rmx, char *errbuf)
{
  int i,k,s;
  int status;

  if (( status = validate_column_zero(rmx, -eslINFINITY, errbuf)) != eslOK) return status;  

  /* Row i=0 */
  for (k = 1; k <= rmx->M; k++)
    for (s = 0; s <= p7R_NSCELLS; s++)
      if ( (status = validate_mainstate(rmx, 0, k, s, -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_E,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_J,     -eslINFINITY, errbuf)) != eslOK) return status;
  if ( ( status = validate_special(rmx, 0, p7R_C,     -eslINFINITY, errbuf)) != eslOK) return status;

  /* Row i=1 has some additional expected -inf's, different from the remaining rows 2..L */
  for (k = 1; k <= rmx->M; k++) {
    if ( (status = validate_mainstate(rmx, 1, k, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
    if ( (status = validate_mainstate(rmx, 1, k, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
  }

  /* Rows 2..L, plus the row i=1 cells we didn't just check */
  for (i = 1; i <= rmx->L; i++)
    {
      if ((status =    validate_mainstate(rmx, i,      1, p7R_DL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i,      1, p7R_DG, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i, rmx->M, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
      if ((status =    validate_mainstate(rmx, i, rmx->M, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
    }

  /* Row i=L has some additional expected -inf's, different from the remaining rows */
  for (k = 1; k <= rmx->M; k++) {
    if ( (status = validate_mainstate(rmx, rmx->L, k, p7R_IL, -eslINFINITY, errbuf)) != eslOK) return status;
    if ( (status = validate_mainstate(rmx, rmx->L, k, p7R_IG, -eslINFINITY, errbuf)) != eslOK) return status;
  }
  if ( (status = validate_special  (rmx, rmx->L,    p7R_N,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_J,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_B,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_L,  -eslINFINITY, errbuf)) != eslOK) return status;
  if ( (status = validate_special  (rmx, rmx->L,    p7R_G,  -eslINFINITY, errbuf)) != eslOK) return status;

  /* and throughout, CC/JJ are only for decoding; Alignment matrix has them as -inf */
  for (i = 0; i <= rmx->L; i++)
    {
      if ( ( status = validate_special(rmx, 0, p7R_JJ,    -eslINFINITY, errbuf)) != eslOK) return status;
      if ( ( status = validate_special(rmx, 0, p7R_CC,    -eslINFINITY, errbuf)) != eslOK) return status;
    }

  return eslOK;
}


/* Function:  p7_refmx_Validate()
 * Synopsis:  Validates a reference DP matrix.
 *
 * Purpose:   Validates the internals of the
 *            DP matrix <rmx>. Returns <eslOK> if
 *            it passes. Returns <eslFAIL> if it fails,
 *            and sets <errbuf> to contain an explanation,
 *            if caller provided an <errbuf>
 *
 * Args:      rmx    - reference DP matrix to validate.
 *            errbuf - allocated space for err msg, or NULL
 *
 * Returns:   <eslOK> on success; <eslFAIL> on failure, and
 *            <errbuf> (if provided) contains explanation.
 *
 * Xref:      See notes in <p7_refmx.h> for an explanation
 *            of the particular patterns we're checking.
 */
int
p7_refmx_Validate(P7_REFMX *rmx, char *errbuf)
{
  int status;

  if ( (status = validate_dimensions(rmx, errbuf)) != eslOK) return status;
  if ( (status = validate_no_nan    (rmx, errbuf)) != eslOK) return status;
  
  switch (rmx->type) {
  case p7R_FORWARD:   if ( (status = validate_forward    (rmx,               errbuf)) != eslOK) return status;  break;
  case p7R_VITERBI:   if ( (status = validate_forward    (rmx,               errbuf)) != eslOK) return status;  break; // Viterbi has same pattern as Forward
  case p7R_BACKWARD:  if ( (status = validate_backward   (rmx,               errbuf)) != eslOK) return status;  break;
  case p7R_DECODING:  if ( (status = validate_decoding   (rmx,               errbuf)) != eslOK) return status;  break;     
  case p7R_ALIGNMENT: if ( (status = validate_alignment  (rmx,               errbuf)) != eslOK) return status;  break; 
  case p7R_UNSET:     if ( (status = validate_column_zero(rmx, -eslINFINITY, errbuf)) != eslOK) return status;  break;
  default:            ESL_FAIL(eslFAIL, errbuf, "no such reference DP algorithm type %d", rmx->type);
  }
  return eslOK;
}
/*------------------ end, validation ----------------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
