/* H4_SPARSEMX: DP matrix for production implementation
 *
 * Contents:
 *   1. H4_SPARSEMX object
 *   2. Debugging tools for H4_SPARSEMX
 *   3. Validation of H4_SPARSEMX
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "h4_path.h"
#include "h4_refmx.h"
#include "h4_sparsemask.h"
#include "h4_sparsemx.h"

/*****************************************************************
 * 1. H4_SPARSEMX object
 *****************************************************************/

/* Function:  h4_sparsemx_Create()
 * Synopsis:  Create a new sparse DP matrix.
 *
 * Purpose:   Create a new sparse DP matrix, defined by the sparse DP
 *            mask <sm>.
 *            
 *            If <sm> is <NULL>, allocate for a reasonable-sized
 *            matrix. Matrix won't be usable until caller calls
 *            <h4_sparsemx_Reinit()>. This style works well in loops
 *            over sequences, where you can create the object once
 *            before any sequences are known, then
 *            <_Reinit()/_Reuse()> it in the loop.
 */
H4_SPARSEMX *
h4_sparsemx_Create(const H4_SPARSEMASK *sm)
{
  H4_SPARSEMX *sx            = NULL;
  int64_t      default_ncell = 4096; /* if <sm> is NULL, what <dalloc> should be (stored sparsified cells) */
  int          default_nx    = 256;  /* if <sm> is NULL, what <xalloc> should be (stored rows of specials) */
  int          status;

  ESL_ALLOC(sx, sizeof(H4_SPARSEMX));
  sx->dp   = NULL;
  sx->xmx  = NULL;
  sx->sm   = sm;
  sx->type = h4S_UNSET;

  /* We must avoid zero-sized mallocs. If there are no rows or cells, alloc the default sizes */
  sx->dalloc = ( (sm && sm->ncells)         ? sm->ncells       : default_ncell);
  sx->xalloc = ( (sm && (sm->nrow + sm->S)) ? sm->nrow + sm->S : default_nx);

  ESL_ALLOC(sx->dp,  sizeof(float) * h4S_NSCELLS * sx->dalloc);
  ESL_ALLOC(sx->xmx, sizeof(float) * h4S_NXCELLS * sx->xalloc);
  return sx;

 ERROR:
  h4_sparsemx_Destroy(sx);
  return NULL;
}

/* Function:  h4_sparsemx_Reinit()
 * Synopsis:  Reinitialize, reallocate sparse DP matrix for new calculation.
 *
 * Purpose:   Reinitialize an existing sparse matrix <sx> to use the
 *            sparse mask <sm>. Equivalent to a call to
 *            <h4_sparsemx_Create()> except we reuse as much memory as
 *            possible in the preexisting <sx>, to minimize
 *            realloc/free calls.
 *
 * Returns:   <eslOK> on success, and <sx> is now ready for sparse DP
 *            calculations defined by the mask <sm>.
 *
 * Throws:    <eslEMEM> on (re-)allocation failure.
 */
int
h4_sparsemx_Reinit(H4_SPARSEMX *sx, const H4_SPARSEMASK *sm)
{
  int64_t dalloc_req = sm->ncells;
  int     xalloc_req = sm->nrow + sm->S;
  int     status;

  if (sx->dalloc < dalloc_req) {
    ESL_REALLOC(sx->dp, sizeof(float) * h4S_NSCELLS * dalloc_req);
    sx->dalloc = dalloc_req;
  }
  if (sx->xalloc < xalloc_req) {
    ESL_REALLOC(sx->xmx, sizeof(float) * h4S_NXCELLS * xalloc_req);
    sx->xalloc = xalloc_req;
  }
  sx->sm   = sm;
  sx->type = h4S_UNSET;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  h4_sparsemx_SetValues()
 * Synopsis:  Sets all cells in sparse matrix to zero.
 */
int
h4_sparsemx_SetValues(H4_SPARSEMX *sx, float val)
{
  esl_vec_FSet(sx->dp,   sx->sm->ncells*h4S_NSCELLS,          val);
  esl_vec_FSet(sx->xmx, (sx->sm->nrow+sx->sm->S)*h4S_NXCELLS, val);
  return eslOK;
}


/* Function:  h4_sparsemx_Sizeof()
 * Synopsis:  Returns current allocated size of a H4_SPARSEMX, in bytes.
 *
 * Purpose:   Returns allocated size of <sx>, in bytes.  Does not
 *            include the size of its <H4_SPARSEMASK>.
 */
size_t
h4_sparsemx_Sizeof(const H4_SPARSEMX *sx)
{
  size_t n = sizeof(H4_SPARSEMX);
  n += sizeof(float) * h4S_NSCELLS * sx->dalloc;
  n += sizeof(float) * h4S_NXCELLS * sx->xalloc;
  return n;
}

/* Function:  h4_sparsemx_MinSizeof()
 * Synopsis:  Returns minimal required allocation size of a H4_SPARSEMX, in bytes.
 *
 * Purpose:   Calculate and return the minimum required size, in 
 *            bytes, of a sparse DP matrix, for a profile/sequence
 *            comparison using the sparse mask in <sm>.
 *            
 *            Taking a <H4_SPARSEMASK> as the arg, not <H4_SPARSEMX>,
 *            is not a typo. Does not require having an actual DP
 *            matrix allocated.  We use this function when
 *            planning/profiling memory allocation strategies.
 */
size_t
h4_sparsemx_MinSizeof(const H4_SPARSEMASK *sm)
{
  size_t n = sizeof(H4_SPARSEMX);
  n += sizeof(float) * h4S_NSCELLS * sm->ncells;          // dp[]
  n += sizeof(float) * h4S_NXCELLS * (sm->nrow + sm->S);  // xmx[]; for each seg ia..ib, ia-1..ib has special cells
  return n;
}

/* Function:  h4_sparsemx_Destroy()
 * Synopsis:  Free a sparse DP matrix.
 */
void
h4_sparsemx_Destroy(H4_SPARSEMX *sx)
{
  if (sx) {
    free(sx->dp);
    free(sx->xmx);
    /* sx->sm is a reference copy. caller remains responsible for it. */
    free(sx);
  }
}
/*----------- end, H4_SPARSEMX implementation -------------------*/



/*****************************************************************
 * 2. Debugging tools
 *****************************************************************/

char *
h4_sparsemx_DecodeState(int type)
{
  switch (type) {
  case h4S_ML: return "ML";
  case h4S_MG: return "MG";
  case h4S_IL: return "IL";
  case h4S_IG: return "IG";
  case h4S_DL: return "DL";
  case h4S_DG: return "DG";
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such H4_SPARSEMX main state code %d\n", type);
  return NULL;
}

char *
h4_sparsemx_DecodeSpecial(int type)
{
  switch (type) {
  case h4S_E:  return "E";
  case h4S_N:  return "N";
  case h4S_J:  return "J";
  case h4S_B:  return "B";
  case h4S_L:  return "L";
  case h4S_G:  return "G";
  case h4S_C:  return "C";
  case h4S_JJ: return "JJ";
  case h4S_CC: return "CC";
  default:     break;
  }
  esl_exception(eslEINVAL, FALSE, __FILE__, __LINE__, "no such H4_SPARSEMX special state type code %d\n", type);
  return NULL;
}


int
h4_sparsemx_Dump(FILE *ofp, H4_SPARSEMX *sx)
{
  return h4_sparsemx_DumpWindow(ofp, sx, 0, sx->sm->L, 0, sx->sm->M);
}

int
h4_sparsemx_DumpWindow(FILE *ofp, const H4_SPARSEMX *sx, int ia, int ib, int ka, int kb)
{
  const H4_SPARSEMASK *sm  = sx->sm;
  float         *dpc = sx->dp;
  float         *xc  = sx->xmx;
  int width          = 9;
  int precision      = 4;
  int i,k,x,z;

  /* Header */
  fprintf(ofp, "       ");
  for (k = ka; k <= kb;         k++) fprintf(ofp, "%*d ", width, k);
  for (x = 0;  x < h4S_NXCELLS; x++) fprintf(ofp, "%*s ", width, h4_sparsemx_DecodeSpecial(x));
  fprintf(ofp, "\n");

  fprintf(ofp, "       ");
  for (k = ka; k <= kb;        k++) fprintf(ofp, "%*.*s ", width, width, "----------");
  for (x = 0; x < h4S_NXCELLS; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");

  /* Skipping ahead in matrix, over rows we're not dumping: */
  for (i = 1; i < ia; i++) 
    if (sm->n[i]) {
      if (sm->n[i-1] == 0) xc += h4S_NXCELLS; // skip an extra chunk of specials (ia-1) before each segment start on ia
      dpc += sm->n[i] * h4S_NSCELLS;          // skip over rows we're not dumping
      xc  += h4S_NXCELLS;                     // skip specials on sparsified rows
    }

  for (i = ia; i <= ib; i++)
    {
      /* If current row has no cells...  */
      if (sm->n[i] == 0) {
        if (i > 1     && sm->n[i-1] > 0) fputs("...\n\n", ofp);  // ... if last row did, then we ended a segment. Print an ellipsis.
        if (i < sm->L && sm->n[i+1] > 0)                         // if next row does, then we're about to start a segment; print specials on an ia-1 row 
          {
            fprintf(ofp, "%3d -- ", i);
            for (k = ka; k <= kb;         k++) fprintf(ofp, "%*s ", width, ".....");
            for (x = 0;  x < h4S_NXCELLS; x++) fprintf(ofp, "%*.*f ", width, precision, xc[x]);
            fputs("\n\n", ofp);
            xc += h4S_NXCELLS;    
          }
        continue;                       
      }

      fprintf(ofp, "%3d ML ", i);
      for (z = 0, k = ka; k <= kb; k++) { 
        while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
        if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*h4S_NSCELLS + h4S_ML));
        else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      for (x = 0; x < h4S_NXCELLS; x++) fprintf(ofp, "%*.*f ", width, precision, xc[x]);
      fputc('\n', ofp);

      fprintf(ofp, "%3d MG ", i);
      for (z = 0, k = ka; k <= kb; k++) { 
        while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
        if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*h4S_NSCELLS + h4S_MG));
        else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d IL ", i);
      for (z = 0, k = ka; k <= kb; k++) { 
        while (z < sm->n[i] && sm->k[i][z] < k)  z++;
        if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*h4S_NSCELLS + h4S_IL));
        else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d IG ", i);
      for (z = 0, k = ka; k <= kb; k++) { 
        while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
        if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*h4S_NSCELLS + h4S_IG));
        else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d DL ", i);
      for (z = 0, k = ka; k <= kb; k++) { 
        while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
        if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*h4S_NSCELLS + h4S_DL));
        else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputc('\n', ofp);

      fprintf(ofp, "%3d DG ", i);
      for (z = 0, k = ka; k <= kb; k++) { 
        while (z < sm->n[i] && sm->k[i][z] < k)  z++; 
        if    (z < sm->n[i] && sm->k[i][z] == k) fprintf(ofp, "%*.*f ", width, precision,  *(dpc + z*h4S_NSCELLS + h4S_DG));
        else                                     fprintf(ofp, "%*s ",   width, ".....");
      }
      fputs("\n\n", ofp);

      dpc += sm->n[i] * h4S_NSCELLS;
      xc  += h4S_NXCELLS;
    }
  return eslOK;
}

/* Function:  h4_sparsemx_CompareReference()
 * Synopsis:  Test sparse DP matrix for tolerable equality to reference
 *
 * Purpose:   Compare all the values in sparse DP matrix <sx> to the
 *            corresponding values in reference DP matrix <rx> for
 *            equality within the absolute epsilon <tol>, using
 *            <esl_FCompare()> calls. Return <eslOK> if comparison
 *            succeeds; return <eslFAIL> otherwise.
 *            
 *            In general, this is only going to succeed if <sx> was
 *            calculated as a 'full' sparse matrix, with all cells
 *            marked, as in <h4_sparsemask_AddAll()>. For truly sparse
 *            calculations, you will more likely want to use
 *            <h4_sparsemx_CompareReferenceAsBound()>, which treats
 *            the reference calculation's matrix values as an upper
 *            bound on the corresponding sparse cell values.
 *            
 *            The comparison tests two different traversal mechanisms
 *            for sparse main cells.  Should be completely redundant,
 *            unless we've corrupted the structure.  The routine is
 *            intended for debugging, not speed.
 *            
 *            As also noted in <h4_refmx_Compare()>, absolute
 *            difference comparison is preferred over relative
 *            differences. Numerical error accumulation in DP scales
 *            more with the number of terms than their magnitude. DP
 *            cells with values close to zero (and hence small
 *            absolute differences) may reasonably have large relative
 *            differences.
 *            
 *            When debugging, you can set a breakpoint at <esl_fail()> 
 *            to stop the debugger where a comparison failed.
 *            
 *            We assume that <h4S_NSCELLS> and <h4S_NXCELLS> main and
 *            special elements are in the same order and number in
 *            <h4R_NSCELLS> and <h4R_NXCELLS>.
 *
 * Args:      sx  - sparse DP matrix
 *            rx  - reference DP matrix
 *            tol - absolute floating point comparison tolerance
 *                 
 * Returns:   <eslOK> on success.
 *            <eslFAIL> if comparison fails.
 */
int
h4_sparsemx_CompareReference(const H4_SPARSEMX *sx, const H4_REFMX *rx, float tol)
{
  char                 msg[] = "comparison of H4_SPARSEMX to H4_REFMX failed";
  const H4_SPARSEMASK *sm    = sx->sm;
  const float   *dpc, *dpc2;
  const float   *xc,  *xc2;
  int            g,i,s,z;
  int            ia,ib;
  
  if (sx->type != rx->type) ESL_FAIL(eslFAIL, NULL, msg);
  if (sm->M    != rx->M)    ESL_FAIL(eslFAIL, NULL, msg);
  if (sm->L    != rx->L)    ESL_FAIL(eslFAIL, NULL, msg);

  /* First traversal way: sm->dp[] is just one big array */
  dpc = sx->dp;
  xc  = sx->xmx;
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1])         /* ia-1 specials at a segment start */
        {
          for (s = 0; s < h4S_NXCELLS; xc++, s++) 
            if (esl_FCompare(*xc, H4R_XMX(rx,i-1,s), /*rtol=*/0.0, tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
        }
      for (z = 0; z < sm->n[i]; z++)       /* sparse cells */
        {
          for (s = 0; s < h4S_NSCELLS; dpc++, s++) 
            if (esl_FCompare(*dpc, H4R_MX(rx,i,sm->k[i][z],s), /*rtol=*/0.0, tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
        }
      if (sm->n[i])       /* specials */
        {
          for (s = 0; s < h4S_NXCELLS; xc++, s++) 
            if (esl_FCompare(*xc, H4R_XMX(rx,i,s), /*rtol=*/0.0, tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
        }
    }

  /* Second way: "segments" */
  dpc2 = sx->dp;
  xc2  = sx->xmx;
  for (g = 1; g <= sm->S; g++)
    {
      ia = sm->seg[g].ia;
      ib = sm->seg[g].ib;

      for (s = 0; s < h4S_NXCELLS; xc2++, s++)        /* ia-1 specials at segment start */
        if (esl_FCompare(*xc2, H4R_XMX(rx,ia-1,s), /*rtol=*/0.0, tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);

      for (i = ia; i <= ib; i++) 
        {
          for (z = 0; z < sm->n[i]; z++)          /* sparse main cells */
            {
              for (s = 0; s < h4S_NSCELLS; dpc2++, s++) 
                if (esl_FCompare(*dpc2, H4R_MX(rx,i,sm->k[i][z],s), /*rtol=*/0.0, tol) == eslFAIL)  ESL_FAIL(eslFAIL, NULL, msg);
            }
          for (s = 0; s < h4S_NXCELLS; xc2++, s++)        /* specials */
            if (esl_FCompare(*xc2, H4R_XMX(rx,i,s), /*rtol=*/0.0, tol) == eslFAIL) ESL_FAIL(eslFAIL, NULL, msg);
        }
    }
  
  /* Both ways must reach the same end */
  if (dpc != dpc2) ESL_FAIL(eslFAIL, NULL, msg);
  if (xc  != xc2)  ESL_FAIL(eslFAIL, NULL, msg);
  return eslOK;
}

/* Function:  h4_sparsemx_CompareReferenceAsBound()
 * Synopsis:  Test sparse DP matrix cell values are upper-bounded by reference
 *
 * Purpose:   Check that all values in sparse matrix <sx> are bounded
 *            (less than or equal to) the corresponding values in reference
 *            matrix <rx>, allowing for a small amount <tol> of floating-point
 *            roundoff error accumulation. That is, <v_s <= v_r+tol> for all
 *            matrix cell values <v>. Return <eslOK> if comparison succeeds;
 *            return <eslFAIL> otherwise.
 *            
 *            This is a variant of <p7_sparsemx_CompareReference()>;
 *            see additional documentation there, with the exception
 *            that this routine only uses one traversal mechanism in
 *            the sparse matrix.
 *
 * Args:      sx  - sparse DP matrix
 *            rx  - reference DP matrix
 *            tol - absolute floating point comparison tolerance
 *                 
 * Returns:   <eslOK> on success.
 *            <eslFAIL> if comparison fails.
 */
int
h4_sparsemx_CompareReferenceAsBound(const H4_SPARSEMX *sx, const H4_REFMX *rx, float tol)
{
  char                 msg[] = "failed comparison of H4_SPARSEMX to upper bound in H4_REFMX"; 
  const H4_SPARSEMASK *sm    = sx->sm;
  const float         *dpc   = sx->dp;
  const float         *xc    = sx->xmx;
  int   i,s,z;

  if (sx->type != rx->type) ESL_FAIL(eslFAIL, NULL, msg);
  if (sm->M    != rx->M)    ESL_FAIL(eslFAIL, NULL, msg);
  if (sm->L    != rx->L)    ESL_FAIL(eslFAIL, NULL, msg);
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {    /* ia-1 specials at a segment start */
        for (s = 0; s < h4S_NXCELLS; xc++, s++) 
          if (*xc > H4R_XMX(rx,i-1,s)+tol) 
            ESL_FAIL(eslFAIL, NULL, msg);
      }
      for (z = 0; z < sm->n[i]; z++)    /* sparse cells */
        for (s = 0; s < h4S_NSCELLS; dpc++, s++) 
          if (*dpc > H4R_MX(rx,i,sm->k[i][z],s)+tol) 
            ESL_FAIL(eslFAIL, NULL, msg);
  
      if (sm->n[i]) {    /* specials */
        for (s = 0; s < h4S_NXCELLS; xc++, s++) 
          if (*xc > H4R_XMX(rx,i,s)+tol) 
            ESL_FAIL(eslFAIL, NULL, msg);
      }
    }
  return eslOK;
}


/* Function:  h4_sparsemx_CompareDecoding()
 * Synopsis:  Compare exact and approximate posterior decoding matrices.
 *
 * Purpose:   Compare exact sparse decoding matrix <sxe> (calculated by 
 *            <h4_sparse_Decoding()> to an approximate one <sxa> 
 *            (calculated by sampling lots of stochastic traces).
 *            Make sure that no nonzero value in <sxe> differs by
 *            more than absolute difference <tol> in <sxa>, and make
 *            sure that an exact zero value in <sxe> is also zero
 *            in <sxa>. Return <eslOK> if these tests succeed; <eslFAIL>
 *            if they do not.
 *            
 *            This comparison is used in the main unit test of
 *            posterior decoding. See
 *            <sparse_dp.c::utest_approx_decoding()>.
 *            
 * Args:      sxe  - exact posterior decoding matrix, from h4_sparse_Decoding
 *            sxa  - approximate decoding mx, from stochastic trace ensemble
 *            tol  - absolute difference to tolerate per cell value
 *
 * Returns:   <eslOK> if all comparisons pass. <eslFAIL> if not.
 */
int
h4_sparsemx_CompareDecoding(const H4_SPARSEMX *sxe, const H4_SPARSEMX *sxa, float tol)
{
  char                 msg[] = "failed comparison of exact and sampled sparse decoding matrices";
  const H4_SPARSEMASK *sm  = sxe->sm;
  const float         *dpe = sxe->dp;
  const float         *dpa = sxa->dp;
  const float         *xce = sxe->xmx;
  const float         *xca = sxa->xmx;
  int   g,i,s,z;

  if (sxe->type != h4S_DECODING) ESL_FAIL(eslFAIL, NULL, msg);
  if (sxa->type != h4S_DECODING) ESL_FAIL(eslFAIL, NULL, msg);
  if (sm->M     != sxa->sm->M)   ESL_FAIL(eslFAIL, NULL, msg);
  if (sm->L     != sxa->sm->L)   ESL_FAIL(eslFAIL, NULL, msg);

  /* main sparse cells */
  for (g = 1; g <= sm->S; g++)
    for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)    // main rows from ia..ib
      for (z = 0; z < sm->n[i]; z++) 
        for (s = 0; s < h4S_NSCELLS; dpe++, dpa++, s++)
          {
            if (*dpe == 0.0) { if (*dpa != 0.0) ESL_FAIL(eslFAIL, NULL, msg); }
            else             { if (esl_FCompare(*dpe, *dpa, 0., tol) != eslOK) ESL_FAIL(eslFAIL, NULL, msg); }
          }

  /* specials */
  for (g = 1; g <= sm->S; g++)
    for (i = sm->seg[g].ia-1; i <= sm->seg[g].ib; i++)  // specials run ia-1..ib
      for (s = 0; s < h4S_NXCELLS; xce++, xca++, s++)
        {
          if (*xce == 0.0) { if (*xca != 0.0) ESL_FAIL(eslFAIL, NULL, msg); }
          else             { if (esl_FCompare(*xce, *xca, 0., tol) != eslOK) ESL_FAIL(eslFAIL, NULL, msg); }
        }         
  return eslOK;
}


/* Function:  h4_sparsemx_CountPath()
 * Synopsis:  Count one (stochastic) path into accumulating decoding matrix.
 *
 * Purpose:   Used when we're approximating posterior decoding by path
 *            sampling.
 */
int
h4_sparsemx_CountPath(const H4_PATH *pi, H4_SPARSEMX *sxd)
{
  static int sttype[h4P_NST] = { -1, -1, h4S_N, h4S_B, h4S_G, h4S_MG, h4S_IG, h4S_DG, h4S_L, h4S_ML, h4S_IL, h4S_DL, h4S_E, h4S_J, h4S_C, -1 }; /* sttype[] translates trace idx to DP matrix idx*/
  const H4_SPARSEMASK *sm  = sxd->sm;
  float               *dpc = sxd->dp;   // ptr that steps thru stored main supercells i,k
  float               *xc  = sxd->xmx;  // ptr that steps thru stored special rows, including ia-1 seg edges
  int   i = 0;                          // current position in sequence
  int   k = 0;                          // current position in profile
  int   y = 0;                          // index in sparse cell list on a row
  int   z;                              // index in trace positions
  int   r;                              // index in path runs

  if (sxd->type == h4S_UNSET) sxd->type = h4S_DECODING; // first time 
  ESL_DASSERT1(( sxd->type == h4S_DECODING ));
  
  for (z = 0; z < pi->Z; z++)
    {
      if (h4_path_IsX(pi->st[z])) // N|J|C
        {
          if (pi->st[z] == h4P_J || pi->st[z] == h4P_C) xc[h4S_E] += 1.0;  // X->E->{JC}. We know i is stored, for first J|C
          if (sm->n[i] || (i < sm->L && sm->n[i+1]))    xc[sttype[(int) pi->st[z]]] += 1.0;  //  ... but not necessarily for first N at i=0.

          //          printf("trying to count a %s at i=%d\n", h4_path_DecodeStatetype(pi->st[z]), i);

          for (r = 1; r < pi->rle[z]; r++)   // any remainder from r=1 are NN|CC|JJ emissions that advance i
            {
              while (y < sm->n[i]) { y++; dpc += h4S_NSCELLS; } // skip remainder of prv sparse row
              if (sm->n[i] || sm->n[i+1])  xc += h4S_NXCELLS;   // advance specials to new row i
              i++;                                              // advance i
              y = 0;                                            

              //printf("trying to count a %s at i=%d\n", h4_path_DecodeStatetype(pi->st[z]), i);

              if (sm->n[i] || (i < sm->L && sm->n[i+1])) {
                xc[sttype[(int) pi->st[z]]] += 1.0;
                if      (pi->st[z] == h4P_J) xc[h4S_JJ] += 1.0;
                else if (pi->st[z] == h4P_C) xc[h4S_CC] += 1.0;
              }
            }
        }
      else if (h4_path_IsB(pi->st[z])) // G|L
        {
          if (sm->n[i] || sm->n[i+1]) { xc[h4S_B] += 1.0; xc[sttype[(int) pi->st[z]]] += 1.0; }
          else ESL_EXCEPTION(eslEINCONCEIVABLE, "in sparse_CountPath, a G|L state must be a stored row");

          //          printf("trying to count a %s at i=%d\n", h4_path_DecodeStatetype(pi->st[z]), i);

          k    = pi->rle[z]-1;  // k=0 for G, or local entry-1 for L
          dpc -= h4S_NSCELLS*y; // back up to start of current row i, in case of G->D1..Dk-> entry
          y    = 0;
        }
      else // MG|IG|DG|ML|IL|DL
        {
          for (r = 0; r < pi->rle[z]; r++)
            {
              if (h4_path_IsM(pi->st[z]) || h4_path_IsI(pi->st[z]))  // do we need to advance i?
                {
                  while (y < sm->n[i]) { y++; dpc += h4S_NSCELLS; }  // skip remainder of prev sparse row 
                  if (sm->n[i] || sm->n[i+1])  xc += h4S_NXCELLS;    // increment specials. note, i+1 is safe here; i<L at this point, because tr->i[z] <= L and i < tr->i[z] 
                  i++;
                  y = 0;
                }

              if (h4_path_IsM(pi->st[z]) || h4_path_IsD(pi->st[z])) k++; // do we need to advance k?

              //              printf("trying to count a %s at i=%d k=%d\n", h4_path_DecodeStatetype(pi->st[z]), i,k);

              while (y < sm->n[i] && sm->k[i][y]  < k) { y++; dpc += h4S_NSCELLS; }           // try to find sparse cell for i,k
              if    (y < sm->n[i] && sm->k[i][y] == k) dpc[sttype[(int) pi->st[z]]] += 1.0;   // did we find it? then increment +1
              else if (pi->st[z] != h4P_DG)
                ESL_EXCEPTION(eslEINCONCEIVABLE, "in sparse_CountPath, M|I|D states (except wing-retracted DG's) must be stored");
            }
        }
    } // end loop over states z in path

  return eslOK;
}
/*------------------ end, debugging tools -----------------------*/



/*****************************************************************
 * 7. Validation of a H4_SPARSEMX
 *****************************************************************/
/* also a debugging tool, but in its own section because it's
 * fiddly and complicated
 */

static int
validate_dimensions(const H4_SPARSEMX *sx, char *errbuf)
{
  const H4_SPARSEMASK *sm  = sx->sm;
  int   g      = 0;
  int   r      = 0;
  int   ncells = 0;
  int   i;

  if ( sm->M <= 0)                  ESL_FAIL(eslFAIL, errbuf, "nonpositive M");
  if ( sm->L <= 0)                  ESL_FAIL(eslFAIL, errbuf, "nonpositive L");
  if ( sm->V <= 0)                  ESL_FAIL(eslFAIL, errbuf, "nonpositive V");
  if ( sm->Q != H4_Q(sm->M,sm->V))  ESL_FAIL(eslFAIL, errbuf, "bad Q");          

  for (r=0, g=0, i = 1; i <= sm->L; i++) {
    if (sm->n[i] && !sm->n[i-1]) g++; /* segment count */
    if (sm->n[i])                r++; /* sparse row count */
    ncells += sm->n[i];
  }
  if (g      != sm->S)          ESL_FAIL(eslFAIL, errbuf, "S (nseg) is wrong");
  if (r      != sm->nrow)       ESL_FAIL(eslFAIL, errbuf, "nrow is wrong");
  if (ncells != sm->ncells)     ESL_FAIL(eslFAIL, errbuf, "ncells is wrong");

  if (sm->L+1    > sm->ralloc)  ESL_FAIL(eslFAIL, errbuf, "k[] row allocation too small");
  if (sm->ncells > sm->kalloc)  ESL_FAIL(eslFAIL, errbuf, "kmem[] cell allocation too small");
  if (sm->S+2    > sm->salloc)  ESL_FAIL(eslFAIL, errbuf, "seg[] segment allocation too small");
  return eslOK;
}

static int
validate_no_nan(const H4_SPARSEMX *sx, char *errbuf)
{
  const H4_SPARSEMASK *sm  = sx->sm;
  float         *dpc = sx->dp;
  float         *xc  = sx->xmx;
  int            i,k,z,s;

  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1])          /* ia-1 specials at a segment start */
        {
          for (s = 0; s < h4S_NXCELLS; s++) {
            if (isnan(*xc)) ESL_FAIL(eslFAIL, errbuf, "nan at i=%d, %s", i, h4_sparsemx_DecodeSpecial(s));
            xc++;
          }
        }
      for (z = 0; z < sm->n[i]; z++)       /* sparse main cells */
        {
          k = sm->k[i][z];
          for (s = 0; s < h4S_NSCELLS; s++) {
            if (isnan(*dpc)) ESL_FAIL(eslFAIL, errbuf, "nan at i=%d, k=%d, %s", i, k, h4_sparsemx_DecodeState(s));
            dpc++;
          }
        }

      if (sm->n[i])                       /* specials on sparse row */
        {
          for (s = 0; s < h4S_NXCELLS; s++) {
            if (isnan(*xc)) ESL_FAIL(eslFAIL, errbuf, "nan at i=%d, %s", i, h4_sparsemx_DecodeSpecial(s));
            xc++;
          }
        }
    }
  return eslOK;
}

static int
validate_fwdvit(const H4_SPARSEMX *sx, char *errbuf)
{
  const H4_SPARSEMASK *sm  = sx->sm;
  float         *dpc = sx->dp;
  float         *xc  = sx->xmx;
  int            i,z;

  /* Check special cases prohibited in the first ia-1 presegment specials: */
  if ( xc[h4S_J] !=  -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "first J not -inf");
  if ( xc[h4S_C] !=  -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "first C not -inf");
  
  /* Sweep, checking for (the most easily spotchecked) prohibited values (must be -inf) */
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {       /* ia-1 specials at a segment start */
        if ( xc[h4S_E]  != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "E seg start for ia=%d not -inf", i);
        if ( xc[h4S_JJ] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "JJ seg start for ia=%d not -inf", i);
        if ( xc[h4S_CC] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "CC seg start for ia=%d not -inf", i);
        xc += h4S_NXCELLS;
      }
      for (z = 0; z < sm->n[i]; z++)       /* sparse main cells */
        {
          /* if k-1 supercell doesn't exist, can't reach D's */
          if ((z == 0 || sm->k[i][z] != sm->k[i][z-1]+1) && dpc[h4S_DL] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "first DL on i=%d not -inf", i);
          if ((z == 0 || sm->k[i][z] != sm->k[i][z-1]+1) && dpc[h4S_DG] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "first DG on i=%d not -inf", i);
          if (   sm->k[i][z] == sm->M                    && dpc[h4S_IL] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "IL on i=%d,k=M not -inf", i);
          if (   sm->k[i][z] == sm->M                    && dpc[h4S_IG] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "IG on i=%d,k=M not -inf", i);
          dpc += h4S_NSCELLS;
          /* there are other conditions where I(i,k) values must be zero but this is more tedious to check */
        }
      if (sm->n[i]) {                     
        if ( xc[h4S_JJ] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "JJ at i=%d not -inf", i);
        if ( xc[h4S_CC] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "CC at i=%d not -inf", i);
        xc += h4S_NXCELLS;
      }
    }
  return eslOK;
}

static int
validate_backward(const H4_SPARSEMX *sx, char *errbuf)
{
  const H4_SPARSEMASK *sm = sx->sm;
  float         *dpc    = sx->dp  + (sm->ncells-1)*h4S_NSCELLS;         // last supercell in dp  
  float         *xc     = sx->xmx + (sm->nrow + sm->S - 1)*h4S_NXCELLS; // last supercell in xmx 
  int            last_n = 0;
  int            i,z;

  /* Backward sweep; many of our prohibits are on ib segment-end rows */
  /* first: some special cases on absolute final stored row ib */
  if (xc[h4S_N] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "N on last row not 0");
  if (xc[h4S_J] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "J on last row not 0");
  /* sweep: */
  for (i = sm->L; i >= 1; i--)
    {
      if (sm->n[i]) {                   /* specials on stored row i */
        if (               xc[h4S_JJ] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "JJ on row i=%d not -inf", i);
        if (               xc[h4S_CC] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "CC on row i=%d not -inf", i);
        if (last_n == 0 && xc[h4S_B]  != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "B on end-seg row ib=%d not -inf", i);
        if (last_n == 0 && xc[h4S_L]  != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "L on end-seg row ib=%d not -inf", i);
        if (last_n == 0 && xc[h4S_G]  != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "G on end-seg row ib=%d not -inf", i);
        xc -= h4S_NXCELLS;
      }
      for (z = sm->n[i]-1; z >= 0; z--) /* sparse main cells */
        {
          if (sm->k[i][z] == sm->M && dpc[h4S_IL] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "IL on i=%d,k=M not -inf", i);
          if (sm->k[i][z] == sm->M && dpc[h4S_IG] != -eslINFINITY) ESL_FAIL(eslFAIL, errbuf, "IG on i=%d,k=M not -inf", i);
          dpc -= h4S_NSCELLS;
        }
      if (sm->n[i] && sm->n[i-1] == 0) xc -= h4S_NXCELLS; /* specials on ia-1 row before a start-segment */
      last_n = sm->n[i];
    }
  return eslOK;
}

static int
is_prob(float val, float tol)
{
  if (val < 0. || val > 1.+tol) return FALSE; 
  return TRUE;
}

static int
validate_decoding(const H4_SPARSEMX *sx, char *errbuf)
{
  const H4_SPARSEMASK *sm  = sx->sm;
  const float         *dpc = sx->dp;
  const float         *xc;
  int            i,z,s;
  float          rowsum;
  float          tol = 0.001;

  /* Check special cases prohibited in the first ia-1 presegment specials: */
  /* Note that the N on the first ia-1 row is *not* necessarily 1, because ia-1 is not necessarily row 0 */
  xc = sx->xmx; // first supercell in xmx
  if ( esl_FCompare(xc[h4S_J],  0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first J not 0");
  if ( esl_FCompare(xc[h4S_C],  0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first C not 0");
  if ( esl_FCompare(xc[h4S_JJ], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first JJ not 0");
  if ( esl_FCompare(xc[h4S_CC], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first CC not 0");

  /* Check special cases in very last ib row specials */
  xc = sx->xmx + (sm->nrow + sm->S - 1)*h4S_NXCELLS; // last supercell in xmx 
  if (esl_FCompare(xc[h4S_N], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "N on last row not 0");
  if (esl_FCompare(xc[h4S_J], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "J on last row not 0");
  if (esl_FCompare(xc[h4S_B], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "B on last row not 0");
  if (esl_FCompare(xc[h4S_L], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "L on last row not 0");
  if (esl_FCompare(xc[h4S_G], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "G on last row not 0");

  /* Sweep, checking for (the most easily spotchecked) prohibited values (must be 0.0's) */
  xc = sx->xmx;
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {       /* ia-1 specials at a segment start */
        if ( esl_FCompare(xc[h4S_E], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "E seg start for ia=%d not 0", i);
        if ( xc[h4S_J]+tol < xc[h4S_JJ])                               ESL_FAIL(eslFAIL, errbuf, "JJ>J at seg start for ia=%d ", i);
        if ( xc[h4S_C]+tol < xc[h4S_CC])                               ESL_FAIL(eslFAIL, errbuf, "CC>C at seg start for ia=%d ", i);
        xc += h4S_NXCELLS;
      }
      for (z = 0; z < sm->n[i]; z++)       /* sparse main cells */
        {
          /* if k-1 supercell doesn't exist, can't reach DL. But all DGk are reachable, because of wing-retracted entry/exit */
          if ((z == 0 || sm->k[i][z] != sm->k[i][z-1]+1) && esl_FCompare(dpc[h4S_DL], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "first DL on i=%d not 0", i);
          if (   sm->k[i][z] == sm->M                    && esl_FCompare(dpc[h4S_IL], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "IL on i=%d,M not 0", i);
          if (   sm->k[i][z] == sm->M                    && esl_FCompare(dpc[h4S_IG], 0.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "IG on i=%d,M not 0", i);
          dpc += h4S_NSCELLS;
          /* there are other conditions where I(i,k) values must be zero but this is more tedious to check */
        }
      if (sm->n[i]) {
        if ( xc[h4S_J]+tol < xc[h4S_JJ]) ESL_FAIL(eslFAIL, errbuf, "JJ>J at i=%d ", i);
        if ( xc[h4S_C]+tol < xc[h4S_CC]) ESL_FAIL(eslFAIL, errbuf, "CC>C at i=%d ", i);
        xc += h4S_NXCELLS;
      }
    }

  /* Sweep again; check all values are probabilities, 0<=x<=1.
   */
  dpc = sx->dp;
  xc  = sx->xmx;
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) {       /* ia-1 specials at a segment start */
        for (s = 0; s < h4S_NXCELLS; xc++, s++) 
          if (! is_prob(*xc, tol)) ESL_FAIL(eslFAIL, errbuf, "bad decode prob %f for %s, seg start, ia=%d\n", *xc, h4_sparsemx_DecodeSpecial(s), i); 
      }
      for (z = 0; z < sm->n[i]; z++)       /* sparse main cells */
        for (s = 0; s < h4S_NSCELLS; dpc++, s++) 
          if (! is_prob(*dpc, tol)) ESL_FAIL(eslFAIL, errbuf, "bad decode prob %f at i=%d,k=%d,%s", *dpc, i,sm->k[i][z], h4_sparsemx_DecodeState(s));
      if (sm->n[i]) {                      /* specials on sparse row */
        for (s = 0; s < h4S_NXCELLS; xc++, s++) 
          if (! is_prob(*xc, tol))  ESL_FAIL(eslFAIL, errbuf, "bad decode prob %f at i=%d,%s", *xc, i, h4_sparsemx_DecodeSpecial(s)); 
      }
    }

  /* Sweep again; each residue i must be emitted by something,
   * so the sum of all emitting states across a row is 1.0.
   */
  dpc = sx->dp;
  xc  = sx->xmx;
  for (i = 1; i <= sm->L; i++)
    {
      if (sm->n[i] && !sm->n[i-1]) /* ia-1 specials at a segment start: check xc only (on row ia-1) before doing row i. */
        {
          rowsum = xc[h4S_N] + xc[h4S_JJ] + xc[h4S_CC]; /* this'd even work for nonemitting row ia-1=0, because there xc[N]=1 and JJ/CC=0 */
          if (esl_FCompare(rowsum, 1.0, /*rtol=*/0.0, tol) != eslOK) ESL_FAIL(eslFAIL, errbuf, "sparse decode mx init-seg row %d sums to %.4f not 1", i, rowsum);
          xc += h4S_NXCELLS;
        }

      if (sm->n[i])
        {
          for (rowsum = 0., z = 0; z < sm->n[i]; z++, dpc += h4S_NSCELLS) // sparse main cells 
            rowsum += dpc[h4S_ML] + dpc[h4S_MG] + dpc[h4S_IL] + dpc[h4S_IG];
          rowsum += xc[h4S_N] + xc[h4S_JJ] + xc[h4S_CC];
          if (esl_FCompare(rowsum, 1.0, /*rtol=*/0.0, tol) != eslOK)  ESL_FAIL(eslFAIL, errbuf, "sparse decode mx row %d sums to %.4f not 1", i, rowsum);
          xc += h4S_NXCELLS;
        }
    }

  return eslOK;
}


/* Function:  h4_sparsemx_Validate()
 * Synopsis:  Validate a sparse DP matrix.
 *
 * Purpose:   Validate the contents of sparse DP matrix <sx>.
 *            Return <eslOK> if it passes. Return <eslFAIL> if
 *            it fails, and set <errbuf> to contain an 
 *            explanation, if caller provides a non-<NULL>
 *            <errbuf>.
 *
 * Args:      sx      - sparse DP matrix to validate
 *            errbuf  - char[eslERRBUFSIZE] space for error msg, or NULL.
 *
 * Returns:   <eslOK> on success.
 *            <eslFAIL> on failure, with an error message in <errbuf>
 *            if <errbuf> was provided.
 */
int
h4_sparsemx_Validate(const H4_SPARSEMX *sx, char *errbuf)
{
  int status;

  if (errbuf) errbuf[0] = '\0';

  if ( (status = validate_dimensions(sx, errbuf)) != eslOK) return status;
  if ( (status = validate_no_nan    (sx, errbuf)) != eslOK) return status;

  switch (sx->type) {
  case h4S_UNSET:      ESL_FAIL(eslFAIL, errbuf, "validating an unset sparse DP matrix? probably not what you meant");
  case h4S_FORWARD:    if ( (status = validate_fwdvit  (sx, errbuf)) != eslOK) return status; break;
  case h4S_BACKWARD:   if ( (status = validate_backward(sx, errbuf)) != eslOK) return status; break;
  case h4S_DECODING:   if ( (status = validate_decoding(sx, errbuf)) != eslOK) return status; break;
  case h4S_VITERBI:    if ( (status = validate_fwdvit  (sx, errbuf)) != eslOK) return status; break;
  default:             ESL_FAIL(eslFAIL, errbuf, "no such sparse DP matrix type %d", sx->type);
  }
  return eslOK;
}
/*----------------- end, H4_SPARSEMX validation -----------------*/



