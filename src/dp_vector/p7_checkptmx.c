/* P7_CHECKPTMX: checkpointed, striped vector DP matrix.
 * 
 * Independent of vector ISA. (Do not add any ISA-specific code.)
 * See p7_checkptmx.md for implementation notes.
 * 
 * Contents:
 *    1. The P7_CHECKPTMX object
 *    2. Debugging, development routines.
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alloc.h"

#include "dp_vector/simdvec.h"
#include "dp_vector/p7_checkptmx.h"


/*****************************************************************
 * 1. The <P7_CHECKPTMX> object
 *****************************************************************/

/* minimum_rows()
 * 
 * Calculate the minimum number of rows needed to checkpoint a forward
 * matrix for a sequence of length <L>, exclusive of other constant
 * row overhead (R0, consisting of two backwards rows and the fwd[0]
 * initial row).
 * 
 * This value is a double. If it has a fraction, a partial checkpoint
 * block ("b" region) is needed, as in this typical code:
 *    Rbc  = minimum_rows(L);
 *    Rc   = floor(Rbc);
 *    Rb   = (Rbc > Rc ? 1 : 0);
 *    minR = (int) ceil(Rbc);    // or, Rc+Rb
 */
double 
minimum_rows(int L)
{
  return (sqrt(9. + 8. * (double) L) - 3.) / 2.;
}


/* checkpointed_rows()
 * 
 * Given L and maxR, solve for the number of checkpointed rows (Rb+Rc)
 * we need. The value is a double, and if it has a fractional part,
 * then a partial checkpoint block is needed, Rb==1.
 * 
 * This equation is obtained by solving 
 *     L = Ra + (Rbc +2)(Rbc+1)/2 - 1
 * for Rbc (i.e. Rb+Rc) using the quadratic equation,
 * after substitution Ra = (maxR - Rbc - R0) to get
 * Rbc in terms of L,maxR.
 */
double
checkpointed_rows(int L, int R)
{
  return (sqrt(1. + 8. * (double) (L - R)) - 1.) / 2.;
}


/* A "full" matrix is easy: Ra = La = L, using Ra+R0 <= maxR rows total. */
static void
set_full(P7_CHECKPTMX *ox, int L)
{
  ox->Ra     = L;
  ox->Rb     = 0;
  ox->Rc     = 0;
  ox->La     = L;
  ox->Lb     = 0;
  ox->Lc     = 0;
}

/* If we can fit a checkpointed matrix into maxR rows, 
 * then the trick is to use all maxR rows, making the
 * "a" region (all rows kept) as large as possible, to
 * minimize computation. This means solving a little
 * quadratic equation for Rb+Rc given L and maxR: see
 * <checkpointed_rows()> for that solution.
 */
void
set_checkpointed(P7_CHECKPTMX *ox, int L, int R)
{
  double Rbc = checkpointed_rows(L, R - ox->R0);
  double Rc  = floor(Rbc);

  ox->Rc     = (int) Rc;
  ox->Rb     = (Rbc > Rc ? 1 : 0);
  ox->Ra     = R - ox->Rb - ox->Rc - ox->R0;
  ox->Lc     = ((ox->Rc + 2) * (ox->Rc + 1)) / 2 - 1;
  ox->La     = ox->Ra;
  ox->Lb     = L - ox->La - ox->Lc;
}   

/* If we can't fit in maxR rows, then we checkpoint
 * the entire matrix; R0+Ra+Rb+Rc > maxR.
 */
void
set_redlined(P7_CHECKPTMX *ox, int L, double minR)
{
  double Rc = floor(minR);

  ox->Rc     = (int) Rc;
  ox->Rb     = (minR > Rc ? 1 : 0);
  ox->Ra     = 0;
  ox->Lc     = ((ox->Rc + 2) * (ox->Rc + 1)) / 2 - 1;
  ox->La     = 0;
  ox->Lb     = L - ox->La - ox->Lc;
}

/* set_row_layout()
 *
 * Determines the size of the a,b,c regions ("all", "between",
 * "checkpointed") of rows in the DP matrix. 
 *
 * Caller must have already set <ox->allocW> and <ox->R0>; they are
 * needed here.
 * 
 * Upon return, we've set the R{0abc} and L{abc} fields in the <ox>
 * structure.
 * 
 * <maxR> is the maximum number of rows the caller *wants* to use,
 * either because a <ramlimit>'ed allocation fits that number of rows,
 * or because an existing matrix has that number of valid rows.  We
 * will exceed this for one comparison if absolutely necessary, but
 * the next <_Reuse()> call will bring the allocation back down.
 * 
 * So there's three possibilities:
 *  1. A full matrix fits into our recommended max memory use.
 *  2. A checkpointed matrix fits into our recommended memory.
 *     Make as much of the matrix uncheckpointed as we can,
 *     using every row in maxR.
 *  3. We can't satisfy the recommended max memory, even fully
 *     checkpointed. Make a fully checkpointed matrix, in which
 *     R0+Ra+Rb+Rc will exceed maxR, and caller will have to 
 *     allocate ("redlined").
 */
static void
set_row_layout(P7_CHECKPTMX *ox, int allocL, int maxR)
{
  double Rbc      = minimum_rows(allocL);               
  int    minR_chk = ox->R0 + (int) ceil(Rbc);	          /* min # rows we need for checkpointing          */
  int    minR_all = ox->R0 + allocL;		          /* min # rows we need for full matrix            */

  if      (minR_all <= maxR) set_full        (ox, allocL);
  else if (minR_chk <= maxR) set_checkpointed(ox, allocL, maxR);
  else                       set_redlined    (ox, allocL, Rbc);
}




/* Function:  p7_checkptmx_Create()
 * Synopsis:  Allocate a new <P7_CHECKPTMX> object.
 *
 * Purpose:   Allocate a new <P7_CHECKPTMX> checkpointed, striped vector
 *            DP matrix sufficient for the Forward/Backward local
 *            decoding calculation for a query model
 *            of up to length <M> and a target sequence of up to
 *            length <L>.
 *            
 *            Try to keep the allocation within <redline> bytes in
 *            memory.  For example, <redline=ESL_MBYTES(128)>, sets a
 *            recommended memory limit of 128 MiB. Allocation can
 *            exceed this, if even a fully checkpointed <MxL>
 *            comparison requires it -- but in this case, any
 *            subsequent <p7_checkptmx_Reinit()> call that attempts to
 *            reuse the matrix will try to reallocate it back
 *            downwards to the <redline>.
 *            
 *            Choice of <redline> should take into account how many
 *            parallel threads there are, because each one will likely
 *            have its own <P7_CHECKPTMX> allocation.
 *            
 *            By design spec, <M> and <L> are $\leq$ 100K.
 *
 * Args:      M        - query profile size, consensus positions (<=100000)
 *            L        - target sequence length, residues (<=100000)
 *            ramlimit - recommended memory limit, bytes
 *
 * Returns:   ptr to new <P7_CHECKPTMX> object on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_CHECKPTMX *
p7_checkptmx_Create(int M, int L, int64_t redline)
{
  P7_CHECKPTMX *ox = NULL;
  int           maxR;
  int           r;
  int           status;
  
  /* Validity of integer variable ranges may depend on design spec:                  */
  ESL_DASSERT1( (M <= 100000) );       /* design spec says  model length M <= 100000 */
  ESL_DASSERT1( (L <= 100000) );       /*           ... and   seq length L <= 100000 */
  ESL_DASSERT1( (L >  0) );
  ESL_DASSERT1( (M >  0) );

  /* Level 1 allocation: the structure itself */
  ESL_ALLOC(ox, sizeof(P7_CHECKPTMX));
  ox->dp_mem  = NULL;
  ox->dpf     = NULL;

  /* Set checkpointed row layout: allocR, R{abc}, L{abc} fields */
  ox->R0          = 3;	  // fwd[0]; bck[prv,cur] 

  /* Calculate allocW, the allocated width (in cells, floats) of each vector row */
  /*                       maxQf        *  3 MDI cells * floats/vector    */
  ox->allocW      = P7_Q(M, p7_VMAX_FB) * p7C_NSCELLS  * p7_VMAX_FB;        
  ox->allocW     += ESL_UPROUND(p7C_NXCELLS, p7_VMAX_FB);    // specials are scalar. maintain vector mem alignment
  ox->redline     = redline / sizeof(float);                 // redline, allocW, allocN are all in <float> units, internally.                             

  maxR            = ox->redline / ox->allocW;
  set_row_layout(ox, L, maxR);
  ox->allocR      = ox->R0 + ox->Ra + ox->Rb + ox->Rc;
  ox->validR      = ox->allocR;

  ESL_DASSERT1(( ox->allocW % p7_VMAX_FB == 0 )); /* verify that concat'ed rows will stay aligned */

  /* Level 2 allocations: row pointers and dp cell memory */
  ox->allocN = ox->allocR * ox->allocW; 
  ox->dp_mem = esl_alloc_aligned(sizeof(float) * ox->allocN, p7_VALIGN);  // dp_mem is aligned vector memory
  if (ox->dp_mem == NULL) { status = eslEMEM; goto ERROR; }

  ESL_ALLOC( ox->dpf, sizeof(float **) * ox->allocR);    
  for (r = 0; r < ox->validR; r++)
    ox->dpf[r] = ox->dp_mem + (r * ox->allocW);  
  
#if eslDEBUGLEVEL > 0
  ox->do_dumping     = FALSE;
  ox->dfp            = NULL;
  ox->dump_maxpfx    = 5;	
  ox->dump_width     = 9;
  ox->dump_precision = 4;
  ox->dump_flags     = p7_DEFAULT;
  ox->fwd            = NULL;
  ox->bck            = NULL;
  ox->pp             = NULL;
  ox->bcksc          = 0.0f;
#endif

  ox->M  = 0;
  ox->L  = 0;
  ox->R  = 0; 
  ox->Q  = 0;
  ox->Vf = 0;  // until striped data are valid. Gets set by vector DP code.
  return ox;

 ERROR:
  p7_checkptmx_Destroy(ox);
  return NULL;
}

/* Function:  p7_checkptmx_Reinit()
 * Synopsis:  Reinitialize checkpointed DP matrix for a new comparison.
 *
 * Purpose:   Given an existing checkpointed matrix structure <ox>,
 *            and the dimensions <M> and <L> of a new comparison,
 *            reinitialize <ox> for that comparison, reallocating
 *            if needed.
 *
 *            Essentially equivalent to _Create(), but while reusing
 *            previous space and minimizing expensive reallocations.
 *            Should only be called in a context similar to a
 *            _Create() call. For efficiency, old data are not copied
 *            to new internal allocations, so any existing data may be
 *            lost.
 *
 *            Usually <ox> only grows. However, if <ox> is redlined
 *            (over its recommended allocation) and the new problem
 *            size <M,L> can fit in the preset recommended allocation,
 *            then <ox> is reallocated down to the smaller recommended
 *            size.
 *            
 * Args:      ox    - existing checkpointed matrix
 *            M     - new query profile length
 *            L     - new target sequence length         
 * 
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if an allocation fails. The state of <ox> is
 *            now undefined, and the caller should not use it. 
 */
int
p7_checkptmx_Reinit(P7_CHECKPTMX *ox, int M, int L)
{
  int     minR_chk      = (int) ceil(minimum_rows(L)) + ox->R0; // minimum number of DP rows needed  
  int     reset_dp_ptrs = FALSE;
  int     maxR;
  int64_t W;			/* minimum row width needed, bytes */
  int     r;
  int     status;

  /* Validity of integer variable ranges may depend on design spec:                  */
  ESL_DASSERT1( (M <= 100000) );       /* design spec says, model length M <= 100000 */
  ESL_DASSERT1( (L <= 100000) );       /*           ... and,  seq length L <= 100000 */
  ESL_DASSERT1( (L >  0) );
  ESL_DASSERT1( (M >  0) );

  /* Reset data-specific stuff for a new comparison */
  ox->M  = 0;
  ox->L  = 0;
  ox->R  = 0;
  ox->Q  = 0;
  ox->Vf = 0;

  /* If we're debugging and we have stored copies of any matrices,
   * reset and realloc them too.  Must do this first, because we have an early 
   * successful exit condition coming below.
   *  (SRE TODO: these should be _Reinit() calls, combining Reuse and GrowTo)
   */
#if eslDEBUGLEVEL > 0
  if (ox->fwd && (status = p7_refmx_Reuse(ox->fwd)) != eslOK) return status;
  if (ox->bck && (status = p7_refmx_Reuse(ox->bck)) != eslOK) return status;
  if (ox->pp  && (status = p7_refmx_Reuse(ox->pp))  != eslOK) return status;
  ox->bcksc = 0.0f;
#endif
#if eslDEBUGLEVEL > 0
  if (ox->fwd && (status = p7_refmx_GrowTo(ox->fwd, M, L)) != eslOK) goto ERROR;
  if (ox->bck && (status = p7_refmx_GrowTo(ox->bck, M, L)) != eslOK) goto ERROR;
  if (ox->pp  && (status = p7_refmx_GrowTo(ox->pp,  M, L)) != eslOK) goto ERROR;
#endif
 
  /* Calculate W, minimum row width needed by NEW allocation, in floats (cells) */
  /*         maxQf        *  3 MDI cells * cells/vector                         */
  W  = P7_Q(M,p7_VMAX_FB) * p7C_NSCELLS  * p7_VMAX_FB;       // vector part, aligned for any vector ISA
  W += ESL_UPROUND(p7C_NXCELLS, p7_VMAX_FB);                 // specials are scalar. maintain vector mem alignment

  /* Are the current allocations satisfactory? */
  if (W <= ox->allocW && ox->allocN <= ox->redline)
    {
      if      (L + ox->R0 <= ox->validR) { set_full        (ox, L);             return eslOK; }
      else if (minR_chk   <= ox->validR) { set_checkpointed(ox, L, ox->validR); return eslOK; }
    }

  /* Do individual matrix rows need to be expanded, by resetting ox->dpf ptrs into dp_mem? */
  if ( W > ox->allocW) 
    {
      ox->allocW    = W;
      ox->validR    = ox->allocN / ox->allocW; // validR is still guaranteed to be <= allocR after this
      reset_dp_ptrs = TRUE;
    }

  /* Does matrix dp_mem need reallocation, either up or down? */
  maxR  = ox->allocN / ox->allocW;                              // max rows if we use up to the current allocation size.      
  if ( (ox->allocN > ox->redline && minR_chk <= maxR) ||        // we were redlined, and current alloc will work: so downsize 
       minR_chk > ox->validR)				        // not enough memory for needed rows: so upsize                   
    {
      set_row_layout(ox, L, maxR); 
      ox->validR = ox->R0 + ox->Ra + ox->Rb + ox->Rc;      // this may be > allocR now; we'll reallocate dp[] next, if so 
      ox->allocN = ox->validR * ox->allocW;

      esl_alloc_free(ox->dp_mem);
      ox->dp_mem = esl_alloc_aligned(sizeof(float) * ox->allocN, p7_VALIGN);
      if (ox->dp_mem == NULL) { status = eslEMEM; goto ERROR; }

      reset_dp_ptrs = TRUE;
    }
  else  // current validR will suffice, either full or checkpointed; but we still need to calculate a layout 
    {
      if   (L+ox->R0 <= ox->validR) set_full(ox, L); 
      else                          set_checkpointed(ox, L, ox->validR);
    }
  
  /* Does the array of row ptrs need reallocation? */
  if (ox->validR > ox->allocR)
    {
      ESL_REALLOC(ox->dpf, sizeof(float *) * ox->validR);
      ox->allocR    = ox->validR;
      reset_dp_ptrs = TRUE;
    }

  /* Do the row ptrs need to be reset? */
  if (reset_dp_ptrs)
    {
      for (r = 0; r < ox->validR; r++)
	ox->dpf[r] = ox->dp_mem + (r * ox->allocW);
    }

  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_checkptmx_Sizeof()
 * Synopsis:  Returns size of checkpointed vector DP matrix, in bytes.
 * 
 * Purpose:   Returns the size of the checkpointed vector DP matrix
 *            in bytes. 
 *            
 *            If code has been compiled in debugging mode, the
 *            returned size includes a negligible amount of extra
 *            space for debugging fields in the structure (about 5
 *            ints, 4 pointers, and a float - around 56 bytes). The
 *            returned size does not include the use of any full
 *            Forward, Backward, or decoding matrices in the debugging
 *            part of the structure. This is because when we're in
 *            debugging mode asking about memory usage, we're usually
 *            interested in the estimated usage of the production
 *            code, because we're optimizing some parameter choices
 *            for example.
 */
size_t
p7_checkptmx_Sizeof(const P7_CHECKPTMX *ox)
{
  size_t n = sizeof(P7_CHECKPTMX);
  n += ox->allocN  * sizeof(float);     // neglects alignment overhead of up to p7_VALIGN bytes
  n += ox->allocR  * sizeof(float *);	  
  return n;
}

/* Function:  p7_checkptmx_MinSizeof()
 * Synopsis:  Returns minimum required size of a <P7_CHECKPTMX>, in bytes.
 *
 * Purpose:   Calculate and return the minimal required size, in bytes,
 *            of a checkpointed f/b matrix, for a comparison of a profile
 *            of length <M> to a sequence of length <L>.
 *            
 *            Does not require having an actual DP matrix allocated.
 *            We use this function when planning/profiling memory
 *            allocation strategies.
 */
size_t
p7_checkptmx_MinSizeof(int M, int L)
{
  size_t n    = sizeof(P7_CHECKPTMX);
  int    Q    = P7_Q(M,p7_VMAX_FB);               // number of vectors needed
  int    minR = 3 + (int) ceil(minimum_rows(L));  // 3 = Ra, 2 rows for backwards, 1 for fwd[0]

  n += minR * (Q * p7C_NSCELLS * p7_VMAX_FB);         // each row: Q supervectors. Each supervector: p7C_NSCELLS=3 vectors MDI. Each vector: p7_VMAX_FB cells
  n += minR * (ESL_UPROUND(p7C_NXCELLS, p7_VMAX_FB)); // dp_mem, specials: maintaining vector memory alignment 
  n += minR * sizeof(float *);                        // dpf[] row ptrs
  return n;
}




/* Function:  p7_checkptmx_Destroy()
 * Synopsis:  Frees a <P7_CHECKPTMX>.
 *
 * Purpose:   Free the <P7_CHECKPTMX> <ox>. <ox> may be <NULL>,
 *            or incompletely allocated.
 */
void
p7_checkptmx_Destroy(P7_CHECKPTMX *ox)
{
  if (ox)
    {
      if (ox->dp_mem) esl_alloc_free(ox->dp_mem);
      if (ox->dpf)    free(ox->dpf);

#if eslDEBUGLEVEL > 0
      if (ox->fwd)    p7_refmx_Destroy(ox->fwd);
      if (ox->bck)    p7_refmx_Destroy(ox->bck);
      if (ox->pp)     p7_refmx_Destroy(ox->pp);
#endif
    }
}
/*--------------- end, P7_CHECKPTMX object -----------------------*/



/*****************************************************************
 * 2. Debugging, development routines
 *****************************************************************/
#if eslDEBUGLEVEL > 0

/* checkptmx_DecodeX()
 * Converts a special X statecode to a string.
 * Used below in _DumpFBHeader()
 */
static char *
checkptmx_DecodeX(enum p7c_xcells_e xcode)
{
  switch (xcode) {
  case p7C_E:     return "E"; 
  case p7C_N:     return "N"; 
  case p7C_JJ:    return "JJ"; 
  case p7C_J:     return "J"; 
  case p7C_B:     return "B"; 
  case p7C_CC:    return "CC"; 
  case p7C_C:     return "C"; 
  case p7C_SCALE: return "SCALE"; 
  default:        return "?";
  }
}


/* checkptmx_get_val()
 * Access a score in a striped vector row.
 * Used below in _DumpFBRow() 
 */
static float
checkptmx_get_val(P7_CHECKPTMX *ox, float *dpc, int k, enum p7c_scells_e s, int logify)
{
  if (k == 0) return 0.0;
  int q = P7_Q_FROM_K(k,ox->Q);
  int z = P7_Z_FROM_K(k,ox->Q);
  if (logify) return esl_logf(dpc[(q*p7C_NSCELLS + s)*ox->Vf + z]);
  else        return          dpc[(q*p7C_NSCELLS + s)*ox->Vf + z];
}



/* Function:  p7_checkptmx_SetDumpMode()
 * Synopsis:  Toggle dump mode flag in a P7_CHECKPTMX.
 *
 * Purpose:   Toggles whether DP matrix rows will be dumped for examination
 *            during <p7_ForwardFilter()>, <p7_BackwardFilter()>.
 *            Dumping has to be done row by row, on the fly during the
 *            DP calculations, not afterwards, because the DP routines
 *            are memory efficient (checkpointed) and they do not save
 *            all their rows.
 * 
 *            To turn on dumping, <truefalse> is <TRUE>, and <dfp> is an
 *            open <FILE> pointer for dump output (such as <stderr>).
 *            
 *            To turn off dumping, <truefalse> is <FALSE>, and <dfp>
 *            should be <NULL>.
 *            
 *            For this to have effect, <eslDEBUGLEVEL> compile-time
 *            flag must be nonzero. If it is not, no dumping will occur,
 *            and this call is a no-op.
 *
 * Args:      ox        - DP matrix object to set
 *            dfp       - open FILE * for debugging output, or NULL
 *            truefalse - TRUE to set dump mode, or FALSE to unset
 *
 * Returns:   <eslOK> on success.
 */
int
p7_checkptmx_SetDumpMode(P7_CHECKPTMX *ox, FILE *dfp, int true_or_false)
{
  ox->do_dumping    = true_or_false;
  ox->dfp           = dfp;
  return eslOK;
}


/* Function:  p7_checkptmx_DumpFBHeader()
 * Synopsis:  Prints the header of the fwd/bck dumps.
 *
 * Purpose:   Print the header that accompanies <p7_checkptmx_DumpFBRow()>.
 */
int
p7_checkptmx_DumpFBHeader(P7_CHECKPTMX *ox)
{
  int maxpfx = ox->dump_maxpfx;
  int width  = ox->dump_width;
  int M      = ox->M;
  int k;

  fprintf(ox->dfp, "%*s", maxpfx, "");
  fprintf(ox->dfp, "      ");
  for (k = 0; k <= M;          k++) fprintf(ox->dfp, " %*d", width, k);
  for (k = 0; k < p7C_NXCELLS; k++) fprintf(ox->dfp, " %*s", width, checkptmx_DecodeX(k));
  fputc('\n', ox->dfp);

  fprintf(ox->dfp, "%*s", maxpfx, "");
  fprintf(ox->dfp, "      ");
  for (k = 0; k <= M+p7C_NXCELLS;  k++) fprintf(ox->dfp, " %*s", width, "--------");
  fputc('\n', ox->dfp);

  return eslOK;
}


/* Function:  p7_checkptmx_DumpFBRow()
 * Synopsis:  Dump one row from fwd or bck version of the matrix.
 *
 * Purpose:   Dump current row <dpc> of forward or backward calculations from
 *            DP matrix <ox> for diagnostics. The index <rowi> is used
 *            as a row label, along with an additional free-text label
 *            <pfx>.  (The checkpointed backward implementation
 *            interleaves backward row calculations with recalculated
 *            fwd rows, both of which it is dumping; they need to be
 *            labeled something like "fwd" and "bck" to distinguish
 *            them in the debugging dump.)
 *
 *            Independent of vector ISA. A vectorized caller passes
 *            the striped row <dpc> cast to <float *>, and we access
 *            it here by translating to normal coords k=1..M.
 */
int
p7_checkptmx_DumpFBRow(P7_CHECKPTMX *ox, int rowi, float *dpc, char *pfx)
{
  float *xc        = dpc + ox->Vf * ox->Q * p7C_NSCELLS;
  int    logify    = (ox->dump_flags & p7_SHOW_LOG) ? TRUE : FALSE;
  int    maxpfx    = ox->dump_maxpfx;
  int    width     = ox->dump_width;
  int    precision = ox->dump_precision;
  int    k,z;

  /* Line 1. M cells, unstriped */
  fprintf(ox->dfp, "%*s %3d M", maxpfx, pfx, rowi);
  for (k = 0; k <= ox->M; k++) 
    fprintf(ox->dfp, " %*.*f", width, precision, checkptmx_get_val(ox, dpc, k, p7C_M, logify));

  /* Line 1 end: specials */
  for (z = 0; z < p7C_NXCELLS; z++)
    fprintf(ox->dfp, " %*.*f", width, precision, (logify ? esl_logf(xc[z]) : xc[z]));
  fputc('\n', ox->dfp);

  /* Line 2: I cells, unstriped */
  fprintf(ox->dfp, "%*s %3d I", maxpfx, pfx, rowi);
  for (k = 0; k <= ox->M; k++) 
    fprintf(ox->dfp, " %*.*f", width, precision, checkptmx_get_val(ox, dpc, k, p7C_I, logify));
  fputc('\n', ox->dfp);

  /* Line 3. D cells, unstriped */
  fprintf(ox->dfp, "%*s %3d D", maxpfx, pfx, rowi);
  for (k = 0; k <= ox->M; k++)
    fprintf(ox->dfp, " %*.*f", width, precision, checkptmx_get_val(ox, dpc, k, p7C_D, logify));
  fputc('\n', ox->dfp);
  fputc('\n', ox->dfp);

  return eslOK;
}
#endif // eslDEBUGLEVEL
/*---------------- end, debugging -------------------------------*/


