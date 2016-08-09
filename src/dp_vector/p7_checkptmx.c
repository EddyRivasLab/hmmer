/* Implementation of P7_CHECKPTMX: checkpointed, striped vector DP matrix.
 * 
 * Contents:
 *    1. API for the P7_CHECKPTMX object
 *    2. Debugging, development routines.
 *    3. Internal routines.
 *    4. Copyright and license information.
*    Many of the functions in this file are now dispatch wrappers that just call the appropriate SIMD version of the function
*   in cases where the correcct version is known, you can achieve (slightly) better performance by calling it directly
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"

#include "dp_vector/simdvec.h"
#include "dp_vector/p7_checkptmx.h"
#include "hardware/hardware.h"

/*****************************************************************
 * 1. API for the <P7_CHECKPTMX> object
 *****************************************************************/

/* Function:  p7_checkptmx_Create()
 * Synopsis:  Allocate a new <P7_CHECKPTMX> object.
 *
 * Purpose:   Allocate a new <P7_CHECKPTMX> checkpointed, striped vector
 *            DP matrix sufficient for the Forward/Backward local
 *            decoding calculation for a query model
 *            of up to length <M> and a target sequence of up to
 *            length <L>.
 *            
 *            Try to keep the allocation within <ramlimit> bytes in
 *            memory.  For example, <ramlimit=ESL_MBYTES(128)>, sets a
 *            recommended memory limit of 128 MiB. Allocation can
 *            exceed this, if even a fully checkpointed <MxL>
 *            comparison requires it -- but in this case, any
 *            subsequent <p7_checkptmx_GrowTo()> call that attempts to
 *            reuse the matrix will try to reallocated it back
 *            downwards to the <ramlimit>.
 *            
 *            Choice of <ramlimit> should take into account how many
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
p7_checkptmx_Create(int M, int L, int64_t ramlimit, SIMD_TYPE simd)
{
switch(simd){
    case SSE:
      return p7_checkptmx_Create_sse(M, L, ramlimit);
      break;
    case AVX:
      return p7_checkptmx_Create_avx(M, L, ramlimit);
      break;
    case AVX512:
      return p7_checkptmx_Create_avx512(M, L, ramlimit);
      break;
    case NEON: case NEON64:
      return p7_checkptmx_Create_neon(M, L, ramlimit);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_checkptmx_Create");  
  }
}

/* Function:  p7_checkptmx_GrowTo()
 * Synopsis:  Resize checkpointed DP matrix for new seq/model comparison.
 *
 * Purpose:   Given an existing checkpointed matrix structure <ox>,
 *            and the dimensions <M> and <L> of a new comparison,
 *            reallocate and reinitialize <ox>.
 *
 *            Essentially the same as free'ing the previous matrix and
 *            creating a new one -- but minimizes expensive memory
 *            allocation/reallocation calls.
 *            
 *            Usually <ox> only grows. The exception is if <ox> is
 *            redlined (over its recommended allocation) and the new
 *            problem size <M,L> can fit in the preset recommended
 *            allocation, then <ox> is reallocated down to the smaller
 *            recommended size.
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
p7_checkptmx_GrowTo(P7_CHECKPTMX *ox, int M, int L)
{
switch(ox->simd){
    case SSE:
      return p7_checkptmx_GrowTo_sse(ox, M, L);
      break;
    case AVX:
      return p7_checkptmx_GrowTo_avx(ox, M, L);
      break;
    case AVX512:
      return p7_checkptmx_GrowTo_avx512(ox, M, L);
      break;
    case NEON: case NEON64:
      return p7_checkptmx_GrowTo_neon(ox, M, L);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_checkptmx_Growto");  
  }
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
switch(ox->simd){
    case SSE:
      return p7_checkptmx_Sizeof_sse(ox);
      break;
    case AVX:
      return p7_checkptmx_Sizeof_avx(ox);
      break;
    case AVX512:
      return p7_checkptmx_Sizeof_avx512(ox);
      break;
    case NEON: case NEON64:
      return p7_checkptmx_Sizeof_neon(ox);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_checkptmx_Sxzeof");  
  }
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
p7_checkptmx_MinSizeof(int M, int L, SIMD_TYPE simd)
{
  switch(simd){
    case SSE:
      return p7_checkptmx_MinSizeof_sse(M, L);
      break;
    case AVX:
      return p7_checkptmx_MinSizeof_avx(M, L);
      break;
    case AVX512:
      return p7_checkptmx_MinSizeof_avx512(M, L);
      break;
    case NEON: case NEON64:
      return p7_checkptmx_MinSizeof_neon(M, L);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_checkptmx_MinSizeof");  
  }
}


/* Function:  p7_checkptmx_Reuse()
 * Synopsis:  Recycle a checkpointed vector DP matrix.
 *
 * Purpose:   Resets the checkpointed vector DP matrix <ox> for reuse,
 *            minimizing free/malloc wastefulness. All information
 *            specific to the DP problem we just computed is
 *            reinitialized. All allocations (and information about
 *            those allocations) are preserved.
 *            
 *            Caller will still need to call <p7_checkptmx_GrowTo()>
 *            before each new DP, to be sure that the allocations are
 *            sufficient, and checkpointed rows are laid out.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_checkptmx_Reuse(P7_CHECKPTMX *ox)
{

 switch(ox->simd){
    case SSE:
      return p7_checkptmx_Reuse_sse(ox);
      break;
    case AVX:
      return p7_checkptmx_Reuse_avx(ox);
      break;
    case AVX512:
      return p7_checkptmx_Reuse_avx512(ox);
      break;
    case NEON: case NEON64:
      return p7_checkptmx_Reuse_neon(ox);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_checkptmx_Reuse");  
  } 
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
 switch(ox->simd){
    case SSE:
      p7_checkptmx_Destroy_sse(ox);
      break;
    case AVX:
      p7_checkptmx_Destroy_avx(ox);
      break;
    case AVX512:
      p7_checkptmx_Destroy_avx512(ox);
      break;
    case NEON: case NEON64:
      p7_checkptmx_Destroy_neon(ox);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_checkptmx_Destroy");  
  }
}
/*--------------- end, P7_CHECKPTMX object -----------------------*/



/*****************************************************************
 * 2. Debugging, development routines
 *****************************************************************/

/* Function:  p7_checkptmx_SetDumpMode()
 * Synopsis:  Toggle dump mode flag in a P7_CHECKPTMX.
 *
 * Purpose:   Toggles whether DP matrix rows will be dumped for examination
 *            during <p7_ForwardFilter()>, <p7_BackwardFilter()>.
 *            Dumping has to be done row by row, on the fly during the
 *            DP calculations, not afterwards, because these routines
 *            are memory efficient (checkpointed) and they do not save
 *            all their rows.
 * 
 *            To turn on dumping, <truefalse> is <TRUE>, and <dfp> is an
 *            open <FILE> pointer for dump output (such as <stderr>).
 *            
 *            To turn off dumping, <truefalse> is <FALSE>, and <dfp>
 *            should be <NULL>.
 *            
 *            For this to have effect, <p7_DEBUGGING> compile-time
 *            flag must be on. If it is not, no dumping will occur,
 *            and this call is a no-op.
 *
 * Args:      ox        - DP matrix object to set
 *            dfp       - open FILE * for debugging output, or NULL
 *            truefalse - TRUE to set dump mode, or FALSE to unset
 *
 * Returns:   <eslOK> on success.
 */
int
p7_checkptmx_SetDumpMode(P7_CHECKPTMX *ox, FILE *dfp, int truefalse)
{
#ifdef p7_DEBUGGING
  ox->do_dumping    = truefalse;
  ox->dfp           = dfp;
#endif
  return eslOK;
}


#ifdef p7_DEBUGGING

/* Function:  p7_checkptmx_DecodeX()
 * Synopsis:  Convert a special X statecode to a string.
 */
char *
p7_checkptmx_DecodeX(enum p7c_xcells_e xcode)
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

/* Function:  p7_checkptmx_DumpFBHeader()
 * Synopsis:  Prints the header of the fwd/bck dumps.
 *
 * Purpose:   Print the header that accompanies 
 *            <p7_checkptmx_DumpFBRow()>.
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
  for (k = 0; k < p7C_NXCELLS; k++) fprintf(ox->dfp, " %*s", width, p7_checkptmx_DecodeX(k));
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
 */
int
p7_checkptmx_DumpFBRow(P7_CHECKPTMX *ox, int rowi, __m128 *dpc, char *pfx)
{
  switch(ox->simd){
    case SSE:
      return p7_checkptmx_DumpFBRow_sse(ox, rowi, dpc, pfx);
      break;
    case AVX:
      return p7_checkptmx_DumpFBRow_avx(ox, rowi, dpc, pfx);
      break;
    case AVX512:
      return p7_checkptmx_DumpFBRow_avx512(ox, rowi, dpc, pfx);
      break;
    case NEON: case NEON64:
      return p7_checkptmx_DumpFBRow_neon(ox, rowi, dpc, pfx);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_checkptmx_DumpFBRow");  
  }
}

#endif /*p7_DEBUGGING*/
/*---------------- end, debugging -------------------------------*/




/*****************************************************************
 * 3. Internal routines
 *****************************************************************/

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
void
set_row_layout(P7_CHECKPTMX *ox, int allocL, int maxR)
{
  double Rbc      = minimum_rows(allocL);               
  int    minR_chk = ox->R0 + (int) ceil(Rbc);	          /* min # rows we need for checkpointing          */
  int    minR_all = ox->R0 + allocL;		          /* min # rows we need for full matrix            */

  if      (minR_all <= maxR) set_full        (ox, allocL);
  else if (minR_chk <= maxR) set_checkpointed(ox, allocL, maxR);
  else                       set_redlined    (ox, allocL, Rbc);
}

/* A "full" matrix is easy: Ra = La = L, using Ra+R0 <= maxR rows total. */
void
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

/* minimum_rows()
 * 
 * Calculate the minimum number of rows needed to checkpoint a
 * forward matrix for a sequence of length L, exclusive of 
 * other constant row overhead (R0: two backwards rows, fwd[0]
 * initial row).
 * 
 * This value is a double; if it has a fraction, a partial checkpoint
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
 * Given L and maxR, solve for the number of checkpointed
 * rows (Rb+Rc) we need. The value is a double; if it has
 * a fractional part, then a partial checkpoint block is
 * needed, Rb==1.
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
/*----------------- end, internals ------------------------------*/


/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
