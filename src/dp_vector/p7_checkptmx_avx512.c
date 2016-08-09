
/* Implementation of P7_CHECKPTMX: checkpointed, striped vector DP matrix.
 * 
 * Contents:
 *    1. API for the P7_CHECKPTMX object
 *    2. Debugging, development routines.
 *    3. Internal routines.
 *    4. Copyright and license information.
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>

#include "easel.h"

#include "dp_vector/simdvec.h"
#include "dp_vector/p7_checkptmx.h"
#ifdef HAVE_AVX512
 #include "immintrin.h"
 #endif
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
p7_checkptmx_Create_avx512(int M, int L, int64_t ramlimit)
{
#ifdef HAVE_AVX512  
  P7_CHECKPTMX *ox = NULL;
  int          maxR;
  int          r;
  int          status;
  
  /* Validity of integer variable ranges may depend on design spec:                  */
  ESL_DASSERT1( (M <= 100000) );       /* design spec says, model length M <= 100000 */
  ESL_DASSERT1( (L <= 100000) );       /*           ... and,  seq length L <= 100000 */
  ESL_DASSERT1( (L >  0) );
  ESL_DASSERT1( (M >  0) );

  /* Level 1 allocation: the structure itself */
  ESL_ALLOC(ox, sizeof(P7_CHECKPTMX));
  ox->simd = AVX512;

  ox->dp_mem_AVX_512  = NULL;
  ox->dpf_AVX_512     = NULL;

  /* Set checkpointed row layout: allocR, R{abc}, L{abc} fields */
  ox->R0          = 3;	                                                   /* fwd[0]; bck[prv,cur] */

  // determine number of rows based on SSE vector sizes to make checking easier
  int64_t base_row_length = sizeof(float) * P7_NVF(M) * p7C_NSCELLS * p7_VNF;    /* accounts for main vector part of the row       */ 
  base_row_length     += ESL_UPROUND(sizeof(float) * p7C_NXCELLS, p7_VALIGN);  /* plus specials (must maintain memory alignment) */
  ox->ramlimit    = ramlimit;
  maxR            = (int) (ox->ramlimit / base_row_length); 
  set_row_layout(ox, L, maxR);

  ox->allocW_AVX_512      = sizeof(float) * P7_NVF_AVX_512(M) * p7C_NSCELLS * p7_VNF_AVX_512;    /* accounts for main vector part of the row       */ 
  ox->allocW_AVX_512     += ESL_UPROUND(sizeof(float) * p7C_NXCELLS, p7_VALIGN_AVX_512);  /* plus specials (must maintain memory alignment) */
  ox->allocR_AVX_512      = ox->R0 + ox->Ra + ox->Rb + ox->Rc;
  ox->validR_AVX_512      = ox->allocR_AVX_512;
  ESL_DASSERT1( (ox->allocW_AVX_512 % p7_VALIGN_AVX_512 == 0) ); /* verify alignment */

  /* Level 2 allocations: row pointers and dp cell memory */
  ox->nalloc_AVX_512 = ox->allocR_AVX_512 * ox->allocW_AVX_512;
  ESL_ALLOC( ox->dp_mem_AVX_512, ox->nalloc_AVX_512 + (p7_VALIGN_AVX_512-1));    /* (p7_VALIGN-1) because we'll hand-align memory */
  ESL_ALLOC( ox->dpf_AVX_512,    sizeof(float *) * ox->allocR_AVX_512);  
  // Static analyzers may complain about the above.
  // sizeof(float *) is correct, even though ox->dpf is char **.
  // ox->dpf will be cast to __m128 SIMD vector in DP code.

  ox->dpf_AVX_512[0] = (char *) ( ((uintptr_t) ox->dp_mem_AVX_512 + p7_VALIGN_AVX_512 - 1) & p7_VALIMASK_AVX_512); /* hand memory alignment */
  for (r = 1; r < ox->validR_AVX_512; r++)
    ox->dpf_AVX_512[r] = ox->dpf_AVX_512[0] + r * ox->allocW_AVX_512;

#ifdef p7_DEBUGGING
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
 
  ox->R_AVX_512  = 0;
  ox->Qf_AVX_512 = 0;
  
  return ox;

 ERROR:
  p7_checkptmx_Destroy(ox);
  return NULL;
#endif //HAVE_AVX512
#ifndef HAVE_AVX512
  return NULL;
#endif      
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
p7_checkptmx_GrowTo_avx512(P7_CHECKPTMX *ox, int M, int L)
{
#ifdef HAVE_AVX512  
  int     minR_chk      = (int) ceil(minimum_rows(L)) + ox->R0; /* minimum number of DP rows needed  */
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

  /* If we're debugging and we have stored copies of any matrices,
   * grow them too.  Must do this first, because we have an early exit
   * condition coming below.
   */
#ifdef p7_DEBUGGING
  if (ox->fwd && (status = p7_refmx_GrowTo(ox->fwd, M, L)) != eslOK) goto ERROR;
  if (ox->bck && (status = p7_refmx_GrowTo(ox->bck, M, L)) != eslOK) goto ERROR;
  if (ox->pp  && (status = p7_refmx_GrowTo(ox->pp,  M, L)) != eslOK) goto ERROR;
#endif

 

  /* Are current allocations satisfactory ? */
  /* When more than one p7_build variable is defined, this counts on the fact that we compute the number of 
  rows to allocate based on SSE data widths and force AVX, AVX2 to use the same number.  With multiple
  p7_build variables defined, this introduces redundant computation, but that's ok because that should only 
  happen when we're checking ISAs against each other */

  if (W <= ox->allocW_AVX_512 && ox->nalloc_AVX_512 <= ox->ramlimit)
    {
      if      (L + ox->R0 <= ox->validR_AVX_512) { set_full        (ox, L);             return eslOK; }
      else if (minR_chk   <= ox->validR_AVX_512) { set_checkpointed(ox, L, ox->validR_AVX_512); return eslOK; }
    }

  reset_dp_ptrs = FALSE; // Reset this just to prevent any wierdness with multiple ISA builds
  W  = sizeof(float) * P7_NVF_AVX_512(M) * p7C_NSCELLS * p7_VNF_AVX_512;     /* vector part of row (MDI)     */
  W += ESL_UPROUND(sizeof(float) * p7C_NXCELLS, p7_VALIGN_AVX_512);  /* float part of row (specials); must maintain p7_VALIGN-byte alignment */
  /* Do individual matrix rows need to expand? */
  if ( W > ox->allocW_AVX_512) 
    {
      ox->allocW_AVX_512    = W;
      ox->validR_AVX_512    = (int) (ox->nalloc_AVX_512 / ox->allocW_AVX_512); /* validR must be <= allocR */
      reset_dp_ptrs = TRUE;
    }

  /* Does matrix dp_mem need reallocation, either up or down? */
  maxR  = (int) (ox->nalloc_AVX_512 / ox->allocW_AVX_512);                      /* max rows if we use up to the recommended allocation size.      */
  if ( (ox->nalloc_AVX_512 > ox->ramlimit && minR_chk <= maxR) ||       /* we were redlined, and recommended alloc will work: so downsize */
       minR_chk > ox->validR_AVX_512)               /* not enough memory for needed rows: so upsize                   */
    {
      set_row_layout(ox, L, maxR);  // Skip this if we did a layout for SSE to force SSE and AVX_512 to use the same number of rows 
      ox->validR_AVX_512 = ox->R0 + ox->Ra + ox->Rb + ox->Rc;   /* this may be > allocR now; we'll reallocate dp[] next, if so     */
      ox->nalloc_AVX_512 = ox->validR_AVX_512 * ox->allocW_AVX_512;
      ESL_REALLOC(ox->dp_mem_AVX_512, ox->nalloc_AVX_512 + (p7_VALIGN_AVX_512-1)); /* (p7_VALIGN-1) because we will manually align dpf ptrs into dp_mem */
      reset_dp_ptrs = TRUE;
    }
  else  /* current validR will suffice, either full or checkpointed; we still need to calculate a layout */
    {
      if   (L+ox->R0 <= ox->validR_AVX_512) set_full(ox, L); 
      else                          set_checkpointed(ox, L, ox->validR_AVX_512);
    }
  
  /* Does the array of row ptrs need reallocation? */
  if (ox->validR_AVX_512 > ox->allocR_AVX_512)
    {
      ESL_REALLOC(ox->dpf_AVX_512, sizeof(float *) * ox->validR_AVX_512);
      ox->allocR_AVX_512    = ox->validR_AVX_512;
      reset_dp_ptrs = TRUE;
    }

  /* Do the row ptrs need to be reset? */
  if (reset_dp_ptrs)
    {
      ox->dpf_AVX_512[0] = (char *) ( ( (uintptr_t) ox->dp_mem_AVX_512 + p7_VALIGN_AVX_512 - 1) & p7_VALIMASK_AVX_512); /* vectors must be aligned on p7_VALIGN-byte boundary */
      for (r = 1; r < ox->validR_AVX_512; r++)
  ox->dpf_AVX_512[r] = ox->dpf_AVX_512[0] + (r * ox->allocW_AVX_512);
    }

  return eslOK;

 ERROR:
  return status;

#endif //HAVE_AVX512
#ifndef HAVE_AVX512
  return eslENORESULT;
#endif      
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
p7_checkptmx_Sizeof_avx512(const P7_CHECKPTMX *ox)
{
  size_t n = sizeof(P7_CHECKPTMX);
 #ifdef HAVE_AVX512
  n += ox->nalloc_AVX_512 + (p7_VALIGN_AVX_512-1);            /* +63 because of manual alignment */
  n += ox->allocR_AVX_512  * sizeof(float *);   
#endif  
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
p7_checkptmx_MinSizeof_avx512(int M, int L)
{
  size_t n    = sizeof(P7_CHECKPTMX);
#ifdef HAVE_AVX512  
  int    minR = 3 + (int) ceil(minimum_rows(L));  // 3 = Ra, 2 rows for backwards, 1 for fwd[0]
  n += p7_VALIGN_AVX_512-1;                                                  // dp_mem_AVX_512 has to be hand-aligned for vectors
  n += minR * (sizeof(float) * p7_VNF_AVX_512 * P7_NVF_AVX_512(M) * p7C_NSCELLS);            
  // dp_mem_AVX, main: QR supercells; each has p7C_NSCELLS=3 cells, MID; each cell is __m512 vector 16 of floats (p7_VNF_AVX_512 = 16 * float)
  n += minR * (ESL_UPROUND(sizeof(float) * p7C_NXCELLS, p7_VALIGN_AVX_512)); // dp_mem_AVX_512, specials: maintaining vector memory alignment 
  n += minR * sizeof(float *);                                       // dpf[] row ptrs
#endif 
  return n;
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
p7_checkptmx_Reuse_avx512(P7_CHECKPTMX *ox)
{
#ifdef HAVE_AVX512  
#ifdef p7_DEBUGGING
  int status;
#endif

  ox->M  = 0;
  ox->L  = 0;
  ox->R_AVX_512  = 0;
  ox->Qf_AVX_512 = 0;


#ifdef p7_DEBUGGING
  if (ox->fwd && (status = p7_refmx_Reuse(ox->fwd)) != eslOK) return status;
  if (ox->bck && (status = p7_refmx_Reuse(ox->bck)) != eslOK) return status;
  if (ox->pp  && (status = p7_refmx_Reuse(ox->pp))  != eslOK) return status;
  ox->bcksc = 0.0f;
#endif

  return eslOK;
#endif //HAVE_AVX512
#ifndef HAVE_AVX512
  return eslENORESULT;
#endif      
}


/* Function:  p7_checkptmx_Destroy()
 * Synopsis:  Frees a <P7_CHECKPTMX>.
 *
 * Purpose:   Free the <P7_CHECKPTMX> <ox>. <ox> may be <NULL>,
 *            or incompletely allocated.
 */
void
p7_checkptmx_Destroy_avx512(P7_CHECKPTMX *ox)
{
 if (ox) {

#ifdef HAVE_AVX512  
   if (ox->dp_mem_AVX_512) free(ox->dp_mem_AVX_512);
   if (ox->dpf_AVX_512)    free(ox->dpf_AVX_512);
#endif

#ifdef p7_DEBUGGING
   if (ox->fwd)    p7_refmx_Destroy(ox->fwd);
   if (ox->bck)    p7_refmx_Destroy(ox->bck);
   if (ox->pp)     p7_refmx_Destroy(ox->pp);
#endif
   free(ox);
 }
}
/*--------------- end, P7_CHECKPTMX object -----------------------*/



/*****************************************************************
 * 2. Debugging, development routines
 *****************************************************************/
#ifdef p7_DEBUGGING
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
p7_checkptmx_DumpFBRow_avx512(P7_CHECKPTMX *ox, int rowi, __m512 *dpc, char *pfx)
{
#ifdef HAVE_AVX512
  union { __m512 v; float x[p7_VNF_AVX_512]; } u;
  float *v         = NULL;		/*  */
  int    Q         = ox->Qf;
  int    M         = ox->M;
  float *xc        = (float *) (dpc + Q*p7C_NSCELLS);
  int    logify    = (ox->dump_flags & p7_SHOW_LOG) ? TRUE : FALSE;
  int    maxpfx    = ox->dump_maxpfx;
  int    width     = ox->dump_width;
  int    precision = ox->dump_precision;
  int    k,q,z;
  int    status;

  ESL_ALLOC(v, sizeof(float) * ( (Q*p7_VNF_AVX_512) + 1));
  v[0] = 0.;

  /* Line 1. M cells: unpack, unstripe, print */
  for (q = 0; q < Q; q++) {
    u.v = P7C_MQ(dpc, q);
    for (z = 0; z < p7_VNF_AVX_512; z++) v[q+Q*z+1] = u.x[z];
  }
  fprintf(ox->dfp, "%*s %3d M", maxpfx, pfx, rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, " %*.*f", width, precision, (logify ? esl_logf(v[k]) : v[k]));
  /* a static analyzer may complain about v[k] being uninitialized
   * if it isn't smart enough to see that M,Q are linked.
   */

  /* Line 1 end: Specials */
  for (z = 0; z < p7C_NXCELLS; z++)
    fprintf(ox->dfp, " %*.*f", width, precision, (logify ? esl_logf(xc[z]) : xc[z]));
  fputc('\n', ox->dfp);

  /* Line 2: I cells: unpack, unstripe, print */
  for (q = 0; q < Q; q++) {
    u.v = P7C_IQ(dpc, q);
    for (z = 0; z < p7_VNF_AVX_512; z++) v[q+Q*z+1] = u.x[z];
  }
  fprintf(ox->dfp, "%*s %3d I", maxpfx, pfx, rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, " %*.*f", width, precision, (logify ? esl_logf(v[k]) : v[k]));
  fputc('\n', ox->dfp);

  /* Line 3. D cells: unpack, unstripe, print */
  for (q = 0; q < Q; q++) {
    u.v = P7C_DQ(dpc, q);
    for (z = 0; z < p7_VNF_AVX_512; z++) v[q+Q*z+1] = u.x[z];
  }
  fprintf(ox->dfp, "%*s %3d D", maxpfx, pfx, rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, " %*.*f", width, precision, (logify ? esl_logf(v[k]) : v[k]));
  fputc('\n', ox->dfp);
  fputc('\n', ox->dfp);

  free(v);
  return eslOK;

 ERROR:
  if (v) free(v);
  return status;
#endif //HAVE_AVX512
#ifndef HAVE_AVX512
  return eslENORESULT
#endif    
}

#endif /*p7_DEBUGGING*/
/*---------------- end, debugging -------------------------------*/

/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
