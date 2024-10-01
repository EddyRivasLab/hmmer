/* SSE implementation of an optimized profile structure.
 * 
 * Contents:
 *   1. The P7_OMX structure: a dynamic programming matrix
 *   2. Debugging dumps of P7_OMX structures
 * 
 * See also:
 *   p7_omx.ai - figure illustrating the layout of a P7_OMX.
 *
 * SRE, Sun Nov 25 11:26:48 2007 [Casa de Gatos]
 */
#include <p7_config.h>

#include <stdio.h>
#include <math.h>
#include <float.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"
#include "esl_sse.h"
#include "esl_cpu.h"
#include "hmmer.h"
#include "impl_avx.h"

/*****************************************************************
 * 1. The P7_OMX structure: a dynamic programming matrix
 *****************************************************************/

/* Function:  p7_omx_Create()
 * Synopsis:  Create an optimized dynamic programming matrix.
 * Incept:    SRE, Tue Nov 27 08:48:20 2007 [Janelia]
 *
 * Purpose:   Allocates a reusable, resizeable <P7_OMX> for models up to
 *            size <allocM> and target sequences up to length
 *            <allocL/allocXL>, for use by any of the various optimized
 *            DP routines.
 *            
 *            To allocate the very memory-efficient one-row matrix
 *            used by *Filter() and *Score() functions that only
 *            calculate scores, <allocM=M>, <allocL=0>, and
 *            <allocXL=0>.
 *            
 *            To allocate the reasonably memory-efficient linear
 *            arrays used by *Parser() functions that only keep
 *            special (X) state scores, <allocM=M>, <allocL=0>,
 *            and <allocXL=L>.
 *            
 *            To allocate a complete matrix suitable for functions
 *            that need the whole DP matrix for traceback, sampling,
 *            posterior decoding, or reestimation, <allocM=M> and
 *            <allocL=allocXL=L>.
 *
 * Returns:   a pointer to the new <P7_OMX>.
 *
 * Throws:    <NULL> on allocation failure.
 */

static P7_OMX *p7_omx_Create_Dispatcher(int allocM, int allocL, int allocXL){
#ifdef P7_TEST_ALL_SIMD
  p7_omx_Create = p7_omx_Create_test_all;
  return p7_omx_Create_test_all_simd(allocM, allocL, allocXL);
#endif

#ifdef P7_TEST_SSE_AVX
  p7_omx_Create = p7_omx_Create_test_sse_avx;
  return p7_omx_Create_test_sse_avx(allocM, allocL, allocXL); 
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
  if (esl_cpu_has_avx512())
    {
      p7_omx_Create = p7_omx_Create_avx512;
      return p7_omx_Create_avx512(allocM, allocL, allocXL);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      p7_omx_Create = p7_omx_Create_avx;
      return p7_omx_Create_avx(allocM, allocL, allocXL);
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse())
    {
      p7_omx_Create = p7_omx_Create_sse;
      return p7_omx_Create_sse(allocM, allocL, allocXL);
    }
#endif

  p7_Die("p7_omx_Create_dispatcher found no vector implementation - that shouldn't happen.");
  return NULL;
}

P7_OMX *
(* p7_omx_Create)(int allocM, int allocL, int allocXL)=p7_omx_Create_Dispatcher;

/* Function:  p7_omx_GrowTo()
 * Synopsis:  Assure that a DP matrix is big enough.
 * Incept:    SRE, Thu Dec 20 09:27:07 2007 [Janelia]
 *
 * Purpose:   Assures that an optimized DP matrix <ox> is allocated for
 *            a model up to <allocM> in length; if not, reallocate to
 *            make it so.
 *            
 *            Because the optimized matrix is one-row, only the model
 *            length matters; the target sequence length isn't
 *            relevant.
 *
 * Returns:   <eslOK> on success, and <gx> may be reallocated upon
 *            return; any data that may have been in <gx> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <gx> must be assumed to be invalidated.
 */

static int p7_omx_GrowTo_Dispatcher(P7_OMX *ox, int allocM, int allocL, int allocXL){
  #ifdef P7_TEST_ALL_SIMD
  p7_omx_GrowTo = p7_omx_GrowTo_test_all;
  return p7_omx_GrowTo_test_all_simd(ox, allocM, allocL, allocXL);
#endif

#ifdef P7_TEST_SSE_AVX
  p7_omx_GrowTo = p7_omx_GrowTo_test_sse_avx;
  return p7_omx_GrowTo_test_sse_avx(ox, allocM, allocL, allocXL); 
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
  if (esl_cpu_has_avx512())
    {
      p7_omx_GrowTo = p7_omx_GrowTo_avx512;
      return p7_omx_GrowTo_avx512(ox, allocM, allocL, allocXL);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      p7_omx_GrowTo = p7_omx_GrowTo_avx;
      return p7_omx_GrowTo_avx(ox, allocM, allocL, allocXL);
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse())
    {
      p7_omx_GrowTo = p7_omx_GrowTo_sse;
      return p7_omx_GrowTo_sse(ox, allocM, allocL, allocXL);
    }
#endif

  p7_Die("p7_omx_GrowTo_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;
}

int (* p7_omx_GrowTo)(P7_OMX *ox, int allocM, int allocL, int allocXL) = p7_omx_GrowTo_Dispatcher;

/* Function:  p7_omx_FDeconvert()
 * Synopsis:  Convert an optimized DP matrix to generic one.
 * Incept:    SRE, Tue Aug 19 17:58:13 2008 [Janelia]
 *
 * Purpose:   Convert the 32-bit float values in optimized DP matrix
 *            <ox> to a generic one <gx>. Caller provides <gx> with sufficient
 *            space to hold the <ox->M> by <ox->L> matrix.
 *            
 *            This function is used to gain access to the
 *            somewhat more powerful debugging and display
 *            tools available for generic DP matrices.
 */

static int p7_omx_FDeconvert_Dispatcher(P7_OMX *ox, P7_GMX *gx){
  #ifdef P7_TEST_ALL_SIMD
  p7_omx_FDeconvert = p7_omx_Fdeconvert_test_all;
  return p7_omx_FDeconvert_test_all_simd(ox, gx);
#endif

#ifdef P7_TEST_SSE_AVX
  p7_omx_FDeconvert = p7_omx_FDeconvert_test_sse_avx;
  return p7_omx_FDeconvert_test_sse_avx(ox, gx); 
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
  if (esl_cpu_has_avx512())
    {
      p7_omx_FDeconvert = p7_omx_FDeconvert_avx512;
      return p7_omx_FDeconvert_avx512(ox, gx);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      p7_omx_FDeconvert = p7_omx_FDeconvert_avx;
      return p7_omx_FDeconvert_avx(ox, gx);
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse())
    {
      p7_omx_FDeconvert = p7_omx_FDeconvert_sse;
      return p7_omx_FDeconvert_sse(ox, gx);
    }
#endif

  p7_Die("p7_omx_FDeconvert_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;
}

int (*p7_omx_FDeconvert)(P7_OMX *ox, P7_GMX *gx) = p7_omx_FDeconvert_Dispatcher;


/* Function:  p7_omx_Reuse()
 * Synopsis:  Recycle an optimized DP matrix.
 * Incept:    SRE, Wed Oct 22 11:31:00 2008 [Janelia]
 *
 * Purpose:   Recycles <ox> for re-use.
 *
 * Returns:   <eslOK> on success.
 */
// Doesn't need separate versions for different SIMD ISAs
int
p7_omx_Reuse(P7_OMX *ox)
{
  ox->M              = 0;
  ox->L              = 0;
  ox->totscale       = 0.0;
  ox->has_own_scales = TRUE;	/* default assumes a Forward matrix, with its own scale factors */
#if eslDEBUGLEVEL > 0
  ox->debugging      = FALSE;
  ox->dfp            = NULL;
#endif
  return eslOK;
}




/* Function:  p7_omx_Destroy()
 * Synopsis:  Frees an optimized DP matrix.
 * Incept:    SRE, Tue Nov 27 09:11:42 2007 [Janelia]
 *
 * Purpose:   Frees optimized DP matrix <ox>.
 *
 * Returns:   (void)
 */

static void p7_omx_Destroy_Dispatcher(P7_OMX *ox){
#ifdef P7_TEST_ALL_SIMD
  p7_omx_Destroy = p7_omx_Destroy_test_all;
  return p7_omx_Destroy_test_all_simd(ox);
#endif

#ifdef P7_TEST_SSE_AVX
  p7_omx_Destroy = p7_omx_Destroy_test_sse_avx;
  return p7_omx_Destroy_test_sse_avx(ox); 
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
  if (esl_cpu_has_avx512())
    {
      p7_omx_Destroy = p7_omx_Destroy_avx512;
      return p7_omx_Destroy_avx512(ox);
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      p7_omx_Destroy = p7_omx_Destroy_avx;
      return p7_omx_Destroy_avx(ox);
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse())
    {
      p7_omx_Destroy = p7_omx_Destroy_sse;
      return p7_omx_Destroy_sse(ox);
    }
#endif

  p7_Die("P7_omx_Destroy_dispatcher found no vector implementation - that shouldn't happen.");
}


void
(* p7_omx_Destroy)(P7_OMX *ox) = p7_omx_Destroy_Dispatcher;

/*------------------- end, P7_OMX structure ---------------------*/



/*****************************************************************
 * 2. Debugging dumps of P7_OMX structures
 *****************************************************************/
/* Because the P7_OMX may be a one-row DP matrix, we can't just run a
 * DP calculation and then dump a whole matrix; we have to dump each
 * row one at a time, as the DP calculation is progressing. Thus we
 * need to call the dump from *within* some DP routines. We'd rather not
 * have anything like this in production code - not even a flag check.
 * So, we use a compile-time debugging idiom, with conditionally
 * compiled debugging code that's added to the DP routines to check a
 * debugging flag in the P7_OMX structure; if it's up, we dump a row.
 *
 * Therefore, the externally exposed API call is p7_omx_SetDumpMode(),
 * rather than the dumping routine itself; and all p7_omx_SetDumpMode()
 * does is sets the debugging flag in <ox>.
 */

/* Function:  p7_omx_SetDumpMode()
 * Synopsis:  Set an optimized DP matrix to be dumped for debugging.
 * Incept:    SRE, Thu Dec 13 10:24:38 2007 [Janelia]
 *
 * Purpose:   Sets debugging mode for DP matrix <ox>.  If <truefalse>
 *            flag is <TRUE>, then whenever a dynamic programming
 *            calculation is run, dump DP matrix <ox> to stream <fp>
 *            for diagnostics.
 *            
 *            When the dump mode is on, the DP routine itself actually
 *            does the dumping, because it has to dump after every row
 *            is calculated. (We're doing an optimized one-row
 *            calculation.)
 *            
 *            If the code has not been compiled with the
 *            <eslDEBUGLEVEL> flag set nonzero, this function is a no-op.
 *
 * Args:      fp        - output stream for diagnostics (stdout, perhaps)
 *            ox        - DP matrix to set debugging mode
 *            truefalse - TRUE to set dumping, FALSE to unset
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      J2/62.
 */
// Doesn't need separate routines for different SIMD ISAs
int
p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse)
{
#if eslDEBUGLEVEL > 0
  ox->debugging = truefalse;
  ox->dfp       = fp;
#endif
  return eslOK;
}


/* Function:  p7_omx_DumpMFRow()
 * Synopsis:  Dump one row from MSV uchar version of a DP matrix.
 * Incept:    SRE, Wed Jul 30 16:47:26 2008 [Janelia]
 *
 * Purpose:   Dump current row of uchar part of DP matrix <ox> for diagnostics,
 *            and include the values of specials <xE>, etc. The index <rowi> for
 *            the current row is used as a row label. This routine has to be 
 *            specialized for the layout of the MSVFilter() row, because it's
 *            all match scores dp[0..q..Q-1], rather than triplets of M,D,I.
 * 
 *            If <rowi> is 0, print a header first too.
 * 
 *            The output format is coordinated with <p7_gmx_Dump()> to
 *            facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. 
 */
// SIMD implementations in p7_omx_sse.c, p7_omx_avx.c, and p7_omx_avx512.c


/* Function:  p7_omx_DumpVFRow()
 * Synopsis:  Dump current row of ViterbiFilter (int16) part of <ox> matrix.
 * Incept:    SRE, Wed Jul 30 16:43:21 2008 [Janelia]
 *
 * Purpose:   Dump current row of ViterbiFilter (int16) part of DP
 *            matrix <ox> for diagnostics, and include the values of
 *            specials <xE>, etc. The index <rowi> for the current row
 *            is used as a row label.
 *
 *            If <rowi> is 0, print a header first too.
 * 
 *            The output format is coordinated with <p7_gmx_Dump()> to
 *            facilitate comparison to a known answer.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */

// SIMD implementations in p7_omx_sse.c, p7_omx_avx.c, and p7_omx_avx512.c

/* Function:  p7_omx_DumpFBRow()
 * Synopsis:  Dump one row from float part of a DP matrix.
 * Incept:    SRE, Wed Jul 30 16:45:16 2008 [Janelia]
 *
 * Purpose:   Dump current row of Forward/Backward (float) part of DP
 *	      matrix <ox> for diagnostics, and include the values of
 *	      specials <xE>, etc. The index <rowi> for the current row
 *	      is used as a row label. 
 *
 *            The output format of the floats is controlled by
 *	      <width>, <precision>; 8,5 is good for pspace, 5,2 is
 *	      fine for lspace.
 * 	       								       
 * 	      If <rowi> is 0, print a header first too.			       
 * 	       								       
 * 	      If <logify> is TRUE, then scores are printed as log(score); this is 
 * 	      useful for comparing DP with pspace scores with other DP matrices   
 * 	      (like generic P7_GMX ones) that use log-odds scores.		       
 * 	       								       
 * 	      The output format is coordinated with <p7_gmx_Dump()> to	       
 * 	      facilitate comparison to a known answer.                            
 * 
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.  
 */

// SIMD implementations in p7_omx_sse.c, p7_omx_avx.c, and p7_omx_avx512.c


 
/*------------- end, debugging dumps of P7_OMX ------------------*/


