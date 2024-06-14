/* Routines for the P7_OPROFILE structure:  
 * a search profile in an optimized implementation.
 * 
 * Contents:
 *   1. The P7_OPROFILE object: allocation, initialization, destruction.
 *   2. Conversion from generic P7_PROFILE to optimized P7_OPROFILE
 *   3. Conversion from optimized P7_OPROFILE to compact score arrays
 *   4. Debugging and development utilities.
 *   5. Benchmark driver.
 *   6. Unit tests.
 *   7. Test driver.
 *   8. Example.
 */
#include <p7_config.h>

#include <stdio.h>
#include <string.h>
#include <math.h>		/* roundf() */


#include "easel.h"
#include "esl_random.h"
#include "esl_cpu.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:03:19 2007 [Casa de Gatos]
 *
 * Purpose:   Allocate for profiles of up to <allocM> nodes for digital alphabet <abc>.
 *
 * Throws:    <NULL> on allocation error.
 */

static P7_OPROFILE *p7_oprofile_Create_Dispatcher(int allocM, const ESL_ALPHABET *abc);

P7_OPROFILE * (*p7_oprofile_Create)(int allocM, const ESL_ALPHABET *abc) = p7_oprofile_Create_Dispatcher;

P7_OPROFILE * p7_oprofile_Create_Dispatcher(int allocM, const ESL_ALPHABET *abc){
#ifdef P7_TEST_ALL_SIMD
 p7_oprofile_Create = p7_oprofile_Create_test_all;
 return p7_oprofile_Create_test_all(allocM, abc);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_Create = p7_oprofile_Create_test_sse_avx;
 return p7_oprofile_Create_test_sse_avx(allocM, abc);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_Create = p7_oprofile_Create_avx512;
     return p7_oprofile_Create_avx512(allocM, abc);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_Create = p7_oprofile_Create_avx;
     return p7_oprofile_Create_avx(allocM, abc); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_Create = p7_oprofile_Create_sse;
      return p7_oprofile_Create_sse(allocM, abc);
    }
#endif

  p7_Die("p7_oprofile_Create_dispatcher found no vector implementation - that shouldn't happen.");
  return(NULL);
}


/* Function:  p7_oprofile_IsLocal()
 * Synopsis:  Returns TRUE if profile is in local alignment mode.
 * Incept:    SRE, Sat Aug 16 08:46:00 2008 [Janelia]
 */
//This function does not need separate versions for different x86 SIMD ISAs
int
p7_oprofile_IsLocal(const P7_OPROFILE *om)
{
  if (om->mode == p7_LOCAL || om->mode == p7_UNILOCAL) return TRUE;
  return FALSE;
}



/* Function:  p7_oprofile_Destroy()
 * Synopsis:  Frees an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:22:21 2007 [Casa de Gatos]
 */
static void p7_oprofile_Destroy_Dispatcher(P7_OPROFILE *om);

void (*p7_oprofile_Destroy)(P7_OPROFILE *om) = p7_oprofile_Destroy_Dispatcher;

void p7_oprofile_Destroy_Dispatcher(P7_OPROFILE *om)
{
#ifdef P7_TEST_ALL_SIMD
 p7_oprofile_Destroy = p7_oprofile_Destroy_test_all;
 p7_oprofile_Destroy_test_all(om);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_Destroy = p7_oprofile_Destroy_test_sse_avx;
 p7_oprofile_Destroy_test_sse_avx(om);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_Destroy = p7_oprofile_Destroy_avx512;
     p7_oprofile_Destroy_avx512(om);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_Destroy = p7_oprofile_Destroy_avx;
     p7_oprofile_Destroy_avx(om); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_Destroy = p7_oprofile_Destroy_sse;
      p7_oprofile_Destroy_sse(om);
    }
#endif

  p7_Die("p7_oprofile_Destroy_dispatcher found no vector implementation - that shouldn't happen.");
}


/* TODO: this is not following the _Copy interface guidelines; it's a _Clone */
/* Function:  p7_oprofile_Copy()
 * Synopsis:  Perform a full copy of a P7_Oprofile structure
 * Incept:    SRE, Sun Nov 25 12:03:19 2007 [Casa de Gatos]
 *
 * Purpose:   Creates a copy of the input data structure.
 *
 * Throws:    <NULL> on allocation error.
 */
static P7_OPROFILE * p7_oprofile_Copy_Dispatcher(P7_OPROFILE *om1);

P7_OPROFILE * (*p7_oprofile_Copy)(P7_OPROFILE *om1)= p7_oprofile_Copy_Dispatcher;

P7_OPROFILE * p7_oprofile_Copy_Dispatcher(P7_OPROFILE *om1){
#ifdef P7_TEST_ALL_SIMD
 p7_oprofile_Copy = p7_oprofile_Copy_test_all;
 return p7_oprofile_Copy_test_all(om1);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_Copy = p7_oprofile_Copy_test_sse_avx;
 return p7_oprofile_Copy_test_sse_avx(om1);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_Copy = p7_oprofile_Copy_avx512;
     return p7_oprofile_Copy_avx512(om1);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_Copy = p7_oprofile_Copy_avx;
     return p7_oprofile_Copy_avx(om1); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_Copy = p7_oprofile_Copy_sse;
      return p7_oprofile_Copy_sse(om1);
    }
#endif

  p7_Die("p7_oprofile_Copy_dispatcher found no vector implementation - that shouldn't happen.");
  return NULL;
}

/* Function:  p7_oprofile_Clone()
 * Synopsis:  Allocate a cloned copy of  the optimized profile structure.  All
 *            allocated memory from the original profile is not reallocated.
 *            The cloned copy will point to the same memory as the original.
 * Incept:    SRE, Sun Nov 25 12:03:19 2007 [Casa de Gatos]
 *
 * Purpose:   Quick copy of an optimized profile used in mutiple threads.
 *
 * Throws:    <NULL> on allocation error.
 */

//This function does not require variants for each of the x86 SIMD ISAs
P7_OPROFILE *
p7_oprofile_Clone(const P7_OPROFILE *om1)
{
  int           status;

  P7_OPROFILE  *om2  = NULL;

  ESL_ALLOC(om2, sizeof(P7_OPROFILE));
  memcpy(om2, om1, sizeof(P7_OPROFILE));

  om2->clone  = 1;

  return om2;

 ERROR:
  p7_oprofile_Destroy(om2);
  return NULL;
}


/* Function:  p7_oprofile_UpdateFwdEmissionScores()
 * Synopsis:  Update the Forward/Backward part of the optimized profile
 *            match emissions to account for new background distribution.
 *
 * Purpose:   This implementation re-orders the loops used to access/modify
 *            the rfv array relative to how it's accessed for example in
 *            fb_conversion(), to minimize the required size of sc_arr.
 *
 * Args:      om              - optimized profile to be updated.
 *            bg              - the new bg distribution
 *            fwd_emissions   - precomputed Fwd (float) residue emission
 *                              probabilities in serial order (gathered from
 *                              the optimized striped <om> with
 *                              p7_oprofile_GetFwdEmissionArray() ).
 *            sc_arr            Preallocated array of at least Kp*4 floats
 */

static int p7_oprofile_UpdateFwdEmissionScores_Dispatcher(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);


int
(* p7_oprofile_UpdateFwdEmissionScores)(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr) = p7_oprofile_UpdateFwdEmissionScores_Dispatcher;

static int p7_oprofile_UpdateFwdEmissionScores_Dispatcher(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr){

#ifdef P7_TEST_ALL_SIMD
 p7_oprofile_UpdateFwdEmissionScores = p7_oprofile_UpdateFwdEmissionScores_test_all;
 return p7_oprofile_UpdateFwdEmissionScores_test_all(om, bg, fwd_emissions, sc_arr);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_UpdateFwdEmissionScores = p7_oprofile_UpdateFwdEmissionScores_test_sse_avx;
 return p7_oprofile_UpdateFwdEmissionScores_test_sse_avx(om, bg, fwd_emissions, sc_arr);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_UpdateFwdEmissionScores = p7_oprofile_UpdateFwdEmissionScores_avx512;
     return p7_oprofile_UpdateFwdEmissionScores_avx512(om, bg, fwd_emissions, sc_arr);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_UpdateFwdEmissionScores = p7_oprofile_UpdateFwdEmissionScores_avx;
     return p7_oprofile_UpdateFwdEmissionScores_avx(om, bg, fwd_emissions, sc_arr); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_UpdateFwdEmissionScores = p7_oprofile_UpdateFwdEmissionScores_sse;
      return p7_oprofile_UpdateFwdEmissionScores_sse(om, bg, fwd_emissions, sc_arr);
    }
#endif

  p7_Die("p7_oprofile_UpdateFwdEmissionScores_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;
}

/* Function:  p7_oprofile_UpdateVitEmissionScores()
 * Synopsis:  Update the Viterbi part of the optimized profile match
 *            emissions to account for new background distribution.
 *.
 * Purpose:   This implementation re-orders the loops used to access/modify
 *            the rmv array relative to how it's accessed for example in
 *            vf_conversion(), to minimize the required size of sc_arr.
 *
 * Args:      om              - optimized profile to be updated.
 *            bg              - the new bg distribution
 *            fwd_emissions   - precomputed Fwd (float) residue emission
 *                              probabilities in serial order (gathered from
 *                              the optimized striped <om> with
 *                              p7_oprofile_GetFwdEmissionArray() ).
 *            sc_arr            Preallocated array of at least Kp*8 floats
 */
static int p7_oprofile_UpdateVitEmissionScores_Dispatcher(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);

int (*p7_oprofile_UpdateVitEmissionScores)(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr) = p7_oprofile_UpdateVitEmissionScores_Dispatcher;

static int p7_oprofile_UpdateVitEmissionScores_Dispatcher(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr){

#ifdef P7_TEST_ALL_SIMD
 p7_oprofile_UpdateVitEmissionScores = p7_oprofile_UpdateVitEmissionScores_test_all;
 return p7_oprofile_UpdateVitEmissionScores_test_all(om, bg, fwd_emissions, sc_arr);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_UpdateVitEmissionScores = p7_oprofile_UpdateVitEmissionScores_test_sse_avx;
 return p7_oprofile_UpdateVitEmissionScores_test_sse_avx(om, bg, fwd_emissions, sc_arr);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_UpdateVitEmissionScores = p7_oprofile_UpdateVitEmissionScores_avx512;
     return p7_oprofile_UpdateVitEmissionScores_avx512(om, bg, fwd_emissions, sc_arr);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_UpdateVitEmissionScores = p7_oprofile_UpdateVitEmissionScores_avx;
     return p7_oprofile_UpdateVitEmissionScores_avx(om, bg, fwd_emissions, sc_arr); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_UpdateVitEmissionScores = p7_oprofile_UpdateVitEmissionScores_sse;
      return p7_oprofile_UpdateVitEmissionScores_sse(om, bg, fwd_emissions, sc_arr);
    }
#endif

  p7_Die("p7_oprofile_UpdateVitEmissionScores_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;
}


/* Function:  p7_oprofile_UpdateMSVEmissionScores()
 * Synopsis:  Update the MSV part of the optimized profile match
 *            emissions to account for new background distribution.
 *.
 * Purpose:   This implementation re-orders the loops used to access/modify
 *            the rbv array relative to how it's accessed for example in
 *            mf_conversion(), to minimize the required size of sc_arr.
 *
 * Args:      om              - optimized profile to be updated.
 *            bg              - the new bg distribution
 *            fwd_emissions   - precomputed Fwd (float) residue emission
 *                              probabilities in serial order (gathered from
 *                              the optimized striped <om> with
 *                              p7_oprofile_GetFwdEmissionArray() ).
 *            sc_arr            Preallocated array of at least Kp*16 floats
 */
static int
p7_oprofile_UpdateMSVEmissionScores_Dispatcher(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);

int
(*p7_oprofile_UpdateMSVEmissionScores)(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr) = p7_oprofile_UpdateMSVEmissionScores_Dispatcher;

static int p7_oprofile_UpdateMSVEmissionScores_Dispatcher(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr){

#ifdef P7_TEST_ALL_SIMD
 p7_oprofile_UpdateMSVEmissionScores = p7_oprofile_UpdateMSVEmissionScores_test_all;
 return p7_oprofile_UpdateMSVEmissionScores_test_all(om, bg, fwd_emissions, sc_arr);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_UpdateMSVEmissionScores = p7_oprofile_UpdateMSVEmissionScores_test_sse_avx;
 return p7_oprofile_UpdateMSVEmissionScores_test_sse_avx(om, bg, fwd_emissions, sc_arr);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_UpdateMSVEmissionScores = p7_oprofile_UpdateMSVEmissionScores_avx512;
     return p7_oprofile_UpdateMSVEmissionScores_avx512(om, bg, fwd_emissions, sc_arr);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_UpdateMSVEmissionScores = p7_oprofile_UpdateMSVEmissionScores_avx;
     return p7_oprofile_UpdateMSVEmissionScores_avx(om, bg, fwd_emissions, sc_arr); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_UpdateMSVEmissionScores = p7_oprofile_UpdateMSVEmissionScores_sse;
      return p7_oprofile_UpdateMSVEmissionScores_sse(om, bg, fwd_emissions, sc_arr);
    }
#endif

  p7_Die("p7_oprofile_UpdateMSVEmissionScores_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;
}

/*----------------- end, P7_OPROFILE structure ------------------*/



/*****************************************************************
 * 2. Conversion from generic P7_PROFILE to optimized P7_OPROFILE
 *****************************************************************/
/* unbiased_byteify()
 * Convert original transition score to a rounded uchar cost
 * Transition scores for MSVFilter get this treatment.
 * e.g. a score of -2.1, with scale 3.0, becomes a cost of 6.
 * (A cost of +255 is our -infinity "prohibited event")
 */
static uint8_t 
unbiased_byteify(P7_OPROFILE *om, float sc)
{
  uint8_t b;

  sc  = -1.0f * roundf(om->scale_b * sc);       /* ugh. sc is now an integer cost represented in a float...    */
  b   = (sc > 255.) ? 255 : (uint8_t) sc;	/* and now we cast and saturate it to an unsigned char cost... */
  return b;
}
 
/* wordify()
 * Converts log probability score to a rounded signed 16-bit integer cost.
 * Both emissions and transitions for ViterbiFilter get this treatment.
 * No bias term needed, because we use signed words. 
 *   e.g. a score of +3.2, with scale 500.0, becomes +1600.
 */
static int16_t 
wordify(P7_OPROFILE *om, float sc)
{
  sc  = roundf(om->scale_w * sc);
  if      (sc >=  32767.0) return  32767;
  else if (sc <= -32768.0) return -32768;
  else return (int16_t) sc;
}

/* Function:  p7_oprofile_Convert()
 * Synopsis:  Converts standard profile to an optimized one.
 * Incept:    SRE, Mon Nov 26 07:38:57 2007 [Janelia]
 *
 * Purpose:   Convert a standard profile <gm> to an optimized profile <om>,
 *            where <om> has already been allocated for a profile of at 
 *            least <gm->M> nodes and the same emission alphabet <gm->abc>.
 *
 * Args:      gm - profile to optimize
 *            om - allocated optimized profile for holding the result.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <gm>, <om> aren't compatible. 
 *            <eslEMEM> on allocation failure.
 */

static int p7_oprofile_Convert_Dispatcher(const P7_PROFILE *gm, P7_OPROFILE *om);

int (*p7_oprofile_Convert)(const P7_PROFILE *gm, P7_OPROFILE *om) = p7_oprofile_Convert_Dispatcher;

static int p7_oprofile_Convert_Dispatcher(const P7_PROFILE *gm, P7_OPROFILE *om){


#ifdef P7_TEST_ALL_SIMD
 p7_oprofile_Convert = p7_oprofile_Convert_test_all;
 return p7_oprofile_Convert_test_all(gm, om);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_Convert = p7_oprofile_Convert_test_sse_avx;
 return p7_oprofile_Convert_test_sse_avx(gm, om);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_Convert = p7_oprofile_Convert_avx512;
     return p7_oprofile_Convert_avx512(gm, om);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_Convert = p7_oprofile_Convert_avx;
     return p7_oprofile_Convert_avx(gm, om); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_Convert = p7_oprofile_Convert_sse;
      return p7_oprofile_Convert_sse(gm, om);
    }
#endif

  p7_Die("p7_oprofile_Convert_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;
}


/* Function:  p7_oprofile_ReconfigLength()
 * Synopsis:  Set the target sequence length of a model.
 * Incept:    SRE, Thu Dec 20 09:56:40 2007 [Janelia]
 *
 * Purpose:   Given an already configured model <om>, quickly reset its
 *            expected length distribution for a new mean target sequence
 *            length of <L>. 
 *            
 *            This doesn't affect the length distribution of the null
 *            model. That must also be reset, using <p7_bg_SetLength()>.
 *            
 *            We want this routine to run as fast as possible, because
 *            this call is in the critical path: it must be called at
 *            each new target sequence in a database search.
 *
 * Returns:   <eslOK> on success. Costs/scores for N,C,J transitions are set
 *            here.
 */
int
p7_oprofile_ReconfigLength(P7_OPROFILE *om, int L)
{
  int status;
  if ((status = p7_oprofile_ReconfigMSVLength (om, L)) != eslOK) return status;
  if ((status = p7_oprofile_ReconfigRestLength(om, L)) != eslOK) return status;
  return eslOK;
}

/* Function:  p7_oprofile_ReconfigMSVLength()
 * Synopsis:  Set the target sequence length of the MSVFilter part of the model.
 * Incept:    SRE, Tue Dec 16 13:39:17 2008 [Janelia]
 *
 * Purpose:   Given an  already configured model <om>, quickly reset its
 *            expected length distribution for a new mean target sequence
 *            length of <L>, only for the part of the model that's used
 *            for the accelerated MSV filter.
 *            
 *            The acceleration pipeline uses this to defer reconfiguring the
 *            length distribution of the main model, mostly because hmmscan
 *            reads the model in two pieces, MSV part first, then the rest.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_oprofile_ReconfigMSVLength(P7_OPROFILE *om, int L)
{
  om->tjb_b = unbiased_byteify(om, logf(3.0f / (float) (L+3)));
  return eslOK;
}

/* Function:  p7_oprofile_ReconfigRestLength()
 * Synopsis:  Set the target sequence length of the main profile.
 * Incept:    SRE, Tue Dec 16 13:41:30 2008 [Janelia]
 *
 * Purpose:   Given an  already configured model <om>, quickly reset its
 *            expected length distribution for a new mean target sequence
 *            length of <L>, for everything except the MSV filter part
 *            of the model.
 *            
 *            Calling <p7_oprofile_ReconfigMSVLength()> then
 *            <p7_oprofile_ReconfigRestLength()> is equivalent to
 *            just calling <p7_oprofile_ReconfigLength()>. The two
 *            part version is used in the acceleration pipeline.
 *
 * Returns:   <eslOK> on success.           
 */
int
p7_oprofile_ReconfigRestLength(P7_OPROFILE *om, int L)
{
  float pmove, ploop;
  
  pmove = (2.0f + om->nj) / ((float) L + 2.0f + om->nj); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1.0f - pmove;

  /* ForwardFilter() parameters: pspace floats */
  om->xf[p7O_N][p7O_LOOP] =  om->xf[p7O_C][p7O_LOOP] = om->xf[p7O_J][p7O_LOOP] = ploop;
  om->xf[p7O_N][p7O_MOVE] =  om->xf[p7O_C][p7O_MOVE] = om->xf[p7O_J][p7O_MOVE] = pmove;

  /* ViterbiFilter() parameters: lspace signed 16-bit ints */
  om->xw[p7O_N][p7O_MOVE] =  om->xw[p7O_C][p7O_MOVE] = om->xw[p7O_J][p7O_MOVE] = wordify(om, logf(pmove));
  /* om->xw[p7O_N][p7O_LOOP] =  om->xw[p7O_C][p7O_LOOP] = om->xw[p7O_J][p7O_LOOP] = wordify(om, logf(ploop)); */ /* 3nat approx in force: these stay 0 */
  /* om->ncj_roundoff        = (om->scale_w * logf(ploop)) - om->xw[p7O_N][p7O_LOOP];                         */ /* and this does too                  */

  om->L = L;
  return eslOK;
}


/* Function:  p7_oprofile_ReconfigMultihit()
 * Synopsis:  Quickly reconfig model into multihit mode for target length <L>.
 * Incept:    SRE, Thu Aug 21 10:04:07 2008 [Janelia]
 *
 * Purpose:   Given a profile <om> that's already been configured once,
 *            quickly reconfigure it into a multihit mode for target 
 *            length <L>. 
 *            
 *            This gets called in domain definition, when we need to
 *            flip the model in and out of unihit mode to
 *            process individual domains.
 *            
 * Note:      You can't just flip uni/multi mode alone, because that
 *            parameterization also affects target length
 *            modeling. You need to make sure uni vs. multi choice is
 *            made before the length model is set, and you need to
 *            make sure the length model is recalculated if you change
 *            the uni/multi mode. Hence, these functions call
 *            <p7_oprofile_ReconfigLength()>.
 */
int
p7_oprofile_ReconfigMultihit(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = 0.5;
  om->xf[p7O_E][p7O_LOOP] = 0.5;
  om->nj = 1.0f;

  om->xw[p7O_E][p7O_MOVE] = wordify(om, -eslCONST_LOG2);
  om->xw[p7O_E][p7O_LOOP] = wordify(om, -eslCONST_LOG2);

  return p7_oprofile_ReconfigLength(om, L);
}

/* Function:  p7_oprofile_ReconfigUnihit()
 * Synopsis:  Quickly reconfig model into unihit mode for target length <L>.
 * Incept:    SRE, Thu Aug 21 10:10:32 2008 [Janelia]
 *
 * Purpose:   Given a profile <om> that's already been configured once,
 *            quickly reconfigure it into a unihit mode for target 
 *            length <L>. 
 *            
 *            This gets called in domain definition, when we need to
 *            flip the model in and out of unihit <L=0> mode to
 *            process individual domains.
 */
int
p7_oprofile_ReconfigUnihit(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = 1.0f;
  om->xf[p7O_E][p7O_LOOP] = 0.0f;
  om->nj = 0.0f;

  om->xw[p7O_E][p7O_MOVE] = 0;
  om->xw[p7O_E][p7O_LOOP] = -32768;

  return p7_oprofile_ReconfigLength(om, L);
}
/*------------ end, conversions to P7_OPROFILE ------------------*/

/*******************************************************************
*   3. Conversion from optimized P7_OPROFILE to compact score arrays
 *******************************************************************/

/* Function:  p7_oprofile_GetFwdTransitionArray()
 * Synopsis:  Retrieve full 32-bit float transition probabilities from an
 *            optimized profile into a flat array
 *
 * Purpose:   Extract an array of <type> (e.g. p7O_II) transition probabilities
 *            from the underlying <om> profile. In SIMD implementations,
 *            these are striped and interleaved, making them difficult to
 *            directly access.
 *
 * Args:      <om>   - optimized profile, containing transition information
 *            <type> - transition type (e.g. p7O_II)
 *            <arr>  - preallocated array into which floats will be placed
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
static int
p7_oprofile_GetFwdTransitionArray_Dispatcher(const P7_OPROFILE *om, int type, float *arr );

int
(* p7_oprofile_GetFwdTransitionArray)(const P7_OPROFILE *om, int type, float *arr ) = p7_oprofile_GetFwdTransitionArray_Dispatcher;

static int
p7_oprofile_GetFwdTransitionArray_Dispatcher(const P7_OPROFILE *om, int type, float *arr ){

  #ifdef P7_TEST_ALL_SIMD
 p7_oprofile_GetFwdTransitionArray = p7_oprofile_GetFwdTransitionArray_test_all;
 return p7_oprofile_GetFwdTransitionArray_test_all(om, type, arr);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_GetFwdTransitionArray = p7_oprofile_GetFwdTransitionArray_test_sse_avx;
 return p7_oprofile_GetFwdTransitionArray_test_sse_avx(om, type, arr);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_GetFwdTransitionArray = p7_oprofile_GetFwdTransitionArray_avx512;
     return p7_oprofile_GetFwdTransitionArray_avx512(om, type, arr);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_GetFwdTransitionArray = p7_oprofile_GetFwdTransitionArray_avx;
     return p7_oprofile_GetFwdTransitionArray_avx(om, type, arr); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_GetFwdTransitionArray = p7_oprofile_GetFwdTransitionArray_sse;
      return p7_oprofile_GetFwdTransitionArray_sse(om, type, arr);
    }
#endif

  p7_Die("p7_oprofile_GetFwdTransitionArray_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;
}


/* Function:  p7_oprofile_GetSSVEmissionScoreArray()
 * Synopsis:  Retrieve MSV residue emission scores from an optimized
 *            profile into an array
 *
 * Purpose:   Extract an implicitly 2D array of 8-bit int SSV residue
 *            emission scores from an optimized profile <om>. <arr> must
 *            be allocated by the calling function to be of size
 *            ( om->abc->Kp * ( om->M  + 1 )), and indexing into the array
 *            is done as  [om->abc->Kp * i +  c ] for character c at
 *            position i.
 *
 *            In SIMD implementations, the residue scores are striped
 *            and interleaved, making them somewhat difficult to
 *            directly access. Faster access is desired, for example,
 *            in SSV back-tracking of a high-scoring diagonal
 *
 * Args:      <om>   - optimized profile, containing transition information
 *            <arr>  - preallocated array into which scores will be placed
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */

static int
p7_oprofile_GetSSVEmissionScoreArray_Dispatcher(const P7_OPROFILE *om, uint8_t *arr );

int (* p7_oprofile_GetSSVEmissionScoreArray)(const P7_OPROFILE *om, uint8_t *arr ) = p7_oprofile_GetSSVEmissionScoreArray_Dispatcher;

static int
p7_oprofile_GetSSVEmissionScoreArray_Dispatcher(const P7_OPROFILE *om, uint8_t *arr )
{
  
  #ifdef P7_TEST_ALL_SIMD
 p7_oprofile_GetSSVEmissionScoreArray = p7_oprofile_GetSSVEmissionScoreArray_test_all;
 return p7_oprofile_GetSSVEmissionScoreArray_test_all(om,  arr);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_GetSSVEmissionScoreArray = p7_oprofile_GetSSVEmissionScoreArray_test_sse_avx;
 return p7_oprofile_GetSSVEmissionScoreArray_test_sse_avx(om,  arr);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_GetSSVEmissionScoreArray = p7_oprofile_GetSSVEmissionScoreArray_avx512;
     return p7_oprofile_GetSSVEmissionScoreArray_avx512(om, arr);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_GetSSVEmissionScoreArray = p7_oprofile_GetSSVEmissionScoreArray_avx;
     return p7_oprofile_GetSSVEmissionScoreArray_avx(om,  arr); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_GetSSVEmissionScoreArray = p7_oprofile_GetSSVEmissionScoreArray_sse;
      return p7_oprofile_GetSSVEmissionScoreArray_sse(om, arr);
    }
#endif

  p7_Die("p7_oprofile_GetSSVEmissionScoreArray_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;
}


/* Function:  p7_oprofile_GetFwdEmissionScoreArray()
 * Synopsis:  Retrieve Fwd (float) residue emission scores from an optimized
 *            profile into an array
 *
 * Purpose:   Extract an implicitly 2D array of 32-bit float Fwd residue
 *            emission scores from an optimized profile <om>. <arr> must
 *            be allocated by the calling function to be of size
 *            ( om->abc->Kp * ( om->M  + 1 )), and indexing into the array
 *            is done as  [om->abc->Kp * i +  c ] for character c at
 *            position i.
 *
 *            In SIMD implementations, the residue scores are striped
 *            and interleaved, making them somewhat difficult to
 *            directly access.
 *
 * Args:      <om>   - optimized profile, containing transition information
 *            <arr>  - preallocated array into which scores will be placed
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Note:      [SRE 2024/0107-h3-iss320] This is inelegant, if not bugged.
 *            It's accessing slots in vectors that aren't valid scores,
 *            because 4*nq >= K, so you can't just take all the slots
 *            in the striped vectors and try to store them somewhere.
 *            I think it's only working because of order of operations:
 *            it improperly stores a few values then overwrites them.
 */
static int
p7_oprofile_GetFwdEmissionScoreArray_Dispatcher(const P7_OPROFILE *om, float *arr );

int (* p7_oprofile_GetFwdEmissionScoreArray)(const P7_OPROFILE *om, float *arr ) = p7_oprofile_GetFwdEmissionScoreArray_Dispatcher;

static int
p7_oprofile_GetFwdEmissionScoreArray_Dispatcher(const P7_OPROFILE *om, float *arr )
{
  
  #ifdef P7_TEST_ALL_SIMD
 p7_oprofile_GetFwdEmissionScoreArray = p7_oprofile_GetFwdEmissionScoreArray_test_all;
 return p7_oprofile_GetFwdEmissionScoreArray_test_all(om,  arr);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_GetFwdEmissionScoreArray = p7_oprofile_GetFwdEmissionScoreArray_test_sse_avx;
 return p7_oprofile_GetFwdEmissionScoreArray_test_sse_avx(om,  arr);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_GetFwdEmissionScoreArray = p7_oprofile_GetFwdEmissionScoreArray_avx512;
     return p7_oprofile_GetFwdEmissionScoreArray_avx512(om, arr);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_GetFwdEmissionScoreArray = p7_oprofile_GetFwdEmissionScoreArray_avx;
     return p7_oprofile_GetFwdEmissionScoreArray_avx(om,  arr); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_GetFwdEmissionScoreArray = p7_oprofile_GetFwdEmissionScoreArray_sse;
      return p7_oprofile_GetFwdEmissionScoreArray_sse(om, arr);
    }
#endif

  p7_Die("p7_oprofile_GetFwdEmissionScoreArray_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;
}


/* Function:  p7_oprofile_GetFwdEmissionArray()
 * Synopsis:  Retrieve Fwd (float) residue emission values from an optimized
 *            profile into an array
 *
 * Purpose:   Extract an implicitly 2D array of 32-bit float Fwd residue
 *            emission values from an optimized profile <om>, converting
 *            back to emission values based on the background. <arr> must
 *            be allocated by the calling function to be of size
 *            ( om->abc->Kp * ( om->M  + 1 )), and indexing into the array
 *            is done as  [om->abc->Kp * i +  c ] for character c at
 *            position i.
 *
 *            In SIMD implementations, the residue scores are striped
 *            and interleaved, making them somewhat difficult to
 *            directly access.
 *
 * Args:      <om>   - optimized profile, containing transition information
 *            <bg>   - background frequencies
 *            <arr>  - preallocated array into which scores will be placed
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Note:      see above comment for p7_oprofile_GetFwdEmissionScoreArray.
 *            This function also appears to be problematic.
 *            (Why are we even doing this. If we want unstriped probabilities,
 *            use the HMM, not the vectorized profile, right?)
 */

static int
p7_oprofile_GetFwdEmissionArray_Dispatcher(const P7_OPROFILE *om, P7_BG *bg, float *arr );

int (* p7_oprofile_GetFwdEmissionArray)(const P7_OPROFILE *om, P7_BG *bg, float *arr ) = p7_oprofile_GetFwdEmissionArray_Dispatcher;

static int
p7_oprofile_GetFwdEmissionArray_Dispatcher(const P7_OPROFILE *om, P7_BG *bg, float *arr )
{
  
  #ifdef P7_TEST_ALL_SIMD
 p7_oprofile_GetFwdEmissionArray = p7_oprofile_GetFwdEmissionArray_test_all;
 return p7_oprofile_GetFwdEmissionArray_test_all(om, bg, arr);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_GetFwdEmissionArray = p7_oprofile_GetFwdEmissionArray_test_sse_avx;
 return p7_oprofile_GetFwdEmissionArray_test_sse_avx(om, bg, arr);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_GetFwdEmissionArray = p7_oprofile_GetFwdEmissionArray_avx512;
     return p7_oprofile_GetFwdEmissionArray_avx512(om, bg, arr);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_GetFwdEmissionArray = p7_oprofile_GetFwdEmissionArray_avx;
     return p7_oprofile_GetFwdEmissionArray_avx(om, bg, arr); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_GetFwdEmissionArray = p7_oprofile_GetFwdEmissionArray_sse;
      return p7_oprofile_GetFwdEmissionArray_sse(om, bg, arr);
    }
#endif

  p7_Die("p7_oprofile_GetFwdEmissionArray_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;
}



/*------------ end, conversions from P7_OPROFILE ------------------*/


/*****************************************************************
 * 4. Debugging and development utilities.
 *****************************************************************/


/* oprofile_dump_mf()
 * 
 * Dump the MSVFilter part of a profile <om> to <stdout>.
 */
static int
oprofile_dump_mf(FILE *fp, const P7_OPROFILE *om)
{
  int     M   = om->M;		/* length of the query                                          */
  int     nq  = p7O_NQB(M);     /* segment length; total # of striped vectors needed            */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* counter over nodes 1..M                                      */
  int     z;			/* counter within elements of one SIMD minivector               */
  union { __m128i v; uint8_t i[16]; } tmp; /* used to align and read simd minivectors           */

  /* Header (rearranged column numbers, in the vectors)  */
  fprintf(fp, "     ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < 16; z++) 
	if (k+z*nq <= M) fprintf(fp, "%4d ", k+z*nq);
	else             fprintf(fp, "%4s ", "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");

  /* Table of residue emissions */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 

      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->rbv[x][q]);
	  for (z = 0; z < 16; z++) fprintf(fp, "%4d ", tmp.i[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");
    }
  fprintf(fp, "\n");

  fprintf(fp, "t_EC,EJ:    %4d\n",  om->tec_b);
  fprintf(fp, "t_NB,JB,CT: %4d\n",  om->tjb_b);
  fprintf(fp, "t_BMk:      %4d\n",  om->tbm_b);
  fprintf(fp, "scale:      %.2f\n", om->scale_b);
  fprintf(fp, "base:       %4d\n",  om->base_b);
  fprintf(fp, "bias:       %4d\n",  om->bias_b);
  fprintf(fp, "Q:          %4d\n",  nq);  
  fprintf(fp, "M:          %4d\n",  M);  
  return eslOK;
}



/* oprofile_dump_vf()
 * 
 * Dump the ViterbiFilter part of a profile <om> to <stdout>.
 */
static int
oprofile_dump_vf(FILE *fp, const P7_OPROFILE *om)
{
  int     M   = om->M;		/* length of the query                                          */
  int     nq  = p7O_NQW(M);     /* segment length; total # of striped vectors needed            */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = p7O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     j;			/* counter in interleaved vector arrays in the profile          */
  union { __m128i v; int16_t i[8]; } tmp; /* used to align and read simd minivectors           */

  /* Emission score header (rearranged column numbers, in the vectors)  */
  fprintf(fp, "     ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < 8; z++) 
	if (k+z*nq <= M) fprintf(fp, "%6d ", k+z*nq);
	else             fprintf(fp, "%6s ", "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");

  /* Table of residue emissions */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 

      /* Match emission scores (insert emissions are assumed zero by design) */
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->rwv[x][q]);
	  for (z = 0; z < 8; z++) fprintf(fp, "%6d ", tmp.i[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");
    }
  fprintf(fp, "\n");

  /* Transitions */
  for (t = p7O_BM; t <= p7O_II; t++)
    {
      switch (t) {
      case p7O_BM: fprintf(fp, "\ntBM: "); break;
      case p7O_MM: fprintf(fp, "\ntMM: "); break;
      case p7O_IM: fprintf(fp, "\ntIM: "); break;
      case p7O_DM: fprintf(fp, "\ntDM: "); break;
      case p7O_MD: fprintf(fp, "\ntMD: "); break;
      case p7O_MI: fprintf(fp, "\ntMI: "); break;
      case p7O_II: fprintf(fp, "\ntII: "); break;
      }

      for (k = 1, q = 0; q < nq; q++, k++)
	{
	  switch (t) {
	  case p7O_BM: kb = k;                 break; 
	  case p7O_MM: kb = (1 + (nq+k-2)) % nq; break; /* MM, DM, IM quads rotated by +1  */
	  case p7O_IM: kb = (1 + (nq+k-2)) % nq; break;  
	  case p7O_DM: kb = (1 + (nq+k-2)) % nq; break;  
	  case p7O_MD: kb = k;                 break; /* the remaining ones are straight up  */
	  case p7O_MI: kb = k;                 break; 
	  case p7O_II: kb = k;                 break; 
	  }
	  fprintf(fp, "[ ");
	  for (z = 0; z < 8; z++) 
	    if (kb+z*nq <= M) fprintf(fp, "%6d ", kb+z*nq);
	    else              fprintf(fp, "%6s ", "xx");
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n     ");	  
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  _mm_store_si128(&tmp.v, om->twv[q*7 + t]);
	  for (z = 0; z < 8; z++) fprintf(fp, "%6d ", tmp.i[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");	  
    }

  /* DD transitions */
  fprintf(fp, "\ntDD: ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < 8; z++) 
	if (k+z*nq <= M) fprintf(fp, "%6d ", k+z*nq);
	else             fprintf(fp, "%6s ", "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n     ");	  
  for (j = nq*7, q = 0; q < nq; q++, j++)
    {
      fprintf(fp, "[ ");
      _mm_store_si128(&tmp.v, om->twv[j]);
      for (z = 0; z < 8; z++) fprintf(fp, "%6d ", tmp.i[z]);
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");	  

  fprintf(fp, "E->C: %6d    E->J: %6d\n", om->xw[p7O_E][p7O_MOVE], om->xw[p7O_E][p7O_LOOP]);
  fprintf(fp, "N->B: %6d    N->N: %6d\n", om->xw[p7O_N][p7O_MOVE], om->xw[p7O_N][p7O_LOOP]);
  fprintf(fp, "J->B: %6d    J->J: %6d\n", om->xw[p7O_J][p7O_MOVE], om->xw[p7O_J][p7O_LOOP]);
  fprintf(fp, "C->T: %6d    C->C: %6d\n", om->xw[p7O_C][p7O_MOVE], om->xw[p7O_C][p7O_LOOP]);

  fprintf(fp, "scale: %6.2f\n", om->scale_w);
  fprintf(fp, "base:  %6d\n",   om->base_w);
  fprintf(fp, "bound: %6d\n",   om->ddbound_w);
  fprintf(fp, "Q:     %6d\n",   nq);  
  fprintf(fp, "M:     %6d\n",   M);  
  return eslOK;
}


/* oprofile_dump_fb()
 * 
 * Dump the Forward/Backward part of a profile <om> to <stdout>.
 * <width>, <precision> control the floating point output:
 *  8,5 is a reasonable choice for prob space,
 *  5,2 is reasonable for log space.
 */
static int
oprofile_dump_fb(FILE *fp, const P7_OPROFILE *om, int width, int precision)
{
  int     M   = om->M;		/* length of the query                                          */
  int     nq  = p7O_NQF(M);     /* segment length; total # of striped vectors needed            */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = p7O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     j;			/* counter in interleaved vector arrays in the profile          */
  union { __m128 v; float x[4]; } tmp; /* used to align and read simd minivectors               */

  /* Residue emissions */
  for (x = 0; x < om->abc->Kp; x++)
    {
      fprintf(fp, "(%c): ", om->abc->sym[x]); 
      for (k =1, q = 0; q < nq; q++, k++)
	{
	  fprintf(fp, "[ ");
	  for (z = 0; z < 4; z++) 
	    if (k+z*nq <= M) fprintf(fp, "%*d ", width, k+z*nq);
	    else             fprintf(fp, "%*s ", width, "xx");
	  fprintf(fp, "]");
	}
      fprintf(fp, "\nmat: ");
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  tmp.v = om->rfv[x][q];
	  for (z = 0; z < 4; z++) fprintf(fp, "%*.*f ", width, precision, tmp.x[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n\n");
    }

  /* Transitions */
  for (t = p7O_BM; t <= p7O_II; t++)
    {
      switch (t) {
      case p7O_BM: fprintf(fp, "\ntBM: "); break;
      case p7O_MM: fprintf(fp, "\ntMM: "); break;
      case p7O_IM: fprintf(fp, "\ntIM: "); break;
      case p7O_DM: fprintf(fp, "\ntDM: "); break;
      case p7O_MD: fprintf(fp, "\ntMD: "); break;
      case p7O_MI: fprintf(fp, "\ntMI: "); break;
      case p7O_II: fprintf(fp, "\ntII: "); break;
      }
      for (k = 1, q = 0; q < nq; q++, k++)
	{
	  switch (t) {
	  case p7O_MM:/* MM, DM, IM quads rotated by +1  */
	  case p7O_IM:
	  case p7O_DM:
		  kb = (1 + (nq+k-2)) % nq;
		  break;
	  case p7O_BM:/* the remaining ones are straight up  */
	  case p7O_MD:
	  case p7O_MI:
	  case p7O_II:
		  kb = k;
		  break;
	  }
	  fprintf(fp, "[ ");
	  for (z = 0; z < 4; z++) 
	    if (kb+z*nq <= M) fprintf(fp, "%*d ", width, kb+z*nq);
	    else              fprintf(fp, "%*s ", width, "xx");
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n     ");	  
      for (q = 0; q < nq; q++)
	{
	  fprintf(fp, "[ ");
	  tmp.v = om->tfv[q*7 + t];
	  for (z = 0; z < 4; z++) fprintf(fp, "%*.*f ", width, precision, tmp.x[z]);
	  fprintf(fp, "]");
	}
      fprintf(fp, "\n");	  
    }

  /* DD transitions */
  fprintf(fp, "\ntDD: ");
  for (k =1, q = 0; q < nq; q++, k++)
    {
      fprintf(fp, "[ ");
      for (z = 0; z < 4; z++) 
	if (k+z*nq <= M) fprintf(fp, "%*d ", width, k+z*nq);
	else             fprintf(fp, "%*s ", width, "xx");
      fprintf(fp, "]");
    }
  fprintf(fp, "\n     ");	  
  for (j = nq*7, q = 0; q < nq; q++, j++)
    {
      fprintf(fp, "[ ");
      tmp.v = om->tfv[j];
      for (z = 0; z < 4; z++) fprintf(fp, "%*.*f ", width, precision, tmp.x[z]);
      fprintf(fp, "]");
    }
  fprintf(fp, "\n");	  
  
  /* Specials */
  fprintf(fp, "E->C: %*.*f    E->J: %*.*f\n", width, precision, om->xf[p7O_E][p7O_MOVE], width, precision, om->xf[p7O_E][p7O_LOOP]);
  fprintf(fp, "N->B: %*.*f    N->N: %*.*f\n", width, precision, om->xf[p7O_N][p7O_MOVE], width, precision, om->xf[p7O_N][p7O_LOOP]);
  fprintf(fp, "J->B: %*.*f    J->J: %*.*f\n", width, precision, om->xf[p7O_J][p7O_MOVE], width, precision, om->xf[p7O_J][p7O_LOOP]);
  fprintf(fp, "C->T: %*.*f    C->C: %*.*f\n", width, precision, om->xf[p7O_C][p7O_MOVE], width, precision, om->xf[p7O_C][p7O_LOOP]);
  fprintf(fp, "Q:     %d\n",   nq);  
  fprintf(fp, "M:     %d\n",   M);  
  return eslOK;
}


/* Function:  p7_oprofile_Dump()
 * Synopsis:  Dump internals of a <P7_OPROFILE>
 * Incept:    SRE, Thu Dec 13 08:49:30 2007 [Janelia]
 *
 * Purpose:   Dump the internals of <P7_OPROFILE> structure <om>
 *            to stream <fp>; generally for testing or debugging
 *            purposes.
 *
 * Args:      fp   - output stream (often stdout)
 *            om   - optimized profile to dump
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om)
{
  int status;

  fprintf(fp, "Dump of a <P7_OPROFILE> ::\n");

  fprintf(fp, "\n  -- float part, odds ratios for Forward/Backward:\n");
  if ((status = oprofile_dump_fb(fp, om, 8, 5)) != eslOK) return status;

  fprintf(fp, "\n  -- sword part, log odds for ViterbiFilter(): \n");
  if ((status = oprofile_dump_vf(fp, om))       != eslOK) return status;

  fprintf(fp, "\n  -- uchar part, log odds for MSVFilter(): \n");
  if ((status = oprofile_dump_mf(fp, om))       != eslOK) return status;

  return eslOK;
}


/* Function:  p7_oprofile_Sample()
 * Synopsis:  Sample a random profile.
 * Incept:    SRE, Wed Jul 30 13:11:52 2008 [Janelia]
 *
 * Purpose:   Sample a random profile of <M> nodes for alphabet <abc>,
 *            using <r> as the source of random numbers. Parameterize
 *            it for generation of target sequences of mean length
 *            <L>. Calculate its log-odds scores using background
 *            model <bg>.
 *            
 * Args:      r       - random number generator
 *            abc     - emission alphabet 
 *            bg      - background frequency model
 *            M       - size of sampled profile, in nodes
 *            L       - configured target seq mean length
 *            opt_hmm - optRETURN: sampled HMM
 *            opt_gm  - optRETURN: sampled normal profile
 *            opt_om  - RETURN: optimized profile
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
		   P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om)
{
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_OPROFILE    *om   = NULL;
  int             status;

  if ((gm = p7_profile_Create (M, abc)) == NULL)  { status = eslEMEM; goto ERROR; }
  if ((om = p7_oprofile_Create(M, abc)) == NULL)  { status = eslEMEM; goto ERROR; }

  if ((status = p7_hmm_Sample(r, M, abc, &hmm))             != eslOK) goto ERROR;
  if ((status = p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL)) != eslOK) goto ERROR;
  if ((status = p7_oprofile_Convert(gm, om))                != eslOK) goto ERROR;
  if ((status = p7_oprofile_ReconfigLength(om, L))          != eslOK) goto ERROR;

  if (opt_hmm != NULL) *opt_hmm = hmm; else p7_hmm_Destroy(hmm);
  if (opt_gm  != NULL) *opt_gm  = gm;  else p7_profile_Destroy(gm);
  *ret_om = om;
  return eslOK;

 ERROR:
  if (opt_hmm != NULL) *opt_hmm = NULL;
  if (opt_gm  != NULL) *opt_gm  = NULL;
  *ret_om = NULL;
  return status;
}


/* Function:  p7_oprofile_Compare()
 * Synopsis:  Compare two optimized profiles for equality.
 * Incept:    SRE, Wed Jan 21 13:29:10 2009 [Janelia]
 *
 * Purpose:   Compare the contents of <om1> and <om2>; return 
 *            <eslOK> if they are effectively identical profiles,
 *            or <eslFAIL> if not.
 * 
 *            Floating point comparisons are done to a tolerance
 *            of <tol> using <esl_FCompare_old()>.
 *            
 *            If a comparison fails, an informative error message is
 *            left in <errmsg> to indicate why.
 *            
 *            Internal allocation sizes are not compared, only the
 *            data.
 *            
 * Args:      om1    - one optimized profile to compare
 *            om2    - the other
 *            tol    - floating point comparison tolerance; see <esl_FCompare_old()>
 *            errmsg - ptr to array of at least <eslERRBUFSIZE> characters.
 *            
 * Returns:   <eslOK> on effective equality;  <eslFAIL> on difference.
 */

static int p7_oprofile_compare_Dispatcher(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);

int
(*p7_oprofile_Compare)(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg) = p7_oprofile_compare_Dispatcher;

static int p7_oprofile_compare_Dispatcher(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg){
  
#ifdef P7_TEST_ALL_SIMD
 p7_oprofile_Compare = p7_oprofile_Compare_test_all;
 return p7_oprofile_Compare_test_all(om1, om2, tol, errmsg);
 
#endif

#ifdef P7_TEST_SSE_AVX
 p7_oprofile_Compare = p7_oprofile_Compare_test_sse_avx;
 return p7_oprofile_Compare_test_sse_avx(om1, om2, tol, errmsg);  
#endif

#ifdef eslENABLE_AVX512  // Fastest first.
 if (esl_cpu_has_avx512())
   {
     p7_oprofile_Compare = p7_oprofile_Compare_avx512;
     return p7_oprofile_Compare_avx512(om1, om2, tol, errmsg);
   }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
     p7_oprofile_Compare = p7_oprofile_Compare_avx;
     return p7_oprofile_Compare_avx(om1, om2, tol, errmsg); 
    }
#endif

#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse4())
    {
      p7_oprofile_Compare = p7_oprofile_Compare_sse;
      return p7_oprofile_Compare_sse(om1, om2, tol, errmsg);
    }
#endif

  p7_Die("p7_oprofile_Compare_dispatcher found no vector implementation - that shouldn't happen.");
  return eslFAIL;

}

/* Function:  p7_profile_SameAsMF()
 * Synopsis:  Set a generic profile's scores to give MSV scores.
 * Incept:    SRE, Wed Jul 30 13:42:49 2008 [Janelia]
 *
 * Purpose:   Set a generic profile's scores so that the normal <dp_generic> DP 
 *            algorithms will give the same score as <p7_MSVFilter()>:
 *            all t_MM scores = 0; all other core transitions = -inf;
 *            multihit local mode; all <t_BMk> entries uniformly <log 2/(M(M+1))>;
 *            <tCC, tNN, tJJ> scores 0; total approximated later as -3;
 *            rounded in the same way as the 8-bit limited precision.
 *
 * Returns:   <eslOK> on success.
 */
// Does not require separate versions for different SIMD ISAs
int
p7_profile_SameAsMF(const P7_OPROFILE *om, P7_PROFILE *gm)
{
  int    k,x;
  float  tbm = roundf(om->scale_b * (log(2.0f / ((float) gm->M * (float) (gm->M+1)))));

  /* Transitions */
  esl_vec_FSet(gm->tsc, p7P_NTRANS * gm->M, -eslINFINITY);
  for (k = 1; k <  gm->M; k++) p7P_TSC(gm, k, p7P_MM) = 0.0f;
  for (k = 0; k <  gm->M; k++) p7P_TSC(gm, k, p7P_BM) = tbm;
  
  /* Emissions */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 0; k <= gm->M; k++)
      {
	gm->rsc[x][k*2]   = (gm->rsc[x][k*2] <= -eslINFINITY) ? -eslINFINITY : roundf(om->scale_b * gm->rsc[x][k*2]);
	gm->rsc[x][k*2+1] = 0;	/* insert score: VF makes it zero no matter what. */
      }	

   /* Specials */
  for (k = 0; k < p7P_NXSTATES; k++)
    for (x = 0; x < p7P_NXTRANS; x++)
      gm->xsc[k][x] = (gm->xsc[k][x] <= -eslINFINITY) ? -eslINFINITY : roundf(om->scale_b * gm->xsc[k][x]);

  /* NN, CC, JJ hardcoded 0 in limited precision */
  gm->xsc[p7P_N][p7P_LOOP] =  gm->xsc[p7P_J][p7P_LOOP] =  gm->xsc[p7P_C][p7P_LOOP] = 0;

  return eslOK;
}


/* Function:  p7_profile_SameAsVF()
 * Synopsis:  Round a generic profile to match ViterbiFilter scores.
 * Incept:    SRE, Wed Jul 30 13:37:48 2008 [Janelia]
 *
 * Purpose:   Round all the scores in a generic (lspace) <P7_PROFILE> <gm> in
 *            exactly the same way that the scores in the
 *            <P7_OPROFILE> <om> were rounded. Then we can test that two profiles
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
 * Args:      <om>  - optimized profile, containing scale information.
 *            <gm>  - generic profile that <om> was built from.          
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
// Does not require separate versions for different SIMD ISAs
int
p7_profile_SameAsVF(const P7_OPROFILE *om, P7_PROFILE *gm)
{
  int k;
  int x;

  /* Transitions */
  /* <= -eslINFINITY test is used solely to silence compiler. really testing == -eslINFINITY */
  for (x = 0; x < gm->M*p7P_NTRANS; x++)
    gm->tsc[x] = (gm->tsc[x] <= -eslINFINITY) ? -eslINFINITY : roundf(om->scale_w * gm->tsc[x]);
  
  /* Enforce the rule that no II can be 0; max of -1 */
  for (x = p7P_II; x < gm->M*p7P_NTRANS; x += p7P_NTRANS) 
    if (gm->tsc[x] == 0.0) gm->tsc[x] = -1.0;

  /* Emissions */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 0; k <= gm->M; k++)
      {
	gm->rsc[x][k*2]   = (gm->rsc[x][k*2]   <= -eslINFINITY) ? -eslINFINITY : roundf(om->scale_w * gm->rsc[x][k*2]);
	gm->rsc[x][k*2+1] = 0.0;	/* insert score: VF makes it zero no matter what. */
      }	

  /* Specials */
  for (k = 0; k < p7P_NXSTATES; k++)
    for (x = 0; x < p7P_NXTRANS; x++)
      gm->xsc[k][x] = (gm->xsc[k][x] <= -eslINFINITY) ? -eslINFINITY : roundf(om->scale_w * gm->xsc[k][x]);

  /* 3nat approximation: NN, CC, JJ hardcoded 0 in limited precision */
  gm->xsc[p7P_N][p7P_LOOP] =  gm->xsc[p7P_J][p7P_LOOP] =  gm->xsc[p7P_C][p7P_LOOP] = 0.0;

  return eslOK;
}
/*------------ end, P7_OPROFILE debugging tools  ----------------*/



/*****************************************************************
 * 5. Benchmark driver.
 *****************************************************************/

#ifdef p7OPROFILE_BENCHMARK
/* Timing profile conversion.
   gcc -o benchmark-oprofile -std=gnu99 -g -Wall -msse2 -I.. -L.. -I../../easel -L../../easel -Dp7OPROFILE_BENCHMARK\
      p7_oprofile.c -lhmmer -leasel -lm 
   icc -o benchmark-oprofile -O3 -static -I.. -L.. -I../../easel -L../../easel -Dp7OPROFILE_BENCHMARK p7_oprofile.c -lhmmer -leasel -lm 
   ./benchmark-sse <hmmfile>         runs benchmark
 */
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_avx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-L",        eslARG_INT,    "400", NULL, NULL,  NULL,  NULL, NULL, "length of target sequence",                        0 },
  { "-N",        eslARG_INT, "100000", NULL, NULL,  NULL,  NULL, NULL, "number of conversions to time",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the generic implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  int             i;

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)           != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    p7_oprofile_Convert(gm, om);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M = %d\n", gm->M);

  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7OPROFILE_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/




  
/*****************************************************************
 * 6. Unit tests
 *****************************************************************/
#ifdef p7OPROFILE_TESTDRIVE


#endif /*p7OPROFILE_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/




/*****************************************************************
 * 7. Test driver
 *****************************************************************/
#ifdef p7OPROFILE_TESTDRIVE


#endif /*p7OPROFILE_TESTDRIVE*/
/*------------------- end, test driver --------------------------*/


/*****************************************************************
 * 8. Example
 *****************************************************************/
#ifdef p7OPROFILE_EXAMPLE
/* gcc -std=gnu99 -g -Wall -Dp7OPROFILE_EXAMPLE -I.. -I../../easel -L.. -L../../easel -o p7_oprofile_example p7_oprofile.c -lhmmer -leasel -lm
 * ./p7_oprofile_example <hmmfile>
 */
#include <p7_config.h>
#include <stdlib.h>
#include "easel.h"
#include "hmmer.h"

int
main(int argc, char **argv)
{
  char         *hmmfile = argv[1];
  ESL_ALPHABET *abc     = NULL;
  P7_HMMFILE   *hfp     = NULL;
  P7_HMM       *hmm     = NULL;
  P7_BG        *bg      = NULL;
  P7_PROFILE   *gm      = NULL;
  P7_OPROFILE  *om1     = NULL;
  P7_OPROFILE  *om2     = NULL;
  int           status;
  char          errbuf[eslERRBUFSIZE];

  status = p7_hmmfile_Open(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);  

  status = p7_hmmfile_Read(hfp, &abc, &hmm);
  if      (status == eslEFORMAT)   p7_Fail("Bad file format in HMM file %s:\n%s\n",          hfp->fname, hfp->errbuf);
  else if (status == eslEINCOMPAT) p7_Fail("HMM in %s is not in the expected %s alphabet\n", hfp->fname, esl_abc_DecodeType(abc->type));
  else if (status == eslEOF)       p7_Fail("Empty HMM file %s? No HMM data found.\n",        hfp->fname);
  else if (status != eslOK)        p7_Fail("Unexpected error in reading HMMs from %s\n",     hfp->fname);

  bg  = p7_bg_Create(abc);
  gm  = p7_profile_Create(hmm->M, abc);   
  om1 = p7_oprofile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, 400, p7_LOCAL);
  p7_oprofile_Convert(gm, om1);
  
  p7_oprofile_Dump(stdout, om1);

  om2 = p7_oprofile_Copy(om1);
  if (p7_oprofile_Compare(om1, om2, 0.001f, errbuf) != eslOK)    printf ("ERROR %s\n", errbuf);

  p7_oprofile_Destroy(om1);
  p7_oprofile_Destroy(om2);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  return eslOK;
}
#endif /*p7OPROFILE_EXAMPLE*/
/*----------------------- end, example --------------------------*/


