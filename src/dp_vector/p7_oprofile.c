/* Routines for the P7_OPROFILE structure:  
 * a search profile in an optimized implementation.
 * 
 * Contents:
 *   1. The P7_OPROFILE object: allocation, initialization, destruction.
 *   2. Conversion from generic P7_PROFILE to optimized P7_OPROFILE
 *   3. Conversion from optimized P7_OPROFILE to compact score arrays
 *   4. Debugging and development utilities.
 *   5. Benchmark driver.
 *   6. Example.
 *   7. Copyright and license information.
 *
 *  Change since HMMER3:  Many of the routines in this file are now front-ends that check the SIMD architecture
 *  that the program is running on and call a routine that uses the appropriate ISA.  See p7_oprofile_sse.c, p7_oprofile_avx.c, etc.
 *  Note that you will get (slightly) better performance if you call the appropriate SIMD routine directly when possible.  For example,
 *  if you call any of the SIMD oprofile routines from within the SSE version of a filter, you can save time by just calling the SSE version
 *  of the function.
 */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* roundf() */

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "base/p7_bg.h"
#include "base/p7_hmm.h"
#include "base/p7_profile.h"

#include "build/modelsample.h"
#include "search/modelconfig.h"

#include "hardware/hardware.h"
#include "dp_vector/p7_oprofile.h"

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
P7_OPROFILE *
p7_oprofile_Create(int allocM, const ESL_ALPHABET *abc, SIMD_TYPE simd)
{
  switch(simd){
    case SSE:
      return p7_oprofile_Create_sse(allocM, abc);
      break;
    case AVX:
      return p7_oprofile_Create_avx(allocM, abc);
      break;
    case AVX512:
      return p7_oprofile_Create_avx512(allocM, abc);
      break;
    case NEON:
      return p7_oprofile_Create_neon(allocM, abc);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_oprofile_Create");  
  }
}

/* Function:  p7_oprofile_IsLocal()
 * Synopsis:  Returns TRUE if profile is in local alignment mode.
 * Incept:    SRE, Sat Aug 16 08:46:00 2008 [Janelia]
 */
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
void
p7_oprofile_Destroy(P7_OPROFILE *om)
{
 switch(om->simd){
    case SSE:
      p7_oprofile_Destroy_sse(om);
      break;
    case AVX:
      p7_oprofile_Destroy_avx(om);
      break;
    case AVX512:
      p7_oprofile_Destroy_avx512(om);
      break;
    case NEON:
      p7_oprofile_Destroy_neon(om);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_oprofile_Destroy");  
  }
}

/* Function:  p7_oprofile_Sizeof()
 * Synopsis:  Return the allocated size of a <P7_OPROFILE>.
 * Incept:    SRE, Wed Mar  2 10:09:21 2011 [Janelia]
 *
 * Purpose:   Returns the allocated size of a <P7_OPROFILE>,
 *            in bytes.
 *            
 *            Very roughly, M*284 bytes, for a model of length M; 60KB
 *            for a typical model; 30MB for a design limit M=100K
 *            model.
 */
size_t
p7_oprofile_Sizeof(const P7_OPROFILE *om)
{
 switch(om->simd){
    case SSE:
      return p7_oprofile_Sizeof_sse(om);
      break;
    case AVX:
      return p7_oprofile_Sizeof_avx(om);
      break;
    case AVX512:
      return p7_oprofile_Sizeof_avx512(om);
      break;
    case NEON:
      return p7_oprofile_Sizeof_neon(om);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_oprofile_Sizeof");  
  }
}


/* Function:  p7_oprofile_Clone()
 * Synopsis:  Create a new copy of an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:03:19 2007 [Casa de Gatos]
 *
 * Purpose:   Create a newly allocated copy of <om1> and return a ptr
 *            to it.
 *            
 * Throws:    <NULL> on allocation error.
 */
P7_OPROFILE *
p7_oprofile_Clone(const P7_OPROFILE *om1)
{
  switch(om1->simd){
    case SSE:
      return p7_oprofile_Clone_sse(om1);
      break;
    case AVX:
      return p7_oprofile_Clone_avx(om1);
      break;
    case AVX512:
      return p7_oprofile_Clone_avx512(om1);
      break;
    case NEON:
      return p7_oprofile_Clone_neon(om1);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_oprofile_Clone");  
  }
}

/* Function:  p7_oprofile_Shadow()
 * Synopsis:  Create a shadow of an optimized profile, for use in multithreading
 *
 * Synopsis:  Allocate a cloned copy of the optimized profile structure.  All
 *            allocated memory from the original profile is not reallocated.
 *            The cloned copy will point to the same memory as the original.
 *
 * Purpose:   Allocate only the shell of a new <P7_OPROFILE>, and memcpy()
 *            the contents of <om1> into it.
 *            
 *            This gets used in multithreading. It's a hack. Almost
 *            all of the data in a profile is constant during a
 *            search, so threads can share pointers to it.  The
 *            exception is the length modeling ENJC transition costs,
 *            which have to be reconfigured for every new target seq;
 *            we need that data to be thread-specific and threadsafe.
 *            It happens that because the length model params are
 *            runtime arrays in a <P7_OPROFILE>, if we memcpy() the
 *            contents, we get reference copies of pointers to all the
 *            dynamic-allocated data (that we treat as constant), but
 *            actual copies of the static arrays (that we need to
 *            change).
 *            
 *            Caller still frees the shadow with
 *            <p7_oprofile_Destroy()>; that routine can tell the
 *            difference between a real (fully allocated) profile and
 *            a shadow, using the <om->is_shadow> flag.
 *
 * Returns:   Pointer to the new shadow. Caller frees with <p7_oprofile_Destroy()>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_OPROFILE *
p7_oprofile_Shadow(const P7_OPROFILE *om1)
{
  P7_OPROFILE  *om2  = NULL;
  int           status;

  ESL_ALLOC(om2, sizeof(P7_OPROFILE));
  memcpy(om2, om1, sizeof(P7_OPROFILE));
  om2->is_shadow  = TRUE;
  return om2;

 ERROR:
  p7_oprofile_Destroy(om2);
  return NULL;
}




/*----------------- end, P7_OPROFILE structure ------------------*/



/*****************************************************************
 * 2. Conversion from generic P7_PROFILE to optimized P7_OPROFILE
 *****************************************************************/

/* biased_byteify()
 * Converts original log-odds residue score to a rounded biased uchar cost.
 * Match emission scores for MSVFilter get this treatment.
 * e.g. a score of +3.2, with scale 3.0 and bias 12, becomes 2.
 *    3.2*3 = 9.6; rounded = 10; bias-10 = 2.
 * When used, we add the bias, then subtract this cost.
 * (A cost of +255 is our -infinity "prohibited event")
 */
uint8_t
biased_byteify(P7_OPROFILE *om, float sc)
{
  uint8_t b;
  sc  = -1.0f * roundf(om->scale_b * sc);
  int32_t q = round(sc);
  uint8_t b1 = (uint8_t) q;                              /* ugh. sc is now an integer cost represented in a float...           */
  b   = (sc > 255 - om->bias_b) ? 255 : b1 + om->bias_b; /* and now we cast, saturate, and bias it to an unsigned char cost... */
  return b;
}
 
/* unbiased_byteify()
 * Convert original transition score to a rounded uchar cost
 * Transition scores for MSVFilter get this treatment.
 * e.g. a score of -2.1, with scale 3.0, becomes a cost of 6.
 * (A cost of +255 is our -infinity "prohibited event")
 */
uint8_t 
unbiased_byteify(P7_OPROFILE *om, float sc)
{
  uint8_t b;
  sc  = -1.0f * roundf(om->scale_b * sc);       /* ugh. sc is now an integer cost represented in a float...    */
  uint32_t q = round(sc);
  b   = (sc > 255.) ? 255 : (uint8_t) q;        /* and now we cast and saturate it to an unsigned char cost... */
  return b;
}
 
/* wordify()
 * Converts log probability score to a rounded signed 16-bit integer cost.
 * Both emissions and transitions for ViterbiFilter get this treatment.
 * No bias term needed, because we use signed words. 
 *   e.g. a score of +3.2, with scale 500.0, becomes +1600.
 */
int16_t 
wordify(P7_OPROFILE *om, float sc)
{
  sc  = roundf(om->scale_w * sc);
  if      (sc >=  32767.0) return  32767;
  else if (sc <= -32768.0) return -32768;
  else return (int16_t) sc;
}


/* sf_conversion():
 * Author: Bjarne Knudsen
 * 
 * Generates the SSVFilter() parts of the profile <om> scores
 * from the completed MSV score.  This includes calculating 
 * special versions of the match scores for using the the
 * ssv filter.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 * 
 * This function should probably never be called.  It should always be possible
 * to call the correct SIMD version directly, but it's included for completeness
 *
 */
static int
sf_conversion(P7_OPROFILE *om)
{
 switch(om->simd){
    case SSE:
      sf_conversion_sse(om);
      break;
    case AVX:
      sf_conversion_avx(om);
      break;
    case AVX512:
      sf_conversion_avx512(om);
      break;
    case NEON:
      sf_conversion_neon(om);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to sf_conversion");  
  }
}

/* mf_conversion(): 
 * 
 * This builds the MSVFilter() parts of the profile <om>, scores
 * in lspace uchars (16-way parallel), by rescaling, rounding, and
 * casting the scores in <gm>.
 * 
 * Returns <eslOK> on success;
 * throws <eslEINVAL> if <om> hasn't been allocated properly.
 * 
 * This function should probably never be called.  It should always be possible
 * to call the correct SIMD version directly, but it's included for completeness
 *
 */
static int
mf_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  switch(om->simd){
    case SSE:
      mf_conversion_sse(gm, om);
      break;
    case AVX:
      mf_conversion_avx(gm, om);
      break;
    case AVX512:
      mf_conversion_avx512(gm, om);
      break;
    case NEON:
      mf_conversion_neon(gm, om);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to mf_conversion");  
  }
}


/* vf_conversion(): 
 * 
 * This builds the ViterbiFilter() parts of the profile <om>, scores
 * in lspace swords (8-way parallel), by rescaling, rounding, and
 * casting the scores in <gm>.
 * 
 * Returns <eslOK> on success;
 * throws <eslEINVAL> if <om> hasn't been allocated properly.
 */
static int
vf_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
 switch(om->simd){
    case SSE:
      vf_conversion_sse(gm, om);
      break;
    case AVX:
      vf_conversion_avx(gm, om);
      break;
    case AVX512:
      vf_conversion_avx512(gm, om);
      break;
    case NEON:
      vf_conversion_neon(gm, om);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to vf_conversion");  
  }
}


/* fb_conversion()
 * This builds the Forward/Backward part of the optimized profile <om>,
 * where we use odds ratios (not log-odds scores).
 */
static int
fb_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  switch(om->simd){
    case SSE:
      fb_conversion_sse(gm, om);
      break;
    case AVX:
      fb_conversion_avx(gm, om);
      break;
    case AVX512:
      fb_conversion_avx512(gm, om);
      break;
    case NEON:
      fb_conversion_neon(gm, om);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to fb_conversion");  
  }
}


/* Function:  p7_oprofile_Convert()
 * Synopsis:  Converts standard profile to an optimized one.
 * Incept:    SRE, Mon Nov 26 07:38:57 2007 [Janelia]
 *
 * Purpose:   Convert a standard profile <gm> to an optimized profile <om>,
 *            where <om> has already been allocated for a profile of at 
 *            least <gm->M> nodes and the same emission alphabet <gm->abc>.
 *            
 *            Retain the length model and uni/multihit config of <gm>.
 *            Set <om> to the appropriate local mode, ignoring whether
 *            <gm> was glocal, dual-mode, or local. Optimized
 *            profiles are local only, not dual-mode local/glocal.
 *            
 *            Usually, <gm> would be expected to be in dual-mode
 *            multihit configuration, in our production code.  The
 *            <om> comes out in local multihit config, with the same
 *            length model.
 *            
 *            <om> cannot be a "shadow" (created by
 *            <p7_oprofile_Shadow()>); it must be a real allocation
 *            for the <P7_OPROFILE>.
 *
 * Args:      gm - profile to optimize
 *            om - allocated optimized profile for holding the result.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_oprofile_Convert(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int status, z;

  ESL_DASSERT1(( ! om->is_shadow ));
  ESL_DASSERT1(( gm->abc->type == om->abc->type));
  ESL_DASSERT1(( gm->M         <= om->allocM));

  if      (gm->nj == 0.0) om->mode = p7_UNILOCAL;
  else if (gm->nj == 1.0) om->mode = p7_LOCAL;
  else    ESL_EXCEPTION(eslEINVAL, "oprofile must be unilocal or local");

  om->L          = gm->L;
  om->M          = gm->M;
  om->nj         = gm->nj;
  om->max_length = gm->max_length;

  switch(om->simd){
      case SSE:
        if ((status =  mf_conversion_sse(gm, om)) != eslOK) return status;   /* MSVFilter()'s information     */
        if ((status =  vf_conversion_sse(gm, om)) != eslOK) return status;   /* ViterbiFilter()'s information */
        if ((status =  fb_conversion_sse(gm, om)) != eslOK) return status;   /* ForwardFilter()'s information */
        break;
      case AVX:
        if ((status =  mf_conversion_avx(gm, om)) != eslOK) return status;   /* MSVFilter()'s information     */
        if ((status =  vf_conversion_avx(gm, om)) != eslOK) return status;   /* ViterbiFilter()'s information */
        if ((status =  fb_conversion_avx(gm, om)) != eslOK) return status;   /* ForwardFilter()'s information */
        break;
      case AVX512:
        if ((status =  mf_conversion_avx512(gm, om)) != eslOK) return status;   /* MSVFilter()'s information     */
        if ((status =  vf_conversion_avx512(gm, om)) != eslOK) return status;   /* ViterbiFilter()'s information */
        if ((status =  fb_conversion_avx512(gm, om)) != eslOK) return status;   /* ForwardFilter()'s information */
        break;
    case NEON:
        if ((status =  mf_conversion_neon(gm, om)) != eslOK) return status;   /* MSVFilter()'s information     */
        if ((status =  vf_conversion_neon(gm, om)) != eslOK) return status;   /* ViterbiFilter()'s information */
        if ((status =  fb_conversion_neon(gm, om)) != eslOK) return status;   /* ForwardFilter()'s information */
        break;
      default:
        p7_Fail("Unrecognized SIMD type passed to p7_oprofile_Convert");  
  }

  if (om->name != NULL) free(om->name);
  if (om->acc  != NULL) free(om->acc);
  if (om->desc != NULL) free(om->desc);
  if ((status = esl_strdup(gm->name, -1, &(om->name))) != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->acc,  -1, &(om->acc)))  != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->desc, -1, &(om->desc))) != eslOK) goto ERROR;
  strcpy(om->rf,        gm->rf);
  strcpy(om->mm,        gm->mm);
  strcpy(om->cs,        gm->cs);
  strcpy(om->consensus, gm->consensus);
  for (z = 0; z < p7_NEVPARAM; z++) om->evparam[z] = gm->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) om->cutoff[z]  = gm->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) om->compo[z]   = gm->compo[z];

  return eslOK;

 ERROR:
  return status;
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
  om->L     = L;
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
 * 3. Conversion from optimized P7_OPROFILE to compact score arrays
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
int
p7_oprofile_GetFwdTransitionArray(const P7_OPROFILE *om, int type, float *arr )
{
  switch(om->simd){
    case SSE:
      return p7_oprofile_GetFwdTransitionArray_sse(om, type, arr);
      break;
    case AVX:
      return p7_oprofile_GetFwdTransitionArray_avx(om, type, arr);
      break;
    case AVX512:
      return p7_oprofile_GetFwdTransitionArray_avx512(om, type, arr);
      break;
    case NEON:
      return p7_oprofile_GetFwdTransitionArray_neon(om, type, arr);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_oprofile_GetFwdTransitionArray");  
  }
  return eslOK;

}

/* Function:  p7_oprofile_GetMSVEmissionScoreArray()
 * Synopsis:  Retrieve MSV residue emission scores from an optimized
 *            profile into an array
 *
 * Purpose:   Extract an implicitly 2D array of 8-bit int MSV residue
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
int
p7_oprofile_GetMSVEmissionScoreArray(const P7_OPROFILE *om, uint8_t *arr )
{
  switch(om->simd){
    case SSE:
      return p7_oprofile_GetMSVEmissionScoreArray_sse(om, arr);
      break;
    case AVX:
      return p7_oprofile_GetMSVEmissionScoreArray_avx(om, arr);
      break;
    case AVX512:
      return p7_oprofile_GetMSVEmissionScoreArray_avx512(om, arr);
      break;
    case NEON:
      return p7_oprofile_GetMSVEmissionScoreArray_neon(om, arr);
      break; 
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_oprofile_GetMSVEmissionScoreArray");  
  }
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
 */
int
p7_oprofile_GetFwdEmissionScoreArray(const P7_OPROFILE *om, float *arr )
{
  switch(om->simd){
    case SSE:
      return p7_oprofile_GetFwdEmissionScoreArray_sse(om, arr);
      break;
    case AVX:
      return p7_oprofile_GetFwdEmissionScoreArray_avx(om, arr);
      break;
    case AVX512:
      return p7_oprofile_GetFwdEmissionScoreArray_avx512(om, arr);
      break;
    case NEON:
      return p7_oprofile_GetFwdEmissionScoreArray_neon(om, arr);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_oprofile_GetFwdEmissionScoreArray");  
  }
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
 */
int
p7_oprofile_GetFwdEmissionArray(const P7_OPROFILE *om, P7_BG *bg, float *arr )
{
  switch(om->simd){
    case SSE:
      return p7_oprofile_GetFwdEmissionArray_sse(om, bg, arr);
      break;
    case AVX:
      return p7_oprofile_GetFwdEmissionArray_avx(om, bg, arr);
      break;
    case AVX512:
      return p7_oprofile_GetFwdEmissionArray_avx512(om, bg, arr);
      break;
    case NEON:
      return p7_oprofile_GetFwdEmissionArray_neon(om, bg, arr);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_oprofile_GetFwdEmissionArray");  
  }
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
 switch(om->simd){
    case SSE:
      return oprofile_dump_mf_sse(fp, om);
      break;
    case AVX:
      return oprofile_dump_mf_avx(fp, om);
      break;
    case AVX512:
      return oprofile_dump_mf_sse(fp, om);
      break;
    case NEON:
      return oprofile_dump_mf_neon(fp, om);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to oprofile_dump_mf");
  }
}

/* oprofile_dump_vf()
 * 
 * Dump the ViterbiFilter part of a profile <om> to <stdout>.
 */
static int
oprofile_dump_vf(FILE *fp, const P7_OPROFILE *om)
{
  switch(om->simd){
    case SSE:
      return oprofile_dump_vf_sse(fp, om);
      break;
    case AVX:
      return oprofile_dump_vf_avx(fp, om);
      break;
    case AVX512:
      return oprofile_dump_vf_sse(fp, om);
      break;
    case NEON:
      return oprofile_dump_vf_neon(fp, om);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to oprofile_dump_vf");
  }
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
  switch(om->simd){
    case SSE:
      return oprofile_dump_fb_sse(fp, om, width, precision);
      break;
    case AVX:
      return oprofile_dump_fb_avx(fp, om, width, precision);
      break;
    case AVX512:
      return oprofile_dump_fb_sse(fp, om, width, precision);
      break;
    case NEON:
      return oprofile_dump_fb_neon(fp, om, width, precision);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to oprofile_dump_fb");
  }
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
 *
 * Purpose:   Sample a random profile of <M> nodes for alphabet <abc>,
 *            using <r> as the source of random numbers. Parameterize
 *            it for generation of target sequences of mean length
 *            <L>. Calculate its log-odds scores using background
 *            model <bg>.
 *            
 *            Caller may optionally obtain the corresponding hmm by
 *            passing a non-<NULL> <opt_hmm>, and/or the corresponding
 *            profile by passing a non-<NULL> <opt_gm>. If the <gm> is
 *            obtained, it is configured for local-only mode and for a
 *            target length of <L>, so that its scores will match the
 *            <om> (as closely as roundoff allows).
 *            
 * Args:      r       - random number generator
 *            abc     - emission alphabet 
 *            bg      - background frequency model
 *            M       - size of sampled profile, in nodes
 *            L       - configured target seq mean length
 *            opt_hmm - optRETURN: sampled HMM
 *            opt_gm  - optRETURN: sampled normal profile, (local,L) mode
 *            opt_om  - RETURN: optimized profile, length config'ed to L
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
  P7_HARDWARE *hw;
  if ((hw = p7_hardware_Create ()) == NULL)  { status = eslEMEM; goto ERROR; }

  if ((gm = p7_profile_Create (M, abc)) == NULL)  { status = eslEMEM; goto ERROR; }
  if ((om = p7_oprofile_Create(M, abc, hw->simd)) == NULL)  { status = eslEMEM; goto ERROR; }

  if ((status = p7_modelsample(r, M, abc, &hmm))        != eslOK) goto ERROR;
  if ((status = p7_profile_ConfigLocal(gm, hmm, bg, L)) != eslOK) goto ERROR;
  if ((status = p7_oprofile_Convert(gm, om))            != eslOK) goto ERROR;
  if ((status = p7_oprofile_ReconfigLength(om, L))      != eslOK) goto ERROR;

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
 *            of <tol> using <esl_FCompare()>.
 *            
 *            If a comparison fails, an informative error message is
 *            left in <errmsg> to indicate why.
 *            
 *            Internal allocation sizes are not compared, only the
 *            data.
 *            
 * Args:      om1    - one optimized profile to compare
 *            om2    - the other
 *            tol    - floating point comparison tolerance; see <esl_FCompare()>
 *            errmsg - ptr to array of at least <eslERRBUFSIZE> characters.
 *            
 * Returns:   <eslOK> on effective equality;  <eslFAIL> on difference.
 */
int
p7_oprofile_Compare(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg)
{
 
  switch(om1->simd){
    case SSE:
      return p7_oprofile_Compare_sse(om1, om2, tol, errmsg);
      break;
    case AVX:
      return p7_oprofile_Compare_avx(om1, om2, tol, errmsg);
      break;
    case AVX512:
      return p7_oprofile_Compare_avx512(om1, om2, tol, errmsg);
      break;
    case NEON:
      return p7_oprofile_Compare_neon(om1, om2, tol, errmsg);
      break;
    default:
      p7_Fail("Unrecognized SIMD type passed to p7_oprofile_Compare");
  }
}


/* Function:  p7_profile_SameAsMF()
 * Synopsis:  Set a generic profile's scores to give MSV scores.
 *
 * Purpose:   Set a generic profile's scores so that the reference Viterbi
 *            implementation will give the same score as <p7_MSVFilter()>.
 *            All t_MM scores = 0; all other core transitions = -inf;
 *            multihit local mode; all <t_BMk> entries uniformly <log 2/(M(M+1))>;
 *            <tCC, tNN, tJJ> scores 0; total approximated later as -3;
 *            rounded in the same way as the 8-bit limited precision.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_profile_SameAsMF(const P7_OPROFILE *om, P7_PROFILE *gm)
{
  int    k,x;
  float  tbm = roundf(om->scale_b * (log(2.0f / ((float) gm->M * (float) (gm->M+1)))));

  /* Transitions */
  esl_vec_FSet(gm->tsc, p7P_NTRANS * gm->M, -eslINFINITY);
  for (k = 1; k <  gm->M; k++) P7P_TSC(gm, k, p7P_MM)  = 0.0f;
  for (k = 0; k <  gm->M; k++) P7P_TSC(gm, k, p7P_LM) = tbm;
  
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
   ./benchmark-sse <hmmfile>      
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

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

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, L);
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
 * 6. Example
 *****************************************************************/
#ifdef p7OPROFILE_EXAMPLE
/* 
 * ./p7_oprofile_example <hmmfile>
 */
#include "p7_config.h"

#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
   /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                            0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "example main() for p7_oprofile.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  char         *hmmfile = esl_opt_GetArg(go, 1);
  ESL_ALPHABET *abc     = NULL;
  P7_HMMFILE   *hfp     = NULL;
  P7_HMM       *hmm     = NULL;
  P7_BG        *bg      = NULL;
  P7_PROFILE   *gm      = NULL;
  P7_OPROFILE  *om1     = NULL;
  P7_OPROFILE  *om2     = NULL;
  int           status;
  char          errbuf[eslERRBUFSIZE];

  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
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

  p7_profile_ConfigLocal(gm, hmm, bg, 400);
  p7_oprofile_Convert(gm, om1);
  
  p7_oprofile_Dump(stdout, om1);

  om2 = p7_oprofile_Clone(om1);
  if (p7_oprofile_Compare(om1, om2, 0.001f, errbuf) != eslOK)    printf ("ERROR %s\n", errbuf);

  p7_oprofile_Destroy(om1);
  p7_oprofile_Destroy(om2);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7OPROFILE_EXAMPLE*/
/*----------------------- end, example --------------------------*/



/*****************************************************************
 * @LICENSE@
 *   
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
