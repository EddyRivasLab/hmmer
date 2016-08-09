/* The MSV filter implementation; ARM NEON64 version.
 * 
 * Ported from Intel/SSE: Tyler Camp (University of Texas, Austin)
 * Refer to Intel/SSE version for general notes.
 * 
 * Contents:
 *   1. p7_MSVFilter() implementation
 *   2. Benchmark driver
 *   3. Unit tests
 *   4. Test driver
 *   5. Example
 */


#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <arm_neon.h>

#include "easel.h"
#include "esl_gumbel.h"

#include "esl_neon64.h"

#include "base/p7_hmmwindow.h"
#include "search/p7_pipeline.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"
#include "dp_vector/ssvfilter.h"
#include "dp_vector/msvfilter.h"


/*****************************************************************
 * 1. The p7_MSVFilter() DP implementation.
 *****************************************************************/
 
/* Function:  p7_MSVFilter()
 * Synopsis:  Calculates MSV score, vewy vewy fast, in limited precision.
 *
 * Purpose:   Calculates an approximation of the MSV score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and the one-row DP matrix <ox>. Return the 
 *            estimated MSV score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow.
 *            
 *            <ox> will be resized if needed. It's fine if it was
 *            just <_Reuse()'d> from a previous, smaller profile.
 *            
 *            The model may be in any mode, because only its match
 *            emission scores will be used. The MSV filter inherently
 *            assumes a multihit local mode, and uses its own special
 *            state transition scores, not the scores in the profile.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - filter DP matrix (one row)
 *            ret_sc  - RETURN: MSV score (in nats)          
 *                      
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score overflows the limited range; in
 *            this case, this is a high-scoring hit.
 *            <ox> may have been resized.
 *
 * Throws:    <eslEMEML> if <ox> reallocation fails.
 */
int
p7_MSVFilter_neon64(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
  #ifdef HAVE_NEON64
  register esl_neon_128i_t mpv;            /* previous row values                                       */
  register esl_neon_128i_t xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register esl_neon_128i_t xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register esl_neon_128i_t sv;		   /* temp storage of 1 curr row value in progress              */
  register esl_neon_128i_t biasv;	   /* emission bias in a vector                                 */
  uint8_t  xJ;                             /* special states' scores                                    */
  int i;			           /* counter over sequence positions 1..L                      */
  int q;			           /* counter over vectors 0..nq-1                              */
  int Q        = P7_NVB(om->M);            /* segment length: # of vectors                              */
  esl_neon_128i_t *dp;   	           /* the dp row memory                                         */
  esl_neon_128i_t *rsc;			   /* will point at om->rbv[x] for residue x[i]                 */
  esl_neon_128i_t xJv;                     /* vector for states score                                   */
  esl_neon_128i_t tjbmv;                   /* vector for cost of moving {JN}->B->M                      */
  esl_neon_128i_t tecv;                    /* vector for E->C  cost                                     */
  esl_neon_128i_t basev;                   /* offset for scores                                         */
  esl_neon_128i_t ceilingv;                /* saturated simd value used to test for overflow            */
  esl_neon_128i_t tempv;                   /* work vector                                               */
  int             status;

  /* Contract checks */
  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  /* note however that it makes no sense to run MSV w/ a model in glocal mode                            */

  /* Try highly optimized Knudsen SSV filter first. 
   * Note that SSV doesn't use any main memory (from <ox>) at all! 
   */
  if (( status = p7_SSVFilter_neon64(dsq, L, om, ret_sc)) != eslENORESULT) return status;

  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of MSV filter matrix failed");
  dp = ox->dp;			/* This MUST be set AFTER the GrowTo() call. */

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_MSVFILTER;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.  */
  biasv.s8x16 = vdupq_n_s8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
  for (q = 0; q < Q; q++) dp[q].s32x4 = vdupq_n_s32(0);

  /* saturate simd register for overflow test */
  ceilingv.u8x16 = vceqq_s8(biasv.s8x16, biasv.s8x16);
  basev.s8x16    = vdupq_n_s8((int8_t) om->base_b);
  tjbmv.s8x16    = vdupq_n_s8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv.s8x16     = vdupq_n_s8((int8_t) om->tec_b);
  xJv.u8x16      = vqsubq_u8(biasv.u8x16, biasv.u8x16);
  xBv.u8x16      = vqsubq_u8(basev.u8x16, tjbmv.u8x16);

#ifdef p7_DEBUGGING
  if (ox->do_dumping)
    {
      uint8_t xB;
      xB = vgetq_lane_s16(xBv.s16x8, 0);
      xJ = vgetq_lane_s16(xJv.s16x8, 0);
      p7_filtermx_DumpMFRow(ox, 0, 0, 0, xJ, xB, xJ);
    }
#endif


  for (i = 1; i <= L; i++)
    {
      rsc = om->rbv[dsq[i]];
      xEv.s32x4 = vdupq_n_s32(0);      

      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically, which is our -infinity.
       */
      mpv.s8x16 = vextq_s8(xEv.s8x16, dp[Q-1].s8x16, 15);
      for (q = 0; q < Q; q++)
	{
        /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
	  sv.u8x16   = vmaxq_u8(mpv.u8x16, xBv.u8x16);
	  sv.u8x16   = vqaddq_u8(sv.u8x16, biasv.u8x16);
	  sv.u8x16   = vqsubq_u8(sv.u8x16, (*rsc).u8x16);   rsc++;
	  xEv.u8x16  = vmaxq_u8(xEv.u8x16, sv.u8x16);

	  mpv   = dp[q];   	  /* Load {MDI}(i-1,q) into mpv */
	  dp[q] = sv;       	  /* Do delayed store of M(i,q) now that memory is usable */
	}

      /* test for the overflow condition */
      tempv.u8x16 = vqaddq_u8(xEv.u8x16, biasv.u8x16);
      tempv.u8x16 = vceqq_s8(tempv.s8x16, ceilingv.s8x16);
      union { esl_neon_128i_t v; uint64_t i[2]; uint8_t k[16]; }s;
      s.v = tempv;
      uint64_t mask = s.i[0] | s.i[1];
      /* Now the "special" states, which start from Mk->E (->C, ->J->B)
       * Use shuffles instead of shifts so when the last max has completed,
       * the last four elements of the simd register will contain the
       * max value.  Then the last shuffle will broadcast the max value
       * to all simd elements.
       */
	  
      xEv.u8x16 = vdupq_n_u8(esl_neon64_hmax_u8(xEv));

      /* immediately detect overflow */
      if (mask != 0) { *ret_sc = eslINFINITY; return eslERANGE; }

      xEv.u8x16 = vqsubq_u8(xEv.u8x16, tecv.u8x16);
      xJv.u8x16 = vmaxq_u8(xJv.u8x16,xEv.u8x16);
      xBv.u8x16 = vmaxq_u8(basev.u8x16, xJv.u8x16);
      xBv.u8x16 = vqsubq_u8(xBv.u8x16, tjbmv.u8x16);
	  
#ifdef p7_DEBUGGING
      if (ox->do_dumping)
	{
	  uint8_t xB, xE;
	  xB = vgetq_lane_s16(xBv.s16x8, 0);
	  xE = vgetq_lane_s16(xEv.s16x8, 0);
	  xJ = vgetq_lane_s16(xJv.s16x8, 0);
	  p7_filtermx_DumpMFRow(ox, i, xE, 0, xJ, xB, xJ);
	}
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  xJ = (uint8_t) vgetq_lane_s16(xJv.s16x8, 0);
  *ret_sc = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
  return eslOK;
  #endif /* HAVE_NEON64 */
  
  #ifndef HAVE_NEON64
    return 0;
  #endif
}
/*------------------ end, p7_MSVFilter() ------------------------*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/

