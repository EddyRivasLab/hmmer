/* The MSV filter implementation; SSE version.
 * 
 * A "filter" is a one-row, O(M), DP implementation that calculates an
 * approximated nat score (i.e. in limited precision - here, uchar)
 * and may have limited numeric range. It will return <eslERANGE> if
 * its numeric range is exceeded, in which case the caller will have
 * to obtain the score by another (probably slower) method.
 * 
 * Contents:
 *   1. p7_MSVFilter() implementation
 *   2. Benchmark driver
 *   3. Unit tests
 *   4. Test driver
 *   5. Example
 * 
 * SRE, Sun Nov 25 11:26:48 2007 [Casa de Gatos]
 */
#include <p7_config.h>

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#include <immintrin.h>    /* AVX2*/

#include "easel.h"
#include "esl_sse.h"
#include "esl_avx.h"
#include "esl_gumbel.h"

#include "hmmer.h"
#include "impl_avx.h"

/*****************************************************************
 * 1. The p7_MSVFilter() DP implementation.
 *****************************************************************/
 
/* Function:  p7_MSVFilter()
 * Synopsis:  Calculates MSV score, vewy vewy fast, in limited precision.
 * Incept:    SRE, Wed Dec 26 15:12:25 2007 [Janelia]
 *
 * Purpose:   Calculates an approximation of the MSV score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and a preallocated one-row DP matrix <ox>. Return the 
 *            estimated MSV score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow.
 *            
 *            The model may be in any mode, because only its match
 *            emission scores will be used. The MSV filter inherently
 *            assumes a multihit local mode, and uses its own special
 *            state transition scores, not the scores in the profile.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: MSV score (in nats)          
 *                      
 * Note:      We misuse the matrix <ox> here, using only a third of the
 *            first dp row, accessing it as <dp[0..Q-1]> rather than
 *            in triplets via <{MDI}MX(q)> macros, since we only need
 *            to store M state values. We know that if <ox> was big
 *            enough for normal DP calculations, it must be big enough
 *            to hold the MSVFilter calculation.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score overflows the limited range; in
 *            this case, this is a high-scoring hit.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_MSVFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m256i mpv;            /* previous row values                                       */
  register __m256i xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m256i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m256i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m256i biasv;	   /* emission bias in a vector                                 */
  uint8_t  xJ;                     /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQB_AVX(om->M);   /* segment length: # of vectors                              */
  __m256i *dp  = ox->dpb_avx[0];	   /* we're going to use dp[0][0..q..Q-1], not {MDI}MX(q) macros*/
  __m256i *rsc;			   /* will point at om->rbv[x] for residue x[i]                 */

  __m256i xJv;                     /* vector for states score                                   */
  __m256i tjbmv;                   /* vector for cost of moving from either J or N through B to an M state */
  __m256i tecv;                    /* vector for E->C  cost                                     */
  __m256i basev;                   /* offset for scores                                         */
  __m256i ceilingv;                /* saturateed simd value used to test for overflow           */
  __m256i tempv;                   /* work vector                                               */

  int cmp;
  int status = eslOK;

  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ16_avx)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M   = om->M;

  /* Try highly optimized ssv filter first */
  status = p7_SSVFilter(dsq, L, om, ret_sc);

  // uncomment when done testing MSV
  //if (status != eslENORESULT) return status;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.
   */
  biasv = _mm256_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
  for (q = 0; q < Q; q++) dp[q] = _mm256_setzero_si256();
  xJ   = 0;

  /* saturate simd register for overflow test */
  ceilingv = _mm256_cmpeq_epi8(biasv, biasv);
  basev = _mm256_set1_epi8((int8_t) om->base_b);

  tjbmv = _mm256_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv = _mm256_set1_epi8((int8_t) om->tec_b);

  xJv = _mm256_subs_epu8(biasv, biasv);
  xBv = _mm256_subs_epu8(basev, tjbmv);

#if eslDEBUGLEVEL > 0
  if (ox->debugging)
  {
      uint8_t xB;
      xB = _mm256_extract_epi16(xBv, 0);
      xJ = _mm256_extract_epi16(xJv, 0);
      p7_omx_DumpMFRow(ox, 0, 0, 0, xJ, xB, xJ);
  }
#endif

  for (i = 1; i <= L; i++)
  {
      rsc = om->rbv_avx[dsq[i]];
      xEv = _mm256_setzero_si256();      

      
      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12.
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros is our -infinity, which is why we reuse xEv as neginfmask for this call
       */
       __m256i dp_temp = dp[Q -1];
      mpv = esl_avx_rightshift_int8(dp_temp, xEv); 
      
      for (q = 0; q < Q; q++)
      {
        /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
        sv   = _mm256_max_epu8(mpv, xBv);
        sv   = _mm256_adds_epu8(sv, biasv);
        sv   = _mm256_subs_epu8(sv, *rsc);   rsc++;
        xEv  = _mm256_max_epu8(xEv, sv);

        mpv   = dp[q];   	  /* Load {MDI}(i-1,q) into mpv */
        dp[q] = sv;       	  /* Do delayed store of M(i,q) now that memory is usable */
      }

      /* test for the overflow condition */
      tempv = _mm256_adds_epu8(xEv, biasv);
      tempv = _mm256_cmpeq_epi8(tempv, ceilingv);
      cmp = _mm256_movemask_epi8(tempv);

      /* Now the "special" states, which start from Mk->E (->C, ->J->B)
       * Use shuffles instead of shifts so when the last max has completed,
       * the last four elements of the simd register will contain the
       * max value.  Then the last shuffle will broadcast the max value
       * to all simd elements.
       */

      xEv = _mm256_set1_epi8(esl_avx_hmax_epu8(xEv));  // broadcast the max byt\e from original xEv_AVX
      // to all bytes of xEv_AVX

      /* immediately detect overflow */
      if (cmp != 0x0000)
      {
        *ret_sc = eslINFINITY;
        return eslERANGE;
      }

      xEv = _mm256_subs_epu8(xEv, tecv);
      xJv = _mm256_max_epu8(xJv,xEv);
      
      xBv = _mm256_max_epu8(basev, xJv);
      xBv = _mm256_subs_epu8(xBv, tjbmv);
	  
#if eslDEBUGLEVEL > 0
      if (ox->debugging)
      {
        uint8_t xB, xE;
        xB = _mm256_extract_epi16(xBv, 0);
        xE = _mm256_extract_epi16(xEv, 0);
        xJ = _mm256_extract_epi16(xJv, 0);
        p7_omx_DumpMFRow(ox, i, xE, 0, xJ, xB, xJ);
      }
#endif
  } /* end loop over sequence residues 1..L */

  xJ = (uint8_t) _mm256_extract_epi16(xJv, 0);

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  *ret_sc = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */

  return eslOK;
}
/*------------------ end, p7_MSVFilter() ------------------------*/

