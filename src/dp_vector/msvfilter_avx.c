/* The MSV filter implementation; AVX version.
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
 *   6. Copyright and license information
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#ifdef HAVE_AVX2
  #include <immintrin.h>  /* AVX2 */
  #include "esl_avx.h"
#endif
#include "easel.h"
#include "esl_sse.h"
#include "esl_gumbel.h"

#include "base/p7_hmmwindow.h"
#include "search/p7_pipeline.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"
#include "dp_vector/ssvfilter.h"
#include "dp_vector/msvfilter.h"
#include "x86intrin.h"

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
p7_MSVFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
#ifdef HAVE_AVX2 
  uint8_t  xJ;                     /* special states' scores                  */
  register __m256i mpv_AVX;            /* previous row values                                       */
  register __m256i xEv_AVX;      /* E state: keeps max for Mk->E as we go                     */
  register __m256i xBv_AVX;      /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m256i sv_AVX;       /* temp storage of 1 curr row value in progress              */
  register __m256i biasv_AVX;    /* emission bias in a vector                                 */
  __m256i *dp_AVX;               /* the dp row memory                                         */
  __m256i *rsc_AVX;        /* will point at om->rbv[x] for residue x[i]                 */
  __m256i xJv_AVX;                     /* vector for states score                                   */
  __m256i tjbmv_AVX;                   /* vector for cost of moving {JN}->B->M                      */
  __m256i tecv_AVX;                    /* vector for E->C  cost                                     */
  __m256i basev_AVX;                   /* offset for scores                                         */
  __m256i ceilingv_AVX;                /* saturated simd value used to test for overflow            */
  __m256i tempv_AVX;                   /* work vector                                               */
  int Q_AVX        = P7_NVB_AVX(om->M);    /* segment length: # of vectors                              */
  int q_AVX;         /* counter over vectors 0..nq-1                              */

  int i;         /* counter over sequence positions 1..L                      */
  int     cmp;
  int     status;
//printf("Starting MSVFilter\n");
  /* Contract checks */
  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  /* note however that it makes no sense to run MSV w/ a model in glocal mode                            */

  /* Try highly optimized Knudsen SSV filter first. 
   * Note that SSV doesn't use any main memory (from <ox>) at all! 
  */
 
   status = p7_SSVFilter_avx(dsq, L, om, ret_sc);
   if (status != eslENORESULT) return status;

  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of MSV filter matrix failed");

  dp_AVX = ox->dp_AVX;  /* ditto this */

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_MSVFILTER;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.  */

  biasv_AVX = _mm256_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */

  for (q_AVX = 0; q_AVX < Q_AVX; q_AVX++) dp_AVX[q_AVX] = _mm256_setzero_si256(); 
  /* saturate simd register for overflow test */

  ceilingv_AVX = _mm256_cmpeq_epi8(biasv_AVX, biasv_AVX);
  basev_AVX    = _mm256_set1_epi8((int8_t) om->base_b);
  tjbmv_AVX    = _mm256_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv_AVX     = _mm256_set1_epi8((int8_t) om->tec_b);
  xJv_AVX      = _mm256_subs_epu8(biasv_AVX, biasv_AVX);
  xBv_AVX      = _mm256_subs_epu8(basev_AVX, tjbmv_AVX);

#ifdef p7_DEBUGGING
  if (ox->do_dumping)
    {
      uint8_t xB;
      xB = _mm_extract_epi16(xBv, 0);
      xJ = _mm_extract_epi16(xJv, 0);
      p7_filtermx_DumpMFRow(ox, 0, 0, 0, xJ, xB, xJ);
    }
#endif


  for (i = 1; i <= L; i++)  /* Outer loop over residues*/
    {

      rsc_AVX = om->rbv_AVX[dsq[i]];
      xEv_AVX = _mm256_setzero_si256();      

      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically, which is our -infinity.
       */
       __m256i dp_temp_AVX = dp_AVX[Q_AVX -1];
      mpv_AVX = esl_avx_leftshift_one(dp_temp_AVX);
      
      for (q_AVX = 0; q_AVX < Q_AVX; q_AVX++)
      {
        /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
        sv_AVX   = _mm256_max_epu8(mpv_AVX, xBv_AVX);
        sv_AVX   = _mm256_adds_epu8(sv_AVX, biasv_AVX);
        sv_AVX   = _mm256_subs_epu8(sv_AVX, *rsc_AVX);   rsc_AVX++;
        xEv_AVX  = _mm256_max_epu8(xEv_AVX, sv_AVX);

        mpv_AVX   = dp_AVX[q_AVX];      /* Load {MDI}(i-1,q) into mpv */
        dp_AVX[q_AVX] = sv_AVX;           /* Do delayed store of M(i,q) now that memory is usable */
      }

      /* test for the overflow condition */
      tempv_AVX = _mm256_adds_epu8(xEv_AVX, biasv_AVX);
      tempv_AVX = _mm256_cmpeq_epi8(tempv_AVX, ceilingv_AVX);
      cmp   = _mm256_movemask_epi8(tempv_AVX);

      /* Now the "special" states, which start from Mk->E (->C, ->J->B)
       * Use shuffles instead of shifts so when the last max has completed,
       * the last four elements of the simd register will contain the
       * max value.  Then the last shuffle will broadcast the max value
       * to all simd elements.
       */

      xEv_AVX = _mm256_set1_epi8(esl_avx_hmax_epu8(xEv_AVX));  // broadcast the max byte from original xEv_AVX
      // to all bytes of xEv_AVX

      /* immediately detect overflow */
      if (cmp != 0x0000) {
   //            MSV_end_time = __rdtsc();
   //     MSV_time += (MSV_end_time - MSV_start_time);
        *ret_sc = eslINFINITY; return eslERANGE; }

      xEv_AVX = _mm256_subs_epu8(xEv_AVX, tecv_AVX);
      xJv_AVX = _mm256_max_epu8(xJv_AVX,xEv_AVX);
      xBv_AVX = _mm256_max_epu8(basev_AVX, xJv_AVX);
      xBv_AVX = _mm256_subs_epu8(xBv_AVX, tjbmv_AVX);

#ifdef p7_DEBUGGING
      if (ox->do_dumping)
	{
	  uint8_t xB, xE;
	  xB = _mm_extract_epi16(xBv, 0);
	  xE = _mm_extract_epi16(xEv, 0);
	  xJ = _mm_extract_epi16(xJv, 0);
	  p7_filtermx_DumpMFRow(ox, i, xE, 0, xJ, xB, xJ);
	}
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  xJ =  _mm256_extract_epi8(xJv_AVX, 0);

  *ret_sc = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
    /*     MSV_end_time = __rdtsc();
        MSV_time += (MSV_end_time - MSV_start_time); */
  return eslOK;
  #endif
  #ifndef HAVE_AVX2
  return eslENORESULT;  // Stub so we have something to link if we build without AVX2 support
  #endif
  }

/*------------------ end, p7_MSVFilter() ------------------------*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/



