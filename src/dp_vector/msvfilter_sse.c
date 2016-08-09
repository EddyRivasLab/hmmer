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
 *   6. Copyright and license information
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>
#if p7_CPU_ARCH == x86
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
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


//uint64_t SSV_time, MSV_time;
uint64_t full_MSV_calls;
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
p7_MSVFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
  #ifdef HAVE_SSE2 // only build this function if we can compile SSE, otherwise make it a null function so that we don't have
  // missing function problems with the library 

  // Don't reuse names to allow use of multiple ISAs at same time for debugging
  register __m128i mpv;            /* previous row values                                       */
  register __m128i xEv;      /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;      /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i sv;       /* temp storage of 1 curr row value in progress              */
  register __m128i biasv;    /* emission bias in a vector                                 */
  __m128i *dp;               /* the dp row memory                                         */
  __m128i *rsc;        /* will point at om->rbv[x] for residue x[i]                 */
  __m128i xJv;                     /* vector for states score                                   */
  __m128i tjbmv;                   /* vector for cost of moving {JN}->B->M                      */
  __m128i tecv;                    /* vector for E->C  cost                                     */
  __m128i basev;                   /* offset for scores                                         */
  __m128i ceilingv;                /* saturated simd value used to test for overflow            */
  __m128i tempv;                   /* work vector                                               */
  int Q        = P7_NVB(om->M);    /* segment length: # of vectors                              */
  int q;         /* counter over vectors 0..nq-1                              */

 uint8_t  xJ;                     /* special states' scores                  */

 // extern uint64_t full_MSV_calls, SSV_time, MSV_time;

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
 //  printf("Calling SSVFilter\n");
//  if (( status = p7_SSVFilter(dsq, L, om, ret_sc)) != eslENORESULT) return status;
//   uint64_t SSV_start_time = __rdtsc();
   status = p7_SSVFilter_sse(dsq, L, om, ret_sc);
//   uint64_t SSV_end_time = __rdtsc();
 //  SSV_time += SSV_end_time - SSV_start_time;
   if (status != eslENORESULT) return status;
//full_MSV_calls++;
 // uint64_t MSV_start_time = __rdtsc();
 // uint64_t MSV_end_time;
  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of MSV filter matrix failed");

  dp = ox->dp;			/* This MUST be set AFTER the GrowTo() call. */

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_MSVFILTER;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.  */
  biasv = _mm_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
  for (q = 0; q < Q; q++) dp[q] = _mm_setzero_si128();
  /* saturate simd register for overflow test */
  ceilingv = _mm_cmpeq_epi8(biasv, biasv);
  basev    = _mm_set1_epi8((int8_t) om->base_b);
  tjbmv    = _mm_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv     = _mm_set1_epi8((int8_t) om->tec_b);
  xJv      = _mm_subs_epu8(biasv, biasv);
  xBv      = _mm_subs_epu8(basev, tjbmv);
  
//printf("Biases set\n");
  

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

      rsc = om->rbv[dsq[i]];
      xEv = _mm_setzero_si128();      

      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically, which is our -infinity.
       */
      mpv = _mm_slli_si128(dp[Q-1], 1);   
      for (q = 0; q < Q; q++)
      {
        /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
        sv   = _mm_max_epu8(mpv, xBv);
        sv   = _mm_adds_epu8(sv, biasv);
        sv   = _mm_subs_epu8(sv, *rsc);   rsc++;
        xEv  = _mm_max_epu8(xEv, sv);

        mpv   = dp[q];   	  /* Load {MDI}(i-1,q) into mpv */
        dp[q] = sv;       	  /* Do delayed store of M(i,q) now that memory is usable */
      }

      /* test for the overflow condition */
      tempv = _mm_adds_epu8(xEv, biasv);
      tempv = _mm_cmpeq_epi8(tempv, ceilingv);
      cmp   = _mm_movemask_epi8(tempv);

      /* Now the "special" states, which start from Mk->E (->C, ->J->B)
       * Use shuffles instead of shifts so when the last max has completed,
       * the last four elements of the simd register will contain the
       * max value.  Then the last shuffle will broadcast the max value
       * to all simd elements.
       */
      tempv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(2, 3, 0, 1));
      xEv   = _mm_max_epu8(xEv, tempv);
      tempv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(0, 1, 2, 3));
      xEv   = _mm_max_epu8(xEv, tempv);
      tempv = _mm_shufflelo_epi16(xEv, _MM_SHUFFLE(2, 3, 0, 1));
      xEv   = _mm_max_epu8(xEv, tempv);
      tempv = _mm_srli_si128(xEv, 1);
      xEv   = _mm_max_epu8(xEv, tempv);
      xEv   = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(0, 0, 0, 0));

      /* immediately detect overflow */
      if (cmp != 0x0000) {
    //    MSV_end_time = __rdtsc();
    //    MSV_time += (MSV_end_time - MSV_start_time);
       *ret_sc = eslINFINITY; return eslERANGE; }

      xEv = _mm_subs_epu8(xEv, tecv);
      xJv = _mm_max_epu8(xJv,xEv);
      xBv = _mm_max_epu8(basev, xJv);
      xBv = _mm_subs_epu8(xBv, tjbmv);

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

  xJ = (uint8_t) _mm_extract_epi16(xJv, 0);

  *ret_sc = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
    /*     MSV_end_time = __rdtsc();
        MSV_time += (MSV_end_time - MSV_start_time); */
  return eslOK;
  #endif //HAVE_SSE2
#ifndef HAVE_SSE2
  return 0;  // put this here so that compiler doesn't complain about reaching the end of function without a return value
#endif
  }

/*------------------ end, p7_MSVFilter_SSE() ------------------------*/



/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/



