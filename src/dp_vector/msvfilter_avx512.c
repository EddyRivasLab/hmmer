/* The MSV filter implementation; AVX-512 version.
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
 */
#include "p7_config.h"
#ifdef eslENABLE_AVX512

#include <stdio.h>
#include <math.h>

#include <x86intrin.h>  

#include "easel.h"
#include "esl_gumbel.h"

#include "base/p7_hmmwindow.h"
#include "search/p7_pipeline.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"
#include "dp_vector/ssvfilter.h"
#include "dp_vector/msvfilter.h"



/*****************************************************************
 * 1. The p7_MSVFilter() DP implementation.
 *****************************************************************/
 
/* Function:  p7_MSVFilter_avx512()
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
p7_MSVFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
  uint8_t  xJ;                     /* special states' scores                  */
  register __m512i mpv_AVX_512;    /* previous row values                                       */
  register __m512i xEv_AVX_512;    /* E state: keeps max for Mk->E as we go                     */
  register __m512i xBv_AVX_512;    /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m512i sv_AVX_512;     /* temp storage of 1 curr row value in progress              */
  register __m512i biasv_AVX_512;  /* emission bias in a vector                                 */
  __m512i *dp_AVX_512;             /* the dp row memory                                         */
  __m512i *rsc_AVX_512;            /* will point at om->rbv[x] for residue x[i]                 */
  __m512i xJv_AVX_512;             /* vector for states score                                   */
  __m512i tjbmv_AVX_512;           /* vector for cost of moving {JN}->B->M                      */
  __m512i tecv_AVX_512;            /* vector for E->C  cost                                     */
  __m512i basev_AVX_512;           /* offset for scores                                         */
  __m512i ceilingv_AVX_512;        /* saturated simd value used to test for overflow            */
  __m512i tempv_AVX_512;           /* work vector                                               */
  int Q_AVX_512        = P7_NVB_AVX_512(om->M);    /* segment length: # of vectors                              */
  int q_AVX_512;                   /* counter over vectors 0..nq-1                              */
  int i;                           /* counter over sequence positions 1..L                      */
  int     cmp;
  int     status;

  /* Contract checks */
  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  /* note however that it makes no sense to run MSV w/ a model in glocal mode                            */

  /* Try highly optimized Knudsen SSV filter first. 
   * Note that SSV doesn't use any main memory (from <ox>) at all! 
  */
  status = p7_SSVFilter_avx512(dsq, L, om, ret_sc);
  if (status != eslENORESULT) return status;

  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of MSV filter matrix failed");
  dp_AVX_512 = ox->dp_AVX_512; /* and this */

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_MSVFILTER;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.  */
  biasv_AVX_512 = _mm512_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
  for (q_AVX_512 = 0; q_AVX_512 < Q_AVX_512; q_AVX_512++) dp_AVX_512[q_AVX_512] = _mm512_setzero_si512();
  /* saturate simd register for overflow test */
  // ceilingv_AVX_512 = _mm512_cmpeq_epi8(biasv_AVX_512, biasv_AVX_512);
  ceilingv_AVX_512 = _mm512_set1_epi8(0xff);
  basev_AVX_512    = _mm512_set1_epi8((int8_t) om->base_b);
  tjbmv_AVX_512    = _mm512_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv_AVX_512     = _mm512_set1_epi8((int8_t) om->tec_b);
  xJv_AVX_512      = _mm512_subs_epu8(biasv_AVX_512, biasv_AVX_512);
  xBv_AVX_512      = _mm512_subs_epu8(basev_AVX_512, tjbmv_AVX_512);

#if eslDEBUGLEVEL > 0
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
      rsc_AVX_512 = om->rbv_AVX_512[dsq[i]];
      xEv_AVX_512 = _mm512_setzero_si512();      

      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically, which is our -infinity.
       */
       __m512i dp_temp_AVX_512 = dp_AVX_512[Q_AVX_512 -1];
 
       // left_shift dp_temp by 128 bits by shuffling and then inzerting zeroes at the low end
       __m512i temp_mask_AVX_512 = _mm512_shuffle_i32x4(dp_temp_AVX_512, dp_temp_AVX_512, 0x90);
       __m128i zero128 = _mm_setzero_si128();
       temp_mask_AVX_512 = _mm512_inserti32x4(temp_mask_AVX_512, zero128, 0);

       //now do the same merge and right-shift trick we used with AVX to create a left-shift by one byte
       mpv_AVX_512 = _mm512_alignr_epi8(dp_temp_AVX_512, temp_mask_AVX_512,15);

       for (q_AVX_512 = 0; q_AVX_512 < Q_AVX_512; q_AVX_512++)
	 {
	   /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
	   sv_AVX_512   = _mm512_max_epu8(mpv_AVX_512, xBv_AVX_512);
	   sv_AVX_512   = _mm512_adds_epu8(sv_AVX_512, biasv_AVX_512);
	   sv_AVX_512   = _mm512_subs_epu8(sv_AVX_512, *rsc_AVX_512);   rsc_AVX_512++;
	   xEv_AVX_512  = _mm512_max_epu8(xEv_AVX_512, sv_AVX_512);

	   mpv_AVX_512   = dp_AVX_512[q_AVX_512];      /* Load {MDI}(i-1,q) into mpv */
	   dp_AVX_512[q_AVX_512] = sv_AVX_512;           /* Do delayed store of M(i,q) now that memory is usable */
	 }

       /* test for the overflow condition */
       tempv_AVX_512 = _mm512_adds_epu8(xEv_AVX_512, biasv_AVX_512);
       __mmask64 compare_mask = _mm512_cmpeq_epi8_mask(tempv_AVX_512, ceilingv_AVX_512);

      /* Now the "special" states, which start from Mk->E (->C, ->J->B)
       * Use shuffles instead of shifts so when the last max has completed,
       * the last four elements of the simd register will contain the
       * max value.  Then the last shuffle will broadcast the max value
       * to all simd elements.
       */

       // Have to do horizontal max as binary reduction 
      __m256i temp3 = _mm512_extracti32x8_epi32(xEv_AVX_512, 0);
      __m256i temp4 = _mm512_extracti32x8_epi32(xEv_AVX_512, 1);
      temp3 = _mm256_max_epu8(temp3, temp4);
      temp4 = _mm256_permute2x128_si256(temp3, temp3, 0x01);
      // Swap the 128-bit halves from xEv_AVX into temp1

      temp3 = _mm256_max_epu8(temp4, temp3); // each 8-bit field in temp3 now has the max of the
      //corresponding fields in the high and low halves of Dmaxv_AVX

      temp4 = _mm256_shuffle_epi32(temp3, 0x4e);  // Swap the 64-bit halves of each 128-bit half of temp3
      temp3 = _mm256_max_epu8(temp4, temp3);  // Each 64-bit quantity in temp2 now has the max of the corresponding
      // 8-bit fields from the 64-bit quarters of Dmaxv_AVX

      temp4 = _mm256_shuffle_epi32(temp3, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp3
      temp3 = _mm256_max_epu8(temp4, temp3);  // Each 32-bit quantity in temp2 now has the max of the corresponding
      // 8 bit fields from the 32-bit eighths of Dmaxv_AVX

      temp4 = _mm256_shufflelo_epi16(temp3, 0xb1); // bottom 32-bits of temp4 now contain the swapped 16-bit halves
      // of the low 32 bits of temp3
      temp3 = _mm256_max_epu8(temp4, temp3);  //bottom 16 bits of temp3 now contain the max of the 16-bit fields of Dmaxv_AVX

      uint8_t temp_stash2 = _mm256_extract_epi8(temp3, 1);
      temp4 = _mm256_insert_epi8(temp3, temp_stash2, 0);  // low byte of temp4 now has byte 2 of temp3
      temp3 = _mm256_max_epu8(temp4, temp3);  //bottom 16 bits of temp3 now contain the max of the 16-bit fields of Dmaxv_AVX
      temp_stash2 = _mm256_extract_epi8(temp3, 0);  // get low byte of temp3
      

      xEv_AVX_512 = _mm512_set1_epi8(temp_stash2);  // broadcast the max byte from original xEv_AVX
        // to all bytes of xEv_AVX_512

      /* immediately detect overflow */
      if (compare_mask) { *ret_sc = eslINFINITY; return eslERANGE; }

      xEv_AVX_512 = _mm512_subs_epu8(xEv_AVX_512, tecv_AVX_512);
      xJv_AVX_512 = _mm512_max_epu8(xJv_AVX_512,xEv_AVX_512);
      xBv_AVX_512 = _mm512_max_epu8(basev_AVX_512, xJv_AVX_512);
      xBv_AVX_512 = _mm512_subs_epu8(xBv_AVX_512, tjbmv_AVX_512);
      //     printf("Ending AVX_512 code in loop\n");

#if eslDEBUGLEVEL > 0
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

  //need to do this in two steps because AVX-512 doesn't have extract byte
  __m256i xJtemp = _mm512_extracti32x8_epi32(xJv_AVX_512, 0);
  uint8_t xJ_AVX_512 = _mm256_extract_epi8(xJtemp, 0);

  xJ =  xJ_AVX_512;

  *ret_sc = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
    /*     MSV_end_time = __rdtsc();
	   MSV_time += (MSV_end_time - MSV_start_time); */
  return eslOK;
}
/*------------------ end, p7_MSVFilter_avx512() ------------------------*/




#else // ! eslENABLE_AVX512

/* Standard compiler-pleasing mantra for an #ifdef'd-out, empty code file. */
void p7_msvfilter_avx512_silence_hack(void) { return; }
#if defined p7MSVFILTER_AVX512_TESTDRIVE || p7MSVFILTER_AVX512_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_AVX512 or not




