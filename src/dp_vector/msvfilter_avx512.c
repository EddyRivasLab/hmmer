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
#if p7_CPU_ARCH == intel
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#include "x86intrin.h"
#ifdef HAVE_AVX512
  #include <immintrin.h>  /* AVX-512 */
#endif
#endif /* x86 */
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
p7_MSVFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
#ifdef HAVE_AVX512 
   uint8_t  xJ;                     /* special states' scores                  */
  register __m512i mpv_AVX_512;            /* previous row values                                       */
  register __m512i xEv_AVX_512;      /* E state: keeps max for Mk->E as we go                     */
  register __m512i xBv_AVX_512;      /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m512i sv_AVX_512;       /* temp storage of 1 curr row value in progress              */
  register __m512i biasv_AVX_512;    /* emission bias in a vector                                 */
  __m512i *dp_AVX_512;               /* the dp row memory                                         */
  __m512i *rsc_AVX_512;        /* will point at om->rbv[x] for residue x[i]                 */
  __m512i xJv_AVX_512;                     /* vector for states score                                   */
  __m512i tjbmv_AVX_512;                   /* vector for cost of moving {JN}->B->M                      */
  __m512i tecv_AVX_512;                    /* vector for E->C  cost                                     */
  __m512i basev_AVX_512;                   /* offset for scores                                         */
  __m512i ceilingv_AVX_512;                /* saturated simd value used to test for overflow            */
  __m512i tempv_AVX_512;                   /* work vector                                               */
  int Q_AVX_512        = P7_NVB_AVX_512(om->M);    /* segment length: # of vectors                              */
  int q_AVX_512;         /* counter over vectors 0..nq-1                              */

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
  #endif  //HAVE_AVX512
  #ifndef HAVE_AVX512
  return eslENORESULT;  // Stub so we have something to link if we build without AVX512 support
  #endif
  }

/*------------------ end, p7_MSVFilter() ------------------------*/


/* Function:  p7_SSVFilter_longtarget()
 * Synopsis:  Finds windows with SSV scores above some threshold (vewy vewy fast, in limited precision)
 *
 * Purpose:   Calculates an approximation of the SSV (single ungapped diagonal)
 *            score for regions of sequence <dsq> of length <L> residues, using
 *            optimized profile <om>, and a preallocated one-row DP matrix <ox>,
 *            and captures the positions at which such regions exceed the score
 *            required to be significant in the eyes of the calling function,
 *            which depends on the <bg> and <p> (usually p=0.02 for nhmmer).
 *            Note that this variant performs only SSV computations, never
 *            passing through the J state - the score required to pass SSV at
 *            the default threshold (or less restrictive) is sufficient to
 *            pass MSV in essentially all DNA models we've tested.
 *
 *            Above-threshold diagonals are captured into a preallocated list
 *            <windowlist>. Rather than simply capturing positions at which a
 *            score threshold is reached, this function establishes windows
 *            around those high-scoring positions, using scores in <msvdata>.
 *            These windows can be merged by the calling function.
 *
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues
 *            om      - optimized profile
 *            ox      - DP matrix
 *            msvdata    - compact representation of substitution scores, for backtracking diagonals
 *            bg         - the background model, required for translating a P-value threshold into a score threshold
 *            P          - p-value below which a region is captured as being above threshold
 *            windowlist - preallocated container for all hits (resized if necessary)
 *
 *
 * Note:      We misuse the matrix <ox> here, using only a third of the
 *            first dp row, accessing it as <dp[0..Q-1]> rather than
 *            in triplets via <{MDI}MX(q)> macros, since we only need
 *            to store M state values. We know that if <ox> was big
 *            enough for normal DP calculations, it must be big enough
 *            to hold the MSVFilter calculation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_SSVFilter_longtarget_avx512(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_FILTERMX *ox, const P7_SCOREDATA *msvdata,
                        P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist)
{
  register __m128i mpv;            /* previous row values                                       */
  register __m128i xEv;            /* E state: keeps max for Mk->E for a single iteration       */
  register __m128i xBv;            /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i sv;             /* temp storage of 1 curr row value in progress              */
  register __m128i biasv;          /* emission bias in a vector                                 */
  int i;                           /* counter over sequence positions 1..L                      */
  int q;                           /* counter over vectors 0..nq-1                              */
  int Q        = P7_NVB(om->M);    /* segment length: # of vectors                              */
  __m128i *dp  = ox->dp;           /* one DP row of memory                                      */
  __m128i *rsc;                    /* will point at om->rbv[x] for residue x[i]                 */
  __m128i tjbmv;                   /* vector for J->B move cost + B->M move costs               */
  __m128i basev;                   /* offset for scores                                         */
  __m128i ceilingv;                /* saturated simd value used to test for overflow            */
  __m128i tempv;                   /* work vector                                               */
  int cmp;
  int k;
  int n;
  int end;
  int rem_sc;
  int start;
  int target_end;
  int target_start;
  int max_end;
  int max_sc;
  int sc;
  int pos_since_max;
  float ret_sc;
  union { __m128i v; uint8_t b[16]; } u;
  int   status;


  /*
   * Computing the score required to let P meet the F1 prob threshold
   * In original code, converting from a scaled int MSV
   * score S (the score getting to state E) to a probability goes like this:
   *  usc =  S - om->tec_b - om->tjb_b - om->base_b;
   *  usc /= om->scale_b;
   *  usc -= 3.0;
   *  P = f ( (usc - nullsc) / eslCONST_LOG2 , mu, lambda)
   * and we're computing the threshold usc, so reverse it:
   *  (usc - nullsc) /  eslCONST_LOG2 = inv_f( P, mu, lambda)
   *  usc = nullsc + eslCONST_LOG2 * inv_f( P, mu, lambda)
   *  usc += 3
   *  usc *= om->scale_b
   *  S = usc + om->tec_b + om->tjb_b + om->base_b
   *
   *  Here, I compute threshold with length model based on max_length.  Doesn't
   *  matter much - in any case, both the bg and om models will change with roughly
   *  1 bit for each doubling of the length model, so they offset.
   */
  float nullsc;
  __m128i sc_threshv;
  uint8_t sc_thresh;
  float invP = esl_gumbel_invsurv(P, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);

  p7_bg_SetLength(bg, om->max_length);
  p7_oprofile_ReconfigMSVLength(om, om->max_length);
  p7_bg_NullOne  (bg, dsq, om->max_length, &nullsc);

  sc_thresh = (int) ceil( ( ( nullsc  + (invP * eslCONST_LOG2) + 3.0 )  * om->scale_b ) + om->base_b +  om->tec_b  + om->tjb_b );
  sc_threshv = _mm_set1_epi8((int8_t) 255 - sc_thresh);


  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of SSV filter matrix failed");

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_SSVFILTER;

  /* Initialization. In offset unsigned  arithmetic, -infinity is 0, and 0 is om->base.
   */
  biasv = _mm_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
  ceilingv = _mm_cmpeq_epi8(biasv, biasv);
  for (q = 0; q < Q; q++) dp[q] = _mm_setzero_si128();

  basev = _mm_set1_epi8((int8_t) om->base_b);
  tjbmv = _mm_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);

  xBv = _mm_subs_epu8(basev, tjbmv);
  for (i = 1; i <= L; i++) {
    rsc = om->rbv[dsq[i]];
    xEv = _mm_setzero_si128();

          /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12.
           * Because ia32 is littlendian, this means a left bit shift.
           * Zeros shift on automatically, which is our -infinity.
           */
          mpv = _mm_slli_si128(dp[Q-1], 1);
          for (q = 0; q < Q; q++) {
                  /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
                  sv   = _mm_max_epu8(mpv, xBv);
                  sv   = _mm_adds_epu8(sv, biasv);
                  sv   = _mm_subs_epu8(sv, *rsc);   rsc++;
                  xEv  = _mm_max_epu8(xEv, sv);

                  mpv   = dp[q];          /* Load {MDI}(i-1,q) into mpv */
                  dp[q] = sv;             /* Do delayed store of M(i,q) now that memory is usable */
          }

          /* test if the pthresh significance threshold has been reached;
           * note: don't use _mm_cmpgt_epi8, because it's a signed comparison, which won't work on uint8s */
          tempv = _mm_adds_epu8(xEv, sc_threshv);
          tempv = _mm_cmpeq_epi8(tempv, ceilingv);
          cmp = _mm_movemask_epi8(tempv);

          if (cmp != 0) {  //hit pthresh, so add position to list and reset values

            //figure out which model state hit threshold
            end = -1;
            rem_sc = -1;
            for (q = 0; q < Q; q++) {  /// Unpack and unstripe, so we can find the state that exceeded pthresh
          u.v = dp[q];
          for (k = 0; k < 16; k++) { // unstripe
            //(q+Q*k+1) is the model position k at which the xE score is found
            if (u.b[k] >= sc_thresh && u.b[k] > rem_sc && (q+Q*k+1) <= om->M) {
              end = (q+Q*k+1);
              rem_sc = u.b[k];
            }
          }
          dp[q] = _mm_set1_epi8(0); // while we're here ... this will cause values to get reset to xB in next dp iteration
            }

            //recover the diagonal that hit threshold
            start = end;
            target_end = target_start = i;
            sc = rem_sc;
            while (rem_sc > om->base_b - om->tjb_b - om->tbm_b) {
              rem_sc -= om->bias_b -  msvdata->msv_scores[start*om->abc->Kp + dsq[target_start]];
              --start;
              --target_start;
            }
            start++;
            target_start++;


            //extend diagonal further with single diagonal extension
            k = end+1;
            n = target_end+1;
            max_end = target_end;
            max_sc = sc;
            pos_since_max = 0;
            while (k<om->M && n<=L) {
              sc += om->bias_b -  msvdata->msv_scores[k*om->abc->Kp + dsq[n]];

              if (sc >= max_sc) {
                max_sc = sc;
                max_end = n;
                pos_since_max=0;
              } else {
                pos_since_max++;
                if (pos_since_max == 5)
                  break;
              }
              k++;
              n++;
            }

            end  +=  (max_end - target_end);
            k    +=  (max_end - target_end);
      target_end = max_end;

      ret_sc = ((float) (max_sc - om->tjb_b) - (float) om->base_b);
      ret_sc /= om->scale_b;
      ret_sc -= 3.0; // that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ

      p7_hmmwindow_new(windowlist, 0, target_start, k, end, end-start+1 , ret_sc, p7_NOCOMPLEMENT);

      i = target_end; // skip forward
          }


  } /* end loop over sequence residues 1..L */
  return eslOK;
#endif /* HAVE_AVX512 */
#ifndef HAVE_AVX512
  return 0;
#endif
}
/*------------------ end, p7_SSVFilter_longtarget() ------------------------*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/



