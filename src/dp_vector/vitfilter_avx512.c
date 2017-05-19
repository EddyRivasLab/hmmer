/* Viterbi filter implementation; AVX-512 version.
 * 
 * This is a SIMD vectorized, striped, interleaved, one-row, reduced
 * precision (epi16) implementation of the Viterbi algorithm.
 * 
 * It calculates a close approximation of the Viterbi score, in
 * limited precision (signed words: 16 bits) and range. It may overflow on
 * high scoring sequences, but this indicates that the sequence is a
 * high-scoring hit worth examining more closely anyway.  It will not
 * underflow, in local alignment mode.
 */
#include "p7_config.h"
#ifdef eslENABLE_AVX512

#include <stdio.h>
#include <math.h>

#include <x86intrin.h>

#include "easel.h"
#include "esl_avx512.h"

#include "esl_gumbel.h"

#include "base/p7_hmmwindow.h"
#include "search/p7_pipeline.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"


/*****************************************************************
 * 1. Viterbi filter implementation.
 *****************************************************************/

/* Function:  p7_ViterbiFilter_avx512()
 * Synopsis:  Calculates Viterbi score, vewy vewy fast, in limited precision.
 *
 * Purpose:   Calculates an approximation of the Viterbi score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and a preallocated one-row DP matrix <ox>. Return the 
 *            estimated Viterbi score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow. 
 *            
 *            <ox> will be resized if needed. It's fine if it was just
 *            <_Reuse()'d> from a previous, smaller profile comparison.
 *            
 *            The model must be in a local alignment mode; other modes
 *            cannot provide the necessary guarantee of no underflow.
 *            
 *            This is a striped SIMD Viterbi implementation using Intel
 *            SSE/SSE2 integer intrinsics \citep{Farrar07}, in reduced
 *            precision (signed words, 16 bits).
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: Viterbi score (in nats)          
 *
 * Returns:   <eslOK> on success;
 *            <eslERANGE> if the score overflows; in this case
 *            <*ret_sc> is <eslINFINITY>, and the sequence can 
 *            be treated as a high-scoring hit.
 *            <ox> may be reallocated.
 *
 * Xref:      [Farrar07] for ideas behind striped SIMD DP.
 *            J2/46-47 for layout of HMMER's striped SIMD DP.
 *            J2/50 for single row DP.
 *            J2/60 for reduced precision (epu8)
 *            J2/65 for initial benchmarking
 *            J2/66 for precision maximization
 *            J4/138-140 for reimplementation in 16-bit precision
 *            J9/110-111 for reimplementation with P7_FILTERMX, memory share w/ checkpointed DP matrix
 *            J10/101 for separating P7_FILTERMX from P7_CHECKPTMX again: don't share these
 */
int
p7_ViterbiFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{
  int i;                                                  // counter over sequence positions 1..L
  register __m512i mpv_AVX_512, dpv_AVX_512, ipv_AVX_512; // previous row values
  register __m512i sv_AVX_512;                            // temp storage of 1 curr row value in progress
  register __m512i dcv_AVX_512;                           // delayed storage of D(i,q+1)
  register __m512i xEv_AVX_512;                           // E state: keeps max for Mk->E as we go
  register __m512i xBv_AVX_512;                           // B state: splatted vector of B[i-1] for B->Mk calculations 
  register __m512i Dmaxv_AVX_512;                         // keeps track of maximum D cell on row                      
  int16_t  xE_AVX_512, xB_AVX_512, xC_AVX_512, xJ_AVX_512, xN_AVX_512;     // special states' scores
  int16_t  Dmax_AVX_512;                                  // maximum D cell score on row
  int q_AVX_512;                                          // counter over vectors 0..nq-1                           
  int Q_AVX_512        = P7_NVW_AVX_512(om->M);           // segment length: # of vectors 
  __m512i *dp_AVX_512;
  __m512i *rsc_AVX_512;                                   // will point at om->ru[x] for residue x[i]
  __m512i *tsc_AVX_512;                                   // will point into (and step thru) om->tu 
  __m512i negInfv_AVX_512;
  int     status;

  /* Contract checks */
  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  /* note however that ViterbiFilter numerics are only guaranteed for local alignment, not glocal        */

  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of Vit filter matrix failed");

  dp_AVX_512 = ox->dp_AVX_512;           /* enables MMX_AVX_512f(), IMX_AVX_512f(), DMX_AVX_512f() access macros. Must be set AFTER the GrowTo, because ox->dp may get moved */
  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */

  ox->M    = om->M;
  ox->type = p7F_VITFILTER;

  /* -infinity is -32768 */
  // Want -32768 in the low 16 bits of negInfv_AVX
  __m256i dummy1 = _mm256_setzero_si256();
  dummy1 = _mm256_insert_epi16(dummy1, -32768, 0);
  negInfv_AVX_512 = _mm512_setzero_si512();
  negInfv_AVX_512 = _mm512_inserti64x4(negInfv_AVX_512, dummy1, 0);
  /* Initialization. In unsigned arithmetic, -infinity is -32768
   */
  for (q_AVX_512 = 0; q_AVX_512 < Q_AVX_512; q_AVX_512++)
    MMX_AVX_512f(q_AVX_512) = IMX_AVX_512f(q_AVX_512) = DMX_AVX_512f(q_AVX_512) = _mm512_set1_epi16(-32768);
  xN_AVX_512   = om->base_w;
  xB_AVX_512   = xN_AVX_512 + om->xw[p7O_N][p7O_MOVE];
  xJ_AVX_512   = -32768;
  xC_AVX_512   = -32768;
  xE_AVX_512   = -32768;

#if eslDEBUGLEVEL > 0
  if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, 0, xE, 0, xJ, xB, xC); /* first 0 is <rowi>: do header. second 0 is xN: always 0 here. */
#endif

  for (i = 1; i <= L; i++)
    {

      rsc_AVX_512   = om->rwv_AVX_512[dsq[i]];
      tsc_AVX_512   = om->twv_AVX_512;
      dcv_AVX_512   = _mm512_set1_epi16(-32768);      /* "-infinity" */
      xEv_AVX_512   = _mm512_set1_epi16(-32768);     
      Dmaxv_AVX_512 = _mm512_set1_epi16(-32768);     
      xBv_AVX_512   = _mm512_set1_epi16(xB_AVX_512);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically; replace it with -32768.
       */
    
      mpv_AVX_512 = MMX_AVX_512f(Q_AVX_512-1);  
     
      mpv_AVX_512 = esl_avx512_leftshift_two(mpv_AVX_512);

      mpv_AVX_512 = _mm512_or_si512(mpv_AVX_512, negInfv_AVX_512);
      
      dpv_AVX_512 = DMX_AVX_512f(Q_AVX_512-1);  
      
      dpv_AVX_512 = esl_avx512_leftshift_two(dpv_AVX_512);

      dpv_AVX_512 = _mm512_or_si512(dpv_AVX_512, negInfv_AVX_512);
      
      ipv_AVX_512 = IMX_AVX_512f(Q_AVX_512-1);  
         //left-shift macro
      ipv_AVX_512 = esl_avx512_leftshift_two(ipv_AVX_512);
      ipv_AVX_512 = _mm512_or_si512(ipv_AVX_512, negInfv_AVX_512);

      for (q_AVX_512 = 0; q_AVX_512 < Q_AVX_512; q_AVX_512++)
      {
    
        /* Calculate new MMXf(i,q); don't store it yet, hold it in sv. */
        sv_AVX_512   =                    _mm512_adds_epi16(xBv_AVX_512, *tsc_AVX_512);  tsc_AVX_512++;
        sv_AVX_512   = _mm512_max_epi16 (sv_AVX_512, _mm512_adds_epi16(mpv_AVX_512, *tsc_AVX_512)); tsc_AVX_512++;
        sv_AVX_512   = _mm512_max_epi16 (sv_AVX_512, _mm512_adds_epi16(ipv_AVX_512, *tsc_AVX_512)); tsc_AVX_512++;
        sv_AVX_512   = _mm512_max_epi16 (sv_AVX_512, _mm512_adds_epi16(dpv_AVX_512, *tsc_AVX_512)); tsc_AVX_512++;
        sv_AVX_512   = _mm512_adds_epi16(sv_AVX_512, *rsc_AVX_512);                      rsc_AVX_512++;
        xEv_AVX_512  = _mm512_max_epi16(xEv_AVX_512, sv_AVX_512);

        /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
         * {MDI}MX(q) is then the current, not the prev row
         */
        mpv_AVX_512 = MMX_AVX_512f(q_AVX_512);
        dpv_AVX_512 = DMX_AVX_512f(q_AVX_512);
        ipv_AVX_512 = IMX_AVX_512f(q_AVX_512);

        /* Do the delayed stores of {MD}(i,q) now that memory is usable */
        MMX_AVX_512f(q_AVX_512) = sv_AVX_512;
        DMX_AVX_512f(q_AVX_512) = dcv_AVX_512;

        /* Calculate the next D(i,q+1) partially: M->D only;
               * delay storage, holding it in dcv
         */
        dcv_AVX_512   = _mm512_adds_epi16(sv_AVX_512, *tsc_AVX_512);  tsc_AVX_512++;
        Dmaxv_AVX_512 = _mm512_max_epi16(dcv_AVX_512, Dmaxv_AVX_512);

        /* Calculate and store I(i,q) */
        sv_AVX_512     =                    _mm512_adds_epi16(mpv_AVX_512, *tsc_AVX_512);  tsc_AVX_512++;
        IMX_AVX_512f(q_AVX_512)= _mm512_max_epi16 (sv_AVX_512, _mm512_adds_epi16(ipv_AVX_512, *tsc_AVX_512)); tsc_AVX_512++;
      }

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      //Now, find the horizontal max of the 16-bit values in xEv_AVX.  This version is all
      //AVX instructions to avoid the AVX->SSE switch penalty

      __m512i temp1_AVX_512 = _mm512_shuffle_i64x2(xEv_AVX_512, xEv_AVX_512, 0x4e);
      temp1_AVX_512 = _mm512_max_epi16(xEv_AVX_512, temp1_AVX_512);  // get max of corresponding 16-bit quantities in high and 
      // low halves of xEv

      __m256i temp3_AVX = _mm512_extracti64x4_epi64(temp1_AVX_512, 0);  //shift to normal AVX for 16-bit operations
      __m256i temp4_AVX = _mm256_permute2x128_si256(temp3_AVX, temp3_AVX, 0x01);
      // Swap the 128-bit halves from temp3 into temp4

      temp3_AVX = _mm256_max_epi16(temp3_AVX, temp4_AVX); // each 16-bit field in xEv_AVX now has the max of the
      //corresponding fields in the high and low halves of xEv_AVX

      temp4_AVX = _mm256_shuffle_epi32(temp3_AVX, 0x4e);  // Swap the 64-bit halves of each 128-bit half of temp3_AVX
      temp3_AVX = _mm256_max_epi16(temp4_AVX, temp3_AVX);  // Each 64-bit quantity in temp4 now has the max of the corresponding
      // 16-bit fields from the 64-bit eighths of xEv_AVX_512

      temp4_AVX = _mm256_shuffle_epi32(temp3_AVX, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp3_AVX
      temp3_AVX = _mm256_max_epi16(temp4_AVX, temp3_AVX);  // Each 32-bit quantity in temp2 now has the max of the corresponding
      // 16 bit fields from the 32-bit sixteenths of xEv_AVX_512

      temp4_AVX = _mm256_shufflelo_epi16(temp3_AVX, 0xb1); // bottom 32-bits of temp1_AVX now contain the swapped 16-bit halves
      // of the low 32 bits of temp3_AVX
      temp3_AVX = _mm256_max_epi16(temp4_AVX, temp3_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of xEv_AVX

      xE_AVX_512 = _mm256_extract_epi16(temp3_AVX, 0);  // extract that byte into retval_AVX

      if (xE_AVX_512 >= 32767) {*ret_sc = eslINFINITY; return eslERANGE; } /* immediately detect overflow */
      xN_AVX_512 = xN_AVX_512 + om->xw[p7O_N][p7O_LOOP];
      xC_AVX_512 = ESL_MAX(xC_AVX_512 + om->xw[p7O_C][p7O_LOOP], xE_AVX_512 + om->xw[p7O_E][p7O_MOVE]);
      xJ_AVX_512 = ESL_MAX(xJ_AVX_512 + om->xw[p7O_J][p7O_LOOP], xE_AVX_512 + om->xw[p7O_E][p7O_LOOP]);
      xB_AVX_512 = ESL_MAX(xJ_AVX_512 + om->xw[p7O_J][p7O_MOVE], xN_AVX_512 + om->xw[p7O_N][p7O_MOVE]);
      /* and now xB will carry over into next i, and xC carries over after i=L */

      /* Finally the "lazy F" loop (sensu [Farrar07]). We can often
       * prove that we don't need to evaluate any D->D paths at all.
       *
       * The observation is that if we can show that on the next row,
       * B->M(i+1,k) paths always dominate M->D->...->D->M(i+1,k) paths
       * for all k, then we don't need any D->D calculations.
       * 
       * The test condition is:
       *      max_k D(i,k) + max_k ( TDD(k-2) + TDM(k-1) - TBM(k) ) < xB(i)
       * So:
       *   max_k (TDD(k-2) + TDM(k-1) - TBM(k)) is precalc'ed in om->dd_bound;
       *   max_k D(i,k) is why we tracked Dmaxv;
       *   xB(i) was just calculated above.
       */
      
      // compute Dmax_AVX = horizontal 16-bit max of Dmaxv_AVX
      temp1_AVX_512 = _mm512_shuffle_i64x2(Dmaxv_AVX_512, Dmaxv_AVX_512, 0x4e);
      temp1_AVX_512 = _mm512_max_epi16(Dmaxv_AVX_512, temp1_AVX_512);  // get max of corresponding 16-bit quantities in high and 
      // low halves of Dmaxv

      temp3_AVX = _mm512_extracti64x4_epi64(temp1_AVX_512, 0);  //shift to normal AVX for 16-bit operations
      temp4_AVX = _mm256_permute2x128_si256(temp3_AVX, temp3_AVX, 0x01);
      // Swap the 128-bit halves from Dmaxv_AVX into temp1

      temp3_AVX = _mm256_max_epi16(temp4_AVX, temp3_AVX); // each 16-bit field in Dmaxv_AVX now has the max of the
      //corresponding fields in the high and low halves of Dmaxv_AVX

      temp4_AVX = _mm256_shuffle_epi32(temp3_AVX, 0x8e);  // Swap the 64-bit halves of each 128-bit half of temp3_AVX
      temp3_AVX = _mm256_max_epi16(temp3_AVX, temp4_AVX);  // Each 64-bit quantity in temp2 now has the max of the corresponding
      // 16-bit fields from the 64-bit quarters of Dmaxv_AVX

      temp4_AVX = _mm256_shuffle_epi32(temp3_AVX, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp2_AVX
      temp3_AVX = _mm256_max_epi16(temp3_AVX, temp4_AVX);  // Each 32-bit quantity in temp2 now has the max of the corresponding
      // 16 bit fields from the 32-bit eighths of Dmaxv_AVX

      temp4_AVX = _mm256_shufflelo_epi16(temp3_AVX, 0xb1); // bottom 32-bits of temp1_AVX now contain the swapped 16-bit halves
      // of the low 32 bits of temp2_AVX
      temp3_AVX = _mm256_max_epi16(temp3_AVX, temp4_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of Dmaxv_AVX

      Dmax_AVX_512 = _mm256_extract_epi16(temp3_AVX, 0);  // extract that byte into retval_AVX

      if (Dmax_AVX_512 + om->ddbound_w > xB_AVX_512) 
	{
	  /* Now we're obligated to do at least one complete DD path to be sure. */
	  /* dcv has carried through from end of q loop above */
	  dcv_AVX_512 = esl_avx512_leftshift_two(dcv_AVX_512);

	  dcv_AVX_512 = _mm512_or_si512(dcv_AVX_512, negInfv_AVX_512);
	  tsc_AVX_512 = om->twv_AVX_512 + 7*Q_AVX_512;  /* set tsc to start of the DD's */
	  for (q_AVX_512 = 0; q_AVX_512 < Q_AVX_512; q_AVX_512++) 
	    {
	      DMX_AVX_512f(q_AVX_512) = _mm512_max_epi16(dcv_AVX_512, DMX_AVX_512f(q_AVX_512));  
	      dcv_AVX_512     = _mm512_adds_epi16(DMX_AVX_512f(q_AVX_512), *tsc_AVX_512); tsc_AVX_512++;
	    }

	  /* We may have to do up to three more passes; the check
	   * is for whether crossing a segment boundary can improve
	   * our score. 
	   */
	  do {
	    dcv_AVX_512 = esl_avx512_leftshift_two(dcv_AVX_512);
	    dcv_AVX_512 = _mm512_or_si512(dcv_AVX_512, negInfv_AVX_512);
	    tsc_AVX_512 = om->twv_AVX_512 + 7*Q_AVX_512;  /* set tsc to start of the DD's */
	    for (q_AVX_512 = 0; q_AVX_512 < Q_AVX_512; q_AVX_512++) 
	      {
		__mmask32 check2 = _mm512_cmpgt_epi16_mask(dcv_AVX_512, DMX_AVX_512f(q_AVX_512));
		// continue if any of the values were greater than the old ones
		if (check2 == 0) break;
		DMX_AVX_512f(q_AVX_512) = _mm512_max_epi16(dcv_AVX_512, DMX_AVX_512f(q_AVX_512));  
		dcv_AVX_512   = _mm512_adds_epi16(DMX_AVX_512f(q_AVX_512), *tsc_AVX_512); tsc_AVX_512++;
	      }     
	  } while (q_AVX_512 == Q_AVX_512);
	}
      else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
	{
	  dcv_AVX_512 = esl_avx512_leftshift_two(dcv_AVX_512);
	  DMX_AVX_512f(0) = _mm512_or_si512(dcv_AVX_512, negInfv_AVX_512);
	}

#if eslDEBUGLEVEL > 0
      if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, i, xE, 0, xJ, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T */
  if (xC_AVX_512 > -32768)
    {
      *ret_sc = (float) xC_AVX_512 + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
      /* *ret_sc += L * om->ncj_roundoff;  see J4/150 for rationale: superceded by -3.0nat approximation*/
      *ret_sc /= om->scale_w;
      *ret_sc -= 3.0; /* the NN/CC/JJ=0,-3nat approximation: see J5/36. That's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ contrib */
    }
  else  *ret_sc = -eslINFINITY;
  return eslOK;
}


#else // ! eslENABLE_AVX512

/* Standard compiler-pleasing mantra for an #ifdef'd-out, empty code file. */
void p7_vitfilter_avx512_silence_hack(void) { return; }
#if defined p7VITFILTER_AVX512_TESTDRIVE || p7VITFILTER_AVX512_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_AVX512 or not
