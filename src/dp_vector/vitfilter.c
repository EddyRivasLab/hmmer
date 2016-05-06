/* Viterbi filter implementation; SSE version.
 * 
 * This is a SIMD vectorized, striped, interleaved, one-row, reduced
 * precision (epi16) implementation of the Viterbi algorithm.
 * 
 * It calculates a close approximation of the Viterbi score, in
 * limited precision (signed words: 16 bits) and range. It may overflow on
 * high scoring sequences, but this indicates that the sequence is a
 * high-scoring hit worth examining more closely anyway.  It will not
 * underflow, in local alignment mode.
 * 
 * Contents:
 *   1. Viterbi filter implementation.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#ifdef p7_build_AVX512
 #include <immintrin.h>
#endif
#include "easel.h"
#include "esl_sse.h"

#ifdef p7_build_AVX2
 #include <immintrin.h>
 #include "esl_avx.h"
#endif 

#include "esl_gumbel.h"

#include "base/p7_hmmwindow.h"
#include "search/p7_pipeline.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_filtermx.h"


/*****************************************************************
 * 1. Viterbi filter implementation.
 *****************************************************************/

/* Function:  p7_ViterbiFilter()
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
p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox, float *ret_sc)
{

//printf("Starting ViterbiFilter\n");
  int i;        /* counter over sequence positions 1..L                      */
  
#ifdef p7_build_SSE
  register __m128i mpv, dpv, ipv;  /* previous row values                                       */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128i dcv;		   /* delayed storage of D(i,q+1)                               */
  register __m128i xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i Dmaxv;          /* keeps track of maximum D cell on row                      */
  int16_t  xE, xB, xC, xJ, xN;	   /* special states' scores                                    */
  int16_t  Dmax;		   /* maximum D cell score on row                               */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = P7_NVW(om->M);    /* segment length: # of vectors                              */
  __m128i *dp;
  __m128i *rsc;			   /* will point at om->ru[x] for residue x[i]                  */
  __m128i *tsc;			   /* will point into (and step thru) om->tu                    */
  __m128i negInfv;
#endif
#ifdef p7_build_AVX2
  register __m256i mpv_AVX, dpv_AVX, ipv_AVX; /*  prvious row values                                       */
  register __m256i sv_AVX;       /* temp storage of 1 curr row value in progress              */
  register __m256i dcv_AVX;      /* delayed storage of D(i,q+1)                               */
  register __m256i xEv_AVX;      /* E state: keeps max for Mk->E as we go                     */
  register __m256i xBv_AVX;      /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m256i Dmaxv_AVX;          /* keeps track of maximum D cell on row                      */
  int16_t  xE_AVX, xB_AVX, xC_AVX, xJ_AVX, xN_AVX;     /* special states' scores                                    */
  int16_t  Dmax_AVX;       /* maximum D cell score on row                               */
  int q_AVX;         /* counter over vectors 0..nq-1                              */
  int Q_AVX        = P7_NVW_AVX(om->M);    /* segment length: # of vectors                              */
  __m256i *dp_AVX;
  __m256i *rsc_AVX;        /* will point at om->ru[x] for residue x[i]                  */
  __m256i *tsc_AVX;        /* will point into (and step thru) om->tu                    */
  __m256i negInfv_AVX;
#endif
#ifdef p7_build_AVX512
  register __m512i mpv_AVX_512, dpv_AVX_512, ipv_AVX_512; /*  prvious row values                                       */
  register __m512i sv_AVX_512;       /* temp storage of 1 curr row value in progress              */
  register __m512i dcv_AVX_512;      /* delayed storage of D(i,q+1)                               */
  register __m512i xEv_AVX_512;      /* E state: keeps max for Mk->E as we go                     */
  register __m512i xBv_AVX_512;      /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m512i Dmaxv_AVX_512;          /* keeps track of maximum D cell on row                      */
  int16_t  xE_AVX_512, xB_AVX_512, xC_AVX_512, xJ_AVX_512, xN_AVX_512;     /* special states' scores                                    */
  int16_t  Dmax_AVX_512;       /* maximum D cell score on row                               */
  int q_AVX_512;         /* counter over vectors 0..nq-1                              */
  int Q_AVX_512        = P7_NVW_AVX_512(om->M);    /* segment length: # of vectors                              */
  __m512i *dp_AVX_512;
  __m512i *rsc_AVX_512;        /* will point at om->ru[x] for residue x[i]                  */
  __m512i *tsc_AVX_512;        /* will point into (and step thru) om->tu                    */
  __m512i negInfv_AVX_512;
#endif
  int     status;

  /* Contract checks */
  ESL_DASSERT1(( om->mode == p7_LOCAL )); /* Production code assumes multilocal mode w/ length model <L> */
  ESL_DASSERT1(( om->L    == L ));	  /*  ... and it's easy to forget to set <om> that way           */
  ESL_DASSERT1(( om->nj   == 1.0f ));	  /*  ... hence the check                                        */
                                          /*  ... which you can disable, if you're playing w/ config     */
  /* note however that ViterbiFilter numerics are only guaranteed for local alignment, not glocal        */


  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of Vit filter matrix failed");
#ifdef p7_build_SSE
  dp = ox->dp;           /* enables MMXf(), IMXf(), DMXf() access macros. Must be set AFTER the GrowTo, because ox->dp may get moved */
#endif
#ifdef p7_build_AVX2
  dp_AVX = ox->dp_AVX;           /* enables MMX_AVXf(), IMX_AVXf(), DMX_AVXf() access macros. Must be set AFTER the GrowTo, because ox->dp may get moved */
#endif
#ifdef p7_build_AVX512
  dp_AVX_512 = ox->dp_AVX_512;           /* enables MMX_AVX_512f(), IMX_AVX_512f(), DMX_AVX_512f() access macros. Must be set AFTER the GrowTo, because ox->dp may get moved */
#endif
  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_VITFILTER;

#ifdef p7_build_SSE
  /* -infinity is -32768 */
  negInfv = _mm_set1_epi16(-32768);
  negInfv = _mm_srli_si128(negInfv, 14);  /* negInfv = 16-byte vector, 14 0 bytes + 2-byte value=-32768, for an OR operation. */

  /* Initialization. In unsigned arithmetic, -infinity is -32768
   */
  for (q = 0; q < Q; q++)
    MMXf(q) = IMXf(q) = DMXf(q) = _mm_set1_epi16(-32768);
  xN   = om->base_w;
  xB   = xN + om->xw[p7O_N][p7O_MOVE];
  xJ   = -32768;
  xC   = -32768;
  xE   = -32768;
#endif
#ifdef p7_build_AVX2
  /* -infinity is -32768 */
  // Want -32768 in the low 16 bits of negInfv_AVX
  negInfv_AVX = _mm256_setzero_si256();
  negInfv_AVX = _mm256_insert_epi16(negInfv_AVX, -32768, 0);
  /* Initialization. In unsigned arithmetic, -infinity is -32768
   */
  for (q_AVX = 0; q_AVX < Q_AVX; q_AVX++)
    MMX_AVXf(q_AVX) = IMX_AVXf(q_AVX) = DMX_AVXf(q_AVX) = _mm256_set1_epi16(-32768);
  xN_AVX   = om->base_w;
  xB_AVX   = xN_AVX + om->xw[p7O_N][p7O_MOVE];
  xJ_AVX   = -32768;
  xC_AVX   = -32768;
  xE_AVX   = -32768;
#endif
#ifdef p7_build_AVX512
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
#endif

#ifdef p7_DEBUGGING
  if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, 0, xE, 0, xJ, xB, xC); /* first 0 is <rowi>: do header. second 0 is xN: always 0 here. */
#endif

  for (i = 1; i <= L; i++)
    {
#ifdef p7_build_SSE // This is somewhat ugly, but putting SSE, AVX, AVX-512 versions in the same loop makes it easier to
// check each iteration         
      rsc   = om->rwv[dsq[i]];
      tsc   = om->twv;
      dcv   = _mm_set1_epi16(-32768);      /* "-infinity" */
      xEv   = _mm_set1_epi16(-32768);     
      Dmaxv = _mm_set1_epi16(-32768);     
      xBv   = _mm_set1_epi16(xB);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically; replace it with -32768.
       */
      mpv = MMXf(Q-1);  mpv = _mm_slli_si128(mpv, 2);  mpv = _mm_or_si128(mpv, negInfv);
      dpv = DMXf(Q-1);  dpv = _mm_slli_si128(dpv, 2);  dpv = _mm_or_si128(dpv, negInfv);
      ipv = IMXf(Q-1);  ipv = _mm_slli_si128(ipv, 2);  ipv = _mm_or_si128(ipv, negInfv);

      for (q = 0; q < Q; q++)
      {
    
        /* Calculate new MMXf(i,q); don't store it yet, hold it in sv. */
        sv   =                    _mm_adds_epi16(xBv, *tsc);  tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(mpv, *tsc)); tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(ipv, *tsc)); tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(dpv, *tsc)); tsc++;
        sv   = _mm_adds_epi16(sv, *rsc);                      rsc++;
        xEv  = _mm_max_epi16(xEv, sv);

        /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
         * {MDI}MX(q) is then the current, not the prev row
         */
        mpv = MMXf(q);
        dpv = DMXf(q);
        ipv = IMXf(q);

        /* Do the delayed stores of {MD}(i,q) now that memory is usable */
        MMXf(q) = sv;
        DMXf(q) = dcv;

        /* Calculate the next D(i,q+1) partially: M->D only;
               * delay storage, holding it in dcv
         */
        dcv   = _mm_adds_epi16(sv, *tsc);  tsc++;
        Dmaxv = _mm_max_epi16(dcv, Dmaxv);

        /* Calculate and store I(i,q) */
        sv     =                    _mm_adds_epi16(mpv, *tsc);  tsc++;
        IMXf(q)= _mm_max_epi16 (sv, _mm_adds_epi16(ipv, *tsc)); tsc++;
      }

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE = esl_sse_hmax_epi16(xEv);
      if (xE >= 32767) { *ret_sc = eslINFINITY; return eslERANGE; }	/* immediately detect overflow */
      xN = xN + om->xw[p7O_N][p7O_LOOP];
      xC = ESL_MAX(xC + om->xw[p7O_C][p7O_LOOP], xE + om->xw[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ + om->xw[p7O_J][p7O_LOOP], xE + om->xw[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ + om->xw[p7O_J][p7O_MOVE], xN + om->xw[p7O_N][p7O_MOVE]);
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
      Dmax = esl_sse_hmax_epi16(Dmaxv);
      if (Dmax + om->ddbound_w > xB) 
	{
   // printf("vitfilter SSE doing DD calc, Dmax = %i\n", Dmax);
	  /* Now we're obligated to do at least one complete DD path to be sure. */
	  /* dcv has carried through from end of q loop above */
	  dcv = _mm_slli_si128(dcv, 2); 
	  dcv = _mm_or_si128(dcv, negInfv);
	  tsc = om->twv + 7*Q;	/* set tsc to start of the DD's */
	  for (q = 0; q < Q; q++) 
	    {
	      DMXf(q) = _mm_max_epi16(dcv, DMXf(q));	
	      dcv     = _mm_adds_epi16(DMXf(q), *tsc); tsc++;
	    }

	  /* We may have to do up to three more passes; the check
	   * is for whether crossing a segment boundary can improve
	   * our score. 
	   */
	  do {
	    dcv = _mm_slli_si128(dcv, 2);
	    dcv = _mm_or_si128(dcv, negInfv);
	    tsc = om->twv + 7*Q;	/* set tsc to start of the DD's */
	    for (q = 0; q < Q; q++) 
	      {
		if (! esl_sse_any_gt_epi16(dcv, DMXf(q))) break;
		DMXf(q) = _mm_max_epi16(dcv, DMXf(q));	
		dcv     = _mm_adds_epi16(DMXf(q), *tsc);   tsc++;
	      }	    
	  } while (q == Q);
	}
      else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
	{
 //   printf("vitfilter SSE skipping DD calc, Dmax = %i\n", Dmax);
	  dcv = _mm_slli_si128(dcv, 2);
	  DMXf(0) = _mm_or_si128(dcv, negInfv);
	}
#endif

#ifdef p7_build_AVX2 // This is somewhat ugly, but putting SSE, AVX, AVX-512 versions in the same loop makes it easier to
// check each iteration         
      rsc_AVX   = om->rwv_AVX[dsq[i]];
      tsc_AVX   = om->twv_AVX;
      dcv_AVX   = _mm256_set1_epi16(-32768);      /* "-infinity" */
      xEv_AVX   = _mm256_set1_epi16(-32768);     
      Dmaxv_AVX = _mm256_set1_epi16(-32768);     
      xBv_AVX   = _mm256_set1_epi16(xB_AVX);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically; replace it with -32768.
       */
    
      mpv_AVX = MMX_AVXf(Q_AVX-1);  
         //left-shift macro
      __m256i temp_mask_AVX = _mm256_permute2x128_si256(mpv_AVX, mpv_AVX, _MM_SHUFFLE(0,0,3,0) );
      mpv_AVX = _mm256_alignr_epi8(mpv_AVX, temp_mask_AVX,14);
      mpv_AVX = _mm256_or_si256(mpv_AVX, negInfv_AVX);
      
      dpv_AVX = DMX_AVXf(Q_AVX-1);  
         //left-shift macro
      temp_mask_AVX = _mm256_permute2x128_si256(dpv_AVX, dpv_AVX, _MM_SHUFFLE(0,0,3,0) );
      dpv_AVX = _mm256_alignr_epi8(dpv_AVX, temp_mask_AVX,14);
       dpv_AVX = _mm256_or_si256(dpv_AVX, negInfv_AVX);
      
      ipv_AVX = IMX_AVXf(Q_AVX-1);  
         //left-shift macro
      temp_mask_AVX = _mm256_permute2x128_si256(ipv_AVX, ipv_AVX, _MM_SHUFFLE(0,0,3,0) );
      ipv_AVX = _mm256_alignr_epi8(ipv_AVX, temp_mask_AVX,14);
       ipv_AVX = _mm256_or_si256(ipv_AVX, negInfv_AVX);

      for (q_AVX = 0; q_AVX < Q_AVX; q_AVX++)
      {
    
        /* Calculate new MMXf(i,q); don't store it yet, hold it in sv. */
        sv_AVX   =                    _mm256_adds_epi16(xBv_AVX, *tsc_AVX);  tsc_AVX++;
        sv_AVX   = _mm256_max_epi16 (sv_AVX, _mm256_adds_epi16(mpv_AVX, *tsc_AVX)); tsc_AVX++;
        sv_AVX   = _mm256_max_epi16 (sv_AVX, _mm256_adds_epi16(ipv_AVX, *tsc_AVX)); tsc_AVX++;
        sv_AVX   = _mm256_max_epi16 (sv_AVX, _mm256_adds_epi16(dpv_AVX, *tsc_AVX)); tsc_AVX++;
        sv_AVX   = _mm256_adds_epi16(sv_AVX, *rsc_AVX);                      rsc_AVX++;
        xEv_AVX  = _mm256_max_epi16(xEv_AVX, sv_AVX);

        /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
         * {MDI}MX(q) is then the current, not the prev row
         */
        mpv_AVX = MMX_AVXf(q_AVX);
        dpv_AVX = DMX_AVXf(q_AVX);
        ipv_AVX = IMX_AVXf(q_AVX);

        /* Do the delayed stores of {MD}(i,q) now that memory is usable */
        MMX_AVXf(q_AVX) = sv_AVX;
        DMX_AVXf(q_AVX) = dcv_AVX;

        /* Calculate the next D(i,q+1) partially: M->D only;
               * delay storage, holding it in dcv
         */
        dcv_AVX   = _mm256_adds_epi16(sv_AVX, *tsc_AVX);  tsc_AVX++;
        Dmaxv_AVX = _mm256_max_epi16(dcv_AVX, Dmaxv_AVX);

        /* Calculate and store I(i,q) */
        sv_AVX     =                    _mm256_adds_epi16(mpv_AVX, *tsc_AVX);  tsc_AVX++;
        IMX_AVXf(q_AVX)= _mm256_max_epi16 (sv_AVX, _mm256_adds_epi16(ipv_AVX, *tsc_AVX)); tsc_AVX++;
      }

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      //Now, find the horizontal max of the 16-bit values in xEv_AVX.  This version is all
      //AVX instructions to avoid the AVX->SSE switch penalty

      __m256i temp1_AVX = _mm256_permute2x128_si256(xEv_AVX, xEv_AVX, 0x01);
      // Swap the 128-bit halves from xEv_AVX into temp1

      __m256i temp2_AVX = _mm256_max_epi16(temp1_AVX, xEv_AVX); // each 16-bit field in xEv_AVX now has the max of the
      //corresponding fields in the high and low halves of xEv_AVX

      temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0x4e);  // Swap the 64-bit halves of each 128-bit half of temp2_AVX
      temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  // Each 64-bit quantity in temp2 now has the max of the corresponding
      // 16-bit fields from the 64-bit quarters of xEv_AVX

      temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp2_AVX
      temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  // Each 32-bit quantity in temp2 now has the max of the corresponding
      // 16 bit fields from the 32-bit eighths of xEv_AVX

      temp1_AVX = _mm256_shufflelo_epi16(temp2_AVX, 0xb1); // bottom 32-bits of temp1_AVX now contain the swapped 16-bit halves
      // of the low 32 bits of temp2_AVX
      temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of xEv_AVX

      xE_AVX = _mm256_extract_epi16(temp2_AVX, 0);  // extract that byte into retval_AVX

      if (xE_AVX >= 32767) {*ret_sc = eslINFINITY; return eslERANGE; } /* immediately detect overflow */
      xN_AVX = xN_AVX + om->xw[p7O_N][p7O_LOOP];
      xC_AVX = ESL_MAX(xC_AVX + om->xw[p7O_C][p7O_LOOP], xE_AVX + om->xw[p7O_E][p7O_MOVE]);
      xJ_AVX = ESL_MAX(xJ_AVX + om->xw[p7O_J][p7O_LOOP], xE_AVX + om->xw[p7O_E][p7O_LOOP]);
      xB_AVX = ESL_MAX(xJ_AVX + om->xw[p7O_J][p7O_MOVE], xN_AVX + om->xw[p7O_N][p7O_MOVE]);
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

      temp1_AVX = _mm256_permute2x128_si256(Dmaxv_AVX, Dmaxv_AVX, 0x01);
      // Swap the 128-bit halves from Dmaxv_AVX into temp1

      temp2_AVX = _mm256_max_epi16(temp1_AVX, Dmaxv_AVX); // each 16-bit field in Dmaxv_AVX now has the max of the
      //corresponding fields in the high and low halves of Dmaxv_AVX

      temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0x4e);  // Swap the 64-bit halves of each 128-bit half of temp2_AVX
      temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  // Each 64-bit quantity in temp2 now has the max of the corresponding
      // 16-bit fields from the 64-bit quarters of Dmaxv_AVX

      temp1_AVX = _mm256_shuffle_epi32(temp2_AVX, 0xb1);  // Swap the 32-bit halves of each 64-bit quarter of temp2_AVX
      temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  // Each 32-bit quantity in temp2 now has the max of the corresponding
      // 16 bit fields from the 32-bit eighths of Dmaxv_AVX

      temp1_AVX = _mm256_shufflelo_epi16(temp2_AVX, 0xb1); // bottom 32-bits of temp1_AVX now contain the swapped 16-bit halves
      // of the low 32 bits of temp2_AVX
      temp2_AVX = _mm256_max_epi16(temp1_AVX, temp2_AVX);  //bottom 16 bits of temp2_AVX now contain the max of the 16-bit fields of Dmaxv_AVX

      Dmax_AVX = _mm256_extract_epi16(temp2_AVX, 0);  // extract that byte into retval_AVX

      if (Dmax_AVX + om->ddbound_w > xB_AVX) 
  {
//printf("vitfilter AVX doing DD calc, Dmax_AVX = %i\n", Dmax_AVX);
    /* Now we're obligated to do at least one complete DD path to be sure. */
    /* dcv has carried through from end of q loop above */
    //left-shift macro
      __m256i temp_mask_AVX = _mm256_permute2x128_si256(dcv_AVX, dcv_AVX, _MM_SHUFFLE(0,0,3,0) );
      dcv_AVX = _mm256_alignr_epi8(dcv_AVX, temp_mask_AVX,14);
 
    dcv_AVX = _mm256_or_si256(dcv_AVX, negInfv_AVX);
    tsc_AVX = om->twv_AVX + 7*Q_AVX;  /* set tsc to start of the DD's */
    for (q_AVX = 0; q_AVX < Q_AVX; q_AVX++) 
      {
        DMX_AVXf(q_AVX) = _mm256_max_epi16(dcv_AVX, DMX_AVXf(q_AVX));  
        dcv_AVX     = _mm256_adds_epi16(DMX_AVXf(q_AVX), *tsc_AVX); tsc_AVX++;
      }

    /* We may have to do up to three more passes; the check
     * is for whether crossing a segment boundary can improve
     * our score. 
     */
    do {
      //left-shift macro
      __m256i temp_mask_AVX = _mm256_permute2x128_si256(dcv_AVX, dcv_AVX, _MM_SHUFFLE(0,0,3,0) );
      dcv_AVX = _mm256_alignr_epi8(dcv_AVX, temp_mask_AVX,14);
 
      dcv_AVX = _mm256_or_si256(dcv_AVX, negInfv_AVX);
      tsc_AVX = om->twv_AVX + 7*Q_AVX;  /* set tsc to start of the DD's */
      for (q_AVX = 0; q_AVX < Q_AVX; q_AVX++) 
        {
            __m256i check1 = _mm256_cmpgt_epi16(dcv_AVX, DMX_AVXf(q_AVX));
            int cont = _mm256_movemask_epi8(check1); // Continue if any of the 16-bit field s in dcv are greater than 
            // the one we compared to
            if (cont == 0) break;
            DMX_AVXf(q_AVX) = _mm256_max_epi16(dcv_AVX, DMX_AVXf(q_AVX));  
            dcv_AVX   = _mm256_adds_epi16(DMX_AVXf(q_AVX), *tsc_AVX); tsc_AVX++;
         }     
    } while (q_AVX == Q_AVX);
  }
      else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
  {
//printf("vitfilter AVX skipping DD calc Dmax_AVX = %i\n", Dmax_AVX);
    //left-shift macro
    __m256i temp_mask_AVX = _mm256_permute2x128_si256(dcv_AVX, dcv_AVX, _MM_SHUFFLE(0,0,3,0) );
    dcv_AVX = _mm256_alignr_epi8(dcv_AVX, temp_mask_AVX,14);
   
    DMX_AVXf(0) = _mm256_or_si256(dcv_AVX, negInfv_AVX);
  }
#endif

#ifdef p7_build_AVX512 // This is somewhat ugly, but putting SSE, AVX, AVX-512 versions in the same loop makes it easier to
// check each iteration         
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
         //left-shift macro

      // left_shift dp_temp by 128 bits by shuffling and then inzerting zeroes at the low end
       __m512i temp_mask_AVX_512 = _mm512_shuffle_i32x4(mpv_AVX_512, mpv_AVX_512, 0x90);
       __m128i zero128 = _mm_setzero_si128();
       temp_mask_AVX_512 = _mm512_inserti32x4(temp_mask_AVX_512, zero128, 0);

       //now do the same merge and right-shift trick we used with AVX to create a left-shift by two bytes
      mpv_AVX_512 = _mm512_alignr_epi8(mpv_AVX_512, temp_mask_AVX_512,14);

      mpv_AVX_512 = _mm512_or_si512(mpv_AVX_512, negInfv_AVX_512);
      
      dpv_AVX_512 = DMX_AVX_512f(Q_AVX_512-1);  
         //left-shift macro
  
     // left_shift dp_temp by 128 bits by shuffling and then inzerting zeroes at the low end
       temp_mask_AVX_512 = _mm512_shuffle_i32x4(dpv_AVX_512, dpv_AVX_512, 0x90);
       temp_mask_AVX_512 = _mm512_inserti32x4(temp_mask_AVX_512, zero128, 0);

       //now do the same merge and right-shift trick we used with AVX to create a left-shift by two bytes
      dpv_AVX_512 = _mm512_alignr_epi8(dpv_AVX_512, temp_mask_AVX_512,14);
       dpv_AVX_512 = _mm512_or_si512(dpv_AVX_512, negInfv_AVX_512);
      
      ipv_AVX_512 = IMX_AVX_512f(Q_AVX_512-1);  
         //left-shift macro

     // left_shift dp_temp by 128 bits by shuffling and then inzerting zeroes at the low end
       temp_mask_AVX_512 = _mm512_shuffle_i32x4(ipv_AVX_512, ipv_AVX_512, 0x90);
       temp_mask_AVX_512 = _mm512_inserti32x4(temp_mask_AVX_512, zero128, 0);

       //now do the same merge and right-shift trick we used with AVX to create a left-shift by one byte
      ipv_AVX_512 = _mm512_alignr_epi8(ipv_AVX_512, temp_mask_AVX_512,14);
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
    //left-shift macro
//printf("vitfilter AVX-512 doing DD calc,  Dmax_AVX_512 = %i\n", Dmax_AVX_512);
       temp_mask_AVX_512 = _mm512_shuffle_i32x4(dcv_AVX_512, dcv_AVX_512, 0x90);
       temp_mask_AVX_512 = _mm512_inserti32x4(temp_mask_AVX_512, zero128, 0);

       //now do the same merge and right-shift trick we used with AVX to create a left-shift by one byte
      dcv_AVX_512 = _mm512_alignr_epi8(dcv_AVX_512, temp_mask_AVX_512,14);

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
      //left-shift macro
      temp_mask_AVX_512 = _mm512_shuffle_i32x4(dcv_AVX_512, dcv_AVX_512, 0x90);
       temp_mask_AVX_512 = _mm512_inserti32x4(temp_mask_AVX_512, zero128, 0);

       //now do the same merge and right-shift trick we used with AVX to create a left-shift by one byte
      dcv_AVX_512 = _mm512_alignr_epi8(dcv_AVX_512, temp_mask_AVX_512,14);
 
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
    //left-shift macro
//printf("vitfilter AVX-512 skipping DD calc,  Dmax_AVX_512 = %i\n", Dmax_AVX_512);
       temp_mask_AVX_512 = _mm512_shuffle_i32x4(dcv_AVX_512, dcv_AVX_512, 0x90);
       temp_mask_AVX_512 = _mm512_inserti32x4(temp_mask_AVX_512, zero128, 0);

       //now do the same merge and right-shift trick we used with AVX to create a left-shift by one byte
      dcv_AVX_512 = _mm512_alignr_epi8(dcv_AVX_512, temp_mask_AVX_512,14);
    DMX_AVX_512f(0) = _mm512_or_si512(dcv_AVX_512, negInfv_AVX_512);
  }
#endif

#ifdef p7_build_check_AVX2
    int16_t *unstriped, *unstriped_AVX;
    union { __m128i v; int16_t i[8]; } tmp_check;
    union { __m256i v; int16_t i[16]; } tmp_check_AVX;

    unstriped = malloc( sizeof(int16_t) * ((Q*8)+1));  // Yes, these allocates are slow,but this is check code that won't be
    unstriped_AVX = malloc( sizeof(int16_t) * ((Q_AVX*16)+1));  // compiled in production
    int q_temp;
    int z_temp;

    // Code here liberally borrowed from p7_filtermx_DumpVFRow
    // First, check the M's
    for (q_temp = 0; q_temp < Q; q_temp++) {
      tmp_check.v = MMXf(q_temp);
      for (z_temp = 0; z_temp < 8; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.i[z_temp];
    }
    for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
      tmp_check_AVX.v = MMX_AVXf(q_temp);
      for (z_temp = 0; z_temp < 16; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.i[z_temp];
    }
    for(q_temp = 1; q_temp <= om->M; q_temp++){
      if(unstriped[q_temp] != unstriped_AVX[q_temp]){
        printf("SSE vs AVX M mis-match in vitfilter.c at position %d, %i vs %i \n", q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
      }
 /*     else{
        printf("SSE and AVX M match in vitfilter.c at position %d, %i and %i \n", q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
      } */
    }  

    // Now, the I's
    for (q_temp = 0; q_temp < Q; q_temp++) {
      tmp_check.v = IMXf(q_temp);
      for (z_temp = 0; z_temp < 8; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.i[z_temp];
    }
    for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
      tmp_check_AVX.v = IMX_AVXf(q_temp);
      for (z_temp = 0; z_temp < 16; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.i[z_temp];
    }
    for(q_temp = 1; q_temp <= om->M; q_temp++){
      if(unstriped[q_temp] != unstriped_AVX[q_temp]){
        printf("SSE vs AVX I mis-match in vitfilter.c at position %d, %i vs %i \n", q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
      }
 /*     else{
        printf("SSE and AVX D match in vitfilter.c at position %d, %i and %i \n", q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
      } */
    }  

       // and the D's
    for (q_temp = 0; q_temp < Q; q_temp++) {
      tmp_check.v = DMXf(q_temp);
      for (z_temp = 0; z_temp < 8; z_temp++) unstriped[q_temp+Q*z_temp+1] = tmp_check.i[z_temp];
    }
    for (q_temp = 0; q_temp < Q_AVX; q_temp++) {
      tmp_check_AVX.v = DMX_AVXf(q_temp);
      for (z_temp = 0; z_temp < 16; z_temp++) unstriped_AVX[q_temp+Q_AVX*z_temp+1] = tmp_check_AVX.i[z_temp];
    }
    for(q_temp = 1; q_temp <= om->M; q_temp++){
      if(unstriped[q_temp] != unstriped_AVX[q_temp]){
        printf("SSE vs AVX D mis-match in vitfilter.c at position %d, %i vs %i \n", q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
      }
 /*     else{
        printf("SSE and AVX D match in vitfilter.c at position %d, %i and %i \n", q_temp, unstriped[q_temp], unstriped_AVX[q_temp]);
      } */
    }  
    free(unstriped);
    free(unstriped_AVX);
#endif
#ifdef p7_build_check_AVX512
    int16_t *unstriped2, *unstriped2_AVX_512;
    union { __m128i v; int16_t i[8]; } tmp_check2;
    union { __m512i v; int16_t i[32]; } tmp_check2_AVX_512;

    unstriped2 = malloc( sizeof(int16_t) * ((Q*8)+1));  // Yes, these allocates are slow,but this is check code that won't be
    unstriped2_AVX_512 = malloc( sizeof(int16_t) * ((Q_AVX_512*32)+1));  // compiled in production
    int q_temp2;
    int z_temp2;

    // Code here liberally borrowed from p7_filtermx_DumpVFRow
    // First, check the M's
    for (q_temp2 = 0; q_temp2 < Q; q_temp2++) {
      tmp_check2.v = MMXf(q_temp2);
      for (z_temp2 = 0; z_temp2 < 8; z_temp2++) unstriped2[q_temp2+Q*z_temp2+1] = tmp_check2.i[z_temp2];
    }
    for (q_temp2 = 0; q_temp2 < Q_AVX_512; q_temp2++) {
      tmp_check2_AVX_512.v = MMX_AVX_512f(q_temp2);
      for (z_temp2 = 0; z_temp2 < 32; z_temp2++) unstriped2_AVX_512[q_temp2+Q_AVX_512*z_temp2+1] = tmp_check2_AVX_512.i[z_temp2];
    }
    for(q_temp2 = 1; q_temp2 <= om->M; q_temp2++){
      if(unstriped2[q_temp2] != unstriped2_AVX_512[q_temp2]){
        printf("SSE vs AVX-512 M mis-match in vitfilter.c at position %d, %i vs %i \n", q_temp2, unstriped2[q_temp2], unstriped2_AVX_512[q_temp2]);
      }
/*      else{
        printf("SSE and AVX-512 M match in vitfilter.c at position %d, %i and %i \n", q_temp2, unstriped2[q_temp2], unstriped2_AVX_512[q_temp2]);
      } */
    }  

    // Now, the I's
    for (q_temp2 = 0; q_temp2 < Q; q_temp2++) {
      tmp_check2.v = IMXf(q_temp2);
      for (z_temp2 = 0; z_temp2 < 8; z_temp2++) unstriped2[q_temp2+Q*z_temp2+1] = tmp_check2.i[z_temp2];
    }
    for (q_temp2 = 0; q_temp2 < Q_AVX_512; q_temp2++) {
      tmp_check2_AVX_512.v = IMX_AVX_512f(q_temp2);
      for (z_temp2 = 0; z_temp2 < 32; z_temp2++) unstriped2_AVX_512[q_temp2+Q_AVX_512*z_temp2+1] = tmp_check2_AVX_512.i[z_temp2];
    }
    for(q_temp2 = 1; q_temp2 <= om->M; q_temp2++){
      if(unstriped2[q_temp2] != unstriped2_AVX_512[q_temp2]){
        printf("SSE vs AVX-512 I mis-match in vitfilter.c at position %d, %i vs %i \n", q_temp2, unstriped2[q_temp2], unstriped2_AVX_512[q_temp2]);
      }
/*    else{
        printf("SSE and AVX-512 I match in vitfilter.c at position %d, %i and %i \n", q_temp2, unstriped2[q_temp2], unstriped2_AVX_512[q_temp2]);
      } */
    }  

       // and the D's
    for (q_temp2 = 0; q_temp2 < Q; q_temp2++) {
      tmp_check2.v = DMXf(q_temp2);
      for (z_temp2 = 0; z_temp2 < 8; z_temp2++) unstriped2[q_temp2+Q*z_temp2+1] = tmp_check2.i[z_temp2];
    }
    for (q_temp2 = 0; q_temp2 < Q_AVX_512; q_temp2++) {
      tmp_check2_AVX_512.v = DMX_AVX_512f(q_temp2);
      for (z_temp2 = 0; z_temp2 < 32; z_temp2++) unstriped2_AVX_512[q_temp2+Q_AVX_512*z_temp2+1] = tmp_check2_AVX_512.i[z_temp2];
    }
    for(q_temp2 = 1; q_temp2 <= om->M; q_temp2++){
      if(unstriped2[q_temp2] != unstriped2_AVX_512[q_temp2]){
        printf("SSE vs AVX-512 D mis-match in vitfilter.c at position %d, %i vs %i \n", q_temp2, unstriped2[q_temp2], unstriped2_AVX_512[q_temp2]);
      }
 /*   else{
        printf("SSE and AVX D match in vitfilter.c at position %d, %i and %i \n", q_temp2, unstriped2[q_temp2], unstriped2_AVX_512[q_temp2]);
      } */
    }  
    free(unstriped2);
    free(unstriped2_AVX_512);
#endif
#ifdef p7_DEBUGGING
      if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, i, xE, 0, xJ, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */


#ifdef p7_build_check_AVX2
      if (xC != xC_AVX){
        printf("SSE and AVX result miss-match in viterbi filter %i vs %i \n", xC, xC_AVX);
      }
  /*    else{
        printf("AVX and SSE matched in viterbi filter\n");
      } */
#endif

#ifdef p7_build_check_AVX512
      if (xC != xC_AVX_512){
        printf("SSE and AVX-512 result miss-match in viterbi filter %i vs %i \n", xC, xC_AVX_512);
      }
/*    else{
        printf("AVX-512 and SSE matched in viterbi filter\n");
      } */
#endif


#ifdef p7_build_SSE
  /* finally C->T */
  if (xC > -32768)
    {
      *ret_sc = (float) xC + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
      /* *ret_sc += L * om->ncj_roundoff;  see J4/150 for rationale: superceded by -3.0nat approximation*/
      *ret_sc /= om->scale_w;
      *ret_sc -= 3.0; /* the NN/CC/JJ=0,-3nat approximation: see J5/36. That's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ contrib */
    }
  else  *ret_sc = -eslINFINITY;
#endif
#ifdef p7_build_AVX2
  /* finally C->T */
  if (xC_AVX > -32768)
    {
      *ret_sc = (float) xC_AVX + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
      /* *ret_sc += L * om->ncj_roundoff;  see J4/150 for rationale: superceded by -3.0nat approximation*/
      *ret_sc /= om->scale_w;
      *ret_sc -= 3.0; /* the NN/CC/JJ=0,-3nat approximation: see J5/36. That's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ contrib */
    }
  else  *ret_sc = -eslINFINITY;
#endif
#ifdef p7_build_AVX512
  /* finally C->T */
  if (xC_AVX_512 > -32768)
    {
      *ret_sc = (float) xC_AVX_512 + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
      /* *ret_sc += L * om->ncj_roundoff;  see J4/150 for rationale: superceded by -3.0nat approximation*/
      *ret_sc /= om->scale_w;
      *ret_sc -= 3.0; /* the NN/CC/JJ=0,-3nat approximation: see J5/36. That's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ contrib */
    }
  else  *ret_sc = -eslINFINITY;
#endif
  return eslOK;
}
/*---------------- end, p7_ViterbiFilter() ----------------------*/



/* Function:  p7_ViterbiFilter_longtarget()
 * Synopsis:  Finds windows within potentially long sequence blocks with Viterbi
 *            scores above threshold (vewy vewy fast, in limited precision)
 *
 * Purpose:   Calculates an approximation of the Viterbi score for regions
 *            of sequence <dsq>, using optimized profile <om>, and a pre-
 *            allocated one-row DP matrix <ox>, and captures the positions
 *            at which such regions exceed the score required to be
 *            significant in the eyes of the calling function (usually
 *            p=0.001).
 *
 *            The resulting landmarks are converted to subsequence
 *            windows by the calling function
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
 *            filtersc   - null or bias correction, required for translating a P-value threshold into a score threshold
 *            P          - p-value below which a region is captured as being above threshold
 *            windowlist - RETURN: preallocated array of hit windows (start and end of diagonal) for the above-threshold areas
 *
 * Returns:   <eslOK> on success;
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small, or if
 *            profile isn't in a local alignment mode. (Must be in local
 *            alignment mode because that's what helps us guarantee
 *            limited dynamic range.)
 *
 * Xref:      See p7_ViterbiFilter()
 */
int
p7_ViterbiFilter_longtarget(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_FILTERMX *ox,
                            float filtersc, double P, P7_HMM_WINDOWLIST *windowlist)
{
  register __m128i mpv, dpv, ipv; /* previous row values                                       */
  register __m128i sv;            /* temp storage of 1 curr row value in progress              */
  register __m128i dcv;           /* delayed storage of D(i,q+1)                               */
  register __m128i xEv;           /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;           /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i Dmaxv;         /* keeps track of maximum D cell on row                      */
  int16_t  xE, xB, xC, xJ, xN;    /* special states' scores                                    */
  int16_t  Dmax;                  /* maximum D cell score on row                               */
  int i;                          /* counter over sequence positions 1..L                      */
  int q;                          /* counter over vectors 0..nq-1                              */
  int Q        = P7_NVW(om->M);   /* segment length: # of vectors                              */
  __m128i *dp  = ox->dp;          /* using {MDI}MXf(q) macro requires initialization of <dp>   */
  __m128i *rsc;                   /* will point at om->ru[x] for residue x[i]                  */
  __m128i *tsc;                   /* will point into (and step thru) om->tu                    */
  __m128i negInfv;
  int16_t sc_thresh;
  float   invP;
  int z;
  union { __m128i v; int16_t i[8]; } tmp;
  int status;

  windowlist->count = 0;

  /* Contract checks / argument validation */
  ESL_DASSERT1( (om->mode == p7_LOCAL || om->mode == p7_UNILOCAL) ); /* VF numerics only work for local alignment */

  /* Resize the filter mx as needed */
  if (( status = p7_filtermx_GrowTo(ox, om->M))    != eslOK) ESL_EXCEPTION(status, "Reallocation of Vit filter matrix failed");

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  ox->M    = om->M;
  ox->type = p7F_VITFILTER;

/*
 *  In p7_ViterbiFilter, converting from a scaled int Viterbi score
 *  S (aka xE the score getting to state E) to a probability
 *  goes like this:
 *    vsc =  S + om->xw[p7O_E][p7O_MOVE] + om->xw[p7O_C][p7O_MOVE] - om->base_w
 *    ret_sc /= om->scale_w;
 *    vsc -= 3.0;
 *    P  = esl_gumbel_surv((vfsc - filtersc) / eslCONST_LOG2  ,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
 *  and we're computing the threshold vsc, so invert it:
 *    (vsc - filtersc) /  eslCONST_LOG2 = esl_gumbel_invsurv( P, om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA])
 *    vsc = filtersc + eslCONST_LOG2 * esl_gumbel_invsurv( P, om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA])
 *    vsc += 3.0
 *    vsc *= om->scale_w
 *    S = vsc - (float)om->xw[p7O_E][p7O_MOVE] - (float)om->xw[p7O_C][p7O_MOVE] + (float)om->base_w
 */
  invP = esl_gumbel_invsurv(P, om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
  sc_thresh =   (int) ceil ( ( (filtersc + (eslCONST_LOG2 * invP) + 3.0) * om->scale_w )
                - (float)om->xw[p7O_E][p7O_MOVE] - (float)om->xw[p7O_C][p7O_MOVE] + (float)om->base_w );

  /* -infinity is -32768 */
  negInfv = _mm_set1_epi16(-32768);
  negInfv = _mm_srli_si128(negInfv, 14);  /* negInfv = 16-byte vector, 14 0 bytes + 2-byte value=-32768, for an OR operation. */

  /* Initialization. In unsigned arithmetic, -infinity is -32768
   */
  for (q = 0; q < Q; q++)
    MMXf(q) = IMXf(q) = DMXf(q) = _mm_set1_epi16(-32768);
  xN   = om->base_w;
  xB   = xN + om->xw[p7O_N][p7O_MOVE];
  xJ   = -32768;
  xC   = -32768;
  xE   = -32768;

#if p7_DEBUGGING
  if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, 0, xE, 0, xJ, xB, xC); /* first 0 is <rowi>: do header. second 0 is xN: always 0 here. */
#endif


  for (i = 1; i <= L; i++)
  {
      rsc   = om->rwv[dsq[i]];
      tsc   = om->twv;
      dcv   = _mm_set1_epi16(-32768);      /* "-infinity" */
      xEv   = _mm_set1_epi16(-32768);
      Dmaxv = _mm_set1_epi16(-32768);
      xBv   = _mm_set1_epi16(xB);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12.
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically; replace it with -32768.
       */
      mpv = MMXf(Q-1);  mpv = _mm_slli_si128(mpv, 2);  mpv = _mm_or_si128(mpv, negInfv);
      dpv = DMXf(Q-1);  dpv = _mm_slli_si128(dpv, 2);  dpv = _mm_or_si128(dpv, negInfv);
      ipv = IMXf(Q-1);  ipv = _mm_slli_si128(ipv, 2);  ipv = _mm_or_si128(ipv, negInfv);

      for (q = 0; q < Q; q++)
      {
        /* Calculate new MMXf(i,q); don't store it yet, hold it in sv. */
        sv   =                    _mm_adds_epi16(xBv, *tsc);  tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(mpv, *tsc)); tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(ipv, *tsc)); tsc++;
        sv   = _mm_max_epi16 (sv, _mm_adds_epi16(dpv, *tsc)); tsc++;
        sv   = _mm_adds_epi16(sv, *rsc);                      rsc++;
        xEv  = _mm_max_epi16(xEv, sv);

        /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
         * {MDI}MX(q) is then the current, not the prev row
         */
        mpv = MMXf(q);
        dpv = DMXf(q);
        ipv = IMXf(q);

        /* Do the delayed stores of {MD}(i,q) now that memory is usable */
        MMXf(q) = sv;
        DMXf(q) = dcv;

        /* Calculate the next D(i,q+1) partially: M->D only;
               * delay storage, holding it in dcv
         */
        dcv   = _mm_adds_epi16(sv, *tsc);  tsc++;
        Dmaxv = _mm_max_epi16(dcv, Dmaxv);

        /* Calculate and store I(i,q) */
        sv     =                    _mm_adds_epi16(mpv, *tsc);  tsc++;
        IMXf(q)= _mm_max_epi16 (sv, _mm_adds_epi16(ipv, *tsc)); tsc++;
      }

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE = esl_sse_hmax_epi16(xEv);

      if (xE >= sc_thresh) {
        //hit score threshold. Add a window to the list, then reset scores.

        /* Unpack and unstripe, then find the position responsible for the hit */
        for (q = 0; q < Q; q++) {
          tmp.v = MMXf(q);
          for (z = 0; z < 8; z++)  { // unstripe
            if ( tmp.i[z] == xE && (q+Q*z+1) <= om->M) {
              // (q+Q*z+1) is the model position k at which the xE score is found
              p7_hmmwindow_new(windowlist, 0, i, 0, (q+Q*z+1), 1, 0.0, p7_NOCOMPLEMENT );
            }
          }
          MMXf(q) = IMXf(q) = DMXf(q) = _mm_set1_epi16(-32768); //reset score to start search for next vit window.
        }

      } else {


        xN = xN + om->xw[p7O_N][p7O_LOOP];
        xC = ESL_MAX(xC + om->xw[p7O_C][p7O_LOOP], xE + om->xw[p7O_E][p7O_MOVE]);
        xJ = ESL_MAX(xJ + om->xw[p7O_J][p7O_LOOP], xE + om->xw[p7O_E][p7O_LOOP]);
        xB = ESL_MAX(xJ + om->xw[p7O_J][p7O_MOVE], xN + om->xw[p7O_N][p7O_MOVE]);
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
        Dmax = esl_sse_hmax_epi16(Dmaxv);
        if (Dmax + om->ddbound_w > xB)
        {
          /* Now we're obligated to do at least one complete DD path to be sure. */
          /* dcv has carried through from end of q loop above */
          dcv = _mm_slli_si128(dcv, 2);
          dcv = _mm_or_si128(dcv, negInfv);
          tsc = om->twv + 7*Q;  /* set tsc to start of the DD's */
          for (q = 0; q < Q; q++)
          {
            DMXf(q) = _mm_max_epi16(dcv, DMXf(q));
            dcv     = _mm_adds_epi16(DMXf(q), *tsc); tsc++;
          }

          /* We may have to do up to three more passes; the check
           * is for whether crossing a segment boundary can improve
           * our score.
           */
          do {
            dcv = _mm_slli_si128(dcv, 2);
            dcv = _mm_or_si128(dcv, negInfv);
            tsc = om->twv + 7*Q;  /* set tsc to start of the DD's */
            for (q = 0; q < Q; q++)
            {
              if (! esl_sse_any_gt_epi16(dcv, DMXf(q))) break;
              DMXf(q) = _mm_max_epi16(dcv, DMXf(q));
              dcv     = _mm_adds_epi16(DMXf(q), *tsc);   tsc++;
            }
          } while (q == Q);
        }
        else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
        {
          dcv = _mm_slli_si128(dcv, 2);
          DMXf(0) = _mm_or_si128(dcv, negInfv);
        }
      }
#if p7_DEBUGGING
      if (ox->do_dumping) p7_filtermx_DumpVFRow(ox, i, xE, 0, xJ, xB, xC);
#endif
  } /* end loop over sequence residues 1..L */


  return eslOK;

}
/*---------------- end, p7_ViterbiFilter_longtarget() ----------------------*/



/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7VITFILTER_BENCHMARK
/* -c, -x are used for debugging, testing; see msvfilter.c for explanation 
   ./vitfilter_benchmark <hmmfile>          runs benchmark 
   ./vitfilter_benchmark -N100 -c <hmmfile> compare scores to generic impl
   ./vitfilter_benchmark -N100 -x <hmmfile> compare scores to exact emulation
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-x", "compare scores to generic implementation (debug)", 0 }, 
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-c", "equate scores to trusted implementation (debug)",  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for Viterbi filter";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_FILTERMX    *ox      = NULL;
  P7_REFMX       *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc1, sc2;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, L);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  if (esl_opt_GetBoolean(go, "-x")) p7_profile_SameAsVF(om, gm);

  ox = p7_filtermx_Create(om->M);
  gx = p7_refmx_Create(gm->M, L);

  /* Get a baseline time: how long it takes just to generate the sequences */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Run the benchmark */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_ViterbiFilter(dsq, L, om, ox, &sc1);   

      if (esl_opt_GetBoolean(go, "-c")) 
	{
	  p7_ReferenceViterbi(dsq, L, gm, gx, NULL, &sc2); 
	  printf("%.4f %.4f\n", sc1, sc2);  
	}

      if (esl_opt_GetBoolean(go, "-x"))
	{
	  p7_ReferenceViterbi(dsq, L, gm, gx, NULL, &sc2); 
	  sc2 /= om->scale_w;
	  if (om->mode == p7_UNILOCAL)   sc2 -= 2.0; /* that's ~ L \log \frac{L}{L+2}, for our NN,CC,JJ */
	  else if (om->mode == p7_LOCAL) sc2 -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
	  printf("%.4f %.4f\n", sc1, sc2);  
	}
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_filtermx_Destroy(ox);
  p7_refmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7VITFILTER_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/




/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7VITFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_viterbi.h"

/* utest_comparison()
 * 
 * Check against the reference Viterbi, after configuring  
 * a profile such that its scores will match the roundoffs in 
 * the ViterbiFilter -- p7_profile_SameAsVF().
 * 
 * Sample a random model of length <M>, and score <N> random
 * test sequences of length <L>.
 *
 * We assume that we don't accidentally generate a high-scoring random
 * sequence that overflows ViterbiFilter()'s limited range.
 * 
 */
static void
utest_comparison(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_FILTERMX *ox  = p7_filtermx_Create(M);
  P7_REFMX    *gx  = p7_refmx_Create(M, L);
  float sc1, sc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_profile_SameAsVF(om, gm);	/* round and scale the scores in <gm> the same as in <om> */

#if 0
  p7_oprofile_Dump(stdout, om);                   // dumps the optimized profile
  p7_filtermx_SetDumpMode(ox, stdout, TRUE);      // makes the fast DP algorithms dump their matrices
#endif

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_ViterbiFilter   (dsq, L, om, ox,       &sc1);
      p7_ReferenceViterbi(dsq, L, gm, gx, NULL, &sc2);

#if 0
      p7_refmx_Dump(stdout, gx);   // dumps a generic DP matrix
#endif
      
      sc2 /= om->scale_w;
      sc2 -= 3.0;

      if (fabs(sc1-sc2) > 0.001) esl_fatal("viterbi filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);
      
      p7_refmx_Reuse(gx);
      p7_filtermx_Reuse(ox);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_filtermx_Destroy(ox);
  p7_refmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7VITFILTER_TESTDRIVE*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7VITFILTER_TESTDRIVE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the SSE implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter() tests, DNA\n");
  utest_comparison(r, abc, bg, M, L, N);   
  utest_comparison(r, abc, bg, 1, L, 10);  
  utest_comparison(r, abc, bg, M, 1, 10);  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter() tests, protein\n");
  utest_comparison(r, abc, bg, M, L, N); 
  utest_comparison(r, abc, bg, 1, L, 10);
  utest_comparison(r, abc, bg, M, 1, 10);

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);

  fprintf(stderr, "#  status = ok\n");
  return eslOK;
}
#endif /*VITFILTER_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/



/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7VITFILTER_EXAMPLE
/* A minimal example.
   Also useful for debugging on small HMMs and sequences.
   ./vitfilter_example <hmmfile> <seqfile>
 */ 
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in one line awkable format",                0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in profmark format",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of Viterbi filter algorithm";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_FILTERMX    *ox      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           vfraw, nullsc, vfscore;
  double          P;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_profile_ConfigLocal(gm, hmm, bg, sq->n);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  ox = p7_filtermx_Create(gm->M);

  /* Useful to place and compile in for debugging: 
     p7_oprofile_Dump(stdout, om);                      // dumps the optimized profile
     p7_filtermx_SetDumpMode(ox, stdout, TRUE);         // makes the fast DP algorithms dump their matrices
     p7_refmx_Dump(stdout, gx);                         // dumps a generic DP matrix
  */

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_profile_SetLength      (gm, sq->n);
      p7_bg_SetLength(bg,            sq->n);

      p7_ViterbiFilter  (sq->dsq, sq->n, om, ox, &vfraw);
      p7_bg_NullOne (bg, sq->dsq, sq->n, &nullsc);
      vfscore = (vfraw - nullsc) / eslCONST_LOG2;
      P        = esl_gumbel_surv(vfscore,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);

      if (esl_opt_GetBoolean(go, "-1"))
	{
	  printf("%-30s\t%-20s\t%9.2g\t%7.2f\n", sq->name, hmm->name, P, vfscore);
	}
      else if (esl_opt_GetBoolean(go, "-P"))
	{ /* output suitable for direct use in profmark benchmark postprocessors: */
	  printf("%g\t%.2f\t%s\t%s\n", P, vfscore, sq->name, hmm->name);
	}
      else
	{
	  printf("target sequence:      %s\n",        sq->name);
	  printf("vit filter raw score: %.2f nats\n", vfraw);
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", vfscore);
	  printf("P-value:              %g\n",        P);
	}
      
      p7_filtermx_Reuse(ox);
      esl_sq_Reuse(sq);
    }

  /* cleanup */
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_filtermx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7VITFILTER_EXAMPLE*/
/*-------------------- end, example -----------------------------*/


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/

