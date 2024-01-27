/* Viterbi filter; x86 SSE version.
 * 
 * This is a SIMD vectorized, striped, interleaved, one-row, reduced
 * precision (16 bit) implementation of the Viterbi algorithm.
 * 
 * It calculates a close approximation of the Viterbi score, in
 * limited precision (signed words: 16 bits) and range. It may overflow on
 * high scoring sequences, but this indicates that the sequence is a
 * high-scoring hit worth examining more closely anyway.  It will not
 * underflow in local alignment mode.
 */
#include <h4_config.h>

#include <stdio.h>
#include <math.h>

#include "easel.h"

#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_filtermx.h"

#include "simdvec.h"

/* The non-vector includes are deliberately outside the #ifdef.  If
 * SSE is enabled at compile-time, we provide an h4_vitfilter_sse()
 * function that works; otherwise we provide one that issues a fatal
 * error.  This allows driver programs (such as vitfilter_benchmark)
 * to have --sse/--avx/--avx512 options to call
 * vitfilter_{sse,avx,avx512} directly, overriding the CPU dispatcher.
 */
#ifdef eslENABLE_SSE4
#include <x86intrin.h>
#include "esl_sse.h"

/*****************************************************************
 * 1. Viterbi filter implementation.
 *****************************************************************/

/* Function:  h4_vitfilter_sse()
 * Synopsis:  Calculates Viterbi score, vewy vewy fast, in limited precision.
 * See:       vitfilter.c::p7_ViterbiFilter()
 */
int
h4_vitfilter_sse(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc)
{
  register __m128i mpv, ipv, dpv;  // previous row values                                       
  register __m128i mcv, icv, dcv;  // current row values
  register __m128i xEv;            // E state: keeps max for Mk->E as we go                     
  register __m128i xBv;            // B state: splatted vector of B[i-1] for B->Mk calculations 
  int16_t  xE, xB, xC, xJ, xN;     // special states' scores. No L,G because we're all-local in VF                                   
  const int Q = H4_Q(hmm->M, h4_VWIDTH_SSE / sizeof(int16_t));  // segment length: # of vectors        
  const __m128i *rsc;                             // this steps thru hmm->rsc[x]
  const __m128i *tsc;                             // this steps thru hmm->tsc (all but DD's)
  const __m128i *itsc = (__m128i *) hmm->twv[Q];  // this steps thru hmm->tsc DD transitions
  const __m128i neginfmask = _mm_insert_epi16( _mm_setzero_si128(), -32768, 0);  // [ * 0 0 0 ]
  __m128i      *dp;
  int           i;                                // counter over sequence positions 1..L                      
  int           q;                                // counter over vectors 0..nq-1                              
  int           status;

  /* Contract checks */
  ESL_DASSERT1(( hmm->flags & h4_HASVECS )); 
  ESL_DASSERT1(( hmm->V == h4_VWIDTH_SSE ));
  ESL_DASSERT1(( mo->L == L ));         
  /* ViterbiFilter implicitly computes local alignment only; it ignores B->L|G parameters in <mo> */

  /* Resize the filter mx as needed */
  if (( status = h4_filtermx_Reinit(fx, hmm->M)) != eslOK) goto ERROR;
  dp = (__m128i *) fx->dp;   // enables MMXf(), IMXf(), DMXf() access macros. Must be set AFTER _Reinit(), because fx->dp may get moved!

  /* Matrix type and size must be set early, not late: debugging dump functions need this information. */
  fx->M    = hmm->M;
  fx->Vw   = h4_VWIDTH_SSE / sizeof(int16_t);

  /* Initialization. In int16_t, our -infinity is -32768  */
  for (q = 0; q < Q; q++)
    MMXf(q) = IMXf(q) = DMXf(q) = _mm_set1_epi16(-32768);
  xN   = (int16_t) h4_BASE_W;   // VF's base score offset is defined in simdvec.h
  xB   = xN + mo->xw[h4_N][h4_MOVE];
  xJ   = -32768;
  xC   = -32768;
  xE   = -32768;

#if eslDEBUGLEVEL > 0
  if (fx->dfp) h4_filtermx_DumpRow(fx, 0, xE, 0, xJ, xB, xC); /* first 0 is <rowi>: do header. second 0 is xN: always 0 here. */
#endif

  for (i = 1; i <= L; i++)
    {
      rsc   = (__m128i *) hmm->rwv[dsq[i]];
      tsc   = (__m128i *) hmm->twv[0];
      dcv   = _mm_set1_epi16(-32768);
      xEv   = _mm_set1_epi16(-32768);     
      xBv   = _mm_set1_epi16(xB);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12. */
      mpv = esl_sse_rightshift_int16(MMXf(Q-1), neginfmask);
      dpv = esl_sse_rightshift_int16(DMXf(Q-1), neginfmask);
      ipv = esl_sse_rightshift_int16(IMXf(Q-1), neginfmask);

      for (q = 0; q < Q; q++)
        {
          /* Calculate new M(i,q); don't store it yet, hold it in mcv. */
          mcv  =                     _mm_adds_epi16(xBv, *tsc);  tsc++;
          mcv  = _mm_max_epi16 (mcv, _mm_adds_epi16(mpv, *tsc)); tsc++;
          mcv  = _mm_max_epi16 (mcv, _mm_adds_epi16(ipv, *tsc)); tsc++;
          mcv  = _mm_max_epi16 (mcv, _mm_adds_epi16(dpv, *tsc)); tsc++;
          mcv  = _mm_adds_epi16(mcv, *rsc);                      rsc++;
          xEv  = _mm_max_epi16(xEv, mcv);

          /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
           * {MDI}MX(q) is about to become current, not the prev row
           */
          mpv = MMXf(q);
          dpv = DMXf(q);
          ipv = IMXf(q);

          /* Do the delayed stores of {MD}(i,q) now that memory is usable */
          MMXf(q) = mcv;
          DMXf(q) = dcv;

          /* Calculate and store I(i,q) */
          icv =                              _mm_adds_epi16(mpv, *tsc);  tsc++;
          icv =          _mm_max_epi16 (icv, _mm_adds_epi16(ipv, *tsc)); tsc++;
          icv = IMXf(q)= _mm_max_epi16 (icv, _mm_adds_epi16(dpv, *tsc)); tsc++;

          /* Calculate the next D(i,q+1) partially; delay storage, holding it in dcv */
          dcv   =                    _mm_adds_epi16(dcv,  itsc[q]);
          dcv   = _mm_max_epi16(dcv, _mm_adds_epi16(mcv, *tsc)); tsc++;
          dcv   = _mm_max_epi16(dcv, _mm_adds_epi16(icv, *tsc)); tsc++;
        }

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE   = esl_sse_hmax_epi16(xEv);
      if (xE >= 32767)                // Immediately detect overflow; calculate max representative bitscore and return
        {                             // (Saturated max is 32767 when we overflow xE. The >= is overkill.)
          *ret_sc = (float) xE + (float) mo->xw[h4_C][h4_MOVE] - h4_BASE_W;  // E->C costs 0; C->C's are in the 3NAT_APPROX below; C->T cost assessed; then take out the offset
          *ret_sc /= h4_SCALE_W;      //  ... and unscale...
          *ret_sc -= h4_3NAT_APPROX;  // ... assess the NN/CC/JJ=0,-3nat approximation: see J5/36. That's L \log_2 \frac{L}{L+3} for L>>0, for our NN,CC,JJ contrib 
          *ret_sc -= mo->nullsc;      // ... and subtract remaining null model score. Now we have a complete bit score.
          return eslERANGE;
        }            
      xC = ESL_MAX(xC,                         xE + mo->xw[h4_E][h4_MOVE]);    // Assume the 3 nat approximation... NN/CC/JJ transitions = 0. 
      xJ = ESL_MAX(xJ,                         xE + mo->xw[h4_E][h4_LOOP]);    //   xN never changes        (otherwise xN = xN + mo->xw[h4_N][h4_LOOP])                               
      xB = ESL_MAX(xJ + mo->xw[h4_J][h4_MOVE], xN + mo->xw[h4_N][h4_MOVE]);    //   xC, xJ elide transition (otherwise xC = ESL_MAX(xC + mo->xw[h4_C][h4_LOOP]...), and resp. for xJ)
      /* and now xB will carry over into next i, and xC carries over after i=L */

      /* "lazy F" loop [Farrar07]
       * Doing a complete segment before the any_gt test is slightly faster
       * than doing the any_gt test at each q, because the any_gt is expensive.
       */
      dcv  = esl_sse_rightshift_int16(dcv, neginfmask);
      while ( esl_sse_any_gt_epi16( dcv, DMXf(0)))
	{
	  for (q = 0; q < Q; q++)
	    {
	      DMXf(q) = _mm_max_epi16(dcv, DMXf(q));
	      dcv     = _mm_adds_epi16(DMXf(q), itsc[q]);
	    }
	  dcv  = esl_sse_rightshift_int16(dcv, neginfmask);
	}

#if eslDEBUGLEVEL > 0
      if (fx->dfp) h4_filtermx_DumpRow(fx, i, xE, 0, xJ, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T */
  if (xC > -32768)   // -32768 = saturated underflow, which should not happen in theory. (In theory there is no difference between theory and practice but...)
    {
      *ret_sc = (float) xC + (float) mo->xw[h4_C][h4_MOVE] - h4_BASE_W;
      *ret_sc /= h4_SCALE_W;
      *ret_sc -= h4_3NAT_APPROX; // the NN/CC/JJ=0,-3nat approximation in bits: see J5/36. That's L \log_2 \frac{L}{L+3} for L>>0, for our NN,CC,JJ contrib 
      *ret_sc -= mo->nullsc;
    }
  else  *ret_sc = -eslINFINITY;
  return eslOK;

 ERROR:
  *ret_sc = -eslINFINITY;
  return status;
}



#else // ! eslENABLE_SSE4
/* provide a callable function even when we're `./configure --disable-sse` */
int
h4_vitfilter_sse(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_FILTERMX *fx, float *ret_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(hmm); ESL_UNUSED(mo); ESL_UNUSED(fx); ESL_UNUSED(ret_sc);  // shut up, compiler, I know what I'm doing
  esl_fatal("SSE support was not enabled at compile time. Can't use h4_vitfilter_sse().");
  return eslFAIL; // NOTREACHED
}
#endif // eslENABLE_SSE4 or not




