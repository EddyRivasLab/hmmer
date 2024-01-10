/* SSV filter; x86 SSE version.
 * Originally contributed by Bjarne Knudsen (CLC Bio); modified since.
 * 
 * See ssvfilter.md for notes.
 *
 * This file is conditionally compiled when eslENABLE_SSE4 is defined.
 */
#include <h4_config.h>

#include <stdio.h>
#include <math.h>

#include "easel.h"

#include "h4_mode.h"
#include "h4_profile.h"

#include "simdvec.h"

#ifdef eslENABLE_SSE4
#include <x86intrin.h>
#include "esl_sse.h"

/* Note that some ifdefs below has to be changed if these values are
   changed. These values are chosen based on some simple speed
   tests. Apparently, two registers are generally used for something
   else, leaving 14 registers on 64 bit versions and 6 registers on 32
   bit versions. */
#ifdef __x86_64__ /* 64 bit version */
#define  MAX_BANDS 14
#else
#define  MAX_BANDS 6
#endif


#define STEP_SINGLE(sv)                   \
  sv   = _mm_adds_epi8(sv, *rsc); rsc++;  \
  xEv  = _mm_max_epi8(xEv, sv);     

#define LENGTH_CHECK(label) \
  if (i >= L) goto label;

#define NO_CHECK(label)

#define STEP_BANDS_1() \
  STEP_SINGLE(sv00)  

#define STEP_BANDS_2() \
  STEP_BANDS_1()       \
  STEP_SINGLE(sv01)

#define STEP_BANDS_3() \
  STEP_BANDS_2()       \
  STEP_SINGLE(sv02)

#define STEP_BANDS_4() \
  STEP_BANDS_3()       \
  STEP_SINGLE(sv03)

#define STEP_BANDS_5() \
  STEP_BANDS_4()       \
  STEP_SINGLE(sv04)

#define STEP_BANDS_6() \
  STEP_BANDS_5()       \
  STEP_SINGLE(sv05)

#define STEP_BANDS_7() \
  STEP_BANDS_6()       \
  STEP_SINGLE(sv06)

#define STEP_BANDS_8() \
  STEP_BANDS_7()       \
  STEP_SINGLE(sv07)

#define STEP_BANDS_9() \
  STEP_BANDS_8()       \
  STEP_SINGLE(sv08)

#define STEP_BANDS_10() \
  STEP_BANDS_9()        \
  STEP_SINGLE(sv09)

#define STEP_BANDS_11() \
  STEP_BANDS_10()       \
  STEP_SINGLE(sv10)

#define STEP_BANDS_12() \
  STEP_BANDS_11()       \
  STEP_SINGLE(sv11)

#define STEP_BANDS_13() \
  STEP_BANDS_12()       \
  STEP_SINGLE(sv12)

#define STEP_BANDS_14() \
  STEP_BANDS_13()       \
  STEP_SINGLE(sv13)

#define STEP_BANDS_15() \
  STEP_BANDS_14()       \
  STEP_SINGLE(sv14)

#define STEP_BANDS_16() \
  STEP_BANDS_15()       \
  STEP_SINGLE(sv15)

#define STEP_BANDS_17() \
  STEP_BANDS_16()       \
  STEP_SINGLE(sv16)

#define STEP_BANDS_18() \
  STEP_BANDS_17()       \
  STEP_SINGLE(sv17)

#define CONVERT_STEP(step, length_check, label, sv, pos) \
  length_check(label)                                    \
  rsc = (__m128i *) hmm->rbv[dsq[i]] + pos;              \
  step()                                                 \
  sv = esl_sse_rightshift_int8(sv, neginfmask);          \
  i++;

#define CONVERT_1(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv00, Q - 1) 

#define CONVERT_2(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv01, Q - 2)  \
  CONVERT_1(step, LENGTH_CHECK, label)

#define CONVERT_3(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv02, Q - 3)  \
  CONVERT_2(step, LENGTH_CHECK, label)

#define CONVERT_4(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv03, Q - 4)  \
  CONVERT_3(step, LENGTH_CHECK, label)

#define CONVERT_5(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv04, Q - 5)  \
  CONVERT_4(step, LENGTH_CHECK, label)

#define CONVERT_6(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv05, Q - 6)  \
  CONVERT_5(step, LENGTH_CHECK, label)

#define CONVERT_7(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv06, Q - 7)  \
  CONVERT_6(step, LENGTH_CHECK, label)

#define CONVERT_8(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv07, Q - 8)  \
  CONVERT_7(step, LENGTH_CHECK, label)

#define CONVERT_9(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv08, Q - 9)  \
  CONVERT_8(step, LENGTH_CHECK, label)

#define CONVERT_10(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv09, Q - 10) \
  CONVERT_9(step, LENGTH_CHECK, label)

#define CONVERT_11(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv10, Q - 11) \
  CONVERT_10(step, LENGTH_CHECK, label)

#define CONVERT_12(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv11, Q - 12) \
  CONVERT_11(step, LENGTH_CHECK, label)

#define CONVERT_13(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv12, Q - 13) \
  CONVERT_12(step, LENGTH_CHECK, label)

#define CONVERT_14(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv13, Q - 14) \
  CONVERT_13(step, LENGTH_CHECK, label)

#define CONVERT_15(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv14, Q - 15) \
  CONVERT_14(step, LENGTH_CHECK, label)

#define CONVERT_16(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv15, Q - 16) \
  CONVERT_15(step, LENGTH_CHECK, label)

#define CONVERT_17(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv16, Q - 17) \
  CONVERT_16(step, LENGTH_CHECK, label)

#define CONVERT_18(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv17, Q - 18) \
  CONVERT_17(step, LENGTH_CHECK, label)

#define RESET_1()                 \
  register __m128i sv00 = beginv;

#define RESET_2()                 \
  RESET_1()                       \
  register __m128i sv01 = beginv;

#define RESET_3()                 \
  RESET_2()                       \
  register __m128i sv02 = beginv;

#define RESET_4()                 \
  RESET_3()                       \
  register __m128i sv03 = beginv;

#define RESET_5()                 \
  RESET_4()                       \
  register __m128i sv04 = beginv;

#define RESET_6()                 \
  RESET_5()                       \
  register __m128i sv05 = beginv;

#define RESET_7()                 \
  RESET_6()                       \
  register __m128i sv06 = beginv;

#define RESET_8()                 \
  RESET_7()                       \
  register __m128i sv07 = beginv;

#define RESET_9()                 \
  RESET_8()                       \
  register __m128i sv08 = beginv;

#define RESET_10()                \
  RESET_9()                       \
  register __m128i sv09 = beginv;

#define RESET_11()                \
  RESET_10()                      \
  register __m128i sv10 = beginv;

#define RESET_12()                \
  RESET_11()                      \
  register __m128i sv11 = beginv;

#define RESET_13()                \
  RESET_12()                      \
  register __m128i sv12 = beginv;

#define RESET_14()                \
  RESET_13()                      \
  register __m128i sv13 = beginv;

#define RESET_15()                \
  RESET_14()                      \
  register __m128i sv14 = beginv;

#define RESET_16()                \
  RESET_15()                      \
  register __m128i sv15 = beginv;

#define RESET_17()                \
  RESET_16()                      \
  register __m128i sv16 = beginv;

#define RESET_18()                \
  RESET_17()                      \
  register __m128i sv17 = beginv;


/* Original Knudsen CALC() was simpler. (Look back in git for comparison.)
 * Nick Carter applied an additional optimization, using 4:1 loop
 * unrolling to reduce loop overhead, a 13% runtime improvement in the
 * AVX version.
 */
#define CALC(reset, step, convert, width)       \
  __m128i  neginfmask = _mm_insert_epi8( _mm_setzero_si128(), -128, 0); \
  int      Q          = H4_Q(hmm->M, h4_VWIDTH_SSE);                    \
  int      w          = width;                  \
  int      i,i2;                                \
  __m128i *rsc;                                 \
  int      num_iters;                           \
                                                \
  dsq++;                                        \
  reset()                                       \
  if (L <= Q-q-w)  num_iters = L;               \
  else             num_iters = Q-q-w;           \
  i = 0;                                        \
  while (num_iters >= 4) {                      \
    rsc = (__m128i *) hmm->rbv[dsq[i]] + i + q; \
    step()                                      \
    i++;                                        \
    rsc = (__m128i *) hmm->rbv[dsq[i]] + i + q; \
    step()                                      \
    i++;                                        \
    rsc = (__m128i *) hmm->rbv[dsq[i]] + i + q; \
    step()                                      \
    i++;                                        \
    rsc = (__m128i *) hmm->rbv[dsq[i]] + i + q; \
    step()                                      \
    i++;                                        \
    num_iters -= 4;                             \
  }                                             \
  while (num_iters > 0) {                       \
    rsc = (__m128i *) hmm->rbv[dsq[i]] + i + q; \
    step()                                      \
    i++;                                        \
    num_iters--;                                \
  }                                             \
  i = Q - q - w;                                \
  convert(step, LENGTH_CHECK, done1)            \
done1:                                          \
  for (i2 = Q - q; i2 < L - Q; i2 += Q)         \
    {                                           \
      i = 0;                                    \
      num_iters = Q - w;                        \
      while (num_iters >= 4) {                        \
        rsc = (__m128i *) hmm->rbv[dsq[i2 + i]] + i;  \
        step()                                        \
        i++;                                          \
        rsc = (__m128i *) hmm->rbv[dsq[i2 + i]] + i;  \
        step()                                        \
        i++;                                          \
        rsc = (__m128i *) hmm->rbv[dsq[i2 + i]] + i;  \
        step()                                        \
          i++;                                        \
        rsc = (__m128i *) hmm->rbv[dsq[i2 + i]] + i;  \
        step()                                        \
        i++;                                          \
        num_iters-= 4;                                \
      }                                               \
      while (num_iters > 0) {                         \
        rsc = (__m128i *) hmm->rbv[dsq[i2 + i]] + i;  \
        step()                                        \
        i++;                                          \
        num_iters--;                                  \
      }                                               \
      i += i2;                                   \
      convert(step, NO_CHECK, )                  \
    }                                            \
  if ((L-i2) < (Q-w))   num_iters = L-i2;        \
  else                  num_iters = Q - w;       \
  i = 0;                                         \
  while (num_iters >= 4) {                       \
    rsc = (__m128i *) hmm->rbv[dsq[i2 + i]] + i; \
    step()                                       \
    i+= 1;                                       \
    rsc = (__m128i *) hmm->rbv[dsq[i2 + i]] + i; \
    step()                                       \
    i+= 1;                                       \
    rsc = (__m128i *) hmm->rbv[dsq[i2 + i]] + i; \
    step()                                       \
    i+= 1;                                       \
    rsc = (__m128i *) hmm->rbv[dsq[i2 + i]] + i; \
    step()                                       \
    i+= 1;                                       \
    num_iters -= 4;                              \
  }                                              \
  while (num_iters > 0) {                        \
    rsc = (__m128i *) hmm->rbv[dsq[i2 + i]] + i; \
    step()                                       \
    i+= 1;                                       \
    num_iters--;                                 \
  }                                              \
  i+=i2;                                         \
  convert(step, LENGTH_CHECK, done2)             \
done2:                                           \
  return xEv;

static __m128i calc_band_1(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_1, STEP_BANDS_1, CONVERT_1, 1) }
static __m128i calc_band_2(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_2, STEP_BANDS_2, CONVERT_2, 2) }
static __m128i calc_band_3(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_3, STEP_BANDS_3, CONVERT_3, 3) }
static __m128i calc_band_4(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_4, STEP_BANDS_4, CONVERT_4, 4) }
static __m128i calc_band_5(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_5, STEP_BANDS_5, CONVERT_5, 5) }
static __m128i calc_band_6(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_6, STEP_BANDS_6, CONVERT_6, 6) }

#if MAX_BANDS > 6 /* Only include needed functions to limit object file size */
static __m128i calc_band_7 (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_7,  STEP_BANDS_7,  CONVERT_7,  7)  }
static __m128i calc_band_8 (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_8,  STEP_BANDS_8,  CONVERT_8,  8)  }
static __m128i calc_band_9 (const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_9,  STEP_BANDS_9,  CONVERT_9,  9)  }
static __m128i calc_band_10(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_10, STEP_BANDS_10, CONVERT_10, 10) }
static __m128i calc_band_11(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_11, STEP_BANDS_11, CONVERT_11, 11) }
static __m128i calc_band_12(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_12, STEP_BANDS_12, CONVERT_12, 12) }
static __m128i calc_band_13(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_13, STEP_BANDS_13, CONVERT_13, 13) }
static __m128i calc_band_14(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_14, STEP_BANDS_14, CONVERT_14, 14) }
#endif /* MAX_BANDS > 6 */

#if MAX_BANDS > 14
static __m128i calc_band_15(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_15, STEP_BANDS_15, CONVERT_15, 15) }
static __m128i calc_band_16(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_16, STEP_BANDS_16, CONVERT_16, 16) }
static __m128i calc_band_17(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_17, STEP_BANDS_17, CONVERT_17, 17) }
static __m128i calc_band_18(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, int q, __m128i beginv, register __m128i xEv) { CALC(RESET_18, STEP_BANDS_18, CONVERT_18, 18) }
#endif /* MAX_BANDS > 14 */


/* Function:  h4_ssvfilter_sse()
 * Synopsis:  SSV filter; x86 SSE version.
 */
int
h4_ssvfilter_sse(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, float *ret_sc)
{
  int     Q       = H4_Q(hmm->M, h4_VWIDTH_SSE); // segment length: # of vectors                  
  __m128i beginv  = _mm_set1_epi8(-128);         // initialization of each set of diagonals
  __m128i xEv     = beginv; 	                 // E state: keeps max for Mk->E as we go    
  int     last_q  = 0;                           // for saving the last q value to find band width           
  int8_t  xE;                                    // DP max, before conversion to nats
  int     q;			                 // counter over vectors 0..nq-1  
  int     bands;                                 // number of bands (rounds) to use   
  int     i;                                     // counter for bands                     
  /* function pointers for the various number of vectors to use */
  __m128i (*fs[MAX_BANDS + 1]) (const ESL_DSQ *, int, const H4_PROFILE *, int, register __m128i, __m128i)
    = {NULL, calc_band_1,  calc_band_2,  calc_band_3,  calc_band_4,  calc_band_5,  calc_band_6
#if MAX_BANDS > 6
           , calc_band_7,  calc_band_8,  calc_band_9,  calc_band_10, calc_band_11, calc_band_12, calc_band_13, calc_band_14
#endif
#if MAX_BANDS > 14
          , calc_band_15, calc_band_16, calc_band_17, calc_band_18
#endif
  };

  /* Use the highest number of bands but no more than MAX_BANDS */
  bands = (Q + MAX_BANDS - 1) / MAX_BANDS;
  for (i = 0; i < bands; i++) 
    {
      q      = (Q * (i + 1)) / bands;
      xEv    = fs[q-last_q](dsq, L, hmm, last_q, beginv, xEv);
      last_q = q;
    }
  xE = esl_sse_hmax_epi8(xEv);  

  /* SSV will overflow normally on good, high-scoring seqs (xE sticks at 127). 
   * In that case, we set sc to max representable raw bitscore, and return ERANGE.
   */
                       // vvv   Add +128 back onto the diagonal score. DP calculated it from -128 baseline.
  *ret_sc = ((float) xE + 128.) / h4_SCALE_B + hmm->tauBM - h4_2NAT_APPROX;   // 2.0 is the tauNN/tauCC "2 nat approximation". tauBM is B->Mk entry score (unscaled!)
  *ret_sc += 2.0 * log2f(2.0 / (float) (L + 2));                              // tauNB, tauCT moves... also unscaled floats.  
  *ret_sc -= mo->nullsc;
  return (xE == 127 ? eslERANGE : eslOK);
}


#else // ! eslENABLE_SSE4
/* provide callable functions even when we're `./configure --disable-sse` */
int
h4_ssvfilter_sse(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, float *ret_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(hmm); ESL_UNUSED(ret_sc);  
  esl_fatal("SSE support was not enabled at compile time. Can't use h4_ssvfilter_sse().");
  return eslFAIL; // NOTREACHED
}
#endif // eslENABLE_SSE4 or not








