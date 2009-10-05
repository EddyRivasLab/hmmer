#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_sse.h"

#define  MAX_BANDS 14
#ifdef __INTEL_COMPILER
#define  OPT_BANDS 12  /* icc */
#else
#define  OPT_BANDS 9   /* e.g. gcc */
#endif


#define STEP_SINGLE(sv)                         \
  sv   = _mm_max_epu8(sv, beginv);              \
  /*  sv   = _mm_adds_epu8(sv, biasv);   */     \
  sv   = _mm_sub_epi8(sv, *rsc); rsc++;        \
  xEv  = _mm_max_epu8(xEv, sv);


#define LENGTH_CHECK(label)                     \
  if (i == L) goto label;


#define NO_CHECK(label)


#define STEP_BANDS_1()                          \
  STEP_SINGLE(sv00)

#define STEP_BANDS_2()                          \
  STEP_BANDS_1()                                \
  STEP_SINGLE(sv01)

#define STEP_BANDS_3()                          \
  STEP_BANDS_2()                                \
  STEP_SINGLE(sv02)

#define STEP_BANDS_4()                          \
  STEP_BANDS_3()                                \
  STEP_SINGLE(sv03)

#define STEP_BANDS_5()                          \
  STEP_BANDS_4()                                \
  STEP_SINGLE(sv04)

#define STEP_BANDS_6()                          \
  STEP_BANDS_5()                                \
  STEP_SINGLE(sv05)

#define STEP_BANDS_7()                          \
  STEP_BANDS_6()                                \
  STEP_SINGLE(sv06)

#define STEP_BANDS_8()                          \
  STEP_BANDS_7()                                \
  STEP_SINGLE(sv07)

#define STEP_BANDS_9()                          \
  STEP_BANDS_8()                                \
  STEP_SINGLE(sv08)

#define STEP_BANDS_10()                         \
  STEP_BANDS_9()                                \
  STEP_SINGLE(sv09)

#define STEP_BANDS_11()                         \
  STEP_BANDS_10()                               \
  STEP_SINGLE(sv10)

#define STEP_BANDS_12()                         \
  STEP_BANDS_11()                               \
  STEP_SINGLE(sv11)

#define STEP_BANDS_13()                         \
  STEP_BANDS_12()                               \
  STEP_SINGLE(sv12)

#define STEP_BANDS_14()                         \
  STEP_BANDS_13()                               \
  STEP_SINGLE(sv13)


#define CONVERT_STEP(step, LENGTH_CHECK, label, sv, pos)        \
  LENGTH_CHECK(label)                                           \
  rsc = om->sb[dsq[i]] + pos;                                   \
  step();                                                       \
  sv = _mm_slli_si128(sv, 1);                                   \
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


#define RESET_1()                               \
  register __m128i sv00 = beginv;

#define RESET_2()                               \
  RESET_1()                                     \
  register __m128i sv01 = beginv;

#define RESET_3()                               \
  RESET_2()                                     \
  register __m128i sv02 = beginv;

#define RESET_4()                               \
  RESET_3()                                     \
  register __m128i sv03 = beginv;

#define RESET_5()                               \
  RESET_4()                                     \
  register __m128i sv04 = beginv;

#define RESET_6()                               \
  RESET_5()                                     \
  register __m128i sv05 = beginv;

#define RESET_7()                               \
  RESET_6()                                     \
  register __m128i sv06 = beginv;

#define RESET_8()                               \
  RESET_7()                                     \
  register __m128i sv07 = beginv;

#define RESET_9()                               \
  RESET_8()                                     \
  register __m128i sv08 = beginv;

#define RESET_10()                              \
  RESET_9()                                     \
  register __m128i sv09 = beginv;

#define RESET_11()                              \
  RESET_10()                                    \
  register __m128i sv10 = beginv;

#define RESET_12()                              \
  RESET_11()                                    \
  register __m128i sv11 = beginv;

#define RESET_13()                              \
  RESET_12()                                    \
  register __m128i sv12 = beginv;

#define RESET_14()                              \
  RESET_13()                                    \
  register __m128i sv13 = beginv;


#define CALC(reset, step, convert, width)       \
  int i;                                        \
  int i2;                                       \
  int Q        = p7O_NQB(om->M);                \
  __m128i *rsc;                                 \
                                                \
  int w = width;                                \
                                                \
  dsq++;                                        \
                                                \
  reset()                                       \
                                                \
  for (i = 0; i < L && i < Q - q - w; i++)      \
    {                                           \
      rsc = om->sb[dsq[i]] + i + q;             \
      step();                                   \
    }                                           \
                                                \
  i = Q - q - w;                                \
  convert(step, LENGTH_CHECK, done1);           \
done1:                                          \
                                                \
 for (i2 = Q - q; i2 < L - Q; i2 += Q)          \
   {                                            \
     for (i = 0; i < Q - w; i++)                \
       {                                        \
         rsc = om->sb[dsq[i2 + i]] + i;         \
         step();                                \
       }                                        \
                                                \
     i += i2;                                   \
     convert(step, NO_CHECK, );                 \
   }                                            \
                                                \
 for (i = 0; i2 + i < L && i < Q - w; i++)      \
   {                                            \
     rsc = om->sb[dsq[i2 + i]] + i;             \
     step();                                    \
   }                                            \
                                                \
 i+=i2;                                         \
 convert(step, LENGTH_CHECK, done2);            \
done2:                                          \
                                                \
 return xEv;


__m128i
calc_band_1(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_1, STEP_BANDS_1, CONVERT_1, 1)
}

__m128i
calc_band_2(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_2, STEP_BANDS_2, CONVERT_2, 2)
}

__m128i
calc_band_3(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_3, STEP_BANDS_3, CONVERT_3, 3)
}

__m128i
calc_band_4(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_4, STEP_BANDS_4, CONVERT_4, 4)
}

__m128i
calc_band_5(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_5, STEP_BANDS_5, CONVERT_5, 5)
}

__m128i
calc_band_6(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_6, STEP_BANDS_6, CONVERT_6, 6)
}

__m128i
calc_band_7(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_7, STEP_BANDS_7, CONVERT_7, 7)
}

__m128i
calc_band_8(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_8, STEP_BANDS_8, CONVERT_8, 8)
}

__m128i
calc_band_9(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_9, STEP_BANDS_9, CONVERT_9, 9)
}

__m128i
calc_band_10(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_10, STEP_BANDS_10, CONVERT_10, 10)
}

__m128i
calc_band_11(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_11, STEP_BANDS_11, CONVERT_11, 11)
}

__m128i
calc_band_12(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_12, STEP_BANDS_12, CONVERT_12, 12)
}

__m128i
calc_band_13(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_13, STEP_BANDS_13, CONVERT_13, 13)
}

__m128i
calc_band_14(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv)
{
  CALC(RESET_14, STEP_BANDS_14, CONVERT_14, 14)
}


int  get_max_band_width() {
  return  MAX_BANDS;
}


uint8_t
get_xE(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om)
{
  __m128i xEv;		           /* E state: keeps max for Mk->E as we go                     */
  __m128i biasv;                   /* emission bias in a vector                                 */
  __m128i beginv;                  /* begin scores                                              */

  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQB(om->M);   /* segment length: # of vectors                              */

  __m128i (*fs[15]) (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, register __m128i biasv, register __m128i beginv, __m128i xEv) = {NULL, calc_band_1, calc_band_2, calc_band_3, calc_band_4, calc_band_5, calc_band_6, calc_band_7, calc_band_8, calc_band_9, calc_band_10, calc_band_11, calc_band_12, calc_band_13, calc_band_14};

  int w;

  int div = 1;

  int last_q = 0;
  int i;

  w = get_max_band_width();

  biasv = _mm_set1_epi8((int8_t) om->bias_b);
  beginv = _mm_set1_epi8((int8_t) om->base_b - om->tjb_b - om->tbm_b);

  xEv = _mm_setzero_si128();      

  for (;;)  {
    if ((Q + div - 1) / div <= MAX_BANDS && Q - OPT_BANDS * div < OPT_BANDS / 2) {
      break;
    }
    div++;
  }

  for (i = 0; i < div; i++) {
    q = (Q * (i + 1)) / div;

    xEv = fs[q-last_q](dsq, L, om, last_q, biasv, beginv, xEv);

    last_q = q;
  }

  return esl_sse_hmax_epu8(xEv);
}


int
p7_MSVFilter_fast(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  uint8_t  xE;
  uint8_t  xJ;

  int Q        = p7O_NQB(om->M);   /* segment length: # of vectors                              */

  if (L > 1000000) {
    /* the scoring scheme is only guaranteed to work for sequences up to 1,000,000 residues */
    return eslENORESULT;
  }

  xE = get_xE(dsq, L, om);

  if (xE >= 255 - om->bias_b)
    {
      *ret_sc = eslINFINITY;
      return eslERANGE;
    }

  xJ = xE - om->tec_b;

  if (om->base_b < xJ)  return eslENORESULT; /* The J state could have been used, so doubt about score */

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  *ret_sc = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */

  return eslOK;
}
