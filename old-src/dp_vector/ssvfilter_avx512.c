/* SSV filter, x86 AVX-512 vector implementation.
 * Adapted from SSE version by Nick Carter
 * 
 * See ssvfilter.md for notes.
 * See ssvfilter_sse.c for a version with more explanation & commentary.
 *
 * This file is conditionally compiled when eslENABLE_AVX512 is defined.
 */
#include <p7_config.h>

#include <math.h>

#include "easel.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/ssvfilter.h"

#ifdef eslENABLE_AVX512
#include <x86intrin.h>
#include "esl_avx512.h"

#define MAX_BANDS 30

#define STEP_SINGLE(sv)                     \
  sv   = _mm512_adds_epi8(sv, *rsc); rsc++; \
  xEv  = _mm512_max_epi8(xEv, sv);

#define LENGTH_CHECK(label) \
  if (i >= L) goto label;

#define NO_CHECK(label)

#define STEP_BANDS_1()  \
  STEP_SINGLE(sv00) 

#define STEP_BANDS_2()  \
  STEP_BANDS_1()        \
  STEP_SINGLE(sv01)

#define STEP_BANDS_3()  \
  STEP_BANDS_2()        \
  STEP_SINGLE(sv02)

#define STEP_BANDS_4()  \
  STEP_BANDS_3()        \
  STEP_SINGLE(sv03)

#define STEP_BANDS_5()  \
  STEP_BANDS_4()        \
  STEP_SINGLE(sv04)

#define STEP_BANDS_6()  \
  STEP_BANDS_5()        \
  STEP_SINGLE(sv05)

#define STEP_BANDS_7()  \
  STEP_BANDS_6()        \
  STEP_SINGLE(sv06)

#define STEP_BANDS_8()  \
  STEP_BANDS_7()        \
  STEP_SINGLE(sv07)

#define STEP_BANDS_9()  \
  STEP_BANDS_8()        \
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

#define STEP_BANDS_19() \
  STEP_BANDS_18()       \
  STEP_SINGLE(sv18)

#define STEP_BANDS_20() \
  STEP_BANDS_19()       \
  STEP_SINGLE(sv19)

#define STEP_BANDS_21() \
  STEP_BANDS_20()       \
  STEP_SINGLE(sv20)

#define STEP_BANDS_22() \
  STEP_BANDS_21()       \
  STEP_SINGLE(sv21)

#define STEP_BANDS_23() \
  STEP_BANDS_22()       \
  STEP_SINGLE(sv22)

#define STEP_BANDS_24() \
  STEP_BANDS_23()       \
  STEP_SINGLE(sv23)

#define STEP_BANDS_25() \
  STEP_BANDS_24()       \
  STEP_SINGLE(sv24)

#define STEP_BANDS_26() \
  STEP_BANDS_25()       \
  STEP_SINGLE(sv25)

#define STEP_BANDS_27() \
  STEP_BANDS_26()       \
  STEP_SINGLE(sv26)

#define STEP_BANDS_28() \
  STEP_BANDS_27()       \
  STEP_SINGLE(sv27)

#define STEP_BANDS_29() \
  STEP_BANDS_28()       \
  STEP_SINGLE(sv28)

#define STEP_BANDS_30() \
  STEP_BANDS_29()       \
  STEP_SINGLE(sv29)

#define STEP_BANDS_31() \
  STEP_BANDS_30()       \
  STEP_SINGLE(sv30)

#define STEP_BANDS_32() \
  STEP_BANDS_31()       \
  STEP_SINGLE(sv31)

#define CONVERT_STEP(step, length_check, label, sv, pos) \
  length_check(label)                                    \
  rsc = (__m512i *) om->rbv[dsq[i]] + pos;               \
  step()                                                 \
  sv = esl_avx512_rightshift_int8(sv, neginfmask);       \
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
#define CONVERT_19(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv18, Q - 19) \
  CONVERT_18(step, LENGTH_CHECK, label)

#define CONVERT_20(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv19, Q - 20) \
  CONVERT_19(step, LENGTH_CHECK, label)

#define CONVERT_21(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv20, Q - 21) \
  CONVERT_20(step, LENGTH_CHECK, label)

#define CONVERT_22(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv21, Q - 22) \
  CONVERT_21(step, LENGTH_CHECK, label)

#define CONVERT_23(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv22, Q - 23) \
  CONVERT_22(step, LENGTH_CHECK, label)

#define CONVERT_24(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv23, Q - 24) \
  CONVERT_23(step, LENGTH_CHECK, label)

#define CONVERT_25(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv24, Q - 25) \
  CONVERT_24(step, LENGTH_CHECK, label)

#define CONVERT_26(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv25, Q - 26) \
  CONVERT_25(step, LENGTH_CHECK, label)

#define CONVERT_27(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv26, Q - 27) \
  CONVERT_26(step, LENGTH_CHECK, label)

#define CONVERT_28(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv27, Q - 28) \
  CONVERT_27(step, LENGTH_CHECK, label)

#define CONVERT_29(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv28, Q - 29) \
  CONVERT_28(step, LENGTH_CHECK, label)

#define CONVERT_30(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv29, Q - 30) \
  CONVERT_29(step, LENGTH_CHECK, label)

#define CONVERT_31(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv30, Q - 31) \
  CONVERT_30(step, LENGTH_CHECK, label)

#define CONVERT_32(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv31, Q - 32) \
  CONVERT_31(step, LENGTH_CHECK, label)

#define RESET_1()                  \
  register __m512i sv00 = beginv;

#define RESET_2()                  \
  RESET_1()                        \
  register __m512i sv01 = beginv;

#define RESET_3()                  \
  RESET_2()                        \
  register __m512i sv02 = beginv;

#define RESET_4()                  \
  RESET_3()                        \
  register __m512i sv03 = beginv;

#define RESET_5()                  \
  RESET_4()                        \
  register __m512i sv04 = beginv;

#define RESET_6()                  \
   RESET_5()                       \
  register __m512i sv05 = beginv;

#define RESET_7()                  \
  RESET_6()                        \
  register __m512i sv06 = beginv;

#define RESET_8()                  \
  RESET_7()                        \
  register __m512i sv07 = beginv;

#define RESET_9()                  \
  RESET_8()                        \
  register __m512i sv08 = beginv;

#define RESET_10()                 \
  RESET_9()                        \
  register __m512i sv09 = beginv;

#define RESET_11()                 \
  RESET_10()                       \
  register __m512i sv10 = beginv;

#define RESET_12()                 \
  RESET_11()                       \
  register __m512i sv11 = beginv;

#define RESET_13()                 \
  RESET_12()                       \
  register __m512i sv12 = beginv;

#define RESET_14()                 \
  RESET_13()                       \
  register __m512i sv13 = beginv;

#define RESET_15()                 \
  RESET_14()                       \
  register __m512i sv14 = beginv;

#define RESET_16()                 \
  RESET_15()                       \
  register __m512i sv15 = beginv;

#define RESET_17()                 \
  RESET_16()                       \
  register __m512i sv16 = beginv;

#define RESET_18()                 \
  RESET_17()                       \
  register __m512i sv17 = beginv;

#define RESET_19()                \
  RESET_18()                      \
  register __m512i sv18 = beginv;

#define RESET_20()                \
  RESET_19()                      \
  register __m512i sv19 = beginv;

#define RESET_21()                \
  RESET_20()                      \
  register __m512i sv20 = beginv;

#define RESET_22()                \
  RESET_21()                      \
  register __m512i sv21 = beginv;

#define RESET_23()                \
  RESET_22()                      \
  register __m512i sv22 = beginv;

#define RESET_24()                \
  RESET_23()                      \
  register __m512i sv23 = beginv;

#define RESET_25()                \
  RESET_24()                      \
  register __m512i sv24 = beginv;

#define RESET_26()                \
  RESET_25()                      \
  register __m512i sv25 = beginv;

#define RESET_27()                \
  RESET_26()                      \
  register __m512i sv26 = beginv;

#define RESET_28()                \
  RESET_27()                      \
  register __m512i sv27 = beginv;

#define RESET_29()                \
  RESET_28()                      \
  register __m512i sv28 = beginv;

#define RESET_30()                \
  RESET_29()                      \
  register __m512i sv29 = beginv;

#define RESET_31()                \
  RESET_30()                      \
  register __m512i sv30 = beginv;

#define RESET_32()                \
  RESET_31()                      \
  register __m512i sv31 = beginv;

#define CALC(reset, step, convert, width)      \
  __m512i  neginfmask = _mm512_mask_blend_epi8(0x1, _mm512_setzero_si512(), beginv); \
  int      Q          = P7_Q(om->M, p7_VWIDTH_AVX512);                               \
  int      w          = width;                 \
  int      i, i2;                              \
  __m512i *rsc;                                \
  int      num_iters;                          \
                                               \
  dsq++;                                       \
  reset()                                      \
  if (L <= Q-q-w) num_iters = L;               \
  else  	  num_iters = Q-q-w;           \
  i = 0;                                       \
  while (num_iters >=4) {                      \
    rsc = (__m512i *) om->rbv[dsq[i]] + i + q; \
    step()                                     \
    i++;                                       \
    rsc = (__m512i *) om->rbv[dsq[i]] + i + q; \
    step()                                     \
    i++;                                       \
    rsc = (__m512i *) om->rbv[dsq[i]] + i + q; \
    step()                                     \
    i++;                                       \
    rsc = (__m512i *) om->rbv[dsq[i]] + i + q; \
    step()                                     \
    i++;                                       \
    num_iters -= 4;                            \
  }                                            \
  while (num_iters >0) {                       \
    rsc = (__m512i *) om->rbv[dsq[i]] + i + q; \
    step()                                     \
    i++;                                       \
    num_iters--;                               \
  }                                            \
  i = Q - q - w;                               \
  convert(step, LENGTH_CHECK, done1)           \
done1:                                         \
  for (i2 = Q - q; i2 < L - Q; i2 += Q)        \
   {                                           \
     i = 0;                                    \
     num_iters = Q - w;                        \
     while (num_iters >= 4) {                  \
       rsc = (__m512i *) om->rbv[dsq[i2 + i]] + i; \
       step()                                      \
       i++;                                        \
       rsc = (__m512i *) om->rbv[dsq[i2 + i]] + i; \
       step()                                      \
       i++;                                        \
       rsc = (__m512i *) om->rbv[dsq[i2 + i]] + i; \
       step()                                      \
       i++;                                        \
       rsc = (__m512i *) om->rbv[dsq[i2 + i]] + i; \
       step()                                      \
       i++;                                        \
       num_iters-= 4;                              \
     }                                             \
     while (num_iters > 0) {                       \
       rsc = (__m512i *) om->rbv[dsq[i2 + i]] + i; \
       step()                                      \
       i++;                                        \
       num_iters--;                                \
     }                                             \
     i += i2;                                      \
     convert(step, NO_CHECK, )                     \
   }                                               \
  if (L-i2 < Q-w)  num_iters = L-i2;            \
  else 		   num_iters = Q-w;             \
  i = 0;                                        \
  while (num_iters >= 4) {                      \
    rsc = (__m512i *) om->rbv[dsq[i2 + i]] + i; \
    step()                                      \
    i+= 1;                                      \
    rsc = (__m512i *) om->rbv[dsq[i2 + i]] + i; \
    step()                                      \
    i+= 1;                                      \
    rsc = (__m512i *) om->rbv[dsq[i2 + i]] + i; \
    step()                                      \
    i+= 1;                                      \
    rsc = (__m512i *) om->rbv[dsq[i2 + i]] + i; \
    step()                                      \
    i+= 1;                                      \
    num_iters -= 4;                             \
  }                                             \
  while (num_iters > 0) {                       \
    rsc = (__m512i *) om->rbv[dsq[i2 + i]] + i; \
    step()                                      \
    i+= 1;                                      \
    num_iters--;                                \
  }                                             \
  i+=i2;                                        \
  convert(step, LENGTH_CHECK, done2)            \
done2:                                          \
  return xEv;

static __m512i
calc_band_1(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_1, STEP_BANDS_1, CONVERT_1, 1)
}

static __m512i
calc_band_2(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_2, STEP_BANDS_2, CONVERT_2, 2)
}

static __m512i
calc_band_3(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_3, STEP_BANDS_3, CONVERT_3, 3)
}

static __m512i
calc_band_4(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_4, STEP_BANDS_4, CONVERT_4, 4)
}

static __m512i
calc_band_5(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_5, STEP_BANDS_5, CONVERT_5, 5)
}

static __m512i
calc_band_6(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_6, STEP_BANDS_6, CONVERT_6, 6)
}

#if MAX_BANDS > 6 /* Only include needed functions to limit object file size */
static __m512i
calc_band_7(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_7, STEP_BANDS_7, CONVERT_7, 7)
}

static __m512i
calc_band_8(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_8, STEP_BANDS_8, CONVERT_8, 8)
}

static __m512i
calc_band_9(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_9, STEP_BANDS_9, CONVERT_9, 9)
}

static __m512i
calc_band_10(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_10, STEP_BANDS_10, CONVERT_10, 10)
}

static __m512i
calc_band_11(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_11, STEP_BANDS_11, CONVERT_11, 11)
}

static __m512i
calc_band_12(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_12, STEP_BANDS_12, CONVERT_12, 12)
}

static __m512i
calc_band_13(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_13, STEP_BANDS_13, CONVERT_13, 13)
}

static __m512i
calc_band_14(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_14, STEP_BANDS_14, CONVERT_14, 14)
}
#endif /* MAX_BANDS > 6 */
#if MAX_BANDS > 14
static __m512i
calc_band_15(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_15, STEP_BANDS_15, CONVERT_15, 15)
}

static __m512i
calc_band_16(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_16, STEP_BANDS_16, CONVERT_16, 16)
}

static __m512i
calc_band_17(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_17, STEP_BANDS_17, CONVERT_17, 17)
}

static __m512i
calc_band_18(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_18, STEP_BANDS_18, CONVERT_18, 18)
}
#endif /* MAX_BANDS > 14 */

#if MAX_BANDS > 18
static __m512i
calc_band_19(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_19, STEP_BANDS_19, CONVERT_19, 19)
}

static __m512i
calc_band_20(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_20, STEP_BANDS_20, CONVERT_20, 20)
}

static __m512i
calc_band_21(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_21, STEP_BANDS_21, CONVERT_21, 21)
}

static __m512i
calc_band_22(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_22, STEP_BANDS_22, CONVERT_22, 22)
}

static __m512i
calc_band_23(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_23, STEP_BANDS_23, CONVERT_23, 23)
}

static __m512i
calc_band_24(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_24, STEP_BANDS_24, CONVERT_24, 24)
}

static __m512i
calc_band_25(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_25, STEP_BANDS_25, CONVERT_25, 25)
}

static __m512i
calc_band_26(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_26, STEP_BANDS_26, CONVERT_26, 26)
}

static __m512i
calc_band_27(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_27, STEP_BANDS_27, CONVERT_27, 27)
}

static __m512i
calc_band_28(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_28, STEP_BANDS_28, CONVERT_28, 28)
}

static __m512i
calc_band_29(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_29, STEP_BANDS_29, CONVERT_29, 29)
}

static __m512i
calc_band_30(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_30, STEP_BANDS_30, CONVERT_30, 30)
}
#endif
#if MAX_BANDS > 30
static __m512i
calc_band_31(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_31, STEP_BANDS_31, CONVERT_18, 31)
}

static __m512i
calc_band_32(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, register __m512i xEv)
{
  CALC(RESET_32, STEP_BANDS_32, CONVERT_32, 32)
} 

#endif

/* Function:  p7_SSVFilter_avx512()
 * Synopsis:  SSV filter; x86 AVX-512 version. 
 * See:       ssvfilter.c::p7_SSVFilter()
 */
int
p7_SSVFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  int     Q      = P7_Q(om->M, p7_VWIDTH_AVX512); // segment length: # of vectors
  __m512i beginv = _mm512_set1_epi8(-128);        // initialization of each set of diagonals
  __m512i xEv    = beginv;                        // keeps max for Mk->E as we go
  int     last_q = 0;                             // for saving the last q value to find band width 
  int8_t  xE;                                     // DP max, before conversion to float, nats    
  int     q;			                  // counter over vectors 0..nq-1
  int     bands;                                  // the number of bands (rounds) to use
  int     i;                                      // counter for bands                              
  __m512i (*fs[MAX_BANDS + 1]) (const ESL_DSQ *, int, const P7_OPROFILE *, int, register __m512i, __m512i)
    = {NULL, calc_band_1,  calc_band_2,  calc_band_3,  calc_band_4,  calc_band_5,  calc_band_6
#if MAX_BANDS > 6
           , calc_band_7,  calc_band_8,  calc_band_9,  calc_band_10, calc_band_11, calc_band_12, calc_band_13, calc_band_14
#endif
#if MAX_BANDS > 14
          , calc_band_15, calc_band_16, calc_band_17, calc_band_18
#endif
#if MAX_BANDS > 18
          , calc_band_19, calc_band_20, calc_band_21, calc_band_22, calc_band_23, calc_band_24, calc_band_25, calc_band_26
       , calc_band_27, calc_band_28, calc_band_29, calc_band_30
#endif
#if MAX_BANDS > 30
       , calc_band_31, calc_band_32
#endif      
  };

  bands = (Q + MAX_BANDS - 1) / MAX_BANDS;
  for (i = 0; i < bands; i++) 
    {
      q   = (Q * (i + 1)) / bands;
      xEv = fs[q-last_q](dsq, L, om, last_q, beginv, xEv);
      last_q = q;
    }
  xE = esl_avx512_hmax_epi8(xEv);

  if (xE == 127)  // Overflow on high scoring sequences is expected.
    {             // Pass filter, but score is unknown.
      *ret_sc = eslINFINITY;
      return eslERANGE;
    }
  else if (xE > -128)
    {                         //v  Add +128 back onto the diagonal score. DP calculated it from -128 baseline.
      *ret_sc = ((float) xE + 128.) / om->scale_b + om->tauBM - 2.0;   // 2.0 is the tauNN/tauCC "2 nat approximation"
      *ret_sc += 2.0 * logf(2.0 / (float) (L + 2));                    // tauNB, tauCT
      return eslOK;
    }
  else 
    {
      *ret_sc = -eslINFINITY;
      return eslOK;
    }
}

#else // ! eslENABLE_AVX512
/* provide a callable function even when we're `./configure --disable-avx512` */
int
p7_SSVFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  ESL_UNUSED(dsq); ESL_UNUSED(L); ESL_UNUSED(om); ESL_UNUSED(ret_sc);  // shut up, compiler, I know what I'm doing
  esl_fatal("AVX512 support was not enabled at compile time. Can't use p7_SSVFilter_avx512().");
  return eslFAIL; // NOTREACHED
}
#if defined p7SSVFILTER_AVX512_TESTDRIVE || p7SSVFILTER_AVX512_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_AVX512 or not

