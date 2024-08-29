/* The SSV filter implementation; SSE version.
 * 
 * Contents:
 *   1. p7_SSVFilter() implementation (AVX Instructions)
 *
 * Nick Carter, ported from the implementation by 
 * Bjarne Knudsen, CLC Bio
 */

// see ssvfiltter.c for description/documentation

#include <p7_config.h>

#include <math.h>

#include "easel.h"
#include "esl_gumbel.h"
#ifdef eslENABLE_AVX512
#include "esl_avx512.h"
#endif
#include "hmmer.h"
#include "impl_avx.h"
#ifdef eslENABLE_AVX512
#include <x86intrin.h>
#endif
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

#ifdef eslENABLE_AVX512
#define STEP_SINGLE(sv)                         \
  sv   = _mm512_subs_epi8(sv, *rsc); rsc++;        \
  xEv  = _mm512_max_epu8(xEv, sv);


#define LENGTH_CHECK(label)                     \
  if (i >= L) goto label;


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

#define STEP_BANDS_15()                         \
  STEP_BANDS_14()                               \
  STEP_SINGLE(sv14)

#define STEP_BANDS_16()                         \
  STEP_BANDS_15()                               \
  STEP_SINGLE(sv15)

#define STEP_BANDS_17()                         \
  STEP_BANDS_16()                               \
  STEP_SINGLE(sv16)

#define STEP_BANDS_18()                         \
  STEP_BANDS_17()                               \
  STEP_SINGLE(sv17)


#define CONVERT_STEP(step, length_check, label, sv, pos)        \
  length_check(label)                                           \
  rsc = om->sbv_avx512[dsq[i]] + pos;                                   \
  step()                                                        \
  sv = esl_avx512_rightshift_int8(sv, low_byte_128); \
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


#define RESET_1()                               \
  register __m512i sv00 = beginv;

#define RESET_2()                               \
  RESET_1()                                     \
  register __m512i sv01 = beginv;

#define RESET_3()                               \
  RESET_2()                                     \
  register __m512i sv02 = beginv;

#define RESET_4()                               \
  RESET_3()                                     \
  register __m512i sv03 = beginv;

#define RESET_5()                               \
  RESET_4()                                     \
  register __m512i sv04 = beginv;

#define RESET_6()                               \
  RESET_5()                                     \
  register __m512i sv05 = beginv;

#define RESET_7()                               \
  RESET_6()                                     \
  register __m512i sv06 = beginv;

#define RESET_8()                               \
  RESET_7()                                     \
  register __m512i sv07 = beginv;

#define RESET_9()                               \
  RESET_8()                                     \
  register __m512i sv08 = beginv;

#define RESET_10()                              \
  RESET_9()                                     \
  register __m512i sv09 = beginv;

#define RESET_11()                              \
  RESET_10()                                    \
  register __m512i sv10 = beginv;

#define RESET_12()                              \
  RESET_11()                                    \
  register __m512i sv11 = beginv;

#define RESET_13()                              \
  RESET_12()                                    \
  register __m512i sv12 = beginv;

#define RESET_14()                              \
  RESET_13()                                    \
  register __m512i sv13 = beginv;

#define RESET_15()                              \
  RESET_14()                                    \
  register __m512i sv14 = beginv;

#define RESET_16()                              \
  RESET_15()                                    \
  register __m512i sv15 = beginv;

#define RESET_17()                              \
  RESET_16()                                    \
  register __m512i sv16 = beginv;

#define RESET_18()                              \
  RESET_17()                                    \
  register __m512i sv17 = beginv;


#define CALC(reset, step, convert, width)       \
  int i;                                        \
  int i2;                                       \
  int Q        = p7O_NQB_AVX512(om->M);                \
  __m512i *rsc;                                 \
  int w = width;                                \
                                                \
  dsq++;                                        \
                                                \
  reset()                                       \
                                                \
  for (i = 0; i < L && i < Q - q - w; i++)      \
    {                                           \
      rsc = om->sbv_avx512[dsq[i]] + i + q;            \
      step()                                    \
    }                                           \
                                                \
  i = Q - q - w;                                \
  convert(step, LENGTH_CHECK, done1)            \
done1:                                          \
                                                \
 for (i2 = Q - q; i2 < L - Q; i2 += Q)          \
   {                                            \
     for (i = 0; i < Q - w; i++)                \
       {                                        \
         rsc = om->sbv_avx512[dsq[i2 + i]] + i;        \
         step()                                 \
       }                                        \
                                                \
     i += i2;                                   \
     convert(step, NO_CHECK, )                  \
   }                                            \
                                                \
 for (i = 0; i2 + i < L && i < Q - w; i++)      \
   {                                            \
     rsc = om->sbv_avx512[dsq[i2 + i]] + i;            \
     step()                                     \
   }                                            \
                                                \
 i+=i2;                                         \
 convert(step, LENGTH_CHECK, done2)             \
done2:                                          \
                                                \
 return xEv;


#define CALC_UNROLLED(reset, step, convert, width)       \
  int i;                                        \
  int i2;                                       \
  int Q        = p7O_NQB_AVX512(om->M);		\
  __m512i *rsc;                                 \
                                                \
  int w = width;                                \
  dsq++;                                        \
                                                \
  reset()                                       \
  int num_iters; \
  if (L <= Q- q -w){ \
        num_iters = L;  \
  }  \
  else{  \
        num_iters = Q -q -w; \
  } \
  i = 0; \
  while(num_iters >=4){ \
        rsc = om->sbv_avx512[dsq[i]] + i + q;            \
    step()    \
    i++; \
        rsc = om->sbv_avx512[dsq[i]] + i + q;            \
    step()    \
    i++; \
    rsc = om->sbv_avx512[dsq[i]] + i + q;            \
    step()    \
    i++; \
    rsc = om->sbv_avx512[dsq[i]] + i + q;            \
    step()    \
    i++; \
    num_iters -= 4; \
  } \
  while(num_iters >0){ \
          rsc = om->sbv_avx512[dsq[i]] + i + q;            \
      step()                                 \
      i++; \
      num_iters--; \
  }  \
  i = Q - q - w;                                \
  convert(step, LENGTH_CHECK, done1)            \
done1:                                          \
 for (i2 = Q - q; i2 < L - Q; i2 += Q)          \
   {                                            \
        i = 0; \
        num_iters = Q - w; \
        while (num_iters >= 4){  \
       rsc = om->sbv_avx512[dsq[i2 + i]] + i;        \
       step() \
       i++; \
       rsc = om->sbv_avx512[dsq[i2 + i]] + i;        \
       step() \
       i++; \
       rsc = om->sbv_avx512[dsq[i2 + i]] + i;        \
       step() \
       i++; \
       rsc = om->sbv_avx512[dsq[i2 + i]] + i;        \
       step() \
       i++; \
       num_iters-= 4; \
        }\
        while(num_iters > 0){ \
	  rsc = om->sbv_avx512[dsq[i2 + i]] + i;	\
         step()    \
         i++; \
         num_iters--;      \
        } \
                                               \
     i += i2;                                   \
     convert(step, NO_CHECK, ) \
   }                                            \
if((L - i2) < (Q-w)){			   \
        num_iters = L -i2; \
        } \
 else{  \
        num_iters = Q - w; \
 } \
 i = 0; \
 while (num_iters >= 4){ \
     rsc = om->sbv_avx512[dsq[i2 + i]] + i;            \
     step()                                     \
     i+= 1;  \
     rsc = om->sbv_avx512[dsq[i2 + i]] + i;            \
     step()                                     \
     i+= 1;  \
     rsc = om->sbv_avx512[dsq[i2 + i]] + i;            \
     step()                                     \
     i+= 1;  \
     rsc = om->sbv_avx512[dsq[i2 + i]] + i;            \
     step()                                     \
     i+= 1;  \
     num_iters -= 4; \
   }                                            \
   while(num_iters > 0) {  \
         rsc = om->sbv_avx512[dsq[i2 + i]] + i;            \
     step()                                     \
     i+= 1;  \
     num_iters--; \
   } \
 i+=i2;                                         \
 convert(step, LENGTH_CHECK, done2)             \
done2:                                          \
        return xEv;

static __m512i
calc_band_1(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_1, STEP_BANDS_1, CONVERT_1, 1)
}

static __m512i
calc_band_2(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_2, STEP_BANDS_2, CONVERT_2, 2)
}

static __m512i
calc_band_3(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_3, STEP_BANDS_3, CONVERT_3, 3)
}

static __m512i
calc_band_4(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_4, STEP_BANDS_4, CONVERT_4, 4)
}

static __m512i
calc_band_5(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_5, STEP_BANDS_5, CONVERT_5, 5)
}

static __m512i
calc_band_6(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_6, STEP_BANDS_6, CONVERT_6, 6)
}

#if MAX_BANDS > 6 /* Only include needed functions to limit object file size */
static __m512i
calc_band_7(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_7, STEP_BANDS_7, CONVERT_7, 7)
}

static __m512i
calc_band_8(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_8, STEP_BANDS_8, CONVERT_8, 8)
}

static __m512i
calc_band_9(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_9, STEP_BANDS_9, CONVERT_9, 9)
}

static __m512i
calc_band_10(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_10, STEP_BANDS_10, CONVERT_10, 10)
}

static __m512i
calc_band_11(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_11, STEP_BANDS_11, CONVERT_11, 11)
}

static __m512i
calc_band_12(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_12, STEP_BANDS_12, CONVERT_12, 12)
}

static __m512i
calc_band_13(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_13, STEP_BANDS_13, CONVERT_13, 13)
}

static __m512i
calc_band_14(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_14, STEP_BANDS_14, CONVERT_14, 14)
}
#endif /* MAX_BANDS > 6 */
#if MAX_BANDS > 14
static __m512i
calc_band_15(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_15, STEP_BANDS_15, CONVERT_15, 15)
}

static __m512i
calc_band_16(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC(RESET_16, STEP_BANDS_16, CONVERT_16, 16)
}

static __m512i
calc_band_17(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_17, STEP_BANDS_17, CONVERT_17, 17)
}

static __m512i
calc_band_18(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_18, STEP_BANDS_18, CONVERT_18, 18)
}
#endif /* MAX_BANDS > 14 */

static __m512i
calc_band_unrolled_1(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_1, STEP_BANDS_1, CONVERT_1, 1)
}

static __m512i
calc_band_unrolled_2(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_2, STEP_BANDS_2, CONVERT_2, 2)
}

static __m512i
calc_band_unrolled_3(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_3, STEP_BANDS_3, CONVERT_3, 3)
}

static __m512i
calc_band_unrolled_4(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_4, STEP_BANDS_4, CONVERT_4, 4)
}

static __m512i
calc_band_unrolled_5(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_5, STEP_BANDS_5, CONVERT_5, 5)
}

static __m512i
calc_band_unrolled_6(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_6, STEP_BANDS_6, CONVERT_6, 6)
}

#if MAX_BANDS > 6 /* Only include needed functions to limit object file size */
static __m512i
calc_band_unrolled_7(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_7, STEP_BANDS_7, CONVERT_7, 7)
}

static __m512i
calc_band_unrolled_8(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_8, STEP_BANDS_8, CONVERT_8, 8)
}

static __m512i
calc_band_unrolled_9(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_9, STEP_BANDS_9, CONVERT_9, 9)
}

static __m512i
calc_band_unrolled_10(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_10, STEP_BANDS_10, CONVERT_10, 10)
}

static __m512i
calc_band_unrolled_11(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_11, STEP_BANDS_11, CONVERT_11, 11)
}

static __m512i
calc_band_unrolled_12(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_12, STEP_BANDS_12, CONVERT_12, 12)
}

static __m512i
calc_band_unrolled_13(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_13, STEP_BANDS_13, CONVERT_13, 13)
}

static __m512i
calc_band_unrolled_14(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_14, STEP_BANDS_14, CONVERT_14, 14)
}
#endif /* MAX_BANDS > 6 */
#if MAX_BANDS > 14
static __m512i
calc_band_unrolled_15(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_15, STEP_BANDS_15, CONVERT_15, 15)
}

static __m512i
calc_band_unrolled_16(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_16, STEP_BANDS_16, CONVERT_16, 16)
}

static __m512i
calc_band_unrolled_17(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_17, STEP_BANDS_17, CONVERT_17, 17)
}

static __m512i
calc_band_unrolled_18(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m512i beginv, __m512i low_byte_128, register __m512i xEv)
{
  CALC_UNROLLED(RESET_18, STEP_BANDS_18, CONVERT_18, 18)
}
#endif /* MAX_BANDS > 14 */
 

/*****************************************************************
 * 2. p7_SSVFilter() implementation
 *****************************************************************/

static uint8_t
get_xE(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om)
{
  __m512i xEv;		           /* E state: keeps max for Mk->E as we go                     */
  __m512i beginv;                  /* begin scores                                              */

  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQB_AVX512(om->M);   /* segment length: # of vectors                              */

  int bands;                       /* the number of bands (rounds) to use                       */

  int last_q = 0;                  /* for saving the last q value to find band width            */
  int i;                           /* counter for bands                                         */

  /* function pointers for the various number of vectors to use */
  __m512i (*fs[MAX_BANDS + 1]) (const ESL_DSQ *, int, const P7_OPROFILE *, int, __m512i, __m512i, register __m512i)
    = {NULL
       , calc_band_1,  calc_band_2,  calc_band_3,  calc_band_4,  calc_band_5,  calc_band_6
#if MAX_BANDS > 6
       , calc_band_7,  calc_band_8,  calc_band_9,  calc_band_10, calc_band_11, calc_band_12, calc_band_13, calc_band_14
#endif
#if MAX_BANDS > 14
       , calc_band_15, calc_band_16, calc_band_17, calc_band_18
#endif
  };

  beginv =  _mm512_set1_epi8(128);
  xEv    =  beginv;
  __m512i low_byte_128; 
  __mmask64 low_byte_mask = 0x0000000000000001; 
  low_byte_128 = _mm512_setzero_si512();
  low_byte_128 = _mm512_mask_blend_epi8(low_byte_mask, low_byte_128, beginv);
  /* Use the highest number of bands but no more than MAX_BANDS */
  bands = (Q + MAX_BANDS - 1) / MAX_BANDS;

  for (i = 0; i < bands; i++) {
    q = (Q * (i + 1)) / bands;

    xEv = fs[q-last_q](dsq, L, om, last_q, beginv, low_byte_128, xEv);

    last_q = q;
  }

  return esl_avx512_hmax_epi8(xEv);
}

static uint8_t
get_xE_unrolled(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om)
{
  __m512i xEv;		           /* E state: keeps max for Mk->E as we go                     */
  __m512i beginv;                  /* begin scores                                              */

  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQB_AVX512(om->M);   /* segment length: # of vectors                              */

  int bands;                       /* the number of bands (rounds) to use                       */

  int last_q = 0;                  /* for saving the last q value to find band width            */
  int i;                           /* counter for bands                                         */

  /* function pointers for the various number of vectors to use */
  __m512i (*fs[MAX_BANDS + 1]) (const ESL_DSQ *, int, const P7_OPROFILE *, int, __m512i, __m512i, register __m512i)
    = {NULL
       , calc_band_unrolled_1,  calc_band_unrolled_2,  calc_band_unrolled_3,  calc_band_unrolled_4,  calc_band_unrolled_5,  calc_band_unrolled_6
#if MAX_BANDS > 6
       , calc_band_unrolled_7,  calc_band_unrolled_8,  calc_band_unrolled_9,  calc_band_unrolled_10, calc_band_unrolled_11, calc_band_unrolled_12, calc_band_unrolled_13, calc_band_unrolled_14
#endif
#if MAX_BANDS > 14
       , calc_band_unrolled_15, calc_band_unrolled_16, calc_band_unrolled_17, calc_band_unrolled_18
#endif
  };

  beginv =  _mm512_set1_epi8(128);
  xEv    =  beginv;
  __m512i low_byte_128; 
  __mmask64 low_byte_mask = 0x0000000000000001; 
  low_byte_128 = _mm512_setzero_si512();
  low_byte_128 = _mm512_mask_blend_epi8(low_byte_mask, low_byte_128, beginv);
  /* Use the highest number of bands but no more than MAX_BANDS */
  bands = (Q + MAX_BANDS - 1) / MAX_BANDS;

  for (i = 0; i < bands; i++) {
    q = (Q * (i + 1)) / bands;

    xEv = fs[q-last_q](dsq, L, om, last_q, beginv, low_byte_128, xEv);

    last_q = q;
  }

  return esl_avx512_hmax_epi8(xEv);
}


int
p7_SSVFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  /* Use 16 bit values to avoid overflow due to moved baseline */
  uint16_t  xE;
  uint16_t  xJ;

  if (om->tjb_b + om->tbm_b + om->tec_b + om->bias_b >= 127) {
    /* the optimizations are not guaranteed to work under these
       conditions (see comments at start of file) */
    return eslENORESULT;
  }

  xE = get_xE(dsq, L, om);

  if (xE >= 255 - om->bias_b)
    {
      /* We have an overflow. */
      *ret_sc = eslINFINITY;
      if (om->base_b - om->tjb_b - om->tbm_b < 128) 
        {
          /* The original MSV filter may not overflow, so we are not sure our result is correct */
          return eslENORESULT;
        }

      /* We know that the overflow will also occur in the original MSV filter */
      return eslERANGE;
    }

  xE += om->base_b - om->tjb_b - om->tbm_b;
  xE -= 128;

  if (xE >= 255 - om->bias_b)
    {
      /* We know that the result will overflow in the original MSV filter */
      *ret_sc = eslINFINITY;
      return eslERANGE;
    }

  xJ = xE - om->tec_b;

  if (xJ > om->base_b)  return eslENORESULT; /* The J state could have been used, so doubt about score */

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  *ret_sc = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */

  return eslOK;
}


int
p7_SSVFilter_avx512_unrolled(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  /* Use 16 bit values to avoid overflow due to moved baseline */
  uint16_t  xE;
  uint16_t  xJ;

  if (om->tjb_b + om->tbm_b + om->tec_b + om->bias_b >= 127) {
    /* the optimizations are not guaranteed to work under these
       conditions (see comments at start of file) */
    return eslENORESULT;
  }

  xE = get_xE_unrolled(dsq, L, om);

  if (xE >= 255 - om->bias_b)
    {
      /* We have an overflow. */
      *ret_sc = eslINFINITY;
      if (om->base_b - om->tjb_b - om->tbm_b < 128) 
        {
          /* The original MSV filter may not overflow, so we are not sure our result is correct */
          return eslENORESULT;
        }

      /* We know that the overflow will also occur in the original MSV filter */
      return eslERANGE;
    }

  xE += om->base_b - om->tjb_b - om->tbm_b;
  xE -= 128;

  if (xE >= 255 - om->bias_b)
    {
      /* We know that the result will overflow in the original MSV filter */
      *ret_sc = eslINFINITY;
      return eslERANGE;
    }

  xJ = xE - om->tec_b;

  if (xJ > om->base_b)  return eslENORESULT; /* The J state could have been used, so doubt about score */

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  *ret_sc = ((float) (xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */

  return eslOK;
}


// Test code
int p7_SSVFilter_test_all(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc){

  float res_sse, res_avx, res_avx512;
  int res, res2, res3;

  res = p7_SSVFilter_sse(dsq, L, om, &res_sse);
  res2 = p7_SSVFilter_avx(dsq, L, om, &res_avx);

  if(res != res2){
    printf("Error: SSV calls returned different results: %d %d\n", res, res2);
  }

  // ret_sc is undefined if the result is eslENORESULT
  if(res != eslENORESULT && (esl_FCompare(res_sse, res_avx, .01, .01) != eslOK)){
    printf("Error: SSV miss-match.  SSE = %f, AVX = %f\n", res_sse, res_avx);
  }
  res3 = p7_SSVFilter_avx512(dsq, L, om, &res_avx512);
  if(res != res3){
    printf("Error: miss-match in SSV return codes.  SSE= %d, AVX512 = %d\n", res, res3);
  }
  if(res != eslENORESULT && esl_FCompare(res_sse, res_avx512, .01, .01) != eslOK){
    printf("Error: SSV miss-match.  SSE = %f, AVX512 = %f\n", res_sse, res_avx);
  }
  *ret_sc = res_sse;

  return res;
}


int
p7_SSVFilter_longtarget_avx512(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *ssvdata,
                        P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist)
{

  register __m512i mpv;            /* previous row values                                       */
  register __m512i xEv;		   /* E state: keeps max for Mk->E for a single iteration       */
  register __m512i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m512i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m512i biasv;	   /* emission bias in a vector                                 */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQB_AVX512(om->M);   /* segment length: # of vectors                              */
  __m512i *dp  = ox->dpb_avx512[0];	   /* we're going to use dp[0][0..q..Q-1], not {MDI}MX(q) macros*/
  __m512i *rsc;			   /* will point at om->rbv[x] for residue x[i]                 */
  __m512i tjbmv;                   /* vector for J->B move cost + B->M move costs               */
  __m512i basev;                   /* offset for scores                                         */
  __m512i ceilingv;                /* saturated simd value used to test for overflow           */
  __m512i tempv;                   /* work vector        */                                       
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

  union { __m512i v; uint8_t b[64]; } u;

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
  __m512i sc_threshv;
  uint8_t sc_thresh;
  float invP = esl_gumbel_invsurv(P, om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);


  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ16_avx512)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M   = om->M;


  p7_bg_SetLength(bg, om->max_length);
  p7_oprofile_ReconfigMSVLength(om, om->max_length);
  p7_bg_NullOne  (bg, dsq, om->max_length, &nullsc);

  sc_thresh = (int) ceil( ( ( nullsc  + (invP * eslCONST_LOG2) + 3.0 )  * om->scale_b ) + om->base_b +  om->tec_b  + om->tjb_b );
  sc_threshv = _mm512_set1_epi8((int8_t) 255 - sc_thresh);

  /* Initialization. In offset unsigned  arithmetic, -infinity is 0, and 0 is om->base.
   */
  biasv = _mm512_set1_epi8((int8_t) om->bias_b); /* yes, you can set1() an unsigned char vector this way */
  ceilingv = _mm512_set1_epi8(-1);
  for (q = 0; q < Q; q++) dp[q] = _mm512_setzero_si512();

  basev = _mm512_set1_epi8((int8_t) om->base_b);
  tjbmv = _mm512_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);

  xBv = _mm512_subs_epu8(basev, tjbmv);

  for (i = 1; i <= L; i++) 
    {
      rsc = om->rbv_avx512[dsq[i]];
      xEv = _mm512_setzero_si512();

      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12.
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically, which is our -infinity.
       */
      mpv =esl_avx512_rightshift_int8(dp[Q-1], _mm512_setzero_si512());
      
      for (q = 0; q < Q; q++) {
        /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
        sv   = _mm512_max_epu8(mpv, xBv);
        sv   = _mm512_adds_epu8(sv, biasv);
        sv   = _mm512_subs_epu8(sv, *rsc);   rsc++;
        xEv  = _mm512_max_epu8(xEv, sv);
        
        mpv   = dp[q];   	  /* Load {MDI}(i-1,q) into mpv */
        dp[q] = sv;       	  /* Do delayed store of M(i,q) now that memory is usable */
      }

      /* test if the pthresh significance threshold has been reached;
       * note: don't use _mm_cmpgt_epi8, because it's a signed comparison, which won't work on uint8s */
      tempv = _mm512_adds_epu8(xEv, sc_threshv);
      __mmask64 compare_mask = _mm512_cmpeq_epi8_mask(tempv, ceilingv);


      if (compare_mask != 0) {  //hit pthresh, so add position to list and reset values
        //figure out which model state hit threshold
        end = -1;
        rem_sc = -1;
        for (q = 0; q < Q; q++) {  /// Unpack and unstripe, so we can find the state that exceeded pthresh
          u.v = dp[q];
          for (k = 0; k < 64; k++) { // unstripe
            //(q+Q*k+1) is the model position k at which the xE score is found
            if (u.b[k] >= sc_thresh && u.b[k] > rem_sc && (q+Q*k+1) <= om->M) {
              end = (q+Q*k+1);
              rem_sc = u.b[k];
            }
          }
          dp[q] = _mm512_set1_epi8(0); // while we're here ... this will cause values to get reset to xB in next dp iteration
        }

        //recover the diagonal that hit threshold
        start = end;                    // model position
        target_end = target_start = i;  // target position
        sc = rem_sc;
        while (rem_sc > om->base_b - om->tjb_b - om->tbm_b) {
          rem_sc -= om->bias_b -  ssvdata->ssv_scores[start*om->abc->Kp + dsq[target_start]];
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
          sc += om->bias_b -  ssvdata->ssv_scores[k*om->abc->Kp + dsq[n]];
          
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
        //k    +=  (max_end - target_end);
        target_end = max_end;

        ret_sc = ((float) (max_sc - om->tjb_b) - (float) om->base_b);
        ret_sc /= om->scale_b;
        ret_sc -= 3.0; // that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ

        p7_hmmwindow_new(  windowlist,
                           0,                  // sequence_id; used in the FM-based filter, but not here
                           target_start,       // position in the target at which the diagonal starts
                           0,                  // position in the target fm_index at which diagonal starts;  not used here, just in FM-based filter
                           end,                // position in the model at which the diagonal ends
                           end-start+1 ,       // length of diagonal
                           ret_sc,             // score of diagonal
                           p7_NOCOMPLEMENT,    // always p7_NOCOMPLEMENT here;  varies in FM-based filter
                           L
                           );

        i = target_end; // skip forward
      }
    } /* end loop over sequence residues 1..L */
  return eslOK;
}
#endif

//Stubs so that functions exist if compiler can't handle AVX-512
#ifndef eslENABLE_AVX512


int
p7_SSVFilter_longtarget_avx512(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *ssvdata,
                        P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist)
{
  return eslEUNSUPPORTEDISA;
}


int
p7_SSVFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  return eslEUNSUPPORTEDISA;
}

int
p7_SSVFilter_avx512_unrolled(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  return eslEUNSUPPORTEDISA;
}

int
p7_SSVFilter_test_all(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  return eslEUNSUPPORTEDISA;
}


#endif
