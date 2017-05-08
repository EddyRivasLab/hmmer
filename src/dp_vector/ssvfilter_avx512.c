/* SSV filter, x86 AVX-512 vector implementation.
 * Adapted from SSE version by Nick Carter
 * 
 * See ssvfilter.md for notes.
 *
 * This file is conditionally compiled when eslENABLE_AVX512 is defined.
 */
#include "p7_config.h"
#ifdef eslENABLE_AVX512

#include <x86intrin.h>
#include <math.h>

#include "easel.h"
#include "esl_avx512.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/ssvfilter.h"

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

#define STEP_SINGLE_AVX_512(sv_AVX_512)                         \
  sv_AVX_512   = _mm512_subs_epi8(sv_AVX_512, *rsc_AVX_512); rsc_AVX_512++;        \
  xEv_AVX_512  = _mm512_max_epu8(xEv_AVX_512, sv_AVX_512);

#define LENGTH_CHECK_AVX_512(label)                     \
  if (i_AVX_512 >= L) goto label;

#define NO_CHECK(label)

#define STEP_BANDS_1_AVX_512()                          \
  STEP_SINGLE_AVX_512(sv00_AVX_512) 

#define STEP_BANDS_2_AVX_512()                          \
  STEP_BANDS_1_AVX_512()                                \
  STEP_SINGLE_AVX_512(sv01_AVX_512)

#define STEP_BANDS_3_AVX_512()                          \
  STEP_BANDS_2_AVX_512()                                \
  STEP_SINGLE_AVX_512(sv02_AVX_512)

#define STEP_BANDS_4_AVX_512()                          \
  STEP_BANDS_3_AVX_512()                                \
  STEP_SINGLE_AVX_512(sv03_AVX_512)

#define STEP_BANDS_5_AVX_512()                          \
  STEP_BANDS_4_AVX_512()                                \
  STEP_SINGLE_AVX_512(sv04_AVX_512)

#define STEP_BANDS_6_AVX_512()                          \
  STEP_BANDS_5_AVX_512()                                \
  STEP_SINGLE_AVX_512(sv05_AVX_512)

#define STEP_BANDS_7_AVX_512()                          \
  STEP_BANDS_6_AVX_512()                                \
  STEP_SINGLE_AVX_512(sv06_AVX_512)

#define STEP_BANDS_8_AVX_512()                          \
  STEP_BANDS_7_AVX_512()                                \
  STEP_SINGLE_AVX_512(sv07_AVX_512)

#define STEP_BANDS_9_AVX_512()                          \
  STEP_BANDS_8_AVX_512()                                \
  STEP_SINGLE_AVX_512(sv08_AVX_512)

#define STEP_BANDS_10_AVX_512()                         \
  STEP_BANDS_9_AVX_512()                                \
  STEP_SINGLE_AVX_512(sv09_AVX_512)

#define STEP_BANDS_11_AVX_512()                         \
  STEP_BANDS_10_AVX_512()                               \
  STEP_SINGLE_AVX_512(sv10_AVX_512)

#define STEP_BANDS_12_AVX_512()                         \
  STEP_BANDS_11_AVX_512()                               \
  STEP_SINGLE_AVX_512(sv11_AVX_512)

#define STEP_BANDS_13_AVX_512()                         \
  STEP_BANDS_12_AVX_512()                               \
  STEP_SINGLE_AVX_512(sv12_AVX_512)

#define STEP_BANDS_14_AVX_512()                         \
  STEP_BANDS_13_AVX_512()                               \
  STEP_SINGLE_AVX_512(sv13_AVX_512)

#define STEP_BANDS_15_AVX_512()                         \
  STEP_BANDS_14_AVX_512()                               \
  STEP_SINGLE_AVX_512(sv14_AVX_512)

#define STEP_BANDS_16_AVX_512()                         \
  STEP_BANDS_15_AVX_512()                               \
  STEP_SINGLE_AVX_512(sv15_AVX_512)

#define STEP_BANDS_17_AVX_512()                         \
  STEP_BANDS_16_AVX_512()                               \
  STEP_SINGLE_AVX_512(sv16_AVX_512)

#define STEP_BANDS_18_AVX_512()                         \
  STEP_BANDS_17_AVX_512()                               \
  STEP_SINGLE_AVX_512(sv17_AVX_512) 

#define CONVERT_STEP_AVX_512(step, length_check, label, sv_AVX_512, pos)        \
  length_check(label)                                           \
  rsc_AVX_512 = om->sbv_AVX_512[dsq[i_AVX_512]] + pos;                                   \
  step()       \
  sv_AVX_512 = esl_avx512_leftshift_one(sv_AVX_512); \
  sv_AVX_512 = _mm512_or_si512(sv_AVX_512, low_byte_128_AVX_512);                                \
  i_AVX_512++;

#define CONVERT_1_AVX_512(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv00_AVX_512, Q_AVX_512 - 1) 

#define CONVERT_2_AVX_512(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv01_AVX_512, Q_AVX_512 - 2)  \
  CONVERT_1_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_3_AVX_512(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv02_AVX_512, Q_AVX_512 - 3)  \
  CONVERT_2_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_4_AVX_512(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv03_AVX_512, Q_AVX_512 - 4)  \
  CONVERT_3_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_5_AVX_512(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv04_AVX_512, Q_AVX_512 - 5)  \
  CONVERT_4_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_6_AVX_512(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv05_AVX_512, Q_AVX_512 - 6)  \
  CONVERT_5_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_7_AVX_512(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv06_AVX_512, Q_AVX_512 - 7)  \
  CONVERT_6_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_8_AVX_512(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv07_AVX_512, Q_AVX_512 - 8)  \
  CONVERT_7_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_9_AVX_512(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv08_AVX_512, Q_AVX_512 - 9)  \
  CONVERT_8_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_10_AVX_512(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv09_AVX_512, Q_AVX_512 - 10) \
  CONVERT_9_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_11_AVX_512(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv10_AVX_512, Q_AVX_512 - 11) \
  CONVERT_10_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_12_AVX_512(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv11_AVX_512, Q_AVX_512 - 12) \
  CONVERT_11_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_13_AVX_512(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv12_AVX_512, Q_AVX_512 - 13) \
  CONVERT_12_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_14_AVX_512(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv13_AVX_512, Q_AVX_512 - 14) \
  CONVERT_13_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_15_AVX_512(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv14_AVX_512, Q_AVX_512 - 15) \
  CONVERT_14_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_16_AVX_512(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv15_AVX_512, Q_AVX_512 - 16) \
  CONVERT_15_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_17_AVX_512(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv16_AVX_512, Q_AVX_512 - 17) \
  CONVERT_16_AVX_512(step, LENGTH_CHECK, label)

#define CONVERT_18_AVX_512(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX_512(step, LENGTH_CHECK, label, sv17_AVX_512, Q_AVX_512 - 18) \
  CONVERT_17_AVX_512(step, LENGTH_CHECK, label)

#define RESET_1_AVX_512()                               \
  register __m512i sv00_AVX_512 = beginv_AVX_512;

#define RESET_2_AVX_512()                              \
  RESET_1_AVX_512()                                     \
  register __m512i sv01_AVX_512 = beginv_AVX_512;

#define RESET_3_AVX_512()                               \
  RESET_2_AVX_512()                                     \
  register __m512i sv02_AVX_512 = beginv_AVX_512;

#define RESET_4_AVX_512()                               \
  RESET_3_AVX_512()                                     \
  register __m512i sv03_AVX_512 = beginv_AVX_512;

#define RESET_5_AVX_512()                               \
  RESET_4_AVX_512()                                     \
  register __m512i sv04_AVX_512 = beginv_AVX_512;

#define RESET_6_AVX_512()                           \
   RESET_5_AVX_512()                                     \
  register __m512i sv05_AVX_512 = beginv_AVX_512;

#define RESET_7_AVX_512()                               \
  RESET_6_AVX_512()                                     \
  register __m512i sv06_AVX_512 = beginv_AVX_512;

#define RESET_8_AVX_512()                               \
  RESET_7_AVX_512()                                     \
  register __m512i sv07_AVX_512 = beginv_AVX_512;

#define RESET_9_AVX_512()                               \
  RESET_8_AVX_512()                                     \
  register __m512i sv08_AVX_512 = beginv_AVX_512;

#define RESET_10_AVX_512()                               \
  RESET_9_AVX_512()                                     \
  register __m512i sv09_AVX_512 = beginv_AVX_512;

#define RESET_11_AVX_512()                               \
  RESET_10_AVX_512()                                     \
  register __m512i sv10_AVX_512 = beginv_AVX_512;

#define RESET_12_AVX_512()                               \
  RESET_11_AVX_512()                                     \
  register __m512i sv11_AVX_512 = beginv_AVX_512;

#define RESET_13_AVX_512()                               \
  RESET_12_AVX_512()                                     \
  register __m512i sv12_AVX_512 = beginv_AVX_512;

#define RESET_14_AVX_512()                               \
  RESET_13_AVX_512()                                     \
  register __m512i sv13_AVX_512 = beginv_AVX_512;

#define RESET_15_AVX_512()                               \
  RESET_14_AVX_512()                                     \
  register __m512i sv14_AVX_512 = beginv_AVX_512;

#define RESET_16_AVX_512()                               \
  RESET_15_AVX_512()                                     \
  register __m512i sv15_AVX_512 = beginv_AVX_512;

#define RESET_17_AVX_512()                               \
  RESET_16_AVX_512()                                     \
  register __m512i sv16_AVX_512 = beginv_AVX_512;

#define RESET_18_AVX_512()                               \
  RESET_17_AVX_512()                                     \
  register __m512i sv17_AVX_512 = beginv_AVX_512;

#define CALC_AVX_512(reset, step, convert, width)       \
  int i_AVX_512;                                        \
  int i2_AVX_512;                                       \
  int Q_AVX_512        = P7_NVB_AVX_512(om->M);                \
  __m512i *rsc_AVX_512;                                 \
  __m512i low_byte_128_AVX_512; \
  __mmask64 low_byte_mask = 0x0000000000000001; \
  low_byte_128_AVX_512 = _mm512_setzero_si512();                                            \
  low_byte_128_AVX_512 = _mm512_mask_blend_epi8(low_byte_mask, low_byte_128_AVX_512, beginv_AVX_512);                     \
  int w_AVX_512 = width;                                \
uint32_t col_AVX_512;                                              \
  dsq++;                                        \
                                                \
  reset()                                       \
  int num_iters_AVX_512; \
  if (L <= Q_AVX_512- q_AVX_512 -w_AVX_512){ \
  	num_iters_AVX_512 = L;  \
  }  \
  else{  \
  	num_iters_AVX_512 = Q_AVX_512 -q_AVX_512 -w_AVX_512; \
  } \
  i_AVX_512 = 0; \
  while(num_iters_AVX_512 >=4){ \
  	rsc_AVX_512 = om->sbv_AVX_512[dsq[i_AVX_512]] + i_AVX_512 + q_AVX_512;            \
    step()    \
    i_AVX_512++; \
	rsc_AVX_512 = om->sbv_AVX_512[dsq[i_AVX_512]] + i_AVX_512 + q_AVX_512;            \
    step()    \
    i_AVX_512++; \
    	rsc_AVX_512 = om->sbv_AVX_512[dsq[i_AVX_512]] + i_AVX_512 + q_AVX_512;            \
    step()    \
    i_AVX_512++; \
    rsc_AVX_512 = om->sbv_AVX_512[dsq[i_AVX_512]] + i_AVX_512 + q_AVX_512;            \
    step()    \
    i_AVX_512++; \
    num_iters_AVX_512 -= 4; \
  } \
  while(num_iters_AVX_512 >0){ \
  	  rsc_AVX_512 = om->sbv_AVX_512[dsq[i_AVX_512]] + i_AVX_512 + q_AVX_512;            \
      step()                                 \
      i_AVX_512++; \
      num_iters_AVX_512--; \
  }  \
  i_AVX_512 = Q_AVX_512 - q_AVX_512 - w_AVX_512;                                \
  convert(step, LENGTH_CHECK_AVX_512, done1)            \
done1:                                          \
 for (i2_AVX_512 = Q_AVX_512 - q_AVX_512; i2_AVX_512 < L - Q_AVX_512; i2_AVX_512 += Q_AVX_512)          \
   {                                            \
  	i_AVX_512 = 0; \
   	num_iters_AVX_512 = Q_AVX_512 - w_AVX_512; \
   	while (num_iters_AVX_512 >= 4){  \
       rsc_AVX_512 = om->sbv_AVX_512[dsq[i2_AVX_512 + i_AVX_512]] + i_AVX_512;        \
       step() \
       i_AVX_512++; \
       rsc_AVX_512 = om->sbv_AVX_512[dsq[i2_AVX_512 + i_AVX_512]] + i_AVX_512;        \
       step() \
       i_AVX_512++; \
       rsc_AVX_512 = om->sbv_AVX_512[dsq[i2_AVX_512 + i_AVX_512]] + i_AVX_512;        \
       step() \
       i_AVX_512++; \
       rsc_AVX_512 = om->sbv_AVX_512[dsq[i2_AVX_512 + i_AVX_512]] + i_AVX_512;        \
       step() \
       i_AVX_512++; \
       num_iters_AVX_512-= 4; \
   	}\
  	while(num_iters_AVX_512 > 0){ \
  		 rsc_AVX_512 = om->sbv_AVX_512[dsq[i2_AVX_512 + i_AVX_512]] + i_AVX_512;        \
         step()    \
         i_AVX_512++; \
         num_iters_AVX_512--;      \
  	} \
                           \
     i_AVX_512 += i2_AVX_512;                                   \
     convert(step, NO_CHECK, )                  \
   }                                            \
	if((L - i2_AVX_512) < (Q_AVX_512-w_AVX_512)){                                         \
 		num_iters_AVX_512 = L -i2_AVX_512; \
 		} \
 	else{  \
 		num_iters_AVX_512 = Q_AVX_512 - w_AVX_512; \
 	} \
 	i_AVX_512 = 0; \
 	while (num_iters_AVX_512 >= 4){ \
     	rsc_AVX_512 = om->sbv_AVX_512[dsq[i2_AVX_512 + i_AVX_512]] + i_AVX_512;            \
     	step()                                     \
     	i_AVX_512+= 1;  \
     	rsc_AVX_512 = om->sbv_AVX_512[dsq[i2_AVX_512 + i_AVX_512]] + i_AVX_512;            \
     	step()                                     \
     	i_AVX_512+= 1;  \
     	rsc_AVX_512 = om->sbv_AVX_512[dsq[i2_AVX_512 + i_AVX_512]] + i_AVX_512;            \
     	step()                                     \
     	i_AVX_512+= 1;  \
     	rsc_AVX_512 = om->sbv_AVX_512[dsq[i2_AVX_512 + i_AVX_512]] + i_AVX_512;            \
     	step()                                     \
     	i_AVX_512+= 1;  \
     	num_iters_AVX_512 -= 4; \
   	}                                            \
   	while(num_iters_AVX_512 > 0) {  \
   	 	rsc_AVX_512 = om->sbv_AVX_512[dsq[i2_AVX_512 + i_AVX_512]] + i_AVX_512;            \
     	step()                                     \
     	i_AVX_512+= 1;  \
     	num_iters_AVX_512--; \
   	} \
                                       \
 i_AVX_512+=i2_AVX_512;                                         \
 convert(step, LENGTH_CHECK_AVX_512, done2)             \
done2:                                          \
	return xEv_AVX_512;

__m512i
calc_band_1_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_1_AVX_512, STEP_BANDS_1_AVX_512, CONVERT_1_AVX_512, 1)
}

__m512i
calc_band_2_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_2_AVX_512, STEP_BANDS_2_AVX_512, CONVERT_2_AVX_512, 2)
}

__m512i
calc_band_3_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_3_AVX_512, STEP_BANDS_3_AVX_512, CONVERT_3_AVX_512, 3)
}

__m512i
calc_band_4_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_4_AVX_512, STEP_BANDS_4_AVX_512, CONVERT_4_AVX_512, 4)
}

__m512i
calc_band_5_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_5_AVX_512, STEP_BANDS_5_AVX_512, CONVERT_5_AVX_512, 5)
}

__m512i
calc_band_6_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_6_AVX_512, STEP_BANDS_6_AVX_512, CONVERT_6_AVX_512, 6)
}

#if MAX_BANDS > 6 /* Only include needed functions to limit object file size */
__m512i
calc_band_7_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_7_AVX_512, STEP_BANDS_7_AVX_512, CONVERT_7_AVX_512, 7)
}

__m512i
calc_band_8_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_8_AVX_512, STEP_BANDS_8_AVX_512, CONVERT_8_AVX_512, 8)
}

__m512i
calc_band_9_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_9_AVX_512, STEP_BANDS_9_AVX_512, CONVERT_9_AVX_512, 9)
}

__m512i
calc_band_10_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_10_AVX_512, STEP_BANDS_10_AVX_512, CONVERT_10_AVX_512, 10)
}

__m512i
calc_band_11_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_11_AVX_512, STEP_BANDS_11_AVX_512, CONVERT_11_AVX_512, 11)
}

__m512i
calc_band_12_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_12_AVX_512, STEP_BANDS_12_AVX_512, CONVERT_12_AVX_512, 12)
}

__m512i
calc_band_13_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_13_AVX_512, STEP_BANDS_13_AVX_512, CONVERT_13_AVX_512, 13)
}

__m512i
calc_band_14_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_14_AVX_512, STEP_BANDS_14_AVX_512, CONVERT_14_AVX_512, 14)
}
#endif /* MAX_BANDS > 6 */
#if MAX_BANDS > 14
__m512i
calc_band_15_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_15_AVX_512, STEP_BANDS_15_AVX_512, CONVERT_15_AVX_512, 15)
}

__m512i
calc_band_16_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_16_AVX_512, STEP_BANDS_16_AVX_512, CONVERT_16_AVX_512, 16)
}

__m512i
calc_band_17_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_17_AVX_512, STEP_BANDS_17_AVX_512, CONVERT_17_AVX_512, 17)
}

__m512i
calc_band_18_AVX_512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX_512, __m512i beginv_AVX_512, register __m512i xEv_AVX_512)
{
  CALC_AVX_512(RESET_18_AVX_512, STEP_BANDS_18_AVX_512, CONVERT_18_AVX_512, 18)
}
#endif /* MAX_BANDS > 14 */



uint8_t
get_xE_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om)
{
  __m512i xEv_AVX_512;		                  // E state: keeps max for Mk->E as we go   
  __m512i beginv_AVX_512;                         // begin scores 
  uint8_t retval_AVX_512;
  int q_AVX_512;			          // counter over vectors 0..nq-1
  int Q_AVX_512        = P7_NVB_AVX_512(om->M);   // segment length: # of vectors
  int bands_AVX_512;                              // the number of bands (rounds) to use
  beginv_AVX_512 =  _mm512_set1_epi8(128);
  xEv_AVX_512    =  beginv_AVX_512;

  /* function pointers for the various number of vectors to use */
  __m512i (*fs_AVX_512[MAX_BANDS + 1]) (const ESL_DSQ *, int, const P7_OPROFILE *, int, register __m512i, __m512i)
    = {NULL
       , calc_band_1_AVX_512,  calc_band_2_AVX_512,  calc_band_3_AVX_512,  calc_band_4_AVX_512,  calc_band_5_AVX_512,  calc_band_6_AVX_512
#if MAX_BANDS > 6
       , calc_band_7_AVX_512,  calc_band_8_AVX_512,  calc_band_9_AVX_512,  calc_band_10_AVX_512, calc_band_11_AVX_512, calc_band_12_AVX_512,
        calc_band_13_AVX_512, calc_band_14_AVX_512
#endif
#if MAX_BANDS > 14
       , calc_band_15_AVX_512, calc_band_16_AVX_512, calc_band_17_AVX_512, calc_band_18_AVX_512
#endif
  };
  int last_q;                  // for saving the last q value to find band width 
  int i;                       // counter for bands                              

  last_q = 0;
  /* Use the highest number of bands but no more than MAX_BANDS */
  bands_AVX_512 = (Q_AVX_512 + MAX_BANDS - 1) / MAX_BANDS;
  for (i = 0; i < bands_AVX_512; i++) {
    q_AVX_512 = (Q_AVX_512 * (i + 1)) / bands_AVX_512;

    xEv_AVX_512 = fs_AVX_512[q_AVX_512-last_q](dsq, L, om, last_q, beginv_AVX_512, xEv_AVX_512);

    last_q = q_AVX_512;
  }
  retval_AVX_512 = esl_avx512_hmax_epu8(xEv_AVX_512); // assign this here to allow checking vs AVX, AVX-512

  return retval_AVX_512;
}

int
p7_SSVFilter_avx512(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  /* Use 16 bit values to avoid overflow due to moved baseline */
  uint16_t  xE;
  uint16_t  xJ;

  if (om->tjb_b + om->tbm_b + om->tec_b + om->bias_b >= 127) {
    /* the optimizations are not guaranteed to work under these
       conditions (see ssvfilter.md) */
    return eslENORESULT;
  }

  xE = get_xE_avx512(dsq, L, om);

  // SRE TODO REVISIT : All this needs to be rechecked and rethought for H4.
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



#else // ! eslENABLE_AVX512

/* Standard compiler-pleasing mantra for an #ifdef'd-out, empty code file. */
void p7_ssvfilter_avx512_silence_hack(void) { return; }
#if defined p7SSVFILTER_AVX512_TESTDRIVE || p7SSVFILTER_AVX512_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_AVX512 or not

