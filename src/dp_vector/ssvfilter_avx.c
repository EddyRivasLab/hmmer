/* SSV filter, x86 AVX vector implementation.
 * Adapted from SSE version by Nick Carter
 * 
 * See ssvfilter.md for notes.
 *
 * This file is conditionally compiled when eslENABLE_AVX is defined.
 */
#include "p7_config.h"
#ifdef eslENABLE_AVX

#include <x86intrin.h>
#include <math.h>

#include "easel.h"
#include "esl_avx.h"

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

#define STEP_SINGLE_AVX(sv_AVX) \
  sv_AVX   = _mm256_subs_epi8(sv_AVX, *rsc_AVX); rsc_AVX++;        \
  xEv_AVX  = _mm256_max_epu8(xEv_AVX, sv_AVX);	

#define LENGTH_CHECK_AVX(label)                     \
  if (i_AVX >= L) goto label;

#define NO_CHECK(label)

#define STEP_BANDS_1_AVX()                          \
  STEP_SINGLE_AVX(sv00_AVX)

#define STEP_BANDS_2_AVX()                          \
  STEP_BANDS_1_AVX()                                \
  STEP_SINGLE_AVX(sv01_AVX)

#define STEP_BANDS_3_AVX()                          \
  STEP_BANDS_2_AVX()                                \
  STEP_SINGLE_AVX(sv02_AVX)

#define STEP_BANDS_4_AVX()                          \
  STEP_BANDS_3_AVX()                                \
  STEP_SINGLE_AVX(sv03_AVX)

#define STEP_BANDS_5_AVX()                          \
  STEP_BANDS_4_AVX()                                \
  STEP_SINGLE_AVX(sv04_AVX)

#define STEP_BANDS_6_AVX()                          \
  STEP_BANDS_5_AVX()                                \
  STEP_SINGLE_AVX(sv05_AVX)

#define STEP_BANDS_7_AVX()                          \
  STEP_BANDS_6_AVX()                                \
  STEP_SINGLE_AVX(sv06_AVX)

#define STEP_BANDS_8_AVX()                          \
  STEP_BANDS_7_AVX()                                \
  STEP_SINGLE_AVX(sv07_AVX)

#define STEP_BANDS_9_AVX()                          \
  STEP_BANDS_8_AVX()                                \
  STEP_SINGLE_AVX(sv08_AVX)

#define STEP_BANDS_10_AVX()                         \
  STEP_BANDS_9_AVX()                                \
  STEP_SINGLE_AVX(sv09_AVX)

#define STEP_BANDS_11_AVX()                         \
  STEP_BANDS_10_AVX()                               \
  STEP_SINGLE_AVX(sv10_AVX)

#define STEP_BANDS_12_AVX()                         \
  STEP_BANDS_11_AVX( )                               \
  STEP_SINGLE_AVX(sv11_AVX)

#define STEP_BANDS_13_AVX()                         \
  STEP_BANDS_12_AVX()                               \
  STEP_SINGLE_AVX(sv12_AVX)

#define STEP_BANDS_14_AVX()                         \
  STEP_BANDS_13_AVX()                               \
  STEP_SINGLE_AVX(sv13_AVX)

#define STEP_BANDS_15_AVX()                         \
  STEP_BANDS_14_AVX()                               \
  STEP_SINGLE_AVX(sv14_AVX)

#define STEP_BANDS_16_AVX()                         \
  STEP_BANDS_15_AVX()                               \
  STEP_SINGLE_AVX(sv15_AVX)

#define STEP_BANDS_17_AVX()                         \
  STEP_BANDS_16_AVX()                               \
  STEP_SINGLE_AVX(sv16_AVX)

#define STEP_BANDS_18_AVX()                         \
  STEP_BANDS_17_AVX()                               \
  STEP_SINGLE_AVX(sv17_AVX)

#define CONVERT_STEP_AVX(step, length_check, label, sv_AVX, pos)        \
  length_check(label)                                           \
  rsc_AVX = om->sbv_AVX[dsq[i_AVX]] + pos;                                   \
  step()       \
  sv_AVX = esl_avx_leftshift_one(sv_AVX); \
  sv_AVX = _mm256_or_si256(sv_AVX, low_byte_128_AVX);                                \
  i_AVX++;

#define CONVERT_1_AVX(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv00_AVX, Q_AVX - 1)

#define CONVERT_2_AVX(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv01_AVX, Q_AVX - 2)  \
  CONVERT_1_AVX(step, LENGTH_CHECK, label)

#define CONVERT_3_AVX(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv02_AVX, Q_AVX - 3)  \
  CONVERT_2_AVX(step, LENGTH_CHECK, label)

#define CONVERT_4_AVX(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv03_AVX, Q_AVX - 4)  \
  CONVERT_3_AVX(step, LENGTH_CHECK, label)

#define CONVERT_5_AVX(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv04_AVX, Q_AVX - 5)  \
  CONVERT_4_AVX(step, LENGTH_CHECK, label)

#define CONVERT_6_AVX(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv05_AVX, Q_AVX - 6)  \
  CONVERT_5_AVX(step, LENGTH_CHECK, label)

#define CONVERT_7_AVX(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv06_AVX, Q_AVX - 7)  \
  CONVERT_6_AVX(step, LENGTH_CHECK, label)

#define CONVERT_8_AVX(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv07_AVX, Q_AVX - 8)  \
  CONVERT_7_AVX(step, LENGTH_CHECK, label)

#define CONVERT_9_AVX(step, LENGTH_CHECK, label)            \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv08_AVX, Q_AVX - 9)  \
  CONVERT_8_AVX(step, LENGTH_CHECK, label)

#define CONVERT_10_AVX(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv09_AVX, Q_AVX - 10) \
  CONVERT_9_AVX(step, LENGTH_CHECK, label)

#define CONVERT_11_AVX(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv10_AVX, Q_AVX - 11) \
  CONVERT_10_AVX(step, LENGTH_CHECK, label)

#define CONVERT_12_AVX(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv11_AVX, Q_AVX - 12) \
  CONVERT_11_AVX(step, LENGTH_CHECK, label)

#define CONVERT_13_AVX(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv12_AVX, Q_AVX - 13) \
  CONVERT_12_AVX(step, LENGTH_CHECK, label)

#define CONVERT_14_AVX(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv13_AVX, Q_AVX - 14) \
  CONVERT_13_AVX(step, LENGTH_CHECK, label)

#define CONVERT_15_AVX(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv14_AVX, Q_AVX - 15) \
  CONVERT_14_AVX(step, LENGTH_CHECK, label)

#define CONVERT_16_AVX(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv15_AVX, Q_AVX - 16) \
  CONVERT_15_AVX(step, LENGTH_CHECK, label)

#define CONVERT_17_AVX(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv16_AVX, Q_AVX - 17) \
  CONVERT_16_AVX(step, LENGTH_CHECK, label)

#define CONVERT_18_AVX(step, LENGTH_CHECK, label)           \
  CONVERT_STEP_AVX(step, LENGTH_CHECK, label, sv17_AVX, Q_AVX - 18) \
  CONVERT_17_AVX(step, LENGTH_CHECK, label)

#define RESET_1_AVX()                               \
  register __m256i sv00_AVX = beginv_AVX;

#define RESET_2_AVX()                               \
  RESET_1_AVX()                                     \
  register __m256i sv01_AVX = beginv_AVX;

#define RESET_3_AVX()                               \
  RESET_2_AVX()                                     \
  register __m256i sv02_AVX = beginv_AVX;

#define RESET_4_AVX()                               \
  RESET_3_AVX()                                     \
  register __m256i sv03_AVX = beginv_AVX;

#define RESET_5_AVX()                               \
  RESET_4_AVX()                                     \
  register __m256i sv04_AVX = beginv_AVX;

#define RESET_6_AVX()                               \
  RESET_5_AVX()                                     \
  register __m256i sv05_AVX = beginv_AVX;

#define RESET_7_AVX()                               \
  RESET_6_AVX()                                     \
  register __m256i sv06_AVX = beginv_AVX;

#define RESET_8_AVX()                               \
  RESET_7_AVX()                                     \
  register __m256i sv07_AVX = beginv_AVX;

#define RESET_9_AVX()                               \
  RESET_8_AVX()                                     \
  register __m256i sv08_AVX = beginv_AVX;

#define RESET_10_AVX()                              \
  RESET_9_AVX()                                     \
  register __m256i sv09_AVX = beginv_AVX;

#define RESET_11_AVX()                              \
  RESET_10_AVX()                                    \
  register __m256i sv10_AVX = beginv_AVX;

#define RESET_12_AVX()                              \
  RESET_11_AVX()                                    \
  register __m256i sv11_AVX = beginv_AVX;

#define RESET_13_AVX()                              \
  RESET_12_AVX()                                    \
  register __m256i sv12_AVX = beginv_AVX;

#define RESET_14_AVX()                              \
  RESET_13_AVX()                                    \
  register __m256i sv13_AVX = beginv_AVX;

#define RESET_15_AVX()                              \
  RESET_14_AVX()                                    \
  register __m256i sv14_AVX = beginv_AVX;

#define RESET_16_AVX()                              \
  RESET_15_AVX()                                    \
  register __m256i sv15_AVX = beginv_AVX;

#define RESET_17_AVX()                              \
  RESET_16_AVX()                                    \
  register __m256i sv16_AVX = beginv_AVX;

#define RESET_18_AVX()                              \
  RESET_17_AVX()                                    \
  register __m256i sv17_AVX = beginv_AVX;

#define CALC_AVX(reset, step, convert, width)       \
  int i_AVX;                                        \
  int i2_AVX;                                       \
  int Q_AVX        = P7_NVB_AVX(om->M);                \
  __m256i *rsc_AVX;                                 \
  __m256i low_byte_128_AVX = _mm256_setzero_si256();                                            \
  low_byte_128_AVX = _mm256_insert_epi8(low_byte_128_AVX, 128, 0); \
  int w_AVX = width;                                \
  uint32_t col_AVX;                                     \
  dsq++;                                        \
                                                \
  reset()                                       \
  int num_iters_AVX; \
  if (L <= Q_AVX- q_AVX -w_AVX){ \
  	num_iters_AVX = L;  \
  }  \
  else{  \
  	num_iters_AVX = Q_AVX -q_AVX -w_AVX; \
  } \
  i_AVX = 0; \
  while(num_iters_AVX >=4){ \
  	rsc_AVX = om->sbv_AVX[dsq[i_AVX]] + i_AVX + q_AVX;            \
    step()    \
    i_AVX++; \
	rsc_AVX = om->sbv_AVX[dsq[i_AVX]] + i_AVX + q_AVX;            \
    step()    \
    i_AVX++; \
    	rsc_AVX = om->sbv_AVX[dsq[i_AVX]] + i_AVX + q_AVX;            \
    step()    \
    i_AVX++; \
    rsc_AVX = om->sbv_AVX[dsq[i_AVX]] + i_AVX + q_AVX;            \
    step()    \
    i_AVX++; \
    num_iters_AVX -= 4; \
  } \
  while(num_iters_AVX >0){ \
  	  rsc_AVX = om->sbv_AVX[dsq[i_AVX]] + i_AVX + q_AVX;            \
      step()                                 \
      i_AVX++; \
      num_iters_AVX--; \
  }  \
  i_AVX = Q_AVX - q_AVX - w_AVX;                                \
  convert(step, LENGTH_CHECK_AVX, done1_AVX)            \
done1_AVX:                                          \
 for (i2_AVX = Q_AVX - q_AVX; i2_AVX < L - Q_AVX; i2_AVX += Q_AVX)          \
   {                                            \
   	i_AVX = 0; \
   	num_iters_AVX = Q_AVX - w_AVX; \
   	while (num_iters_AVX >= 4){  \
       rsc_AVX = om->sbv_AVX[dsq[i2_AVX + i_AVX]] + i_AVX;        \
       step() \
       i_AVX++; \
       rsc_AVX = om->sbv_AVX[dsq[i2_AVX + i_AVX]] + i_AVX;        \
       step() \
       i_AVX++; \
       rsc_AVX = om->sbv_AVX[dsq[i2_AVX + i_AVX]] + i_AVX;        \
       step() \
       i_AVX++; \
       rsc_AVX = om->sbv_AVX[dsq[i2_AVX + i_AVX]] + i_AVX;        \
       step() \
       i_AVX++; \
       num_iters_AVX-= 4; \
   	}\
  	while(num_iters_AVX > 0){ \
  		 rsc_AVX = om->sbv_AVX[dsq[i2_AVX + i_AVX]] + i_AVX;        \
         step()    \
         i_AVX++; \
         num_iters_AVX--;      \
  	} \
                                                \
     i_AVX += i2_AVX;                                   \
     convert(step, NO_CHECK, ) \
   }                                            \
if((L - i2_AVX) < (Q_AVX-w_AVX)){                                         \
 	num_iters_AVX = L -i2_AVX; \
 	} \
 else{  \
 	num_iters_AVX = Q_AVX - w_AVX; \
 } \
 i_AVX = 0; \
 while (num_iters_AVX >= 4){ \
     rsc_AVX = om->sbv_AVX[dsq[i2_AVX + i_AVX]] + i_AVX;            \
     step()                                     \
     i_AVX+= 1;  \
     rsc_AVX = om->sbv_AVX[dsq[i2_AVX + i_AVX]] + i_AVX;            \
     step()                                     \
     i_AVX+= 1;  \
     rsc_AVX = om->sbv_AVX[dsq[i2_AVX + i_AVX]] + i_AVX;            \
     step()                                     \
     i_AVX+= 1;  \
     rsc_AVX = om->sbv_AVX[dsq[i2_AVX + i_AVX]] + i_AVX;            \
     step()                                     \
     i_AVX+= 1;  \
     num_iters_AVX -= 4; \
   }                                            \
   while(num_iters_AVX > 0) {  \
   	 rsc_AVX = om->sbv_AVX[dsq[i2_AVX + i_AVX]] + i_AVX;            \
     step()                                     \
     i_AVX+= 1;  \
     num_iters_AVX--; \
   } \
 i_AVX+=i2_AVX;                                         \
 convert(step, LENGTH_CHECK_AVX, done2_AVX)             \
done2_AVX:                                          \
                                                \
 return xEv_AVX;

__m256i
calc_band_1_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_1_AVX, STEP_BANDS_1_AVX, CONVERT_1_AVX, 1)
}

__m256i
calc_band_2_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_2_AVX, STEP_BANDS_2_AVX, CONVERT_2_AVX, 2)
}

__m256i
calc_band_3_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_3_AVX, STEP_BANDS_3_AVX, CONVERT_3_AVX, 3)
}

__m256i
calc_band_4_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_4_AVX, STEP_BANDS_4_AVX, CONVERT_4_AVX, 4)
}

__m256i
calc_band_5_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_5_AVX, STEP_BANDS_5_AVX, CONVERT_5_AVX, 5)
}

__m256i
calc_band_6_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_6_AVX, STEP_BANDS_6_AVX, CONVERT_6_AVX, 6)
}


#if MAX_BANDS > 6 /* Only include needed functions to limit object file size */
__m256i
calc_band_7_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_7_AVX, STEP_BANDS_7_AVX, CONVERT_7_AVX, 7)
}

__m256i
calc_band_8_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_8_AVX, STEP_BANDS_8_AVX, CONVERT_8_AVX, 8)
}

__m256i
calc_band_9_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_9_AVX, STEP_BANDS_9_AVX, CONVERT_9_AVX, 9)
}

__m256i
calc_band_10_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_10_AVX, STEP_BANDS_10_AVX, CONVERT_10_AVX, 10)
}

__m256i
calc_band_11_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_11_AVX, STEP_BANDS_11_AVX, CONVERT_11_AVX, 11)
}

__m256i
calc_band_12_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_12_AVX, STEP_BANDS_12_AVX, CONVERT_12_AVX, 12)
}

__m256i
calc_band_13_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_13_AVX, STEP_BANDS_13_AVX, CONVERT_13_AVX, 13)
}

__m256i
calc_band_14_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_14_AVX, STEP_BANDS_14_AVX, CONVERT_14_AVX, 14)
}

#endif /* MAX_BANDS > 6 */
#if MAX_BANDS > 14
__m256i
calc_band_15_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_15_AVX, STEP_BANDS_15_AVX, CONVERT_15_AVX, 15)
}

__m256i
calc_band_16_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_16_AVX, STEP_BANDS_16_AVX, CONVERT_16_AVX, 16)
}

__m256i
calc_band_17_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_17_AVX, STEP_BANDS_17_AVX, CONVERT_17_AVX, 17)
}

__m256i
calc_band_18_AVX(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q_AVX, __m256i beginv_AVX, register __m256i xEv_AVX)
{
  CALC_AVX(RESET_18_AVX, STEP_BANDS_18_AVX, CONVERT_18_AVX, 18)
}
#endif /* MAX_BANDS > 14 */



uint8_t
get_xE_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om)
{
  __m256i xEv_AVX;		          // E state: keeps max for Mk->E as we go  
  __m256i beginv_AVX;                     // begin scores                           
  uint8_t retval_AVX;
  int q_AVX;			          // counter over vectors 0..nq-1                              
  int Q_AVX        = P7_NVB_AVX(om->M);   // segment length: # of vectors 
  int bands_AVX;                          // number of bands (rounds) to use
  beginv_AVX =  _mm256_set1_epi8(128);
  xEv_AVX    =  beginv_AVX;

  /* function pointers for the various number of vectors to use */
  __m256i (*fs_AVX[MAX_BANDS + 1]) (const ESL_DSQ *, int, const P7_OPROFILE *, int, register __m256i, __m256i)
    = {NULL
       , calc_band_1_AVX,  calc_band_2_AVX,  calc_band_3_AVX,  calc_band_4_AVX,  calc_band_5_AVX,  calc_band_6_AVX
#if MAX_BANDS > 6
       , calc_band_7_AVX,  calc_band_8_AVX,  calc_band_9_AVX,  calc_band_10_AVX, calc_band_11_AVX, calc_band_12_AVX
       , calc_band_13_AVX, calc_band_14_AVX
#endif
#if MAX_BANDS > 14
       , calc_band_15_AVX, calc_band_16_AVX, calc_band_17_AVX, calc_band_18_AVX
#endif
  };
  int last_q;                  // for saving the last q value to find band width            
  int i;                       // counter for bands   

  last_q = 0;   // reset in case we also ran SSE code
  /* Use the highest number of bands but no more than MAX_BANDS */
  bands_AVX = (Q_AVX + MAX_BANDS - 1) / MAX_BANDS;
  for (i = 0; i < bands_AVX; i++) {
    q_AVX = (Q_AVX * (i + 1)) / bands_AVX;

    xEv_AVX = fs_AVX[q_AVX-last_q](dsq, L, om, last_q, beginv_AVX, xEv_AVX);
    last_q = q_AVX;
  }

  retval_AVX = esl_avx_hmax_epu8(xEv_AVX);

  return retval_AVX;
}

int
p7_SSVFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  /* Use 16 bit values to avoid overflow due to moved baseline */
  uint16_t  xE;
  uint16_t  xJ;

  if (om->tjb_b + om->tbm_b + om->tec_b + om->bias_b >= 127) {
    /* the optimizations are not guaranteed to work under these
       conditions (see ssvfilter.md) */
    return eslENORESULT;
  }

  xE = get_xE_avx(dsq, L, om);

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



#else // ! eslENABLE_AVX

/* Standard compiler-pleasing mantra for an #ifdef'd-out, empty code file. */
void p7_ssvfilter_avx_silence_hack(void) { return; }
#if defined p7SSVFILTER_AVX_TESTDRIVE || p7SSVFILTER_AVX_EXAMPLE
int main(void) { return 0; }
#endif 
#endif // eslENABLE_AVX or not

