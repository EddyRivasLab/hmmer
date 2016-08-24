/* The SSV filter implementation: ARM NEON version.
 * 
 * Contents:
 *   1. Bjarne's magic macros
 *   2. p7_SSVFilter() implementation
 * 
 * Original author:    Bjarne Knudsen (CLCbio, Aarhus, Denmark)
 * Ported to ARM NEON: Tyler Camp (University of Texas, Austin)
 *
 * For Bjarne's notes, see the Intel SSE implementation.
 */
#include "p7_config.h"

#include <math.h>

#if p7_CPU_ARCH == arm || p7_CPU_ARCH == arm64
#include <arm_neon.h>
#endif /* arm/arm64 arch */

#include "easel.h"
#include "esl_neon.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/ssvfilter.h"

/*****************************************************************
 * 1. Bjarne's magic macros
 *****************************************************************/

/* Note that some ifdefs below has to be changed if these values are
   changed. These values are chosen based on some simple speed
   tests. Apparently, two registers are generally used for something
   else, leaving 14 registers on 64 bit versions and 6 registers on 32
   bit versions. */
#define  MAX_BANDS 6

#ifdef HAVE_NEON
#define STEP_SINGLE(sv)                         \
  sv.s8x16   = vqsubq_s8(sv.s8x16, (*rsc).s8x16); rsc++;        \
  xEv.u8x16  = vmaxq_u8(xEv.u8x16, sv.u8x16);


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
  rsc = om->sbv[dsq[i]] + pos;                                   \
  step()                                                        \
  sv.s8x16 = vextq_s8(beginv.s8x16, sv.s8x16, 15);  \
                             \
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
  register esl_neon_128i_t sv00 = beginv;

#define RESET_2()                               \
  RESET_1()                                     \
  register esl_neon_128i_t sv01 = beginv;

#define RESET_3()                               \
  RESET_2()                                     \
  register esl_neon_128i_t sv02 = beginv;

#define RESET_4()                               \
  RESET_3()                                     \
  register esl_neon_128i_t sv03 = beginv;

#define RESET_5()                               \
  RESET_4()                                     \
  register esl_neon_128i_t sv04 = beginv;

#define RESET_6()                               \
  RESET_5()                                     \
  register esl_neon_128i_t sv05 = beginv;

#define RESET_7()                               \
  RESET_6()                                     \
  register esl_neon_128i_t sv06 = beginv;

#define RESET_8()                               \
  RESET_7()                                     \
  register esl_neon_128i_t sv07 = beginv;

#define RESET_9()                               \
  RESET_8()                                     \
  register esl_neon_128i_t sv08 = beginv;

#define RESET_10()                              \
  RESET_9()                                     \
  register esl_neon_128i_t sv09 = beginv;

#define RESET_11()                              \
  RESET_10()                                    \
  register esl_neon_128i_t sv10 = beginv;

#define RESET_12()                              \
  RESET_11()                                    \
  register esl_neon_128i_t sv11 = beginv;

#define RESET_13()                              \
  RESET_12()                                    \
  register esl_neon_128i_t sv12 = beginv;

#define RESET_14()                              \
  RESET_13()                                    \
  register esl_neon_128i_t sv13 = beginv;

#define RESET_15()                              \
  RESET_14()                                    \
  register esl_neon_128i_t sv14 = beginv;

#define RESET_16()                              \
  RESET_15()                                    \
  register esl_neon_128i_t sv15 = beginv;

#define RESET_17()                              \
  RESET_16()                                    \
  register esl_neon_128i_t sv16 = beginv;

#define RESET_18()                              \
  RESET_17()                                    \
  register esl_neon_128i_t sv17 = beginv;


#define CALC(reset, step, convert, width)       \
  int i;                                        \
  int i2;                                       \
  int Q        = P7_NVB(om->M);                \
  esl_neon_128i_t *rsc;                                 \
                                                \
  int w = width;                                \
                                                \
  dsq++;                                        \
                                                \
  reset()                                       \
                                                \
  for (i = 0; i < L && i < Q - q - w; i++)      \
    {                                           \
      rsc = om->sbv[dsq[i]] + i + q;            \
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
         rsc = om->sbv[dsq[i2 + i]] + i;        \
         step()                                 \
       }                                        \
                                                \
     i += i2;                                   \
     convert(step, NO_CHECK, )                  \
   }                                            \
                                                \
 for (i = 0; i2 + i < L && i < Q - w; i++)      \
   {                                            \
     rsc = om->sbv[dsq[i2 + i]] + i;            \
     step()                                     \
   }                                            \
                                                \
 i+=i2;                                         \
 convert(step, LENGTH_CHECK, done2)             \
done2:                                          \
                                                \
 return xEv;


esl_neon_128i_t
calc_band_neon_1(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_1, STEP_BANDS_1, CONVERT_1, 1)
}

esl_neon_128i_t
calc_band_neon_2(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_2, STEP_BANDS_2, CONVERT_2, 2)
}

esl_neon_128i_t
calc_band_neon_3(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_3, STEP_BANDS_3, CONVERT_3, 3)
}

esl_neon_128i_t
calc_band_neon_4(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_4, STEP_BANDS_4, CONVERT_4, 4)
}

esl_neon_128i_t
calc_band_neon_5(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_5, STEP_BANDS_5, CONVERT_5, 5)
}

esl_neon_128i_t
calc_band_neon_6(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_6, STEP_BANDS_6, CONVERT_6, 6)
}

#if MAX_BANDS > 6 /* Only include needed functions to limit object file size */
esl_neon_128i_t
calc_band_neon_7(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_7, STEP_BANDS_7, CONVERT_7, 7)
}

esl_neon_128i_t
calc_band_neon_8(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_8, STEP_BANDS_8, CONVERT_8, 8)
}

esl_neon_128i_t
calc_band_neon_9(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_9, STEP_BANDS_9, CONVERT_9, 9)
}

esl_neon_128i_t
calc_band_neon_10(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_10, STEP_BANDS_10, CONVERT_10, 10)
}

esl_neon_128i_t
calc_band_neon_11(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_11, STEP_BANDS_11, CONVERT_11, 11)
}

esl_neon_128i_t
calc_band_neon_12(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_12, STEP_BANDS_12, CONVERT_12, 12)
}

esl_neon_128i_t
calc_band_neon_13(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_13, STEP_BANDS_13, CONVERT_13, 13)
}

esl_neon_128i_t
calc_band_neon_14(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_14, STEP_BANDS_14, CONVERT_14, 14)
}
#endif /* MAX_BANDS > 6 */
#if MAX_BANDS > 14
esl_neon_128i_t
calc_band_neon_15(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_15, STEP_BANDS_15, CONVERT_15, 15)
}

esl_neon_128i_t
calc_band_neon_16(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_16, STEP_BANDS_16, CONVERT_16, 16)
}

esl_neon_128i_t
calc_band_neon_17(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_17, STEP_BANDS_17, CONVERT_17, 17)
}

esl_neon_128i_t
calc_band_neon_18(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, esl_neon_128i_t beginv, register esl_neon_128i_t xEv)
{
  CALC(RESET_18, STEP_BANDS_18, CONVERT_18, 18)
}
#endif /* MAX_BANDS > 14 */


/*****************************************************************
 * 2. p7_SSVFilter() implementation
 *****************************************************************/

uint8_t
get_xE_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om)
{
  esl_neon_128i_t xEv;		   /* E state: keeps max for Mk->E as we go                     */
  esl_neon_128i_t beginv;          /* begin scores                                              */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = P7_NVB(om->M);    /* segment length: # of vectors                              */
  int bands;                       /* the number of bands (rounds) to use                       */
  int last_q = 0;                  /* for saving the last q value to find band width            */
  int i;                           /* counter for bands                                         */

  /* function pointers for the various number of vectors to use */
  esl_neon_128i_t (*fs[MAX_BANDS + 1]) (const ESL_DSQ *, int, const P7_OPROFILE *, int, register esl_neon_128i_t, esl_neon_128i_t)
    = {NULL
       , calc_band_neon_1,  calc_band_neon_2,  calc_band_neon_3,  calc_band_neon_4,  calc_band_neon_5,  calc_band_neon_6
#if MAX_BANDS > 6
       , calc_band_neon_7,  calc_band_neon_8,  calc_band_neon_9,  calc_band_neon_10, calc_band_neon_11, calc_band_neon_12, calc_band_neon_13, calc_band_neon_14
#endif
#if MAX_BANDS > 14
       , calc_band_neon_15, calc_band_neon_16, calc_band_neon_17, calc_band_neon_18
#endif
  };

  beginv.s8x16 =  vdupq_n_s8(128);
  xEv    =  beginv;

  /* Use the highest number of bands but no more than MAX_BANDS */
  bands = (Q + MAX_BANDS - 1) / MAX_BANDS;

  for (i = 0; i < bands; i++) 
    {
      q      = (Q * (i + 1)) / bands;
      xEv    = fs[q-last_q](dsq, L, om, last_q, beginv, xEv);
      last_q = q;
    }
  return esl_neon_hmax_u8(xEv);
}


int
p7_SSVFilter_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  /* Use 16 bit values to avoid overflow due to moved baseline */
  uint16_t  xE;
  uint16_t  xJ;

  if (om->tjb_b + om->tbm_b + om->tec_b + om->bias_b >= 127) {
    /* the optimizations are not guaranteed to work under these
       conditions (see comments at start of file) */
    return eslENORESULT;
  }

  xE = get_xE_neon(dsq, L, om);
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
#endif /* HAVE_NEON */
#ifndef HAVE_NEON 
p7_SSVFilter_neon(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  return eslENORESULT;
}
#endif
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
