#ifndef P7_FM_SSE_INCLUDED
#define P7_FM_SSE_INCLUDED

#include "p7_config.h"

#include "esl_getopts.h"

#include "fm/fm.h"


/* Effectively global variables, to be initialized once in fm_initConfig(),
 * then passed around among threads to avoid recomputing them
 *
 * When allocated, must be 16-byte aligned, and all _m128i elements
 * must precede other types
 */
typedef struct {
  /* mask arrays, and 16-byte-offsets into them */
  __m128i *fm_masks_mem;
  __m128i *fm_masks_v;
  __m128i *fm_reverse_masks_mem;
  __m128i *fm_reverse_masks_v;
  __m128i *fm_chars_mem;
  __m128i *fm_chars_v;

  /*various precomputed vectors*/
  __m128i fm_allones_v;
  __m128i fm_zeros_v;
  __m128i fm_neg128_v;
  __m128i fm_m0f;  //00 00 11 11
  __m128i fm_m01;  //01 01 01 01
  __m128i fm_m11;  //00 00 00 11

  /* no non-__m128i- elements above this line */

  /*suffix-array mask and offset values*/
  int maskSA;
  int shiftSA;

  /*counter, to compute FM-index speed*/
  int occCallCnt;

  /*bounding cutoffs*/
  int max_depth;
  int neg_len_limit;
  int consec_pos_req; //6
  float score_ratio_req; //.49
  int msv_length;
  float max_scthreshFM;

  /*pointer to FM-index metadata*/
  FM_METADATA *meta;

} FM_CFG;


//used to convert from a byte array to an __m128i
typedef union {
        uint8_t bytes[16];
        __m128i m128;
        } byte_m128;


/* Gather the sum of all counts in a 16x8-bit element into a single 16-bit
 *  element of the register (the 0th element)
 *
 *  the _mm_sad_epu8  accumulates 8-bit counts into 16-bit counts:
 *      left 8 counts (64-bits) accumulate in counts_v[0],
 *      right 8 counts in counts_v[4]  (the other 6 16-bit ints are 0)
 *  the _mm_shuffle_epi32  flips the 4th int into the 0th slot
 */
#define FM_GATHER_8BIT_COUNTS( in_v, mid_v, out_v  ) do {\
    mid_v = _mm_sad_epu8 (in_v, cfg->fm_zeros_v);\
    tmp_v = _mm_shuffle_epi32(mid_v, _MM_SHUFFLE(1, 1, 1, 2));\
    out_v = _mm_add_epi16(mid_v, tmp_v);\
  } while (0)

// use this instead?
// D = _mm_srli_si128(C, 8);
// E = _mm_add_epi32(C, D);
// sad = mm_cvtsi128_si32(E)


/* Macro for SSE operations to turn 2-bit character values into 2-bit binary
 * (00 or 01) match/mismatch values representing occurrences of a character in a
 * 4-char-per-byte packed BWT.
 *
 * Typically followed by a call to FM_COUNT_SSE_4PACKED, possibly with a
 * mask in between to handle the case where we don't want to add over all
 * positions in the vector
 *
 * tmp_v and tmp2_v are used as temporary vectors throughout, and hold meaningless values
 * at the end
 *
 * xor(in_v, c_v)        : each 2-bit value will be 00 if a match, and non-0 if a mismatch
 * and(in_v, 01010101)   : look at the right bit of each 2-bit value,
 * srli(1)+and()         : look at the left bit of each 2-bit value,
 * or()                  : if either left bit or right bit is non-0, 01, else 00 (match is 00)
 *
 * subs()                : invert, so match is 01, mismatch is 00
 *
 */
#define FM_MATCH_2BIT(in_v, c_v, a_v, b_v, out_v) do {\
    a_v = _mm_xor_si128(in_v, c_v);\
    \
    b_v = _mm_and_si128(a_v, cfg->fm_m01);\
    a_v  = _mm_srli_epi16(a_v, 1);\
    a_v  = _mm_and_si128(a_v, cfg->fm_m01);\
    a_v  = _mm_or_si128(a_v, b_v);\
    \
    out_v  = _mm_subs_epi8(cfg->fm_m01,a_v);\
  } while (0)


/*Macro for SSE operations to count bits produced by FM_MATCH_SSE_4PACKED
 *
 * tmp_v and tmp2_v are used as temporary vectors throughout, and hold meaningless values
 * at the end
 *
 * then add up the 2-bit values:
 * srli(4)+add()         : left 4 bits shifted right, added to right 4 bits
 *
 * srli(2)+and(00000011) : left 2 bits (value 0..2) shifted right, masked, so no other bits active
 * and(00000011)         : right 2 bits (value 0..2) masked so no other bits active
 *
 * final 2 add()s        : tack current counts on to already-tabulated counts.
 */
#define FM_COUNT_2BIT(a_v, b_v, cnts_v) do {\
        b_v = _mm_srli_epi16(a_v, 4);\
        a_v  = _mm_add_epi16(a_v, b_v);\
        \
        b_v = _mm_srli_epi16(a_v, 2);\
        a_v  = _mm_and_si128(a_v,cfg->fm_m11);\
        b_v = _mm_and_si128(b_v,cfg->fm_m11);\
        \
        cnts_v = _mm_add_epi16(cnts_v, a_v);\
        cnts_v = _mm_add_epi16(cnts_v, b_v);\
  } while (0)



/* Macro for SSE operations that turns a vector of 4-bit character values into
 * 2 vectors representing matches. Each byte in the input vector consists of
 * a left half (4 bits) and a right half (4 bits). The 16 left-halves produce
 * one vector, which contains all-1s for bytes in which the left half matches
 * the c_v character (and 0s if it doesn't), while the 16 right-halves produce
 * the other vector, again with each byte either all-1s or all-0s.
 *
 * The expectation is that FM_COUNT_4BIT will be called after this, to
 * turn these binary values into sums over a series of vectors. The macros
 * are split up to allow one end or other to be trimmed in the case that
 * counting is not expected to include the full vector.
 *
 * srli(4)+and() : capture the left 4-bit value   (need the mask because 16-bit shift leaves garbage in left-4-bit chunks)
 * and()         : capture the right 4-bit value
 *
 * cmpeq()x2     : test if both left and right == c.  For each, if ==c , value = 11111111 (-1)
 */
#define FM_MATCH_4BIT(in_v, c_v, out1_v, out2_v) do {\
    out1_v    = _mm_srli_epi16(in_v, 4);\
    out2_v    = _mm_and_si128(in_v, cfg->fm_m0f);\
    out1_v    = _mm_and_si128(out1_v, cfg->fm_m0f);\
    \
    out1_v    = _mm_cmpeq_epi8(out1_v, c_v);\
    out2_v    = _mm_cmpeq_epi8(out2_v, c_v);\
  } while (0)


/* Macro for SSE operations that turns a vector of 4-bit character values into
 * 2 vectors representing matches. Each byte in the input vector consists of
 * a left half (4 bits) and a right half (4 bits). The 16 left-halves produce
 * one vector, which contains all-1s for bytes in which the left half is less than
 * the c_v character (and 0s if it doesn't), while the 16 right-halves produce
 * the other vector, again with each byte either all-1s or all-0s.
 *
 * The expectation is that FM_COUNT_4BIT will be called after this, to
 * turn these binary values into sums over a series of vectors. The macros
 * are split up to allow one end or other to be trimmed in the case that
 * counting is not expected to include the full vector.
 *
 * srli(4)+and() : capture the left 4-bit value   (need the mask because 16-bit shift leaves garbage in left-4-bit chunks)
 * and()         : capture the right 4-bit value
 *
 * cmplt()x2     : test if both left and right < c.  For each, if <c , value = 11111111 (-1)
 */
#define FM_LT_4BIT(in_v, c_v, out1_v, out2_v) do {\
    out1_v    = _mm_srli_epi16(in_v, 4);\
    out2_v    = _mm_and_si128(in_v, cfg->fm_m0f);\
    out1_v    = _mm_and_si128(out1_v, cfg->fm_m0f);\
    \
    out1_v    = _mm_cmplt_epi8(out1_v, c_v);\
    out2_v    = _mm_cmplt_epi8(out2_v, c_v);\
  } while (0)



/* Macro for SSE operations to add occurrence counts to the tally vector counts_v,
 * in the 4-bits-per-character case
 *
 * The expectation is that in[12]_v will contain bytes that are either
 *   00000000  =  0
 *  or
 *   11111111  = -1
 * so subtracting the value of the byte is the same as adding 0 or 1.
 */
#define FM_COUNT_4BIT(in1_v, in2_v, cnts_v) do {\
    cnts_v = _mm_subs_epi8(cnts_v, in1_v);\
    cnts_v = _mm_subs_epi8(cnts_v, in2_v);\
  } while (0)



extern int fm_initConfig      (FM_CFG *cfg, ESL_GETOPTS *go);
extern int fm_destroyConfig   (FM_CFG *cfg );
extern int fm_getOccCount     (const FM_DATA *fm, FM_CFG *cfg, int pos, uint8_t c);
extern int fm_getOccCountLT   (const FM_DATA *fm, FM_CFG *cfg, int pos, uint8_t c, uint32_t *cnteq, uint32_t *cntlt);

#endif /*P7_FM_SSE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
