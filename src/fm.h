
#include "easel.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include <unistd.h> //getopt
#include <getopt.h>
#include <assert.h>


// fm.c
extern void createAlphabet (int alph_type, char **alph, char **inv_alph, int *alph_size, int *alph_bits);
void print_m128 (__m128i in);
void print_m128_rev (__m128i in);
char * commaprint(unsigned long n);
void memset_16aligned(void *space, size_t nbytes);


/*****************************************************************
 * Macros implementing memory allocation conventions
 *****************************************************************/
/* FM_ALLOC(), FM_RALLOC():
 *
 * Allocation and reallocation wrappers.
 * Both require <int status> in scope, and <ERROR:> goto target.
 *
 * Essentially lifted from Easel library
 */


#define FM_MAX_LINE 256

/*
#define PRINTOUT( str, ...) do {\
	 if (out != NULL) \
        fprintf( out, str, __VA_ARGS__);\
     }\
     while (0)
*/

//used to convert from a byte array to an __m128i
typedef union {
        uint8_t bytes[16];
        __m128i m128;
        } byte_m128;


/*****************************************************************
 * BWT: Macros/structs used in building and searching with BWT / FM-index
 *****************************************************************/

/* Structure the 2D occ array into a single array.  "type" is either b or sb.
 * Note that one extra count value is required by RLE, one 4-byte int for
 * each superblock count vector, and one 2-byte short for each block count
 * vector. This is small overhead, even for a small alphabet like dna.
 */
#define FM_OCC_CNT( type, i, c)  ( occCnts_##type[(meta->alph_size)*(i) + (c)]) // good for non-RLE




/* Gather the sum of all counts in a 16x8-bit element into a single 16-bit
 *  element of the register (the 0th element)
 *
 *  the _mm_sad_epu8  accumulates 8-bit counts into 16-bit counts:
 *			left 8 counts (64-bits) accumulate in counts_v[0],
 *			right 8 counts in counts_v[4]  (the other 6 16-bit ints are 0)
 *	the _mm_shuffle_epi32  flips the 4th int into the 0th slot
 */
#define FM_GATHER_8BIT_COUNTS( in_v, mid_v, out_v  ) do {\
		mid_v = _mm_sad_epu8 (in_v, zeros_v);\
		tmp_v = _mm_shuffle_epi32(mid_v, _MM_SHUFFLE(1, 1, 1, 2));\
		out_v = _mm_add_epi16(mid_v, tmp_v);\
	} while (0)


/* Extract the sum of counts from either a collection of 8 16-bit
 * numbers in a m128 element, or 16 8-bit numbers
 */
#define FM_EXTRACT_COUNT(in_v, buffer_v, out_s ) do {\
				FM_EXTRACT_8BIT_COUNTS(in_v,buffer_v,buffer_v);\
				out_s =  _mm_extract_epi16(buffer_v, 0) ;\
		} while (0)




/*Macro for performing SSE operations to count occurrences of a character in a
 * 4-char-per-byte packed BWT.
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
 * then add up the 2-bit values:
 * srli(4)+add()         : left 4 bits shifted right, added to right 4 bits
 *
 * srli(2)+and(00000011) : left 2 bits (value 0..2) shifted right, masked, so no other bits active
 * and(00000011)         : right 2 bits (value 0..2) masked so no other bits active
 *
 * final 2 add()s        : tack current counts on to already-tabulated counts.
 */
#define FM_COUNTS_SSE_4PACKED(in_v, c_v, tmp_v, tmp2_v, cnts_v) do {\
		tmp_v = _mm_xor_si128(in_v, c_v);\
        \
        tmp2_v = _mm_and_si128(tmp_v, m01);\
        tmp_v  = _mm_srli_epi16 (tmp_v, 1);\
        tmp_v  = _mm_and_si128(tmp_v, m01);\
        tmp_v  = _mm_or_si128(tmp_v, tmp2_v);\
        \
        tmp_v  = _mm_subs_epi8(m01,tmp_v);\
        \
        tmp2_v = _mm_srli_epi16(tmp_v, 4);\
        tmp_v  = _mm_add_epi16(tmp_v, tmp2_v);\
        \
        tmp2_v = _mm_srli_epi16(tmp_v, 2);\
        tmp_v  = _mm_and_si128(tmp_v,m11);\
        tmp2_v = _mm_and_si128(tmp2_v,m11);\
        \
        cnts_v = _mm_add_epi16(cnts_v, tmp_v);\
        cnts_v = _mm_add_epi16(cnts_v, tmp2_v);\
	} while (0)



/* Macro for performing SSE operations to count occurrences of a character in
 * a 2-char-per-byte packed BWT.
 *
 * tmp_v and tmp2_v are used as temporary vectors throughout, and hold meaningless
 * values at the end.
 *
 * and()         : capture the right 4-bit value
 * srli(4)+and() : capture the left 4-bit value   (need the mask because 16-bit shift leave garbage in left-4-bit chunks)
 *
 * cmpeq()x2     : test if both left and right == c.  If so, value = 11111111 (-1)
 *
 * subs()x2      : subtracting -1 == adding 1
 */
#define FM_COUNTS_SSE_2PACKED(in_v, c_v, tmp_v, tmp2_v, counts_v) do {\
		tmp_v    = _mm_srli_epi16(BWT_v, 4);\
        tmp2_v   = _mm_and_si128(BWT_v, c_v);\
		tmp_v    = _mm_and_si128(tmp_v, c_v);\
		\
		tmp_v    = _mm_cmpeq_epi8(tmp_v, c_v);\
		tmp2_v   = _mm_cmpeq_epi8(tmp2_v, c_v);\
		\
		counts_v = _mm_subs_epi8(counts_v, tmp_v);\
		counts_v = _mm_subs_epi8(counts_v, tmp2_v);\
	} while (0)


enum fm_alphabettypes_e {
  fm_DNA        = 0,  //acgt,  2 bit
  fm_DNA_full   = 1,  //includes ambiguity codes, 4 bit
  fm_AMINO      = 2,  // 5 bit
};

typedef struct fm_interval_s {
  int   lower;
  int   upper;
} FM_INTERVAL;

typedef struct fm_hit_s {
  int   start;
//  int   length;
} FM_HIT;


typedef struct fm_metadata_s {
	//TODO: these don't need to be ints - uint8_t ?
  int alph_type;
  int alph_size;
  int charBits;
  int N; //length of text
  int L; //bytes used to store BWT/T
  int freq_SA; //frequency with which SA is sampled
  int freq_cnt_sb; //frequency with which full cumulative counts are captured
  int freq_cnt_b; //frequency with which intermittent counts are captured
  int SA_shift;
  int cnt_shift_sb;
  int cnt_shift_b;
  int term_loc; // location in the BWT at which the '$' char is found (replaced in the sequence with 'a')
} FM_METADATA;

typedef struct fm_data_s {
  uint8_t  *T;  //text corresponding to the BWT
  uint8_t  *BWT;
  //byte_m128  *BWT; //BWT, each element of this array is either an array of 16 bytes or a single m128 element
  uint32_t  *SA; // sampled suffix array
  int32_t   *C; //the first position of each letter of the alphabet if all of T is sorted.  (signed, as I use that to keep tract of presence/absence)
  uint32_t  *occCnts_sb;
  uint16_t  *occCnts_b;
} FM_DATA;
