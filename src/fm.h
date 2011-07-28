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

enum fm_status_e {
  fmOK   = 0,
  fmFAIL = 1,
  fmMEM    = 2
};


enum fm_compression_e {
  fmUNPACKED    = 0,
  fmPACKED      = 1,
  fmSTRIPED     = 2
};

enum fm_rle_e {
  fmRLE         = 0,
  fmNORLE       = 1
};


enum fm_sse_e {
  fmSSE     = 0,
  fmDUMMY   = 1
};

#define FM_ALLOC(p, size) do {\
     if (((p) = malloc(size)) == NULL) {\
       status = fmMEM;\
       goto ERROR;\
     }} while (0)

#define FM_RALLOC(p, newsize) do {\
	 void* tmp;\
     if ((p) == NULL) { (tmp) = malloc(newsize);         }\
     else             { (tmp) = realloc((p), (newsize)); }\
     if ((tmp) != NULL) (p) = (tmp);\
     else {\
       status = fmMEM;\
       goto ERROR;\
     }} while (0)

#define FM_FAIL( str, ...) do {\
     fprintf(stderr, str, __VA_ARGS__);\
     exit(fmFAIL); }\
     while (0)

#define FM_MAX_LINE 256

#define PRINTOUT( str, ...) do {\
	 if (out != NULL) \
        fprintf( out, str, __VA_ARGS__);\
     }\
     while (0)


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




/* Extract the sum of all counts in a 16x8-bit element
 *  the _mm_sad_epu8  accumulates 8-bit counts into 16-bit counts:
 *			left 8 counts (64-bits) accumulate in counts_v[0],
 *			right 8 counts in counts_v[4]  (the other 6 16-bit ints are 0)
 *	the _mm_shuffle_epi32  flips the 4th int into the 0th slot
 */
#define FM_EXTRACT_8BIT_COUNTS( in_v, mid_v, out_v  ) do {\
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
/*  16-bit counts
  #define FM_EXTRACT_COUNT(in_v, buffer_v, out_s ) do {\
        buffer_v  = _mm_srli_epi32(in_v,16);\
		in_v      = _mm_add_epi16(in_v, buffer_v);\
		buffer_v  = _mm_shuffle_epi32(in_v, _MM_SHUFFLE(2, 3, 0, 1));\
		in_v      = _mm_add_epi16(in_v, buffer_v);\
		buffer_v  = _mm_shuffle_epi32(in_v, _MM_SHUFFLE(1, 1, 1, 2));\
		in_v      = _mm_add_epi16(in_v, buffer_v);\
		out_s     = _mm_extract_epi16(in_v, 0) ;\
		} while (0)
		*/


/*Macro for performing SSE operations to count occurrences of a character.
 * Doing left and right halves in parallel ekes out a little better performance.
 * Byte values for matching chars are 0xff, or -1, so the value is subtracted from
 * counts_v to increment counts (on a vector-parallel basis). The calling function
 * initializes counts_v to -128, so even though the max of a counts_v element
 * is +127, it's shifted down -128, so effectively max of 255
 */
#define FM_COUNTS_SSE_PACKED() do {\
		tmp_v    = _mm_and_si128(BWT_v, leftmask_v);\
		tmp2_v    = _mm_and_si128(BWT_v, rightmask_v);\
		\
		tmp_v    = _mm_cmpeq_epi8(tmp_v, leftc_tmp_v);\
		tmp2_v   = _mm_cmpeq_epi8(tmp2_v, rightc_tmp_v);\
		\
		counts_v = _mm_subs_epi8(counts_v, tmp_v);\
		counts_v = _mm_subs_epi8(counts_v, tmp2_v);\
	} while (0)





/*"population count", counting the number of bits turned on in b1
* Counts are accumulated in 16-bit blocks, so there's effectively no bound
* on the number of blocks for which sums can be captured (really, it's
* 8*64K =~ 0.5 million).  Later, I'll add up all the 16 bit numbers into a
* 16-bit field, so the real limit is 64K.
*  Similar to the binary divide-and-conquer popcount given in popcount_2() from
*    http://en.wikipedia.org/wiki/Hamming_weight,
*  but removing 2 instructions since the #s are 16 bit
*  See notes in ~/notebook/1220_fmindex/00README:  "Tue Jan 25 09:01:40 EST 2011".
*  Sadly, this takes so long to do that it overwhelms the slick gathering of bits
*  that feeds this part.
*
*  Explanation:
*  		b2 = _mm_srli_epi16(b1,1);  //put count of each 2 bits into those 2 bits
*		b2 = _mm_and_si128(b2,m55);
*		b1 = _mm_sub_epi16(b1,b2);
*
*		b2 = _mm_srli_epi16(b1,2);  //put count of each 4 bits into those 4 bits
*		b2 = _mm_and_si128(b2,m33);
*		b1 = _mm_and_si128(b1,m33);
*		b1 = _mm_add_epi16(b1,b2);
*
*		b2 = _mm_srli_epi16(b1,4);  //put count of each 8 bits into those 8 bits
*		b1 = _mm_add_epi16(b1,b2);
*		b1 = _mm_and_si128(b1,m0f);
*
*		b2 = _mm_srli_epi16(b1,8);  //put count of each 16 bits into those 16 bits
*		b1 = _mm_add_epi16(b1,b2);
*		b1 = _mm_and_si128(b1,m00ff);
*/
#define FM_POP_COUNT_16BIT(in, out) do {\
		b2 = _mm_srli_epi16(in,1);\
		b2 = _mm_and_si128(b2,m55);\
		b1 = _mm_sub_epi16(in,b2);\
		b2 = _mm_srli_epi16(b1,2);\
		b2 = _mm_and_si128(b2,m33);\
		b1 = _mm_and_si128(b1,m33);\
		b1 = _mm_add_epi16(b1,b2);\
		b2 = _mm_srli_epi16(b1,4);\
		b1 = _mm_add_epi16(b1,b2);\
		b1 = _mm_and_si128(b1,m0f);\
		b2 = _mm_srli_epi16(b1,8);\
		b1 = _mm_add_epi16(b1,b2);\
		out = _mm_and_si128(b1,m00ff);\
	} while (0)

#define FM_POP_COUNT_8BIT(in, out) do {\
		b2 = _mm_srli_epi16(in,1);\
		b2 = _mm_and_si128(b2,m55);\
		b1 = _mm_sub_epi16(in,b2);\
		b2 = _mm_srli_epi16(b1,2);\
		b2 = _mm_and_si128(b2,m33);\
		b1 = _mm_and_si128(b1,m33);\
		b1 = _mm_add_epi16(b1,b2);\
		b2 = _mm_srli_epi16(b1,4);\
		b1 = _mm_add_epi16(b1,b2);\
		out = _mm_and_si128(b1,m0f);\
	} while (0)

#define FM_POP_COUNT_4BIT(in, out) do {\
		b2 = _mm_srli_epi16(in,1);\
		b2 = _mm_and_si128(b2,m55);\
		b1 = _mm_sub_epi16(b1,b2);\
		b2 = _mm_srli_epi16(b1,2);\
		b2 = _mm_and_si128(b2,m33);\
		b1 = _mm_and_si128(b1,m33);\
		out = _mm_add_epi16(b1,b2);\
	} while (0)



#define FM_POP_COUNT(bits, counts_in, counts_out) do {\
		FM_POP_COUNT_8BIT(bits, bits);\
		counts_out = _mm_add_epi16(bits,counts_in);\
	} while (0)




#define FM_EXTRACT16(target, source, index) do {\
		switch (index) {\
			case 0 :\
				target = _mm_extract_epi16(source, 0);\
				break;\
			case 1 :\
				target = _mm_extract_epi16(source, 1);\
				break;\
			case 2 :\
				target = _mm_extract_epi16(source, 2);\
				break;\
			case 3 :\
				target = _mm_extract_epi16(source, 3);\
				break;\
			case 4 :\
				target = _mm_extract_epi16(source, 4);\
				break;\
			case 5 :\
				target = _mm_extract_epi16(source, 5);\
				break;\
			case 6 :\
				target = _mm_extract_epi16(source, 6);\
				break;\
			case 7 :\
				target = _mm_extract_epi16(source, 7);\
				break;\
		}\
	} while (0)


enum fm_alphabettypes_e {
  fm_DNA        = 0,  //acgt$,  3 bit
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
  int freq_SA; //frequency with which SA is sampled
  int freq_cnt_sb; //frequency with which full cumulative counts are captured
  int freq_cnt_b; //frequency with which intermittent counts are captured
  int SA_shift;
  int cnt_shift_sb;
  int cnt_shift_b;
} FM_METADATA;

typedef struct fm_data_s {
  uint8_t  *BWT;
  uint32_t  *SA; // sampled suffix array
  uint32_t   *C; //the first position of each letter of the alphabet if all of T is sorted.
  uint32_t  *occCnts_sb;
  uint16_t  *occCnts_b;
} FM_DATA;
