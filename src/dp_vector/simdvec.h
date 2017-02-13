/* Configuration constants and initialization of our SIMD vector environment.
 * 
 * This header file may be included by non-vector code: for example,
 * P7_SPARSEMASK code that is called by the vectorized fwd-bck local
 * decoder.
 * 
 * See simdvec.md for additional notes.
 */
#ifndef p7SIMDVEC_INCLUDED
#define p7SIMDVEC_INCLUDED

/* P7_NV*() macros
 * These are for calculating Q, the number of vectors (of width V bytes)
 * that we need in a DP or profile row for a model of length M. 
 * 
 * Q must be at least 2, or a striped implementation fails.
 * 
 * Data structures hold Q and V information, set by any SIMD vector
 * routine that stores or modifies data. We assume that all SIMD
 * vector routines called in the same process use the same SIMD ISA.
 *
 * Q and V are sufficient to enable (slow) serial access to striped
 * memory in a SIMD ISA-independent fashion, so data structure
 * debugging code is typically ISA-independent.
 * 
 *             V (bytes)
 *             ---------
 * SSE         : 16
 * AVX         : 32
 * AVX-512     : 64
 * NEON        : 16
 * Altivec/VMX : 16
 */
#define P7_NVB(M,V)   ( ESL_MAX(2, ((((M)-1) /  V)    + 1)))  
#define P7_NVW(M,V)   ( ESL_MAX(2, ((((M)-1) / (V/2)) + 1)))  
#define P7_NVF(M,V)   ( ESL_MAX(2, ((((M)-1) / (V/4)) + 1)))  

/* p7_MAXV* constants
 * 
 * Data structures are allocated in an ISA-independent fashion, by
 * allocating for the maximum size required by any runtime-available
 * ISA. For one row, that's something like: 
 *      P7_NVB(M,p7_MAXVB) * p7_MAXVB bytes.
 * or similar for words, floats.
 *
 * Vector memory is aligned on p7_VALIGN byte boundaries.
 */
#if     defined(eslENABLE_AVX512)  // If we support up to 512b vector ISAs, then:
#define p7_MAXVB   64              // max vector width in bytes
#define p7_MAXVW   32              //                ... in words (int16)
#define p7_MAXVF   16              //                ... in floats (or int32)
#define p7_VALIGN  64              // memory alignment (same as MAXVB, less error prone)
#elseif defined(eslENABLE_AVX)     // Or for up to 256b vector ISAs:
#define p7_VALIGN  32
#define p7_MAXVB   32
#define p7_MAXVW   16
#define p7_MAXVF    8
#else                              // Or for 128b vector ISAs (SSE,VMX,NEON) :
#define p7_VALIGN  16
#define p7_MAXVB   16
#define p7_MAXVW    8
#define p7_MAXVF    4
#endif



// Files to convert striped data from one vector size to another
extern void p7_restripe_byte(char *source, char *dest, int length, int source_vector_length, int dest_vector_length);
extern void p7_restripe_short(int16_t *source, int16_t *dest, int length, int source_vector_length, int dest_vector_length);
extern void p7_restripe_float(float *source, float *dest, int length, int source_vector_length, int dest_vector_length);

extern void p7_simdvec_Init(void);

#endif  /*p7SIMDVEC_INCLUDED*/

