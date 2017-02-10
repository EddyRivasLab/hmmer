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

/* In calculating Q, the number of vectors (of width V bytes) we need
 * in a row for a model of length M, we have to make sure there's at
 * least 2, or a striped implementation fails.
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




// Files to convert striped data from one vector size to another
extern void p7_restripe_byte(char *source, char *dest, int length, int source_vector_length, int dest_vector_length);
extern void p7_restripe_short(int16_t *source, int16_t *dest, int length, int source_vector_length, int dest_vector_length);
extern void p7_restripe_float(float *source, float *dest, int length, int source_vector_length, int dest_vector_length);

extern void p7_simdvec_Init(void);

#endif  /*p7SIMDVEC_INCLUDED*/

