/* Configuration constants and initialization of our SIMD vector environment.
 * 
 * This header file is designed to be included by non-vectorized
 * code. For example, P7_SPARSEMASK code that is called by the
 * vectorized fwd-bck local decoder.
 *
 * See simdvec.md for additional notes.
 */
#ifndef p7SIMDVEC_INCLUDED
#define p7SIMDVEC_INCLUDED
#include "easel.h"

#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
/* Actual vector widths (in # of elements) that our SSV, VF, and F/B
 * filters use. These are runtime queries, using cpu dispatching.
 *  
 * Striping and unstriping requires knowing V, so even nonvector code
 * may need these if it's accessing striped data.
 * 
 */
#define P7_V_SSV   (p7_simdvec_Width())     // (Actual) width of SSV score vectors (# of int8)
#define P7_V_VF    (p7_simdvec_Width()/2)   //            ... of VF score vectors (# of int16)
#define P7_V_FB    (p7_simdvec_Width()/4)   //            ... of FB score vectors (# of float)


/* Maximum vector widths (in # of elements).  When we allocate space
 * for vector arrays, we do allocations that will work for *any*
 * vector ISA. For that, we only need to know the maximum vector sizes
 * that are possible, to be able to calculate Q, the number of vectors
 * per DP row.
 *
 * We also need to align memory that's suitable for any vector ISA, so
 * that also depends on the maximum vector size. It's the same number
 * as p7_VMAX_SSV, but we call it p7_VALIGN for "clarity".
 * 
 * And we also need to know the maximum width of a vector in *bytes*,
 * for allocation purposes; we call that p7_VWIDTH, but it's also the
 * same as p7_VALIGN and p7_VMAX_SSV.
 *
 * If we add a new wide (>128b) vector ISA, tests here need to be edited.
 *
 *                     ----------------- V ----------------------
 *                     bits   int8 (SSV)  int16 (VF)  float (F/B)
 *                     ----   ----------  ----------  -----------
 * ARM NEON            128       16           8            4
 * Power Altivec/VMX   128       16           8            4
 * x86 SSE             128       16           8            4
 * x86 AVX             256       32          16            8
 * x86 AVX-512         512       64          32           16
 */
#if     defined(eslENABLE_AVX512)  // If we support up to 512b vector ISAs, then:
#define p7_VALIGN     64           //   memory alignment, in bytes.
#define p7_VWIDTH     64           //   max vector width, in bytes.
#define p7_VMAX_SSV   64           //                 ... in int8's
#define p7_VMAX_VF    32           //                 ... in int16's
#define p7_VMAX_FB    16           //                 ... in floats 
#elif   defined(eslENABLE_AVX)     // Or for supporting up to 256b vector ISAs:
#define p7_VALIGN     32
#define p7_VWIDTH     32
#define p7_VMAX_SSV   32
#define p7_VMAX_VF    16
#define p7_VMAX_FB     8          
#else                              // Or for supporting up to 128b vector ISAs (SSE,VMX,NEON) :
#define p7_VALIGN     16
#define p7_VWIDTH     16
#define p7_VMAX_SSV   16
#define p7_VMAX_VF     8
#define p7_VMAX_FB     4          
#endif


/* In vector ISA-specific code (*_sse.c, for example), rather than
 * hardcoding vector width in bytes, int16's, or floats as 16, 8, 4,
 * we write
 *       p7_VWIDTH_SSE 
 *       p7_VWIDTH_SSE / sizeof(int16_t)
 *       p7_VWIDTH_SSE / sizeof(float)
 * to make it more obvious what the number means.
 */
#define p7_VWIDTH_SSE    16
#define p7_VWIDTH_AVX    32
#define p7_VWIDTH_AVX512 64
#define p7_VWIDTH_NEON   16
#define p7_VWIDTH_VMX    16


/* P7_Q(M,V): calculate length of a DP matrix row, in vectors.
 * 
 * Q must be at least 2, or striped DP implementations fail.
 * Allocations account for slop caused by roundoff to integer number
 * of vectors: i.e. an allocation might be for, in bytes:
 *    P7_Q(M, p7_VMAX_SSV) * p7_VWIDTH
 *    P7_Q(M, p7_VMAX_VF)  * p7_VWIDTH
 *    P7_Q(M, p7_VMAX_FB)  * p7_VWIDTH
 * and this allocation must also be aligned on a p7_VALIGN byte boundary.
 */
#define P7_Q(M,V)  ( ESL_MAX(2, ((((M)-1) / (V)) + 1)) )



/* These macros translate between striped and unstriped coordinate
 * systems. They're conveniences, intended for code that's not
 * performance-critical. See simdvec.md notes on striped coordinates.
 */
#define P7_Q_FROM_K(k,Q)        ( ((k)-1)%(Q) )
#define P7_Z_FROM_K(k,Q)        ( ((k)-1)/(Q) )
#define P7_Y_FROM_K(k,Q,V)      ( (V) * (((k)-1)%(Q)) + ((k)-1)/(Q) )
#define P7_K_FROM_QZ(q,z,Q)     ( (z)*(Q)+((q)+1) )
#define P7_Y_FROM_QZ(q,z,V)     ( (q)*(V)+(z) )
#define P7_K_FROM_Y(y,Q,V)      ( ((y)%(V))*(Q) + (y)/(V) + 1 )
#define P7_Q_FROM_Y(y,V)        ( (y)/(V) )
#define P7_Z_FROM_Y(y,V)        ( (y)%(V) )

#define p7O_EXTRA_SB 17    /* see ssvfilter.c for explanation */


// Files to convert striped data from one vector size to another
extern void p7_restripe_byte(char *source, char *dest, int length, int source_vector_length, int dest_vector_length);
extern void p7_restripe_short(int16_t *source, int16_t *dest, int length, int source_vector_length, int dest_vector_length);
extern void p7_restripe_float(float *source, float *dest, int length, int source_vector_length, int dest_vector_length);

extern void p7_simdvec_Init(void);
extern int  p7_simdvec_Width(void);
#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif  /*p7SIMDVEC_INCLUDED*/

