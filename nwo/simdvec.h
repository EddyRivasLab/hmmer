/* Configuration constants and initialization of our SIMD vector environment.
 * 
 * This header file is designed to be included by non-vectorized code;
 * for example, H4_SPARSEMASK code that gets called by the vectorized
 * fwd-bck local decoder.
 *
 * See simdvec.md for additional notes.
 */
#ifndef h4SIMDVEC_INCLUDED
#define h4SIMDVEC_INCLUDED
#include "h4_config.h"

#include "easel.h"

/* Scale, offsets for reduced precision integer log-odds scores.
 * SSV filter:  8-bit (bytes, B).
 * Vit filter: 16-bit (words, W).
 */
#define h4_SCALE_B      3.0   // SSV is in 1/3 bits
#define h4_SCALE_W    500.0   // VF is in 1/500 bits
#define h4_BASE_W   12000.0   // VF scores are offset 

#define h4_2NAT_APPROX  2.8853901  // for SSV singlehit: NN/CC=0;   asymptotically they score -2nat
#define h4_3NAT_APPROX  4.3280851  // for VF multihit:   NN/CC/JJ=0,  "" -3nat


/* Actual vector widths (in # of elements) that our SSV, VF, and F/B
 * filters use. These are runtime queries, using cpu dispatching.
 * h4_simdvec_width() is a runtime call that returns # of bytes in one
 * vector, in the ISA in use.
 *  
 * Striping and unstriping requires knowing V, so even nonvector code
 * may need these if it's accessing striped data.
 * 
 */
#define H4_V_SSV   (h4_simdvec_width())     // (Actual) width of SSV score vectors (# of int8)
#define H4_V_VF    (h4_simdvec_width()/2)   //            ... of VF score vectors (# of int16)
#define H4_V_FB    (h4_simdvec_width()/4)   //            ... of FB score vectors (# of float)


/* Maximum vector widths (in # of elements).  When we allocate space
 * for vector arrays, we do allocations that will work for *any*
 * vector ISA. For that, we only need to know the maximum vector sizes
 * that are possible, to be able to calculate Q, the number of vectors
 * per DP row.
 *
 * We also need to align memory that's suitable for any vector ISA, so
 * that also depends on the maximum vector size. It's the same number
 * as h4_VMAX_SSV, but we call it h4_VALIGN for "clarity".
 * 
 * And we also need to know the maximum width of a vector in *bytes*,
 * for allocation purposes; we call that h4_VWIDTH, but it's also the
 * same as h4_VALIGN and h4_VMAX_SSV.
 *
 * If we add a new wide (>128b) vector ISA, the ifdef tests here need
 * to be edited.
 *
 *                            ---------------- V ----------------
 *                     bits   int8 (SSV)  int16 (VF)  float (F/B)
 *                     ----   ----------  ----------  -----------
 * ARM NEON            128       16           8            4
 * Power Altivec/VMX   128       16           8            4
 * x86 SSE             128       16           8            4
 * x86 AVX             256       32          16            8
 * x86 AVX-512         512       64          32           16
 */
#if     defined(eslENABLE_AVX512)  // If we support up to 512b vector ISAs, then:
#define h4_VALIGN     64           //   memory alignment, in bytes.
#define h4_VWIDTH     64           //   max vector width, in bytes.
#define h4_VMAX_SSV   64           //                 ... in int8's
#define h4_VMAX_VF    32           //                 ... in int16's
#define h4_VMAX_FB    16           //                 ... in floats 
#elif   defined(eslENABLE_AVX)     // Or for supporting up to 256b vector ISAs:
#define h4_VALIGN     32
#define h4_VWIDTH     32
#define h4_VMAX_SSV   32
#define h4_VMAX_VF    16
#define h4_VMAX_FB     8          
#else                              // Or for supporting up to 128b vector ISAs (SSE,VMX,NEON) :
#define h4_VALIGN     16
#define h4_VWIDTH     16
#define h4_VMAX_SSV   16
#define h4_VMAX_VF     8
#define h4_VMAX_FB     4          
#endif


/* In vector ISA-specific code (*_sse.c, for example), rather than
 * hardcoding vector width in bytes, int16's, or floats as 16, 8, 4,
 * we write
 *       h4_VWIDTH_SSE 
 *       h4_VWIDTH_SSE / sizeof(int16_t)
 *       h4_VWIDTH_SSE / sizeof(float)
 * to make it more obvious what the number means.
 */
#define h4_VWIDTH_SSE    16
#define h4_VWIDTH_AVX    32
#define h4_VWIDTH_AVX512 64
#define h4_VWIDTH_NEON   16
#define h4_VWIDTH_VMX    16


/* H4_Q(M,V): calculate length of a DP matrix row, in vectors.
 * 
 * Q must be at least 2, or striped DP implementations fail.
 * Allocations account for slop caused by roundoff to integer number
 * of vectors: i.e. an allocation might be for, in bytes:
 *    H4_Q(M, h4_VMAX_SSV) * h4_VWIDTH
 *    H4_Q(M, h4_VMAX_VF)  * h4_VWIDTH
 *    H4_Q(M, h4_VMAX_FB)  * h4_VWIDTH
 * and this allocation must also be aligned on a h4_VALIGN byte boundary.
 */
#define H4_Q(M,V)  ( ESL_MAX(2, ((((M)-1) / (V)) + 1)) )

/* These macros translate between striped and unstriped coordinate
 * systems. They're conveniences, intended for code that's not
 * performance-critical. See simdvec.md notes on striped coordinates.
 */
#define H4_Q_FROM_K(k,Q)        ( ((k)-1)%(Q) )
#define H4_Z_FROM_K(k,Q)        ( ((k)-1)/(Q) )
#define H4_Y_FROM_K(k,Q,V)      ( (V) * (((k)-1)%(Q)) + ((k)-1)/(Q) )
#define H4_K_FROM_QZ(q,z,Q)     ( (z)*(Q)+((q)+1) )
#define H4_Y_FROM_QZ(q,z,V)     ( (q)*(V)+(z) )
#define H4_K_FROM_Y(y,Q,V)      ( ((y)%(V))*(Q) + (y)/(V) + 1 )
#define H4_Q_FROM_Y(y,V)        ( (y)/(V) )
#define H4_Z_FROM_Y(y,V)        ( (y)%(V) )

#define h4_EXTRA_SB 17   // see ssvfilter.c for explanation


extern int     h4_simdvec_width  (void);
extern int8_t  h4_simdvec_byteify(float sc);
extern int16_t h4_simdvec_wordify(float sc);

#endif // h4SIMDVEC_INCLUDED
