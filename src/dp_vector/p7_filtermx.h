#ifndef p7FILTERMX_INCLUDED
#define p7FILTERMX_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#if p7_CPU_ARCH == intel 
#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#ifdef HAVE_AVX2
	#include <immintrin.h>  /* AVX2 */
#endif
#ifdef HAVE_AVX512
	#include <immintrin.h>  /* AVX-512 */
#endif
#endif /* intel arch */

#if p7_CPU_ARCH == arm || p7_CPU_ARCH == arm64
#include <arm_neon.h>
#ifdef HAVE_NEON
	#include "esl_neon.h"
#endif
#endif /* arm/arm64 */

#include "hardware/hardware.h"

enum p7f_mxtype_e { p7F_NONE = 0, p7F_SSVFILTER = 1, p7F_MSVFILTER = 2, p7F_VITFILTER = 3 };
                                
typedef struct p7_filtermx_s {  /* MSV needs P7_NVB(M); VF needs 3*P7_NVW(M) __m128i vectors. */
  int      M;			/* current profile size: determines width of <dp> row    */
 SIMD_TYPE simd;   //Which SIMD ISA are we using?
#ifdef HAVE_SSE2
  __m128i *dp;			/* aligned, one row of DP memory: >= 3*P7_NVW(M) vectors */

  void    *dp_mem;		/* unaligned raw memory, where we allocate    */
  int      allocM;		/* <dp_mem> is allocated to hold up to M=allocM */
#endif

#ifdef HAVE_AVX2
	__m256i *dp_AVX;   // add separate versions of each of these to support running the different
	void *dp_mem_AVX;  // vector ISAs simultaneously for testing.
	int allocM_AVX;
#endif

#ifdef HAVE_AVX512
	__m512i *dp_AVX_512;   // add separate versions of each of these to support running the different
	void *dp_mem_AVX_512;  // vector ISAs simultaneously for testing.
	int allocM_AVX_512;
#endif

#ifdef HAVE_NEON
  esl_neon_128i_t *dp;
  void *dp_mem;
  int allocM;
#endif

  enum p7f_mxtype_e type;	/* p7F_NONE | p7F_SSVFILTER | p7F_MSVFILTER | p7F_VITFILTER */

#ifdef p7_DEBUGGING
  int      do_dumping;		/* TRUE if we're dumping each row for diagnostics */
  FILE    *dfp;			/* where we're dumping it (or NULL, if debugging=FALSE) */
#endif
} P7_FILTERMX;


/* Viterbi filter has three cells. MSV/SSV only have match states. */
enum p7f_scells_e { p7F_M = 0, p7F_D = 1, p7F_I = 2 };
#define p7F_NSCELLS 3

/* Viterbi filter uses these macros to improve clarity of accesses */
#ifdef HAVE_SSE2
#define MMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_M])
#define DMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_D])
#define IMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_I])
#endif

/* Viterbi filter uses these macros to improve clarity of accesses */
#ifdef HAVE_NEON
#define MMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_M])
#define DMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_D])
#define IMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_I])
#endif
 
#ifdef HAVE_AVX2
/* Viterbi filter uses these macros to improve clarity of accesses */
#define MMX_AVXf(q)   (dp_AVX[(q) * p7F_NSCELLS + p7F_M])
#define DMX_AVXf(q)   (dp_AVX[(q) * p7F_NSCELLS + p7F_D])
#define IMX_AVXf(q)   (dp_AVX[(q) * p7F_NSCELLS + p7F_I])
#endif

#ifdef HAVE_AVX512
/* Viterbi filter uses these macros to improve clarity of accesses */
#define MMX_AVX_512f(q)   (dp_AVX_512[(q) * p7F_NSCELLS + p7F_M])
#define DMX_AVX_512f(q)   (dp_AVX_512[(q) * p7F_NSCELLS + p7F_D])
#define IMX_AVX_512f(q)   (dp_AVX_512[(q) * p7F_NSCELLS + p7F_I])
#endif

extern P7_FILTERMX *p7_filtermx_Create(int allocM, SIMD_TYPE simd);
extern int          p7_filtermx_GrowTo(P7_FILTERMX *fx, int allocM);
extern size_t       p7_filtermx_Sizeof(const P7_FILTERMX *fx);
extern size_t       p7_filtermx_MinSizeof(int M, SIMD_TYPE simd);
extern int          p7_filtermx_Reuse  (P7_FILTERMX *fx);
extern void         p7_filtermx_Destroy(P7_FILTERMX *fx);

extern P7_FILTERMX *p7_filtermx_Create_sse(int allocM);
extern int          p7_filtermx_GrowTo_sse(P7_FILTERMX *fx, int allocM);
extern size_t       p7_filtermx_Sizeof_sse(const P7_FILTERMX *fx);
extern size_t       p7_filtermx_MinSizeof_sse(int M);
extern void         p7_filtermx_Destroy_sse(P7_FILTERMX *fx);
 
extern P7_FILTERMX *p7_filtermx_Create_avx(int allocM);
extern int          p7_filtermx_GrowTo_avx(P7_FILTERMX *fx, int allocM);
extern size_t       p7_filtermx_Sizeof_avx(const P7_FILTERMX *fx);
extern size_t       p7_filtermx_MinSizeof_avx(int M);
extern void         p7_filtermx_Destroy_avx(P7_FILTERMX *fx);

extern P7_FILTERMX *p7_filtermx_Create_avx512(int allocM);
extern int          p7_filtermx_GrowTo_avx512(P7_FILTERMX *fx, int allocM);
extern size_t       p7_filtermx_Sizeof_avx512(const P7_FILTERMX *fx);
extern size_t       p7_filtermx_MinSizeof_avx512(int M);
extern void         p7_filtermx_Destroy_avx512(P7_FILTERMX *fx);

extern P7_FILTERMX *p7_filtermx_Create_neon(int allocM);
extern int          p7_filtermx_GrowTo_neon(P7_FILTERMX *fx, int allocM);
extern size_t       p7_filtermx_Sizeof_neon(const P7_FILTERMX *fx);
extern size_t       p7_filtermx_MinSizeof_neon(int M);
extern void         p7_filtermx_Destroy_neon(P7_FILTERMX *fx);

extern P7_FILTERMX *p7_filtermx_Create_neon64(int allocM);
extern int          p7_filtermx_GrowTo_neon64(P7_FILTERMX *fx, int allocM);
extern size_t       p7_filtermx_Sizeof_neon64(const P7_FILTERMX *fx);
extern size_t       p7_filtermx_MinSizeof_neon64(int M);
extern void         p7_filtermx_Destroy_neon64(P7_FILTERMX *fx);

extern int p7_filtermx_SetDumpMode(P7_FILTERMX *fx, FILE *dfp, int truefalse);
#ifdef p7_DEBUGGING
extern int p7_filtermx_DumpMFRow(const P7_FILTERMX *fx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
extern int p7_filtermx_DumpMFRow_sse(const P7_FILTERMX *fx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
extern int p7_filtermx_DumpMFRow_avx(const P7_FILTERMX *fx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
extern int p7_filtermx_DumpMFRow_avx512(const P7_FILTERMX *fx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
extern int p7_filtermx_DumpMFRow_neon(const P7_FILTERMX *fx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
extern int p7_filtermx_DumpVFRow(const P7_FILTERMX *fx, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
extern int p7_filtermx_DumpVFRow_sse(const P7_FILTERMX *fx, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
extern int p7_filtermx_DumpVFRow_avx(const P7_FILTERMX *fx, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
extern int p7_filtermx_DumpVFRow_avx512(const P7_FILTERMX *fx, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
extern int p7_filtermx_DumpVFRow_neon(const P7_FILTERMX *fx, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
#endif

#endif /*p7FILTERMX_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/


