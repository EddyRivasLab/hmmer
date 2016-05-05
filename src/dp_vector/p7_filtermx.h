#ifndef p7FILTERMX_INCLUDED
#define p7FILTERMX_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */
#ifdef p7_use_AVX
	#include <immintrin.h>  /* AVX2 */
#endif
#ifdef p7_use_AVX_512
	#include <immintrin.h>  /* AVX-512 */
#endif

enum p7f_mxtype_e { p7F_NONE = 0, p7F_SSVFILTER = 1, p7F_MSVFILTER = 2, p7F_VITFILTER = 3 };
                                
typedef struct p7_filtermx_s {  /* MSV needs P7_NVB(M); VF needs 3*P7_NVW(M) __m128i vectors. */
  int      M;			/* current profile size: determines width of <dp> row    */

//#ifdef p7_use_SSE
  __m128i *dp;			/* aligned, one row of DP memory: >= 3*P7_NVW(M) vectors */

  void    *dp_mem;		/* unaligned raw memory, where we allocate    */
  int      allocM;		/* <dp_mem> is allocated to hold up to M=allocM */
//#endif

#ifdef p7_use_AVX
	__m256i *dp_AVX;   // add separate versions of each of these to support running the different
	void *dp_mem_AVX;  // vector ISAs simultaneously for testing.
	int allocM_AVX;
#endif

#ifdef p7_use_AVX_512
	__m512i *dp_AVX_512;   // add separate versions of each of these to support running the different
	void *dp_mem_AVX_512;  // vector ISAs simultaneously for testing.
	int allocM_AVX_512;
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
#define MMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_M])
#define DMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_D])
#define IMXf(q)   (dp[(q) * p7F_NSCELLS + p7F_I])

#ifdef p7_use_AVX
/* Viterbi filter uses these macros to improve clarity of accesses */
#define MMX_AVXf(q)   (dp_AVX[(q) * p7F_NSCELLS + p7F_M])
#define DMX_AVXf(q)   (dp_AVX[(q) * p7F_NSCELLS + p7F_D])
#define IMX_AVXf(q)   (dp_AVX[(q) * p7F_NSCELLS + p7F_I])
#endif

#ifdef p7_use_AVX_512
/* Viterbi filter uses these macros to improve clarity of accesses */
#define MMX_AVX_512f(q)   (dp_AVX_512[(q) * p7F_NSCELLS + p7F_M])
#define DMX_AVX_512f(q)   (dp_AVX_512[(q) * p7F_NSCELLS + p7F_D])
#define IMX_AVX_512f(q)   (dp_AVX_512[(q) * p7F_NSCELLS + p7F_I])
#endif

extern P7_FILTERMX *p7_filtermx_Create(int allocM);
extern int          p7_filtermx_GrowTo(P7_FILTERMX *fx, int allocM);
extern size_t       p7_filtermx_Sizeof(const P7_FILTERMX *fx);
extern size_t       p7_filtermx_MinSizeof(int M);
extern int          p7_filtermx_Reuse  (P7_FILTERMX *fx);
extern void         p7_filtermx_Destroy(P7_FILTERMX *fx);

extern int p7_filtermx_SetDumpMode(P7_FILTERMX *fx, FILE *dfp, int truefalse);
#ifdef p7_DEBUGGING
extern int p7_filtermx_DumpMFRow(const P7_FILTERMX *fx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
extern int p7_filtermx_DumpVFRow(const P7_FILTERMX *fx, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
#endif

#endif /*p7FILTERMX_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/


