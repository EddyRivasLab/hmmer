#ifndef p7OPROFILE_INCLUDED
#define p7OPROFILE_INCLUDED

#include "p7_config.h"

#if p7_CPU_ARCH == intel 
#include <xmmintrin.h>    /* SSE  */
#include <emmintrin.h>    /* SSE2 */
#ifdef HAVE_AVX2
  #include <immintrin.h>  /* AVX2 */
#endif
#ifdef HAVE_AVX512
  #include <immintrin.h>  /* AVX-512 */
#endif
#ifdef _PMMINTRIN_H_INCLUDED
#include <pmmintrin.h>   /* DENORMAL_MODE */
#endif

#endif /* p7_CPU_ARCH == intel */

#if p7_CPU_ARCH == arm || p7_CPU_ARCH == arm64
//#include <arm_neon.h>
#ifdef HAVE_NEON
	#include "esl_neon.h"
#endif
#endif /* p7_CPU_ARCH == arm/arm64 */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"

#include "base/general.h"
#include "base/p7_bg.h"
#include "base/p7_hmm.h"
#include "base/p7_profile.h"
#include "hardware/hardware.h"
#include "dp_vector/simdvec.h"

/* The OPROFILE is striped [Farrar07] and interleaved, as is the DP matrix.
 * For example, the layout of a profile for an M=14 model (xref J2/46):
 * 
 * rfv[x] : striped blocks of M emissions, starting with q=0
 *                1     11     1      1  
 *             1593   2604   371x   482* 
 *
 *          to get k given q,z:  k = zQ+q+1
 *          to get q,z given k:  q = (k-1)%Q;  z = (k-1)/Q
 *
 *          unused values (marked * above) are 0.0 odds ratio (-inf score)
 *          for x = gap, none, missing: all odds ratios are 0.0
 * 
 * tsc:  grouped in order of accession in DP for 7 transition scores;
 *       starting at q=0 for all but the three transitions to M, which
 *       are rotated by -1 and rightshifted. DD's follow separately, 
 *       starting at q=0.
 *
 *        {     1      1     1     1     1     1     1 }
 *        {  1593   x482  x482  x482  1593  1593  1593 }    
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 * 
 *        {    11      1     1     1    11    11    11 }
 *        {  2604   1593  1593  1593  2604  2604  2604 } 
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *        
 *        {    1      11    11    11    1     1     1  }
 *        {  371x   2604  2604  2604  371x  371x  371x }
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *        
 *        {    1      1     1     1     1     1     1  }
 *        {  482x   371x  371x  371x  482x  482x  482x }
 *        { [tBMk] [tMM] [tIM] [tDM] [tMD] [tMI] [tII] }
 *        
 *        {     1    11    1     1  }
 *        {  1593  2604  371x  482x }
 *        { [TDD] [TDD] [TDD] [TDD] }
 *        
 */

#define p7O_NXSTATES  4    /* special states stored: ENJC                       */
#define p7O_NXTRANS   2    /* special states all have 2 transitions: move, loop */
#define p7O_NTRANS    8    /* 7 core transitions + BMk entry                    */
enum p7o_xstates_e      { p7O_E    = 0, p7O_N    = 1,  p7O_J  = 2,  p7O_C  = 3 };
enum p7o_xtransitions_e { p7O_MOVE = 0, p7O_LOOP = 1 };
enum p7o_tsc_e          { p7O_BM   = 0, p7O_MM   = 1,  p7O_IM = 2,  p7O_DM = 3, p7O_MD   = 4, p7O_MI   = 5,  p7O_II = 6,  p7O_DD = 7 };

typedef struct p7_oprofile_s {
  /* MSVFilter uses scaled, biased uchars: 16x unsigned byte vectors                 */

  SIMD_TYPE simd;   //Which SIMD ISA are we using?

 #ifdef HAVE_SSE2 
  __m128i **rbv;                /* match scores [x][q]: rm, rm[0] are allocated      */
  __m128i **sbv;                /* match scores for ssvfilter         */
  /* Our actual vector mallocs, before we align the memory                           */
  __m128i  *rbv_mem;
  __m128i  *sbv_mem;
#endif
   #ifdef HAVE_AVX2
  __m256i **rbv_AVX;                /* match scores [x][q]: rm, rm[0] are allocated      */
  __m256i **sbv_AVX;                /* match scores for ssvfilter         */
  /* Our actual vector mallocs, before we align the memory                           */
  __m256i  *rbv_mem_AVX;
  __m256i  *sbv_mem_AVX;
  #endif

   #ifdef HAVE_AVX512
  __m512i **rbv_AVX_512;                /* match scores [x][q]: rm, rm[0] are allocated      */
  __m512i **sbv_AVX_512;                /* match scores for ssvfilter         */
  /* Our actual vector mallocs, before we align the memory                           */
  __m512i  *rbv_mem_AVX_512;
  __m512i  *sbv_mem_AVX_512;
  #endif

  #ifdef HAVE_NEON
  esl_neon_128i_t **rbv;        /* match scores [x][q]: rm, rm[0] are allocated      */
  esl_neon_128i_t **sbv;        /* match scores for ssvfilter                        */
  esl_neon_128i_t  *rbv_mem;
  esl_neon_128i_t  *sbv_mem;
  #endif

  uint8_t   tbm_b;              /* constant B->Mk cost:    scaled log 2/M(M+1)       */
  uint8_t   tec_b;              /* constant E->C  cost:    scaled log 0.5            */
  uint8_t   tjb_b;              /* constant NCJ move cost: scaled log 3/(L+3)        */
  float     scale_b;            /* typically 3 / log2: scores scale to 1/3 bits      */
  uint8_t   base_b;             /* typically +190: offset of uchar scores            */
  uint8_t   bias_b;             /* positive bias to emission scores, make them >=0   */


  /* ViterbiFilter uses scaled swords: 8x signed 16-bit integer vectors              */
 #ifdef HAVE_SSE2 
  __m128i **rwv;                /* [x][q]: rw, rw[0] are allocated  [Kp][Q8]         */
  __m128i  *twv;                /* transition score blocks          [8*Q8]           */
 __m128i  *rwv_mem;
  __m128i  *twv_mem;
#endif
#ifdef HAVE_AVX2
 __m256i **rwv_AVX;                /* [x][q]: rw, rw[0] are allocated  [Kp][Q8]         */
  __m256i  *twv_AVX;                /* transition score blocks          [8*Q8]           */
 __m256i  *rwv_mem_AVX;
  __m256i  *twv_mem_AVX;
#endif
#ifdef HAVE_AVX512
 __m512i **rwv_AVX_512;                /* [x][q]: rw, rw[0] are allocated  [Kp][Q8]         */
  __m512i  *twv_AVX_512;                /* transition score blocks          [8*Q8]           */
 __m512i  *rwv_mem_AVX_512;
  __m512i  *twv_mem_AVX_512;
#endif

#ifdef HAVE_NEON
  esl_neon_128i_t **rwv;        /* [x][q]: rw, rw[0] are allocated  [Kp][Q8]         */
  esl_neon_128i_t  *twv;        /* transition score blocks          [8*Q8]           */
  esl_neon_128i_t  *rwv_mem;
  esl_neon_128i_t  *twv_mem;
#endif

  int16_t   xw[p7O_NXSTATES][p7O_NXTRANS]; /* ENJC state transition costs            */
  float     scale_w;            /* score units: typically 500 / log(2), 1/500 bits   */
  int16_t   base_w;             /* offset of sword scores: typically +12000          */
  int16_t   ddbound_w;          /* threshold precalculated for lazy DD evaluation    */
  float     ncj_roundoff;       /* missing precision on NN,CC,JJ after rounding      */

  /* Forward, Backward use IEEE754 single-precision floats: 4x vectors               */
  
  float    xf[p7O_NXSTATES][p7O_NXTRANS]; /* ENJC transition costs                   */

  #ifdef HAVE_SSE2
  __m128 **rfv;                 /* [x][q]:  rf, rf[0] are allocated [Kp][Q4]         */
  __m128  *tfv;                 /* transition probability blocks    [8*Q4]           */
  __m128   *tfv_mem;
  __m128   *rfv_mem;
  #endif
 #ifdef HAVE_AVX2
  __m256 **rfv_AVX;                 /* [x][q]:  rf, rf[0] are allocated [Kp][Q4]         */
  __m256  *tfv_AVX;                 /* transition probability blocks    [8*Q4]           */
  __m256   *tfv_mem_AVX;
  __m256   *rfv_mem_AVX;
#endif
#ifdef HAVE_AVX512
  __m512 **rfv_AVX_512;                 /* [x][q]:  rf, rf[0] are allocated [Kp][Q4]         */
  __m512  *tfv_AVX_512;                 /* transition probability blocks    [8*Q4]           */
  __m512   *tfv_mem_AVX_512;
  __m512   *rfv_mem_AVX_512;
  #endif

#ifdef HAVE_NEON
  esl_neon_128f_t **rfv;        /* [x][q]:  rf, rf[0] are allocated [Kp][Q4]         */
  esl_neon_128f_t  *tfv;        /* transition probability blocks    [8*Q4]           */
  esl_neon_128f_t  *tfv_mem;
  esl_neon_128f_t  *rfv_mem;
#endif

  /* Disk offset information for hmmpfam's fast model retrieval                      */
  off_t  offs[p7_NOFFSETS];     /* p7_{MFP}OFFSET, or -1                             */

  /* Disk offset bookkeeping for h3f:                                                */
  off_t  roff;                  /* record offset (start of record); -1 if none       */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown      */

  /* Information, annotation copied from parent profile:                             */
  char  *name;                  /* unique name of model                              */
  char  *acc;                   /* unique accession of model, or NULL                */
  char  *desc;                  /* brief (1-line) description of model, or NULL      */
  char  *rf;                    /* reference line           1..M; *rf='\0' = unused  */
  char  *mm;                    /* modelmask line           1..M; *mm='\0' = unused  */
  char  *cs;                    /* consensus structure line 1..M, *cs='\0' = unused  */
  char  *consensus;             /* consensus residues for ali display, 1..M          */
  float  evparam[p7_NEVPARAM];  /* parameters for determining E-values, or UNSET     */
  float  cutoff[p7_NCUTOFFS];   /* per-seq/per-dom bit cutoffs, or UNSET             */
  float  compo[p7_MAXABET];     /* per-model HMM filter composition, or UNSET        */
  const ESL_ALPHABET *abc;      /* copy of ptr to alphabet information               */

  /* Information about current configuration, size, allocation                       */
  int    L;                     /* current configured target seq length              */
  int    M;                     /* model length                                      */
  int    max_length;            /* upper bound on emitted sequence length            */
  int    allocM;                /* maximum model length currently allocated for      */
#ifdef HAVE_SSE2  
  int    allocQ4;               /* P7_NVF(allocM): alloc size for tf, rf             */
  int    allocQ8;               /* P7_NVW(allocM): alloc size for tw, rw             */
  int    allocQ16;              /* P7_NVB(allocM): alloc size for rb                 */
#endif
#ifdef HAVE_AVX2  
  int    allocQ4_AVX;               /* P7_NVF_AVX(allocM): alloc size for tf, rf             */
  int    allocQ8_AVX;               /* P7_NVW_AVX(allocM): alloc size for tw, rw             */
  int    allocQ16_AVX;              /* P7_NVB_AVX(allocM): alloc size for rb                 */
#endif
#ifdef HAVE_AVX512  
  int    allocQ4_AVX_512;               /* P7_NVF_AVX_512(allocM): alloc size for tf, rf             */
  int    allocQ8_AVX_512;               /* P7_NVW_AVX_512(allocM): alloc size for tw, rw             */
  int    allocQ16_AVX_512;              /* P7_NVB_AVX_512(allocM): alloc size for rb                 */
#endif
#ifdef HAVE_NEON
  int    allocQ4;               /* P7_NVF(allocM): alloc size for tf, rf             */
  int    allocQ8;               /* P7_NVW(allocM): alloc size for tw, rw             */
  int    allocQ16;              /* P7_NVB(allocM): alloc size for rb                 */
#endif

  int    mode;                  /* currently must be p7_LOCAL                        */
  float  nj;                    /* expected # of J's: 0 or 1, uni vs. multihit       */

  int    is_shadow;             /* TRUE if this profile shadows another, and its ptrs are references */
} P7_OPROFILE;

typedef struct {
  int            count;       /* number of <P7_OPROFILE> objects in the block */
  int            listSize;    /* maximum number elements in the list          */
  P7_OPROFILE  **list;        /* array of <P7_OPROFILE> objects               */
} P7_OM_BLOCK;


// The definition of p7_oprofile_FGetEmission has moved into the p7_alidisplay_<simd> files to make it easier to
// provide SIMD ISA-specific versions of this inlined function.  The original code is retained here for archival purposes
/* retrieve match odds ratio [k][x]
 * this gets used in p7_alidisplay.c, when we're deciding if a residue is conserved or not */
/*static inline float 
p7_oprofile_FGetEmission(const P7_OPROFILE *om, int k, int x)
{
  union { __m128 v; float p[4]; } u;
  int   Q = P7_NVF(om->M);
  int   q = ((k-1) % Q);
  int   r = (k-1)/Q;
  u.v = om->rfv[x][q];
  return u.p[r];
} */

uint8_t unbiased_byteify(P7_OPROFILE *om, float sc);
uint8_t biased_byteify(P7_OPROFILE *om, float sc);
int16_t wordify(P7_OPROFILE *om, float sc);

extern P7_OPROFILE *p7_oprofile_Create(int allocM, const ESL_ALPHABET *abc, SIMD_TYPE simd);
extern void         p7_oprofile_Destroy(P7_OPROFILE *om);
extern int          p7_oprofile_IsLocal(const P7_OPROFILE *om);
extern size_t       p7_oprofile_Sizeof (const P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Clone  (const P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Shadow (const P7_OPROFILE *om);

extern int          p7_oprofile_Convert(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_ReconfigLength    (P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMSVLength (P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigRestLength(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMultihit  (P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigUnihit    (P7_OPROFILE *om, int L);

extern int          p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om);
extern int          p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
				       P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om);
extern int          p7_oprofile_Compare(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
extern int          p7_profile_SameAsMF(const P7_OPROFILE *om, P7_PROFILE *gm);
extern int          p7_profile_SameAsVF(const P7_OPROFILE *om, P7_PROFILE *gm);

extern int          p7_oprofile_GetFwdTransitionArray(const P7_OPROFILE *om, int type, float *arr );
extern int          p7_oprofile_GetMSVEmissionScoreArray(const P7_OPROFILE *om, uint8_t *arr );
extern int          p7_oprofile_GetFwdEmissionScoreArray(const P7_OPROFILE *om, float *arr );
extern int          p7_oprofile_GetFwdEmissionArray(const P7_OPROFILE *om, P7_BG *bg, float *arr );

// SSE versions of SIMD functions
extern P7_OPROFILE *p7_oprofile_Create_sse(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy_sse(P7_OPROFILE *om);
extern size_t       p7_oprofile_Sizeof_sse(const P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Clone_sse(const P7_OPROFILE *om);
extern int          p7_oprofile_Convert_sse(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_Compare_sse(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
extern int          p7_oprofile_GetFwdTransitionArray_sse(const P7_OPROFILE *om, int type, float *arr );
extern int          p7_oprofile_GetMSVEmissionScoreArray_sse(const P7_OPROFILE *om, uint8_t *arr );
extern int          p7_oprofile_GetFwdEmissionScoreArray_sse(const P7_OPROFILE *om, float *arr );
extern int          p7_oprofile_GetFwdEmissionArray_sse(const P7_OPROFILE *om, P7_BG *bg, float *arr );
int sf_conversion_sse(P7_OPROFILE *om);
int mf_conversion_sse(const P7_PROFILE *gm, P7_OPROFILE *om);
int vf_conversion_sse(const P7_PROFILE *gm, P7_OPROFILE *om);
int fb_conversion_sse(const P7_PROFILE *gm, P7_OPROFILE *om);
int oprofile_dump_mf_sse(FILE *fp, const P7_OPROFILE *om);
int oprofile_dump_vf_sse(FILE *fp, const P7_OPROFILE *om);
int oprofile_dump_fb_sse(FILE *fp, const P7_OPROFILE *om, int width, int precision);
int p7_oprofile_Compare_sse(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);

// NEON versions of SIMD functions
extern P7_OPROFILE *p7_oprofile_Create_neon(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy_neon(P7_OPROFILE *om);
extern size_t       p7_oprofile_Sizeof_neon(const P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Clone_neon(const P7_OPROFILE *om);
extern int          p7_oprofile_Convert_neon(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_Compare_neon(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
extern int          p7_oprofile_GetFwdTransitionArray_neon(const P7_OPROFILE *om, int type, float *arr );
extern int          p7_oprofile_GetMSVEmissionScoreArray_neon(const P7_OPROFILE *om, uint8_t *arr );
extern int          p7_oprofile_GetFwdEmissionScoreArray_neon(const P7_OPROFILE *om, float *arr );
extern int          p7_oprofile_GetFwdEmissionArray_neon(const P7_OPROFILE *om, P7_BG *bg, float *arr );
int sf_conversion_neon(P7_OPROFILE *om);
int mf_conversion_neon(const P7_PROFILE *gm, P7_OPROFILE *om);
int vf_conversion_neon(const P7_PROFILE *gm, P7_OPROFILE *om);
int fb_conversion_neon(const P7_PROFILE *gm, P7_OPROFILE *om);
int oprofile_dump_mf_neon(FILE *fp, const P7_OPROFILE *om);
int oprofile_dump_vf_neon(FILE *fp, const P7_OPROFILE *om);
int oprofile_dump_fb_neon(FILE *fp, const P7_OPROFILE *om, int width, int precision);
int p7_oprofile_Compare_neon(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);

// NEON64 versions of SIMD functions
extern P7_OPROFILE *p7_oprofile_Create_neon64(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy_neon64(P7_OPROFILE *om);
extern size_t       p7_oprofile_Sizeof_neon64(const P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Clone_neon64(const P7_OPROFILE *om);
extern int          p7_oprofile_Convert_neon64(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_Compare_neon64(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
extern int          p7_oprofile_GetFwdTransitionArray_neon64(const P7_OPROFILE *om, int type, float *arr );
extern int          p7_oprofile_GetMSVEmissionScoreArray_neon64(const P7_OPROFILE *om, uint8_t *arr );
extern int          p7_oprofile_GetFwdEmissionScoreArray_neon64(const P7_OPROFILE *om, float *arr );
extern int          p7_oprofile_GetFwdEmissionArray_neon64(const P7_OPROFILE *om, P7_BG *bg, float *arr );
int sf_conversion_neon64(P7_OPROFILE *om);
int mf_conversion_neon64(const P7_PROFILE *gm, P7_OPROFILE *om);
int vf_conversion_neon64(const P7_PROFILE *gm, P7_OPROFILE *om);
int fb_conversion_neon64(const P7_PROFILE *gm, P7_OPROFILE *om);
int oprofile_dump_mf_neon64(FILE *fp, const P7_OPROFILE *om);
int oprofile_dump_vf_neon64(FILE *fp, const P7_OPROFILE *om);
int oprofile_dump_fb_neon64(FILE *fp, const P7_OPROFILE *om, int width, int precision);
int p7_oprofile_Compare_neon64(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);



// AVX versions of SIMD functions
extern P7_OPROFILE *p7_oprofile_Create_avx(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy_avx(P7_OPROFILE *om);
extern size_t       p7_oprofile_Sizeof_avx(const P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Clone_avx(const P7_OPROFILE *om);
extern int          p7_oprofile_Convert_avx(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_Compare_avx(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
extern int          p7_oprofile_GetFwdTransitionArray_avx(const P7_OPROFILE *om, int type, float *arr );
extern int          p7_oprofile_GetMSVEmissionScoreArray_avx(const P7_OPROFILE *om, uint8_t *arr );
extern int          p7_oprofile_GetFwdEmissionScoreArray_avx(const P7_OPROFILE *om, float *arr );
extern int          p7_oprofile_GetFwdEmissionArray_avx(const P7_OPROFILE *om, P7_BG *bg, float *arr );
int sf_conversion_avx(P7_OPROFILE *om);
int mf_conversion_avx(const P7_PROFILE *gm, P7_OPROFILE *om);
int vf_conversion_avx(const P7_PROFILE *gm, P7_OPROFILE *om);
int fb_conversion_avx(const P7_PROFILE *gm, P7_OPROFILE *om);
int oprofile_dump_mf_avx(FILE *fp, const P7_OPROFILE *om);
int oprofile_dump_vf_avx(FILE *fp, const P7_OPROFILE *om);
int oprofile_dump_fb_avx(FILE *fp, const P7_OPROFILE *om, int width, int precision);
int p7_oprofile_Compare_avx(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
//AVX-512 versions of SIMD functions
extern P7_OPROFILE *p7_oprofile_Create_avx512(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy_avx512(P7_OPROFILE *om);
extern size_t       p7_oprofile_Sizeof_avx512(const P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Clone_avx512(const P7_OPROFILE *om);
extern int          p7_oprofile_Convert_avx512(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_Compare_avx512(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
extern int          p7_oprofile_GetFwdTransitionArray_avx512(const P7_OPROFILE *om, int type, float *arr );
extern int          p7_oprofile_GetMSVEmissionScoreArray_avx512(const P7_OPROFILE *om, uint8_t *arr );
extern int          p7_oprofile_GetFwdEmissionScoreArray_avx512(const P7_OPROFILE *om, float *arr );
extern int          p7_oprofile_GetFwdEmissionArray_avx512(const P7_OPROFILE *om, P7_BG *bg, float *arr );
int sf_conversion_avx512(P7_OPROFILE *om);
int mf_conversion_avx512(const P7_PROFILE *gm, P7_OPROFILE *om);
int vf_conversion_avx512(const P7_PROFILE *gm, P7_OPROFILE *om);
int fb_conversion_avx512(const P7_PROFILE *gm, P7_OPROFILE *om);
int oprofile_dump_mf_avx512(FILE *fp, const P7_OPROFILE *om);
int oprofile_dump_vf_avx512(FILE *fp, const P7_OPROFILE *om);
int oprofile_dump_fb_avx512(FILE *fp, const P7_OPROFILE *om, int width, int precision);
int p7_oprofile_Compare_avx512(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
#endif /*p7OPROFILE_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
