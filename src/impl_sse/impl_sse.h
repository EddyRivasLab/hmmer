/* SSE implementation of optimized Viterbi and Forward routines:
 * structures, declarations, and macros.
 * 
 * Currently (and this remains in flux as of 14 Dec 07) an optimized
 * implementation is required to provide an MSVFilter(),
 * ViterbiFilter() and a ForwardFilter() implementation. A call to
 * p7_oprofile_Convert() makes an optimized profile that works for
 * all filters.
 * 
 * Any "Filter" returns a score may be an approximation (with
 * characterized or at least characterizable error), and which may
 * have limited upper range, such that high scores are returned as
 * eslINFINITY. Additionally, Filters might only work on local
 * alignment modes, because they are allowed to make assumptions about
 * the range of scores.
 * 
 * Here, MSVFilter() and ViterbiFilter() are 8-bit lspace
 * implementations with limited precision and limited range (max 20
 * bits); ForwardFilter() is a pspace float implementation with
 * correct precision and limited range (max ~127 bits). Both require
 * local mode models.
 * 
 * An optimized implementation may also provide other optimized
 * routines. It provides specialized Convert*() functions for these,
 * which may no-op (if the OPROFILE already suffices), or may
 * overwrite parts of the OPROFILE that Filters or other routines
 * might need. Therefore, after using a "bonus" function, a fresh
 * Convert() will be needed before a Filter() is called again. This
 * API is tentative.
 * 
 * For example, here, ViterbiScore() is a 32-bit lspace float SSE
 * implementation of the Viterbi algorithm.
 *
 * A "Score" function might be an additional target for optimization,
 * for example. A "Score" function returns a correct score with full
 * floating-point precision and range, and works for any mode model.
 * 
 * In the generic implementation, profile scores are 32-bit floating
 * point log-odds scores. In an optimized implementation, internally,
 * profile scores can be of any type, and may be in log space (lspace)
 * or probability space (pspace). (Calculations in probability space
 * are useful in the Forward algorithm, but always limit range.)  A
 * shorthand of "lspace uchar" means log-odds scores stored as
 * unsigned chars, for example; "pspace float" means odds ratios
 * stored as floats.
 * 
 * A note on memory alignment: malloc() is required to return a
 * pointer "suitably aligned so that it may be aligned to a pointer of
 * any type of object" (C99 7.20.3). __m128 vectors are 128-bits wide,
 * so malloc() ought to return a pointer aligned on a 16-byte
 * boundary.  However, this is not the case for glibc, and apparently
 * other system libraries. Google turns up threads of arguments
 * between glibc and gcc developers over whose problem this is; this
 * argument has apparently not been resolved, and is of no help.
 * Here, we manually align the relevant pointers by overallocating in
 * *_mem with malloc, then arithmetically manipulating the address to
 * mask off (~0xf).
 * 
 * SRE, Sun Nov 25 11:23:02 2007
 * SVN $Id$
 */
#ifndef P7_IMPL_SSE_INCLUDED
#define P7_IMPL_SSE_INCLUDED

#include "p7_config.h"

#include "esl_alphabet.h"
#include "esl_random.h"

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

/* In calculating Q, the number of vectors we need in a row, we have
 * to make sure there's at least 2, or a striped implementation fails.
 */
#define p7O_NQU(M)   ( ESL_MAX(2, ((((M)-1) / 16) + 1)))   /* 16 uchars */
#define p7O_NQF(M)   ( ESL_MAX(2, ((((M)-1) / 4)  + 1)))   /* 4 floats  */


/*****************************************************************
 * 1. P7_OPROFILE: an optimized score profile
 *****************************************************************/
/* The SSE OPROFILE normally contains scores for the Filters;
 * it can optionally replace the ForwardFilter's 4 pspace floats
 * with 4 lspace floats, for the SSE ViterbiScore function.
 */

#define p7O_NXSTATES  4
#define p7O_NXTRANS   2
#define p7O_NR        2
#define p7O_NTRANS    8
enum p7o_xstates_e      { p7O_E    = 0, p7O_N    = 1,  p7O_J  = 2,  p7O_C  = 3 };
enum p7o_xtransitions_e { p7O_MOVE = 0, p7O_LOOP = 1 };
enum p7o_rsc_e          { p7O_MSC  = 0, p7O_ISC  = 1 };
enum p7o_tsc_e          { p7O_BM   = 0, p7O_MM   = 1,  p7O_IM = 2,  p7O_DM = 3,
			  p7O_MD   = 4, p7O_MI   = 5,  p7O_II = 6,  p7O_DD = 7 };

/* The OPROFILE is striped [Farrar07] and interleaved, as is the DP matrix.
 * For example, the layout of a profile for an M=14 model (xref J2/46):
 * 
 * rsc[x] : interleaved blocks of M and I emissions, starting with q=0
 *                1      1     11     11     1      1      1      1 
 *             1593   1593   2604   2604   371x   371x   482x   482x
 *            [MMMM] [IIII] [MMMM] [IIII] [MMMM] [IIII] [MMMM] [IIII] 
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
typedef struct p7_oprofile_s {
  /* tu, ru, xu are for ViterbiFilter(): lspace uchars, 16x vectors            */
  __m128i  *tu;	        	/* transition score blocks                     */
  __m128i **ru;     		/* [x][q]:  r16 array and r16[0] are allocated */
  uint8_t   xu[p7O_NXSTATES][p7O_NXTRANS];
  int       allocQ16;		/* how many uchar vectors                      */

  /* info for the MSVFilter() implementation                                   */
  __m128i **rm;     		/* [x][q]:  m16 array and m16[0] are allocated */
  uint8_t   tbm;		/* constant B->Mk cost: scaled log 2/M(M+1)    */
  uint8_t   tec;		/* constant E->C  cost: scaled log 0.5         */
  uint8_t   tjb;		/* constant J->B  cost: scaled log 3/(L+3)     */

  /* info for the ViterbiFilter() implementation                               */
  int       ddbound_u;	        /* used for lazy DD. it must be signed!        */
  float     scale;		/* typically (2,3)*log(2): half or third bits  */
  uint8_t   base;  	        /* typically +127: offset of uchar scores      */
  uint8_t   bias;	        /* positive bias for emission scores           */

  /* tf, rf, xf are for ForwardFilter():    pspace floats, 4x vectors          */
  __m128  *tf;	    		/* transition probability blocks               */
  __m128 **rf;     		/* [x][q]:  rf array and rf[0] are allocated   */
  float    xf[p7O_NXSTATES][p7O_NXTRANS];
  int      allocQ4;		/* how many float vectors                      */

  /* Info specific to the optional ViterbiScore() implementation               */
  float    ddbound_f;		/* full precision bound for lazy DD            */
  int      lspace_f;		/* TRUE if tf,rf are in lspace for Score()     */

  /* Our actual vector mallocs, before we align the memory                     */
  __m128i  *tu_mem;
  __m128i  *ru_mem;
  __m128i  *rm_mem;
  __m128   *tf_mem;
  __m128   *rf_mem;

  /* Information copied from parent profile:                                       */
  char  *name;			/* unique name of model                            */
  float  evparam[p7_NEVPARAM]; 	/* parameters for determining E-values             */
  float  cutoff[p7_NCUTOFFS]; 	/* per-seq/per-domain gather, trust, noise cutoffs */
  float  nj;			/* expected # of J's: 0 or 1, uni vs. multihit     */
  int    mode;			/* p7_LOCAL, for example                           */
  int    allocM;		/* maximum model length currently allocated for    */
  int    M;			/* model length                                    */
  const ESL_ALPHABET *abc;	/* copy of ptr to alphabet information             */
} P7_OPROFILE;


/*****************************************************************
 * 2. P7_OMX: a one-row dynamic programming matrix
 *****************************************************************/

enum p7x_scells_e { p7X_M = 0, p7X_D = 1, p7X_I = 2 };
#define p7X_NSCELLS 3

/* Besides ENJBC states, we may also store a rescaling factor on each row  */
enum p7x_xcells_e { p7X_E = 0, p7X_N = 1, p7X_J = 2, p7X_B = 3, p7X_C = 4, p7X_SCALE = 5 }; 
#define p7X_NXCELLS 6

/* 
 * 
 * dpf[][] 
 *    dpf[i] is  / [M1 M2 M3 M4] [D1 D2 D3 D4] [I1 I2 I3 I4] / [M5 M6 M7 M8] ...
 *    to access M(i,k) for i=0,1..L; k=1..M:  dpf[i][(k-1)/4 + p7X_M].element[(k-1)%4]
 * 
 * xmx[] arrays for individual special states:
 *    xmx[ENJBC] = [0 1 2 3][4 5 6 7]..[L-2 L-1 L x]     XRQ >= (L/4)+1
 *    to access B[i] for example, for i=0..L:   xmx[B][i/4].x[i%4]  (quad i/4; element i%4).
 */  
typedef struct p7_omx_s {
  int       M;			/* current actual model dimension                              */
  int       L;			/* current actual sequence dimension                           */

  /* The main dynamic programming matrix for M,D,I states                                      */
  __m128  **dpf;		/* striped DP matrix for [0,1..L][0..Q-1][MDI], float vectors  */
  __m128i **dpu;		/* striped DP matrix for [0,1..L][0..Q-1][MDI], uchar vectors  */
  void     *dp_mem;		/* DP memory shared by <dpu>, <dpf>                            */
  int       allocR;		/* current allocated # rows in dp{uf}. allocR >= validR >= L+1 */
  int       validR;		/* current # of rows actually pointing at DP memory            */
  int       allocQ4;		/* current set row width in <dpf> quads:   allocQ4*4 >= M      */
  int       allocQ16;		/* current set row width in <dpu> 16-mers: allocQ16*16 >= M    */
  size_t    ncells;		/* current allocation size of <dp_mem>, in accessible cells    */

  /* The X states (for full,parser; or NULL, for scorer                                        */
  float    *xmx[p7X_NXCELLS];	/* [ENJBC][0..(allocXRQ*4)-1]                                  */
  void     *x_mem;		/* X memory before 16-byte alignment                           */
  int       allocXRQ;		/* # of quads allocated in each xmx[] array; XRQ*4 >= L+1      */
  float     totscale;		/* log of the product of all scale factors (0.0 if unscaled)   */

#ifdef p7_DEBUGGING  
  /* Parsers,scorers only hold a row at a time, so to get them to dump full matrix, it
   * must be done during a DP calculation, after each row is calculated 
   */
  int     debugging;		/* TRUE if we're in debugging mode                             */
  FILE   *dfp;			/* output stream for diagnostics                               */
#endif
} P7_OMX;

/* ?MXo(q) access macros work for either uchar or float, so long as you
 * init your "dp" to point to the appropriate array.
 */
#define MMXo(q) (dp[(q) * p7X_NSCELLS + p7X_M])
#define DMXo(q) (dp[(q) * p7X_NSCELLS + p7X_D])
#define IMXo(q) (dp[(q) * p7X_NSCELLS + p7X_I])


/*****************************************************************
 * 3. Declarations of the external API.
 *****************************************************************/

/* p7_omx.c */
extern P7_OMX      *p7_omx_Create(int allocM, int allocL, int allocXL);
extern int          p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL, int allocXL);
extern void         p7_omx_Destroy(P7_OMX *ox);
extern int          p7_omx_DomainPosteriors(P7_OPROFILE *om, P7_OMX *oxf, P7_OMX *oxb, P7_DOMAINDEF *ddef);

extern int          p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse);
extern int          p7_omx_DumpCharRow(P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
extern int          p7_omx_DumpFloatRow(P7_OMX *ox, int logify, int rowi, int width, int precision, float xE, float xN, float xJ, float xB, float xC);
extern int          p7_omx_DumpMSVRow(P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);

/* p7_oprofile.c */
extern P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy(P7_OPROFILE *om);

extern int          p7_oprofile_Convert(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_ReconfigLength(P7_OPROFILE *om, int L);
extern int          p7_oprofile_Logify(P7_OPROFILE *om);
extern int          p7_oprofile_Probify(P7_OPROFILE *om);
extern int          p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
				       P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om);
extern int          p7_oprofile_SameRounding(const P7_OPROFILE *om, P7_PROFILE *gm);
extern int          p7_oprofile_SameMSV(const P7_OPROFILE *om, P7_PROFILE *gm);
extern int          p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om);

/* fbparsers.c */
extern int p7_ForwardParser (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *ret_sc);
extern int p7_BackwardParser(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *ret_sc);

/* fwdfilter.c */
extern int p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

/* msvfilter.c */
extern int p7_MSVFilter    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

/* vitfilter.c */
extern int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);

/* vitscore.c */
extern int p7_ViterbiScore (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);


#endif /* P7_IMPL_SSE_INCLUDED */

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
