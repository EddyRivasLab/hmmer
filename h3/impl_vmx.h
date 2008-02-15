/* :%s/SSE/VMX/g ... in theory */

/* SSE implementation of optimized Viterbi and Forward routines:
 * structures, declarations, and macros.
 * 
 * Currently (this remains in flux as of 14 Dec 07) an optimized
 * implementation is required to provide a ViterbiFilter() and a
 * ForwardFilter() implementation. A call to p7_oprofile_Convert()
 * makes an optimized profile that works for both filters.
 * 
 * A "Filter" returns a score may be an approximation (with
 * characterized or at least characterizable error), and which may
 * have limited upper range, such that high scores are returned as
 * eslINFINITY. Additionally, Filters may only work on local alignment
 * modes, because they are allowed to make assumptions about the range
 * of scores.
 * 
 * Here, ViterbiFilter() is a 8-bit lspace implementation with limited
 * precision and limited range (max 20 bits); ForwardFilter() is a
 * pspace float implementation with correct precision and limited
 * range (max ~127 bits). Both require local mode models.
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
 * SRE, Sun Nov 25 11:23:02 2007
 * SVN $Id$
 */
#ifndef P7_IMPL_VMX_INCLUDED
#define P7_IMPL_VMX_INCLUDED
#include "p7_config.h"

#include <esl_alphabet.h>

#include <altivec.h>

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
 *       are rotated by -1. DD's follow separately, starting at q=0.
 *
 *        {     1     1     1     1      1     1     1 }
 *        {  1593   482x  482x  482x  1593  1593  1593 }    
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
  /* Memory pointers with oversized allocations so we can force boundaries */
  vector unsigned char *tu_mem;
  vector unsigned char *ru_mem;
  vector unsigned char *rm_mem;
  vector float         *tf_mem;
  vector float         *rf_mem;

  /* tu, ru, xu are for ViterbiFilter(): lspace uchars, 16x vectors            */
  vector unsigned char  *tu;	        	/* transition score blocks                     */
  vector unsigned char **ru;     		/* [x][q]:  r16 array and r16[0] are allocated */
  uint8_t                xu[p7O_NXSTATES][p7O_NXTRANS];
  int                    allocQ16;		/* how many uchar vectors                      */

  /* info for the MSPFilter() implementation                                   */
  vector unsigned char **rm;     		/* [x][q]:  m16 array and m16[0] are allocated */
  uint8_t                tbm;		/* constant B->Mk cost: scaled log 2/M(M+1)    */
  uint8_t                tec;		/* constant E->C  cost: scaled log 0.5         */
  uint8_t                tjb;		/* constant J->B  cost: scaled log 3/(L+3)     */

  /* info for the ViterbiFilter() implementation                               */
  int       ddbound_u;	        /* used for lazy DD. it must be signed!        */
  float     scale;		/* typically (2,3)*log(2): half or third bits  */
  uint8_t   base;  	        /* typically +127: offset of uchar scores      */
  uint8_t   bias;	        /* positive bias for emission scores           */

  /* tf, rf, xf are for ForwardFilter():    pspace floats, 4x vectors          */
  vector float  *tf;	    		/* transition probability blocks               */
  vector float **rf;     		/* [x][q]:  rf array and rf[0] are allocated   */
  float          xf[p7O_NXSTATES][p7O_NXTRANS];
  int            allocQ4;		/* how many float vectors                      */

  /* Info specific to the optional ViterbiScore() implementation               */
  float    ddbound_f;		/* full precision bound for lazy DD            */
  int      lspace_f;		/* TRUE if tf,rf are in lspace for Score()     */

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

typedef struct p7_omx_s {
  vector unsigned char *dpu_mem;
  vector unsigned char *dpu;			/* one row of a striped DP matrix for [0..q-1][MDI] for uchars */
  int                   allocQ16;		/* total uchar vectors allocated                               */
  int                   Q16;			/* when omx is in use: how many quads are valid (= p7O_NQU(M)) */

  vector float *dpf_mem;
  vector float *dpf;			/* one row of a striped DP matrix for floats                   */
  int           allocQ4;		/* total float vectors allocated                               */
  int           Q4;			/* when omx is in use: how many quads are valid (= p7O_NQF(M)) */

  int     allocM;		/* current allocation size (redundant with Q16 and Q4, really) */
  int     M;			/* when omx is in use: how big is the query                    */
#ifdef p7_DEBUGGING  
  int     debugging;		/* TRUE if we're in debugging mode                             */
  FILE   *dfp;			/* output stream for diagnostics                               */
#endif
} P7_OMX;


/*****************************************************************
 * 3. Declarations of the external API.
 *****************************************************************/

extern P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy(P7_OPROFILE *om);

extern P7_OMX      *p7_omx_Create(int allocM);
extern int          p7_omx_GrowTo(P7_OMX *ox, int allocM);
extern void         p7_omx_Destroy(P7_OMX *ox);

extern int          p7_oprofile_Dump(FILE *fp, P7_OPROFILE *om);
extern int          p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse);

extern int          p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_ReconfigLength(P7_OPROFILE *om, int L);

extern int p7_MSPFilter    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_ViterbiScore (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);


#endif /* P7_IMPL_VMX_INCLUDED */

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
