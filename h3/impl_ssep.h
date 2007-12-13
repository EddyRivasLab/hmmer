/* SSE implementation of Forward filter
 * structures, declarations, and macros
 * 
 * SRE, Thu Dec 13 08:40:08 2007 [Janelia]
 * SVN $Id$
 */
#ifndef P7_IMPL_SSEP_INCLUDED
#define P7_IMPL_SSEP_INCLUDED

#include <esl_alphabet.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#define p7O_QWIDTH  4		/* each _m128 holds four 32-bit IEEE754 floats */
#define p7O_NQ(M)   ( ESL_MAX(2, ((((M)-1) / p7O_QWIDTH) + 1)))


/*****************************************************************
 * 1. P7_OPROFILE: an optimized score profile
 *****************************************************************/

#define p7O_NTRANS    8

enum p7o_xstates_e {  p7O_E = 0,  p7O_N = 1,   p7O_J = 2,   p7O_C = 3 };
#define p7O_NXSTATES  4

enum p7o_xtransitions_e { p7O_MOVE = 0,  p7O_LOOP = 1};
#define p7O_NXTRANS 2

enum p7o_rsc_e {  p7O_MSC = 0,  p7O_ISC = 1 };
#define p7O_NR 2

enum p7o_tsc_e {
  p7O_BM = 0,
  p7O_MM = 1,
  p7O_IM = 2,
  p7O_DM = 3,
  p7O_MD = 4,
  p7O_MI = 5,
  p7O_II = 6,
  p7O_DD = 7
};
#define p7O_NTRANS 8

typedef struct p7_oprofile_s {
  __m128  *tp;	    		/* transition probability blocks               */
  __m128 **rp;     		/* [x][q]:  rp array and rp[0] are allocated   */
  float    xp[p7O_NXSTATES][p7O_NXTRANS];

  const ESL_ALPHABET *abc;
  int mode;
  int M;
  int allocQ;			/* number of quad blocks we're allocated for */

} P7_OPROFILE;

/*****************************************************************
 * 2. P7_OMX: a (two-row) dynamic programming matrix
 *****************************************************************/

enum p7x_scells_e {
  p7X_M = 0, 
  p7X_D = 1,
  p7X_I = 2 
};
#define p7X_NSCELLS 3

#define MMX(q) (dp[(q) * p7X_NSCELLS + p7X_M])
#define DMX(q) (dp[(q) * p7X_NSCELLS + p7X_D])
#define IMX(q) (dp[(q) * p7X_NSCELLS + p7X_I])

typedef struct p7_omx_s {
  __m128 *dp;			/* one row of a striped DP matrix for [0..q-1][MDI] */

  int     M;			/* when omx is in use: how big is the query  */
  int     Q;			/* when omx is in use: how many quads are valid (= p7O_NQ(M))*/
  int     allocQ;		/* total quads allocated */
#ifdef p7_DEBUGGING  
  int     debugging;		/* TRUE if we're in debugging mode                            */
  FILE   *dfp;			/* output stream for diagnostics                              */
#endif
} P7_OMX;


extern P7_OPROFILE *p7_oprofile_Create(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy(P7_OPROFILE *om);
extern int          p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_Dump(FILE *fp, P7_OPROFILE *om);

extern P7_OMX *p7_omx_Create(int allocM);
extern void    p7_omx_Destroy(P7_OMX *ox);
extern int     p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox);

extern int     p7_ForwardFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);


#endif /* P7_IMPL_SSEP_INCLUDED */

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
