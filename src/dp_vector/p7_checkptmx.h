/* P7_CHECKPTMX: checkpointed, striped vector DP matrix.
 * 
 * Independent of vector ISA. (Do not add any ISA-specific code.)
 * See p7_checkptmx.md for notes.
 */
#ifndef p7CHECKPTMX_INCLUDED
#define p7CHECKPTMX_INCLUDED
#include "p7_config.h"

#include <stdio.h>

#include "dp_reference/p7_refmx.h"
#include "dp_vector/simdvec.h"


#define p7C_NSCELLS 3
enum p7c_scells_e { p7C_M = 0,  p7C_D = 1,  p7C_I = 2 };

#define p7C_NXCELLS 8
enum p7c_xcells_e {
  p7C_E     = 0,
  p7C_N     = 1,
  p7C_JJ    = 2,
  p7C_J     = 3,
  p7C_B     = 4,
  p7C_CC    = 5,
  p7C_C     = 6,
  p7C_SCALE = 7
};

#define P7C_MQ(dp, q)     ((dp)[(q) * p7C_NSCELLS + p7C_M])
#define P7C_DQ(dp, q)     ((dp)[(q) * p7C_NSCELLS + p7C_D])
#define P7C_IQ(dp, q)     ((dp)[(q) * p7C_NSCELLS + p7C_I])

typedef struct p7_checkptmx_s {
  /* Aligned memory allocation (for all rows):                                       */
  int M;	/* current actual query model dimension (consensus positions)        */
  int L;	/* current actual target seq dimension (residues)                    */
  int V;        /* width of striped vectors, in # elements. 0: no striping set yet   */
  int Q;        /* current actual number of fb vectors = P7_NVF(M)                   */
  int R;        /* current actual number of rows (<=Ra+Rb+Rc), excluding R0          */
  float   *dp_mem;     /* raw memory allocation, that dp[] rows point into           */
  int64_t  allocW;     /* alloced width/row, bytes; multiple of p7_MAXVF             */
  int64_t  nalloc;     /* total # of alloc'ed bytes: nalloc >= (validR)(allocW)      */
  int64_t  redline;    /* recommended RAM limit on dp_mem; can temporarily exceed it */

  /* Forward/Backward matrix rows:                                                     */
  float  **dpf;		/* row ptrs, dpf[0.R0-1,R0..R0+R-1], aligned memory            */
  int      allocR;	/* allocated dpf[]. R+R0 <= R0+Ra+Rb+rc <= validR <= allocR    */
  int      validR;	/* # dpf[] rows on valid dp_mem; maybe < allocR after Reinit() */

  /* Checkpointed layout, mapping rows 1..R to residues 1..L:                         */
  int R0;	/* # of extra rows: one for fwd[0] boundary, two for bck[prv,cur]     */
  int Ra;	/* # of rows used in "all" region (uncheckpointed)                    */
  int Rb;	/* # of rows in "between" region (one incomplete checkpoint segment)  */
  int Rc;	/* # of rows in "checkpointed" region                                 */
  int La;	/* residues 1..La are in "all" region                                 */
  int Lb;      	/* residues La+1..La+Lb are in "between" region                       */
  int Lc;	/* residues La+Lb+1..La+Lb+Lc=L are in "checkpointed" region          */

#if eslDEBUGLEVEL > 0
  /* Info for dumping debugging info, conditionally compiled                        */
  int       do_dumping;		/* TRUE if matrix is in dumping mode                */
  FILE     *dfp;		/* open output stream for debug dumps               */
  int       dump_maxpfx;	/* each line prefixed by tag of up to this # chars  */
  int       dump_width;		/* cell values in diagnostic output are fprintf'ed: */
  int       dump_precision;	/*   dfp, "%*.*f", dbg_width, dbg_precision, val    */
  uint32_t  dump_flags;		/* p7_DEFAULT | p7_HIDE_SPECIALS | p7_SHOW_LOG      */

  P7_REFMX *fwd;		/* full Forward matrix, saved for unit test diffs   */
  P7_REFMX *bck;		/* ... full Backward matrix, ditto                  */
  P7_REFMX *pp;			/* ... full posterior probability matrix, ditto     */
  float     bcksc;		/* Backwards score: which we check against Forward  */
#endif // eslDEBUGLEVEL
} P7_CHECKPTMX;


extern P7_CHECKPTMX *p7_checkptmx_Create(int M, int L, int64_t redline);
extern int           p7_checkptmx_Reinit(P7_CHECKPTMX *ox, int M, int L);
extern size_t        p7_checkptmx_Sizeof(const P7_CHECKPTMX *ox);
extern size_t        p7_checkptmx_MinSizeof(int M, int L);
extern void          p7_checkptmx_Destroy(P7_CHECKPTMX *ox);

#if eslDEBUGLEVEL > 0
extern int   p7_checkptmx_SetDumpMode (P7_CHECKPTMX *ox, FILE *dfp, int true_or_false);
extern int   p7_checkptmx_DumpFBHeader(P7_CHECKPTMX *ox);
extern int   p7_checkptmx_DumpFBRow   (P7_CHECKPTMX *ox, int rowi, float *dpc, char *pfx);
#endif  // esl_DEBUGLEVEL

#endif /*p7CHECKPTMX_INCLUDED*/

