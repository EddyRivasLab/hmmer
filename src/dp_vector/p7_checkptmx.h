/* P7_CHECKPTMX is the striped SIMD vector, checkpointed DP matrix
 * used by the vectorized local Forwards/Backwards calculation that
 * builds our sparse mask, for subsequent postprocessing with the more
 * complex glocal/local model.
 * 
 * Contents:
 *    1. The P7_CHECKPTMX object and its access macros
 *    2. Function declarations
 *    3. Notes:
 *       [a] Layout of the matrix, in checkpointed rows
 *       [b] Layout of one row, in vectors and floats
 *    4. Copyright and license information.
 */
#ifndef p7CHECKPTMX_INCLUDED
#define p7CHECKPTMX_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "dp_reference/p7_refmx.h"

#define p7C_NSCELLS 3
enum p7c_scells_e {
  p7C_M     = 0,
  p7C_D     = 1,
  p7C_I     = 2,
};

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
  int M;	/* current actual query model dimension (consensus positions)         */
  int L;	/* current actual target seq dimension (residues)                     */
  int R;	/* current actual number of rows (<=Ra+Rb+Rc), excluding R0           */
  int Qf;	/* current actual number of fb vectors = P7_NVF(M)                    */

  /* Checkpointed layout, mapping rows 1..R to residues 1..L:                         */
  int R0;	/* # of extra rows: one for fwd[0] boundary, two for bck[prv,cur]     */
  int Ra;	/* # of rows used in "all" region (uncheckpointed)                    */
  int Rb;	/* # of rows in "between" region (one incomplete checkpoint segment)  */
  int Rc;	/* # of rows in "checkpointed" region                                 */
  int La;	/* residues 1..La are in "all" region                                 */
  int Lb;      	/* residues La+1..La+Lb are in "between" region                       */
  int Lc;	/* residues La+Lb+1..La+Lb+Lc=L are in "checkpointed" region          */

  /* Raw memory allocation */
  char    *dp_mem;	/* raw memory allocation, that dp[] rows point into           */
  int64_t  allocW;	/* alloced width/row, bytes; multiple of p7_VALIGN            */
  int64_t  nalloc;	/* total # of alloc'ed bytes: nalloc >= (validR)(allocW)      */
  int64_t  ramlimit;	/* recommended RAM limit on dp_mem; can temporarily exceed it */

  /* Forward/Backward matrix rows */
  char   **dpf;		/* row ptrs, dpf[0.R0-1,R0..R0+R-1]; aligned on (p7_VALIGN)-byte boundary  */
  int      allocR;	/* allocated size of dpf[]. R+R0 <= R0+Ra+Rb+rc <= validR <= allocR        */
  int      validR;	/* # of dpf[] rows pointing to valid dp_mem; may be < allocR after GrowTo() */

#ifdef p7_DEBUGGING
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
#endif /*p7_DEBUGGING*/
} P7_CHECKPTMX;


extern P7_CHECKPTMX *p7_checkptmx_Create   (int M, int L, int64_t ramlimit);
extern int           p7_checkptmx_GrowTo   (P7_CHECKPTMX *ox, int M, int L);
extern size_t        p7_checkptmx_Sizeof   (const P7_CHECKPTMX *ox);
extern size_t        p7_checkptmx_MinSizeof(int M, int L);
extern int           p7_checkptmx_Reuse    (P7_CHECKPTMX *ox);
extern void          p7_checkptmx_Destroy  (P7_CHECKPTMX *ox);

extern int           p7_checkptmx_SetDumpMode(P7_CHECKPTMX *ox, FILE *dfp, int truefalse);
#ifdef p7_DEBUGGING
extern char *        p7_checkptmx_DecodeX(enum p7c_xcells_e xcode);
extern int           p7_checkptmx_DumpFBHeader(P7_CHECKPTMX *ox);
extern int           p7_checkptmx_DumpFBRow(P7_CHECKPTMX *ox, int rowi, __m128 *dpc, char *pfx);
#endif /*p7_DEBUGGING*/

/*****************************************************************
 * 3. Notes
 *****************************************************************/

/* [a] Layout of the matrix, in checkpointed rows
 * 
 * One P7_CHECKPTMX data structure is used for both Forward and Backward
 * computations on a target sequence. The Forward calculation is
 * checkpointed. The Backward calculation is linear memory in two
 * rows. The end result is a Forward score and a posterior-decoded set
 * of DP bands. (Additionally, the SSV, MSV, and Viterbi filters use
 * a single row of memory from the structure.)
 *
 * In the diagram below, showing the row layout for the main matrix (MDI states):
 *   O = a checkpointed row; 
 *   x = row that isn't checkpointed;
 *   * = boundary row 0, plus row(s) used for Backwards
 * 
 *   i = index of residues in a target sequence of length L
 *   r = index of rows in the DP matrix, R0+R in total
 *
 *               |------------------------- L -------------------------------|   
 *               |-----La----| |-Lb-| |-------------- Lc --------------------|
 * i =  .  .  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
 *      *  *  *  O  O  O  O  O  x  O  x  x  x  x  O  x  x  x  O  x  x  O  x  O
 * r =  0  1  2  3  4  5  6  7  .  8  .  .  .  .  9  .  .  . 10  .  . 11  . 12
 *      |--R0-|  |-----Ra----| |-Rb-| |-------------- Rc --------------------|
 *               |------------------------- R -------------------------------|   
 *   
 * There are four regions in the rows:
 *    region 0 (R0)                : boundary row 0, and Backwards' two rows
 *    region a ("all"; Ra)         : all rows are kept (no checkpointing)
 *    region b ("between"; Rb)     : partially checkpointed
 *    region c ("checkpointed; Rc) : fully checkpointed
 *   
 * In region a, La = Rb
 * In region b, Rb = 0|1, Lb = 0..Rc+1
 *              more specifically: (Rb=0 && Lb=0) || (Rb=1 && 1 <= Lb <= Rc+1)
 * In region c, Lc = {{Rc+2} \choose {2}}-1 = (Rc+2)(Rc+1)/2 - 1
 * 
 * In this example:
 *    R0 = 3
 *    Ra = 5  La = 5
 *    Rb = 1  La = 2
 *    Rc = 4  Lc = 14
 *                                                             
 * In checkpointed regions, we refer to "blocks", often indexed
 * <b>.  There are Rb+Rc blocks, and each block ends in a checkpointed
 * row. The "width" of each block, often called <w>, decrements from
 * Rc+1 down to 2 in the fully checkpointed region.
 *
 * The reason to mix checkpointing and non-checkpointing is that we
 * use as many rows as we can, given a set memory ceiling, to minimize
 * computation time.
 * 
 * The special states (ENJBC) are kept in xmx for all rows 1..L, not
 * checkpointed.
 */


/* [b] Layout of one row, in striped vectors and floats
 *
 *  [1 5 9 13][1 5 9 13][1 5 9 13] [2 6 10 14][2 6 10 14][2 6 10 14] [3 7 11 x][3 7 11 x][3 7 11 x] [4 8 12 x][4 8 12 x][4 8 12 x] [E N JJ J B CC C SCALE]
 *  |-- M ---||-- D ---||-- I ---| |--- M ---||--- D ---||--- I ---| |-- M ---||-- D ---||-- I ---| |-- M ---||-- D ---||-- I ---| 
 *  |---------- q=0 -------------| |------------ q=1 --------------| |---------- q=2 -------------| |---------- q=3 -------------|
 *  |------------------------------------ P7_NVF(M) * p7C_NSCELLS ---------------------------------------------------------------| |---- p7C_NXCELLS ----|
 *  
 *  Number of elements in a vector = p7_VNF      =  4  (assuming 16-byte wide SIMD vectors; 8, for 32-byte AVX vectors)
 *  Number of vectors on a row     = P7_NVF(M)   =  max( 2, ((M-1) / p7_VNF) + 1)
 *  Number of main states          = p7C_NSCELLS =  3  (e.g. M,I,D)
 *  Number of special state vals   = p7C_NXCELLS =  8  (e.g. E, N, JJ, J, B, CC, C, SCALE)
 *  Total size of row              = sizeof(float) * (P7_NVF(M) * P7C_NSCELLS * p7_VNF + p7C_NXCELLS)
 *
 */

#endif /*p7CHECKPTMX_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
