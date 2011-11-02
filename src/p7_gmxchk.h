/* P7_GMXCHK implementation
 * Checkpointed forward/backward dynamic programming matrix.
 */

/* 
 * One P7_GMXCHK data structure is used for both Forward and
 * Backward computations on a target sequence. The result is
 * a Forward score and a posterior-decoded set of DP bands.
 * 
 * The Forward matrix may be checkpointed. The Backwards matrix is
 * linear-memory with two rows.
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
 *              more specificially: (Rb=0 && Lb=0) || (Rb=1 && 1 <= Lb <= Rc+1)
 * In region c, Lc = {{Rc+2} \choose {2}}-1 = (Rc+2)(Rc+1)/2 - 1
 * 
 * In this example:
 *    R0 = 3
 *    Ra = 5  La = 5
 *    Rb = 1  La = 2
 *    Rc = 4  Lc = 14
 *                                                             
 * In checkpointed regions, we will refer to "blocks", often indexed
 * <b>.  There are Rb+Rc blocks, and each block ends in a checkpointed
 * row. The "width" of each block, often called <w>, decrements from
 * Rc+1 down to 2 in the fully checkpointed region.
 *
 * The reason to mix checkpointing and non-checkpointing is that we
 * use as many rows as we can, given a set memory ceiling, to minimize
 * computation time.
 * 
 * The special states (ENJBC) are kept in xmx for all rows 1..L, just
 * as in a normal (uncheckpointed) P7_GMX.
 */
#ifndef P7_GMXCHK_INCLUDED
#define P7_GMXCHK_INCLUDED

#include "p7_config.h"

typedef struct p7_gmxchk_s {
  int      M;	        /* actual query model dimension of current comparison                 */
  int      L;	        /* actual target sequence dimension of current comparison             */
  int      R;	        /* actual # rows in current fwd matrix (<= Ra+Rb+Rc), excluding R0    */
  
  /* Checkpointed layout, mapping rows 1..R to residues 1..L:                                 */
  int      R0;	        /* # of extra rows: one for fwd[0] boundary, two for bck[prv,cur]     */
  int      Ra;	        /* # of rows used in "all" region (uncheckpointed)                    */
  int      Rb;	        /* # of rows in "between" region (one incomplete checkpoint segment)  */
  int      Rc;	        /* # of rows in "checkpointed" region                                 */
  int      La;	        /* residues 1..La are in "all" region                                 */
  int      Lb;      	/* residues La+1..La+Lb are in "between" region                       */
  int      Lc;	        /* residues La+Lb+1..La+Lb+Lc=L are in "checkpointed" region          */

  int      allocW;	/* allocated width per row, in supercells (M+1 <= allocW)             */

  float   *dp_mem;	/* raw memory allocation, that dp[] rows point into                   */
  int64_t  ncells;	/* total # of alloc'ed supercells: ncells >= (validR)(allocW)         */
  int64_t  ncell_limit;	/* recommended RAM limit on dp_mem; can temporarily exceed it         */

  float  **dp;		/* dp[0..R0-1,R0..R0+R-1][0.1..M][0..p7G_NSCELLS-1]; indexed [r][k*p7G_NSCELLS+s] */
  int      allocR;	/* allocated size of dp[]. R+R0 <= R0+Ra+Rb+Rc <= validR <= allocR                */
  int      validR;	/* # of rows pointing at DP memory; may be < allocR after a GrowTo() call         */ 

  float   *xmx;		/* logically [0.1..L][p7G_NXCELLS-1]; indexed [i*p7G_NXCELLS+s]       */
  int      allocXR;	/* allocated # of rows for special states. L+1 <= allocXR             */
} P7_GMXCHK;

#define MMR(p, k) ((p)[(k)* p7G_NSCELLS + p7G_M])
#define IMR(p, k) ((p)[(k)* p7G_NSCELLS + p7G_I])
#define DMR(p, k) ((p)[(k)* p7G_NSCELLS + p7G_D])

extern P7_GMXCHK *p7_gmxchk_Create (int M, int L, int ramlimit);
extern int        p7_gmxchk_GrowTo (P7_GMXCHK *gxc, int M, int L);
extern size_t     p7_gmxchk_Sizeof (const P7_GMXCHK *gxc);
extern int        p7_gmxchk_Reuse  (P7_GMXCHK *gxc);
extern void       p7_gmxchk_Destroy(P7_GMXCHK *gxc);

extern int        p7_gmxchk_Dump(FILE *ofp, P7_GMXCHK *gxc, int flags);

#endif /*P7_GMXCHK_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * References:
 *    SRE:J8/109-112, Oct 2011: Implementation plan
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
