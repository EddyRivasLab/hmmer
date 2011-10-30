/* P7_GMXCHK implementation
 * Checkpointed forward/backward dynamic programming matrix.
 */


#ifndef P7_GMXCHK_INCLUDED

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


extern P7_GMXCHK *p7_gmxchk_Create(int M, int L, int ramlimit);
extern int        p7_gmxchk_GrowTo(P7_GMXCHK *gxc, int M, int L);
extern void       p7_gmxchk_Destroy(P7_GMXCHK *gxc);

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
