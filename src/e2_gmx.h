/* e2_gmx
 *
 *   
*/
#ifndef E2_GMX_INCLUDED
#define E2_GMX_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_getopts.h"	/* ESL_GETOPTS           */
#include "esl_random.h"

#include "e1_rate.h"
#include "e1_bg.h"
#include "e2.h"
#include "e2_profilesq.h"
#include "e2_trace.h"

#include "hmmer.h"

/* linear-memory equivalents of the indices we'll need for the dp */
#define ID(i,j,L)     ( (i) * ( (L) + 1) + (j) )
#define IDF(i,j,L)    ( (i) * ( (L) + 1) + (j) ) % ( (L) + 3 )
#define IDB(i,Li,j,L) ( ((Li)-(i)) * ( (L) + 1) + ((L)-(j)) ) % ( (L) + 3 )


#define e2G_NSCELLS 13
enum e2g_cells_e {
  e2G_BB = 0,   
  e2G_IB = 1,
  e2G_SS = 2,
  e2G_DS = 3,
  e2G_IS = 4,
  e2G_SD = 5,
  e2G_DD = 6,
  e2G_ID = 7,
  e2G_BI = 8,
  e2G_SI = 9,
  e2G_DI = 10,
  e2G_II = 11,
  e2G_ii = 12
 };

#define e2G_NXCELLS 12
#define e2G_EE  0
#define e2G_N1  1
#define e2G_N2  2
#define e2G_J1  3
#define e2G_J2  4
#define e2G_C1  5
#define e2G_C2  6
#define e2G_NN2 7	/* NN2 (J emission on transition) only needed in decoding matrix */
#define e2G_JJ1 8	/* JJ (J emission on transition) only needed in decoding matrix */
#define e2G_JJ2 9	
#define e2G_CC1 10	/* CC, ditto */
#define e2G_CC2 11	

/* the same data structure gets used in several DP contexts.
 * the <type> field gets set by each algorithm implementation,
 * so p7_refmx_Validate() knows what type of DP matrix it is.
 */
#define E2_UNSET     0
#define E2_FORWARD   1
#define E2_BACKWARD  2
#define E2_DECODING  3
#define E2_ALIGNMENT 4
#define E2_VITERBI   5

/*****************************************************************
 * 6. e2_GMX: a "generic" dynamic programming matrix
 *****************************************************************/

typedef struct e2_gmx_s {
  int  M;	        /* current DP matrix values valid for model of length M     */
  int  Lrow;		/* longer  sequence dimension (1..Lrow)                     */
  int  Lcol;		/* shorter sequence dimension (1..Lcol)                     */
  int  L;		/*  (Lrow+1)*(Lcol+1)                                       */
  int  rowsq;           /* identifies the longer sequence [e2P_SL or e2P_SR]        */

 
  float  **dp;          /* logically [0.1..Lcol+2][0..e2G_NCELLS-1]; indexed [i][s] */
  int      allocL;      /* current allocated # of cols : Lcol+3 <= validL <= allocL */
  int      allocW;      /* width of each dp[i] row, in floats.                      */
  int      validL;	/* # of cols actually pointing at DP memory                 */
  uint64_t ncells;	/* total # of allocated cells in array : ncells >= (validL) */

  float   *dp_mem;      /* matrix memory available. dp[i] rows point into this      */
  int64_t  allocN;      /* # DP cells (floats) allocated. allocN >= allocR*allocW   */

  int      type;        /* E2_UNSET | E2_FORWARD | E2_BACKWARD | E2_DECODING        */
} E2_GMX;

/* Macros below implement indexing idioms for generic DP routines.
 * They require the following setup, for profile <gm> and matrix <gx>:
 *   float const *tsc = gm->tsc;
 *   float      **dp  = gx->dp;
 * and for each row i (target residue x_i in digital seq <dsq>):
 *   float const **sssc = gm->sssc;
 *   float const  **ssc = gm->ssc;
 *   float const  **isc = gm->rsc;
 */


#define BBMX(x)                 (dp[(x)][e2G_BB])
#define IBMX(x)                 (dp[(x)][e2G_IB])
#define SSMX(x)                 (dp[(x)][e2G_SS])
#define DSMX(x)                 (dp[(x)][e2G_DS])
#define ISMX(x)                 (dp[(x)][e2G_IS])
#define SDMX(x)                 (dp[(x)][e2G_SD])
#define DDMX(x)                 (dp[(x)][e2G_DD])
#define IDMX(x)                 (dp[(x)][e2G_ID])
#define BIMX(x)                 (dp[(x)][e2G_BI])
#define SIMX(x)                 (dp[(x)][e2G_SI])
#define DIMX(x)                 (dp[(x)][e2G_DI])
#define IIMX(x)                 (dp[(x)][e2G_II]) /* don't need to use e2G_ii here */

#define BBMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_BB])
#define IBMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_IB])
#define SSMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_SS])
#define DSMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_DS])
#define ISMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_IS])
#define SDMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_SD])
#define DDMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_DD])
#define IDMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_ID])
#define BIMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_BI])
#define SIMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_SI])
#define DIMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_DI])
#define iIMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_II]) /* assign II to emit right */
#define IiMXM(x,k)               (dp[(x)][(k) * e2G_NSCELLS + e2G_ii])

#define E2G_XMX(gx,x,s)  ( (gx)->dp[(x)][ ((gx)->M+1) * e2G_NSCELLS + (s)] )
#define E2G_MX(gx,x,k,s) ( (gx)->dp[(x)][ (k)         * e2G_NSCELLS + (s)] )

#define e2TSC(s)        (tsc[(s)])
#define e2hmmerTSC(k,s) (tsc[(k)][(s)])

/* Flags that control E2_GMX debugging dumps */
#define e2_SHOW_LOG      (1<<0)

/* e2_gmx.c */
extern E2_GMX *e2_gmx_Create (int M, int allocL1, int allocL2);
extern int     e2_gmx_GrowTo (E2_GMX *gx, int M, int allocL1, int allocL2);
extern size_t  e2_gmx_Sizeof (E2_GMX *gx);
extern int     e2_gmx_Reuse  (E2_GMX *gx);
extern void    e2_gmx_Destroy(E2_GMX *gx);
extern int     e2_gmx_Compare(E2_GMX *gx1, E2_GMX *gx2, float tolerance);
extern int     e2_gmx_Dump(FILE *fp, E2_GMX *gx, int flags);
extern int     e2_gmx_DumpWindow(FILE *ofp, E2_GMX *gx, int istart, int iend, int jstart, int jend, int kstart, int kend, int flags);
extern int     e2_gmx_NTFromTag(enum e2g_cells_e *ret_e2cell, char *tag);
extern int     e2_gmx_NTtag(enum e2g_cells_e e2cell, char **ret_tag);
extern char   *e2_gmx_DecodeSpecial(int type);
extern char   *e2_gmx_DecodeState(enum e2g_cells_e e2cell);

extern int     e2_gmx_Validate(E2_GMX *gx, char *errbuf);

#endif
