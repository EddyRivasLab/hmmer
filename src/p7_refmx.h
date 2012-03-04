/* P7_REFMX :
 *  dynamic programming matrix for dual-mode (local/glocal) alignment;
 *  full (quadratic memory, not banded, not checkpointed);
 *  "reference implementation", baseline for comparison to the
 *  production code.
 *           
 * This structure is not part of H3's main code path.  It is used as
 * an example and for testing. Production code uses two other
 * Forward/Backward implementations: a vectorized banding filter (see
 * impl_{xxx}/fwdfilter.c) and a banded local/glocal alignment (see
 * banded_fwdback.c).
 * 
 */          
#ifndef p7REFMX_INCLUDED
#define p7REFMX_INCLUDED

#include "hmmer.h"

#define p7R_NSCELLS 6
#define p7R_ML 0
#define p7R_MG 1
#define p7R_IL 2
#define p7R_IG 3
#define p7R_DL 4
#define p7R_DG 5

#define p7R_NXCELLS 7
#define p7R_E  0
#define p7R_N  1
#define p7R_J  2
#define p7R_B  3
#define p7R_L  4
#define p7R_G  5
#define p7R_C  6


/* Layout of each row dp[i] of the P7_REFMX dynamic programming matrix:
 * dp[i]:   [ML MG IL IG DL DG] [ML MG IL IG DL DG] [ML MG IL IG DL DG]  ...  [ML MG IL IG DL DG]  [E  N  J  B  L  G  C]
 *     k:   |------- 0 -------| |------- 1 -------| |------- 2 -------|  ...  |------- M -------|  
 *          |--------------------------------- (M+1)*p7R_NSCELLS -------------------------------|  |--- p7R_NXCELLS --|
 * Initializations: * = -inf, . = calculated value, 0 = 0:
 *     0:    *  *  *  *  *  *    *  *  *  *  *  *    *  *  *  *  *  *          *  *  *  *  *  *     *  0  *  .  .  .  *
 *     i:    *  *  *  *  *  *    .  .  .  .  *  .    .  .  .  .  .  .          .  .  *  *  .  .     .  .  .  .  .  .  .
 * Access:
 *  Row dp[r]:                     gxd->dp_mem+(r*allocW) = dpc
 *  Main state s at node k={0..M}: dpc[k*p7R_NSCELLS+s]   
 *  Special state s={ENJBLGC}:     dpc[(M+1)*p7R_NSCELLS+s]
 */
typedef struct p7_gmxd_s {
  int      M;	     /* current DP matrix values valid for model of length M   */
  int      L;	     /* current DP matrix values valid for seq of length L     */

  float   *dp_mem;   /* matrix memory available. dp[i] rows point into this    */
  int64_t  allocN;   /* # DP cells (floats) allocated. allocN >= allocR*allocW */

  float  **dp;	     /* dp[i] rows of matrix. 0..i..L; L+1 <= validR <= allocR */
  int      allocR;   /* # of allocated rows, dp[]                              */
  int      allocW;   /* width of each dp[i] row, in floats.                    */
  int      validR;   /* # of dp[] ptrs validly placed in dp_mem                */
} P7_REFMX;


/* Usually we access the matrix values by stepping pointers thru,
 * exploiting detailed knowledge of their order. Sometimes, either for
 * code clarity or robustness against layout changes, it's worth
 * having access macros, though we can expect these to be relatively
 * expensive to evaluate:
 */
#define P7R_XMX(gxd,i,s)  ( (gxd)->dp[(i)][ ( ((gxd)->M) +1) * p7R_NSCELLS + (s)] )
#define P7R_MX(gxd,i,k,s) ( (gxd)->dp[(i)][ (k)              * p7R_NSCELLS + (s)] )

/* from p7_refmx.c */
extern P7_REFMX *p7_refmx_Create(int M, int L);
extern int       p7_refmx_GrowTo (P7_REFMX *gxd, int M, int L);
extern char *    p7_refmx_DecodeSpecial(int type);
extern int       p7_refmx_Reuse  (P7_REFMX *gxd);
extern void      p7_refmx_Destroy(P7_REFMX *gxd);

extern int       p7_refmx_Dump(FILE *ofp, P7_REFMX *gxd);
extern int       p7_refmx_DumpWindow(FILE *ofp, P7_REFMX *gxd, int istart, int iend, int kstart, int kend);
extern int       p7_refmx_DumpCSV(FILE *fp, P7_REFMX *pp, int istart, int iend, int kstart, int kend);

/* from reference_fwdback.c */
extern int      p7_ReferenceForward (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, float *opt_sc);
extern int      p7_ReferenceBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, float *opt_sc);
extern int      p7_ReferenceDecoding(const P7_PROFILE *gm, const P7_REFMX *fwd, P7_REFMX *bck, P7_REFMX *pp);
extern int      p7_ReferenceAlignMEA(const P7_PROFILE *gm, const P7_REFMX *pp,  P7_REFMX *rmx, P7_TRACE *tr);

/* from reference_viterbi.c */
extern int      p7_ReferenceViterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, P7_TRACE *opt_tr, float *opt_sc);

#endif /*p7REFMX_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
