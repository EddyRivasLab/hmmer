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

/* the same data structure gets used in several DP contexts.
 * the <type> field gets set by each algorithm implementation,
 * so p7_refmx_Validate() knows what type of DP matrix it is.
 */
#define p7R_UNSET    0
#define p7R_FORWARD  1
#define p7R_BACKWARD 2
#define p7R_DECODING 3
#define p7R_MEA      4
#define p7R_VITERBI  5


/* Layout of each row dp[i] of the P7_REFMX dynamic programming matrix:
 * dp[i]:   [ML MG IL IG DL DG] [ML MG IL IG DL DG] [ML MG IL IG DL DG]  ...  [ML MG IL IG DL DG]  [E  N  J  B  L  G  C]
 *     k:   |------- 0 -------| |------- 1 -------| |------- 2 -------|  ...  |------- M -------|  
 *          |--------------------------------- (M+1)*p7R_NSCELLS -------------------------------|  |--- p7R_NXCELLS --|
 * The Validate() routine checks the following pattern: where * = -inf, . = calculated value, 0 = 0:
 * Forward:
 *     0:    *  *  *  *  *  *    *  *  *  *  *  *    *  *  *  *  *  *          *  *  *  *  *  *     *  0  *  .  .  .  *     
 *     1:    *  *  *  *  *  *    .  .  *  *  *  *    .  .  *  *  .  .          .  .  *  *  .  .     .  .  .  .  .  .  . 
 *  2..L:    *  *  *  *  *  *    .  .  .  .  *  *    .  .  .  .  .  .          .  .  *  *  .  .     .  .  .  .  .  .  .
 * Backward:
 *      0:   *  *  *  *  *  *    *  *  *  *  *  *    *  *  *  *  *  *          *  *  *  *  *  *     *  .  *  .  .  .  *
 * 1..L-1:   *  *  *  *  *  *    .  .  .  .  .  .    .  .  .  .  .  .          .  .  *  *  .  .     .  .  .  .  .  .  . 
 *      L:   *  *  *  *  *  *    .  .  *  *  .  .    .  .  *  *  .  .          .  .  *  *  .  .     .  *  *  *  *  *  .  
 * Decoding:
 *      0:   0  0  0  0  0  0    0  0  0  0  0  0    0  0  0  0  0  0          0  0  0  0  0  0     0  0  0  0  0  0  0
 *      1:   0  0  0  0  0  0    .  .  0  0  0  0    .  .  0  0  .  .          .  .  0  0  .  .     0  .  0  0  0  0  0
 * 2..L-1:   0  0  0  0  0  0    .  .  .  .  0  0    .  .  .  .  .  .          .  .  0  0  .  .     0  .  .  0  0  0  .
 *      L:   0  0  0  0  0  0    .  .  0  0  0  0    .  .  0  0  .  .          .  .  0  0  .  .     0  0  0  0  0  0  .
 *
 * rationale:
 *   k=0 columns are only present for indexing k=1..M conveniently
 *   i=0 row is Forward's initialization condition: only S->N->B->{LG} path prefix is possible, and S->N is 1.0
 *   i=0 row is Backward's termination condition: unneeded for posterior decoding; if we need Backwards score, we need N->B->{LG}-> path
 *   DL1,DG1 states removed by entry transition distributions (uniform entry, wing retraction)
 *   DL1 value is valid in Backward because it can be reached (via D->E local exit) but isn't ever used; saves having to special case its nonexistence.
 *   DG1 value is valid in Backward because we intentionally leave D1->{DM} distribution in the P7_PROFILE, for use outside DP algorithms;
 *     in p7_trace_Score() for example. Forward's initialization of DG1 to -inf is sufficient to make DG1 unused in Decoding.
 *   ILm,IGm state never exists.
 *   at i=L, no IL/IG state is possible, because any IL/IG must be followed by at least one more M state and therefore at least one more residue.
 *     IL,IG values at i=L allowed in Forward because they can be reached, but cannot be extended; saves having to special case their nonexistence.
 *   similar for i=1; IL/IG state must be preceded by at least one M state and therefore at least one residue.
 *     IL,IG values at i=1 allowed in Backward because they can be reached, but not extended; saves special casing.
 * 
 * Access:
 *  Row dp[r]:                     gxd->dp_mem+(r*allocW) = dpc
 *  Main state s at node k={0..M}: dpc[k*p7R_NSCELLS+s]   
 *  Special state s={ENJBLGC}:     dpc[(M+1)*p7R_NSCELLS+s]
 */
typedef struct p7_refmx_s {
  int      M;	     /* current DP matrix values valid for model of length M   */
  int      L;	     /* current DP matrix values valid for seq of length L     */

  float   *dp_mem;   /* matrix memory available. dp[i] rows point into this    */
  int64_t  allocN;   /* # DP cells (floats) allocated. allocN >= allocR*allocW */

  float  **dp;	     /* dp[i] rows of matrix. 0..i..L; L+1 <= validR <= allocR */
  int      allocR;   /* # of allocated rows, dp[]                              */
  int      allocW;   /* width of each dp[i] row, in floats.                    */
  int      validR;   /* # of dp[] ptrs validly placed in dp_mem                */

  int      type;     /* p7R_UNSET | p7R_FORWARD | p7R_BACKWARD | p7R_DECODING  */
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
extern int       p7_refmx_Reuse  (P7_REFMX *gxd);
extern void      p7_refmx_Destroy(P7_REFMX *gxd);

extern char *p7_refmx_DecodeSpecial(int type);
extern char *p7_refmx_DecodeState(int type);
extern int   p7_refmx_Dump(FILE *ofp, P7_REFMX *gxd);
extern int   p7_refmx_DumpWindow(FILE *ofp, P7_REFMX *gxd, int istart, int iend, int kstart, int kend);
extern int   p7_refmx_DumpCSV(FILE *fp, P7_REFMX *pp, int istart, int iend, int kstart, int kend);

extern int   p7_refmx_Validate(P7_REFMX *rmx, char *errbuf);

/* from reference_fwdback.c */
extern int p7_ReferenceForward (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, float *opt_sc);
extern int p7_ReferenceBackward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, float *opt_sc);
extern int p7_ReferenceDecoding(const P7_PROFILE *gm, const P7_REFMX *fwd, P7_REFMX *bck, P7_REFMX *pp);
extern int p7_ReferenceAlignMEA(const P7_PROFILE *gm, const P7_REFMX *pp,  P7_REFMX *rmx, P7_TRACE *tr);

/* from reference_viterbi.c */
extern int p7_ReferenceViterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, P7_TRACE *opt_tr, float *opt_sc);

#endif /*p7REFMX_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
