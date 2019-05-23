/* H4_REFMX: DP matrix for reference implementations 
 * 
 * Dynamic programming matrix for the reference versions of Forward,
 * Backward, decoding, and alignment calculations. The reference
 * implementation is used for testing and debugging, and to provide a
 * baseline for comparison to our production code.
 * 
 * The reference DP algorithms use quadratic memory (not banded, not
 * checkpointed), with values in standard floats (not vectors).
 * 
 * Contents:
 *    1. The H4_REFMX object and its access macros.
 *    2. Function declarations.
 */          
#ifndef h4REFMX_INCLUDED
#define h4REFMX_INCLUDED

#include "h4_config.h"

/*****************************************************************
 * 1. The H4_REFMX object and its access macros
 *****************************************************************/

typedef struct h4_refmx_s {
  int      M;	     /* current DP matrix values valid for model of length M   */
  int      L;	     /* current DP matrix values valid for seq of length L     */

  float   *dp_mem;   /* matrix memory available. dp[i] rows point into this    */
  int64_t  allocN;   /* # DP cells (floats) allocated. allocN >= allocR*allocW */

  float  **dp;	     /* dp[i] rows of matrix. 0..i..L; L+1 <= validR <= allocR */
  int      allocR;   /* # of allocated rows, dp[]                              */
  int      allocW;   /* width of each dp[i] row, in floats.                    */
  int      validR;   /* # of dp[] ptrs validly placed in dp_mem                */

  int      type;     /* h4R_UNSET | h4R_FORWARD | h4R_BACKWARD | h4R_DECODING  */
} H4_REFMX;

/* Indices for main states in each <dp[i]> supercell.
 * Some code has hardcoded an assumption of this order of the indices.
 */
#define h4R_NSCELLS 6
#define h4R_ML 0
#define h4R_MG 1
#define h4R_IL 2
#define h4R_IG 3
#define h4R_DL 4
#define h4R_DG 5

/* Codes/indices for special states ENJBLGC, and emission-on-transition JJ/CC 
 * Some code (reference_fwdback for example) has hardcoded this order of indices.
 */
#define h4R_NXCELLS 9
#define h4R_E  0
#define h4R_N  1
#define h4R_J  2
#define h4R_B  3
#define h4R_L  4
#define h4R_G  5
#define h4R_C  6
#define h4R_JJ 7	/* JJ (J emission on transition) only needed in decoding matrix */
#define h4R_CC 8	/* CC, ditto */

/* The same data structure gets used in several DP contexts.
 * The <type> field gets set by each algorithm implementation,
 * so p7_refmx_Validate() knows what type of DP matrix it is.
 * 
 * Some of these codes must be sync'ed with h4_sparsemx.h. Some unit
 * tests compare reference, sparse matrices, including their
 * type codes.
 */
#define h4R_UNSET            0   // = h4S_UNSET
#define h4R_FORWARD          1   // = h4S_FORWARD
#define h4R_BACKWARD         2   // = h4S_BACKWARD
#define h4R_DECODING         3   // = h4S_DECODING
#define h4R_VITERBI          4   // = h4S_VITERBI
#define h4R_AEC_ALIGN        5   // = h4S_AEC_ALIGN
#define h4R_ASC_FWD_UP       6
#define h4R_ASC_FWD_DOWN     7
#define h4R_ASC_BCK_UP       8
#define h4R_ASC_BCK_DOWN     9
#define h4R_ASC_DECODE_UP   10
#define h4R_ASC_DECODE_DOWN 11

/* Usually we access the matrix values by stepping pointers thru,
 * exploiting detailed knowledge of their order. Sometimes, either for
 * code clarity or robustness against layout changes, it's worth
 * having access macros, though we can expect these to be relatively
 * expensive to evaluate:
 */
#define H4R_XMX(rmx,i,y)  ( (rmx)->dp[(i)][ ( (rmx)->M +1) * h4R_NSCELLS + (y)] )
#define H4R_MX(rmx,i,k,y) ( (rmx)->dp[(i)][ (k)            * h4R_NSCELLS + (y)] )

/*****************************************************************
 * 2. Function declarations
 *****************************************************************/

extern H4_REFMX *h4_refmx_Create   (int M, int L);
extern int       h4_refmx_GrowTo   (H4_REFMX *rx, int M, int L);
extern int       h4_refmx_SetValues(H4_REFMX *rx, float val);
extern int       h4_refmx_SetType  (H4_REFMX *rx, int M, int L, int type);
extern int       h4_refmx_Scale    (H4_REFMX *rx, float scale);
extern int       h4_refmx_Reuse    (H4_REFMX *rx);
extern void      h4_refmx_Destroy  (H4_REFMX *rx);

extern char     *h4_refmx_DecodeSpecial(int type);
extern char     *h4_refmx_DecodeState  (int type);
extern int       h4_refmx_Dump      (FILE *ofp, H4_REFMX *rx);
extern int       h4_refmx_DumpWindow(FILE *ofp, H4_REFMX *rx, int istart, int iend, int kstart, int kend);
extern int       h4_refmx_CountPath(const H4_PATH *pi, H4_REFMX *rxd);
extern int       h4_refmx_CompareDecoding(const H4_REFMX *ppe, const H4_REFMX *ppa, float a_tol);
extern int       h4_refmx_Validate(H4_REFMX *rmx, char *errbuf);
#endif /*h4REFMX_INCLUDED*/
