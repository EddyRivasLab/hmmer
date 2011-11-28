/* P7_GMXD : dynamic programming matrix for dual-mode (local/glocal) alignment;
 *           full (quadratic memory, not banded, not checkpointed);
 *           "generic" (standard C code, not vectorized).
 *           
 * This structure is not part of H3's main code path.
 * It is used as an example, and for testing.
 * The main code path uses a more complicated variant, P7_GMXB,
 * which is a banded DP dual-mode matrix; see p7_gmxb.[ch].
 */          

#define p7GD_NSCELLS 6
#define p7GD_ML 0
#define p7GD_MG 1
#define p7GD_IL 2
#define p7GD_IG 3
#define p7GD_DL 4
#define p7GD_DG 5

#define p7GD_NXCELLS 7
#define p7GD_E  0
#define p7GD_N  1
#define p7GD_J  2
#define p7GD_B  3
#define p7GD_L  4
#define p7GD_G  5
#define p7GD_C  6

typedef struct p7_gmxd_s {
  int      M;
  int      L; 

  float   *dp_mem;
  int64_t  allocN;

  float  **dp;
  int      allocR;
  int      allocW;
  int      validR;
} P7_GMXD;


/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
