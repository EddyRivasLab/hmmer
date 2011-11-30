/* P7_GMXD : dynamic programming matrix for dual-mode (local/glocal) alignment;
 *           full (quadratic memory, not banded, not checkpointed);
 *           "generic" (standard C code, not vectorized).
 *           
 * This structure is not part of H3's main code path.
 * It is used as an example, and for testing.
 * The main code path uses a more complicated variant, P7_GMXB,
 * which is a banded DP dual-mode matrix; see p7_gmxb.[ch].
 */          
#ifndef p7GMXD_INCLUDED
#define p7GMXD_INCLUDED

#include "hmmer.h"

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

/* using these macros requires some variable initialization:
 *    float **dp = gm->dp
 *    int     M  = gm->M
 */
#define P7_GMXD_XMX(i,s) (dp[(i)][ (M+1) * p7GD_NSCELLS + (s)];



/* from p7_gmxd.c */
extern P7_GMXD *p7_gmxd_Create(int M, int L);
extern int      p7_gmxd_GrowTo (P7_GMXD *gxd, int M, int L);
extern int      p7_gmxd_Reuse  (P7_GMXD *gxd);
extern void     p7_gmxd_Destroy(P7_GMXD *gxd);

extern int      p7_gmxd_Dump(FILE *ofp, P7_GMXD *gxd);
extern int      p7_gmxd_DumpWindow(FILE *ofp, P7_GMXD *gxd, int istart, int iend, int kstart, int kend);

/* from generic_fwdback_dual.c */
extern int      p7_GForwardDual(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXD *gxd, float *opt_sc);



#endif /*p7GMXD_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
