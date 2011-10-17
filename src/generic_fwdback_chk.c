/* Forward/Backward, checkpointed: generic (non-SIMD) demonstration version.
 * 
 * Contents:
 *    x. 
 *    x. Copyright and license information.
 *    
 * References:
 *    SRE J8/109-112, Oct 2011: implementation plan.
 */

#include "p7_config.h"

#include "easel.h"

#include "hmmer.h"
#include "p7_gmxchk.h"

#define MMXR(p, k) ((p)[(k)* p7G_NSCELLS * p7G_M])
#define IMXR(p, k) ((p)[(k)* p7G_NSCELLS * p7G_I])
#define DMXR(p, k) ((p)[(k)* p7G_NSCELLS * p7G_D])

#define FORWARD_ROW(dpp, dpc) \
  for (k = 1; k < M; k++)\
    {\
      /* match state */\
      sc = p7_FLogsum(p7_FLogsum(MMR(dpp,k-1)   + TSC(p7P_MM,k-1), \
				 IMR(dpp,k-1)   + TSC(p7P_IM,k-1)),\
		      p7_FLogsum(XMX(i-1,p7G_B) + TSC(p7P_BM,k-1),\
				 DMR(dpp,k-1)   + TSC(p7P_DM,k-1)));\
      MMR(dpc,k) = sc + MSC(k);\
      /* insert state */\
      sc = p7_FLogsum(MMR(dpp,k) + TSC(p7P_MI,k),\
		      IMR(dpp,k) + TSC(p7P_II,k));\
      IMR(dpc,k) = sc + ISC(k);\
      /* delete state */\
      DMR(dpc,k) = p7_FLogsum(MMR(dpc,k-1) + TSC(p7P_MD,k-1),\
			      DMR(dpc,k-1) + TSC(p7P_DD,k-1));\
      /* E state update */\
      XMX(i,p7G_E) = p7_FLogsum(p7_FLogsum(MMR(dpc,k) + esc,\
					   DMR(dpc,k) + esc),\
				XMX(i,p7G_E));\
    }\
  /* unrolled match state M_M */\
  sc = p7_FLogsum(p7_FLogsum(MMR(dpp,M-1)   + TSC(p7P_MM,M-1), \
			     IMR(dpp,M-1)   + TSC(p7P_IM,M-1)),\
		  p7_FLogsum(XMX(i-1,p7G_B) + TSC(p7P_BM,M-1),\
			     DMR(dpp,M-1)   + TSC(p7P_DM,M-1)));\
  MMR(dpc,M) = sc + MSC(M);\
  IMR(dpc,M) = -eslINFINITY;\
  /* unrolled delete state D_M */\
  DMR(dpc,M) = p7_FLogsum(MMR(dpc,M-1) + TSC(p7P_MD,M-1),\
			  DMR(dpc,M-1) + TSC(p7P_DD,M-1));


#define GET_SAVE_ROW(gxc,dpc) \
  gxc->R++;\
  (dpc) = gxc->dp[gxc->R0+gxc->R];

#define GET_TEMP_ROW(gxc,dpc) \
  (dpc) = gxc->dp[0];




int
p7_GForwardCheckpointed(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXCHK *gxc, float *opt_sc)
{
  float *dpc;			/* ptr to current row  */
  float *dpp;			/* ptr to previous row */
  float *xmx = gxc->xmx;
  int    M   = gm->M;
  int    W;			/* size of a segment of DP rows ended w/ one checkpointed row, inclusive */
  int    i,k;

  /* Initialization of the zero row, fwd[0] */
  dpc = gxc->dp[gxc->R0-1];	/* i.e., dp[2], the reserved fwd[0] row; 0,1 are reserved for backwards calculation */
  XMX(0,p7G_N) = 0;
  XMX(0,p7G_B) = gm->xsc[p7P_N][p7P_MOVE];                /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;  /* need seq to get here */
  for (k = 0; k <= M; k++)
    MMXR(dpc,k) = IMXR(dpc,k) = DMXR(dpc,k) = -eslINFINITY;
  dpp = dpc;
  
  /* Phase one: "a" region: uncheckpointed rows of matrix */
  for (i = 1; i <= gxc->La; i++)
    {
      GET_SAVE_ROW(gxc, dpc);
      FORWARD_ROW(dpp, dpc);
      dpp = dpc;
    }

  /* Phase two: "b" and "c" regions: partially and fully checkpointed */
  /* i= gxc->La+1, from previous loop's end */
  for (r = gxc->Rb + gxc->Rc, W = gxc->Lb; i++; i <= L)
    {
      W--;
      if (!W) { 		/* we're on the checkpoint: we'll save this row */
	GET_SAVE_ROW(gxc, dpc);
	r--;
	W=r+1;			/* next segment */
      }
      else { GET_TMP_ROW(gxc, dpc); }
    }
  
  

  
  

}

int
p7_GBackwardCheckpointed_PosteriorBands(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXCHK *gxc, P7_GBANDS *bnd);
{



}


      

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
