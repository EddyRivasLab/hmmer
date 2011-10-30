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

#define MMR(p, k) ((p)[(k)* p7G_NSCELLS * p7G_M])
#define IMR(p, k) ((p)[(k)* p7G_NSCELLS * p7G_I])
#define DMR(p, k) ((p)[(k)* p7G_NSCELLS * p7G_D])

/*  dpp - pointer into previous row
 *  dpc - pointer into current row
 *  i   - index of current row
 *  tsc - points  at gm->tsc
 *  xmx - points at gxc->xmx
 *  rsc - points at gm->rsc[dsq[i]]
 *  esc - 0 or -eslINFINITY for local vs. glocal
 *  gm  - the model
 *  M   - model length
 *  
 *  int k
 *  float sc
 */
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
			  DMR(dpc,M-1) + TSC(p7P_DD,M-1));\
  /* unrolled E state update */\
  XMX(i,p7G_E) = p7_FLogsum(p7_FLogsum(MMR(dpc,M),\
	 		               DMR(dpc,M)),\
				       XMX(i,p7G_E));\
  /* J state */\
  XMX(i,p7G_J) = p7_FLogsum(XMX(i-1,p7G_J) + gm->xsc[p7P_J][p7P_LOOP],\
	                    XMX(i,  p7G_E) + gm->xsc[p7P_E][p7P_LOOP]);\
  /* C state */\
  XMX(i,p7G_C) = p7_FLogsum(XMX(i-1,p7G_C) + gm->xsc[p7P_C][p7P_LOOP],\
 			    XMX(i,  p7G_E) + gm->xsc[p7P_E][p7P_MOVE]);\
  /* N state */\
 XMX(i,p7G_N) = XMX(i-1,p7G_N) + gm->xsc[p7P_N][p7P_LOOP];\
 /* B state */\
 XMX(i,p7G_B) = p7_FLogsum(XMX(i,  p7G_N) + gm->xsc[p7P_N][p7P_MOVE],\
			   XMX(i,  p7G_J) + gm->xsc[p7P_J][p7P_MOVE]);


int
p7_GForwardCheckpointed(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXCHK *gxc, float *opt_sc)
{
  float const *tsc = gm->tsc;
  float       *xmx = gxc->xmx;
  float *dpc;			/* ptr to current row  */
  float *dpp;			/* ptr to previous row */
  float  sc;
  int    M   = gm->M;
  int    W;			/* size of a segment of DP rows ended w/ one checkpointed row, inclusive */
  int    r;			/* counts down over segment number (Rb+Rc..1); each segment has r+1 rows */
  int    i,k;
  float  esc = p7_profile_IsLocal(gm) ? 0. : -eslINFINITY;

  /* Initialization of the zero row, fwd[0] */
  dpc = gxc->dp[gxc->R0-1];	/* i.e., dp[2], the reserved fwd[0] row; 0,1 are reserved for tmp space and backwards calculation */
  XMX(0,p7G_N) = 0;
  XMX(0,p7G_B) = gm->xsc[p7P_N][p7P_MOVE];                    /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;  /* need seq to get here */
  for (k = 0; k <= M; k++)
    MMR(dpc,k) = IMR(dpc,k) = DMR(dpc,k) = -eslINFINITY;
  dpp = dpc;
  
  /* Phase one: "a" region: uncheckpointed rows of matrix */
  for (i = 1; i <= gxc->La; i++)
    {
      float const *rsc = gm->rsc[dsq[i]];
      dpc = gxc->dp[gxc->R0+gxc->R]; gxc->R++; /* idiomatic for "get next saved row" */
      FORWARD_ROW(dpp, dpc);
      dpp = dpc;
    }

  /* Phase two: "b" and "c" regions: partially and fully checkpointed */
  /* i= gxc->La+1, from previous loop's end */
  for (r = gxc->Rb + gxc->Rc, W = gxc->Lb; i <= L; i++)
    {
      float const *rsc = gm->rsc[dsq[i]];
      W--;
      if (!W)
	{ 		                           /* we're on the last row in segment: this row is saved    */
	  dpc = gxc->dp[gxc->R0+gxc->R]; gxc->R++; /* idiomatic for "get next saved row"                     */
	  r--;					   /* decrement segment number counter; last segment is r=1  */
	  W=r+1;			           /* next segment has this many rows, ending in a saved row */
	}
      else
	dpc =  gxc->dp[i%2];	/* idiomatic for "get next temp row", 0 or 1; i%2 is to make sure dpp != dpc  */
      
      FORWARD_ROW(dpp, dpc);
      dpp = dpc;
    }
  
  gxc->M = M;
  gxc->L = L;
  gxc->R = gxc->Ra + gxc->Rb + gxc->Rc;
  if (opt_sc) *opt_sc = XMX(L,p7G_C) + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}

/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_CHK_EXAMPLE
/* 
   gcc -g -O2 -o generic_fwdback_chk_example -Dp7GENERIC_FWDBACK_CHK_EXAMPLE -I. -I../easel -L. -L../easel generic_fwdback_chk.c -lhmmer -leasel -lm
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_gmxchk.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of checkpointed Forward/Backward, generic implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *fwd     = NULL;
  P7_GMXCHK      *fwdc    = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fsc;
  float           nullsc;
  int             status;

  /* Initialize log-sum calculator */
  p7_FLogsumInit();

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);

  /* Allocate matrices */
  fwd  = p7_gmx_Create(gm->M, sq->n);
  fwdc = p7_gmxchk_Create(gm->M, sq->n, 32);

  while ( (status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
      else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

      /* Resize the DP matrices if necessary */
      p7_gmx_GrowTo(fwd, gm->M, sq->n);
      p7_gmx_GrowTo(bck, gm->M, sq->n);

      /* Set the profile and null model's target length models */
      p7_bg_SetLength(bg,   sq->n);
      p7_ReconfigLength(gm, sq->n);

      /* Run Forward, Backward */
      p7_GForward (sq->dsq, sq->n, gm, fwd, &fsc);
      p7_GBackward(sq->dsq, sq->n, gm, bck, &bsc);

      //p7_gmx_Dump(stdout, fwd, p7_DEFAULT);

      /* Those scores are partial log-odds likelihoods in nats.
       * Subtract off the rest of the null model, convert to bits.
       */
      p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

      printf("%-30s   %10.4f %10.4f   %10.4f %10.4f\n", 
	     sq->name, 
	     fsc, bsc, 
	     (fsc - nullsc) / eslCONST_LOG2, (bsc - nullsc) / eslCONST_LOG2);

      p7_gmx_Reuse(fwd);
      p7_gmx_Reuse(bck);
      esl_sq_Reuse(sq);
    }

  /* Cleanup */
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_FWDBACK_EXAMPLE*/


/*------------------- end, example driver -----------------------*/
      

/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
