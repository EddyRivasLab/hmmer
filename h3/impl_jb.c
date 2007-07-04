/* Optimized implementation of HMMER2 from Jeremy Buhler (Washington
 * University, St. Louis); adapted to HMMER3.
 * 
 * SRE, Wed Jun 27 07:57:26 2007 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"
#include "impl_jb.h"

/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Wed Jun 27 08:07:01 2007 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
P7_OPROFILE *
p7_oprofile_Create(int M, const ESL_ALPHABET *abc)
{
  int          status;
  P7_OPROFILE *om = NULL;
  int          x;

  ESL_ALLOC(om, sizeof(P7_OPROFILE));
  om->tsc = NULL;
  om->rsc = NULL;

  ESL_ALLOC(om->tsc,   sizeof(int)   * M  * p7X_NT);
  ESL_ALLOC(om->rsc,   sizeof(int *) * abc->Kp);                     
  om->rsc[0] = NULL;

  ESL_ALLOC(om->rsc[0],sizeof(int)   * abc->Kp * (M+1) * p7X_NR);
  for (x = 1; x < abc->Kp; x++)
    om->rsc[x] = om->rsc[0] + (x * (M+1) * p7X_NR);

  om->M     = M;
  om->abc_r = abc;

  /* Initializations of unused model transitions 
   * tsc[0][0..6] are unused: set to -infty. [7=TBM is used, though]. These values are touched by DP algs.
   */
  for (x = 0; x < p7X_NT; x++) om->tsc[x] = p7_IMPOSSIBLE;
  return om;

 ERROR:
  p7_oprofile_Destroy(om);
  return NULL;
}

/* Function:  p7_oprofile_Destroy()
 * Synopsis:  Free an optimized profile structure.
 * Incept:    SRE, Wed Jun 27 08:07:08 2007 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
void
p7_oprofile_Destroy(P7_OPROFILE *om)
{
  if (om == NULL) return;
  
  if (om->rsc != NULL) {
    if (om->rsc[0] != NULL) free(om->rsc[0]);
    free(om->rsc);
  }
  if (om->tsc != NULL) free(om->tsc);
  free(om);
  return;
}


/*****************************************************************
 * 2. The P7_OMX structure: a dynamic programming matrix.
 *****************************************************************/

/* Function:  p7_omx_Create()
 * Synopsis:  Create an optimized dynamic programming matrix.
 * Incept:    SRE, Wed Jun 27 08:07:29 2007 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
P7_OMX *
p7_omx_Create(int allocM, int allocL)
{
  int      status;
  P7_OMX  *ox;
  int      i;

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->dp  = NULL;
  ox->xmx = NULL;

  ESL_ALLOC(ox->dp,  sizeof(int *) * (allocL+1));
  ESL_ALLOC(ox->xmx, sizeof(int)   * (allocL+1) * p7X_NXCELLS);
  ox->dp[0] = NULL;
  
  ESL_ALLOC(ox->dp[0], sizeof(int) * (allocL+1) * (allocM+1) * p7X_NSCELLS);
  for (i = 1; i <= allocL; i++)
    ox->dp[i] = ox->dp[0] + i*(allocM+1)*p7X_NSCELLS;
  
  ox->nrows  = allocL+1;
  ox->ncells = (allocM+1)*(allocL+1);
  ox->M      = allocM;
  ox->L      = allocL;
  return ox;

 ERROR:
  p7_omx_Destroy(ox);
  return NULL;
}

/* Function:  p7_omx_Destroy()
 * Synopsis:  Frees an optimized DP matrix.
 * Incept:    SRE, Wed Jun 27 08:07:37 2007 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
void
p7_omx_Destroy(P7_OMX *ox)
{
  if (ox == NULL) return;
  if (ox->dp != NULL) {
    if (ox->dp[0] != NULL) free(ox->dp[0]);
    free(ox->dp);
  }
  if (ox->xmx != NULL) free(ox->xmx);
  free(ox);
}


/*****************************************************************
 * Configuration of/conversion to an optimized profile.
 *****************************************************************/

/* Function:  p7_oprofile_Convert()
 * Synopsis:  Convert generic profile to optimized profile.
 * Incept:    SRE, Wed Jun 27 08:09:02 2007 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
p7_oprofile_Convert(P7_PROFILE *gm, P7_OPROFILE *om)
{
  int  k, x;
  
  /* Transition scores */
  for (k = 1 ; k < gm->M; k++) {
    int *tsc = om->tsc + k*p7X_NT;
    
    tsc[p7X_TMM] = gm->tsc[p7_TMM][k];
    tsc[p7X_TMI] = gm->tsc[p7_TMI][k];
    tsc[p7X_TMD] = gm->tsc[p7_TMD][k];
    tsc[p7X_TIM] = gm->tsc[p7_TIM][k];
    tsc[p7X_TII] = gm->tsc[p7_TII][k];
    tsc[p7X_TDM] = gm->tsc[p7_TDM][k];
    tsc[p7X_TDD] = gm->tsc[p7_TDD][k];
  }
  
  /* Begin scores stored specially off-by-one: [k-1][TBM] is the score for entering at Mk */
  for (k = 0; k < gm->M; k++) {
    int *tsc = om->tsc + k*p7X_NT;
    tsc[p7X_BSC] = gm->bsc[k+1];
  }

  /* Match scores */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1; k <= gm->M; k++)
      {
	int *rsc = om->rsc[x] + k*p7X_NR;
	rsc[p7X_MSC] = gm->msc[x][k];
      }

  /* Insert scores */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1; k < gm->M; k++)
      {
	int *rsc = om->rsc[x] + k*p7X_NR;
	rsc[p7X_ISC] = gm->isc[x][k];
      }

  om->xsc[p7X_XTE][p7X_LOOP] = gm->xsc[p7_XTE][p7_LOOP];  
  om->xsc[p7X_XTE][p7X_MOVE] = gm->xsc[p7_XTE][p7_MOVE];
  om->xsc[p7X_XTN][p7X_LOOP] = gm->xsc[p7_XTN][p7_LOOP];
  om->xsc[p7X_XTN][p7X_MOVE] = gm->xsc[p7_XTN][p7_MOVE];
  om->xsc[p7X_XTC][p7X_LOOP] = gm->xsc[p7_XTC][p7_LOOP];
  om->xsc[p7X_XTC][p7X_MOVE] = gm->xsc[p7_XTC][p7_MOVE];
  om->xsc[p7X_XTJ][p7X_LOOP] = gm->xsc[p7_XTJ][p7_LOOP];
  om->xsc[p7X_XTJ][p7X_MOVE] = gm->xsc[p7_XTJ][p7_MOVE];

  om->abc_r = gm->abc;
  return eslOK;
}

/*****************************************************************
 * 4. Dynamic programming algorithm implementations.
 *****************************************************************/

#define MMX(i,k) (dp[(i)][(k) * p7X_NSCELLS + p7X_XMM])
#define IMX(i,k) (dp[(i)][(k) * p7X_NSCELLS + p7X_XMI])
#define DMX(i,k) (dp[(i)][(k) * p7X_NSCELLS + p7X_XMD])
#define XMX(i,s) (xmx[(i) * p7X_NXCELLS + (s)])

#define TSC(s,k) (tsc[(k) * p7X_NT + (s)])
#define BSC(k)   (tsc[(k) * p7X_NT + p7X_BSC])
#define MSC(k)   (rsc[(k) * p7X_NR + p7X_MSC])
#define ISC(k)   (rsc[(k) * p7X_NR + p7X_ISC])

/* Function:  p7_Viterbi()
 * Synopsis:  Viterbi algorithm.
 * Incept:    SRE, Wed Jun 27 08:49:53 2007 [Janelia]
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
p7_Viterbi(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, int *ret_sc)
{
  int const *tsc = om->tsc;
  int      **dp  = ox->dp;
  int       *xmx = ox->xmx;
  int        M   = om->M;
  int        i,k;

  /* Initialization of the zero row.  */
  XMX(0,p7X_XMN) = 0;                                                /* S->N, p=1            */
  XMX(0,p7X_XMB) = om->xsc[p7X_XTN][p7X_MOVE];                       /* S->N->B, no N-tail   */
  XMX(0,p7X_XME) = XMX(0,p7X_XMC) = XMX(0,p7X_XMJ) = p7_IMPOSSIBLE;  /* need seq to get here */
  for (k = 0; k <= om->M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = p7_IMPOSSIBLE;               /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   */
  for (i = 1; i <= L; i++) 
    {
      int const * rsc = om->rsc[dsq[i]];
      int sc;
	  
      MMX(i,0) = IMX(i,0) = DMX(i,0) = p7_IMPOSSIBLE;
      XMX(i,p7X_XME) = p7_IMPOSSIBLE;
      
      for (k = 1; k < M; k++) 
	{
	  /* match state */
	  sc =             MMX(i-1,k-1)  + TSC(p7X_TMM,k-1);
	  sc = ESL_MAX(sc, IMX(i-1,k-1)  + TSC(p7X_TIM,k-1));
	  sc = ESL_MAX(sc, DMX(i-1,k-1)  + TSC(p7X_TDM,k-1));
	  sc = ESL_MAX(sc, XMX(i-1,p7X_XMB) + BSC(k-1)); /* remember, bsc is slid off-by-one */
	  MMX(i,k) = ESL_MAX(sc + MSC(k), p7_IMPOSSIBLE);

	  /* E state update */
	  XMX(i,p7X_XME) = ESL_MAX(XMX(i,p7X_XME), MMX(i,k));
	  
	  /* insert state */
	  sc =             MMX(i-1,k) + TSC(p7X_TMI,k);
	  sc = ESL_MAX(sc, IMX(i-1,k) + TSC(p7X_TII,k));
	  IMX(i,k) = ESL_MAX(sc + ISC(k), p7_IMPOSSIBLE);
	  
	  /* delete state */
	  sc =             MMX(i,k-1) + TSC(p7X_TMD,k-1);
	  sc = ESL_MAX(sc, DMX(i,k-1) + TSC(p7X_TDD,k-1));
	  DMX(i,k) = ESL_MAX(sc, p7_IMPOSSIBLE);
	}
 
      /* unrolled match state M*/
      sc =             MMX(i-1,M-1)  + TSC(p7X_TMM,M-1);
      sc = ESL_MAX(sc, IMX(i-1,M-1)  + TSC(p7X_TIM,M-1));
      sc = ESL_MAX(sc, DMX(i-1,M-1)  + TSC(p7X_TDM,M-1));
      sc = ESL_MAX(sc, XMX(i-1,p7X_XMB) + BSC(M-1));
      MMX(i,M) = ESL_MAX(sc + MSC(M), p7_IMPOSSIBLE);
      
      /* E state update */
      XMX(i,p7X_XME) = ESL_MAX(XMX(i,p7X_XME), MMX(i,M));
 
      /* Now the special states. Order is important here.
       * remember, N, C and J emissions are zero score by definition.
       */
      
      /* J state */
      sc             =             XMX(i-1,p7X_XMJ) + om->xsc[p7X_XTJ][p7X_LOOP];
      sc             = ESL_MAX(sc, XMX(i,  p7X_XME) + om->xsc[p7X_XTE][p7X_LOOP]);
      XMX(i,p7X_XMJ) = ESL_MAX(sc, p7_IMPOSSIBLE);
      
      /* C state */
      sc             =             XMX(i-1,p7X_XMC) + om->xsc[p7X_XTC][p7X_LOOP];
      sc             = ESL_MAX(sc, XMX(i,  p7X_XME) + om->xsc[p7X_XTE][p7X_MOVE]);
      XMX(i,p7X_XMC) = ESL_MAX(sc, p7_IMPOSSIBLE);
      
      /* N state */
      XMX(i,p7X_XMN) = ESL_MAX(XMX(i-1,p7X_XMN) + om->xsc[p7X_XTN][p7X_LOOP], p7_IMPOSSIBLE);
      
      /* B state */
      sc             =             XMX(i,p7X_XMN) + om->xsc[p7X_XTN][p7X_MOVE];
      sc             = ESL_MAX(sc, XMX(i,p7X_XMJ) + om->xsc[p7X_XTJ][p7X_MOVE]);
      XMX(i,p7X_XMB) = ESL_MAX(sc, p7_IMPOSSIBLE);
    }
  
  /* T state (not stored) */
  *ret_sc = XMX(L,p7X_XMC) + om->xsc[p7X_XTC][p7X_MOVE];
  return eslOK;
}


/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef p7IMPL_JB_BENCHMARK
/* gcc -o benchmark-jb -g -O2 -I. -L. -I../easel -L../easel -Dp7IMPL_JB_BENCHMARK impl_jb.c -lhmmer -leasel -lm
 * ./benchmark-jb <hmmfile>
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-v",        eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual scores",             0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the Buhler implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_RANDOMNESS *r       = esl_randomness_Create(42);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  int             sc;
  int             nullsc;
  double          bitscore;
  

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, p7_UNILOCAL);
  p7_ReconfigLength(gm, L);

  om = p7_oprofile_Create(hmm->M, abc);
  p7_oprofile_Convert(gm, om);

  ox = p7_omx_Create(om->M, L);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, bg->f, abc->K, L, dsq);

      p7_Viterbi(dsq, L, om, ox, &sc);
      p7_bg_NullOne(bg, dsq, L, &nullsc);
      bitscore = p7_SILO2Bitscore(sc - nullsc);
      
      if (esl_opt_GetBoolean(go, "-v")) printf("%.2f bits  (%d SILO)\n", bitscore, sc); 

    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  free(dsq);
  p7_omx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7IMPL_FP_BENCHMARK*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
