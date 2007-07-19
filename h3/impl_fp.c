/* The test "floating-point" implementation.
 * 
 * Based on the Buhler implementation.
 * 
 * SRE, Fri Jun 22 14:29:54 2007 [Janelia]
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
#include "impl_fp.h"


/*****************************************************************
 * 1. The P7_OPROFILE structure: a score profile.
 *****************************************************************/

/* Function:  p7_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Fri Jun 22 14:31:28 2007 [Janelia]
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

  ESL_ALLOC(om->tsc,   sizeof(float)   * M  * p7X_NT);
  ESL_ALLOC(om->rsc,   sizeof(float *) * abc->Kp);                     
  om->rsc[0] = NULL;

  ESL_ALLOC(om->rsc[0],sizeof(float)   * abc->Kp * (M+1) * p7X_NR);
  for (x = 1; x < abc->Kp; x++)
    om->rsc[x] = om->rsc[0] + x*(M+1)*p7X_NR;

  om->M     = M;
  om->abc_r = abc;

  /* Initializations of unused model transitions 
   * tsc[0][0..6] are unused: set to -infty. [7=TBM is used, though]. These values are touched by DP algs.
   */
  for (x = 0; x < p7X_NT; x++) om->tsc[x] = -eslINFINITY;
  return om;


 ERROR:
  p7_oprofile_Destroy(om);
  return NULL;
}

/* Function:  p7_oprofile_Destroy()
 * Synopsis:  Free an optimized profile structure.
 * Incept:    SRE, Fri Jun 22 14:43:30 2007 [Janelia]
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
 * Incept:    SRE, Fri Jun 22 14:50:27 2007 [Janelia]
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

  ESL_ALLOC(ox->dp,  sizeof(float *) * (allocL+1));
  ESL_ALLOC(ox->xmx, sizeof(float)   * (allocL+1) * p7X_NXCELLS);
  ox->dp[0] = NULL;
  
  ESL_ALLOC(ox->dp[0], sizeof(float) * (allocL+1) * (allocM+1) * p7X_NSCELLS);
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
 * Incept:    SRE, Fri Jun 22 14:59:07 2007 [Janelia]
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
 * 3. Configuration of an optimized profile.
 *****************************************************************/

/* Function:  p7_oprofile_Config()
 * Synopsis:  Configure a score profile, from an HMM.
 * Incept:    SRE, Fri Jun 22 15:04:09 2007 [Janelia]
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
p7_oprofile_Config(const P7_HMM *hmm, const P7_BG *bg, P7_OPROFILE *om, int mode)
{
  int    status;
  float *occ = NULL;
  int    k, x;
  float  Z;
  float  *tp, *rp;
  float   sc[p7_MAXCODE];

  /* E state loop/move probabilities: nonzero for MOVE allows loops/multihits
   * N,C,J transitions are set later by length config 
   */
  if (mode == p7_LOCAL || mode == p7_GLOCAL) {
    om->xsc[p7X_XTE][p7X_MOVE] = log(0.5);  
    om->xsc[p7X_XTE][p7X_LOOP] = log(0.5);  
  } else {
    om->xsc[p7X_XTE][p7X_MOVE]  = 0.;  
    om->xsc[p7X_XTE][p7X_LOOP]  = -eslINFINITY;  
  }

  /* Precalculations for begin probabilities, which
   * are occ[k] /( \sum_i occ[i] * (M-k+1))
   */
  ESL_ALLOC(occ, sizeof(float) * (hmm->M+1));
  if ((status = p7_hmm_CalculateOccupancy(hmm, occ)) != eslOK) goto ERROR;
  for (Z = 0., k = 1; k <= hmm->M; k++) 
    Z += occ[k] * (float) (hmm->M-k+1);

  /* Transition probabilities. */
  for (k = 1; k < om->M; k++) {
    tp = om->tsc + k * p7X_NT;
    
    tp[p7X_TMM] = log(hmm->t[k][p7H_MM]);
    tp[p7X_TMI] = log(hmm->t[k][p7H_MI]);
    tp[p7X_TMD] = log(hmm->t[k][p7H_MD]);
    tp[p7X_TIM] = log(hmm->t[k][p7H_IM]);
    tp[p7X_TII] = log(hmm->t[k][p7H_II]);
    tp[p7X_TDM] = log(hmm->t[k][p7H_DM]);
    tp[p7X_TDD] = log(hmm->t[k][p7H_DD]);
    tp[p7X_BSC] = log(occ[k+1] / Z);        /* TBM is stored specially off-by-one: [k-1][TBM] is the score for entering at Mk  */
  }
  tp = om->tsc;			/* [0][BSC] */
  tp[p7X_BSC] = log(occ[1] / Z);

  /* Match emission probabilities. */
  sc[om->abc_r->K]     = -eslINFINITY; /* gap character */
  sc[om->abc_r->Kp-1]  = -eslINFINITY; /* missing data character */
  for (k = 1; k <= om->M; k++) {
    for (x = 0; x < om->abc_r->K; x++) 
      sc[x] = log(hmm->mat[k][x] / bg->f[x]);
    esl_abc_FExpectScVec(om->abc_r, sc, bg->f); 
    for (x = 0; x < om->abc_r->Kp; x++) {
      rp = om->rsc[x] + k*p7X_NR;
      rp[p7X_MSC] = sc[x];
    }
  }

  /* Insert emission probabilities. */
  sc[om->abc_r->K]     = -eslINFINITY; /* gap character */
  sc[om->abc_r->Kp-1]  = -eslINFINITY; /* missing data character */
  for (k = 1; k < om->M; k++) {
    for (x = 0; x < om->abc_r->K; x++) 
      sc[x] = log(hmm->ins[k][x] / bg->f[x]); 
    esl_abc_FExpectScVec(om->abc_r, sc, bg->f); 
    for (x = 0; x < om->abc_r->Kp; x++) {
      rp = om->rsc[x] + k*p7X_NR;
      rp[p7X_ISC] = sc[x];
    }
  }    

  /* Remaining specials, [N-CJ][MOVE | LOOP] are set by ReconfigLength()
   */
  free(occ);
  return eslOK;

 ERROR:
  if (occ != NULL) free(occ);
  return status;
}
                

/* Function:  p7_oprofile_ReconfigLength()
 * Synopsis:  Set the target sequence length distribution of a model.
 * Incept:    SRE, Mon Jun 25 09:20:44 2007 [Janelia]
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
p7_oprofile_ReconfigLength(P7_OPROFILE *om, int L)
{
  float nj;
  float pmove, ploop;
  float loopsc, movesc;
  
  /* Figure out the expected number of uses of the J state.
   */
  nj = 0.;			/* FIXME: this assumes loop p = 0: unilocal only. */

  /* configure N,J,C transitions so they bear L/(2+nj) of the total
   * unannotated sequence length L. 
   */
  pmove = (2. + nj) / ((float) L + 2. + nj); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1. - pmove;
  movesc = log(pmove);
  loopsc = log(ploop);
  om->xsc[p7X_XTN][p7X_LOOP] = loopsc;
  om->xsc[p7X_XTN][p7X_MOVE] = movesc;
  om->xsc[p7X_XTC][p7X_LOOP] = loopsc;
  om->xsc[p7X_XTC][p7X_MOVE] = movesc;
  om->xsc[p7X_XTJ][p7X_LOOP] = loopsc;
  om->xsc[p7X_XTJ][p7X_MOVE] = movesc;
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
#define MSC(k)   (rsc[(k) * p7X_NR + p7X_MSC])
#define ISC(k)   (rsc[(k) * p7X_NR + p7X_ISC])

/* Function:  p7_Viterbi()
 * Synopsis:  Viterbi algorithm.
 * Incept:    SRE, Mon Jun 25 09:23:51 2007 [Janelia]
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
p7_Viterbi(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  float const *tsc  = om->tsc;
  float      **dp   = ox->dp;
  float       *xmx  = ox->xmx;
  int          M    = om->M;
  int     i,k;


  /* Initialization of the zero row.  */
  XMX(0,p7X_XMN) = 0;                                             /* S->N, p=1            */
  XMX(0,p7X_XMB) = om->xsc[p7X_XTN][p7X_MOVE];                    /* S->N->B, no N-tail   */
  XMX(0,p7X_XME) = XMX(0,p7X_XMC) = XMX(0,p7X_XMJ) = -eslINFINITY;  /* need seq to get here */
  for (k = 0; k <= om->M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = -eslINFINITY;               /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = impossible for all eight transitions (no node 0)
   *    I_M is wastefully calculated (doesn't exist)
   *    D_M is wastefully calculated (provably can't appear in a Viterbi path)
   */
  for (i = 1; i <= L; i++) 
    {
      float const *rsc = om->rsc[dsq[i]];
      float sc;

      MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;
      XMX(i,p7X_XME) = -eslINFINITY;
    
      for (k = 1; k < om->M; k++) 
	{
  	  /* match state */
	  sc =             MMX(i-1,k-1)     + TSC(p7X_TMM,k-1);
	  sc = ESL_MAX(sc, IMX(i-1,k-1)     + TSC(p7X_TIM,k-1));
	  sc = ESL_MAX(sc, DMX(i-1,k-1)     + TSC(p7X_TDM,k-1));
	  sc = ESL_MAX(sc, XMX(i-1,p7X_XMB) + TSC(p7X_BSC,k-1));
	  MMX(i,k) = sc + MSC(k); 
	  
	  /* E state update */
	  XMX(i,p7X_XME) = ESL_MAX(XMX(i,p7X_XME), MMX(i,k));

	  /* insert state */
	  sc =             MMX(i-1,k) + TSC(p7X_TMI,k);
	  sc = ESL_MAX(sc, IMX(i-1,k) + TSC(p7X_TII,k));
	  IMX(i,k) = sc + ISC(k);
	  
	  /* delete state */
	  sc =                    MMX(i,k-1) + TSC(p7X_TMD,k-1);
	  DMX(i,k) =  ESL_MAX(sc, DMX(i,k-1) + TSC(p7X_TDD,k-1));
	}

      /* Unrolled match state M. */
      sc =             MMX(i-1,M-1)     + TSC(p7X_TMM,M-1);
      sc = ESL_MAX(sc, IMX(i-1,M-1)     + TSC(p7X_TIM,M-1));
      sc = ESL_MAX(sc, DMX(i-1,M-1 )    + TSC(p7X_TDM,M-1));
      sc = ESL_MAX(sc, XMX(i-1,p7X_XMB) + TSC(p7X_BSC,M-1));
      MMX(i,M) = sc + MSC(M);
      
      /* E state update */
      XMX(i,p7X_XME) = ESL_MAX(XMX(i,p7X_XME), MMX(i,M));
   
      /* Now the special states. Order is important here.
       * remember, N, C and J emissions are zero score by definition.
       */
      /* J state */
      sc             =             XMX(i-1,p7X_XMJ) + om->xsc[p7X_XTJ][p7X_LOOP];
      XMX(i,p7X_XMJ) = ESL_MAX(sc, XMX(i,  p7X_XME) + om->xsc[p7X_XTE][p7X_LOOP]);
      
      /* C state */
      sc =                         XMX(i-1,p7X_XMC) + om->xsc[p7X_XTC][p7X_LOOP];
      XMX(i,p7X_XMC) = ESL_MAX(sc, XMX(i,  p7X_XME) + om->xsc[p7X_XTE][p7X_MOVE]);
      
      /* N state */
      XMX(i,p7X_XMN) = XMX(i-1,p7X_XMN) + om->xsc[p7X_XTN][p7X_LOOP];
      
      /* B state */
      sc             =             XMX(i,p7X_XMN) + om->xsc[p7X_XTN][p7X_MOVE];
      XMX(i,p7X_XMB) = ESL_MAX(sc, XMX(i,p7X_XMJ) + om->xsc[p7X_XTJ][p7X_MOVE]);
    }
  
  /* T state (not stored) */
  *ret_sc = XMX(L,p7X_XMC) + om->xsc[p7X_XTC][p7X_MOVE];
  return eslOK;
}




/*****************************************************************
 * Benchmark driver
 *****************************************************************/
#ifdef p7IMPL_FP_BENCHMARK
/* gcc -o benchmark-fp -g -O2 -I. -L. -I../easel -L../easel -Dp7IMPL_FP_BENCHMARK impl_fp.c -lhmmer -leasel -lm
 * ./benchmark-fp <hmmfile>
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
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual scores",             0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the floating-point test implementation";

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
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc, nullsc;
  float           bitscore;
  double          x;
  

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  om = p7_oprofile_Create(hmm->M, abc);
  p7_oprofile_Config(hmm, bg, om, p7_UNILOCAL);
  p7_oprofile_ReconfigLength(om, L);

  ox = p7_omx_Create(om->M, L);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, bg->f, abc->K, L, dsq);

      p7_Viterbi(dsq, L, om, ox, &sc);
      p7_bg_NullOne(bg, dsq, L, &nullsc);
      bitscore = (sc - nullsc) / eslCONST_LOG2;
      
      if (esl_opt_GetBoolean(go, "-v")) printf("%.4f bits  (%.4f raw)\n", bitscore, sc); 
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  free(dsq);
  p7_omx_Destroy(ox);
  p7_oprofile_Destroy(om);
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
