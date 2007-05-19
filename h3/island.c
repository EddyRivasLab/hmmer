/* Implementation of the "island method" for estimating extreme value
 * distribution parameters. 
 * 
 * Contents:
 *    1. The exposed API.
 *    2. Statistics collection driver.
 *    3. Copyright stamp.
 * 
 * References: 
 *  [Altschul01] SF Altschul, R Bundschuh, R Olsen, T Hwa. The estimation of
 *               statistical parameters for local alignment score
 *               distributions. Nucleic Acids Research 29:351-361, 2001.
 *               
 *  [Olsen99]    R Olsen, R Bundschuh, T Hwa. Rapid assessment of extremal
 *               statistics for gapped local alignment. Proc Int Conf
 *               Intell Syst Mol Biol 1999; 211-222.            
 *               
 * 
 * Original implementation: STL9/p75; SRE, Thu Apr 21 09:57:21 2005
 * Adaptation to H3:        J1/22;    SRE, Fri Mar 23 10:34:53 2007 [Janelia]
 * SVN $Id$
 */

#include "p7_config.h"

#include <easel.h>
#include <esl_histogram.h>

#include "hmmer.h"


/* Function: p7_island_Viterbi()
 * Date:     SRE, Thu Apr 21 11:23:36 2005 [St. Louis]
 *
 * Purpose:  An implementation of the Olsen/Bundschuh/Hwa island
 *           algorithm \citep{OlsenHwa99}, extended to Plan 7 profile HMMs.
 *           
 *           Every S->N^n->B->M entry is counted as an independent island
 *           beginning; every M->E->C^n->T is counted as an
 *           independent island ending. But an island may still
 *           contain multiple hits to the model, passing through the J
 *           state. 
 *
 *           Caller provides an allocated DP matrix <mx> with at least
 *           <L=4> rows and <M=gm->M> columns. Internally, we use
 *           two-row, rolling pointer, linear-memory alignment. The
 *           reason we need 4 rows instead of 2 is because we have two
 *           different types of storage: scores and island tags.
 *
 *           Caller provides a histogram <h> to collect the island
 *           scores. A normal histogram will collect only binned data;
 *           a full histogram will collect the individual data points
 *           as well.
 *
 * Args:     dsq          - sequence in digitized form
 *           L            - length of dsq
 *           gm           - score profile
 *           mx           - DP matrix allocated for at least L=4 rows, M=gm->M columns.
 *           h            - histogram that will collect island scores.
 *
 * Returns:  <eslOK> on success.
 * 
 * Xref:     STL9/75.
 */
int
p7_island_Viterbi(ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_GMX *mx, ESL_HISTOGRAM *h)
{
  int  status;
  int  **xmx, **mmx, **dmx, **imx; /* convenience ptrs to score matrix   */  
  int  **xi,  **mi,  **di,  **ii;  /* convenience ptrs to island tags    */

  int   *I = NULL;	/* island scores, 0..inum-1 */
  void  *tmp;
  int    ialloc;	/* current allocation of I */
  int    inum;		/* current # of islands, in I[0..inum-1] */
  
  int    sc;			/* integer score of optimal alignment  */
  int    i,k;		        /* index for seq, model position */
  int    cur, prv;		/* indices for rolling dp matrix */
  int    nullsc;
  float  bitscore;

  /* Reach into the DP matrix memory, set our convenience ptrs.
   */
  xmx = mx->xmx;		/* use only two rows, xmx[0..1] */
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;
  xi  = mx->xmx+2;		/* yes, this works. xi[0] now points to mx->xmx[2] row */
  mi  = mx->mmx+2;
  ii  = mx->imx+2;
  di  = mx->dmx+2;
  
  /* Initial allocation of the island score vector.
   */
  ialloc = L;
  inum   = 0;
  ESL_ALLOC(I, sizeof(int) * ialloc);

  /* Initialization of the zero row.
   */
  xmx[0][p7_XMN] = 0;		                                     /* S->N, p=1            */
  xmx[0][p7_XMB] = gm->xsc[p7_XTN][p7_MOVE];                        /* S->N->B, no N-tail   */
  xmx[0][p7_XME] = xmx[0][p7_XMC] = xmx[0][p7_XMJ] = p7_IMPOSSIBLE;  /* need seq to get here */
  for (k = 0; k <= gm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = p7_IMPOSSIBLE;               /* need seq to get here too*/
  xi[0][p7_XMB]  = -1;

  /* Recursion. Done as a pull, with rolling index trick.
   *    
   * Notes on start position and island tag propagation.
   *  - In the path B->E, we propagate the i that B was aligned to 
   *    in the optimal alignment, via mtr, dtr, and itr. 
   *  - When we reach an E, we record the i of the B it started from in etr.
   *  - In a looping path E->J...->B or terminal path E->C...->T, we propagate
   *    the i that E was aligned to in the optimal alignment via xtr[][XMC]
   *    and xtr[][XMJ].
   *  - When we enter B, we record the i of the best previous E, or 0 if there
   *    isn't one, in btr.
   */
  for (i = 1; i <= L; i++) {
    cur = i % 2;
    prv = !cur;
    mmx[cur][0] = imx[cur][0] = dmx[cur][0] = p7_IMPOSSIBLE;

    for (k = 1; k <= gm->M; k++) {
				/* match state transitions*/
      mmx[cur][k] = p7_IMPOSSIBLE;
      if ((sc = mmx[prv][k-1] + gm->tsc[p7_TMM][k-1]) > mmx[cur][k]) { 
	mmx[cur][k] = sc; 
	mi[cur][k]  = mi[prv][k-1]; /* propagate the island tag */
      }
      if ((sc = imx[prv][k-1] + gm->tsc[p7_TIM][k-1]) > mmx[cur][k]) {
	mmx[cur][k] = sc;
	mi[cur][k]  = ii[prv][k-1];
      }
      if ((sc = dmx[prv][k-1] + gm->tsc[p7_TDM][k-1]) > mmx[cur][k]) {
	mmx[cur][k] = sc;
	mi[cur][k]  = di[prv][k-1];
      }
      if ((sc = xmx[prv][p7_XMB] + gm->bsc[k]) > mmx[cur][k]) {
	mmx[cur][k] = sc; 

	if (xi[prv][p7_XMB] == -1) {
	  mi[cur][k]  = inum;	/* start a new island, if B->M is the best path. */

	  I[inum] = p7_IMPOSSIBLE;
	  inum++;			
	  if (inum == ialloc) { /* reallocation needed? */
	    ialloc *= 2;
	    ESL_RALLOC(I, tmp, sizeof(int) * ialloc);
	  }
	} else {
	  mi[cur][k] = xi[prv][p7_XMB]; /* propagating a loop */
	}
      }

				/* match emissions */
      if (mmx[cur][k] != p7_IMPOSSIBLE && gm->msc[dsq[i]][k] != p7_IMPOSSIBLE)
	mmx[cur][k] += gm->msc[dsq[i]][k];
      else
	mmx[cur][k] = p7_IMPOSSIBLE;

				/* delete state */
      dmx[cur][k] = p7_IMPOSSIBLE;
      if ((sc = mmx[cur][k-1] + gm->tsc[p7_TMD][k-1]) > dmx[cur][k]) {
	dmx[cur][k] = sc; 
	di[cur][k]  = mi[cur][k-1];
      }
      if ((sc = dmx[cur][k-1] + gm->tsc[p7_TDD][k-1]) > dmx[cur][k]) {
	dmx[cur][k] = sc; 
	di[cur][k]  = di[cur][k-1];
      }

				/* insert state */
      if (k < gm->M) {
	imx[cur][k] = p7_IMPOSSIBLE;
	if ((sc = mmx[prv][k] + gm->tsc[p7_TMI][k]) > imx[cur][k]) {
	  imx[cur][k] = sc; 
	  ii[cur][k]  = mi[prv][k];
	}
	if ((sc = imx[prv][k] + gm->tsc[p7_TII][k]) > imx[cur][k]) {
	  imx[cur][k] = sc; 
	  ii[cur][k]  = ii[prv][k];
	}

	if (imx[cur][k] != p7_IMPOSSIBLE && gm->isc[dsq[i]][k] != p7_IMPOSSIBLE)
	  imx[cur][k] += gm->isc[dsq[i]][k];
	else
	  imx[cur][k] = p7_IMPOSSIBLE;
      }
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition.
     */
				/* N state */
    xmx[cur][p7_XMN] = p7_IMPOSSIBLE;
    if ((sc = xmx[prv][p7_XMN] + gm->xsc[p7_XTN][p7_LOOP]) > p7_IMPOSSIBLE)
      xmx[cur][p7_XMN] = sc;
				/* E state */
    xmx[cur][p7_XME] = p7_IMPOSSIBLE;
    for (k = 1; k <= gm->M; k++)
      {
	if ((sc =  mmx[cur][k] + gm->esc[k]) > xmx[cur][p7_XME]) {
	  xmx[cur][p7_XME] = sc;
	  xi[cur][p7_XME]  = mi[cur][k];
	}
	/* calculate what island sc would be, if we ended it here. */
	sc = sc  + gm->xsc[p7_XTE][p7_MOVE] + (L-i)* gm->xsc[p7_XTC][p7_LOOP] + gm->xsc[p7_XTC][p7_MOVE];
	if (sc > I[mi[cur][k]])	  
	  I[mi[cur][k]]  = sc;
      }
				/* J state */
    xmx[cur][p7_XMJ] = p7_IMPOSSIBLE;
    if ((sc = xmx[prv][p7_XMJ] + gm->xsc[p7_XTJ][p7_LOOP]) > p7_IMPOSSIBLE) {
      xmx[cur][p7_XMJ] = sc; 
      xi[cur][p7_XMJ]  = xi[prv][p7_XMJ];
    }
    if ((sc = xmx[cur][p7_XME] + gm->xsc[p7_XTE][p7_LOOP]) > xmx[cur][p7_XMJ]) {
      xmx[cur][p7_XMJ] = sc;
      xi[cur][p7_XMJ]  = xi[cur][p7_XME];
    }
				/* B state */
    xmx[cur][p7_XMB] = p7_IMPOSSIBLE;
    if ((sc = xmx[cur][p7_XMN] + gm->xsc[p7_XTN][p7_MOVE]) > p7_IMPOSSIBLE) {
      xmx[cur][p7_XMB] = sc; 
      xi[cur][p7_XMB]  = -1;       /* if coming from N, we are islandless */
    }
    if ((sc = xmx[cur][p7_XMJ] + gm->xsc[p7_XTJ][p7_MOVE]) > xmx[cur][p7_XMB]) {
      xmx[cur][p7_XMB] = sc; 
      xi[cur][p7_XMB]  = xi[cur][p7_XMJ]; /* if from J, then propagate island tag */
    }
				/* C state */
    xmx[cur][p7_XMC] = p7_IMPOSSIBLE;
    if ((sc = xmx[prv][p7_XMC] + gm->xsc[p7_XTC][p7_LOOP]) > p7_IMPOSSIBLE)
      xmx[cur][p7_XMC] = sc; 
    if ((sc = xmx[cur][p7_XME] + gm->xsc[p7_XTE][p7_MOVE]) > xmx[cur][p7_XMC])
      xmx[cur][p7_XMC] = sc; 
  }
				/* T state (not stored) */
  sc = xmx[cur][p7_XMC] + gm->xsc[p7_XTC][p7_MOVE];
  /* sc is the overall optimal score, but we don't do anything with it here. */


  /* Count the maximal island scores into the histogram.
   */
  p7_bg_NullOne(gm->bg, dsq, L, &nullsc);
  for (i = 0; i < inum; i++)
    {
      if (I[i] > p7_IMPOSSIBLE) {
	bitscore = p7_SILO2Bitscore(I[i] - nullsc);
	esl_histogram_Add(h, bitscore);
      }
    }
  free(I);
  return eslOK;

 ERROR:
  if (I != NULL) free(I);
  return status;
}


/*****************************************************************
 * 2. Statistics collection driver.
 *****************************************************************/


#ifdef p7ISLAND_STATS
/* gcc -g -O2 -Dp7ISLAND_STATS -I. -I../easel -L. -L../easel -o statprog island.c -lhmmer -leasel -lm
 * ./statprog
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_histogram.h"
#include "esl_getopts.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"
#include "esl_dmatrix.h"
#include "esl_ratematrix.h"
#include "esl_exponential.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   1 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0",     NULL,      NULL,    NULL, "length of random target seqs",           1 },
  { "-N",        eslARG_INT,      "1", NULL, "n>0",     NULL,      NULL,    NULL, "number of random target seqs",           1 },
  { "-2",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "configure profile in old HMMER2 style",  1 },

  { "-m",        eslARG_INFILE,  NULL, NULL, NULL,      NULL,      NULL,    NULL, "input HMM from file <f>",                2 },
  { "-S",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "sample a sequence and make HMM from it", 2 },

  { "-M",        eslARG_INT,     "50", NULL, "n>0",     NULL,      NULL,    "-m", "length of a sampled query HMM or seq",         3 },
  { "-t",        eslARG_REAL,   "2.0", NULL, "x>0",     NULL,      "-S",    NULL, "branch length, for parameterizing seq-based query", 3 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "statprog [-options]";

int
main(int argc, char **argv)
{
  int status;
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);
  ESL_RANDOMNESS *r   = esl_randomness_CreateTimeseeded();
  ESL_HISTOGRAM  *h   = esl_histogram_CreateFull(-50., 50., 0.5);
  ESL_DSQ        *dsq = NULL;
  P7_BG          *bg  = p7_bg_Create(abc);
  P7_HMM         *hmm = NULL;
  P7_PROFILE     *gm  = NULL;
  P7_GMX         *vmx = NULL;

  int nseq;			/* counter over sequences */

  char            *hmmfile;            /* file to read HMM(s) from                */
  int              do_seqquery;        /* TRUE to make model from random sequence */
  int              L;		       /* length of generated seqs                */
  int              M;		       /* length of a sampled HMM                 */
  int              N;		       /* number of seqs to generate              */
  double           t;		       /* branch length for seq query model       */
  int              do_h2;              /* TRUE to use H2 exit/entry configuration */

  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  ESL_GETOPTS    *go  = esl_getopts_Create(options);

  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_GetBoolean(go, "-h")) {
    puts(usage);
    puts("\ngeneral options are:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 2 = indentation; 80=textwidth*/
    puts("\nalternatives for random query model :");
    esl_opt_DisplayHelp(stdout, go, 2, 2, 80); /* 2 = indentation; 80=textwidth*/
    puts("\ncontrolling character of random query model :");
    esl_opt_DisplayHelp(stdout, go, 3, 2, 80); /* 2 = indentation; 80=textwidth*/
    return eslOK;
  }

  hmmfile     = esl_opt_GetString (go, "-m");
  do_seqquery = esl_opt_GetBoolean(go, "-S");
  L           = esl_opt_GetInteger(go, "-L");
  M           = esl_opt_GetInteger(go, "-M");
  N           = esl_opt_GetInteger(go, "-N");
  do_h2       = esl_opt_GetBoolean(go, "-2");
  t           = esl_opt_GetReal   (go, "-t");

  if (esl_opt_ArgNumber(go) != 0) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return eslFAIL;
  }
  esl_getopts_Destroy(go);


  /*****************************************************************
   * Create the model we're going to use.
   *****************************************************************/

  if (hmmfile != NULL)		/* Read in an HMM from file... */
    {
      P7_HMMFILE      *hfp     = NULL;    

      status = p7_hmmfile_Open(hmmfile, NULL, &hfp);
      if (status == eslENOTFOUND) esl_fatal("Failed to open hmm file %s for reading.\n", hmmfile);
      else if (status != eslOK)   esl_fatal("Unexpected error in opening hmm file %s.\n", hmmfile);

      if ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslOK) {
	if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
	else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s", hmmfile);
	else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets", hmmfile);
	else                             esl_fatal("Unexpected error in reading HMMs");
      }
      M = hmm->M;
      p7_hmmfile_Close(hfp);
    }
  else if (do_seqquery)
    {
      double   pi[20];
      ESL_DMATRIX *Q     = esl_dmatrix_Create(abc->K,abc->K);
      ESL_DMATRIX *P     = esl_dmatrix_Create(abc->K,abc->K);
      ESL_DSQ     *query = malloc(sizeof(ESL_DSQ) * (M+2));

      esl_vec_F2D(bg->f, 20, pi);
      esl_rmx_SetWAG(Q, pi);
      esl_dmx_Exp(Q, t, P);

      printf("expected score   of WAG at t=%f  is %f bits\n", t, esl_rmx_ExpectedScore  (P, pi));
      printf("relative entropy of WAG at t=%f  is %f bits\n", t, esl_rmx_RelativeEntropy(P, pi));

      esl_rnd_xfIID(r, bg->f, abc->K, M, query);
      p7_Seqmodel(abc, query, M, P, bg->f, 0.05, 0.5, 0.05, 0.2, &hmm); /* tmi, tii, tmd, tdd */
      
      esl_dmatrix_Destroy(Q);
      esl_dmatrix_Destroy(P);
      free(dsq);
    }
  else				/* or generate a simulated HMM. */
    {
      p7_hmm_Sample(r, M, abc, &hmm);
    }


  /*****************************************************************
   * Remaining initializations.
   *****************************************************************/

  vmx     = p7_gmx_Create(hmm->M, 4); /* 4 rows needed for island method. */
  gm      = p7_profile_Create(hmm->M, abc);
  hmm->bg = bg;
  hmm->gm = gm;
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));

  if (do_h2) {
    p7_H2_ProfileConfig(hmm, gm, p7_UNILOCAL);
  } else {
    p7_ProfileConfig(hmm, gm, p7_UNILOCAL);
    p7_ReconfigLength(hmm->gm,  L);
  }

  /*****************************************************************
   * Run the island method on random target sequences.
   *****************************************************************/

  for (nseq = 1; nseq <= N; nseq++) 
    {
      esl_rnd_xfIID(r, hmm->bg->f, abc->K, L, dsq);
      p7_island_Viterbi(dsq, L, gm, vmx, h);
    }
  
  /*****************************************************************
   * Output to xmgrace file.
   *****************************************************************/

  double *xv;
  int     n;
  double  mu, lambda;
  double  param[2];

  esl_histogram_GetTailByMass(h, 0.01, &xv, &n, NULL);
  esl_exp_FitComplete(xv, n, &mu, &lambda);
  param[0] = mu;
  /*param[1] = lambda; */
  param[1] = 0.693; 

  esl_histogram_SetExpectedTail(h, mu, 0.01, &esl_exp_generic_cdf, &param);
  esl_histogram_PlotSurvival(stdout, h);

  fprintf(stderr, "lambda = %f\n", lambda);

  free(dsq);
  p7_gmx_Destroy(vmx);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  return 0;

 ERROR:
  return status;
}


#endif /* p7ISLAND_STATS */




/************************************************************
 * @LICENSE@
 ************************************************************/

