/* Forward/Backward, checkpointed: generic (non-SIMD) demonstration version.
 * 
 * Contents:
 *    1. Forwards:  checkpointed fill, Forwards nat score.
 *    2. Backwards: linear-memory back pass, recovering posterior-decoded bands
 *    3. Benchmark driver
 *    4. Example main().
 *    5. Copyright and license information.
 */
#include "p7_config.h"

#include "easel.h"

#include "hmmer.h"
#include "p7_gmxchk.h"

/*****************************************************************
 *= 1. Forwards: checkpointed fill, Forwards nat score
 *****************************************************************/

/* forward_row()
 * 
 * inlined function, because we have to call the row calculation in
 * several places in a checkpointed calculation, including in the Backwards()
 * routine.
 * 
 */
static inline void
forward_row(const ESL_DSQ *dsq, const P7_PROFILE *gm, P7_GMXCHK *gxc, const float *dpp, float *dpc, int i)
{
  float const *tsc = gm->tsc;
  float const *rsc = gm->rsc[dsq[i]];
  float       *xmx = gxc->xmx;
  int            M = gm->M;
  float        esc = p7_profile_IsLocal(gm) ? 0. : -eslINFINITY;
  int            k;
  float         sc;

  MMR(dpc,0) = IMR(dpc,0) = DMR(dpc,0) = -eslINFINITY;
  XMX(i,p7G_E) = -eslINFINITY;
  for (k = 1; k < M; k++)
    {
      /* match state */
      sc = p7_FLogsum(p7_FLogsum(MMR(dpp,k-1)   + TSC(p7P_MM,k-1), 
				 IMR(dpp,k-1)   + TSC(p7P_IM,k-1)),
		      p7_FLogsum(XMX(i-1,p7G_B) + TSC(p7P_BM,k-1),
				 DMR(dpp,k-1)   + TSC(p7P_DM,k-1)));
      MMR(dpc,k) = sc + MSC(k);

      /* insert state */
      sc = p7_FLogsum(MMR(dpp,k) + TSC(p7P_MI,k),
		      IMR(dpp,k) + TSC(p7P_II,k));
      IMR(dpc,k) = sc + ISC(k);

      /* delete state */
      DMR(dpc,k) = p7_FLogsum(MMR(dpc,k-1) + TSC(p7P_MD,k-1),
			      DMR(dpc,k-1) + TSC(p7P_DD,k-1));

      /* E state update */
      XMX(i,p7G_E) = p7_FLogsum(p7_FLogsum(MMR(dpc,k) + esc,
					   DMR(dpc,k) + esc),
				XMX(i,p7G_E));
    }

  /* unrolled match state M_M */
  sc = p7_FLogsum(p7_FLogsum(MMR(dpp,M-1)   + TSC(p7P_MM,M-1), 
			     IMR(dpp,M-1)   + TSC(p7P_IM,M-1)),
		  p7_FLogsum(XMX(i-1,p7G_B) + TSC(p7P_BM,M-1),
			     DMR(dpp,M-1)   + TSC(p7P_DM,M-1)));
  MMR(dpc,M) = sc + MSC(M);
  IMR(dpc,M) = -eslINFINITY;

  /* unrolled delete state D_M */
  DMR(dpc,M) = p7_FLogsum(MMR(dpc,M-1) + TSC(p7P_MD,M-1),
			  DMR(dpc,M-1) + TSC(p7P_DD,M-1));

  /* unrolled E state update */
  XMX(i,p7G_E) = p7_FLogsum(p7_FLogsum(MMR(dpc,M),
	 		               DMR(dpc,M)),
				       XMX(i,p7G_E));

  /* J state */
  XMX(i,p7G_J) = p7_FLogsum(XMX(i-1,p7G_J) + gm->xsc[p7P_J][p7P_LOOP],
	                    XMX(i,  p7G_E) + gm->xsc[p7P_E][p7P_LOOP]);

  /* C state */
  XMX(i,p7G_C) = p7_FLogsum(XMX(i-1,p7G_C) + gm->xsc[p7P_C][p7P_LOOP],
 			    XMX(i,  p7G_E) + gm->xsc[p7P_E][p7P_MOVE]);
  /* N state */
 XMX(i,p7G_N) = XMX(i-1,p7G_N) + gm->xsc[p7P_N][p7P_LOOP];

 /* B state */
 XMX(i,p7G_B) = p7_FLogsum(XMX(i,  p7G_N) + gm->xsc[p7P_N][p7P_MOVE],
			   XMX(i,  p7G_J) + gm->xsc[p7P_J][p7P_MOVE]);
}
  

/* Function:  p7_GForwardCheckpointed()
 * Synopsis:  Forward pass in a checkpointed generic DP matrix
 *
 * Purpose:   Compute the Forward pass of a comparison of model <gm>
 *            against digital sequence <dsq> of length <L>, resulting
 *            in a filled DP matrix <gxc> and (optionally)
 *            a Forward score <opt_sc> in nats.
 *            
 *            The caller has already allocated and laid out <gxc>
 *            appropriately for the <gm->M> by <L> comparison, either
 *            with <p7_gmxchk_Create()> or <p7_gmxchk_GrowTo()>.
 *            
 *            The caller has also already configured the length model
 *            in <gm> for the target sequence length <L>, for example
 *            by calling <p7_ReconfigLength()>.
 *            
 * Args:      dsq    : digital sequence target
 *            L      : length of <dsq> in residues
 *            gm     : query profile
 *            gxc    : checkpointed DP matrix to fill
 *            opt_sc : optRETURN: Forward raw lod score, in nats
 *
 * Returns:   <eslOK> on success. <gxc> contains checkpointed DP
 *            Forward matrix, ready for backwards pass. <opt_sc>
 *            is the raw lod score in nats.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_GForwardCheckpointed(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXCHK *gxc, float *opt_sc)
{
  float       *xmx = gxc->xmx;
  float *dpc;			/* ptr to current row  */
  float *dpp;			/* ptr to previous row */
  int    M   = gm->M;
  int    w;			/* size of a segment of DP rows ended w/ one checkpointed row, inclusive */
  int    b;			/* counts down over segment number (Rb+Rc..1); each segment has r+1 rows */
  int    i,k;

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
      dpc = gxc->dp[gxc->R0+gxc->R]; gxc->R++; /* idiomatic for "get next saved row" */
      forward_row(dsq, gm, gxc, dpp, dpc, i);
      dpp = dpc;
    }

  /* Phase two: "b" and "c" regions: partially and fully checkpointed */
  /* i= gxc->La+1, from previous loop's end */
  for (b = gxc->Rb + gxc->Rc, w = (gxc->Rb ? gxc->Lb : gxc->Rc+1); i <= L; i++)
    {
      if (! (--w))
	{ 		                           /* we're on the last row in segment: this row is saved    */
	  dpc = gxc->dp[gxc->R0+gxc->R]; gxc->R++; /* idiomatic for "get next saved row"                     */
	  w=b;  			           /* next segment has this many rows, ending in a saved row */
	  b--;					   /* decrement segment number counter; last segment is r=1  */
	}
      else
	dpc =  gxc->dp[i%2];	/* idiomatic for "get next temp row", 0 or 1; i%2 is to make sure dpp != dpc  */
      
      forward_row(dsq, gm, gxc, dpp, dpc, i);
      dpp = dpc;
    }
  
  gxc->M = M;
  gxc->L = L;
  gxc->R = gxc->Ra + gxc->Rb + gxc->Rc;
  if (opt_sc) *opt_sc = XMX(L,p7G_C) + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}
/*-------------------- end, forwards ----------------------------*/


/*****************************************************************
 *= 2. Backwards: linear-memory back pass, recovering posterior-decoded bands
 *****************************************************************/

#if 0				/* WORK IN PROGRESS */
int
p7_GBackwardCheckpointed(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXCHK *gxc)
{
  float *bck;

  /* We have to handle the first block b=1 as a special case (rows L-1, L)
   * because row L is a special case.
   */
  
  /* Initialize backwards row L */
  


  /* init i=L row for backwards, in dpc */
  dpp = dpc;

  for (i = L, b = 1; b <= gxc->Rb+gxc->Rc; b++)
    {
      w = (b <= gxc->Rc ? b+1 : gxc->Lb);

      /* current row i (r=R0+R-1) ends a block and is checkpointed. */
      gxc->R--;
      fwd = gxc->dp[gxc->R0+gxc->R]; /* pop last forward row off "stack" */
      bck = gxc->dp[i%2];
      backward_row(bck, dpp);
      posterior_decode(fwd, bck, &band);
      
      /* compute Forwards from last checkpoint */
      dpp = gxc->dp[gxc->R0+gxc->R-1];
      for (i2 = i-w+1; i2 <= i-1; i2++)
	{
	  fwd = gxc->dp[gxc->R0+gxc->R]; gxc->R++; /* push new forward row on "stack" */
	  forward_row(dsq, gm, gxc, dpp, fwd, i);
	  dpp = fwd;	  
	}

      /* now compute Backwards over the block we just calculated */
      dpp = bck;
      for (i2 = i-1; i2 >= i-w+1; i2--)
	{
	  gxc->R--;
	  fwd = gxc->dp[gxc->R0+gxc->R]; /* pop last forward row off "stack" */
	  bck = gxc->dp[i2%2];

	  backward_row(bck, dpp);
	  posterior_decode(fwd, bck);

	  dpp = bck;
	}

    }
  /* now i=La */

  for (; i >= 1; i--)
    {
      gxc->R--; 
      fwd = gxc->dp[gxc->R0+gxc->R];
      bck = gxc->dp[i%2];

      backward_row(bck, i);
      posterior_decoding(fwd, bck);

      dpp = bck;
    }
#endif /*WORK IN PROGRESS*/

/*--------------------- end, backwards --------------------------*/



/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_CHK_BENCHMARK
/*
   gcc -g -O2      -o generic_fwdback_chk_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_FWDBACK_CHK_BENCHMARK generic_fwdback_chk.c -lhmmer -leasel -lm
   icc -O3 -static -o generic_fwdback_chk_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_FWDBACK_CHK_BENCHMARK generic_fwdback_chk.c -lhmmer -leasel -lm
   ./generic_fwdback_chk_benchmark <hmmfile>
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "p7_gmxchk.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,   "2000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for checkpointed generic Forward/Backward";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMXCHK      *gxc     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL);

  gxc = p7_gmxchk_Create(gm->M, L, 32);

  /* Baseline time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_GForwardCheckpointed (dsq, L, gm, gxc, &sc);
      
      p7_gmxchk_Reuse(gxc);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_gmxchk_Destroy(gxc);
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
#endif /*p7GENERIC_FWDBACK_CHK_BENCHMARK*/
/*----------------- end, benchmark ------------------------------*/


/*****************************************************************
 * 4. Example main()
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
  float           fsc, fsc2;
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
  if      (status == eslENOTFOUND) p7_Fail("No such file, or file open failed");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, 400, p7_LOCAL);

  /* Allocate matrices */
  fwd  = p7_gmx_Create(gm->M, 400);
  fwdc = p7_gmxchk_Create(gm->M, 400, 32);

  while ( (status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
      else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

      /* Resize the DP matrices if necessary */
      p7_gmx_GrowTo   (fwd,  gm->M, sq->n);
      p7_gmxchk_GrowTo(fwdc, gm->M, sq->n);

      //printf("Allocation: %ld\n", p7_gmxchk_Sizeof(fwdc));

      /* Set the profile and null model's target length models */
      p7_bg_SetLength(bg,   sq->n);
      p7_ReconfigLength(gm, sq->n);

      /* Run Forward in both modes */
      p7_GForward            (sq->dsq, sq->n, gm, fwd,  &fsc);
      p7_GForwardCheckpointed(sq->dsq, sq->n, gm, fwdc, &fsc2);

      /* Dump the DP matrices. (Voluminous; only small examples are reasonable) */
      //p7_gmx_Dump(stdout,    fwd,  p7_DEFAULT);
      //p7_gmxchk_Dump(stdout, fwdc, p7_DEFAULT);

      /* Those scores are partial log-odds likelihoods in nats.
       * Subtract off the rest of the null model, convert to bits.
       */
      p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

      printf("%-30s   %10.4f %10.4f   %10.4f %10.4f\n", 
	     sq->name, 
	     fsc, fsc2, 
	     (fsc - nullsc) / eslCONST_LOG2, (fsc2 - nullsc) / eslCONST_LOG2);

      p7_gmx_Reuse(fwd);
      p7_gmxchk_Reuse(fwdc);
      esl_sq_Reuse(sq);
    }

  /* Cleanup */
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_gmx_Destroy(fwd);
  p7_gmxchk_Destroy(fwdc);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_FWDBACK_CHK_EXAMPLE*/
/*------------------- end, example main() -----------------------*/
      

/*****************************************************************
 * @LICENSE@
 *    
 * References:
 *    SRE J8/109-112, Oct 2011: implementation plan.
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
