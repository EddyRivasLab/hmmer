/* Forward/Backward implementation variant:
 *   dual-mode, local/glocal alignment;
 *   quadratic memory (simplest variant; not banded, not checkpointed);
 *   "generic" (standard C code; not striped/vectorized);
 *   using P7_GMXD DP matrix structure.
 *   
 * Contents:  
 *   1. Forward, Backward implementations.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 */


#include "p7_config.h"

#include "easel.h"

#include "hmmer.h"
#include "p7_gmxd.h"


/*****************************************************************
 * 1. Forward, Backward implementations
 *****************************************************************/

/* Function:  
 * Synopsis:  
 *
 * Purpose:   
 * 
 *            This function uses <p7_FLogsum()>. Caller must have
 *            initialized its static lookup table with a
 *            <p7_FLogsumInit()> call.
 *            
 *            
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Notes:     This function makes assumptions about the order
 *            of the state indices in p7_gmxd.h:
 *              main states: ML MG IL IG DL DG
 *              specials:    E N J B L G C
 *              
 *            If L=0, the score is -infinity, by construction; HMMER
 *            profiles generate sequences of L>=1.
 */
int
p7_GForwardDual(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXD *gxd, float *opt_sc)
{
  float *dpc, *dpp;
  const float *tsc;		/* ptr for stepping thru profile's transition parameters */
  const float *rsc;		/* ptr for stepping thru profile's emission parameters   */
  int    M          = gm->M;
  int    i, k, s;
  float  mlv, mgv;	      /* ML,MG cell values on current row   */
  float  dlv, dgv; 	      /* pushed-ahead DL,DG cell k+1 values */
  float  xE, xN, xJ, xB, xL, xG;
  
  /* Initialization of the zero row. */
  dpc = gxd->dp[0];
  for (s = 0; s < (M+1) * p7GD_NSCELLS; s++)
    *dpc++ = -eslINFINITY; 	                               // all M,I,D; k=0..M
  *dpc++ = -eslINFINITY;	                               // E
  *dpc++ = 0.0;			                               // N
  *dpc++ = -eslINFINITY;                                       // J
  *dpc++ = gm->xsc[p7P_N][p7P_MOVE];                           // B
  *dpc++ = xL = gm->xsc[p7P_N][p7P_MOVE] + gm->xsc[p7P_B][0];  // L
  *dpc++ = xG = gm->xsc[p7P_N][p7P_MOVE] + gm->xsc[p7P_B][1];  // G
  *dpc        = -eslINFINITY;                                  // C
  /* *dpc must end and stay on C state, to handle L=0 case where the recursion body below doesn't run */

  /* Main DP recursion */
  for (i = 1; i <= L; i++)
    {
      /* Initialization for a new row */
      rsc = gm->rsc[dsq[i]] + p7P_NR;	/* this ptr steps through the row's emission scores 1..M. skip k=0 */
      tsc = gm->tsc;			/* this ptr steps through profile's transition scores 0..M         */

      dpp = gxd->dp[i-1];               /* previous row dpp is already set, and at k=0 */
      dpc = gxd->dp[i];                 /* current DP row, skip k=0, start at k=1.  */
      for (s = 0; s < p7GD_NSCELLS; s++) *dpc++ = -eslINFINITY;

      dlv = dgv = -eslINFINITY;
      xE  =       -eslINFINITY;

      /* Main inner loop of the recursion */
      for (k = 1; k < M; k++)
	{
	  /* match states MLk, MGk */
	  mlv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7GD_ML) + *(tsc + p7P_MM),
						       *(dpp+p7GD_IL) + *(tsc + p7P_IM)),
					    p7_FLogsum(*(dpp+p7GD_DL) + *(tsc + p7P_DM),
						       xL             + *(tsc + p7P_BLM)));

	  mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7GD_MG) + *(tsc + p7P_MM),
						       *(dpp+p7GD_IG) + *(tsc + p7P_IM)),
 					    p7_FLogsum(*(dpp+p7GD_DG) + *(tsc + p7P_DM),
						       xG             + *(tsc + p7P_BGM)));

	  rsc++;                /* rsc advances to insert score for position k */
	  tsc += p7P_NTRANS;    /* tsc advances to transitions in states k     */
	  dpp += p7GD_NSCELLS;	/* dpp advances to cells for states k          */

	  /* Insert state calculations ILl, IGk. */
	  *dpc++ = *rsc + p7_FLogsum( *(dpp + p7GD_ML) + *(tsc + p7P_MI), *(dpp + p7GD_IL) + *(tsc + p7P_II));
	  *dpc++ = *rsc + p7_FLogsum( *(dpp + p7GD_MG) + *(tsc + p7P_MI), *(dpp + p7GD_IG) + *(tsc + p7P_II));
	  rsc++;		/* rsc advances to next match state emission   */

	  /* E state update; local paths only, DLk->E, MLk->E */
	  xE  = p7_FLogsum( p7_FLogsum(mlv, dlv), xE);

	  /* Delete state, deferred storage trick */
	  *dpc++ = dlv;
	  *dpc++ = dgv;
	  dlv = p7_FLogsum( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD));
	  dgv = p7_FLogsum( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD));
	}

      /* k=M node is unrolled and handled separately. No I state, and glocal exits. */
      mlv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7GD_ML) + *(tsc + p7P_MM),
						   *(dpp+p7GD_IL) + *(tsc + p7P_IM)),
					p7_FLogsum(*(dpp+p7GD_DL) + *(tsc + p7P_DM),
						   xL             + *(tsc + p7P_BLM)));

      mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7GD_MG) + *(tsc + p7P_MM),
						   *(dpp+p7GD_IG) + *(tsc + p7P_IM)),
                    				   *(dpp+p7GD_DG) + *(tsc + p7P_DM));
      dpp  += p7GD_NSCELLS; 

      /* I_M state doesn't exist      */
      *dpc++ = -eslINFINITY;	/* IL */
      *dpc++ = -eslINFINITY;	/* IG */

      /* E state update now includes glocal exits */
      xE  = p7_FLogsum( xE, p7_FLogsum( p7_FLogsum(mlv, dlv), p7_FLogsum(mgv, dgv)));
      
      /* D_M state: deferred storage only */
      *dpc++ = dlv;
      *dpc++ = dgv;
    
      /* row i is now finished, and dpc[] is positioned exactly on first special state, E */
      dpp += p7GD_NSCELLS;    /* now dpp[] is also positioned exactly on first special, E */
      
      *dpc++ = xE;		/* E */
      *dpc++ = xN = *(dpp + p7GD_N) + gm->xsc[p7P_N][p7P_LOOP]; /* N */
      *dpc++ = xJ = p7_FLogsum( *(dpp + p7GD_J) + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_LOOP]); /* J */
      *dpc++ = xB = p7_FLogsum(             xN  + gm->xsc[p7P_N][p7P_MOVE],  xJ + gm->xsc[p7P_J][p7P_MOVE]); /* B */
      *dpc++ = xL = xB  + gm->xsc[p7P_B][0]; /* L */
      *dpc++ = xG = xB  + gm->xsc[p7P_B][1]; /* G */
      *dpc        = p7_FLogsum( *(dpp + p7GD_C) + gm->xsc[p7P_C][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_MOVE]); /* C */
    }
  /* Done with all rows i. As we leave, dpc is still sitting on the xC value for i=L ... including even the L=0 case */
  
  if (opt_sc) *opt_sc = *dpc + gm->xsc[p7P_C][p7P_MOVE];
  gxd->M = M;
  gxd->L = L;
  return eslOK;
}


/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_DUAL_BENCHMARK
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "p7_gmxd.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,   "2000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Backward",                        0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Forward",                         0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic dual-mode Forward/Backward";

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
  P7_GMXD        *fwd     = NULL;
  P7_GMXD        *bck     = NULL;
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
  p7_profile_Config(gm, hmm, bg);
  p7_profile_SetLength(gm, L);

  fwd = p7_gmxd_Create(gm->M, L);
  bck = p7_gmxd_Create(gm->M, L);

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
      if (! esl_opt_GetBoolean(go, "-B"))  p7_GForwardDual (dsq, L, gm, fwd, &sc);
      //if (! esl_opt_GetBoolean(go, "-F"))  p7_GBackward(dsq, L, gm, bck, NULL);

      p7_gmxd_Reuse(fwd);
      //p7_gmxd_Reuse(bck);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_gmxd_Destroy(bck);
  p7_gmxd_Destroy(fwd);
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
#endif /*p7GENERIC_FWDBACK_DUAL_BENCHMARK*/
/*----------------- end, benchmark ------------------------------*/



/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_DUAL_TESTDRIVE
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* The p7_GForward() function only evaluates local alignment,
 * regardless of configuration of glocal/local in <gm>.  p7_GForward()
 * with a model configured in dual-mode should give the same score as
 * p7_GForwardDual() with a model configured in local-only mode.
 *
 * Make two profiles from a sampled <hmm>, one local-only and one dual-mode (both multihit, L=L).
 * Generate <nseq> iid random sequences of length <L>, using <bg> frequencies and the seeded <rng>.
 * Score each sequence, using p7_GForward(dual) and p7_GForwardDual(local). 
 * Check that raw nat scores match (within an absolute floating point tolerance).
 * Also, check that average bit score (e.g. expected score on random seqs) is nonpositive.
 */
static void
utest_compare_local(ESL_GETOPTS *go, ESL_RANDOMNESS *rng)
{
  char          msg[]  = "generic_fwdback_dual : compare-local unit test failed";
  ESL_DSQ      *dsq    = NULL;
  ESL_ALPHABET *abc    = NULL;
  P7_HMM       *hmm    = NULL;
  P7_PROFILE   *gmd    = NULL;
  P7_PROFILE   *gml    = NULL;
  P7_GMX       *gx     = NULL;
  P7_GMXD      *gxd    = NULL;
  P7_BG        *bg     = NULL;
  //int           M      = 100;
  //int           L      = 200;
  int           M      = 10;
  int           L      = 10;
  int           nseq   = 20;
  float         avg_sc = 0.0;
  float         sc1, sc2, nullsc;
  int           idx;
  char          errbuf[eslERRBUFSIZE];

  if ((abc = esl_alphabet_Create(eslAMINO))      == NULL)  esl_fatal(msg);
  if ( p7_hmm_Sample(rng, M, abc, &hmm)          != eslOK) esl_fatal(msg);
  if (p7_hmm_Validate (hmm, errbuf, 0.0001)      != eslOK) esl_fatal("bad hmm: %s", errbuf);

  if ((bg = p7_bg_Create(abc))                   == NULL)  esl_fatal(msg);
  if ( p7_bg_SetLength(bg, L)                    != eslOK) esl_fatal(msg);                 

  if (( gmd = p7_profile_Create(hmm->M, abc) )   == NULL)  esl_fatal(msg);
  if (( gml = p7_profile_Create(hmm->M, abc) )   == NULL)  esl_fatal(msg);

  if ( p7_profile_Config(gmd, hmm, bg)           != eslOK)  esl_fatal(msg); /* gmd is dual-mode, multihit, L=L */
  if ( p7_profile_SetLength(gmd, L)              != eslOK)  esl_fatal(msg);

  if ( p7_profile_ConfigLocal(gml, hmm, bg, L)   != eslOK)  esl_fatal(msg); /* gml is local-mode, multihit, L=L */

  if ( p7_profile_Validate(gmd,  errbuf, 0.0001) != eslOK) esl_fatal("bad profile: %s", errbuf);
  if ( p7_profile_Validate(gml,  errbuf, 0.0001) != eslOK) esl_fatal("bad profile: %s", errbuf);

  if (( dsq = malloc(sizeof(ESL_DSQ) * (L+2)))   == NULL)  esl_fatal(msg);
  if (( gx  = p7_gmx_Create(hmm->M, L))          == NULL)  esl_fatal(msg);
  if (( gxd = p7_gmxd_Create(hmm->M, L))         == NULL)  esl_fatal(msg);

  for (idx = 0; idx < nseq; idx++)
    {
      if ( esl_rsq_xfIID(rng, bg->f, abc->K, L, dsq) != eslOK) esl_fatal(msg);
      if ( p7_GForward    (dsq, L, gmd, gx,  &sc1)   != eslOK) esl_fatal(msg);
      if ( p7_GForwardDual(dsq, L, gml, gxd, &sc2)   != eslOK) esl_fatal(msg);
      if ( p7_bg_NullOne  (bg, dsq, L, &nullsc)      != eslOK) esl_fatal(msg);

      if (fabs(sc1-sc2) > 0.0001) esl_fatal(msg);

      avg_sc += (sc2 - nullsc) / eslCONST_LOG2; /* bit conversion is for consistency; unnecessary here because we're only going to check for nonpositive value */
    }
  
  avg_sc /= (float) nseq;
  if (avg_sc > 0.0) esl_fatal(msg);
  
  p7_gmxd_Destroy(gxd);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gmd);
  p7_profile_Destroy(gml);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  free(dsq);
} 


/* The "duality" test uses the fact that for unihit models, the
 * dual-mode score should be equal to FLogsum(local_score +
 * glocal_score) - 1 bit.
 */
static void
utest_duality(ESL_GETOPTS *go, ESL_RANDOMNESS *rng)
{
  char          msg[]  = "generic_fwdback_dual : duality unit test failed";
  ESL_DSQ      *dsq    = NULL;
  ESL_ALPHABET *abc    = NULL;
  P7_HMM       *hmm    = NULL;
  P7_BG        *bg     = NULL;
  P7_PROFILE   *gmd    = NULL;
  P7_PROFILE   *gml    = NULL;
  P7_PROFILE   *gmg    = NULL;
  P7_GMXD      *gxd    = NULL;
  int           M      = 100;
  int           L      = 200;
  int           nseq   = 20;
  float         dual_sc, local_sc, glocal_sc, combined_sc;
  int           idx;

  if ((abc = esl_alphabet_Create(eslAMINO))   == NULL)  esl_fatal(msg);
  if ( p7_hmm_Sample(rng, M, abc, &hmm)       != eslOK) esl_fatal(msg);
  if ((bg = p7_bg_Create(abc))                == NULL)  esl_fatal(msg);

  if (( gmd = p7_profile_Create(hmm->M, abc) )   == NULL)  esl_fatal(msg);
  if (( gml = p7_profile_Create(hmm->M, abc) )   == NULL)  esl_fatal(msg);
  if (( gmg = p7_profile_Create(hmm->M, abc) )   == NULL)  esl_fatal(msg);

  if ( p7_profile_ConfigCustom(gmd, hmm, bg, L, 0.0, 0.5)   != eslOK) esl_fatal(msg); /* unihit, dual mode        */
  if ( p7_profile_ConfigCustom(gml, hmm, bg, L, 0.0, 0.0)   != eslOK) esl_fatal(msg); /* unihit, local-only mode  */
  if ( p7_profile_ConfigCustom(gmg, hmm, bg, L, 0.0, 1.0)   != eslOK) esl_fatal(msg); /* unihit, glocal-only mode */

  if (( dsq = malloc(sizeof(ESL_DSQ) * (L+2))) == NULL)  esl_fatal(msg);
  if (( gxd = p7_gmxd_Create(hmm->M, L))       == NULL)  esl_fatal(msg);

  for (idx = 0; idx < nseq; idx++)
    {
      if ( esl_rsq_xfIID(rng, bg->f, abc->K, L, dsq) != eslOK) esl_fatal(msg);

      if ( p7_GForwardDual(dsq, L, gmd, gxd, &dual_sc)   != eslOK) esl_fatal(msg);
      if ( p7_gmxd_Reuse(gxd)                            != eslOK) esl_fatal(msg);

      if ( p7_GForwardDual(dsq, L, gml, gxd, &local_sc)  != eslOK) esl_fatal(msg);
      if ( p7_gmxd_Reuse(gxd)                            != eslOK) esl_fatal(msg);

      if ( p7_GForwardDual(dsq, L, gmg, gxd, &glocal_sc) != eslOK) esl_fatal(msg);
      if ( p7_gmxd_Reuse(gxd)                            != eslOK) esl_fatal(msg);

      combined_sc = p7_FLogsum(local_sc, glocal_sc) - eslCONST_LOG2;

      if (fabs(dual_sc-combined_sc) > 0.001)  esl_fatal(msg);
    }
  
  p7_gmxd_Destroy(gxd);
  p7_profile_Destroy(gmg);
  p7_profile_Destroy(gml);
  p7_profile_Destroy(gmd);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  free(dsq);
}


/* The "enumeration" test samples a random enumerable HMM. This HMM
 * has tII transitions all zero, so the generated sequence space
 * ranges from L=0..2M+1 for the HMM. It uses this to create an
 * enumerable profile, by using a unihit L=0 configuration; this
 * profile generates all sequences of lengths L=1..2M-1. (The
 * differences from the HMM are 1) the I0 and Im states are normalized
 * away, and 2) the B->DDD->E mute path that generates zero residues
 * is also normalized away.)
 * 
 * The profile is configured in dual local/glocal mode, so that the
 * test will exercise all paths (except II transitions) in dual-mode
 * DP calculations.
 * 
 * Then the test enumerates all those sequences, scores them with
 * p7_GForwardDual(), obtains P(seq | profile) from the score, and
 * sums P(seq | profile) over the enumerable space of profiles.  This
 * sum should be 1.0, within floating point tolerance.
 * 
 * All P(seq | profile) terms need to be >> DBL_EPSILON for the
 * summation to work. This means M must be very small -- perhaps on
 * the order of ~5. Small M also helps the enumeration run quickly.
 * Even a short M suffices to detect most conceivable failure modes in
 * a DP implementation.
 * 
 * To speed up the enumeration we use a tiny alphabet, <eslCOINS>.
 * Incidentally, this also helps us test this rarely-used Easel
 * alphabet, and whether HMMER can deal with non-bio alphabets.
 *
 * The enumeration test in generic_fwdback.c is similar, but uses
 * a different enumeration: p7_hmm_SampleEnumerable() instead of
 * p7_hmm_SampleEnumerable2(). p7_hmm_SampleEnumerable() sets all
 * transitions to insert to 0, so it enumerates a smaller seq space of
 * L=0..M (no inserts are possible at all.)
 */
static void
utest_enumeration(ESL_GETOPTS *go, ESL_RANDOMNESS *rng)
{
  char          msg[] = "enumeration unit test failed";
  ESL_ALPHABET *abc   = esl_alphabet_Create(eslCOINS);
  P7_HMM       *hmm   = NULL;
  P7_BG        *bg    = NULL;
  ESL_DSQ      *dsq   = NULL;
  P7_PROFILE   *gm    = NULL;
  P7_GMXD      *gxd   = NULL;
  int           M     = 8;
  int           maxL  = 2*M-1;	
  int           i, L;
  float         fsc;
  float         bg_ll;
  double        fp;
  double        total_p = 0.0;

  if ( p7_hmm_SampleEnumerable2(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if (( bg = p7_bg_Create(abc))                    == NULL)  esl_fatal(msg);
  if (( gm = p7_profile_Create(hmm->M, abc))       == NULL)  esl_fatal(msg);
  
                                         /* L, nj,  pglocal:  L=0 unihit dual-mode */
  if ( p7_profile_ConfigCustom(gm, hmm, bg, 0, 0.0, 0.5) != eslOK) esl_fatal(msg);

  if (( dsq = malloc(sizeof(ESL_DSQ) * (maxL+3))) == NULL)  esl_fatal(msg); /* 1..2*M-1, +2 for sentinels at 0, 2*M, +1 for the maxL+1 test seq */
  if (( gxd = p7_gmxd_Create(hmm->M, maxL+1))     == NULL)  esl_fatal(msg); /* +1 because of the maxL+1 final test */

  /* L=0 included just to test that an L=0 sequence does indeed get a score of -inf, as it should */
  for (L = 0; L <= maxL; L++)
    {
      /* initialize dsq[1..L] at "0000..." */
      dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;
      for (i = 1; i <= L; i++) dsq[i] = 0;

      /* enumerate and score all sequences of length L */
      while (1)	
	{
	  if ( p7_GForwardDual(dsq, L, gm, gxd, &fsc) != eslOK) esl_fatal(msg);
	  
	  /* calculate bg log likelihood component of the scores */
	  for (bg_ll = 0., i = 1; i <= L; i++)  bg_ll += log(bg->f[dsq[i]]);
	  
	  /* convert to probability P(seq|model), adding the bg LL back to the LLR */
	  fp =  exp(fsc + bg_ll);
	  total_p += fp;

	  /* Increment dsq to next seq, like a reversed odometer; works for any alphabet */
	  for (i = 1; i <= L; i++) 
	    if (dsq[i] < abc->K-1) { dsq[i]++; break; } else { dsq[i] = 0; }
	  if (i > L) break;	/* we're done enumerating sequences */

	  p7_gmxd_Reuse(gxd);
	}
    }

  /* That sum is subject to significant numerical error because of
   * discretization error in FLogsum(); don't expect it to be too close.
   */
  if (total_p < 0.999 || total_p > 1.001) esl_fatal(msg);

  /* And any sequence of length L > maxL should get score -infinity. */
  if ( esl_rsq_xfIID(rng, bg->f, abc->K, maxL+1, dsq) != eslOK) esl_fatal(msg);
  if ( p7_GForwardDual(dsq, maxL+1, gm, gxd, &fsc)    != eslOK) esl_fatal(msg);
  if ( fsc != -eslINFINITY) esl_fatal(msg);                                    

  p7_gmxd_Destroy(gxd);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  free(dsq);
}
#endif /*p7GENERIC_FWDBACK_DUAL_TESTDRIVE*/

/*----------------- end, unit tests -----------------------------*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_DUAL_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"

#include "hmmer.h"
#include "p7_gmxd.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be very verbose",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for the generic Forward/Backward dual-mode implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));

  p7_FLogsumInit();

  utest_compare_local(go, r);
  utest_duality      (go, r);
  utest_enumeration  (go, r);

  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}


#endif /*p7GENERIC_FWDBACK_DUAL_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/



/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_DUAL_EXAMPLE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_gmxd.h"

#define STYLES     "--fs,--sw,--ls,--s"	               /* Exclusive choice for alignment mode     */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "--fs",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "multihit local alignment",                         0 },
  { "--sw",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "unihit local alignment",                           0 },
  { "--ls",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "multihit glocal alignment",                        0 },
  { "--s",       eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "unihit glocal alignment",                          0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "very verbose debugging output: inc. DP matrix",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of Forward/Backward, generic dual local/glocal implementation";

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
  P7_GMXD        *fwd     = NULL;
  P7_GMXD        *bck     = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fsc, bsc;
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

  /* Now reconfig the models however we were asked to */
  if      (esl_opt_GetBoolean(go, "--fs"))  p7_profile_ConfigLocal    (gm, hmm, bg, sq->n);
  else if (esl_opt_GetBoolean(go, "--sw"))  p7_profile_ConfigUnilocal (gm, hmm, bg, sq->n);
  else if (esl_opt_GetBoolean(go, "--ls"))  p7_profile_ConfigGlocal   (gm, hmm, bg, sq->n);
  else if (esl_opt_GetBoolean(go, "--s"))   p7_profile_ConfigUniglocal(gm, hmm, bg, sq->n);
  else                                      p7_profile_Config         (gm, hmm, bg);

  /* Allocate matrices */
  fwd = p7_gmxd_Create(gm->M, sq->n);

  printf("%-30s   %-10s %-10s   %-10s %-10s\n", "# seq name",      "fwd (raw)",   "bck (raw) ",  "fwd (bits)",  "bck (bits)");
  printf("%-30s   %10s %10s   %10s %10s\n",     "#--------------", "----------",  "----------",  "----------",  "----------");

  while ( (status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
      else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

      /* Resize the DP matrices if necessary */
      p7_gmxd_GrowTo(fwd, gm->M, sq->n);

      /* Set the profile and null model's target length models */
      p7_bg_SetLength     (bg, sq->n);
      p7_profile_SetLength(gm, sq->n);

      //p7_profile_Dump(stdout, gm);

      /* Run Forward, Backward */
      p7_GForwardDual (sq->dsq, sq->n, gm, fwd, &fsc);
      bsc = 0.0;

      p7_gmxd_Dump(stdout, fwd);

      if (esl_opt_GetBoolean(go, "--vv")) p7_gmxd_Dump(stdout, fwd);

      /* Those scores are partial log-odds likelihoods in nats.
       * Subtract off the rest of the null model, convert to bits.
       */
      p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

      printf("%-30s   %10.4f %10.4f   %10.4f %10.4f\n", 
	     sq->name, 
	     fsc, bsc, 
	     (fsc - nullsc) / eslCONST_LOG2, (bsc - nullsc) / eslCONST_LOG2);

      p7_gmxd_Reuse(fwd);
      //p7_gmxd_Reuse(bck);
      esl_sq_Reuse(sq);
    }

  /* Cleanup */
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_gmxd_Destroy(fwd);
  p7_gmxd_Destroy(bck);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_FWDBACK_DUAL_EXAMPLE*/
/*-------------------- end, example -----------------------------*/


/*****************************************************************
 * @LICENSE@
 *
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
