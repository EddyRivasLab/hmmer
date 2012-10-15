/* Reference implementation of Viterbi scoring and alignment.
 *   dual-mode (local/glocal)
 *   quadratic memory (not banded, not checkpointed)
 *   standard C code (not striped, not vectorized)
 *   
 * The reference implementation is for testing and debugging, and to
 * provide an example simpler than our production DP code, which
 * layers on some more complicated techniques (banding, vectorization,
 * checkpointing).
 * 
 * Contents:
 *   1. Viterbi DP fill.
 *   2. Viterbi optimal traceback.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 *   6. Copyright and license information.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_trace.h"
#include "dp_reference/reference_viterbi.h"


/*****************************************************************
 * 1. Viterbi DP fill.
 *****************************************************************/

/* Function:  p7_ReferenceViterbi()
 * Synopsis:  Reference implementation of the Viterbi algorithm.
 *
 * Purpose:   Given a target sequence <dsq> of length <L> residues, a
 *            query profile <gm>, and a DP matrix <rmx>;
 *            do the Viterbi optimal alignment algorithm.
 *            Return the Viterbi score in nats in <*opt_sc>, if
 *            caller provides it. 
 *            Return the Viterbi optimal trace in <*opt_tr>, if
 *            caller provides an allocated trace structure.
 *            
 *            <rmx> will be reallocated if needed, so it can be reused
 *            from a previous calculation, even a smaller one.
 *            
 * Args:      dsq     - digital target sequence 1..L
 *            L       - length of <dsq> in residues
 *            gm      - query profile
 *            rmx     - DP matrix
 *            opt_tr  - optRETURN: trace structure for Viterbi alignment, or NULL
 *            opt_sc  - optRETURN: Viterbi raw score in nats, or NULL
 *
 * Returns:   <eslOK> on success. <rmx> contains the Viterbi matrix.
 * 
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_ReferenceViterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_REFMX *rmx, P7_TRACE *opt_tr, float *opt_sc)
{
  int    M          = gm->M;
  float *dpc, *dpp;
  const float *tsc;		/* ptr for stepping thru profile's transition parameters */
  const float *rsc;		/* ptr for stepping thru profile's emission parameters   */
  int    i, k, s;
  float  mlv, mgv;	      /* ML,MG cell values on current row   */
  float  dlv, dgv; 	      /* pushed-ahead DL,DG cell k+1 values */
  float  xE, xL, xG;
  float  vsc;
  int    status;

  /* contract checks / arg validation */
  ESL_DASSERT1( ( gm->L == L || gm->L == 0) ); /* length model in profile is either L (usually) or 0 (some unit tests) */

  /* reallocation, if needed */
  if ( (status = p7_refmx_GrowTo(rmx, gm->M, L)) != eslOK) return status;
  rmx->M    = M;
  rmx->L    = L;
  rmx->type = p7R_VITERBI;


  /* Initialization of the zero row. */
  dpc = rmx->dp[0];
  for (s = 0; s < (M+1) * p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; // all M,I,D; k=0..M
  dpc[p7R_E]      = -eslINFINITY;	                               
  dpc[p7R_N]      = 0.0f;			                       
  dpc[p7R_J]      = -eslINFINITY;                                   
  dpc[p7R_B]      = gm->xsc[p7P_N][p7P_MOVE];                       
  dpc[p7R_L] = xL = gm->xsc[p7P_N][p7P_MOVE] + gm->xsc[p7P_B][0];
  dpc[p7R_G] = xG = gm->xsc[p7P_N][p7P_MOVE] + gm->xsc[p7P_B][1];
  dpc[p7R_C]      = -eslINFINITY;                              
  dpc[p7R_JJ]     = -eslINFINITY;                              
  dpc[p7R_CC]     = -eslINFINITY;                              
  /* *dpc must be on specials of 0 row, to handle L=0 case where the recursion body below doesn't run */

  /* Main DP recursion */
  for (i = 1; i <= L; i++)
    {
      /* Initialization for a new row */
      rsc = gm->rsc[dsq[i]] + p7P_NR;	/* this ptr steps through the row's emission scores 1..M. skip k=0 */
      tsc = gm->tsc;			/* this ptr steps through profile's transition scores 0..M         */
      dpp = rmx->dp[i-1];               /* previous row dpp is already set, and at k=0 */
      dpc = rmx->dp[i];                 /* current DP row, skip k=0, start at k=1.  */
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;

      dlv = dgv = -eslINFINITY;
      xE  =       -eslINFINITY;

      /* Main inner loop of the recursion */
      for (k = 1; k < M; k++)
	{
	  /* match states MLk, MGk */
	  mlv = dpc[p7R_ML] = *rsc + ESL_MAX( ESL_MAX(dpp[p7R_ML] + tsc[p7P_MM],
						      dpp[p7R_IL] + tsc[p7P_IM]),
					      ESL_MAX(dpp[p7R_DL] + tsc[p7P_DM],
						      xL          + tsc[p7P_LM]));
	  mgv = dpc[p7R_MG] = *rsc + ESL_MAX( ESL_MAX(dpp[p7R_MG] + tsc[p7P_MM],
					              dpp[p7R_IG] + tsc[p7P_IM]),
 					      ESL_MAX(dpp[p7R_DG] + tsc[p7P_DM],
					   	      xG          + tsc[p7P_GM]));

	  rsc++;                /* rsc advances to insert score for position k */
	  tsc += p7P_NTRANS;    /* tsc advances to transitions in states k     */
	  dpp += p7R_NSCELLS;	/* dpp advances to cells for states k          */

	  /* Insert state calculations ILk, IGk. */
	  dpc[p7R_IL] = *rsc + ESL_MAX( dpp[p7R_ML] + tsc[p7P_MI], dpp[p7R_IL] + tsc[p7P_II]);
	  dpc[p7R_IG] = *rsc + ESL_MAX( dpp[p7R_MG] + tsc[p7P_MI], dpp[p7R_IG] + tsc[p7P_II]);
	  rsc++;		/* rsc advances to next match state emission   */

	  /* E state update; local paths only, transition prob 1.0 in implicit probability model, Dk->E can never be optimal */
	  xE  = ESL_MAX( mlv, xE);

	  /* Delete state, deferred storage trick */
	  dpc[p7R_DL] = dlv;
	  dpc[p7R_DG] = dgv;
	  dlv = ESL_MAX( mlv + tsc[p7P_MD], dlv + tsc[p7P_DD]);
	  dgv = ESL_MAX( mgv + tsc[p7P_MD], dgv + tsc[p7P_DD]);

	  dpc += p7R_NSCELLS;
	}

      /* k=M node is unrolled and handled separately. No I state, and glocal exits. */
      mlv = dpc[p7R_ML] = *rsc + ESL_MAX( ESL_MAX(dpp[p7R_ML] + tsc[p7P_MM],
						  dpp[p7R_IL] + tsc[p7P_IM]),
					  ESL_MAX(dpp[p7R_DL] + tsc[p7P_DM],
						  xL          + tsc[p7P_LM]));
      mgv = dpc[p7R_MG] = *rsc + ESL_MAX( ESL_MAX(dpp[p7R_MG] + tsc[p7P_MM],
						  dpp[p7R_IG] + tsc[p7P_IM]),
					  ESL_MAX(dpp[p7R_DG] + tsc[p7P_DM],
						  xG          + tsc[p7P_GM]));
      dpp  += p7R_NSCELLS; 

      /* I_M state doesn't exist */
      dpc[p7R_IL] = -eslINFINITY;
      dpc[p7R_IG] = -eslINFINITY;

      /* E state update now includes glocal exits: transition prob 1.0 from MG_m; DLk->E still can't be optimal, but DGk->E can */
      xE  = ESL_MAX( ESL_MAX( mgv, dgv),
		     ESL_MAX( xE,  mlv));
      
      /* D_M state: deferred storage only */
      dpc[p7R_DL] = dlv;
      dpc[p7R_DG] = dgv;
    
      /* row i is now finished; position both dpc, dpp on specials */
      dpc += p7R_NSCELLS;
      dpp += p7R_NSCELLS;   
      
      dpc[p7R_E]      = xE;	
      dpc[p7R_N]      = dpp[p7R_N] + gm->xsc[p7P_N][p7P_LOOP];
      dpc[p7R_J]      = ESL_MAX( dpp[p7R_J] + gm->xsc[p7P_J][p7P_LOOP],  dpc[p7R_E] + gm->xsc[p7P_E][p7P_LOOP]); 
      dpc[p7R_B]      = ESL_MAX( dpc[p7R_N] + gm->xsc[p7P_N][p7P_MOVE],  dpc[p7R_J] + gm->xsc[p7P_J][p7P_MOVE]); 
      dpc[p7R_L] = xL = dpc[p7R_B] + gm->xsc[p7P_B][0];
      dpc[p7R_G] = xG = dpc[p7R_B] + gm->xsc[p7P_B][1];
      dpc[p7R_C]      = ESL_MAX( dpp[p7R_C] + gm->xsc[p7P_C][p7P_LOOP],  dpc[p7R_E] + gm->xsc[p7P_E][p7P_MOVE]); 
      dpc[p7R_JJ]     = -eslINFINITY;
      dpc[p7R_CC]     = -eslINFINITY;
    }
  /* Done with all rows i. As we leave, dpc is still sitting on the specials for i=L ... including even the L=0 case */
  
  vsc =  dpc[p7R_C] + gm->xsc[p7P_C][p7P_MOVE];
  if (opt_sc) *opt_sc = vsc;
  if (opt_tr && vsc != -eslINFINITY) return p7_reference_trace_Viterbi(gm, rmx, opt_tr);  // if no paths are possible at all: leave tr->N=0, our convention for impossible trace
  else                               return eslOK;
}
/*------------------ end, viterbi DP fill -----------------------*/


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7REFERENCE_VITERBI_TESTDRIVE

#include "esl_random.h"
#include "esl_randomseq.h"

/* "randomseq" test:
 * Does Viterbi on query <gm> against random sequences. 
 * Tests:
 *   1. Viterbi score equals score of Viterbi trace.
 *   2. Viterbi traces pass trace validation
 */
static void
utest_randomseq(ESL_RANDOMNESS *rng, P7_PROFILE *gm, P7_BG *bg, int nseq, int L)
{
  char      msg[]   = "viterbi randomseq test failed";
  ESL_DSQ   *dsq    = NULL;
  P7_REFMX  *rmx    = NULL;
  P7_TRACE  *tr     = NULL;
  int        idx;
  float      sc1, sc2;
  char       errbuf[eslERRBUFSIZE];

  if (( dsq = malloc(sizeof(ESL_DSQ) * (L+2))) == NULL) esl_fatal(msg);
  if (( tr  = p7_trace_Create())               == NULL) esl_fatal(msg);
  if (( rmx = p7_refmx_Create(gm->M, L))       == NULL) esl_fatal(msg);
  
  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rsq_xfIID(rng, bg->f, gm->abc->K, L, dsq)  != eslOK) esl_fatal(msg);
      if (p7_ReferenceViterbi(dsq, L, gm, rmx, tr, &sc1) != eslOK) esl_fatal(msg);
      if (p7_trace_Validate(tr, gm->abc, dsq, errbuf)    != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
      if (p7_trace_Score(tr, dsq, gm, &sc2)              != eslOK) esl_fatal(msg);
      if (esl_FCompareAbs(sc1, sc2, 1e-6)                != eslOK) esl_fatal(msg);
      if (! isfinite(sc1))                                         esl_fatal(msg);
      if (! isfinite(sc2))                                         esl_fatal(msg);

      p7_trace_Reuse(tr);
      p7_refmx_Reuse(rmx);
    }
  p7_refmx_Destroy(rmx);
  p7_trace_Destroy(tr);
  free(dsq);
}

/* "generation" test
 * Sample sequences emitted from the profile ("true positives").
 * Test:
 *    1. Score of emission trace <= score of Viterbi trace
 *    2. Viterbi score == score of Viterbi trace
 *    3. Viterbi traces pass trace validation
 */
static void
utest_generation(ESL_RANDOMNESS *rng, P7_HMM *hmm, P7_PROFILE *gm, P7_BG *bg, int nseq, int L)
{
  char      msg[]   = "viterbi generation test failed";
  ESL_SQ    *sq     = esl_sq_CreateDigital(gm->abc);
  P7_REFMX  *rmx    = p7_refmx_Create(gm->M, 400);
  P7_TRACE  *tr     = p7_trace_Create();
  int        idx;
  float      sc1, sc2, sc3;
  char       errbuf[eslERRBUFSIZE];

  for (idx = 0; idx < nseq; idx++)
    {
      p7_profile_SetLength(gm, L);
      do {
	esl_sq_Reuse(sq);
	if (p7_ProfileEmit(rng, hmm, gm, bg, sq, tr)             != eslOK) esl_fatal(msg);
      } while (sq->n > (gm->M+gm->L) * 3); /* keep sampled sequence length down to something arbitrarily reasonable */
      if (p7_trace_Validate(tr, gm->abc, sq->dsq, errbuf)        != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
      if (p7_trace_Score   (tr, sq->dsq, gm, &sc1)               != eslOK) esl_fatal(msg);
      p7_trace_Reuse(tr);

      p7_profile_SetLength(gm, sq->n);
      if (p7_ReferenceViterbi(sq->dsq, sq->n, gm, rmx, tr, &sc2) != eslOK) esl_fatal(msg);
      if (p7_trace_Validate(tr, gm->abc, sq->dsq, errbuf)        != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
      if (p7_trace_Score(tr, sq->dsq, gm, &sc3)                  != eslOK) esl_fatal(msg);

      if (! isfinite(sc1))                            esl_fatal(msg);
      if (! isfinite(sc2))                            esl_fatal(msg);
      if (! isfinite(sc3))                            esl_fatal(msg);
      if (sc1 > sc2)                                  esl_fatal(msg); /* score of generated trace should be <= Viterbi score */
      if (esl_FCompareAbs(sc2, sc3, 0.0001) != eslOK) esl_fatal(msg);

      esl_sq_Reuse(sq);
      p7_trace_Reuse(tr);
      p7_refmx_Reuse(rmx);
    }
  p7_refmx_Destroy(rmx);
  p7_trace_Destroy(tr);
  esl_sq_Destroy(sq);
}
#endif /*p7REFERENCE_VITERBI_TESTDRIVE*/
/*---------------- end, unit tests ------------------------------*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7REFERENCE_VITERBI_TESTDRIVE

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 }, 
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for reference Viterbi implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go           = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  int             M            = esl_opt_GetInteger(go, "-M");
  int             L            = esl_opt_GetInteger(go, "-L");
  int             N            = esl_opt_GetInteger(go, "-N");
  ESL_RANDOMNESS *rng          = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc          = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg           = p7_bg_Create(abc);
  P7_HMM         *hmm          = NULL;
  P7_PROFILE     *gm           = p7_profile_Create(M, abc);
  
  p7_hmm_Sample(rng, M, abc, &hmm);
  p7_profile_Config   (gm, hmm, bg);   
  p7_profile_SetLength(gm, L);

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_randomseq (rng,      gm, bg, N, L);
  utest_generation(rng, hmm, gm, bg, N, L);

  fprintf(stderr, "#  status = ok\n");
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
}
#endif /*p7REFERENCE_VITERBI_TESTDRIVE*/
/*-------------- end, test driver -------------------------------*/



/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7REFERENCE_VITERBI_EXAMPLE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "p7_refmx.h"

#define STYLES     "--fs,--sw,--ls,--s"	               /* Exclusive choice for alignment mode     */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "--fs",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "multihit local alignment",                         0 },
  { "--sw",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "unihit local alignment",                           0 },
  { "--ls",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "multihit glocal alignment",                        0 },
  { "--s",       eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "unihit glocal alignment",                          0 },
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump Viterbi matrix for inspection",               0 },
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "dump optimal Viterbi trace for inspection",        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of using the Viterbi reference implementation";

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
  P7_REFMX       *vit     = p7_refmx_Create(200, 400); /* will grow as needed */
  P7_TRACE       *tr      = p7_trace_Create();
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           vsc, nullsc, dvsc;
  int             d;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  if      (esl_opt_GetBoolean(go, "--fs"))  p7_profile_ConfigLocal    (gm, hmm, bg, 400);
  else if (esl_opt_GetBoolean(go, "--sw"))  p7_profile_ConfigUnilocal (gm, hmm, bg, 400);
  else if (esl_opt_GetBoolean(go, "--ls"))  p7_profile_ConfigGlocal   (gm, hmm, bg, 400);
  else if (esl_opt_GetBoolean(go, "--s"))   p7_profile_ConfigUniglocal(gm, hmm, bg, 400);
  else                                      p7_profile_Config         (gm, hmm, bg);

  while ( (status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
      else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

      /* Set the profile and null model's target length models */
      p7_bg_SetLength     (bg, sq->n);
      p7_profile_SetLength(gm, sq->n);

      /* Run Viterbi - get raw score and optimal trace */
      p7_ReferenceViterbi(sq->dsq, sq->n, gm, vit, tr, &vsc);    
      if (esl_opt_GetBoolean(go, "-D")) p7_refmx_Dump(stdout, vit);
      if (esl_opt_GetBoolean(go, "-T")) p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);
      
      /* Those scores are partial log-odds likelihoods in nats.
       * Subtract off the rest of the null model, convert to bits.
       */
      p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

      p7_trace_Index(tr);
      for (d = 0; d < tr->ndom; d++)
	{
	  p7_trace_ScoreDomain(tr, sq->dsq, gm, d, &dvsc);

	  printf("%-20s %-30s %10.2f %2d/%-2d %5d %5d %10.2f\n",
		 gm->name, sq->name,
		 (vsc - nullsc) / eslCONST_LOG2,
		 d+1, tr->ndom,
		 tr->sqfrom[d], tr->sqto[d],
		 (dvsc - nullsc) / eslCONST_LOG2);
	}

#if 0
      printf("%-30s   %10.4f nats (raw)    %10.4f bits\n", 
	     sq->name,
	     vsc, 
	     (vsc - nullsc) / eslCONST_LOG2);
#endif

      p7_refmx_Reuse(vit);
      p7_trace_Reuse(tr);
      esl_sq_Reuse(sq);
    }

  /* Cleanup */
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(vit);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_VITERBI_EXAMPLE*/



/*****************************************************************
 * @LICENSE@
 * 
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/


