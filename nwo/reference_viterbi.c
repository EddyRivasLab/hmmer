/* Reference implementation of Viterbi scoring and alignment.
 *    - quadratic memory (not banded, not checkpointed) 
 *    - standard C code (not striped, not vectorized)
 *   
 * The reference implementation is for testing and debugging. It also
 * provides an example simpler than our production DP code, which
 * layers on some more complicated techniques (banding, vectorization,
 * checkpointing).
 * 
 * Contents:
 *   1. Viterbi DP fill.
 *   2. Stats driver: numerical error characterization.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_sq.h"
#include "esl_vectorops.h"

#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_path.h"
#include "h4_refmx.h"

#include "reference_trace.h"
#include "reference_viterbi.h"




/*****************************************************************
 * 1. Viterbi DP fill.
 *****************************************************************/

/* Function:  h4_ReferenceViterbi()
 * Synopsis:  Reference implementation of the Viterbi algorithm.
 *
 * Purpose:   Given a target sequence <dsq> of length <L> residues, a
 *            query profile <hmm> in mode <mo>, and a DP matrix <rmx>; do the
 *            Viterbi optimal alignment algorithm.  Return the Viterbi
 *            score in nats in <*opt_sc>, if caller provides it.
 *            Return the Viterbi optimal trace in <*opt_pi>, if caller
 *            provides an allocated trace structure.
 *            
 *            <rmx> will be reallocated if needed, so it can be reused
 *            from a previous calculation, even a smaller one.
 *            
 * Args:      dsq     - digital target sequence 1..L
 *            L       - length of <dsq> in residues
 *            hmm     - query profile
 *            mo      - profile's mode (algorithm-dependent parameters)
 *            rmx     - DP matrix
 *            opt_pi  - optRETURN: path structure for Viterbi alignment, or NULL
 *            opt_sc  - optRETURN: Viterbi raw score in nats, or NULL
 *
 * Returns:   <eslOK> on success. <rmx> contains the Viterbi matrix.
 * 
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
h4_ReferenceViterbi(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx, H4_PATH *opt_pi, float *opt_sc)
{
  int    M = hmm->M;
  float *dpc, *dpp;
  const float *tsc;	      // ptr for stepping thru profile's transition parameters 
  const float *rsc;	      // ptr for stepping thru profile's emission parameters   
  int    i, k, s;
  float  dlv, dgv; 	      // pushed-ahead DL,DG cell k+1 values 
  float  xE, xL, xG;
  float  vsc;
  int    status;

  /* contract checks / arg validation */
  ESL_DASSERT1( ( mo->L == L || mo->L == 0) ); /* length model is either L (usually) or 0 (some unit tests) */

  /* reallocation, if needed */
  if ( (status = h4_refmx_GrowTo(rmx, hmm->M, L)) != eslOK) return status;
  rmx->M    = M;
  rmx->L    = L;
  rmx->type = h4R_VITERBI;


  /* Initialization of the zero row. */
  dpc = rmx->dp[0];
  for (s = 0; s < (M+1) * h4R_NSCELLS; s++) *dpc++ = -eslINFINITY; // all M,I,D; k=0..M
  dpc[h4R_E]      = -eslINFINITY;	                               
  dpc[h4R_N]      = 0.0f;			                       
  dpc[h4R_J]      = -eslINFINITY;                                   
  dpc[h4R_B]      = mo->xsc[h4_N][h4_MOVE];                       
  dpc[h4R_L] = xL = mo->xsc[h4_N][h4_MOVE] + mo->xsc[h4_B][0];
  dpc[h4R_G] = xG = mo->xsc[h4_N][h4_MOVE] + mo->xsc[h4_B][1];
  dpc[h4R_C]      = -eslINFINITY;                              
  dpc[h4R_JJ]     = -eslINFINITY;                              
  dpc[h4R_CC]     = -eslINFINITY;                              
  /* *dpc is (and must) be on specials of 0 row, to handle L=0 case where the recursion body below doesn't run */

  /* Main DP recursion */
  for (i = 1; i <= L; i++)
    {
      /* Initialization for a new row */
      rsc = hmm->rsc[dsq[i]] + 1;	/* this ptr steps through the row's emission scores 1..M. skip k=0 */
      tsc = hmm->tsc[0];		/* this ptr steps through profile's transition scores 0..M. tsc[0] can be treated as one long vector. */
      dpp = rmx->dp[i-1];               /* previous row dpp is already calculated; set to k=0 */
      dpc = rmx->dp[i];                 /* current DP row, skip k=0, start at k=1.  */
      for (s = 0; s < h4R_NSCELLS; s++) *dpc++ = -eslINFINITY;

      dlv = dgv = -eslINFINITY;
      xE  =       -eslINFINITY;

      /* Main inner loop of the recursion */
      for (k = 1; k < M; k++)
	{
	  /* match states MLk, MGk */
	  dpc[h4R_ML] = *rsc + ESL_MAX( ESL_MAX(dpp[h4R_ML] + tsc[h4_MM],
	 				        dpp[h4R_IL] + tsc[h4_IM]),
					ESL_MAX(dpp[h4R_DL] + tsc[h4_DM],
						xL          + tsc[h4_LM]));
	  dpc[h4R_MG] = *rsc + ESL_MAX( ESL_MAX(dpp[h4R_MG] + tsc[h4_MM],
						dpp[h4R_IG] + tsc[h4_IM]),
					ESL_MAX(dpp[h4R_DG] + tsc[h4_DM],
						xG          + tsc[h4_GM]));

	  rsc++;                /* rsc advances to next score (for position k+1) */
	  tsc += h4_NTSC;       /* tsc advances to transitions in states k       */
	  dpp += h4R_NSCELLS;	/* dpp advances to cells for states k            */

	  /* Insert state calculations ILk, IGk. */
	  dpc[h4R_IL] = ESL_MAX( ESL_MAX( dpp[h4R_ML] + tsc[h4_MI],
					  dpp[h4R_IL] + tsc[h4_II]),
				          dpp[h4R_DL] + tsc[h4_DI]);
	  dpc[h4R_IG] = ESL_MAX( ESL_MAX( dpp[h4R_MG] + tsc[h4_MI],
					  dpp[h4R_IG] + tsc[h4_II]),
				          dpp[h4R_DG] + tsc[h4_DI]);

	  /* E state update; local paths only, transition prob 1.0 in implicit probability model, Dk->E can never be optimal */
	  xE  = ESL_MAX( dpc[h4R_ML], xE);

	  /* Delete state, deferred storage trick */
	  dpc[h4R_DL] = dlv;
	  dpc[h4R_DG] = dgv;
	  dlv = ESL_MAX( ESL_MAX (dpc[h4R_ML] + tsc[h4_MD],
				  dpc[h4R_IL] + tsc[h4_ID]),
				  dpc[h4R_DL] + tsc[h4_DD]);
	  dgv = ESL_MAX( ESL_MAX (dpc[h4R_MG] + tsc[h4_MD],
				  dpc[h4R_IG] + tsc[h4_ID]),
				  dpc[h4R_DG] + tsc[h4_DD]);
	  dpc += h4R_NSCELLS;
	}

      /* k=M node is unrolled and handled separately. No I state, and glocal exits. */
      dpc[h4R_ML] = *rsc + ESL_MAX( ESL_MAX(dpp[h4R_ML] + tsc[h4_MM],
					    dpp[h4R_IL] + tsc[h4_IM]),
				    ESL_MAX(dpp[h4R_DL] + tsc[h4_DM],
					    xL          + tsc[h4_LM]));
      dpc[h4R_MG] = *rsc + ESL_MAX( ESL_MAX(dpp[h4R_MG] + tsc[h4_MM],
					    dpp[h4R_IG] + tsc[h4_IM]),
				    ESL_MAX(dpp[h4R_DG] + tsc[h4_DM],
					    xG          + tsc[h4_GM]));
      dpp  += h4R_NSCELLS; 

      /* I_M state doesn't exist */
      dpc[h4R_IL] = -eslINFINITY;
      dpc[h4R_IG] = -eslINFINITY;

      /* E state update now includes glocal exits: transition prob 1.0 from MG_m; DLk->E still can't be optimal, but DGk->E can */
      xE  = ESL_MAX( ESL_MAX( dpc[h4R_MG], dgv),
		     ESL_MAX( dpc[h4R_ML], xE));
      
      /* D_M state: deferred storage only */
      dpc[h4R_DL] = dlv;
      dpc[h4R_DG] = dgv;
    
      /* row i is now finished; position both dpc, dpp on specials */
      dpc += h4R_NSCELLS;
      dpp += h4R_NSCELLS;   
      
      dpc[h4R_E]      = xE;	
      dpc[h4R_N]      = dpp[h4R_N] + mo->xsc[h4_N][h4_LOOP];
      dpc[h4R_J]      = ESL_MAX( dpp[h4R_J] + mo->xsc[h4_J][h4_LOOP],  dpc[h4R_E] + mo->xsc[h4_E][h4_LOOP]); 
      dpc[h4R_B]      = ESL_MAX( dpc[h4R_N] + mo->xsc[h4_N][h4_MOVE],  dpc[h4R_J] + mo->xsc[h4_J][h4_MOVE]); 
      dpc[h4R_L] = xL = dpc[h4R_B] + mo->xsc[h4_B][0];
      dpc[h4R_G] = xG = dpc[h4R_B] + mo->xsc[h4_B][1];
      dpc[h4R_C]      = ESL_MAX( dpp[h4R_C] + mo->xsc[h4_C][h4_LOOP],  dpc[h4R_E] + mo->xsc[h4_E][h4_MOVE]); 
      dpc[h4R_JJ]     = -eslINFINITY;
      dpc[h4R_CC]     = -eslINFINITY;
    }
  /* Done with all rows i. As we leave, dpc is still sitting on the specials for i=L ... including even the L=0 case */
  
  vsc =  dpc[h4R_C] + mo->xsc[h4_C][h4_MOVE];
  if (opt_sc) *opt_sc = vsc;
  if (opt_pi && vsc != -eslINFINITY) return h4_reference_trace_Viterbi(hmm, mo, rmx, opt_pi);  // if no paths are possible at all: leave pi->N=0, our convention for impossible trace
  else                               return eslOK;
}
/*------------------ end, viterbi DP fill -----------------------*/





/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef h4REFERENCE_VITERBI_EXAMPLE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"

#include "general.h"
#include "reference_viterbi.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                        docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help",             0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show HMMER version info",     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of using the Viterbi reference implementation";


int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  H4_REFMX       *rmx     = h4_refmx_Create(100, 100);
  H4_PATH        *pi      = h4_path_Create();
  float           vsc;
  int             status;

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile");
  h4_hmmfile_Close(hfp);

  if ( esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, /*env=*/NULL, &sqfp) != eslOK) esl_fatal("couldn't open sequence file %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  h4_profile_Config(hmm);
  mo = h4_mode_Create();

  while (( status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      h4_mode_SetLength(mo, sq->n);

      h4_ReferenceViterbi(sq->dsq, sq->n, hmm, mo, rmx, pi, &vsc);

      printf("%s %.3f\n", sq->name, vsc);

      //h4_refmx_Dump(stdout, rmx);
      h4_path_Dump(stdout, pi);

      h4_refmx_Reuse(rmx);
      h4_path_Reuse(pi);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  h4_path_Destroy(pi);
  h4_refmx_Destroy(rmx);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
}

#endif // h4REFERENCE_VITERBI_EXAMPLE
