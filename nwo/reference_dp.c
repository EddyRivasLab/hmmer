/* Reference implementation of basic dynamic programming algorithms
 * for dual-mode glocal/local model.
 *
 * All reference implementation code is for testing. It is not used in
 * HMMER's main executables. Sparse DP code is the production version.
 *
 * Contents:
 *   1. Viterbi
 *   2. Forward
 *   3. Backward
 *   4. Posterior decoding
 *   5. Tracebacks (Viterbi and stochastic)
 *    ... a. Selection functions for Viterbi traceback
 *    ... b. Selection functions for stochastic traceback
 *    ... c. Traceback engine
 *    ... d. API calls (wrappers around engine)
 *   6. Benchmark driver
 *   7. "Brute" unit test
 *   8. Other unit tests
 *   9. Test driver
 *  10. Example
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_random.h"      // stochastic tracebacks
#include "esl_vectorops.h"

#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_refmx.h"
#include "logsum.h"

#include "reference_dp.h"


/*****************************************************************
 * 1. Viterbi
 *****************************************************************/

/* Function:  h4_reference_Viterbi()
 * Synopsis:  Reference implementation of the Viterbi algorithm.
 *
 * Purpose:   Given a target sequence <dsq> of length <L> residues, a
 *            query profile <hmm> in mode <mo>, and a DP
 *            matrix <rmx>; do the Viterbi optimal alignment
 *            algorithm.  Return the Viterbi score in bits in
 *            <*opt_sc>, if caller provides it.  Return the Viterbi
 *            optimal trace in <*opt_pi>, if caller provides an
 *            allocated trace structure.
 *            
 *            <rmx> is allocated by the caller but can be any size; it
 *            will be reallocated if needed, so it can be reused from
 *            a previous calculation, even a smaller one.
 *            
 * Args:      dsq     - digital target sequence 1..L
 *            L       - length of <dsq> in residues
 *            hmm     - query profile
 *            mo      - profile's mode (algorithm-dependent parameters)
 *            rmx     - DP matrix
 *            opt_pi  - optRETURN: path structure for Viterbi alignment, or NULL
 *            opt_sc  - optRETURN: Viterbi raw score in bits, or NULL
 *
 * Returns:   <eslOK> on success. <rmx> contains the Viterbi matrix.
 * 
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
h4_reference_Viterbi(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx, H4_PATH *opt_pi, float *opt_sc)
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
  if (( status = h4_refmx_GrowTo (rmx, M, L))              != eslOK) goto ERROR;
  if (( status = h4_refmx_SetType(rmx, M, L, h4R_VITERBI)) != eslOK) goto ERROR;

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
				 ESL_MAX( dpp[h4R_DG] + tsc[h4_DI],
					  xG          + tsc[h4_GI]));

	  /* E state update; local paths only, transition prob 1.0 in implicit probability model
           * You might think that Dk->E can never be optimal, and in Plan7 models this was so, 
           * but in Plan9 models there are edge cases where Ik->Dk+1->E is an optimal path.
	   */
	  xE  = ESL_MAX(xE, ESL_MAX(dpc[h4R_ML], dlv));

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

      /* E state update now includes glocal exits: transition prob 1.0 from MG_m
       * DLk->E and DGk->E both need to be included here; DGk->E is common, DLk->E is more weirdo-edge-casey */
      xE  = ESL_MAX( ESL_MAX( ESL_MAX( dpc[h4R_MG],
				       dgv),
		              ESL_MAX( dpc[h4R_ML],
				       dlv)),
   		                       xE);
      
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
  if (opt_pi && vsc != -eslINFINITY) return h4_reference_ViterbiTrace(hmm, mo, rmx, opt_pi);  // if no paths are possible at all: leave pi->N=0, our convention for impossible trace
  else                               return eslOK;

 ERROR:
  if (opt_pi) h4_path_Reuse(opt_pi);
  if (opt_sc) *opt_sc = -eslINFINITY;
  return status;
}
/*------------------ end, viterbi -------------------------------*/



/*****************************************************************
 * 2. Forward
 *****************************************************************/

/* Function:  h4_reference_Forward()
 * Synopsis:  Reference implementation of the Forward algorithm.
 *
 * Purpose:   The Forward algorithm, comparing profile <hmm> in mode
 *            <mo> to target sequence <dsq> of length <L>. Caller
 *            provides an allocated <H4_REFMX> DP matrix <rmx> of
 *            any size; this will be reallocated if needed, so it can be
 *            reused from any previous calculation, even a smaller
 *            one.
 *
 *            Caller also has initialized with a <h4_logsum_Init()>
 *            call, because this function will use <h4_logsum()>.
 *
 *            Upon successful return, the raw Forward score (in bits)
 *            is optionally returned in <*opt_sc>, and the DP matrix
 *            <rmx> is filled.
 *
 * Args:      dsq    : digital target sequence of length <L>
 *            L      : length of the target sequence
 *            hmm    : query profile
 *            mo     : algorithm-dependent parameters, the alignment "mode"
 *            rmx    : RESULT: DP matrix
 *            opt_sc : optRETURN: raw Forward score in bits
 *
 * Returns:   <eslOK> on success. <rmx> contains the Forward matrix;
 *            its internals may have been reallocated.
 *
 * Throws:    <eslEMEM> if reallocation is attempted and fails.
 *
 * Notes:     This function makes assumptions about the order
 *            of the state indices in h4_refmx.h:
 *              main states: ML MG IL IG DL DG
 *              specials:    E N J B L G C JJ CC
 *
 *            If L=0, the score is -infinity, by construction; HMMER
 *            profiles generate sequences of L>=1.
 */
int
h4_reference_Forward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx, float *opt_sc)
{
  int    M = hmm->M;
  float *dpc, *dpp;
  const float *tsc;	      // ptr for stepping thru profile's transition parameters 
  const float *rsc;	      // ptr for stepping thru profile's emission parameters   
  int    i, k, s;
  float  dlv, dgv; 	      // pushed-ahead DL,DG cell k+1 values
  float  xE, xL, xG;
  int    status;

  /* contract checks / arg validation */
  ESL_DASSERT1( ( mo->L == L || mo->L == 0) ); /* length model is either L (usually) or 0 (some unit tests) */

  /* reallocation, if needed */
  if (( status = h4_refmx_GrowTo (rmx, M, L))              != eslOK) goto ERROR;
  if (( status = h4_refmx_SetType(rmx, M, L, h4R_FORWARD)) != eslOK) goto ERROR; 

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
  /* *dpc is on the specials' subrow; must be there to make L=0 case work, where loop below is skipped */

  /* Main DP recursion */
  for (i = 1; i <= L; i++)
    {
      /* Initialization for a new row */
      rsc = hmm->rsc[dsq[i]] + 1;       /* this ptr steps through the row's emission scores 1..M. skip k=0 */
      tsc = hmm->tsc[0];                /* this ptr steps through profile's transition scores 0..M         */

      dpp = rmx->dp[i-1];               /* previous row dpp is already set, and at k=0 */
      dpc = rmx->dp[i];                 /* current DP row, skip k=0, start at k=1.  */
      for (s = 0; s < h4R_NSCELLS; s++) *dpc++ = -eslINFINITY;

      dlv = dgv = -eslINFINITY;
      xE  =       -eslINFINITY;

      /* Main inner loop of the recursion */
      for (k = 1; k < M; k++)
        {
          /* match states MLk, MGk */
          dpc[h4R_ML] = *rsc + h4_logsum( h4_logsum(dpp[h4R_ML] + tsc[h4_MM],
                                                    dpp[h4R_IL] + tsc[h4_IM]),
                                          h4_logsum(dpp[h4R_DL] + tsc[h4_DM],
                                                    xL          + tsc[h4_LM]));
          dpc[h4R_MG] = *rsc + h4_logsum( h4_logsum(dpp[h4R_MG] + tsc[h4_MM],
                                                    dpp[h4R_IG] + tsc[h4_IM]),
                                          h4_logsum(dpp[h4R_DG] + tsc[h4_DM],
                                                    xG          + tsc[h4_GM]));

          rsc++;                // rsc advances to next score (for position k+1)
          tsc += h4_NTSC;       // tsc advances to transitions in states k
          dpp += h4R_NSCELLS;   // dpp advances to cells for states k

          /* Insert state calculations ILk, IGk. */
          dpc[h4R_IL] = h4_logsum( h4_logsum( dpp[h4R_ML] + tsc[h4_MI],
                                              dpp[h4R_IL] + tsc[h4_II]),
                                              dpp[h4R_DL] + tsc[h4_DI]);
          dpc[h4R_IG] = h4_logsum( h4_logsum( dpp[h4R_MG] + tsc[h4_MI],
                                              dpp[h4R_IG] + tsc[h4_II]),
                                   h4_logsum( dpp[h4R_DG] + tsc[h4_DI],
					      xG          + tsc[h4_GI]));

          /* E state update; local paths only, DLk->E, MLk->E; transition prob 1.0 in implicit probability model */
          xE  = h4_logsum( xE, h4_logsum( dpc[h4R_ML], dlv ));  // dlv because dpc[h4R_DL] hasn't been stored yet, that's next

          /* Delete state, deferred storage trick */
          dpc[h4R_DL] = dlv;
          dpc[h4R_DG] = dgv;
          dlv = h4_logsum( h4_logsum(dpc[h4R_ML] + tsc[h4_MD],
                                     dpc[h4R_IL] + tsc[h4_ID]),
                                     dpc[h4R_DL] + tsc[h4_DD]);
          dgv = h4_logsum( h4_logsum(dpc[h4R_MG] + tsc[h4_MD],
                                     dpc[h4R_IG] + tsc[h4_ID]),
                                     dpc[h4R_DG] + tsc[h4_DD]);

          dpc += h4R_NSCELLS;
        }

      /* k=M node is unrolled and handled separately. No I state, and glocal exits. */
      dpc[h4R_ML] = *rsc + h4_logsum( h4_logsum(dpp[h4R_ML] + tsc[h4_MM],
                                                dpp[h4R_IL] + tsc[h4_IM]),
                                      h4_logsum(dpp[h4R_DL] + tsc[h4_DM],
                                                xL          + tsc[h4_LM]));
      dpc[h4R_MG] = *rsc + h4_logsum( h4_logsum(dpp[h4R_MG] + tsc[h4_MM],
                                                dpp[h4R_IG] + tsc[h4_IM]),
                                      h4_logsum(dpp[h4R_DG] + tsc[h4_DM],
                                                xG          + tsc[h4_GM]));
      dpp  += h4R_NSCELLS; 

      /* I_M state doesn't exist */
      dpc[h4R_IL] = -eslINFINITY;
      dpc[h4R_IG] = -eslINFINITY;


      /* E state update now includes glocal exits: transition prob 1.0 from MG_m, DG_m */
      xE  = h4_logsum( h4_logsum( h4_logsum(dpc[h4R_MG],
                                            dgv),
                                  h4_logsum(dpc[h4R_ML],
                                            dlv)),
                                            xE);

      /* D_M state: deferred storage only */
      dpc[h4R_DL] = dlv;
      dpc[h4R_DG] = dgv;
    
      /* row i is now finished; position both dpc, dpp on specials */
      dpc += h4R_NSCELLS;
      dpp += h4R_NSCELLS;   
      
      dpc[h4R_E]      = xE;		
      dpc[h4R_N]      =            dpp[h4R_N] + mo->xsc[h4_N][h4_LOOP]; 
      dpc[h4R_J]      = h4_logsum( dpp[h4R_J] + mo->xsc[h4_J][h4_LOOP], dpc[h4R_E] + mo->xsc[h4_E][h4_LOOP]); 
      dpc[h4R_B]      = h4_logsum( dpc[h4R_N] + mo->xsc[h4_N][h4_MOVE], dpc[h4R_J] + mo->xsc[h4_J][h4_MOVE]); 
      dpc[h4R_L] = xL =            dpc[h4R_B] + mo->xsc[h4_B][0]; 
      dpc[h4R_G] = xG =            dpc[h4R_B] + mo->xsc[h4_B][1]; 
      dpc[h4R_C]      = h4_logsum( dpp[h4R_C] + mo->xsc[h4_C][h4_LOOP], dpc[h4R_E] + mo->xsc[h4_E][h4_MOVE]);
      dpc[h4R_JJ]     = -eslINFINITY;                                                                           
      dpc[h4R_CC]     = -eslINFINITY;                                                                           
    }
  /* Done with all rows i. As we leave, dpc is still sitting on the final special value for i=L ... including even the L=0 case */
  
  if (opt_sc) *opt_sc = dpc[h4R_C] + mo->xsc[h4_C][h4_MOVE]; /* C->T */
  return eslOK;

 ERROR:
  if (opt_sc) *opt_sc = -eslINFINITY;
  return status;
}
/*-----------  end, Forwards implementation ---------------------*/


/*****************************************************************
 * 3. Backward
 *****************************************************************/

/* Function:  h4_reference_Backward()
 * Synopsis:  Reference implementation of the Backward algorithm
 *
 * Purpose:   The Backward algorithm, comparing profile <hmm> in mode
 *            <mo> to target sequence <dsq> of length <L>. Caller
 *            provides an allocated <P7_REFMX> DP matrix <rmx>; this
 *            matrix will be reallocated if needed, so it can be
 *            reused from a previous calculation, including a smaller
 *            one.
 *            
 *            Caller must have initialized with a <h4_logsum_Init()>
 *            call, because this function will use <h4_logsum()>.
 *            
 *            Upon successful return, the raw Backward score (in bits)
 *            is optionally returned in <*opt_sc>, and the DP matrix
 *            <rmx> is filled in.
 *
 * Args:      dsq    : digital target sequence of length <L>
 *            L      : length of the target sequence
 *            hmm    : query profile 
 *            mo     : algorithm-dependent parameters, the alignment "mode"
 *            rmx    : allocated DP matrix
 *            opt_sc : optRETURN: raw Backward score in bits
 *
 * Returns:   <eslOK> on success. <rmx> contains the Backward matrix;
 *            its internals may have been reallocated.
 *            
 * Throws:    <eslEMEM> if reallocation is attempted and fails.           
 *
 * Notes:     Order of evaluation in the code is carefully
 *            arranged to guarantee that dpc,dpn could be pointing
 *            into the same row of memory in a memory-efficient
 *            one-row DP implementation... even though this particular
 *            function, working with a <rmx>, knows it has all rows
 *            <0,1..L>. This is so this code could be cribbed for a
 *            one-row implementation.
 */
int
h4_reference_Backward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx, float *opt_sc)
{
  float *dpc;                   // ptr into current DP row, rmx->dp[i]                     
  float *dpn;                   // ptr into next DP row, rmx->dp[i+1]                       
  const float *rsn;             // ptr to next row's residue score x_{i+1} vector in <hmm> 
  const float *tsc;             // ptr to model transition score vector hmm->tsc[]         
  float dgc, dlc;               // DG,DL tmp values on current row, [i,?,DG], [i,?,DL]     
  float mgc, mlc;               // MG,ML tmp values on current row i
  float mgn, mln;               // MG,ML (+ emission) tmp values from next row i+1
  float ign, iln;               // IG/IL tmp values from next row i+1
  float xE, xG, xL;             // tmp vars for special state values
  int   i;                      // counter over sequence positions 1..L 
  int   k;                      // counter over model positions 1..M    
  const int M  = hmm->M;
  int   status;

 /* contract checks / arg validation */
  ESL_DASSERT1( ( mo->L == L || mo->L == 0) ); /* length model in profile is either L (usually) or 0 (some unit tests) */

  /* reallocation, if needed */
  if (( status = h4_refmx_GrowTo (rmx, M, L))               != eslOK) goto ERROR;
  if (( status = h4_refmx_SetType(rmx, M, L, h4R_BACKWARD)) != eslOK) goto ERROR;

  /* Initialize row L, starting with specials: */
  dpc = rmx->dp[L] + (M+1)*h4R_NSCELLS;        // position <dpc> on row L+1 specials
  dpc[h4R_CC]      = -eslINFINITY;                 
  dpc[h4R_JJ]      = -eslINFINITY;                 
  dpc[h4R_C]       = mo->xsc[h4_C][h4_MOVE];   // C <- T
  dpc[h4R_G]       = -eslINFINITY;                 
  dpc[h4R_L]       = -eslINFINITY;                 
  dpc[h4R_B]       = -eslINFINITY;                 
  dpc[h4R_J]       = -eslINFINITY;                 
  dpc[h4R_N]       = -eslINFINITY;                 
  dpc[h4R_E]  = xE = dpc[h4R_C] + mo->xsc[h4_E][h4_MOVE]; 
  dpc -= h4R_NSCELLS;                   // dpc now on M

  /* init supercell M on row L: */
  dpc[h4R_ML]       = xE;               // ML: M_M->E (transition prob 1.0)  
  dpc[h4R_MG]       = xE;               //   (same for MG)
  dpc[h4R_IL]       = -eslINFINITY;     // no ILm
  dpc[h4R_IG]       = -eslINFINITY;     // no IGm
  dpc[h4R_DL] = dlc = xE;               // DL: D_M->E (transition prob 1.0)  
  dpc[h4R_DG] = dgc = xE;               //   (same for DG)  
  dpc -= h4R_NSCELLS;                   // dpc now on M-1

  /* initialize main cells [k=1..M-1] on row L, which has no pulls
   * from I or M because there's no next row. 
   */
  tsc  = hmm->tsc[M-1];                // position <tsc> on node M-1    
  for (k = M-1; k >= 1; k--)
    {
      dpc[h4R_ML]       = h4_logsum(xE, dlc + tsc[h4_MD]); // ML -> DL or E              
      dpc[h4R_MG]       =               dgc + tsc[h4_MD];  // MG -> DG
      dpc[h4R_IL]       =               dlc + tsc[h4_ID];  // IL -> DL
      dpc[h4R_IG]       =               dgc + tsc[h4_ID];  // IG -> DG 
      dpc[h4R_DL] = dlc = h4_logsum(xE, dlc + tsc[h4_DD]); // DL: Dk->Dk+1 or Dk->E 
      dpc[h4R_DG] = dgc =               dgc + tsc[h4_DD];  // DG: only D->D path is possible 
      tsc -= h4_NTSC; 
      dpc -= h4R_NSCELLS;
    }
  /* k=0 cells are already -inf */

  /* The main recursion over rows i=L-1 down to 0.
   * On row i=0, we're going to break out early from the code below, after we
   * do the specials but before the main states.
   */
  for (i = L-1; i >= 0; i--)
    {
      /* Initialization of xL/xG by pulling from next row i+1.
       * Done by sweeping forward through k=1..M; with care taken to
       * leave ptrs sitting where they need to be for backwards sweeps.
       */
      xL  = xG = -eslINFINITY;
      tsc = hmm->tsc[0];
      dpn = rmx->dp[i+1] + h4R_NSCELLS;
      rsn = hmm->rsc[dsq[i+1]];
      for (k = 1; k <= M; k++)
        {
          rsn++;
          xL   = h4_logsum(xL, dpn[h4R_ML] + *rsn + tsc[h4_LM]);
          xG   = h4_logsum(xG, dpn[h4R_MG] + *rsn + tsc[h4_GM]);
          tsc += h4_NTSC;                                        // LMk and GMk entries are stored off-by-one, hence the delayed bump of tsc
          xG   = h4_logsum(xG, dpn[h4R_IG] + tsc[h4_GI]);
          dpn += h4R_NSCELLS; 
        }
      // <tsc>, <rsn> at M; and <dpn> at start of next row's specials (i.e. "M+1")

      /* Calculation of the special states on row i. */
      dpc = rmx->dp[i] + (M+1)*h4R_NSCELLS;                               // put dpc on start of current row's specials        
      dpc[h4R_CC]      = -eslINFINITY;                                    // CC unused (only used in decoding) 
      dpc[h4R_JJ]      = -eslINFINITY;                                    // JJ unused too 
      dpc[h4R_C]       = dpn[h4R_C] + mo->xsc[h4_C][h4_LOOP];             // C = C<-C 
      dpc[h4R_G]       = xG;                                              // G was just calculated above
      dpc[h4R_L]       = xL;                                              //  ... and L too
      dpc[h4R_B]       = h4_logsum(dpc[h4R_G] + mo->xsc[h4_B][h4_MOVE],   // B<-G 
                                   dpc[h4R_L] + mo->xsc[h4_B][h4_LOOP]);  // B<-L 
      dpc[h4R_J]       = h4_logsum(dpn[h4R_J] + mo->xsc[h4_J][h4_LOOP],   // J<-J 
                                   dpc[h4R_B] + mo->xsc[h4_J][h4_MOVE]);  // J<-B 
      dpc[h4R_N]       = h4_logsum(dpn[h4R_N] + mo->xsc[h4_N][h4_LOOP],   // N<-N 
                                   dpc[h4R_B] + mo->xsc[h4_N][h4_MOVE]);  // N<-B 
      dpc[h4R_E]  = xE = h4_logsum(dpc[h4R_C] + mo->xsc[h4_E][h4_MOVE],
                                   dpc[h4R_J] + mo->xsc[h4_E][h4_LOOP]);
      /* Forgive the weirdness, but this is the actual end of the
       * iteration over i.  When we reach row i=0, we need to
       * calculate the specials, exactly like we did for other rows;
       * but we don't calculate the rest of row i=0. Break out here.
       * Decided this is better than duplicating code.
       * Need dpc to remain on specials (M+1) so break before decrement.
       */
      if (i == 0) break;
      dpc -= h4R_NSCELLS;       /* dpc now on [i,M] supercell   */
      dpn -= h4R_NSCELLS;       /* dpn now on [i+1,M] supercell */


      /* Initialization of k=M supercell */
      dpc[h4R_ML]       = xE;           // MLm->E only
      dpc[h4R_MG]       = xE;           // MGm->E only
      dpc[h4R_IL]       = -eslINFINITY; // no ILm state
      dpc[h4R_IG]       = -eslINFINITY; // no IGm state
      dpc[h4R_DL] = dlc = xE;           // DLm->E, and store dlc i,k for next k loop iteration
      dpc[h4R_DG] = dgc = xE;           // DGm->E, and "" 

      mgn = *rsn + dpn[h4R_MG];         // pick up MG(i+1,k=M) + s(i+1,k=M) for next k loop iteration
      mln = *rsn + dpn[h4R_ML];
      rsn--;                            // now rsn on i+1,k=M-1
      tsc -= h4_NTSC;                   // now tsc on M-1
      dpc -= h4R_NSCELLS;               // now dpc on i,M-1   
      dpn -= h4R_NSCELLS;               // now dpn on i+1,M-1

      /* The main recursion over model positions k=M-1 down to 1. */
      for (k = M-1; k >= 1; k--)
        {
          ign = dpn[h4R_IG];                              // dpn is on k. pick up IG|IL from [i+1] so we can write into dpc[IG] 
          iln = dpn[h4R_IL];                              //   ... insert emission is zero
                                                          // tsc is on k.
          mlc =  h4_logsum( h4_logsum(mln + tsc[h4_MM],   // mln is i+1,k+1 plus i+1,k+1 emission
                                      iln + tsc[h4_MI]),  // iln is i+1,k just picked up
                            h4_logsum(dlc + tsc[h4_MD],   // dlc is i,k+1
                                      xE));               // ML->E trans = 1.0  
          mgc =  h4_logsum( h4_logsum(mgn + tsc[h4_MM],   // mgn is i+1,k+1 plus i+1,k+1 emission
                                      ign + tsc[h4_MI]),  // ign is i+1,k just picked up
                                      dgc + tsc[h4_MD]);  // dgc is i,k+1

          dpc[h4R_IL] = h4_logsum( h4_logsum(mln + tsc[h4_IM],           // safe to store in dpc[IG|IL] because we picked up ign|iln already
                                             iln + tsc[h4_II]),
                                             dlc + tsc[h4_ID]);
          dpc[h4R_IG] = h4_logsum( h4_logsum(mgn + tsc[h4_IM],     
                                             ign + tsc[h4_II]),
                                             dgc + tsc[h4_ID]);

          dpc[h4R_DL] = dlc = h4_logsum( h4_logsum(mln + tsc[h4_DM],     // dlc|dgc i,k picked up here, used in next loop iteration
                                                   iln + tsc[h4_DI]),
                                         h4_logsum(dlc + tsc[h4_DD],
                                                   xE));
          dpc[h4R_DG] = dgc = h4_logsum( h4_logsum(mgn + tsc[h4_DM],   
                                                   ign + tsc[h4_DI]),
                                                   dgc + tsc[h4_DD]);
          tsc -= h4_NTSC;               // now tsc is on k-1 supercell

                                        // rsn on i+1,k
          mgn = *rsn + dpn[h4R_MG];     // mgn i+1,k + emission i+1,k picked up here, used in next loop iteration
          mln = *rsn + dpn[h4R_ML];     
          rsn--;                        // rsn now on i+1,k-1
          dpc[h4R_MG] = mgc;            // delayed store of [i,k,MG] value enables dpc,dpn to point into same single row 
          dpc[h4R_ML] = mlc;

          dpn -= h4R_NSCELLS;           // dpn now on i+1,k-1 supercell
          dpc -= h4R_NSCELLS;           // dpc now on   i,k-1 supercell
          /* as we loop around now and decrement k:
           *   dpn [i+1,k-1] => [i+1,k] 
           *   dpc [i,k-1]   => [i,k] 
           *   tsc [k-1]     => [k]
           *   rsn [i+1,k-1] => [i+1,k]
           *   rsc [i,k-1]   => [i,k]
           *   dgc,dlc [i,k] => [i,k+1] 
           *   mgn,mln [i+1,k] => [i+1,k+1] including emission of i+1,k+1
           */
        } // end of loop over model positions k 
      // k=0 cells are -inf 
      // xG,xL values are now ready for deferred storage on next row above 
    } /* end of loop over rows i. */

  /* We just broke out of the iteration above on row i=0, having just
   * calculated the i=0 special states. Only N,B,G,L states are
   * ultimately reachable on i=0, but we calculate C,J,E anyway, to
   * match sparse implementation exactly. Decoding will see -inf for
   * these cells in Fwd matrix, and decode them to impossible. 
   */
  /* for complete cleanliness: set all the main states on row 0 to -inf */
  esl_vec_FSet(rmx->dp[0], (M+1)*h4R_NSCELLS, -eslINFINITY);

  if (opt_sc) *opt_sc = dpc[h4R_N]; // assumes dpc is sitting on i=0 specials (M+1)
  return eslOK;

 ERROR:
  if (opt_sc) *opt_sc = -eslINFINITY;
  return status;
}
/*-------------- end, backwards implementation ------------------*/



/*****************************************************************
 * 4. Posterior decoding
 *****************************************************************/

/* Function:  h4_reference_Decoding()
 * Synopsis:  Reference implementation of posterior decoding.
 *
 * Purpose:   Given previously calculated Forward and Backward matrices
 *            <fwd> and <bck>, for a comparison of query profile <hmm>
 *            in mode <mo> to digital sequence <dsq> of length <L>,
 *            perform posterior decoding for all states. 
 *
 *            The resulting posterior decoding matrix is left in <pp>,
 *            which has been provided by the caller. The caller can
 *            provide any allocated <pp> regardless of size; it will
 *            be reallocated here as needed. The caller can also
 *            provide <bck> itself, in which case <bck> is overwritten
 *            with posterior decoding probabilities.
 *
 *            <pp(i,y)> is the probability that state y is used on row
 *            i (except for DG states). Values for DG states are
 *            expected counts, not probabilities, due to the "mute
 *            partial cycle" issue. When <y> is an emitting state, it
 *            <pp(i,y> is the probability that state <y> emits residue
 *            <x_i>. The sum over emitting y's on any given row i=1..L
 *            is 1.0.
 *            
 * Args:      dsq - digital sequence, 1..L
 *            L   - length of <dsq>
 *            hmm - query profile
 *            mo  - alignment mode for the profile/seq comparison
 *            fwd - Forward matrix, provided by caller
 *            bck - Backward matrix, provided by caller
 *            pp  - RESULT: posterior decoding matrix
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> if a reallocation fails.
 * 
 * Notes:     For more information on the "mute partial cycle" issue, see
 *            reference_dp.md; and see the
 *            `utest_mute_partial_cycle()` unit test that exercises
 *            the issue.
 */
int
h4_reference_Decoding(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_REFMX *fwd, const H4_REFMX *bck, H4_REFMX *pp)
{
  const float *rsc;		// ptr to current row's residue score x_i vector in <hmm> 
  const int    M = hmm->M;
  float     xJ, xC, xG;         // tmp vars for holding J/C/G vals from previous row iteration (prv row i-1)
  float    *fwdp;               // ptr that will step through curr Forward values on row i 
  float    *bckp;               //  ... ditto for Backward
  float    *ppp;                //  ... and for decoding. bckp and ppp may be pointing to same memory; care is taken to access bckp before writing ppp
  float     denom;              // sum of pp for all emitting states on row i; should be 1.0; used for renormalization of roundoff error accumulation.
  float     sc;                 // overall Forward/Backward score, used for normalization
  int       i,k;                // indices for rows/residues, columns/profile positions
  float     delta;		// additions to DGk's, resulting from unfolding wing-retracted entry/exit paths 
  int       status;

  /* Contract checks / argument validation */
  ESL_DASSERT1( (fwd->type == h4R_FORWARD) );
  ESL_DASSERT1( (bck->type == h4R_BACKWARD) );
  ESL_DASSERT1( (fwd->L == bck->L && fwd->M == bck->M && fwd->M == hmm->M) );

  /* Reallocation, if needed. Be careful not to reallocate pp when it's the same as bck. */
  if (bck != pp && ( status = h4_refmx_GrowTo(pp, M, L))   != eslOK) goto ERROR;
  if (( status = h4_refmx_SetType(pp, M, L, h4R_DECODING)) != eslOK) goto ERROR;

  /* On row 0, all main states are 0; initialize them so. (leave k=0 col at -inf)
   * N is 1.0 by definition.
   * Some specials (BLG) can be reached, calculate w/ fwdp,bckp on specials.
   */
  esl_vec_FSet(pp->dp[0]+h4R_NSCELLS, M*h4R_NSCELLS, 0.0);

  ppp  = pp->dp[0]  + (h4R_NSCELLS * (M+1)); // fwdp, bckp,ppp on first special of row 0 
  fwdp = fwd->dp[0] + (h4R_NSCELLS * (M+1)); 
  bckp = bck->dp[0] + (h4R_NSCELLS * (M+1)); 
  sc   = bckp[h4R_N];                        // bck[0][N] is the tot score. Use it to normalize 
  ppp[h4R_E]  = 0.0f;
  ppp[h4R_N]  = 1.0f;
  ppp[h4R_J]  = 0.0f;
  ppp[h4R_B]  = exp2f(fwdp[h4R_B] + bckp[h4R_B] - sc);
  ppp[h4R_L]  = exp2f(fwdp[h4R_L] + bckp[h4R_L] - sc);
  ppp[h4R_G]  = exp2f(fwdp[h4R_G] + bckp[h4R_G] - sc);
  ppp[h4R_C]  = 0.0f;
  ppp[h4R_JJ] = 0.0f;
  ppp[h4R_CC] = 0.0f;
  
  /* xJ/xC/xG hold the previous row i-1 values from forward matrix;
   * needed for decoding emit-on-transition CC/JJ, and G->D1..Dk-1->Mk entry
   */
  xJ = fwdp[h4R_J];     /* i.e. -inf */
  xC = fwdp[h4R_C];	/* i.e. -inf */
  xG = fwdp[h4R_G];

  /* main recursion */
  for (i = 1; i <= L; i++)
    {
      /* Wing retraction contributes posterior probability terms to DG's on _previous_ row i-1.
       * Calculation uses bckp[MG|IG] terms on current row i, so it must come before ppp[MG|IG] are written.
       * Calculation decodes the _transition_, so it needs an <*rsc> emission term.
       */
      bckp  = bck->dp[i]  + M*h4R_NSCELLS;
      ppp   = pp->dp[i-1] + M*h4R_NSCELLS;
      rsc   = hmm->rsc[dsq[i]] + M;
      delta = 0.;
      for (k = M; k >= 1; k--, rsc--, ppp-=h4R_NSCELLS)
	{
	  ppp[h4R_DG] += delta;                                                         // bump prev row i-1 DGk with G->DG1...DGk->{IGk,MGk+1} entry
	  delta       += exp2f(xG + hmm->tsc[k-1][h4_GM] + *rsc + bckp[h4R_MG] - sc);   // G->Mk path will contribute pp delta to D1..Dk-1
	  bckp        -= h4R_NSCELLS;                                                   // bckp now on i,k-1
	  delta       += exp2f(xG + hmm->tsc[k-1][h4_GI]        + bckp[h4R_IG] - sc);   // G->Ik-1 path will contribute pp delta to D1..Dk-1 too
	}  // bckp, rsc now on i,k=0
      

      /* Main decoding. */
      fwdp  = fwd->dp[i] + h4R_NSCELLS;               // fwdp on i,k=1
      bckp += h4R_NSCELLS;                            // bckp on i,k=1
      ppp   = pp->dp[i]  + h4R_NSCELLS;               // ppp  on i,k=1
      denom = 0.0;
      for (k = 1; k <= M; k++, ppp+= h4R_NSCELLS, fwdp+= h4R_NSCELLS, bckp += h4R_NSCELLS )
	{
	  ppp[h4R_ML] = exp2f(fwdp[h4R_ML] + bckp[h4R_ML] - sc); denom += ppp[h4R_ML]; 
	  ppp[h4R_MG] = exp2f(fwdp[h4R_MG] + bckp[h4R_MG] - sc); denom += ppp[h4R_MG]; 
	  ppp[h4R_IL] = exp2f(fwdp[h4R_IL] + bckp[h4R_IL] - sc); denom += ppp[h4R_IL];  // at k=M IL=0.0, because fwd/bck are -inf
	  ppp[h4R_IG] = exp2f(fwdp[h4R_IG] + bckp[h4R_IG] - sc); denom += ppp[h4R_IG];  //   ... (ditto for IG)
	  ppp[h4R_DL] = exp2f(fwdp[h4R_DL] + bckp[h4R_DL] - sc);                  
	  ppp[h4R_DG] = exp2f(fwdp[h4R_DG] + bckp[h4R_DG] - sc);                  
	} // at exit, ppp/fwdp/bckp are at start of specials
      

      /* [ E N J B L G C JJ CC ] */
      /* JJ, CC must be done first; they access bckp[J,C], which ppp[J,C] calc will clobber */
      ppp[h4R_JJ] = exp2f(xJ + mo->xsc[h4_J][h4_LOOP] + bckp[h4R_J] - sc); xJ = fwdp[h4R_J]; denom += ppp[h4R_JJ];
      ppp[h4R_CC] = exp2f(xC + mo->xsc[h4_C][h4_LOOP] + bckp[h4R_C] - sc); xC = fwdp[h4R_C]; denom += ppp[h4R_CC];
      ppp[h4R_E]  = exp2f(fwdp[h4R_E] + bckp[h4R_E] - sc);                   
      ppp[h4R_N]  = exp2f(fwdp[h4R_N] + bckp[h4R_N] - sc);                                   denom += ppp[h4R_N]; /* only NN is possible for i>=1, so N=NN */
      ppp[h4R_J]  = exp2f(fwdp[h4R_J] + bckp[h4R_J] - sc);
      ppp[h4R_B]  = exp2f(fwdp[h4R_B] + bckp[h4R_B] - sc);
      ppp[h4R_L]  = exp2f(fwdp[h4R_L] + bckp[h4R_L] - sc);
      ppp[h4R_G]  = exp2f(fwdp[h4R_G] + bckp[h4R_G] - sc);                xG = fwdp[h4R_G];
      ppp[h4R_C]  = exp2f(fwdp[h4R_C] + bckp[h4R_C] - sc);


      /* Renormalization.
       * 
       * Roundoff error accumulation in F/B is significant. For large
       * target seqs, it isn't unusual to have a whole nat of
       * difference in overall fwd vs. bck raw score; for example, fn3
       * vs. TITIN_HUMAN. Since we only use the bck score (at i=0) for
       * normalization above, pp's would have very large systematic
       * error. 
       * 
       * To squash this error in production code, we renormalize.
       * Because renormalization can hide real errors, we don't
       * renormalize when we've compiled the code for unit testing.
       * Default unit tests don't run large enough M/L to create a
       * lot of error accumulation.
       */
#ifndef h4REFERENCE_DP_TESTDRIVE
      esl_vec_FScale(pp->dp[i], (M+1)*h4R_NSCELLS + h4R_NXCELLS, 1.0/denom);
#endif
    }
  return eslOK;

 ERROR:
  return status;
}
/*-------------- end, posterior decoding ------------------------*/





/*****************************************************************
 * 5. Tracebacks (Viterbi and stochastic)
 *****************************************************************/

/* Traceback routines for Viterbi, MEG alignment, and stochastic
 * sampling, for the reference implementation. All tracebacks use the
 * same machinery (the reference_traceback_engine()), using different
 * select*() functions.
 */


/********* 5a. Selection functions for Viterbi tracebacks *******/

static inline int8_t
v_select_ml(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxv, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[4] = { h4P_ML, h4P_IL, h4P_DL, h4P_L };
  float  path[4];

  path[0] = H4R_MX (rxv, i-1, k-1, h4R_ML) + hmm->tsc[k-1][h4_MM];
  path[1] = H4R_MX (rxv, i-1, k-1, h4R_IL) + hmm->tsc[k-1][h4_IM];
  path[2] = H4R_MX (rxv, i-1, k-1, h4R_DL) + hmm->tsc[k-1][h4_DM];
  path[3] = H4R_XMX(rxv, i-1,      h4R_L)  + hmm->tsc[k-1][h4_LM];
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int8_t
v_select_mg(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxv, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[4] = { h4P_MG, h4P_IG, h4P_DG, h4P_G };
  float  path[4];

  path[0] = H4R_MX (rxv, i-1, k-1, h4R_MG) + hmm->tsc[k-1][h4_MM];
  path[1] = H4R_MX (rxv, i-1, k-1, h4R_IG) + hmm->tsc[k-1][h4_IM];
  path[2] = H4R_MX (rxv, i-1, k-1, h4R_DG) + hmm->tsc[k-1][h4_DM];
  path[3] = H4R_XMX(rxv, i-1,      h4R_G)  + hmm->tsc[k-1][h4_GM];
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int8_t
v_select_il(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxv, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[3] = { h4P_ML, h4P_IL, h4P_DL };
  float  path[3];

  path[0] = H4R_MX(rxv, i-1, k, h4R_ML) + hmm->tsc[k][h4_MI];
  path[1] = H4R_MX(rxv, i-1, k, h4R_IL) + hmm->tsc[k][h4_II];
  path[2] = H4R_MX(rxv, i-1, k, h4R_DL) + hmm->tsc[k][h4_DI];
  return state[esl_vec_FArgMax(path, 3)];
}

static inline int8_t
v_select_ig(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxv, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[4] = { h4P_MG, h4P_IG, h4P_DG, h4P_G };
  float  path[4];

  path[0] = H4R_MX (rxv, i-1, k, h4R_MG) + hmm->tsc[k][h4_MI];
  path[1] = H4R_MX (rxv, i-1, k, h4R_IG) + hmm->tsc[k][h4_II];
  path[2] = H4R_MX (rxv, i-1, k, h4R_DG) + hmm->tsc[k][h4_DI];
  path[3] = H4R_XMX(rxv, i-1,    h4R_G)  + hmm->tsc[k][h4_GI];
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int8_t
v_select_dl(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxv, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[3] = { h4P_ML, h4P_IL, h4P_DL };
  float  path[3];

  path[0] = H4R_MX(rxv, i, k-1, h4R_ML) + hmm->tsc[k-1][h4_MD];
  path[1] = H4R_MX(rxv, i, k-1, h4R_IL) + hmm->tsc[k-1][h4_ID];
  path[2] = H4R_MX(rxv, i, k-1, h4R_DL) + hmm->tsc[k-1][h4_DD];
  return state[esl_vec_FArgMax(path, 3)];
}

static inline int8_t
v_select_dg(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxv, int i, int k)
{
  ESL_UNUSED(r);
  int8_t state[3] = { h4P_MG, h4P_IG, h4P_DG };
  float  path[3];

  path[0] = H4R_MX(rxv, i, k-1, h4R_MG) + hmm->tsc[k-1][h4_MD];
  path[1] = H4R_MX(rxv, i, k-1, h4R_IG) + hmm->tsc[k-1][h4_ID];
  path[2] = H4R_MX(rxv, i, k-1, h4R_DG) + hmm->tsc[k-1][h4_DD];
  return state[esl_vec_FArgMax(path, 3)];
}

static inline int8_t
v_select_e(ESL_RANDOMNESS *r, float *wrk, const H4_PROFILE *hmm, const H4_REFMX *rxv, int i, int *ret_k)
{
  ESL_UNUSED(r);
  ESL_UNUSED(wrk);
  float  max  = -eslINFINITY;
  int8_t smax = -1;
  int    kmax = -1;
  int    k;

  for (k = 1; k <= hmm->M; k++)
    {
      if (H4R_MX(rxv, i, k, h4R_ML) > max) { max = H4R_MX(rxv, i, k, h4R_ML); smax = h4P_ML; kmax = k; }
      if (H4R_MX(rxv, i, k, h4R_DL) > max) { max = H4R_MX(rxv, i, k, h4R_DL); smax = h4P_DL; kmax = k; }
    }

  if (H4R_MX(rxv, i, hmm->M, h4R_MG) > max) { max = H4R_MX(rxv, i, hmm->M, h4R_MG); smax = h4P_MG; kmax = hmm->M; }
  if (H4R_MX(rxv, i, hmm->M, h4R_DG) > max) {                                       smax = h4P_DG; kmax = hmm->M; }

  *ret_k = kmax;
  return smax;
}

static inline int8_t
v_select_j(ESL_RANDOMNESS *r, const H4_MODE *mo, const H4_REFMX *rxv, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = H4R_XMX(rxv, i-1, h4R_J) + mo->xsc[h4_J][h4_LOOP];
  path[1] = H4R_XMX(rxv, i,   h4R_E) + mo->xsc[h4_E][h4_LOOP];
  return ( (path[0] > path[1]) ? h4P_J : h4P_E);
}

static inline int8_t
v_select_b( ESL_RANDOMNESS *r, const H4_MODE *mo, const H4_REFMX *rxv, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = H4R_XMX(rxv, i, h4R_N) + mo->xsc[h4_N][h4_MOVE];
  path[1] = H4R_XMX(rxv, i, h4R_J) + mo->xsc[h4_J][h4_MOVE];
  return ( (path[0] > path[1]) ? h4P_N : h4P_J);
}

static inline int8_t
v_select_c(ESL_RANDOMNESS *r, const H4_MODE *mo, const H4_REFMX *rxv, int i)
{
  ESL_UNUSED(r);
  float path[2];
  path[0] = H4R_XMX(rxv, i-1, h4R_C) + mo->xsc[h4_C][h4_LOOP];
  path[1] = H4R_XMX(rxv, i,   h4R_E) + mo->xsc[h4_E][h4_MOVE];
  return ( (path[0] > path[1]) ? h4P_C : h4P_E);
}



/****** 5b. Selection functions for stochastic tracebacks ********/

static inline int8_t
sto_select_ml(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxf, int i, int k)
{
  int8_t state[4] = { h4P_ML, h4P_IL, h4P_DL, h4P_L };
  float  path[4];

  path[0] = H4R_MX (rxf, i-1, k-1, h4R_ML) + hmm->tsc[k-1][h4_MM];
  path[1] = H4R_MX (rxf, i-1, k-1, h4R_IL) + hmm->tsc[k-1][h4_IM];
  path[2] = H4R_MX (rxf, i-1, k-1, h4R_DL) + hmm->tsc[k-1][h4_DM];
  path[3] = H4R_XMX(rxf, i-1,      h4R_L)  + hmm->tsc[k-1][h4_LM];
  esl_vec_FLog2Norm(path,4);
  return state[esl_rnd_FChoose(r, path, 4)];
}

static inline int8_t
sto_select_mg(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxf, int i, int k)
{
  int8_t state[4] = { h4P_MG, h4P_IG, h4P_DG, h4P_G };
  float  path[4];

  path[0] = H4R_MX (rxf, i-1, k-1, h4R_MG) + hmm->tsc[k-1][h4_MM];
  path[1] = H4R_MX (rxf, i-1, k-1, h4R_IG) + hmm->tsc[k-1][h4_IM];
  path[2] = H4R_MX (rxf, i-1, k-1, h4R_DG) + hmm->tsc[k-1][h4_DM];
  path[3] = H4R_XMX(rxf, i-1,      h4R_G)  + hmm->tsc[k-1][h4_GM];
  esl_vec_FLog2Norm(path,4);
  return state[esl_rnd_FChoose(r, path, 4)];
}

static inline int8_t
sto_select_il(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxf, int i, int k)
{
  int8_t state[3] = { h4P_ML, h4P_IL, h4P_DL };
  float  path[3];

  path[0] = H4R_MX(rxf, i-1, k, h4R_ML) + hmm->tsc[k][h4_MI];
  path[1] = H4R_MX(rxf, i-1, k, h4R_IL) + hmm->tsc[k][h4_II];
  path[2] = H4R_MX(rxf, i-1, k, h4R_DL) + hmm->tsc[k][h4_DI];
  esl_vec_FLog2Norm(path, 3);
  return state[esl_rnd_FChoose(r, path, 3)];
}

static inline int8_t
sto_select_ig(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxf, int i, int k)
{
  int8_t state[4] = { h4P_MG, h4P_IG, h4P_DG, h4P_G };
  float  path[4];

  path[0] = H4R_MX (rxf, i-1, k, h4R_MG) + hmm->tsc[k][h4_MI];
  path[1] = H4R_MX (rxf, i-1, k, h4R_IG) + hmm->tsc[k][h4_II];
  path[2] = H4R_MX (rxf, i-1, k, h4R_DG) + hmm->tsc[k][h4_DI];
  path[3] = H4R_XMX(rxf, i-1,    h4R_G)  + hmm->tsc[k][h4_GI];
  esl_vec_FLog2Norm(path, 4);
  return state[esl_rnd_FChoose(r, path, 4)];
}

static inline int8_t
sto_select_dl(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxf, int i, int k)
{
  int8_t state[3] = { h4P_ML, h4P_IL, h4P_DL };
  float  path[3];

  path[0] = H4R_MX(rxf, i, k-1, h4R_ML) + hmm->tsc[k-1][h4_MD];
  path[1] = H4R_MX(rxf, i, k-1, h4R_IL) + hmm->tsc[k-1][h4_ID];
  path[2] = H4R_MX(rxf, i, k-1, h4R_DL) + hmm->tsc[k-1][h4_DD];
  esl_vec_FLog2Norm(path, 3);
  return state[esl_rnd_FChoose(r, path, 3)];
}

static inline int8_t
sto_select_dg(ESL_RANDOMNESS *r, const H4_PROFILE *hmm, const H4_REFMX *rxf, int i, int k)
{
  int8_t state[3] = { h4P_MG, h4P_IG, h4P_DG };
  float  path[3];

  path[0] = H4R_MX(rxf, i, k-1, h4R_MG) + hmm->tsc[k-1][h4_MD];
  path[1] = H4R_MX(rxf, i, k-1, h4R_IG) + hmm->tsc[k-1][h4_ID];
  path[2] = H4R_MX(rxf, i, k-1, h4R_DG) + hmm->tsc[k-1][h4_DD];
  esl_vec_FLog2Norm(path, 3);
  return state[esl_rnd_FChoose(r, path, 3)];
}

static inline int8_t
sto_select_e(ESL_RANDOMNESS *r, float *wrk, const H4_PROFILE *hmm, const H4_REFMX *rxf, int i, int *ret_k)
{
  int   M = hmm->M;
  int   k;

  wrk[0]                             = H4R_MX(rxf, i, M, h4R_DG); // 0       : DGm->E exit           
  for (k = 1; k <= M; k++)  wrk[k]   = H4R_MX(rxf, i, k, h4R_DL); // 1..M    : DLk->E exits for k=1..M 
  for (k = 1; k <= M; k++)  wrk[k+M] = H4R_MX(rxf, i, k, h4R_ML); // M+1..2M : MLk->E exits for k=1..M 
  wrk[2*M+1]                         = H4R_MX(rxf, i, M, h4R_MG); // 2M+1    : MGk->E exit 

  esl_vec_FLog2Norm(wrk, 2*M+2);
  k = esl_rnd_FChoose(r, wrk, 2*M+2);
  if      (k == 2*M+1) { *ret_k = M;   return h4P_MG; }
  else if (k == 0)     { *ret_k = M;   return h4P_DG; }
  else if (k <= M)     { *ret_k = k;   return h4P_DL; }
  else                 { *ret_k = k-M; return h4P_ML; }
}

static inline int8_t
sto_select_j(ESL_RANDOMNESS *r, const H4_MODE *mo, const H4_REFMX *rxf, int i)
{
  float path[2];
  path[0] = H4R_XMX(rxf, i-1, h4R_J) + mo->xsc[h4_J][h4_LOOP];
  path[1] = H4R_XMX(rxf, i,   h4R_E) + mo->xsc[h4_E][h4_LOOP];
  esl_vec_FLog2Norm(path, 2);
  return (  (esl_rnd_FChoose(r,path,2) == 0) ? h4P_J : h4P_E);
}

static inline int8_t
sto_select_b( ESL_RANDOMNESS *r, const H4_MODE *mo, const H4_REFMX *rxf, int i)
{
  float path[2];
  path[0] = H4R_XMX(rxf, i, h4R_N) + mo->xsc[h4_N][h4_MOVE];
  path[1] = H4R_XMX(rxf, i, h4R_J) + mo->xsc[h4_J][h4_MOVE];
  esl_vec_FLog2Norm(path, 2);
  return (  (esl_rnd_FChoose(r,path,2) == 0) ? h4P_N : h4P_J);
}

static inline int8_t
sto_select_c(ESL_RANDOMNESS *r, const H4_MODE *mo, const H4_REFMX *rxf, int i)
{
  float path[2];
  path[0] = H4R_XMX(rxf, i-1, h4R_C) + mo->xsc[h4_C][h4_LOOP];
  path[1] = H4R_XMX(rxf, i,   h4R_E) + mo->xsc[h4_E][h4_MOVE];
  esl_vec_FLog2Norm(path, 2);
  return (  (esl_rnd_FChoose(r,path,2) == 0) ? h4P_C : h4P_E);
}



/***************** 5c. Traceback engine **************************/

static int 
reference_trace_engine(ESL_RANDOMNESS *rng, float *wrk, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_REFMX *rmx, H4_PATH *pi)
{
  int    i    = rmx->L;
  int    k    = 0;
  int8_t sprv = h4P_C;
  int8_t scur;
  int    status;

  int8_t (*select_ml)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_mg)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_il)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_ig)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_dl)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_dg)(ESL_RANDOMNESS *,          const H4_PROFILE *, const H4_REFMX *, int, int);
  int8_t (*select_e) (ESL_RANDOMNESS *, float *, const H4_PROFILE *, const H4_REFMX *, int, int *);
  int8_t (*select_j) (ESL_RANDOMNESS *,          const H4_MODE *,    const H4_REFMX *, int);
  int8_t (*select_b) (ESL_RANDOMNESS *,          const H4_MODE *,    const H4_REFMX *, int);
  int8_t (*select_c) (ESL_RANDOMNESS *,          const H4_MODE *,    const H4_REFMX *, int);
  
  /* if the target sequence is impossible: leave trace alone, with tr->Z=0, our convention for an impossible trace */
  if (H4R_XMX(rmx, rmx->L, h4R_C) == -eslINFINITY)  return eslOK; 

  if (rmx->type == h4R_VITERBI)	/* configure for Viterbi optimal traceback */
    {
      select_ml = &v_select_ml;      select_mg = &v_select_mg;
      select_il = &v_select_il;      select_ig = &v_select_ig;
      select_dl = &v_select_dl;      select_dg = &v_select_dg;
      select_j  = &v_select_j;       select_c  = &v_select_c;
      select_e  = &v_select_e;       select_b  = &v_select_b;
    }
  else if (rmx->type == h4R_FORWARD) /* configure for stochastic traceback */
    {
      select_ml = &sto_select_ml;      select_mg = &sto_select_mg;
      select_il = &sto_select_il;      select_ig = &sto_select_ig;
      select_dl = &sto_select_dl;      select_dg = &sto_select_dg;
      select_j  = &sto_select_j;       select_c  = &sto_select_c;
      select_e  = &sto_select_e;       select_b  = &sto_select_b;
    }

  if ((status = h4_path_Append(pi, h4P_C)) != eslOK) return status;
  while (sprv != h4P_S)
    {
      switch (sprv) {
      case h4P_ML: scur = (*select_ml)(rng,      hmm, rmx, i, k); k--; i--; break;
      case h4P_MG: scur = (*select_mg)(rng,      hmm, rmx, i, k); k--; i--; break;
      case h4P_IL: scur = (*select_il)(rng,      hmm, rmx, i, k);      i--; break;
      case h4P_IG: scur = (*select_ig)(rng,      hmm, rmx, i, k);      i--; break;
      case h4P_DL: scur = (*select_dl)(rng,      hmm, rmx, i, k); k--;      break;
      case h4P_DG: scur = (*select_dg)(rng,      hmm, rmx, i, k); k--;      break;
      case h4P_E:  scur = (*select_e) (rng, wrk, hmm, rmx, i, &k);          break;
      case h4P_N:  scur = (i==0 ? h4P_S : h4P_N);                           break;
      case h4P_J:  scur = (*select_j) (rng, mo, rmx, i);                    break;
      case h4P_B:  scur = (*select_b) (rng, mo, rmx, i);                    break;
      case h4P_L:  scur = h4P_B;                                            break;
      case h4P_G:  scur = h4P_B;                                            break;
      case h4P_C:  scur = (*select_c) (rng, mo, rmx, i);                    break;
      default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
      }

      /* A glocal B->G->{M|I}k wing-retraction entry: unfold it 
       */
      if (scur == h4P_G) {
	while (k--) { if ( (status = h4_path_Append(pi, h4P_DG)) != eslOK) return status; }
      }

      if (scur == h4P_L) status = h4_path_AppendElement(pi, scur, k+1); // if we just backed into an L state, record the k of L->Mk in its rle[].
      else               status = h4_path_Append       (pi, scur);
      if (status != eslOK) return status;

      if ( (scur == h4P_N || scur == h4P_J || scur == h4P_C) && scur == sprv) i--;
      sprv = scur;
    }
  
  return h4_path_Reverse(pi);
}

/****** 5d. Exposed API calls: wrappers around the engine ********/


/* Function:  h4_reference_ViterbiTrace()
 * Synopsis:  Traceback an optimal path from a Viterbi DP matrix
 * Incept:    SRE, Tue 21 May 2019
 *
 * Purpose:   Given Viterbi matrix <rxv> that was computed for a comparison of
 *            profile <hmm> in mode <mo> to a sequence, get the optimal 
 *            traceback path and return it in <pi>.
 *            
 *            Caller provides an allocated and empty <pi>.
 */
int
h4_reference_ViterbiTrace(const H4_PROFILE *hmm, const H4_MODE *mo, const H4_REFMX *rxv, H4_PATH *pi)
{
  ESL_DASSERT1(( pi->Z == 0 ));   // caller has create'd or reuse'd the path
  ESL_DASSERT1(( hmm->M == rxv->M ));
  ESL_DASSERT1(( rxv->type == h4R_VITERBI ));

  return reference_trace_engine(NULL, NULL, hmm, mo, rxv, pi);
}


/* Function:  h4_reference_StochasticTrace()
 * Synopsis:  Sample a path from posterior distribution from a Forward DP matrix
 * Incept:    SRE, Tue 21 May 2019
 *
 * Purpose:   Given Forward matrix <rxf> that was computed for a
 *            comparison of profile <hmm> in mode <mo> to a sequence,
 *            use random number generator <rng> to do a stochastic
 *            traceback to sample a path from the posterior
 *            distribution; return that path in <pi>.
 *            
 *            Requires a tmp workspace vector, at least 2M+2 floats
 *            long, that the caller may optionally provide in <wrk>.
 *            <wrk> follows the Easel "bypass" idiom. If caller passes
 *            <NULL>, workspace will be internally allocated and
 *            managed. If caller passes <&wrk> for <float *wrk =
 *            NULL>, then <*wrk> is allocated here for 2M+2 floats,
 *            and returned; caller can keep passing it to subsequent
 *            calls. Or caller can allocate <float *wrk> for at least
 *            2M+2 floats itself and pass <&wrk>.
 *
 * Args:      rng      - random number generator to use
 *            wrk_byp  - NULL, or ptr to <float *> tmp workspace of at least 2M+2 floats
 *            hmm      - profile HMM
 *            mo       - alignment mode
 *            rxf      - forward matrix
 *            pi       - RESULT: sampled path
 */
int
h4_reference_StochasticTrace(ESL_RANDOMNESS *rng, float **wrk_byp, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_REFMX *rxf, H4_PATH *pi)
{
  float *wrk = NULL;
  int    status;

  ESL_DASSERT1(( pi->Z == 0 ));   // caller has create'd or reuse'd the path
  ESL_DASSERT1(( hmm->M == rxf->M ));
  ESL_DASSERT1(( rxf->type == h4R_FORWARD ));

  if      (esl_byp_IsInternal(wrk_byp)) { ESL_ALLOC  (wrk,      (2*hmm->M+2) * sizeof(float));                 }
  else if (esl_byp_IsReturned(wrk_byp)) { ESL_ALLOC  (*wrk_byp, (2*hmm->M+2) * sizeof(float)); wrk = *wrk_byp; }
  else if (esl_byp_IsProvided(wrk_byp)) { ESL_REALLOC(*wrk_byp, (2*hmm->M+2) * sizeof(float)); wrk = *wrk_byp; }

  status = reference_trace_engine(rng, wrk, hmm, mo, rxf, pi);

 ERROR: // also normal exit:
  if  (esl_byp_IsInternal(wrk_byp)) { free(wrk); }
  return status;
}


/*---------------- end, tracebacks ------------------------------*/

/*****************************************************************
 * 6. Benchmark driver
 *****************************************************************/
#ifdef h4REFERENCE_DP_BENCHMARK

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "h4_profile.h"
#include "h4_refmx.h"

#include "general.h"
#include "logsum.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,   "2000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "--version", eslARG_NONE,    NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version/release information",         0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 1, argc, argv,
						"benchmark driver for generic dual-mode Forward/Backward DP",
						"[-options] <hmmfile>");
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  H4_REFMX       *rx      = h4_refmx_Create(100, 100);  // will be resized as needed.
  H4_PATH        *pi      = h4_path_Create();
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, F_time, B_time, V_time;
  double          F_speed, B_speed, V_speed;

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile");
  h4_hmmfile_Close(hfp);
  h4_mode_SetLength(mo, L);

  /* Baseline time for seq generation alone */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(rng, hmm->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Forward benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(rng, hmm->f, abc->K, L, dsq);
      h4_reference_Forward (dsq, L, hmm, mo, rx, &sc);
      h4_refmx_Reuse(rx);
    }
  esl_stopwatch_Stop(w);
  F_time  = w->user - base_time;
  F_speed = (double) N * (double) L * (double) hmm->M * 1e-6 / F_time;

  /* Backward */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(rng, hmm->f, abc->K, L, dsq);
      h4_reference_Backward (dsq, L, hmm, mo, rx, &sc);
      h4_refmx_Reuse(rx);
    }
  esl_stopwatch_Stop(w);
  B_time  = w->user - base_time;
  B_speed = (double) N * (double) L * (double) hmm->M * 1e-6 / B_time;

  /* Viterbi */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(rng, hmm->f, abc->K, L, dsq);
      h4_reference_Viterbi (dsq, L, hmm, mo, rx, pi, &sc);
      h4_refmx_Reuse(rx);
      h4_path_Reuse(pi);
    }
  esl_stopwatch_Stop(w);
  V_time  = w->user - base_time;
  V_speed = (double) N * (double) L * (double) hmm->M * 1e-6 / V_time;

  printf("# %s (M= %d)\n", "(name here)", hmm->M);
  printf("# Reference Forward:  %8.1f Mc/s\n", F_speed);
  printf("# Reference Backward: %8.1f Mc/s\n", B_speed);
  printf("# Reference Viterbi:  %8.1f Mc/s\n", V_speed);

  free(dsq);
  h4_refmx_Destroy(rx);
  h4_path_Destroy(pi);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*h4REFERENCE_DP_BENCHMARK*/
/*----------------- end, benchmark ------------------------------*/





/*****************************************************************
 * 7. "Brute" unit test
 *****************************************************************/
#ifdef h4REFERENCE_DP_TESTDRIVE
#include "esl_sq.h"
#include "esl_stack.h"

#include "h4_path.h"
#include "emit.h"
#include "logsum.h"
#include "modelsample.h"
#include "standardize.h"

#include "reference_dp.h"

/* The "brute" test, elaborate and savage, enumerates all possible
 * paths for a comparison of a profile of length M and sequence of
 * length L and calculates the Viterbi and Forward/Backward score by
 * brute force maximization and summation over all paths.
 * 
 * Obviously this is only feasible for small M and L, because the
 * number of paths explodes, but even small M and L suffice to
 * exercise every possible transition. 
 *  
 *  M   L     npaths
 *  -   -     ------
 *  1   1          2
 *  2   2         54
 *  3   3       2239  
 *  3   4      30587  *default test
 *  4   4     153306   
 *  4   5    3128767
 *   
 * Larger than that gets crazy. 
 * 
 * We can't expect exact agreement, because of roundoff error
 * (especially in the approximate table-driven logsumexp), and I
 * haven't tried to analyze the expected error (which seems
 * complicated; I have analyzed the expected error of one logsumexp
 * operation, but each path has its own number of them). But just
 * checking for "close" agreement is an excellent bug smoker.
 * 
 * The H4 brute test is stronger than the H3 version. The H3 version
 * enumerated all paths for just one specific small example HMM and
 * sequence.  The H4 version can systematically enumerate all paths
 * for *any* profile/sequence comparison, using a tree traversal
 * algorithm in enumerate_paths(), allowing it to randomly sample
 * profiles and sequences for testing.
 * 
 * The test requires a bit of setup (support for a partial path
 * structure wrapped around H4_PATH, and enumerate_paths() itself),
 * so I broke it out into its own long section in the code.
 */


/* struct partial_path_s
 * 
 * enumerate_paths() is going to use a depth-first tree traversal to
 * enumerate all possible paths for a profile/seq comparison.  This
 * structure, wrapped around a H4_PATH, lets us keep track of an
 * incomplete path as a node in that tree traversal, and when it's
 * been completed.
 */
struct partial_path_s {
  H4_PATH *pi;
  int      k;      // model position of last added pi->st, for MDI; else 0.
  int      L;      // current length of seq this partial path generates

  int      modM;   // length of test profile 
  int      seqL;   // length of test sequence
};

/* declare the support functions for the structure.  enumerate_paths()
 * and the utest itself also lurk below, but there's no need to
 * declare them early.
 */
static struct partial_path_s *create_partial_path(int M, int L);
static int  extended_partial_path(const struct partial_path_s *ppi, int8_t st, struct partial_path_s **ret_np);
static void destroy_partial_path(struct partial_path_s *ppi);


/* create_partial_path()
 * 
 * Create and initialize a new partial path for a brute test of a
 * profile of length <M> against a sequence of length <L>.
 * 
 * This creates the root of a tree of possible paths that we'll
 * enumerate by depth first traversal.
 * 
 * Returns the new partial path structure, or NULL on allocation
 * failure.
 */
static struct partial_path_s *
create_partial_path(int M, int L)
{
  struct partial_path_s *ppi = NULL;
  int status;
  ESL_ALLOC(ppi, sizeof(struct partial_path_s));

  if (( ppi->pi = h4_path_Create() ) == NULL) goto ERROR;  // path starts with Z=0 and nothing in it; h4P_S state is implicit in H4_PATH, not explicit, so check pi->Z == 0 instead 
  ppi->k    = 0;
  ppi->L    = 0;
  ppi->modM = M;
  ppi->seqL = L;
  return ppi;

 ERROR:
  destroy_partial_path(ppi);
  return NULL;
}


/* extended_partial_path()
 * 
 * Try to extend path <ppi> by one state <st>.
 *
 * If successful, return the (newly allocated) extended path in
 * <*ret_np>, and return <eslOK>. (That is, create an entirely new
 * partial path structure, because we'll also be extending <ppi> in
 * other ways too.)
 *
 * However, if:
 *   - adding <st> would exceed test sequence length, or
 *   - adding <st> would cause k position to exceed test profile length, or
 *   - <st> is DG_m and would complete an all-delete "mute" path thru the profile
 * then don't. Return <eslFAIL>, and <*ret_np> is NULL.
 */
static int
extended_partial_path(const struct partial_path_s *ppi, int8_t st, struct partial_path_s **ret_np)
{
  struct partial_path_s *np  = NULL;
  int     k   = ppi->k;       
  int     L   = ppi->L;
  int8_t  prv = ppi->pi->st[ppi->pi->Z-1];
  int     status;

  /* Figure out what adding <st> would do to curr k, L.
   */
  switch (st) {
  case h4P_N:  case h4P_J: case h4P_C: k = 0; if (st == prv) L++; break; // st==prv: i.e. NN, CC, JJ emission-on-transition
  case h4P_G:  case h4P_L:             k = 0;                     break;
  case h4P_MG: case h4P_ML:            k++;   L++;                break;
  case h4P_IG: case h4P_IL:                   L++;                break;
  case h4P_DG: case h4P_DL:            k++;                       break;
  default: esl_fatal("inconceivable");
  }

  /* If k or L ran off end of test profile or test seq, return NULL */
  if (k > ppi->modM || L > ppi->seqL)
    { status = eslFAIL; goto ERROR; }

  /* If we reached DG_m, check for mute cycle; reject such a path */
  if (k == ppi->modM && st == h4P_DG && (prv == h4P_G || (prv == h4P_DG && ppi->pi->rle[ppi->pi->Z-1] == ppi->modM-1)))
    { status = eslFAIL; goto ERROR; }

  /* Otherwise extend the path by <st> */
  ESL_ALLOC(np, sizeof(struct partial_path_s));
  if ((np->pi = h4_path_Clone(ppi->pi) ) == NULL) { status = eslEMEM; goto ERROR; }
  np->k    = k;
  np->L    = L;
  np->modM = ppi->modM;
  np->seqL = ppi->seqL;
  if (( status = h4_path_Append(np->pi, st)) != eslOK) goto ERROR; 

  *ret_np = np;
  return eslOK;

 ERROR:
  *ret_np = NULL;
  return status;
}
  
/* destroy_partial_path()
 * Frees a partial path.
 */
static void
destroy_partial_path(struct partial_path_s *ppi)
{
  if (ppi)
    {
      h4_path_Destroy(ppi->pi);
      free(ppi);
    }
}


/* enumerate_paths()
 * 
 * Enumerate the set of all possible state paths for a test profile of
 * length <M> aligned to a test profile of length <L>. Return the
 * array of paths in <***ret_arrpi>, and their number in <*ret_n>.
 *
 * The number of paths explodes quickly, so <M> and <L> need to be
 * small.
 */
static int
enumerate_paths(int M, int L, H4_PATH ***ret_arrpi, int *ret_n)
{
  struct partial_path_s *ppi = NULL;    // parent partial path that we'll extend
  struct partial_path_s *np  = NULL;    // child partial path that we've extended
  ESL_STACK *pda    = NULL;             // pushdown stack of partial paths, enabling depth first tree traversal
  H4_PATH  **arrpi  = NULL;             // growing array of complete paths, [0..npath-1]
  int        npath  = 0;                // number of paths in <arrpi>
  int        nalloc = 256;              // current allocation of <arrpi>
  int8_t     st;                        // a state type
  int        k;                         // a profile position index for MID state
  int        status;

  ESL_ALLOC(arrpi, sizeof(H4_PATH *) * nalloc);
  if ((    ppi = create_partial_path(M, L))          == NULL) { status = eslEMEM; goto ERROR; }  // S isn't explicit in paths; path starts with N instead
  if (( status = h4_path_Append(ppi->pi, h4P_N))     != eslOK) goto ERROR;
  if ((    pda = esl_stack_PCreate())                == NULL) { status = eslEMEM; goto ERROR; }  
  if (( status = esl_stack_PPush(pda, (void *) ppi)) != eslOK) goto ERROR;

  while ((status = esl_stack_PPop(pda, (void **) &ppi)) == eslOK)
    {
      //printf("Popped. Working on:\n");
      //h4_path_Dump(stdout, ppi->pi);

      st = ppi->pi->st[ppi->pi->Z-1];   // pull this out to a tmp var to clarify a bit when we extend NN, JJ below
      switch (st) {
	
      case h4P_N:
      case h4P_J:
	if ( extended_partial_path(ppi, st,    &np) == eslOK)  esl_stack_PPush(pda, np); 
	if ( extended_partial_path(ppi, h4P_G, &np) == eslOK)  esl_stack_PPush(pda, np);  // H4_PATH doesn't store B explicitly; go to G|L next
	if ( extended_partial_path(ppi, h4P_L, &np) == eslOK)  esl_stack_PPush(pda, np);
	break;

      case h4P_G:
	if ( extended_partial_path(ppi, h4P_MG, &np) == eslOK) esl_stack_PPush(pda, np);
	if ( extended_partial_path(ppi, h4P_DG, &np) == eslOK) esl_stack_PPush(pda, np);
	break;

      case h4P_MG:
      case h4P_IG:
      case h4P_DG:
	if (ppi->k == ppi->modM) // this catches {MD}{LG}_m states, where we can only go to E next
	  {
	    if ( extended_partial_path(ppi, h4P_C, &np) == eslOK) esl_stack_PPush(pda, np);  // E isn't explicit in path; go directly to C|J
	    if ( extended_partial_path(ppi, h4P_J, &np) == eslOK) esl_stack_PPush(pda, np);
	  }
	else
	  {
	    if ( extended_partial_path(ppi, h4P_MG, &np) == eslOK) esl_stack_PPush(pda, np);
	    if ( extended_partial_path(ppi, h4P_IG, &np) == eslOK) esl_stack_PPush(pda, np);
	    if ( extended_partial_path(ppi, h4P_DG, &np) == eslOK) esl_stack_PPush(pda, np);
	  }
	break;

      case h4P_L:
	for (k = 1; k <= M; k++)
	  if ( extended_partial_path(ppi, h4P_ML, &np) == eslOK)
	    {
	      np->pi->rle[np->pi->Z-2] = k;  // set rle[] field for the L state to k of LMk entry
	      np->k = k;                     // also override what extended_partial_path() set k to (which was 1)
	      esl_stack_PPush(pda, np);
	    }
	break;

      case h4P_ML:
      case h4P_DL:
	if ( extended_partial_path(ppi, h4P_J, &np) == eslOK)  esl_stack_PPush(pda, np);
	if ( extended_partial_path(ppi, h4P_C, &np) == eslOK)  esl_stack_PPush(pda, np);
      case h4P_IL: // flowthru is deliberate. IL doesn't transit to E.
	if ( extended_partial_path(ppi, h4P_ML, &np) == eslOK) esl_stack_PPush(pda, np);
	if ( extended_partial_path(ppi, h4P_IL, &np) == eslOK) esl_stack_PPush(pda, np);
	if ( extended_partial_path(ppi, h4P_DL, &np) == eslOK) esl_stack_PPush(pda, np);
	break;
	
      case h4P_C:
	if ( extended_partial_path(ppi, h4P_C, &np) == eslOK)  esl_stack_PPush(pda, np);
	// we've reached the T state (which isn't explicit in the
	// path), so <ppi> is now a complete path. If it's the right
	// length, store it. It's safe to pull it right out of <ppi>,
	// leaving <ppi> as a shell structure that we will soon free.
	if (ppi->L == L)
	  {
	    if (npath == nalloc) {
	      ESL_REALLOC(arrpi, sizeof(H4_PATH *) * nalloc * 2);
	      nalloc *= 2;
	    }
	    arrpi[npath] = ppi->pi;
	    ppi->pi = NULL;
	    npath++;
	  }
	break;

      default: esl_fatal("inconceivable");
      } // end of the big switch statement
      
      destroy_partial_path(ppi);
    } // end of the Pop(), the main loop of the tree traversal.
  if (status != eslEOD) goto ERROR; // i.e. the _Pop() actually failed, instead of stack just being empty

  esl_stack_Destroy(pda);
  *ret_arrpi = arrpi;
  *ret_n     = npath;
  return eslOK;

 ERROR:
  if (arrpi) {
    int i;
    for (i = 0; i < npath; i++) h4_path_Destroy(arrpi[i]);
    free(arrpi);
  }
  esl_stack_Destroy(pda);
  *ret_arrpi = NULL;
  *ret_n     = 0;
  return status;
}


/* utest_brute()
 * 
 * The brute force unit test itself.
 * 
 * Characterizing error:
 *     ./reference_dp_utest --diag brute --bruteN 10000 > foo
 *     avg -f3   # viterbi absolute error
 *     avg -f7   # fwd abs error
 *
 * viterbi absolute error:
 *    default:       -5e-10 +- 1e-7     10sd => 1e-6  but that's close to range so => 1e-5
 *    exact logsum:   2e-9  +- 1e-7      ""
 * fwd absolute error:
 *    default:       2.7e-6 +- 0.0003   10sd => 0.003  
 *    exact logsum:  1.7e-7 +- 3e-7     10sd => 3e-6  but close to range so => 1e-5
 */
static void
utest_brute(FILE *diagfp, ESL_RANDOMNESS *rng, int M, int L, int ntrials)
{
  char            msg[]    = "reference_fwdback brute unit test failed";
  ESL_ALPHABET   *abc      = esl_alphabet_Create(eslAMINO);
  ESL_SQ         *sq       = esl_sq_CreateDigital(abc);
  H4_PROFILE     *hmm      = NULL;
  H4_MODE        *mo       = h4_mode_Create();
  H4_PATH        *true_pi  = h4_path_Create();
  H4_PATH        *vpi      = h4_path_Create();
  H4_PATH       **arrpi    = NULL;
  H4_REFMX       *rmx      = h4_refmx_Create(M, L); 
  float          *pathsc   = NULL;
  float           vtol_abs = 1e-5;
  float           ftol_abs = (h4_logsum_IsSlowExact() ? 1e-5 : 0.003);
  int   npath;
  float vsc, fsc, bsc, brute_vsc, brute_fsc;
  int   i;
  int   status;

  /* path enumeration is independent of model, params, seq */
  if ( enumerate_paths(M, L, &arrpi, &npath) != eslOK) esl_fatal(msg);
  ESL_ALLOC(pathsc, sizeof(float) * npath);

  /* all paths should validate */
  for (i = 0; i < npath; i++)
    if ( h4_path_Validate(arrpi[i], M, L, NULL) != eslOK) esl_fatal(msg);

  if (diagfp) {
    fprintf(diagfp, "# npath    = %d\n", npath);
    fprintf(diagfp, "# RNG seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));
    esl_dataheader(diagfp,
		   12, "bru_vsc", 12, "vsc", 12, "vit abserr", 12, "vit relerr",
		   12, "bru_fsc", 12, "fsc", 12, "fwd abserr", 12, "fwd relerr", 0);
  }

  while (ntrials--)
    {
      /* sample a random model of length M; half the time, zeropepper it for extra pathology */
      if (esl_rnd_Roll(rng, 2) == 0) h4_modelsample_zeropeppered(rng, abc, M, &hmm);
      else                           h4_modelsample(rng, abc, M, &hmm);

      /* configure the model and alignment mode; sample a length L sequence from model */
      if ( h4_mode_SetLength(mo, L) != eslOK) esl_fatal(msg);
      do {
	if (h4_emit(rng, hmm, mo, sq, true_pi) != eslOK) esl_fatal(msg);
      } while (sq->n != L);

      //h4_profile_Dump(stdout, hmm);
      //h4_mode_Dump(stdout, mo);

      /* run Forward, Backward, Viterbi */
      if (h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rmx,       &fsc) != eslOK) {esl_fatal(msg);}   h4_refmx_Reuse(rmx);  // {..} are to shut compiler up about misleading indentation
      if (h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rmx,       &bsc) != eslOK) {esl_fatal(msg);}   h4_refmx_Reuse(rmx);
      if (h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rmx,  vpi, &vsc) != eslOK) {esl_fatal(msg);}

      /* Backward and Forward scores should match */
      if (fabs(fsc - bsc) > ftol_abs) esl_fatal(msg);

      /* Score all brute paths with this seq/profile  */
      for (i = 0; i < npath; i++)
	if ( h4_path_Score(arrpi[i], sq->dsq, hmm, mo, &(pathsc[i])) != eslOK) esl_fatal(msg);

      //for (i = 0; i < npath; i++) {
      //  printf("%d : %9.5f ", i, pathsc[i]);
      //  h4_path_DumpCigar(stdout, arrpi[i]);
      //}

      /* Get brute force vsc, fsc by maximization, summation over paths */
      esl_vec_FSortIncreasing(pathsc, npath);        // After this sort, we've messed up the indexing of individual paths
      brute_vsc = pathsc[npath-1];
      brute_fsc = esl_vec_FLog2Sum(pathsc, npath);   // Sorting scores from small to large increases numerical accuracy 

      if (diagfp)
	fprintf(diagfp, "%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n",
		brute_vsc, vsc, brute_vsc - vsc, fabs((brute_vsc - vsc) / vsc), 
		brute_fsc, fsc, brute_fsc - fsc, fabs((brute_fsc - fsc) / fsc));
      else
	{
	  /* The main test: reference DP calculations are acceptably close to brute force calculations */
	  if (fabs(brute_vsc - vsc) > vtol_abs) esl_fatal(msg);
	  if (fabs(brute_fsc - fsc) > ftol_abs) esl_fatal(msg);
	}

      /* Deliberately not testing that the max-scoring brute path matches the Viterbi path,
       * because ties and/or roundoff error can select a different path with a near-identical 
       * score.
       */
      h4_path_Reuse(vpi);
      h4_path_Reuse(true_pi);
      h4_profile_Destroy(hmm);
      esl_sq_Reuse(sq);
      h4_refmx_Reuse(rmx);
    }

  free(pathsc);
  h4_refmx_Destroy(rmx);
  for (i = 0; i < npath; i++) h4_path_Destroy(arrpi[i]);
  free(arrpi);
  h4_path_Destroy(vpi);
  h4_path_Destroy(true_pi);
  h4_mode_Destroy(mo);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  return;

 ERROR: 
  esl_fatal("allocation failed");
}
#endif // h4REFERENCE_DP_TESTDRIVE
/*----------------- end, brute unit test ------------------------*/


/*****************************************************************
 * 8. Other unit tests
 *****************************************************************/
#ifdef h4REFERENCE_DP_TESTDRIVE

/* utest_generation()
 * 
 * The "generation" test samples (emits) a sequence from a randomly
 * generated profile HMM, then compares them by DP. It tests that:
 *    1. Vit score >= score of actual emitted path
 *    2. Fwd, Bck scores >= Vit score
 *    3. Fwd = Bck score
 *    4. bitscores have positive expectation (because they're generated positives)
 *    5. all DP matrices and paths Validate() (obey documented conventions)  
 *    6. (overwrite:) Decoding gets exactly the same result whether
 *       overwriting Backwards matrix or using dedicated Decoding matrix
 *    7. (rowsum:) in a decoding matrix, on rows 1..L, sum over emitting
 *       states (ML|MG|IL|IG|N|CC|JJ) is 1.0
 *    8. (colsum:) for a uniglocal model (only), in a decoding matrix,
 *       each consensus column k=1..M is visited once and only once in
 *       a Mk or Dk: \sum_i Mk+Dk = 1.0 for each k.  \sum_i = 1 for
 *       E,B,G specials as well.  (Mute partial cycle issue isn't in
 *       play because model is uniglocal for this test.)
 *
 * The colsum test only gets done if <do_colsum> is TRUE, since it
 * only works for a uniglocal model. Other tests work for any model or
 * length configuration (and when <do_colsum> is FALSE, the
 * configuration is randomized.)
 *
 * Can fail stochastically, in theory, but tolerances are set to ~10
 * sd's to make that unlikely. Nonetheless, best to run with fixed RNG
 * seed in production code.
 * 
 * To collect data on error magnitudes:
 * (--diag generation is similar but doesn't give colsum stats)
 *     ./reference_dp_utest --diag gencol     -N 10000 > foo   
 *
 *               ------- default ---------       ----- exact logsum -----
 *                mean      sd      ~10sd      mean     sd        ~10sd
 *  fsc-bsc     -3.8e-6   0.0005    0.005        -7e-7   8e-6     0.0001
 *  row_maxerr   0.0004   0.0002    0.002        3.8e-6  4e-6     0.00004
 *  col_maxerr   0.0004   0.0002    0.002        2.4e-6  3e-6     0.00003
 *
 *  Chose abs tol of 0.005, 0.0001 based on this.
 */
static void
utest_generation(FILE *diagfp, ESL_RANDOMNESS *rng, int alphatype, int M, int L, int nseq, int do_colsum)
{
  char            msg[]    = "reference dp generation unit test failed";
  ESL_ALPHABET   *abc      = esl_alphabet_Create(alphatype);
  ESL_SQ         *sq       = esl_sq_CreateDigital(abc);
  H4_PROFILE     *hmm      = NULL;
  H4_MODE        *mo       = h4_mode_Create();
  H4_REFMX       *rxv      = h4_refmx_Create(M, L);
  H4_REFMX       *rxf      = h4_refmx_Create(M, L);
  H4_REFMX       *rxb      = h4_refmx_Create(M, L);
  H4_REFMX       *rxd      = h4_refmx_Create(M, L);
  H4_PATH        *epi      = h4_path_Create();
  H4_PATH        *vpi      = h4_path_Create();
  float           a_tol    = ( h4_logsum_IsSlowExact() ? 0.0001 : 0.005);
  float           avgsc    = 0.0f;
  float          *colsum   = malloc(sizeof(float) * (M+1));
  float           xcolsum[h4R_NXCELLS];
  float           rowsum;
  float           max_colerr, max_rowerr;
  float           psc, vsc, fsc, bsc;
  float          *dpc;
  int             i,k,y,idx;
  char            errbuf[eslERRBUFSIZE];

  if (diagfp) esl_dataheader(diagfp, 12, "vsc-psc", 12, "fsc-vsc", 12, "fsc-bsc", 12, "max_rowerr", 12, "max_colerr", 12, "bitscore", 0);

  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);

  if (! do_colsum) // choose a random configuration:
    {
      switch (esl_rnd_Roll(rng, 6)) {
      case 0: h4_mode_SetDefault(mo);   break;
      case 1: h4_mode_SetLocal(mo);     break;
      case 2: h4_mode_SetGlocal(mo);    break;
      case 3: h4_mode_SetUnihit(mo);    break;
      case 4: h4_mode_SetUnilocal(mo);  break;
      case 5: h4_mode_SetUniglocal(mo); break;
      }
    }
  else h4_mode_SetUniglocal(mo);

  for (idx = 0; idx < nseq; idx++)
    {
      if ( h4_mode_SetLength(mo, L)                   != eslOK) esl_fatal(msg);                              // set length model to L for emission
      do { if ( h4_emit(rng, hmm, mo, sq, epi)        != eslOK) esl_fatal(msg); } while (sq->n > (M+L)*3);   // keep seq length reasonable
      if ( h4_mode_SetLength(mo, sq->n)               != eslOK) esl_fatal(msg);   // reset length model to this sampled seq's length before scoring or DP 
      if ( h4_path_Score(epi, sq->dsq, hmm, mo, &psc) != eslOK) esl_fatal(msg);   // emitted path score. Viterbi can improve on this.

      if ( h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc)  != eslOK) esl_fatal(msg);
      if ( h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxb,      &bsc)  != eslOK) esl_fatal(msg);
      if ( h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxv, vpi, &vsc)  != eslOK) esl_fatal(msg);
      if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxb, rxd)   != eslOK) esl_fatal(msg);
      avgsc += fsc - mo->nullsc; // avg score test (positive expectation, for generated seqs);

      if ( h4_refmx_Validate(rxf, errbuf)          != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
      if ( h4_refmx_Validate(rxb, errbuf)          != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
      if ( h4_refmx_Validate(rxv, errbuf)          != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
      if ( h4_refmx_Validate(rxd, errbuf)          != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
      if ( h4_path_Validate(vpi, M, sq->n, errbuf) != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
      if ( h4_path_Validate(epi, M, sq->n, errbuf) != eslOK) esl_fatal("%s\n  %s", msg, errbuf);

      /* overwrite test */
      if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxb, rxb) != eslOK) esl_fatal(msg);
      if ( h4_refmx_Validate(rxb, errbuf)                                != eslOK) esl_fatal("%s\n  %s", msg, errbuf);
      if ( h4_refmx_CompareDecoding(rxd, rxb, 0.)                        != eslOK) esl_fatal(msg);  // 0. = exact equality expected.

      /* colsum test */
      if (do_colsum)
	{
	  esl_vec_FSet(colsum,  M+1,         0.);
	  esl_vec_FSet(xcolsum, h4R_NXCELLS, 0.);
	  for (i = 0; i <= sq->n; i++)
	    {
	      for (k = 0, dpc = rxd->dp[i]; k <= M; k++, dpc += h4R_NSCELLS)   // colsum[0] will be -inf, but we aren't going to check there.
		colsum[k] += dpc[h4R_ML] + dpc[h4R_MG] + dpc[h4R_DL] + dpc[h4R_DG];
	      for (y = 0; y < h4R_NXCELLS; y++)
		xcolsum[y] += dpc[y];
	    }
	  max_colerr = 0.;
	  for (k = 1; k <= M; k++) max_colerr = ESL_MAX(max_colerr, fabs(colsum[k]  - 1.0));
	  max_colerr = ESL_MAX(max_colerr, fabs(xcolsum[h4R_E] - 1.0));
	  max_colerr = ESL_MAX(max_colerr, fabs(xcolsum[h4R_B] - 1.0));
	  max_colerr = ESL_MAX(max_colerr, fabs(xcolsum[h4R_G] - 1.0));
	}

      /* rowsum test */
      max_rowerr = 0.;
      for (i = 1; i <= sq->n; i++)
	{
	  for (rowsum = 0., k = 1, dpc = rxd->dp[i] + h4R_NSCELLS; k <= M; k++, dpc += h4R_NSCELLS)
	    rowsum += dpc[h4R_ML] + dpc[h4R_MG] + dpc[h4R_IL] + dpc[h4R_IG];
	  rowsum += dpc[h4R_N] + dpc[h4R_JJ] + dpc[h4R_CC];
	  max_rowerr = ESL_MAX(max_rowerr, fabs(rowsum - 1.0));
	}

      if (!diagfp)
	{ 
	  if (! isfinite(vsc) || ! isfinite(fsc) || ! isfinite(bsc))        esl_fatal(msg);
	  if (psc > vsc || vsc > fsc)                                       esl_fatal(msg); // Viterbi >= emitted path, and fwd >= vit
	  if (esl_FCompareNew(fsc, bsc, 0., a_tol) != eslOK)                esl_fatal(msg); // fwd = bck 
	  if (do_colsum && ( ! isfinite(max_colerr) || max_colerr > a_tol)) esl_fatal(msg); // colsum test
	  if (max_rowerr > a_tol)                                           esl_fatal(msg); // rowsum test
	}
      else // non-NULL diagfp dumps data we can use to estimate chance failure rates
	fprintf(diagfp, "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", vsc-psc, fsc-vsc, fsc-bsc, max_rowerr, max_colerr, fsc-mo->nullsc);

      h4_refmx_Reuse(rxd); h4_refmx_Reuse(rxb); h4_refmx_Reuse(rxf); h4_refmx_Reuse(rxv);
      h4_path_Reuse(vpi);  h4_path_Reuse(epi);
    }
  if (avgsc < 0.0f) esl_fatal(msg);  // positive expectation

  free(colsum);
  h4_path_Destroy(epi);
  h4_path_Destroy(vpi);
  h4_refmx_Destroy(rxd);
  h4_refmx_Destroy(rxb);
  h4_refmx_Destroy(rxf);
  h4_refmx_Destroy(rxv);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_sq_Destroy(sq);
}



/* utest_duality()
 * 
 * The "duality" unit test uses the fact that for unihit models, the
 * dual-mode score should be equal to logsum(local_score +
 * glocal_score) - 1 bit.
 *                 
 * Unihit, because in dual multihit mode there's paths with both
 * glocal and local subpaths, that glocal-only and local-only can't
 * see.
 * 
 * Characterizing error magnitude:
 *     ./reference_dp_utest --diag duality -N 10000 > foo
 *     avg -f5 foo
 *  
 *      default:  -5e-7 +- 0.0001   10sd => 0.001    but range is -0.001..0.001 so => 0.002
 * exact logsum:  -3e-7 +- 1e-6          => 0.00001  but range is up to 6e-6 so    => 0.00002
 */
static void
utest_duality(FILE *diagfp, ESL_RANDOMNESS *rng, int alphatype, int M, int L, int nseq)
{
  char          msg[]  = "reference_dp duality unit test failed";
  ESL_DSQ      *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_ALPHABET *abc    = esl_alphabet_Create(alphatype);
  H4_PROFILE   *hmm    = NULL;
  H4_MODE      *mo     = h4_mode_Create();
  H4_REFMX     *rxf    = h4_refmx_Create(M,L);
  float         dual_sc, local_sc, glocal_sc, combined_sc;
  float         tol    = ( h4_logsum_IsSlowExact() ? 0.00002 : 0.002 );
  int           idx;

  if (diagfp)
    esl_dataheader(diagfp, 12, "dual_sc", 12, "local_sc", 12, "glocal_sc", 12, "combined_sc",
		   20, "dual-combined", 0);

  if ( h4_modelsample(rng, abc, M, &hmm)      != eslOK) esl_fatal(msg);

  for (idx = 0; idx < nseq; idx++)
    {
      if ( esl_rsq_xfIID(rng, hmm->f, abc->K, L, dsq) != eslOK) esl_fatal(msg);
      
      if ( h4_mode_SetUnihit(mo)    != eslOK) esl_fatal(msg);
      if ( h4_mode_SetLength(mo, L) != eslOK) esl_fatal(msg);
      if ( h4_reference_Forward (dsq, L, hmm, mo, rxf, &dual_sc) != eslOK) esl_fatal(msg);
      if ( h4_refmx_Reuse(rxf)      != eslOK) esl_fatal(msg);

      if ( h4_mode_SetUnilocal(mo)  != eslOK) esl_fatal(msg);
      if ( h4_mode_SetLength(mo, L) != eslOK) esl_fatal(msg);
      if ( h4_reference_Forward (dsq, L, hmm, mo, rxf, &local_sc) != eslOK) esl_fatal(msg);
      if ( h4_refmx_Reuse(rxf)      != eslOK) esl_fatal(msg);

      if ( h4_mode_SetUniglocal(mo) != eslOK) esl_fatal(msg);
      if ( h4_mode_SetLength(mo, L) != eslOK) esl_fatal(msg);
      if ( h4_reference_Forward (dsq, L, hmm, mo, rxf, &glocal_sc) != eslOK) esl_fatal(msg);
      if ( h4_refmx_Reuse(rxf)      != eslOK) esl_fatal(msg);

      combined_sc = h4_logsum(local_sc, glocal_sc) - 1.;

      if (diagfp)
	fprintf(diagfp, "%12.6f %12.6f %12.6f %12.6f %20g\n",
		dual_sc, local_sc, glocal_sc, combined_sc, dual_sc - combined_sc);
      else
	{ if (esl_FCompareAbs(dual_sc, combined_sc, tol) != eslOK) esl_fatal(msg);  }
    }
  
  h4_refmx_Destroy(rxf);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  free(dsq);
}


/* utest_enumeration()
 * 
 * The "enumeration" unit test makes sure that Forward scores sum to
 * one over all possible sequences, by creating a restricted profile
 * (called an "enumerable" profile) that only generates sequences of
 * length 1..2M-1. 
 * 
 * The enumerable profile has no nonzero II transitions, and the
 * alignment mode is unihit with length L=0.  The mode is still
 * dual-mode glocal/local, so all transitions except NN|CC|JJ|II will
 * be tested.
 * 
 * All P(seq | profile) terms must be >> DBL_EPSILON, so M needs to be
 * very small (<10). Small M also helps this expensive test run in
 * reasonable time. 
 * 
 * Using a tiny alphabet helps limit the size of the enumerated
 * sequence space, so we use eslCOINS.
 * 
 * Characterizing error magnitude:
 *     ./reference_dp_utest --diag enumeration -N 10000 > foo
 *     avg foo
 * 
 *      default: -2e-7  +- 5e-5   10sd => 0.001
 * exact logsum: -8e-10 +- 9e8         => 1e-6 but range gets close so => 1e-5
 */
static void
utest_enumeration(FILE *diagfp, ESL_RANDOMNESS *rng, int M)
{
  char          msg[] = "reference_dp enumeration unit test failed";
  ESL_ALPHABET *abc   = esl_alphabet_Create(eslCOINS); // using the eslCOINS alphabet helps keep the enumeration's size down.
  H4_PROFILE   *hmm   = NULL;
  H4_MODE      *mo    = h4_mode_Create();
  int           maxL  = 2*M-1;	
  ESL_DSQ      *dsq   = malloc(sizeof(ESL_DSQ) * (maxL+3));  // +3 = +2 for sentinels; +1 for maxL+1 test seq with P=0.
  H4_REFMX     *rxf   = h4_refmx_Create(M, maxL+1);          // +1 for maxL+1 test seq
  int           i, L;
  float         fsc, mutesc;
  float         bg_ll;
  double        fp;
  double        total_p = 0.0;
  double        tol     = ( h4_logsum_IsSlowExact() ? 0.00001 : 0.001 );

  if ( h4_modelsample_Enumerable2(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  if ( h4_mode_SetUnihit(mo)                         != eslOK) esl_fatal(msg);
  if ( h4_mode_SetLength(mo, 0)                      != eslOK) esl_fatal(msg);

  /* By design, the profile neglects a small amount of probability
   * mass for the L=0 case (a mute path through G->D1...Dm->E) because
   * it removes the mute path using wing retraction. To test that the
   * total sum is 1.0, we have to add that small mass back. 
   */
  if ( h4_profile_MutePathScore(hmm, &mutesc)        != eslOK) esl_fatal(msg);
  total_p = exp2(mutesc) * mo->pglocal;
  //printf("mute = %g\n", total_p);

  /* L=0 is included deliberately, to test that L=0 seq gets -inf score as it should */
  for (L = 0; L <= maxL; L++)
    {
      /* initialize dsq[1..L] at "0000..." */
      dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;  // DP routines don't need these sentinels, actually.
      for (i = 1; i <= L; i++) dsq[i] = 0;

      /* enumerate and score all sequences of length L */
      while (1)	
	{
	  if ( h4_reference_Forward(dsq, L, hmm, mo, rxf, &fsc) != eslOK) esl_fatal(msg);
	  
	  /* calculate bg log likelihood component of the scores */
	  for (bg_ll = 0., i = 1; i <= L; i++)  bg_ll += log2f(hmm->f[dsq[i]]);
	  
	  /* convert to probability P(seq|model), adding the bg LL back to the LLR */
	  fp =  exp2(fsc + bg_ll);
	  //printf("fp = %g\n", fp);
	  total_p += fp;

	  /* Increment dsq to next seq, like a reversed odometer; works for any alphabet! */
	  for (i = 1; i <= L; i++) 
	    if (dsq[i] < abc->K-1) { dsq[i]++; break; } else { dsq[i] = 0; }
	  if (i > L) break;	/* we're done enumerating sequences */

	  h4_refmx_Reuse(rxf);
	}
    }

  if (!diagfp)
    {
      if (fabs(total_p - 1.0) > tol) esl_fatal(msg);
    }
  else
    fprintf(diagfp, "%g\n", 1.0 - total_p);

  /* And before we go ... any sequence of length L > maxL should get score -infty */
  if ( esl_rsq_xfIID(rng, hmm->f, abc->K, maxL+1, dsq)        != eslOK) esl_fatal(msg);
  if ( h4_reference_Forward(dsq, maxL+1, hmm, mo, rxf, &fsc)  != eslOK) esl_fatal(msg);
  if ( fsc != -eslINFINITY) esl_fatal(msg);

  h4_refmx_Destroy(rxf);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  free(dsq);
}


/* utest_singlepath()
 * 
 * The "singlepath" unit test creates a special profile with only a single
 * possible parse (including its emitted residues): it samples an HMM
 * whose transitions and emissions can only be 1.0 or 0.0, and it
 * creates a uniglocal L=0 profile from it. Then it samples that 
 * path (and sequence), and runs Forward, Backward, and Viterbi on it.
 * It tests:
 *    1. Generated trace is the same as the Viterbi trace.
 *    2. Generated trace score is the same as Viterbi score.
 *    3. Forward score = Backward score.
 *    4. Forward score = score of generated trace
 *
 * Characterizing error:
 *    ./reference_dp_utest --diag singlepath -N 10000 > foo    
 *    avg -f1   # psc-vsc. should be 0   (not sure it's guaranteed tho)
 *    avg -f2   # fsc-bsc. 
 *    avg -f3   # fsc-psc. should be 0   (ditto)
 *
 * logsum has no effect in a singlepath test, so no difference is
 * expected in exact logsum vs. default.
 * 
 * h4_path_Score() is written carefully to sum in exactly the
 * same order of evaluation as Vit/Fwd, so we can expect exact
 * equality - but not sure it is absolutely guaranteed.
 * 
 *  psc-vsc  default          0
 *           exact logsum     0
 *  fsc-bsc  default          2e-7 +- 9e-6   10sd => 0.0001
 *           exact logsum    -4e-8 +- 9e-6        => 0.0001
 *  fsc-psc  default          0
 *           exact logsum     0             
 */
static void
utest_singlepath(FILE *diagfp, ESL_RANDOMNESS *rng, int alphatype, int M)
{
  char          msg[] = "reference_dp singlepath unit test failed";
  ESL_ALPHABET *abc   = esl_alphabet_Create(alphatype);
  H4_PROFILE   *hmm   = NULL;
  H4_MODE      *mo    = h4_mode_Create();
  ESL_SQ       *sq    = esl_sq_CreateDigital(abc);
  H4_PATH      *pi0   = h4_path_Create();            // true (generated) path
  H4_PATH      *vpi   = h4_path_Create();            // Viterbi path
  H4_PATH      *pi    = h4_path_Create();
  H4_REFMX     *rxv   = h4_refmx_Create(M, 100);
  H4_REFMX     *rxf   = h4_refmx_Create(M, 100);
  H4_REFMX     *rxd   = h4_refmx_Create(M, 100);
  H4_REFMX     *rxa   = NULL;
  float        *wrk   = NULL;
  float         psc, vsc, fsc, bsc;
  float         a_tol = 0.0001;
  int           N     = 10;
  char          errbuf[eslERRBUFSIZE];
  int           idx;

  /* Create a profile that has only a single possible path (including
   * emissions) thru it; requires configuring in uniglocal mode w/ L=0
   */
  if ( h4_modelsample_SinglePath(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);
  if ( h4_mode_SetUniglocal(mo)                     != eslOK) esl_fatal(msg);
  if ( h4_mode_SetLength(mo, 0)                     != eslOK) esl_fatal(msg);

  /* Sample that sequence and path */
  if ( h4_emit(rng, hmm, mo, sq, pi0)               != eslOK) esl_fatal(msg);
  
  /* Run DP routines, collect scores that should all match  */
  if (h4_path_Score(pi0, sq->dsq, hmm, mo, &psc)                      != eslOK) esl_fatal(msg);
  if (h4_reference_Viterbi (sq->dsq, sq->n, hmm, mo, rxv, vpi, &vsc)  != eslOK) esl_fatal(msg);  
  if (h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,      &fsc)  != eslOK) esl_fatal(msg);
  if (h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxd,      &bsc)  != eslOK) esl_fatal(msg);
  if (h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxd, rxd)   != eslOK) esl_fatal(msg);

  /* Approximate decoding by stochastic traceback.  Since there's only
   * one possible path, the <rxa> matrix will have only 1|0 values.
   */
  if ( (rxa = h4_refmx_Create(M, sq->n))               == NULL)  esl_fatal(msg);
  if ( h4_refmx_SetType  (rxa, M, sq->n, h4R_DECODING) != eslOK) esl_fatal(msg);      
  if ( h4_refmx_SetValues(rxa, 0.0)                    != eslOK) esl_fatal(msg);
  for (idx = 0; idx < N; idx++)
    {
      if ( h4_reference_StochasticTrace(rng, &wrk, hmm, mo, rxf, pi) != eslOK) esl_fatal(msg);
      if ( h4_refmx_CountPath(pi, rxa)                               != eslOK) esl_fatal(msg);
      if ( h4_path_Compare(pi0, pi)                                  != eslOK) esl_fatal(msg); 
      if ( h4_path_Reuse(pi)                                         != eslOK) esl_fatal(msg);
    }
  h4_refmx_Scale(rxa, 1./(float)N);

  if (diagfp)
    fprintf(diagfp, "%20g %20g %20g\n", psc-vsc, fsc-bsc, fsc-psc);

  if (h4_path_Compare(pi0, vpi)                 != eslOK) esl_fatal(msg); 
  if (psc != vsc)                                         esl_fatal(msg);  
  if (psc != fsc)                                         esl_fatal(msg);
  if (esl_FCompareNew(fsc, bsc, 0., a_tol)      != eslOK) esl_fatal(msg);
  if (h4_refmx_CompareDecoding(rxd, rxa, a_tol) != eslOK) esl_fatal(msg);
  if (h4_refmx_Validate(rxv, errbuf)            != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
  if (h4_refmx_Validate(rxf, errbuf)            != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
  if (h4_refmx_Validate(rxd, errbuf)            != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
  if (h4_refmx_Validate(rxa, errbuf)            != eslOK) esl_fatal("%s:\n%s", msg, errbuf);

  free(wrk);
  h4_refmx_Destroy(rxa);
  h4_refmx_Destroy(rxf);
  h4_refmx_Destroy(rxd);
  h4_refmx_Destroy(rxv);
  esl_sq_Destroy(sq);
  h4_path_Destroy(pi);
  h4_path_Destroy(pi0);
  h4_path_Destroy(vpi);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_alphabet_Destroy(abc);
}

/* utest_approx_decoding()
 * 
 * The "approx-decoding" test compares exact posterior decoding to a
 * stochastic approximation (by taking a large sample of stochastic
 * tracebacks). It does this for a randomly sampled profile, put into
 * a random alignment mode, for either an L=L or L=0 length config,
 * compared to one sequence (sampled from the profile, so we're sure
 * that it can align to the profile with nonzero probability).
 * 
 * Only one sequence, because the stochastic tracebacks are expensive.
 * 
 * Tests that the two decoding matrix values are identical (within some
 * absolute numerical error tolerance). 
 * 
 * Characterizing error (only -N100 because this is expensive)
 *    ./reference_dp_utest  --diag decoding -N 100 > foo
 *    avg foo
 *    
 *      default: 0.003 +- 0.001    10sd = 0.01
 * exact logsum: 0.003 +- 0.0008   10sd = 0.01   
 *    
 * default, exact logsum show little difference because here the
 * absolute difference is driven mostly by sampling error, from the
 * finite number of stochastic traces.
 */
static void
utest_approx_decoding(FILE *diagfp, ESL_RANDOMNESS *rng, int alphatype, int M, int L)
{
  char          msg[] = "reference_dp approx_decoding unit test failed";
  ESL_ALPHABET *abc   = esl_alphabet_Create(alphatype);
  H4_PROFILE   *hmm   = NULL;
  H4_MODE      *mo    = h4_mode_Create();
  ESL_SQ       *sq    = esl_sq_CreateDigital(abc);
  H4_REFMX     *rxf   = h4_refmx_Create(M, 100);
  H4_REFMX     *rxd   = h4_refmx_Create(M, 100);  // also used for Backward.
  H4_REFMX     *rxa   = NULL;
  H4_PATH      *pi    = h4_path_Create();
  float        *wrk   = NULL;                     // tmp space needed by stochastic traceback
  int           N     = 100000;
  float         a_tol = 0.01;
  int           idx;
  char          errbuf[eslERRBUFSIZE];

  /* Sample a random profile HMM */
  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);

  /* Select a random alignment mode */
  switch (esl_rnd_Roll(rng, 6)) {
  case 0: h4_mode_SetDefault(mo);   break;
  case 1: h4_mode_SetLocal(mo);     break;
  case 2: h4_mode_SetGlocal(mo);    break;
  case 3: h4_mode_SetUnihit(mo);    break;
  case 4: h4_mode_SetUnilocal(mo);  break;
  case 5: h4_mode_SetUniglocal(mo); break;
  }
 
  /* randomly set mode's length to L or 0 for generating a sequence */
  switch (esl_rnd_Roll(rng, 2)) {
  case 0: h4_mode_SetLength(mo, 0); break;
  case 1: h4_mode_SetLength(mo, L); break;
  }

  /* generate (sample) a sequence from the profile */
  do {   
    esl_sq_Reuse(sq);
    h4_emit(rng, hmm, mo, sq, NULL);
  } while (sq->n > L * 3); // keep sequence length from getting ridiculous; long seqs do have higher abs error per cell 

  /* for DP, length model has to be sq->n or 0 -- original L doesn't work. */
  if (mo->L == L) h4_mode_SetLength(mo, sq->n);

  /* posterior decoding, by DP (exact - we hope) */
  if ( h4_reference_Forward (sq->dsq, sq->n, hmm, mo, rxf,  NULL)     != eslOK) esl_fatal(msg);
  if ( h4_reference_Backward(sq->dsq, sq->n, hmm, mo, rxd,  NULL)     != eslOK) esl_fatal(msg);
  if ( h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, rxf, rxd, rxd)  != eslOK) esl_fatal(msg);

  /* Approximate decoding by stochastic traceback  */
  if ( (rxa = h4_refmx_Create(M, sq->n))               == NULL)  esl_fatal(msg);
  if ( h4_refmx_SetType  (rxa, M, sq->n, h4R_DECODING) != eslOK) esl_fatal(msg);      
  if ( h4_refmx_SetValues(rxa, 0.0)                    != eslOK) esl_fatal(msg);
  for (idx = 0; idx < N; idx++)
    {
      if ( h4_reference_StochasticTrace(rng, &wrk, hmm, mo, rxf, pi) != eslOK) esl_fatal(msg);
      if ( h4_refmx_CountPath(pi, rxa)                               != eslOK) esl_fatal(msg);
      if ( h4_path_Reuse(pi)                                         != eslOK) esl_fatal(msg);
    }
  h4_refmx_Scale(rxa, 1./(float)N);

  if (diagfp)
    {
      int    i,y;
      float *dpe, *dpa;
      float  max_err = 0.0;
      for (i = 0; i <= sq->n; i++)
	for (y = h4R_NSCELLS, dpe = rxd->dp[i]+h4R_NSCELLS, dpa = rxa->dp[i]+h4R_NSCELLS; y < (M+1)*h4R_NSCELLS + h4R_NXCELLS; y++, dpe++, dpa++)
	  max_err = ESL_MAX(max_err, fabs(*dpe - *dpa));
      printf("%g\n", max_err);
    }
  else	  
    {
      /* Tests */
      if ( h4_refmx_CompareDecoding(rxd, rxa, a_tol) != eslOK) esl_fatal(msg);
      if ( h4_refmx_Validate(rxd, errbuf)            != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
      if ( h4_refmx_Validate(rxa, errbuf)            != eslOK) esl_fatal("%s:\n%s", msg, errbuf);
    }
  
  free(wrk);
  h4_path_Destroy(pi);
  h4_refmx_Destroy(rxa);
  h4_refmx_Destroy(rxd);
  h4_refmx_Destroy(rxf);
  esl_sq_Destroy(sq);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
}


/* utest_mute_partial_cycle()
 * 
 * The "mute partial cycle" issue is why posterior decoding values for
 * DGk states are expected counts, not probabilities. See
 * reference_dp.md for details.
 * 
 * This example is contrived to maximize the problem. The test profile
 * (with multihit glocal L=0 mode) allows two different glocal paths
 * for a sequence of length 6:
 *  
 *                ___========___________________________========
 *  S N B G M1 M2 M3 D4 D5 D6 D7 D8 D9 E J B G D1 D2 D3 D4 D5 D6 M7 M8 M9 E C T
 *                ___========___
 *  S N B G M1 M2 M3 D4 D5 D6 M7 M8 M9 E C T
 *  
 *  The overlined regions are on the same DP row (i=3). The double
 *  overlined regions show how the first path uses D4/D5/D6 twice on
 *  that row. After the posterior decoding code adds the contribution
 *  of G->D1..D6->M7 left wing retraction, posterior decoding values
 *  of DG4..DG6 states can give values > 1.
 *  
 *  To maximize the issue, we contrive to maximize the probability of
 *  the first path.  The second path has a strict subset of the
 *  transition/emission probabilities as the first model, so the
 *  second path cannot have less probability than the first.  So we
 *  want to maximize the probability of the transitions that only the
 *  first path uses. Thus:
 *  
 *    B->G   = 1      i.e. glocal-only
 *    G->D1  ~ 1      can't be exactly 1, because we need G->M1 > 0 too
 *    D1->D2 = 1       .. (ditto for D2,D3,D7,D8-> DD transitions)
 *    D6->D7 ~ 1      again can't be exactly 1, because D6->M7 > 0 too
 *    E->J   = large  Let it be the standard multihit 0.5.
 *    J->B   = 1      i.e. use a L=0 length model
 *    
 *  Thus the tEJ=0.5 term for standard multihit alignment sets an
 *  upper bound of 0.5 for the probability of path 1 relative to path
 *  2. Thus the maximum decoded DGk expected count is 1.33 ( 1.0*1 +
 *  0.5*2 / ( 1.0 + 0.5). (In this test, you'll see 1.3289 as the
 *  decoded "probabilities" for DG4,5,6 on row i=3.) A validation that
 *  looks for probabilities exceeding 1.0 will fail on such DGk's.
 *
 *  This isn't much of a test, per se... more a memorialization of the
 *  issue.
 */
static void
utest_mute_partial_cycle(void)
{
  char          msg[] = "reference_dp mute partial cycle unit test failed";
  ESL_ALPHABET *abc   = esl_alphabet_Create(eslAMINO);
  ESL_DSQ       dsq[] = { eslDSQ_SENTINEL, 1, 2, 3, 7, 8, 9, eslDSQ_SENTINEL }; // e.g. "ACDHIK", digitized
  int           L     = 6;
  int           M     = 9;
  H4_PROFILE   *hmm   = h4_profile_Create(abc, M);
  H4_MODE      *mo    = h4_mode_Create();
  H4_REFMX     *rxf   = h4_refmx_Create(M, L);
  H4_REFMX     *rxd   = h4_refmx_Create(M, L);
  char          errbuf[eslERRBUFSIZE];
  int           k;

  /* Manually create contrived example profile */
  for (k = 1; k <= M; k++)   hmm->e[k][k] = 1.0;   // state 1 = A, 2 = C, 3 = D... : ACDEFGHIK
  for (k = 1; k < M;  k++) { hmm->t[k][h4_TMM] = hmm->t[k][h4_TIM] = hmm->t[k][h4_TDD] = 1.0; }
  hmm->t[3][h4_TMM] = 0.0;
  hmm->t[3][h4_TMD] = 1.0; 
  hmm->t[0][h4_TMD] = hmm->t[6][h4_TDD] = 0.99;
  hmm->t[0][h4_TMM] = hmm->t[6][h4_TDM] = 0.01;
  hmm->flags |= h4_HASPROBS;

  /* Configure it, in multihit glocal L=0 mode */
  if ( h4_standardize(hmm)      != eslOK) esl_fatal(msg);
  if ( h4_mode_SetGlocal(mo)    != eslOK) esl_fatal(msg);
  if ( h4_mode_SetLength(mo, 0) != eslOK) esl_fatal(msg);

  /* posterior decoding of L=6 seq ACDHIK, which aligns as ACD---HIK to model */
  if ( h4_reference_Forward (dsq, L, hmm, mo, rxf, NULL)     != eslOK) esl_fatal(msg);
  if ( h4_reference_Backward(dsq, L, hmm, mo, rxd, NULL)     != eslOK) esl_fatal(msg);
  if ( h4_reference_Decoding(dsq, L, hmm, mo, rxf, rxd, rxd) != eslOK) esl_fatal(msg);

  if ( h4_refmx_Validate(rxd, errbuf) != eslOK) esl_fatal(msg);                                  // Validate() knows DG's can be >1.
  if ( esl_FCompareNew( 1.3289, H4R_MX(rxd, 3, 4, h4R_DG), 0., 0.0001) != eslOK) esl_fatal(msg); // Individually check the ones that are.
  if ( esl_FCompareNew( 1.3289, H4R_MX(rxd, 3, 5, h4R_DG), 0., 0.0001) != eslOK) esl_fatal(msg);
  if ( esl_FCompareNew( 1.3289, H4R_MX(rxd, 3, 6, h4R_DG), 0., 0.0001) != eslOK) esl_fatal(msg);

  //h4_refmx_Dump(stdout, rxd);

  h4_refmx_Destroy(rxd);
  h4_refmx_Destroy(rxf);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_alphabet_Destroy(abc);
}
#endif //h4REFERENCE_DP_TESTDRIVE
/*----------------- end, other unit tests -----------------------*/



/*****************************************************************
 * 9. Test driver
 *****************************************************************/
#ifdef h4REFERENCE_DP_TESTDRIVE

#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                          docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help summary",               0 },
  { "-s",         eslARG_INT,     "0", NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed",      0 },
  { "-L",         eslARG_INT,    "50", NULL, NULL,  NULL,  NULL, NULL, "set test seq length",                   0 },
  { "-M",         eslARG_INT,    "40", NULL, NULL,  NULL,  NULL, NULL, "set test profile length",               0 },
  { "-N",         eslARG_INT,   "100", NULL, NULL,  NULL,  NULL, NULL, "set number of test seqs to run",        0 },
  { "--bruteL",   eslARG_INT,     "4", NULL, NULL,  NULL,  NULL, NULL, "set brute test seq length",             0 },  // brute test is expensive and requires tiny model/seq
  { "--bruteM",   eslARG_INT,     "3", NULL, NULL,  NULL,  NULL, NULL, "set brute test model length",           0 },
  { "--bruteN",   eslARG_INT,    "10", NULL, NULL,  NULL,  NULL, NULL, "number of brute test trials",           0 },
  { "--enumM",    eslARG_INT,     "5", NULL, NULL,  NULL,  NULL, NULL, "set enumeration test model length",     0 },
  { "--diag",     eslARG_STRING, NULL, NULL, NULL,  NULL,  NULL, NULL, "dump data on one utest's failure rate", 0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv, "test driver for reference_dp", "[-options]");
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  int M          = esl_opt_GetInteger(go, "-M");
  int L          = esl_opt_GetInteger(go, "-L");
  int nseq       = esl_opt_GetInteger(go, "-N");
  int enumM      = esl_opt_GetInteger(go, "--enumM");
  int bruteM     = esl_opt_GetInteger(go, "--bruteM");
  int bruteL     = esl_opt_GetInteger(go, "--bruteL");
  int bruteN     = esl_opt_GetInteger(go, "--bruteN");

  if (esl_opt_IsOn(go, "--diag"))  // --diag is for studying error distributions, to set tolerances
    {                              // passing an open <FILE> as first arg to unit test toggles diagnostics mode
      char *which = esl_opt_GetString(go, "--diag");

      if      (strcmp(which, "generation")  == 0) utest_generation     (stdout, rng, eslAMINO, M, L, nseq, /*do_colsum=*/ FALSE);
      else if (strcmp(which, "gencol")      == 0) utest_generation     (stdout, rng, eslAMINO, M, L, nseq, /*do_colsum=*/ TRUE);
      else if (strcmp(which, "duality")     == 0) utest_duality        (stdout, rng, eslAMINO, M, L, nseq);
      else if (strcmp(which, "enumeration") == 0) { while (nseq--) utest_enumeration    (stdout, rng, enumM);          }
      else if (strcmp(which, "singlepath")  == 0) { while (nseq--) utest_singlepath     (stdout, rng, eslAMINO, M);    }
      else if (strcmp(which, "decoding")    == 0) { while (nseq--) utest_approx_decoding(stdout, rng, eslAMINO, M, L); }
      else if (strcmp(which, "brute")       == 0) utest_brute          (stdout, rng, bruteM, bruteL, bruteN);
      else esl_fatal("--diag takes: randomseq, generation, duality, enumeration, singlepath, decoding, brute");
    }
  else
    {
      fprintf(stderr, "## %s\n", argv[0]);
      fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng)); 

      utest_generation     (NULL, rng, eslAMINO, M, L, nseq, FALSE);     // 0.6s
      utest_generation     (NULL, rng, eslAMINO, M, L, nseq, TRUE);     // 0.6s
      utest_duality        (NULL, rng, eslAMINO, M, L, nseq);     // 0.5s
      utest_enumeration    (NULL, rng, enumM);                    // 0.1s
      utest_singlepath     (NULL, rng, eslAMINO, M);              // 0.01s
      utest_brute          (NULL, rng, bruteM, bruteL, bruteN);   // 0.3s
      utest_approx_decoding(NULL, rng, eslAMINO, M, L);           // 1.7s
      utest_mute_partial_cycle();                             // 0.01s

      fprintf(stderr, "#  status   = ok\n");
    }

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif // h4REFERENCE_DP_TESTDRIVE


/*****************************************************************
 * 10. Example
 *****************************************************************/
#ifdef h4REFERENCE_DP_EXAMPLE

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_hmmfile.h"
#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_refmx.h"

#include "general.h"
#include "logsum.h"

#include "reference_dp.h"


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                        docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help",             0 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show HMMER version info",     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of using the Forward/Backward reference implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = h4_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  H4_HMMFILE     *hfp     = NULL;
  H4_PROFILE     *hmm     = NULL;
  H4_MODE        *mo      = h4_mode_Create();
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  H4_PATH        *vpi     = h4_path_Create();
  H4_REFMX       *vit     = h4_refmx_Create(100, 100);
  H4_REFMX       *fwd     = h4_refmx_Create(100, 100);
  H4_REFMX       *bck     = h4_refmx_Create(100, 100);
  H4_REFMX       *pp      = h4_refmx_Create(100, 100);
  float           vsc, fsc, bsc;
  char            errbuf[eslERRBUFSIZE];
  int             status;

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);

  if ( esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, /*env=*/NULL, &sqfp) != eslOK) esl_fatal("couldn't open sequence file %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      h4_mode_SetLength(mo, sq->n);

      h4_reference_Viterbi(sq->dsq, sq->n, hmm, mo, vit, vpi, &vsc);
      h4_reference_Forward(sq->dsq, sq->n, hmm, mo, fwd, &fsc);
      //h4_refmx_Dump(stdout, fwd);
      printf("%s vit %.6f\n", sq->name, vsc);
      printf("%s fwd %.6f\n", sq->name, fsc);

      h4_reference_Backward(sq->dsq, sq->n, hmm, mo, bck, &bsc);
      //h4_refmx_Dump(stdout, bck);
      printf("%s bck %.6f\n", sq->name, bsc);

      h4_reference_Decoding(sq->dsq, sq->n, hmm, mo, fwd, bck, pp);
      h4_refmx_Dump(stdout, pp);
      if (h4_refmx_Validate(pp, errbuf) != eslOK) esl_fatal(errbuf);

      h4_refmx_Reuse(vit);
      h4_refmx_Reuse(fwd);
      h4_refmx_Reuse(bck);
      h4_refmx_Reuse(pp);
      h4_path_Reuse(vpi);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  h4_refmx_Destroy(pp);
  h4_refmx_Destroy(bck);
  h4_refmx_Destroy(fwd);
  h4_refmx_Destroy(vit);
  h4_path_Destroy(vpi);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
}  
#endif //h4REFERENCE_DP_EXAMPLE
/*-------------------- end, example ---------------------------- */
