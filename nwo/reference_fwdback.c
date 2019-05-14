/* Reference implementation of Forward and Backward algorithms
 * for dual-mode glocal/local model.
 *
 * All reference implementation code is for testing. It is not used in
 * HMMER's main executables. Sparse DP code is the production version.
 *
 * Contents:
 *   1. Forward
 *   2. Backward
 *   3. Benchmark driver
 *   4. "Brute" unit test
 *   5. Other unit tests
 *   6. Test driver
 *   7. Example
 */
#include "h4_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "h4_mode.h"
#include "h4_profile.h"
#include "h4_refmx.h"
#include "logsum.h"

#include "reference_fwdback.h"


/*****************************************************************
 * 1. Forward
 *****************************************************************/

/* Function:  h4_ReferenceForward()
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
h4_ReferenceForward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx, float *opt_sc)
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
  if ( (status = h4_refmx_GrowTo(rmx, hmm->M, L)) != eslOK) return status;
  rmx->M    = M;
  rmx->L    = L;
  rmx->type = h4R_FORWARD;

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
}
/*-----------  end, Forwards implementation ---------------------*/


/*****************************************************************
 * 2. Backward
 *****************************************************************/

/* Function:  h4_ReferenceBackward()
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
h4_ReferenceBackward(const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, H4_REFMX *rmx, float *opt_sc)
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
  if ( (status = h4_refmx_GrowTo(rmx, M, L)) != eslOK) return status;
  rmx->M    = M;
  rmx->L    = L;
  rmx->type = h4R_BACKWARD;

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
}
/*-------------- end, backwards implementation ------------------*/


/*****************************************************************
 * x. Benchmark driver
 *****************************************************************/
#ifdef h4REFERENCE_FWDBACK_BENCHMARK

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
#include "reference_fwdback.h"
#include "reference_viterbi.h"

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
  ESL_RANDOMNESS *rng     = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
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

  h4_logsum_Init();

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile");
  h4_hmmfile_Close(hfp);
  h4_profile_Config(hmm);
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
      h4_ReferenceForward (dsq, L, hmm, mo, rx, &sc);
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
      h4_ReferenceBackward (dsq, L, hmm, mo, rx, &sc);
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
      h4_ReferenceViterbi (dsq, L, hmm, mo, rx, pi, &sc);
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
#endif /*h4REFERENCE_FWDBACK_BENCHMARK*/
/*----------------- end, benchmark ------------------------------*/





/*****************************************************************
 * x. "Brute" unit test
 *****************************************************************/
#ifdef h4REFERENCE_FWDBACK_TESTDRIVE

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
#include "esl_sq.h"
#include "esl_stack.h"

#include "h4_path.h"
#include "emit.h"
#include "logsum.h"
#include "modelsample.h"
#include "reference_viterbi.h"

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
 */
static void
utest_brute(ESL_RANDOMNESS *rng, int M, int L, int ntrials, int be_verbose)
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

  while (ntrials--)
    {
      /* sample a random model of length M; half the time, zeropepper it for extra pathology */
      if (esl_rnd_Roll(rng, 2) == 0) h4_modelsample_zeropeppered(rng, abc, M, &hmm);
      else                           h4_modelsample(rng, abc, M, &hmm);

      /* configure the model and alignment mode; sample a length L sequence from model */
      if ( h4_profile_Config(hmm)   != eslOK) esl_fatal(msg);
      if ( h4_mode_SetLength(mo, L) != eslOK) esl_fatal(msg);
      do {
	if (h4_emit(rng, hmm, mo, sq, true_pi) != eslOK) esl_fatal(msg);
      } while (sq->n != L);

      if (be_verbose) h4_profile_Dump(stdout, hmm);
      if (be_verbose) h4_mode_Dump(stdout, mo);

      /* run Forward, Backward, Viterbi */
      if (h4_ReferenceForward (sq->dsq, sq->n, hmm, mo, rmx,       &fsc) != eslOK) esl_fatal(msg);
      if (be_verbose) h4_refmx_Dump(stdout, rmx);
      h4_refmx_Reuse(rmx);
      if (h4_ReferenceBackward(sq->dsq, sq->n, hmm, mo, rmx,       &bsc) != eslOK) esl_fatal(msg);
      h4_refmx_Reuse(rmx);
      if (h4_ReferenceViterbi (sq->dsq, sq->n, hmm, mo, rmx,  vpi, &vsc) != eslOK) esl_fatal(msg);

      /* Backward and Forward scores should match */
      if (fabs(fsc - bsc) > ftol_abs) esl_fatal(msg);

      /* Score all brute paths with this seq/profile  */
      for (i = 0; i < npath; i++)
	if ( h4_path_Score(arrpi[i], sq->dsq, hmm, mo, &(pathsc[i])) != eslOK) esl_fatal(msg);

      /* dump all paths and their scores before we sort scores for accurate summation */
      if (be_verbose)
	{
	  printf("number of paths  = %d\n", npath);
	  for (i = 0; i < npath; i++)
	    {
	      printf("%d : %9.5f ", i, pathsc[i]);
	      h4_path_DumpCigar(stdout, arrpi[i]);
	    }
	}

      /* Get brute force vsc, fsc by maximization, summation over paths */
      esl_vec_FSortIncreasing(pathsc, npath);        // After this sort, we've messed up the indexing of individual paths
      brute_vsc = pathsc[npath-1];
      brute_fsc = esl_vec_FLog2Sum(pathsc, npath);   // Sorting scores from small to large increases numerical accuracy 

      if (be_verbose)
	{
	  printf("brute_vsc = %10.6f  vsc = %10.5f  abs_err = %12.9f  rel_err = %10.6f  seed = %" PRIu32 "\n",
		 brute_vsc, vsc, brute_vsc - vsc, fabs((brute_vsc - vsc) / vsc), esl_randomness_GetSeed(rng));
	  printf("brute_fsc = %10.6f  fsc = %10.5f  abs_err = %12.9f  rel_err = %10.6f  seed = %" PRIu32 "\n",
		 brute_fsc, fsc, brute_fsc - fsc, fabs((brute_fsc - fsc) / fsc), esl_randomness_GetSeed(rng));
	}


      /* The main test: reference DP calculations are acceptably close to brute force calculations */
      if (fabs(brute_vsc - vsc) > vtol_abs) esl_fatal(msg);
      if (fabs(brute_fsc - fsc) > ftol_abs) esl_fatal(msg);

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
#endif // h4REFERENCE_FWDBACK_TESTDRIVE

/*****************************************************************
 * x. Test driver
 *****************************************************************/
#ifdef h4REFERENCE_FWDBACK_TESTDRIVE

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                          docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help summary",             0 },
  { "-L",         eslARG_INT,     "4", NULL, NULL,  NULL,  NULL, NULL, "set test seq length",                 0 },
  { "-M",         eslARG_INT,     "3", NULL, NULL,  NULL,  NULL, NULL, "set test model length",               0 },
  { "--nbrute",   eslARG_INT,    "10", NULL, NULL,  NULL,  NULL, NULL, "number of brute test trials",         0 },
  { "-v",         eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "verbose: dump debugging diagnostics", 0 },
  { "--seed",     eslARG_INT,     "0", NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed",    0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv,
					    "test driver for reference_fwdback",
					    "[-options]");
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  int bruteM     = esl_opt_GetInteger(go, "-M");
  int bruteL     = esl_opt_GetInteger(go, "-L");
  int bruteN     = esl_opt_GetInteger(go, "--nbrute");
  int be_verbose = esl_opt_GetBoolean(go, "-v");

  h4_logsum_Init();

  esl_fprintf(stderr, "## %s\n", argv[0]);
  esl_fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_brute(rng, bruteM, bruteL, bruteN, be_verbose);

  esl_fprintf(stderr, "#  status   = ok\n");

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}

#endif // h4REFERENCE_FWDBACK_TESTDRIVE


/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef h4REFERENCE_FWDBACK_EXAMPLE

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
#include "reference_fwdback.h"

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
  H4_MODE        *mo      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  H4_REFMX       *rmx     = h4_refmx_Create(100, 100);
  float           fsc, bsc;
  int             status;

  h4_logsum_Init();

  if ( h4_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) esl_fatal("couldn't open profile file %s", hmmfile);
  if ( h4_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) esl_fatal("failed to read profile from file %s", hmmfile);
  h4_hmmfile_Close(hfp);

  if ( esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, /*env=*/NULL, &sqfp) != eslOK) esl_fatal("couldn't open sequence file %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  h4_profile_Config(hmm);
  mo = h4_mode_Create();  // default mode = dual-mode multihit glocal/local

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      h4_mode_SetLength(mo, sq->n);

      h4_ReferenceForward(sq->dsq, sq->n, hmm, mo, rmx, &fsc);
      //h4_refmx_Dump(stdout, rmx);
      h4_refmx_Reuse(rmx);
      printf("%s fwd %.6f\n", sq->name, fsc);

      h4_ReferenceBackward(sq->dsq, sq->n, hmm, mo, rmx, &bsc);
      //h4_refmx_Dump(stdout, rmx);
      printf("%s bck %.6f\n", sq->name, bsc);

      h4_refmx_Reuse(rmx);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  h4_refmx_Destroy(rmx);
  h4_profile_Destroy(hmm);
  h4_mode_Destroy(mo);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
}  
#endif //h4REFERENCE_FWDBACK_EXAMPLE
/*-------------------- end, example ---------------------------- */
