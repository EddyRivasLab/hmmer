/* Forward/Backward implementation variant:
 *   dual-mode alignment (local/glocal);
 *   quadratic memory (simplest variant; not banded, not checkpointed);
 *   "generic" (standard C code; not striped/vectorized);
 *   using P7_GMXD DP matrix structure.
 *   
 * Contents:  
 *   1. Forwards.
 *   2. Backwards.
 *   3. Posterior decoding.
 *   4. Gamma-centroid MGE alignment.
 *   5. Benchmark driver.
 *   6. Unit tests.
 *   7. Test driver.
 *   8. Example.
 *   9. Copyright and license information.
 */


#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_gmxd.h"


/*****************************************************************
 * 1. Forward implementation
 *****************************************************************/

/* Function:  p7_GForwardDual()
 * Synopsis:  Forward, dual-mode, quadratic memory, generic profile
 *
 * Purpose:   The Forward algorithm, comparing profile <gm> to target
 *            sequence <dsq> of length <L>. Caller provides an
 *            allocated <P7_GMXD> DP matrix <gxd>, sized for an
 *            <gm->M> by <L> problem. 
 *            
 *            Caller also has initialized with a <p7_FLogsumInit()>
 *            call; this function will use <p7_FLogsum()>.
 *            
 *            Upon successful return, the raw Forward score (in nats)
 *            is optionally returned in <*opt_sc>, and the DP matrix
 *            <gxd> is filled in.
 *
 * Args:      dsq    : digital target sequence of length <L>
 *            L      : length of the target sequence
 *            gm     : query profile 
 *            gxd    : allocated DP matrix
 *            opt_sc : optRETURN: raw Forward score in nats
 *
 * Returns:   <eslOK> on success.
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
						       xL             + *(tsc + p7P_LM)));

	  mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7GD_MG) + *(tsc + p7P_MM),
						       *(dpp+p7GD_IG) + *(tsc + p7P_IM)),
 					    p7_FLogsum(*(dpp+p7GD_DG) + *(tsc + p7P_DM),
						       xG             + *(tsc + p7P_GM)));

	  rsc++;                /* rsc advances to insert score for position k */
	  tsc += p7P_NTRANS;    /* tsc advances to transitions in states k     */
	  dpp += p7GD_NSCELLS;	/* dpp advances to cells for states k          */

	  /* Insert state calculations ILk, IGk. */
	  *dpc++ = *rsc + p7_FLogsum( *(dpp + p7GD_ML) + *(tsc + p7P_MI), *(dpp + p7GD_IL) + *(tsc + p7P_II));
	  *dpc++ = *rsc + p7_FLogsum( *(dpp + p7GD_MG) + *(tsc + p7P_MI), *(dpp + p7GD_IG) + *(tsc + p7P_II));
	  rsc++;		/* rsc advances to next match state emission   */

	  /* E state update; local paths only, DLk->E, MLk->E; transition prob 1.0 in implicit probability model */
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
						   xL             + *(tsc + p7P_LM)));

      mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7GD_MG) + *(tsc + p7P_MM),
						   *(dpp+p7GD_IG) + *(tsc + p7P_IM)),
                    				   *(dpp+p7GD_DG) + *(tsc + p7P_DM));
      dpp  += p7GD_NSCELLS; 

      /* I_M state doesn't exist      */
      *dpc++ = -eslINFINITY;	/* IL */
      *dpc++ = -eslINFINITY;	/* IG */

      /* E state update now includes glocal exits: transition prob 1.0 from MG_m, DG_m */
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
/*-----------  end, Forwards implementation ---------------------*/



/*****************************************************************
 * 2. Backwards implementation
 *****************************************************************/

/* Function:  p7_GBackwardDual()
 * Synopsis:  Backward, dual-mode, quadratic memory, generic profile
 *
 * Purpose:   The Backward algorithm, comparing profile <gm> to target
 *            sequence <dsq> of length <L>. Caller provides an
 *            allocated <P7_GMXD> DP matrix <gxd>, sized for an
 *            <gm->M> by <L> problem. 
 *            
 *            Caller also has initialized with a <p7_FLogsumInit()>
 *            call; this function will use <p7_FLogsum()>.
 *            
 *            Upon successful return, the raw Backward score (in nats)
 *            is optionally returned in <*opt_sc>, and the DP matrix
 *            <gxd> is filled in.
 *
 * Args:      dsq    : digital target sequence of length <L>
 *            L      : length of the target sequence
 *            gm     : query profile 
 *            gxd    : allocated DP matrix
 *            opt_sc : optRETURN: raw Backward score in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Notes:     In <gm->rsc>, assumes p7P_NR = 2 and order [M I]
 *            In <gm->tsc>, does not make assumptions about p7P_NTRANS or order of values
 *            in <gxd->dp[i]>, assumes p7GD_NSCELLS=6 in order [ ML MG IL IG DL DG]
 *                             assumes p7GD_NXCELLS=7 in order [ E N J B L G C ]
 *                             
 *            Order of evaluation in the code is pretty carefully
 *            arranged to guarantee that dpc,dpn could be pointing
 *            into the same row of memory in a memory-efficient
 *            one-row DP implementation... even though this particular
 *            function, working with a <gxd>, knows it has rows
 *            <0,1..L>. This is so this code could be cribbed for a
 *            one-row implementation.
 */
int
p7_GBackwardDual(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXD *gxd, float *opt_sc)
{
  float *dpc;			/* ptr into current DP row, gxd->dp[i]                    */
  float *dpn;	            	/* ptr into next DP row, gxd->dp[i+1]                     */   // dpc, dpn could point to same row, in a single-row implementation.
  const float *rsc;		/* ptr to current row's residue score x_i vector in <gm>  */
  const float *rsn;		/* ptr to next row's residue score x_{i+1} vector in <gm> */
  const float *tsc;		/* ptr to model transition score vector gm->tsc[]         */
  float dgc, dlc;	        /* DG,DL tmp values on current row, [i,?,DG], [i,?,DL]    */
  float mgc, mlc;
  float mgn, mln;
  float ign, iln;
  float xN, xJ, xC, xE, xG, xL, xB;	/* temp vars for special state values                     */
  int   i;			/* counter over sequence positions 1..L */
  int   k;			/* counter over model positions 1..M    */
  const int M  = gm->M;

  /* Initialize row L. */
  /* Specials are in order ENJBLGC: step backwards thru them */
  dpc  = gxd->dp[L] + (M+1)*p7GD_NSCELLS + p7GD_C;
  rsc  = gm->rsc[dsq[L]] + M*p7P_NR;

  *dpc-- = xC = gm->xsc[p7P_C][p7P_MOVE];      /* C : C<-T */
  *dpc--      = -eslINFINITY;                  /* G : impossible w/o residues after it */
  *dpc--      = -eslINFINITY;	               /* L : ditto, impossible */
  *dpc--      = -eslINFINITY;	               /* B : ditto, impossible */
  *dpc--      = -eslINFINITY;	               /* J : ditto, impossible */
  *dpc--      = -eslINFINITY;	               /* N : ditto, impossible */
  *dpc-- = xE = xC + gm->xsc[p7P_E][p7P_MOVE]; /* E: E<-C<-T, no tail */
  /* dpc is now sitting on [M][DG] */
  
  /* dpc main cells for k=M*/
  tsc = gm->tsc + (M-1)*p7P_NTRANS;
  
  xG   = xE + *rsc + *(tsc + p7P_GM);
  xL   = xE + *rsc + *(tsc + p7P_LM);
  rsc -= p7P_NR;

  *dpc-- = dgc = xE;		/* DG: D_M->E (transition prob 1.0)  */
  *dpc-- = dlc = xE;		/* DL: ditto */
  *dpc-- = -eslINFINITY;	/* IG_M: no such state: always init'ed to -inf */
  *dpc-- = -eslINFINITY;	/* IL_M: no such state: always init'ed to -inf */
  *dpc-- = xE;			/* MG: M_M->E (transition prob 1.0)  */
  *dpc-- = xE;			/* ML: ditto */
  /* dpc is now sitting on [M-1][DG] */

  /* initialize main cells [k=1..M-1] on row i=L*/
  for (k = M-1; k >= 1; k--)
    {
      mgc =                 dgc + *(tsc + p7P_MD);
      mlc =  p7_FLogsum(xE, dlc + *(tsc + p7P_MD));

      xG   = p7_FLogsum(xG, mgc + *rsc + *(tsc + p7P_GM - p7P_NTRANS)); /* off-by-one: tGMk stored as [k-1,GM] */
      xL   = p7_FLogsum(xL, mlc + *rsc + *(tsc + p7P_LM - p7P_NTRANS));
      rsc -= p7P_NR;

      *dpc-- = dgc =                dgc + *(tsc + p7P_DD);  /* DG: only D->D path is possible */
      *dpc-- = dlc = p7_FLogsum(xE, dlc + *(tsc + p7P_DD)); /* DL: Dk->Dk+1 or Dk->E */
      *dpc--       = -eslINFINITY;  	                    /* IG impossible w/o residues following it */
      *dpc--       = -eslINFINITY;	                    /* IL, ditto */
      *dpc--       = mgc;
      *dpc--       = mlc;
      tsc -= p7P_NTRANS;
    }

  /* k=0 cells are -inf */


  /* The main recursion over rows i=L-1 down to 1. (residues x_{L-1} down to x_1) */
  for (i = L-1; i >= 1; i--)
    {
                                                        /* ...xG,xL inherited from previous loop...               */
      rsn = gm->rsc[dsq[i+1]] + M * p7P_NR;        	/* residue x_{i+1} scores in *next* row:  start at end, M */
      rsc = gm->rsc[dsq[i]]   + M * p7P_NR;	        /* residue x_{i} scores in *current* row: start at end, M */
      dpc = gxd->dp[i]   + (M+1)*p7GD_NSCELLS + p7GD_C;	/* dpc is on [i,(M),C]   : end of current row's specials  */
      dpn = gxd->dp[i+1] + (M+1)*p7GD_NSCELLS;	        /* dpn is on [i+1,(M),0] : start of next row's specials   */
      tsc = gm->tsc + M * p7P_NTRANS;		        /* tsc is on t[M,0]: vector [MM IM DM LM GM MD DD MI II]  */

      /* Calculation of the special states. */
      /* dpc is on dp[i][C] special, will now step backwards thru [E N J B L G C ] */
      *dpc-- = xC = *(dpn + p7GD_C) + gm->xsc[p7P_C][p7P_LOOP]; /* C = C<-C */

      *dpc-- = xG;     /* G was calculated during prev row (G->Mk wing unfolded) */
      *dpc-- = xL;     /* L was calculated during prev row */

      *dpc-- = xB = p7_FLogsum( xG + gm->xsc[p7P_B][1],    /* B<-G */
				xL + gm->xsc[p7P_B][0]);   /* B<-L */
      
      *dpc-- = xJ = p7_FLogsum( *(dpn + p7GD_J) + gm->xsc[p7P_J][p7P_LOOP],   /* J<-J */
				xB              + gm->xsc[p7P_J][p7P_MOVE]);  /* J<-B */
      
      *dpc--      = p7_FLogsum( *(dpn + p7GD_N) + gm->xsc[p7P_N][p7P_LOOP],  /* N<-N */
				xB              + gm->xsc[p7P_N][p7P_MOVE]); /* N<-B */

      *dpc-- = xE = p7_FLogsum( xC + gm->xsc[p7P_E][p7P_MOVE],
				xJ + gm->xsc[p7P_E][p7P_LOOP]);
      dpn -= 5;		      /* backs dpn up to [i+1,M,MG], skipping [IL IG DL DG] at k=M */
      

      /* Initialization of the k=M states */
      /* dpc on [i,k=M,DG], init at k=M, step back thru [ ML MG IL IG DL DG] */
      /* dpn on [i+1,k=M,MG] */
      mgn = *rsn + *dpn--;	/* pick up MG(i+1,k=M) + s(x_i+1,k=M, M) */
      mln = *rsn + *dpn--;	/* pick up ML(i+1,k=M) + s(x_i+1,k=M, M) */
      rsn--;			/* rsn now on s(x_i+1, k=M-1, I)         */

      xG     = xE + *rsc + *(tsc + p7P_GM - p7P_NTRANS); /* t[k-1][GM] is G->Mk wing-folded entry, recall off-by-one storage   */
      xL     = xE + *rsc + *(tsc + p7P_LM - p7P_NTRANS); /* t[k-1][LM] is L->Mk uniform local entry */
      rsc -= p7P_NR;		/* rsc now on s[x_{i},M-1,M] */
      tsc -= p7P_NTRANS;	/* tsc now on t[M-1,0]       */

      *dpc-- = dgc = xE;		/* DGm->E */
      *dpc-- = dlc = xE;		/* DLm->E */
      *dpc--       = -eslINFINITY;	/* IGm nonexistent */
      *dpc--       = -eslINFINITY;	/* ILm nonexistent */
      *dpc--       = xE;		/* MGm->E */
      *dpc--       = xE;		/* MLm->E */
      /* dpc on [i,M-1,DG]; dpn on [i+1,M-1,DG] */


      /* The main recursion over model positions k=M-1 down to 1. */
      for (k = M-1; k >= 1; k--)
	{
                             	    /* rsn is on residue score [x_{i+1},k,I]    */
	                            /* dpn is on [i+1,k,DG]                     */
	  dpn -= 2;	   	    /* skip DG/DL values on next row            */
	  ign = *dpn--;		    /* pick up IG value from dp[i+1]            */ // if inserts had nonzero score: + *rsn 
	  iln = *dpn--;		    /* pick up IL value                         */ // if inserts had nonzero score: + *rsn
	  rsn--;		    /* skip residue score for I (zero)          */ 
	                            /* dpn is now sitting on dp[i+1,k,MG]       */

                                                                 /* tsc is on tsc[k,0] */
	  mgc =  p7_FLogsum( p7_FLogsum(mgn + *(tsc + p7P_MM),   /* mgn = [i+1,k+1,MG] */
					ign + *(tsc + p7P_MI)),  /* ign = [i+1,k,  IG] */
 			                dgc + *(tsc + p7P_MD));  /* dgc = [i,  k+1,DG] */

	  mlc =  p7_FLogsum( p7_FLogsum(mln + *(tsc + p7P_MM),   /* mln = [i+1,k+1,ML] */
					iln + *(tsc + p7P_MI)),  /* iln = [i+1,k,  IL] */
			     p7_FLogsum(dlc + *(tsc + p7P_MD),   /* dlc = [i,  k+1,DL] */
					xE));                    /* ML->E trans = 1.0  */

	  xG   = p7_FLogsum(xG, mgc + *rsc + *(tsc + p7P_GM - p7P_NTRANS)); /* t[k-1][GM] is G->Mk wing-retracted glocal entry */
	  xL   = p7_FLogsum(xL, mlc + *rsc + *(tsc + p7P_LM - p7P_NTRANS)); /* t[k-1][LM] is L->Mk uniform local entry         */
	  rsc -= p7P_NR;				       /* rsc now on s[x_i, k-1, M] */

	  /* dpc is on [i,k,DG] and will now step backwards thru: [ ML MG IL IG DL DG ] */
	  *dpc-- = dgc = p7_FLogsum( mgn + *(tsc + p7P_DM),   /* dgc picked up for next loop of k */
				     dgc + *(tsc + p7P_DD));
	  *dpc-- = dlc = p7_FLogsum( p7_FLogsum( mln + *(tsc + p7P_DM),   /* dlc picked up for next loop of k */
						 dlc + *(tsc + p7P_DD)),
				     xE);

	  *dpc-- = p7_FLogsum( mgn + *(tsc + p7P_IM),
			       ign + *(tsc + p7P_II));
	  *dpc-- = p7_FLogsum( mln + *(tsc + p7P_IM),
			       iln + *(tsc + p7P_II));

                                /* recall that dpn is on dp[i+1][k,MG]    */
	  mgn = *rsn + *dpn--;	/* pick up M[i+1,k]; add score[x_i+1,k,M] */
	  mln = *rsn + *dpn--;
	  rsn--;		/* rsn is now on score[i+1,k-1,I] */

	  *dpc-- = mgc;		/* delayed store of [i,k,MG] value enables dpc,dpn to point into same single row */
	  *dpc-- = mlc;

	  tsc -= p7P_NTRANS;

	  /* as we loop around now and decrement k:
           *   dpn is on [i+1,k-1,DG] which becomes [i+1,k,DG] 
           *   dpc is on [i,k-1,DG]   which becomes [i,k,DG] 
	   *   tsc is on tsc[k-1,0]   which becomes tsc[k,0]
	   *   rsn is on s[i+1,k-1,I] which becomes s[i+1,k,I]
	   *   rsc is on s[i,  k-1,M] which becomes s[i,k,M]
	   *   dgc is [i,k,DG],   which becomes [i,k+1,DG] value  (and analog. for dlc,DL)
	   *   mgn is [i+1,k,MG], which becomes [i+1,k+1,MG] value (and analog. for ML)
	   */
	} /* end of loop over model positions k */

      /* k=0 cells are -inf */

      /* xG,xL values are now ready for next row */
    } /* end of loop over rows i. */
  /* now on row i=0. Only N,B,G,L states are reachable on this initial row. G,L values are already done. */
  
  dpc = gxd->dp[0] + (M+1)*p7GD_NSCELLS + p7GD_C;	/* dpc is on [0,(M),C] : end of row 0 specials  */
  dpn = gxd->dp[1] + (M+1)*p7GD_NSCELLS;	        /* dpn is on [1,(M),0] : start of row 1 specials   */

  *dpc--      = -eslINFINITY;                                           /* C */
  *dpc--      = xG;                                                     /* G */
  *dpc--      = xL;                                                     /* L */
  *dpc-- = xB = p7_FLogsum( xG + gm->xsc[p7P_B][1],                     /* B */
			    xL + gm->xsc[p7P_B][0]);   
  *dpc--      = -eslINFINITY;                                           /* J */
  *dpc-- = xN = p7_FLogsum( *(dpn + p7GD_N) + gm->xsc[p7P_N][p7P_LOOP],	/* N */
			    xB              + gm->xsc[p7P_N][p7P_MOVE]); 
  *dpc--      = -eslINFINITY;                                           /* E */

  gxd->M = M;
  gxd->L = L;
  if (opt_sc) *opt_sc = xN;
  return eslOK;
}
/*-------------- end, backwards implementation ------------------*/


/*****************************************************************
 * 3. Posterior decoding
 *****************************************************************/

int
p7_GDecodingDual(const P7_PROFILE *gm, const P7_GMXD *fwd, P7_GMXD *bck, P7_GMXD *pp)
{
  const int L = fwd->L;
  const int M = fwd->M;
  float xN, xJ, xC;
  float    *fwdp;
  float    *bckp;
  float    *ppp;
  float     denom;
  float     sc = P7_GMXD_XMX(fwd,L,p7GD_C) + gm->xsc[p7P_C][p7P_MOVE];
  int       i;
  int       k;
  int       x;

  ppp  = pp->dp[0];
  for (x = 0; x < p7GD_NSCELLS * (M+1) + p7GD_NXCELLS; x++) *ppp++ = 0.0;
  xN = 0.0;
  xJ = xC = -eslINFINITY;

  for (i = 1; i <= L; i++)
    {
      fwdp  = fwd->dp[i] + p7GD_NSCELLS;
      bckp  = bck->dp[i] + p7GD_NSCELLS;
      ppp   = pp->dp[i];
      denom = 0.0;
      for (x = 0; x < p7GD_NSCELLS; x++) *ppp++ = 0.0;

      for (k = 1; k < M; k++)
	{
	  /* [ ML MG IL IG DL DG] */
	  *ppp = expf(*fwdp + *bckp - sc); denom += *ppp++; fwdp++;  bckp++;  /* ML */
	  *ppp = expf(*fwdp + *bckp - sc); denom += *ppp++; fwdp++;  bckp++;  /* MG */
	  *ppp = expf(*fwdp + *bckp - sc); denom += *ppp++; fwdp++;  bckp++;  /* IL */
	  *ppp = expf(*fwdp + *bckp - sc); denom += *ppp++; fwdp+=3; bckp+=3; /* IG; skip DL,DG */
	  *ppp++ = 0.0;		/* DL */
	  *ppp++ = 0.0;		/* DG */
	}
      *ppp = expf(*fwdp + *bckp - sc); denom += *ppp++; fwdp++;  bckp++;  /* ML_m */
      *ppp = expf(*fwdp + *bckp - sc); denom += *ppp++; fwdp++;  bckp++;  /* MG_m */
      *ppp++ = 0.0;
      *ppp++ = 0.0;
      *ppp++ = 0.0; 
      *ppp++ = 0.0; fwdp += 5; bckp += 5; /* +4, +1 to skip E; fwd,bck now on N */

      /* [ E N J B L G C ] */
      *ppp++ = 0.0;		/* E */
      *ppp   = expf(xN + *bckp + gm->xsc[p7P_N][p7P_LOOP] - sc); denom += *ppp++; xN = *fwdp++;  bckp++;  /* N */
      *ppp   = expf(xJ + *bckp + gm->xsc[p7P_J][p7P_LOOP] - sc); denom += *ppp++; xJ = *fwdp++;  bckp++; /* J; fwdp, bckp advance to C */
      *ppp++ = 0.0;		/* B */
      *ppp++ = 0.0;		/* L */
      *ppp++ = 0.0; 	        /* G */
      fwdp += 3; bckp += 3;
      *ppp   = expf(xC + *bckp + gm->xsc[p7P_C][p7P_LOOP] - sc); denom += *ppp;   xC = *fwdp;  /* C */
      /* note delayed store for N/J/C, because it's fwd[i-1][X] + bck[i-1][X] + tXX */

      denom = 1.0 / denom;	/* multiplication faster than division... */
      ppp = pp->dp[i] + p7GD_NSCELLS;
      for (k = 1; k < M; k++) { 
	*ppp++ *= denom;
	*ppp++ *= denom;
	*ppp++ *= denom;
	*ppp   *= denom;
	ppp += 3;
      }
      *ppp++ *= denom;
      *ppp++ *= denom;
      ppp += 5;
      
      *(ppp + p7GD_N) *= denom;
      *(ppp + p7GD_J) *= denom;
      *(ppp + p7GD_C) *= denom;
    }
  
  pp->M = M;
  pp->L = L;
  return eslOK;
}
/*-------------- end, posterior decoding ------------------------*/



/*****************************************************************
 * 4. Gamma-centroid MGE alignment: fill
 *****************************************************************/

#define P7_DELTAT(val, tsc) ( ((tsc) == -eslINFINITY) ? -eslINFINITY : (val))

/* Function:  p7_GCentroidAlignDual()
 * Synopsis:  Gamma centroid MGE alignment against a dual mode model
 *
 * Purpose:   
 *
 * Args:      
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      [HamadaAsai11] for details of gamma-centroid estimators
 *            [Kall05] for more on why the delta function on transitions is needed
 */
int
p7_GCentroidAlignDual(float gamma, const P7_PROFILE *gm, const P7_GMXD *pp, P7_GMXD *gxd)
{
  float       *dpp;		     /* ptr into previous DP row in <gxd> */
  float       *dpc;		     /* ptr into current DP row in <gxd> */
  const float *tsc;		     /* ptr into transition scores in <gm> */
  const float *ppp;		     /* ptr into decoded posterior probs in <pp> */
  const float  gam1 = gamma + 1.0;   /* residue score is  pp(x_i) * (gamma + 1.0) - 1.0 */
  const int    L    = pp->L;
  const int    M    = pp->M;
  float        dlv, dgv;
  float        mlv, mgv;
  float        xE,xG,xL,xN,xJ,xB;
  int          i,k,s;

  /* Initialization of the zero row. 
   * Main states    [ ML MG IL IG DL DG] 0..M; then special states [E N J B L G C] 
   * Gamma-centroid values can be negative, so "impossible" is -eslINFINITY.
   * P7_DELTAT() causes impossible transitions to result in impossible-scoring paths
   */
  dpc = gxd->dp[0];
  for (s = 0; s < (M+1) * p7GD_NSCELLS; s++) *dpc++ = -eslINFINITY;   /* all M,I,D; k=0..M */
  *dpc++ = -eslINFINITY;	               /* E */
  *dpc++ = 0.0;			               /* N */
  *dpc++ = -eslINFINITY;                       /* J */
  *dpc++ = 0.0;				       /* B  assumes tNB > 0.0, which must be true */
  *dpc++ = P7_DELTAT( 0.0, gm->xsc[p7P_B][0]); /* L */
  *dpc++ = P7_DELTAT( 0.0, gm->xsc[p7P_B][1]); /* G */
  *dpc   = -eslINFINITY;		       /* C */

  /* Main DP recursion for rows 1..L */
  for (i = 1; i <= L; i++)
    {
      ppp = pp->dp[i] + p7GD_NSCELLS; /* positioned at pp[i,k=1,ML]     */
      tsc = gm->tsc;		      /* model on k=0 transitions (k-1 as we enter loop) */
      dpp = gxd->dp[i-1];	      /* prev row on k=0 (k-1 as we enter loop) */
      dpc = gxd->dp[i];		      
      for (s = 0; s < p7GD_NSCELLS; s++) *dpc++ = -eslINFINITY;

      dlv = dgv = xE = -eslINFINITY;

      for (k = 1; k < M; k++)
	{ /* main states [ ML MG IL IG DL DG] */ 	 
	  /* ML calculation */   
	  mlv = *dpc++ = (*ppp++ * gam1 - 1.) + 
	    ESL_MAX( ESL_MAX( P7_DELTAT(*(dpp + p7GD_ML), *(tsc + p7P_MM)),
			      P7_DELTAT(*(dpp + p7GD_IL), *(tsc + p7P_IM))),
		     ESL_MAX( P7_DELTAT(*(dpp + p7GD_DL), *(tsc + p7P_DM)),
			      P7_DELTAT( xL,              *(tsc + p7P_LM))));

	  /* MG calculation */
	  mgv = *dpc++ = (*ppp++ * gam1 - 1.) + 
	    ESL_MAX( ESL_MAX( P7_DELTAT(*(dpp + p7GD_MG), *(tsc + p7P_MM)),
			      P7_DELTAT(*(dpp + p7GD_IG), *(tsc + p7P_IM))),
		     ESL_MAX( P7_DELTAT(*(dpp + p7GD_DG), *(tsc + p7P_DM)),
			      P7_DELTAT( xG,              *(tsc + p7P_GM))));

	  tsc += p7P_NTRANS;	/* transition scores now on gm->tsc[k] transitions */
	  dpp += p7GD_NSCELLS;	/* prev row now on dp[i-1][k] cells       */

	  /* IL/IG calculation. */
	  *dpc++ = (*ppp++ * gam1 - 1.) +
	    ESL_MAX( P7_DELTAT(*(dpp + p7GD_ML), *(tsc + p7P_MI)),
		     P7_DELTAT(*(dpp + p7GD_IL), *(tsc + p7P_II)));

	  *dpc++ = (*ppp++ * gam1 - 1.) +
	    ESL_MAX( P7_DELTAT(*(dpp + p7GD_MG), *(tsc + p7P_MI)),
		     P7_DELTAT(*(dpp + p7GD_IG), *(tsc + p7P_II)));

	  /* E state update with ML->E and DL->E local exits k<M */
	  xE = ESL_MAX(xE, ESL_MAX(mlv, dlv));

	  /* Delete states: delayed storage trick */
	  *dpc++ = dlv;
	  *dpc++ = dgv;

	  dlv = ESL_MAX( P7_DELTAT( mlv, *(tsc + p7P_MD)),
			 P7_DELTAT( dlv, *(tsc + p7P_DD)));
	  dgv = ESL_MAX( P7_DELTAT( mgv, *(tsc + p7P_MD)),
			 P7_DELTAT( dgv, *(tsc + p7P_DD)));

	  ppp += 2;		/* skip pp past DL/DG, which are zero and unused */
	}

      /* k=M is unrolled: no I state, and glocal exits. */
      mlv = *dpc++ = (*ppp++ * gam1 - 1.) + 
	ESL_MAX( ESL_MAX( P7_DELTAT(*(dpp + p7GD_ML), *(tsc + p7P_MM)),
			  P7_DELTAT(*(dpp + p7GD_IL), *(tsc + p7P_IM))),
		 ESL_MAX( P7_DELTAT(*(dpp + p7GD_DL), *(tsc + p7P_DM)),
			  P7_DELTAT( xL,              *(tsc + p7P_LM))));
      mgv = *dpc++ = (*ppp++ * gam1 - 1.) + 
	ESL_MAX( ESL_MAX( P7_DELTAT(*(dpp + p7GD_MG), *(tsc + p7P_MM)),
			  P7_DELTAT(*(dpp + p7GD_IG), *(tsc + p7P_IM))),
		 ESL_MAX( P7_DELTAT(*(dpp + p7GD_DG), *(tsc + p7P_DM)),
			  P7_DELTAT( xG,              *(tsc + p7P_GM))));
      *dpc++ = 0.0;	/* IL */
      *dpc++ = 0.0;	/* IG */
      *dpc++ = dlv;
      *dpc++ = dgv;
      ppp   += 4;

      /* Special states [E N J B L G C] */
      /* row i main states finished; dpc, ppp are now on first special state, E */
      dpp += 2*p7GD_NSCELLS;	/* and now dpp is too */
      
      *dpc++ = xE = ESL_MAX( xE, ESL_MAX( ESL_MAX(mlv, dlv),    /* E */
			    	          ESL_MAX(mgv, dgv)));

      *dpc++ = xN =          P7_DELTAT( *(dpp + p7GD_N) + *(ppp + p7GD_N) * gam1 - 1., gm->xsc[p7P_N][p7P_LOOP]);
      
      *dpc++ = xJ = ESL_MAX( P7_DELTAT( *(dpp + p7GD_J) + *(ppp + p7GD_J) * gam1 - 1., gm->xsc[p7P_J][p7P_LOOP]),
			     P7_DELTAT( xE,                                            gm->xsc[p7P_E][p7P_LOOP]));
      
      *dpc++ = xB = ESL_MAX( P7_DELTAT( xN, gm->xsc[p7P_N][p7P_MOVE]),
			     P7_DELTAT( xJ, gm->xsc[p7P_J][p7P_MOVE]));
      
      *dpc++ = xL = P7_DELTAT(xB, gm->xsc[p7P_B][0]);
      *dpc++ = xG = P7_DELTAT(xB, gm->xsc[p7P_B][1]);
      *dpc        = ESL_MAX( P7_DELTAT( *(dpp + p7GD_C) + *(ppp + p7GD_C) * gam1 - 1., gm->xsc[p7P_C][p7P_LOOP]), 
			     P7_DELTAT( xE,                                            gm->xsc[p7P_E][p7P_MOVE]));
    }
  return eslOK;
}
/*--------------- end, MGE alignment DP -------------------------*/


/*****************************************************************
 * 5. Gamma-centroid alignment: traceback
 *****************************************************************/


static inline int
select_mg(const P7_PROFILE *gm, const P7_GMXD *gxd, int i, int k)
{
  int   state[4] = { p7T_MG, p7T_IG, p7T_DG, p7T_G };
  float path[4];

  path[0] = P7_DELTAT( P7_GMXD_MX(gxd, i-1, k-1, p7GD_MG), P7P_TSC(gm, k-1, p7P_MM));
  path[1] = P7_DELTAT( P7_GMXD_MX(gxd, i-1, k-1, p7GD_IG), P7P_TSC(gm, k-1, p7P_IM));
  path[2] = P7_DELTAT( P7_GMXD_MX(gxd, i-1, k-1, p7GD_DG), P7P_TSC(gm, k-1, p7P_DM));
  path[3] = P7_DELTAT( P7_GMXD_XMX(gxd, i-1, p7GD_G),      P7P_TSC(gm, k-1, p7P_GM));
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
select_ml(const P7_PROFILE *gm, const P7_GMXD *gxd, int i, int k)
{
  int   state[4] = { p7T_ML, p7T_IL, p7T_DL, p7T_L };
  float path[4];

  path[0] = P7_DELTAT( P7_GMXD_MX(gxd, i-1, k-1, p7GD_ML), P7P_TSC(gm, k-1, p7P_MM));
  path[1] = P7_DELTAT( P7_GMXD_MX(gxd, i-1, k-1, p7GD_IL), P7P_TSC(gm, k-1, p7P_IM));
  path[2] = P7_DELTAT( P7_GMXD_MX(gxd, i-1, k-1, p7GD_DL), P7P_TSC(gm, k-1, p7P_DM));
  path[3] = P7_DELTAT( P7_GMXD_XMX(gxd, i-1, p7GD_L),      P7P_TSC(gm, k-1, p7P_LM));
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
select_ig(const P7_PROFILE *gm, const P7_GMXD *gxd, int i, int k)
{
  float path[2];

  path[0] = P7_DELTAT( P7_GMXD_MX(gxd, i-1, k, p7GD_MG), P7P_TSC(gm, k, p7P_MI));
  path[1] = P7_DELTAT( P7_GMXD_MX(gxd, i-1, k, p7GD_IG), P7P_TSC(gm, k, p7P_II));
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_IG);
}

static inline int
select_il(const P7_PROFILE *gm, const P7_GMXD *gxd, int i, int k)
{
  float path[2];

  path[0] = P7_DELTAT( P7_GMXD_MX(gxd, i-1, k, p7GD_ML), P7P_TSC(gm, k, p7P_MI));
  path[1] = P7_DELTAT( P7_GMXD_MX(gxd, i-1, k, p7GD_IL), P7P_TSC(gm, k, p7P_II));
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_IL);
}

static inline int
select_dg(const P7_PROFILE *gm, const P7_GMXD *gxd, int i, int k)
{
  float path[2];

  path[0] = P7_DELTAT( P7_GMXD_MX(gxd, i, k-1, p7GD_MG), P7P_TSC(gm, k-1, p7P_MD));
  path[1] = P7_DELTAT( P7_GMXD_MX(gxd, i, k-1, p7GD_DG), P7P_TSC(gm, k-1, p7P_DD));
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_DG);
}

static inline int
select_dl(const P7_PROFILE *gm, const P7_GMXD *gxd, int i, int k)
{
  float path[2];

  path[0] = P7_DELTAT( P7_GMXD_MX(gxd, i, k-1, p7GD_ML), P7P_TSC(gm, k-1, p7P_MD));
  path[1] = P7_DELTAT( P7_GMXD_MX(gxd, i, k-1, p7GD_DL), P7P_TSC(gm, k-1, p7P_DD));
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_DL);
}

static inline int
select_e(const P7_PROFILE *gm, const P7_GMXD *gxd, int i, int *ret_k)
{
  float *dpc = gxd->dp[i] + p7GD_NSCELLS; /* position at k=1, ML */
  float max  = -eslINFINITY;
  int   smax = -1;
  int   kmax = -1;
  int   k;

  for (k = 1; k < gm->M; k++)
    { /* [ML MG IL IG DL DG] */
      if (*dpc > max) { max = *dpc; smax = p7T_ML; kmax = k; }   dpc+=4;
      if (*dpc > max) { max = *dpc; smax = p7T_DL; kmax = k; }   dpc+=2;
    }
  if (*dpc > max) { max = *dpc; smax = p7T_ML; kmax = k; }       dpc++;
  if (*dpc > max) { max = *dpc; smax = p7T_MG; kmax = k; }       dpc+=3;
  if (*dpc > max) { max = *dpc; smax = p7T_DL; kmax = k; }       dpc++;
  if (*dpc > max) { max = *dpc; smax = p7T_DG; kmax = k; }       

  *ret_k = kmax;
  return smax;
}

static inline int
select_n(int i)
{
  return ((i==0) ? p7T_S : p7T_N);
}

static inline int
select_j(const P7_PROFILE *gm, const P7_GMXD *gxd, int i, const P7_GMXD *pp, float gamma)
{
  float path[2];
  float rsc = P7_GMXD_XMX(pp, i, p7GD_J) * (gamma + 1.) - 1.0;

  path[0] = P7_DELTAT ( rsc + P7_GMXD_XMX(gxd, i-1, p7GD_J), gm->xsc[p7P_J][p7P_LOOP]);
  path[1] = P7_DELTAT (       P7_GMXD_XMX(gxd, i,   p7GD_E), gm->xsc[p7P_E][p7P_LOOP]);
  return ( (path[0] > path[1]) ? p7T_J : p7T_E);
}

static inline int
select_b( const P7_PROFILE *gm, const P7_GMXD *gxd, int i)
{
  float path[2];

  path[0] = P7_DELTAT( P7_GMXD_XMX(gxd, i, p7GD_N), gm->xsc[p7P_N][p7P_MOVE]);
  path[1] = P7_DELTAT( P7_GMXD_XMX(gxd, i, p7GD_J), gm->xsc[p7P_J][p7P_MOVE]);
  return ( (path[0] > path[1]) ? p7T_N : p7T_J);
}

static inline int
select_c(const P7_PROFILE *gm, const P7_GMXD *gxd, int i, const P7_GMXD *pp, float gamma)
{
  float path[2];
  float rsc = P7_GMXD_XMX(pp, i, p7GD_C) * (gamma + 1.) - 1.0;

  path[0] = P7_DELTAT ( rsc + P7_GMXD_XMX(gxd, i-1, p7GD_C), gm->xsc[p7P_C][p7P_LOOP]);
  path[1] = P7_DELTAT (       P7_GMXD_XMX(gxd, i,   p7GD_E), gm->xsc[p7P_E][p7P_MOVE]);
  return ( (path[0] > path[1]) ? p7T_C : p7T_E);
}


int
p7_GCentroidTrace(float gamma, const P7_PROFILE *gm, const P7_GMXD *pp, const P7_GMXD *gxd, P7_TRACE *tr)
{
  int   i    = gxd->L;
  int   k    = 0;
  int   sprv = p7T_C;
  int   scur;
  float ppv;
  int   status;

  if ((status = p7_trace_AppendWithPP(tr, p7T_T, k, i, 0.0)) != eslOK) return status;
  if ((status = p7_trace_AppendWithPP(tr, p7T_C, k, i, 0.0)) != eslOK) return status;

  while (sprv != p7T_S)
    {
      switch (sprv) {
      case p7T_ML: scur = select_ml(gm, gxd, i, k); k--; i--; break;
      case p7T_MG: scur = select_mg(gm, gxd, i, k); k--; i--; break;
      case p7T_IL: scur = select_il(gm, gxd, i, k);      i--; break;
      case p7T_IG: scur = select_ig(gm, gxd, i, k);      i--; break;
      case p7T_DL: scur = select_dl(gm, gxd, i, k); k--;      break;
      case p7T_DG: scur = select_dg(gm, gxd, i, k); k--;      break;
      case p7T_E:  scur = select_e (gm, gxd, i, &k);          break;
      case p7T_N:  scur = select_n (         i    );          break;
      case p7T_J:  scur = select_j (gm, gxd, i, pp, gamma);   break;
      case p7T_B:  scur = select_b (gm, gxd, i);              break;
      case p7T_L:  scur = p7T_B;                              break;
      case p7T_G:  scur = p7T_B;                              break;
      case p7T_C:  scur = select_c (gm, gxd, i, pp, gamma);   break;
      default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in traceback");
      }

      switch (scur) {		/* posterior probs of residues are marginalized over glocal/local path choice */
      case p7T_ML:  ppv = P7_GMXD_MX(pp, i, k, p7GD_ML) + P7_GMXD_MX(pp, i, k, p7GD_MG); break; 
      case p7T_MG:  ppv = P7_GMXD_MX(pp, i, k, p7GD_ML) + P7_GMXD_MX(pp, i, k, p7GD_MG); break;
      case p7T_IL:  ppv = P7_GMXD_MX(pp, i, k, p7GD_IL) + P7_GMXD_MX(pp, i, k, p7GD_IG); break;
      case p7T_IG:  ppv = P7_GMXD_MX(pp, i, k, p7GD_IL) + P7_GMXD_MX(pp, i, k, p7GD_IG); break;
      case p7T_N:   ppv = (sprv==scur ? P7_GMXD_XMX(pp, i, p7GD_N) : 0.0); break;
      case p7T_J:   ppv = (sprv==scur ? P7_GMXD_XMX(pp, i, p7GD_J) : 0.0); break;
      case p7T_C:   ppv = (sprv==scur ? P7_GMXD_XMX(pp, i, p7GD_C) : 0.0); break;
      default:      ppv = 0.0;  break;
      }

      /* A glocal B->G->Mk wing-retraction entry: unfold it */
      if (scur == p7T_G) {
	while (k > 1) {
	  if ( (status = p7_trace_AppendWithPP(tr, p7T_DG, k-1, i, 0.0)) != eslOK) return status;
	  k--;
	}
      }

      if ( (status = p7_trace_AppendWithPP(tr, scur, k, i, ppv)) != eslOK) return status;

      /* For NCJ, we had to defer i decrement. */
      if ( (scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;
      sprv = scur;
    }
  
  tr->M = gxd->M;
  tr->L = gxd->L;
  return p7_trace_Reverse(tr);
}
/*-------------- end, MGE alignment traceback -------------------*/

/*****************************************************************
 * 6. Benchmark driver.
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
      if (! esl_opt_GetBoolean(go, "-F"))  p7_GBackwardDual(dsq, L, gm, bck, NULL);

      p7_gmxd_Reuse(fwd);
      p7_gmxd_Reuse(bck);
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
 * 6. Unit tests
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
 * 7. Test driver
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
 * 8. Example
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
  { "--csv",  eslARG_OUTFILE,   NULL,  NULL, NULL,   NULL,  NULL, NULL, "output CSV file of decoding matrix, for heat map", 0 },
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
  char           *csvfile = esl_opt_GetString(go, "--csv");
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMXD        *fwd     = NULL;
  P7_GMXD        *bck     = NULL;
  P7_TRACE       *tr      = NULL;
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
  bck = p7_gmxd_Create(gm->M, sq->n);
  tr  = p7_trace_CreateWithPP();

  printf("%-30s   %-10s %-10s   %-10s %-10s\n", "# seq name",      "fwd (raw)",   "bck (raw) ",  "fwd (bits)",  "bck (bits)");
  printf("%-30s   %10s %10s   %10s %10s\n",     "#--------------", "----------",  "----------",  "----------",  "----------");

  while ( (status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
      else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

      /* Resize the DP matrices if necessary */
      p7_gmxd_GrowTo(fwd, gm->M, sq->n);
      p7_gmxd_GrowTo(bck, gm->M, sq->n);

      /* Set the profile and null model's target length models */
      p7_bg_SetLength     (bg, sq->n);
      p7_profile_SetLength(gm, sq->n);

      //p7_profile_Dump(stdout, gm);

      /* Run Forward, Backward */
      p7_GForwardDual (sq->dsq, sq->n, gm, fwd, &fsc);

      if (esl_opt_GetBoolean(go, "--vv")) p7_gmxd_Dump(stdout, fwd);

      p7_GBackwardDual(sq->dsq, sq->n, gm, bck, &bsc);
      p7_GDecodingDual(gm, fwd, bck, bck);                   /* <bck> now becomes the pp matrix */

      // quickie check added in a hurry, for a lab mtg stat
      printf("Total cells: %" PRId64 "\n", (int64_t) sq->n * (int64_t) gm->M);
      {
	int   i,k;
	float val;
	int64_t ncells = 0;
	float pthresh = 0.001;
	P7_GMXD *pp = bck;

	for (i = 1; i <= sq->n; i++)
	  for (k = 1; k <= gm->M; k++)
	    {
	      val =  
		pp->dp[i][k * p7GD_NSCELLS + p7GD_ML] + pp->dp[i][k * p7GD_NSCELLS + p7GD_MG] + 
		pp->dp[i][k * p7GD_NSCELLS + p7GD_IL] + pp->dp[i][k * p7GD_NSCELLS + p7GD_IG] + 
		pp->dp[i][k * p7GD_NSCELLS + p7GD_DL] + pp->dp[i][k * p7GD_NSCELLS + p7GD_DG];
	      if (val >= pthresh) ncells++;
	    }
	printf("cells over thresh: %" PRId64 "\n", ncells);
      }


      if (csvfile) {
	FILE *csvfp = fopen(csvfile, "w");
	p7_gmxd_DumpCSV(csvfp, bck, 1, sq->n, 1, gm->M);
	fclose(csvfp);
      }

      if (esl_opt_GetBoolean(go, "--vv")) p7_gmxd_Dump(stdout, bck);

      p7_GCentroidAlignDual(/*gamma=*/0.5, gm, bck, fwd);    /* <fwd> is now the centroid DP alignment fill */

      if (esl_opt_GetBoolean(go, "--vv")) p7_gmxd_Dump(stdout, fwd);

      p7_GCentroidTrace(/*gamma=*/0.5, gm, bck, fwd, tr);

      p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);

      /* Those scores are partial log-odds likelihoods in nats.
       * Subtract off the rest of the null model, convert to bits.
       */
      p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

      printf("%-30s   %10.4f %10.4f   %10.4f %10.4f\n", 
	     sq->name, 
	     fsc, bsc, 
	     (fsc - nullsc) / eslCONST_LOG2, (bsc - nullsc) / eslCONST_LOG2);

      p7_gmxd_Reuse(fwd);
      p7_gmxd_Reuse(bck);
      p7_trace_Reuse(tr);
      esl_sq_Reuse(sq);
    }

  /* Cleanup */
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
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
