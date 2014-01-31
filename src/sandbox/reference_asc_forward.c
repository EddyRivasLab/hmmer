
#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_coords2.h"

#include "misc/logsum.h"

#include "dp_reference/p7_refmx.h"

#include "sandbox/reference_asc_forward.h"


/* Function:  p7_ReferenceASCForward()
 * Synopsis:  Calculate an anchor-set-constrained (ASC) Forward score.
 *
 * Purpose:   The ASC Forward algorithm. Given digital sequence <dsq> of
 *            length <L>; profile <gm> to compare it to; and an anchor
 *            set <anch> for <D> domains -- calculate ASC Forward score,
 *            and return it in <*opt_sc>.
 *            
 *            Caller provides two DP matrices <mxu> and <mxd>. They
 *            can be of any allocated size; they will be reallocated
 *            if needed. Upon return, <mxu> and <mxd> contain the ASC
 *            Forward DP matrices for the UP and DOWN sectors of the
 *            calculation, respectively.
 *
 *            Caller must have initialized at least once (per program
 *            invocation) with a <p7_FLogsumInit()> call, because this
 *            function uses <p7_FLogsum()>.
 *            
 *            The two coords in <anch>, <anch[].n1> and <anch[].n2>,
 *            are assigned to (i,k) pairs (in that order). The anchors
 *            in <anch> must be sorted in order of increasing sequence
 *            position <i>.
 *            
 *            <anch> and <D> might be data in a <P7_COORDS2> list
 *            management container: for example, for <P7_COORDS2 *dom>,
 *            you would pass <dom->arr> and <dom->n>.
 *
 * Args:      dsq    : digital target sequence 1..L
 *            L      : length of <dsq>
 *            gm     : profile
 *            anch   : array of (i,k) anchors defining <dsq>'s domain structure
 *            D      : number of anchors in <anch> array -- # of domains
 *            mxu    : UP matrix
 *            mxd    : DOWN matrix
 *            opt_sc : optRETURN - ASC Forward lod score, in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_ReferenceASCForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_COORD2 *anch, int D,
		       P7_REFMX *mxu, P7_REFMX *mxd, float *opt_sc)
{
  const float *tsc;		/* ptr for stepping thru profile's transition parameters */
  const float *rsc;		/* ptr for stepping thru profile's emission parameters   */
  float *dpp;                   /* ptr into main states ({MID}{LG}) for prev row... */
  float *dpc;			/*   ... and a current row. */
  float *xp;			/* ptr into specials (ENJBLGC) for a previous row... */
  float *xc;			/*   ... and a current row. */
  int   d;			/* counter over domains: 0..D-1 */
  int   i;			/* counter over sequence positions (rows): 0,1..L */
  int   k;			/* counter over model positions (columns): 0,1..M */
  int   s;			/* counter over states */
  float mlv, mgv;		/* tmp variables for ML, MG scores... */
  float dlv, dgv;		/*   ... and for DL, DG scores */
  float xE;			/* tmp var for accumulating E score for a row */
  int   iend;			/* tmp var for row index for end of a DOWN sector */
  int   M = gm->M;		/* for a bit more clarity, less dereference clutter */
  int   status;

  /* contract checks / arg validation */
  ESL_DASSERT1( ( gm->L == L || gm->L == 0) ); /* length model in profile is either L (usually) or 0 (some unit tests) */

  /* reallocation, if needed */
  if ( (status = p7_refmx_GrowTo(mxu, gm->M, L)) != eslOK) return status;
  if ( (status = p7_refmx_GrowTo(mxd, gm->M, L)) != eslOK) return status;
  mxu->M    = mxd->M    = M;
  mxu->L    = mxd->L    = L;
  mxu->type = mxd->type = p7R_FORWARD;

  /* Initialize i=0..anch[0].i-1 specials; specials are stored in down matrix */
  for (i = 0; i < anch[0].n1; i++)
    {
      xc = mxd->dp[i] + (M+1) * p7R_NSCELLS; 
      xc[p7R_E]  = -eslINFINITY;
      xc[p7R_N]  = gm->xsc[p7P_N][p7P_LOOP] * i;
      xc[p7R_J]  = -eslINFINITY;
      xc[p7R_B]  = xc[p7R_N] + gm->xsc[p7P_N][p7P_MOVE];
      xc[p7R_L]  = xc[p7R_B] + gm->xsc[p7P_B][0]; 
      xc[p7R_G]  = xc[p7R_B] + gm->xsc[p7P_B][1]; 
      xc[p7R_C]  = -eslINFINITY;
      xc[p7R_JJ] = -eslINFINITY;
      xc[p7R_CC] = -eslINFINITY;
    }

  /* Iterate over domains d=0..D-1: */
  for (d = 0; d < D; d++)
    {
      /*****************************************************************
       * Part 1. UP matrix sector for domain d
       *    In an UP sector, we can enter the model, but not exit;
       *    so we make use of L,G state values from a previous row.
       *
       *    The UP sector includes:
       *       i = anch[d-1].i+1 to anch[d].i  (1..anch[d].i for d==0)
       *       k = 1 to anch[d].k
       *
       *    We'll initialize the top row; then do the remaining rows.
       *****************************************************************/

      /* Initialization of top row, anch[d-1].i+1 (or, 1 for d=0) */
      i   = (d == 0 ? 1 : anch[d-1].n1 + 1);                   // i is always our current row index, 0.1..L 
      rsc = gm->rsc[dsq[i]] + p7P_NR;	                       // Start <rsc> at k=1 on row i 
      tsc = gm->tsc;			                       // Start <tsc> at k=0, where off-by-one {LG}->M transition scores are
      xp  = mxd->dp[i-1] + (M+1) * p7R_NSCELLS;                // <xp> set to the specials of the previous row, i-1
      dpc = mxu->dp[i];                                        // Start <dpc> on k=0 of row i...
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; //   ... initialize that cell, and now <dpc> is on k=1.
      dlv = dgv = -eslINFINITY;

      for (k = 1; k <= anch[d].n2; k++)
	{
	  *dpc++ = mlv = *rsc + xp[p7R_L] + *(tsc + p7P_LM);  // ML.  On top row, only L->ML entry is possible.
	  *dpc++ = mgv = *rsc + xp[p7R_G] + *(tsc + p7P_GM);  // MG.          ... only G->MG.
	  *dpc++ = -eslINFINITY;                              // IL.  Can't be reached on top row.
	  *dpc++ = -eslINFINITY;                              // IG.  Ditto.
	  *dpc++ = dlv;                                       // DL.  Normal calculation, reachable by {MD}->D transitions.
	  *dpc++ = dgv;                                       // DG.

	  tsc += p7P_NTRANS;  
	  rsc += 2;	

	  dlv = p7_FLogsum( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD));
	  dgv = p7_FLogsum( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD));
	}

      /* Now we recurse for remaining rows, down to anch[d].i.
       */
      for (i = i+1 ; i <= anch[d].n1; i++)
	{
	  rsc = gm->rsc[dsq[i]] + p7P_NR;
	  tsc = gm->tsc;
	  xp  = mxd->dp[i-1] + (M+1) * p7R_NSCELLS;
	  dpp = mxu->dp[i-1];
	  dpc = mxu->dp[i];
	  for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; 
	  dlv = dgv = -eslINFINITY;

	  for (k = 1; k <= anch[d].n2; k++)
	    {
	      mlv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_ML) + *(tsc + p7P_MM),
							   *(dpp+p7R_IL) + *(tsc + p7P_IM)),
						p7_FLogsum(*(dpp+p7R_DL) + *(tsc + p7P_DM),
							      xp[p7R_L]  + *(tsc + p7P_LM)));
	      mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_MG) + *(tsc + p7P_MM),
							   *(dpp+p7R_IG) + *(tsc + p7P_IM)),
						p7_FLogsum(*(dpp+p7R_DG) + *(tsc + p7P_DM),
							      xp[p7R_G]  + *(tsc + p7P_GM)));

	      rsc++;              
	      tsc += p7P_NTRANS;  
	      dpp += p7R_NSCELLS;	

	      *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_ML) + *(tsc + p7P_MI), *(dpp + p7R_IL) + *(tsc + p7P_II));
	      *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_MG) + *(tsc + p7P_MI), *(dpp + p7R_IG) + *(tsc + p7P_II));
	      rsc++;	

	      *dpc++ = dlv;
	      *dpc++ = dgv;
	      dlv = p7_FLogsum( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD));
	      dgv = p7_FLogsum( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD));
	    }
	}

      /* The very last cell we calculated was the anchor cell for
       * domain d in the Up matrix, mxu.  dpc just stepped past it, so
       * we can step dpc back to see it again. 
       */
      dpp = dpc - p7R_NSCELLS;

      /*****************************************************************
       * Part 2. DOWN matrix sector for domain d
       *    In a DOWN sector, we can exit the model, but not enter,
       *    so we collect xE on each row,
       *    and use it to set the specials for that row.
       *    
       *    The DOWN sector includes:
       *       i = anch[d].i to anch[d+1].i-1  (anch[d].i to L, for last domain d=D-1)
       *       k = anch[d].k to M
       *       
       *    Again we'll initialize the top row, then do the remaining rows.
       *    The top row starts with the anchor cell, 
       *    which we copy from UP sector calculation we just finished.
       *****************************************************************/


      /* Start with anch[d].k-1 on first row, and set all cells to -inf*/
      i   = anch[d].n1;
      dpc = mxd->dp[i] + (anch[d].n2-1) * p7R_NSCELLS;
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;
      /* Then copy the anchor cell from mxu, in which <dpp> was already set, above. */
      for (s = 0; s < p7R_NSCELLS; s++)	*dpc++ = dpp[s];

      /* Initialization of the rest of the top row from k=anch.k+1 to M 
       * For the D cells on this row: we can do M->D and D->D from the
       * anchor (first) cell, but after that, since we can't enter and
       * and all remaining M's are -inf, only D->D propagates along
       * the row.
       */
                                                // No need to initialize <dpc>, we did that when we copied anchor cell.
      tsc = gm->tsc + anch[d].n2 * p7P_NTRANS;  // Start <tsc> on anch.k, i.e. k-1 relative to start of calculation

      dlv = p7_FLogsum( dpp[p7R_ML] + *(tsc + p7P_MD), dpp[p7R_DL] + *(tsc + p7P_DD));  // DL for anch[d].k + 1 cell, delayed store...
      dgv = p7_FLogsum( dpp[p7R_MG] + *(tsc + p7P_MD), dpp[p7R_DG] + *(tsc + p7P_DD));  // DG, ditto.

      /* xE initialization is not -inf here, because it counts exits
       * from the anchor cell we copied.  It has to watch out for the
       * glocal exit case when the anchor cell (unusually) sits on M
       * itself. And unlike the rest of the cells on the row, M may be
       * finite, so it includes M->E transitions.
       */
      xE  = (anch[d].n2 == M ?
	     p7_FLogsum( p7_FLogsum(dpp[p7R_ML], dpp[p7R_MG]), p7_FLogsum(dpp[p7R_DL], dpp[p7R_DG])) :
	     p7_FLogsum( dpp[p7R_ML], dpp[p7R_DL] ));
 
      for (k = anch[d].n2+1; k <= M; k++)
	{
	  *dpc++ = -eslINFINITY; // ML. No entry, and unreachable from other cells too. 
	  *dpc++ = -eslINFINITY; // MG. Ditto.
	  *dpc++ = -eslINFINITY; // IL. Not reachable on top row. 
	  *dpc++ = -eslINFINITY; // IG. Ditto.
	  *dpc++ = dlv;          // DL. Customary delayed store of prev calculation.
	  *dpc++ = dgv;          // DG. Ditto.

	  tsc   += p7P_NTRANS;
	  
	  xE  = (k == M ?                                  // Glocal exit included if k==M.
		 p7_FLogsum( xE, p7_FLogsum( dlv, dgv)) :  // We know all non-anchor-cell M's are -inf on top row, so 
		 p7_FLogsum( xE, dlv));			   //  we don't include M in these sums.
	  
	  dlv    = dlv + (*tsc + p7P_DD);
	  dgv    = dgv + (*tsc + p7P_DD);
	}

      /* dpc now sits on the start of the specials, in mxd */
      xc = dpc;
      xc[p7R_E]  = xE;
      xc[p7R_N]  = -eslINFINITY; 
      xc[p7R_J]  = (d == D-1 ? -eslINFINITY : xc[p7R_E] + gm->xsc[p7P_E][p7P_LOOP]);
      xc[p7R_B]  = xc[p7R_J] + gm->xsc[p7P_J][p7P_MOVE];
      xc[p7R_L]  = xc[p7R_B] + gm->xsc[p7P_B][0]; 
      xc[p7R_G]  = xc[p7R_B] + gm->xsc[p7P_B][1]; 
      xc[p7R_C]  = (d == D-1 ? xc[p7R_E] + gm->xsc[p7P_E][p7P_LOOP] : -eslINFINITY);
      xc[p7R_JJ] = -eslINFINITY;
      xc[p7R_CC] = -eslINFINITY;

      /* Now we can do the remaining rows in the Down sector of domain d. */
      iend = (d < D-1 ? anch[d+1].n1 : L);
      for (i = i+1 ; i < iend; i++)
	{
	  rsc = gm->rsc[dsq[i]] + anch[d].n2     * p7P_NR;         // Start <rsc> on (x_i, anchor_k, MAT) */
	  tsc = gm->tsc         + (anch[d].n2-1) * p7P_NTRANS;	   // Start <tsc> on (anchor_k-1), to pick up LMk,GMk entries 
	  dpp = mxd->dp[i-1]    + (anch[d].n2-1) * p7R_NSCELLS;	   // Start <dpp> on (i-1, anchor_k-1) 
	  dpc = mxd->dp[i]      + (anch[d].n2-1) * p7R_NSCELLS;	   // Start <dpc> on (i, anchor_k-1)... 
	  for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; //  ... and initialize the k-1 cells to -inf... 
                                                           	   //  ... so, now dpc is on anchor_k.
	  dlv = dgv = xE = -eslINFINITY;

  	  for (k = anch[d].n2; k <= M; k++) 
	    {				  
	      mlv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_ML) + *(tsc + p7P_MM),
							   *(dpp+p7R_IL) + *(tsc + p7P_IM)),
						           *(dpp+p7R_DL) + *(tsc + p7P_DM));
	      mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_MG) + *(tsc + p7P_MM),
							   *(dpp+p7R_IG) + *(tsc + p7P_IM)),
					 	           *(dpp+p7R_DG) + *(tsc + p7P_DM));

	      rsc++;                // rsc advances to insert score for position k 
	      tsc += p7P_NTRANS;    // tsc advances to transitions in states k     
	      dpp += p7R_NSCELLS;   // dpp advances to cells for states k          

	      *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_ML) + *(tsc + p7P_MI), *(dpp + p7R_IL) + *(tsc + p7P_II));
	      *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_MG) + *(tsc + p7P_MI), *(dpp + p7R_IG) + *(tsc + p7P_II));
	      rsc++;		    // rsc advances to next match state emission  

	      xE  = (k == M ?
		     p7_FLogsum( xE, p7_FLogsum( p7_FLogsum(mlv, dlv), p7_FLogsum(mgv, dgv))) : // k=M includes glocal exits  
		     p7_FLogsum( xE, p7_FLogsum(mlv, dlv)));                                    // k<M allows local exit only 

	      *dpc++ = dlv;                                                    // DL. Customary delayed store.
	      *dpc++ = dgv;                                                    //   ... ditto for DG store.
	      dlv = p7_FLogsum( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD)); // Precalculation of DL for next k.
	      dgv = p7_FLogsum( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD)); //   ... ditto for DG calculation.
	    }

	  /*****************************************************************
	   *  Having finished and stored the DOWN calculation on row i, with value xE,
           *  we can calculate and store the specials - also in the DOWN matrix.
           *    dpc[] is already on the special state storage.
	   *****************************************************************/

	  xc = dpc;
	  xp = dpp + p7R_NSCELLS;	
	  xc[p7R_E]  = xE;		
	  xc[p7R_N]  = -eslINFINITY; 
	  xc[p7R_J]  = (d == D-1 ? -eslINFINITY : p7_FLogsum( xp[p7R_J] + gm->xsc[p7P_J][p7P_LOOP], xc[p7R_E] + gm->xsc[p7P_E][p7P_LOOP]));
	  xc[p7R_B]  = xc[p7R_J] + gm->xsc[p7P_J][p7P_MOVE]; 
	  xc[p7R_L]  = xc[p7R_B] + gm->xsc[p7P_B][0]; 
	  xc[p7R_G]  = xc[p7R_B] + gm->xsc[p7P_B][1]; 
	  xc[p7R_C]  = (d == D-1 ? p7_FLogsum( xp[p7R_C] + gm->xsc[p7P_C][p7P_LOOP], xc[p7R_E] + gm->xsc[p7P_E][p7P_LOOP]) : -eslINFINITY);
	  xc[p7R_JJ] = -eslINFINITY;                                                                           
	  xc[p7R_CC] = -eslINFINITY;       
	} /* end loop over rows i of DOWN sector for domain d */

    } /* end loop over domains d=0..D-1; DP calculation complete. */

  /* As we leave the DP recursion, <xc> is still sitting on the
   * special states for the last row L... even for the edge case
   * of D==0 (and the edge case L=0 which must also have D==0).
   */
  if (opt_sc) *opt_sc = xc[p7R_C] + gm->xsc[p7P_C][p7P_MOVE]; /* C->T */
  return eslOK;
}



/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7REFERENCE_ASC_FORWARD_EXAMPLE
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of ASC Forward reference implementation";

int
compute_anchors_from_trace(const P7_REFMX *pp, const P7_TRACE *tr, P7_COORD2 **ret_anch, int *ret_n)
{
  P7_COORD2 *anch;
  int   n;
  float ppv;
  int   d       = 0;
  float best_pp = -1.;
  float *dpc;
  int    z, i,k,s;
  int    status;

  n = tr->ndom;
  ESL_ALLOC(anch, sizeof(P7_COORD2) * n);

  for (z = 0; z < tr->N; z++)
    {
      if (tr->i[z]) i = tr->i[z]; /* keeps track of last i emitted, for when a D state is best */

      if (p7_trace_IsMain(tr->st[z]))
	{
	  k = tr->k[z];
	  
	  dpc = pp->dp[i] + k * p7R_NSCELLS;
	  ppv = 0.0;
	  for (s = 0; s < p7R_NSCELLS; s++) ppv += dpc[s];
	  if (ppv > best_pp)
	    {
	      anch[d].n1 = i;
	      anch[d].n2 = k;
	    }
	}
      else if (tr->st[z] == p7T_E)
	{
	  d++;
	  best_pp = -1.;
	}
    }
  
  *ret_anch = anch;
  *ret_n    = n;
  return eslOK;

 ERROR:
  if (anch) free(anch);
  return status;

}


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
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_REFMX       *vit     = p7_refmx_Create(100, 100);
  P7_REFMX       *fwd     = p7_refmx_Create(100, 100);
  P7_REFMX       *bck     = p7_refmx_Create(100, 100);
  P7_REFMX       *pp      = p7_refmx_Create(100, 100);
  P7_REFMX       *mxu     = p7_refmx_Create(100, 100);
  P7_REFMX       *mxd     = p7_refmx_Create(100, 100);
  P7_TRACE       *tr      = p7_trace_Create();
  P7_COORD2      *anch    = NULL;
  int             D;
  int             d;
  float           fsc, vsc, asc_sc;
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
  p7_profile_Config(gm, hmm, bg);

  while ( (status = esl_sqio_Read(sqfp, sq)) != eslEOF)
    {
      if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
      else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

      /* Set the profile and null model's target length models */
      p7_bg_SetLength     (bg, sq->n);
      p7_profile_SetLength(gm, sq->n);

      p7_ReferenceViterbi (sq->dsq, sq->n, gm, vit, tr, &vsc);
      p7_ReferenceForward (sq->dsq, sq->n, gm, fwd, &fsc);   
      p7_ReferenceBackward(sq->dsq, sq->n, gm, bck, NULL);   
      p7_ReferenceDecoding(sq->dsq, sq->n, gm, fwd, bck, pp);   

      p7_trace_Index(tr);

      //p7_trace_Dump(stdout, tr);

      compute_anchors_from_trace(pp, tr, &anch, &D);

      //      printf("%d anchors/domains in trace:\n", D);
      //for (d = 0; d < D; d++)
      // printf("  domain %3d: anchor at i = %d, k = %d\n", d, anch[d].n1, anch[d].n2);



      p7_ReferenceASCForward(sq->dsq, sq->n, gm, anch, D, mxu, mxd, &asc_sc);

      //p7_refmx_Dump(stdout, mxu);
      //p7_refmx_Dump(stdout, mxd);

      printf("%-30s   %10.4f %10.4f %10.4f %10.4g\n", 
	     sq->name, 
	     vsc, 
	     asc_sc,
	     fsc, 
	     exp(asc_sc - fsc));

      p7_refmx_Reuse(vit);
      p7_refmx_Reuse(fwd);
      p7_refmx_Reuse(mxu);
      p7_refmx_Reuse(mxd);
      p7_trace_Reuse(tr);
      esl_sq_Reuse(sq);
      free(anch);
    }

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(mxu);
  p7_refmx_Destroy(mxd);
  p7_refmx_Destroy(pp);
  p7_refmx_Destroy(bck);
  p7_refmx_Destroy(fwd);
  p7_refmx_Destroy(vit);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_MPL_FWD_EXAMPLE*/




/*****************************************************************
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
