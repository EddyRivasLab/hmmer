/* Reference implementation of MPL Forward, a restricted Forward
 * algorithm for calculating the likelihood of a domain "labeling" of
 * a target sequence.
 * 
 * MPL = most probable labeling.
 * 
 * All reference implementation code is for testing. It is not used
 * in HMMER's main executables. Sparse DP code is the production
 * version.
 * 
 * This code derives from the standard Forward implementation in
 * reference_fwdback.c.
 * 
 * Contents:
 *   1. MPL Forward.
 */

#include <p7_config.h>

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_coords2.h"

#include "misc/logsum.h"

#include "dp_reference/p7_refmx.h"

#include "sandbox/reference_mpl_fwd.h"


/*****************************************************************
 * 1. MPL Forward
 *****************************************************************/

/* Function:  p7_ReferenceMPLForward()
 * Synopsis:  MPL-restricted Forward algorithm; reference implementation
 *
 * Purpose:   Compare profile <gm> to target sequence <dsq> of length
 *            <L>, using an MPL-restricted Forward algorithm that sums
 *            over those paths consistent with <ndom> homology regions
 *            (domains) with start/end coords given by the list in
 *            <dom>.
 *            
 *            Return the MPL-restricted Forward lod score in
 *            <*opt_sc>, in nats.
 *            
 *            Caller provides an allocated <P7_REFMX> DP matrix <rmx>.
 *            This DP matrix will be filled in by the MPL-restricted
 *            Forward calculation. This matrix will be reallocated as
 *            needed. Caller does not need to worry about DP matrix
 *            size, so the matrix can be a smaller or bigger one,
 *            reused from a previous calculation of any reference
 *            algorithm (Forward, Backward, etc; restricted or not).
 *            
 *            Caller must have initialized at least once (per program
 *            invocation) with a <p7_FLogsumInit()> call, because this
 *            function uses <p7_FLogsum()>.
 *            
 *            The domain start/end positions in <dom> must be sorted
 *            and nonoverlapping.
 *            
 *            <dom> and <ndom> might be data in a <P7_COORDS2> list
 *            management container: for example, for <P7_COORDS2 *c2>,
 *            you would pass <c2->arr> and <c2->n>.
 *
 * Args:      dsq    : digital target sequence dsq[1..L] 
 *            L      : length of <dsq> in residues
 *            gm     : query profile
 *            dom    : dom[0..ndom-1].{start,end} coords of each domain
 *            ndom   : # of domains in dom[0..ndom-1]; >=1 
 *            rmx    : RESULT: filled DP matrix
 *            opt_sc : optRETURN: raw MPL Forward score, nats
 *
 * Returns:   <eslOK> on success, and <rmx> contains the MPL-restricted
 *            Forward matrix; its internals may have been reallocated.
 *
 * Throws:    <eslEMEM> if reallocation of <rmx> is attempted and fails.
 *            Now <*opt_sc> (if provided) is 0, and <rmx> is left the
 *            way the caller provided it.
 *
 * Notes:     Makes assumptions about the order of state indices in p7_refmx.h:
 *              main states: ML MG IL IG DL DG
 *              specials:    E N J B L G C JJ CC
 *              
 *            If L=0, the score is -infinity, by construction; HMMER
 *            profiles generate sequences of L>=1.
 *            
 * Xref:      Derived from reference_fwdback.c::p7_ReferenceForward().           
 */
int
p7_ReferenceMPLForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, 
		       P7_COORD2 *dom, int ndom,
		       P7_REFMX *rmx, float *opt_sc)
{
  float *dpc, *dpp;   /* ptr into current, previous DP matrix row              */
  const float *tsc;   /* ptr for stepping thru profile's transition parameters */
  const float *rsc;   /* ptr for stepping thru profile's emission parameters   */
  int    M = gm->M;   /* just for clearer code below                           */
  int    d = 0;       /* domain # as we step thru dom[d], 0..ndom-1.           */
  int    i, k, s;     /* seq pos 1..L; profile pos 0,1..M; state idx           */
  float  mlv, mgv;    /* ML,MG cell values on current row                      */
  float  dlv, dgv;    /* pushed-ahead DL,DG cell k+1 values                    */
  float  xE, xL, xG;	      
  int    status;

  /* Contract checks, arg validation */
  ESL_DASSERT1( (gm->L == L || gm->L == 0) );  /* length model in profile is either L (usually) or 0 (some unit tests) */
  ESL_DASSERT1( (ndom > 0) );		       /* has to be at least one domain */
  
  /* Reallocation, if needed */
  if ( (status = p7_refmx_GrowTo(rmx, gm->M, L)) != eslOK) return status;
  rmx->M    = M;
  rmx->L    = L;
  rmx->type = p7R_FORWARD;	/* I don't think we have to worry (yet) about marking difference between regular vs. MPL */

  /* Initialization of row 0 */
  dpc = rmx->dp[0];
  for (s = 0; s < (M+1) * p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; // all M,I,D; k=0..M
  dpc[p7R_E]      = -eslINFINITY;	                           
  dpc[p7R_N]      = 0.0;			                   
  dpc[p7R_J]      = -eslINFINITY;                                  
  dpc[p7R_B]      = (dom[d].n1 == 1) ? gm->xsc[p7P_N][p7P_MOVE] : -eslINFINITY; /* MPL restriction. B only before a domain starts */
  dpc[p7R_L] = xL = gm->xsc[p7P_N][p7P_MOVE] + gm->xsc[p7P_B][0];  
  dpc[p7R_G] = xG = gm->xsc[p7P_N][p7P_MOVE] + gm->xsc[p7P_B][1];  
  dpc[p7R_C]      = -eslINFINITY;                                  
  dpc[p7R_JJ]     = -eslINFINITY; /* JJ,CC values aren't used in Forward calculations; they'll all be -inf */
  dpc[p7R_CC]     = -eslINFINITY;                             
  /* *dpc is on the specials' subrow; must be there to make L=0 case work, where loop below is skipped */

  /* note about <d>:
   *    <d> = current domain index if <i> is in a domain;
   *          else it's the index of the *next* domain.
   *          When i > the last domain's end, d == ndom.
   */

  for (i = 1; i <= L; i++)
    {
      /* Initialization for a new row */
      rsc = gm->rsc[dsq[i]] + p7P_NR;	/* this ptr steps through profile's emission scores 1..M for this row; skipping k=0 */
      tsc = gm->tsc;			/* this ptr steps through profile's transition scores 0..M                          */

      dpp = rmx->dp[i-1];	        /* prev DP row: start at k=0.                          */
      dpc = rmx->dp[i];                 /* current DP row, init k=0 to -inf and start at k=1.  */
      for (s = 0; s < p7R_NSCELLS; s++) *dpc++ = -eslINFINITY;

      dlv = dgv = -eslINFINITY;
      xE  =       -eslINFINITY;

      /* MPL restriction: case 1.  
       *    We are within a domain (and not on its first row): i > start,
       *    and i <= end.  Here, we know that both the prev and curr row
       *    have valid values for all M,I,D.  Thus we can do the standard
       *    Fwd calculation with two modifications:
       *      1. no entries from L, G: because that can only happen in i==start.
       *      2. We collect xE, but we will only store it if i==end.
       */
      if (d < ndom && i > dom[d].n1 && i <= dom[d].n2)
	{
	  /* Main inner loop of the standard Forward recursion. */
	  for (k = 1; k < M; k++)
	    {
	      /* match states MLk, MGk */
	      mlv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_ML) + *(tsc + p7P_MM),
							   *(dpp+p7R_IL) + *(tsc + p7P_IM)),
						           *(dpp+p7R_DL) + *(tsc + p7P_DM)); /* MPL: no entry from L */
	      mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_MG) + *(tsc + p7P_MM),
							   *(dpp+p7R_IG) + *(tsc + p7P_IM)),
				  		           *(dpp+p7R_DG) + *(tsc + p7P_DM)); /* MPL: no entry from G */

	      rsc++;                /* rsc advances to insert score for position k */
	      tsc += p7P_NTRANS;    /* tsc advances to transitions in states k     */
	      dpp += p7R_NSCELLS;   /* dpp advances to cells for states k          */

	      /* Insert state calculations ILk, IGk. */
	      *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_ML) + *(tsc + p7P_MI), *(dpp + p7R_IL) + *(tsc + p7P_II));
	      *dpc++ = *rsc + p7_FLogsum( *(dpp + p7R_MG) + *(tsc + p7P_MI), *(dpp + p7R_IG) + *(tsc + p7P_II));
	      rsc++;		/* rsc advances to next match state emission   */

	      /* E state update; local paths only, DLk->E, MLk->E; transition prob 1.0 in implicit probability model 
	       */
	      xE  = p7_FLogsum( p7_FLogsum(mlv, dlv), xE); /* MPL: xE will be unused unless i==end */

	      /* Delete state, deferred storage trick */
	      *dpc++ = dlv;
	      *dpc++ = dgv;
	      dlv = p7_FLogsum( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD));
	      dgv = p7_FLogsum( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD));
	    } /* end loop over k up to M-1 *dpc is now ready to do k==M. *dpp is on M-1.  */

	  /* k=M node is unrolled and handled separately. No I state, and glocal exits. */
	  mlv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_ML) + *(tsc + p7P_MM),
						       *(dpp+p7R_IL) + *(tsc + p7P_IM)),
					               *(dpp+p7R_DL) + *(tsc + p7P_DM));
	  mgv = *dpc++ = *rsc + p7_FLogsum( p7_FLogsum(*(dpp+p7R_MG) + *(tsc + p7P_MM),
						       *(dpp+p7R_IG) + *(tsc + p7P_IM)),
					               *(dpp+p7R_DG) + *(tsc + p7P_DM));
	  dpp  += p7R_NSCELLS; 

	  /* I_M state doesn't exist      */
	  *dpc++ = -eslINFINITY;    /* IL */
	  *dpc++ = -eslINFINITY;    /* IG */

	  /* E state update now includes glocal exits: transition prob 1.0 from MG_m, DG_m */
	  xE  = p7_FLogsum( xE, p7_FLogsum( p7_FLogsum(mlv, dlv), p7_FLogsum(mgv, dgv)));
      
	  /* D_M state: deferred storage only */
	  *dpc++ = dlv;
	  *dpc++ = dgv;
    
	  /* row i is now finished, and dpc[] is positioned exactly on first special state, E */
	  dpp += p7R_NSCELLS;    /* now dpp[] is also positioned exactly on first special,  E */
	}

      /* MPL restriction: case 2.
       *    We're on the first row of domain d.
       *    We can reach match states only on {LG}-> entry transitions.
       *    We can't reach an insert state yet.
       *    We can reach delete states as usual.
       *    dpp main states are unused: we'll advance it once, immediately to its specials.
       */
      else if (d < ndom && i == dom[d].n1)
	{
	  for (k = 1; k < M; k++)
	    {
	      mlv = *dpc++ = *rsc + xL + *(tsc + p7P_LM); /* MPL: L->MLk only */
	      mgv = *dpc++ = *rsc + xG + *(tsc + p7P_GM); /* MPL: G->MGk only */
	  
	      rsc += 2;		    /* MPL: +2 because we immediately skip I emission */
	      tsc += p7P_NTRANS;    /* tsc advances to transitions in states k        */
	      /* MPL: no advance on dpp */

	      *dpc++ = -eslINFINITY; /* ILk */
	      *dpc++ = -eslINFINITY; /* IGk */

	      /* E state update; local paths only, DLk->E, MLk->E; transition prob 1.0 in implicit probability model */
	      xE  = p7_FLogsum( p7_FLogsum(mlv, dlv), xE);

	      /* Delete state, deferred storage trick */
	      *dpc++ = dlv;
	      *dpc++ = dgv;
	      dlv = p7_FLogsum( mlv + *(tsc + p7P_MD), dlv + *(tsc + p7P_DD));
	      dgv = p7_FLogsum( mgv + *(tsc + p7P_MD), dgv + *(tsc + p7P_DD));
	    }
	  /* k=M node unrolled and handled separately */
	  *dpc++ = mlv = *rsc + xL + *(tsc + p7P_LM); /* MPL: L->MLk only */
	  *dpc++ = mgv = *rsc + xG + *(tsc + p7P_GM); /* MPL: G->MGk only */
	  *dpc++       = -eslINFINITY; /* ILk */
	  *dpc++       = -eslINFINITY; /* IGk */
	  
	  /* E state update now includes glocal exits: transition prob 1.0 from MG_m, DG_m */
	  xE  = p7_FLogsum( xE, p7_FLogsum( p7_FLogsum(mlv, dlv), p7_FLogsum(mgv, dgv)));
   
	  /* D_M state: deferred storage only */
	  *dpc++ = dlv;
	  *dpc++ = dgv;

	  /* row i is now finished, and dpc[] is positioned exactly on first special state, E */
	  dpp += (M+1) * p7R_NSCELLS; /* MPL: now dpp[] is also positioned exactly on E */
	}

      /* MPL restriction, case 3:
       *    We're not in a domain. d = index of our *next* domain. (So i < dom[d].n1)
       *    Here, no main state is possible.
       *    We set everything to -inf but the only reason to do that is prettification, so 
       *    the matrix looks good in debugging/development. The MPL Forward code doesn't 
       *    depend on this initialization. (Like Forward itself, it *does* depend on 
       *    the cells in k=0 being init'ed to -inf but we already did that.)
       */
      else 
	{
	  for (s = 0; s < M * p7R_NSCELLS; s++) *dpc++ = -eslINFINITY; /* see note as above: we only do this for prettification */
	  dpp += (M+1) * p7R_NSCELLS; /* and now both dpc[] and dpp[] are positioned exactly on E */
	}
  
      /* Specials. For all three MPL cases -- we can now calculate the specials at the end of *dpc row */
      dpc[p7R_E]      = (d < ndom && i == dom[d].n2)  ? xE                                    : -eslINFINITY;
      dpc[p7R_N]      = (d == 0   && i < dom[0].n1) ? dpp[p7R_N] + gm->xsc[p7P_N][p7P_LOOP] : -eslINFINITY; /* note: yes, that's dom[0]! */

      if      (d < ndom && i < dom[d].n1) dpc[p7R_J] = dpp[p7R_J] + gm->xsc[p7P_J][p7P_LOOP];
      else if (d < ndom && i == dom[d].n2)  dpc[p7R_J] = dpc[p7R_E] + gm->xsc[p7P_E][p7P_LOOP]; 
      else                                   dpc[p7R_J] = -eslINFINITY;

      if      (d == ndom-1 &&  i == dom[d].n2) dpc[p7R_C] = dpc[p7R_E] + gm->xsc[p7P_E][p7P_MOVE];
      else if (d >= ndom)                       dpc[p7R_C] = dpp[p7R_C] + gm->xsc[p7P_C][p7P_LOOP];
      else                                      dpc[p7R_C] = -eslINFINITY;

      /* Advance d to next domain, before doing the B/L/G states, which need to test for whether we're on row before next domain start */
      if (d < ndom && i == dom[d].n2) d++;	
	
      if (d < ndom && i+1 == dom[d].n1)
	{
	  dpc[p7R_B]      = p7_FLogsum( dpc[p7R_N] + gm->xsc[p7P_N][p7P_MOVE], dpc[p7R_J] + gm->xsc[p7P_J][p7P_MOVE]); 
	  dpc[p7R_L] = xL =             dpc[p7R_B] + gm->xsc[p7P_B][0]; 
	  dpc[p7R_G] = xG =             dpc[p7R_B] + gm->xsc[p7P_B][1]; 
	}
      else
	{
	  dpc[p7R_B]      = -eslINFINITY;
	  dpc[p7R_L] = xL = -eslINFINITY;
	  dpc[p7R_G] = xG = -eslINFINITY;
	}
    } /* done with row i */
  /* Done with all rows i. As we leave, dpc is still sitting on the final special value for i=L ... including even the L=0 case */
 
  if (opt_sc) *opt_sc = dpc[p7R_C] + gm->xsc[p7P_C][p7P_MOVE]; /* C->T */
  return eslOK;
}

/*****************************************************************
 * x. Example
 *****************************************************************/
#ifdef p7REFERENCE_MPL_FWD_EXAMPLE
#include <p7_config.h>

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
static char banner[] = "example of MPL Forward reference implementation";

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
  P7_REFMX       *fwd     = p7_refmx_Create(100, 100);
  P7_REFMX       *vit     = p7_refmx_Create(100, 100);
  P7_REFMX       *mpl     = p7_refmx_Create(100, 100);
  P7_TRACE       *tr      = p7_trace_Create();
  P7_COORD2      *dom     = NULL;
  float           fsc, vsc, mplsc;
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

      /* Create a coord2 list, crudely. We can do better later. */
      p7_trace_Index(tr);
      ESL_ALLOC(dom, sizeof(P7_COORD2) * tr->ndom);
      for (d = 0; d < tr->ndom; d++)
	{
	  dom[d].n1 = tr->sqfrom[d];
	  dom[d].n2 = tr->sqto[d];
	}

      p7_ReferenceMPLForward(sq->dsq, sq->n, gm, dom, tr->ndom, mpl, &mplsc);

      printf("%-30s   %10.4f %10.4f %10.4f %10.4g\n", 
	     sq->name, 
	     vsc, 
	     mplsc,
	     fsc, 
	     exp(mplsc - fsc));

      p7_refmx_Reuse(vit);
      p7_refmx_Reuse(fwd);
      p7_refmx_Reuse(mpl);
      p7_trace_Reuse(tr);
      esl_sq_Reuse(sq);
      free(dom);
    }

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_refmx_Destroy(mpl);
  p7_refmx_Destroy(fwd);
  p7_refmx_Destroy(vit);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;

 ERROR:
  return status;
}
#endif /*p7REFERENCE_MPL_FWD_EXAMPLE*/







