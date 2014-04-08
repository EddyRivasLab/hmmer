/* Reference implementation of anchor set constrained (ASC) 
 * posterior decoding.
 * 
 * All reference implementation code is for development and
 * testing. It is not used in HMMER's main executables. Production
 * code uses sparse dynamic programming.
 * 
 * Contents:
 *    1. ASC Decoding
 *    x. Example
 *    x. Copyright and license information
 */

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_coords2.h"

#include "misc/logsum.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_asc_fwdback.h"
#include "dp_reference/reference_asc_decoding.h"

/*****************************************************************
 * 1. ASC Decoding
 *****************************************************************/

/* Function:  p7_ReferenceASCDecoding()
 * Synopsis:  Anchor-set-constrained (ASC) posterior decoding.
 *
 * Purpose:   The anchor set constrained (ASC) posterior decoding
 *            algorithm.  Given digital sequence <dsq> of length <L>,
 *            profile <gm> to compare it to, an anchor set <anch> for
 *            <D> domains, the already calculated ASC Forward UP and
 *            DOWN matrices <afu> and <afd> and ASC Backward UP and
 *            DOWN matrices <abu> and <abd>; do posterior decoding, to
 *            create the decoding UP and DOWN matrices <apu> and
 *            <apd>.
 *            
 *            The caller may provide empty matrices for <apu> and
 *            <apd>.  They can be of any allocated size, and they will
 *            be reallocated as needed.
 *            
 *            Alternatively, the caller can overwrite the ASC Backward 
 *            matrices by passing <abu> for <apu>, and <abd> for
 *            <apd>. That is, the call would look like <(... afu, afd,
 *            abu, abd, abu, abd)>. 
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
 * Args:      dsq  : digital target sequence 1..L
 *            L    : length of <dsq>
 *            gm   : profile
 *            anch : array of (i,k) anchors defining <dsq>'s domain structure
 *            D    : number of anchors in <anch> array = # of domains in <dsq>
 *            afu  : ASC Forward UP matrix
 *            afd  : ASC Forward DOWN matrix
 *            abu  : ASC Backward UP matrix
 *            abd  : ASC Backward DOWN matrix
 *            apu  : RESULT : ASC Decoding UP matrix   (can be <abu>)
 *            apd  : RESULT : ASC Decoding DOWN matrix (can be <abd>)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_ReferenceASCDecoding(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_COORD2 *anch, int D, 
			const P7_REFMX *afu, const P7_REFMX *afd, P7_REFMX *abu, P7_REFMX *abd, P7_REFMX *apu, P7_REFMX *apd)
{
  const float *tsc = gm->tsc;	/* activates TSC() convenience macro, used in G->Mk wing unfolding */
  const float *rsc;		/* ptr to current row's residue emission scores in <gm>            */
  int          d, i, k, s;	/* indices for D domains, L residues, M nodes, and 6 states        */
  int          iend;		/* tmp var for start or end of a chunk in UP or DOWN matrices      */
  float        totsc;		/* overall Backward (=Forward) score, our normalization factor     */
  float        denom;		/* sum of pp for all emitting states for this residue, for renorm  */
  float       *fwdp;		/* ptr into row of Forward DP matrix, UP or DOWN (<afu> or <afd>)  */
  float       *bckp;		/* ptr into row of Backward DP matrix, UP or DOWN (<abu>, <abd>)   */
  float       *ppp;		/* ptr into row of Decoding DP matrix, <apu> or <apd>              */
  float        delta;		/* piece of probability alloted to Dj in G->D1..Dk-1->Mk wing      */
  float        xJ, xC;		/* used to keep J,C fwd scores from prev row, for decoding JJ/CC   */
  float        xG;		/* for clarity, tmp var, G(i-1) pulled out for wing unfolding      */
  const int    M = gm->M;	/* for clarity, pull out model's size                              */
  int          status;

  /* Contract checks / argument validation */
  ESL_DASSERT1( (afu->type == p7R_ASC_FWD_UP)   );
  ESL_DASSERT1( (afd->type == p7R_ASC_FWD_DOWN) );
  ESL_DASSERT1( (abu->type == p7R_ASC_BCK_UP)   );
  ESL_DASSERT1( (abd->type == p7R_ASC_BCK_DOWN) );
  ESL_DASSERT1( (afu->L == L && afd->L == L && abu->L == L && abd->L == L) );
  ESL_DASSERT1( (afu->M == M && afd->M == M && abu->M == M && abd->M == M) );
  
  /* Reallocation, if needed. 
   * Caller is allowed to overwrite abu -> apu, abd -> apd
   */
  if ( apu != abu && ( status = p7_refmx_GrowTo(apu, M, L)) != eslOK) return status;
  if ( apd != abd && ( status = p7_refmx_GrowTo(apd, M, L)) != eslOK) return status;
  apu->L = apd->L = L;
  apu->M = apd->M = M;
  apu->type = p7R_ASC_DECODE_UP;
  apd->type = p7R_ASC_DECODE_DOWN;


  /* Initialize specials on rows 1..anch[0].i-1 
   * We've above the first anchor, so only S->N->B->LG is possible in specials. 
   * We pick up totsc from row 0 of backwards.
   */
  for (i = 0; i < anch[0].n1; i++)
    {
      fwdp  = afd->dp[i] + (M+1) * p7R_NSCELLS;
      bckp  = abd->dp[i] + (M+1) * p7R_NSCELLS;
      ppp   = apd->dp[i] + (M+1) * p7R_NSCELLS;
      denom = 0.0;

      ppp[p7R_JJ] = 0.0;	
      ppp[p7R_CC] = 0.0; 
      ppp[p7R_E]  = 0.0; 
      ppp[p7R_N]  = (i == 0 ? 1.0 : expf(fwdp[p7R_N] + bckp[p7R_N] - totsc));
      ppp[p7R_J]  = 0.0;  
      ppp[p7R_B]  = expf(fwdp[p7R_B] + bckp[p7R_B] - totsc);
      ppp[p7R_L]  = expf(fwdp[p7R_L] + bckp[p7R_L] - totsc);
      ppp[p7R_G]  = expf(fwdp[p7R_G] + bckp[p7R_G] - totsc); 
      ppp[p7R_C]  = 0.0;

      if (i == 0) totsc = bckp[p7R_N]; 
      else        *(apu->dp[i]) = ppp[p7R_N];  // that's a hack. We stash denom (here, p7R_N) in k=0,ML slot of UP, which is unused
    }


  for (d = 0; d < D; d++)
    {
      /* UP matrix 
       */
      iend = (d = 0 ? 1 : anch[d-1].n1 + 1);
      for (i = iend; i < anch[d].n1; i++)
	{
	  /* Wing retraction of G->D1..Dk-1->MGk paths.
	   * 
	   * In Forward/Backward we used G->MGk directly, but in
           * decoding we need the pp's of the D's. Each Dj in the
           * PREVIOUS row gets an added correction <delta>, which is
           * the sum of all G->Mk paths that run through it, j<k.
           * This step is the only reason we need <rsc>, hence <dsq>,
           * in this function.
	   */
	  bckp  = abu->dp[i]      + (anch[d].n2-1) * p7R_NSCELLS;  // bckp on i, k0-1
	  ppp   = apu->dp[i-1]    + (anch[d].n2-1) * p7R_NSCELLS;  // ppp starts on i+1, k0-1 (PREVIOUS row)
	  rsc   = gm->rsc[dsq[i]] + (anch[d].n2-1) * p7P_NR;
    	  xG    = *(afd->dp[i-1] + (M+1) * p7R_NSCELLS + p7R_G);   // I don't see any good way of avoiding this reach out into memory
	  delta = 0.0;
    
	  for (k = anch[d].n2-1; k >= 1; k--)
	    {
	      ppp[p7R_DG] += delta;
	      delta       += expf(xG + TSC(p7P_GM, k-1) + *rsc + bckp[p7R_MG] - totsc);

	      ppp  -= p7R_NSCELLS;
	      bckp -= p7R_NSCELLS;
	      rsc  -= p7P_NR;
	    }

	  fwdp  = afu->dp[i] + p7R_NSCELLS;
	  bckp  = abu->dp[i] + p7R_NSCELLS;
	  ppp   = apu->dp[i];
	  denom = *ppp;		/* pick up what we stashed earlier; we're now going to finish row i and renormalize if needed */
	  for (s = 0; s < p7R_NSCELLS; s++) *ppp++ = 0.0;

	  /* Main decoding recursion:
	   * [ ML MG IL IG DL DG ] 
	   */
	  for (k = 1; k < anch[d].n2; k++)
	    {
	      ppp[p7R_ML] = expf(fwdp[p7R_ML] + bckp[p7R_ML] - totsc); denom += ppp[p7R_ML];
	      ppp[p7R_MG] = expf(fwdp[p7R_MG] + bckp[p7R_MG] - totsc); denom += ppp[p7R_MG];
	      ppp[p7R_IL] = expf(fwdp[p7R_IL] + bckp[p7R_IL] - totsc); denom += ppp[p7R_IL];
	      ppp[p7R_IG] = expf(fwdp[p7R_IG] + bckp[p7R_IG] - totsc); denom += ppp[p7R_IG];
	      ppp[p7R_DL] = expf(fwdp[p7R_DL] + bckp[p7R_DL] - totsc);
	      ppp[p7R_DG] = expf(fwdp[p7R_DG] + bckp[p7R_DG] - totsc);

	      fwdp += p7R_NSCELLS;
	      bckp += p7R_NSCELLS;
	      ppp  += p7R_NSCELLS;
	    }

#ifndef p7REFERENCE_ASC_DECODING_TESTDRIVE
	  denom = 1.0 / denom;		    // multiplication may be faster than division
	  ppp   = apu->dp[i] + p7R_NSCELLS;                                          // that's k=1 in UP
	  for (s = 0; s < (anch[d].n2-1) * p7R_NSCELLS; s++) *ppp++ *= denom;        // UP matrix row i renormalized
	  ppp   = apd->dp[i] + (anch[d].n2) * p7R_NSCELLS;                           // that's k0 in DOWN
	  for (s = 0; s < (M - anch[d].n2 + 1) * p7R_NSCELLS; s++) *ppp++ *= denom;  // DOWN matrix row i renormalized
	  for (s = 0; s < p7R_NXCELLS; s++) ppp[s] *= denom;
#endif
	} /* end loop over i's in UP chunk for domain d */



      /* DOWN matrix */
      xJ = xC = -eslINFINITY;
      iend = (d == D-1 ? L+1 : anch[d+1].n1);
      for (i = anch[d].n1; i < iend; i++)
	{
	  fwdp  = afd->dp[i] + anch[d].n2 * p7R_NSCELLS;
	  bckp  = abd->dp[i] + anch[d].n2 * p7R_NSCELLS;
	  ppp   = apd->dp[i] + anch[d].n2 * p7R_NSCELLS;
	  denom = 0.0;
	  
	  for (k = anch[d].n2; k <= M; k++)
	    {
	      ppp[p7R_ML] = expf(fwdp[p7R_ML] + bckp[p7R_ML] - totsc); denom += ppp[p7R_ML];
	      ppp[p7R_MG] = expf(fwdp[p7R_MG] + bckp[p7R_MG] - totsc); denom += ppp[p7R_MG];
	      ppp[p7R_IL] = expf(fwdp[p7R_IL] + bckp[p7R_IL] - totsc); denom += ppp[p7R_IL];
	      ppp[p7R_IG] = expf(fwdp[p7R_IG] + bckp[p7R_IG] - totsc); denom += ppp[p7R_IG];
	      ppp[p7R_DL] = expf(fwdp[p7R_DL] + bckp[p7R_DL] - totsc);
	      ppp[p7R_DG] = expf(fwdp[p7R_DG] + bckp[p7R_DG] - totsc);
	    }
	  /* fwdp, bckp, ppp now all sit at M+1, start of specials */

	  ppp[p7R_JJ] = (d == D-1 ? 0.0 : expf(xJ + gm->xsc[p7P_J][p7P_LOOP] + bckp[p7R_J] - totsc)); xJ = fwdp[p7R_J]; denom += ppp[p7R_JJ];
	  ppp[p7R_CC] = (d  < D-1 ? 0.0 : expf(xC + gm->xsc[p7P_C][p7P_LOOP] + bckp[p7R_C] - totsc)); xC = fwdp[p7R_C]; denom += ppp[p7R_CC];
	  ppp[p7R_E]  = expf(fwdp[p7R_E] + bckp[p7R_E] - totsc); 
	  ppp[p7R_N]  = 0.0;
	  ppp[p7R_J]  = (d == D-1 ? 0.0 : expf(fwdp[p7R_J] + bckp[p7R_J] - totsc));
	  ppp[p7R_B]  = expf(fwdp[p7R_B] + bckp[p7R_B] - totsc); 
	  ppp[p7R_L]  = expf(fwdp[p7R_L] + bckp[p7R_L] - totsc);
	  ppp[p7R_G]  = expf(fwdp[p7R_G] + bckp[p7R_G] - totsc);  
	  ppp[p7R_C]  = (d  < D-1 ? 0.0 : expf(fwdp[p7R_C] + bckp[p7R_C] - totsc));

	  
	  if (d < D-1) *(apu->dp[i]) = denom;  // hack: stash denom, which we'll pick up again when we do the next UP matrix.
#ifndef p7REFERENCE_ASC_DECODING_TESTDRIVE
	  if (d == D-1) 
	    { // UP matrices only go through anch[D-1].i-1, so for last DOWN chunk, we need to renormalize.
	      denom = 1.0 / denom;
	      ppp = apd->dp[i] + (anch[d].n2) * p7R_NSCELLS;                           // that's k0 in DOWN
	      for (s = 0; s < (M - anch[d].n2 + 1) * p7R_NSCELLS; s++) *ppp++ *= denom; 
	      for (s = 0; s < p7R_NXCELLS; s++) ppp[s] *= denom;
	    }
#endif
	} /* end loop over i's in DOWN chunk for domain d*/

    } /* end loop over domains d=0..D-1 */
  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 *
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
