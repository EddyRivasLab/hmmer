/* Reference implementation of anchor set constrained (ASC) 
 * posterior decoding.
 * 
 * All reference implementation code is for development and
 * testing. It is not used in HMMER's main executables. Production
 * code uses sparse dynamic programming.
 * 
 * Contents:
 *    1. ASC Decoding
 *    2. Unit tests
 *    3. Test driver
 *    5. Copyright and license information
 */

#include "p7_config.h"

#include "easel.h"

#include "base/p7_profile.h"
#include "base/p7_anchors.h"

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
 *
 * Args:      dsq  : digital target sequence 1..L
 *            L    : length of <dsq>
 *            gm   : profile
 *            anch : array of (i0,k0) anchors defining <dsq>'s domain structure; (0) 1..D (D+1) with sentinels
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
p7_ReferenceASCDecoding(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_ANCHOR *anch, int D, 
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


  /* Initialize specials on rows 1..anch[1].i-1 
   * We've above the first anchor, so only S->N->B->LG is possible in specials. 
   * We pick up totsc from row 0 of backwards.
   */
  for (i = 0; i < anch[1].i0; i++)
    {
      fwdp  = afd->dp[i] + (M+1) * p7R_NSCELLS;
      bckp  = abd->dp[i] + (M+1) * p7R_NSCELLS;
      ppp   = apd->dp[i] + (M+1) * p7R_NSCELLS;
      denom = 0.0;

      if (i == 0) totsc = bckp[p7R_N]; 

      ppp[p7R_JJ] = 0.0;	
      ppp[p7R_CC] = 0.0; 
      ppp[p7R_E]  = 0.0; 
      ppp[p7R_N]  = (i == 0 ? 1.0 : expf(fwdp[p7R_N] + bckp[p7R_N] - totsc));
      ppp[p7R_J]  = 0.0;  
      ppp[p7R_B]  = expf(fwdp[p7R_B] + bckp[p7R_B] - totsc);
      ppp[p7R_L]  = expf(fwdp[p7R_L] + bckp[p7R_L] - totsc);
      ppp[p7R_G]  = expf(fwdp[p7R_G] + bckp[p7R_G] - totsc); 
      ppp[p7R_C]  = 0.0;

      *(apu->dp[i]) = ppp[p7R_N];  // that's a hack. We stash denom (here, p7R_N) in k=0,ML slot of UP, which is unused
    }


  for (d = 1; d <= D; d++)
    {
      /* UP matrix */

      /* Row iend-1 (0, for example, but also the top edge of each UP matrix) is a boundary case:
       * only reachable by wing retraction. Set it all to zero.
       */
      iend = anch[d-1].i0+1;
      for (s = 0; s < anch[d].k0*p7R_NSCELLS; s++) apu->dp[iend-1][s] = 0.;

      for (i = iend; i < anch[d].i0; i++)
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
	  bckp  = abu->dp[i]      + (anch[d].k0-1) * p7R_NSCELLS;  // bckp on i, k0-1
	  ppp   = apu->dp[i-1]    + (anch[d].k0-1) * p7R_NSCELLS;  // ppp starts on i+1, k0-1 (PREVIOUS row)
	  rsc   = gm->rsc[dsq[i]] + (anch[d].k0-1) * p7P_NR;
    	  xG    = *(afd->dp[i-1] + (M+1) * p7R_NSCELLS + p7R_G);   // I don't see any good way of avoiding this reach out into memory
	  delta = 0.0;
    
	  for (k = anch[d].k0-1; k >= 1; k--)
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
	  denom = *ppp;	 /* pick up what we stashed earlier; we're now going to finish row i and renormalize if needed */
	  for (s = 0; s < p7R_NSCELLS; s++) *ppp++ = 0.0;

	  /* Main decoding recursion:
	   * [ ML MG IL IG DL DG ] 
	   */
	  for (k = 1; k < anch[d].k0; k++)
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

	  //#ifndef p7REFERENCE_ASC_DECODING_TESTDRIVE
	  denom = 1.0 / denom;		                                             // multiplication may be faster than division
	  ppp   = apu->dp[i] + p7R_NSCELLS;                                          // that's k=1 in UP
	  for (s = 0; s < (anch[d].k0-1) * p7R_NSCELLS; s++) *ppp++ *= denom;        // UP matrix row i renormalized
	  ppp   = apd->dp[i] + (anch[d].k0) * p7R_NSCELLS;                           // that's k0 in DOWN
	  for (s = 0; s < (M - anch[d].k0 + 1) * p7R_NSCELLS; s++) *ppp++ *= denom;  // DOWN matrix row i renormalized
	  for (s = 0; s < p7R_NXCELLS; s++) ppp[s] *= denom;
	  //#endif
	} /* end loop over i's in UP chunk for domain d */

      /* One last wing-retracted row:
       * On last row of each UP chunk (anch[d].i-1 row), those DG's
       * are reachable on a G(i)->D1...Dk-1->MG(i+1,k) path that we
       * need to unfold; and that MG(i+1,k) is the anchor itself,
       * whose value is in the DOWN matrix.
       */
      i     = anch[d].i0;	
      k     = anch[d].k0;
      xG    = *(afd->dp[i-1] + (M+1) * p7R_NSCELLS + p7R_G); // ugly, sorry. Reach out and get the xG value on row i.
      bckp  = abd->dp[i]      + k * p7R_NSCELLS;             // *bckp is the anchor cell i0,k0 (in DOWN)
      rsc   = gm->rsc[dsq[i]] + k * p7P_NR;
      delta = expf(xG + TSC(p7P_GM, k-1) + *rsc + bckp[p7R_MG] - totsc);  // k-1 in TSC because GMk is stored off by one, in -1 slot
      ppp   = apu->dp[i-1] + p7R_NSCELLS;                    // ppp starts on i-1,k=1 cells
      for (k = 1; k < anch[d].k0; k++) 
	{
	  ppp[p7R_DG] += delta;   // unlike the main wing retraction above, there's only one path going through these 
	  ppp += p7R_NSCELLS;     //   DGk's (the G->Mk0 entry into the anchor cell), so there's no retro-accumulation of delta
	}
      

      /* DOWN matrix */
      xJ = xC = -eslINFINITY;
      for (i = anch[d].i0; i < anch[d+1].i0; i++)
	{
	  fwdp  = afd->dp[i] + anch[d].k0 * p7R_NSCELLS;
	  bckp  = abd->dp[i] + anch[d].k0 * p7R_NSCELLS;
	  ppp   = apd->dp[i] + anch[d].k0 * p7R_NSCELLS;
	  denom = 0.0;
	  
	  for (k = anch[d].k0; k <= M; k++)
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
	  /* fwdp, bckp, ppp now all sit at M+1, start of specials */

	  ppp[p7R_JJ] = (d == D ? 0.0 : expf(xJ + gm->xsc[p7P_J][p7P_LOOP] + bckp[p7R_J] - totsc)); xJ = fwdp[p7R_J]; denom += ppp[p7R_JJ];
	  ppp[p7R_CC] = (d  < D ? 0.0 : expf(xC + gm->xsc[p7P_C][p7P_LOOP] + bckp[p7R_C] - totsc)); xC = fwdp[p7R_C]; denom += ppp[p7R_CC];
	  ppp[p7R_E]  = expf(fwdp[p7R_E] + bckp[p7R_E] - totsc); 
	  ppp[p7R_N]  = 0.0;
	  ppp[p7R_J]  = (d == D ? 0.0 : expf(fwdp[p7R_J] + bckp[p7R_J] - totsc));
	  ppp[p7R_B]  = expf(fwdp[p7R_B] + bckp[p7R_B] - totsc); 
	  ppp[p7R_L]  = expf(fwdp[p7R_L] + bckp[p7R_L] - totsc);
	  ppp[p7R_G]  = expf(fwdp[p7R_G] + bckp[p7R_G] - totsc);  
	  ppp[p7R_C]  = (d  < D ? 0.0 : expf(fwdp[p7R_C] + bckp[p7R_C] - totsc));

	  
	  /* Even with forced renormalization, I don't think you can
	   * do much better than the error inherent in the default
	   * FLogsum() calculation. For example, the calculation of
	   * ppp[C] above can give a pp of > 1.0, even if the sum of
	   * all emitters on the row is denom=1.0, because of FLogsum
	   * table approximation error in fwdp[C].
	   */
	  if (d < D) *(apu->dp[i]) = denom;  // hack: stash denom, which we'll pick up again when we do the next UP matrix.
	  //#ifndef p7REFERENCE_ASC_DECODING_TESTDRIVE
	  if (d == D) 
	    { // UP matrices only go through anch[D].i-1, so for last DOWN chunk, we need to renormalize.
	      denom = 1.0 / denom;
	      ppp = apd->dp[i] + (anch[d].k0) * p7R_NSCELLS;                           // that's k0 in DOWN
	      for (s = 0; s < (M - anch[d].k0 + 1) * p7R_NSCELLS; s++) *ppp++ *= denom; 
	      for (s = 0; s < p7R_NXCELLS; s++) ppp[s] *= denom;
	    }
	  //#endif
	} /* end loop over i's in DOWN chunk for domain d*/

    } /* end loop over domains d=0..D-1 */
  return eslOK;
}
/*------------------ end, ASC decoding --------------------------*/


/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef p7REFERENCE_ASC_DECODING_TESTDRIVE
#include "hmmer.h"

/* The unit tests are essentially the same as the unit tests for
 * Forward/Backward, but carried through to the additional step of
 * posterior decoding. The general idea is to create contrived
 * profile/seq/anchorset comparisons where we can compare standard DP
 * matrices and/or results to ASC matrices and/or results. 
 * 
 * Most of the tests (all but "singlemulti") create comparisons where
 * the standard and the ASC calculations have identical path
 * ensembles, by forcing the standard DP calculation to use the anchor
 * cells even without constraint by the ASC algorithm. Forcing this
 * requires unusual models, but there are various ways to do it, and
 * different ways allows us to test different aspects of the model, as
 * summarized below (also see SRE:J13/59):
 *                                                                     same   
 *                      multiple   N/C/J     I             multiple   std,ASC
 *  Unit test            paths?    emits?  emits?  Local?  domains?   ensemble?
 *  ---------------     --------  -------  ------  ------  --------   ---------
 *  singlepath              -        -      1/0       -        -        YES
 *  singlesingle            -       YES     YES       -        -        YES
 *  singlemulti             -       YES     YES       -       YES         -
 *  multisingle           YES       YES     YES       -        -        YES
 *  multipath_local       YES        -       -       YES       -          -
 *  multimulti            YES        -       -        -       YES         -  
 * 
 * Single pathed tests verify that the decoding matrix has pp=1.0
 * along the trace, 0.0 elsewhere. (The singlemulti test is special:
 * the ASC Decoding matrix has only one nonzero path, the trace, but
 * standard DP has an ensemble.) Ensemble tests compare standard to
 * ASC matrices. 
 * 
 * When comparing ASC matrices to standard Decoding matrices
 * cell-by-cell, if there is an ensemble, we have to marginalize
 * UP+DOWN before doing the comparison.
 */




/* ascmatrix_compare_std()
 * Compare a standard Decoding matrix to ASC UP/DOWN matrices,
 * cell by cell, for equality.
 * 
 * This version of ascmatrix_compare() is similar to the one in 
 * reference_asc_fwdback.c, but is specialized to decoding matrices. 
 * 
 * Key observation: For all but the singlemulti test, standard and
 * ASC algorithms recursively enumerate exactly the same set of paths,
 * so we can directly compare the ASC calculations to the standard
 * calculations, even at the level of individual DP matrix cells. But to
 * do this comparison, ASC Decoding values must be marginalized over
 * UP/DOWN sectors (sum UP + DOWN values) before comparing to a
 * standard DP cell value.
 * 
 * Return <eslOK> if all valid values in the standard Decoding DP
 * matrix <rxd> matches the values in the ASC UP/DOWN matrices
 * <apu>/<apd>; <eslFAIL> if they don't.
 */
static int
ascmatrix_compare_std(P7_REFMX *rxd, P7_REFMX *apu, P7_REFMX *apd, P7_ANCHOR *anch, int D, float epsilon)
{
  int M = rxd->M;
  int L = rxd->L;
  int d, i, k, s;
  float ascval;
  int killmenow = FALSE;
#ifdef p7_DEBUGGING
  killmenow = TRUE;
#endif

  /* contract check, argument validation */
  ESL_DASSERT1( (apu->M == M && apu->L == L));
  ESL_DASSERT1( (apd->M == M && apd->L == L));
  ESL_DASSERT1( (rxd->type == p7R_DECODING  && apu->type == p7R_ASC_DECODE_UP && apd->type == p7R_ASC_DECODE_DOWN) );

  /* We use an access pattern that lets us check UP and DOWN sectors
   * for the same i,k,s,d cell at the same time
   */
  for (d = 1, i = 0; i <= L; i++)
    {      
      if (i == anch[d].i0) d++;            // d = index of current UP matrix. d-1 = DOWN matrix;  d will be D+1 for final chunk
      for (k = 1; k <= M; k++)             //   ... other way to think about it: d = index of next anchor we will reach
	for (s = 0; s < p7R_NSCELLS; s++)
	  {
	    ascval = 0.0;
	    if (k >= anch[d-1].k0) ascval += P7R_MX(apd,i,k,s);   // sentinel k0(0)   = M+1, so no k gets evaluated for DOWN(d-1=0) 
	    if (k <  anch[d].k0)   ascval += P7R_MX(apu,i,k,s);   // sentinel k0(D+1) = 0,   so no k gets evaluated for UP(d=D+1)
	    
	    if (esl_FCompareAbs( ascval, P7R_MX(rxd,i,k,s), epsilon) != eslOK)
	      { if (killmenow) abort(); return eslFAIL; }
	  }
      for (s = 0; s < p7R_NXCELLS; s++)
	if (esl_FCompareAbs( P7R_XMX(apd,i,s), P7R_XMX(rxd,i,s), epsilon) != eslOK)
	  { if (killmenow) abort(); return eslFAIL; }
    }
  return eslOK;
}

/* Same as above, except for use with --diag: collect the largest
 * absolute difference seen for any cell in the compared matrices;
 * i.e. maxdiff = max fabs(diff(i,k,s))
 */
static float 
ascmatrix_maxdiff_std(P7_REFMX *rxd, P7_REFMX *apu, P7_REFMX *apd, P7_ANCHOR *anch, int D)
{
  int   d,i,k,s;
  float ascval;
  float diff;
  float maxdiff = 0.;

  for (d = 1, i = 0; i <= rxd->L; i++)
    {
      if (i == anch[d].i0) d++;
      for (k = 1; k <= rxd->M; k++)
	for (s = 0; s < p7R_NSCELLS; s++)
	  {
	    ascval = 0.0;	
	    if (k >= anch[d-1].k0) ascval += P7R_MX(apd,i,k,s); 
	    if (k <  anch[d].k0)   ascval += P7R_MX(apu,i,k,s); 

	    diff = ascval - P7R_MX(rxd,i,k,s);  
	    if (fabs(diff) > fabs(maxdiff)) maxdiff = diff;
	  }
      for (s = 0; s < p7R_NXCELLS; s++)
	{
	  diff = P7R_XMX(apd,i,s) - P7R_XMX(rxd,i,s);
	  if (fabs(diff) > fabs(maxdiff)) maxdiff = diff;
	}
    }
  return maxdiff;
}



/* ascmatrix_compare_asc()
 */
static int
ascmatrix_compare_asc(P7_REFMX *apu1, P7_REFMX *apd1, P7_REFMX *apu2, P7_REFMX *apd2, P7_ANCHOR *anch, int D, float epsilon)
{
  int M = apu1->M;
  int L = apu1->L;
  int d, i, k, s;
  int killmenow = FALSE;
#ifdef p7_DEBUGGING
  killmenow = TRUE;
#endif

  /* contract check, argument validation */
  ESL_DASSERT1(( apu1->M == apd1->M && apu1->M == apu2->M && apu1->M == apd2->M));
  ESL_DASSERT1(( apu1->L == apd1->L && apu1->L == apu2->L && apu1->L == apd2->L));
  ESL_DASSERT1(( apu1->type == p7R_ASC_DECODE_UP && apd1->type == p7R_ASC_DECODE_DOWN ));
  ESL_DASSERT1(( apu2->type == p7R_ASC_DECODE_UP && apd2->type == p7R_ASC_DECODE_DOWN ));

  for (d = 1, i = 0; i <= L; i++)
    {
      if (i == anch[d].i0) d++;             // d=1..D+1, with D+1 in final leg; think of d as "index of next anchor we will reach"

      /* DOWN row, if one exists for this i */
      for (k = anch[d-1].k0; k <= M; k++)   // sentinel k0(0) = M+1, so no k gets evaluated at d=1 for nonexistent DOWN(0)
	for (s = 0; s < p7R_NSCELLS; s++)   //   ... i.e., first leg has only an UP(1) matrix.
	  if (esl_FCompareAbs( P7R_MX(apd1,i,k,s), P7R_MX(apd2,i,k,s), epsilon) != eslOK)
	    { if (killmenow) abort(); return eslFAIL; }

      /* specials exist for all rows. */
      for (s = 0; s < p7R_NXCELLS; s++)
	if (esl_FCompareAbs( P7R_XMX(apd1,i,s), P7R_XMX(apd2,i,s), epsilon) != eslOK)
	  { if (killmenow) abort(); return eslFAIL; }

      /* UP row, if one exists for this i */
      for (k = 1; k < anch[d].k0; k++)   // sentinel k0(D+1) = 0, so no k gets evaluated at d=D+1 for nonexistent UP(D+1)
	for (s = 0; s < p7R_NSCELLS; s++)
	  if (esl_FCompareAbs( P7R_MX(apu1,i,k,s), P7R_MX(apu2,i,k,s), epsilon) != eslOK)
	    { if (killmenow) abort(); return eslFAIL; }
    }
  return eslOK;
}

/* refmx_trace_embed()
 * 
 * Create a posterior decoding matrix for a single path <tr>. Return
 * the new standard decoding matrix in <ret_ppt>. This matrix has values
 * of 0.0 or 1.0.
 * 
 * Used by both ascmatrix_trace_compare() and ascmatrix_trace_maxdiff().
 */
static int
refmx_trace_embed(P7_TRACE *tr, P7_REFMX **ret_ppt)
{
  static int sttype[p7T_NSTATETYPES] = { -1, p7R_ML, p7R_MG, p7R_IL, p7R_IG, p7R_DL, p7R_DG, -1, p7R_N, p7R_B, p7R_L, p7R_G, p7R_E, p7R_C, p7R_J, -1 }; /* sttype[] translates trace idx to DP matrix idx*/
  P7_REFMX *ppt = p7_refmx_Create(tr->M, tr->L);
  int       i,z;

  p7_refmx_SetType  (ppt, tr->M, tr->L, p7R_DECODING);
  p7_refmx_SetValues(ppt, 0.0);
  for (i = 0, z = 1; z < tr->N-1; z++)
    {
      if (tr->i[z]) i = tr->i[z];
      if (p7_trace_IsMain(tr->st[z]))  
	P7R_MX (ppt, i, tr->k[z], sttype[(int)tr->st[z]]) = 1.0;
      else {
	P7R_XMX(ppt, i, sttype[(int)tr->st[z]]) = 1.0;
	if (tr->st[z] == p7T_C && tr->i[z]) P7R_XMX(ppt, i, p7R_CC) = 1.0;
	if (tr->st[z] == p7T_J && tr->i[z]) P7R_XMX(ppt, i, p7R_JJ) = 1.0;
      }
    }
  *ret_ppt = ppt;
  return eslOK;
}


/* ascmatrix_trace_compare()
 * 
 * For a single pathed test, the cells in a posterior decoding matrix
 * have values 1.0 or 0.0; specifically, 1.0 wherever the single path
 * visits a cell. Test this, given the single path <tr> and
 * ASC posterior decoding matrix <apu/apd>.
 * 
 * Returns <eslOK> if the comparison is correct; <eslFAIL> if not.
 * 
 * Shares features w/ p7_refmx_CountTrace(). For example, it's
 * easy to mess up the conversion from trace state types (p7T types)
 * to matrix cell types (p7R types), and easy to forget to do the
 * CC/JJ cells in the decoding matrix.
 */
static int
ascmatrix_trace_compare(P7_TRACE *tr, P7_REFMX *apu, P7_REFMX *apd, P7_ANCHOR *anch, int D, float epsilon)
{
  P7_REFMX *ppt  = NULL;
  int       M    = apu->M;
  int       L    = apu->L;
  int       status;

  /* Contract checks */
  ESL_DASSERT1(( apd->M == M && tr->M == M));
  ESL_DASSERT1(( apd->L == L && tr->L == L));
  ESL_DASSERT1(( apu->type == p7R_ASC_DECODE_UP && apd->type == p7R_ASC_DECODE_DOWN ));

  refmx_trace_embed(tr, &ppt);

  //printf("### Trace Matrix:\n"); p7_refmx_Dump(stdout, ppt);
      
  /* Now we can just compare this std DP matrix to the ASC matrices */
  status = ascmatrix_compare_std(ppt, apu, apd, anch, D, epsilon);
  
  p7_refmx_Destroy(ppt);
  return status;
}

/* ascmatrix_trace_maxdiff()
 * 
 * Same as above, but for --diag: return the maximum absolute
 * difference between the actual values in <apu/apd> and the expected
 * 0.0/1.0 for decoded cells in a single pathed test, for the path
 * <tr>.
 */
static float 
ascmatrix_trace_maxdiff(P7_TRACE *tr, P7_REFMX *apu, P7_REFMX *apd, P7_ANCHOR *anch, int D)
{
  P7_REFMX *ppt  = NULL;
  float     maxdiff;

  refmx_trace_embed(tr, &ppt);
  maxdiff = ascmatrix_maxdiff_std(ppt, apu, apd, anch, D);

  p7_refmx_Destroy(ppt);
  return maxdiff;
}

/* ascmatrix_validate()
 *
 * Validate an ASC Decoding matrix UP/DOWN pair.  Testing here is less
 * thorough than p7_refmx_Validate().  Here we're just looking for
 * obviously bad values: values that are not probabilities.
 *
 * Throws <eslFAIL> if validation fails.
 */
static int
ascmatrix_validate(P7_REFMX *apu, P7_REFMX *apd, P7_ANCHOR *anch, int D, float epsilon, char *errbuf)
{
  int L = apu->L;
  int M = apu->M;
  int d,i,k,s;
  float val;
  
  for (d = 1, i = 0; i <= L; i++)
    {
      if (i == anch[d].i0) d++;

      /* DOWN sector for domain d, if one exists for this i */
      for (k = anch[d-1].k0; k <= M; k++)
	for (s = 0; s < p7R_NSCELLS; s++)
	  {
	    val = P7R_MX(apd,i,k,s);
	    if (isnan(val) || val < 0. || val > 1.0+epsilon) 
	      ESL_FAIL(eslFAIL, errbuf, "bad DOWN value %f at i=%d, k=%d, %s", val, i, k, p7_refmx_DecodeState(s));
	  }	

      /* specials exist on all rows, stored in <apd> */
      for (s = 0; s < p7R_NXCELLS; s++)
	{
	  val = P7R_XMX(apd,i,s);
	  if (isnan(val) || val < 0. || val > 1.0+epsilon) 
	    ESL_FAIL(eslFAIL, errbuf, "bad special value at i=%d, %s", i, p7_refmx_DecodeSpecial(s));
	}	

      /* UP sector for domain d */
      for (k = 1; k < anch[d].k0; k++)
	for (s = 0; s < p7R_NSCELLS; s++)
	  {
	    val = P7R_MX(apu,i,k,s);
	    if (isnan(val) || val < 0. || val > 1.0+epsilon) 
	      ESL_FAIL(eslFAIL, errbuf, "bad UP value %f at i=%d, k=%d, %s", val, i, k, p7_refmx_DecodeState(s));
	    }	
    }
  return eslOK;
}

/* ascmatrix_max_over_one()
 * 
 * Like ascmatrix_validate(), except for --diag: look for values
 * in the decoding matrix > 1.0, and return the maximum excursion
 * seen (i.e. max(val-1.0) for val>1.0).
 */
static float
ascmatrix_max_over_one(P7_REFMX *apu, P7_REFMX *apd, P7_ANCHOR *anch, int D)
{
  int d, i, k, s;
  float max = 0.;
  
  for (d = 1, i = 0; i <= apu->L; i++)
    {
      if (i == anch[d].i0) d++;

      for (k = anch[d-1].k0; k <= apu->M; k++)
	for (s = 0; s < p7R_NSCELLS; s++)
	  if (P7R_MX(apd,i,k,s) > max) max = P7R_MX(apd,i,k,s);

      for (s = 0; s < p7R_NXCELLS; s++)
	if (P7R_XMX(apd,i,s) > max) max = P7R_XMX(apd,i,s);

      for (k = 1; k < anch[d].k0; k++)
	for (s = 0; s < p7R_NSCELLS; s++)
	  if (P7R_MX(apu,i,k,s) > max) max = P7R_MX(apu,i,k,s);
    }
  
  return (max > 1.0 ? max-1.0 : 0.0);
}


/* "overwrite" test
 * 
 * Same idea as in reference_decoding.c. Verify an important feature in 
 * the API: we're allowed to overwrite the input Backwards matrices with
 * the new decoding matrix, thus saving matrix allocations. But this
 * can be tricky to get right in implementing Decoding. 
 * 
 * The unit tests samples a random profile/sequence comparison,
 * decodes it with and without overwriting, and verifies that the
 * resulting ASC decoding matrices are valid and exactly identical.
 * 
 *****************************************************************
 * May fail stochastically. SRE:J11/128 27-Jun-14
 *  ./reference_asc_decoding_utest --diag overwrite -N 10000
 * Field 1: max_over_one
 * 
 * Though the main test is for exact equality of decoding matrices
 * produced in the two alternative call styles, the test also 
 * calls Validate() on the decoding matrices, and a validate call
 * checks for values > 1.0. This can happen because of roundoff
 * and FLogsum() error.
 * 
 * Default: mean 2e-5 sd 4e-5 max 0.0004 => 0.01
 * Exact:   mean 9e-8 sd 1e-7 max 1e-6   => 0.0001
 */
static void
utest_overwrite(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding overwrite unit test failed";
  ESL_SQ     *sq        = esl_sq_CreateDigital(abc);
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = p7_profile_Create(M, abc);
  P7_ANCHORS *anch      = p7_anchors_Create();
  P7_TRACE   *vtr       = p7_trace_Create();
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxd       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  P7_REFMX   *apu       = p7_refmx_Create(M, 20);
  P7_REFMX   *apd       = p7_refmx_Create(M, 20);
  float       tolerance = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  char        errbuf[eslERRBUFSIZE];

  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(failmsg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(failmsg);
  
  /* Emit sequence from model, using a length model of <L=M>;
   * restrict the emitted sequence length to 6M, arbitrarily, to 
   * keep it down to something reasonable.
   */
  if ( p7_profile_SetLength(gm, M)    != eslOK) esl_fatal(failmsg);
  do {
    esl_sq_Reuse(sq);
    if (p7_ProfileEmit(rng, hmm, gm, bg, sq, /*tr=*/NULL) != eslOK) esl_fatal(failmsg);
  } while (sq->n > M*6);
  if (p7_profile_SetLength(gm, sq->n) != eslOK) esl_fatal(failmsg);

  /* We need an anchor set. It doesn't have to be good, just valid. 
   * To use SetFromTrace(), though, we need a standard Decoding matrix,
   * and we may as well get a Viterbi trace too.
   */
  if ( p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxv, vtr, NULL) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      NULL) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceBackward(sq->dsq, sq->n, gm, rxd,      NULL) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxd,  rxd) != eslOK) esl_fatal(failmsg);
  if ( p7_reference_anchors_SetFromTrace(rxd, vtr, anch)        != eslOK) esl_fatal(failmsg);

  /* F/B, then decode both ways */
  if ( p7_ReferenceASCForward (sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, NULL)               != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch->a, anch->D, abu, abd, NULL)               != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, abu, abd, apu, apd) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch->a, anch->D, afu, afd, abu, abd, abu, abd) != eslOK) esl_fatal(failmsg);

  //printf("### no overwrite: ASC Up:\n");     p7_refmx_Dump(stdout, apu);
  //printf("### no overwrite: ASC Down:\n");   p7_refmx_Dump(stdout, apd);
  //printf("### with overwrite: ASC Up:\n");   p7_refmx_Dump(stdout, abu);
  //printf("### with overwrite: ASC Down:\n"); p7_refmx_Dump(stdout, abd);

  /* Now both apu/apd and abu/abd should contain valid and equal ASC Decoding matrices */
  if (!diagfp) 
    {
      if ( ascmatrix_validate(apu, apd, anch->a, anch->D, tolerance, errbuf)  != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
      if ( ascmatrix_validate(abu, abd, anch->a, anch->D, tolerance, errbuf)  != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
      if ( ascmatrix_compare_asc(apu, apd, abu, abd, anch->a, anch->D, 0.0)   != eslOK) esl_fatal(failmsg);  // expect exact match, so tolerance = 0.
    }
  else
    fprintf(diagfp, "%20g\n", ascmatrix_max_over_one(apu, apd, anch->a, anch->D));

  p7_anchors_Destroy(anch);
  p7_trace_Destroy(vtr);
  p7_refmx_Destroy(apu); p7_refmx_Destroy(apd);
  p7_refmx_Destroy(abu); p7_refmx_Destroy(abd);
  p7_refmx_Destroy(afu); p7_refmx_Destroy(afd);
  p7_refmx_Destroy(rxv); p7_refmx_Destroy(rxf); p7_refmx_Destroy(rxd);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_sq_Destroy(sq);
}


/* "singlepath" test.
 * 
 * p7_modelsample_SinglePathed() gives us a profile HMM contrived
 * such that there is only one possible path; the profile must be
 * uniglocal L=0 to do this, and thus lacks N/C/J emissions, local
 * alignment, and multihit.
 * 
 * Since only one path is possible, every cell along this path decodes
 * to pp=1.0, and all other cells to 0.0. 
 * 
 * We test:
 *    1. All scores are equal (V, F, B, ASC F, ASC B)
 *    2. Vit trace = generated trace
 *    3. ASC Decoding matrix has 1.0 for all cells in trace
 *    4. ASC Decoding matrix = standard Decoding mx, cell by cell.
 *    
 *****************************************************************
 * May fail stochastically. Analysis: SRE: J13/126.
 *    ./reference_asc_decoding_utest --diag singlepath -N 10000
 * 1: max_over_one  (in Validation, how much can a cell's P > 1.0)
 * 2: trace_maxdiff (in trace_compare, how much can a cell differ)
 * 3: maxdiff_std   (in compare_std, how much can a cell differ)
 * 
 * We expect all of these to be exactly zero. The singlepath test
 * creates a model that must generate a unique sequence with 
 * P=1, and that sets all emission/transition probs to 1.0 or 0.0,
 * which can be exactly represented both in probspace and logspace.
 *
 * For tests that duplicate fwdback_utest, see analysis over there.
 * 
 * Both default and exact logsum do indeed show all zero difference. 
 */
static void
utest_singlepath(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding singlepath unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = p7_profile_Create(M, bg->abc);
  ESL_SQ     *sq        = esl_sq_CreateDigital(bg->abc);
  P7_TRACE   *tr        = p7_trace_Create();
  P7_TRACE   *vtr       = p7_trace_Create();
  P7_ANCHOR  *anch      = malloc(sizeof(P7_ANCHOR)* 3);   // utest is uniglocal, so it will have D=1 and only one anchor; +2 for sentinels
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);          // 20 is arbitrary; all DP routines will reallocate as needed.
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxd       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  P7_REFMX   *apu       = p7_refmx_Create(M, 20);
  P7_REFMX   *apd       = p7_refmx_Create(M, 20);
  int         D;
  int         nm, z, which;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       fbtol     = 0.0001;                // floating pt tolerance on fwdback tests
  float       dectol    = 0.;                    // zero floating pt tolerance (demand exact equality) on decoding tests
  char        errbuf[eslERRBUFSIZE];

  /* Sample single-pathed HMM, create uniglocal L=0 profile */
  if ( p7_modelsample_SinglePathed(rng, M, bg->abc, &hmm) != eslOK) esl_fatal(failmsg);
  if ( p7_profile_ConfigUniglocal(gm, hmm, bg, 0)         != eslOK) esl_fatal(failmsg); 

  /* Emit sequence and trace. Get the trace score. */
  if (p7_ProfileEmit(rng, hmm, gm, bg, sq, tr) != eslOK) esl_fatal(failmsg);
  if (p7_trace_Index(tr)                       != eslOK) esl_fatal(failmsg);
  if (p7_trace_Score(tr, sq->dsq, gm, &sc)     != eslOK) esl_fatal(failmsg);

  /* Create anchor "set". We know it's uniglocal, so D=1. 
   * Choose anchor by randomly choosing a match state in the trace.
   */
  D  = 1;
  for (nm = 0, z = 1; z < tr->N; z++) 
    if (p7_trace_IsM(tr->st[z])) nm++;
  which = 1 + esl_rnd_Roll(rng, nm);
  for (z = 1; z < tr->N; z++) 
    if (p7_trace_IsM(tr->st[z])) 
      {
	which--;
	if (which == 0) {
	  anch[D].i0 = tr->i[z];
	  anch[D].k0 = tr->k[z];
	  break;
	}
      }
  p7_anchor_SetSentinels(anch, D, tr->L, tr->M);

  /* Run the DP routines, both standard and ASC */
  if ( p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxv, vtr, &vsc) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceBackward(sq->dsq, sq->n, gm, rxb,      &bsc) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxb,  rxd) != eslOK) esl_fatal(failmsg);

  if ( p7_ReferenceASCForward (sq->dsq, sq->n, gm, anch, D, afu, afd, &asc_f)             != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch, D, abu, abd, &asc_b)             != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch, D, afu, afd, abu, abd, apu, apd) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Vit:\n");      p7_refmx_Dump(stdout, rxv);
  //printf("### Reference Fwd:\n");      p7_refmx_Dump(stdout, rxf);
  //printf("### Reference Bck:\n");      p7_refmx_Dump(stdout, rxb);
  //printf("### Reference Decoding:\n"); p7_refmx_Dump(stdout, rxd);
  //printf("### ASC Fwd UP:\n");         p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");       p7_refmx_Dump(stdout, afd);
  //printf("### ASC Bck UP:\n");         p7_refmx_Dump(stdout, abu);
  //printf("### ASC Bck DOWN:\n");       p7_refmx_Dump(stdout, abd);
  //printf("### ASC Decode UP:\n");      p7_refmx_Dump(stdout, apu);
  //printf("### ASC Decode DOWN:\n");    p7_refmx_Dump(stdout, apd);
  //p7_trace_DumpAnnotated(stdout, tr,  gm, sq->dsq);
  //p7_trace_DumpAnnotated(stdout, vtr, gm, sq->dsq);

  /* Tests.
   * Some tests are redundant with what we do in asc_fwdback,
   * but that's ok.
   */
  if (!diagfp)
    {
      if (esl_FCompareAbs(sc, vsc,   fbtol) != eslOK) esl_fatal(failmsg);  // generated trace score = Viterbi score
      if (esl_FCompareAbs(sc, fsc,   fbtol) != eslOK) esl_fatal(failmsg);  //  ... = Forward score
      if (esl_FCompareAbs(sc, bsc,   fbtol) != eslOK) esl_fatal(failmsg);  //  ... = Backward score
      if (esl_FCompareAbs(sc, asc_f, fbtol) != eslOK) esl_fatal(failmsg);  //  ... = ASC Forward score
      if (esl_FCompareAbs(sc, asc_b, fbtol) != eslOK) esl_fatal(failmsg);  //  ... = ASC Backward score
      if (p7_trace_Compare(tr, vtr, 0)      != eslOK) esl_fatal(failmsg);  // generated trace = Viterbi trace

      if (ascmatrix_validate(apu, apd, anch, D, dectol, errbuf)  != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
      if (ascmatrix_trace_compare(tr, apu, apd, anch, D, dectol) != eslOK) esl_fatal(failmsg);
      if (ascmatrix_compare_std (rxd, apu, apd, anch, D, dectol) != eslOK) esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g %20g\n",
	    ascmatrix_max_over_one(apu, apd, anch, D),
	    ascmatrix_trace_maxdiff(tr, apu, apd, anch, D),
	    ascmatrix_maxdiff_std(rxd, apu, apd, anch, D));

  p7_refmx_Destroy(apu);   p7_refmx_Destroy(apd);
  p7_refmx_Destroy(abu);   p7_refmx_Destroy(abd);
  p7_refmx_Destroy(afu);   p7_refmx_Destroy(afd);
  p7_refmx_Destroy(rxv);
  p7_refmx_Destroy(rxf);
  p7_refmx_Destroy(rxb);
  p7_refmx_Destroy(rxd);
  p7_trace_Destroy(tr);
  p7_trace_Destroy(vtr);
  free(anch);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
}


/* "singlesingle" test
 *
 * The "singlesingle" test (for single path, single domain) uses
 * p7_modelsample_SinglePathedSeq() to sample a contrived profile/seq
 * pair that has only one possible path with probability 1 (that is,
 * P(\pi | x, H) = 1) for this sequence. This must be a glocal single domain
 * path. Choose any anchor A in that domain. 
 * 
 * Then:
 *   1. Viterbi = Fwd = Bck = ASC Fwd = ASC Bck scores.
 *   2. Viterbi trace = trace that generated the test sequence. 
 *   3. ASC Decoding matrix has 1.0 for all cells in trace
 *   4. ASC Decoding matrix = standard Decoding mx, cell by cell.
 *
 *****************************************************************
 * May fail stochastically. Analysis: SRE:J13/126 26-Jun-14
 *  ./reference_asc_decoding_utest --diag singlesingle -N 10000
 * 1: max_over_one  (in Validation, how much can a cell's P > 1.0)
 * 2: trace_maxdiff (in trace_compare, how much can a cell differ)
 * 3: maxdiff_std   (in compare_std, how much can a cell differ)
 * 
 * Unlike the singlepath test, now we have a model with various
 * emission/transition probabilities, subject to roundoff error.
 * maxdiff_std is exactly zero because order of operations is
 * identical (I believe), but I haven't proven that that is
 * platform-independent. 
 * 
 * logsum error has no effect in a single path test, so we set
 * same tolerance for default and exact logsum cases.
 *
 * For tests that duplicate fwdback_utest, see analysis over there,
 * which determined tolerance of 0.001.
 *
 * Default:
 *   max_over_one:  mean 2e-7, sd 7e-7, range up to 2.2e-5  => 0.0001
 *   trace_maxdiff: mean 5e-8, sd 1e-6, range -3e-5 .. 2e-5 => 0.0001
 *   maxdiff_std:   zero
 * Exact is similar.  
 */
static void
utest_singlesingle(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding singlesingle unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_TRACE   *vtr       = p7_trace_Create();
  P7_ANCHOR  *anch      = NULL;
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxb       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxd       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  P7_REFMX   *apu       = p7_refmx_Create(M, 20);
  P7_REFMX   *apd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, vsc, fsc, bsc, asc_f, asc_b;
  float       fbtol      = 0.001;
  float       dectol     = 0.0001;
  char        errbuf[eslERRBUFSIZE];

  if ( p7_modelsample_SinglePathedSeq(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc) != eslOK) esl_fatal(failmsg);

  if ( p7_ReferenceViterbi (dsq, L, gm, rxv, vtr, &vsc) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceForward (dsq, L, gm, rxf,      &fsc) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceBackward(dsq, L, gm, rxb,      &bsc) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceDecoding(dsq, L, gm, rxf, rxb,  rxd) != eslOK) esl_fatal(failmsg);

  if ( p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd, &asc_f)             != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd, &asc_b)             != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, abu, abd, apu, apd) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Vit:\n");   p7_refmx_Dump(stdout, rxv);
  //printf("### Reference Fwd:\n");   p7_refmx_Dump(stdout, rxf);
  //printf("### ASC Fwd UP:\n");      p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");    p7_refmx_Dump(stdout, afd);
  //printf("### ASC Decode UP:\n");   p7_refmx_Dump(stdout, apu);
  //printf("### ASC Decode DOWN:\n"); p7_refmx_Dump(stdout, apd);
  //p7_trace_DumpAnnotated(stdout, tr,  gm, dsq);
  //p7_trace_DumpAnnotated(stdout, vtr, gm, dsq);

  if (!diagfp) 
    {
      if (esl_FCompareAbs(sc, vsc,   fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(sc, fsc,   fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(sc, bsc,   fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(sc, asc_f, fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(sc, asc_b, fbtol) != eslOK) esl_fatal(failmsg);
      if (p7_trace_Compare(tr, vtr, 0)      != eslOK) esl_fatal(failmsg);

      if (ascmatrix_validate(apu, apd, anch, D, dectol, errbuf)  != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
      if (ascmatrix_trace_compare(tr, apu, apd, anch, D, dectol) != eslOK) esl_fatal(failmsg);
      if (ascmatrix_compare_std (rxd, apu, apd, anch, D, dectol) != eslOK) esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g %20g\n",
	    ascmatrix_max_over_one(apu, apd, anch, D),
	    ascmatrix_trace_maxdiff(tr, apu, apd, anch, D),
	    ascmatrix_maxdiff_std(rxd, apu, apd, anch, D));

  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(abu);  p7_refmx_Destroy(abd);
  p7_refmx_Destroy(apu);  p7_refmx_Destroy(apd);
  p7_refmx_Destroy(rxd);  p7_refmx_Destroy(rxb);
  p7_refmx_Destroy(rxf);  p7_refmx_Destroy(rxv);
  p7_trace_Destroy(vtr);  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}

/* "singlemulti" test
 * 
 * The "singlemulti" test (single path, multiple domains) samples
 * a contrived profile/seq/anchorset triplet that has only one 
 * path when anchor set constrained; that is, P(\pi | x, A, H) = 1.
 * This must be a glocal path, and may contain multiple domains.
 * 
 * Then:
 *   1. Trace score = ASC Fwd score = ASC Bck score.
 *   2. ASC Decoding matrix has 1.0 for all cells in trace
 *   
 *****************************************************************
 * May fail stochastically. Analysis: SRE:J13/126 26-Jun-14
 *  ./reference_asc_decoding_utest --diag singlemulti -N 10000
 * Field 1: max_over_one  (cells in Validation)
 * Field 2: trace_maxdiff (cells in trace_compare)
 * 
 * For fb tests redone here, use tolerance from reference_asc_fwdbacK: 0.002.
 * 
 * Single pathed test, so default vs. exact logsum has no effect; single tolerance.
 * 
 * Default:
 *    max_over_one: mean 1.5e-6 sd 9e-6 range up to 0.0006   => 0.01
 *    trace_maxdiff:  mean 1e-7 sd 1e-5 range -0.001..0.0003 => 0.01
 * Exact is similar.
 * This seems excessively fat tailed, but I think it's because the test
 * can generate quite long test sequences on occasion.
 */
static void
utest_singlemulti(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding singlemulti unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_ANCHOR  *anch      = NULL;
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  P7_REFMX   *apu       = p7_refmx_Create(M, 20);
  P7_REFMX   *apd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, asc_f, asc_b;
  float       fbtol     = 0.002;
  float       dectol    = 0.01;
  char        errbuf[eslERRBUFSIZE];
  
  if ( p7_modelsample_SinglePathedASC(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc) != eslOK) esl_fatal(failmsg);

  if ( p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd, &asc_f)             != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd, &asc_b)             != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, abu, abd, apu, apd) != eslOK) esl_fatal(failmsg);

  //p7_trace_DumpAnnotated(stdout, tr, gm, dsq);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);
  //printf("### ASC Decode UP:\n");    p7_refmx_Dump(stdout, apu);
  //printf("### ASC Decode DOWN:\n");  p7_refmx_Dump(stdout, apd);
  //p7_trace_DumpAnnotated(stdout, tr, gm, dsq);

  if (!diagfp)
    {
      if (esl_FCompareAbs(sc, asc_f, fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(sc, asc_b, fbtol) != eslOK) esl_fatal(failmsg);

      if (ascmatrix_validate(apu, apd, anch, D, dectol, errbuf)  != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
      if (ascmatrix_trace_compare(tr, apu, apd, anch, D, dectol) != eslOK) esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g\n",
	    ascmatrix_max_over_one(apu, apd, anch, D),
	    ascmatrix_trace_maxdiff(tr, apu, apd, anch, D));	    

  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(abu);  p7_refmx_Destroy(abd);
  p7_refmx_Destroy(apu);  p7_refmx_Destroy(apd);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}

/* "multisingle" test
 *
 * The "multisingle" test (multiple path, single domain) samples
 * a contrived profile/seq/anchor triplet for which all paths must
 * pass through the anchor. Thus, standard and ASC F/B give the same
 * results. 
 *
 * To make this work, the contrived model is uniglocal, with a chosen
 * match state k0 that emits an anchor residue X with probability 1,
 * and the contrived sequence has a single X at the anchor position
 * (thus forcing all paths to pass through it, even without an ASC
 * algorithm).
 *
 * Then: 
 *    1. Fwd score = ASC Fwd score = Bck score = ASC Bck score.
 *    2. ASC Decoding matrix = standard Decoding mx, cell by cell.
 *    
 *****************************************************************
 * May fail stochastically. Analysis SRE:J13/126 26-Jun-14
 *  ./reference_asc_decoding_utest --diag multisingle -N 10000
 * Field 1 = max_over_one (ASC Matrix cells in Validate())
 * Field 2 = maxdiff_std  (reference decoding cells compared to ASC decoding)
 * 
 * Also repeats fwdback tests; analysis in reference_asc_fwdback 
 * gave fbtol = 0.001 | 0.002.
 * 
 * Default:
 *  max_over_one: mean 1e-5 sd 3e-5 max 0.0003           => 0.01
 *  maxdiff_std:  mean 8e-9 sd 1e-7 range -1e-7 .. 1e-7  => 1e-5
 *
 * Exact:
 *  max_over_one: mean 4e-8 sd 9e-8 max 1e-6             => 1e-4
 *  maxdiff_std:  mean 1e-8 sd 1e-7 range -2e-7 .. 3e-7  => 1e-5
 */
static void
utest_multisingle(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding multisingle unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_ANCHOR  *anch      = NULL;
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxd       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *apu       = p7_refmx_Create(M, 20);
  P7_REFMX   *apd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, fsc, bsc, asc_f, asc_b;
  float       fbtol     = ( p7_logsum_IsSlowExact() ? 0.001 : 0.002 );
  float       onetol    = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  float       dectol    = 1e-5;
  char        errbuf[eslERRBUFSIZE];
  
  if ( p7_modelsample_AnchoredUni(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc) != eslOK) esl_fatal(failmsg);

  if ( p7_ReferenceForward (dsq, L, gm, rxf,       &fsc) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceBackward(dsq, L, gm, rxd,       &bsc) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceDecoding(dsq, L, gm, rxf, rxd,  rxd)  != eslOK) esl_fatal(failmsg);

  if ( p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd,  &asc_f)            != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCBackward(dsq, L, gm, anch, D, apu, apd,  &asc_b)            != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Fwd:\n");    p7_refmx_Dump(stdout, rxf);
  //printf("### Reference Decode:\n");    p7_refmx_Dump(stdout, rxd);
  //printf("### ASC Fwd UP:\n");       p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");     p7_refmx_Dump(stdout, afd);
  //printf("### ASC Decode UP:\n");    p7_refmx_Dump(stdout, apu);
  //printf("### ASC Decode DOWN:\n");  p7_refmx_Dump(stdout, apd);
  //printf("FWD = BCK = ASC_FWD = ASC_BCK = %.2f\n", fsc);

  if (!diagfp) 
    {
      if (esl_FCompareAbs(fsc, bsc,   fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(fsc, asc_f, fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(bsc, asc_b, fbtol) != eslOK) esl_fatal(failmsg);

      if (ascmatrix_validate(apu, apd, anch, D, onetol, errbuf) != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
      if (ascmatrix_compare_std(rxd, apu, apd, anch, D, dectol) != eslOK) esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g\n",
	    ascmatrix_max_over_one(apu, apd, anch, D),
	    ascmatrix_maxdiff_std(rxd, apu, apd, anch, D));

  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(apu);  p7_refmx_Destroy(apd);
  p7_refmx_Destroy(rxf);  p7_refmx_Destroy(rxd);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}

/* "multipath_local" test.
 * 
 * Contrives profile/seq/anchorset for which all paths must pass
 * through the anchorset.  In this version (alone amongst the unit
 * tests), local alignment is allowed, in a dual-mode (local/glocal)
 * model.
 * 
 * In order to allow local alignment transitions, while enforcing the
 * anchorset, only M states can generate residues. N/C/J/I emission is
 * prevented by setting L=0 length model and all tMI = 0. Thus, model
 * must be unihit L=0 dual-mode, and anchor set has a single D=1
 * domain.
 * 
 * Anchor 'set' is enforced by choosing special residue X; only anchor state
 * Mk0 can generate this residue with nonzero probability; target sequence
 * has exactly 1 X residue, at the anchor.
 * 
 * As in fwdbck, we test for:
 *    1. Fwd score = ASC Fwd score
 *    2. Bck score = ASC Bck score
 * and specifically for Decoding, we test:
 *    3. ASC DP UP/DOWN matrices pass validation.
 *    4. When marginalized UP+DOWN, ASC decoding matrices ==
 *       standard decoding, cell by cell.   
 *       
 *****************************************************************
 * May fail stochastically. SRE:J13/126 26-Jun-14
 *  ./reference_asc_decoding_utest --diag multipath_local -N 10000
 * Field 1: max_over_one
 * Field 2: maxdiff_std
 * 
 * For F/B tests, borrows tolerance analysis from reference_asc_fwdback:
 * 0.0001 | 0.002.
 *
 * maxdiff_std gives me zeros empirically, but I don't understand why
 * it should, so I set a nonzero threshold anyway.
 *   
 * Default: 
 *   max_over_one: mean 3e-5 sd 6e-5 max 0.0004 => 0.01
 *   maxdiff_std:  all zero (!)                 => 1e-6
 * 
 * Exact: 
 *   max_over_one: mean 6e-8 sd 1e-7 max 2e-6   => 0.0001
 *   maxdiff_std:  all zero                     
 */
static void
utest_multipath_local(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding multipath_local unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_ANCHOR  *anch      = NULL;
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxd       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *apu       = p7_refmx_Create(M, 20);   // will use as both Bck and Decode
  P7_REFMX   *apd       = p7_refmx_Create(M, 20);   // ditto
  int         D;
  float       sc, fsc, bsc, asc_f, asc_b;
  float       fbtol     = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.002 );
  float       onetol    = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01  );
  float       dectol    = 1e-6;
  char        errbuf[eslERRBUFSIZE];
  
  if ( p7_modelsample_AnchoredLocal(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc) != eslOK) esl_fatal(failmsg);

  if ( p7_ReferenceForward    (dsq, L, gm, rxf,      &fsc)   != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceBackward   (dsq, L, gm, rxd,      &bsc)   != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceDecoding   (dsq, L, gm, rxf, rxd, rxd)    != eslOK) esl_fatal(failmsg);

  if ( p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd,  &asc_f)            != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCBackward(dsq, L, gm, anch, D, apu, apd,  &asc_b)            != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Fwd:\n"); p7_refmx_Dump(stdout, rxf);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);

  if (!diagfp)
    {
      if (esl_FCompareAbs(fsc, bsc,   fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(fsc, asc_f, fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(bsc, asc_b, fbtol) != eslOK) esl_fatal(failmsg);

      if (ascmatrix_validate(apu, apd, anch, D, onetol, errbuf) != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
      if (ascmatrix_compare_std(rxd, apu, apd, anch, D, dectol) != eslOK) esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g\n",
	    ascmatrix_max_over_one(apu, apd, anch, D),
	    ascmatrix_maxdiff_std(rxd, apu, apd, anch, D));    

  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);
  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(apu);
  p7_refmx_Destroy(apd);
  p7_refmx_Destroy(rxd);
  p7_refmx_Destroy(rxf);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}


/* "multimulti" test
 * 
 * To get multiple domains, a different contrivance is needed. Only M
 * states can generate residues. N/C/J/I emission is prevented by
 * setting L=0 length model and all tMI=0. Model is multihit glocal.
 * Choose a special residue X; only the anchor state Mk0 can generate
 * this residue with nonzero probability. The target sequence has
 * exactly D X residues, one per domain. The profile has to be in
 * glocal-only mode, because local alignment paths may dip in and out
 * of the target sequence anywhere, not just on X residues.
 * 
 * Now:
 *    1. Fwd score = Bck score = ASC Fwd score = ASC Bck score
 *    2. ASC Decoding matrices pass validation.
 *    3. Std, ASC decoding matrices are equal cell-by-cell.
 *    
 *****************************************************************
 * May follow stochastically. SRE J13:126 26-Jun-14
 *  ./reference_asc_decoding_utest --diag multimulti -N 10000
 * Field 1: max_over_one
 * Field 2: maxdiff_std
 * 
 * For f/b comparisons, fbtol from reference_asc_fwdback: 0.0001 : 0.01
 *    
 * Default:
 *    max_over_one: mean 5e-5 sd 0.0001 max 0.001              => 0.01
 *    maxdiff_std:  mean -9e-7 sd 0.0001 range -0.001 .. 0.001 => 0.01
 * Exact:
 *    max_over_one: mean 3e-7 sd 7e-7 max 2e-4                 => 0.001
 *    maxdiff_std:  mean 1e-8 sd 1e-6 range -2e-5 .. 2e-5      => 0.001
 */
static void
utest_multimulti(FILE *diagfp, ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding multimulti unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_ANCHOR  *anch      = NULL;
  P7_REFMX   *rxf       = p7_refmx_Create(M, 20);
  P7_REFMX   *rxd       = p7_refmx_Create(M, 20);
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *apu       = p7_refmx_Create(M, 20);
  P7_REFMX   *apd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, fsc, bsc, asc_f, asc_b;
  float       fbtol     = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  float       dectol    = ( p7_logsum_IsSlowExact() ? 0.001  : 0.01 );
  char        errbuf[eslERRBUFSIZE];
  
  if ( p7_modelsample_AnchoredMulti(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc) != eslOK) esl_fatal(failmsg);

  if ( p7_ReferenceForward (dsq, L, gm, rxf,     &fsc) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceBackward(dsq, L, gm, rxd,     &bsc) != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceDecoding(dsq, L, gm, rxf, rxd, rxd) != eslOK) esl_fatal(failmsg);

  if ( p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd,  &asc_f)            != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCBackward(dsq, L, gm, anch, D, apu, apd,  &asc_b)            != eslOK) esl_fatal(failmsg);
  if ( p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, apu, apd, apu, apd) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Fwd:\n"); p7_refmx_Dump(stdout, rxf);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);

  if (!diagfp)
    {
      if (esl_FCompareAbs(fsc, bsc,   fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(fsc, asc_f, fbtol) != eslOK) esl_fatal(failmsg);
      if (esl_FCompareAbs(bsc, asc_b, fbtol) != eslOK) esl_fatal(failmsg);

      if (ascmatrix_validate(apu, apd, anch, D, dectol, errbuf) != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
      if (ascmatrix_compare_std(rxd, apu, apd, anch, D, dectol) != eslOK) esl_fatal(failmsg);
    }
  else
    fprintf(diagfp, "%20g %20g\n",
	    ascmatrix_max_over_one(apu, apd, anch, D),
	    ascmatrix_maxdiff_std(rxd, apu, apd, anch, D));    

  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(apu);  p7_refmx_Destroy(apd);
  p7_refmx_Destroy(rxf);  p7_refmx_Destroy(rxd);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
}
#endif /*p7REFERENCE_ASC_DECODING_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/



/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7REFERENCE_ASC_DECODING_TESTDRIVE

#include "p7_config.h"

#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

#if p7_DEVELOPMENT
#define p7_RNGSEED "0"
#else
#define p7_RNGSEED "42"
#endif

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",     eslARG_NONE,      FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",     eslARG_INT,  p7_RNGSEED, NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-N",     eslARG_INT,       "100", NULL, NULL,  NULL,  NULL, NULL, "number of times to run utest for --diag",        0 },
  { "--diag", eslARG_STRING,     NULL, NULL, NULL,  NULL,  NULL, NULL, "dump data on a utest's chance failure rate",     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for reference ASC posterior decoding";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng  = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  int             M    = 10;
  int             N    = esl_opt_GetInteger(go, "-N");

  if (esl_opt_IsOn(go, "--diag"))
    {  // --diag is for studying chance failure rates, error distributions
      char *which = esl_opt_GetString(go, "--diag");

      if      (strcmp(which, "overwrite")       == 0) while (N--) utest_overwrite      (stdout, rng, M, abc);
      else if (strcmp(which, "singlepath")      == 0) while (N--) utest_singlepath     (stdout, rng, M, abc);
      else if (strcmp(which, "singlesingle")    == 0) while (N--) utest_singlesingle   (stdout, rng, M, abc);
      else if (strcmp(which, "singlemulti")     == 0) while (N--) utest_singlemulti    (stdout, rng, M, abc);
      else if (strcmp(which, "multisingle")     == 0) while (N--) utest_multisingle    (stdout, rng, M, abc);
      else if (strcmp(which, "multipath_local") == 0) while (N--) utest_multipath_local(stdout, rng, M, abc);
      else if (strcmp(which, "multimulti")      == 0) while (N--) utest_multimulti     (stdout, rng, M, abc);
      else esl_fatal("--diag takes: overwrite, singlepath, singlesingle, singlemulti, multisingle, multipath_local, multimulti");
    }
  else // Running the unit tests is what we usually do:
    {
      fprintf(stderr, "## %s\n", argv[0]);
      fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

      utest_overwrite      (NULL, rng, M, abc);  
      utest_singlepath     (NULL, rng, M, abc);
      utest_singlesingle   (NULL, rng, M, abc);
      utest_singlemulti    (NULL, rng, M, abc);
      utest_multisingle    (NULL, rng, M, abc);
      utest_multipath_local(NULL, rng, M, abc);
      utest_multimulti     (NULL, rng, M, abc);

      fprintf(stderr, "#  status = ok\n");
    }

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7REFERENCE_ASC_DECODING_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/





/*****************************************************************
 * @LICENSE@
 *
 * SVN $URL$
 * SVN $Id$
 *****************************************************************/
