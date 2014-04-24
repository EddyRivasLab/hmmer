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
#if eslDEBUGLEVEL >= 1
  p7_refmx_Zero(apu, M, L);
  p7_refmx_Zero(apd, M, L);
#endif
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


  for (d = 0; d < D; d++)
    {
      /* UP matrix */
      iend = (d == 0 ? 1 : anch[d-1].n1+1);
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
	  denom = *ppp;	 /* pick up what we stashed earlier; we're now going to finish row i and renormalize if needed */
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

	  //#ifndef p7REFERENCE_ASC_DECODING_TESTDRIVE
	  denom = 1.0 / denom;		                                             // multiplication may be faster than division
	  ppp   = apu->dp[i] + p7R_NSCELLS;                                          // that's k=1 in UP
	  for (s = 0; s < (anch[d].n2-1) * p7R_NSCELLS; s++) *ppp++ *= denom;        // UP matrix row i renormalized
	  ppp   = apd->dp[i] + (anch[d].n2) * p7R_NSCELLS;                           // that's k0 in DOWN
	  for (s = 0; s < (M - anch[d].n2 + 1) * p7R_NSCELLS; s++) *ppp++ *= denom;  // DOWN matrix row i renormalized
	  for (s = 0; s < p7R_NXCELLS; s++) ppp[s] *= denom;
	  //#endif
	} /* end loop over i's in UP chunk for domain d */

      /* One last wing-retracted row:
       * On last row of each UP chunk (anch[d].i-1 row), those DG's
       * are reachable on a G(i)->D1...Dk-1->MG(i+1,k) path that we
       * need to unfold; and that MG(i+1,k) is the anchor itself,
       * whose value is in the DOWN matrix.
       */
      i     = anch[d].n1;	
      k     = anch[d].n2;
      xG    = *(afd->dp[i-1] + (M+1) * p7R_NSCELLS + p7R_G); // ugly, sorry. Reach out and get the xG value on row i.
      bckp  = abd->dp[i]      + k * p7R_NSCELLS;             // *bckp is the anchor cell i0,k0 (in DOWN)
      rsc   = gm->rsc[dsq[i]] + k * p7P_NR;
      delta = expf(xG + TSC(p7P_GM, k-1) + *rsc + bckp[p7R_MG] - totsc);  // k-1 in TSC because GMk is stored off by one, in -1 slot
      ppp   = apu->dp[i-1] + p7R_NSCELLS;                    // ppp starts on i-1,k=1 cells
      for (k = 1; k < anch[d].n2; k++) 
	{
	  ppp[p7R_DG] += delta;   // unlike the main wing retraction above, there's only one path going through these 
	  ppp += p7R_NSCELLS;     //   DGk's (the G->Mk0 entry into the anchor cell), so there's no retro-accumulation of delta
	}
      

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

	      fwdp += p7R_NSCELLS;
	      bckp += p7R_NSCELLS;
	      ppp  += p7R_NSCELLS;
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
	  //#ifndef p7REFERENCE_ASC_DECODING_TESTDRIVE
	  if (d == D-1) 
	    { // UP matrices only go through anch[D-1].i-1, so for last DOWN chunk, we need to renormalize.
	      denom = 1.0 / denom;
	      ppp = apd->dp[i] + (anch[d].n2) * p7R_NSCELLS;                           // that's k0 in DOWN
	      for (s = 0; s < (M - anch[d].n2 + 1) * p7R_NSCELLS; s++) *ppp++ *= denom; 
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

/* General note:
 * The unit tests are essentially the same as the unit tests
 * for Forward/Backward, but carried through to the additional
 * step of posterior decoding. Single pathed tests verify 
 * that the decoding matrix has pp=1.0 along the trace, 0.0
 * elsewhere. Ensemble tests compare standard to ASC matrices.
 */

/* ascmatrix_validate()
 *
 * Validate an ASC decoding matrix UP/DOWN pair.
 * Testing here is less thorough than p7_refmx_Validate().
 * Here we're just looking for obviously bad values.
 */
static int
ascmatrix_validate(P7_REFMX *apu, P7_REFMX *apd, P7_COORD2 *anch, int D, float epsilon, char *errbuf)
{
  int L = apu->L;
  int M = apu->M;
  int d,iz,i,k,s;
  float val;
  
  /* i=0..anch[0].i-1 specials are a boundary case... */
  iz = (D == 0 ? 1 : anch[0].n1);
  for (i = 0; i < iz; i++)
    for (s = 0; s < p7R_NXCELLS; s++)
      {
	val = P7R_XMX(apd,i,s);
	if (isnan(val) || val < 0. || val > 1.0+epsilon) 
	  ESL_FAIL(eslFAIL, errbuf, "bad value at i=%d, %s", i, p7_refmx_DecodeSpecial(s));
      }	

  for (d = 0; d < D; d++)
    {  
      /* UP sector for domain d */
      iz = (d == 0 ? 1  : anch[d-1].n1+1);
      for (i = iz; i < anch[d].n1; i++)
	for (k = 1; k < anch[d].n2; k++)
	  for (s = 0; s < p7R_NSCELLS; s++)
	    {
	      val = P7R_MX(apd,i,k,s);
	      if (isnan(val) || val < 0. || val > 1.0+epsilon) 
		ESL_FAIL(eslFAIL, errbuf, "bad value at i=%d, k=%d, %s", i, k, p7_refmx_DecodeState(s));
	    }	
      
      /* DOWN sector for domain d */
      iz = ( d < D-1 ? anch[d+1].n1 : L+1);
      for (i = anch[d].n1; i < iz; i++)
	{
	  for (k = anch[d].n2; k <= M; k++)
	    for (s = 0; s < p7R_NSCELLS; s++)
	    {
	      val = P7R_MX(apd,i,k,s);
	      if (isnan(val) || val < 0. || val > 1.0+epsilon) 
		ESL_FAIL(eslFAIL, errbuf, "bad value at i=%d, k=%d, %s", i, k, p7_refmx_DecodeState(s));
	    }	

	  /* Specials are in the DOWN matrix */
	  for (s = 0; s < p7R_NXCELLS; s++)
	    {
	      val = P7R_XMX(apd,i,s);
	      if (isnan(val) || val < 0. || val > 1.0+epsilon) 
		ESL_FAIL(eslFAIL, errbuf, "bad value at i=%d, %s", i, p7_refmx_DecodeSpecial(s));
	    }	
	}
    }
  return eslOK;
}


/* ascmatrix_compare()
 * 
 * This version of ascmatrix_compare() is very similar to
 * the one in reference_asc_fwdback.c, but is specialized
 * to decoding matrices. Because Fwd and Bck calculations
 * score prefixes and suffixes, respectively, there are
 * some expected differences when comparing ASC to Std
 * F/B matrices (i.e. there are path prefixes/suffixes that
 * can have finite probability in standard DP, but not
 * in ASC, but even in standard DP those paths cannot 
 * be completed with nonzero probability).
 * 
 * Returns <eslOK> if standard DP matrices match the
 * ASC matrices, <eslFAIL> if they don't.
 */
static int
ascmatrix_compare(P7_REFMX *std, P7_REFMX *apu, P7_REFMX *apd, P7_COORD2 *anch, int D, float epsilon)
{
  int M         = std->M;
  int L         = std->L;
  int killmenow = FALSE;
  int d,i,k,s,iz;
#ifdef p7_DEBUGGING
  killmenow = TRUE;
#endif

  ESL_DASSERT1( (apu->M == M && apu->L == L));
  ESL_DASSERT1( (apd->M == M && apd->L == L));
  ESL_DASSERT1( (std->type == p7R_DECODING  && apu->type == p7R_ASC_DECODE_UP && apd->type == p7R_ASC_DECODE_DOWN) );

  /* i=0..anch[0].i-1 specials are a boundary case... */
  iz = (D == 0 ? 1 : anch[0].n1);
  for (i = 0; i < iz; i++)
    for (s = 0; s < p7R_NXCELLS; s++)
      if (esl_FCompareAbs( P7R_XMX(std,i,s), P7R_XMX(apd,i,s), epsilon) == eslFAIL) 
	{ if (killmenow) abort(); return eslFAIL; }

  for (d = 0; d < D; d++)
    {  
      /* UP sector for domain d */
      iz = (d == 0 ? 1  : anch[d-1].n1+1);
      for (i = iz; i < anch[d].n1; i++)
	for (k = 1; k < anch[d].n2; k++)
	  for (s = 0; s < p7R_NSCELLS; s++)
	    if (esl_FCompareAbs( P7R_MX(std,i,k,s), P7R_MX(apu,i,k,s), epsilon) == eslFAIL)
	      { if (killmenow) abort(); return eslFAIL; }
      
      /* DOWN sector for domain d */
      iz = ( d < D-1 ? anch[d+1].n1 : L+1);
      for (i = anch[d].n1; i < iz; i++)
	{
	  for (k = anch[d].n2; k <= M; k++)
	    for (s = 0; s < p7R_NSCELLS; s++)
	      if (esl_FCompareAbs( P7R_MX(std,i,k,s), P7R_MX(apd,i,k,s), epsilon) == eslFAIL) 
		{ if (killmenow) abort(); return eslFAIL; }

	  /* Specials are in the DOWN matrix */
	  for (s = 0; s < p7R_NXCELLS; s++)
	    if (esl_FCompareAbs( P7R_XMX(std,i,s), P7R_XMX(apd,i,s), epsilon) == eslFAIL) 
	      { if (killmenow) abort(); return eslFAIL; }
	}
    }
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
ascmatrix_trace_compare(P7_TRACE *tr, P7_REFMX *apu, P7_REFMX *apd, P7_COORD2 *anch, int D, float epsilon)
{
  static int sttype[p7T_NSTATETYPES] = { -1, p7R_ML, p7R_MG, p7R_IL, p7R_IG, p7R_DL, p7R_DG, -1, p7R_N, p7R_B, p7R_L, p7R_G, p7R_E, p7R_C, p7R_J, -1 }; /* sttype[] translates trace idx to DP matrix idx*/
  P7_REFMX *ppt  = NULL;
  int       M    = apu->M;
  int       L    = apu->L;
  int       z, i;
  int       status;

  /* Contract checks */
  ESL_DASSERT1(( apd->M == M && tr->M == M));
  ESL_DASSERT1(( apd->L == L && tr->L == L));
  ESL_DASSERT1(( apu->type == p7R_ASC_DECODE_UP && apd->type == p7R_ASC_DECODE_DOWN ));

  /* Make a standard decoding matrix with 1.0 whereever
   * the trace visits a cell, and 0.0 elsewhere
   */
  ppt = p7_refmx_Create(M, L);
  p7_refmx_Zero(ppt, M, L);     // This also sets M,L, and type=decoding in <ppt>.
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
      
  /* Now we can just compare this std DP matrix to the ASC matrices */
  status = ascmatrix_compare(ppt, apu, apd, anch, D, epsilon);
  
  p7_refmx_Destroy(ppt);
  return status;
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
 */
static void
utest_singlepath(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding singlepath unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = p7_profile_Create(M, bg->abc);
  ESL_SQ     *sq        = esl_sq_CreateDigital(bg->abc);
  P7_TRACE   *tr        = p7_trace_Create();
  P7_TRACE   *vtr       = p7_trace_Create();
  P7_COORD2  *anch      = malloc(sizeof(P7_COORD2));   // utest is uniglocal, so it will have D=1 and only one anchor.
  P7_REFMX   *rxv       = p7_refmx_Create(M, 20);      // 20 is arbitrary; all DP routines will reallocate as needed.
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
  float       epsilon   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  char        errbuf[eslERRBUFSIZE];
  int         status;

  /* Sample single-pathed HMM, create uniglocal L=0 profile */
  if ((status = p7_modelsample_SinglePathed(rng, M, bg->abc, &hmm)) != eslOK) esl_fatal(failmsg);
  if ( p7_profile_ConfigUniglocal(gm, hmm, bg, 0)                   != eslOK) esl_fatal(failmsg); 

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
	  anch[0].n1 = tr->i[z];
	  anch[0].n2 = tr->k[z];
	  break;
	}
      }

  /* Run the DP routines, both standard and ASC */
  if ((status = p7_ReferenceViterbi (sq->dsq, sq->n, gm, rxv, vtr, &vsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceForward (sq->dsq, sq->n, gm, rxf,      &fsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceBackward(sq->dsq, sq->n, gm, rxb,      &bsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceDecoding(sq->dsq, sq->n, gm, rxf, rxb,  rxd)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceASCForward (sq->dsq, sq->n, gm, anch, D, afu, afd, &asc_f))             != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(sq->dsq, sq->n, gm, anch, D, abu, abd, &asc_b))             != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCDecoding(sq->dsq, sq->n, gm, anch, D, afu, afd, abu, abd, apu, apd)) != eslOK) esl_fatal(failmsg);

  /* Dumping, for debugging. 
   * Uncomment what you need, for what problems you have.
   */
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
  if (esl_FCompare(sc, vsc,   epsilon) != eslOK) esl_fatal(failmsg);  // generated trace score = Viterbi score
  if (esl_FCompare(sc, fsc,   epsilon) != eslOK) esl_fatal(failmsg);  //  ... = Forward score
  if (esl_FCompare(sc, bsc,   epsilon) != eslOK) esl_fatal(failmsg);  //  ... = Backward score
  if (esl_FCompare(sc, asc_f, epsilon) != eslOK) esl_fatal(failmsg);  //  ... = ASC Forward score
  if (esl_FCompare(sc, asc_b, epsilon) != eslOK) esl_fatal(failmsg);  //  ... = ASC Backward score
  if (p7_trace_Compare(tr, vtr, 0)     != eslOK) esl_fatal(failmsg);  // generated trace = Viterbi trace

  if (ascmatrix_validate(apu, apd, anch, D, epsilon, errbuf)  != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
  if (ascmatrix_trace_compare(tr, apu, apd, anch, D, epsilon) != eslOK) esl_fatal(failmsg);
  if (ascmatrix_compare     (rxd, apu, apd, anch, D, epsilon) != eslOK) esl_fatal(failmsg);

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
 * P(\pi | x, H) = 1). This must be a glocal single domain
 * path. Choose any anchor A in that domain. 
 * 
 * Then:
 *   1. Viterbi = Fwd = Bck = ASC Fwd = ASC Bck scores.
 *   2. Viterbi trace = trace that generated the test sequence. 
 *   3. ASC Decoding matrix has 1.0 for all cells in trace
 *   4. ASC Decoding matrix = standard Decoding mx, cell by cell.
 */
static void
utest_singlesingle(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding singlesingle unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_TRACE   *vtr       = p7_trace_Create();
  P7_COORD2  *anch      = NULL;
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
  float       epsilon   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  char        errbuf[eslERRBUFSIZE];
  int         status;

  if ((status = p7_modelsample_SinglePathedSeq(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceViterbi (dsq, L, gm, rxv, vtr, &vsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceForward (dsq, L, gm, rxf,      &fsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceBackward(dsq, L, gm, rxb,      &bsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceDecoding(dsq, L, gm, rxf, rxb,  rxd)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd, &asc_f))             != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd, &asc_b))             != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, abu, abd, apu, apd)) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Vit:\n"); p7_refmx_Dump(stdout, rxv);
  //printf("### Reference Fwd:\n"); p7_refmx_Dump(stdout, rxf);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);
  //p7_trace_DumpAnnotated(stdout, tr,  gm, dsq);
  //p7_trace_DumpAnnotated(stdout, vtr, gm, dsq);

  if (esl_FCompare(sc, vsc, epsilon)   != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(sc, fsc, epsilon)   != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(sc, bsc, epsilon)   != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(sc, asc_f, epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(sc, asc_b, epsilon) != eslOK) esl_fatal(failmsg);
  if (p7_trace_Compare(tr, vtr, 0)     != eslOK) esl_fatal(failmsg);

  if (ascmatrix_validate(apu, apd, anch, D, epsilon, errbuf)  != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
  if (ascmatrix_trace_compare(tr, apu, apd, anch, D, epsilon) != eslOK) esl_fatal(failmsg);
  if (ascmatrix_compare     (rxd, apu, apd, anch, D, epsilon) != eslOK) esl_fatal(failmsg);

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
 */
static void
utest_singlemulti(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding singlemulti unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_COORD2  *anch      = NULL;
  P7_REFMX   *afu       = p7_refmx_Create(M, 20);
  P7_REFMX   *afd       = p7_refmx_Create(M, 20);
  P7_REFMX   *abu       = p7_refmx_Create(M, 20);
  P7_REFMX   *abd       = p7_refmx_Create(M, 20);
  P7_REFMX   *apu       = p7_refmx_Create(M, 20);
  P7_REFMX   *apd       = p7_refmx_Create(M, 20);
  int         D;
  float       sc, asc_f, asc_b;
  float       epsilon   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  char        errbuf[eslERRBUFSIZE];
  int         status;
  
  if ((status = p7_modelsample_SinglePathedASC(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd, &asc_f))             != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd, &asc_b))             != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, abu, abd, apu, apd)) != eslOK) esl_fatal(failmsg);

  //p7_trace_DumpAnnotated(stdout, tr, gm, dsq);
  //printf("### ASC Fwd UP:\n");    p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");  p7_refmx_Dump(stdout, afd);

  if (esl_FCompare(sc, asc_f, epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(sc, asc_b, epsilon) != eslOK) esl_fatal(failmsg);

  if (ascmatrix_validate(apu, apd, anch, D, epsilon, errbuf)  != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
  if (ascmatrix_trace_compare(tr, apu, apd, anch, D, epsilon) != eslOK) esl_fatal(failmsg);

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

/* utest_multisingle()
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
 *    2. Viterbi score >= sampled trace score.
 *    3. Fwd/Bck scores >= Viterbi score.
 *    4. ASC Decoding matrix = standard Decoding mx, cell by cell.
 */
static void
utest_multisingle(ESL_RANDOMNESS *rng, int M, const ESL_ALPHABET *abc)
{
  char        failmsg[] = "reference_asc_decoding multisingle unit test failed";
  P7_BG      *bg        = p7_bg_Create(abc);
  P7_HMM     *hmm       = NULL;
  P7_PROFILE *gm        = NULL;
  ESL_DSQ    *dsq       = NULL; 
  int         L;
  P7_TRACE   *tr        = NULL;
  P7_COORD2  *anch      = NULL;
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
  float       epsilon   = ( p7_logsum_IsSlowExact() ? 0.0001 : 0.01 );
  char        errbuf[eslERRBUFSIZE];
  int         status;
  
  if ((status = p7_modelsample_AnchoredUni(rng, M, bg, &hmm, &gm, &dsq, &L, &tr, &anch, &D, &sc)) != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceViterbi (dsq, L, gm, rxv, NULL, &vsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceForward (dsq, L, gm, rxf,       &fsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceBackward(dsq, L, gm, rxb,       &bsc)) != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceDecoding(dsq, L, gm, rxf, rxb,  rxd))  != eslOK) esl_fatal(failmsg);

  if ((status = p7_ReferenceASCForward (dsq, L, gm, anch, D, afu, afd,  &asc_f))            != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCBackward(dsq, L, gm, anch, D, abu, abd,  &asc_b))            != eslOK) esl_fatal(failmsg);
  if ((status = p7_ReferenceASCDecoding(dsq, L, gm, anch, D, afu, afd, abu, abd, apu, apd)) != eslOK) esl_fatal(failmsg);

  //printf("### Reference Fwd:\n");    p7_refmx_Dump(stdout, rxf);
  printf("### Reference Decode:\n");    p7_refmx_Dump(stdout, rxd);
  //printf("### ASC Fwd UP:\n");       p7_refmx_Dump(stdout, afu);
  //printf("### ASC Fwd DOWN:\n");     p7_refmx_Dump(stdout, afd);
  printf("### ASC Decode UP:\n");    p7_refmx_Dump(stdout, apu);
  printf("### ASC Decode DOWN:\n");  p7_refmx_Dump(stdout, apd);
  //printf("FWD = BCK = ASC_FWD = ASC_BCK = %.2f\n", fsc);

  if (esl_FCompare(fsc, bsc,   epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(fsc, asc_f, epsilon) != eslOK) esl_fatal(failmsg);
  if (esl_FCompare(bsc, asc_b, epsilon) != eslOK) esl_fatal(failmsg);
  if (sc  > vsc+epsilon)                          esl_fatal(failmsg);
  if (vsc > fsc+epsilon)                          esl_fatal(failmsg);

  if (ascmatrix_validate(apu, apd, anch, D, epsilon, errbuf) != eslOK) esl_fatal("%s:\n%s", failmsg, errbuf);
  if (ascmatrix_compare(rxd, apu, apd, anch, D, epsilon)     != eslOK) esl_fatal(failmsg);

  free(anch);
  free(dsq);
  p7_refmx_Destroy(afu);  p7_refmx_Destroy(afd);
  p7_refmx_Destroy(abu);  p7_refmx_Destroy(abd);
  p7_refmx_Destroy(apu);  p7_refmx_Destroy(apd);
  p7_refmx_Destroy(rxf);  p7_refmx_Destroy(rxb);
  p7_refmx_Destroy(rxv);  p7_refmx_Destroy(rxd);
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

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
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

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_singlepath     (rng, M, abc);
  utest_singlesingle   (rng, M, abc);
  utest_singlemulti    (rng, M, abc);
  utest_multisingle    (rng, M, abc);
  //utest_multipath_local(rng, M, abc);
  //utest_multimulti     (rng, M, abc);

  fprintf(stderr, "#  status = ok\n");

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
