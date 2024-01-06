/* Sparse dynamic programming. Posterior decoding implementation.
 * 
 * Contents:
 *   1. Sparse posterior decoding
 *   2. Unit tests
 *   3. Test driver
 */
#include <p7_config.h>

#include <math.h>

#include "easel.h"

#include "base/p7_profile.h"

#include "misc/logsum.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/sparse_decoding.h"


/*****************************************************************
 * 1. Sparse Decoding
 *****************************************************************/

/* Function:  p7_SparseDecoding()
 * Synopsis:  Posterior decoding algorithm, in sparse DP.
 *
 * Purpose:   Given Forward matrix <sxf> and Backward matrix <sxb>, that
 *            have been computed for a comparison of profile <gm> against
 *            digital sequence <dsq> of length <L>;
 *            fill in the posterior decoding matrix <sxd>.
 *            
 *            <sxb> and <sxd> can point to the same structure, in
 *            which case posterior decoding overwrites the Backward
 *            matrix; a trick that the caller might use to save
 *            allocating a new matrix. (Can't do the same with
 *            Forward, because of a detail involving CC/JJ transition
 *            decoding.)
 * 
 *            If <sxd> is an independent matrix (i.e. not overwriting
 *            <sxb>), it will be reinitialized (and possibly reallocated)
 *            to use the same sparse mask as <sxf> and <sxb>. 
 *            
 *            {M/I}(i,k) is the prob that state {MI}k generated
 *            residue <x_i>.  D(i,k) is the prob that a state path
 *            used state Dk after already generating residue <x_i>
 *            with another M state.  {BLG}(i) is the probability of
 *            being in a {BLG} state just as we start a new domain at
 *            <x_i+1>. E(i) is the probability of being in the end
 *            state on row <i>, having ended the domain on this row.
 *            These all arise from standard posterior decoding eqn's.
 *            
 *            Watch out for {NJC}, though, the states that emit on
 *            transition; these are decoded in a slightly nonstandard
 *            way. {NN/CC/JJ}(i) is the probability that we emitted
 *            residue <i> on an {NN/CC/JJ} transition. {NCJ}(i) is the
 *            probability that the state path uses {NCJ} on (i), which
 *            is a sum of emitting <i> and an initial mute transition
 *            {S->N,E->C,E->J}.
 *
 * Args:      dsq  - digital sequence, 1..L
 *            L    - length of <dsq>
 *            gm   - profile
 *            sxf  - Forward matrix, computed by caller
 *            sxb  - Backward matrix, computed by caller
 *            sxd  - space for result - the posterior decoding matrix
 *            
 * Returns:   <eslOK> on success, and <sxd> is filled by the posterior decoding
 *            calculation.
 *
 * Throws:    <eslEMEM> if a reallocation of <sxd> fails.
 */
int
p7_SparseDecoding(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_SPARSEMX *sxf, P7_SPARSEMX *sxb, P7_SPARSEMX *sxd)
{
  const P7_SPARSEMASK *sm = sxf->sm;
  const float   *dpf      = sxf->dp;
  float         *dpb      = sxb->dp;
  float         *dpd;
  const float   *xf       = sxf->xmx;
  float         *xb       = sxb->xmx;
  float         *xd;
  float         *dpd2;
  const float   *rsc;		      /* residue scores on current row i; enables MSC(k) macro shorthand */
  float          totsc  = xf[p7S_N] + xb[p7S_N];  // remember, the first ia-1 is not necessarily row 0, so xb[N] alone does not suffice.
  float          norm;                // except for numerical roundoff error accumulation, so we also explicitly renormalize each row.
  float          xC,xJ,xG;            // Forward values remembered from prev row i-1.
  int            i,k,y,z;
  float          delta;		     /* additions to DGk's, resulting from unfolding the wing-retracted entry/exit paths */
  const float   *tsc    = gm->tsc;   /* activates the TSC() macro */
#ifndef p7SPARSE_DECODING_TESTDRIVE  /* see note on renormalization further below. */
  int            x;
#endif
  int            status;

  /* Contract validation (argument checks) */
  ESL_DASSERT1 ( (sxf->type == p7S_FORWARD) );
  ESL_DASSERT1 ( (sxb->type == p7S_BACKWARD) );
  ESL_DASSERT1 ( (sxf->sm == sxb->sm) );
  ESL_DASSERT1 ( (sm->L   == L) );
  ESL_DASSERT1 ( (sm->M   == gm->M) );

  /* Assure that <sxd> is allocated large enough (we might be reusing it),
   * either because we're overwriting <sxb> (which we know is large enough),
   * or because we reinitialize it (if it isn't already using <sm>). 
   */
  if (sxd != sxb && (status = p7_sparsemx_Reinit(sxd, sm)) != eslOK) return status;
  sxd->type = p7S_DECODING;

  xJ = xC = -eslINFINITY;
  dpd = sxd->dp;
  xd  = sxd->xmx;
  for (i = 1; i <= L; i++)
    {
      rsc = gm->rsc[dsq[i]];	/* MSC(k), ISC(k) macros now get residue scores for this row */

      /* i=ia-1, initialization of a sparse segment; special storage */
      if (sm->n[i] && !sm->n[i-1])
	{
	  norm       = 0.0f;
	  xd[p7S_E]  = 0.0f;                                
	  xd[p7S_N]  = expf(xf[p7S_N] + xb[p7S_N] - totsc);                  norm += xd[p7S_N];
	  xd[p7S_J]  = expf(xf[p7S_J] + xb[p7S_J] - totsc);  
	  xd[p7S_B]  = expf(xf[p7S_B] + xb[p7S_B] - totsc);
	  xd[p7S_L]  = expf(xf[p7S_L] + xb[p7S_L] - totsc);
	  xd[p7S_G]  = expf(xf[p7S_G] + xb[p7S_G] - totsc);  xG = xf[p7S_G];
	  xd[p7S_C]  = expf(xf[p7S_C] + xb[p7S_C] - totsc);
	  xd[p7S_JJ] = xd[p7S_J];                            xJ = xf[p7S_J]; norm += xd[p7S_JJ];
	  xd[p7S_CC] = xd[p7S_C];                            xC = xf[p7S_C]; norm += xd[p7S_CC];
	  
#ifndef p7SPARSE_DECODING_TESTDRIVE /* see note on renormalization further below. */
	  norm = 1.0f/norm;
	  for (x = 0; x < p7S_NXCELLS; x++) xd[x] *= norm;
#endif

	  xf += p7S_NXCELLS;
	  xb += p7S_NXCELLS;
	  xd += p7S_NXCELLS;
	}

      /* For each sparse cell z (k=k[i][z]) on row: 
       * DG is special, because we have to unfold the wing-retracted entries and exits. 
       * Exits: unfolding right wing retractions 
       *
       * Caution: Although we're unfolding wings correctly, by design,
       * we're only storing that probability mass in i,k cells in the
       * sparse mask. The sparse mask was created by local-only
       * decoding in ForwardFilter(). Suppose there's a high
       * probability G->DDDD->Mk entry, for example. The D's won't
       * necessarily be in the mask! Thus, sparse decoding does NOT
       * give you a complete decoding of D's. I don't think we care;
       * but I'm leaving a bread crumb here, just in case.
       */
      dpd2 = dpd;
      for (delta = 0.0f, z = 0; z < sm->n[i]; z++)
	{
	  k = sm->k[i][z];
	  dpd[p7S_DG] = delta + expf(dpf[p7S_DG] + dpb[p7S_DG] - totsc);   /* because we might be overwriting backwards mx with decoding, this is the last time we can access dpb[p7S_DG] values on the row! */
	  if (z == sm->n[i]-1 || sm->k[i][z+1] != k+1) /* No cell to our right? then {MD}Gk->Dk+1..Dm->E path delta added to all subsequent stored Dk+1..m */
	    delta += expf(p7_FLogsum( dpf[p7S_MG] + TSC(p7P_MD, k),  dpf[p7S_DG] + TSC(p7P_DD,k)) 
			  + TSC(p7P_DGE,k) + xb[p7S_E] - totsc);
	  dpd += p7S_NSCELLS; dpf += p7S_NSCELLS; dpb += p7S_NSCELLS;
	}
      /* Entries: unfolding left wing retractions. The DG's for a G->D..D->Mk path are on the PREVIOUS ROW from the Mk! */
      /* All Mk on current row contributes a delta; and each delta applies to all 1..k-1 on prev row. Sparsity on both rows makes this tricky */
      /* and remember, we need the residue score for e(x_i,MGk) as well as the backwards score MGk,i */
      for (delta = 0.0f, y = sm->n[i-1]-1, z = sm->n[i]-1; z >= 0; z--)
	{
	  k = sm->k[i][z];
	  while (y >= 0 && sm->k[i-1][y] >= k) { dpd2 -= p7S_NSCELLS; dpd2[p7S_DG] += delta; y--; }

	  dpb -= p7S_NSCELLS;
	  delta       += expf(xG + TSC(p7P_GM, k-1) + MSC(k) + dpb[p7S_MG] - totsc); /* G->D1..Dk-1->Mk entry path, added to all stored D1..Dk-1 */
	}
      while (y >= 0) { dpd2 -= p7S_NSCELLS; dpd2[p7S_DG] += delta; y--; } /* dpd2 now sits on first sparse cell of prev row i-1. */
      /* note that dpb ran up and back down; dpf and dpd only ran up, and need to be run back down */
      dpf -= sm->n[i]*p7S_NSCELLS;
      dpd -= sm->n[i]*p7S_NSCELLS;

      norm = 0.0;
      for (z = 0; z < sm->n[i]; z++, dpd += p7S_NSCELLS, dpf += p7S_NSCELLS, dpb += p7S_NSCELLS)
	{
	  dpd[p7S_ML] = expf(dpf[p7S_ML] + dpb[p7S_ML] - totsc); norm += dpd[p7S_ML];
	  dpd[p7S_MG] = expf(dpf[p7S_MG] + dpb[p7S_MG] - totsc); norm += dpd[p7S_MG];
	  dpd[p7S_IL] = expf(dpf[p7S_IL] + dpb[p7S_IL] - totsc); norm += dpd[p7S_IL];
	  dpd[p7S_IG] = expf(dpf[p7S_IG] + dpb[p7S_IG] - totsc); norm += dpd[p7S_IG];
	  dpd[p7S_DL] = expf(dpf[p7S_DL] + dpb[p7S_DL] - totsc);                       // nonemitters don't count toward normalization
	  // DG was already done above.
	}
      
      /* specials on each stored row */
      if (sm->n[i]) {
	xd[p7S_JJ] = expf(  xJ      + xb[p7S_J] + gm->xsc[p7P_J][p7P_LOOP] - totsc); xJ = xf[p7S_J]; norm += xd[p7S_JJ];  // JJ,CC calculations must come before J,C; they depend on xb[J,C], which we may overwrite when we calc J,C
	xd[p7S_CC] = expf(  xC      + xb[p7S_C] + gm->xsc[p7P_C][p7P_LOOP] - totsc); xC = xf[p7S_C]; norm += xd[p7S_CC];
	xd[p7S_E]  = expf(xf[p7S_E] + xb[p7S_E] - totsc);                                  
	xd[p7S_N]  = expf(xf[p7S_N] + xb[p7S_N] - totsc);                                            norm += xd[p7S_N];
	xd[p7S_J]  = expf(xf[p7S_J] + xb[p7S_J] - totsc);  
	xd[p7S_B]  = expf(xf[p7S_B] + xb[p7S_B] - totsc);                                  
	xd[p7S_L]  = expf(xf[p7S_L] + xb[p7S_L] - totsc);                                  
	xd[p7S_G]  = expf(xf[p7S_G] + xb[p7S_G] - totsc);                            xG = xf[p7S_G];                              
	xd[p7S_C]  = expf(xf[p7S_C] + xb[p7S_C] - totsc);                   

	/* Renormalization.
	 * 
	 * Roundoff error accumulation in F/B is significant. For large
	 * target seqs, it isn't unusual to have a whole nat of
	 * difference in overall fwd vs. bck raw score; for example, fn3
	 * vs. TITIN_HUMAN. Since we only use the bck score (at i=0) for
	 * normalization above, pp's would have very large systematic
	 * error. 
	 * 
	 * To squash this error in production code, we renormalize rows.
	 * Because renormalization can hide real errors, we don't
	 * renormalize when we've compiled the code for unit testing.
	 * Default unit tests don't run large enough M/L to create a
	 * lot of error accumulation.
	 *
	 * (A similar note appears in reference_decoding.)
	 */
#ifndef p7SPARSE_DECODING_TESTDRIVE	
	norm = 1.0f/norm;
	dpd -= sm->n[i]*p7S_NSCELLS;  /* back up to start of row again */
	for (x = 0; x < sm->n[i]*p7S_NSCELLS; x++) *dpd++ *= norm;  
	for (x = 0; x < p7S_NXCELLS;          x++) xd[x]  *= norm;
#endif
	  
	xd += p7S_NXCELLS;
	xf += p7S_NXCELLS;
	xb += p7S_NXCELLS;
      }
    }
  return eslOK;
}
/*--------------- end, posterior decoding -----------------------*/



/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef p7SPARSE_DECODING_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "base/p7_bg.h"

#include "build/modelsample.h"
#include "search/modelconfig.h"
#include "misc/emit.h"

#include "dp_vector/p7_oprofile.h"
#include "dp_vector/p7_checkptmx.h"
#include "dp_vector/fwdfilter.h"

#include "dp_sparse/sparse_fwdback.h"
#include "dp_sparse/sparse_trace.h"

/* The "overwrite" utest verifies an important wrinkle in the API:
 * we're allowed to overwrite the input Backwards matrix with the new
 * posterior decoding matrix, thus saving a matrix allocation.  The
 * utest samples a random profile, compares it against a generated
 * sequence, decodes it both with and without overwriting, and
 * verifies that the resulting decoding matrices are valid and
 * identical.
 */
static void
utest_overwrite(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L)
{
  char           msg[] = "reference_decoding: overwrite unit test failed";
  ESL_SQ        *sq    = esl_sq_CreateDigital(abc);
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = p7_profile_Create(M, abc);
  P7_OPROFILE   *om    = p7_oprofile_Create(M, abc);
  P7_CHECKPTMX  *ox    = NULL;
  P7_SPARSEMASK *sm    = NULL;
  P7_SPARSEMX   *fwd   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *bck   = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *pp    = p7_sparsemx_Create(NULL);
  float         tol   = 0.0f;	/* exact match is expected! */
  char          errbuf[eslERRBUFSIZE];

  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);
  if ( p7_oprofile_Convert(gm, om)       != eslOK) esl_fatal(msg);
  
  /* Emit sequence from model, using a length model of <L>;
   * restrict the emitted sequence length to 3x (L+M), arbitrarily, to 
   * keep it down to something reasonable.
   */
  if ( p7_profile_SetLength(gm, L)      != eslOK) esl_fatal(msg);
  do {
    esl_sq_Reuse(sq);
    if (p7_ProfileEmit(rng, hmm, gm, bg, sq, /*tr=*/NULL) != eslOK) esl_fatal(msg);
  } while (sq->n > (gm->M+gm->L) * 3); 

  /* Allocate matrices, set length models */
  if ( p7_profile_SetLength(gm, sq->n)        != eslOK)  esl_fatal(msg);
  if ( p7_oprofile_ReconfigLength(om, sq->n)  != eslOK)  esl_fatal(msg);

  /* F/B filter to calculate the sparse mask */
  if ( (ox = p7_checkptmx_Create(gm->M, sq->n, ESL_MBYTES(32)))          == NULL)  esl_fatal(msg);
  if ( (sm = p7_sparsemask_Create(gm->M, sq->n))                         == NULL)  esl_fatal(msg);
  if ( p7_ForwardFilter (sq->dsq, sq->n, om, ox, /*fsc=*/NULL)           != eslOK) esl_fatal(msg);
  if ( p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH) != eslOK) esl_fatal(msg);

  /* F/B, then decode both ways */
  if (p7_SparseForward (sq->dsq, sq->n, gm, sm, fwd, /*fsc=*/NULL)  != eslOK) esl_fatal(msg);
  if (p7_SparseBackward(sq->dsq, sq->n, gm, sm, bck, /*bsc=*/NULL)  != eslOK) esl_fatal(msg);
  if (p7_SparseDecoding(sq->dsq, sq->n, gm, fwd, bck, pp)           != eslOK) esl_fatal(msg); /* <pp> is decoded independently      */
  if (p7_SparseDecoding(sq->dsq, sq->n, gm, fwd, bck, bck)          != eslOK) esl_fatal(msg); /* <bck> is overwritten with decoding */

  /* Tests. */
  if (p7_sparsemx_Compare(pp, bck, tol) != eslOK) esl_fatal(msg);
  if (p7_sparsemx_Validate(pp,  errbuf) != eslOK) esl_fatal("%s:\n%s\n", msg, errbuf);
  if (p7_sparsemx_Validate(bck, errbuf) != eslOK) esl_fatal("%s:\n%s\n", msg, errbuf);

  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
  p7_checkptmx_Destroy(ox);
  p7_sparsemask_Destroy(sm);
  p7_sparsemx_Destroy(pp);
  p7_sparsemx_Destroy(bck);
  p7_sparsemx_Destroy(fwd);
}

/* The "rowsum" utest verifies that the sum of posterior probs for all
 * emitters (M,I,NN,CC,JJ) on any row 1..L is ~1.0, for comparisons of
 * a randomly sampled profile to <N> randomly sampled homologs.
 * 
 * TODO: sparse_asc does this better, by making rowsum test part of
 * Validate().
 */
static void
utest_rowsum(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char           msg[] = "reference_decoding: rowsum unit test failed";
  ESL_SQ        *sq    = esl_sq_CreateDigital(abc);
  P7_HMM        *hmm   = NULL;
  P7_PROFILE    *gm    = p7_profile_Create(M, abc);
  P7_OPROFILE   *om    = p7_oprofile_Create(M, abc);
  P7_CHECKPTMX  *ox    = p7_checkptmx_Create(M, L, ESL_MBYTES(32));
  P7_SPARSEMASK *sm    = p7_sparsemask_Create(M, L);
  P7_SPARSEMX   *sxf   = p7_sparsemx_Create(sm);
  P7_SPARSEMX   *sxb   = p7_sparsemx_Create(sm);
  P7_SPARSEMX   *sxd   = p7_sparsemx_Create(sm);
  float          tol   = ( p7_logsum_IsSlowExact() ? 0.001 : 0.01); /* tuned to the test's default <M>,<L>,<N> */
  int            idx;
  int            i,z;
  float          rowsum;
  float          bsc, fsc;
  const float   *dpc;
  const float   *xc;
  float          max = 0.;

  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);
  if ( p7_oprofile_Convert(gm, om)       != eslOK) esl_fatal(msg);
  
  for (idx = 0; idx < N; idx++)
    {
      /* Emit sequence from model, using a length model of <L>;
       * restrict the emitted sequence length to 3x (L+M), arbitrarily, to 
       * keep it down to something reasonable.
       */
      if ( p7_profile_SetLength(gm, L)      != eslOK) esl_fatal(msg);
      do {
	esl_sq_Reuse(sq);
	if (p7_ProfileEmit(rng, hmm, gm, bg, sq, /*tr=*/NULL) != eslOK) esl_fatal(msg);
      } while (sq->n > (gm->M+gm->L) * 3); 
      if ( p7_profile_SetLength(gm,       sq->n)   != eslOK) esl_fatal(msg);
      if ( p7_oprofile_ReconfigLength(om, sq->n)   != eslOK) esl_fatal(msg);

      /* F/B filter to calculate the sparse mask */
      if ( p7_checkptmx_Reinit (ox,  gm->M, sq->n)                           != eslOK) esl_fatal(msg);
      if ( p7_ForwardFilter (sq->dsq, sq->n, om, ox, /*fsc=*/NULL)           != eslOK) esl_fatal(msg);
      if ( p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH) != eslOK) esl_fatal(msg); // sparse mask is _Reinit'ed here.

      /* F/B, then decode in place */
      if (p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf, &fsc)     != eslOK) esl_fatal(msg);
      if (p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxb, &bsc)     != eslOK) esl_fatal(msg);
      if (p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxb, sxd)     != eslOK) esl_fatal(msg); 

      /* Collect row sums from <pp> */
      dpc = sxd->dp;
      xc  = sxd->xmx;
      for (i = 1; i <= sq->n; i++)
	{
	  if (sm->n[i] && !sm->n[i-1]) /* ia-1 specials at a segment start. */
	    {
	      rowsum = xc[p7S_N] + xc[p7S_JJ] + xc[p7S_CC]; /* this even works for nonemitting row ia-1=0, because there xc[N]=1 and JJ/CC=0 */
	      max    = ESL_MAX(rowsum-1.0, max);
	      if (esl_FCompare(rowsum, 1.0, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg);
	      xc += p7S_NXCELLS;
	    }

	  for (rowsum = 0., z = 0; z < sm->n[i]; z++, dpc += p7R_NSCELLS) /* sparse main cells */
	    rowsum += dpc[p7S_ML] + dpc[p7S_MG] + dpc[p7S_IL] + dpc[p7S_IG];
	  
	  if (sm->n[i]) {
	    rowsum += xc[p7S_N] + xc[p7S_JJ] + xc[p7S_CC];
	    max    = ESL_MAX(rowsum-1.0, max);
	    if (esl_FCompare(rowsum, 1.0, /*rtol=*/0.0, tol) != eslOK) esl_fatal(msg);
	    xc += p7S_NXCELLS;
	  }
	}
      
      p7_sparsemask_Reuse(sm);
      p7_sparsemx_Reuse(sxf);
      p7_sparsemx_Reuse(sxb);
      p7_sparsemx_Reuse(sxd);
    }

  //printf("max = %.4f\n", max);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
  p7_checkptmx_Destroy(ox);
  p7_sparsemask_Destroy(sm);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxf);
}


/* The 'approx-decoding' utest compares exact posterior decoding (via
 * p7_SparseDecoding()) to stochastic approximation (by a large
 * ensemble of stochastic tracebacks). It does this for a randomly
 * sampled profile HMM of length <M> compared against one homologous
 * (generated) sequence. (Only one of them, not <N>, because the
 * stochastic tracebacks are computationally expensive.)
 * 
 * Tests:
 * 1. The two decoding approaches give identical matrices, within 
 *    a given sampling error tolerance. (Additionally, cells that
 *    are exactly zero in exact posterior decoding must not be
 *    visited in any stochastic trace.) All this is checked by 
 *    <p7_sparsemx_CompareDecoding()>.
 * 2. The two decoding matrices both Validate().
 * 
 */
static void
utest_approx_decoding(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int ntr)
{
  char           msg[]  = "sparse fwdback, approx-decoding unit test failed";
  P7_HMM        *hmm    = NULL;
  P7_PROFILE    *gm     = p7_profile_Create(M, abc);
  ESL_SQ        *sq     = esl_sq_CreateDigital(abc);       /* space for generated (homologous) target seqs              */
  P7_OPROFILE   *om     = p7_oprofile_Create(M, abc);
  P7_CHECKPTMX  *ox     = NULL;
  P7_TRACE      *tr     = NULL;
  P7_SPARSEMASK *sm     = NULL;
  P7_SPARSEMX   *sxf    = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxb    = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxd    = p7_sparsemx_Create(NULL);
  P7_SPARSEMX   *sxs    = NULL;
  float         *wrk    = NULL;	/* reusable scratch workspace needed by stochastic trace */
  int            idx;
  float          tol    = 0.02;	         /* with utest's defaults, max diff will be ~0.004 or so; exact v. approx logsum seems to make no difference   */

  /* Sample a profile. 
   * Config as usual: multihit dual-mode local/glocal, so all paths in it are valid.
   */
  if ( p7_modelsample(rng, M, abc, &hmm) != eslOK) esl_fatal(msg);
  if ( p7_profile_Config(gm, hmm, bg)    != eslOK) esl_fatal(msg);
  if ( p7_oprofile_Convert(gm, om)       != eslOK) esl_fatal(msg);

  /* Generate (sample) a sequence from the profile */
  if ( p7_profile_SetLength(gm, L)       != eslOK) esl_fatal(msg);   /* config to generate mean length of L */
  do {
    esl_sq_Reuse(sq);
    p7_ProfileEmit(rng, hmm, gm, bg, sq, NULL);
  } while (sq->n > L * 3); /* keep sequence length from getting ridiculous; long seqs do have higher abs error per cell */
  if ( p7_profile_SetLength(gm, sq->n)       != eslOK) esl_fatal(msg);
  if ( p7_oprofile_ReconfigLength(om, sq->n) != eslOK) esl_fatal(msg);

  /* Fwd/Bck local filter to calculate the sparse mask */
  if (  (ox = p7_checkptmx_Create(M, sq->n, ESL_MBYTES(32)))             == NULL)  esl_fatal(msg);
  if (  (sm = p7_sparsemask_Create(M, sq->n))                            == NULL)  esl_fatal(msg);
  if ( p7_checkptmx_Reinit(ox, M, sq->n)                                 != eslOK) esl_fatal(msg);
  if ( p7_ForwardFilter (sq->dsq, sq->n, om, ox, /*fsc=*/NULL)           != eslOK) esl_fatal(msg);
  if ( p7_BackwardFilter(sq->dsq, sq->n, om, ox, sm, p7_SPARSIFY_THRESH) != eslOK) esl_fatal(msg);

  /* Sparse DP calculations, and exact posterior decoding */
  if ( p7_SparseForward   (sq->dsq, sq->n, gm, sm, sxf, NULL) != eslOK) esl_fatal(msg);
  if ( p7_SparseBackward  (sq->dsq, sq->n, gm, sm, sxd, NULL) != eslOK) esl_fatal(msg); /* Backwards mx temporarily in sxd... */
  if ( p7_SparseDecoding  (sq->dsq, sq->n, gm, sxf, sxd, sxd) != eslOK) esl_fatal(msg); /* followed by in-place decoding      */

  /* Approximate decoding by stochastic traceback  */
  if ( (sxs = p7_sparsemx_Create(sm)) == NULL) esl_fatal(msg);
  if ( (tr  = p7_trace_Create())      == NULL) esl_fatal(msg);
  if ( p7_sparsemx_Zero(sxs)         != eslOK) esl_fatal(msg);
  for (idx = 0; idx < ntr; idx++)
    {
      if ( p7_sparse_trace_Stochastic(rng, &wrk, gm, sxf, tr) != eslOK) esl_fatal(msg);

      //p7_trace_DumpAnnotated(stdout, tr, gm, sq->dsq);

      if ( p7_sparsemx_CountTrace(tr, sxs)                    != eslOK) esl_fatal(msg);
      if ( p7_trace_Reuse(tr)                                 != eslOK) esl_fatal(msg);
    }
  esl_vec_FScale(sxs->dp,   sxs->sm->ncells*p7S_NSCELLS,           1./(float)ntr);
  esl_vec_FScale(sxs->xmx, (sxs->sm->nrow+sxs->sm->S)*p7S_NXCELLS, 1./(float)ntr);

  //  p7_sparsemx_Dump(stdout, sxd);
  //  p7_sparsemx_Dump(stdout, sxs);

  /* Tests */
  if ( p7_sparsemx_CompareDecoding(sxd, sxs, tol) != eslOK) esl_fatal(msg);
  if ( p7_sparsemx_Validate(sxd, NULL)            != eslOK) esl_fatal(msg);
  if ( p7_sparsemx_Validate(sxs, NULL)            != eslOK) esl_fatal(msg);
  
  if (wrk) free(wrk);
  p7_sparsemx_Destroy(sxs);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemx_Destroy(sxb);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemask_Destroy(sm);
  p7_trace_Destroy(tr);
  p7_checkptmx_Destroy(ox);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  esl_sq_Destroy(sq);
}
#endif /*p7SPARSE_DECODING_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/



/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7SPARSE_DECODING_TESTDRIVE
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

/* Some tests involve stochastic sampling. Use a fixed RNG seed, to avoid
 * rare statistical excursions triggering rare failures at their expected
 * frequency.
 */
static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,     "20", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  { "-n",        eslARG_INT, "100000", NULL, NULL,  NULL,  NULL, NULL, "number of stochastic traces in approx_decoding", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for sparse posterior decoding, with dual-mode profile";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg   = p7_bg_Create(abc);
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");
  int             ntr  = esl_opt_GetInteger(go, "-n");

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  utest_overwrite      (r, abc, bg, M, L);
  utest_rowsum         (r, abc, bg, M, L, N);
  utest_approx_decoding(r, abc, bg, M, L, ntr);

  fprintf(stderr, "#  status = ok\n");

  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7SPARSE_DECODING_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/

