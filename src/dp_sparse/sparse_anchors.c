/* Most probable anchor set (MPAS) algorithm 
 * Sparse implementation: production code.
 * 
 * Contents:
 *     1. Sparse MPAS algorithm, fast version (segmental divide & conquer)
 *     2. Sparse MPAS algorithm, fallback version (global)
 *     3. Making trial anchor sets
 *     4. Footnotes 
 *     5. Statistics driver
 *     6. Example
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_random.h"

#include "base/p7_anchorhash.h"
#include "base/p7_anchors.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_sparse/p7_sparsemx.h"
#include "dp_sparse/p7_spascmx.h"
#include "dp_sparse/sparse_asc_fwdback.h"
#include "dp_sparse/sparse_trace.h"

#include "misc/logsum.h"

#include "search/p7_mpas.h"     // P7_MPAS_PARAMS and P7_MPAS_STATS

#include "dp_sparse/sparse_anchors.h"


/*****************************************************************
 * 1. MPAS algorithm, segmental divide & conquer version
 *****************************************************************/

static int segment_forward_score(const P7_PROFILE *gm, const P7_SPARSEMX *sxf, int g, float *dpc, float *xc, 
				 float **opt_dpn, float **opt_xn, float *ret_p0, float *ret_sF);
static int sparse_anchors_CatenateFromSegTrace(const P7_SPARSEMX *sxd, const P7_TRACE *tr, int g, float *dpc, P7_ANCHORS *anch, int D0, float **opt_dpn);



/* Function:  p7_sparse_Anchors()
 * Synopsis:  Most probable anchor set (MPAS) algorithm: sparse segmental version.
 *
 * Purpose:   Find the most probable anchor set for comparison of query
 *            profile <gm> to digital sequence <dsq> of length <L>.
 *            
 *            This is a segmental divide and conquer version of MPAS.
 *            It differs from <p7_sparse_AnchorsGlobal()> by building the
 *            optimal anchor set sequentially, one segment at a time,
 *            exploiting things we can prove from sparsity in the DP
 *            matrix. 
 *            
 *            Rarely, no individual segment wants to have an anchor,
 *            yielding a D=0 solution overall, which isn't possible; a
 *            profile must pass through the homology model at least
 *            once. In this case, the routine fails over to the global
 *            version of MPAS to find the correct solution.
 *
 *            MPAS algorithm uses stochastic path sampling, so caller
 *            provides random number generator <rng>.
 *            
 *            Caller has already done a so-called "First Analysis" on
 *            the sequence: (unconstrained) sparse Viterbi,
 *            Forward/Backward, and posterior Decoding. From this,
 *            caller gives us the Viterbi and Forward raw nat scores <vsc> and
 *            <fsc>; the sparse Forward and Decoding matrices <sxf>
 *            and <sxd>; and <vanch>, the anchor set implied by the
 *            Viterbi path. These remain constant during the MPAS
 *            calculation and are returned unchanged.
 *            
 *            Caller provides three empty workspaces. All of these
 *            will be reallocated as needed, so they can be of any
 *            initial allocation size. <tr> is a trace structure used
 *            for stochastic path samples. <byp_wrk> is an array using
 *            the Easel 'bypass' idiom, needed by the stochastic
 *            traceback algorithm. <ah> is an anchorset hash table,
 *            used for fast checking of which anchor set suffixes
 *            we've already calculated scores for. 
 *            
 *            For retrieving the result, caller provides an empty
 *            sparse matrix <asf>, an empty anchor set <anch>, and
 *            <ret_asc> for the ASC score. Upon return, <asf>
 *            contains the sparse ASC Forward matrix for the optimal
 *            anchor set; <anch> is that anchor set; and <asc_sc>
 *            is the raw ASC Forward score in nats.
 *            
 *            Optionally, caller may override default parameters by
 *            passing in a <prm> structure. Pass <NULL> to use
 *            defaults.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
p7_sparse_Anchors(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,
		  float vsc, float fsc, const P7_SPARSEMX *sxf, const P7_SPARSEMX *sxd, const P7_ANCHORS *vanch,
		  P7_TRACE *tr, float **byp_wrk, P7_ANCHORHASH *ah,
		  P7_SPARSEMX *asf, P7_ANCHORS *anch, float *ret_asc,
		  P7_MPAS_PARAMS *prm)
{
  int     max_iterations  = (prm ? prm->max_iterations      : 1000);
  float   loglossthresh   = (prm ? log(prm->loss_threshold) : log (0.001));
  const P7_SPARSEMASK *sm = sxf->sm; // sparse mask: sxf, sxd, and asf all have links to the same mask
  float  *dpc_f = sxf->dp;           // strides thru <sxf> segment by segment: ptr to next sparse DP cell at ia(g) (=ia(g)-1, because no main cells at ia(g)-1)
  float  *xc_f  = sxf->xmx;          //  ... ditto. ptr to next special DP cell at ia(g)-1.
  float  *dpn_f = NULL;              // will get set (by segment_forward_score) to next main DP cell > ib(g) cells; i.e. ia(g+1)
  float  *xn_f  = NULL;              //  ... ditto; next special DP cell > ib(g), i.e. ia(g+1)-1
  float  *dpc_d = sxd->dp;           // keeps state in <sxd> decoding matrix; ptr to next DP cell at ia(g)
  float  *dpn_d = NULL;              // will get set (by CatenateFromSegTrace) to next main DP cell > ib(g); i.e. ia(g+1)
  float  *dpc_a = NULL;              // keeps current state in <asf->dp>
  float  *xc_a  = NULL;              //                ... and <asf->xmx>
  float  *dpn_a = NULL;              // gets set to next state of <asf->dp> 
  float  *xn_a  = NULL;              //                   ... and <asf->xmx>
  int     D0    = 0;                 // <anch> prefix optimized so far, for segments 1..g-1
  int     D0v   = 0;                 // <vanch> prefix used so far, for segments 1..g-1
  float   xNc   = 0.;                // xN score at ib(g-1), last row of prev segment; starts at 0.
  float   xJc   = -eslINFINITY;      //   ... and xJ
  float   xCc   = -eslINFINITY;      //   ... and xC
  float   xNn, xJn, xCn;             // xN,xJ,xC scores at ib(g) of this segment, after a sparse_asc_ForwardSeg() calculation.
  float   sF;			     // segment Forward score for current segment
  float   p0;                        // segment score of the empty path from ia(g)-1..ib(g); sF sum includes this path, sF = logsum(deltaS + p0)
  int     g;                         // counter over segments, 1..sm->S
  int     Dg;                        // number of anchors in <anch> for this segment g; total D in anch = D0 prefix + Dg suffix
  int     iteration;                 // counter over iterations. 0=empty set; 1=Viterbi set; 2 and up are sets from sampled paths
  int     i;                         // counter over sequence residues 1..L
  float   asc, best_ascprob;
  int32_t keyidx, asc_keyidx, best_keyidx;
  float   best_asc;
  int     ngap;
  int     status;

  /* Contract checks and argument validation */

  ESL_DASSERT1(( gm->xsc[p7P_N][p7P_LOOP] == gm->xsc[p7P_J][p7P_LOOP] ));  // segment Fwd score depends on assumption tNN=tJJ=tCC
  ESL_DASSERT1(( gm->xsc[p7P_N][p7P_LOOP] == gm->xsc[p7P_C][p7P_LOOP] ));
  ESL_DASSERT1(( sxf->type == p7S_FORWARD  ));
  ESL_DASSERT1(( sxd->type == p7S_DECODING ));
  ESL_DASSERT1(( sm->L == L && sm->M == gm->M ));
  ESL_DASSERT1(( vanch->D  >  0 ));                                        // <vanch> is anchorset implied by Viterbi path
  ESL_DASSERT1(( anch->D   == 0 ));                                        // <anch> is provided empty
  ESL_DASSERT1(( tr->N     == 0 ));                                        //   ... so it <tr> 
  ESL_DASSERT1(( ah->nkeys == 0 ));                                        //   ... so is <ah>

  /* If we're lucky, we can immediately prove that the Viterbi-implied
   * anchor set has to be the MPAS, in which case we're immediately
   * done.
   */
  if (exp(vsc - fsc) > 0.5)
    {
      p7_anchors_Copy(vanch, anch);
      p7_sparse_asc_Forward(dsq, L, gm, anch->a, anch->D, sm, asf, ret_asc);  // Determines overall <asc>; fills ASC Fwd matrix <asf>. <asf> is resized by sparse_asc_Forward as needed.
      return eslOK;
    }
  /* Otherwise, on to the segmental MPAS algorithm. */



  /* Reallocate <asf> sufficient to work with any anchor set. 
   * We may revisit this allocation strategy. This allocates
   * the upper bound, which is 2x the sparse mask in main states,
   * and 1x in specials. Average-case usage in sparse ASC
   * for typical anchor sets is more like ~0.9x/0.9x, I 
   * believe, so we're overallocating by ~2x. That's not too
   * much of a big deal because we're likely going to see
   * a larger profile/seq comparison in a future call where
   * we're reusing this DP matrix. It only becomes a serious
   * concern when we're seeing an outlier comparison size, like
   * titin x titin.
   */
  p7_spascmx_Resize(asf, sm, NULL, 0);
  asf->type = p7S_ASC_FWD;
  dpc_a     = asf->dp;
  xc_a      = asf->xmx;

  /* anchor set in <anch> must be valid even when empty:
   * must set its sentinels at start. They will be copied
   * and moved as needed, as <anch> grows.
   */
  p7_anchor_SetSentinels(anch->a, 0, L, gm->M);



  /*****************************************************************
   * Start of the segmental MPAS algorithm.
   * We optimize each segment's anchor set separately.
   * Overall MPAS is the concantenation of these.
   * At any time, <anch> is a current guess at the MPAS for
   * segments 1..g, consisting of the MPAS prefix 1..g-1 and
   * a current guess at a possible optimum for g.
   */
  for (g = 1; g <= sm->S; g++)
    {
      /* Get the "segment forward score" for this segment.
       * sF is sum of all paths (including D=0) from ia(g)-1 to ib(g)
       * for any starting special X ({N|J}) to any ending special Y
       * ({N|J|C}). Assumes tNN=tJJ=tCC. sF is the denominator for
       * calculating posterior prob of an anchor set.
       * 
       * This is also the best place to figure out <dpn_f>, <xn_f>:
       * location of DP ptrs after this segment, next position >ib(g).
       */
      segment_forward_score(gm, sxf, g, dpc_f, xc_f, &dpn_f, &xn_f, &p0, &sF);
      dpc_f = dpn_f - sm->n[sm->seg[g].ib] * p7S_NSCELLS;  // Now dpc_f is on start of row ib(g) - where stochastic tracebacks need it to be
      xc_f  = xn_f  - p7S_NXCELLS;                         //  ... ditto for xc.
      dpn_d = NULL;                                        // <dpn_d> in <sxd> will only get set properly if we do at least one stochastic trace ... which we may not do.


      /*****************************************************************
       * Start of the main MPAS optimization loop for one segment g
       */
      iteration = 0;
      best_asc  = -eslINFINITY;
      while (1)
	{
	  ESL_DASSERT1(( xc_a != NULL && dpc_a != NULL));

	  /* First step is to propose an anchor set suffix for 
	   * segment <g>, and append it to the MPAS solution so far
	   * in <anch>. There's three ways to propose this suffix:
	   *   1. Anchor set is empty.
           *   2. Anchor set is the Viterbi-implied set.
	   *   3. Do a stochastic segment traceback and get anchor set from that.
           *
           * <anch> will have D anchors total:
           *    1..D0 are the MPAS solution for segments 1..g-1
           *    D0+1..D0+Dg is the proposed solution for segment g
	   */
	  if      (iteration == 0)                                            // 1st guess: empty anchor set for this segment; D=D0. <anch> stays as it is.
	    {                                                                 
	      Dg = 0;              
	    }
	  else if (iteration == 1)                                            // 2nd guess: Viterbi path for this segment yields optimal anchor set.
	    {                                                                 // (Even if V segment score alone is enough to prove that it's the seg's MPAS, we do still need ASC Fwd matrix across the segment.)
	      for (Dg = 0, i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)        // Figure out how many anchors are in the Viterbi path for this segment
		if (i == vanch->a[D0v+Dg+1].i0) Dg++;
	      p7_anchors_Catenate(anch, D0, vanch->a + D0v, Dg);	      // anch [ 1..D0 ] gets vanch [D0v+1..D0v+Dg+1] appended to it; <anch> now has D = D0+Dg
	      D0v += Dg;
	    }
	  else                                                                // Else: all remaining tries take a stochastic trace for this segment.
	    {
	      p7_sparse_trace_StochasticSeg(rng, byp_wrk, gm, sxf, p0, sF, g, dpc_f, xc_f, tr);   // returned tr->N can be 0 here, in the special case of the empty path being sampled.
	      sparse_anchors_CatenateFromSegTrace(sxd, tr, g, dpc_d, anch, D0, &dpn_d);           // dpn_d moves to start of seg g+1 in <sxd>
	      Dg = anch->D - D0;                                                                  // Dg may be 0, in the case of the empty trace tr->N=0
	    }
	    
	  /* Step 2: Now use the anchor hash table to see if we've
           * already seen this anchor suffix. If we've already seen
	   * it, we don't need to waste time calculating its ASC 
	   * score. <eslOK> status means it's new.
	   * <eslEDUP> means we've already seen it. Dg=0 (empty path,
           * no anchors) is special cased for a little extra efficiency:
           * we know that we see it on iteration 0, so thereafter Dg=0
	   * must be a dup.
	   */
	  if (Dg > 0)   status = p7_anchorhash_Store(ah, anch, D0, &keyidx);
	  else        { status = (iteration == 0 ? eslOK : eslEDUP);         keyidx = -1; }
	  

	  /* Step 3. If it's the first time we've seen this proposed
	   * anchor set for this segment, then we need to calculate
	   * its ASC Forward score. We also need to update the state
	   * information that we'll use to advance to the next segment
	   * g+1 when we're all done with g. 
	   */
	  if (status == eslOK)
	    {            
	      ngap  = (g == 1 ? sm->seg[g].ia - 1 : sm->seg[g].ia - sm->seg[g-1].ib - 1); // Length of intersegment space ib(g-1)+1..ia(g)-1. First segment must be special cased because sentinel ib(0) is -1, not 0
	      if (Dg == 0) 
		{ 
		  xNn   = xNc + p0 + (ngap ? (float) ngap * gm->xsc[p7P_N][p7P_LOOP] : 0.);   // Update state of xC/xN/xJ by adding score for crossing intersegment space,
		  xJn   = xJc + p0 + (ngap ? (float) ngap * gm->xsc[p7P_J][p7P_LOOP] : 0.);   //   and score for crossing segment <g>. Numerical error must match other sparse DP 
		  xCn   = xCc + p0 + (ngap ? (float) ngap * gm->xsc[p7P_C][p7P_LOOP] : 0.);   //   algorithms: interseg is done by multiplication, within-seg by summation
		  dpn_a = dpc_a; 
		  xn_a  = xc_a;
		  asc   = p0; 
		}
	      else 
		p7_sparse_asc_ForwardSeg(dsq, L, gm, anch->a, anch->D, D0+1, sm, ngap, sm->seg[g].ia, sm->seg[g].ib, 
					 asf,                // <asf> has already been reallocated. We must not Reuse() it; it contains a growing prefix we need to keep!
					 xNc, xJc, xCc, dpc_a, xc_a,
					 &xNn, &xJn, &xCn, /*opt_d:*/NULL, &dpn_a, &xn_a, &asc);
	      asc_keyidx = keyidx;     // asc_keyidx remembers which <keyidx> corresponds to <asf>
	      
	      if (asc > best_asc) {
		best_asc     = asc;
		best_ascprob = exp(asc - sF);   // sF segment score, not fsc overall score
		best_keyidx  = keyidx;
	      }
	    }
	  
	  p7_trace_Reuse(tr);

	  /* Termination conditions */
	  if (best_ascprob >= 0.5) break;
	  if ((iteration - 1) * log(1.0 - best_ascprob) < loglossthresh) break;
	  if (iteration == max_iterations) break;
     
	  iteration++;
	} // end of the while (1) loop that finds MPAS for one segment <g>

      if (keyidx != best_keyidx)                               // ... then <anch> does not currently have the optimal suffix; fetch it back.
	{
	  if   (best_keyidx == -1)                                  // if best_keyidx is -1, that's the Dg=0 empty set; truncate back to D0
	    {
	      anch->D = D0;                                         // ... so set D to what we had up to the previous seg g
	      p7_anchor_SetSentinels(anch->a, anch->D, L, gm->M);   // ... remember to reset sentinels! if last anch was for Dg>0
	    }
	  else p7_anchorhash_Get(ah, best_keyidx, D0, anch);   // else, get a suffix Dg>0. _Get resets sentinels
	}
      if (asc_keyidx != best_keyidx)                           // ... then <asf> does not currently have the DP matrix for the optimal suffix; recalculate it
	{                                                      //     also, the state information associated <asf> needs to be recalculated.

	  if (best_keyidx == -1)                               // Is the optimum the empty set? Then we don't need to run an ASC Forward calc on segment g;
	    {                                                  //   instead, we know all the relevant information to update state and set asc.
	      asc   = p0;
	      ngap  = (g == 1 ? sm->seg[g].ia - 1 : sm->seg[g].ia - sm->seg[g-1].ib - 1); // Length of intersegment space ib(g-1)+1..ia(g)-1. First segment must be special cased because sentinel ib(0) is -1, not 0
	      xNn   = xNc + p0 + (ngap ? (float) ngap * gm->xsc[p7P_N][p7P_LOOP] : 0.);   // Update state of xC/xN/xJ by adding score for crossing intersegment space,
	      xJn   = xJc + p0 + (ngap ? (float) ngap * gm->xsc[p7P_J][p7P_LOOP] : 0.);   //   and score for crossing segment <g>. Numerical error must match other sparse DP 
	      xCn   = xCc + p0 + (ngap ? (float) ngap * gm->xsc[p7P_C][p7P_LOOP] : 0.);   //   algorithms: interseg is done by multiplication, within-seg by summation
	      dpn_a = dpc_a;                                                              // Update the DP ptrs, so they point to where they were for g; they don't
	      xn_a  = xc_a;                                                               //   move, because there's no storage for segment g if we use the empty path p0
	    }
	  else  // but if the optimum is Dg>0, then we do need an ASC Fwd calculation for this segment.
	    {
	      ngap  = (g == 1 ? sm->seg[g].ia - 1 : sm->seg[g].ia - sm->seg[g-1].ib - 1);
	      p7_sparse_asc_ForwardSeg(dsq, L, gm, anch->a, anch->D, D0+1, sm, ngap, sm->seg[g].ia, sm->seg[g].ib, asf, 
				       xNc, xJc, xCc, dpc_a, xc_a, 
				       &xNn, &xJn, &xCn, /*opt_d:*/NULL, &dpn_a, &xn_a, &asc);
	    }
	}

      /* We've finished segment <g>;
       *   <anch> contains the optimal anchor set for 1..g;
       *   <asf> contains the sparse ASC matrix for <anch> for segments 1..g.
       * Now bump the state information along, initializing it for the next segment g+1.
       */
      dpc_f = dpn_f;
      xc_f  = xn_f;

      if (dpn_d)	// <dpn> gets set if we did a stochastic trace
	dpc_d = dpn_d;
      else              //  ... but we may not have, in which case we have to bump it manually.
	for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++) 
	  dpc_d += sm->n[i] * p7S_NSCELLS;

      dpc_a = dpn_a;
      xc_a  = xn_a;
      xNc   = xNn;
      xJc   = xJn;
      xCc   = xCn;
      D0    = anch->D;
    } // end loop over all segments <g>

  /* It's possible for segmental D&C MPAS to find a D=0 solution,
   * where no individual segment wants to have any anchors.  A D=0
   * solution is not allowed, because the profile must have at
   * least one pass thru the homology model. Therefore D=0 is a
   * special case failure mode; we have to fail over to the
   * slower, global version of MPAS.
   */
  if (anch->D == 0)
    {
      p7_trace_Reuse(tr);
      p7_anchorhash_Reuse(ah);
      p7_sparsemx_Reuse(asf);
      p7_anchors_Reuse(anch);
      
      return p7_sparse_AnchorsGlobal(rng, dsq, L, gm,
				     vsc, fsc, sxf, sxd, vanch,
				     tr, byp_wrk, ah, 
				     asf, anch, ret_asc,
				     prm);
    }
  else
    {
      if (L > sm->seg[sm->S].ib) 
	xCc +=  (float) (L - sm->seg[sm->S].ib) * gm->xsc[p7P_C][p7P_LOOP];
      *ret_asc = xCc +  gm->xsc[p7P_C][p7P_MOVE];
      return eslOK;
    }
}



/*****************************************************************
 * 2. MPAS algorithm, complete sequence version
 *****************************************************************/

/* Function:  p7_sparse_AnchorsGlobal()
 * Synopsis:  The most probable anchor set (MPAS) algorithm; sparse, complete-seq version.
 *
 * Purpose:   Find the most probable anchor set for comparison of query profile
 *            <gm> to digital sequence <dsq> of length <L>. 
 *
 *            This version of sparse MPAS solves for the optimum
 *            anchor set of the complete sequence as a whole.
 *            Normally we use <p7_sparse_Anchors()>, which uses a
 *            faster divide and conquer algorithm, but that algorithm
 *            has a (rare) failure mode in which it can produce a
 *            <D=0> 'solution', in which case we need to fail over to
 *            this algorithm.
 *            
 *            The MPAS algorithm depends on stochastic path sampling,
 *            so caller provides a random number generator <rng>.
 *            
 *            Caller has already done what we call "First Analysis" on
 *            the sequence: standard Viterbi, Forward/Backward, and
 *            posterior Decoding. From this, caller gives us the
 *            Viterbi and Forward raw nat scores <vsc> and <fsc>; the
 *            sparse Forward and Decoding matrices <sxf> and <sxd>;
 *            and <vanch>, the anchorset implied by the Viterbi path.
 *            These remain constant during the MPAS calculations and
 *            are returned unchanged.
 *            
 *            Caller provides three empty workspaces. All are
 *            reallocated as needed, so they can be of any initial
 *            allocation size. <tr> is a trace structure used for
 *            stochastic path samples. <byp_wrk> is a "bypass"
 *            workspace that the stochastic trace algorithm knows how
 *            to deal with.  <ah> is an empty hash table that we need
 *            for fast checking of which anchor sets we've already
 *            calculated scores for.
 *            
 *            For retrieving the result, caller provides an empty
 *            sparse matrix <asf>, an empty anchor set <anch>, and
 *            <ret_asc> for the ASC score.  Upon return, <asf>
 *            contains the sparse ASC Forward matrix for the optimal
 *            anchor set; <anch> is that anchor set; and <asc_sc> is
 *            the raw ASC Forward score in nats.
 *            
 *            Optionally, caller may override default parameters by
 *            passing in a <prm> structure. Pass <NULL> to use
 *            defaults.
 *
 * Args:      rng     : random number generator
 *            dsq     : target sequence, digital (1..L)
 *            L       : length of <dsq>
 *            gm      : query profile, dual-mode (1..M)
 *            vsc     : sparse Viterbi score (raw, nats)
 *            fsc     : sparse Forward score
 *            sxf     : sparse Forward matrix calculated for <gm>/<dsq> comparison
 *            sxd     : sparse Decoding matrix calculated for <gm>/<dsq> comparison
 *            vanch   : anchorset implied by the Viterbi trace
 *            tr      : empty traceback; used as space for sampled traces.
 *            byp_wrk : BYPASS: workspace array for stochastic tracebacks
 *            ah      : empty hash table for alternative <anch> solutions
 *            asf     : empty sparse matrix, to be used for sparse ASC Forward calculations
 *            anch    : empty anchor data structure to hold the result, the MPAS
 *            ret_asc : score of most probable anchor set (raw, nats)
 *            prm     : OPTIONAL: non-default control parameters
 *
 * Returns:   <eslOK> on success, and:
 *            anch    : contains most probable anchor set
 *            asf     : contains sparse ASC Forward matrix for MPAS solution
 *            ret_asc : ASC forward score of it, raw (nats)
 *            
 */
int
p7_sparse_AnchorsGlobal(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,
			float vsc, float fsc, const P7_SPARSEMX *sxf, const P7_SPARSEMX *sxd, const P7_ANCHORS *vanch,
			P7_TRACE *tr, float **byp_wrk, P7_ANCHORHASH *ah,
			P7_SPARSEMX *asf, P7_ANCHORS *anch, float *ret_asc,
			P7_MPAS_PARAMS *prm)
{
  int     max_iterations = (prm ? prm->max_iterations      : 1000);
  float   loglossthresh  = (prm ? log(prm->loss_threshold) : log (0.001));
  float   best_asc       = -eslINFINITY;
  int     iteration      = 0;
  float   asc, best_ascprob;
  int32_t keyidx, asc_keyidx, best_keyidx;
  int     status;

  /* Contract checks, argument validation */
  ESL_DASSERT1(( sxf->type == p7S_FORWARD  ));
  ESL_DASSERT1(( sxd->type == p7S_DECODING ));
  ESL_DASSERT1(( sxf->sm->L == L && sxf->sm->M == gm->M ));
  ESL_DASSERT1(( sxd->sm->L == L && sxd->sm->M == gm->M ));
  ESL_DASSERT1(( vanch->D  >  0 ));                                        // <vanch> is anchorset implied by Viterbi path
  ESL_DASSERT1(( anch->D   == 0 ));                                        // <anch> is provided empty
  ESL_DASSERT1(( tr->N     == 0 ));                                        //   ... so it <tr> 
  ESL_DASSERT1(( ah->nkeys == 0 ));                                        //   ... so is <ah>

  /* anchorset implied by Viterbi path is our initial guess */
  p7_anchors_Copy(vanch, anch);

  while (1)
    {
      status = p7_anchorhash_Store(ah, anch, 0, &keyidx);

      if (status == eslOK)
	{
	  p7_sparsemx_Reuse(asf);
	  p7_sparse_asc_Forward(dsq, L, gm, anch->a, anch->D, sxf->sm, asf, &asc);
	  asc_keyidx = keyidx;   // asc_keyidx remembers which <keyidx> corresponds to <asf>

	  if (asc > best_asc) {
	    best_asc     = asc;
	    best_ascprob = exp(asc - fsc);
	    best_keyidx  = keyidx;
	  }

	} 

      if (best_ascprob >= 0.5)                                 break;
      if (iteration * log(1.0 - best_ascprob) < loglossthresh) break;
      if (iteration == max_iterations)                         break;

      p7_trace_Reuse(tr);
      p7_anchors_Reuse(anch);
      p7_sparse_trace_Stochastic(rng, byp_wrk, gm, sxf, tr);
      p7_sparse_anchors_SetFromTrace(sxd, tr, anch);

      iteration++;
    }

  if (keyidx     != best_keyidx) p7_anchorhash_Get(ah, best_keyidx, 0, anch);
  if (asc_keyidx != best_keyidx) p7_sparse_asc_Forward(dsq, L, gm, anch->a, anch->D, sxf->sm, asf, &asc);

  *ret_asc = asc;
  return eslOK;
}





/*****************************************************************
 * 3. Making trial anchor sets
 *****************************************************************/

/* Function:  p7_sparse_anchors_SetFromTrace()
 * Synopsis:  Select an anchor set from a path.
 * 
 * Purpose:   Select an anchor set from a path, by choosing the highest
 *            posterior probability match i,k cell (marginalizing over 
 *            local/glocal). 
 *            
 *            Caller provides a sparse posterior decoding matrix <sxd>
 *            and the path <tr>. Resulting anchor set is stored in
 *            <anch>, which will be reallocated as needed. Caller 
 *            provides existing but empty <anch> object.
 *            
 *            The trace must be compatible with the sparse mask in
 *            <sxd>: every i,k used in the trace must be in the sparse
 *            mask.  (This is automatically true if the trace comes
 *            from sparse viterbi or sparse stochastic traceback with
 *            the same mask.)
 *
 * Args:      sxd  : sparse decoding matrix
 *            tr   : path
 *            anch : resulting anchor set selected from <tr>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 *
 * Note:      Similar to sparse_anchors_SetFromSegmentTrace().           
 */
int
p7_sparse_anchors_SetFromTrace(const P7_SPARSEMX *sxd, const P7_TRACE *tr, P7_ANCHORS *anch)
{
  const P7_SPARSEMASK *sm  = sxd->sm;
  float               *dpc = sxd->dp;      // <dpc> is on start of sparse rows i
  float               *dpc2;               // <dpc2> steps across sparse cells k on a given row i
  int   z;                                 // counter in trace 0..N-1
  int   i;                                 // sequence position 1..L
  int   y;                                 // counter in sparse cells on a row, 0..sm->n[i]-1
  float ppv;                               // posterior probability for M(i,k)
  float best_ppv = -1.;                    // best ppv seen so far for domain anch->D+1
  int   status;

  ESL_DASSERT1(( sxd->type == p7S_DECODING ));
  ESL_DASSERT1(( anch->D   == 0 ));

  for (z = 0; z < tr->N; z++)
    {
      i = tr->i[z];

      if (p7_trace_IsM(tr->st[z]))
	{
	  y    = 0;
	  dpc2 = dpc;
	  while (y < sm->n[i] && sm->k[i][y] < tr->k[z]) { y++; dpc2 += p7S_NSCELLS; }
	  ESL_DASSERT1(( sm->k[i][y] == tr->k[z] ));  // all i,k in trace must be present in sparse mask
	  ppv = dpc2[p7S_ML] + dpc2[p7S_MG];          // posterior probability at match(i,k) marginalized over local/glocal
	  //printf("                 i %d k %d ppv %.4f offset %d\n", tr->i[z], tr->k[z], ppv, (dpc2 - sxd->dp) / 6);
	  if (ppv > best_ppv) 
	    {
	      anch->a[anch->D+1].i0 = tr->i[z];        
	      anch->a[anch->D+1].k0 = tr->k[z];
	      best_ppv = ppv;
	    }
	}
      else if (tr->st[z] == p7T_E)
	{
	  anch->D++;
	  if ((status = p7_anchors_Resize(anch, anch->D+1)) != eslOK) goto ERROR;
	  best_ppv = -1.;
	}
  
      if (i) dpc += p7S_NSCELLS * sm->n[i];
    }
  p7_anchor_SetSentinels(anch->a, anch->D, tr->L, tr->M);
  return eslOK;

 ERROR:
  return status;
}



/* sparse_anchors_CatenateFromSegTrace()
 *
 * Purpose:   Given a segment trace <tr> for segment <g> in sparse
 *            decoding matrix <sxd>; for every domain in the segment
 *            trace, choose the best anchor <i0,k0> by choosing the
 *            match state (ML+MG, marginalized posterior) with highest
 *            posterior probability. Store the anchors in order in
 *            <anch>, an allocated and empty structure provided by the
 *            caller.
 *            
 *            The caller may also provide <dpc>, optionally. This is a
 *            pointer into the sparse DP matrix <sxd>, to the start of
 *            the main matrix cells for row <ia>. The caller may know
 *            this location, if it is traversing the segments
 *            sequentially g=1..G. Otherwise, caller passes
 *            <dpc=NULL>, and it will be set internally the slow way.
 *            
 *            Upon return, in addition to <anch> being complete with
 *            <anch->D> anchors, we also optionally return a pointer <opt_dpn>:
 *            this is a ptr to the next segment's main matrix row ia.
 *            Caller can use this to sequentially traverse segments.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 *
 * Xref:      See footnote [4] for more info on the "segment trace" hack
 *            of P7_TRACE structure.
 */
static int
sparse_anchors_CatenateFromSegTrace(const P7_SPARSEMX *sxd, const P7_TRACE *tr, int g, float *dpc, P7_ANCHORS *anch, int D0, float **opt_dpn)
{
  const P7_SPARSEMASK *sm = sxd->sm;
  int           *kc;
  float         *dpc2;
  float          ppv;
  float          best_ppv = -1.;
  int            i;
  int            z;
  int            y;
  int            status;

  ESL_DASSERT1(( sxd->type == p7S_DECODING ));

  /* If we don't have <dpc>, find it the slow way */
  if (dpc == NULL)
    {
      dpc = sxd->dp;
      for (i = 1; i < sm->seg[g].ia; i++) 
	dpc += p7S_NSCELLS * sm->n[i];
    }

  anch->D = D0;                 // That's sufficient to remove the previous suffix.

  for (z = 0; z < tr->N; z++)   // If N=0 (empty trace, empty path, no anchors), then this loop no-ops, which is correct.
    {
      if (p7_trace_IsM(tr->st[z]))
	{
	  i    = tr->i[z];
	  y    = 0; 
	  kc   = sm->k[i];
	  dpc2 = dpc;
	  while (y < sm->n[i] && kc[y]  < tr->k[z]) { y++; dpc2 += p7S_NSCELLS; }
	  if    (y < sm->n[i] && kc[y] == tr->k[z]) {
	    ppv = dpc2[p7S_ML] + dpc2[p7S_MG];
	    //printf("                 i %d k %d ppv %.4f offset %d\n", tr->i[z], tr->k[z], ppv, (dpc2 - sxd->dp) / 6);
	    if (ppv > best_ppv)

	      {
		anch->a[anch->D+1].i0 = tr->i[z];  
		anch->a[anch->D+1].k0 = tr->k[z];
		best_ppv = ppv;
		//printf("new best anchor. i %d k %d ppv %.4f\n", tr->i[z], tr->k[z], ppv);
	      }
	  }
	}
      else if (tr->st[z] == p7T_E)
	{
	  anch->D++;
	  if ((status = p7_anchors_Resize(anch, anch->D+1)) != eslOK) goto ERROR;
	  best_ppv = -1.;
	}

      if (tr->i[z]) {
	//printf("%4d  offset %d; advancing by %d\n", tr->i[z], (dpc - sxd->dp) / 6, sm->n[tr->i[z]]);
	dpc += p7S_NSCELLS * sm->n[tr->i[z]];
      }
    }

  p7_anchor_SetSentinels(anch->a, anch->D, tr->L, tr->M);
  *opt_dpn = dpc;
  return eslOK;

 ERROR:
  return status;
}


/* segment_forward_score()
 *
 * Purpose:   See footnote [1]. The "segment forward score" is the score
 *            of all paths that account for the residues in segment <g>,
 *            starting from N|J|C state at ia-1 and ending at N|J|C state
 *            at ib. Caller provides the profile <gm> and the sparse Forward
 *            matrix <sxf>. 
 *
 *            This calculation requires that t_NN = t_JJ = t_CC in the
 *            profile.
 *            
 *            Segment Forward score for segment <g> is returned thru
 *            <ret_sF>. The empty path component of that sum is
 *            returned thru <p0>.

 *            If caller knows the current "state" of a sequential
 *            segment-based traversal of the <sxf> matrix... i.e., if
 *            it knows valid ptrs <dpc> and <xc> to the main and
 *            special rows for <seg[g].ia-1> in <sxf>, it should
 *            provide these. Else, pass <NULL>, and the needed
 *            pointers will be determined by a somewhat more laborious
 *            means. (Sparse matrices are designed for sequential, not
 *            random access.)
 *            
 *            To support sequential access patterns, the next state
 *            ptrs (pointing to the start of the next segment,
 *            i.e. just past the end of this segment) are optionally
 *            returned thru <opt_dpn>, <opt_xn>.
 *            
 *            <p0> is calculated as \sum_n t_XX for n residues in the
 *            segment, not n * t_XX, in order to match numerical error
 *            in other DP routines.
 *
 * Returns:   <eslOK> on success.
 */
static int
segment_forward_score(const P7_PROFILE *gm, const P7_SPARSEMX *sxf, int g, float *dpc, float *xc, 
		      float **opt_dpn, float **opt_xn, float *ret_p0, float *ret_sF)

{
  const P7_SPARSEMASK *sm  = sxf->sm;
  float Jia1;			// Forward score in <sx> at J(ia-1)       
  float Nia1;			// Forward score in <sx> at N(ia-1)       
  float Jib;			// Forward score in <sx> at J(ib)         
  int   Ls;			// Length of segment g in residues        
  float p0;			// Score term for D=0 path from ia-1..ib  
  float deltaC;			// Intermediate term in calculating sF(g); see footnote [2] 
  float deltaS;			//   ... another intermediate term; see footnote [2]        
  float term2;			//   ... another intermediate term; see footnote [2]        
  float *dpn;
  int   g2,i;
  
  /* Contract checks.
   * Using segmented MPAS algorithm optimization requires t_NN = t_JJ = t_CC;
   * see footnote [1].
   */
  ESL_DASSERT1(( gm->xsc[p7P_N][p7P_LOOP] == gm->xsc[p7P_J][p7P_LOOP] ));
  ESL_DASSERT1(( gm->xsc[p7P_N][p7P_LOOP] == gm->xsc[p7P_C][p7P_LOOP] ));
  ESL_DASSERT1(( sxf->type == p7S_FORWARD ));

  /* If caller doesn't have <dpc>,<xc> positioned already on ia-1, 
   * it passes dpc=xc=NULL, and we set them the slow way.
   */
  if (xc == NULL)
    {
      xc = sxf->xmx;
      for (g2 = 1; g2 < g; g++)
	xc += p7S_NXCELLS * (sm->seg[g2].ib - sm->seg[g2].ia+2);  // ia-1..ib = ib-ia+2 to skip segment
    }
  if (dpc == NULL)
    {
      dpc = sxf->dp;
      for (g2 = 1; g2 < g; g++)
	for (i = sm->seg[g2].ia; i <= sm->seg[g2].ib; i++)
	  dpc += p7S_NSCELLS * sm->n[i];
    }

  /* p0, alas, must be determined by summation, not multiplication,
   * because we need to match summation roundoff error accumulation
   * in other DP routines.
   * 
   * While we're at it, figure out where <dpn> needs to be: ptr
   * past the last row ib(g), i.e. at start of next segment ia(g+1).
   */
  p0  = 0.;
  dpn = dpc;
  for (i = sm->seg[g].ia; i <= sm->seg[g].ib; i++)
    {
      p0  += gm->xsc[p7P_J][p7P_LOOP];
      dpn += sm->n[i] * p7S_NSCELLS;
    }

  /* Finally, we can do the algebra to calculate deltaS
   * (sum of all paths with D>=1 domains) from the sparse
   * Forward scores
   */
  Ls     = sm->seg[g].ib - sm->seg[g].ia + 1;
  Jia1   = xc[p7S_J];
  Nia1   = xc[p7S_N];
  xc    += p7S_NXCELLS * Ls;  // now <xc> is on <ib>
  Jib    = xc[p7S_J];
  deltaC = Jib - Jia1 - p0;
	 
  /* See footnote [3], numerical issues with calculating deltaS from deltaC: */
  if      (deltaC < eslSMALLX1)            term2 = log(deltaC);
  else if (deltaC > -1. * log(eslSMALLX1)) term2 = -1. * exp(-1. * deltaC);
  else                                     term2 = log(1. - exp(-1. * deltaC));

  deltaS  = Jib + term2 - p7_FLogsum(Nia1, Jia1);

  if (opt_dpn) *opt_dpn = dpn;
  if (opt_xn)  *opt_xn  = xc + p7S_NXCELLS;  // because xc was on ib(g), and xn needs to be ib(g)+1 i.e. ia(g+1)-1
  *ret_p0 = p0;
  *ret_sF = p7_FLogsum(deltaS, p0);
  return eslOK;
}


/*****************************************************************
 * 4. Footnotes
 *****************************************************************
 *
 *
 * [1] SEGMENT FORWARD SCORES
 * 
 * Given a segment ia..ib, define the "segment Forward score" s_F as
 * the log-sum of all paths that get us from the N|J|C state at ia-1
 * to the N|J|C state at ib (having generated residues ia..ib),
 * including the no-homology path (N...N, J...J, or C...C).
 * 
 * s_F is independent of the identity of the starting and ending state
 * (N|J|C at ia-1, ib) in the special case of t_NN = t_JJ = t_CC,
 * which is the default configuration of HMMER search profiles.  In
 * this special case, s_F can be calculated independently for each
 * segment, because s_F is independent of what happens in other
 * segments. This allows the MPAS algorithm to solve a MPAS on each
 * independent segment, one at a time, which can be a considerable
 * speedup.
 * 
 * If t_NN, t_JJ, t_CC are not equal, segment Forward scores are not
 * independent, and this partitioning cannot be applied correctly. In
 * that case, the complete-sequence MPAS algorithm must be used.
 * 
 * Another complication: it is possible for every individual segment
 * to have an MPAS solution of "no anchor" (because we consider the
 * empty NNN/JJJ/CCC path across it), but an anchor set for the
 * overall sequence must have D>=1 anchors. If the segmental MPAS
 * algorithm leads to a D=0 solution, again we must fall back on a
 * complete-sequence MPAS algorithm.
 * 
 *
 * 
 *****************************************************************
 ***************************************************************** 
 * 
 * [2] CALCULATING SEGMENT FORWARD SCORE FROM SPARSE FWD MATRIX
 * 
 * The segment may contain D>=0 domains (anchors). Let t_XX = t_NN =
 * t_CC = t_JJ.
 * 
 * If there is no domain in the segment, there's only a single path
 * with log probability (ib-ia+1) t_XX. Call this score p0.
 *
 * If there is one or more domain in the segment, define "deltaS" as
 * the log-sum of all paths that get from state X(ia-1) to Y(ib),
 * where X can be N or J (depending on whether ia..ib contains the
 * first domain in the sequence or not) and Y can be J or C (depending
 * on whether ia..ib contains the last domain in the sequence of not).
 * 
 * Note that deltaS is exclusive of the empty D=0 path, but contains
 * all other paths that get us from X(ia-1) to Y(ib), so:
 * 
 *    s_F = logsum(deltaS, p0)
 *    
 * Now let's see how to get deltaS from an existing Forward matrix.
 *    
 * Because t_EJ = t_EC and t_JJ = t_CC in the special case of
 * interest, we have C(j) = J(j) for all j in a standard Forward
 * matrix, so we can use either; let's arbitrarily use J().  The
 * summed score in J(ib) includes:
 *    1. Paths that went thru N(ia-1) and contain >=1 domain  = N(ia-1) + deltaS
 *    2. Paths that went thru J(ia-1) and contain >=1 domain  = J(ia-1) + deltaS
 *    3. Paths that went thru J(ia-1) and contain 0 domains.  = J(ia-1) + p0.
 *    
 * Thus:              
 *     J(ib) = log[ exp(J(ia-1) + p0) +
 *                  exp(N(ia-1) + deltaS) +
 *                  exp(J(ia-1) + deltaS)]
 *
 * Given a sparse Forward matrix, ia and ib for a sparse segment, and
 * t_XX, we know N(ia-1), J(ia-1), J(ib). So we can solve for deltaS:
 * 
 *     deltaS = J(ib) + log[ 1 - exp(J(ia-1) + p0 - J(ib)) ] - log[ exp(N(ia-1)) + exp(J(ia-1)) ]
 *                               
 *            = J(ib) + log[ 1 - exp(-deltaC)] -  logsum(N(ia-1), J(ia-1))
 * 
 * We've rewritten the middle term in terms of the difference "deltaC":
 * 
 *    deltaC = J(ib) - J(ia-1) - p0
 *    
 * That is: if the only possible path were the empty path p0, we'd
 * reach J(ib) with score J(ia-1) + p0, and deltaC would be 0. Other
 * paths contribute nonnegative log (or log-odds) probability to the
 * difference deltaC.  Thus, deltaC must be >= 0.
 * 
 * xref: SRE:J10/43
 * 
 * 
 * 
 *****************************************************************
 ***************************************************************** 
 *    
 * [3] NUMERICAL ISSUES WITH CALCULATING deltaS FROM deltaC
 * 
 * In the small limit deltaC -> 0, log[1 - exp(-deltaC)] -> -inf, and deltaS -> -inf.
 * We use lim{x->0} 1-e^-x = x:
 *    if deltaC < eslSMALLX1   :   log[1 - exp(-deltaC)] = log(deltaC)      -> -inf
 *    
 * In the large limit deltaC -> inf, exp(-deltaC) -> 0, and we can use
 * lim{x->0} log(1-x) = -x:
 *    if exp(-deltaC) < eslSMALLX1 :  log[1 - exp(deltaC)] = -exp(-deltaC)  -> 0 from below
 *    
 * To see what's going on, suppose we know that we started in J(ia-1)
 * (so we don't have to think about the last term and the summed
 * uncertainty over N|J at ia-1):
 * 
 * If deltaC = 0, that means that only the empty path gets us from
 * ia-1..ib, and deltaS = -inf.
 * 
 * If deltaC = 1000, that means that the contribution of the empty
 * path is negligible; the middle term asymptotes to zero (from below; i.e. term2 is
 * negative), and deltaS converges toward just the difference J(ib) - J(ia-1).
 *                                              
 * xref: SRE:J10/44
 * 
 *****************************************************************
 ***************************************************************** 
 *
 * [4] SEGMENT TRACES
 * 
 * A "segment trace" is a hack on a P7_TRACE structure to allow us to
 * deal with a subpath over a segment ia..ib. 
 * 
 * tr->st[0]   is p7T_N | p7T_J    (instead of p7T_S)
 * tr->st[N-1] is p7T_J            (instead of p7T_T)
 * 
 * tr->i[z] are ia..ib (or 0)      (instead of 1..L)
 *   Note that <i> coords remain relative to the original 1..L sequence.
 *   Same is true for the anchor set we construct from a segment trace.
 *   Because a segment trace can only emit ia..ib, all of which
 *   are (of course) in a segment, we also know sm->n[i] > 0 for
 *   all residues i=ia..ib.
 *   
 * tr->k[z], if nonzero, must be in the sparse list sm->k[tr->i[z]].
 * 
 * A segment trace may have zero domains (unlike a standard P7_TRACE), 
 * in which case:
 *    tr->N = Ls+1 = ib-ia+2
 *    tr->st[z] = p7T_J for z = 0..Ls
 *    tr->i[z]  = {0,ia..ib}
 *    tr->k[z]  = 0 for all z
 *
 *
 *****************************************************************
 *****************************************************************
 * 
 * [5] NONINDEPENDENCE OF SEGMENTAL MPAS OPTIMA
 *     see "A glitch in segmental MPAS algorithm", SRE:J13/146.
 *     
 * The segment forward score is defined as logsum(deltaS, p0): the log
 * sum of all paths from X(ia-1) to Y(ib) that have D>0 (that have
 * zero or more passes through the homology model).
 * 
 * For almost all choices of X and Y, the ensemble deltaS+p0 is the
 * correct ensemble, and the 
 * 
 * The sole exception is when there are no domains upstream or
 * downstream.  Now X must be N and Y must be C, and we cannot use the
 * empty path p0; the valid path ensemble is only deltaS. 
 * 
 * We can still run MPAS by sampling traces (thus anchor sets) on the
 * deltaS+p0 ensemble. If we obtain a D>0 optimal solution on the
 * incorrect deltaS+p0 ensemble, that's fine, because that solution
 * must also be the optimum in the correct deltaS ensemble. A problem
 * only arises if we get a D=0 optimal MPAS solution on every segment.
 * Then we have to fail over to a complete-sequence MPAS, because
 * there must be at least one anchor somewhere. We don't expect this
 * case to arise: for a target sequence to make it as far as the MPAS
 * calculation, probability mass has already convinced us that there
 * are one or more probable domains in the sequence.
 * 
 * xref SRE:J13/146.
 * 
 *****************************************************************
 ***************************************************************** 
 *
 * [6] SAMPLING D=0 EMPTY TRACES FROM THE sF (deltaS + p0) ENSEMBLE
 * 
 * Stochastic segment trace starts in J(ib) and traces back to either
 * N(ia-1) or J(ia-1). That samples N->J and J->J paths (equiv to
 * other choices of X,Y endpoints such as N->C and J->C) but does not
 * fully account for the probability mass of the empty path p0 because
 * we don't sample N->N (resp. C->C). 
 * 
 * So the empty path must be handled specially. We sample traces by a
 * two-step process, handling the deltaS and p0 components of sF
 * separately.
 * 
 * First, we ask if we use the empty path. The probability of using
 * the empty path is exp(p0 - sF).
 *                            
 * Else, we must sample from the deltaS ensemble (D>0). A traceback 
 * from J(ib) samples from an ensemble of deltaS + the part of the p0 
 * "ensemble" that uses J->J. So we can use rejection sampling:
 * sample, reject sample if D=0, iterate until we get D>0 trace.
 * 
 * Danger: what happens if p0 is large, but by chance we don't sample
 * the empty path in step 1, and we make it to the rejection sampling
 * stage in step 2? We could be rejection sampling for an arbitrarily
 * long time. Fortunately, we don't have to worry about this danger.
 * If p0 > deltaS, then we already know that the optimal MPAS solution
 * for this segment is D=0 (no anchors), and there's no need to run
 * the MPAS sampling optimization in the first place! We only need to
 * sample segment traces in the case where p0 <= deltaS, which means a
 * minimum 50% probability of getting D>0 sample.
 */

/*****************************************************************
 * 5. Statistics driver
 *****************************************************************/
#ifdef p7SPARSE_ANCHORS_STATS

#include "p7_config.h"

#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type           default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",          eslARG_NONE,   FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",                   0 },
  { "-s",          eslARG_INT,      "0",  NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                          0 },
  { "-G",          eslARG_NONE,   FALSE,  NULL, NULL,   NULL,  NULL, NULL, "run 'global' version instead of segmental D&C",          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "statistics collection on sparse MPAS algorithm";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create( esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  P7_HMMFILE     *hfp     = NULL;
  ESL_ALPHABET   *abc     = NULL;
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_SQ         *sq      = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  ESL_SQFILE     *sqfp    = NULL;
  P7_BG          *bg      = NULL;
  P7_HMM         *hmm     = NULL;
  P7_PROFILE     *gm      = NULL;           /* profile in H4's standard dual-mode local/glocal */
  P7_OPROFILE    *om      = NULL;
  P7_FILTERMX    *fx      = p7_filtermx_Create(100);
  P7_CHECKPTMX   *cx      = p7_checkptmx_Create(100, 100, ESL_MBYTES(p7_SPARSIFY_RAMLIMIT));
  P7_SPARSEMASK  *sm      = p7_sparsemask_Create(100, 100);
  P7_ANCHORS     *anch    = p7_anchors_Create();
  P7_ANCHORS     *vanch   = p7_anchors_Create();
  P7_ANCHORHASH  *ah      = p7_anchorhash_Create();
  P7_SPARSEMX    *sxf     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *sxd     = p7_sparsemx_Create(NULL);
  P7_SPARSEMX    *asf     = p7_sparsemx_Create(NULL);
  P7_TRACE       *tr      = p7_trace_Create();
  float          *wrk     = NULL;
  float           nullsc, vsc, fsc, asc;
  clock_t         start_c, end_c;
  clock_t         total_c = 0;
  clock_t         init_c  = clock();
  int             status;

  /* Read in one HMM. Set alphabet to whatever the HMM's alphabet is. */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Configure vector, dual-mode, and local-only profiles from HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create (hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_oprofile_Convert(gm, om);
  
  /* Open sequence file */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* For each target sequence... */
  while (( status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_bg_SetLength           (bg,  sq->n);
      p7_profile_SetLength      (gm,  sq->n);
      p7_oprofile_ReconfigLength(om,  sq->n);

      if (( status = p7_pipeline_AccelerationFilter(sq->dsq, sq->n, om, bg, fx, cx, sm)) == eslOK)
	{
	  printf("%-15s  %-30s  ", gm->name, sq->name);

	  p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

	  p7_SparseViterbi (sq->dsq, sq->n, gm, sm, sxf, tr, &vsc);
	  p7_SparseForward (sq->dsq, sq->n, gm, sm, sxf, &fsc);
	  p7_SparseBackward(sq->dsq, sq->n, gm, sm, sxd, NULL);   
	  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxd, sxd);   


	  start_c = clock();

	  p7_sparse_anchors_SetFromTrace(sxd, tr, vanch);
	  p7_trace_Reuse(tr);

	  if (esl_opt_GetBoolean(go, "-G"))
	    p7_sparse_AnchorsGlobal(rng, sq->dsq, sq->n, gm, 
				    vsc, fsc, sxf, sxd, vanch,
				    tr, &wrk, ah, 
				    asf, anch, &asc,
				    NULL);
	  else
	    p7_sparse_Anchors(rng, sq->dsq, sq->n, gm, 
			      vsc, fsc, sxf, sxd, vanch, 
			      tr, &wrk, ah,
			      asf, anch, &asc, NULL);

	  end_c   = clock();
	  total_c += end_c - start_c;
	     
	  p7_sparsemx_Reuse(sxf);
	  p7_sparsemx_Reuse(sxd);
	  p7_sparsemx_Reuse(asf);
	  p7_trace_Reuse(tr);
	  p7_anchorhash_Reuse(ah);
	  p7_anchors_Reuse(anch);
	  p7_anchors_Reuse(vanch);
	}

      p7_sparsemask_Reuse(sm);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));
  else if (status != eslEOF)     p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);

  printf("# Total time in p7_sparse_Anchors = %.4f sec\n", (double) total_c / (double) CLOCKS_PER_SEC);
  printf("# Total time overall              = %.4f sec\n", (double) (clock() - init_c) / (double) CLOCKS_PER_SEC);

  if (wrk) free(wrk);
  p7_anchors_Destroy(vanch);
  p7_anchors_Destroy(anch);
  p7_anchorhash_Destroy(ah);
  p7_trace_Destroy(tr);
  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(sxf);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(cx);
  p7_filtermx_Destroy(fx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7SPARSE_ANCHORS_STATS*/

/*----------------- end, statistics driver ----------------------*/




/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef p7SPARSE_ANCHORS_EXAMPLE
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
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "include all cells in sparse mx",                   0 },
  { "-s",        eslARG_INT,      "0", NULL, NULL,   NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-G",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "run 'global' version instead of segmental D&C",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of running sparse MPAS (most probable anchor set) algorithm";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng     = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  P7_CHECKPTMX   *cx      = NULL;
  P7_SPARSEMASK  *sm      = NULL;
  P7_ANCHORS     *anch    = p7_anchors_Create();
  P7_ANCHORS     *vanch   = p7_anchors_Create();
  P7_ANCHORHASH  *ah      = p7_anchorhash_Create();
  P7_SPARSEMX    *sxf     = NULL;
  P7_SPARSEMX    *sxd     = NULL;
  P7_SPARSEMX    *asf     = NULL;
  P7_TRACE       *tr      = NULL;
  float          *wrk     = NULL;
  P7_MPAS_PARAMS  prm;
  float           fsc, vsc, asc, nullsc;
  int             status;

  /* Customize parameters for MPAS */
  prm.max_iterations = 1000;
  prm.loss_threshold = 0.001;
  prm.nmax_sampling  = FALSE;
  prm.be_verbose     = TRUE;

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
 
  /* Read one sequence */
  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslOK)      p7_Fail("Unexpected error %d reading sequence file %s", status, sqfp->filename);
  esl_sqfile_Close(sqfp);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  om = p7_oprofile_Create(hmm->M, abc);
  p7_profile_Config(gm, hmm, bg);
  p7_oprofile_Convert(gm, om);

  /* Set the profile and null model's target length models */
  p7_bg_SetLength     (bg, sq->n);
  p7_profile_SetLength(gm, sq->n);
  p7_oprofile_ReconfigLength(om, sq->n);

  /* Make a sparse mask, <sm> */
  cx = p7_checkptmx_Create(hmm->M, sq->n, ESL_MBYTES(32));
  sm = p7_sparsemask_Create(gm->M, sq->n);
  if (esl_opt_GetBoolean(go, "-a")) 
    p7_sparsemask_AddAll(sm);
  else {
    p7_ForwardFilter (sq->dsq, sq->n, om, cx, /*fsc=*/NULL);
    p7_BackwardFilter(sq->dsq, sq->n, om, cx, sm, p7_SPARSIFY_THRESH);
  }

  /* Allocate DP matrices, tracebacks */
  sxf  = p7_sparsemx_Create(sm);
  sxd  = p7_sparsemx_Create(sm);
  asf  = p7_sparsemx_Create(sm);
  tr   = p7_trace_Create();

  /* First pass analysis */
  p7_SparseViterbi (sq->dsq, sq->n, gm, sm,  sxf,  tr, &vsc);
  p7_SparseForward (sq->dsq, sq->n, gm, sm,  sxf,      &fsc);
  p7_SparseBackward(sq->dsq, sq->n, gm, sm,  sxd,      NULL);
  p7_SparseDecoding(sq->dsq, sq->n, gm, sxf, sxd,      sxd);
  p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc);

  //p7_trace_DumpAnnotated(stdout, vtr, gm, sq->dsq);

  /* Do it */

  p7_sparse_anchors_SetFromTrace(sxd, tr, vanch);
  p7_trace_Reuse(tr);

  if (esl_opt_GetBoolean(go, "-G")) 
    p7_sparse_AnchorsGlobal(rng, sq->dsq, sq->n, gm, 
			    vsc, fsc, sxf, sxd, vanch, 
			    tr, &wrk, ah, 
			    asf, anch, &asc, 
			    &prm);
  else
    p7_sparse_Anchors(rng, sq->dsq, sq->n, gm, 
		      vsc, fsc, sxf, sxd, vanch,
		      tr, &wrk, ah,
		      asf, anch, &asc, &prm);

  printf("Viterbi score = %.2f nats\n", vsc);
  printf("Forward score = %.2f nats\n", fsc);
  printf("Optimal anchor set:  ");
  p7_anchors_DumpOneLine(stdout, anch);
  printf("ASC Forward   = %.2f nats\n", asc);
  printf("ASC prob      = %.4f\n", exp(asc - fsc));

  if (wrk) free(wrk);
  p7_trace_Destroy(tr);
  p7_sparsemx_Destroy(asf);
  p7_sparsemx_Destroy(sxd);
  p7_sparsemx_Destroy(sxf);
  p7_anchorhash_Destroy(ah);
  p7_anchors_Destroy(anch);
  p7_anchors_Destroy(vanch);
  p7_sparsemask_Destroy(sm);
  p7_checkptmx_Destroy(cx);
  esl_sq_Destroy(sq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7SPARSE_ANCHORS_EXAMPLE*/
