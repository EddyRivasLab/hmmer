/* Most probable anchor set (MPAS) algorithm 
 * Sparse implementation: production code.
 * 
 * Contents:
 *     1. MPAS algorithm
 *     x. Notes.
 *     x. Copyright and license information.
 */

#include "search/p7_mpas.h"     // P7_MPAS_PARAMS and P7_MPAS_STATS

static int segment_forward_score(const P7_PROFILE *gm, const P7_SPARSEMX *sx, int g, float *xc, float *ret_sF);
static int sparse_anchors_SetFromSegmentTrace(const P7_SPARSEMX *sxd, const P7_TRACE *tr, int g, float *dpc, P7_ANCHORS *anch, float *ret_dpn);

/*****************************************************************
 * 1. MPAS algorithm
 *****************************************************************/


int
p7_sparse_Anchors(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,
		  const P7_SPARSEMX *sxf, const P7_SPARSEMX *sxd,
		  float vsc, float fsc,
		  P7_TRACE *tr, float **byp_wrk, P7_ANCHORHASH *ah,
		  P7_SPARSEMX *afu, P7_SPARSEMX *afd, P7_ANCHORS *anch, P7_ANCHORS *tmpa, float *ret_asc,
		  P7_MPAS_PARAMS *prm, P7_MPAS_STATS *stats)
{
  
  
  if (exp(vsc - fsc) > 0.5)
    {
      sparse_anchors_SetFromTrace(sxd, tr, anch);
      return eslOK;
    }
  
  

  for (g = 1; g < sm->S; g++)
    {
      
      /* First guess: D=0, no anchor in segment; 
       */

      


    }
  




}









/* sparse_anchors_SetFromSegmentTrace()
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
 *            <anch->D> anchors, we also return a pointer <ret_dpn>:
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
sparse_anchors_SetFromSegmentTrace(const P7_SPARSEMX *sxd, const P7_TRACE *tr, int g, float *dpc, P7_ANCHORS *anch, float *ret_dpn)
{
  P7_SPARSEMASK *sm = sxd->sm;
  const float   *dpc2;
  float          ppv;
  float          best_ppv = -1.;
  int            z;
  int            y;
  int            status;

  p7_anchors_Reuse(anch);
  anch->D = 0;               // incremented when/if we see an E state

  ESL_DASSERT1(( sxd->type == p7S_DECODING ));

  if (dpc == NULL)
    {
      dpc = sxd->dp;
      for (i = 1; i < sm->seg[g].ia; i++) 
	dpc += p7_NSCELLS * sm->n[i];
    }

  for (z = 0; z < tr->N; z++)
    {
      if (p7_trace_IsM(tr->st[z]))
	{
	  i    = tr->i[z];
	  y    = 0; 
	  kc   = sm->k[i];
	  dpc2 = dpc;
	  while (y < sm->n[i] && kc[y]  < k) { y++; dpc2 += p7S_NSCELLS; }
	  if    (y < sm->n[i] && kc[y] == k) {
	    ppv = dpc2[p7S_ML] + dpc2[p7S_MG];
	    if (ppv > best_ppv)
	      {
		anch->a[anch->D+1].i0 = tr->i[z];  // D+1 because anchors are numbered 1..D, but D here is a counter that starts at 0
		anch->a[anch->D+1].k0 = tr->k[z];
	      }
	  }
	}
      else if (tr->st[z] == p7T_E)
	{
	  anch->D++;
	  best_ppv = -1.;
	  if ((status = p7_anchors_Grow(anch)) != eslOK) goto ERROR;
	}

      if (tr->i[z])
	dpc += p7_NSCELLS * sm->n[i];
    }

  p7_anchor_SetSentinels(anch->a, anch->D, tr->L, tr->M);
  *ret_dpn = dpc;
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
 *            matrix <sx>. 
 *
 *            If caller has a valid ptr <xc> to the special row for
 *            <seg[g].ia-1> in <sx>, it should provide it. Else, pass
 *            <xc=NULL>, and the needed special values will be looked
 *            up by a somewhat more laborious means. (Sparse matrices
 *            are designed for sequential, not random access.)
 *            
 *            Calculation requires that t_NN = t_JJ = t_CC in the
 *            profile.
 *            
 *            Segment Forward score for segment <g> is returned thru
 *            <ret_sF>.
 *
 * Returns:   <eslOK> on success.
 */
static int
segment_forward_score(const P7_PROFILE *gm, const P7_SPARSEMX *sx, int g, float *xc, float *ret_sF)
{
  float Jia1;			// Forward score in <sx> at J(ia-1)       
  float Nia1;			// Forward score in <sx> at N(ia-1)       
  float Jib;			// Forward score in <sx> at J(ib)         
  int   Ls;			// Length of segment g in residues        
  float p0;			// Score term for D=0 path from ia-1..ib  
  float deltaC;			// Intermediate term in calculating sF(g); see footnote [2] 
  float deltaS;			//   ... another intermediate term; see footnote [2]        
  float term2;			//   ... another intermediate term; see footnote [2]        
  
  /* Contract checks.
   * Using segmented MPAS algorithm optimization requires t_NN = t_JJ = t_CC;
   * see footnote [1].
   */
  ESL_DASSERT1(( gm->xsc[p7P_N][p7P_N] == gm->xsc[p7P_J][p7P_J] ));
  ESL_DASSERT1(( gm->xsc[p7P_N][p7P_N] == gm->xsc[p7P_C][p7P_C] ));
  ESL_DASSERT1(( sx->type == p7S_FORWARD ));

  /* If caller doesn't have <xc> positioned already on ia-1, 
   * it passes xc=NULL, and we set it the slow way.
   */
  if (xc == NULL)
    {
      int g2;
      xc = sx->xmx;
      for (g2 = 1; g2 < g; g++)
	xc += p7S_NXCELLS * (sx->sm->seg[g].ib- sx->sm->seg[g].ia+2);  // ia-1..ib = ib-ia+2 to skip segment
    }
  
  Ls     = sx->sm->seg[g].ib- sx->sm->seg[g].ia+1;
  Jia1   = xc[p7S_J];
  Nia1   = xc[p7S_N];
  xc    += p7S_NXCELLS * Ls;  // now <xc> is on <ib>
  Jib    = xc[p7S_J];
  p0     = Ls * gm->xsc[p7P_J][p7P_J];
  deltaC = Jib - Jia1 - p0;
	 
  /* See footnote [3], numerical issues with calculating deltaS from deltaC: */
  if      (deltaC < eslSMALLX1)            term2 = log(deltaC);
  else if (deltaC > -1. * log(eslSMALLX1)) term2 = -1. * exp(-1. * deltaC);
  else                                     term2 = log(1. - exp(-1. * deltaC));

  deltaS  = Jib + term2 - p7_FLogsum(Nia1, Jia1);
  *ret_sF = p7_FLogsum(deltaS, p0);
  return eslOK;
}


/*****************************************************************
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
 * @LICENSE@
 * 
 * SVN $Id$
 * SVN $URL$
 ***************************************************************** 
