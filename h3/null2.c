/* "null2" model: biased composition correction
 * 
 * See p7_domaindef.c -- null2 correction of per-seq and per-domain
 * scores is essentially embedded within p7_domaindef's logic, but we
 * split it out to a separate file because it's so important.
 * 
 * SRE, Thu Feb 28 09:51:27 2008 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"
#include "easel.h"
#include "hmmer.h"



/* Function:  p7_Null2Corrections()
 * Synopsis:  Calculate biased composition corrections
 * Incept:    SRE, Thu Feb 28 09:52:28 2008 [Janelia]
 *
 * Purpose:   Calculate biased composition corrections for per-domain
 *            and per-sequence scores, for scoring an envelope $i..j$
 *            with model <gm>. The envelope $i..j$ is defined by
 *            passing a <dsq> pointer to residue $i-1$, and an
 *            envelope length <Ld>, where <Ld> $ = j-i+1$.
 *            
 *            (Why not pass the original full-length <dsq> and $i,j$
 *            envelope coords, instead of the offset <dsq> pointer?
 *            Because the caller has calculated small Forward and
 *            Backward matrices for the offset <dsq>.)
 *            
 *            Caller provides calculated Forward and Backward matrices
 *            in <fwd> and <bck>, for the chunk of digital sequence in
 *            the envelope: e.g. the same offset <dsq+i-1> pointer,
 *            over length <Ld>. The null2 correction uses these to
 *            calculate a log odds emission score vector as the
 *            posterior weighted average over all emission vectors
 *            used in alignments (inclusive of N,C,J) in the envelope
 *            <i>..<j>.
 *            
 *            The caller also provides <noverlap>, the number of
 *            residues that overlap with a previous envelope. This is
 *            used to avoid overcounting null2 corrections against the
 *            per-sequence score. If no preceding envelope overlaps
 *            <i>..<j>, <noverlap> is 0.
 *            
 *            The calculation returns two corrections. The per-domain
 *            correction <*opt_domcorrection> is the score (in nats)
 *            that the caller should subtract from the domain envelope
 *            score for <i>..<j>. The per-sequence correction
 *            <*opt_seqcorrection> is the score (in nats) that the
 *            caller should subtract from the overall per-sequence
 *            score.  
 *            
 *            If the envelope <i..j> is independent (nonoverlapping)
 *            with any previous domain, these two corrections are the
 *            same number; when the envelope overlaps with a previous
 *            envelope, the per-seq correction only accounts for the
 *            last <Ld - noverlap> residues of the envelope, even
 *            though the null model is calculated over the entire
 *            envelope. This isn't well principled; it's just a way of
 *            making sure each residue only gets corrected once in the
 *            per-sequence score.
 *
 * Args:      gm                - profile, in any mode, target length model set to <L>
 *            dsq               - offset digital seq being scored; <dsq+i-1>; <1..Ld>
 *            Ld                - length of domain envelope
 *            noverlap          - number of residues that overlap with a previous envelope
 *            fwd               - Forward matrix,  filled with <gm> against domain envelope in <dsq>
 *            bck               - Backward matrix, filled with <gm> against domain envelope in <dsq>
 *            opt_domcorrection - optRETURN: per-domain score correction (in nats) to subtract
 *            opt_seqcorrection - optRETURN: per-seq score correction (in nats) to subtract
 *
 * Returns:   <eslOK> on success, and <*opt_domcorrection>, <*opt_seqcorrection> contain
 *            the corrections.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
p7_Null2Corrections(const P7_PROFILE *gm, const ESL_DSQ *dsq, int Ld, int noverlap,
		    const P7_GMX *fwd, const P7_GMX *bck,
		    float *opt_domcorrection, float *opt_seqcorrection)
{
  int    M          = gm->M;
  float  overall_sc = bck->xmx[p7G_N];
  float  null2[p7_MAXABET];
  float  domsc, seqsc;
  float  omega;
  int    x;			/* over symbols 0..K-1                       */
  int    i;			/* over offset envelope dsq positions 1..Ld  */
  int    k;			/* over model M states 1..M, I states 1..M-1 */

  /* Calculate null2's log odds emission probabilities,
   * by taking posterior weighted sum over all emission vectors used in paths explaining the domain.
   * This is dog-slow; a point for future optimization
   */
  for (x = 0; x < gm->abc->K; x++)
    {
      null2[x] = -eslINFINITY;
      for (i = 1; i <= Ld; i++) 
	{
	  for (k = 1; k < M; k++)
	    {
	      null2[x] = p7_FLogsum(null2[x], 
				    p7P_MSC(gm, k, x) + 
				    fwd->dp[i][k*p7G_NSCELLS + p7G_M] + bck->dp[i][k*p7G_NSCELLS + p7G_M] - overall_sc);
	      null2[x] = p7_FLogsum(null2[x], 
				    p7P_ISC(gm, k, x) + 
				    fwd->dp[i][k*p7G_NSCELLS + p7G_I] + bck->dp[i][k*p7G_NSCELLS + p7G_I] - overall_sc);
	    }	  
	  null2[x] = p7_FLogsum(null2[x], 
				p7P_MSC(gm, M, x) + 
				fwd->dp[i][M*p7G_NSCELLS + p7G_M] + bck->dp[i][M*p7G_NSCELLS + p7G_M] - overall_sc);

	  /* N,C,J contribute log odds scores of 0 (by def'n), weighted by posterior decoding */
	  null2[x] = p7_FLogsum(null2[x],  fwd->xmx[p7G_NXCELLS*(i-1) + p7G_N] + bck->xmx[p7G_NXCELLS*i + p7G_N] + gm->xsc[p7P_N][p7P_LOOP] - overall_sc);
	  null2[x] = p7_FLogsum(null2[x],  fwd->xmx[p7G_NXCELLS*(i-1) + p7G_J] + bck->xmx[p7G_NXCELLS*i + p7G_J] + gm->xsc[p7P_J][p7P_LOOP] - overall_sc);
	  null2[x] = p7_FLogsum(null2[x],  fwd->xmx[p7G_NXCELLS*(i-1) + p7G_C] + bck->xmx[p7G_NXCELLS*i + p7G_C] + gm->xsc[p7P_C][p7P_LOOP] - overall_sc);
	}
      null2[x] -= log(Ld);
    }

  /* Tote up the null2 path scores for the envelope */
  noverlap = (noverlap > Ld ? Ld : noverlap); /* just checking. */
  domsc = seqsc = 0.0;
  for (i = 1; i <= noverlap; i++) domsc += null2[dsq[i]];
  for ( ;     i <= Ld;       i++) seqsc += null2[dsq[i]];
  domsc += seqsc;

  /* Sum over the null1, null2 paths for this domain */
  omega = 1.0f / 256.0f;
  domsc = -1. * p7_FLogsum( log (1. - omega), log(omega) + domsc);
  seqsc = -1. * p7_FLogsum( log (1. - omega), log(omega) + seqsc);

#if 0
  /* Add in the corrections to s_d, the per-domain score */
  nsc += log (3.0 / ((float) L + 3.));          /* or N->B score, if <gm> were in L=L mode */
  nsc -= eslCONST_LOG2;		                /* or E->J/E->C score, if <gm> were in multihit mode */
  nsc -= (float) Ld * log( (float) L / (float) (L+3)); /* or Ld * NN/CC/JJ score, if <gm> were in L=L mode */
#endif

  if (opt_domcorrection != NULL) *opt_domcorrection = domsc;
  if (opt_seqcorrection != NULL) *opt_seqcorrection = seqsc;
  return eslOK;
}

