/* 
 * 
 * In terms of control flow, see p7_domaindef.c -- null2.c
 * (calculation of the null2 correction, and per-domain scores)
 * is essentially embedded within p7_domaindef's logic.
 * 
 * SRE, Thu Feb 28 09:51:27 2008 [Janelia]
 * SVN $Id$
 */
#include "p7_config.h"
#include "easel.h"
#include "hmmer.h"



/* Function:  
 * Synopsis:  
 * Incept:    SRE, Thu Feb 28 09:52:28 2008 [Janelia]
 *
 * Purpose:   
 *
 * Args:      dsq - digital seq pointing to the start-1 position of
 *                  domain to be score-corrected. This is typically
 *                  dsq + (i-1) for an original (longer) target seq
 *                  with a domain i..j we're rescoring.
 *             Ld - length of the domain. This is typically j-i+1 for
 *                  one domain i..j in a larger target sequence.
 *             L  - the original total length of the sequence.
 *            fwd - Forward matrix, filled with <gm> against domain in <dsq>
 *            bck - Backward matrix, filled with <gm> against domain in <dsq>
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
p7_Null2(const P7_PROFILE *gm, const ESL_DSQ *dsq, int Ld, int L, const P7_GMX *fwd, const P7_GMX *bck)
{
  int    M          = gm->M;
  float  overall_sc = bck->xmx[p7G_N];
  int x;			/* over symbols 0..K-1     */
  int i;			/* over dsq positions 1..L */
  int k;			/* over model M states 1..M, I states 1..M-1 */

  float        null2[p7_MAXABET];
  float        nsc;
  float        omega;

  /* Calculate null2's log odds emission probabilities,
   * by taking posterior weighted sum over all emission vectors used in paths explaining the domain.
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

  /* Tote up the null2 path score for the domain */
  nsc = 0.0;
  for (i = 1; i <= Ld; i++) nsc += null2[dsq[i]];

  /* Sum over the null1, null2 paths for this domain */
  omega = 1.0f / 256.0f;
  nsc   = -1. * p7_FLogsum( log (1. - omega), log(omega) + nsc);

  printf("nsc = %.2f    (the null2 correction; i.e. after logsum)\n\n", nsc);

  /* Add in the corrections to s_d, the per-domain score */
  nsc += log (3.0 / ((float) L + 3.));          /* or N->B score, if <gm> were in L=L mode */
  nsc -= eslCONST_LOG2;		                /* or E->J/E->C score, if <gm> were in multihit mode */
  nsc -= (float) Ld * log( (float) L / (float) (L+3)); /* or Ld * NN/CC/JJ score, if <gm> were in L=L mode */

  printf("total correction to raw domain score: %.2f\n\n", nsc);
}

