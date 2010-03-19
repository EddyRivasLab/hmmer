/* Miyazawa's (1994) partition function
 *
 * Changes to the original algorithm:
 *
 *              - local instead of global
 *              - no I->D and D->I transitions allowed
 *              - multiple gaps that lead to ambiguous parsing not allowed
 *
 * Simple implementation:
 *
 *              - no sse code involved
 *              - not using HMMER3 MMX, DMX and IMX matrices
 *              - keeping whole dp matrices, which is wasteful
 *
 *  Returns:   zscore
 */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_scorematrix.h"

#include "hmmer.h"

/*****************************************************************
 * 1. Miyazawa implementation.
 *****************************************************************/

/* Function:  p7_GMiyazawa()
 * Synopsis:  The Miyazawa's partition function algorithm.
 * Incept:    SC, Wed Mar 17 11:37:34 EDT 2010
 *
 * Purpose:   The standard Smith-Waterman dynamic programming algorithm.
 *
 *            Given a digital target sequence <dsq> of length <L>, a
 *            query sequence and score system <sm>, and DP matrix <gx>
 *            allocated for at least <L> by <sm->n> cells; calculate the
 *            partition function (sum over all paths) by Miyazawa; return
 *            the Miyazawa score in <ret_sc>, and the Miyazawa matrix is
 *            in <gx>.
 *
 *            The Miyazawa lod score is returned in nats. The
 *            caller needs convert to bits. REALLY???
 *
 * Args:      dsq    - target sequence in digitized form, 1..L
 *            L      - length of dsq
 *            sm     - query sequence and score system
 *            gx     - DP matrix with room for an MxL alignment
 *            opt_sc - optRETURN: Viterbi lod score in nats
 *
 * Return:   <eslOK> on success.
 */

int
p7_GMiyazawa(const ESL_DSQ *dsq, int L, P7_SCORESYS *sm, P7_GMX *gx, float *ret_sc)
{
	float      **dp   = gx->dp;
	int          M    = sm->n;
	int          i;                    /* index over rows (target)   */
	int 				 k;						         /* index over columns (query) */
	double  slambda   = sm->slambda;
	double     **MSC  = sm->Q->mx;
	double    lopen   = sm->lopen;
	double  lextend   = sm->lextend;
	float        sc   = 0;

	/* Initialization */
	for (k = 0; k <= M; k++)
	{
		MMX(0,k) = IMX(0,k) = DMX(0,k) = -eslINFINITY;   /* first row */
	}

	/* Recursion */
	for (i = 1; i <= L; i++)              /* loop over target sequence */
	{
		/* Initialization */
		MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY; /* first column */

		for (k = 1; k <= M; k++)            /* loop over query sequence */
		{
			/* match */
			sc = p7_FLogsum(p7_FLogsum(0, MMX(i-1,k-1)), p7_FLogsum(IMX(i-1,k-1), DMX(i-1,k-1)));
			MMX(i,k) = sc + MSC[dsq[i]][sm->dsq[k]];

			/* insert (in the target) */
			IMX(i,k) = p7_FLogsum(MMX(i-1,k) - lopen, IMX(i-1,k) - lextend);

			/* delete (in the target) */
			DMX(i,k) = p7_FLogsum(MMX(i,k-1) - lopen, DMX(i,k-1) - lextend);

		} /* end loop over query sequence */

	} /* end loop over target sequence */

	*ret_sc = MMX(L,M); // IS THIS CORRECT?

	return eslOK;

}
