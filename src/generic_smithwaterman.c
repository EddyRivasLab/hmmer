/* Smith-Waterman algorithm (non-sse version) */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_scorematrix.h"

#include "hmmer.h"

/*****************************************************************
 * 1. Smith-Waterman implementation.
 *****************************************************************/

/* Function:  p7_GSmithWaterman()
 * Synopsis:  The Smith-Waterman algorithm.
 * Incept:    SC, Fri Mar 12 10:30:19 EST 2010
 *
 * Purpose:   The standard Smith-Waterman dynamic programming algorithm.
 *
 *            Given a digital target sequence <dsq> of length <L>, a
 *            query sequence and score system <sm>, and DP matrix <gx>
 *            allocated for at least <L> by <sm->n> cells; calculate the
 *            maximum scoring path by Smith-Waterman; return the Smith-
 *            Waterman score in <ret_sc>, and the Smith-Waterman matrix
 *            is in <gx>.
 *
 *            The caller may then retrieve the Smith-Waterman path by
 *            calling <NOT IMPLEMENTED>.
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
p7_GSmithWaterman(const ESL_DSQ *dsq, int L, P7_SCORESYS *sm, P7_GMX *gx, float *ret_sc)
{
	float      **dp   = gx->dp;
	int          M    = sm->n;
	int          i;                    /* index over rows (target)   */
	int 				 k;						         /* index over columns (query) */
	double    sopen   = sm->sopen;
	double  sextend   = sm->sextend;
  int      **MSC    = sm->S->s;
	float        sc   = 0;
	float       osc   = 0;

	/* Initialization */
	for (k = 0; k <= M; k++)
	{
		MMX(0,k) = 0; 											  /* first row */
		IMX(0,k) = DMX(0,k) = -eslINFINITY;   /* first row */
	}

	/* Recursion */
	for (i = 1; i <= L; i++)              /* loop over target sequence */
	{
		/* Initialization */
		MMX(i,0) = 0;                       /* first column */
		IMX(i,0) = DMX(i,0) = -eslINFINITY; /* first column */

		for (k = 1; k <= M; k++)            /* loop over query sequence */
		{
			/* match */
			sc = ESL_MAX(MMX(i-1,k-1), IMX(i-1,k-1));
			sc = ESL_MAX(sc, DMX(i-1,k-1));

			sc = sc + MSC[dsq[i]][sm->dsq[k]];      /* s[target][query] */

			if (sc > 0) MMX(i,k) = sc; else MMX(i,k) = 0;

			/* optimal score (in MMX) */
			if (MMX(i,k) >= osc) osc = MMX(i,k);

			/* insert (in the target) */
			IMX(i,k) = ESL_MAX((MMX(i-1,k) - sopen), (IMX(i-1,k) - sextend));

			/* delete (in the target) */
			DMX(i,k) = ESL_MAX((MMX(i,k-1) - sopen), (DMX(i,k-1) - sextend));

		} /* end loop over query sequence */

	} /* end loop over target sequence */

	*ret_sc = osc;

	return eslOK;
}
