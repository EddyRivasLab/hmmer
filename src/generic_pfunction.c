/*
 * Miyazawa's (1994) partition function
 *
 * Changes to the original algorithm:
 *
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"

int
pfunction(ESL_DSQ *qsq,  ESL_DSQ *dbsq, int N, int M, double popen, double pextend, double lambda, ESL_SCOREMATRIX *SMX, double **ret_zscore)
{

	int status;
	int i,j;      /* indexes through query and target sequences */

	double zscore;

	/* Don't need the full matrices to compute the zscore (thanks gotoh!) */

	/* ALLOCATION */


	/* INITIALIZATION */

	/* RECURSION */

	/* TERMINATION */

	*ret_zscore = zscore;

	return eslOK;

}
