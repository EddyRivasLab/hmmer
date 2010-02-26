/*
 * Miyazawa's (1994) partition function
 *
 * Changes to the original algorithm:
 *
 */

#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_scorematrix.h"

#include "hmmer.h"

int
pfunction(ESL_DSQ *qsq, int N, ESL_DSQ *dbsq, int M, double popen, double pextend, double lambda, ESL_SCOREMATRIX *SMX, double *ret_zscore)
{

//	int status;
	int i,j;      /* indexes through query and target sequences */
	double zscore;

	ESL_DMATRIX *MMX = NULL;
	ESL_DMATRIX *IMX = NULL;
	ESL_DMATRIX *JMX = NULL;  /* Caller should provide these matrices in the future */

	/* ALLOCATION */ /* This is wasteful 'cos don't need the full matrices to compute the zscore (thanks gotoh!) */
	MMX = esl_dmatrix_Create(N+1,M+1);
	if (MMX == NULL) esl_fatal("MMX not allocated!!!\n");

	IMX = esl_dmatrix_Create(N+1,M+1);
	if (IMX == NULL) esl_fatal("IMX not allocated!!!\n");

	JMX = esl_dmatrix_Create(N+1,M+1);
	if (JMX == NULL) esl_fatal("JMX not allocated!!!\n");

	/* INITIALIZATION */
	MMX->mx[0][0] = 1;
	IMX->mx[0][0] = 0;
	JMX->mx[0][0] = 0;

	for (i = 1; i <= N; i++)  MMX->mx[i][0] = exp(-lambda * popen) * pow(exp(-lambda * pextend), (double)i-1);  /* affine gap cost */
	for (j = 1; j <= M; j++)  MMX->mx[0][j] = exp(-lambda * popen) * pow(exp(-lambda * pextend), (double)j-1);  /* affine gap cost */

	for (i = 1; i <= N; i++)  IMX->mx[i][0] = exp(-lambda * popen) * pow(exp(-lambda * pextend), (double)i-1);  /* affine gap cost */
	for (j = 1; j <= M; j++)  IMX->mx[0][j] = 0;                                                                /* so that it starts with gap open in i=1 */

	for (i = 1; i <= N; i++)  JMX->mx[i][0] = 0;                                                                /* so that it starts with gap open in j=1 */
	for (j = 1; j <= M; j++)  JMX->mx[0][j] = exp(-lambda * popen) * pow(exp(-lambda * pextend), (double)j-1);  /* affine gap cost */

	/* RECURSION */
	for (i = 1; i <= N; i++)            /* X sequence */
		for (j = 1; j <= M; j++)        /* Y sequence */
		{

		/* IMX */
			IMX->mx[i][j] = MMX->mx[i-1][j] * exp(-lambda * popen) + IMX->mx[i-1][j] * exp(-lambda * pextend);

		/* JMX */
			JMX->mx[i][j] = MMX->mx[i][j-1] * exp(-lambda * popen) + JMX->mx[i][j-1] * exp(-lambda * pextend);

		/* MMX */
			MMX->mx[i][j] = (MMX->mx[i-1][j-1] + IMX->mx[i-1][j-1] + JMX->mx[i-1][j-1]) * exp(lambda * SMX->s[qsq[i]][dbsq[j]]); /* CHECK SMX INDEXES!!! */

		} /* END RECURSION */

	/* TERMINATION */
	zscore = MMX->mx[N][M] + IMX->mx[N][M] + JMX->mx[N][M];

	*ret_zscore = zscore;

	/* CLEANUP */
	esl_dmatrix_Destroy(MMX); MMX = NULL;
	esl_dmatrix_Destroy(IMX); IMX = NULL;
	esl_dmatrix_Destroy(JMX); JMX = NULL;

	return eslOK;

}
