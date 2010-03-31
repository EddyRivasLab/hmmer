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

#include "easel.h"
#include "esl_alphabet.h"
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
 *            sm     - query sequence (of length M = sm->n) and score system
 *            gx     - DP matrix with room for an LxM alignment
 *            ret_sc - RETURN: Miyazawa score
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
	double     **MSC  = sm->Q->mx;
	double    lopen   = sm->lopen;
	double  lextend   = sm->lextend;
	float        sc   = 0;
	float         Z   = 0;

	/* Log sum initialization*/
	p7_FLogsumInit();

	/* Initialization */
	MMX(0,0) = IMX(0,0) = DMX(0,0) = 0;

	for (k = 1; k <= M; k++)
	{
		MMX(0,k) = IMX(0,k) = DMX(0,k) = -eslINFINITY; /* first row */
	}

	/* Recursion */
	for (i = 1; i <= L; i++)              /* loop over target sequence */
	{
		/* Initialization */
		MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY; /* first column */

		for (k = 1; k <= M; k++)            /* loop over query sequence */
		{
			/* match */
			sc = p7_FLogsum(p7_FLogsum(0, MMX(i-1,k-1)), p7_FLogsum(IMX(i-1,k-1), DMX(i-1,k-1))); /* e^0 = 1 starts a new local alignment in each match cell */
			MMX(i,k) = sc + (float)MSC[dsq[i]][sm->dsq[k]];

			/* insert (in the target) */
			IMX(i,k) = p7_FLogsum(MMX(i-1,k) - lopen, IMX(i-1,k) - lextend);

			/* delete (in the target) */
			DMX(i,k) = p7_FLogsum(MMX(i,k-1) - lopen, DMX(i,k-1) - lextend);

			Z = p7_FLogsum(Z, MMX(i,k));                                                          /* we sum the score of each local alignment terminating in each match cell */

		} /* end loop over query sequence */

	} /* end loop over target sequence */

	*ret_sc = Z;

	return eslOK;
}
/*------------- end: miyazawa ------------------*/

/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
/*----------------- end, benchmark ------------------------------*/

/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_MIYAZAWA_TESTDRIVE
#include <string.h>
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* Optimal local alignment
 * with Miyazawa scores
 */
int
miyazawa_optimal(const ESL_DSQ *dsq, int L, P7_SCORESYS *sm, P7_GMX *gx, float *ret_sc)
{
	float      **dp   = gx->dp;
	int          M    = sm->n;
	int          i;                    /* index over rows (target)   */
	int 				 k;						         /* index over columns (query) */
	double    lopen   = sm->lopen;
	double  lextend   = sm->lextend;
	double     **MSC  = sm->Q->mx;
	float       sc;
	float      fsc    = -eslINFINITY;

	/* Initialization */
	for (k = 0; k <= M; k++)
	{
		MMX(0,k) = 0; 											  /* first row */
		IMX(0,k) = DMX(0,k) = -eslINFINITY;   /* first row. Because this is an optimal alignment algorithm */
		                                      /* DMX(0,k) can also be sopen + (k-1)sextend                 */
	}

	/* Recursion */
	for (i = 1; i <= L; i++)              /* loop over target sequence */
	{
		/* Initialization */
		MMX(i,0) = 0;                       /* first column */
		IMX(i,0) = DMX(i,0) = -eslINFINITY; /* first column. Because this is an optimal alignment algorithm */
		                                    /* IMX(i,0) can also be sopen + (i-1)sextend                    */

		for (k = 1; k <= M; k++)            /* loop over query sequence */
		{
			/* match */
			sc = ESL_MAX(0, MMX(i-1,k-1));
			sc = ESL_MAX(sc, IMX(i-1,k-1));
			sc = ESL_MAX(sc, DMX(i-1,k-1));
			MMX(i,k) = sc + (float)MSC[dsq[i]][sm->dsq[k]];

			/* optimal score (in MMX) */
			if (MMX(i,k) > fsc) fsc = MMX(i,k);

			/* insert (in the target) */
			IMX(i,k) = ESL_MAX((MMX(i-1,k) - lopen), (IMX(i-1,k) - lextend));

			/* delete (in the target) */
			DMX(i,k) = ESL_MAX((MMX(i,k-1) - lopen), (DMX(i,k-1) - lextend));

		} /* end loop over query sequence */

	} /* end loop over target sequence */

	*ret_sc = fsc;

	return eslOK;
}

/* Miyazawa is hard to validate.
 * We do know that the total Miyazawa score is >= optimal Miyazawa score
 */
static void
utest_miyazawa(ESL_GETOPTS *go, ESL_RANDOMNESS *r, const double *p, P7_SCORESYS *sm, int nseq, int L)
{

	ESL_DSQ  *dsq     = NULL;       /* target sequence                              */
	P7_GMX   *ogx     = NULL;       /* Dynamic matrix for optimal scores            */
	P7_GMX   *pgx     = NULL;       /* Dynamic matrix for partition funciton scores */
	int       sstatus;              /* search status            */
	float    osc;                   /* optimal score            */
	float    psc;                   /* partition function score */
  int       i;

	/* Create dynamic programming matrices */
	if ((ogx = p7_gmx_Create(sm->n, L)) == NULL) esl_fatal("Dynamic programming matrix allocation failure");
	if ((pgx = p7_gmx_Create(sm->n, L)) == NULL) esl_fatal("Dynamic programming matrix allocation failure");

	/* Create target sequence object */
	if ((dsq = malloc(sizeof(ESL_DSQ) * (L+2))) == NULL) esl_fatal("allocation failed");

  /* Loop over target sequences */
  for (i = 0; i < nseq; i++)
  {
  	/* Create target sequence */
  	if (esl_rsq_xIID(r, p, sm->abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");

  	/* Compute local optimal alignment with miyazawa scores */
  	sstatus = miyazawa_optimal(dsq, L, sm, ogx, &osc);
  	if (sstatus != eslOK)  esl_fatal ("Failed to compute the optimal score!\n");

  	/* Compute partition function with miyazawa scores */
  	sstatus = p7_GMiyazawa(dsq, L, sm, pgx, &psc);
  	if (sstatus != eslOK)  esl_fatal ("Failed to compute the partition function score!\n");

  	printf("Scores: %f %f\n", osc, psc);
  	if (psc < osc) esl_fatal("Partition function score can't be less than the optimal score");
  }

  p7_gmx_Destroy(ogx);
  p7_gmx_Destroy(pgx);
  free(dsq);
  return;
}
#endif
/*------------------------- end, unit tests ---------------------*/

/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
/* gcc -g -Wall -Dp7GENERIC_MIYAZAWA_TESTDRIVE -I. -I../easel -L. -L../easel -o generic_miyazawa_utest generic_miyazawa.c -lhmmer -leasel -lm
 */
#ifdef p7GENERIC_MIYAZAWA_TESTDRIVE
#include "easel.h"
#include "esl_getopts.h"
#include "esl_scorematrix.h"

#include "p7_config.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be very verbose",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options]";
static char banner[] = "unit test driver for the generic Miyazawa implementation";

int
main(int argc, char **argv)
{
	ESL_GETOPTS    *go   = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
	ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
	P7_SCORESYS    *sm   = NULL;
	double sopen         = 11;
	double sextend       = 1;
  int             M    = 400;   /* query length  */
  int             L    = 1600;   /* target length */
  int             nseq = 1000;
  int             i, k, l;
  double         *p;
  int             status;
  float     expected_s;
  float     expected_m;

  /* Create score system (query + scores) */
  ESL_ALLOC(sm, sizeof(*sm));
  if ((sm->abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal("failed to create alphabet");

  /* Allocate random query sequence */
  sm->n = M;
  if ((sm->dsq = malloc(sizeof(ESL_DSQ) * (M+2))) == NULL) esl_fatal("allocation failed");

  /* Create background frequency vector */
  if ((p = malloc(sizeof(double) * (20))) == NULL) esl_fatal("allocation failed");
  esl_composition_SW50(p);           /* set background frequencies to Swiss-Prot 50.8 */

  /* Create random query sequence */
  if (esl_rsq_xIID(r, p, sm->abc->K, sm->n, sm->dsq) != eslOK) esl_fatal("query sequence generation failed");

  /* Set score system */
 	status = p7_builder_SetMiyScoreSystem(sm, NULL, NULL, sopen, sextend); /* will default to BLOSUM62 scores */
 	if (status != eslOK) esl_fatal("Failed to set single seq score system:\n%s\n", sm->errbuf);

  /* Run utest */
	utest_miyazawa(go, r, p, sm, nseq, L);

//	/* Expected S/W score */
//	for(k = 0; k <= sm->abc->K-1; k++)
//		for(l = 0; l <= sm->abc->K-1; l++)
//			expected_s += p[k] * p[l] * sm->S->s[k][l];
//
//	printf("Expected S/W: %f\n", expected_s);
//
//	/* Expected Miyazawa score */
//	for(k = 0; k <= sm->abc->K-1; k++)
//		for(l = 0; l <= sm->abc->K-1; l++)
//			expected_m += p[k] * p[l] * sm->Q->mx[k][l];

//	printf("Expected Miy: %f\n", expected_m);

	/* Clean up */
  esl_alphabet_Destroy(sm->abc);
  free(sm->dsq);
  esl_scorematrix_Destroy(sm->S);
  if (sm->Q != NULL) free(sm->Q);
  free(sm);
  esl_randomness_Destroy(r);

  return 0;

  ERROR:
  printf("allocation failed");
  exit(1);
}
#endif /* p7GENERIC_MIYAZAWA_TESTDRIVE */
/*-------------------- end, test driver -------------------------*/


/*****************************************************************
 * 5. Example driver
 *****************************************************************/
/*-------------------- end, example -----------------------------*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/


