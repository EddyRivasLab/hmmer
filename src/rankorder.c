/* rank order statistics
 *
 * Compute pvalues using rank order statistics
 * for the alternative score systems in shmmer
 *
 * SC, Mon Apr  5 16:02:15 EDT 2010 [Janelia]
 */

#include "p7_config.h"

#include "easel.h"
#include "hmmer.h"

#define LEFT 0
#define RIGHT p7_RANKORDER_LIMIT

static int binary_search(float random_scores[], int low, int high, float sc);

extern double
rank_order(FILE *rsfp, float sc)
{
	static float random_scores[p7_RANKORDER_LIMIT]; /* scores from random seqs. Array order is from high (left) to low (right) scores
																									 * use static array to keep data out of the stack
																									 * A dynamic (growable) array is a better option here
																									 * We do not need to set a limit and the space is only
																									 * allocated if we do use rank order statistics
																									 */
	int          num;                               /* number of random sequence scores >= sc */
	double       P;                                 /* P-value of a score                     */
	int          i = 0;

	/* read forward scores from random sequences */
	while (fscanf(rsfp, " %f", &random_scores[i]) != EOF) ++i;

	num = binary_search(random_scores, LEFT, RIGHT - 1, sc);
	if (num == -eslFAIL) esl_fatal("Failed to compute rank order statistics");

	/* Calculate P-value */
	P = (double)(num + 1) / RIGHT; /* P(sc >= t) */

	return P;
}

static int
binary_search(float random_scores[], int left, int right, float sc)
{

	/* Check whether our score is greater or lower
	 * than any score from random sequences
	 */
	if (sc > random_scores[left])       return left;   /* left = 0                         , P-value < 1/length(random_scores) */
	else if (sc < random_scores[right]) return right;  /* right = length(random_scores) - 1, P-value = 1                       */

	while (left <= right)
	{
		int mid = left + (right - left)/2; /* or (left + right)/2 */

		/* Check whether our score is between two
		 * consecutive scores from random sequences
		 */
		if      (sc < random_scores[mid] && sc > random_scores[mid+1]) return mid;   /* sc between array values */
		else if (sc > random_scores[mid] && sc < random_scores[mid-1]) return mid-1; /* sc between array values */

		/* Check whether we need to keep partitioning
		 * the scores from random sequences
		 */
		else if (sc < random_scores[mid]) left  = mid;
		else if (sc > random_scores[mid]) right = mid;
		else return mid;                               /* sc = random_scores[mid] */
	}

	return -eslFAIL;
}
