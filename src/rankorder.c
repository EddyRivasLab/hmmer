/* rank order statistics
 *
 * Compute pvalues using rank order statistics
 * for the alternative score systems in shmmer
 *
 *  1. Rank order function.
 *  2. Private functions.
 *
 * SC, Mon Apr  5 16:02:15 EDT 2010 [Janelia] [The Man of Stone]
 * $Id$
 */

#include "p7_config.h"

#include "easel.h"
#include "hmmer.h"

#define LEFT 0                                    /* first array index of random scores */
#define RIGHT p7_RANKORDER_LIMIT - 1              /* last array index of random scores  */

static int binary_search(float random_scores[], int low, int high, float sc);

/*****************************************************************
 * 1. Rank order function
 *****************************************************************/

/* Function:  rank_order()
 * Incept:    SC, Mon Apr  5 16:02:15 EDT 2010 [Janelia]
 *
 * Purpose:   Compute rank order p-value for a score <sc>.
 *
 * Args:      random_scores - array of scores from random sequences
 *						sc            - score to compute p-value
 *
 * Returns:   p-value on success;
 *
 * Throws:    Nothing
 */
extern double
rank_order(float random_scores[], float sc)
{
	int          num;                               /* number of random sequence scores >= sc */
	double       P = 0;                             /* P-value of a score                     */

	num = binary_search(random_scores, LEFT, RIGHT, sc);
	if (num == -eslFAIL) esl_fatal("Failed to compute rank order statistics");

	/* Calculate P-value */
	P = (double)(num + 1) / p7_RANKORDER_LIMIT; /* P(sc >= t) */

	return P;
}

/*****************************************************************
 * 1. Private functions
 *****************************************************************/

/* Function:  binary_search()
 * Incept:    SC, Mon Apr  5 16:02:15 EDT 2010 [Janelia]
 *
 * Purpose:   Compute the number of random scores in array
 * 						<random_scores> at or above a given score <sc>.
 *
 * Args:      random_scores - array of scores from random sequences
 * 						left          - lowest array index (usually 0)
 * 						right         - highest array index;
 * 													  set to p7_RANKORDER_LIMIT - 1
 *						sc            - score in need of a p-value
 *
 * Returns:   index of score in <random_scores> at or above <sc> on success;
 * 						<-eslFAIL> on failure
 *
 * Throws:    Nothing
 */
static int
binary_search(float random_scores[], int left, int right, float sc)
{

	int k;

	/* Check whether our score is greater or lower
	 * than any score from random sequences
	 * Return array index
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
		else if (sc < random_scores[mid]) left  = mid + 1;
		else if (sc > random_scores[mid]) right = mid - 1;
		else if (sc == random_scores[mid] && sc == random_scores[mid+1])  /* sc is within a patch of duplicate scores */
		{
			for (k = 2; k <= right; k++)
				if (sc > random_scores[mid + k]) return (mid + k - 1); /* find the rightmost duplicate score */
		}
		else return mid;                                           /* sc = random_scores[mid] */
	}

	return -eslFAIL;
}
