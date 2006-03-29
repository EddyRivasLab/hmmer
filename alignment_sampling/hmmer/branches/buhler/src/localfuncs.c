#include "plan7.h"
#include "structs.h"
#include "funcs.h"

/* Function:  ViterbiSpaceOK()
 * Incept:    SRE, Wed Oct  1 12:53:13 2003 [St. Louis]
 *
 * Purpose:   Returns TRUE if the existing matrix allocation
 *            is already sufficient to hold the requested MxN, or if
 *            the matrix can be expanded in M and/or N without
 *            exceeding RAMLIMIT megabytes. 
 *            
 *            This gets called anytime we're about to do Viterbi().
 *            If it returns FALSE, we switch into the appropriate
 *            small-memory alternative: P7SmallViterbi() or
 *            P7WeeViterbi().
 *            
 *            Checking the DP problem size against P7ViterbiSize()
 *            is not enough, because the DP matrix may be already
 *            allocated in MxN. For example, if we're already
 *            allocated to maxM,maxN of 1447x979, and we try to
 *            do a problem of MxN=12x125000, P7ViterbiSize() may
 *            think that's fine - but when we resize, we can only
 *            grow, so we'll expand to 1447x125000, which is 
 *            likely over the RAMLIMIT. [bug#h26; xref SLT7 p.122]
 *
 * Args:      L  - length of sequence
 *            M  - length of HMM
 *            mx - an allocated model
 *
 * Returns:   TRUE if we can run Viterbi(); FALSE if we need
 *            to use a small memory variant.
 *
 * Xref:      STL7 p.122.
 */
int   ViterbiSpaceOK(int L, int M, cust_dpmatrix_s *mx){

}

/*
 * Function: DispatchViterbi()
 * Date:     CRS, Wed 18 Aug 2005 [J. Buhler's student, St. Louis]
 *
 * Purpose:  Determines the appropriate Viterbi algorithm to call,
 *           based on the values of the parameters provided.  
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - reused dp matrix
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted,
                      see below
 *           need_trace - true if traceback is needed, false otherwise
 *
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float DispatchViterbi(unsigned char *dsq, int L, struct plan7_s *hmm,
		      cust_dpmatrix_s *mx, struct p7trace_s **ret_tr,
		      int need_trace)
{
  
}


/* Function: Viterbi()
 *
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           Identical to Forward() except that max's
 *           replace sum's.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - reused DP matrix
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
float Viterbi(unsigned char *dsq, int L, struct plan7_s *hmm, 
	      cust_dpmatrix_s *mx, struct p7trace_s **ret_tr)
{

}


/* Function: Forward()
 *
 * Purpose:  The Forward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float Forward(unsigned char *dsq, int L, struct plan7_s *hmm, 
	      cust_dpmatrix_s **ret_mx)
{
  
}


/* Function: Backward()
 * 
 * Purpose:  The Backward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ILogsum(); this is a slow but robust approach.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float Backward(unsigned char *dsq, int L, struct plan7_s *hmm,	
	       cust_dpmatrix_s **ret_mx)
{
  
}
