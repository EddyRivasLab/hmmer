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
 *            This gets called anytime we're about to do P7Viterbi().
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
 * Args:      
 *
 * Returns:   TRUE if we can run P7Viterbi(); FALSE if we need
 *            to use a small memory variant.
 *
 * Xref:      STL7 p.122.
 */
extern int   ViterbiSpaceOK(int L, int M, cust_dpmatrix_s *mx){

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
extern float Backward(unsigned char *dsq, int L, struct plan7_s *hmm,	
		      cust_dpmatrix_s **ret_mx){

}

/* Function: P7Forward()
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
extern float Forward(unsigned char *dsq, int L, struct plan7_s *hmm, 
		     cust_dpmatrix_s **ret_mx){

}

/* Function: Viterbi()
 *
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           Identical to Forward() except that max's
 *           replace sum's. 
 *           
 *           This is the slower, more understandable version
 *           of P7Viterbi(). The default version in fastfuncs.c
 *           is portably optimized and more difficult to understand;
 *           the ALTIVEC version in altivecfuncs.c is vectorized
 *           with Altivec-specific code, and is pretty opaque.
 *           
 *           This function is not enabled by default; it is only
 *           activated by -DSLOW at compile time.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - reused DP matrix
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score
 */
extern float Viterbi(unsigned char *dsq, int L, struct plan7_s *hmm, 
		     cust_dpmatrix_s *mx, struct p7trace_s **ret_tr){

}

/* Function: ViterbiTrace()
 * Date:     SRE, Sat Aug 23 10:30:11 1997 (St. Louis Lambert Field) 
 *
 * Purpose:  Traceback of a Viterbi matrix: i.e. retrieval 
 *           of optimum alignment.
 *           
 * Args:     hmm    - hmm, log odds form, used to make mx
 *           dsq    - sequence aligned to (digital form) 1..N  
 *           N      - length of seq
 *           mx     - the matrix to trace back in, N x hmm->M
 *           ret_tr - RETURN: traceback.
 *           
 * Return:   (void)
 *           ret_tr is allocated here. Free using P7FreeTrace().
 */
extern void  ViterbiTrace(struct plan7_s *hmm, unsigned char *dsq, int N,
			  cust_dpmatrix_s *mx, struct p7trace_s **ret_tr){

}
