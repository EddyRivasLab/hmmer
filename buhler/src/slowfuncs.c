#include "structs.h"
#include "funcs.h"

/*
 * We have to have to provide a declaration for this function for the
 * the code to compile, but since we don't optimize anything, we jut need to
 * call the default P7Logoddsify().  We can't call that function directly, 
 * howwever, since the memory for it probably hasn't been allocated yet.
 * Thus, we just dispatch back to  REQUIRE_P7LOGODDS, which will allocate 
 * the memory for us, and then call P7Logoddsify().  We inline this
 * function to avoid paying a performance penalty for making a useless
 * function call. - CRS 20 June 2005
 */
inline void
Logoddsify(struct plan7_s *hmm){
  REQUIRE_P7LOGODDS(hmm);
}

/*
 * This is the "slow" Viterbi algorithm for non-optimized versions
 * of the application.  We just dispatch back to the P7Viterbi() 
 * algorithm defined in core_algorihms.c.  We have to provide this
 * shell here, however, to support customizability.  We make it inlined
 * however, so that we don't pay a performance penalty for the useless
 * function call. - CRS 20 June 2005
 */
inline float
Viterbi(unsigned char *dsq, int L, struct plan7_s *hmm, 
	  struct dpmatrix_s *mx, struct p7trace_s **ret_tr)
{
  P7Viterbi(dsq, L, hmm, mx, ret_tr);
}

