#include "structs.h"
#include "funcs.h"

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

