#include "structs.h"
#include "funcs.h"

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
		      int need_trace){
  float sc = 0.0;

  if (ViterbiSpaceOK(L, hmm->M, mx))
    {
      SQD_DPRINTF1(("   ... using normal Viterbi(); size ~%d MB\n", 
		    P7ViterbiSize(L, hmm->M)));
      sc = Viterbi(dsq, L, hmm, mx, ret_tr);
    }
  else
    {
      SQD_DPRINTF1(("   ... using P7SmallViterbi(); size ~%d MB\n",
		    P7ViterbiSize(L, hmm->M)));
      sc = P7SmallViterbi(dsq, L, hmm, mx, ret_tr);
    }
  
  return sc;
}
