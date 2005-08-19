#include "structs.h"

/*
 * Function: AllocLogoddsShell()
 * Date:     CRS, 10 June 2005
 * 
 * Purpose:  Allocates the memory needed by the shell of the logodds_s
 *           structure.  Called from within AllocPlan7Shell().
 *
 * Args:     hmm - the profile hmm containing the logodds_s structure
 *
 * Returns:  (void)
 *           Memory is allocated for the shell of hmm->lom.  Freed 
 *           from FreeLogodds(), which is called from within FreePlan7().
 */
inline void AllocLogoddsShell(struct plan7_s *hmm){ 
  /*
   * Uncomment the following line if your implementation does not use
   * the logodds_s structure.
   */
  /* hmm->lom = NULL; */
}

/*
 * Function: AllocLogoddsBody()
 * Date:     CRS, 10 June 2005
 *
 * Purpose:  Allocates the memory needed by the logodds_s strcuture.
 *           Called from within AllocPlan7Body().
 *
 * Args:     hmm - the profile hmm containing the logodds_s structure
 *
 * Returns:  (void)
 *           Memory is allocated for  hmm->lom.  Freed from FreeLogodds(), 
 *           which is called from within FreePlan7().
 */
inline void AllocLogoddsBody(struct plan7_s *hmm){

}

/*
 * Function: FreeLogodds()
 * Date:     CRS, 10 June 2005
 *
 * Purpose:  Frees ALL memory (shell and body) used by the logodds_s
 *           structure.  Called from FreePlan7().
 * 
 * Args:     hmm - the profile hmm containing the logodds_s structure
 *
 * Returns:  (void)
 *           Memory used by hmm->lom is freed.  
 */
inline void FreeLogodds(struct plan7_s *hmm){

}

/*
 * Function: FillCustomLogodds()
 * Date:     CRS, 13 July 2005
 *
 * Purpose:  Fills the custom_s logodds structure using the data from the
 *           profile hmm.  One can assume all non-customized values of the 
 *           hmm have been filled in before this function is called.  No
 *           memory should be allocated here; that should be done in
 *           AllocLogoddsBody().
 *
 * Args:     hmm - the profile hmm containing the logodds_s structure
 *
 * Returns:  (void)
 */
inline void FillCustomLogodds(struct plan7_s *hmm){

}

/* Function: CreateDPMatrix()
 *
 * Purpose:  Create a dynamic programming matrix for standard Forward,
 *           Backward, or Viterbi, with scores kept as scaled log-odds
 *           integers. Keeps 2D arrays compact in RAM in an attempt 
 *           to maximize cache hits. 
 *           
 *           The mx structure can be dynamically grown, if a new
 *           HMM or seq exceeds the currently allocated size. Dynamic
 *           growing is more efficient than an alloc/free of a whole
 *           matrix for every new target. The ResizePlan7Matrix()
 *           call does this reallocation, if needed. Here, in the
 *           creation step, we set up some pads - to inform the resizing
 *           call how much to overallocate when it realloc's. 
 *           
 * Args:     N     - N+1 rows are allocated, for sequence.  
 *           M     - size of model in nodes
 *           padN  - over-realloc in seq/row dimension, or 0
 *           padM  - over-realloc in HMM/column dimension, or 0
 *                 
 * Return:   mx
 *           mx is allocated here. Caller frees with FreeDPMatrix(mx).
 */
inline struct dpmatrix_s *CreateDPMatrix(int N, int M, int padN, int padM){

}

/* Function: ResizeDPMatrix()
 * 
 * Purpose:  Reallocate a dynamic programming matrix, if necessary,
 *           for a problem of NxM: sequence length N, model size M.
 *           (N=1 for small memory score-only variants; we allocate
 *           N+1 rows in the DP matrix.) 
 *           
 *           We know (because of the way hmmsearch and hmmpfam are coded)
 *           that only one of the two dimensions is going to change
 *           in size after the first call to ResizePlan7Matrix();
 *           that is, for hmmsearch, we have one HMM of fixed size M
 *           and our target sequences may grow in N; for hmmpfam,
 *           we have one sequence of fixed size N and our target models
 *           may grow in M. What we have to watch out for is P7SmallViterbi()
 *           working on a divide and conquer problem and passing us N < maxN,
 *           M > maxM; we should definitely *not* reallocate a smaller N.
 *           Since we know that only one dimension is going to grow,
 *           we aren't scared of reallocating to maxN,maxM. (If both
 *           M and N could grow, we would be more worried.)
 *
 *           Returns individual ptrs to the four matrix components
 *           as a convenience.
 *           
 * Args:     mx    - an already allocated model to grow.
 *           N     - seq length to allocate for; N+1 rows
 *           M     - size of model
 *           xmx, mmx, imx, dmx 
 *                 - RETURN: ptrs to four mx components as a convenience
 *                   
 * Return:   (void)
 *           mx is (re)allocated here.
 */
inline void ResizeDPMatrix(struct dpmatrix_s *mx, int N, int M){

}

/* Function: FreeDPMatrix()
 *
 * Purpose:  Free a dynamic programming matrix allocated by CreateDPMatrix().
 * 
 * Return:   (void)
 */
inline void FreeDPMatrix(struct dpmatrix_s *mx){

}
