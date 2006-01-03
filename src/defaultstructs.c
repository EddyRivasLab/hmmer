#include "structs.h"
#include "funcs.h"

/*
 * Note:  The default implementation doesn't use any custom data
 *        structures, so we these function bodies are empty.  We
 *        still have to provide them, however, for the program to
 *        compile successfully.  - CRS 20 June 2005
 */
void AllocLogoddsShell(struct plan7_s *hmm){ hmm->lom = NULL; }
void AllocLogoddsBody(struct plan7_s *hmm){}
void FreeLogodds(struct plan7_s *hmm){}
void FillCustomLogodds(struct plan7_s *hmm){}


/*
 * Note: Since the cust_dpmatrix_s structure doesn't differ from the 
 *       default dpmatrix_s structure, we just dispatch to the 
 *       functions for the default dpmatrix_s structure.  We still
 *       have to provide these definitions, however, to get things
 *       to compile under the architecture. - CRS 15 Aug 2005
 */
cust_dpmatrix_s *
CreateDPMatrix(int N, int M, int padN, int padM)
{
  return (cust_dpmatrix_s*)CreatePlan7Matrix(N, M, padN, padM);
}

void
ResizeDPMatrix(cust_dpmatrix_s *mx, int N, int M)
{
  ResizePlan7Matrix(mx, N, M, NULL, NULL, NULL, NULL);
  return;
}

void
FreeDPMatrix(cust_dpmatrix_s *mx)
{
  FreePlan7Matrix(mx);
}
