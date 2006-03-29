/*
 * JDBSTRUCTS.C
 * Structure management code for JB implementation of HMMer for x86
 *
 * For information, please contact jbuhler@cse.wustl.edu.
 */

#include "config.h"
#include "plan7.h"
#include "structs.h"

/*
 * Function: AllocLogoddsShell()
 * Date:     CRS, 10 June 2005 [J. Buhler's student, St. Louis]
 * 
 * Purpose:  Allocates the memory needed by the shell of the logodds_s
 *           structure.  Called from within AllocPlan7Shell().
 *
 * Args:     hmm - the profile hmm that contains the logodds_s structure
 *
 * Returns:  (void)
 *           Memory is allocated for the shell of hmm->lom.  Freed 
 *           from FreeLogodds(), which is called from within FreePlan7().
 */
inline void AllocLogoddsShell(struct plan7_s *hmm)
{ 
  hmm->lom = (struct logodds_s *) malloc(sizeof(struct logodds_s));
}

/*
 * Function: AllocLogoddsBody()
 * Date:     CRS, 10 June 2005 [J. Buhler's student, St. Louis]
 *
 * Purpose:  Allocates the memory needed by the logodds_s strcuture.
 *           Called from within AllocPlan7Body().
 *
 * Args:     hmm - the profile hmm that contains the logodds_s structure
 *
 * Returns:  (void)
 *           Memory is allocated for  hmm->lom.  Freed from FreeLogodds(), 
 *           which is called from within FreePlan7().
 */
/* JB: user must know about MAXCODE, as well as hmm->M, to do this right */
inline void AllocLogoddsBody(struct plan7_s *hmm)
{
  unsigned int a;
  int *miscBase;
  
  struct logodds_s *lom = hmm->lom;
  
  lom->tsc = malloc((hmm->M + 1) * N_TMOVES * sizeof(int));
  
  lom->misc = malloc(MAXCODE * sizeof(int *));
  miscBase = malloc(MAXCODE * (hmm->M + 1) * N_MISC * sizeof(int));
  for (a = 0; a < MAXCODE; a++)
    {
      lom->misc[a] = miscBase + a * (hmm->M + 1) * N_MISC;
    }
}

/*
 * Function: FreeLogodds()
 * Date:     CRS, 10 June 2005 [J. Buhler's student, St. Louis]
 *
 * Purpose:  Frees ALL memory (shell and body) used by the logodds_s
 *           structure.  Called from FreePlan7().
 * 
 * Args:     hmm - the profile hmm that contains the logodds_s structure
 *
 * Returns:  (void)
 *           Memory used by hmm->lom is freed.  
 */
inline void FreeLogodds(struct plan7_s *hmm)
{
  struct logodds_s *lom = hmm->lom;
  
  free(lom->tsc);
  
  {
    int *miscBase = lom->misc[0];
    free(miscBase);
  }
  
  free(lom->misc);
  free(lom);
}


/*
 * Function: FillCustomLogodds()
 * Date:     CRS, 13 July 2005 [J. Buhler's student, St. Louis]
 *
 * Purpose:  Fills the custom logodds_s structure using the data from the
 *           profile hmm.  One can assume all non-customized values of the 
 *           hmm have been filled in before this function is called.  No
 *           memory should be allocated here; that should be done in
 *           AllocLogoddsBody().
 *
 * Args:     hmm - the profile hmm that contains the logodds_s structure
 *
 * Returns:  (void)
 */
inline void FillCustomLogodds(struct plan7_s *hmm)
{
  struct logodds_s *lom = hmm->lom;
  unsigned int k, a;
  
  /* we only use values 0..M-1; 0 values are always -INFTY */
  for (k = 0; k < hmm->M; k++)
    {
      int *tsck = lom->tsc + k * N_TMOVES;
      tsck[_TMM] = hmm->tsc[TMM][k];
      tsck[_TMI] = hmm->tsc[TMI][k];
      tsck[_TMD] = hmm->tsc[TMD][k];
      tsck[_TDM] = hmm->tsc[TDM][k];
      tsck[_TDD] = hmm->tsc[TDD][k];
      tsck[_TIM] = hmm->tsc[TIM][k];
      tsck[_TII] = hmm->tsc[TII][k];
    }
  
  /* we only use values 1..M */
  for (k = 1; k <= hmm->M; k++)
    {
      int *tsck = lom->tsc + k * N_TMOVES;
      
      tsck[_BSC] = hmm->bsc[k];
      tsck[_ESC] = hmm->esc[k];
    }
  
  for (a = 0; a < MAXCODE; a++)
    {
      /* we only use values 1..M */
      for (k = 1; k <= hmm->M; k++)
        {
          int *misck = lom->misc[a] + k * N_MISC;
          
          misck[_MSC] = hmm->msc[a][k];
          misck[_ISC] = hmm->isc[a][k];
        }
    }
  
  lom->xsc[_XTE][_LOOP] = hmm->xsc[XTE][LOOP];  
  lom->xsc[_XTE][_MOVE] = hmm->xsc[XTE][MOVE];
  lom->xsc[_XTN][_LOOP] = hmm->xsc[XTN][LOOP];
  lom->xsc[_XTN][_MOVE] = hmm->xsc[XTN][MOVE];
  lom->xsc[_XTC][_LOOP] = hmm->xsc[XTC][LOOP];
  lom->xsc[_XTC][_MOVE] = hmm->xsc[XTC][MOVE];
  lom->xsc[_XTJ][_LOOP] = hmm->xsc[XTJ][LOOP];
  lom->xsc[_XTJ][_MOVE] = hmm->xsc[XTJ][MOVE];
}


/* Function: CreateDPMatrix()
 *
 * Purpose:  Create a dynamic programming matrix for standard Forward,
 *           Backward, or Viterbi, with scores kept as scaled log-odds
 *           integers. Keeps 2D arrays compact in RAM in an attempt 
 *           to maximize cache hits. 
 *           
 *           The structure can be dynamically grown, if a new
 *           HMM or seq exceeds the currently allocated size. Dynamic
 *           growing is more efficient than an alloc/free of a whole
 *           matrix for every new target. The ResizeDPMatrix()
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
inline cust_dpmatrix_s *CreateDPMatrix(int N, int M, int padN, int padM)
{
  cust_dpmatrix_s *mx = malloc(sizeof(struct dpmatrix_s));
  
  mx->maxM = M;
  mx->maxN = N;
  mx->padM = padM;
  mx->padN = padN;
  
  mx->dp  = malloc((mx->maxN + 1) * sizeof(int *));
  mx->xmx = malloc((mx->maxN + 1) * N_XSTATES * sizeof(int));
  
  {
    int *dpbase  = malloc((mx->maxN + 1) * (mx->maxM + 1) * 
                          N_MSTATES * sizeof(int));
    unsigned int i;
    
    for (i = 0; i <= mx->maxN; i++)
      {
        mx->dp[i]  = dpbase + i * (mx->maxM + 1) * N_MSTATES;
      }
  }
  
  return mx;
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
 *           in size after the first call to ResizeDPMatrix();
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
 * Args:     mx    - an already allocated model to grow.
 *           N     - seq length to allocate for; N+1 rows
 *           M     - size of model
 *                   
 * Return:   (void)
 *           mx is (re)allocated here.
 */
inline void ResizeDPMatrix(cust_dpmatrix_s *mx, int N, int M)
{
  if (M > mx->maxM || N > mx->maxN)
    {
      int *dpbase = mx->dp[0];
      unsigned int i;
      
      free(dpbase);
      
      if (M > mx->maxM) 
        mx->maxM = M + mx->padM;
      
      if (N > mx->maxN) 
        {
          mx->maxN = N + mx->padN;
          
          free(mx->dp);
          free(mx->xmx);
          
          mx->dp  = malloc((mx->maxN + 1) * sizeof(int *));
          mx->xmx = malloc((mx->maxN + 1) * N_XSTATES * sizeof(int));
        }
      
      dpbase = malloc((mx->maxN + 1) * (mx->maxM + 1) * 
                      N_MSTATES * sizeof(int));
      
      for (i = 0; i <= mx->maxN; i++)
        {
          mx->dp[i]  = dpbase + i * (mx->maxM + 1) * N_MSTATES;
        }
    }
}


/* Function: FreeDPMatrix()
 *
 * Purpose:  Free the dynamic programming matrix allocated by CreateDPMatrix().
 * 
 * Return:   (void)
 */
inline void FreeDPMatrix(cust_dpmatrix_s *mx)
{
  free(mx->dp[0]); /* block ptrs */
  free(mx->xmx);  
  
  free(mx->dp);    /* ptr arrays */
  
  free(mx);
}
