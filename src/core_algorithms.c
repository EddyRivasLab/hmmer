/************************************************************
 * @LICENSE@
 ************************************************************/

/* core_algorithms.c
 * SRE, Mon Nov 11 15:58:52 1996
 * CVS $Id$
 * 
 * Simple and robust "research" implementations of Forward, Backward,
 * and Viterbi for Plan7. For optimized replacements for some of these functions,
 * see fast_algorithms.c.
 */
#include "structs.h"
#include "config.h"
#include "funcs.h"
#include "squid.h"

#include <string.h>
#include <assert.h>

static float get_wee_midpt(struct plan7_s *hmm, char *dsq, int L, 
			   int k1, char t1, int s1,
			   int k3, char t3, int s3,
			   int *ret_k2, char *ret_t2, int *ret_s2);


#ifndef ALTIVEC
/* Function: CreatePlan7Matrix()
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
 *           call how much to overallocate when it realloc's. If a pad
 *           is zero, we will not resize in that dimension.
 *           
 * Args:     N     - N+1 rows are allocated, for sequence.  
 *           M     - size of model in nodes
 *           padN  - over-realloc in seq/row dimension, or 0
 *           padM  - over-realloc in HMM/column dimension, or 0
 *                 
 * Return:   mx
 *           mx is allocated here. Caller frees with FreePlan7Matrix(mx).
 */
struct dpmatrix_s *
CreatePlan7Matrix(int N, int M, int padN, int padM)
{
  struct dpmatrix_s *mx;
  int i;

  mx          = (struct dpmatrix_s *) MallocOrDie (sizeof(struct dpmatrix_s));
  mx->xmx     = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->mmx     = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->imx     = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->dmx     = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->xmx_mem = (void *) MallocOrDie (sizeof(int) * ((N+1)*5));
  mx->mmx_mem = (void *) MallocOrDie (sizeof(int) * ((N+1)*(M+2)));
  mx->imx_mem = (void *) MallocOrDie (sizeof(int) * ((N+1)*(M+2)));
  mx->dmx_mem = (void *) MallocOrDie (sizeof(int) * ((N+1)*(M+2)));

  /* The indirect assignment below looks wasteful; it's actually
   * used for aligning data on 16-byte boundaries as a cache 
   * optimization in the fast altivec implementation
   */
  mx->xmx[0] = (int *) mx->xmx_mem;
  mx->mmx[0] = (int *) mx->mmx_mem;
  mx->imx[0] = (int *) mx->imx_mem;
  mx->dmx[0] = (int *) mx->dmx_mem;
  for (i = 1; i <= N; i++)
    {
      mx->xmx[i] = mx->xmx[0] + (i*5); 
      mx->mmx[i] = mx->mmx[0] + (i*(M+2));
      mx->imx[i] = mx->imx[0] + (i*(M+2));
      mx->dmx[i] = mx->dmx[0] + (i*(M+2));
    }

  if (padM > 0 && padN > 0)  Die("there's trouble with RAMLIMIT if you grow in both M and N.");
  mx->maxN = N;
  mx->maxM = M;
  mx->padN = padN;
  mx->padM = padM;
  
  return mx;
}
#endif /*ALTIVEC*/

#ifndef ALTIVEC
/* Function: ResizePlan7Matrix()
 * 
 * Purpose:  Reallocate a dynamic programming matrix, if necessary,
 *           for a problem of NxM: sequence length N, model size M.
 *           (N=1 for small memory score-only variants; we allocate
 *           N+1 rows in the DP matrix.) 
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
void
ResizePlan7Matrix(struct dpmatrix_s *mx, int N, int M, 
		  int ***xmx, int ***mmx, int ***imx, int ***dmx)
{
  int i;

  if (N <= mx->maxN && M <= mx->maxM) goto DONE;
  
  if (N > mx->maxN) {
    N          += mx->padN; 
    mx->maxN    = N; 
    mx->xmx     = (int **) ReallocOrDie (mx->xmx, sizeof(int *) * (N+1));
    mx->mmx     = (int **) ReallocOrDie (mx->mmx, sizeof(int *) * (N+1));
    mx->imx     = (int **) ReallocOrDie (mx->imx, sizeof(int *) * (N+1));
    mx->dmx     = (int **) ReallocOrDie (mx->dmx, sizeof(int *) * (N+1));
  }

  if (M > mx->maxM) {
    M += mx->padM; 
    mx->maxM = M; 
  }

  mx->xmx_mem = (void *) ReallocOrDie (mx->xmx_mem, sizeof(int) * ((N+1)*5));
  mx->mmx_mem = (void *) ReallocOrDie (mx->mmx_mem, sizeof(int) * ((N+1)*(M+2)));
  mx->imx_mem = (void *) ReallocOrDie (mx->imx_mem, sizeof(int) * ((N+1)*(M+2)));
  mx->dmx_mem = (void *) ReallocOrDie (mx->dmx_mem, sizeof(int) * ((N+1)*(M+2)));

  mx->xmx[0] = (int *) mx->xmx_mem;
  mx->mmx[0] = (int *) mx->mmx_mem;
  mx->imx[0] = (int *) mx->imx_mem;
  mx->dmx[0] = (int *) mx->dmx_mem;

  for (i = 1; i <= N; i++)
    {
      mx->xmx[i] = mx->xmx[0] + (i*5); 
      mx->mmx[i] = mx->mmx[0] + (i*(M+2));
      mx->imx[i] = mx->imx[0] + (i*(M+2));
      mx->dmx[i] = mx->dmx[0] + (i*(M+2));
    }

 DONE:
  if (xmx != NULL) *xmx = mx->xmx;
  if (mmx != NULL) *mmx = mx->mmx;
  if (imx != NULL) *imx = mx->imx;
  if (dmx != NULL) *dmx = mx->dmx;
}
#endif /*ALTIVEC*/

/* Function: AllocPlan7Matrix()
 * Date:     SRE, Tue Nov 19 07:14:47 2002 [St. Louis]
 *
 * Purpose:  Used to be the main allocator for dp matrices; we used to
 *           allocate, calculate, free. But this spent a lot of time
 *           in malloc(). Replaced with Create..() and Resize..() to
 *           allow matrix reuse in P7Viterbi(), the main alignment 
 *           engine. But matrices are alloc'ed by other alignment engines
 *           too, ones that are less frequently called and less 
 *           important to optimization of cpu performance. Instead of
 *           tracking changes through them, for now, provide
 *           an Alloc...() call with the same API that's just a wrapper.
 *
 * Args:     rows  - generally L+1, or 2; # of DP rows in seq dimension to alloc
 *           M     - size of model, in nodes
 *           xmx, mmx, imx, dmx 
 *                 - RETURN: ptrs to four mx components as a convenience
 *
 * Returns:  mx
 *           Caller free's w/ FreePlan7Matrix()
 */
struct dpmatrix_s *
AllocPlan7Matrix(int rows, int M, int ***xmx, int ***mmx, int ***imx, int ***dmx)
{
  struct dpmatrix_s *mx;
  mx = CreatePlan7Matrix(rows-1, M, 0, 0);
  if (xmx != NULL) *xmx = mx->xmx;
  if (mmx != NULL) *mmx = mx->mmx;
  if (imx != NULL) *imx = mx->imx;
  if (dmx != NULL) *dmx = mx->dmx;
  return mx;
}


/* Function: FreePlan7Matrix()
 * 
 * Purpose:  Free a dynamic programming matrix allocated by CreatePlan7Matrix().
 * 
 * Return:   (void)
 */
void
FreePlan7Matrix(struct dpmatrix_s *mx)
{
  free (mx->xmx_mem);
  free (mx->mmx_mem);
  free (mx->imx_mem);
  free (mx->dmx_mem);
  free (mx->xmx);
  free (mx->mmx);
  free (mx->imx);
  free (mx->dmx);
  free (mx);
}

/* Function: AllocShadowMatrix()
 * 
 * Purpose:  Allocate a dynamic programming traceback pointer matrix for 
 *           a Viterbi algorithm. 
 *           
 * Args:     rows  - number of rows to allocate; typically L+1
 *           M     - size of model
 *           xtb, mtb, itb, dtb 
 *                 - RETURN: ptrs to four mx components as a convenience
 *                 
 * Return:   mx
 *           mx is allocated here. Caller frees with FreeDPMatrix(mx).
 */

struct dpshadow_s *
AllocShadowMatrix(int rows, int M, char ***xtb, char ***mtb, char ***itb, char ***dtb)
{
  struct dpshadow_s *tb;
  int i;

  tb         = (struct dpshadow_s *) MallocOrDie (sizeof(struct dpshadow_s));
  tb->xtb    = (char **) MallocOrDie (sizeof(char *) * rows);
  tb->mtb    = (char **) MallocOrDie (sizeof(char *) * rows);
  tb->itb    = (char **) MallocOrDie (sizeof(char *) * rows);
  tb->dtb    = (char **) MallocOrDie (sizeof(char *) * rows);
  tb->esrc   = (int *)   MallocOrDie (sizeof(int)  * rows);
  tb->xtb[0] = (char *)  MallocOrDie (sizeof(char) * (rows*5));
  tb->mtb[0] = (char *)  MallocOrDie (sizeof(char) * (rows*(M+2)));
  tb->itb[0] = (char *)  MallocOrDie (sizeof(char) * (rows*(M+2)));
  tb->dtb[0] = (char *)  MallocOrDie (sizeof(char) * (rows*(M+2)));
  for (i = 1; i < rows; i++)
    {
      tb->xtb[i] = tb->xtb[0] + (i*5); 
      tb->mtb[i] = tb->mtb[0] + (i*(M+2));
      tb->itb[i] = tb->itb[0] + (i*(M+2));
      tb->dtb[i] = tb->dtb[0] + (i*(M+2));
    }

  if (xtb != NULL) *xtb = tb->xtb;
  if (mtb != NULL) *mtb = tb->mtb;
  if (itb != NULL) *itb = tb->itb;
  if (dtb != NULL) *dtb = tb->dtb;
  return tb;
}

/* Function: FreeShadowMatrix()
 * 
 * Purpose:  Free a dynamic programming matrix allocated by AllocShadowMatrix().
 * 
 * Return:   (void)
 */
void
FreeShadowMatrix(struct dpshadow_s *tb)
{
  free (tb->xtb[0]);
  free (tb->mtb[0]);
  free (tb->itb[0]);
  free (tb->dtb[0]);
  free (tb->esrc);
  free (tb->xtb);
  free (tb->mtb);
  free (tb->itb);
  free (tb->dtb);
  free (tb);
}

/* Function: P7ViterbiSize()
 * Date:     SRE, Fri Mar  6 15:13:20 1998 [St. Louis]
 *
 * Purpose:  Returns the ballpark predicted memory requirement for a 
 *           P7Viterbi() alignment, in MB. 
 *           
 *           Currently L must fit in an int (< 2 GB), but we have
 *           to deal with LM > 2 GB - e.g. watch out for overflow, do
 *           the whole calculation in floating point. Bug here detected
 *           in 2.1.1 by David Harper, Sanger Centre.
 *                                    
 * Args:     L  - length of sequence
 *           M  - length of HMM       
 *
 * Returns:  # of MB
 */
int
P7ViterbiSize(int L, int M)
{
  float Mbytes;

  /* We're excessively precise here, but it doesn't cost
   * us anything to be pedantic. The four terms are:
   *   1. the matrix structure itself;
   *   2. the O(NM) main matrix (this dominates!)
   *   3. ptrs into the rows of the matrix
   *   4. storage for 5 special states. (xmx)
   */
  Mbytes =  (float) sizeof(struct dpmatrix_s); 
  Mbytes += 3. * (float) (L+1) * (float) (M+2) * (float) sizeof(int); 
  Mbytes += 4. * (float) (L+1) * (float) sizeof(int *); 
  Mbytes += 5. * (float) (L+1) * (float) sizeof(int);
  Mbytes /= 1048576.;
  return (int) Mbytes;
}

/* Function: P7SmallViterbiSize()
 * Date:     SRE, Fri Mar  6 15:20:04 1998 [St. Louis]
 *
 * Purpose:  Returns the ballpark predicted memory requirement for
 *           a P7SmallViterbi() alignment, in MB. 
 *           
 *           P7SmallViterbi() is a wrapper, calling both P7ParsingViterbi()
 *           and P7WeeViterbi(). P7ParsingViterbi() typically dominates
 *           the memory requirement, so the value returned
 *           is the P7ParsingViterbi() number.
 *     
 *           We don't (yet) worry about overflow issues like we did with
 *           P7ViterbiSize(). We'll have many other 32-bit int issues in the
 *           code if we overflow here.
 *
 * Args:     L - length of sequence
 *           M - length of HMM   
 *               
 * Returns:  # of MB
 */
int
P7SmallViterbiSize(int L, int M)
{
  return ((2 * sizeof(struct dpmatrix_s) +
	   12 * (M+2) * sizeof(int)      + /* 2 matrices w/ 2 rows */ 
	   16 * sizeof(int *)            + /* ptrs into rows of matrix */
	   20 * sizeof(int)              + /* 5 special states     */
	   2  * (L+1) * sizeof(int))       /* traceback indices    */      
	  / 1000000);
}


/* Function: P7WeeViterbiSize()
 * Date:     SRE, Fri Mar  6 15:40:42 1998 [St. Louis]
 *
 * Purpose:  Returns the ballpark predicted memory requirement for
 *           a P7WeeViterbi() alignment, in MB. 
 *
 * Args:     L - length of sequence
 *           M - length of HMM   
 *               
 * Returns:  # of MB 
 */
int
P7WeeViterbiSize(int L, int M)
{
  return ((2 * sizeof(struct dpmatrix_s) +
	   12 * (M+2) * sizeof(int)      + /* 2 matrices w/ 2 rows */ 
	   16 * sizeof(int *)            + /* ptrs into rows of matrix */
	   20 * sizeof(int)              + /* 5 special states      */
	   2  * (L+2) * sizeof(int) +      /* stacks for starts/ends (overkill) */      
	   (L+2) * sizeof(int) +           /* k assignments to seq positions */
	   (L+2) * sizeof(char))           /* state assignments to seq pos */
	  / 1000000);
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
float
P7Forward(char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s **ret_mx)
{
  struct dpmatrix_s *mx;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;

  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 
  mx = AllocPlan7Matrix(L+1, hmm->M, &xmx, &mmx, &imx, &dmx);

  /* Initialization of the zero row.
   * Note that xmx[i][stN] = 0 by definition for all i,
   *    and xmx[i][stT] = xmx[i][stC], so neither stN nor stT need
   *    to be calculated in DP matrices.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 1; i <= L; i++)
    {
      mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;
      for (k = 1; k < hmm->M; k++)
	{
	  mmx[i][k]  = ILogsum(ILogsum(mmx[i-1][k-1] + hmm->tsc[TMM][k-1],
				     imx[i-1][k-1] + hmm->tsc[TIM][k-1]),
			      ILogsum(xmx[i-1][XMB] + hmm->bsc[k],
				     dmx[i-1][k-1] + hmm->tsc[TDM][k-1]));
	  mmx[i][k] += hmm->msc[(int) dsq[i]][k];

	  dmx[i][k]  = ILogsum(mmx[i][k-1] + hmm->tsc[TMD][k-1],
			      dmx[i][k-1] + hmm->tsc[TDD][k-1]);
	  imx[i][k]  = ILogsum(mmx[i-1][k] + hmm->tsc[TMI][k],
			      imx[i-1][k] + hmm->tsc[TII][k]);
	  imx[i][k] += hmm->isc[(int) dsq[i]][k];
	}
      mmx[i][hmm->M] = ILogsum(ILogsum(mmx[i-1][hmm->M-1] + hmm->tsc[TMM][hmm->M-1],
				   imx[i-1][hmm->M-1] + hmm->tsc[TIM][hmm->M-1]),
			       ILogsum(xmx[i-1][XMB] + hmm->bsc[hmm->M-1],
				   dmx[i-1][hmm->M-1] + hmm->tsc[TDM][hmm->M-1]));
      mmx[i][hmm->M] += hmm->msc[(int) dsq[i]][hmm->M];

      /* Now the special states.
       * remember, C and J emissions are zero score by definition
       */
      xmx[i][XMN] = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP];

      xmx[i][XME] = -INFTY;
      for (k = 1; k <= hmm->M; k++)
	xmx[i][XME] = ILogsum(xmx[i][XME], mmx[i][k] + hmm->esc[k]);

      xmx[i][XMJ] = ILogsum(xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP],
			   xmx[i][XME]   + hmm->xsc[XTE][LOOP]);

      xmx[i][XMB] = ILogsum(xmx[i][XMN] + hmm->xsc[XTN][MOVE],
			    xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]);

      xmx[i][XMC] = ILogsum(xmx[i-1][XMC] + hmm->xsc[XTC][LOOP],
			    xmx[i][XME] + hmm->xsc[XTE][MOVE]);
    }
			    
  sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

  if (ret_mx != NULL) *ret_mx = mx;
  else                FreePlan7Matrix(mx);

  return Scorify(sc);		/* the total Forward score. */
}

      
#ifdef SLOW
/* Function: P7Viterbi()
 * 
 * Purpose:  The Viterbi dynamic programming algorithm. 
 *           Identical to Forward() except that max's
 *           replace sum's. 
 *           
 *           This is the slower, more understandable version
 *           of P7Viterbi(). The default version in fast_algorithms.c
 *           is portably optimized and more difficult to understand;
 *           the ALTIVEC version in fast_algorithms.c is vectorized
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
float
P7Viterbi(char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s *mx, struct p7trace_s **ret_tr)
{
  struct p7trace_s  *tr;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;

  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 
  ResizePlan7Matrix(mx, L, hmm->M, &xmx, &mmx, &imx, &dmx);

  /* Initialization of the zero row.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 1; i <= L; i++) {
    mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;

    for (k = 1; k <= hmm->M; k++) {
				/* match state */
      mmx[i][k]  = -INFTY;
      if ((sc = mmx[i-1][k-1] + hmm->tsc[TMM][k-1]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = imx[i-1][k-1] + hmm->tsc[TIM][k-1]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = xmx[i-1][XMB] + hmm->bsc[k]) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = dmx[i-1][k-1] + hmm->tsc[TDM][k-1]) > mmx[i][k])
	mmx[i][k] = sc;
      if (hmm->msc[(int) dsq[i]][k] != -INFTY) mmx[i][k] += hmm->msc[(int) dsq[i]][k];
      else                                     mmx[i][k] = -INFTY;

				/* delete state */
      dmx[i][k] = -INFTY;
      if ((sc = mmx[i][k-1] + hmm->tsc[TMD][k-1]) > dmx[i][k])
	dmx[i][k] = sc;
      if ((sc = dmx[i][k-1] + hmm->tsc[TDD][k-1]) > dmx[i][k])
	dmx[i][k] = sc;

				/* insert state */
      if (k < hmm->M) {
	imx[i][k] = -INFTY;
	if ((sc = mmx[i-1][k] + hmm->tsc[TMI][k]) > imx[i][k])
	  imx[i][k] = sc;
	if ((sc = imx[i-1][k] + hmm->tsc[TII][k]) > imx[i][k])
	  imx[i][k] = sc;
	if (hmm->isc[(int)dsq[i]][k] != -INFTY) imx[i][k] += hmm->isc[(int) dsq[i]][k];
	else                                    imx[i][k] = -INFTY;   
      }
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
				/* N state */
    xmx[i][XMN] = -INFTY;
    if ((sc = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[i][XMN] = sc;

				/* E state */
    xmx[i][XME] = -INFTY;
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  mmx[i][k] + hmm->esc[k]) > xmx[i][XME])
	xmx[i][XME] = sc;
				/* J state */
    xmx[i][XMJ] = -INFTY;
    if ((sc = xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP]) > -INFTY)
      xmx[i][XMJ] = sc;
    if ((sc = xmx[i][XME]   + hmm->xsc[XTE][LOOP]) > xmx[i][XMJ])
      xmx[i][XMJ] = sc;

				/* B state */
    xmx[i][XMB] = -INFTY;
    if ((sc = xmx[i][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
      xmx[i][XMB] = sc;
    if ((sc = xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[i][XMB])
      xmx[i][XMB] = sc;

				/* C state */
    xmx[i][XMC] = -INFTY;
    if ((sc = xmx[i-1][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[i][XMC] = sc;
    if ((sc = xmx[i][XME] + hmm->xsc[XTE][MOVE]) > xmx[i][XMC])
      xmx[i][XMC] = sc;
  }
				/* T state (not stored) */
  sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

  if (ret_tr != NULL) {
    P7ViterbiTrace(hmm, dsq, L, mx, &tr);
    *ret_tr = tr;
  }

  return Scorify(sc);		/* the total Viterbi score. */
}
#endif /*SLOW*/

/* Function: P7ViterbiTrace()
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
void
P7ViterbiTrace(struct plan7_s *hmm, char *dsq, int N,
	       struct dpmatrix_s *mx, struct p7trace_s **ret_tr)
{
  struct p7trace_s *tr;
  int curralloc;		/* current allocated length of trace */
  int tpos;			/* position in trace */
  int i;			/* position in seq (1..N) */
  int k;			/* position in model (1..M) */
  int **xmx, **mmx, **imx, **dmx;
  int sc;			/* temp var for pre-emission score */

  /* Overallocate for the trace.
   * S-N-B- ... - E-C-T  : 6 states + N is minimum trace;
   * add N more as buffer.             
   */
  curralloc = N * 2 + 6; 
  P7AllocTrace(curralloc, &tr);

  xmx = mx->xmx;
  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;

  /* Initialization of trace
   * We do it back to front; ReverseTrace() is called later.
   */
  tr->statetype[0] = STT;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tr->statetype[1] = STC;
  tr->nodeidx[1]   = 0;
  tr->pos[1]       = 0;
  tpos = 2;
  i    = N;			/* current i (seq pos) we're trying to assign */

  /* Traceback
   */
  while (tr->statetype[tpos-1] != STS) {
    switch (tr->statetype[tpos-1]) {
    case STM:			/* M connects from i-1,k-1, or B */
      sc = mmx[i+1][k+1] - hmm->msc[(int) dsq[i+1]][k+1];
      if (sc <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (sc == xmx[i][XMB] + hmm->bsc[k+1])
	{
				/* Check for wing unfolding */
	  if (Prob2Score(hmm->begin[k+1], hmm->p1) + 1 * INTSCALE <= hmm->bsc[k+1])
	    while (k > 0)
	      {
		tr->statetype[tpos] = STD;
		tr->nodeidx[tpos]   = k--;
		tr->pos[tpos]       = 0;
		tpos++;
		if (tpos == curralloc) 
		  {				/* grow trace if necessary  */
		    curralloc += N;
		    P7ReallocTrace(tr, curralloc);
		  }
	      }

	  tr->statetype[tpos] = STB;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	}
      else if (sc == mmx[i][k] + hmm->tsc[TMM][k])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == imx[i][k] + hmm->tsc[TIM][k])
	{
	  tr->statetype[tpos] = STI;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == dmx[i][k] + hmm->tsc[TDM][k])
	{
	  tr->statetype[tpos] = STD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	}
      else
	Die("traceback failed");
      break;

    case STD:			/* D connects from M,D */
      if (dmx[i][k+1] <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (dmx[i][k+1] == mmx[i][k] + hmm->tsc[TMD][k])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (dmx[i][k+1] == dmx[i][k] + hmm->tsc[TDD][k]) 
	{
	  tr->statetype[tpos] = STD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	}
      else Die("traceback failed");
      break;

    case STI:			/* I connects from M,I */
      sc = imx[i+1][k] - hmm->isc[(int) dsq[i+1]][k];
      if (sc <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (sc == mmx[i][k] + hmm->tsc[TMI][k])
	{
	  tr->statetype[tpos] = STM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == imx[i][k] + hmm->tsc[TII][k])
	{
	  tr->statetype[tpos] = STI;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = i--;
	}
      else Die("traceback failed");
      break;

    case STN:			/* N connects from S, N */
      if (i == 0 && xmx[i][XMN] == 0)
	{
	  tr->statetype[tpos] = STS;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	}
      else if (i > 0 && xmx[i+1][XMN] == xmx[i][XMN] + hmm->xsc[XTN][LOOP])
	{
	  tr->statetype[tpos] = STN;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;    /* note convention adherence:  */
	  tr->pos[tpos-1]     = i--;  /* first N doesn't emit        */
	}
      else Die("traceback failed");
      break;

    case STB:			/* B connects from N, J */
      if (xmx[i][XMB] <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (xmx[i][XMB] == xmx[i][XMN] + hmm->xsc[XTN][MOVE])
	{
	  tr->statetype[tpos] = STN;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	}
      else if (xmx[i][XMB] == xmx[i][XMJ] + hmm->xsc[XTJ][MOVE])
	{
	  tr->statetype[tpos] = STJ;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	}

      else Die("traceback failed");
      break;

    case STE:			/* E connects from any M state. k set here */
      if (xmx[i][XME] <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      for (k = hmm->M; k >= 1; k--)
	if (xmx[i][XME] == mmx[i][k] + hmm->esc[k])
	  {
				/* check for wing unfolding */
	    if (Prob2Score(hmm->end[k], 1.) + 1*INTSCALE <=  hmm->esc[k])
	      {
		int dk;		/* need a tmp k while moving thru delete wing */
		for (dk = hmm->M; dk > k; dk--)
		  {
		    tr->statetype[tpos] = STD;
		    tr->nodeidx[tpos]   = dk;
		    tr->pos[tpos]       = 0;
		    tpos++;
		    if (tpos == curralloc) 
		      {				/* grow trace if necessary  */
			curralloc += N;
			P7ReallocTrace(tr, curralloc);
		      }
		  }
	      }

	    tr->statetype[tpos] = STM;
	    tr->nodeidx[tpos]   = k--;
	    tr->pos[tpos]       = i--;
	    break;
	  }
      if (k < 0) Die("traceback failed");
      break;

    case STC:			/* C comes from C, E */
      if (xmx[i][XMC] <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (xmx[i][XMC] == xmx[i-1][XMC] + hmm->xsc[XTC][LOOP])
	{
	  tr->statetype[tpos] = STC;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;    /* note convention adherence: */
	  tr->pos[tpos-1]     = i--;  /* first C doesn't emit       */
	}
      else if (xmx[i][XMC] == xmx[i][XME] + hmm->xsc[XTE][MOVE])
	{
	  tr->statetype[tpos] = STE;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0; /* E is a nonemitter */
	}
      
      else Die("Traceback failed.");
      break;

    case STJ:			/* J connects from E, J */
      if (xmx[i][XMJ] <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (xmx[i][XMJ] == xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP])
	{
	  tr->statetype[tpos] = STJ;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;    /* note convention adherence: */
	  tr->pos[tpos-1]     = i--;  /* first J doesn't emit       */
	}
      else if (xmx[i][XMJ] == xmx[i][XME] + hmm->xsc[XTE][LOOP])
	{
	  tr->statetype[tpos] = STE;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0; /* E is a nonemitter */
	}

      else Die("Traceback failed.");
      break;

    default:
      Die("traceback failed");

    } /* end switch over statetype[tpos-1] */
    
    tpos++;
    if (tpos == curralloc) 
      {				/* grow trace if necessary  */
	curralloc += N;
	P7ReallocTrace(tr, curralloc);
      }

  } /* end traceback, at S state; tpos == tlen now */
  tr->tlen = tpos;
  P7ReverseTrace(tr);
  *ret_tr = tr;
}


/* Function: P7SmallViterbi()
 * Date:     SRE, Fri Mar  6 15:29:41 1998 [St. Louis]
 *
 * Purpose:  Wrapper function, for linear memory alignment
 *           with same arguments as P7Viterbi(). 
 *           
 *           Calls P7ParsingViterbi to break the sequence
 *           into fragments. Then, based on size of fragments,
 *           calls either P7Viterbi() or P7WeeViterbi() to 
 *           get traces for them. Finally, assembles all these
 *           traces together to produce an overall optimal
 *           trace for the sequence.
 *           
 *           If the trace isn't needed for some reason,
 *           all we do is call P7ParsingViterbi.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           mx     - DP matrix, growable
 *           ret_tr - RETURN: traceback; pass NULL if it's not wanted
 *
 * Returns:  Score of optimal alignment in bits.
 */
float
P7SmallViterbi(char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s *mx, struct p7trace_s **ret_tr)
{
  struct p7trace_s *ctr;        /* collapsed trace of optimal parse */
  struct p7trace_s *tr;         /* full trace of optimal alignment */
  struct p7trace_s **tarr;      /* trace array */   
  int   ndom;			/* number of subsequences */
  int   i;			/* counter over domains   */
  int   pos;			/* position in sequence */
  int   tpos;			/* position in trace */
  int   tlen;			/* length of full trace   */
  int   sqlen;			/* length of a subsequence */
  int   totlen;                 /* length of L matched by model (as opposed to N/C/J) */
  float sc;			/* score of optimal alignment */
  int   t2;			/* position in a subtrace */
  
  /* Step 1. Call P7ParsingViterbi to calculate an optimal parse
   *         of the sequence into single-hit subsequences; this parse
   *         is returned in a "collapsed" trace
   */
  sc = P7ParsingViterbi(dsq, L, hmm, &ctr);

  /* If we don't want full trace, we're done;
   * also, if parsing viterbi returned a NULL trace we're done. */
  if (ctr == NULL || ret_tr == NULL)
    {
      P7FreeTrace(ctr);
      return sc;
    }
  
  /* Step 2. Call either P7Viterbi or P7WeeViterbi on each subsequence
   *         to recover a full traceback of each, collecting them in 
   *         an array. 
   */
  ndom = ctr->tlen/2 - 1;
  tarr = MallocOrDie(sizeof(struct p7trace_s *) * ndom);
  tlen = totlen = 0;
  for (i = 0; i < ndom; i++)
    {
      sqlen = ctr->pos[i*2+2] - ctr->pos[i*2+1];   /* length of subseq */

      if (P7ViterbiSize(sqlen, hmm->M) > RAMLIMIT)
	P7WeeViterbi(dsq + ctr->pos[i*2+1], sqlen, hmm, &(tarr[i]));
      else
	P7Viterbi(dsq + ctr->pos[i*2+1], sqlen, hmm, mx, &(tarr[i]));

      tlen  += tarr[i]->tlen - 4; /* not counting S->N,...,C->T */
      totlen += sqlen;
    }

  /* Step 3. Compose the subtraces into one big final trace.
   *         This is wasteful because we're going to TraceDecompose()
   *         it again in both hmmsearch and hmmpfam to look at
   *         individual domains; but we do it anyway so the P7SmallViterbi
   *         interface looks exactly like the P7Viterbi interface. Maybe
   *         long traces shouldn't include all the N/J/C states anyway,
   *         since they're unambiguously implied.
   */

  /* Calculate total trace len and alloc;
   * nonemitting SNCT + nonemitting J's + emitting NJC
   */
  tlen += 4 + (ndom-1) + (L-totlen);
  P7AllocTrace(tlen, &tr);
  tr->tlen = tlen;

  /* Add N-terminal trace framework
   */
  tr->statetype[0] = STS;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tr->statetype[1] = STN;
  tr->nodeidx[1]   = 0;
  tr->pos[1]       = 0;
  tpos = 2;
				/* add implied N's */
  for (pos = 1; pos <= ctr->pos[1]; pos++)
    {
      tr->statetype[tpos] = STN;
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = pos;
      tpos++;
    }

  /* Add each subseq trace in, with its appropriate 
   * sequence offset derived from the collapsed trace
   */
  for (i = 0; i < ndom; i++)
    {				/* skip SN, CT framework at ends */
      for (t2 = 2; t2 < tarr[i]->tlen-2; t2++)
	{
	  tr->statetype[tpos] = tarr[i]->statetype[t2];
	  tr->nodeidx[tpos]   = tarr[i]->nodeidx[t2];
	  if (tarr[i]->pos[t2] > 0) 
	    tr->pos[tpos]       = tarr[i]->pos[t2] + ctr->pos[i*2+1];
	  else
	    tr->pos[tpos]       = 0;
	  tpos++;
	}
				/* add nonemitting J or C */
      tr->statetype[tpos] = (i == ndom-1) ? STC : STJ;
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = 0;
      tpos++;
				/* add implied emitting J's */
      if (i != ndom-1)
	for (pos = ctr->pos[i*2+2]+1; pos <= ctr->pos[(i+1)*2+1]; pos++)
	  {
	    tr->statetype[tpos] = STJ;
	    tr->nodeidx[tpos]   = 0;
	    tr->pos[tpos]       = pos;
	    tpos++; 
	  }
    }

				/* add implied C's */
  for (pos = ctr->pos[ndom*2]+1; pos <= L; pos++)
    {
      tr->statetype[tpos] = STC;
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = pos;
      tpos++;
    }
				/* add terminal T */
  tr->statetype[tpos] = STT;
  tr->nodeidx[tpos]   = 0;
  tr->pos[tpos]       = 0;
  tpos++;
	    
  for (i = 0; i < ndom; i++) P7FreeTrace(tarr[i]);
  free(tarr);
  P7FreeTrace(ctr);

  *ret_tr = tr;
  return sc;
}




/* Function: P7ParsingViterbi()
 * Date:     SRE, Wed Mar  4 14:07:31 1998 [St. Louis]
 *
 * Purpose:  The "hmmfs" linear-memory algorithm for finding
 *           the optimal alignment of a very long sequence to
 *           a looping, multihit (e.g. Plan7) model, parsing it into
 *           a series of nonoverlapping subsequences that match  
 *           the model once. Other algorithms (e.g. P7Viterbi()
 *           or P7WeeViterbi()) are applied subsequently to
 *           these subsequences to recover complete alignments.
 *           
 *           The hmmfs algorithm appears briefly in [Durbin98],
 *           but is otherwise unpublished.
 *           
 *           The traceback structure returned is special: a
 *           "collapsed" trace S->B->E->...->B->E->T, where
 *           stateidx is unused and pos is used to indicate the
 *           position of B and E in the sequence. The matched
 *           subsequence is B_pos+1...E_pos. The number of
 *           matches in the trace is (tlen/2)-1.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model (log odds scores ready)
 *           ret_tr - RETURN: a collapsed traceback.      
 *
 * Returns:  Score of the optimal Viterbi alignment, in bits.
 */
float
P7ParsingViterbi(char *dsq, int L, struct plan7_s *hmm, struct p7trace_s **ret_tr)
{
  struct dpmatrix_s *mx;        /* two rows of score matrix */
  struct dpmatrix_s *tmx;       /* two rows of misused score matrix: traceback ptrs */
  struct p7trace_s  *tr;        /* RETURN: collapsed traceback */
  int  **xmx, **mmx, **dmx, **imx; /* convenience ptrs to score matrix */  
  int  **xtr, **mtr, **dtr, **itr; /* convenience ptrs to traceback pointers */
  int   *btr, *etr;             /* O(L) trace ptrs for B, E state pts in seq */    
  int    sc;			/* integer score of optimal alignment  */
  int    i,k,tpos;		/* index for seq, model, trace position */
  int    cur, prv;		/* indices for rolling dp matrix */
  int    curralloc;		/* size of allocation for tr */


  /* Alloc a DP matrix and traceback pointers, two rows each, O(M).
   * Alloc two O(L) arrays to trace back through the sequence thru B and E.
   */
  mx  = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  tmx = AllocPlan7Matrix(2, hmm->M, &xtr, &mtr, &itr, &dtr);
  btr = MallocOrDie(sizeof(int) * (L+1));
  etr = MallocOrDie(sizeof(int) * (L+1));

  /* Initialization of the zero row.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  btr[0]      = 0;
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  etr[0]      = -1; 
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Recursion. Done as a pull. Rolling index trick. Trace ptr propagation trick.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   *    
   * Notes on traceback pointer propagation.
   *  - In the path B->E, we propagate the i that B was aligned to in the optimal
   *    alignment, via mtr, dtr, and itr. 
   *  - When we reach an E, we record the i of the B it started from in etr.
   *  - In a looping path E->J...->B or terminal path E->C...->T, we propagate
   *    the i that E was aligned to in the optimal alignment via xtr[][XMC]
   *    and xtr[][XMJ].
   *  - When we enter B, we record the i of the best previous E, or 0 if there
   *    isn't one, in btr.
   */
  for (i = 1; i <= L; i++) {
    cur = i % 2;
    prv = !cur;
    
    mmx[cur][0] = imx[cur][0] = dmx[cur][0] = -INFTY;

    for (k = 1; k <= hmm->M; k++) {
				/* match state */
      mmx[cur][k] = -INFTY;
      if ((sc = mmx[prv][k-1] + hmm->tsc[TMM][k-1]) > -INFTY)
	{ mmx[cur][k] = sc; mtr[cur][k] = mtr[prv][k-1]; }
      if ((sc = imx[prv][k-1] + hmm->tsc[TIM][k-1]) > mmx[cur][k])
	{ mmx[cur][k] = sc; mtr[cur][k] = itr[prv][k-1]; }
      if ((sc = xmx[prv][XMB] + hmm->bsc[k]) > mmx[cur][k])
	{ mmx[cur][k] = sc; mtr[cur][k] = i-1; }
      if ((sc = dmx[prv][k-1] + hmm->tsc[TDM][k-1]) > mmx[cur][k])
	{ mmx[cur][k] = sc; mtr[cur][k] = dtr[prv][k-1]; }
      if (hmm->msc[(int) dsq[i]][k] != -INFTY)
	mmx[cur][k] += hmm->msc[(int) dsq[i]][k];
      else
	mmx[cur][k] = -INFTY;

				/* delete state */
      dmx[cur][k] = -INFTY;
      if ((sc = mmx[cur][k-1] + hmm->tsc[TMD][k-1]) > -INFTY)
	{ dmx[cur][k] = sc; dtr[cur][k] = mtr[cur][k-1]; }
      if ((sc = dmx[cur][k-1] + hmm->tsc[TDD][k-1]) > dmx[cur][k])
	{ dmx[cur][k] = sc; dtr[cur][k] = dtr[cur][k-1]; }

				/* insert state */
      if (k < hmm->M) {
	imx[cur][k] = -INFTY;
	if ((sc = mmx[prv][k] + hmm->tsc[TMI][k]) > -INFTY)
	  { imx[cur][k] = sc; itr[cur][k] = mtr[prv][k]; }
	if ((sc = imx[prv][k] + hmm->tsc[TII][k]) > imx[cur][k])
	  { imx[cur][k] = sc; itr[cur][k] = itr[prv][k]; }
	if (hmm->isc[(int) dsq[i]][k] != -INFTY)
	  imx[cur][k] += hmm->isc[(int) dsq[i]][k];
	else
	  imx[cur][k] = -INFTY;
      }
    }

    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
				/* N state */
    xmx[cur][XMN] = -INFTY;
    if ((sc = xmx[prv][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[cur][XMN] = sc;
				/* E state */
    xmx[cur][XME] = -INFTY;
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  mmx[cur][k] + hmm->esc[k]) > xmx[cur][XME])
	{ xmx[cur][XME] = sc; etr[i] = mtr[cur][k]; }
				/* J state */
    xmx[cur][XMJ] = -INFTY;
    if ((sc = xmx[prv][XMJ] + hmm->xsc[XTJ][LOOP]) > -INFTY)
      { xmx[cur][XMJ] = sc; xtr[cur][XMJ] = xtr[prv][XMJ]; }
    if ((sc = xmx[cur][XME]   + hmm->xsc[XTE][LOOP]) > xmx[cur][XMJ])
      { xmx[cur][XMJ] = sc; xtr[cur][XMJ] = i; }
				/* B state */
    xmx[cur][XMB] = -INFTY;
    if ((sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
      { xmx[cur][XMB] = sc; btr[i] = 0; }
    if ((sc = xmx[cur][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[cur][XMB])
      { xmx[cur][XMB] = sc; btr[i] = xtr[cur][XMJ]; }
				/* C state */
    xmx[cur][XMC] = -INFTY;
    if ((sc = xmx[prv][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      { xmx[cur][XMC] = sc; xtr[cur][XMC] = xtr[prv][XMC]; }
    if ((sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > xmx[cur][XMC])
      { xmx[cur][XMC] = sc; xtr[cur][XMC] = i; }
  }
				/* T state (not stored) */
  sc = xmx[cur][XMC] + hmm->xsc[XTC][MOVE];

  /*****************************************************************
   * Collapsed traceback stage. 
   * xtr[L%2][XMC] contains the position j of the previous E
   * etr[j]        contains the position i of the previous B
   * btr[i]        contains the position j of the previous E, or 0
   * continue until btr[i] = 0.
   *****************************************************************/

  curralloc = 2;		/* minimum: no hits */
  P7AllocTrace(curralloc, &tr);

  /* Init of collapsed trace. Back to front; we ReverseTrace() later.
   */
  tpos = 0;
  tr->statetype[tpos] = STT;
  tr->pos[tpos]       = 0;
  i                   = xtr[L%2][XMC];
  while (i > 0)
    {
      curralloc += 2;
      P7ReallocTrace(tr, curralloc);

      tpos++;
      tr->statetype[tpos] = STE;
      tr->pos[tpos] = i;
      i = etr[i];

      tpos++;
      tr->statetype[tpos] = STB;
      tr->pos[tpos] = i;
      i = btr[i];
    }

  tpos++;
  tr->statetype[tpos] = STS;
  tr->pos[tpos]       = 0;
  tr->tlen = tpos + 1;
  P7ReverseTrace(tr);
  
  FreePlan7Matrix(mx);
  FreePlan7Matrix(tmx);
  free(btr);
  free(etr);

  *ret_tr = tr;
  return Scorify(sc);
}

/* Function: P7WeeViterbi()
 * Date:     SRE, Wed Mar  4 08:24:04 1998 [St. Louis]
 *
 * Purpose:  Hirschberg/Myers/Miller linear memory alignment.
 *           See [Hirschberg75,MyM-88a] for the idea of the algorithm.
 *           Adapted to HMM implementation. 
 *           
 *           Requires that you /know/ that there's only
 *           one hit to the model in the sequence: either
 *           because you're forcing single-hit, or you've
 *           previously called P7ParsingViterbi to parse
 *           the sequence into single-hit segments. The reason
 *           for this is that a cyclic model (a la Plan7)
 *           defeats the nice divide and conquer trick.
 *           (I think some trickery with propagated trace pointers
 *           could get around this but haven't explored it.)
 *           This is implemented by ignoring transitions
 *           to/from J state.
 *
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_tr - RETURN: traceback.
 *
 * Returns:  Score of the optimal Viterbi alignment.
 */
float
P7WeeViterbi(char *dsq, int L, struct plan7_s *hmm, struct p7trace_s **ret_tr)
{
  struct p7trace_s *tr;         /* RETURN: traceback */
  int          *kassign;        /* 0..L+1, alignment of seq positions to model nodes */
  char         *tassign;        /* 0..L+1, alignment of seq positions to state types */
  int          *endlist;        /* stack of end points on sequence to work on */
  int          *startlist;      /* stack of start points on sequence to work on */
  int          lpos;            /* position in endlist, startlist */
  int          k1, k2, k3;	/* start, mid, end in model      */
  char         t1, t2, t3;	/* start, mid, end in state type */
  int          s1, s2, s3;	/* start, mid, end in sequence   */
  float        sc;		/* score of segment optimal alignment */
  float        ret_sc;		/* optimal score over complete seq */
  int          tlen;		/* length needed for trace */
  int          i, k, tpos;	/* index in sequence, model, trace */


  /* Initialize.
   */
  kassign   = MallocOrDie (sizeof(int) * (L+1));
  tassign   = MallocOrDie (sizeof(char)* (L+1));
  endlist   = MallocOrDie (sizeof(int) * (L+1));
  startlist = MallocOrDie (sizeof(int) * (L+1));

  lpos = 0;
  startlist[lpos] = 1;
  endlist[lpos]   = L;
  kassign[1]      = 1;
  kassign[L]      = hmm->M;
  tassign[1]      = STS;	/* temporary boundary condition! will become N or M */
  tassign[L]      = STT;	/* temporary boundary condition! will become M or C */

  /* Recursive divide-and-conquer alignment.
   */
  while (lpos >= 0)
    {
				/* Pop a segment off the stack */
      s1 = startlist[lpos];
      k1 = kassign[s1];
      t1 = tassign[s1];
      s3 = endlist[lpos];
      k3 = kassign[s3];
      t3 = tassign[s3];
      lpos--;
				/* find optimal midpoint of segment */
      sc = get_wee_midpt(hmm, dsq, L, k1, t1, s1, k3, t3, s3, &k2, &t2, &s2);
      kassign[s2] = k2;
      tassign[s2] = t2;
                               /* score is valid on first pass */
      if (t1 == STS && t3 == STT) ret_sc = sc;

				/* push N-terminal segment on stack */
      if (t2 != STN && (s2 - s1 > 1 || (s2 - s1 == 1 && t1 == STS)))
	{
	  lpos++;
	  startlist[lpos] = s1;
	  endlist[lpos]   = s2;
	}
				/* push C-terminal segment on stack */
      if (t2 != STC && (s3 - s2 > 1 || (s3 - s2 == 1 && t3 == STT)))
	{
          lpos++;
          startlist[lpos] = s2;
          endlist[lpos]   = s3;
	}

      if (t2 == STN)
	{			/* if we see STN midpoint, we know the whole N-term is STN */
	  for (; s2 >= s1; s2--) {
	    kassign[s2] = 1;
	    tassign[s2] = STN;
	  }
	}
      if (t2 == STC)
	{			/* if we see STC midpoint, we know whole C-term is STC */
	  for (; s2 <= s3; s2++) {
	    kassign[s2] = hmm->M;
	    tassign[s2] = STC;
	  }
	}
    }

  /*****************************************************************
   * Construct a traceback structure from kassign/tassign by interpolating
   * necessary states. 
   * Trace allocation is as follows. We clearly need L emitting states.
   * We also need nonemitting states as follows:
   * STS,STN,STB,STE,STC,STT = 6
   * STD: count k2-k1-1 in kassign M->M's
   * Also, count N->M's and M->C's (potential wing unfoldings)...
   *   ...and be careful to check wing unfoldings when there aren't
   *      any emitting N or C flanks! (bugfix, 2.1.1b)
   *****************************************************************/ 

  tlen = L + 6;
  for (i = 1; i < L; i++)
    {
      if (tassign[i] == STM && tassign[i+1] == STM)
	tlen += kassign[i+1] - kassign[i] - 1;
      if (tassign[i] == STN && tassign[i+1] == STM)
	tlen += kassign[i+1] - 1;
      if (tassign[i] == STM && tassign[i+1] == STC)
	tlen += hmm->M - kassign[i];
    }
  if (tassign[1] == STM) tlen += kassign[1] - 1; 
  if (tassign[L] == STM) tlen += hmm->M - kassign[L];
  P7AllocTrace(tlen, &tr);

  tr->statetype[0] = STS;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tr->statetype[1] = STN;
  tr->nodeidx[1]   = 0;
  tr->pos[1]       = 0;
  tpos = 2;
  
  for (i = 1; i <= L; i++)
    {
      switch(tassign[i]) {
      case STM:
				/* check for first match state */
	if (tr->statetype[tpos-1] == STN) {
	  tr->statetype[tpos] = STB;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0;
	  tpos++;
				/* check for wing unfolding */
	  if (Prob2Score(hmm->begin[kassign[i]], hmm->p1) + INTSCALE <= hmm->bsc[kassign[i]])
	    for (k = 1; k < kassign[i]; k++) {
	      tr->statetype[tpos] = STD;
	      tr->nodeidx[tpos]   = k;
	      tr->pos[tpos]       = 0;
	      tpos++;
	    }
	}
				/* do the match state itself */
	tr->statetype[tpos] = STM;
	tr->nodeidx[tpos]   = kassign[i];
	tr->pos[tpos]       = i;
	tpos++;
				/* do any deletes necessary 'til next match */
	if (i < L && tassign[i+1] == STM && kassign[i+1] - kassign[i] > 1)
	  for (k = kassign[i] + 1; k < kassign[i+1]; k++)
	    {
	      tr->statetype[tpos] = STD;
	      tr->nodeidx[tpos]   = k;
	      tr->pos[tpos]       = 0;
	      tpos++;
	    }
				/* check for last match state */
	if (i == L || tassign[i+1] == STC) {
				/* check for wing unfolding */
	  if (Prob2Score(hmm->end[kassign[i-1]], 1.) + INTSCALE <=  hmm->esc[kassign[i-1]])
	    for (k = kassign[i]+1; k <= hmm->M; k++)
	      {
		tr->statetype[tpos] = STD;
		tr->nodeidx[tpos]   = k;
		tr->pos[tpos]       = 0;
		tpos++;
	      }
				/* add on the end state */
	  tr->statetype[tpos] = STE;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0; 
	  tpos++;
				/* and a nonemitting C state */
	  tr->statetype[tpos] = STC;
	  tr->nodeidx[tpos]   = 0;
	  tr->pos[tpos]       = 0; 
	  tpos++;
	}
	break;
	
      case STI:
	tr->statetype[tpos] = STI;
	tr->nodeidx[tpos]   = kassign[i];
	tr->pos[tpos]       = i;
	tpos++;
	break;

      case STN:
	tr->statetype[tpos] = STN;
	tr->nodeidx[tpos]   = 0;
	tr->pos[tpos]       = i;
	tpos++;
	break;

      case STC:
	tr->statetype[tpos] = STC;
	tr->nodeidx[tpos]   = 0;
	tr->pos[tpos]       = i; 
	tpos++;
	break;

      default: Die("Bogus state %s", Statetype(tassign[i]));
      }
    }
				/* terminate the trace */
  tr->statetype[tpos] = STT;
  tr->nodeidx[tpos]   = 0;
  tr->pos[tpos]       = 0; 
  tr->tlen = tpos+1;

  *ret_tr = tr;

  free(kassign);
  free(tassign);
  free(startlist);
  free(endlist);
  return ret_sc;
}


/* Function: Plan7ESTViterbi()
 * 
 * Purpose:  Frameshift-tolerant alignment of protein model to cDNA EST.
 *           
 * 
 */
float
Plan7ESTViterbi(char *dsq, int L, struct plan7_s *hmm, struct dpmatrix_s **ret_mx)
{
  struct dpmatrix_s *mx;
  int **xmx;
  int **mmx;
  int **imx;
  int **dmx;
  int   i,k;
  int   sc;
  int   codon;
  
  /* Allocate a DP matrix with 0..L rows, 0..M+1 columns.
   */ 
  mx = AllocPlan7Matrix(L+1, hmm->M, &xmx, &mmx, &imx, &dmx);

  /* Initialization of the zero row (DNA sequence of length 0)
   * Note that xmx[i][stN] = 0 by definition for all i,
   *    and xmx[i][stT] = xmx[i][stC], so neither stN nor stT need
   *    to be calculated in DP matrices.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */

  /* Initialization of the first row (DNA sequence of length 1);
   * only N state can make this nucleotide.
   */
  xmx[1][XMN] = xmx[0][XMN] + hmm->xsc[XTN][LOOP];
  xmx[1][XMB] = xmx[1][XMN] + hmm->xsc[XTN][MOVE]; 
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need 2 nt to get here */
  for (k = 0; k <= hmm->M; k++)
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need 2 nt to get into model */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 2; i <= L; i++) {
    mmx[i][0] = imx[i][0] = dmx[i][0] = -INFTY;

				/* crude calculation of lookup value for codon */
    if (i > 2) {
      if (dsq[i-2] < 4 && dsq[i-1] < 4 && dsq[i] < 4)
	codon = dsq[i-2] * 16 + dsq[i-1] * 4 + dsq[i];
      else
	codon = 64;		/* ambiguous codon; punt */
    }

    for (k = 1; k <= hmm->M; k++) {
				/* match state */
      if (i > 2) {
	mmx[i][k]  = mmx[i-3][k-1] + hmm->tsc[TMM][k-1];
	if ((sc = imx[i-3][k-1] + hmm->tsc[TIM][k-1]) > mmx[i][k])
	  mmx[i][k] = sc;
	if ((sc = xmx[i-3][XMB] + hmm->bsc[k]) > mmx[i][k])
	  mmx[i][k] = sc;
	if ((sc = dmx[i-3][k-1] + hmm->tsc[TDM][k-1]) > mmx[i][k])
	  mmx[i][k] = sc;
	mmx[i][k] += hmm->dnam[codon][k];
      }
				/* -1 frameshifts into match state */
      if ((sc = mmx[i-2][k-1] + hmm->tsc[TMM][k-1] + hmm->dna2) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = imx[i-2][k-1] + hmm->tsc[TIM][k-1] + hmm->dna2) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = xmx[i-2][XMB] + hmm->bsc[k] + hmm->dna2) > mmx[i][k])
	mmx[i][k] = sc;
      if ((sc = dmx[i-2][k-1] + hmm->tsc[TDM][k-1] + hmm->dna2) > mmx[i][k])
	mmx[i][k] = sc;
      
				/* +1 frameshifts into match state */
      if (i > 3) {
	if ((sc = mmx[i-4][k-1] + hmm->tsc[TMM][k-1] + hmm->dna4) > mmx[i][k])
	  mmx[i][k] = sc;
	if ((sc = imx[i-4][k-1] + hmm->tsc[TIM][k-1] + hmm->dna4) > mmx[i][k])
	  mmx[i][k] = sc;
	if ((sc = xmx[i-4][XMB] + hmm->bsc[k] + hmm->dna4) > mmx[i][k])
	  mmx[i][k] = sc;
	if ((sc = dmx[i-4][k-1] + hmm->tsc[TDM][k-1] + hmm->dna4) > mmx[i][k])
	  mmx[i][k] = sc;
      }
      				/* delete state */
      dmx[i][k]  = mmx[i][k-1] + hmm->tsc[TMD][k-1];
      if ((sc = dmx[i][k-1] + hmm->tsc[TDD][k-1]) > dmx[i][k])
	dmx[i][k] = sc;

				/* insert state */
      if (i > 2) {
	imx[i][k] = mmx[i-3][k] + hmm->tsc[TMI][k];
	if ((sc = imx[i-3][k] + hmm->tsc[TII][k]) > imx[i][k])
	  imx[i][k] = sc;
	imx[i][k] += hmm->dnai[codon][k];
      }

				/* -1 frameshifts into insert state */
      if ((sc = mmx[i-2][k] + hmm->tsc[TMI][k] + hmm->dna2) > imx[i][k])
	imx[i][k] = sc;
      if ((sc = imx[i-2][k] + hmm->tsc[TII][k] + hmm->dna2) > imx[i][k])
	imx[i][k] = sc;

				/* +1 frameshifts into insert state */
      if (i > 4) {
	if ((sc = mmx[i-4][k] + hmm->tsc[TMI][k] + hmm->dna4) > imx[i][k])
	  imx[i][k] = sc;
	if ((sc = imx[i-4][k] + hmm->tsc[TII][k] + hmm->dna4) > imx[i][k])
	  imx[i][k] = sc;
      }
    }
    /* Now the special states. Order is important here.
     * remember, C and J emissions are zero score by definition,
     */
				/* N state: +1 nucleotide */
    xmx[i][XMN] = xmx[i-1][XMN] + hmm->xsc[XTN][LOOP];
                                /* E state: collect from M's, and last D  */
    xmx[i][XME] = dmx[i][hmm->M];    /* transition prob from last D = 1.0 */
    for (k = 1; k <= hmm->M; k++)
      if ((sc =  mmx[i][k] + hmm->esc[k]) > xmx[i][XME])
        xmx[i][XME] = sc;
                                /* J state: +1 nucleotide */
    xmx[i][XMJ] = xmx[i-1][XMJ] + hmm->xsc[XTJ][LOOP];
    if ((sc = xmx[i][XME]   + hmm->xsc[XTE][LOOP]) > xmx[i][XMJ])
      xmx[i][XMJ] = sc;
                                /* B state: collect from N,J */
    xmx[i][XMB] = xmx[i][XMN] + hmm->xsc[XTN][MOVE];
    if ((sc = xmx[i][XMJ] + hmm->xsc[XTJ][MOVE]) > xmx[i][XMB])
      xmx[i][XMB] = sc;
				/* C state: +1 nucleotide */
    xmx[i][XMC] = xmx[i-1][XMC] + hmm->xsc[XTC][LOOP];
    if ((sc = xmx[i][XME] + hmm->xsc[XTE][MOVE]) > xmx[i][XMC])
      xmx[i][XMC] = sc;
  }

  sc = xmx[L][XMC] + hmm->xsc[XTC][MOVE];

  if (ret_mx != NULL) *ret_mx = mx;
  else                FreePlan7Matrix(mx);

  return Scorify(sc);            /* the total Viterbi score. */
}



/* Function: get_wee_midpt()
 * Date:     SRE, Wed Mar  4 08:27:11 1998 [St. Louis]
 *
 * Purpose:  The heart of the divide and conquer algorithm
 *           for P7WeeViterbi(). This function is called
 *           recursively to find successive optimal midpoints
 *           in the alignment matrix. See P7WeeViterbi() for
 *           further comments on the assumptions of this algorithm.
 *           
 * Args:     hmm   - the model, set up for integer scores
 *           dsq   - the sequence, digitized
 *           L     - length of the sequence
 *           k1    - model node to start with, 1..M
 *           t1    - state type to start with, STM | STI | STN | STC; STS to start
 *           s1    - sequence position to start with, 1..L; 1 to start
 *           k3    - model node to end with, 1..M
 *           t3    - state type to end with, STM | STI | STN | STC; STT to start
 *           s3    - sequence position to end with, 1..L; L to start
 *          ret_k2 - RETURN: optimal midpoint, node position in model
 *          ret_t2 - RETURN: optimal midpoint, state type 
 *          ret_s2 - RETURN: optimal midpoint, sequence position
 *
 * Returns: score of optimal alignment, in bits. 
 */
static float
get_wee_midpt(struct plan7_s *hmm, char *dsq, int L, 
	      int k1, char t1, int s1,
	      int k3, char t3, int s3,
	      int *ret_k2, char *ret_t2, int *ret_s2)
{
  struct dpmatrix_s *fwd;
  struct dpmatrix_s *bck;
  int        **xmx;             /* convenience ptr into special states */
  int        **mmx;             /* convenience ptr into match states   */
  int        **imx;             /* convenience ptr into insert states  */
  int        **dmx;             /* convenience ptr into delete states  */
  int          k2;
  char         t2;
  int          s2;
  int          cur, prv, nxt;	/* current, previous, next row index (0 or 1)*/
  int          i,k;		/* indices for seq, model */
  int          sc;		/* integer score */
  int          max;		/* maximum integer score */
  int          start;		/* s1 to start at (need, for STS special case) */

 
  /* Choose our midpoint.
   * Special cases: s1, s3 adjacent and t1 == STS: s2 = s1
   *                s1, s3 adjacent and t3 == STT: s2 = s3
   *                (where we must replace STS, STT eventually)
   */
  s2 = s1 + (s3-s1) / 2;
  if (s3-s1 == 1 && t1 == STS) s2 = s1;
  if (s3-s1 == 1 && t3 == STT) s2 = s3;

  /* STS is a special case. STS aligns to row zero by convention,
   * but we'll be passed s1=1, t1=STS. We have to init on row
   * zero then start DP on row 1.
   */
  start = (t1 == STS) ? 0 : s1;

  /* Allocate our forward two rows.
   * Initialize row zero.
   */
  fwd = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  cur = start%2;
  xmx[cur][XMN] = xmx[cur][XMB] = -INFTY;
  xmx[cur][XME] = xmx[cur][XMC] = -INFTY;  
  for (k = k1; k <= k3; k++)
    mmx[cur][k] = imx[cur][k] = dmx[cur][k] = -INFTY;      

  /* Where to put our zero for our start point...
   * (only possible to start on an emitting state; J disallowed)
   */
  switch (t1) {
  case STM: mmx[cur][k1]  = 0; break;
  case STI: imx[cur][k1]  = 0; break;
  case STN: xmx[cur][XMN] = 0; break;
  case STC: xmx[cur][XMC] = 0; break;
  case STS: xmx[cur][XMN] = 0; break;
  default:  Die("you can't init get_wee_midpt with a %s\n", Statetype(t1));
  }

  /* Still initializing.
   * Deal with pulling horizontal matrix moves in initial row.
   * These are any transitions to nonemitters:
   *    STM-> E, D
   *    STI-> none
   *    STN-> B
   *    STC-> (T, but we never observe this in the forward pass of a d&c)
   *    STE-> C
   *    STS-> (N, already implied by setting xmx[cur][XMN] = 0)
   *    STB-> M
   */ 
  if (t1 == STM)
    {
      for (k = k1+1; k <= k3; k++)
	{				/* transits into STD */
	  dmx[cur][k] = -INFTY;
	  if ((sc = mmx[cur][k-1] + hmm->tsc[TMD][k-1]) > -INFTY)
	    dmx[cur][k] = sc;
	  if ((sc = dmx[cur][k-1] + hmm->tsc[TDD][k-1]) > dmx[cur][k])
	    dmx[cur][k] = sc;
	}
				/* transit into STE */
      xmx[cur][XME] = -INFTY;
      if ((sc = mmx[cur][k1] + hmm->esc[k1]) > -INFTY)
	xmx[cur][XME] = sc;
    }
				/* transit into STB from STN */
  xmx[cur][XMB] = -INFTY;
  if ((sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
    xmx[cur][XMB] = sc;
				/* transit into STC from STE */
  xmx[cur][XMC] = -INFTY;
  if ((sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > -INFTY)
    xmx[cur][XMC] = sc;
  
  /* Done initializing.
   * Start recursive DP; sweep forward to chosen s2 midpoint. Done as a pull.
   */
  for (i = start+1; i <= s2; i++) {
    cur = i % 2;
    prv = !cur;

    mmx[cur][k1] = imx[cur][k1] = dmx[cur][k1] = -INFTY;

    /* Insert state in column k1, and B->M transition in k1.
     */
    if (k1 < hmm->M) {
      imx[cur][k1] = -INFTY;
      if ((sc = mmx[prv][k1] + hmm->tsc[TMI][k1]) > -INFTY)
	imx[cur][k1] = sc;
      if ((sc = imx[prv][k1] + hmm->tsc[TII][k1]) > imx[cur][k1])
	imx[cur][k1] = sc;
      if (hmm->isc[(int) dsq[i]][k1] != -INFTY)
	imx[cur][k1] += hmm->isc[(int) dsq[i]][k1];
      else
	imx[cur][k1] = -INFTY;
    }
    if ((sc = xmx[prv][XMB] + hmm->bsc[k1]) > -INFTY)
      mmx[cur][k1] = sc;
    if (hmm->msc[(int) dsq[i]][k1] != -INFTY)
      mmx[cur][k1] += hmm->msc[(int) dsq[i]][k1];
    else
      mmx[cur][k1] = -INFTY;

    /* Main chunk of recursion across model positions
     */
    for (k = k1+1; k <= k3; k++) {
				/* match state */
      mmx[cur][k]  = -INFTY;
      if ((sc = mmx[prv][k-1] + hmm->tsc[TMM][k-1]) > -INFTY)
	mmx[cur][k] = sc;
      if ((sc = imx[prv][k-1] + hmm->tsc[TIM][k-1]) > mmx[cur][k])
	mmx[cur][k] = sc;
      if ((sc = xmx[prv][XMB] + hmm->bsc[k]) > mmx[cur][k])
	mmx[cur][k] = sc;
      if ((sc = dmx[prv][k-1] + hmm->tsc[TDM][k-1]) > mmx[cur][k])
	mmx[cur][k] = sc;
      if (hmm->msc[(int) dsq[i]][k] != -INFTY)
	mmx[cur][k] += hmm->msc[(int) dsq[i]][k];
      else
	mmx[cur][k] = -INFTY;

				/* delete state */
      dmx[cur][k] = -INFTY;
      if (k < hmm->M) {
	if ((sc = mmx[cur][k-1] + hmm->tsc[TMD][k-1]) > -INFTY)
	  dmx[cur][k] = sc;
	if ((sc = dmx[cur][k-1] + hmm->tsc[TDD][k-1]) > dmx[cur][k])
	  dmx[cur][k] = sc;
      }

				/* insert state */
      imx[cur][k] = -INFTY;
      if (k < hmm->M) {
	if ((sc = mmx[prv][k] + hmm->tsc[TMI][k]) > -INFTY)
	  imx[cur][k] = sc;
	if ((sc = imx[prv][k] + hmm->tsc[TII][k]) > imx[cur][k])
	  imx[cur][k] = sc;
	if (hmm->isc[(int) dsq[i]][k] != -INFTY)
	  imx[cur][k] += hmm->isc[(int) dsq[i]][k];
	else
	  imx[cur][k] = -INFTY;
      }
    }
				/* N state */
    xmx[cur][XMN] = -INFTY;
    if ((sc = xmx[prv][XMN] + hmm->xsc[XTN][LOOP]) > -INFTY)
      xmx[cur][XMN] = sc;
				/* E state */
    xmx[cur][XME] = -INFTY;
    for (k = k1; k <= k3 && k <= hmm->M; k++)
      if ((sc =  mmx[cur][k] + hmm->esc[k]) > xmx[cur][XME])
	xmx[cur][XME] = sc;
				/* B state */
    xmx[cur][XMB] = -INFTY;
    if ((sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > -INFTY)
      xmx[cur][XMB] = sc;
				/* C state */
    xmx[cur][XMC] = -INFTY;
    if ((sc = xmx[prv][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[cur][XMC] = sc;
    if ((sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > xmx[cur][XMC])
      xmx[cur][XMC] = sc;
  }

  /* Row s2%2 in fwd matrix now contains valid scores from s1 (start) to s2,
   * with J transitions disallowed (no cycles through model). 
   */

  /*****************************************************************
   * Backwards pass.
   *****************************************************************/ 

  /* Allocate our backwards two rows. Init last row.
   */
  bck = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  nxt = s3%2;
  xmx[nxt][XMN] = xmx[nxt][XMB] = -INFTY;
  xmx[nxt][XME] = xmx[nxt][XMC] = -INFTY;  
  for (k = k1; k <= k3 + 1; k++)
    mmx[nxt][k] = imx[nxt][k] = dmx[nxt][k] = -INFTY;      
  cur = !nxt;
  mmx[cur][k3+1] = imx[cur][k3+1] = dmx[cur][k3+1] = -INFTY;      

  /* Where to put the zero for our end point on last row.
   */
  switch (t3) {
  case STM: mmx[nxt][k3]  = 0; break;
  case STI: imx[nxt][k3]  = 0; break;
  case STN: xmx[nxt][XMN] = 0; break;
  case STC: xmx[nxt][XMC] = 0; break;   /* must be an emitting C */
  case STT: xmx[nxt][XMC] = hmm->xsc[XTC][MOVE];  break; /* C->T implied */
  default:  Die("you can't init get_wee_midpt with a %s\n", Statetype(t3));
  }

  /* Still initializing.
   * In the case t3==STT, there are a few horizontal moves possible 
   * on row s3, because STT isn't an emitter. All other states are 
   * emitters, so their connections have to be to the previous row s3-1.
   */
  if (t3 == STT) 
    {				/* E->C */
      xmx[nxt][XME] = xmx[nxt][XMC] + hmm->xsc[XTE][MOVE];
				/* M->E */
      for (k = k3; k >= k1; k--) {
	mmx[nxt][k] = xmx[nxt][XME] + hmm->esc[k];
	if (s3 != s2)
	  mmx[nxt][k] += hmm->msc[(int)dsq[s3]][k];
      }
    }

  /* Start recursive DP; sweep backwards to chosen s2 midpoint.
   * Done as a pull. M, I scores at current row do /not/ include
   * emission scores. Be careful of integer underflow.
   */
  for (i = s3-1; i >= s2; i--) {
				/* note i < L, so i+1 is always a legal index */
    cur = i%2;
    nxt = !cur;
				/* C pulls from C (T is special cased) */
    xmx[cur][XMC] = -INFTY;
    if ((sc = xmx[nxt][XMC] + hmm->xsc[XTC][LOOP]) > -INFTY)
      xmx[cur][XMC] = sc;
				/* B pulls from M's */
    xmx[cur][XMB] = -INFTY;
    for (k = k1; k <= k3; k++)
      if ((sc = mmx[nxt][k] + hmm->bsc[k]) > xmx[cur][XMB])
	xmx[cur][XMB] = sc;
				/* E pulls from C (J disallowed) */
    xmx[cur][XME] = -INFTY;
    if ((sc = xmx[cur][XMC] + hmm->xsc[XTE][MOVE]) > -INFTY)
      xmx[cur][XME] = sc;
				/* N pulls from B, N */
    xmx[cur][XMN] = -INFTY;
    if ((sc = xmx[cur][XMB] + hmm->xsc[XTN][MOVE]) > -INFTY)
      xmx[cur][XMN] = sc;
    if ((sc = xmx[nxt][XMN] + hmm->xsc[XTN][LOOP]) > xmx[cur][XMN])
      xmx[cur][XMN] = sc;

    /* Main recursion across model
     */
    for (k = k3; k >= k1; k--)  {
				/* special case k == M */
      if (k == hmm->M) {
	mmx[cur][k] = xmx[cur][XME]; /* p=1 transition to E by definition */
	dmx[cur][k] = -INFTY;	/* doesn't exist */
	imx[cur][k] = -INFTY;	/* doesn't exist */
	if (i != s2)
	  mmx[cur][k] += hmm->msc[(int)dsq[i]][k];
	continue;		
      }    	/* below this k < M, so k+1 is a legal index */

				/* pull into match state */
      mmx[cur][k] = -INFTY;
      if ((sc = xmx[cur][XME] + hmm->esc[k]) > -INFTY)
	mmx[cur][k] = sc; 
      if ((sc = mmx[nxt][k+1] + hmm->tsc[TMM][k]) > mmx[cur][k])
	mmx[cur][k] = sc; 
      if ((sc = imx[nxt][k] + hmm->tsc[TMI][k]) > mmx[cur][k])
	mmx[cur][k] = sc; 
      if ((sc = dmx[cur][k+1] + hmm->tsc[TMD][k]) > mmx[cur][k])
	mmx[cur][k] = sc;
      if (i != s2) 
	mmx[cur][k] += hmm->msc[(int)dsq[i]][k];

				/* pull into delete state */
      dmx[cur][k] = -INFTY;
      if ((sc = mmx[nxt][k+1] + hmm->tsc[TDM][k]) > -INFTY)
	dmx[cur][k] = sc;
      if ((sc = dmx[cur][k+1] + hmm->tsc[TDD][k]) > dmx[cur][k])
	dmx[cur][k] = sc;
				/* pull into insert state */
      imx[cur][k] = -INFTY;
      if ((sc = mmx[nxt][k+1] + hmm->tsc[TIM][k]) > -INFTY)
	imx[cur][k] = sc;
      if ((sc = imx[nxt][k] + hmm->tsc[TII][k]) > imx[cur][k])
	imx[cur][k] = sc;
      if (i != s2)
	imx[cur][k] += hmm->isc[(int)dsq[i]][k];
      
    }
  }
   
  /*****************************************************************
   * DP complete; we have both forward and backward passes. Now we
   * look across the s2 row and find the optimal emitting state.
   *****************************************************************/  

  cur = s2%2;
  max = -INFTY;
  for (k = k1; k <= k3; k++)
    {
      if ((sc = fwd->mmx[cur][k] + bck->mmx[cur][k]) > max)
	{ k2 = k; t2 = STM; max = sc; }
      if ((sc = fwd->imx[cur][k] + bck->imx[cur][k]) > max)
	{ k2 = k; t2 = STI; max = sc; }
    }
  if ((sc = fwd->xmx[cur][XMN] + bck->xmx[cur][XMN]) > max)
    { k2 = 1;        t2 = STN; max = sc; }
  if ((sc = fwd->xmx[cur][XMC] + bck->xmx[cur][XMC]) > max)
    { k2 = hmm->M;   t2 = STC; max = sc; }

  /*****************************************************************
   * Garbage collection, return.
   *****************************************************************/
  
  FreePlan7Matrix(fwd);
  FreePlan7Matrix(bck);
  *ret_k2 = k2;
  *ret_t2 = t2;
  *ret_s2 = s2;
  return Scorify(max);
}


/* Function: P7ViterbiAlignAlignment()
 * Date:     SRE, Sat Jul  4 13:39:00 1998 [St. Louis]
 *
 * Purpose:  Align a multiple alignment to an HMM without
 *           changing the multiple alignment itself.
 *           Adapted from P7Viterbi().
 *           
 *           Heuristic; not a guaranteed optimal alignment.
 *           Guaranteeing an optimal alignment appears difficult.
 *           [cryptic note to myself:] In paths connecting to I* metastates,
 *           recursion breaks down; if there is a gap in the
 *           previous column for a given seq, we can't determine what state the
 *           I* metastate corresponds to for this sequence, unless we
 *           look back in the DP matrix. The lookback would either involve
 *           recursing back to the previous M* metastate (giving a
 *           O(MN^2) algorithm instead of O(MN)) or expanding the I*
 *           metastate into 3^nseq separate I* metastates to keep track
 *           of which of three states each seq is in. Since the second
 *           option blows up exponentially w/ nseq, it is not attractive.
 *           If the first option were used, the correct algorithm would be related to 
 *           modelmakers.c:Maxmodelmaker(), but somewhat more difficult.
 *           
 *           The heuristic approach here is to calculate a "consensus"
 *           sequence from the alignment, and align the consensus to the HMM.
 *           Some hackery is employed, weighting transitions and emissions
 *           to make things work (re: con and mocc arrays).
 *
 * Args:     aseq  - aligned sequences
 *           ainfo - info for aseqs (includes alen, nseq, wgt)
 *           hmm   - model to align to        
 *
 * Returns:  Traceback. Caller must free with P7FreeTrace().
 *           pos[] contains alignment columns, indexed 1..alen.
 *           statetype[] contains metastates M*, etc. as STM, etc.
 */
struct p7trace_s *
P7ViterbiAlignAlignment(MSA *msa, struct plan7_s *hmm)
{
  struct dpmatrix_s *mx;        /* Viterbi calculation lattice (two rows) */
  struct dpshadow_s *tb;        /* shadow matrix of traceback pointers */
  struct p7trace_s  *tr;        /* RETURN: traceback */
  int  **xmx, **mmx, **imx, **dmx;
  char **xtb, **mtb, **itb, **dtb;
  float **con;                  /* [1..alen][0..Alphabet_size-1], consensus counts */
  float  *mocc;                 /* fractional occupancy of a column; used to weight transitions */
  int     i;			/* counter for columns */
  int     k;			/* counter for model positions */
  int     idx;			/* counter for seqs */
  int     sym;			/* counter for alphabet symbols */
  int     sc;			/* temp variable for holding score */
  float   denom;		/* total weight of seqs; used to "normalize" counts */
  int     cur, prv;

  /* The "consensus" is a counts matrix, [1..alen][0..Alphabet_size-1].
   * Gaps are not counted explicitly, but columns with lots of gaps get
   * less total weight because they have fewer counts.
   */
				/* allocation */
  con  = MallocOrDie(sizeof(float *) * (msa->alen+1));
  mocc = MallocOrDie(sizeof(float)   * (msa->alen+1));
  for (i = 1; i <= msa->alen; i++) {
    con[i] = MallocOrDie(sizeof(float) * Alphabet_size);
    FSet(con[i], Alphabet_size, 0.0);
  }
  mocc[0] = -9999.;
				/* initialization */
				/* note: aseq is off by one, 0..alen-1 */
				/* "normalized" to have a max total count of 1 per col */
  denom = FSum(msa->wgt, msa->nseq);
  for (i = 1; i <= msa->alen; i++)
    {
      for (idx = 0; idx < msa->nseq; idx++)
	if (! isgap(msa->aseq[idx][i-1]))
	  P7CountSymbol(con[i], SYMIDX(msa->aseq[idx][i-1]), msa->wgt[idx]);
      FScale(con[i], Alphabet_size, 1./denom);
      mocc[i] = FSum(con[i], Alphabet_size);
    }
  
  /* Allocate a DP matrix with 2 rows, 0..M columns,
   * and a shadow matrix with 0,1..alen rows, 0..M columns.
   */ 
  mx = AllocPlan7Matrix(2, hmm->M, &xmx, &mmx, &imx, &dmx);
  tb = AllocShadowMatrix(msa->alen+1, hmm->M, &xtb, &mtb, &itb, &dtb);

  /* Initialization of the zero row.
   */
  xmx[0][XMN] = 0;		                     /* S->N, p=1            */
  xtb[0][XMN] = STS;
  xmx[0][XMB] = hmm->xsc[XTN][MOVE];                 /* S->N->B, no N-tail   */
  xtb[0][XMB] = STN;
  xmx[0][XME] = xmx[0][XMC] = xmx[0][XMJ] = -INFTY;  /* need seq to get here */
  tb->esrc[0] = 0;
  xtb[0][XMC] = xtb[0][XMJ] = STBOGUS;  
  for (k = 0; k <= hmm->M; k++) {
    mmx[0][k] = imx[0][k] = dmx[0][k] = -INFTY;      /* need seq to get here */
    mtb[0][k] = itb[0][k] = dtb[0][k] = STBOGUS;
  }
  
  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   *    D_M and I_M are wastefully calculated (they don't exist)
   */
  for (i = 1; i <= msa->alen; i++) {
    cur = i % 2;
    prv = ! cur;

    mmx[cur][0] = imx[cur][0] = dmx[cur][0] = -INFTY;
    mtb[i][0]   = itb[i][0]   = dtb[i][0]   = STBOGUS;

    for (k = 1; k <= hmm->M; k++) {
				/* match state */
      mmx[cur][k]  = -INFTY;
      mtb[i][k]    = STBOGUS;
      if (mmx[prv][k-1] > -INFTY && hmm->tsc[TMM][k-1] > -INFTY &&
	  (sc = mmx[prv][k-1] + hmm->tsc[TMM][k-1]) > mmx[cur][k])
	{ mmx[cur][k] = sc; mtb[i][k] = STM; }
      if (imx[prv][k-1] > -INFTY && hmm->tsc[TIM][k-1] > -INFTY &&
	  (sc = imx[prv][k-1] + hmm->tsc[TIM][k-1] * mocc[i-1]) > mmx[cur][k])
	{ mmx[cur][k] = sc; mtb[i][k] = STI; }
      if ((sc = xmx[prv][XMB] + hmm->bsc[k]) > mmx[cur][k])
	{ mmx[cur][k] = sc; mtb[i][k] = STB; }
      if (dmx[prv][k-1] > -INFTY && hmm->tsc[TDM][k-1] > -INFTY &&
	  (sc = dmx[prv][k-1] + hmm->tsc[TDM][k-1]) > mmx[cur][k])
	{ mmx[cur][k] = sc; mtb[i][k] = STD; }
				/* average over "consensus" sequence */
      for (sym = 0; sym < Alphabet_size; sym++)
	{
	  if (con[i][sym] > 0 && hmm->msc[sym][k] == -INFTY) { mmx[cur][k] = -INFTY; break; }
	  mmx[cur][k] += hmm->msc[sym][k] * con[i][sym];
	}

				/* delete state */
      dmx[cur][k] = -INFTY;
      dtb[i][k]   = STBOGUS;
      if (mmx[cur][k-1] > -INFTY && hmm->tsc[TMD][k-1] > -INFTY &&
	  (sc = mmx[cur][k-1] + hmm->tsc[TMD][k-1]) > dmx[cur][k])
	{ dmx[cur][k] = sc; dtb[i][k] = STM; }
      if (dmx[cur][k-1] > -INFTY && hmm->tsc[TDD][k-1] > -INFTY &&
	  (sc = dmx[cur][k-1] + hmm->tsc[TDD][k-1]) > dmx[cur][k])
	{ dmx[cur][k] = sc; dtb[i][k] = STD; }

				/* insert state */
      if (k < hmm->M) {
	imx[cur][k] = -INFTY;
	itb[i][k]   = STBOGUS;
	if (mmx[prv][k] > -INFTY && hmm->tsc[TMI][k] > -INFTY &&
	    (sc = mmx[prv][k] + hmm->tsc[TMI][k] * mocc[i]) > imx[cur][k])
	  { imx[cur][k] = sc; itb[i][k] = STM; }
	if (imx[prv][k] > -INFTY && hmm->tsc[TII][k] > -INFTY &&
	    (sc = imx[prv][k] + hmm->tsc[TII][k] * mocc[i-1] * mocc[i]) > imx[cur][k])
	  { imx[cur][k] = sc; itb[i][k] = STI; }
				/* average over "consensus" sequence */
	for (sym = 0; sym < Alphabet_size; sym++)
	  {
	    if (con[i][sym] > 0 && hmm->isc[sym][k] == -INFTY) { imx[cur][k] = -INFTY; break; }
	    imx[cur][k] += hmm->isc[sym][k] * con[i][sym];
	  }
      }
    }
    
    /* Now the special states. Order is important here.
     * remember, N, C, and J emissions are zero score by definition.
     */
				/* N state */
    xmx[cur][XMN] = -INFTY;
    xtb[i][XMN]   = STBOGUS;
    if (xmx[prv][XMN] > -INFTY && hmm->xsc[XTN][LOOP] > -INFTY &&
	(sc = xmx[prv][XMN] + hmm->xsc[XTN][LOOP] * mocc[i]) > -INFTY)
      { xmx[cur][XMN] = sc; xtb[i][XMN] = STN; }
				/* E state */
    xmx[cur][XME] = -INFTY;
    xtb[i][XME]   = STBOGUS;
    for (k = 1; k <= hmm->M; k++)
      if (mmx[cur][k] > -INFTY && hmm->esc[k] > -INFTY &&
	  (sc =  mmx[cur][k] + hmm->esc[k]) > xmx[cur][XME])
	{ xmx[cur][XME] = sc; tb->esrc[i] = k; }

				/* we don't check J state */
				/* B state; don't connect from J */
    xmx[cur][XMB] = -INFTY;
    xtb[i][XMB]   = STBOGUS;
    if (xmx[cur][XMN] > -INFTY && hmm->xsc[XTN][MOVE] > -INFTY &&
	(sc = xmx[cur][XMN] + hmm->xsc[XTN][MOVE]) > xmx[cur][XMB])
      { xmx[cur][XMB] = sc; xtb[i][XMB] = STN; }

				/* C state */
    xmx[cur][XMC] = -INFTY;
    xtb[i][XMC]   = STBOGUS;
    if (xmx[prv][XMC] > -INFTY && hmm->xsc[XTC][LOOP] > -INFTY &&
	(sc = xmx[prv][XMC] + hmm->xsc[XTC][LOOP] * mocc[i]) > -INFTY)
      { xmx[cur][XMC] = sc; xtb[i][XMC] = STC; }
    if (xmx[cur][XME] > -INFTY && hmm->xsc[XTE][MOVE] > -INFTY &&
	(sc = xmx[cur][XME] + hmm->xsc[XTE][MOVE]) > xmx[cur][XMC])
      { xmx[cur][XMC] = sc; xtb[i][XMC] = STE; }
  }
				/* T state (not stored in mx) */
  sc = xmx[msa->alen%2][XMC] + hmm->xsc[XTC][MOVE];

				/* do the traceback */
  tr = ShadowTrace(tb, hmm, msa->alen);
				/* cleanup and return */
  FreePlan7Matrix(mx);
  FreeShadowMatrix(tb);
  for (i = 1; i <= msa->alen; i++)
    free(con[i]);
  free(con);
  free(mocc);

  return tr;
}



/* Function: ShadowTrace()
 * Date:     SRE, Sun Jul  5 11:38:24 1998 [St. Louis]
 *
 * Purpose:  Given a shadow matrix, trace it back, and return
 *           the trace.
 *
 * Args:     tb  - shadow matrix of traceback pointers
 *           hmm - the model (needed for figuring out wing unfolding)
 *           L   - sequence length
 *
 * Returns:  traceback. Caller must free w/ P7FreeTrace().
 */
struct p7trace_s *
ShadowTrace(struct dpshadow_s *tb, struct plan7_s *hmm, int L)
{
  struct p7trace_s *tr;
  int curralloc;		/* current allocated length of trace */
  int tpos;			/* position in trace */
  int i;			/* position in seq (1..N) */
  int k;			/* position in model (1..M) */
  char nxtstate;        	/* next state to assign in traceback */

  /* Overallocate for the trace.
   * S-N-B- ... - E-C-T  : 6 states + L is minimum trace;
   * add L more as buffer.             
   */
  curralloc = L * 2 + 6; 
  P7AllocTrace(curralloc, &tr);

  /* Initialization of trace
   * We do it back to front; ReverseTrace() is called later.
   */
  tr->statetype[0] = STT;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tpos     = 1;
  i        = L;			/* current i (seq pos) we're trying to assign   */
  k        = 0;			/* current k (model pos) we're trying to assign */
  nxtstate = STC;		/* assign the C state first, for C->T */

  /* Traceback
   */
  while (nxtstate != STS) {
    switch (nxtstate) {
    case STM:
      tr->statetype[tpos] = STM;
      nxtstate            = tb->mtb[i][k];
      tr->nodeidx[tpos]   = k--;
      tr->pos[tpos]       = i--;
      tpos++;
      break;

    case STI:
      tr->statetype[tpos] = STI;
      nxtstate            = tb->itb[i][k];
      tr->nodeidx[tpos]   = k;
      tr->pos[tpos]       = i--;
      tpos++;
      break;
      
    case STD:
      tr->statetype[tpos] = STD;
      nxtstate            = tb->dtb[i][k];
      tr->nodeidx[tpos]   = k--;
      tr->pos[tpos]       = 0;
      tpos++;
      break;

    case STN:	
      tr->statetype[tpos] = STN;
      nxtstate            = tb->xtb[i][XMN];
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]  = (nxtstate == STN) ? i-- : 0; /* N->N; 2nd one emits. */
      tpos++;
      break;

    case STB:
				/* Check for wing unfolding */
      if (Prob2Score(hmm->begin[k+1], hmm->p1) + 1 * INTSCALE <= hmm->bsc[k+1])
	while (k > 0)
	  {
	    tr->statetype[tpos] = STD;
	    tr->nodeidx[tpos]   = k--;
	    tr->pos[tpos]       = 0;
	    tpos++;
	    if (tpos == curralloc) 
	      {				/* grow trace if necessary  */
		curralloc += L;
		P7ReallocTrace(tr, curralloc);
	      }
	  }

      tr->statetype[tpos] = STB;
      nxtstate            = tb->xtb[i][XMB];
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = 0;
      tpos++;
      break;

    case STJ:
      tr->statetype[tpos] = STJ;
      nxtstate            = tb->xtb[i][XMJ];
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]  = (nxtstate == STJ) ? i-- : 0; /* J->J; 2nd one emits. */
      tpos++;
      break;

    case STE:			
      tr->statetype[tpos] = STE;
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]       = 0;
      k                   = tb->esrc[i];
      nxtstate            = STM;
      tpos++;
				/* check for wing unfolding */
      if (Prob2Score(hmm->end[k], 1.) + 1*INTSCALE <=  hmm->esc[k])
	{
	  int dk;		/* need a tmp k while moving thru delete wing */
	  for (dk = hmm->M; dk > k; dk--)
	    {
	      tr->statetype[tpos] = STD;
	      tr->nodeidx[tpos]   = dk;
	      tr->pos[tpos]       = 0;
	      tpos++;
	      if (tpos == curralloc) 
		{				/* grow trace if necessary  */
		  curralloc += L;
		  P7ReallocTrace(tr, curralloc);
		}
	    }
	}
      break;

    case STC:	
      tr->statetype[tpos] = STC;
      nxtstate            = tb->xtb[i][XMC];
      tr->nodeidx[tpos]   = 0;
      tr->pos[tpos]  = (nxtstate == STC) ? i-- : 0; /* C->C; 2nd one emits. */
      tpos++;
      break;

    default:
      Die("HMMER: Bad state (%s) in ShadowTrace()\n", Statetype(nxtstate));

    } /* end switch over nxtstate */
    
    if (tpos == curralloc) 
      {				/* grow trace if necessary  */
	curralloc += L;
	P7ReallocTrace(tr, curralloc);
      }

  } /* end traceback, just before assigning S state */
  
  tr->statetype[tpos] = STS;
  tr->nodeidx[tpos]   = 0;
  tr->pos[tpos]       = 0;
  tr->tlen            = tpos + 1;

  P7ReverseTrace(tr);
  return tr;
}



/* Function: PostprocessSignificantHit()
 * Date:     SRE, Wed Dec 20 12:11:01 2000 [StL]
 * 
 * Purpose:  Add a significant hit to per-seq and per-domain hit
 *           lists, after postprocessing the scores appropriately,
 *           and making sure per-domain scores add up to the per-seq 
 *           score.
 *  
 *          [doesn't really belong in core_algorithms.c, because
 *           it's more of a hack than an algorithm, but on the other
 *           hand it's now part of the core of how HMMER scores
 *           things. Maybe there should be a core_hacks.c.]
 *
 *           Given: active hit lists for per-seq and per-domain
 *           scores (e.g. hmmpfam and hmmsearch, collating their
 *           results), and a new hit that's significant enough
 *           that it may need to be reported in final output.
 *           
 *           Breaks the traceback into individual domain traces;
 *           scores each one of them, then applies null2 correction
 *           for biased composition. Recalculates the per-seq score
 *           as the sum of the per-domain scores. Stores the hits
 *           in the lists, for eventual sorting and output by the
 *           caller.
 * 
 * Notes:    In principle we've got the score, and a pvalue, and a traceback
 *           by doing the Viterbi algorithm, right? What else is left
 *           to do? Well, in practice, life is more complicated, because
 *           of the trace-dependent null2 score correction.
 *           
 *           After a null2 score correction is carried out on
 *           each domain (the default) the number of detected domains
 *           with scores > 0 may have decreased. We want the
 *           global (per-seq) hit list to have the recalculated number of 
 *           domains, not necessarily what Viterbi gave us. 
 *           
 *           Also, since we want the global score to be the sum of
 *           the individual domains, but the null2 correction is
 *           applied to each domain individually, we have to calculate
 *           an adjusted global score. (To do otherwise invites
 *           subtle inconsistencies; xref bug 2.)
 *
 *           We don't have final evalues, so we may put a few
 *           more hits into the hit lists than we end up reporting.
 *           The main output routine is responsible for final
 *           enforcement of the thresholds.
 *          
 *           This routine is NOT THREADSAFE. When multithreaded,
 *           with using shared ghit/dhit output buffers, calls to
 *           PostprocessSignificantHit() need to be protected.
 *           
 * Args:     ghit     - an active list of per-seq (global) hits
 *           dhit     - an active list of per-domain hits
 *           tr       - the significant HMM/seq traceback we'll report on
 *           hmm      - ptr to the HMM
 *           dsq      - digitized sequence (1..L)
 *           L        - length of dsq
 *           seqname  - name of sequence (same as targname, in hmmsearch)
 *           seqacc   - seq's accession (or NULL)
 *           seqdesc  - seq's description (or NULL)
 *           do_forward  - TRUE if we've already calculated final per-seq score
 *           sc_override - per-seq score to use if do_forward is TRUE
 *           do_null2    - TRUE to apply the null2 scoring correction
 *           thresh      - contains the threshold/cutoff information.
 *           hmmpfam_mode - TRUE if called by hmmpfam, else assumes hmmsearch;
 *                          affects how the lists' sort keys are set.
 *
 * Returns:  the recalculated per-seq score (or sc_override),
 *           as appropriate, for subsequent storage in the histogram
 */
float
PostprocessSignificantHit(struct tophit_s    *ghit, 
			  struct tophit_s    *dhit,
			  struct p7trace_s   *tr,
			  struct plan7_s     *hmm,
			  char               *dsq,
			  int                 L,
			  char               *seqname,
			  char               *seqacc,
			  char               *seqdesc,
			  int                 do_forward,
			  float               sc_override,
			  int                 do_null2,
			  struct threshold_s *thresh,
			  int                 hmmpfam_mode)
{
  struct p7trace_s **tarr;      /* array of per-domain traces */
  struct fancyali_s *ali;       /* alignment of a domain      */ 
  int ntr;			/* number of domain traces from Viterbi */
  int tidx;			/* index for traces (0..ntr-1) */
  int ndom;			/* # of domains accepted in sequence */
  int didx;			/* index for domains (1..ndom) */
  int k1, k2;			/* start, stop coord in model */
  int i1, i2;			/* start, stop in sequence    */
  float   whole_sc;		/* whole sequence score = \sum domain scores */
  float  *score;                /* array of raw scores for each domain */
  int    *usedomain;            /* TRUE if this domain is accepted */ 
  double  whole_pval;
  double  pvalue;
  double  sortkey;

  /* Special case: rarely, the alignment was totally impossible
   * and tr is NULL.
   */
  if (tr == NULL) return whole_sc;

  /* Break the trace into one or more individual domains.
   */
  TraceDecompose(tr, &tarr, &ntr);
  if (ntr == 0) Die("TraceDecompose() screwup"); /* "can't happen" (!) */

  /* Rescore each domain, apply null2 correction if asked.
   * Mark positive-scoring ones (we'll definitely report those),
   * and include their score in the whole sequence score.
   */
  score     = MallocOrDie(sizeof(float) * ntr);
  usedomain = MallocOrDie(sizeof(int)   * ntr);
  ndom      = 0;
  whole_sc  = 0.;
  for (tidx = 0; tidx < ntr; tidx++)
    {
      score[tidx]  = P7TraceScore(hmm, dsq, tarr[tidx]);
      if (do_null2) score[tidx] -= TraceScoreCorrection(hmm, tarr[tidx], dsq);
      if (score[tidx] > 0.0) {
	usedomain[tidx] = TRUE;
	ndom++;
	whole_sc += score[tidx];
      } else 
	usedomain[tidx] = FALSE;
    }

  /* Make sure at least one positive scoring domain is in
   * the trace. If not, invoke "weak single domain" rules:
   * we will always report at least one domain per sequence, even
   * if it has a negative score. (HMMER's Plan7 architecture can report
   * one negative scoring domain but not more.)
   */
  if (ndom == 0) {
    tidx            = FArgMax(score, ntr);
    usedomain[tidx] = TRUE;
    whole_sc        = score[tidx];
    ndom            = 1;
  }

  /* Implement --do_forward: override the trace-dependent sum-of-domain
   * whole score, use the P7Forward() score that the called passed
   * us instead. This is a hack; null2 is trace-dependent and
   * thus undefined for P7Forward() scoring; see commentary in hmmpfam.c. 
   */
  if (do_forward) whole_sc = sc_override;

  /* Go through and put all the accepted domains into the hit list.
   */
  whole_pval = PValue(hmm, whole_sc);
  for (tidx = 0, didx = 1; tidx < ntr; tidx++) {
    if (! usedomain[tidx]) continue;

    TraceSimpleBounds(tarr[tidx], &i1, &i2, &k1, &k2);
    pvalue = PValue(hmm, score[tidx]); 

    if (pvalue <= thresh->domE && score[tidx] >= thresh->domT) {
      ali     = CreateFancyAli(tarr[tidx], hmm, dsq, seqname);

      if (hmmpfam_mode) 
	sortkey = -1.*(double)i1; /* hmmpfam: sort on position in seq    */
      else
	sortkey = score[tidx];	  /* hmmsearch: sort on E (monotonic w/ sc) */

      RegisterHit(dhit, sortkey,
		  pvalue,     score[tidx], 
		  whole_pval, whole_sc,
		  hmmpfam_mode ? hmm->name : seqname,
		  hmmpfam_mode ? hmm->acc  : seqacc,
		  hmmpfam_mode ? hmm->desc : seqdesc,
		  i1,i2, L, 
		  k1,k2, hmm->M, 
		  didx,ndom,ali);
    }
    didx++;
  }

  /* Now register the global hit, with the domain-derived score.
   */

  /* sorting: 
   * hmmpfam has to worry that score and E-value are not monotonic
   * when multiple HMMs (with different EVD parameters) are potential
   * targets. Therefore in hmmpfam_mode we apply a weird hack
   * to sort primarily on E-value, but on score
   * for really good hits with E=0.0... works because we can
   * assume 100000. > -log(DBL_MIN).
   * hmmsearch simply sorts on score (which for a single HMM, we
   * know is monotonic with E-value).
   */  
  if (hmmpfam_mode)
    sortkey = (whole_pval > 0.0) ? -1.*log(whole_pval) : 100000. + whole_sc;
  else
    sortkey = whole_sc;

  /* Note: we've recalculated whole_sc and it may have decreased
   *       after the null2 correction was applied. For Pfam GA, TC,
   *       or NC cutoffs, we have to be sure that everything on the
   *       hitlist is correct (the hmmpfam output routine assumes it, 
   *       otherwise it would have to reload each HMM to get its
   *       cutoffs). In all other cases, though, we don't care if
   *       the hit list has a bit too many things on it, because the
   *       output routine in hmmsearch or hmmpfam will check against
   *       the cutoffs. Hence we only need to check against globT
   *       (it may be set by GA, TC, or NC) but not globE.
   *                 - SRE, CSHL genome mtg May 2001
   */ 
  if (whole_sc >= thresh->globT) {
    RegisterHit(ghit, sortkey,
		whole_pval, whole_sc,
		0., 0.,	                  /* no mother seq */
		hmmpfam_mode ? hmm->name : seqname,
		hmmpfam_mode ? hmm->acc  : seqacc,
		hmmpfam_mode ? hmm->desc : seqdesc,
		0,0,0,                	  /* seq positions  */
		0,0,0,	                  /* HMM positions  */
		0, ndom,	          /* # domains info    */
		NULL);	                  /* alignment info */
  }

  /* Clean up and return.
   */
  for (tidx = 0; tidx < ntr; tidx++)
    P7FreeTrace(tarr[tidx]);
  free(tarr);
  free(score);
  free(usedomain);
  return whole_sc;
}
