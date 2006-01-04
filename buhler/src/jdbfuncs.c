/*
 * JDBFUNCS.C
 * Implementation of the core HMMer algorithms for the JDB data structures.
 * This implementation is made possible by Christopher Swope's infrastructure
 * hacking on HMMER.
 *
 * For information, please contact jbuhler@cse.wustl.edu.
 */

#include "config.h"

#include "plan7.h"
#include "structs.h"
#include "funcs.h"

/* local prototype */
static void  ViterbiTrace(struct plan7_s *hmm, unsigned char *dsq, int N,
			  cust_dpmatrix_s *mx, struct p7trace_s **ret_tr);



/* Function:  ViterbiSpaceOK()
 *
 * Purpose:   Returns TRUE if the existing matrix allocation
 *            is already sufficient to hold the requested MxN, or if
 *            the matrix can be expanded in M and/or N without
 *            exceeding RAMLIMIT megabytes. 
 *
 * JB: It's just silly that we have to reimplement this, but
 * it's the only way not to expose the resizing mechanism
 * (the padding and max stuff) to the core code.  It would
 * be preferable to have ViterbiSize take mx and do
 * the right thing, but there are places where it is called
 * without reference to a matrix.  This can be worked around,
 * but we haven't done so yet.
 */
int ViterbiSpaceOK(int L, int M, cust_dpmatrix_s *mx)
{
  int newM;
  int newN;
  
  if (M <= mx->maxM && L <= mx->maxN) return TRUE;
  
  if (M > mx->maxM) newM = M + mx->padM; else newM = mx->maxM;
  if (L > mx->maxN) newN = L + mx->padN; else newN = mx->maxN;
  
  if (ViterbiSize(newN, newM) <= RAMLIMIT)
    return TRUE;
  else
    return FALSE;
}


/*
 * ViterbiSize()
 * Given size of input sequence (L) and motif (M), compute how much space will
 * be needed to allocate the DP matrix in the Viterbi algorithm.
 */
int ViterbiSize(int L, int M)
{
  float Mbytes;
  
  Mbytes  = (float) sizeof(cust_dpmatrix_s);
  Mbytes += N_MSTATES * (float) (L+1) * (float) (M+1) * (float) sizeof(int);
  Mbytes += N_XSTATES * (float) (L+1) * (float) sizeof(int);
  Mbytes += (float) (L+1) * (float) sizeof(int *);
  Mbytes /= 1048576.0;
  return (int) Mbytes;
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
  float sc = 0.0;
  
  if (ViterbiSpaceOK(L, hmm->M, mx))
    {
      SQD_DPRINTF1(("   ... using normal P7Viterbi(); size ~%d MB\n", 
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

/***********************************************************************
 * Local sum-of-logs implementation for Forward and Backward. This
 * version uses a more precise table than the one in mathsupport.c and
 * moves the firsttime test out of the high-traffice ilogsum function
 * into a separate initialization procedure.
 *
 * Inlining ilogsum is possible but does not appear helpful.
 ***********************************************************************/

static int ilogsum_lookup[LOGSUM_TBL];

static void init_ilogsum(void)
{
  static int firsttime = TRUE;
  if (!firsttime)  return;
  
  firsttime = FALSE;
    
  {
    const double ln2 = log(2.0);
    double iln2      = 1.0 / ln2;
    int i;
    
    for (i = 0; i < LOGSUM_TBL; i++) 
      ilogsum_lookup[i] = rint(INTSCALE * iln2 *
			       (log(1.+exp(ln2 * (double) -i/INTSCALE))));
  }
}

static int ilogsum(int p1, int p2)
{
  int diff = p1-p2;
  int result;
  
  if (diff >= LOGSUM_TBL)
    { result = p1; }
  else if (diff <= -LOGSUM_TBL)
    { result = p2; }
  else if (diff > 0)
    { result = p1 + ilogsum_lookup[diff]; }
  else
    { result = p2 + ilogsum_lookup[-diff]; }
  
  return result;
} 

/***********************************************************************
 * Portably optimized Viterbi, Forward, and Backward algorithms, plus
 * an unoptimized ViterbiTrace.  The main wins in this implementation
 * come from:
 *   - rearranging data structures to reduce the number
 *     of registers needed in the inner loop
 *   - eliminating branches from the inner loop, by
 *     unrolling the Mth iteration in Viterbi and by
 *     replacing a bunch of "if" tests with MAX.
 *   - exposing opportunities for hoisting and strength reduction
 *     to the compiler.
 *
 * These Viterbi and Forward implementations were found to produce
 * very good inner loop code on an AMD Opteron, and decent code
 * on a Pentium 4, with gcc 3.4.4.  On x86 and x86-64, one may expect 
 * speedups of around 2x for Viterbi over the "fast" implementation
 * from HMMer 2.3.2.  The Forward implementation is also faster,
 * but somewhat less so.
 **********************************************************************/

#undef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define MMX(i,k) (dp[(i)][(k) * N_MSTATES + _MMX])
#define IMX(i,k) (dp[(i)][(k) * N_MSTATES + _IMX])
#define DMX(i,k) (dp[(i)][(k) * N_MSTATES + _DMX])
#define XMX(i,s) (xmx[(i) * N_XSTATES + (s)])

#define TSC(s,k) (tsc[(k) * N_TMOVES + (s)])
#define BSC(k)   (tsc[(k) * N_TMOVES + _BSC])
#define ESC(k)   (tsc[(k) * N_TMOVES + _ESC])
#define MSC(k)   (misc[(k) * N_MISC + _MSC])
#define ISC(k)   (misc[(k) * N_MISC + _ISC])

float
Viterbi(unsigned char *dsq, int L, struct plan7_s *hmm, cust_dpmatrix_s *mx, 
	struct p7trace_s **ret_tr)
{
  struct p7trace_s *tr;
  const struct logodds_s *lom = hmm->lom;
  int const *tsc = lom->tsc;
  
  int M = hmm->M;
  int **dp;
  int *xmx;
  int final_score;
  int i, k;
  
  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 

  ResizeDPMatrix(mx, L, M);
  dp = mx->dp;
  xmx = mx->xmx;
  
  /* Initialization of the zero row.
   * Note that xmx[i][stN] = 0 by definition for all i,
   *    and xmx[i][stT] = xmx[i][stC], so neither stN nor stT need
   *    to be calculated in DP matrices.
   */
  XMX(0,_XMN) = 0;		                     /* S->N, p=1            */
  XMX(0,_XMB) = lom->xsc[_XTN][_MOVE];               /* S->N->B, no N-tail   */
  XMX(0,_XME) = XMX(0,_XMC) = XMX(0,_XMJ) = -INFTY;  /* need seq to get here */
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = -INFTY;         /* need seq to get here */
  
  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   */
  for (i = 1; i <= L; i++) 
    {
      int const * misc = lom->misc[dsq[i]];
      int sc;
	  
      MMX(i,0) = IMX(i,0) = DMX(i,0) = -INFTY;
       XMX(i,_XME) = -INFTY;
      
      for (k = 1; k < M; k++) 
	{
	  /* match state */
	  sc =         MMX(i-1,k-1)  + TSC(_TMM,k-1);
	  sc = MAX(sc, IMX(i-1,k-1)  + TSC(_TIM,k-1));
	  sc = MAX(sc, DMX(i-1,k-1)  + TSC(_TDM,k-1));
	  sc = MAX(sc, XMX(i-1,_XMB) + BSC(k));
	  
	  MMX(i,k) = MAX(sc + MSC(k), -INFTY);
	  
	  /* E state update */
	  XMX(i,_XME) = MAX(XMX(i,_XME), MMX(i,k) + ESC(k));
	  
	  /* insert state */
	  sc =         MMX(i-1,k) + TSC(_TMI,k);
	  sc = MAX(sc, IMX(i-1,k) + TSC(_TII,k));
	  
	  IMX(i,k) = MAX(sc + ISC(k), -INFTY);
	  
	  /* delete state */
	  sc =         MMX(i,k-1) + TSC(_TMD,k-1);
	  sc = MAX(sc, DMX(i,k-1) + TSC(_TDD,k-1));
	  
	  DMX(i,k) = MAX(sc, -INFTY);
	}
      
      /* match state */
      sc =         MMX(i-1,M-1)  + TSC(_TMM,M-1);
      sc = MAX(sc, IMX(i-1,M-1)  + TSC(_TIM,M-1));
      sc = MAX(sc, DMX(i-1,M-1)  + TSC(_TDM,M-1));
      sc = MAX(sc, XMX(i-1,_XMB) + BSC(M));
      
      MMX(i,M) = MAX(sc + MSC(M), -INFTY);
      
      /* E state update */
      XMX(i,_XME) = MAX(XMX(i,_XME), MMX(i,M) + ESC(M));
      
      /* Now the special states. Order is important here.
       * remember, N, C and J emissions are zero score by definition.
       */
      
      /* J state */
      sc          =         XMX(i-1,_XMJ) + lom->xsc[_XTJ][_LOOP];
      sc          = MAX(sc, XMX(i,_XME)   + lom->xsc[_XTE][_LOOP]);
      XMX(i,_XMJ) = MAX(sc, -INFTY);
      
      /* C state */
      sc =                  XMX(i-1,_XMC) + lom->xsc[_XTC][_LOOP];
      sc =          MAX(sc, XMX(i,_XME)   + lom->xsc[_XTE][_MOVE]);
      XMX(i,_XMC) = MAX(sc, -INFTY);
      
      /* N state */
      XMX(i,_XMN) = MAX(XMX(i-1,_XMN) + lom->xsc[_XTN][_LOOP], -INFTY);
      
      /* B state */
      sc          =         XMX(i,_XMN) + lom->xsc[_XTN][_MOVE];
      sc          = MAX(sc, XMX(i,_XMJ) + lom->xsc[_XTJ][_MOVE]);
      XMX(i,_XMB) = MAX(sc, -INFTY);
    }
  
  /* T state (not stored) */
  final_score = XMX(L,_XMC) + lom->xsc[_XTC][_MOVE];
  printf("Viterbi: %d\n", final_score);
  
  if (ret_tr != NULL) 
    {
      ViterbiTrace(hmm, dsq, L, mx, &tr);
      *ret_tr = tr;
    }
  
  return Scorify(final_score);		/* the total Viterbi score. */
}



/* Function: Forward()
 *
 * Purpose:  The Forward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ilogsum(); this is a slow but robust approach.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
 float 
 Forward(unsigned char *dsq, int L, struct plan7_s *hmm, 
	 cust_dpmatrix_s **ret_mx)
{
  cust_dpmatrix_s *mx;
  const struct logodds_s *lom = hmm->lom;
  int const *tsc = lom->tsc;
  
  int M = hmm->M;
  int **dp;
  int *xmx;
  int final_score;
  int i, k;  
  
  init_ilogsum();
  
  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 
  mx = CreateDPMatrix(L, M, 0, 0);
  dp = mx->dp;
  xmx = mx->xmx;
  
  /* Initialization of the zero row.
   * Note that xmx[i][stN] = 0 by definition for all i,
   *    and xmx[i][stT] = xmx[i][stC], so neither stN nor stT need
   *    to be calculated in DP matrices.
   */
  XMX(0,_XMN) = 0;                                   /* S->N, p=1            */
  XMX(0,_XMB) = lom->xsc[_XTN][_MOVE];               /* S->N->B, no N-tail   */
  XMX(0,_XME) = XMX(0,_XMC) = XMX(0,_XMJ) = -INFTY;  /* need seq to get here */
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = -INFTY;         /* need seq to get here */

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = -INFTY for all eight transitions (no node 0)
   */
  for (i = 1; i <= L; i++)
    {
      int const * misc = lom->misc[dsq[i]];
      int sc;
      
      MMX(i,0) = IMX(i,0) = DMX(i,0) = -INFTY;
      XMX(i,_XME) = -INFTY;
      
      for (k = 1; k < M; k++)
        {
          sc =
	    ilogsum(ilogsum(MMX(i-1,k-1) + TSC(_TMM,k-1),
			    IMX(i-1,k-1) + TSC(_TIM,k-1)),
		    ilogsum(XMX(i-1,_XMB) + BSC(k),
			    DMX(i-1,k-1) + TSC(_TDM,k-1)));
	  
	  MMX(i,k) = MAX(sc + MSC(k), -INFTY);
	  
	  /* update E state */
	  XMX(i,_XME) = ilogsum(XMX(i,_XME), MMX(i,k) + ESC(k));
	  
	  sc = 
	    ilogsum(MMX(i-1,k) + TSC(_TMI,k),
		    IMX(i-1,k) + TSC(_TII,k));
	  
	  IMX(i,k) = MAX(sc + ISC(k), -INFTY);
	  
          sc =
	    ilogsum(MMX(i,k-1) + TSC(_TMD,k-1),
		    DMX(i,k-1) + TSC(_TDD,k-1));
         
	  DMX(i,k) = MAX(sc, -INFTY);
	}
      
      sc =
	ilogsum(ilogsum(MMX(i-1,M-1) + TSC(_TMM,M-1),
			IMX(i-1,M-1) + TSC(_TIM,M-1)),
		ilogsum(XMX(i-1,_XMB) + BSC(M),
			DMX(i-1,M-1) + TSC(_TDM,M-1)));
      
      MMX(i,M) = MAX(sc + MSC(M), -INFTY);
      
      /* update E state */
      XMX(i,_XME) = ilogsum(XMX(i,_XME), MMX(i,M) + ESC(M));
      
      /* Now the special states.
       * remember, N, C, and J emissions are zero score by definition
       */
      
      /* J state */
      XMX(i,_XMJ) = ilogsum(XMX(i-1,_XMJ) + lom->xsc[_XTJ][_LOOP],
			    XMX(i,_XME)   + lom->xsc[_XTE][_LOOP]);
      
      /* C state */
      XMX(i,_XMC) = ilogsum(XMX(i-1,_XMC) + lom->xsc[_XTC][_LOOP],
			    XMX(i,_XME)   + lom->xsc[_XTE][_MOVE]);

      /* N state */
      XMX(i,_XMN) = XMX(i-1,_XMN) + lom->xsc[_XTN][_LOOP];
      
      /* B state */
      XMX(i,_XMB) = ilogsum(XMX(i,_XMN) + lom->xsc[_XTN][_MOVE],
			    XMX(i,_XMJ) + lom->xsc[_XTJ][_MOVE]);
    }
  
  /* T state (not stored) */
  final_score = XMX(L,_XMC) + lom->xsc[_XTC][_MOVE];
  printf("Forward: %d\n", final_score);
  
  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeDPMatrix(mx);
  
  return Scorify(final_score);           /* the total Forward score. */
}



/* Function: Backward()
 * 
 * Purpose:  The Backward dynamic programming algorithm.
 *           The scaling issue is dealt with by working in log space
 *           and calling ilogsum(); this is a slow but robust approach.
 *           
 * Args:     dsq    - sequence in digitized form
 *           L      - length of dsq
 *           hmm    - the model
 *           ret_mx - RETURN: dp matrix; pass NULL if it's not wanted
 *           
 * Return:   log P(S|M)/P(S|R), as a bit score.
 */
float 
Backward(unsigned char *dsq, int L, struct plan7_s *hmm,	
	 cust_dpmatrix_s **ret_mx)
{
  cust_dpmatrix_s *mx;
  const struct logodds_s *lom = hmm->lom;
  int *tsc = lom->tsc;                 /* cannot be const because of hack */
  
  int M = hmm->M;
  int **dp;
  int *xmx;
  int final_score;
  int i, k;
  
  init_ilogsum();
  
  /* Allocate a DP matrix with 0..L rows, 0..M-1 columns.
   */ 
  mx = CreateDPMatrix(L+1, hmm->M, 0, 0);
  dp = mx->dp;
  xmx = mx->xmx;
  
  /* Initialization of the L row.
   * Note that xmx[i][stS] = xmx[i][stN] by definition for all i,
   *    so stS need not be calculated in backward DP matrices.
   */
  XMX(L,_XMC) = lom->xsc[_XTC][_MOVE];                /* C<-T */
  XMX(L,_XME) = XMX(L,_XMC) + lom->xsc[_XTE][_MOVE];  /* E<-C, no C-tail */
  XMX(L,_XMJ) = XMX(L,_XMB) = XMX(L,_XMN) = -INFTY;   /* need seq to get out from here */
  
  {
    int const * misc = lom->misc[dsq[L]];
    
    for (k = M; k >= 1; k--) 
      {
	MMX(L,k) = XMX(L,_XME) + ESC(k);                /* M<-E ...                      */
	MMX(L,k) += MSC(k);                             /* ... + emitted match symbol    */
	IMX(L,k) = DMX(L,k) = -INFTY;                   /* need seq to get out from here */
      }
  }
  
  /* Recursion. Done as a pull.
   * Note slightly wasteful boundary conditions:
   *    M_M precalculated, D_M set to -INFTY,
   *    D_1 wastefully calculated.
   * Scores for transitions to D_M also have to be hacked to -INFTY,
   * as Plan7Logoddsify does not do this for us (I think? - ihh).
   */
  TSC(_TDD,M-1) = TSC(_TMD,M-1) = -INFTY;  /* no D_M state -- HACK -- should be in Plan7Logoddsify */
  for (i = L-1; i >= 0; i--)
    {
      /* Do the special states first.
       * remember, C, N and J emissions are zero score by definition
       */
      XMX(i,_XMC) = XMX(i+1,_XMC) + lom->xsc[_XTC][_LOOP];
      
      XMX(i,_XMB) = -INFTY;
      /* The following section has been hacked to fit a bug in core_algorithms.c
       * The "correct" code is:
       * for (k = M; k >= 1; k--)
       *   XMX(i,_XMB) = ilogsum(XMX(i,_XMB), MMX(i+1,k) + BSC(k))
       *
       * The following code gives the same results as core_algorithms.c:
       */
      XMX(i,_XMB) = ilogsum(XMX(i,_XMB), MMX(i+1,M) + BSC(M-1));
      
      for (k = M-1; k >= 1; k--)
	XMX(i,_XMB) = ilogsum(XMX(i,_XMB), MMX(i+1,k) + BSC(k));
      
      XMX(i,_XMJ) = ilogsum(XMX(i,_XMB)   + lom->xsc[_XTJ][_MOVE],
			    XMX(i+1,_XMJ) + lom->xsc[_XTJ][_LOOP]);
      
      XMX(i,_XME) = ilogsum(XMX(i,_XMC) + lom->xsc[_XTE][_MOVE],
			    XMX(i,_XMJ) + lom->xsc[_XTE][_LOOP]);
      
      XMX(i,_XMN) = ilogsum(XMX(i,_XMB)   + lom->xsc[_XTN][_MOVE],
			    XMX(i+1,_XMN) + lom->xsc[_XTN][_LOOP]);
      
      /* Now the main states. Note the boundary conditions at M.
       */
      if (i > 0) 
	{
	  int const * misc = lom->misc[dsq[i]];
	  int sc;
	  
	  MMX(i,M) = XMX(i,_XME) + ESC(M) + MSC(M);
	  DMX(i,M) = -INFTY;
	  
	  for (k = M-1; k >= 1; k--)
	    {
	      sc = ilogsum(ilogsum(XMX(i,_XME)  + ESC(k),
				   MMX(i+1,k+1) + TSC(_TMM,k)),
			   ilogsum(IMX(i+1,k) + TSC(_TMI,k),
				   DMX(i,k+1) + TSC(_TMD,k)));
	      MMX(i,k) = sc + MSC(k);
	      
	      sc = ilogsum(IMX(i+1,k)   + TSC(_TII,k),
			   MMX(i+1,k+1) + TSC(_TIM,k));
	      
	      IMX(i,k) = sc + ISC(k);
	      
	      sc = ilogsum(DMX(i,k+1)   + TSC(_TDD,k),
			   MMX(i+1,k+1) + TSC(_TDM,k));
	      
	      DMX(i,k) = sc;
	    }
	}
    }
  
  final_score = XMX(0,_XMN);
  printf("Backward: %d\n", final_score);
  
  if (ret_mx != NULL) *ret_mx = mx;
  else                FreeDPMatrix(mx);
  
  return Scorify(final_score);		/* the total Backward score. */
}


/***************************************************************************/

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
static void  ViterbiTrace(struct plan7_s *hmm, unsigned char *dsq, int N,
			  cust_dpmatrix_s *mx, struct p7trace_s **ret_tr)
{
  struct p7trace_s *tr;
  int curralloc;                /* current allocated length of trace */
  int tpos;                     /* position in trace */
  int i;                        /* position in seq (1..N) */
  int k;                        /* position in model (1..M) */
  int **dp;
  int *xmx;
  int sc;                       /* temp var for pre-emission score */
  
  /* Overallocate for the trace.
   * S-N-B- ... - E-C-T  : 6 states + N is minimum trace;
   * add N more as buffer.             
   */
  curralloc = N * 2 + 6; 
  P7AllocTrace(curralloc, &tr);
  
  dp  = mx->dp;
  xmx = mx->xmx;
  
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
  i    = N;        /* current i (seq pos) we're trying to assign */

  /* Traceback
   */
  while (tr->statetype[tpos-1] != STS) {
    switch (tr->statetype[tpos-1]) {
    case STM:                   /* M connects from i-1,k-1, or B */
      sc = MMX(i+1,k+1) - hmm->msc[dsq[i+1]][k+1];
      if (sc <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (sc == XMX(i,_XMB) + hmm->bsc[k+1])
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
                  {                             /* grow trace if necessary  */
                    curralloc += N;
                    P7ReallocTrace(tr, curralloc);
                  }
              }

          tr->statetype[tpos] = STB;
          tr->nodeidx[tpos]   = 0;
          tr->pos[tpos]       = 0;
        }
      else if (sc == MMX(i,k) + hmm->tsc[TMM][k])
        {
          tr->statetype[tpos] = STM;
          tr->nodeidx[tpos]   = k--;
          tr->pos[tpos]       = i--;
        }
      else if (sc == IMX(i,k) + hmm->tsc[TIM][k])
        {
          tr->statetype[tpos] = STI;
          tr->nodeidx[tpos]   = k;
          tr->pos[tpos]       = i--;
        }
      else if (sc == DMX(i,k) + hmm->tsc[TDM][k])
        {
          tr->statetype[tpos] = STD;
          tr->nodeidx[tpos]   = k--;
          tr->pos[tpos]       = 0;
        }
      else
        Die("traceback failed");
      break;

    case STD:                   /* D connects from M,D */
      if (DMX(i,k+1) <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (DMX(i,k+1) == MMX(i,k) + hmm->tsc[TMD][k])
        {
          tr->statetype[tpos] = STM;
          tr->nodeidx[tpos]   = k--;
          tr->pos[tpos]       = i--;
        }
      else if (DMX(i,k+1) == DMX(i,k) + hmm->tsc[TDD][k]) 
        {
          tr->statetype[tpos] = STD;
          tr->nodeidx[tpos]   = k--;
          tr->pos[tpos]       = 0;
        }
      else Die("traceback failed");
      break;

    case STI:                   /* I connects from M,I */
      sc = IMX(i+1,k) - hmm->isc[dsq[i+1]][k];
      if (sc <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (sc == MMX(i,k) + hmm->tsc[TMI][k])
        {
          tr->statetype[tpos] = STM;
          tr->nodeidx[tpos]   = k--;
          tr->pos[tpos]       = i--;
        }
      else if (sc == IMX(i,k) + hmm->tsc[TII][k])
        {
          tr->statetype[tpos] = STI;
          tr->nodeidx[tpos]   = k;
          tr->pos[tpos]       = i--;
        }
      else Die("traceback failed");
      break;

    case STN:                   /* N connects from S, N */
      if (i == 0 && XMX(i,_XMN) == 0)
        {
          tr->statetype[tpos] = STS;
          tr->nodeidx[tpos]   = 0;
          tr->pos[tpos]       = 0;
        }
      else if (i > 0 && XMX(i+1,_XMN) == XMX(i,_XMN) + hmm->xsc[XTN][LOOP])
        {
          tr->statetype[tpos] = STN;
          tr->nodeidx[tpos]   = 0;
          tr->pos[tpos]       = 0;    /* note convention adherence:  */
          tr->pos[tpos-1]     = i--;  /* first N doesn't emit        */
        }
      else Die("traceback failed");
      break;

    case STB:                   /* B connects from N, J */
      if (XMX(i,_XMB) <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (XMX(i,_XMB) == XMX(i,_XMN) + hmm->xsc[XTN][MOVE])
        {
          tr->statetype[tpos] = STN;
          tr->nodeidx[tpos]   = 0;
          tr->pos[tpos]       = 0;
        }
      else if (XMX(i,_XMB) == XMX(i,_XMJ) + hmm->xsc[XTJ][MOVE])
        {
          tr->statetype[tpos] = STJ;
          tr->nodeidx[tpos]   = 0;
          tr->pos[tpos]       = 0;
        }

      else Die("traceback failed");
      break;

    case STE:                   /* E connects from any M state. k set here */
      if (XMX(i,_XME) <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      for (k = hmm->M; k >= 1; k--)
        if (XMX(i,_XME) == MMX(i,k) + hmm->esc[k])
          {
            /* check for wing unfolding */
            if (Prob2Score(hmm->end[k], 1.) + 1*INTSCALE <=  hmm->esc[k])
              {
                int dk;         /* need a tmp k while moving thru delete wing */
                for (dk = hmm->M; dk > k; dk--)
                  {
                    tr->statetype[tpos] = STD;
                    tr->nodeidx[tpos]   = dk;
                    tr->pos[tpos]       = 0;
                    tpos++;
                    if (tpos == curralloc) 
                      {                         /* grow trace if necessary  */
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

    case STC:                   /* C comes from C, E */
      if (XMX(i,_XMC) <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (XMX(i,_XMC) == XMX(i-1,_XMC) + hmm->xsc[XTC][LOOP])
        {
          tr->statetype[tpos] = STC;
          tr->nodeidx[tpos]   = 0;
          tr->pos[tpos]       = 0;    /* note convention adherence: */
          tr->pos[tpos-1]     = i--;  /* first C doesn't emit       */
        }
      else if (XMX(i,_XMC) == XMX(i,_XME) + hmm->xsc[XTE][MOVE])
        {
          tr->statetype[tpos] = STE;
          tr->nodeidx[tpos]   = 0;
          tr->pos[tpos]       = 0; /* E is a nonemitter */
        }
      
      else Die("Traceback failed.");
      break;

    case STJ:                   /* J connects from E, J */
      if (XMX(i,_XMJ) <= -INFTY) { P7FreeTrace(tr); *ret_tr = NULL; return; }
      else if (XMX(i,_XMJ) == XMX(i-1,_XMJ) + hmm->xsc[XTJ][LOOP])
        {
          tr->statetype[tpos] = STJ;
          tr->nodeidx[tpos]   = 0;
          tr->pos[tpos]       = 0;    /* note convention adherence: */
          tr->pos[tpos-1]     = i--;  /* first J doesn't emit       */
        }
      else if (XMX(i,_XMJ) == XMX(i,_XME) + hmm->xsc[XTE][LOOP])
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
      {                         /* grow trace if necessary  */
        curralloc += N;
        P7ReallocTrace(tr, curralloc);
      }

  } /* end traceback, at S state; tpos == tlen now */
  tr->tlen = tpos;
  P7ReverseTrace(tr);
  *ret_tr = tr;
}
