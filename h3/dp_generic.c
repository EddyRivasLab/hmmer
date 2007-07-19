/* HMMER's standard implementation of the DP algorithms.
 * 
 * This implementation is modified from an optimized implementation
 * contributed by Jeremy D. Buhler (Washington University in
 * St. Louis).
 * 
 * Relative to the implementation in HMMER2, Jeremy rearranged data
 * structures to reduce the number of registers needed in the inner
 * loop; eliminated branches from the inner loop by unrolling the Mth
 * iteration in Viterbi and by replacing a bunch of "if" tests with
 * MAX; and exposed opportunities for hoisting and strength reduction
 * to the compiler. (The preceding sentence is nearly verbatim from
 * Jeremy's notes.) I then uplifted the JB code to H3, most notably 
 * involving conversion from H2's scaled integers to H3's floating
 * point calculations.
 * 
 * Contents:
 *     1. HMM alignment algorithms.
 *     2. Traceback algorithms. 
 *     3. Benchmark driver.
 *     4. Unit tests.
 *     5. Test driver.
 * 
 * SRE, Tue Jan 30 10:49:43 2007 [at Einstein's in St. Louis]
 * SVN $Id$
 */

#include "p7_config.h"

#include <easel.h>
#include <esl_alphabet.h>

#include "hmmer.h"


/*****************************************************************
 * 1. HMM alignment algorithms
 *****************************************************************/

#define MMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_M])
#define IMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_I])
#define DMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_D])
#define XMX(i,s) (xmx[(i) * p7G_NXCELLS + (s)])

#define TSC(s,k) (tsc[(k) * p7P_NTRANS + (s)])
#define MSC(k)   (rsc[(k) * p7P_NR     + p7P_MSC])
#define ISC(k)   (rsc[(k) * p7P_NR     + p7P_ISC])


/* Function:  p7_GViterbi()
 * Synopsis:  The Viterbi algorithm.
 * Incept:    SRE, Tue Jan 30 10:50:53 2007 [Einstein's, St. Louis]
 * 
 * Purpose:   The standard Viterbi dynamic programming algorithm. 
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <gm>, and DP matrix <gx> allocated for at least <L>
 *            by <gm->M> cells; calculate the maximum scoring path by
 *            Viterbi; return the Viterbi score in <ret_sc>, and the
 *            Viterbi matrix is in <gx>.
 *            
 *            The caller may then retrieve the Viterbi path by calling
 *            <p7_GTrace()>.
 *           
 *            The Viterbi lod score is returned in nats. The caller
 *            needs to subtract a null model lod score, then convert
 *            to bits.
 *           
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm     - profile. Does not need to contain any
 *                     reference pointers (alphabet, HMM, or null model)
 *            gx     - DP matrix with room for an MxL alignment
 *            ret_sc - RETURN: Viterbi lod score in nats
 *           
 * Return:   <eslOK> on success.
 */
int
p7_GViterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx, float *ret_sc)
{
  float const *tsc  = gm->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          M    = gm->M;
  int          i,k;
  float        esc  = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;

  /* Initialization of the zero row.  */
  XMX(0,p7G_N) = 0;                                           /* S->N, p=1            */
  XMX(0,p7G_B) = gm->xsc[p7P_N][p7P_MOVE];                    /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;  /* need seq to get here */
  for (k = 0; k <= gm->M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = -eslINFINITY;            /* need seq to get here */


  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = impossible for all eight transitions (no node 0)
   *    I_M is wastefully calculated (doesn't exist)
   *    D_M is wastefully calculated (provably can't appear in a Viterbi path)
   */
  for (i = 1; i <= L; i++) 
    {
      float const *rsc = gm->rsc[dsq[i]];
      float sc;

      MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;
      DMX(i,M) = -eslINFINITY;	/* for Viterbi (not Forward) paths thru D_M can't win */
      XMX(i,p7G_E) = -eslINFINITY;
    
      for (k = 1; k < gm->M; k++) 
	{
  	  /* match state */
	  sc       = ESL_MAX(    MMX(i-1,k-1)   + TSC(p7P_MM,k-1), 
				 IMX(i-1,k-1)   + TSC(p7P_IM,k-1));
	  sc       = ESL_MAX(sc, DMX(i-1,k-1)   + TSC(p7P_DM,k-1));
	  sc       = ESL_MAX(sc, XMX(i-1,p7G_B) + TSC(p7P_BM,k-1));
	  MMX(i,k) = sc + MSC(k);

	  /* E state update */
	  XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), MMX(i,k) + esc);
	  /* in Viterbi alignments, Dk->E can't win in local mode (and
	   * isn't possible in glocal mode), so don't bother
	   * looking. */

	  /* insert state */
	  sc = ESL_MAX(MMX(i-1,k) + TSC(p7P_MI,k),
		       IMX(i-1,k) + TSC(p7P_II,k));
	  IMX(i,k) = sc + ISC(k);
	  
	  /* delete state */
	  DMX(i,k) =  ESL_MAX(MMX(i,k-1) + TSC(p7P_MD,k-1),
			      DMX(i,k-1) + TSC(p7P_DD,k-1));
	}

      /* Unrolled match state M. */
      sc       = ESL_MAX(    MMX(i-1,M-1)   + TSC(p7P_MM,M-1),
			     IMX(i-1,M-1)   + TSC(p7P_IM,M-1));
      sc       = ESL_MAX(sc, DMX(i-1,M-1 )  + TSC(p7P_DM,M-1));
      sc       = ESL_MAX(sc, XMX(i-1,p7G_B) + TSC(p7P_BM,M-1));
      MMX(i,M) = sc + MSC(M);
      
      /* Unrolled delete state D_M 
       * (Unlike internal Dk->E transitions that can never appear in 
       * Viterbi alignments, D_M->E is possible in glocal mode.)
       */
      DMX(i,M) = ESL_MAX(MMX(i,M-1) + TSC(p7P_MD,M-1),
			 DMX(i,M-1) + TSC(p7P_DD,M-1));

      /* E state update; transition from M_M scores 0 by def'n */
      sc  =          ESL_MAX(XMX(i,p7G_E), MMX(i,M));
      XMX(i,p7G_E) = ESL_MAX(sc,           DMX(i,M));
   
      /* Now the special states. Order is important here.
       * remember, N, C and J emissions are zero score by definition.
       */
      /* J state */
      sc           =             XMX(i-1,p7G_J) + gm->xsc[p7P_J][p7P_LOOP];   /* J->J */
      XMX(i,p7G_J) = ESL_MAX(sc, XMX(i,  p7G_E) + gm->xsc[p7P_E][p7P_LOOP]);  /* E->J is E's "loop" */
      
      /* C state */
      sc           =             XMX(i-1,p7G_C) + gm->xsc[p7P_C][p7P_LOOP];
      XMX(i,p7G_C) = ESL_MAX(sc, XMX(i,  p7G_E) + gm->xsc[p7P_E][p7P_MOVE]);
      
      /* N state */
      XMX(i,p7G_N) = XMX(i-1,p7G_N) + gm->xsc[p7P_N][p7P_LOOP];
      
      /* B state */
      sc           =             XMX(i,p7G_N) + gm->xsc[p7P_N][p7P_MOVE];   /* N->B is N's move */
      XMX(i,p7G_B) = ESL_MAX(sc, XMX(i,p7G_J) + gm->xsc[p7P_J][p7P_MOVE]);  /* J->B is J's move */
    }
  
  /* T state (not stored) */
  *ret_sc = XMX(L,p7G_C) + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}




/* Function:  p7_GForward()
 * Synopsis:  The Forward algorithm.
 * Incept:    SRE, Mon Apr 16 13:57:35 2007 [Janelia]
 *
 * Purpose:   The Forward dynamic programming algorithm. 
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <gm>, and DP matrix <gm> allocated for at least <gm->M>
 *            by <L> cells; calculate the probability of the sequence
 *            given the model using the Forward algorithm; return the
 *            Forward matrix in <mx>, and the Forward score in <ret_sc>.
 *           
 *            The Forward score is in lod score form.  To convert to a
 *            bitscore, the caller needs to subtract a null model lod
 *            score, then convert to bits.
 *           
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm     - profile. Does not need to contain any
 *                     reference pointers (alphabet, HMM, or null model)
 *            mx     - DP matrix with room for an MxL alignment
 *            ret_sc - RETURN: Forward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
p7_GForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx, float *ret_sc)
{
  float const *tsc  = gm->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx; 			    
  int          M    = gm->M;
  int          i, k;  
  float        esc  = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;

  /* Initialization of the zero row, and the lookup table of the log
   * sum routine.
   */
  XMX(0,p7G_N) = 0;                                           /* S->N, p=1            */
  XMX(0,p7G_B) = gm->xsc[p7P_N][p7P_MOVE];                    /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;  /* need seq to get here */
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = -eslINFINITY;            /* need seq to get here */
  p7_FLogsumInit();

  /* Recursion. Done as a pull.
   * Note some slightly wasteful boundary conditions:  
   *    tsc[0] = impossible for all eight transitions (no node 0)
   *    D_1 is wastefully calculated (doesn't exist)
   */
  for (i = 1; i <= L; i++) 
    {
      float const *rsc = gm->rsc[dsq[i]];
      float sc;

      MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;
      XMX(i, p7G_E) = -eslINFINITY;

      for (k = 1; k < M; k++)
	{
	  /* match state */
	  sc = p7_FLogsum(p7_FLogsum(MMX(i-1,k-1)   + TSC(p7P_MM,k-1), 
				     IMX(i-1,k-1)   + TSC(p7P_IM,k-1)),
			  p7_FLogsum(XMX(i-1,p7G_B) + TSC(p7P_BM,k-1),
				     DMX(i-1,k-1)   + TSC(p7P_DM,k-1)));
	  MMX(i,k) = sc + MSC(k);

	  /* insert state */
	  sc = p7_FLogsum(MMX(i-1,k) + TSC(p7P_MI,k),
			  IMX(i-1,k) + TSC(p7P_II,k));
	  IMX(i,k) = sc + ISC(k);

	  /* delete state */
	  DMX(i,k) = p7_FLogsum(MMX(i,k-1) + TSC(p7P_MD,k-1),
				DMX(i,k-1) + TSC(p7P_DD,k-1));

	  /* E state update */
	  XMX(i,p7G_E) = p7_FLogsum(p7_FLogsum(MMX(i,k) + esc,
					       DMX(i,k) + esc),
				               XMX(i,p7G_E));
	}
      /* unrolled match state M_M */
      sc = p7_FLogsum(p7_FLogsum(MMX(i-1,M-1)   + TSC(p7P_MM,M-1), 
				 IMX(i-1,M-1)   + TSC(p7P_IM,M-1)),
		      p7_FLogsum(XMX(i-1,p7G_B) + TSC(p7P_BM,M-1),
				 DMX(i-1,M-1)   + TSC(p7P_DM,M-1)));
      MMX(i,M) = sc + MSC(M);

      /* unrolled delete state D_M */
      DMX(i,M) = p7_FLogsum(MMX(i,M-1) + TSC(p7P_MD,M-1),
			    DMX(i,M-1) + TSC(p7P_DD,M-1));

      /* unrolled E state update */
      XMX(i,p7G_E) = p7_FLogsum(p7_FLogsum(MMX(i,M),
					   DMX(i,M)),
					   XMX(i,p7G_E));

      /* J state */
      XMX(i,p7G_J) = p7_FLogsum(XMX(i-1,p7G_J) + gm->xsc[p7P_J][p7P_LOOP],
				XMX(i,  p7G_E) + gm->xsc[p7P_E][p7P_LOOP]);
      /* C state */
      XMX(i,p7G_C) = p7_FLogsum(XMX(i-1,p7G_C) + gm->xsc[p7P_C][p7P_LOOP],
				XMX(i,  p7G_E) + gm->xsc[p7P_E][p7P_MOVE]);
      /* N state */
      XMX(i,p7G_N) = XMX(i-1,p7G_N) + gm->xsc[p7P_N][p7P_LOOP];

      /* B state */
      XMX(i,p7G_B) = p7_FLogsum(XMX(i,  p7G_N) + gm->xsc[p7P_N][p7P_MOVE],
				XMX(i,  p7G_J) + gm->xsc[p7P_J][p7P_MOVE]);
    }

  *ret_sc = XMX(L,p7G_C) + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}


/* Function:  p7_GHybrid()
 * Synopsis:  The "hybrid" algorithm.
 * Incept:    SRE, Sat May 19 10:01:46 2007 [Janelia]
 *
 * Purpose:   The profile HMM version of the Hwa "hybrid" alignment
 *            algorithm \citep{YuHwa02}. The "hybrid" score is the
 *            maximum score in the Forward matrix. 
 *            
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <gm>, and DP matrix <mx> allocated for at least <gm->M>
 *            by <L> cells; calculate the probability of the sequence
 *            given the model using the Forward algorithm; return
 *            the calculated Forward matrix in <mx>, and optionally
 *            return the Forward score in <opt_fwdscore> and/or the
 *            Hybrid score in <opt_hybscore>.
 *           
 *            This is implemented as a wrapper around <p7_GForward()>.
 *            The Forward matrix and the Forward score obtained from
 *            this routine are identical to what <p7_GForward()> would
 *            return.
 *           
 *            The scores are returned in lod form.  To convert to a
 *            bitscore, the caller needs to subtract a null model lod
 *            score, then convert to bits.
 *           
 * Args:      dsq          - sequence in digitized form, 1..L
 *            L            - length of dsq
 *            gm           - profile. Does not need to contain any
 *                           reference pointers (alphabet, HMM, or null model)
 *            gx           - DP matrix with room for an MxL alignment
 *            opt_fwdscore - optRETURN: Forward lod score in nats.
 *            opt_hybscore - optRETURN: Hybrid lod score in nats.
 *
 * Returns:   <eslOK> on success, and results are in <mx>, <opt_fwdscore>,
 *            and <opt_hybscore>.
 */
int
p7_GHybrid(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx, float *opt_fwdscore, float *opt_hybscore)
{
  float   F    = -eslINFINITY;
  float   H    = -eslINFINITY;
  float **dp   = gx->dp;
  int     i,k;
  int     status;

  if ((status = p7_GForward(dsq, L, gm, gx, &F)) != eslOK)  goto ERROR;
  for (i = 1; i <= L; i++)
    for (k = 1 ; k <= gm->M; k++)
      H = ESL_MAX(H, MMX(i,k));
  
  if (opt_fwdscore != NULL) *opt_fwdscore = F;
  if (opt_hybscore != NULL) *opt_hybscore = H;
  return eslOK;

 ERROR:
  if (opt_fwdscore != NULL) *opt_fwdscore = 0;
  if (opt_hybscore != NULL) *opt_hybscore = 0;
  return status;
}

/*****************************************************************
 * 2. Traceback
 *****************************************************************/

/* Function: p7_GTrace()
 * Incept:   SRE, Thu Feb  1 10:25:56 2007 [UA 8018 St. Louis to Dulles]
 * 
 * Purpose:  Traceback of a Viterbi matrix: retrieval 
 *           of optimum alignment.
 *           
 *           This function is currently implemented as a
 *           reconstruction traceback, rather than using a shadow
 *           matrix. Because H3 uses floating point scores, and we
 *           can't compare floats for equality, we have to compare
 *           floats for near-equality and therefore, formally, we can
 *           only guarantee a near-optimal traceback. However, the
 *           score difference from true optimal will be negligible.
 *           
 * Args:     dsq    - sequence aligned to (digital form) 1..L 
 *           L      - length of dsq gm    
 *           gm     - profile model; does not need any ref ptrs
 *           mx     - the matrix to trace, L x M
 *           tr     - storage for the recovered traceback.
 *           
 * Return:   <eslOK> on success.
 *           <eslFAIL> if even the optimal path has zero probability;
 *           in this case, the trace is set blank (<tr->N = 0>).
 */
int
p7_GTrace(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr)
{
  int     status;
  int     i;			/* position in seq (1..L) */
  int     k;			/* position in model (1..M) */
  int     M   = gm->M;
  float   sc;
  float **dp  = gx->dp;
  float  *xmx = gx->xmx;
  float   tol = 1e-5;
  float   esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
  float const *tsc  = gm->tsc;


  if ((status = p7_trace_Reuse(tr)) != eslOK) goto ERROR;

  /* Initialization.
   * (back to front. ReverseTrace() called later.)
   */
  if ((status = p7_trace_Append(tr, p7T_T, 0, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_Append(tr, p7T_C, 0, 0)) != eslOK) goto ERROR;
  i    = L;			/* next position to explain in seq */

  /* Traceback
   */
  while (tr->st[tr->N-1] != p7T_S) {
    float const *rsc = gm->rsc[dsq[i]];

    switch (tr->st[tr->N-1]) {
    case p7T_C:		/* C(i) comes from C(i-1) or E(i) */
      if   (XMX(i,p7G_C) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible C reached at i=%d", i);

      if (esl_FCompare(XMX(i, p7G_C), XMX(i-1, p7G_C) + gm->xsc[p7P_C][p7P_LOOP], tol) == eslOK) {
	tr->i[tr->N-1]    = i--;  /* first C doesn't emit: subsequent ones do */
	status = p7_trace_Append(tr, p7T_C, 0, 0);
      } else if (esl_FCompare(XMX(i, p7G_C), XMX(i, p7G_E) + gm->xsc[p7P_E][p7P_MOVE], tol) == eslOK) 
	status = p7_trace_Append(tr, p7T_E, 0, 0);
      else ESL_XEXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:		/* E connects from any M state. k set here */
      if (XMX(i, p7G_E) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if (esl_FCompare(XMX(i, p7G_E), MMX(i,M), tol) == eslOK) { k = M; status = p7_trace_Append(tr, p7T_M, M, i); }
      else {
	for (k = M-1; k >= 1; k--)
	  if (esl_FCompare(XMX(i, p7G_E), MMX(i,k) + esc, tol) == eslOK)
	    { status = p7_trace_Append(tr, p7T_M, k, i); break; }
	if (k < 0) ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      }
      break;

    case p7T_M:			/* M connects from i-1,k-1, or B */
      sc = MMX(i,k) - MSC(k);
      if (sc == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      if      (esl_FCompare(sc, XMX(i-1,p7G_B) + TSC(p7P_BM, k-1), tol) == eslOK) status = p7_trace_Append(tr, p7T_B, 0,   0);
      else if (esl_FCompare(sc, MMX(i-1,k-1)   + TSC(p7P_MM, k-1), tol) == eslOK) status = p7_trace_Append(tr, p7T_M, k-1, i-1);
      else if (esl_FCompare(sc, IMX(i-1,k-1)   + TSC(p7P_IM, k-1), tol) == eslOK) status = p7_trace_Append(tr, p7T_I, k-1, i-1);
      else if (esl_FCompare(sc, DMX(i-1,k-1)   + TSC(p7P_DM, k-1), tol) == eslOK) status = p7_trace_Append(tr, p7T_D, k-1, 0);
      else ESL_XEXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
      if (status != eslOK) goto ERROR;
      k--; 
      i--;
      break;

    case p7T_D:			/* D connects from M,D at i,k-1 */
      if (DMX(i, k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);

      if      (esl_FCompare(DMX(i,k), MMX(i, k-1) + TSC(p7P_MD, k-1), tol) == eslOK) status = p7_trace_Append(tr, p7T_M, k-1, i);
      else if (esl_FCompare(DMX(i,k), DMX(i, k-1) + TSC(p7P_DD, k-1), tol) == eslOK) status = p7_trace_Append(tr, p7T_D, k-1, 0);
      else ESL_XEXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      if (status != eslOK) goto ERROR;
      k--;
      break;

    case p7T_I:			/* I connects from M,I at i-1,k*/
      sc = IMX(i,k) - ISC(k);
      if (sc == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      if      (esl_FCompare(sc, MMX(i-1,k) + TSC(p7P_MI, k), tol) == eslOK) status = p7_trace_Append(tr, p7T_M, k, i-1);
      else if (esl_FCompare(sc, IMX(i-1,k) + TSC(p7P_II, k), tol) == eslOK) status = p7_trace_Append(tr, p7T_I, k, i-1);
      else ESL_XEXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      if (status != eslOK) goto ERROR;
      i--;
      break;

    case p7T_N:			/* N connects from S, N */
      if (XMX(i, p7G_N) <= p7_IMPOSSIBLE) ESL_XEXCEPTION(eslFAIL, "impossible N reached at i=%d", i);

      if      (i == 0 && esl_FCompare(XMX(i, p7G_N), 0.0f, tol) == eslOK)
	status = p7_trace_Append(tr, p7T_S, 0, 0);
      else if (i > 0  && esl_FCompare(XMX(i,p7G_N), XMX(i-1, p7G_N) + gm->xsc[p7P_N][p7P_LOOP], tol) == eslOK) {
	tr->i[tr->N-1] = i--;
	status = p7_trace_Append(tr, p7T_N, 0, 0);
      } else ESL_XEXCEPTION(eslFAIL, "N at i=%d couldn't be traced", i);
      if (status != eslOK) goto ERROR;
      break;

    case p7T_B:			/* B connects from N, J */
      if (XMX(i,p7G_B) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if (esl_FCompare(XMX(i,p7G_B), XMX(i, p7G_N) + gm->xsc[p7P_N][p7P_MOVE], tol)  == eslOK)
	status = p7_trace_Append(tr, p7T_N, 0, 0);
      else if (esl_FCompare(XMX(i,p7G_B),  XMX(i, p7G_J) + gm->xsc[p7P_J][p7P_MOVE], tol) == eslOK)
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      else  ESL_XEXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:			/* J connects from E(i) or J(i-1) */
      if (XMX(i,p7G_J) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible J reached at i=%d", i);

      if (esl_FCompare(XMX(i,p7G_J), XMX(i-1,p7G_J) + gm->xsc[p7P_J][p7P_LOOP], tol) == eslOK) {
	tr->i[tr->N-1] = i--;
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      } else if (esl_FCompare(XMX(i,p7G_J), XMX(i,p7G_E) + gm->xsc[p7P_E][p7P_LOOP], tol) == eslOK) 
	status = p7_trace_Append(tr, p7T_E, 0, 0);
      else  ESL_XEXCEPTION(eslFAIL, "J at i=%d couldn't be traced", i);
      break;

    default: ESL_XEXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    if (status != eslOK) goto ERROR;
  } /* end traceback, at S state */

  if ((status = p7_trace_Reverse(tr)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  return status;
}




/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7DP_GENERIC_BENCHMARK
/* gcc -o benchmark-generic -g -O2 -I. -L. -I../easel -L../easel -Dp7DP_GENERIC_BENCHMARK dp_generic.c -lhmmer -leasel -lm
 * ./benchmark-generic <hmmfile>
 */
/* As of Tue Jul 17 13:14:43 2007
 * 61 Mc/s for Viterbi, 8.6 Mc/s for Forward.
 * (gcc -g -O2, 3.2GHz Xeon, N=50K, L=400, M=72 RRM_1 model)
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-f",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "use the Forward algorithm",                      0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose: show individual scores",             0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the generic implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = NULL;
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *gx      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  float           nullsc;
  float           bitscore;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)     != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);

  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL);

  gx = p7_gmx_Create(gm->M, L);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rnd_xfIID(r, bg->f, abc->K, L, dsq);

      if (esl_opt_GetBoolean(go, "-f")) p7_GForward(dsq, L, gm, gx, &sc);
      else                              p7_GViterbi(dsq, L, gm, gx, &sc);

      p7_bg_NullOne(bg, dsq, L, &nullsc);
      bitscore = (sc - nullsc) / eslCONST_LOG2;
      
      if (esl_opt_GetBoolean(go, "-v")) printf("%.4f bits  (%.4f raw)\n", bitscore, sc); 
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  free(dsq);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7DP_GENERIC_BENCHMARK*/


/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7DP_GENERIC_TESTDRIVE
#include <string.h>

#include "esl_getopts.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

/* The "basic" utest is a minimal driver for making a small DNA profile and a small DNA sequence,
 * then running Viterbi and Forward. It's useful for dumping DP matrices and profiles for debugging.
 */
static void
utest_basic(ESL_GETOPTS *go)
{
  char           *query= "# STOCKHOLM 1.0\n\nseq1 GAATTC\nseq2 GAATTC\n//\n";
  int             fmt  = eslMSAFILE_STOCKHOLM;
  char           *targ = "GAATTC";
  ESL_ALPHABET   *abc  = NULL;
  ESL_MSA        *msa  = NULL;
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_BG          *bg   = NULL;
  P7_DPRIOR      *pri  = NULL;	
  ESL_DSQ        *dsq  = NULL;
  P7_GMX         *gx   = NULL;
  P7_TRACE        *tr  = NULL;
  int             L    = strlen(targ);
  float           vsc, vsc2, fsc;

  if ((abc = esl_alphabet_Create(eslDNA))          == NULL)  esl_fatal("failed to create alphabet");
  if ((pri = p7_dprior_CreateNucleic())            == NULL)  esl_fatal("failed to create prior");
  if ((msa = esl_msa_CreateFromString(query, fmt)) == NULL)  esl_fatal("failed to create MSA");
  if (esl_msa_Digitize(abc, msa)                   != eslOK) esl_fatal("failed to digitize MSA");
  if (p7_Fastmodelmaker(msa, 0.5, &hmm, NULL)      != eslOK) esl_fatal("failed to create GAATTC model");
  if (p7_ParameterEstimation(hmm, pri)             != eslOK) esl_fatal("failed to parameterize GAATTC model");
  if ((bg = p7_bg_Create(abc))                     == NULL)  esl_fatal("failed to create DNA null model");
  if ((gm = p7_profile_Create(hmm->M, abc))        == NULL)  esl_fatal("failed to create GAATTC profile");
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL)!= eslOK) esl_fatal("failed to config profile");
  if (p7_profile_Validate(gm, 0.0001)              != eslOK) esl_fatal("whoops, profile is bad!");
  if (esl_abc_CreateDsq(abc, targ, &dsq)           != eslOK) esl_fatal("failed to create GAATTC digital sequence");
  if ((gx = p7_gmx_Create(gm->M, L))               == NULL)  esl_fatal("failed to create DP matrix");
  if ((tr = p7_trace_Create(L))                    == NULL)  esl_fatal("trace creation failed");

  p7_GViterbi   (dsq, L, gm, gx, &vsc);
  if (esl_opt_GetBoolean(go, "-v")) printf("Viterbi score: %.4f\n", vsc);
  if (esl_opt_GetBoolean(go, "-v")) p7_gmx_Dump(stdout, gx);

  p7_GTrace     (dsq, L, gm, gx, tr);
  p7_trace_Score(tr, dsq, gm, &vsc2);
  if (esl_opt_GetBoolean(go, "-v")) p7_trace_Dump(stdout, tr, gm, dsq);
  
  if (esl_FCompare(vsc, vsc2, 1e-5) != eslOK)  esl_fatal("trace score and Viterbi score don't agree.");

  p7_GForward   (dsq, L, gm, gx, &fsc);
  if (esl_opt_GetBoolean(go, "-v")) printf("Forward score: %.4f\n", fsc);
  if (esl_opt_GetBoolean(go, "-v")) p7_gmx_Dump(stdout, gx);

  p7_trace_Destroy(tr);
  p7_gmx_Destroy(gx);
  free(dsq);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_msa_Destroy(msa);
  p7_dprior_Destroy(pri);
  esl_alphabet_Destroy(abc);
  return;
}

/* Viterbi validation is done by comparing the returned score
 * to the score of the optimal trace. Not foolproof, but catches
 * many kinds of errors.
 * 
 * Another check is that the average score should be <= 0,
 * since the random sequences are drawn from the null model.
 */ 
static void
utest_viterbi(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, P7_PROFILE *gm, int nseq, int L)
{
  float     avg_sc = 0.;
  char      errbuf[eslERRBUFSIZE];
  ESL_DSQ  *dsq = NULL;
  P7_GMX   *gx  = NULL;
  P7_TRACE *tr  = NULL;
  int       idx;
  float     sc1, sc2;

  if ((dsq    = malloc(sizeof(ESL_DSQ) *(L+2))) == NULL)  esl_fatal("malloc failed");
  if ((tr     = p7_trace_Create(L))             == NULL)  esl_fatal("trace creation failed");
  if ((gx     = p7_gmx_Create(gm->M, L))        == NULL)  esl_fatal("matrix creation failed");

  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rnd_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");
      if (p7_GViterbi(dsq, L, gm, gx, &sc1)       != eslOK) esl_fatal("viterbi failed");
      if (p7_GTrace  (dsq, L, gm, gx, tr)         != eslOK) esl_fatal("trace failed");
      if (p7_trace_Validate(tr, abc, dsq, errbuf) != eslOK) esl_fatal("trace invalid:\n%s", errbuf);
      if (p7_trace_Score(tr, dsq, gm, &sc2)       != eslOK) esl_fatal("trace score failed");
      if (sc1 != sc2) esl_fatal("Trace score not equal to Viterbi score");
      if (p7_bg_NullOne(bg, dsq, L, &sc2)         != eslOK) esl_fatal("null score failed");

      avg_sc += (sc1 - sc2);

      if (esl_opt_GetBoolean(go, "--vv"))       
	printf("utest_viterbi: Viterbi score: %.4f (null %.4f) (total so far: %.4f)\n", sc1, sc2, avg_sc);
    }

  avg_sc /= (float) nseq;
  if (avg_sc > 0.) esl_fatal("Viterbi scores have positive expectation (%f nats)", avg_sc);

  p7_gmx_Destroy(gx);
  p7_trace_Destroy(tr);
  free(dsq);
  return;
}


/* Forward is harder to validate. 
 * We do know that the Forward score is >= Viterbi.
 * We also know that the expected score on random seqs is <= 0 (not
 * exactly - we'd have to sample the random length from the background
 * model too, not just use a fixed L - but it's close enough to
 * being true to be a useful test.)
 */
static void
utest_forward(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, P7_PROFILE *gm, int nseq, int L)
{
  float     avg_sc;
  ESL_DSQ  *dsq = NULL;
  P7_GMX   *gx  = NULL;
  P7_TRACE *tr  = NULL;
  int       idx;
  float     vsc, fsc, nullsc;

  if ((dsq    = malloc(sizeof(ESL_DSQ) *(L+2))) == NULL)  esl_fatal("malloc failed");
  if ((gx     = p7_gmx_Create(gm->M, L))        == NULL)  esl_fatal("matrix creation failed");

  avg_sc = 0.;
  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rnd_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");
      if (p7_GViterbi(dsq, L, gm, gx, &vsc)       != eslOK) esl_fatal("viterbi failed");
      if (p7_GForward(dsq, L, gm, gx, &fsc)       != eslOK) esl_fatal("forward failed");
      if (fsc < vsc) esl_fatal("Forward score can't be less than Viterbi score");
      if (p7_bg_NullOne(bg, dsq, L, &nullsc)      != eslOK) esl_fatal("null score failed");

      avg_sc += fsc - nullsc;

      if (esl_opt_GetBoolean(go, "--vv")) 
	printf("utest_forward: Forward score: %.4f (total so far: %.4f)\n", fsc, avg_sc);
    }

  avg_sc /= (float) nseq;
  if (avg_sc > 0.) esl_fatal("Forward scores have positive expectation (%f nats)", avg_sc);

  p7_gmx_Destroy(gx);
  p7_trace_Destroy(tr);
  free(dsq);
  return;
}

/* The "generation" test scores sequences generated by the same profile.
 * Each Viterbi and Forward score should be >= the trace score of the emitted seq.
 * The expectation of Forward scores should be positive.
 */
static void
utest_generation(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc,
		 P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, int nseq)
{
  ESL_SQ   *sq = esl_sq_CreateDigital(abc);
  P7_GMX   *gx = p7_gmx_Create(gm->M, 100);
  P7_TRACE *tr = p7_trace_Create(256);
  float     vsc, fsc, nullsc, tracesc;
  float     avg_fsc;
  int       idx;

  avg_fsc = 0.0;
  for (idx = 0; idx < nseq; idx++)
    {
      if (p7_ProfileEmit(r, hmm, gm, bg, sq, tr)     != eslOK) esl_fatal("profile emission failed");

      if (p7_gmx_GrowTo(gx, gm->M, sq->n)            != eslOK) esl_fatal("failed to reallocate gmx");
      if (p7_GViterbi(sq->dsq, sq->n, gm, gx, &vsc)  != eslOK) esl_fatal("viterbi failed");
      if (p7_GForward(sq->dsq, sq->n, gm, gx, &fsc)  != eslOK) esl_fatal("forward failed");
      if (p7_trace_Score(tr, sq->dsq, gm, &tracesc)  != eslOK) esl_fatal("trace score failed");
      if (p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc) != eslOK) esl_fatal("null score failed");

      if (vsc < tracesc) esl_fatal("viterbi score is less than trace");
      if (fsc < tracesc) esl_fatal("forward score is less than trace");

      if (esl_opt_GetBoolean(go, "--vv")) 
	printf("generated:  len=%d v=%8.4f  f=%8.4f  t=%8.4f\n", sq->n, vsc, fsc, tracesc);
      
      avg_fsc += (fsc - nullsc);
    }
  
  avg_fsc /= (float) nseq;
  if (avg_fsc < 0.) esl_fatal("generation: Forward scores have negative expectation (%f nats)", avg_fsc);

  p7_gmx_Destroy(gx);
  p7_trace_Destroy(tr);
  esl_sq_Destroy(sq);
}

/* The "enumeration" test samples a random enumerable HMM (transitions to insert are 0,
 * so the generated seq space only includes seqs of L<=M). 
 *
 * The test scores all seqs of length <=M by both Viterbi and Forward, verifies that 
 * the two scores are identical, and verifies that the sum of all the probabilities is
 * 1.0. It also verifies that the score of a sequence of length M+1 is indeed -infinity.
 * 
 * Because this function is going to work in unscaled probabilities, adding them up,
 * all P(seq) terms must be >> DBL_EPSILON.  That means M must be small; on the order 
 * of <= 10. 
 */
static void
utest_enumeration(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, int M)
{
  char            errbuf[eslERRBUFSIZE];
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_BG          *bg   = NULL;
  ESL_DSQ        *dsq  = NULL;
  P7_GMX         *gx   = NULL;
  float  vsc, fsc;
  float  bg_ll;   		/* log P(seq | bg) */
  double vp, fp;		/* P(seq,\pi | model) and P(seq | model) */
  int L;
  int i;
  double total_p;
  char   *seq;
    
  /* Sample an enumerable HMM & profile of length M.
   */
  if (p7_hmm_SampleEnumerable(r, M, abc, &hmm)      != eslOK) esl_fatal("failed to sample an enumerable HMM");
  if ((bg = p7_bg_Create(abc))                      == NULL)  esl_fatal("failed to create null model");
  if ((gm = p7_profile_Create(hmm->M, abc))         == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig(hmm, bg, gm, 0, p7_UNILOCAL) != eslOK) esl_fatal("failed to config profile");
  if (p7_hmm_Validate    (hmm,     0.0001, errbuf)  != eslOK) esl_fatal("whoops, HMM is bad!");
  if (p7_profile_Validate(gm, 0.0001)               != eslOK) esl_fatal("whoops, profile is bad!");

  if (  (dsq = malloc(sizeof(ESL_DSQ) * (M+3)))     == NULL)  esl_fatal("allocation failed");
  if (  (seq = malloc(sizeof(char)    * (M+2)))     == NULL)  esl_fatal("allocation failed");
  if ((gx     = p7_gmx_Create(hmm->M, M+3))         == NULL)  esl_fatal("matrix creation failed");

  /* Enumerate all sequences of length L <= M
   */
  total_p = 0;
  for (L = 0; L <= M; L++)
    {
      /* Initialize dsq of length L at 0000... */
      dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;
      for (i = 1; i <= L; i++) dsq[i] = 0;

      while (1) 		/* enumeration of seqs of length L*/
	{
	  if (p7_GViterbi(dsq, L, gm, gx, &vsc)  != eslOK) esl_fatal("viterbi failed");
	  if (p7_GForward(dsq, L, gm, gx, &fsc)  != eslOK) esl_fatal("forward failed");
 
	  /* calculate bg log likelihood component of the scores */
	  for (bg_ll = 0., i = 1; i <= L; i++)  bg_ll += log(bg->f[dsq[i]]);
	  
	  /* convert to probabilities, adding the bg LL back to the LLR */
	  vp =  exp(vsc + bg_ll);
	  fp =  exp(fsc + bg_ll);

	  if (esl_opt_GetBoolean(go, "--vv")) {
	    esl_abc_Textize(abc, dsq, L, seq);
	    printf("probability of sequence: %10s   %16g  (lod v=%8.4f f=%8.4f)\n", seq, fp, vsc, fsc);
	  }
	  total_p += fp;

	  /* Increment dsq like a reversed odometer */
	  for (i = 1; i <= L; i++) 
	    if (dsq[i] < abc->K-1) { dsq[i]++; break; } else { dsq[i] = 0; }
	  if (i > L) break;	/* we're done enumerating sequences */
	}
    }

  /* That sum is subject to significant numerical error because of
   * discretization error in FLogsum(); don't expect it to be too close.
   */
  if (total_p < 0.999 || total_p > 1.001) esl_fatal("Enumeration unit test failed: total Forward p isn't near 1.0 (%g)", total_p);
  if (esl_opt_GetBoolean(go, "-v")) {
    printf("enumeration test: total p is %g\n", total_p);
  }
  
  p7_gmx_Destroy(gx);
  p7_bg_Destroy(bg);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  free(dsq);
  free(seq);
}
#endif /*p7DP_GENERIC_TESTDRIVE*/




/*****************************************************************
 * 5. Test driver.
 *****************************************************************/
/* gcc -g -Wall -Dp7DP_GENERIC_TESTDRIVE -I. -I../easel -L. -L../easel -o dp_generic_utest dp_generic.c -lhmmer -leasel -lm
 */
#ifdef p7DP_GENERIC_TESTDRIVE
#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"

#include "p7_config.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-r",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "set random number seed randomly",                0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be very verbose",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for the generic implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  char            errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS *r    = NULL;
  ESL_ALPHABET   *abc  = NULL;
  P7_HMM         *hmm  = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_BG          *bg   = NULL;
  int             M    = 100;
  int             L    = 200;
  int             nseq = 20;

  if (esl_opt_GetBoolean(go, "-r"))  r = esl_randomness_CreateTimeseeded();
  else                               r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  utest_basic(go);
  
  if ((abc = esl_alphabet_Create(eslAMINO))         == NULL)  esl_fatal("failed to create alphabet");
  if (p7_hmm_Sample(r, M, abc, &hmm)                != eslOK) esl_fatal("failed to sample an HMM");
  if ((bg = p7_bg_Create(abc))                      == NULL)  esl_fatal("failed to create null model");
  if ((gm = p7_profile_Create(hmm->M, abc))         == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL) != eslOK) esl_fatal("failed to config profile");
  if (p7_hmm_Validate    (hmm,     0.0001, errbuf)  != eslOK) esl_fatal("whoops, HMM is bad!");
  if (p7_profile_Validate(gm, 0.0001)               != eslOK) esl_fatal("whoops, profile is bad!");

  utest_viterbi    (go, r, abc, bg, gm, nseq, L);
  utest_forward    (go, r, abc, bg, gm, nseq, L);
  utest_generation (go, r, abc, gm, hmm, bg, nseq);
  utest_enumeration(go, r, abc, 4);	/* can't go much higher than 5; enumeration test is cpu-intensive. */
  
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}

#endif /*p7DP_GENERIC_TESTDRIVE*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
