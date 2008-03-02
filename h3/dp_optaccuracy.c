/* Implementations of posterior decoding and optimal accuracy alignment.
 * 
 * Contents:
 * 
 * SRE, Fri Feb 29 12:48:46 2008 [Janelia]
 * SVN $Id$
 */

#include "p7_config.h"

#include <float.h>

#include "easel.h"
#include "hmmer.h"



#define MMX(i,k)      (dp[(i)][(k) * p7G_NSCELLS + p7G_M])
#define IMX(i,k)      (dp[(i)][(k) * p7G_NSCELLS + p7G_I])
#define DMX(i,k)      (dp[(i)][(k) * p7G_NSCELLS + p7G_D])
#define XMX(i,s)      (xmx[(i) * p7G_NXCELLS + (s)])
#define TSCDELTA(s,k) ( (tsc[(k) * p7P_NTRANS + (s)] == -eslINFINITY) ? FLT_MIN : 1.0)

/* The TSCDELTA is used to make impossible paths impossible in the
 * optimal accuracy decoding algorithm; see Kall et al (2005). What we
 * want to do is multiply by a Kronecker delta that's 1 when the
 * transition probability is finite, and 0 when it's zero (when the
 * log prob is -eslINFINITY). But we can't do that easily, when we're
 * in log space, because 0 * -eslINFINITY = NaN. Instead, we use a
 * tiny number (FLT_MIN, ~1e-37).
 * 
 * A side concern is that we don't want to put a bunch of if-else
 * branches in the code; compilers should be able to generate more
 * efficient code from the TSCDELTA() construction.
 */


/* Function:  p7_PosteriorDecoding()
 * Synopsis:  Calculate posterior decoding probabilities, using forward/backward.
 * Incept:    SRE, Fri Feb 29 10:16:21 2008 [Janelia]
 *
 * Purpose:   Calculates a posterior decoding of the residues <1..L> in
 *            a target sequence of length <L>, given profile <gm>, a
 *            filled Forward matrix <fwd>, and a Backward matrix <bck>.
 *            
 *            The resulting posterior decoding is stored in a DP matrix <pp>,
 *            provided by the caller.
 *            
 *            Each residue <i> must have been emitted by match state
 *            <1..M>, insert state <1..M-1>, or an NN, CC, or JJ loop
 *            transition.  For <dp = pp->dp>, <xmx = pp->xmx>,
 *            <MMX(i,k)> is the probability that match <k> emitted
 *            residue <i>; <IMX(i,k)> is the probability that insert
 *            <k> emitted residue <i>; <XMX(i,N)>,<XMX(i,C)>,
 *            <XMX(i,J)> are the probabilities that residue <i> was
 *            emitted on an NN, CC, or JJ transition. The sum over all
 *            these possibilities for a given residue <i> is 1.0.
 *            
 *            The caller may pass the Backward matrix <bck> as <pp>,
 *            in which case <bck> will be overwritten with
 *            <pp>. However, the caller may \emph{not} overwrite <fwd>
 *            this way; an <(i-1)> dependency in the calculation of
 *            NN, CC, JJ transitions prevents this.
 *
 * Args:      L    - length of the target sequence
 *            gm   - profile (must be the same that was used to fill <fwd>, <bck>).
 *            fwd  - filled Forward matrix 
 *            bck  - filled Backward matrix
 *            pp   - RESULT: posterior decoding matrix.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_PosteriorDecoding(int L, const P7_PROFILE *gm, const P7_GMX *fwd, P7_GMX *bck, P7_GMX *pp)
{
  float      **dp   = pp->dp;
  float       *xmx  = pp->xmx;
  int          M    = gm->M;
  int          i,k;
  float        overall_sc = bck->xmx[p7G_NXCELLS*0 + p7G_N];
  
  XMX(0, p7G_E) = 0.0;
  XMX(0, p7G_N) = 0.0;		
  XMX(0, p7G_J) = 0.0;		
  XMX(0, p7G_B) = 0.0;
  XMX(0, p7G_C) = 0.0;
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = 0.0;
  
  for (i = 1; i <= L; i++)
    {
      MMX(i,0) = IMX(i,0) = DMX(i,0) = 0.0;
      for (k = 1; k < M; k++)
	{
	  MMX(i,k) = expf(fwd->dp[i][k*p7G_NSCELLS + p7G_M] + bck->dp[i][k*p7G_NSCELLS + p7G_M] - overall_sc);
	  IMX(i,k) = expf(fwd->dp[i][k*p7G_NSCELLS + p7G_I] + bck->dp[i][k*p7G_NSCELLS + p7G_I] - overall_sc);
	  DMX(i,k) = 0.;
	}
      MMX(i,M)     = expf(fwd->dp[i][M*p7G_NSCELLS + p7G_M] + bck->dp[i][M*p7G_NSCELLS + p7G_M] - overall_sc);
      DMX(i,M)     = 0.;
      
      /* order doesn't matter.  note that this whole function is trivially simd parallel */
      XMX(i,p7G_E) = 0.;
      XMX(i,p7G_N) = expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_N] + bck->xmx[p7G_NXCELLS*i + p7G_N] + gm->xsc[p7P_N][p7P_LOOP] - overall_sc);
      XMX(i,p7G_J) = expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_J] + bck->xmx[p7G_NXCELLS*i + p7G_J] + gm->xsc[p7P_J][p7P_LOOP] - overall_sc);
      XMX(i,p7G_B) = 0.;
      XMX(i,p7G_C) = expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_C] + bck->xmx[p7G_NXCELLS*i + p7G_C] + gm->xsc[p7P_C][p7P_LOOP] - overall_sc);
    }
  return eslOK;
}


/* Function:  p7_OptimalAccuracyDP()
 * Synopsis:  Optimal accuracy decoding: fill. 
 * Incept:    SRE, Fri Feb 29 11:56:49 2008 [Janelia]
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm \citep{Kall05}.
 *            
 *            Caller provides the posterior decoding matrix <pp>,
 *            which was calculated by Forward/Backward on a target sequence
 *            of length <L> using the query model <gm>.
 *            
 *            Caller also provides a DP matrix <gx>, allocated for the
 *            <gm->M> by <L> comparison. The routine fills this in
 *            with OA scores.
 *            
 * Args:      L     - length of original target sequence 
 *            gm    - query profile      
 *            pp    - posterior decoding matrix created by <p7_PosteriorDecoding()>
 *            gx    - RESULT: caller provided DP matrix for <gm->M> by <L> 
 *            ret_e - RETURN: expected number of correctly decoded positions 
 *
 * Returns:   <eslOK> on success, and <*ret_e> contains the final OA
 *            score, which is the expected number of correctly decoded
 *            positions in the target sequence (up to <L>).
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_OptimalAccuracyDP(int L, const P7_PROFILE *gm, const P7_GMX *pp, P7_GMX *gx, float *ret_e)
{
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float const *tsc  = gm->tsc;
  int          i,k;
  int          M    = gm->M;
  float        esc  = p7_profile_IsLocal(gm) ? 1.0 : 0.0;
  float        t1, t2;

  /* Initialization of the zero row (i=0; no residues to account for.  */
  XMX(0,p7G_N) = 0.;                                          /* S->N, p=1            */
  XMX(0,p7G_B) = 0.;                                          /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;  /* need seq to get here */
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = -eslINFINITY;            /* need seq to get here */

  for (i = 1; i <= L; i++)
    {
      MMX(i,0) = IMX(i,0) = DMX(i,0) = XMX(i,p7G_E) = -eslINFINITY;

      for (k = 1; k < M; k++)
	{
	  MMX(i,k)     = ESL_MAX(ESL_MAX(TSCDELTA(p7P_MM, k-1) * (MMX(i-1,k-1)  + pp->dp[i][k*p7G_NSCELLS + p7G_M]),
					 TSCDELTA(p7P_IM, k-1) * (IMX(i-1,k-1)  + pp->dp[i][k*p7G_NSCELLS + p7G_M])),
				 ESL_MAX(TSCDELTA(p7P_DM, k-1) * (DMX(i-1,k-1)  + pp->dp[i][k*p7G_NSCELLS + p7G_M]),
					 TSCDELTA(p7P_BM, k-1) * (XMX(i-1,p7G_B)+ pp->dp[i][k*p7G_NSCELLS + p7G_M])));

	  XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), 
				 esc * MMX(i,k));

	  IMX(i,k)     = ESL_MAX(TSCDELTA(p7P_MI, k) * (MMX(i-1,k) + pp->dp[i][k*p7G_NSCELLS + p7G_I]),
				 TSCDELTA(p7P_II, k) * (IMX(i-1,k) + pp->dp[i][k*p7G_NSCELLS + p7G_I]));

	  DMX(i,k)     = ESL_MAX(TSCDELTA(p7P_MD, k-1) * MMX(i,k-1),
				 TSCDELTA(p7P_DD, k-1) * DMX(i,k-1));
	} 

      /* last node (k=M) is unrolled; it has no I state, and it has a p=1.0 {MD}->E transition even in local mode */
      MMX(i,M)     = ESL_MAX(ESL_MAX(TSCDELTA(p7P_MM, M-1) * (MMX(i-1,M-1)  + pp->dp[i][M*p7G_NSCELLS + p7G_M]),
				     TSCDELTA(p7P_IM, M-1) * (IMX(i-1,M-1)  + pp->dp[i][M*p7G_NSCELLS + p7G_M])),
			     ESL_MAX(TSCDELTA(p7P_DM, M-1) * (DMX(i-1,M-1)  + pp->dp[i][M*p7G_NSCELLS + p7G_M]),
				     TSCDELTA(p7P_BM, M-1) * (XMX(i-1,p7G_B)+ pp->dp[i][M*p7G_NSCELLS + p7G_M])));

      DMX(i,M)     = ESL_MAX(TSCDELTA(p7P_MD, M-1) * MMX(i,M-1),
			     TSCDELTA(p7P_DD, M-1) * DMX(i,M-1));

      /* note: we calculated XMX before DMX in the loop, because we probably had MMX(i,k) in a register. 
       * but now we can't do that, because XMX depends on DMX
       */
      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), ESL_MAX(MMX(i,M), DMX(i, M)));

      /* now the special states; it's important that E is already done, and B is done after N,J */
      t1 = ( (gm->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
      t2 = ( (gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
      XMX(i, p7G_J) = ESL_MAX( t1 * (XMX(i-1,p7G_J) + pp->xmx[i*p7G_NXCELLS + p7G_J]),
			       t2 * XMX(i,  p7G_E));

      t1 = ( (gm->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
      t2 = ( (gm->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
      XMX(i,p7G_C) = ESL_MAX( t1 * (XMX(i-1,p7G_C) + pp->xmx[i*p7G_NXCELLS + p7G_C]),
			      t2 * XMX(i,  p7G_E));
      
      t1 = ( (gm->xsc[p7P_N][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
      XMX(i,p7G_N) = t1 *  (XMX(i-1,p7G_N) + pp->xmx[i*p7G_NXCELLS + p7G_N]);

      t1 = ( (gm->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
      t2 = ( (gm->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
      XMX(i,p7G_B) = ESL_MAX( t1 * XMX(i,  p7G_N), 
			      t2 * XMX(i,  p7G_J));
    }
  
  *ret_e = XMX(L,p7G_C);
  return eslOK;
}


/* Function:  p7_OATrace()
 * Synopsis:  Optimal accuracy decoding: traceback.
 * Incept:    SRE, Fri Feb 29 12:59:11 2008 [Janelia]
 *
 * Purpose:   The traceback stage of the optimal accuracy decoding algorithm
 *            \citep{Kall05}.
 *            
 *            Caller provides the OA DP matrix <gx> that was just
 *            calculated by <p7_OptimalAccuracyDP()>, as well as the
 *            posterior decoding matrix <pp>, which was calculated by
 *            Forward/Backward on a target sequence of length <L>
 *            using the query model <gm>.
 *            
 *            The resulting optimal accuracy decoding traceback is put
 *            in a caller-provided traceback structure <tr>. This
 *            structure will be reallocated internally if necessary.
 *
 * Args:      L     - length of original target sequence 
 *            gm    - query profile      
 *            pp    - posterior decoding matrix created by <p7_PosteriorDecoding()>
 *            gx    - OA DP matrix calculated by  <p7_OptimalAccuracyDP()>
 *            tr    - RESULT: OA traceback
 *
 * Returns:   <eslOK> on success, and <tr> contains the OA traceback.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_OATrace(int L, const P7_PROFILE *gm, const P7_GMX *pp, const P7_GMX *gx, P7_TRACE *tr)
{
  int     status;
  int     i;			/* position in seq (1..L) */
  int     k;			/* position in model (1..M) */
  int           M   = gm->M;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float const *tsc  = gm->tsc;
  float   tol = 1e-4;
  float   t1, t2;

  if ((status = p7_trace_Reuse(tr)) != eslOK) goto ERROR;

  /* Initialization. (back to front. ReverseTrace() called later.)  */
  if ((status = p7_trace_Append(tr, p7T_T, 0, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_Append(tr, p7T_C, 0, 0)) != eslOK) goto ERROR;
  i    = L;			/* next position to explain in seq */

  /* Traceback  */
  while (tr->st[tr->N-1] != p7T_S) {

    switch (tr->st[tr->N-1]) {
    case p7T_C:		/* C(i) comes from C(i-1) or E(i) */
      if   (XMX(i,p7G_C) < 0.0f) ESL_XEXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      t1 = ( (gm->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
      t2 = ( (gm->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);

      if (esl_FCompare(XMX(i, p7G_C), t1 * (XMX(i-1, p7G_C) + pp->xmx[i*p7G_NXCELLS + p7G_C]), tol) == eslOK) {
	tr->i[tr->N-1]    = i--;  /* first C doesn't emit: subsequent ones do */
	status = p7_trace_Append(tr, p7T_C, 0, 0);
      } else if (esl_FCompare(XMX(i, p7G_C), t2 * XMX(i, p7G_E), tol) == eslOK) 
	status = p7_trace_Append(tr, p7T_E, 0, 0);
      else ESL_XEXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:		/* E connects from any M state, or D_M. k set here */
      if (XMX(i, p7G_E) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      
      if      (esl_FCompare(XMX(i, p7G_E), MMX(i,M), tol) == eslOK) { k = M; status = p7_trace_Append(tr, p7T_M, M, i); }
      else if (esl_FCompare(XMX(i, p7G_E), DMX(i,M), tol) == eslOK) { k = M; status = p7_trace_Append(tr, p7T_D, M, i); }
      else if (p7_profile_IsLocal(gm)) /* local exits allowed? */
	{
	  for (k = M-1; k >= 1; k--)
	    if (esl_FCompare(XMX(i, p7G_E), MMX(i,k), tol) == eslOK) { status = p7_trace_Append(tr, p7T_M, k, i); break; }
	  if (k <= 0) ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      }
      else ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      break;

    case p7T_M:			/* M connects from i-1,k-1, or B */
      if (MMX(i,k) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      if      (esl_FCompare(MMX(i,k), TSCDELTA(p7P_MM, k-1) * (MMX(i-1,k-1)   + pp->dp[i][k*p7G_NSCELLS + p7G_M]), tol) == eslOK) status = p7_trace_Append(tr, p7T_M, k-1, i-1);
      else if (esl_FCompare(MMX(i,k), TSCDELTA(p7P_IM, k-1) * (IMX(i-1,k-1)   + pp->dp[i][k*p7G_NSCELLS + p7G_M]), tol) == eslOK) status = p7_trace_Append(tr, p7T_I, k-1, i-1);
      else if (esl_FCompare(MMX(i,k), TSCDELTA(p7P_DM, k-1) * (DMX(i-1,k-1)   + pp->dp[i][k*p7G_NSCELLS + p7G_M]), tol) == eslOK) status = p7_trace_Append(tr, p7T_D, k-1, 0);
      else if (esl_FCompare(MMX(i,k), TSCDELTA(p7P_BM, k-1) * (XMX(i-1,p7G_B) + pp->dp[i][k*p7G_NSCELLS + p7G_M]), tol) == eslOK) status = p7_trace_Append(tr, p7T_B, 0,   0);
      else ESL_XEXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

      if (status != eslOK) goto ERROR;
      k--; i--;
      break;

    case p7T_D:			/* D connects from M,D at i,k-1 */
      if (DMX(i, k) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);

      if      (esl_FCompare(DMX(i,k), TSCDELTA(p7P_MD, k-1) * MMX(i, k-1), tol) == eslOK) status = p7_trace_Append(tr, p7T_M, k-1, i);
      else if (esl_FCompare(DMX(i,k), TSCDELTA(p7P_DD, k-1) * DMX(i, k-1), tol) == eslOK) status = p7_trace_Append(tr, p7T_D, k-1, 0);
      else ESL_XEXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      if (status != eslOK) goto ERROR;
      k--;
      break;

    case p7T_I:			/* I connects from M,I at i-1,k*/
      if (IMX(i,k) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      if      (esl_FCompare(IMX(i,k), TSCDELTA(p7P_MI, k) * (MMX(i-1,k) + pp->dp[i][k*p7G_NSCELLS + p7G_I]), tol) == eslOK) status = p7_trace_Append(tr, p7T_M, k, i-1);
      else if (esl_FCompare(IMX(i,k), TSCDELTA(p7P_II, k) * (IMX(i-1,k) + pp->dp[i][k*p7G_NSCELLS + p7G_I]), tol) == eslOK) status = p7_trace_Append(tr, p7T_I, k, i-1);
      else ESL_XEXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      if (status != eslOK) goto ERROR;
      i--;
      break;

    case p7T_N:			/* N connects from S, N */
      t1 = ( (gm->xsc[p7P_N][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);

      if (i == 0) status = p7_trace_Append(tr, p7T_S, 0, 0);
      else if (esl_FCompare(XMX(i,p7G_N), t1 * (XMX(i-1, p7G_N) + pp->xmx[i*p7G_NXCELLS + p7G_N]), tol) == eslOK)
	{
	  tr->i[tr->N-1] = i--;
	  status = p7_trace_Append(tr, p7T_N, 0, 0);
	} 
      else ESL_XEXCEPTION(eslFAIL, "N at i=%d couldn't be traced", i);
      break;

    case p7T_B:			/* B connects from N, J */
      if (XMX(i,p7G_B) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      t1 = ( (gm->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
      t2 = ( (gm->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);

      if      (esl_FCompare(XMX(i,p7G_B), t1 * XMX(i, p7G_N), tol) == eslOK) status = p7_trace_Append(tr, p7T_N, 0, 0);
      else if (esl_FCompare(XMX(i,p7G_B), t2 * XMX(i, p7G_J), tol) == eslOK) status = p7_trace_Append(tr, p7T_J, 0, 0);
      else    ESL_XEXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:			/* J connects from E(i) or J(i-1) */
      if (XMX(i,p7G_J) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible J reached at i=%d", i);

      t1 = ( (gm->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
      t2 = ( (gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);

      if (esl_FCompare(XMX(i,p7G_J), t1 * (XMX(i-1,p7G_J) + pp->xmx[i*p7G_NXCELLS + p7G_J]), tol) == eslOK) {
	tr->i[tr->N-1] = i--;
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      } else if (esl_FCompare(XMX(i,p7G_J), t2 * XMX(i,p7G_E), tol) == eslOK) 
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
