
#include "p7_config.h"
#include "easel.h"
#include "hmmer.h"

#define MMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_M])
#define IMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_I])
#define DMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_D])
#define XMX(i,s) (xmx[(i) * p7G_NXCELLS + (s)])

int
p7_PostProb(int L, const P7_PROFILE *gm, P7_GMX *fwd, P7_GMX *bck, P7_GMX *pp)
{
  float      **dp   = pp->dp;
  float       *xmx  = pp->xmx;
  int          M    = gm->M;
  int          i,k;
  float        overall_sc = bck->xmx[p7G_NXCELLS*0 + p7G_N];
  
  XMX(0, p7G_N) = 1.0;		
  XMX(0, p7G_B) = expf(fwd->xmx[p7G_NXCELLS*0 + p7G_B] + bck->xmx[p7G_NXCELLS*0 + p7G_B] - overall_sc);
  XMX(0, p7G_E)  = 0.0;
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = 0.0;
  
  for (i = 1; i <= L; i++)
    {
      MMX(i,0) = IMX(i,0) = DMX(i,0) = 0.0;
      for (k = 1; k < M; k++)
	{
	  MMX(i,k) = expf(fwd->dp[i][k*p7G_NSCELLS + p7G_M] + bck->dp[i][k*p7G_NSCELLS + p7G_M] - overall_sc);
	  DMX(i,k) = 0.;
	  IMX(i,k) = expf(fwd->dp[i][k*p7G_NSCELLS + p7G_I] + bck->dp[i][k*p7G_NSCELLS + p7G_I] - overall_sc);
	}
      MMX(i,M)     = expf(fwd->dp[i][M*p7G_NSCELLS + p7G_M] + bck->dp[i][M*p7G_NSCELLS + p7G_M] - overall_sc);
      DMX(i,M)     = 0.;
      
      /* order doesn't matter. note that this whole function is trivially simd parallel */
      XMX(i,p7G_E) = 0.;
      XMX(i,p7G_N) = expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_N] + bck->xmx[p7G_NXCELLS*i + p7G_N] + gm->xsc[p7P_N][p7P_LOOP] - overall_sc);
      XMX(i,p7G_J) = expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_J] + bck->xmx[p7G_NXCELLS*i + p7G_J] + gm->xsc[p7P_J][p7P_LOOP] - overall_sc);
      XMX(i,p7G_B) = 0.;
      XMX(i,p7G_C) = expf(fwd->xmx[p7G_NXCELLS*(i-1) + p7G_C] + bck->xmx[p7G_NXCELLS*i + p7G_C] + gm->xsc[p7P_C][p7P_LOOP] - overall_sc);
    }
  return eslOK;
}


int
p7_OptimalAccuracyDP(int L, const P7_PROFILE *gm, const P7_GMX *pp, P7_GMX *gx, float *ret_e)
{
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          i,k;
  int          M    = gm->M;
  float        esc  = p7_profile_IsLocal(gm) ? 1.0 : 0.0;

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
	  MMX(i,k)     = ESL_MAX(ESL_MAX(MMX(i-1,k-1), IMX(i-1,k-1)),
				 ESL_MAX(DMX(i-1,k-1), XMX(i-1,p7G_B))) + pp->dp[i][k*p7G_NSCELLS + p7G_M];
	  XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), MMX(i,k) * esc);
	  IMX(i,k)     = ESL_MAX(MMX(i-1,k), IMX(i-1,k)) + pp->dp[i][k*p7G_NSCELLS + p7G_I];
	  DMX(i,k)     = ESL_MAX(MMX(i,k-1), DMX(i,k-1));
	} 

      /* last node (k=M) is unrolled; it has no I state, and it has a p=1.0 {MD}->E transition even in local mode */
      MMX(i,M)     = ESL_MAX(ESL_MAX(MMX(i-1,M-1), IMX(i-1,M-1)),
			     ESL_MAX(DMX(i-1,M-1), XMX(i-1,p7G_B))) + pp->dp[i][M*p7G_NSCELLS + p7G_M];
      DMX(i,M)     = ESL_MAX(MMX(i,M-1), DMX(i,M-1));
      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), ESL_MAX(MMX(i,M), DMX(i,M)));

      /* now the special states; it's important that E is already done, and B is done after N,J */
      XMX(i,p7G_J) = ESL_MAX( XMX(i-1,p7G_J) + pp->xmx[i*p7G_NXCELLS + p7G_J], XMX(i,  p7G_E));
      XMX(i,p7G_C) = ESL_MAX( XMX(i-1,p7G_C) + pp->xmx[i*p7G_NXCELLS + p7G_C], XMX(i,  p7G_E));
      XMX(i,p7G_N) =          XMX(i-1,p7G_N) + pp->xmx[i*p7G_NXCELLS + p7G_N];
      XMX(i,p7G_B) = ESL_MAX( XMX(i,  p7G_N),                                  XMX(i,  p7G_J));
    }
  
  *ret_e = XMX(L,p7G_C);
  return eslOK;
}

int
p7_OATrace(int L, const P7_PROFILE *gm, const P7_GMX *pp, const P7_GMX *gx, P7_TRACE *tr)
{
  int     status;
  int     i;			/* position in seq (1..L) */
  int     k;			/* position in model (1..M) */
  int     M   = gm->M;
  float **dp  = gx->dp;
  float  *xmx = gx->xmx;
  float   tol = 1e-4;
  float   esc = p7_profile_IsLocal(gm) ? 1.0f : 0.0f;

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

      if (esl_FCompare(XMX(i, p7G_C), XMX(i-1, p7G_C) + pp->xmx[i*p7G_NXCELLS + p7G_C], tol) == eslOK) {
	tr->i[tr->N-1]    = i--;  /* first C doesn't emit: subsequent ones do */
	status = p7_trace_Append(tr, p7T_C, 0, 0);
      } else if (esl_FCompare(XMX(i, p7G_C), XMX(i, p7G_E), tol) == eslOK) 
	status = p7_trace_Append(tr, p7T_E, 0, 0);
      else ESL_XEXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:		/* E connects from any M state. k set here */
      if (XMX(i, p7G_E) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if (esl_FCompare(XMX(i, p7G_E), MMX(i,M), tol) == eslOK) { k = M; status = p7_trace_Append(tr, p7T_M, M, i); }
      else {
	for (k = M-1; k >= 1; k--)
	  if (esl_FCompare(XMX(i, p7G_E), MMX(i,k) * esc, tol) == eslOK)
	    { status = p7_trace_Append(tr, p7T_M, k, i); break; }
	if (k < 0) ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      }
      break;

    case p7T_M:			/* M connects from i-1,k-1, or B */
      if (MMX(i,k) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      if      (esl_FCompare(MMX(i,k), MMX(i-1,k-1)   + pp->dp[i][k*p7G_NSCELLS + p7G_M], tol) == eslOK) status = p7_trace_Append(tr, p7T_M, k-1, i-1);
      else if (esl_FCompare(MMX(i,k), IMX(i-1,k-1)   + pp->dp[i][k*p7G_NSCELLS + p7G_M], tol) == eslOK) status = p7_trace_Append(tr, p7T_I, k-1, i-1);
      else if (esl_FCompare(MMX(i,k), DMX(i-1,k-1)   + pp->dp[i][k*p7G_NSCELLS + p7G_M], tol) == eslOK) status = p7_trace_Append(tr, p7T_D, k-1, 0);
      else if (esl_FCompare(MMX(i,k), XMX(i-1,p7G_B) + pp->dp[i][k*p7G_NSCELLS + p7G_M], tol) == eslOK) status = p7_trace_Append(tr, p7T_B, 0,   0);
      else ESL_XEXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

      if (status != eslOK) goto ERROR;
      k--; 
      i--;
      break;

    case p7T_D:			/* D connects from M,D at i,k-1 */
      if (DMX(i, k) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);

      if      (esl_FCompare(DMX(i,k), MMX(i, k-1), tol) == eslOK) status = p7_trace_Append(tr, p7T_M, k-1, i);
      else if (esl_FCompare(DMX(i,k), DMX(i, k-1), tol) == eslOK) status = p7_trace_Append(tr, p7T_D, k-1, 0);
      else ESL_XEXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      if (status != eslOK) goto ERROR;
      k--;
      break;

    case p7T_I:			/* I connects from M,I at i-1,k*/
      if (IMX(i,k) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      if      (esl_FCompare(IMX(i,k), MMX(i-1,k) + pp->dp[i][k*p7G_NSCELLS + p7G_I], tol) == eslOK) status = p7_trace_Append(tr, p7T_M, k, i-1);
      else if (esl_FCompare(IMX(i,k), IMX(i-1,k) + pp->dp[i][k*p7G_NSCELLS + p7G_I], tol) == eslOK) status = p7_trace_Append(tr, p7T_I, k, i-1);
      else ESL_XEXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      if (status != eslOK) goto ERROR;
      i--;
      break;

    case p7T_N:			/* N connects from S, N */
      if (i == 0) status = p7_trace_Append(tr, p7T_S, 0, 0);
      else if (esl_FCompare(XMX(i,p7G_N), XMX(i-1, p7G_N) + pp->xmx[i*p7G_NXCELLS + p7G_N], tol) == eslOK)
	{
	  tr->i[tr->N-1] = i--;
	  status = p7_trace_Append(tr, p7T_N, 0, 0);
	} 
      else ESL_XEXCEPTION(eslFAIL, "N at i=%d couldn't be traced", i);
      break;

    case p7T_B:			/* B connects from N, J */
      if (XMX(i,p7G_B) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if (esl_FCompare(XMX(i,p7G_B), XMX(i, p7G_N), tol)  == eslOK)
	status = p7_trace_Append(tr, p7T_N, 0, 0);
      else if (esl_FCompare(XMX(i,p7G_B),  XMX(i, p7G_J), tol) == eslOK)
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      else  ESL_XEXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:			/* J connects from E(i) or J(i-1) */
      if (XMX(i,p7G_J) < 0.0) ESL_XEXCEPTION(eslFAIL, "impossible J reached at i=%d", i);

      if (esl_FCompare(XMX(i,p7G_J), XMX(i-1,p7G_J) + pp->xmx[i*p7G_NXCELLS + p7G_J], tol) == eslOK) {
	tr->i[tr->N-1] = i--;
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      } else if (esl_FCompare(XMX(i,p7G_J), XMX(i,p7G_E), tol) == eslOK) 
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
