/* Optimal accuracy alignment; SSE version.
 * 
 * Contents:
 *   1. Optimal accuracy alignment, DP fill
 *   2. OA traceback
 *   3. Benchmark driver
 *   4. Unit tests
 *   5. Test driver
 *   6. Example
 * 
 * SRE, Mon Aug 18 20:01:01 2008 [Casa de Gatos]
 */
#include <p7_config.h>

#include <float.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"


/*****************************************************************
 * 1. Optimal accuracy alignment, DP fill
 *****************************************************************/

/* Function:  p7_OptimalAccuracy()
 * Synopsis:  DP fill of an optimal accuracy alignment calculation.
 * Incept:    SRE, Mon Aug 18 11:04:48 2008 [Janelia]
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm \citep{Kall05}.
 *            
 *            Caller provides the posterior decoding matrix <pp>,
 *            which was calculated by Forward/Backward on a target sequence
 *            of length <pp->L> using the query model <om>.
 *            
 *            Caller also provides a DP matrix <ox>, allocated for a full
 *            <om->M> by <L> comparison. The routine fills this in
 *            with OA scores.
 *  
 * Args:      om    - query profile      
 *            pp    - posterior decoding matrix created by <p7_GPosteriorDecoding()>
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
p7_OptimalAccuracy_avx512(const P7_OPROFILE *om, const P7_OMX *pp, P7_OMX *ox, float *ret_e)
{
  register __m128 mpv, dpv, ipv;   /* previous row values                                       */
  register __m128 sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128 xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128 xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128 dcv;
  float  *xmx = ox->xmx;
  __m128 *dpc = ox->dpf[0];        /* current row, for use in {MDI}MO(dpp,q) access macro       */
  __m128 *dpp;                     /* previous row, for use in {MDI}MO(dpp,q) access macro      */
  __m128 *ppp;			   /* quads in the <pp> posterior probability matrix            */
  __m128 *tp;			   /* quads in the <om->tfv> transition scores                  */
  __m128 zerov = _mm_setzero_ps();
  __m128 infv  = _mm_set1_ps(-eslINFINITY);
  int M = om->M;
  int Q = p7O_NQF(M);
  int q;
  int j;
  int i;
  float t1, t2;

  ox->M = om->M;
  ox->L = pp->L;
  for (q = 0; q < Q; q++) MMO(dpc, q) = IMO(dpc,q) = DMO(dpc,q) = infv;
  XMXo(0, p7X_E)    = -eslINFINITY;
  XMXo(0, p7X_N)    = 0.;
  XMXo(0, p7X_J)    = -eslINFINITY;
  XMXo(0, p7X_B)    = 0.;
  XMXo(0, p7X_C)    = -eslINFINITY;

  for (i = 1; i <= pp->L; i++)
    {
      dpp = dpc;		/* previous DP row in OA matrix */
      dpc = ox->dpf[i];   	/* current DP row in OA matrix  */
      ppp = pp->dpf[i];		/* current row in the posterior probabilities per position */
      tp  = om->tfv;		/* transition probabilities */
      dcv = infv;
      xEv = infv;
      xBv = _mm_set1_ps(XMXo(i-1, p7X_B));

      mpv = esl_sse_rightshift_ps(MMO(dpp,Q-1), infv);  /* Right shifts by 4 bytes. 4,8,12,x becomes x,4,8,12. */
      dpv = esl_sse_rightshift_ps(DMO(dpp,Q-1), infv);
      ipv = esl_sse_rightshift_ps(IMO(dpp,Q-1), infv);
      for (q = 0; q < Q; q++)
	{
	  sv  =                _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), xBv);  tp++;
	  sv  = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), mpv)); tp++;
	  sv  = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), ipv)); tp++;
	  sv  = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), dpv)); tp++;
	  sv  = _mm_add_ps(sv, *ppp);                                      ppp += 2;
	  xEv = _mm_max_ps(xEv, sv);
	  
	  mpv = MMO(dpp,q);
	  dpv = DMO(dpp,q);
	  ipv = IMO(dpp,q);

	  MMO(dpc,q) = sv;
	  DMO(dpc,q) = dcv;

	  dcv = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), sv); tp++;

	  sv         =                _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), mpv);   tp++;
	  sv         = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), ipv));  tp++;
	  IMO(dpc,q) = _mm_add_ps(sv, *ppp);                                       ppp++;
	}
      
      /* dcv has carried through from end of q loop above; store it 
       * in first pass, we add M->D and D->D path into DMX
       */
      dcv = esl_sse_rightshift_ps(dcv, infv); 
      tp  = om->tfv + 7*Q;	/* set tp to start of the DD's */
      for (q = 0; q < Q; q++)
	{
	  DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
	  dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), DMO(dpc,q));   tp++;
	}

      /* fully serialized D->D; can optimize later */
      for (j = 1; j < 4; j++)
	{
	  dcv = esl_sse_rightshift_ps(dcv, infv);
	  tp  = om->tfv + 7*Q;	
	  for (q = 0; q < Q; q++)
	    {
	      DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
	      dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), dcv);   tp++;
	    }
	}

      /* D->E paths */
      for (q = 0; q < Q; q++) xEv = _mm_max_ps(xEv, DMO(dpc,q));
      
      /* Specials */
      esl_sse_hmax_ps(xEv, &(XMXo(i,p7X_E)));
      
      t1 = ( (om->xf[p7O_J][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[(i-1)*p7X_NXCELLS+p7X_J] + pp->xmx[i*p7X_NXCELLS+p7X_J]);
      t2 = ( (om->xf[p7O_E][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[   i *p7X_NXCELLS+p7X_E]);
      ox->xmx[i*p7X_NXCELLS+p7X_J] = ESL_MAX(t1, t2);

      t1 = ( (om->xf[p7O_C][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] + pp->xmx[i*p7X_NXCELLS+p7X_C]);
      t2 = ( (om->xf[p7O_E][p7O_MOVE] == 0.0) ? 0.0 : ox->xmx[   i *p7X_NXCELLS+p7X_E]);
      ox->xmx[i*p7X_NXCELLS+p7X_C] = ESL_MAX(t1, t2);
      
      ox->xmx[i*p7X_NXCELLS+p7X_N] = ((om->xf[p7O_N][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[(i-1)*p7X_NXCELLS+p7X_N] + pp->xmx[i*p7X_NXCELLS+p7X_N]);
      
      t1 = ( (om->xf[p7O_N][p7O_MOVE] == 0.0) ? 0.0 : ox->xmx[i*p7X_NXCELLS+p7X_N]);
      t2 = ( (om->xf[p7O_J][p7O_MOVE] == 0.0) ? 0.0 : ox->xmx[i*p7X_NXCELLS+p7X_J]);
      ox->xmx[i*p7X_NXCELLS+p7X_B] = ESL_MAX(t1, t2);
    }

  *ret_e = ox->xmx[pp->L*p7X_NXCELLS+p7X_C];
  return eslOK;
}
/*------------------- end, OA DP fill ---------------------------*/





/*****************************************************************
 * 2. OA traceback
 *****************************************************************/

static inline float get_postprob(const P7_OMX *pp, int scur, int sprv, int k, int i);

static inline int select_m(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int select_d(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int select_i(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int select_n(int i);
static inline int select_c(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, int i);
static inline int select_j(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, int i);
static inline int select_e(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int *ret_k);
static inline int select_b(const P7_OPROFILE *om,                   const P7_OMX *ox, int i);


/* Function:  p7_OATrace()
 * Synopsis:  Optimal accuracy decoding: traceback.
 * Incept:    SRE, Mon Aug 18 13:53:33 2008 [Janelia]
 *
 * Purpose:   The traceback stage of the optimal accuracy decoding algorithm
 *            \citep{Kall05}.
 *            
 *            Caller provides the OA DP matrix <ox> that was just
 *            calculated by <p7_OptimalAccuracyDP()>, as well as the
 *            posterior decoding matrix <pp>, which was calculated by
 *            Forward/Backward on a target sequence using the query
 *            model <gm>. Because the calculation depends only on
 *            <pp>, the target sequence itself need not be provided.
 *            
 *            The resulting optimal accuracy decoding traceback is put
 *            in a caller-provided traceback structure <tr>, which the
 *            caller has allocated for optional posterior probability
 *            annotation on residues (with <p7_trace_CreateWithPP()>,
 *            generally). This structure will be reallocated
 *            internally if necessary.
 *
 * Args:      om  - profile
 *            pp  - posterior probability matrix
 *            ox  - OA matrix to trace, LxM
 *            tr  - storage for the recovered traceback
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if the trace <tr> isn't empty (needs to be Reuse()'d).
 */
int
p7_OATrace_avx512(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr)
{
  int   i   = ox->L;		/* position in sequence 1..L */
  int   k   = 0;		/* position in model 1..M */
  int   s0, s1;			/* choice of a state */
  float postprob;
  int   status;			
  
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  if ((status = p7_trace_AppendWithPP(tr, p7T_T, k, i, 0.0)) != eslOK) return status;
  if ((status = p7_trace_AppendWithPP(tr, p7T_C, k, i, 0.0)) != eslOK) return status;

  s0 = tr->st[tr->N-1];
  while (s0 != p7T_S)
    {
      switch (s0) {
      case p7T_M: s1 = select_m(om,     ox, i, k);  k--; i--; break;
      case p7T_D: s1 = select_d(om,     ox, i, k);  k--;      break;
      case p7T_I: s1 = select_i(om,     ox, i, k);       i--; break;
      case p7T_N: s1 = select_n(i);                           break;
      case p7T_C: s1 = select_c(om, pp, ox, i);               break;
      case p7T_J: s1 = select_j(om, pp, ox, i);               break;
      case p7T_E: s1 = select_e(om,     ox, i, &k);           break;
      case p7T_B: s1 = select_b(om,     ox, i);               break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (s1 == -1) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");

      postprob = get_postprob(pp, s1, s0, k, i);
      if ((status = p7_trace_AppendWithPP(tr, s1, k, i, postprob)) != eslOK) return status;

      if ( (s1 == p7T_N || s1 == p7T_J || s1 == p7T_C) && s1 == s0) i--;
      s0 = s1;
    } /* end traceback, at S state */
  tr->M = om->M;
  tr->L = ox->L;
  return p7_trace_Reverse(tr);
}

static inline float
get_postprob(const P7_OMX *pp, int scur, int sprv, int k, int i)
{
  int     Q     = p7O_NQF(pp->M);
  int     q     = (k-1) % Q;		/* (q,r) is position of the current DP cell M(i,k) */
  int     r     = (k-1) / Q;
  union { __m128 v; float p[4]; } u;

  switch (scur) {
  case p7T_M: u.v = MMO(pp->dpf[i], q); return u.p[r]; 
  case p7T_I: u.v = IMO(pp->dpf[i], q); return u.p[r]; 
  case p7T_N: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_N];
  case p7T_C: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_C];
  case p7T_J: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_J];
  default:    return 0.0;
  }
}

/* M(i,k) is reached from B(i-1), M(i-1,k-1), D(i-1,k-1), or I(i-1,k-1). */
static inline int
select_m(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;		/* (q,r) is position of the current DP cell M(i,k) */
  int     r     = (k-1) / Q;
  __m128 *tp    = om->tfv + 7*q;       	/* *tp now at start of transitions to cur cell M(i,k) */
  __m128  xBv   = _mm_set1_ps(ox->xmx[(i-1)*p7X_NXCELLS+p7X_B]);
  __m128  mpv, dpv, ipv;
  union { __m128 v; float p[4]; } u, tv;
  float   path[4];
  int     state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
  
  if (q > 0) {
    mpv = ox->dpf[i-1][(q-1)*3 + p7X_M];
    dpv = ox->dpf[i-1][(q-1)*3 + p7X_D];
    ipv = ox->dpf[i-1][(q-1)*3 + p7X_I];
  } else {
    mpv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_M]);
    dpv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_D]);
    ipv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_I]);
  }	  

  /* paths are numbered so that most desirable choice in case of tie is first. */
  u.v = xBv;  tv.v = *tp;  path[3] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
  u.v = mpv;  tv.v = *tp;  path[0] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
  u.v = ipv;  tv.v = *tp;  path[1] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
  u.v = dpv;  tv.v = *tp;  path[2] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  
  return state[esl_vec_FArgMax(path, 4)];
}


/* D(i,k) is reached from M(i, k-1) or D(i,k-1). */
static inline int
select_d(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;		/* (q,r) is position of the current DP cell D(i,k) */
  int     r     = (k-1) / Q;
  union { __m128 v; float p[4]; } mpv, dpv, tmdv, tddv;
  float   path[2];

  if (q > 0) {
    mpv.v  = ox->dpf[i][(q-1)*3 + p7X_M];
    dpv.v  = ox->dpf[i][(q-1)*3 + p7X_D];
    tmdv.v = om->tfv[7*(q-1) + p7O_MD];
    tddv.v = om->tfv[7*Q + (q-1)];
  } else {
    mpv.v  = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*3 + p7X_M]);
    dpv.v  = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*3 + p7X_D]);
    tmdv.v = esl_sse_rightshiftz_float(om->tfv[7*(Q-1) + p7O_MD]);
    tddv.v = esl_sse_rightshiftz_float(om->tfv[8*Q-1]);
  }	  

  path[0] = ((tmdv.p[r] == 0.0) ? -eslINFINITY : mpv.p[r]);
  path[1] = ((tddv.p[r] == 0.0) ? -eslINFINITY : dpv.p[r]);
  return  ((path[0] >= path[1]) ? p7T_M : p7T_D);
}

/* I(i,k) is reached from M(i-1, k) or I(i-1,k). */
static inline int
select_i(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q    = p7O_NQF(ox->M);
  int     q    = (k-1) % Q;		/* (q,r) is position of the current DP cell D(i,k) */
  int     r    = (k-1) / Q;
  __m128 *tp   = om->tfv + 7*q + p7O_MI;
  union { __m128 v; float p[4]; } tv, mpv, ipv;
  float   path[2];

  mpv.v = ox->dpf[i-1][q*3 + p7X_M]; tv.v = *tp;  path[0] = ((tv.p[r] == 0.0) ? -eslINFINITY : mpv.p[r]);  tp++;
  ipv.v = ox->dpf[i-1][q*3 + p7X_I]; tv.v = *tp;  path[1] = ((tv.p[r] == 0.0) ? -eslINFINITY : ipv.p[r]);  
  return  ((path[0] >= path[1]) ? p7T_M : p7T_I);
}

/* N(i) must come from N(i-1) for i>0; else it comes from S */
static inline int
select_n(int i)
{
  return ((i==0) ? p7T_S : p7T_N);
}

/* C(i) is reached from E(i) or C(i-1). */
static inline int
select_c(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, int i)
{
  float path[2];
  path[0] = ( (om->xf[p7O_C][p7O_LOOP] == 0.0) ? -eslINFINITY : ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] + pp->xmx[i*p7X_NXCELLS+p7X_C]);
  path[1] = ( (om->xf[p7O_E][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[   i *p7X_NXCELLS+p7X_E]);
  return  ((path[0] > path[1]) ? p7T_C : p7T_E);
}

/* J(i) is reached from E(i) or J(i-1). */
static inline int
select_j(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, int i)
{
  float path[2];

  path[0] = ( (om->xf[p7O_J][p7O_LOOP] == 0.0) ? -eslINFINITY : ox->xmx[(i-1)*p7X_NXCELLS+p7X_J] + pp->xmx[i*p7X_NXCELLS+p7X_J]);
  path[1] = ( (om->xf[p7O_E][p7O_LOOP] == 0.0) ? -eslINFINITY : ox->xmx[   i *p7X_NXCELLS+p7X_E]);
  return  ((path[0] > path[1]) ? p7T_J : p7T_E);
}
 
/* E(i) is reached from any M(i, k=1..M) or D(i, k=2..M). */
/* This assumes all M_k->E, D_k->E are 1.0 */
static inline int
select_e(const P7_OPROFILE *om, const P7_OMX *ox, int i, int *ret_k)
{
  int     Q     = p7O_NQF(ox->M);
  __m128 *dp    = ox->dpf[i];
  union { __m128 v; float p[4]; } u;
  float  max   = -eslINFINITY;
  int    smax, kmax;
  int    q,r;

  /* precedence rules in case of ties here are a little tricky: M beats D: note the >= max!  */
  for (q = 0; q < Q; q++)
    {
      u.v   = *dp; dp++;  for (r = 0; r < 4; r++) if (u.p[r] >= max) { max = u.p[r]; smax = p7T_M; kmax = r*Q + q + 1; }
      u.v   = *dp; dp+=2; for (r = 0; r < 4; r++) if (u.p[r] > max)  { max = u.p[r]; smax = p7T_D; kmax = r*Q + q + 1; }
    }
  *ret_k = kmax;
  return smax;
}


/* B(i) is reached from N(i) or J(i). */
static inline int
select_b(const P7_OPROFILE *om, const P7_OMX *ox, int i)
{
  float path[2];

  path[0] = ( (om->xf[p7O_N][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_N]);
  path[1] = ( (om->xf[p7O_J][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_J]);
  return  ((path[0] > path[1]) ? p7T_N : p7T_J);
}
/*---------------------- end, OA traceback ----------------------*/


