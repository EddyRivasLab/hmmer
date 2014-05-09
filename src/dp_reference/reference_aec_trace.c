/* Reference implementation of alignment tracebacks for anchor/envelope
 * constrained (AEC) alignment.
 * 
 * All reference implementation code is for development and
 * testing. It is not used in HMMER's main executables. Production
 * code uses sparse dynamic programming.
 * 
 * Style here is taken from reference_trace.c. The choice functions
 * are abstracted, so the same traceback engine works on any
 * optimization criterion for the AEC path. (Even though the only
 * currently implemented optimization is MEG.)
 * 
 * Someday we want to replace P7_TRACE with a more compressed data
 * structure.
 * 
 * Contents:
 *    1. Choice selection for MEG paths
 *    2. Traceback engine
 *    3. Exposed API - wrappers around the engine
 *    4. Copyright and license information.
 */

#include "easel.h"
#include "esl_vectorops.h"

#include "base/p7_envelopes.h"
#include "base/p7_profile.h"
#include "base/p7_trace.h"

#include "dp_reference/p7_refmx.h"
#include "dp_reference/reference_aec_trace.h"


/*****************************************************************
 * 1. Choice selection functions for MEG traces
 *****************************************************************/

static inline int
aec_meg_select_ml(const P7_PROFILE *gm, const P7_REFMX *mx, int i, int k)
{
  int   state[4] = { p7T_ML, p7T_IL, p7T_DL, p7T_L };
  float path[4];

  path[0] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_ML), P7P_TSC(gm, k-1, p7P_MM));
  path[1] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_IL), P7P_TSC(gm, k-1, p7P_IM));
  path[2] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_DL), P7P_TSC(gm, k-1, p7P_DM));
  path[3] = P7_DELTAT( P7R_XMX(mx, i-1, p7R_L),      P7P_TSC(gm, k-1, p7P_LM));
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
aec_meg_select_mg(const P7_PROFILE *gm, const P7_REFMX *mx, int i, int k)
{
  int   state[4] = { p7T_MG, p7T_IG, p7T_DG, p7T_G };
  float path[4];

  path[0] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_MG), P7P_TSC(gm, k-1, p7P_MM));
  path[1] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_IG), P7P_TSC(gm, k-1, p7P_IM));
  path[2] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_DG), P7P_TSC(gm, k-1, p7P_DM));
  path[3] = P7_DELTAT( P7R_XMX(mx, i-1, p7R_L),      P7P_TSC(gm, k-1, p7P_LM));
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
aec_meg_select_il(const P7_PROFILE *gm, const P7_REFMX *mx, int i, int k)
{
  float path[2];

  path[0] = P7_DELTAT( P7R_MX(mx, i-1, k, p7R_ML), P7P_TSC(gm, k, p7P_MI));
  path[1] = P7_DELTAT( P7R_MX(mx, i-1, k, p7R_IL), P7P_TSC(gm, k, p7P_II));
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_IL);
}

static inline int
aec_meg_select_ig(const P7_PROFILE *gm, const P7_REFMX *mx, int i, int k)
{
  float path[2];

  path[0] = P7_DELTAT( P7R_MX(mx, i-1, k, p7R_MG), P7P_TSC(gm, k, p7P_MI));
  path[1] = P7_DELTAT( P7R_MX(mx, i-1, k, p7R_IG), P7P_TSC(gm, k, p7P_II));
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_IG);
}

static inline int
aec_meg_select_dl(const P7_PROFILE *gm, const P7_REFMX *mx, int i, int k)
{
  float path[2];

  path[0] = P7_DELTAT( P7R_MX(mx, i, k-1, p7R_ML), P7P_TSC(gm, k-1, p7P_MD));
  path[1] = P7_DELTAT( P7R_MX(mx, i, k-1, p7R_DL), P7P_TSC(gm, k-1, p7P_DD));
  return ( (path[0] >= path[1]) ? p7T_ML : p7T_DL);
}

static inline int
aec_meg_select_dg(const P7_PROFILE *gm, const P7_REFMX *mx, int i, int k)
{
  float path[2];

  path[0] = P7_DELTAT( P7R_MX(mx, i, k-1, p7R_MG), P7P_TSC(gm, k-1, p7P_MD));
  path[1] = P7_DELTAT( P7R_MX(mx, i, k-1, p7R_DG), P7P_TSC(gm, k-1, p7P_DD));
  return ( (path[0] >= path[1]) ? p7T_MG : p7T_DG);
}

static inline int
aec_meg_select_e(const P7_PROFILE *gm, const P7_REFMX *mx, int i, int *ret_k, const P7_ENVELOPE *env)
{
  float *dpc = mx->dp[i] + env->k0*p7R_NSCELLS; /* position at k0 */
  float max  = -eslINFINITY;
  int   smax = -1;
  int   kmax = -1;
  int   k;

  if ( env->flags & p7E_IS_GLOCAL)
    {
      dpc  = mx->dp[i] + gm->M * p7R_NSCELLS;
      kmax = gm->M;
      smax = P7R_MX(mx, i, gm->M, p7R_MG) >= P7R_MX(mx, i, gm->M, p7R_DG) ? p7T_MG : p7T_DG;
    }
  else
    {
      for (k = env->k0; k <= gm->M; k++)
	{ 
	  if (dpc[p7R_ML] > max) { max = dpc[p7R_ML]; smax = p7T_ML; kmax = k; }   
	  if (dpc[p7R_DL] > max) { max = dpc[p7R_DL]; smax = p7T_DL; kmax = k; }   
	  dpc += p7R_NSCELLS;
	}
    }

  *ret_k = kmax;
  return smax;
}

/* The J and C selection functions need <apd>, because MEG optimization
 * (using a label gain function) uses posterior probabilities J(i), C(i).
 */
static inline int
aec_meg_select_j(const P7_PROFILE *gm, const P7_REFMX *mx, int i, const P7_REFMX *apd)
{
  float path[2];

  path[0] = P7_DELTAT ( P7R_XMX(mx, i-1, p7R_J), gm->xsc[p7P_J][p7P_LOOP]) + P7R_XMX(apd, i, p7R_JJ);
  path[1] = P7_DELTAT ( P7R_XMX(mx, i,   p7R_E), gm->xsc[p7P_E][p7P_LOOP]);
  return ( (path[0] > path[1]) ? p7T_J : p7T_E);
}
  
static inline int
aec_meg_select_c(const P7_PROFILE *gm, const P7_REFMX *mx, int i, const P7_REFMX *apd)
{
  float path[2];

  path[0] = P7_DELTAT ( P7R_XMX(mx, i-1, p7R_C), gm->xsc[p7P_C][p7P_LOOP]) + P7R_XMX(apd, i, p7R_CC);
  path[1] = P7_DELTAT ( P7R_XMX(mx, i,   p7R_E), gm->xsc[p7P_E][p7P_MOVE]);
  return ( (path[0] > path[1]) ? p7T_C : p7T_E);
}

/* We don't need a select_b() function for AEC, because AEC knows what domain
 * it's in at any given time, so it knows whether a B traces back to J vs. N
 */
/*------------ end, MEG selection functions ---------------------*/




/*****************************************************************
 * 2. Traceback engine.
 *****************************************************************/

/* reference_aec_trace_engine()
 * 
 * Adapted from reference_trace.c. 
 */
static int
reference_aec_trace_engine(const P7_PROFILE *gm, P7_ENVELOPES *env, const P7_REFMX *apd, const P7_REFMX *mx, P7_TRACE *tr)
{
  int i    = mx->L;
  int k    = 0;
  int d;
  int sprv = p7T_C;
  int scur;
  int status;
  /* abstracting the per-cell traceback selection functions lets us 
   * use one engine for any optimization criterion:
   */
  int (*select_ml)(const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_mg)(const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_il)(const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_ig)(const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_dl)(const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_dg)(const P7_PROFILE *, const P7_REFMX *, int, int);
  int (*select_e) (const P7_PROFILE *, const P7_REFMX *, int, int *, const P7_ENVELOPE *);
  int (*select_j) (const P7_PROFILE *, const P7_REFMX *, int, const P7_REFMX *);
  int (*select_c) (const P7_PROFILE *, const P7_REFMX *, int, const P7_REFMX *);


  /* Contract checks, argument validation */
  ESL_DASSERT1(( gm->M == mx->M && gm->M == apd->M ));
  ESL_DASSERT1(( mx->L == apd->L ));
  ESL_DASSERT1(( apd->type == p7R_ASC_DECODE_DOWN ));

  /* Currently we only need traceback for MEG alignments, but as in
   * reference_trace.c, design allows for generalization to other
   * optimization criteria.
   */
  if (mx->type == p7R_AEC_MEG_ALIGN)
    {
      select_ml = &aec_meg_select_ml;   select_mg = &aec_meg_select_mg;
      select_il = &aec_meg_select_il;   select_ig = &aec_meg_select_ig;
      select_dl = &aec_meg_select_dl;   select_dg = &aec_meg_select_dg;
      select_j  = &aec_meg_select_j;    select_c  = &aec_meg_select_c;
      select_e  = &aec_meg_select_e;
    }
  else esl_fatal("can't trace that matrix type; not AEC?");


  /* Now the main recursion.
   * The code below could be simpler. The standard engine
   * (reference_trace.c) will work fine with just one exception: we
   * must be careful on rows i==ia that we don't use select_{ml,mg},
   * because we cannot look at row ia-1; it might be the ib row of a
   * prev domain d-1. However, the code below also breaks out the
   * flanking and interdomain N/C/J chunks, because I think it will
   * eventually prove advantageous to use a new trace structure (instead
   * of P7_TRACE) that does not bother to represent N/C/J paths
   * explicitly. 
   */
  if   ((status = p7_trace_Append(tr, p7T_T, 0, i)) != eslOK) return status;
  if   ((status = p7_trace_Append(tr, p7T_C, 0, i)) != eslOK) return status; 
  for (i = mx->L; i > env->arr[env->n-1].ib; i--)
    if ((status = p7_trace_Append(tr, p7T_C, 0, i)) != eslOK) return status;
  /* Now i=ib for last domain D-1, in the C state. */

  for (d = env->n-1; d >= 0; d--)
    {
      /* now i=ib for this domain, and we've already added C|J for ib.
       * sprv is C|J. Work backward until we get to the B state.
       */
      while (sprv != p7T_B)
	{
	  /* on rows i=ia, you cannot use select_{ml,mg}, because you cannot
	   * look at row i-1; that could be an ib row of prev domain d-1
	   */
	  if      (i == env->arr[d].ia && sprv == p7T_ML) { scur = p7T_L; k--; i--; } // k-- because we have to act like select_{ml,mg}; ali endpoints depend on it
	  else if (i == env->arr[d].ia && sprv == p7T_MG) { scur = p7T_G; k--; i--; }
	  else {
	    switch (sprv) {
	    case p7T_ML: scur = (*select_ml)(gm, mx, i, k);        k--; i--; break;
	    case p7T_MG: scur = (*select_mg)(gm, mx, i, k);        k--; i--; break;
	    case p7T_IL: scur = (*select_il)(gm, mx, i, k);             i--; break;
	    case p7T_IG: scur = (*select_ig)(gm, mx, i, k);             i--; break;
	    case p7T_DL: scur = (*select_dl)(gm, mx, i, k);        k--;      break;
	    case p7T_DG: scur = (*select_dg)(gm, mx, i, k);        k--;      break;
	    case p7T_E:  scur = (*select_e) (gm, mx, i, &k, &(env->arr[d])); break;
	    case p7T_J:  scur = (*select_j) (gm, mx, i, apd);                break;
	    case p7T_L:  scur = p7T_B;                                       break;
	    case p7T_G:  scur = p7T_B;                                       break;
	    case p7T_C:  scur = (*select_c) (gm, mx, i, apd);                break;
	    default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in aec traceback");
	    }
	  }

	  /* Unfold the left wing */
	  if (scur == p7T_G) {  
	    while (k > 0) {
	      if  ( (status = p7_trace_Append(tr, p7T_DG, k, i)) != eslOK) return status;
	      k--;
	    }
	  }

	  /* Record start, end points of alignment in <env>.
	   * This must follow left wing unfolding, to get glocal ka=1 right
           * ka setting is a bit subtle. Takes advantage of the fact that switch() 
           * above has decremented k; also, that glocal wing unfolding leaves k=0.
	   */
	  if (sprv == p7T_E) { env->arr[d].alib = i;   env->arr[d].kb   = k;   }
	  if (scur == p7T_B) { env->arr[d].alia = i+1; env->arr[d].ka   = k+1; }

	  if ( (status = p7_trace_Append(tr, scur, k, i)) != eslOK) return status;
	  /* For NCJ, we had to defer i decrement. */
	  if ( (scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;
	  sprv = scur;
	}

      /* Now, sprv = B, i == -1 from first i of alignment; >= ia-1 of envelope.
       * Residues i back to ib+1 for prev domain (or 1) are generated by J (or N).
       * This might be no residues; it's possible i==ib (or i==0) already.
       * We always add at least one J/N state, for i==ib (or i==0) itself.
       * In reversed traces (like what we have now), the first J we add is 
       * a nonemitter.
       */
      if (d > 0)
	{
	  if ((status = p7_trace_Append(tr, p7T_J, 0, i)) != eslOK) return status;  // always at least one J added, as a nonemitter
	  for ( ; i > env->arr[d-1].ib; i--)
	    if ((status = p7_trace_Append(tr, p7T_J, 0, i)) != eslOK) return status;
	  /* Now we're on row ib[d-1], in the J state. */
	  sprv = p7T_J;
	}
    }
  
  if ((status = p7_trace_Append(tr, p7T_N, 0, i)) != eslOK) return status;
  for ( ; i >= 1; i--)
    if ((status = p7_trace_Append(tr, p7T_N, 0, i)) != eslOK) return status;
  if   ((status = p7_trace_Append(tr, p7T_S, 0, 0)) != eslOK) return status;

  tr->M = mx->M;
  tr->L = mx->L;
  return p7_trace_Reverse(tr);
}
/*------------- end, AEC traceback engine -----------------------*/




/*****************************************************************
 * 3. Exposed API - wrappers around the engine
 *****************************************************************/

int
p7_reference_aec_trace_MEG(const P7_PROFILE *gm, P7_ENVELOPES *env, const P7_REFMX *apd, const P7_REFMX *mx, P7_TRACE *tr)
{
  return (reference_aec_trace_engine(gm, env, apd, mx, tr));
}

/*-------------- end, exposed API wrappers ----------------------*/





/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
