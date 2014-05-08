/* Reference implementation of alignment tracebacks for anchor set
 * constrained alignments.
 * 
 * All reference implementation code is for development and
 * testing. It is not used in HMMER's main executables. Production
 * code uses sparse dynamic programming.
 * 
 * Contents:
 *    1.
 *    x. Copyright and license information.
 */

#include "easel.h"
#include "esl_vectorops.h"


/*****************************************************************
 * 1. Choice selection functions for MEG traces
 *****************************************************************/

static inline int
meg_select_m(const P7_PROFILE *mx, int i, int k)
{
  ESL_UNUSED(rng);
  int   state[4] = { p7T_ML, p7T_IL, p7T_DL, p7T_B };
  float path[4];

  path[0] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_ML), P7P_TSC(gm, k-1, p7P_MM));
  path[1] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_IL), P7P_TSC(gm, k-1, p7P_IM));
  path[2] = P7_DELTAT( P7R_MX(mx, i-1, k-1, p7R_DL), P7P_TSC(gm, k-1, p7P_DM));
  path[3] = P7_DELTAT( P7R_XMX(mx, i-1, p7R_B),      P7P_TSC(gm, k-1, p7P_GM));
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
meg_select_c(ESL_RANDOMNESS *rng, const P7_PROFILE *gm,...)
{
  
  if (env->flags & p7E_IS_GLOCAL)
    {

    }

  for (k = env->arr[d].k0; k <= M; k++)
    { 
      if (dpc[p7R_ML] > max) { max = dpc[p7R_ML]; smax = p7T_
  



static int
reference_asc_trace_engine()
{
  int i    = mx->L;
  int k    = 0;
  int sprv = p7T_C;
  int scur;
  int status;

  if ((status = p7_trace_Append(tr, p7T_T, k, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, k, i)) != eslOK) return status;

  for (d = env->n-1; d >= 0; d--)
    {

      /* From previous ia-1 (L) to ib : we know they're C/J */
      scur = sprv; /* C|J; tricky */
      k    = 0;
      for ( ; i >= env[d].ib; i--)
	if ((status = p7_trace_Append(tr, scur,  k, i)) != eslOK) return status;

      /* now i=ib for this domain, and we've already added C|J for ib.
       * sprv is C|J.
       */
      while (sprv != p7T_B)
	{
	  switch (sprv) {
	    
	  case p7T_C : 
	  default: ESL_EXCEPTION(eslEINCONCEIVABLE, "lost in asc traceback");
	  }

	  if (scur == p7T_B) {
	    if (env[d].flags & p7E_IS_GLOCAL) {
	      while (k) {
		if  ( (status = p7_trace_Append(tr, p7T_DG, k, i)) != eslOK) return status;
		k--;
	      }
	      if    ( (status = p7_trace_Append(tr, p7T_G, k, i)) != eslOK) return status;
	    }
	    else if ( (status = p7_trace_Append(tr, p7T_L, k, i)) != eslOK) return status;
	  }
	}

	      



      
      /* now you're  */
      



    }

}


/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
