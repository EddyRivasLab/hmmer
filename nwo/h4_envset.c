/* H4_ENVSET
 *
 * Contents:
 *   1. the H4_ENVSET object: a set of envelopes
 * 
 */
#include "h4_config.h"

#include "easel.h"

#include "h4_anchorset.h"
#include "h4_envset.h"


/*****************************************************************
 * 1. H4_ENVSET object
 *****************************************************************/


/* Function:  h4_envset_Create()
 * Synopsis:  Create a new <H4_ENVSET>
 * Incept:    SRE, Thu 25 Feb 2021
 *
 * Purpose:   Returns ptr to a new <H4_ENVSET> object for defining <D> domains, for a
 *            comparison of a sequence of length <L> to a profile of length <M>.
 *
 *            <D> can be 0, to create an initially empty envelope set.
 *
 *            If <D> is 0, <L> and <M> can also be passed as 0, if they aren't known
 *            yet. In this case, caller will need to call <h4_envset_SetSentinels()>
 *            when it starts using the structure for a particular seq/profile
 *            comparison.
 */
H4_ENVSET *
h4_envset_Create(int D, int L, int M)
{
  H4_ENVSET *env             = NULL;
  int        default_nalloc  = 8;     // i.e. for up to D=6, +2 sentinels
  int        default_redline = 64;    // i.e. for up to D=62
  int        status;

  ESL_DASSERT1(( ( D >= 0 && L > 0 && M > 0 ) ||
                 ( D == 0 && L == 0 && M == 0 ) ));

  ESL_ALLOC(env, sizeof(H4_ENVSET));
  env->e        = NULL;
  env->D        = D;
  env->L        = L;
  env->M        = M;
  env->nalloc   = esl_resize(ESL_MAX(D+2, default_nalloc), 0, default_redline);
  env->nredline = default_redline;

  ESL_ALLOC(env->e, sizeof(H4_ENVELOPE) * env->nalloc);
  return env;

 ERROR:
  h4_envset_Destroy(env);
  return NULL;
}


/* Function:  h4_envset_Resize()
 * Synopsis:  Reallocate an H4_ENVSET to hold at least D domains.
 * Incept:    SRE, Wed 17 Mar 2021
 */
int
h4_envset_Resize(H4_ENVSET *env, int D)
{
  int nalloc = esl_resize(D+2, env->nalloc, env->nredline);  // +2 for sentinels
  int status;

  if (nalloc != env->nalloc) {
    ESL_REALLOC(env->e, sizeof(H4_ENVELOPE) * nalloc);
    env->nalloc = nalloc;
  }
  return eslOK;

 ERROR:
  return status;
}


/* Function:  h4_envset_SetSentinels()
 * Synopsis:  Set the sentinels in an envelope set
 * Incept:    SRE, Sun 11 Jul 2021 [Josh Ritter, Another New World]
 */
int
h4_envset_SetSentinels(H4_ENVSET *env, int D, int L, int M)
{
  env->e[0].oa = env->e[0].ia = env->e[0].i0 = env->e[0].ib = env->e[0].ob = 0;
  env->e[0].ka = env->e[0].k0 = env->e[0].kb = M+1;

  env->e[D+1].oa = env->e[D+1].ia = env->e[D+1].i0 = env->e[D+1].ib = env->e[D+1].ob = L+1;
  env->e[D+1].ka = env->e[D+1].k0 = env->e[D+1].kb = 0;

  env->e[0].env_sc   = env->e[D+1].env_sc   = 0.;
  env->e[0].null2_sc = env->e[D+1].null2_sc = 0.;  
  env->e[0].flags    = env->e[D+1].flags    = 0;

  env->D = D;
  env->L = L;
  env->M = M;
  return eslOK;
}


/* Function:  h4_envset_CopyFromAnchorset()
 * Synopsis:  Use an anchor set to initialize an H4_ENVSET
 * Incept:    SRE, Wed 17 Mar 2021
 *
 * Purpose:   Use anchorset <anch> to initialize <env>, including sentinels at d=0 and
 *            d=D+1.  <env> must already be allocated for at least D domains; caller
 *            should use <h4_envset_Resize()> first if you need to.
 */
int
h4_envset_CopyFromAnchorset(const H4_ANCHORSET *anch, H4_ENVSET *env)
{
  int d;
  int status;

  if ((status = h4_anchorset_GetSentinels(anch, &(env->L), &(env->M))) != eslOK) return status;

  for (d = 1; d <= anch->D; d++)  
    {
      env->e[d].i0 = anch->a[d].i0;
      env->e[d].k0 = anch->a[d].k0;

      env->e[d].ia = env->e[d].ib = 0;
      env->e[d].ka = env->e[d].kb = 0;
      env->e[d].oa = env->e[d].ob = 0;

      env->e[d].env_sc   = 0.;
      env->e[d].null2_sc = 0.;
      env->e[d].flags    = 0;
    }
  h4_envset_SetSentinels(env, anch->D, env->L, env->M);
  return eslOK;
}


/* Function:  h4_envset_Dump()
 * Synopsis:  Dump contents of an H4_ENVSET for inspection
 * Incept:    SRE, Wed 17 Mar 2021
 */
int
h4_envset_Dump(FILE *ofp, const H4_ENVSET *env)
{
  int d;

  esl_dataheader(ofp, 5, "dom",
                 5, "oa", 5, "ia", 5, "i0", 5, "ib", 5, "ob", 5, "ka", 5, "k0", 5, "kb", 
		 6, "env_sc", 
		 3, "app", 3, "glo",
		 0);

  for (d = 1; d <= env->D; d++)
    {
      fprintf(ofp, "%-5d %5d %5d %5d %5d %5d %5d %5d %5d %6.2f %3s %3s\n",
              d,
              env->e[d].oa, env->e[d].ia, env->e[d].i0, env->e[d].ib, env->e[d].ob,
              env->e[d].ka, env->e[d].k0, env->e[d].kb,
              env->e[d].env_sc,
              (env->e[d].flags & h4E_ENVSC_APPROX) ? "YES" : "no",
              (env->e[d].flags & h4E_IS_GLOCAL)    ? "YES" : "no");
    }
              
  fprintf(ofp, "# Total domains = %d\n",    env->D);
  fprintf(ofp, "# M,L           = %d,%d\n", env->M, env->L);
  fprintf(ofp, "# nalloc        = %d\n",    env->nalloc);
  fprintf(ofp, "# nredline      = %d\n",    env->nredline);
  return eslOK;
}


/* Function:  h4_envset_Destroy()
 * Synopsis:  Free an H4_ENVSET
 * Incept:    SRE, Wed 17 Mar 2021
 */
void
h4_envset_Destroy(H4_ENVSET *env)
{
  if (env) {
    free(env->e);
    free(env);
  }
}
