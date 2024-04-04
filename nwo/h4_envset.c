/* H4_ENVSET
 *
 * Contents:
 *   1. The H4_ENVSET object: a set of envelopes
 *   2. Internal (static) functions
 *   3. Unit tests
 *   4. Test driver
 */
#include <h4_config.h>

#include <string.h>

#include "easel.h"

#include "h4_anchorset.h"
#include "h4_envset.h"
#include "h4_mode.h"
#include "h4_path.h"
#include "h4_pathidx.h"
#include "h4_profile.h"

static int build_envpath(const H4_PATH *pi, const H4_PATHIDX *pidx, int whichd, H4_PATH *envpi);

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
 *
 * Purpose:   Reallocate <env> to hold at least <D> domains.
 *
 *            Domains are numbered 1..D, and the <env> structure has
 *            sentinels at 0 and D+1.
 *
 *            The resize obeys standard Easel conventions for a
 *            growing/shrinking structure resize, with a redline
 *            (<env->redline>). 
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


/* Function:  h4_envset_SetFromPath()
 * Synopsis:  Construct envelope set using path (rather than by ensemble inference)
 * Incept:    SRE, Mon 11 Mar 2024 
 *
 * Purpose:   Given path <pi> and anchorset <anch> for a comparison of a
 *            sequence <dsq> of length <L> to profile <hmm> in mode
 *            <mo>, construct an envelope set using that path as a
 *            point estimate: assigns ia..ib and ka..kb to domain
 *            endpoints, outer envelope coords oa=ia and ob=ib (since
 *            there's no uncertainty to represent), and set envelope
 *            scores to partial path bit scores for each domain. Store the
 *            results in <env>, which the caller provides in any
 *            allocation. <env> is reused and internally reallocated
 *            as needed.
 *            
 *            This function is used in experiments comparing ensemble
 *            inference to Viterbi inference. In these experiments,
 *            <pi> would typically be a Viterbi optimal path. We would
 *            not use this function in normal code. Normally we
 *            determine envelopes by ensemble inference.
 *
 *            Caller typically obtains <anch> by a prior call to
 *            <h4_reference_mpas_path2anchors()> or
 *            <h4_sparse_mpas_path2anchors()>, depending on whether it's
 *            using reference or sparse DP.
 *
 *            As a very special case, for unit testing only, <anch>
 *            can be passed as NULL. In this case, anchor coords will
 *            be set to -1/-1. This saves us from having to do DP for
 *            posterior decoding in the envset unit tests below, only
 *            to set anchor coords that aren't part of what we test.
 *
 *            Because we only use this code in experiments, not
 *            production, it is not time critical. To calculate
 *            "envelope scores" for partial paths, we make a path
 *            specific for each individual domain and call
 *            <h4_path_Score()> on it. This may not be the fastest way
 *            to do things, but it seemed to require the least fiddly
 *            code.
 *
 * Args:      pi   : path to use to define envelopes. Typically a Viterbi path.
 *            dsq  : digital target sequence
 *            L    : length of <dsq>
 *            hmm  : profile HMM
 *            mo   : mode for profile/seq comparison that gave <pi> as a result
 *            anch : anchor set defined from <pi> with h4_{reference,sparse}_mpas_path2anchors(), depending on DP impl. (Can be NULL, only for unit tests)
 *            env  : RESULT: envelope set.  Caller provides this in any allocation; reused/resized as needed.
 * 
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEINCONCEIVABLE> on coding errors.
 *            On errors, the state of <env> is undefined; don't use it.
 *
 * Xref:      H15/110.
 *
 * Notes:     If you just need the bounds of each domain in a (Viterbi)
 *            path, use the simpler H4_PATHIDX. The H4_ENVSET
 *            structure is for H4-style ensemble decoding of domain
 *            locations.
 */
int
h4_envset_SetFromPath(const H4_PATH *pi, const ESL_DSQ *dsq, int L, const H4_PROFILE *hmm, const H4_MODE *mo, const H4_ANCHORSET *anch, H4_ENVSET *env)
{
  H4_PATH    *envpi = NULL;  // To get envelope scores, we construct single-domain paths in this structure, then call h4_path_Score()
  H4_PATHIDX *pidx  = NULL;
  int do_envpaths;
  int d;
  int status;

  if (( status = h4_pathidx_Build(pi, &pidx)) != eslOK) return status;

  if (( status = h4_envset_Resize(env, pidx->D))                  != eslOK) return status;  // eslEMEM
  if (( status = h4_envset_SetSentinels(env, pidx->D, L, hmm->M)) != eslOK) return eslEINCONCEIVABLE; 

  /* If there's more than 1 domain in <pi> we'll construct
   * single-domain paths to score to obtain env scores. If there's
   * only 1 domain in the path to begin with, we can just use <pi>
   * itself.
   */
  do_envpaths = (pidx->D > 1 ? TRUE : FALSE);
  if (do_envpaths) {
    if (( envpi  = h4_path_Create())             == NULL)  return eslEMEM;   // A structure big enough to hold original path 
    if (( status = h4_path_Resize(envpi, pi->Z)) != eslOK) return status;    // is big enough to hold any single envelope.
  }

  for (d = 1; d <= pidx->D; d++)
    {
      env->e[d].ia = env->e[d].oa = pidx->ia[d];
      env->e[d].ib = env->e[d].ob = pidx->ib[d];
      env->e[d].ka = pidx->ka[d];
      env->e[d].kb = pidx->kb[d];

      env->e[d].i0 = (anch ? anch->a[d].i0 : -1);   // <anch> can be NULL in a special case used for envset unit testing.
      env->e[d].k0 = (anch ? anch->a[d].k0 : -1);

      if (do_envpaths) {   
        if (( status = build_envpath(pi, pidx, d, envpi))                       != eslOK) return status;
        if (( status = h4_path_Score(envpi, dsq, hmm, mo, &(env->e[d].env_sc))) != eslOK) return status;
      } else {  
        if (( status = h4_path_Score(pi, dsq, hmm, mo, &(env->e[d].env_sc)))    != eslOK) return status;
      }
      
      env->e[d].null2_sc = 0.0;       // TK TK
      env->e[d].flags    = (pidx->is_glocal[d] ? h4E_IS_GLOCAL : 0);
    }

  h4_pathidx_Destroy(pidx);
  h4_path_Destroy(envpi);
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
/*******************  end, H4_ENVSET *****************************/


/*****************************************************************
 * 2. Internal (static) functions
 *****************************************************************/


/* build_envpath()
 *
 * Used for calculating envelope score for a domain, in the off-kilter
 * case where we're using Viterbi inference and a single Viterbi path
 * to define envelopes. We only do this for experiments evaluating how
 * well H4's probabilistic inference works, by comparison to a Viterbi
 * baseline. See h4_envset_SetFromPath().
 *
 * Of a variety of ways of obtaining this that I considered, it seems
 * like the simplest is to construct that path for each domain (in
 * <envpi>), then score it with h4_path_Score().  Writing something
 * like a h4_path_PartialScore() variant of h4_path_Score() itself
 * might be slightly more compute-efficient but seemed much more
 * fiddly in terms of extra code.
 */
static int
build_envpath(const H4_PATH *pi, const H4_PATHIDX *pidx, int whichd, H4_PATH *envpi)
{
  int za = pidx->za[whichd];
  int zb = pidx->zb[whichd];
  int ia = pidx->ia[whichd];
  int ib = pidx->ib[whichd];
  int status;

  if ((status = h4_path_Resize(envpi, zb-za+3)) != eslOK) return status; // N + chunk + C

  envpi->st[0]  = h4P_N;
  envpi->rle[0] = ia;
  memcpy(envpi->st+1,  pi->st+za,  (zb-za+1)*sizeof(int8_t));
  memcpy(envpi->rle+1, pi->rle+za, (zb-za+1)*sizeof(int));
  envpi->st[zb-za+2]  = h4P_C;
  envpi->rle[zb-za+2] = pidx->L - ib + 1;    // L-ib nonhomologous residues + 1 extra C state
  envpi->Z = zb-za+3;              
  return eslOK;
}
/*******************  end, internals *****************************/



/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef h4ENVSET_TESTDRIVE

#include "esl_alphabet.h"
#include "esl_sq.h"

#include "emit.h"
#include "modelsample.h"

static void
utest_SetFromPath(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc)
{
  char          msg[]   = "h4_envset SetFromPath unit test failed";
  ESL_SQ       *sq      = esl_sq_CreateDigital(abc);
  H4_PROFILE   *hmm     = NULL;
  H4_MODE      *mo      = h4_mode_Create();   // default: dual glocal/local
  H4_PATH      *pi      = h4_path_Create();
  H4_PATHIDX   *pidx    = NULL;
  H4_PATH      *envpi   = h4_path_Create();
  H4_ENVSET    *env     = NULL;
  int           M       = 20;
  int           ntrials = 100;
  int           d;
  float         totsc,sc,sc_offset;

  if ( h4_mode_SetLength(mo, 10)         != eslOK) esl_fatal(msg);
  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);

  while (ntrials--)
    {
      if ( h4_emit(rng, hmm, mo, sq, pi)    != eslOK) esl_fatal(msg);
      if ( h4_pathidx_Build(pi, &pidx)      != eslOK) esl_fatal(msg);

      if (( env = h4_envset_Create(pidx->D, sq->n, hmm->M)) == NULL) esl_fatal(msg);

      if ( h4_envset_SetFromPath(pi, sq->dsq, sq->n, hmm, mo, /*anch=*/NULL, env) != eslOK) esl_fatal(msg);  // anch=NULL skips copying of anchor i0/k0 coords into <env>
        
      /* If there's only a single domain in the emitted path, then
       * SetFromPath() didn't call build_envpath() to calculate envelope score;
       * it just scored the path itself. So here we can test build_envpath() 
       * and make sure it gets the same result - indeed, the same path exactly.
       */
      if (pidx->D == 1)
        {
          if ( build_envpath(pi, pidx, /*d=*/1, envpi) != eslOK) esl_fatal(msg);
          if ( h4_path_Compare(pi, envpi)              != eslOK) esl_fatal(msg);
        }
      /* If there's more than one domain, our envelope score test
       * is more arcane: we can calculate the sequence score algebraically
       * from the sum of envelope scores [xref reference_envelopes.md].
       */
      else
        {
          h4_path_Score(pi, sq->dsq, hmm, mo, &sc);                  // sequence score
          totsc = 0.0;
          for (d = 1; d <= pidx->D; d++) totsc += env->e[d].env_sc;  // sum of envelope scores

          sc_offset = (float) (pidx->D - 1) *                        // calculated difference (assumes default multihit <mo>)
            ((float) sq->n * mo->xsc[h4_C][h4_LOOP] + mo->xsc[h4_C][h4_MOVE] - mo->nullsc);

          if ( esl_FCompare(sc, totsc-sc_offset, 0.0, 0.01) != eslOK) esl_fatal(msg);   // allow +/- 0.01 bit absolute difference from roundoff error
        }

      h4_envset_Destroy(env);
      h4_pathidx_Destroy(pidx);
    }

  h4_path_Destroy(envpi);
  h4_path_Destroy(pi);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_sq_Destroy(sq);
}


#endif // h4ENVSET_TESTDRIVE
/*******************  end, unit tests ****************************/



/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef h4ENVSET_TESTDRIVE

#include <h4_config.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "general.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                          docgroup*/
  { "-h",         eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help summary",             0 },
  { "-s",         eslARG_INT,     "0", NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed",    0 },
  { "--version",  eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show HMMER version number",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = h4_CreateDefaultApp(options, 0, argc, argv, "test driver for H4_PATHIDX", "[-options]");
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  fprintf(stderr, "## %s\n", argv[0]);
  fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_SetFromPath(rng, abc);

  fprintf(stderr, "#  status   = ok\n");

  esl_randomness_Destroy(rng);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  exit(0);
}
#endif // h4ENVSET_TESTDRIVE
/*******************  end, test driver ***************************/
