/* H4_PATHIDX - index of the locations of domains in a path
 *
 * See h4_pathidx.md for notes.
 *
 * Contents:
 *    1. The H4_PATHIDX object
 *    2. Debugging & development tools
 *    3. Unit tests
 *    4. Test driver
 */
#include <h4_config.h>

#include "easel.h"

#include "h4_path.h"
#include "h4_pathidx.h"


/*****************************************************************
 * 1. H4_PATHIDX 
 *****************************************************************/

/* Function:  h4_pathidx_Build()
 * Synopsis:  Build an index of domain locations in a path.
 * Incept:    SRE, Wed 13 Mar 2024
 *
 * Purpose:   Given path <pi>, build an index of the locations of each
 *            of its domains; return the index in <*ret_pidx>.
 *
 *            In the edge case of a zero-length homology path [xref
 *            h4_path.md], returns an empty index with D=0; only the L
 *            field (target sequence length) is valid.
 *
 * Returns:   <eslOK> on success; <*ret_pidx> is the new index.
 *
 * Throws:    <eslEMEM> on allocation failure; now <*ret_pidx> is NULL.
 */
int
h4_pathidx_Build(const H4_PATH *pi, H4_PATHIDX **ret_pidx)
{
  H4_PATHIDX *pidx = NULL;
  int         D    = h4_path_GetDomainCount(pi);
  int         d,z,i,k;
  int         is_zerolen;
  int         status;

  ESL_ALLOC(pidx, sizeof(H4_PATHIDX));
  pidx->ia = NULL;

  ESL_ALLOC(pidx->ia, sizeof(int) * 7*(D+1));
  pidx->ib        = pidx->ia +   (D+1);
  pidx->ka        = pidx->ia + 2*(D+1);
  pidx->kb        = pidx->ia + 3*(D+1);
  pidx->za        = pidx->ia + 4*(D+1);
  pidx->zb        = pidx->ia + 5*(D+1);
  pidx->is_glocal = pidx->ia + 6*(D+1);

  pidx->ia[0] = pidx->ib[0] = 0;  // d=0 is unused, because we index as d=1..D
  pidx->ka[0] = pidx->kb[0] = 0;
  pidx->za[0] = pidx->zb[0] = 0;
  pidx->is_glocal[0]        = FALSE;
  pidx->D                   = D;

  for (z = 0, d = 0, i = 0; z < pi->Z; z++)
    {
      if (h4_path_IsB(pi->st[z])) // L|G. Start of a domain (barring zero-length homology, anyway)
        {
          if (pi->rle[z]) {
            pidx->ia[++d]      = i+1;             // i = current position in seq. The domain is going to start at i+1.
            pidx->ka[d]        = k = pi->rle[z];  // k=1 for G; k=ka for L. This is the *next* k we visit.
            pidx->za[d]        = z;
            pidx->is_glocal[d] = (pi->st[z] == h4P_G ? TRUE : FALSE );
            is_zerolen         = FALSE;
          } else
            is_zerolen         = TRUE;
        }
      else if (pi->st[z] == h4P_C || pi->st[z] == h4P_J)
        {
          if (! is_zerolen) {
            pidx->ib[d] = i;              // i = last i we generated, so that was ib
            pidx->kb[d] = k-1;            // since k=next k we use, k-1 is the last one we used.
            pidx->zb[d] = z-1;         
          }
          i += pi->rle[z]-1;            // the -1 is because C|J emit on transition.
                                        // we still count i even for a zero-length homology.
        }
      else if (pi->st[z] == h4P_N)     { i += pi->rle[z]-1;                   }
      else if (h4_path_IsM(pi->st[z])) { i += pi->rle[z];    k += pi->rle[z]; }
      else if (h4_path_IsI(pi->st[z])) { i += pi->rle[z];                     }
      else if (h4_path_IsD(pi->st[z])) {                     k += pi->rle[z]; }
    }
  pidx->L = i;

  ESL_DASSERT1(( d == D ));
  *ret_pidx = pidx;
  return eslOK;

 ERROR:
  *ret_pidx = NULL;
  return status;
}

/* Function:  h4_pathidx_Destroy()
 * Synopsis:  Free an H4_PATHIDX
 * Incept:    SRE, Thu 14 Mar 2024
 */
void
h4_pathidx_Destroy(H4_PATHIDX *pidx)
{
  if (pidx) {
    free(pidx->ia);
    free(pidx);
  }
}

/*****************************************************************
 * 2. Development & debugging
 *****************************************************************/

/* Function:  h4_pathidx_Validate()
 * Synopsis:  Validate an H4_PATHIDX object
 * Incept:    SRE, Wed 13 Mar 2024
 *
 * Purpose:   Validate the internal data in a path index structure <pidx>,
 *            that indexes domain coords in path <pi>. 
 *
 * Args:      
 *
 * Returns:   <eslOK> if path index looks fine.
 *            <eslFAIL> if a problem is detected and <errbuf> is set to say why.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */

int
h4_pathidx_Validate(const H4_PATHIDX *pidx, const H4_PATH *pi, char *errbuf)
{
  int M = 0;             // profile length unknown until we see a glocal domain; but then we know and can do some extra checks
  int d;
  int8_t st1, st2, st3;

  if (pidx == NULL) ESL_FAIL(eslFAIL, errbuf, "path index is NULL");
  if (pi   == NULL) ESL_FAIL(eslFAIL, errbuf, "path is NULL");

  /* Have a pass over the domains and see if one is glocal. If so, we know M */
  for (d = 1; d <= pidx->D; d++)
    if (pidx->is_glocal[d]) { M = pidx->kb[d]; break; }

  /* Now do some checks that domain boundaries make sense */
  for (d = 1; d <= pidx->D; d++)
    {
      st1 = pi->st[pidx->za[d]];     // domain should start with L|G.   (And that must be followed by Mk|Dk, but if that isn't so, it's a path problem)
      st2 = pi->st[pidx->zb[d]];     //  ... and end with Mk|Dk
      st3 = pi->st[pidx->zb[d]+1];   //  ... and that Mk|Dk should be the last before J|C

      if (! h4_path_IsB(st1))                       ESL_FAIL(eslFAIL, errbuf, "domain %d doesn't start with L|G",    d);
      if (! h4_path_IsM(st2) && ! h4_path_IsD(st2)) ESL_FAIL(eslFAIL, errbuf, "domain %d doesn't end with Mk|Dk",    d);
      if (st3 != h4P_J && st3 != h4P_C)             ESL_FAIL(eslFAIL, errbuf, "domain %d end isn't followed by J|C", d);

      if (pidx->is_glocal[d]) {
        if (pidx->ka[d] != 1 || pidx->kb[d] != M)   ESL_FAIL(eslFAIL, errbuf, "domain %d glocal but ka..kb not 1..M", d);
      } else {
        if (M > 0 && pidx->kb[d] > M)               ESL_FAIL(eslFAIL, errbuf, "domain %d (local) has kb > M", d);
      }
    }

  if (errbuf) errbuf[0] = '\0';
  return eslOK;
}



/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef h4PATHIDX_TESTDRIVE

#include "esl_random.h"

#include "esl_alphabet.h"
#include "esl_sq.h"

#include "h4_mode.h"
#include "h4_path.h"
#include "h4_pathidx.h"
#include "h4_profile.h"

#include "emit.h"
#include "modelsample.h"

static void
utest_sanity(ESL_RANDOMNESS *rng, ESL_ALPHABET *abc)
{
  char          msg[]   = "h4_pathidx sanity unit test failed";
  ESL_SQ       *sq      = esl_sq_CreateDigital(abc);
  H4_PROFILE   *hmm     = NULL;
  H4_MODE      *mo      = h4_mode_Create();   // default: dual glocal/local
  H4_PATH      *pi      = h4_path_Create();
  H4_PATHIDX   *pidx    = NULL;
  int           M       = 20;
  int           ntrials = 100;
  char          errbuf[eslERRBUFSIZE];

  if ( h4_mode_SetLength(mo, 10)         != eslOK) esl_fatal(msg);
  if ( h4_modelsample(rng, abc, M, &hmm) != eslOK) esl_fatal(msg);

  while (ntrials--)
    {
      if ( h4_emit(rng, hmm, mo, sq, pi) != eslOK) esl_fatal(msg);
      if ( h4_pathidx_Build(pi, &pidx)   != eslOK) esl_fatal(msg);

      if (sq->n != pidx->L) esl_fatal(msg);
      if ( h4_pathidx_Validate(pidx, pi, errbuf) != eslOK) esl_fatal("%s:\n  %s", msg, errbuf);

      h4_pathidx_Destroy(pidx);
    }

  h4_path_Destroy(pi);
  h4_mode_Destroy(mo);
  h4_profile_Destroy(hmm);
  esl_sq_Destroy(sq);
}
#endif // h4PATHIDX_TESTDRIVE




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef h4PATHIDX_TESTDRIVE

#include <h4_config.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "general.h"
#include "h4_path.h"

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

  utest_sanity(rng, abc);

  fprintf(stderr, "#  status   = ok\n");

  esl_randomness_Destroy(rng);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  exit(0);
}

#endif // h4PATHIDX_TESTDRIVE
