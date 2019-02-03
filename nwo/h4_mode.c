/* H4_MODE: "algorithm-dependent" parameters of HMMER4 profiles
 * 
 * Contents:
 *   1. The H4_MODE object.
 *   2. Unit tests
 *   3. Test driver
 *
 * See also:
 *   h4_profile : the profile itself; "model-dependent" parameters and metadata.
 */
#include "h4_config.h"

#include "easel.h"

#include "h4_mode.h"


/*****************************************************************
 * 1. The H4_MODE object
 *****************************************************************/

/* Function:  h4_mode_Create()
 * Synopsis:  Creates H4_MODE, algorithm-dependent profile params.
 * Incept:    SRE, Wed 30 Jan 2019
 *
 * Purpose:   Create an <H4_MODE> object with default values for
 *            algorithm-dependent profile parameters, and return
 *            a ptr to it.
 *            
 *            The default alignment mode is multihit dual-mode
 *            (glocal/local).  To set a different standard mode, see
 *            <h4_mode_SetGlocal()> and friends. To set a
 *            custom-parameterized mode, see <h4_mode_SetCustom()>.
 *            
 *            For search statistics to work well, you need to set the
 *            length parameterization to the actual length of each
 *            individual target sequence before scoring or aligning to
 *            it. See <h4_mode_SetLength()>. The object is initialized
 *            with <L=400> as a "typical" protein length.
 *
 * Returns:   ptr to new <H4_MODE>.
 *
 * Throws:    <NULL> on allocation failure.
 */
H4_MODE *
h4_mode_Create(void)
{
  H4_MODE *mo = NULL;
  int      status;

  ESL_ALLOC(mo, sizeof(H4_MODE));
  h4_mode_SetCustom(mo, -1, -1., -1.);  // requests defaults for L, nj, pglocal 
  return mo;

 ERROR:
  return NULL;
}


/* Function:  h4_mode_SetCustom()
 * Synopsis:  Set algorithm-dependent parameters for a search profile.
 * Incept:    SRE, Thu 31 Jan 2019
 *
 * Purpose:   Set the algorithm-dependent parameters in the <H4_MODE>
 *            structure that accompanies a profile/seq comparison:
 *            target sequence length <L>, expected additional number
 *            of domains per alignment <nj>, and glocal (vs. local)
 *            alignment probability <pglocal>.
 *
 *            Passing <L=-1>, <nj=-1.0>, and/or <pglocal=-1.0> sets
 *            default values, which are currently <L=400>, <nj=1.0>,
 *            <pglocal=0.5>.
 *
 * Returns:   <eslOK> on success.
 */
int
h4_mode_SetCustom(H4_MODE *mo, int L, float nj, float pglocal)
{
  /* defaults. */
  if (L       == -1)  L       = 400;
  if (nj      == -1.) nj      = 1.;
  if (pglocal == -1.) pglocal = 0.5;

  /* contract checks */
  ESL_DASSERT1(( L >= 0 && L <= 100000 ));
  ESL_DASSERT1(( nj >= 0. ));
  ESL_DASSERT1(( pglocal >= 0. && pglocal <= 1. ));

  mo->xsc[h4_E][h4_LOOP] = esl_logf(nj / (nj + 1.));
  mo->xsc[h4_E][h4_MOVE] = esl_logf(1. / (nj + 1.));

  mo->xsc[h4_B][h4_LOOP] = 1. - pglocal;
  mo->xsc[h4_B][h4_MOVE] = pglocal;

  h4_mode_SetLength(mo, L);

  mo->nj      = nj;
  mo->pglocal = pglocal;
  return eslOK;
}


/* Function:  h4_mode_SetDefault(), h4_mode_SetLocal(), h4_mode_SetGlocal(), h4_mode_SetUnilocal(), h4_mode_SetUniglocal(), h4_mode_SetGlobal()
 * Synopsis:  Set alignment-dependent parameters to a standard mode
 * Incept:    SRE, Thu 31 Jan 2019
 *
 * Purpose:   Set the algorithm parameters in the <H4_MODE> structure to
 *            one of several standard modes:
 *            
 *            |  function                | mode                          |  nj | pglocal | L   | 
 *            |--------------------------|-------------------------------|-----|---------|-----|
 *            | <h4_mode_SetDefault()>   | default multihit glocal/local | 1.0 | 0.5     | 400 |
 *            | <h4_mode_SetLocal()>     | multihit local-only           | 1.0 |   0     | 400 |
 *            | <h4_mode_SetGlocal()>    | multihit glocal-only          | 1.0 | 1.0     | 400 |
 *            | <h4_mode_SetUnilocal()>  | unihit local-only             |   0 |   0     | 400 |
 *            | <h4_mode_SetUniglocal()> | unihit glocal-only            |   0 | 1.0     | 400 |
 *            | <h4_mode_SetGlobal()>    | single global alignment       |   0 | 1.0     |   0 |
 *
 *            You'll need to reset the <L> parameter to each actual
 *            target sequence length with <h4_mode_SetLength()>, so
 *            the default <L> parameter setting is just a placeholder;
 *            except for <_SetGlobal()>, where the <L=0> setting
 *            prevents any nonhomologous flanking residues at all. 
 *
 * Returns:   <eslOK> on success.
 */
int h4_mode_SetDefault   (H4_MODE *mo){ return h4_mode_SetCustom(mo, 400,  1.0, 0.5); }
int h4_mode_SetLocal     (H4_MODE *mo){ return h4_mode_SetCustom(mo, 400,  1.0, 0.0); }
int h4_mode_SetGlocal    (H4_MODE *mo){ return h4_mode_SetCustom(mo, 400,  1.0, 1.0); }
int h4_mode_SetUnilocal  (H4_MODE *mo){ return h4_mode_SetCustom(mo, 400,  0.0, 0.0); }
int h4_mode_SetUniglocal (H4_MODE *mo){ return h4_mode_SetCustom(mo, 400,  0.0, 1.0); }
int h4_mode_SetGlobal    (H4_MODE *mo){ return h4_mode_SetCustom(mo,  0.,  0.0, 1.0); }


/* Function:  h4_mode_SetLength()
 * Synopsis:  Set the length parameterization for a profile search.
 * Incept:    SRE, Fri 01 Feb 2019
 *
 * Purpose:   For each individual target sequence, you use
 *            <h4_mode_SetLength()> to set the <L> parameterization
 *            in the alignment-dependent parameters.
 *            
 *            This gets called _outside_ the DP algorithms, because
 *            the caller might want to configure <L=0> for example.
 *            (Don't be tempted to put <_SetLength()> calls inside DP
 *            algorithms, even if it looks like that's a good idea
 *            since the DP algorithms have target seq length <L> as an
 *            argument.)
 */
int
h4_mode_SetLength(H4_MODE *mo, int L)
{
  ESL_DASSERT1(( L >= 0 && L <= 100000 ));

  mo->xsc[h4_N][h4_LOOP] = mo->xsc[h4_J][h4_LOOP] = mo->xsc[h4_C][h4_LOOP] = esl_logf( (float) L    / ( (float) L + 2. + mo->nj));
  mo->xsc[h4_N][h4_MOVE] = mo->xsc[h4_J][h4_MOVE] = mo->xsc[h4_C][h4_MOVE] = esl_logf((2. + mo->nj) / ( (float) L + 2. + mo->nj));
  mo->L = L;
  return eslOK;
}

/* Function:  h4_mode_Destroy()
 * Synopsis:  Frees an <H4_MODE>
 */
void
h4_mode_Destroy(H4_MODE *mo)
{
  free(mo);
}


/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef h4MODE_TESTDRIVE

static void
utest_minimal(void)
{
  char     msg[] = "h4_mode minimal unit test failed";
  H4_MODE *mo    = NULL;

  if ((mo = h4_mode_Create())               == NULL)  esl_fatal(msg);
  if ( h4_mode_SetCustom(mo, 400, 1.0, 0.5) != eslOK) esl_fatal(msg);
  if ( h4_mode_SetDefault(mo)               != eslOK) esl_fatal(msg);
  if ( h4_mode_SetLocal(mo)                 != eslOK) esl_fatal(msg);
  if ( h4_mode_SetGlocal(mo)                != eslOK) esl_fatal(msg);
  if ( h4_mode_SetUnilocal(mo)              != eslOK) esl_fatal(msg);
  if ( h4_mode_SetUniglocal(mo)             != eslOK) esl_fatal(msg);
  if ( h4_mode_SetGlobal(mo)                != eslOK) esl_fatal(msg);
  if ( h4_mode_SetLength(mo, 400)           != eslOK) esl_fatal(msg);
  h4_mode_Destroy(mo);
}
#endif /* h4MODE_TESTDRIVE */

/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef h4MODE_TESTDRIVE

#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                             docgroup*/
  { "-h",  eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for h4_mode";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
 
  fprintf(stderr, "## %s\n", argv[0]);

  utest_minimal();

  fprintf(stderr, "#  status = ok\n");
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /* h4MODE_TESTDRIVE */



