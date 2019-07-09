/* H4_MODE: "algorithm-dependent" parameters of HMMER4 profiles
 * 
 * Contents:
 *   1. The H4_MODE object
 *   2. Debugging and development tools
 *   3. Unit tests
 *   4. Test driver
 *
 * See also:
 *   h4_profile : the profile itself; "model-dependent" parameters and metadata.
 */
#include "h4_config.h"

#include "easel.h"

#include "simdvec.h"
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


/* Function:  h4_mode_Clone()
 * Synopsis:  Clone an H4_MODE
 * Incept:    SRE, Sun 07 Jul 2019
 */
H4_MODE *
h4_mode_Clone(const H4_MODE *mo)
{
  H4_MODE *m2 = NULL;

  if (( m2 = h4_mode_Create()) == NULL)  goto ERROR;
  if ( h4_mode_Copy(mo, m2)    != eslOK) goto ERROR;
  return m2;

 ERROR:
  h4_mode_Destroy(m2);
  return NULL;
}


/* Function:  h4_mode_Copy()
 * Synopsis:  Copy an H4_MODE.
 * Incept:    SRE, Sun 07 Jul 2019
 */
int
h4_mode_Copy(const H4_MODE *mo, H4_MODE *m2)
{
  int s,t;

  for (s = 0; s < h4_NX; s++)
    for (t = 0; t < h4_NXT; t++)
      {
	m2->xsc[s][t] = mo->xsc[s][t];
	m2->xw[s][t]  = mo->xw[s][t];
      }

  m2->nullsc  = mo->nullsc;
  m2->L       = mo->L;
  m2->nj      = mo->nj;
  m2->pglocal = mo->pglocal;
  return eslOK;
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
  int s,t;

  /* defaults. */
  if (L       == -1)  L       = 400;
  if (nj      == -1.) nj      = 1.;
  if (pglocal == -1.) pglocal = 0.5;

  /* contract checks */
  ESL_DASSERT1(( L >= 0 && L <= 100000 ));
  ESL_DASSERT1(( nj >= 0. ));
  ESL_DASSERT1(( pglocal >= 0. && pglocal <= 1. ));

  mo->xsc[h4_E][h4_LOOP] = esl_log2f(nj / (nj + 1.));
  mo->xsc[h4_E][h4_MOVE] = esl_log2f(1. / (nj + 1.));

  mo->xsc[h4_B][h4_LOOP] = esl_log2f(1. - pglocal);
  mo->xsc[h4_B][h4_MOVE] = esl_log2f(pglocal);

  mo->nj      = nj;
  mo->pglocal = pglocal;
  h4_mode_SetLength(mo, L); // mo->nj must be set before calling SetLength().

  /* Vectorized 16bit scaled scores.
   * See notes on the 3 nat approximation, for why N/J/C transitions are hardcoded zero 
   * VF assumes local alignment. It doesn't even use the B transitions, but set them anyway. 
   */
  for (s = 0; s < h4_NX; s++)
    for (t = 0; t < h4_NXT; t++)
      mo->xw[s][t] = h4_simdvec_wordify(mo->xsc[s][t]);
  mo->xw[h4_N][h4_LOOP] = mo->xw[h4_J][h4_LOOP] = mo->xw[h4_C][h4_LOOP] = 0;  
  mo->xw[h4_B][h4_LOOP]  = 0.;         // B->L (VF and other filters are local-only)
  mo->xw[h4_B][h4_MOVE]  = -32768;     // B->G 

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
 *            | <h4_mode_SetUnihit()>    | unihit glocal/local           |   0 | 0.5     | 400 |
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
int h4_mode_SetUnihit    (H4_MODE *mo){ return h4_mode_SetCustom(mo, 400,  0.0, 0.5); }
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
  float p    = (float) L / (float) (L+1);
  float q    = 1.0f / (float) (L+1);
  mo->nullsc = (float) L * log2f(p) + log2f(q);

  mo->xsc[h4_N][h4_LOOP] = mo->xsc[h4_J][h4_LOOP] = mo->xsc[h4_C][h4_LOOP] = esl_log2f( (float) L    / ( (float) L + 2. + mo->nj));
  mo->xsc[h4_N][h4_MOVE] = mo->xsc[h4_J][h4_MOVE] = mo->xsc[h4_C][h4_MOVE] = esl_log2f((2. + mo->nj) / ( (float) L + 2. + mo->nj));

  mo->xw[h4_N][h4_LOOP]  = mo->xw[h4_J][h4_LOOP]  = mo->xw[h4_C][h4_LOOP]  = 0;
  mo->xw[h4_N][h4_MOVE]  = mo->xw[h4_J][h4_MOVE]  = mo->xw[h4_C][h4_MOVE]  = h4_simdvec_wordify(mo->xsc[h4_N][h4_MOVE]);

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
 * 2. Debugging and development tools
 *****************************************************************/

/* Function:  h4_mode_Dump()
 * Synopsis:  Dump contents of an H4_PATH for inspection.
 */
int
h4_mode_Dump(FILE *fp, const H4_MODE *mo)
{
  printf("E->J: %9.5f    E->C: %9.5f\n", mo->xsc[0][0], mo->xsc[0][1]);
  printf("N->N: %9.5f    N->B: %9.5f\n", mo->xsc[1][0], mo->xsc[1][1]);
  printf("J->J: %9.5f    J->B: %9.5f\n", mo->xsc[2][0], mo->xsc[2][1]); 
  printf("C->C: %9.5f    C->T: %9.5f\n", mo->xsc[3][0], mo->xsc[3][1]); 
  printf("B->L: %9.5f    B->G: %9.5f\n", mo->xsc[4][0], mo->xsc[4][1]); 
  return eslOK;
}


/* Function:  h4_mode_SameAsVF()
 * Synopsis:  Scale and round to match VF.
 * Incept:    SRE, Sun 07 Jul 2019
 *
 * Purpose:   Make a copy of <mo> with standard <xsc> scores set to
 *            match the scaling and rounding of the Viterbi filter:
 *            scaled by <h4_SCALE_W> and rounded to nearest integer.
 *            
 *            With <h4_profile_SameAsVF()>, a unit test (for example)
 *            can obtain a (local-only alignment) model that will give
 *            identical results for reference Viterbi compared to the
 *            Viterbi filter, so long as the comparison doesn't
 *            overflow in VF. See <h4_profile_SameAsVF()> for further
 *            notes on how to do this.
 *            
 *            The new, scaled/rounded <H4_MODE> is returned in
 *            <*ret_xmo>. Caller frees, with <h4_mode_Destroy()>.  You
 *            can't call any reparameterization on <xmo>, including
 *            <h4_mode_SetLength()>, because such routines would set
 *            standard (unscaled, unrounded) scores.
 */
int
h4_mode_SameAsVF(const H4_MODE *mo, H4_MODE **ret_xmo)
{
  H4_MODE *xmo = NULL;
  int      s,t;
  int      status;

  if ((xmo = h4_mode_Clone(mo)) == NULL) { status = eslEMEM; goto ERROR; }
  
  for (s = 0; s < h4_NX; s++)
    for (t = 0; t < h4_NXT; t++)
      xmo->xsc[s][t] = (float) xmo->xw[s][t];

  *ret_xmo = xmo;
  return eslOK;

 ERROR:
  h4_mode_Destroy(xmo);
  *ret_xmo = NULL;
  return status;
}



/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef h4MODE_TESTDRIVE

static void
utest_sanity(void)
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
 * 4. Test driver
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

  utest_sanity();

  fprintf(stderr, "#  status = ok\n");
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /* h4MODE_TESTDRIVE */



