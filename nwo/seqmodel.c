/* Building profile HMMs from single sequence queries.
 *
 * Contents:
 *   1. h4_seqmodel()
 *   2. Example.
 */
#include <h4_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_composition.h"
#include "esl_dmatrix.h"
#include "esl_scorematrix.h"
#include "esl_vectorops.h"

#include "h4_profile.h"

#include "standardize.h"
#include "vectorize.h"

/*****************************************************************
 * 1. h4_seqmodel()
 *****************************************************************/

/* Function:  h4_seqmodel()
 * Synopsis:  Make a profile from a single sequence
 * Incept:    SRE, Fri 23 Jul 2021 [Johnny Cash, One]
 *
 * Purpose:   Create a profile from a single query sequence <dsq> of
 *            length <L> in alphabet <abc>; return it in
 *            <*ret_hmm>.
 *            
 *            This is currently a draft implementation with no choice
 *            of scoring system. The scoring system is BLOSUM62,
 *            reversed to inferred conditional probabilities using the
 *            BLOSUM62 background frequencies. Insert/delete
 *            transitions are set using a gap-open probability of 0.02
 *            (tMI, tMD) and a gap-extend probability of 0.40 (tII,
 *            tDD). ID and DI transitions are set to 0.0. All this
 *            needs to be optimized against benchmarks in the future.
 *
 *            The <h4_SINGLE> flag is raised on the <hmm>. Model
 *            configuration (h4_standardize()) detects this flag.
 *            <L->Mk> local entry transitions use a match state
 *            occupancy term for profiles, but for single sequence
 *            queries, that term is assumed 1.0 for all positions.
 *
 * Args:      dsq     - single query sequence
 *            L       - length of <dsq>
 *            abc     - digital alphabet of <dsq>
 *            ret_hmm - RESULT: new HMM, allocated here
 *            errbuf  - optional: space for informative mesg on failure (or NULL)
 *
 * Returns:   <eslOK> on success. <*ret_hmm> points to the new profile.
 *            <errbuf>, if provided, is an empty string.
 *
 * Throws:    <eslEMEM> on allocation failure
 *            <eslEINVAL> if score matrix can't be reverse-engineered to probabilities
 *            <eslENOHALT> if reverse-engineering failed to converge
 */
int
h4_seqmodel(const ESL_DSQ *dsq, int L, const ESL_ALPHABET *abc, H4_PROFILE **ret_hmm, char *errbuf)
{
  H4_PROFILE      *hmm = NULL;
  ESL_SCOREMATRIX *S   = esl_scorematrix_Create(abc);
  double          *fa  = NULL;
  ESL_DMATRIX     *P   = NULL;
  double           lambda;
  float            popen   = 0.02;   // TODO: optimize, and provide for setting
  float            pextend = 0.4;    //  (ditto)
  int              a,b,k;
  int              status;

  if (errbuf) errbuf[0] = '\0';

  if ((    hmm = h4_profile_Create(abc, L))          == NULL)  { status = eslEMEM; goto ERROR; }
  if (( status = esl_scorematrix_Set("BLOSUM62", S)) != eslOK) goto ERROR;

  /* Use BLOSUM62 bg frequencies (since we're using the BLOSUM62 matrix),
   * and marginalize degeneracies 
   */
  ESL_ALLOC(fa, sizeof(double) * abc->Kp);
  esl_composition_BL62(fa);
  for (a = abc->K+1; a < abc->Kp-3; a++)
    {
      fa[a] = 0.;
      for (b = 0; b < abc->K; b++)
        if (abc->degen[a][b]) fa[a] += fa[b];
    }
  fa[abc->K]    = 0.;  // gap, -
  fa[abc->Kp-3] = 1.0; // any, X or N
  fa[abc->Kp-2] = 0.0; // none, *
  fa[abc->Kp-1] = 0.0; // missing, ~
  esl_vec_D2F(fa, abc->K, hmm->f);
      
  /* Backcalculate P(a,b) joints from the score matrix,
   * and convert to conditionals P(b|a) = P(ab) / P(a) 
   */
  status = esl_scorematrix_ProbifyGivenBG(S, fa, fa, &lambda, &P); 
  if      (status == eslEINVAL)  ESL_XEXCEPTION(eslEINVAL,  errbuf, "built-in score matrix BLOSUM62 has no valid solution for lambda");
  else if (status == eslENOHALT) ESL_XEXCEPTION(eslENOHALT, errbuf, "failed to solve score matrix BLOSUM62 for lambda");
  else if (status != eslOK)      ESL_XEXCEPTION(status,     errbuf, "unexpected error in solving score matrix BLOSUM62 for probability parameters");

  if (( status = esl_scorematrix_JointToConditionalOnQuery(abc, P)) != eslOK)
    ESL_XEXCEPTION(status, errbuf, "joint to conditional failed");
  // now P->mx[a][b] = P(b | a), for residue a in the query, b in the target seq

  /* Parameterize the profile model of targets, given the query seq */
  for (k = 1; k <= L; k++)
    esl_vec_D2F(P->mx[(int) dsq[k]], abc->K, hmm->e[k]);

  hmm->t[0][h4_TMM] = 1.0 - popen;
  hmm->t[0][h4_TMD] = popen;
  for (k = 1; k < L; k++)
    {
      hmm->t[k][h4_TMM] = 1.0 - 2*popen;
      hmm->t[k][h4_TMI] = popen;
      hmm->t[k][h4_TMD] = popen;

      hmm->t[k][h4_TIM] = 1.0 - pextend;
      hmm->t[k][h4_TII] = pextend;
      hmm->t[k][h4_TID] = 0.0;

      hmm->t[k][h4_TDM] = 1.0 - pextend;
      hmm->t[k][h4_TDI] = 0.0;
      hmm->t[k][h4_TDD] = pextend;
    }
  // h4_profile_Create already set remaining fixed boundary conditions on t[0], t[M] and e[0].

  hmm->flags |= (h4_SINGLE | h4_HASPROBS);
  if (( status = h4_standardize(hmm))            != eslOK) goto ERROR;
  if (( status = h4_vectorize(hmm))              != eslOK) goto ERROR;

  free(fa);
  esl_dmatrix_Destroy(P);
  esl_scorematrix_Destroy(S);
  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  free(fa);
  esl_dmatrix_Destroy(P);
  esl_scorematrix_Destroy(S);
  h4_profile_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}


/*****************************************************************
 * 2. Example
 *****************************************************************/
#ifdef h4SEQMODEL_EXAMPLE

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "h4_profile.h"

#include "general.h"
#include "seqmodel.h"

static ESL_OPTIONS options[] = {
  /* name          type       default  env  range toggles reqs incomp  help         docgroup */
 { "-h",          eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help",         0 },
 { "--dna",       eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "use DNA alphabet",        0 },
 { "--rna",       eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "use RNA alphabet",        0 },
 { "--amino",     eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "use amino alphabet",      0 },
 { "--informat",  eslARG_STRING, NULL, NULL, NULL, NULL, NULL, NULL, "set input format",        0 },
 { "--version",   eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show HMMER version info", 0 },
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <seqfile>";
static char banner[] = "example of building profiles from single sequences";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go        = h4_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *seqfile   = esl_opt_GetArg(go, 1);
  int             infmt     = eslSQFILE_UNKNOWN;
  int             alphatype = eslUNKNOWN;
  ESL_SQFILE     *sqfp      = NULL;
  ESL_ALPHABET   *abc       = NULL;
  ESL_SQ         *sq        = NULL;
  H4_PROFILE     *hmm       = NULL;
  char            errbuf[eslERRBUFSIZE];
  int             status;

  if (esl_opt_IsOn(go, "--informat")) {
    if ((infmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat")))==eslSQFILE_UNKNOWN)
      esl_fatal("%s is not a valid input sequence file format for --informat"); 
  }

  status = esl_sqfile_Open(seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) esl_fatal("No such file.");
  else if (status == eslEFORMAT)   esl_fatal("Format couldn't be determined.");
  else if (status != eslOK)        esl_fatal("Open failed, code %d.", status);

  if      (esl_opt_GetBoolean(go, "--rna"))   alphatype = eslRNA;
  else if (esl_opt_GetBoolean(go, "--dna"))   alphatype = eslDNA;
  else if (esl_opt_GetBoolean(go, "--amino")) alphatype = eslAMINO;
  else {
    status = esl_sqfile_GuessAlphabet(sqfp, &alphatype);
    if      (status == eslENOALPHABET)  esl_fatal("Couldn't guess alphabet");
    else if (status == eslEFORMAT)      esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));     
    else if (status == eslENODATA)      esl_fatal("Sequence file empty?");
    else if (status != eslOK)           esl_fatal("Unexpected error guessing alphabet");
  }
  abc = esl_alphabet_Create(alphatype);
  sq  = esl_sq_CreateDigital(abc);
  esl_sqfile_SetDigital(sqfp, abc);

 while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {  
      if ( h4_seqmodel(sq->dsq, sq->n, abc, &hmm, errbuf) != eslOK) esl_fatal("h4_seqmodel failed\n  %s\n",      errbuf);
      if ( h4_profile_Validate(hmm, errbuf)               != eslOK) esl_fatal("model validation failed\n  %s\n", errbuf);

      h4_profile_Dump(stdout, hmm);

      esl_sq_Reuse(sq);
      h4_profile_Destroy(hmm);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed\n  %s", esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d in reading", status);
 
  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif //h4SEQMODEL_EXAMPLE

