/* The Plan7 core HMM data structure.
 *
 * Contents:
 *   1. The P7_HMM object: allocation, initialization, destruction.
 *   2. Unit tests.
 *   3. Test driver.
 *   4. Copyright and license.
 *
 */
#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"
#include "esl_random.h"
#include "esl_dirichlet.h"

#include "hmmer.h"


/*********************************************************************
 *# 1. The P7_MSVDATA object: allocation, initialization, destruction.
 *********************************************************************/

/* Function:  p7_hmm_GetScoreArrays()
 * Synopsis:  Get compact representation of substitution scores and maximal extenstions
 *
 * Purpose:   Extract 8-bit (MSV-style) substitution scores from gm, and
 *            for each position in the model capture the maximum possible
 *            score that can be added to a diagonal's score (in both
 *            directions) by extending lengths 1..10.
 *
 *            The scores will be used in both standard MSV diagonal recovery
 *            and FM-MSV diagonal scoring. The extension scores are used in
 *            FM-MSV's pruning step.
 *
 * Args:      gm         - P7_PROFILE containing scores used to produce MSVDATA contents
 *            data       - where scores and will be stored
 *            do_opt_ext - boolean, TRUE if optimal-extension scores are required (for FM-MSV)
 *            scale      - used to produce 8-bit extracted scores
 *            bias       - used to produce 8-bit extracted scores
 *
 * Returns:   data->scores and possibly ->opt_ext_(fwd|rev) are filled in;
 *            return eslEMEM on allocation failure, eslOK otherwise.
 */
static int
p7_hmm_GetScoreArrays(P7_PROFILE *gm, P7_MSVDATA *data, int do_opt_ext, float scale, uint bias ) {
  int i, j, status;

  //gather values from gm->rsc into a succinct 2D array
  float *max_scores;
  float sc_fwd, sc_rev;

  int K = gm->abc->Kp;

  data->M = gm->M;

  ESL_ALLOC(data->scores, (gm->M + 1) * K * sizeof(float));
  ESL_ALLOC(max_scores, (gm->M + 1) * sizeof(float));

  for (i = 1; i <= gm->M; i++) {
    max_scores[i] = 0;
    for (j=0; j<K; j++) {
      //based on p7_oprofile's biased_byteify()
      float x =  -1.0f * roundf(scale * gm->rsc[j][(i) * p7P_NR     + p7P_M]);
      data->scores[i*K + j] = (x > 255. - bias) ? 255 : (uint8_t) (x + bias);
      if (data->scores[i*K + j] > max_scores[i]) max_scores[i] = data->scores[i*K + j];
    }
  }

  if (do_opt_ext) {
    //for each position in the query, what's the highest possible score achieved by extending X positions, for X=1..10
    ESL_ALLOC(data->opt_ext_fwd, (gm->M + 1) * sizeof(uint8_t*));
    ESL_ALLOC(data->opt_ext_rev, (gm->M + 1) * sizeof(uint8_t*));

    for (i=1; i<=gm->M; i++) {
      ESL_ALLOC(data->opt_ext_fwd[i], 10 * sizeof(uint8_t));
      ESL_ALLOC(data->opt_ext_rev[i], 10 * sizeof(uint8_t));
      sc_fwd = 0;
      sc_rev = 0;
      for (j=0; j<10 && i+j+1<=gm->M; j++) {
        sc_fwd += max_scores[i+j+1];
        data->opt_ext_fwd[i][j] = sc_fwd;

        sc_rev += max_scores[gm->M-i-j];
        data->opt_ext_rev[i][j] = sc_rev;
      }
      for ( ; j<10; j++) { //fill in empty values
        data->opt_ext_fwd[i][j] = data->opt_ext_fwd[i][j-1];
        data->opt_ext_rev[i][j] = data->opt_ext_rev[i][j-1];
      }

    }
  }

  return eslOK;

ERROR:
  return eslEMEM;
}


/* Function:  p7_hmm_MSVDataDestroy()
 *
 * Synopsis:  Destroy a <P7_MSVDATA> object.
 */
void
p7_hmm_MSVDataDestroy(P7_MSVDATA *data )
{
  int i;
  if (data != NULL) {
    if (data->scores != NULL)  free( data->scores);
    if (data->prefix_lengths != NULL) free( data->prefix_lengths);
    if (data->suffix_lengths != NULL) free( data->suffix_lengths);

    if (data->opt_ext_fwd != NULL) {
      for (i=1; i<=data->M; i++)
        free(data->opt_ext_fwd[i]);
      free(data->opt_ext_fwd);
    }
    if (data->opt_ext_rev != NULL) {
      for (i=1; i<=data->M; i++)
        free(data->opt_ext_rev[i]);
      free( data->opt_ext_rev);
    }

    free(data);
  }

}

/* Function:  p7_hmm_MSVDataCreate()
 * Synopsis:  Create a <P7_MSVDATA> model object.
 *
 * Purpose:   Allocate a <P7_MSVDATA> object used in both FM-MSV and
 *            MSV_LongTarget diagonal recovery/extension, then populate
 *            it with data based on the given gm and hmm.
 *
 *            This approach of computing the prefix/suffix length, used
 *            in establishing windows around a seed diagonal, is very fast
 *            by virtue of using a simple closed-form computation of the
 *            length L_i for each position i at which all but
 *            (1-p7_DEFAULT_WINDOW_BETA) of position i's match- and
 *            insert-state emissions are length L_i or shorter.
 *
 * Args:      gm         - P7_PROFILE containing scores used to produce MSVDATA contents
 *            hmm        - P7_HMM containing transitions probabilities used to compute (pre|suf)fix lengths
 *            do_opt_ext - boolean, TRUE if optimal-extension scores are required (for FM-MSV)
 *            scale      - used to produce 8-bit extracted scores
 *            bias       - used to produce 8-bit extracted scores
 *
 * Returns:   a pointer to the new <P7_MSVDATA> object.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_MSVDATA *
p7_hmm_MSVDataCreate(P7_PROFILE *gm, P7_HMM *hmm, int do_opt_ext, float scale, int bias )
{
  P7_MSVDATA *data = NULL;
  int    status;
  int i;
  float sum;

  ESL_ALLOC(data, sizeof(P7_MSVDATA));
  ESL_ALLOC(data->prefix_lengths, (gm->M+1) * sizeof(float));
  ESL_ALLOC(data->suffix_lengths, (gm->M+1) * sizeof(float));

  data->scores         = NULL;
  data->opt_ext_fwd    = NULL;
  data->opt_ext_rev    = NULL;

  p7_hmm_GetScoreArrays(gm, data, do_opt_ext, scale, bias ); /* for FM-index string tree traversal */

  sum = 0;
  for (i=1; i < gm->M; i++) {
    //TODO: this requires that the full HMM be read before performing MSV. If that's
    // a problem for speed, e.g. for hmmscan, I'll need to preprocess these prefix/suffix
    // values, and store them in the hmm file
    data->prefix_lengths[i] = 2 + (int)(log(p7_DEFAULT_WINDOW_BETA / hmm->t[i][p7H_MI] )/log(hmm->t[i][p7H_II]));
    sum += data->prefix_lengths[i];
  }
  for (i=1; i < gm->M; i++)
    data->prefix_lengths[i] /=  sum;

  data->suffix_lengths[gm->M] = data->prefix_lengths[gm->M-1];
  for (i=gm->M - 1; i >= 1; i--)
    data->suffix_lengths[i] = data->suffix_lengths[i+1] + data->prefix_lengths[i-1];
  for (i=2; i < gm->M; i++)
    data->prefix_lengths[i] += data->prefix_lengths[i-1];


  return data;

ERROR:
 p7_hmm_MSVDataDestroy(data);
 return NULL;
}



/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef p7MSVDATA_TESTDRIVE

static void
utest_createMSVData(ESL_GETOPTS *go, ESL_RANDOMNESS *r )
{
  char          msg[]       = "create MSVData unit test failed";
  P7_HMM        *hmm        = NULL;
  ESL_ALPHABET  *abc        = NULL;
  P7_PROFILE    *gm         = NULL;
  P7_OPROFILE   *om         = NULL;
  P7_MSVDATA    *msvdata    = NULL;

  if ( (abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal(msg);

  if (  p7_hmm_Sample(r, 100, abc, &hmm)        != eslOK) esl_fatal(msg);
  if (  (gm = p7_profile_Create (hmm->M, abc))  == NULL ) esl_fatal(msg);
  if (  (om = p7_oprofile_Create(hmm->M, abc))  == NULL ) esl_fatal(msg);

  if (  (msvdata = p7_hmm_MSVDataCreate(gm, hmm, FALSE, om->scale_b, om->bias_b))  == NULL ) esl_fatal(msg);

  p7_hmm_MSVDataDestroy(msvdata);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);

}
#endif /*p7BG_TESTDRIVE*/


/*****************************************************************
 * 3. Test driver
 *****************************************************************/

#ifdef p7MSVDATA_TESTDRIVE
#include "esl_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
   /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                            0},
  {"-s",  eslARG_INT,       "0", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0},
  {"-v",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show verbose commentary/output",                 0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for p7_bg";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go          = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *rng         = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             be_verbose  = esl_opt_GetBoolean(go, "-v");

  if (be_verbose) printf("p7_msvdata unit test: rng seed %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_createMSVData(go, rng);

  esl_randomness_Destroy(rng);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /* p7BG_TESTDRIVE */

/************************************************************
 * @LICENSE@
 *
 * SVN $Id: p7_msvdata.c 3784 2011-12-07 21:51:25Z wheelert $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/src/hmmer/trunk/src/p7_msvdata.c $
 ************************************************************/


