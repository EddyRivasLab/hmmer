/* The Plan7 MSVDATA data structure, which holds a compact representation
 * of substitution scores and maximal extensions, used by nhmmer.
 *
 * Contents:
 *   1. The P7_MSVDATA object: allocation, initialization, destruction.
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
 * Synopsis:  Get compact representation of substitution scores and maximal extensions
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
p7_hmm_GetScoreArrays(P7_OPROFILE *om, P7_MSVDATA *data, int do_opt_ext ) {
  int i, j, status;

  //gather values from gm->rsc into a succinct 2D array
  uint8_t *max_scores;
  float sc_fwd, sc_rev;
  int K = om->abc->Kp;
  data->M = om->M;

  ESL_ALLOC(data->scores, (om->M + 1) * K * sizeof(uint8_t));
  p7_oprofile_GetMSVEmissionArray(om, data->scores);

  if (do_opt_ext) {
    ESL_ALLOC(max_scores, (om->M + 1) * sizeof(float));
    for (i = 1; i <= om->M; i++) {
      max_scores[i] = 0;
      for (j=0; j<K; j++)
        if (data->scores[i*K + j] > max_scores[i]) max_scores[i] = data->scores[i*K + j];
    }

    //for each position in the query, what's the highest possible score achieved by extending X positions, for X=1..10
    ESL_ALLOC(data->opt_ext_fwd, (om->M + 1) * sizeof(uint8_t*));
    ESL_ALLOC(data->opt_ext_rev, (om->M + 1) * sizeof(uint8_t*));

    for (i=1; i<=om->M; i++) {
      ESL_ALLOC(data->opt_ext_fwd[i], 10 * sizeof(uint8_t));
      ESL_ALLOC(data->opt_ext_rev[i], 10 * sizeof(uint8_t));
      sc_fwd = 0;
      sc_rev = 0;
      for (j=0; j<10 && i+j+1<=om->M; j++) {
        sc_fwd += max_scores[i+j+1];
        data->opt_ext_fwd[i][j] = sc_fwd;

        sc_rev += max_scores[om->M-i-j];
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
 * Synopsis:  Create a <P7_MSVDATA> model object, based on MSV-filter
 *            part of profile
 *
 * Purpose:   Allocate a <P7_MSVDATA> object used in both FM-MSV and
 *            MSV_LongTarget diagonal recovery/extension, then populate
 *            it with data based on the given gm.
 *
 *            Once a hit passes the MSV filter, and the prefix/suffix
 *            values of P7_MSVDATA are required, p7_hmm_MSVDataComputeRest()
 *            must be called.
 *
 * Args:      om         - P7_PROFILE containing scores used to produce MSVDATA contents
 *            do_opt_ext - boolean, TRUE if optimal-extension scores are required (for FM-MSV)
 *
 * Returns:   a pointer to the new <P7_MSVDATA> object.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_MSVDATA *
p7_hmm_MSVDataCreate(P7_OPROFILE *om, int do_opt_ext )
{
  P7_MSVDATA *data = NULL;
  int    status;

  ESL_ALLOC(data, sizeof(P7_MSVDATA));

  data->scores         = NULL;
  data->opt_ext_fwd    = NULL;
  data->opt_ext_rev    = NULL;
  data->prefix_lengths = NULL;
  data->suffix_lengths = NULL;

  p7_hmm_GetScoreArrays(om, data, do_opt_ext); /* for FM-index string tree traversal */

  return data;

ERROR:
 p7_hmm_MSVDataDestroy(data);
 return NULL;
}


/* Function:  p7_hmm_MSVDataClone()
 * Synopsis:  Clone a <P7_MSVDATA> model object
 *
 * Purpose:   Allocate a <P7_MSVDATA> object used in both FM-MSV and
 *            MSV_LongTarget diagonal recovery/extension, then
 *            copy data into it from another populated instance
 *
 *            Once a hit passes the MSV filter, and the prefix/suffix
 *            values of P7_MSVDATA are required, p7_hmm_MSVDataComputeRest()
 *            must be called.
 *
 * Args:      src        - P7_MSVDATA upon which clone will be based
 *            Kp          - alphabet size, including degeneracy codes, gaps
 *
 * Returns:   a pointer to the new <P7_MSVDATA> object.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_MSVDATA *
p7_hmm_MSVDataClone(P7_MSVDATA *src, int Kp) {
  P7_MSVDATA *new;
  int status;
  int i;

  if (src == NULL)
    return NULL;

  ESL_ALLOC(new, sizeof(P7_MSVDATA));
  new->M = src->M;
  new->scores         = NULL;
  new->opt_ext_fwd    = NULL;
  new->opt_ext_rev    = NULL;
  new->prefix_lengths = NULL;
  new->suffix_lengths = NULL;

  if (src->scores != NULL) {
    ESL_ALLOC(new->scores, (src->M + 1) * Kp * sizeof(uint8_t));
    memcpy(new->scores, src->scores, (src->M + 1) * Kp * sizeof(uint8_t)  );
  }

  if (src->opt_ext_fwd != NULL) {
     ESL_ALLOC(new->opt_ext_fwd, (src->M + 1) * sizeof(uint8_t*));
     for (i=1; i<=src->M; i++) {
       ESL_ALLOC(new->opt_ext_fwd[i], 10 * sizeof(uint8_t));
       memcpy(new->opt_ext_fwd+i, src->opt_ext_fwd+i, 10 * sizeof(uint8_t));
     }
  }
  if (src->opt_ext_rev != NULL) {
     ESL_ALLOC(new->opt_ext_rev, (src->M + 1) * sizeof(uint8_t*));
     for (i=1; i<=src->M; i++) {
       ESL_ALLOC(new->opt_ext_rev[i], 10 * sizeof(uint8_t));
       memcpy(new->opt_ext_rev+i, src->opt_ext_rev+i, 10 * sizeof(uint8_t));
     }
  }


  if (src->prefix_lengths != NULL) {
     ESL_ALLOC(new->prefix_lengths, (src->M+1) * sizeof(float));
     memcpy(new->prefix_lengths, src->prefix_lengths, (src->M+1) * sizeof(float));
  }
  if (src->suffix_lengths != NULL) {
     ESL_ALLOC(new->suffix_lengths, (src->M+1) * sizeof(float));
     memcpy(new->suffix_lengths, src->suffix_lengths, (src->M+1) * sizeof(float));
  }

  return new;

ERROR:
  return NULL;
}

/* Function:  p7_hmm_MSVDataComputeRest()
 * Synopsis:  Compute prefix_ and suffix_ lengths for a <P7_MSVDATA>
 *            model object that was created by p7_hmm_MSVDataCreate.
 *
 * Purpose:   This approach of computing the prefix/suffix length, used
 *            in establishing windows around a seed diagonal, is fast
 *            becuse it uses a simple closed-form computation of the
 *            length L_i for each position i at which all but
 *            (1-p7_DEFAULT_WINDOW_BETA) of position i's match- and
 *            insert-state emissions are length L_i or shorter.
 *
 * Args:      om         - P7_OPROFILE containing transitions probabilities used to compute (pre|suf)fix lengths
 *            data       - P7_MSVDATA into which the computed prefix/suffix values are placed
 *
 * Returns:   eslEMEM on failure, else eslOK
 *
 * Throws:    <NULL> on allocation failure.
 */
int
p7_hmm_MSVDataComputeRest(P7_OPROFILE *om, P7_MSVDATA *data )
{
  int    status;
  int i;
  float sum;


  float *t_mis;
  float *t_iis;


  ESL_ALLOC(t_mis, sizeof(float) * (om->M+1));
  p7_oprofile_GetFwdTransitionArray(om, p7O_MI, t_mis );

  ESL_ALLOC(t_iis, sizeof(float) * (om->M+1));
  p7_oprofile_GetFwdTransitionArray(om, p7O_II, t_iis );


  ESL_ALLOC(data->prefix_lengths, (om->M+1) * sizeof(float));
  ESL_ALLOC(data->suffix_lengths, (om->M+1) * sizeof(float));


  sum = 0;
  for (i=1; i < om->M; i++) {
    data->prefix_lengths[i] = 2 + (int)(log(p7_DEFAULT_WINDOW_BETA / t_mis[i] )/log(t_iis[i]));
    sum += data->prefix_lengths[i];
  }
  data->prefix_lengths[0] = data->prefix_lengths[om->M] = 0;


  for (i=1; i < om->M; i++)
    data->prefix_lengths[i] /=  sum;

  data->suffix_lengths[om->M] = data->prefix_lengths[om->M-1];
  for (i=om->M - 1; i >= 1; i--)
    data->suffix_lengths[i] = data->suffix_lengths[i+1] + data->prefix_lengths[i-1];
  for (i=2; i < om->M; i++)
    data->prefix_lengths[i] += data->prefix_lengths[i-1];

  free(t_mis);
  free(t_iis);

  return eslOK;

  ERROR:
   if (t_mis) free(t_mis);
   if (t_iis) free(t_iis);
   p7_hmm_MSVDataDestroy(data);
   return eslEMEM;
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

  uint8_t scale = 3.0 / eslCONST_LOG2;                    /* scores in units of third-bits */
  uint8_t bias;
  int x;
  float max = 0.0;

  if ( (abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal(msg);

  if (  p7_hmm_Sample(r, 100, abc, &hmm)        != eslOK) esl_fatal(msg);
  if (  (gm = p7_profile_Create (hmm->M, abc))  == NULL ) esl_fatal(msg);
  if (  (om = p7_oprofile_Create(hmm->M, abc))  == NULL ) esl_fatal(msg);

  for (x = 0; x < gm->abc->K; x++)  max = ESL_MAX(max, esl_vec_FMax(gm->rsc[x], (gm->M+1)*2));
  //based on unbiased_byteify
  max  = -1.0f * roundf(scale * max);
  bias   = (max > 255.) ? 255 : (uint8_t) max;


  if (  (msvdata = p7_hmm_MSVDataCreate(om, FALSE))  == NULL ) esl_fatal(msg);

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


