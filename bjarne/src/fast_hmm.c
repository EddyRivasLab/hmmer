#include <stdio.h>
#include <math.h>

#include "fast_hmm.h"


/* Function:  fast_hmm_Create()
 * Synopsis:  Allocates a new HMM.
 * Incept:    BK, Wed Aug  3 12:00:00 2009 [Benasque]
 *
 * Purpose:   Allocates a new two-state HMM modeling
 *            strings in the alphabet <abc>.
 *
 * Returns:   a pointer to the new HMM.
 *
 * Throws:    <NULL> on allocation failure.
 */
FAST_HMM *
fast_hmm_Create(const ESL_ALPHABET *abc)
{
  int       alpha_size = abc->Kp;
  void     *ptr;
  FAST_HMM *hmm = NULL;
  int       status;

  ESL_ALLOC(ptr, sizeof(FAST_HMM) + 15);
  hmm =  (FAST_HMM*) (((size_t) ptr + 15) & ~15);

  ESL_ALLOC(hmm->prob_ptr, 5 * alpha_size * sizeof(__m128) + 15);

  hmm->emission_probs  =  (__m128*) (((size_t) hmm->prob_ptr + 15) & ~15);
  hmm->combined_probs  =  hmm->emission_probs + 2 * alpha_size;
  hmm->start_probs     =  hmm->emission_probs + 4 * alpha_size;
  hmm->alpha_size      =  alpha_size;
  hmm->up_to_date      =  0;
  hmm->ptr             =  ptr;
  hmm->t_s0            =  0.0;
  hmm->t_s1            =  0.0;
  hmm->t_se            =  0.0;
  hmm->t_00            =  0.0;
  hmm->t_01            =  0.0;
  hmm->t_0e            =  0.0;
  hmm->t_10            =  0.0;
  hmm->t_11            =  0.0;
  hmm->t_1e            =  0.0;

  return hmm;

 ERROR:
  fast_hmm_Destroy(hmm);
  return NULL;
}


/* Function:  esl_hmm_Destroy()
 * Synopsis:  Destroys an HMM.
 * Incept:    BK, Wed Aug  3 12:00:00 2009 [Benasque]
 *
 * Purpose:   Frees an HMM.
 */
void
fast_hmm_Destroy(FAST_HMM *hmm)
{
  if (hmm == NULL) return;

  /* Free the original pointers instead of the aligned pointers */
  if (hmm->prob_ptr != NULL) {
    free(hmm->prob_ptr);
  }
  if (hmm->ptr != NULL) {
    free(hmm->ptr);
  }
}


/* Function:  fast_hmm_SetEmissions()
 * Synopsis:  Set emission probabilities for an HMM.
 * Incept:    BK, Wed Aug  3 12:00:00 2009 [Benasque]
 *
 * Purpose:   Set the emission probabilities for the two states in
 *            the HMM. The emission probabilities are given as
 *            an array for each state, <emissions_0> and
 *            <emissions_1>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
fast_hmm_SetEmissions(FAST_HMM *hmm, float *emissions_0, float *emissions_1)
{
  int   i;
  float factor = 0.0;

  for (i = 0; i < hmm->alpha_size; i++)
    {
      if (emissions_0[i] > factor) factor = emissions_0[i];
      if (emissions_1[i] > factor) factor = emissions_1[i];
    }

  for (i = 0; i < hmm->alpha_size; i++)
    {
      ((float*) (hmm->emission_probs + 2 * i    ))[0] = emissions_0[i] / factor;
      ((float*) (hmm->emission_probs + 2 * i    ))[1] = emissions_1[i] / factor;
      ((float*) (hmm->emission_probs + 2 * i    ))[2] = emissions_0[i] / factor;
      ((float*) (hmm->emission_probs + 2 * i    ))[3] = emissions_1[i] / factor;

      ((float*) (hmm->emission_probs + 2 * i + 1))[0] = emissions_0[i] / factor;
      ((float*) (hmm->emission_probs + 2 * i + 1))[1] = emissions_0[i] / factor;
      ((float*) (hmm->emission_probs + 2 * i + 1))[2] = emissions_1[i] / factor;
      ((float*) (hmm->emission_probs + 2 * i + 1))[3] = emissions_1[i] / factor;
    }

  hmm->emission_adj = logf(factor);

  hmm->up_to_date = 0;

  return eslOK;
}


/* Function:  fast_hmm_SetTransitions()
 * Synopsis:  Set the transition probabilities for an HMM.
 * Incept:    BK, Wed Aug  3 12:00:00 2009 [Benasque]
 *
 * Purpose:   Set the transition probabilities for an HMM.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
fast_hmm_SetTransitions(FAST_HMM *hmm, float t_s0, float t_s1, float t_se, float t_00, float t_01, float t_0e, float t_10, float t_11, float t_1e)
{
  float factor = 0.0;

  /* Do not include t_se in factor since it is not used in the vectors */
  if (t_s0 > factor) factor = t_s0;
  if (t_s1 > factor) factor = t_s1;
  if (t_00 > factor) factor = t_00;
  if (t_01 > factor) factor = t_01;
  if (t_0e > factor) factor = t_0e;
  if (t_10 > factor) factor = t_10;
  if (t_11 > factor) factor = t_11;
  if (t_1e > factor) factor = t_1e;

  hmm->t_s0 = t_s0 / factor;
  hmm->t_s1 = t_s1 / factor;
  hmm->t_se = t_se / factor;
  hmm->t_00 = t_00 / factor;
  hmm->t_01 = t_01 / factor;
  hmm->t_0e = t_0e / factor;
  hmm->t_10 = t_10 / factor;
  hmm->t_11 = t_11 / factor;
  hmm->t_1e = t_1e / factor;

  hmm->transition_adj = logf(factor);

  hmm->up_to_date = 0;

  return eslOK;
}


void
set_probs(FAST_HMM *hmm)
{
  int     i;
  __m128  t_start;
  __m128  t_even;
  __m128  t_odd;

  ((float*) &t_start)[0]  =  hmm->t_s0;
  ((float*) &t_start)[1]  =  hmm->t_s1;
  ((float*) &t_start)[2]  =  hmm->t_s0;
  ((float*) &t_start)[3]  =  hmm->t_s1;

  ((float*) &t_even)[0]  =  hmm->t_00;
  ((float*) &t_even)[1]  =  hmm->t_01;
  ((float*) &t_even)[2]  =  hmm->t_10;
  ((float*) &t_even)[3]  =  hmm->t_11;

  ((float*) &t_odd)[0]  =  hmm->t_00;
  ((float*) &t_odd)[1]  =  hmm->t_10;
  ((float*) &t_odd)[2]  =  hmm->t_01;
  ((float*) &t_odd)[3]  =  hmm->t_11;

  ((float*) &(hmm->end_probs))[0] = hmm->t_0e;
  ((float*) &(hmm->end_probs))[1] = 0;
  ((float*) &(hmm->end_probs))[2] = 0;
  ((float*) &(hmm->end_probs))[3] = hmm->t_1e;

  for (i = 0; i < hmm->alpha_size; i++)
    {
      hmm->start_probs   [    i    ] = _mm_mul_ps(hmm->emission_probs[2 * i    ], t_start);

      hmm->combined_probs[2 * i    ] = _mm_mul_ps(hmm->emission_probs[2 * i    ], t_even );
      hmm->combined_probs[2 * i + 1] = _mm_mul_ps(hmm->emission_probs[2 * i + 1], t_odd  );
    }
}


/* Function:  fast_hmm_GetLogProb()
 * Synopsis:  Get the log of the probability of emitting a given sequence.
 * Incept:    BK, Wed Aug  3 12:00:00 2009 [Benasque]
 *
 * Purpose:   Calculate the probability of an <hmm>
 *            generating a given sequence, <dsq>. Return
 *            the log probability in <p>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
fast_hmm_GetLogProb(const ESL_DSQ *dsq, int L, FAST_HMM *hmm, float* p)
{
  int     i;
  __m128  fwd;
  __m128 *combined_probs  =  hmm->combined_probs;
  __m128  tmp;
  float   factor = 1e30;
  int     adj_count = 0;
  __m128  adj_vector;
  float   res;

  // let the sequence start at zero
  dsq++;

  ((float*) &adj_vector)[0] = factor;
  ((float*) &adj_vector)[1] = factor;
  ((float*) &adj_vector)[2] = factor;
  ((float*) &adj_vector)[3] = factor;

  if (!hmm->up_to_date) set_probs(hmm);

  if (L == 0) {
    *p = logf(hmm->t_se);
    return eslOK;
  }

  fwd = hmm->start_probs[dsq[0]];

  i = 1;

  for (; i + 1 < L; i += 2)
    {
      float  tmpf;

      fwd = _mm_mul_ps(fwd, combined_probs[2 * dsq[i    ] + 1]);
      tmp = _mm_shuffle_ps(fwd, fwd, 0xb1);  // 10110001
      fwd = _mm_add_ps(fwd, tmp);

      fwd = _mm_mul_ps(fwd, combined_probs[2 * dsq[i + 1]    ]);
      tmp = _mm_shuffle_ps(fwd, fwd, 0x4e);  // 01001110
      fwd = _mm_add_ps(fwd, tmp);

      _mm_store_ss(&tmpf, fwd);
      if (tmpf < 1) {
        fwd = _mm_mul_ps(fwd, adj_vector);
        adj_count++;
      }
    }

  if (i < L) {
    fwd = _mm_mul_ps(fwd, combined_probs[2 * dsq[i    ] + 1]);
    tmp = _mm_shuffle_ps(fwd, fwd, 0xb1);  // 10110001
    fwd = _mm_add_ps(fwd, tmp);
  }

  fwd = _mm_mul_ps(fwd, hmm->end_probs);

  tmp = _mm_shuffle_ps(fwd, fwd, 0x1b);  // 00011011
  fwd = _mm_add_ps(fwd, tmp);

  _mm_store_ss(&res, fwd);

  *p = logf(res) - adj_count * logf(factor) + L * hmm->emission_adj + (1 + L) * hmm->transition_adj;

  return eslOK;
}



#ifdef FAST_HMM_EXAMPLE
/*::cexcerpt::fast_hmm_example::begin::*/
/* compile: gcc -Wall -I../easel -L../easel -msse2 -o fast_hmm_example -DFAST_HMM_EXAMPLE fast_hmm.c -leasel -lm
 * run:     ./hmm_example <sequence file>
 */
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "fast_hmm.h"
#include "esl_sq.h"
#include "esl_sqio.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "example of the HMM module";



int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go  = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_ALPHABET *abc = esl_alphabet_Create(eslDICE);
  FAST_HMM     *hmm = fast_hmm_Create(abc);
  int           i;
  float        *fair = NULL;
  float        *loaded = NULL;
  ESL_DSQ      *dsq  = NULL;
  int           status;
  float         log_prob;

  ESL_ALLOC(dsq,    sizeof(ESL_DSQ) * 5);
  ESL_ALLOC(fair,   abc->Kp * sizeof(float));
  ESL_ALLOC(loaded, abc->Kp * sizeof(float));

  for (i = 0; i < abc->K;  i++)     fair[i]   = 1.0f /  abc->K;
  for (     ; i < abc->Kp; i++)     fair[i]   = 0.0f;

  for (i = 0; i < abc->K - 1; i++) loaded[i] = 0.5f / (abc->K - 1);
  for (     ; i < abc->Kp;    i++) loaded[i]   = 0.0f;
  loaded[abc->K-1] = 0.5f;

  fast_hmm_SetEmissions(hmm, fair, loaded);
  
  fast_hmm_SetTransitions(hmm,
                          1.00, 0.00, 0.00,
                          0.96, 0.03, 0.01,
                          0.05, 0.95, 0.00);

  printf("\n");
  printf(" Series       NLL      Probability\n");
  printf(" =======  ==========   ==========\n");

  dsq[1] = 0; dsq[2] = 0; dsq[3] = 0; dsq[4] = 0;

  fast_hmm_GetLogProb(dsq, 4, hmm, &log_prob);

  for (i = 1; i < 5; i++) printf(" %c", abc->sym[dsq[i]]);
  printf("  %10.4f   %10.8f\n", log_prob, expf(log_prob));

  dsq[1] = 5; dsq[2] = 5; dsq[3] = 5; dsq[4] = 5;

  fast_hmm_GetLogProb(dsq, 4, hmm, &log_prob);

  for (i = 1; i < 5; i++) printf(" %c", abc->sym[dsq[i]]);
  printf("  %10.4f   %10.8f\n", log_prob, expf(log_prob));

  printf("\n");

  fast_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  free(dsq);
  free(loaded);
  free(fair);

  return 0;

 ERROR:
  return 1;
}
#endif /*FAST_HMM_EXAMPLE*/
