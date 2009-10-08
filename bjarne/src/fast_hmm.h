/* A fast sse implementation of two state hidden Markov models (discrete; of alphabetic strings)
 * 
 * SRE, Fri Jul 18 08:54:41 2008 [Janelia] 
 * SVN $Id$
 */
#ifndef FAST_HMM_INCLUDED
#define FAST_HMM_INCLUDED

#include <xmmintrin.h>

#include "esl_alphabet.h"

typedef struct {
  __m128  end_probs;            /* A vector of probabilities from general states to end */
  __m128 *emission_probs;       /* A set of vectors of emission probabilities */
  __m128 *combined_probs;       /* A set of vectors of combined transition and emission probabilities */
  __m128 *start_probs;          /* A set of vectors of combined start and emission probabilities */
  float   t_s0;                 /* Individual transition probabilities */
  float   t_s1;
  float   t_se;
  float   t_00;
  float   t_01;
  float   t_0e;
  float   t_10;
  float   t_11;
  float   t_1e;
  void   *prob_ptr;             /* An unaligned memory pointer used for vectors (after aligning) */
  void   *ptr;                  /* An unaligned memory pointer used for this struct (which is then aligned) */
  int     alpha_size;           /* Alphabet size */
  int     up_to_date;           /* Indicates whether the probability vectors are up to date */
  float   emission_adj;         /* An adjustment factor of the emissions to postpone underflow */
  float   transition_adj;       /* An adjustment factor of the transitions to postpone underflow */
} FAST_HMM;


extern FAST_HMM *fast_hmm_Create(const ESL_ALPHABET *abc);
extern void      fast_hmm_Destroy(FAST_HMM *hmm);
extern int       fast_hmm_SetEmissions(FAST_HMM *hmm, float *emissions_0, float *emissions_1);
extern int       fast_hmm_SetTransitions(FAST_HMM *hmm, float t_s0, float t_s1, float t_se, float t_00, float t_01, float t_0e, float t_10, float t_11, float t_1e);
extern int       fast_hmm_GetLogProb(const ESL_DSQ *dsq, int L, FAST_HMM *hmm, float *p);


#endif /*!FAST_HMM_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
