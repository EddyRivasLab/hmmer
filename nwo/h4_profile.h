/* H4_PROFILE: dual-mode local/glocal profile HMM
 */
#ifndef h4PROFILE_INCLUDED
#define h4PROFILE_INCLUDED
#include "h4_config.h"

#include "esl_alphabet.h"

enum h4_transitions_e {
  h4_TMM = 0,
  h4_TMI = 1,
  h4_TMD = 2,
  h4_TIM = 3,
  h4_TII = 4,
  h4_TID = 5, 
  h4_TDM = 6, 
  h4_TDI = 7,
  h4_TDD = 8
};
#define h4_NTRANSITIONS 9


typedef struct {
  int     M;                // model length in nodes (consensus positions)
  float **t;                // transitions.     [0..M][0..8]
  float **e;                // match emissions. [0..M][0..K-1]

  const ESL_ALPHABET *abc;  // reference ptr to alphabet. (only a copy; don't free)
} H4_PROFILE;


extern H4_PROFILE *h4_profile_Create(ESL_ALPHABET *abc, int M);
extern H4_PROFILE *h4_profile_CreateShell(void);
extern int         h4_profile_CreateBody(H4_PROFILE *hmm, const ESL_ALPHABET *abc, int M);
extern void        h4_profile_Destroy(H4_PROFILE *hmm);

#endif /* h4PROFILE_INCLUDED */
