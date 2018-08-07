/* H4_PROFILE: dual-mode local/glocal profile HMM
 * 
 * See h4_profile.md for notes on the H4_PROFILE structure.
 */
#ifndef h4PROFILE_INCLUDED
#define h4PROFILE_INCLUDED
#include "h4_config.h"

#include "esl_alphabet.h"

/* Don't change this order. Code assumes that 
 * 0..2, 3..5, 6..8 are three separate prob distributions for M,I,D,
 * and sometimes hardcodes the +3, +6 offsets.
 */
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


/* H4_PROFILE
 * A HMMER4 profile HMM, dual-mode local/glocal 
 * 
 * Brief summary of fixed edge conditions; also see h4_profile.md, and h4_profile_SetConventions()
 *  t[0] = [ t(G->M1), 0, t(G->D1); 1,0,0; 1,0,0 ]. Only TMM,TMD values are free params.
 *  t[M] = [ 1,0,0; 1,0,0; 1,0,0 ]. None free; these are t(M->E), Im doesn't exist, t(D->E).
 *  e[0] = [ 1,0..0].
 * So, if you're setting free params: 2 values at t[0], t[1..M-1], e[1..M] 
 */
typedef struct {
  int     M;                // model length in nodes (consensus positions)
  float **t;                // transitions.     [0..M][0..8].   [0]: only MM,MD are free. [M]: none free. 
  float **e;                // match emissions. [0..M][0..K-1]  [0]: not free.

  const ESL_ALPHABET *abc;  // reference ptr to alphabet. (only a copy; don't free)
} H4_PROFILE;


extern H4_PROFILE *h4_profile_Create(const ESL_ALPHABET *abc, int M);
extern H4_PROFILE *h4_profile_CreateShell(void);
extern int         h4_profile_CreateBody(H4_PROFILE *hmm, const ESL_ALPHABET *abc, int M);
extern void        h4_profile_Destroy(H4_PROFILE *hmm);

extern int         h4_profile_SetConventions(H4_PROFILE *hmm);
extern int         h4_profile_Renormalize   (H4_PROFILE *hmm);


extern int         h4_profile_Dump(FILE *fp, H4_PROFILE *hmm);
extern int         h4_profile_Validate(const H4_PROFILE *hmm, char *errbuf);
extern int         h4_profile_Compare(const H4_PROFILE *h1, const H4_PROFILE *h2);

#endif /* h4PROFILE_INCLUDED */
