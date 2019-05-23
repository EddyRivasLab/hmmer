/* H4_PROFILE: dual-mode local/glocal profile HMM
 * 
 * See h4_profile.md for notes on the H4_PROFILE structure.
 */
#ifndef h4PROFILE_INCLUDED
#define h4PROFILE_INCLUDED

#include "h4_config.h"

#include "esl_alphabet.h"

/* H4_PROFILE
 * A HMMER4 profile HMM, dual-mode local/glocal, with both probabilities and log-odds scores.
 * 
 * Brief summary of fixed edge conditions; also see h4_profile.md, and h4_profile_SetConventions()
 *  t[0] = [ tGM1, 0, tGD1;  1,0,0;  1,0,0 ].   Only TMM,TMD values are free params.
 *  t[M] = [    1, 0,    0;  1,0,0;  1,0,0 ].   None are free. These are t(M->E), Im doesn't exist, t(D->E).
 *  e[0] = [ 1,0..0].
 * So, if you're setting free params: you set t[0][0,2], t[1..M-1][], and e[1..M][].
 *
 * t,e,tsc,rsc are std Easel 2D arrays that can be manipulated with esl_mat_F*() functions.
 */
typedef struct {
  int     M;                // model length in nodes (consensus positions)
  float **t;                // transition probs.      [0..M][0..h4_NT-1].  [0]: only MM,MD are free. [M]: none free. 
  float **e;                // match emissions probs. [0..M][0..K-1]       [0]: not free.

  float  *f;                // nonhomologous emission probs, for NCJI states. [0..K-1]

  float **tsc;              // transition scores.     [0.1..M][0..h4_NTSC-1].  Sometimes we step through tsc[0] as one big (M+1) x 12 vector.
  float **rsc;              // match emission scores. [0..Kp-1][0.1..M]

  uint32_t flags;           // various boolean flags about this profile

  const ESL_ALPHABET *abc;  // reference ptr to alphabet. (only a copy; don't free)
} H4_PROFILE;


/* H4_PROFILE_CT
 * Profile HMM, counts-collection form.
 * 
 * Profile structure used to collect counts, where we use double
 * instead of float, and we only need the probability model part of
 * the H4_PROFILE structure. Deep alignments (and simulations) can
 * exceed dynamic range of float - especially on II transitions.
 */
typedef struct {
  int    M;
  double **t;
  double **e;

  const ESL_ALPHABET *abc;
} H4_PROFILE_CT;



/* Constants defining fixed sizes of parameter arrays in a profile.  
 *
 * It's not that you can change these easily; rather, these are so you
 * don't see bare numbers like '9' in the code and wonder where it
 * came from and what it's supposed to mean.
 */
#define h4_NT      9     // number of transition probabilities in hmm->t
#define h4_NTSC    13    // number of transition scores in hmm->tsc


/* 9 (h4_NT) transition probabilities, ordered for convenience in
 * normalizing them.
 *
 * Don't change this order. Code assumes that 0..2, 3..5, 6..8 are
 * three separate prob distributions for M,I,D, and sometimes
 * hardcodes the +3, +6 offsets.
 */
#define h4_TMM  0
#define h4_TMI  1
#define h4_TMD  2
#define h4_TIM  3
#define h4_TII  4
#define h4_TID  5 
#define h4_TDM  6 
#define h4_TDI  7
#define h4_TDD  8

/* 13 (h4_NTSC) transition scores, ordered for efficiency in (forward)
 * dynamic programming algorithms. Besides a score for each transition
 * probability, four additional scores are precomputed:
 *   -  L->Mk local entry
 *   -  left wing retracted G->D1..Dk-1->Mk glocal entry
 *   -               ...and G->D1..Dk->Ik glocal entry
 *   -  right wing retracted Dk->...E glocal exit
 */
#define h4_MM   0
#define h4_IM   1
#define h4_DM   2
#define h4_LM   3
#define h4_GM   4
#define h4_MI   5
#define h4_II   6
#define h4_DI   7
#define h4_GI   8
#define h4_MD   9
#define h4_ID   10
#define h4_DD   11
#define h4_DGE  12


/* flags, in flux */
#define h4_HASPROBS (1<<0)
#define h4_HASBITS  (1<<1)
#define h4_SINGLE   (1<<2)


extern H4_PROFILE *h4_profile_Create(const ESL_ALPHABET *abc, int M);
extern H4_PROFILE *h4_profile_CreateShell(void);
extern int         h4_profile_CreateBody(H4_PROFILE *hmm, const ESL_ALPHABET *abc, int M);
extern H4_PROFILE *h4_profile_Clone  (const H4_PROFILE *hmm);
extern int         h4_profile_Copy   (const H4_PROFILE *src, H4_PROFILE *dst);
extern size_t      h4_profile_Sizeof (const H4_PROFILE *hmm);
extern void        h4_profile_Destroy(H4_PROFILE *hmm);

extern H4_PROFILE_CT *h4_profile_ct_Create(const ESL_ALPHABET *abc, int M);
extern void           h4_profile_ct_Destroy(H4_PROFILE_CT *ctm);

extern int         h4_profile_SetConventions(H4_PROFILE *hmm);
extern int         h4_profile_Renormalize   (H4_PROFILE *hmm);
extern int         h4_profile_Occupancy(const H4_PROFILE *hmm, float *mocc, float *iocc, float *opt_mtot, float *opt_itot);

extern int         h4_profile_Config(H4_PROFILE *hmm);

extern int         h4_profile_Dump   (FILE *fp, H4_PROFILE *hmm);
extern int         h4_profile_ct_Dump(FILE *fp, H4_PROFILE_CT *ctm);
extern int         h4_profile_Validate(const H4_PROFILE *hmm, char *errbuf);
extern int         h4_profile_Compare(const H4_PROFILE *h1, const H4_PROFILE *h2);
extern int         h4_profile_MutePathScore(const H4_PROFILE *hmm, float *ret_sc);

#endif /* h4PROFILE_INCLUDED */
