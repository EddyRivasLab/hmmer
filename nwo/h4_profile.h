/* H4_PROFILE: dual-mode local/glocal profile HMM
 * 
 * See h4_profile.md for notes on the H4_PROFILE structure.
 */
#ifndef h4PROFILE_INCLUDED
#define h4PROFILE_INCLUDED

#include "h4_config.h"

#include "esl_alphabet.h"

/* H4_PROFILE
 * A HMMER4 profile HMM, dual-mode local/glocal 
 * 
 * Brief summary of fixed edge conditions; also see h4_profile.md, and h4_profile_SetConventions()
 *  t[0] = [ t(G->M1), 0, t(G->D1); 1,0,0; 1,0,0 ]. Only TMM,TMD values are free params.
 *  t[M] = [ 1,0,0; 1,0,0; 1,0,0 ]. None free; these are t(M->E), Im doesn't exist, t(D->E).
 *  e[0] = [ 1,0..0].
 * So, if you're setting free params: 2 values at t[0], t[1..M-1], e[1..M] 
 *
 * t,e,tsc,rsc are standard Easel 2D arrays, and can be manipulated with esl_mat_F*() functions.
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


/* Constants defining fixed sizes of parameter arrays in a profile.  
 *
 * It's not that you can change these easily; rather, these are so you
 * don't see bare numbers like '9' in the code and wonder where it
 * came from and what it's supposed to mean.
 */
#define h4_NT      9     // number of transition probabilities in hmm->t
#define h4_NTSC    12    // number of transition scores in hmm->tsc


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

/* 12 (h4_TSC) transition scores, ordered for efficiency in dynamic
 * programming algorithms. Besides a score for each transition
 * probability, three scores are precomputed for L->Mk and
 * "wing-retracted" G->...->Mk entries, and for Mk->...->E exits.
 */
#define h4_MM   0
#define h4_IM   1
#define h4_DM   2
#define h4_LM   3
#define h4_GM   4
#define h4_MD   5
#define h4_ID   6
#define h4_DD   7
#define h4_MI   8
#define h4_II   9
#define h4_DI   10
#define h4_DGE  11


/* Flags that can be raised in <hmm->flags>.
 * Which optional annotations are available, for example.
 * 
 * Flags marked with ! may not be changed nor used for other meanings,
 * because they're codes used by HMMER2 (and earlier) that must be
 * preserved for reverse compatibility with old HMMER files.
 * 
 * Why use flags? (So I don't ask this question of myself again:)
 *   1. The way we allocate an HMM, we need to know if we're allocating
 *      M-width annotation fields (RF, CS, CA, MAP) before we read the
 *      annotation from a binary HMM file.
 *   2. Historically, H2 used flags, so we still need to read H2 flags
 *      from H2 files.
 */
#define h4_HASBITS (1<<0)    /* obsolete (was: model has log-odds scores)       !*/
#define h4_DESC    (1<<1)    /* description exists (legacy; xref SRE:J5/114)    !*/
#define h4_RF      (1<<2)    /* #RF annotation available                        !*/
#define h4_CS      (1<<3)    /* #CS annotation available                        !*/
#define h4_XRAY    (1<<4)    /* obsolete (was: structural data available)       !*/
#define h4_HASPROB (1<<5)    /* obsolete (was: model in probability form)       !*/
#define h4_HASDNA  (1<<6)    /* obsolete (was: protein HMM->DNA seq params set) !*/
#define h4_STATS   (1<<7)    /* model has E-value statistics calibrated         !*/
#define h4_MAP     (1<<8)    /* alignment map is available                      !*/
#define h4_ACC     (1<<9)    /* accession is available (legacy; xref SRE:J5/114)!*/
#define h4_GA      (1<<10)   /* gathering thresholds available                  !*/
#define h4_TC      (1<<11)   /* trusted cutoffs available                       !*/
#define h4_NC      (1<<12)   /* noise cutoffs available                         !*/
#define h4_CA      (1<<13)   /* surface accessibilities available               !*/
#define h4_COMPO   (1<<14)   /* model-specific residue composition available     */
#define h4_CHKSUM  (1<<15)   /* model has an alignment checksum                  */
#define h4_CONS    (1<<16)   /* consensus residue line available                 */
#define h4_MMASK   (1<<17)   /* #MM annotation available                        !*/
#define h4_SINGLE  (1<<18)   /* model was from single query w/ Seqmodel()        */



extern H4_PROFILE *h4_profile_Create(const ESL_ALPHABET *abc, int M);
extern H4_PROFILE *h4_profile_CreateShell(void);
extern int         h4_profile_CreateBody(H4_PROFILE *hmm, const ESL_ALPHABET *abc, int M);
extern H4_PROFILE *h4_profile_Clone  (const H4_PROFILE *hmm);
extern int         h4_profile_Copy   (const H4_PROFILE *src, H4_PROFILE *dst);
extern size_t      h4_profile_Sizeof (const H4_PROFILE *hmm);
extern void        h4_profile_Destroy(H4_PROFILE *hmm);

extern int         h4_profile_SetConventions(H4_PROFILE *hmm);
extern int         h4_profile_Renormalize   (H4_PROFILE *hmm);
extern int         h4_profile_CalculateOccupancy(const H4_PROFILE *hmm, float *mocc, float *iocc);

extern int         h4_profile_Config(H4_PROFILE *hmm);

extern int         h4_profile_Dump(FILE *fp, H4_PROFILE *hmm);
extern int         h4_profile_Validate(const H4_PROFILE *hmm, char *errbuf);
extern int         h4_profile_Compare(const H4_PROFILE *h1, const H4_PROFILE *h2);

#endif /* h4PROFILE_INCLUDED */
