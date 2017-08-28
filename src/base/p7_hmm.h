/* P7_HMM: a core model.
 */
#ifndef p7HMM_INCLUDED
#define p7HMM_INCLUDED

#include "p7_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"

#include "base/general.h"
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif
/* Bit flags used in <hmm->flags>: optional annotation in an HMM
 * 
 * Flags marked with ! may not be changed nor used for other meanings,
 * because they're codes used by HMMER2 (and earlier) that must be
 * preserved for reverse compatibility with old HMMER files.
 * 
 * Why use flags? (So I don't ask this question of myself again:)
 *   1. The way we allocate an HMM, we need to know if we're allocating
 *      M-width annotation fields (RF, CS, CA, MAP) before we read the
 *      annotation from the file.
 *   2. Historically, H2 used flags, so we still need to read H2 flags
 *      for backwards compatibility; so we may as well keep using them.
 */
#define p7H_HASBITS (1<<0)    /* obsolete (was: model has log-odds scores)       !*/
#define p7H_DESC    (1<<1)    /* description exists (legacy; xref SRE:J5/114)    !*/
#define p7H_RF      (1<<2)    /* #RF annotation available                        !*/
#define p7H_CS      (1<<3)    /* #CS annotation available                        !*/
#define p7H_XRAY    (1<<4)    /* obsolete (was: structural data available)       !*/
#define p7H_HASPROB (1<<5)    /* obsolete (was: model in probability form)       !*/
#define p7H_HASDNA  (1<<6)    /* obsolete (was: protein HMM->DNA seq params set) !*/
#define p7H_STATS   (1<<7)    /* model has E-value statistics calibrated         !*/
#define p7H_MAP     (1<<8)    /* alignment map is available                      !*/
#define p7H_ACC     (1<<9)    /* accession is available (legacy; xref SRE:J5/114)!*/
#define p7H_GA      (1<<10)   /* gathering thresholds available                  !*/
#define p7H_TC      (1<<11)   /* trusted cutoffs available                       !*/
#define p7H_NC      (1<<12)   /* noise cutoffs available                         !*/
#define p7H_CA      (1<<13)   /* surface accessibilities available               !*/
#define p7H_COMPO   (1<<14)   /* model-specific residue composition available     */
#define p7H_CHKSUM  (1<<15)   /* model has an alignment checksum                  */
#define p7H_CONS    (1<<16)   /* consensus residue line available                 */
#define p7H_MMASK   (1<<17)   /* #MM annotation available                        !*/
#define p7H_SINGLE  (1<<18)   /* model was from single query w/ Seqmodel()        */

/* Indices of Plan7 main model state transitions, hmm->t[k][] */
enum p7h_transitions_e {
  p7H_MM = 0,
  p7H_MI = 1,
  p7H_MD = 2,
  p7H_IM = 3,
  p7H_II = 4,
  p7H_DM = 5,
  p7H_DD = 6 
};
#define p7H_NTRANSITIONS 7

/* How the hmm->t[k] vector is interpreted as separate probability vectors. */
#define P7H_TMAT(hmm, k) ((hmm)->t[k])
#define P7H_TINS(hmm, k) ((hmm)->t[k]+3)
#define P7H_TDEL(hmm, k) ((hmm)->t[k]+5)
#define p7H_NTMAT 3
#define p7H_NTDEL 2
#define p7H_NTINS 2
#define p7H_NTMAX 3	/* max size of a transition prob vector */

#define p7H_II_SAMPLE_MAX 0.99	/* an artificial cap on sampled tII parameters, to avoid inf-length seqs */

/* Some notes:
 *   0. The model might be either in counts or probability form.
 *   1. t[0] is special: t[0][TMM,TMI,TMD] are the begin->M_1,I_0,D_1 entry probabilities,
 *      t[0][TIM,TII] are the I_0 transitions, and delete state 0 doesn't
 *      exist. Therefore D[0] transitions and mat[0] emissions are unused.
 *      To simplify some normalization code, we adopt a convention that these are set
 *      to valid probability distributions: 1.0 for t[0][TDM] and mat[0][0],
 *      and 0 for the rest.
 *                      MM     MI    MD     IM    II     DM   DD
 *              t[0]  B->M1  B->I0  B->D1 I0->M1 I0->I0  (1)  (0)
 *   2. t[M] is also special: there is no next M/D, only E. So TMD and
 *      TDD are 0; TMM, TIM, and TDM are interpreted as transitions to
 *      E; and thus t[M][TDM] must be 1.0 and [TMD] must be 0:
 *                      MM     MI    MD    IM    II    DM   DD
 *              t[M]  Mm->E  Mm->Im   0  Im->E Im->Im   1    0
 */
typedef struct p7_hmm_s {
  /*::cexcerpt::plan7_core::begin::*/
  int     M;                    /* length of the model (# nodes)                           */
  float **t;                    /* transition prob's. t[(0),1..M][0..p7H_NTRANSITIONS-1]   */
  float **mat;                  /* match emissions.  mat[1..M][0..K-1]                     */ 
  float **ins;                  /* insert emissions. ins[1..M][0..K-1]                     */
  /*::cexcerpt::plan7_core::end::*/

  /* Annotation. Everything but <name> is optional. Flags are set when
   * optional values are set. All the char *'s are proper nul-terminated
   * strings, not just arrays. (hmm->map is an int array).
   */
  char    *name;                 /* name of the model                     (mandatory)      */ /* String, \0-terminated   */
  char    *acc;	                 /* accession number of model (Pfam)      (optional: NULL) */ /* String, \0-terminated   */
  char    *desc;                 /* brief (1-line) description of model   (optional: NULL) */ /* String, \0-terminated   */
  char    *rf;                   /* reference line from alignment 1..M    (p7H_RF)         */ /* String; 0=' ', M+1='\0' */
  char    *mm;                   /* model mask line from alignment 1..M   (p7H_MM)         */ /* String; 0=' ', M+1='\0' */
  char    *consensus;	         /* consensus residue line        1..M    (p7H_CONS)       */ /* String; 0=' ', M+1='\0' */
  char    *cs;                   /* consensus structure line      1..M    (p7H_CS)         */ /* String; 0=' ', M+1='\0' */
  char    *ca;	                 /* consensus accessibility line  1..M    (p7H_CA)         */ /* String; 0=' ', M+1='\0' */

  char    *comlog;               /* command line(s) that built model      (optional: NULL) */ /* String, \0-terminated   */
  int      nseq;	         /* number of training sequences          (optional: -1)   */
  float    eff_nseq;             /* effective number of seqs (<= nseq)    (optional: -1)   */
  int	   max_length;           /* upper bound length, all but 1e-7 prob (optional: -1)   */
  char    *ctime;	         /* creation date                         (optional: NULL) */
  int     *map;	                 /* map of alignment cols onto model 1..M (p7H_MAP)        */ /* Array; map[0]=0 */
  uint32_t checksum;             /* checksum of training sequences        (p7H_CHKSUM)     */
  float    evparam[p7_NEVPARAM]; /* E-value params                        (p7H_STATS)      */
  float    cutoff[p7_NCUTOFFS];  /* Pfam score cutoffs                    (p7H_{GA,TC,NC}) */
  float    compo[p7_MAXABET];    /* model bg residue comp                 (p7H_COMPO)      */

  off_t    offset;               /* HMM record offset on disk                              */
  const ESL_ALPHABET *abc;       /* ptr to alphabet info (hmm->abc->K is alphabet size)    */
  int      flags;                /* status flags                                           */
} P7_HMM;


/*      1. The P7_HMM object: allocation, initialization, destruction. */
extern P7_HMM *p7_hmm_Create(int M, const ESL_ALPHABET *abc);
extern P7_HMM *p7_hmm_CreateShell(void);
extern int     p7_hmm_CreateBody(P7_HMM *hmm, int M, const ESL_ALPHABET *abc);
extern void    p7_hmm_Destroy(P7_HMM *hmm);
extern int     p7_hmm_CopyParameters(const P7_HMM *src, P7_HMM *dest);
extern P7_HMM *p7_hmm_Clone(const P7_HMM *hmm);
extern int     p7_hmm_Zero(P7_HMM *hmm);
/*      2. Convenience routines for setting fields in an HMM. */
extern int     p7_hmm_SetName       (P7_HMM *hmm, char *name);
extern int     p7_hmm_SetAccession  (P7_HMM *hmm, char *acc);
extern int     p7_hmm_SetDescription(P7_HMM *hmm, char *desc);
extern int     p7_hmm_AppendComlog  (P7_HMM *hmm, int argc, char **argv);
extern int     p7_hmm_SetCtime      (P7_HMM *hmm);
extern int     p7_hmm_SetComposition(P7_HMM *hmm);
extern int     p7_hmm_SetConsensus  (P7_HMM *hmm, ESL_SQ *sq);
/*      3. Renormalization and rescaling counts in core HMMs. */
extern int     p7_hmm_Scale      (P7_HMM *hmm, double scale);
extern int     p7_hmm_Renormalize(P7_HMM *hmm);
/*      4. Debugging and development code. */
extern int     p7_hmm_Dump(FILE *fp, P7_HMM *hmm);
extern int     p7_hmm_Compare(P7_HMM *h1, P7_HMM *h2, float tol);
extern int     p7_hmm_Validate(P7_HMM *hmm, char *errbuf, float tol);
/*      5. Other routines in the API */
extern int     p7_hmm_CalculateOccupancy(const P7_HMM *hmm, float *mocc, float *iocc);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif /*p7HMM_INCLUDED*/


