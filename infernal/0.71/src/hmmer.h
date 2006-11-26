/* hmmer.h
 *    1. Declaration of structures/objects.
 *        - P7_PROFILE
 *        - P7_EVINFO
 *        - P7_HMM
 *        - P7_TRACE
 *        - P7_GMX
 *        - [P7_OPROFILE]
 *        - [P7_OMX]
 *    2. Definition of constants
 *    3. Declaration of functions.
 * 
 * SRE, Fri Apr 14 09:08:56 2006 [St. Louis]
 * SVN $Id$
 */
#ifndef HMMER_INCLUDED
#define HMMER_INCLUDED

#include "p7_config.h"

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_msa.h>


/*****************************************************************
 * 1. Declarations of structures.
 *****************************************************************/


/*----------------------------------------------------------------
 * P7_PROFILE 
 * Generic score profile.
 */
typedef struct {
  int    mode;         	/* configured algorithm mode (e.g. p7_LOCAL)   */ 
  ESL_ALPHABET *abc;	/* copy of pointer to appropriate alphabet     */
  int    M;		/* number of nodes in the model                */
  int  **tsc;           /* transition scores     [0.6][1.M-1]          */
  int  **msc;           /* match emission scores [0.Kp-1][1.M]         */
  int  **isc;           /* ins emission scores   [0.Kp-1][1.M-1]       */
  int    xsc[4][2];     /* N,E,C,J transitions   [][LOOP,MOVE]         */
  int   *bsc;           /* begin transitions     [1.M]                 */
  int   *esc;		/* end transitions       [1.M]                 */
} P7_PROFILE;




/* ----------------------------------------------------------------
 * P7_EVINFO
 * Statistical parameters for E-value calculations.
 */
typedef struct {
  int    mode;		/* what type of distribution this is          */
  double mu;		/* location param (for gumbel or exponential) */
  double lambda;	/* scale param (gumbel or exponential)        */
  double tailmass;	/* 1.0, or fraction of tail that's fit        */

  int N;		/* this info was calibrated on <N> seqs...    */
  int L;		/*   ...of length <L>.                        */
} P7_EVINFO;

/* Possible distributions (for P7_EVINFO mode): */
#define p7_EXPTAIL 0
#define p7_GUMBEL  1


/*----------------------------------------------------------------
 * P7_HMM
 * The Plan7 HMM structure. 
 * 
 * Includes the core probability model, annotation, a copy of the
 * alphabet pointer; may also contain a profile, an optimized profile,
 * and statistical parameters.
 */
typedef struct {
  /* Alphabet information (includes hmm->abc->K, alphabet size)
   */
  ESL_ALPHABET *abc;

  /* The core model in probability form.
   * p7_HASPROBS flag is raised when these probs are all valid.
   */
  /*::cexcerpt::plan7_core::begin::*/
  int     M;            /* length of the model (# nodes)          */
  float **t;            /* transition prob's. t[(0),1..M-1][0..6] */
  float **mat;          /* match emissions.  mat[1..M][0..K-1]    */ 
  float **ins;          /* insert emissions. ins[1..M-1][0..K-1]  */
  /*::cexcerpt::plan7_core::end::*/

  /* The null model probabilities.
   */
  float  *null;		/* "random sequence" emission prob's     */
  float   p1;           /* null model loop probability           */

  /* The unique states of Plan 7 in probability form.
   * These are the algorithm-dependent, data-independent probabilities.
   */
  float   xt[4][2];     /* N,E,C,J extra states: 2 transitions        +*/
  float  *begin;        /* 1..M B->M state transitions                +*/
  float  *end;          /* 1..M M->E state transitions (!= a dist!)   +*/

  /* The model's log-odds score form: xref modelconfig.c
   * PLAN7_HASBITS flag is up when gm is non-NULL.
   */
  P7_PROFILE  *gm;	/* common profile                                */
  struct p7_oprof_s *om;/* P7_OPROFILE; optional optimized profile       */
  float  lscore;	/* sc for unaligned res when LCORRECT flag is up */

  /* Annotation on the model. A name is mandatory.
   * Other fields are optional; whether they are present is
   * flagged in the stateflags bit array, as indicated in comments below.
   */
  char  *name;        /* name of the model                             */
  char  *acc;	      /* accession number of model (Pfam)      p7_ACC  */
  char  *desc;        /* brief description of model            p7_DESC */ 
  char  *rf;          /* reference line from alignment 0..M    p7_RF   */
  char  *cs;          /* consensus structure line      0..M    p7_CS   */ 
  char  *ca;	      /* consensus accessibility line  0..M    p7_CA   */
  char  *comlog;      /* command line(s) that built model              */
  int    nseq;	      /* number of training sequences                  */
  char  *ctime;	      /* creation date                                 */
  int   *map;	      /* map of alignment cols onto model 1..M p7_MAP  */
  int    checksum;    /* checksum of training sequences                */

  /* Pfam-specific score cutoffs.
   * 
   * ga1, ga2 are valid if PLAN7_GA is set in flags.
   * tc1, tc2 are valid if PLAN7_TC is set in flags.
   * nc1, nc2 are valid if PLAN7_NC is set in flags.
   */
  float  ga1, ga2;	/* per-seq/per-domain gathering thresholds (bits) +*/
  float  tc1, tc2;	/* per-seq/per-domain trusted cutoff (bits)       +*/
  float  nc1, nc2;	/* per-seq/per-domain noise cutoff (bits)         +*/

  /* E-value statistical parameters.
   * If any are non-NULL, a flag is raised in <flags>: 
   * PLAN7_STATS_LV for lvstats, etc.
   */
  P7_EVINFO *lvstats;	/* local Viterbi  */
  P7_EVINFO *lfstats;	/* local Forward  */
  P7_EVINFO *gvstats;	/* glocal Viterbi */
  P7_EVINFO *gfstats;	/* glocal Forward */

  int flags;            /* bit flags indicating state of HMM, valid data +*/
} P7_HMM;


/*----------------------------------------------------------------
 * P7_TRACE
 * Traceback structure for alignments of model to sequence.
 * Element 0 always p7_STS. Element N-1 always p7_STT.
 */
typedef struct {
  int   N;		/* length of traceback                       */
  int   nalloc;		/* allocated length of traceback             */
  char *st;		/* state type code                   [0..N-1]*/
  int  *k;		/* node index; 1..M if M,D,I; else 0 [0..N-1]*/
  int  *i;		/* position in dsq, 1..L; 0 if none  [0..N-1]*/
} P7_TRACE;


/*----------------------------------------------------------------
 * P7_GMX
 * The generic DP matrix.
 */
typedef struct {
  int **xmx;			/* special scores [0.1..N][BECJN]     */
  int **mmx;			/* match scores [0.1..N][0.1..M]      */
  int **imx;			/* insert scores [0.1..N][0.1..M-1.M] */
  int **dmx;			/* delete scores [0.1..N][0.1..M-1.M] */
  
  int N;			/* alloc'ed for seq of length N; N+1 rows */
  int M;			/* alloc'ed for HMM of length M; M+1 cols */

  /* If either pad is 0, we're not growable in that direction.       */
  int padN;			/* extra pad in sequence length/rows */
  int padM;			/* extra pad in HMM length/columns   */
} P7_GMX;


#ifdef    p7_IMPL_REFERENCE
#include "p7_dp_reference.h"	/* the reference implementation.        */
#elif     p7_IMPL_FAST
#include "p7_dp_fast.h"		/* our optimized implementation.        */
#elif     p7_IMPL_ALTIVEC
#include "p7_dp_altivec.h"	/* Erik Lindahl's Altivec for PowerPC   */
#elif     p7_IMPL_BUHLER
#include "p7_dp_buhler.h"	/* Jeremy Buhler's optimization         */
#elif     p7_IMPL_SSE
#include "p7_dp_sse.h"		/* Apple's SSE implementation.          */
#endif


/* Flag codes for plan7->flags.
 * Don't ever change these; we want to be able to read old binary save files
 * that used these values.
 */
#define p7_HASBITS (1<<0)    /* model has log-odds scores            */
#define p7_DESC    (1<<1)    /* description exists                   */
#define p7_RF      (1<<2)    /* #RF annotation available             */
#define p7_CS      (1<<3)    /* #CS annotation available             */
#define p7_XRAY    (1<<4)    /* structural data available            */
#define p7_HASPROB (1<<5)    /* model has probabilities              */
#define p7_HASDNA  (1<<6)	/* protein HMM->DNA seq params set      */
#define p7_STATS   (1<<7)    /* obsolete (2.3 and earlier)           */
#define p7_MAP     (1<<8)	/* alignment map is available           */
#define p7_ACC     (1<<9)	/* accession number is available        */
#define p7_GA      (1<<10)	/* gathering thresholds available       */
#define p7_TC      (1<<11)	/* trusted cutoffs available            */
#define p7_NC      (1<<12)	/* noise cutoffs available              */
#define p7_CA      (1<<13)   /* surface accessibility avail.         */
#define p7_BIMPOSED (1<<14)  /* all entries are B->M_k (not D)       */
#define p7_EIMPOSED (1<<15)  /* all ends are M_k->E (not D)          */
#define p7_LCORRECT (1<<16)  /* require L-dependent score correction */
#define p7_STATS_LV (1<<17)	/* local Viterbi E-val params available */
#define p7_STATS_LF (1<<18)	/* local Forward E-val params available */
#define p7_STATS_GV (1<<19)	/* have glocal Viterbi E-val params     */
#define p7_STATS_GF (1<<20)	/* have glocal Forward E-val params     */


/* Indices for special state types, I: used for dynamic programming xmx[][]
 * mnemonic: eXtra Matrix for B state = XMB
 */
#define p7_XMB 0
#define p7_XME 1
#define p7_XMC 2
#define p7_XMJ 3
#define p7_XMN 4

/* Indices for special state types, II: used for hmm->xt[] indexing
 * mnemonic: eXtra Transition for N state = XTN
 */
#define p7_XTN  0
#define p7_XTE  1
#define p7_XTC  2
#define p7_XTJ  3

/* Indices for Plan7 main model state transitions.
 * Used for indexing hmm->t[k][]
 * mnemonic: Transition from Match to Match = TMM
 */
#define p7_TMM  0
#define p7_TMI  1
#define p7_TMD  2
#define p7_TIM  3
#define p7_TII  4
#define p7_TDM  5
#define p7_TDD  6 

/* Indices for extra state transitions
 * Used for indexing hmm->xt[][].
 */
#define p7_MOVE 0          /* trNB, trEC, trCT, trJB */
#define p7_LOOP 1          /* trNN, trEJ, trCC, trJJ */

/* Plan 7 model state types (esp. used in P7_TRACE structure)
 */
#define p7_BOGUS 0
#define p7_STM   1
#define p7_STD   2
#define p7_STI   3
#define p7_STS   4
#define p7_STN   5
#define p7_STB   6
#define p7_STE   7
#define p7_STC   8
#define p7_STT   9
#define p7_STJ   10     

/* Search modes.
 */
#define p7_LOCAL     1		/* multi-hit local:  "fs" mode   */
#define p7_GLOCAL    2		/* multi-hit glocal: "ls" mode   */
#define p7_UNILOCAL  3		/* one-hit local: "sw" mode      */
#define p7_UNIGLOCAL 4		/* one-hit glocal: "s" mode      */

/* Error codes thrown by HMMER.
 */
#define p7ERR_CONSENSUS   1001    /* modelmakers.c:matassign2hmm(),     M=0, no consensus columns */
#define p7ERR_RF          1002    /* modelmakers.c:p7_Handmodelmaker(), no #=RF annotation        */





/* hmm.c:  the model structure.
 */
extern P7_HMM *p7_hmm_Create(int M, ESL_ALPHABET *abc);
extern P7_HMM *p7_hmm_CreateShell(void);
extern int     p7_hmm_CreateBody(P7_HMM *hmm, int M, ESL_ALPHABET *abc);
extern void    p7_hmm_Destroy(P7_HMM *hmm);
extern int     p7_hmm_ZeroCounts(P7_HMM *hmm);
extern int     p7_hmm_SetName(P7_HMM *hmm, char *name);
extern int     p7_hmm_SetAccession(P7_HMM *hmm, char *acc);
extern int     p7_hmm_SetDescription(P7_HMM *hmm, char *desc);
extern int     p7_hmm_Comlog(P7_HMM *hmm, int argc, char **argv);
extern int     p7_hmm_SetCtime(P7_HMM *hmm);
extern int     p7_hmm_SetNull(P7_HMM *hmm, float *null, int K);
extern int     p7_hmm_Rescale(P7_HMM *hmm, float scale);
extern int     p7_hmm_Renormalize(P7_HMM *hmm);
extern void    p7_hmm_Dump(FILE *fp, P7_HMM *hmm);
extern char   *p7_hmm_Statetype(char st);

/* modelmakers.c:  constructing HMMs from alignments
 */
extern int p7_Handmodelmaker(ESL_MSA *msa, ESL_ALPHABET *abc, char **dsq, 
			     char *isfrag,
			     P7_HMM **ret_hmm, P7_TRACE ***ret_tr);
extern int p7_Fastmodelmaker(ESL_MSA *msa, ESL_ALPHABET *abc, char **dsq,
                             char *isfrag, float symfrac,
			     P7_HMM **ret_hmm, P7_TRACE ***ret_tr);

/* profile.c: the P7_PROFILE structure.
 */
extern P7_PROFILE *p7_profile_Create(int M);
extern void        p7_profile_Destroy(P7_PROFILE *gm);
extern int         p7_profile_GetTScore(P7_PROFILE *gm, 
					char st1, int k1, char st2, int k2,
					int *ret_tsc);

/* trace.c: the P7_TRACE structure.
 */
extern int  p7_trace_Create(int N, P7_TRACE **ret_tr);
extern int  p7_trace_Expand(P7_TRACE *tr);
extern int  p7_trace_ExpandTo(P7_TRACE *tr, int N);
extern void p7_trace_Destroy(P7_TRACE *tr);
extern int  p7_trace_Dump(FILE *fp, P7_TRACE *tr, P7_PROFILE *gm, char *dsq);
extern int  p7_trace_Append(P7_TRACE *tr, char st, int k, int i);
extern int  p7_trace_Reverse(P7_TRACE *tr);
extern int  p7_trace_Count(P7_HMM *hmm, char *dsq, float wt, P7_TRACE *tr, int mode);
extern int  p7_trace_Score(P7_PROFILE *gm, char *dsq, P7_TRACE *tr, int *ret_sc);
extern int  p7_trace_DomainCount(P7_TRACE *tr);
extern int  p7_trace_Decompose(P7_TRACE *otr, P7_TRACE ***ret_tr, int *ret_ntr);
extern int  p7_trace_GetDomainCoords(P7_TRACE *tr, int which,
				     int *ret_i1, int *ret_i2,
				     int *ret_k1, int *ret_k2,
				     int *ret_avlen);





#endif /*HMMER_INCLUDED*/
