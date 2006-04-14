/* plan7.h
 * The Plan7 HMM structure used by HMMER.
 * 
 * SRE, Sat Apr 30 13:05:35 2005
 * SVN $Id$
 */
#ifndef PLAN7_INCLUDED
#define PLAN7_INCLUDED

#include "p7config.h"
#include <esl_alphabet.h>	/* ESL_ALPHABET       */
#include "p7_profile.h"		/* P7_PROFILE, P7_GMX */
#include "p7_evalues.h"		/* P7_EVINFO          */

/* P7_OPROFILE, P7_OMX are provided by one and only one optimized impl: */
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





/* Structure: P7_HMM 
 *    
 * 2. The "configured" model is the scoring form.
 *    It adds the S,N,B; J; E,C,T special states of the Plan7
 *    "algorithm-dependent" probability architecture.
 *    
 *    The S and T states are not explicitly represented anywhere.
 *    S->N is implicitly 1.0. 
 *    
 *    N,E,C,J transitions are set in xt[]. N,C,J emissions are implicitly
 *    equal to the null model emission probabilities, and these states emit
 *    on transition. 
 *    
 *    Entry into the model is controlled by begin[], which is a normalized
 *    prob distribution \sum_{k=1}^{k=M} t(B->M_k) = 1.0. Exit from the model
 *    is controlled by end[], which is not a probability distribution itself;
 *    rather, \sum_{x=MDI} t(M_k->x) + t(M_k->E) = 1.0. 
 *    
 *    A configured model has no D_1 or D_M state. A process called
 *    "wing retraction" is applied by the configuration functions, to
 *    eliminate the mute all-delete B->D_1...D_M->E path. Wing
 *    retraction requires altering the transition probabilities in the
 *    nodes. See modelconfig.c.
 *    
 *    Only two combinations of wing retraction and internal entry/exit
 *    are allowed. In sw/fs mode models, algorithm-dependent entry and
 *    exit probabilities partially trump wing retraction, so begin[]
 *    and end[] are unaffected, though other effects of wing
 *    retraction take place; in an alignment, all B->M_k and M_k->E
 *    begins and exits are interpreted literally as B->M_k and M_k->E
 *    transitions. In ls/s mode models, there is no algorithmic
 *    internal entry/exit probability, but wing retraction contributes
 *    data-dependent probability to begin[] and [end]; in an
 *    alignment, B->M_k and M_k->E internal begins/ends are
 *    interpreted instead as delete paths B->D1->...->D_k-1->M_k or
 *    M_k->D_k+1->...->D_M->E.  That is, internal entry and exit probs
 *    may only be *entirely* algorithm dependent (from configuration)
 *    or *entirely* data-dependent (from wing retraction of the core
 *    model), not a mixture of both. (The reason is that if they're a
 *    sum of the two, you can't correctly extract a Viterbi alignment
 *    from the sum.) When the configuration is using
 *    algorithm-dependent entry, the PLAN7_BIMPOSED flag is raised;
 *    when it is using algorithm-dependent exit, the PLAN7_EIMPOSED
 *    flag is raised.  (xref STL9/79, and modelconfig.c).
 *    
 *    After configuration, the core model is unchanged; xt, begin, end
 *    are set to probabilities; and bsc, esc, tsc, msc, isc, and & xsc
 *    are set to scores. In a configured model, tbd1 is not used (all
 *    entry into the model is controlled by B->M_k begins); D_1 and
 *    D_M don't exist and aren't used; D_M-1->M_M = 1.0 and M_M-1->D_M
 *    = 0.0 (because D_M doesn't exist); and all exits from the model
 *    are controlled by M_k->E ends.
 *    
 *    The PLAN7_HASBITS flag is up when the model is configured into
 *    score form.
 *    
 * hmmbuild creates a core model, then configures it, then saves
 * it. Both the core model and the configured model are saved.
 * 
 * Search programs like hmmpfam and hmmsearch read only the configured
 * score model, usually ignoring the core probability model.
 */
typedef struct {
  /* Alphabet information (includes hmm->abc->K, alphabet size)
   */
  ESL_ALPHABET *abc;

  /* The core model in probability form.
   * P7_HASPROBS flag is raised when these probs are all valid.
   */
  /*::cexcerpt::plan7_core::begin::*/
  int     M;                    /* length of the model (# nodes)          */
  float **t;                    /* transition prob's. t[(0),1..M-1][0..6] */
  float **mat;                  /* match emissions.  mat[1..M][0..K-1]    */ 
  float **ins;                  /* insert emissions. ins[1..M-1][0..K-1]  */
  /*::cexcerpt::plan7_core::end::*/

  /* The null model probabilities.
   */
  float  *null;		         /* "random sequence" emission prob's     */
  float   p1;                    /* null model loop probability           */

  /* The unique states of Plan 7 in probability form.
   * These are the algorithm-dependent, data-independent probabilities.
   */
  float   xt[4][2];             /* N,E,C,J extra states: 2 transitions        +*/
  float  *begin;                /* 1..M B->M state transitions                +*/
  float  *end;                  /* 1..M M->E state transitions (!= a dist!)   +*/

  /* The model's log-odds score form: xref modelconfig.c
   *
   * Note that emission distributions are over possible alphabet
   * symbols, not just the unambiguous protein or DNA alphabet: we
   * precalculate the scores for all IUPAC degenerate symbols we may
   * see.
   *
   * Note the reversed indexing on msc, isc, tsc -- for efficiency
   * reasons. They're not probability vectors, so we can reorder them.
   * 
   * The _mem ptrs are where the real memory is alloc'ed and free'd,
   * as opposed to where it is accessed.  This came in with Erik
   * Lindahl's altivec port; it allows alignment on 16-byte
   * boundaries. In the non-altivec code, this is just a little
   * redundancy; tsc and tsc_mem point to the same thing, for example.
   * 
   * PLAN7_HASBITS flag is up when these scores are valid.
   */
  P7_PROFILE  *gm;
  P7_OPROFILE *om;
  float  lscore;		/* score for each unaligned res when LCORRECT flag is up */

  /* Annotation on the model. A name is mandatory.
   * Other fields are optional; whether they are present is
   * flagged in the stateflags bit array.
   * 
   * desc is only valid if PLAN7_DESC is set in flags.
   *  acc is only valid if PLAN7_ACC is set in flags.
   *   rf is only valid if PLAN7_RF is set in flags.
   *   cs is only valid if PLAN7_CS is set in flags.
   *   ca is only valid if PLAN7_CA is set in flags.
   *  map is only valid if PLAN7_MAP is set in flags.
   */
  char  *name;                  /* name of the model                    +*/
  char  *acc;			/* accession number of model (Pfam)     +*/
  char  *desc;                  /* brief description of model           +*/ 
  char  *rf;                    /* reference line from alignment 0..M   +*/
  char  *cs;                    /* consensus structure line      0..M   +*/ 
  char  *ca;			/* consensus accessibility line  0..M    */
  char  *comlog;		/* command line(s) that built model     +*/
  int    nseq;			/* number of training sequences         +*/
  char  *ctime;			/* creation date                        +*/
  int   *map;			/* map of alignment cols onto model 1..M+*/
  int    checksum;              /* checksum of training sequences       +*/

  /* Pfam-specific score cutoffs.
   * 
   * ga1, ga2 are valid if PLAN7_GA is set in flags.
   * tc1, tc2 are valid if PLAN7_TC is set in flags.
   * nc1, nc2 are valid if PLAN7_NC is set in flags.
   */
  float  ga1, ga2;		/* per-seq/per-domain gathering thresholds (bits) +*/
  float  tc1, tc2;		/* per-seq/per-domain trusted cutoff (bits)       +*/
  float  nc1, nc2;		/* per-seq/per-domain noise cutoff (bits)         +*/

  /* E-value statistical parameters.
   * If any are non-NULL, a flag is raised in <flags>: 
   * PLAN7_STATS_LV for lvstats, etc.
   */
  P7_EVINFO *lvstats;		/* local Viterbi  */
  P7_EVINFO *lfstats;		/* local Forward  */
  P7_EVINFO *gvstats;		/* glocal Viterbi */
  P7_EVINFO *gfstats;		/* glocal Forward */

  int flags;                    /* bit flags indicating state of HMM, valid data +*/
} P7_HMM;

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

/* Functions in plan7.c:  the model structure.
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

/* Functions in modelmakers.c:  constructing HMMs from alignments
 */
extern int p7_Handmodelmaker(ESL_MSA *msa, ESL_ALPHABET *abc, char **dsq, 
			     P7_HMM **ret_hmm, P7_TRACE ***ret_tr);



#endif /* PLAN7_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/

