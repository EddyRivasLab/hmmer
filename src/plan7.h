/* plan7.h
 * The Plan7 HMM structure used by HMMER.
 * 
 * SRE, Sat Apr 30 13:05:35 2005
 * SVN $Id$
 */
#ifndef PLAN7_INCLUDED
#define PLAN7_INCLUDED

/* Plan 7 algorithm modes
 */
enum p7_algmode {
  P7_NO_MODE    = 0,
  P7_LS_MODE    = 1,
  P7_FS_MODE    = 2,
  P7_SW_MODE    = 3,
  P7_S_MODE     = 4,
};

/* Structure: plan7_s
 * 
 * Declaration of a Plan 7 profile-HMM structure.
 * 
 * The model has two forms:
 * 1. The "core" model is the probability form.
 *    It has 1..M nodes, and one additional <tbd1> parameter.
 *    Nodes 1..M-1 have M,D,I states; node M has M,D states. 
 *    
 *    t[1..M-1] are the state transition probs. (t[M] are special,
 *      because this node transits to the end. M_M->E and D_M->E
 *      are implicitly 1.0. All other transitions are implicitly 0.0.)
 *    
 *    mat[1..M] are match emission probs.
 *    ins[1..M-1] are insert emission probs. (No I_M state.)
 *    
 *    Entry into the first node is controlled by <tbd1>: B->M1 is
 *      (1-tbd1), and B->D1 is (tbd1). 
 *    
 *    The PLAN7_HASPROB flag is up when these all correspond to a fully normalized
 *    profile HMM.
 *    
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
struct plan7_s {
  /* The core model in probability form: data-dependent probabilities.
   * Transition probabilities are usually accessed as a
   *   two-D array: hmm->t[k][TMM], for instance. They are allocated
   *   such that they can also be stepped through in 1D by pointer
   *   manipulations, for efficiency in DP algorithms.
   * PLAN7_HASPROBS flag is raised when these probs are all valid.
   */
  int     M;                    /* length of the model (# nodes)        +*/
  float **t;                    /* transition prob's. t[1..M-1][0..6]   +*/
  float **mat;                  /* match emissions.  mat[1..M][0..19]   +*/ 
  float **ins;                  /* insert emissions. ins[1..M-1][0..19] +*/
  float   tbd1;			/* B->D1 prob (data dependent)          +*/
  float   tbm1;			/* B->M1 prob (data dependent)          +*/

  /* The null model probabilities.
   */
  float  null[MAXABET];         /* "random sequence" emission prob's     +*/
  float  p1;                    /* null model loop probability           +*/

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
  enum p7_algmode mode;		/* configured algorithm mode  */
  int  **tsc;                   /* transition scores     [0.6][1.M-1]       +*/
  int  **msc;                   /* match emission scores [0.MAXCODE-1][1.M] +*/
  int  **isc;                   /* ins emission scores [0.MAXCODE-1][1.M-1] +*/
  int    xsc[4][2];             /* N,E,C,J transitions                      +*/
  int   *bsc;                   /* begin transitions     [1.M]              +*/
  int   *esc;			/* end transitions       [1.M]              +*/
  int   *tsc_mem, *msc_mem, *isc_mem, *bsc_mem, *esc_mem;
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

  /* The following are annotations added to support work by Michael Asman, 
   * CGR Stockholm. They are not stored in model files; they are only
   * used in model construction.
   * 
   * #=GC X-PRM (PRT,PRI) annotation is picked up by hmmbuild and interpreted
   * as specifying which mixture Dirichlet component to use. If these flags
   * are non-NULL, the normal mixture Dirichlet code is bypassed, and a
   * single specific Dirichlet is used at each position.
   */
  int   *tpri;                  /* which transition mixture prior to use */ 
  int   *mpri;                  /* which match mixture prior to use */
  int   *ipri;                  /* which insert mixture prior to use */

  /* Pfam-specific score cutoffs.
   * 
   * ga1, ga2 are valid if PLAN7_GA is set in flags.
   * tc1, tc2 are valid if PLAN7_TC is set in flags.
   * nc1, nc2 are valid if PLAN7_NC is set in flags.
   */
  float  ga1, ga2;		/* per-seq/per-domain gathering thresholds (bits) +*/
  float  tc1, tc2;		/* per-seq/per-domain trusted cutoff (bits)       +*/
  float  nc1, nc2;		/* per-seq/per-domain noise cutoff (bits)         +*/

  /* DNA translation scoring parameters
   * For aligning protein Plan7 models to DNA sequence.
   * Lookup value for a codon is calculated by pos1 * 16 + pos2 * 4 + pos3,
   * where 'pos1' is the digitized value of the first nucleotide position;
   * if any of the positions are ambiguous codes, lookup value 64 is used
   * (which will generally have a score of zero)
   * 
   * Only valid if PLAN7_HASDNA is set.
   */
  int  **dnam;                  /* triplet match scores  [0.64][1.M]       -*/
  int  **dnai;                  /* triplet insert scores [0.64][1.M]       -*/
  int    dna2;			/* -1 frameshift, doublet emission, M or I -*/
  int    dna4;			/* +1 frameshift, doublet emission, M or I -*/

  /* P-value and E-value statistical parameters
   * Only valid if PLAN7_STATS is set.
   */
  float  mu;			/* Gumbel mu                                +*/
  float  lambda;		/* Gumbel lambda                            +*/
  float  kappa;			/* mean length of optimal local alignments  +*/
  float  kappa_g;		/* mean length of optimal glocal alignments +*/

  int flags;                    /* bit flags indicating state of HMM, valid data +*/
};

/* Flag codes for plan7->flags.
 */
#define PLAN7_HASBITS (1<<0)    /* model has log-odds scores            */
#define PLAN7_DESC    (1<<1)    /* description exists                   */
#define PLAN7_RF      (1<<2)    /* #RF annotation available             */
#define PLAN7_CS      (1<<3)    /* #CS annotation available             */
#define PLAN7_XRAY    (1<<4)    /* structural data available            */
#define PLAN7_HASPROB (1<<5)    /* model has probabilities              */
#define PLAN7_HASDNA  (1<<6)	/* protein HMM->DNA seq params set      */
#define PLAN7_STATS   (1<<7)	/* EVD parameters are available         */
#define PLAN7_MAP     (1<<8)	/* alignment map is available           */
#define PLAN7_ACC     (1<<9)	/* accession number is available        */
#define PLAN7_GA      (1<<10)	/* gathering thresholds available       */
#define PLAN7_TC      (1<<11)	/* trusted cutoffs available            */
#define PLAN7_NC      (1<<12)	/* noise cutoffs available              */
#define PLAN7_CA      (1<<13)   /* surface accessibility avail.         */
#define PLAN7_BIMPOSED (1<<14)  /* all entries are B->M_k (not D)       */
#define PLAN7_EIMPOSED (1<<15)  /* all ends are M_k->E (not D)          */
#define PLAN7_LCORRECT (1<<16)  /* require L-dependent score correction */

/* Indices for special state types, I: used for dynamic programming xmx[][]
 * mnemonic: eXtra Matrix for B state = XMB
 */
#define XMB 0
#define XME 1
#define XMC 2
#define XMJ 3
#define XMN 4

/* Indices for special state types, II: used for hmm->xt[] indexing
 * mnemonic: eXtra Transition for N state = XTN
 */
#define XTN  0
#define XTE  1
#define XTC  2
#define XTJ  3

/* Indices for Plan7 main model state transitions.
 * Used for indexing hmm->t[k][]
 * mnemonic: Transition from Match to Match = TMM
 */
#define TMM  0
#define TMI  1
#define TMD  2
#define TIM  3
#define TII  4
#define TDM  5
#define TDD  6 

/* Indices for extra state transitions
 * Used for indexing hmm->xt[][].
 */
#define MOVE 0          /* trNB, trEC, trCT, trJB */
#define LOOP 1          /* trNN, trEJ, trCC, trJJ */


#endif /* PLAN7_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/

