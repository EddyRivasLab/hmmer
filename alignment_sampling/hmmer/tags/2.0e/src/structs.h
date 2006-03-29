/************************************************************
 * HMMER - Biological sequence analysis with profile-HMMs
 * Copyright (C) 1992-1998 Washington University School of Medicine
 *
 *   This source code is distributed under the terms of the 
 *   GNU General Public License. See the files COPYING and 
 *   GNULICENSE for details.
 *    
 ************************************************************/

/* structs.h
 * 
 * Data structures used in HMMER.
 * Also, a few miscellaneous macros and global variable declarations.
 * 
 * RCS $Id$
 */

#ifndef STRUCTSH_INCLUDED
#define STRUCTSH_INCLUDED

#include "squid.h"
#include "config.h"

/* Miscellaneous math macros used in the package
 */
#define sreLOG2(x)  ((x) > 0 ? log(x) * 1.44269504 : -9999.)
#define sreEXP2(x)  (exp((x) * 0.69314718 )) 
#define SQR(x)      ((x) * (x))

/* an idiom for determining a symbol's position in the array
 * by pointer arithmetic.
 * does no error checking, so caller must already be damned sure x is
 * valid in the alphabet!
 */
#define SYMIDX(x)   (strchr(Alphabet, (x)) - Alphabet)

/* Worry() is a conditional warning macro, used for debugging;
 * it inserts file name and line number information and
 * calls VerboseWorry(). See debug.c. 
 * 
 * PANIC is called for failures of Std C/POSIX functions,
 * instead of my own functions. It calls perror() and exits
 * abnormally.
 */
#define Worry(x,s)  VerboseWorry(x, __FILE__, __LINE__, s)
#define PANIC       (Panic(__FILE__, __LINE__))

/* The symbol alphabet.
 * Must deal with IUPAC degeneracies. Nondegenerate symbols 
 * come first in Alphabet[], followed by degenerate symbols.
 * Nucleic alphabet also must deal with other common symbols
 * like U (in RNA) and X (often misused for N).     
 * Example: 
 *   Nucleic: "ACGTUNRYMKSWHBVDX"          size=4  iupac=17
 *   Amino:   "ACDEFGHIKLMNPQRSTVWYBZX"    size=20 iupac=23
 *
 * Parts of the code assume that the last symbol is a
 * symbol for an unknown residue, i.e. 'X'.
 * 
 * MAXCODE and MAXABET constants are defined in config.h
 */   
extern char  Alphabet[MAXCODE]; /* "ACDEFGHIKLMNPQRSTVWYBZX" for example */
extern int   Alphabet_type;     /* hmmNUCLEIC or hmmAMINO                */
extern int   Alphabet_size;     /* uniq alphabet size: 4 or 20           */
extern int   Alphabet_iupac;    /* total size of alphabet + IUPAC degen. */
extern char  Degenerate[MAXCODE][MAXABET];
extern int   DegenCount[MAXCODE];
#define hmmNOTSETYET 0
#define hmmNUCLEIC   2		/* compatibility with squid's kDNA   */
#define hmmAMINO     3		/* compatibility with squid's kAmino */

/**********************************************************************
 *
 * Plan7 
 * Implementation of the new Plan7 HMM architecture.
 * Fully probabilistic even for hmmsw, hmmls, and hmmfs;
 * No insert->delete or delete->insert transitions;
 * Improved structure layout.
 * 
 * The strategy is to infiltrate plan7 code into HMMER in
 * an evolutionary rather than revolutionary manner. 
 *
 **********************************************************************/

/* Plan 7 construction strategies.
 */
enum p7_construction {
  P7_MAP_CONSTRUCTION,		/* maximum a posteriori architecture */
  P7_HAND_CONSTRUCTION,		/* hand specified architecture       */
  P7_FAST_CONSTRUCTION		/* fast ad hoc architecture          */
};

/* Plan 7 parameter optimization strategies
 */
enum p7_param {
  P7_MAP_PARAM,			/* standard maximum a posteriori    */
  P7_MD_PARAM,			/* maximum discrimination           */
  P7_MRE_PARAM,			/* maximum relative entropy         */
  P7_WMAP_PARAM			/* ad hoc weighted MAP              */
};

/* Structure: plan7_s
 * 
 * Declaration of a Plan 7 profile-HMM.
 */
struct plan7_s {
  /* Annotation on the model. A name is mandatory.
   * Other fields are optional; whether they are present is
   * flagged in the stateflags bit array.
   * 
   * desc is only valid if PLAN7_DESC is set.
   *   rf is only valid if PLAN7_RF is set.
   *   cs is only valid if PLAN7_CS is set.
   */
  char  *name;                  /* name of the model                    +*/
  char  *desc;                  /* brief description of model           +*/ 
  char  *rf;                    /* reference line from alignment 0..M   +*/
  char  *cs;                    /* consensus structure line      0..M   +*/ 
  char  *comlog;		/* command line(s) that built model     +*/
  int    nseq;			/* number of training sequences         +*/
  char  *ctime;			/* creation date                        +*/

  /* The main model in probability form: data-dependent probabilities.
   * This is the core Krogh/Haussler model.
   * Transition probabilities are usually accessed as a
   *   two-D array: hmm->t[k][TMM], for instance. They are allocated
   *   such that they can also be stepped through in 1D by pointer
   *   manipulations, for efficiency in DP algorithms.
   */
  int     M;                    /* length of the model (# nodes)        +*/
  float **t;                    /* transition prob's. t[1..M-1][0..6]   +*/
  float **mat;                  /* match emissions.  mat[1..M][0..19]   +*/ 
  float **ins;                  /* insert emissions. ins[1..M-1][0..19] +*/
  float   tbd1;			/* B->D1 prob (data dependent)          +*/

  /* The unique states of Plan 7 in probability form.
   * These are the algorithm-dependent, data-independent probabilities.
   */
  float  xt[4][2];              /* N,E,C,J extra states: 2 transitions      +*/
  float *begin;                 /* 1..M B->M state transitions              +*/
  float *end;                   /* 1..M M->E state transitions (!= a dist!) +*/

  /* The null model probabilities.
   */
  float  null[MAXABET];         /* "random sequence" emission prob's     +*/
  float  p1;                    /* null model loop probability           +*/

  /* The model in log-odds score form.
   * These are created from the probabilities by LogoddsifyHMM().
   * By definition, null[] emission scores are all zero.
   * Note that emission distributions are over 26 upper-case letters,
   * not just the unambiguous protein or DNA alphabet: we
   * precalculate the scores for all IUPAC degenerate symbols we
   * may see. Non-IUPAC symbols simply have a -INFTY score.
   * Note the reversed indexing on msc and isc -- for efficiency reasons.
   * 
   * Only valid if PLAN7_HASBITS is set.
   */
  int  **tsc;                   /* transition scores     [1.M-1][0.6]       -*/
  int  **msc;                   /* match emission scores [0.MAXCODE-1][1.M] -*/
  int  **isc;                   /* ins emission scores [0.MAXCODE-1][1.M-1] -*/
  int    xsc[4][2];             /* N,E,C,J transitions                      -*/
  int   *bsc;                   /* begin transitions     [1.M]              -*/
  int   *esc;			/* end transitions       [1.M]              -*/

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
  float  mu;			/* EVD mu       +*/
  float  lambda;		/* EVD lambda   +*/
  float  wonka;			/* EVD fit display fudge factor +*/

  int flags;                    /* bit flags indicating state of HMM    +*/
};

/* Flags for plan7->flags.
 * Note: Some models have scores but no probabilities (for instance,
 *       after reading from an HMM save file). Other models have
 *       probabilities but no scores (for instance, during training
 *       or building). Since it costs time to convert either way,
 *       I use PLAN7_HASBITS and PLAN7_HASPROB flags to defer conversion
 *       until absolutely necessary. This means I have to be careful
 *       about keeping these flags set properly when I fiddle a model. 
 */
#define PLAN7_HASBITS (1<<0)    /* raised if model has log-odds scores      */
#define PLAN7_DESC    (1<<1)    /* raised if description exists             */
#define PLAN7_RF      (1<<2)    /* raised if #RF annotation available       */
#define PLAN7_CS      (1<<3)    /* raised if #CS annotation available       */
#define PLAN7_XRAY    (1<<4)    /* raised if structural data available      */
#define PLAN7_HASPROB (1<<5)    /* raised if model has probabilities        */
#define PLAN7_HASDNA  (1<<6)	/* raised if protein HMM->DNA seq params set*/
#define PLAN7_STATS   (1<<7)	/* raised if EVD parameters are available   */

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

/* Declaration of Plan7 dynamic programming matrix structure.
 */
struct dpmatrix_s {
  int **xmx;			/* special scores [0.1..N][BECJN]     */
  int **mmx;			/* match scores [0.1..N][0.1..M]      */
  int **imx;			/* insert scores [0.1..N][0.1..M-1.M] */
  int **dmx;			/* delete scores [0.1..N][0.1..M-1.M] */
};

/* Structure: HMMFILE
 * 
 * Purpose:   An open HMM file or HMM library. See hmmio.c.
 */
struct hmmfile_s {
  FILE *f;			/* pointer to file opened for reading     */
  int (*parser)(struct hmmfile_s *, struct plan7_s **);  /* parsing function*/
  int   is_binary;		/* TRUE if format is a binary one         */
  int   byteswap;               /* TRUE if binary and byteswapped         */
};
typedef struct hmmfile_s HMMFILE; 


/* Plan 7 model state types
 * used in traceback structure
 */
enum p7stype { STM, STD, STI, STS, STN, STB, STE, STC, STT, STJ, STBOGUS };

/* Structure: p7trace_s
 * 
 * Traceback structure for alignments of model to sequence.
 * Each array in a trace_s is 0..tlen-1.
 * Element 0 is always to STATE_S. Element tlen-1 is always to STATE_T. 
 */
struct p7trace_s {
  int   tlen;                   /* length of traceback                      */
  enum p7stype *statetype;      /* state type used for alignment            */
  int  *nodeidx;                /* index of aligned node, 1..M (if M,D,I)   */
  int  *pos;                    /* position in dsq, 1..L or -1              */ 
};

/* Structure: p7prior_s
 * 
 * Dirichlet priors on HMM parameters.
 */
struct p7prior_s {
  int   strategy;               /* PRI_DCHLET, etc.                          */

  int   tnum;                   /* number of transition Dirichlet mixtures   */
  float tq[MAXDCHLET];          /* probabilities of tnum components          */
  float t[MAXDCHLET][7];        /* transition terms per mix component        */

  int   mnum;                   /* number of mat emission Dirichlet mixtures */
  float mq[MAXDCHLET];          /* probabilities of mnum components          */
  float m[MAXDCHLET][MAXABET];  /* match emission terms per mix component    */

  int   inum;			/* number of insert emission Dirichlet mixes */
  float iq[MAXDCHLET];		/* probabilities of inum components          */
  float i[MAXDCHLET][MAXABET];	/* insert emission terms                     */
};
#define PRI_DCHLET 0            /* simple or mixture Dirichlets   */
#define PRI_PAM    1		/* PAM prior hack                 */


/**********************************************************************
 * Other structures, not having to do with HMMs.
 **********************************************************************/

/* Structure: histogram_s 
 * 
 * Keep a score histogram. 
 * 
 * The main implementation issue here is that the range of
 * scores is unknown, and will go negative. histogram is
 * a 0..max-min array that represents the range min..max.
 * A given score is indexed in histogram array as score-min.
 * The AddToHistogram() function deals with dynamically 
 * resizing the histogram array when necessary.
 */  
struct histogram_s {
  int *histogram;		/* counts of hits                     */
  int  min;			/* elem 0 of histogram == min         */
  int  max;                     /* last elem of histogram == max      */
  int  highscore;		/* highest active elem has this score */
  int  lowscore;		/* lowest active elem has this score  */
  int  lumpsize;		/* when resizing, overalloc by this   */
  int  total;			/* total # of hits counted            */

  float *expect;		/* expected counts of hits            */
  int    fit_type;		/* flag indicating distribution type  */
  float  param[3];		/* parameters used for fits           */
  float  chisq;			/* chi-squared val for goodness of fit*/
  float  chip;			/* P value for chisquared             */
};
#define HISTFIT_NONE     0	/* no fit done yet               */
#define HISTFIT_EVD      1	/* fit type = extreme value dist */
#define HISTFIT_GAUSSIAN 2	/* fit type = Gaussian           */
#define EVD_MU		 0	/* EVD fit parameter mu          */
#define EVD_LAMBDA       1	/* EVD fit parameter lambda      */
#define EVD_WONKA        2      /* EVD fit fudge factor          */
#define GAUSS_MEAN       0	/* Gaussian parameter mean       */
#define GAUSS_SD         1	/* Gaussian parameter std. dev.  */

/* Structure: fancyali_s
 * 
 * Alignment of a hit to an HMM, for printing.
 */
struct fancyali_s {
  char *rfline;                 /* reference coord info                 */
  char *csline;                 /* consensus structure info             */
  char *model;                  /* aligned query consensus sequence     */
  char *mline;                  /* "identities", conservation +'s, etc. */
  char *aseq;                   /* aligned target sequence              */
  int   len;			/* length of strings                    */
  char *query;			/* name of query HMM                    */
  char *target;			/* name of target sequence              */
  int   sqfrom;			/* start position on sequence (1..L)    */
  int   sqto;		        /* end position on sequence   (1..L)    */
};

/* Structure: hit_s
 * 
 * Info about a high-scoring database hit.
 * We keep this info in memory, so we can output a
 * sorted list of high hits at the end.
 *
 * sqfrom and sqto are the coordinates that will be shown
 * in the results, not coords in arrays... therefore, reverse
 * complements have sqfrom > sqto
 */
struct hit_s {
  double sortkey;		/* number to sort by; big is better */
  float  score;			/* score of the hit                 */
  double pvalue;		/* P-value of the hit               */
  float  mothersc;		/* score of whole sequence          */
  double motherp;		/* P-value of whole sequence        */
  char   *name;			/* name of the database seq         */
  char   *desc;			/* description of database seq      */
  int    sqfrom;		/* start position in seq (1..N)     */
  int    sqto;			/* end position in seq (1..N)       */
  int    sqlen;			/* length of sequence (N)           */
  int    hmmfrom;		/* start position in HMM (1..M)     */
  int    hmmto;			/* end position in HMM (1..M)       */
  int    hmmlen;		/* length of HMM (M)                */
  int    domidx;		/* index of this domain             */
  int    ndom;			/* total # of domains in this seq   */
  struct fancyali_s  *ali;	/* ptr to optional alignment info   */
};


/* Structure: tophit_s
 * 
 * Array of high scoring hits, suitable for efficient sorting
 * when we prepare to output results. "hit" list is NULL and 
 * unavailable until after we do a sort.
 */
struct tophit_s {
  struct hit_s **hit;           /* array of ptrs to top scoring hits        */
  struct hit_s  *unsrt;         /* unsorted array                           */
  int            alloc;		/* current allocation size                  */
  int            num;		/* number of hits in list now               */
  int            lump;       	/* allocation lumpsize                      */
};


/**********************************************************
 * BLAST-HMMs.
 *
 * The following structures and definitions are relevant to 
 * the implementation of BLAST algorithms for HMM searching.
 **********************************************************/

/* Struct: hmmword_s
 * 
 * Information about a single HMM neighborhood word.
 */
struct hmmword_s {
  char word[BLAST_MAXWORD];     /* in unambiguous DNA or protein seq  */
  int   len;                    /* length of the word */   
  int   startk;			/* starting node on model */
  int   starty;			/* starts on MATCH, DELETE, or INSERT state */
  int   currk;			/* ending node on model   */
  int   curry;			/* ends on MATCH, DELETE, or INSERT state */
  float logp;			/* current logp for prefix. */
  float bound;			/* logp + bck: best logp attainable, upper bound */
};

/* Struct: wordheap_s
 * 
 * Manages a heap of neighborhood words. Sometimes this structure
 * is really a heap. Sometimes it is just a pushdown stack.
 */
struct wordheap_s {
  int                alloc;	/* max N allocated for          */
  int                N;		/* current number of nodes used */
  struct hmmword_s **node;      /* array of ptrs to words       */ 
};

/* Struct: wordpool_s
 * 
 * Manages allocation of lots of words. Implemented
 * as an array of arrays of words and a pushdown stack of
 * valid pointers. "Freeing" a pointer just pushes it
 * back into the pool.
 */
struct wordpool_s {
  int lev;			/* which level of the pool we're on  */
  int N;			/* index of next ptr to give away    */
  struct hmmword_s **node;      /* 2D array of nodes                 */
  struct hmmword_s **ptr;       /* array of active pointers to give out */
};


/**********************************************************
 * Plan 9: obsolete (?) HMMER1 code
 **********************************************************/

/* We define a "basic" state, which covers the basic match, insert, and
 * delete states from the Haussler paper. Numbers are stored as
 * pre-calculated negative logs.
 */
struct basic_state {
  float t[3];			/* state transitions to +1 M, +0 I, +1 D */
  float p[MAXABET];            	/* symbol emission probabilities         */
};

/* A complete hidden Markov model
 */
struct hmm_struc {
  int    M;			/* length of the model                */
  struct basic_state *ins;      /* insert states 0..M+1               */
  struct basic_state *mat;      /* match 0..M+1; 0 = BEGIN, M+1 = END */
  struct basic_state *del;      /* delete 0..M+1                      */

  float  null[MAXABET];         /* the *suggested* null model         */

  /* Optional annotation on the HMM, taken from alignment
   */
  char  *name;                  /* a name for the HMM                     */
  char  *ref;			/* reference coords and annotation        */
  char  *cs;                    /* consensus structure annotation         */         
  float *xray;	/* Structural annotation: xray[0..M+1][NINPUTS], indexed manually */

  int    flags;			/* flags for what optional info is in HMM */
};

/* Flags for optional info in an HMM structure
 */
#define HMM_REF   (1<<0) 
#define HMM_CS    (1<<1)
#define HMM_XRAY  (1<<2)

/* Array indices for structural inputs in HMM 
 */
#define XRAY_bias  0		/* bias: always 1 */
#define XRAY_E     1		/* P(sheet), 0..1 */
#define XRAY_H     2		/* P(helix), 0..1 */
#define XRAY_SA    3		/* relative side chain solv. access. (0..15) */


/* We have to increase the dynamic range of the calculations to floats
 * when we do simulated annealing.
 * Rather than sacrificing speed in the rest of the code, we define
 * new HMM structures specifically for simulated annealing use. 
 */
struct sa_state_s {
  double t[3];			/* transition probability to M, I, D */
  double p[MAXABET];		/* symbol emission probabilities     */
};

struct sa_hmm_s {
  int    M;			/* length of the model            */
  struct sa_state_s *ins;       /* insert states 0..M+1           */
  struct sa_state_s *mat;       /* match 0..M+1; 0 = impossible   */
  struct sa_state_s *del;       /* delete 0..M+1; 0 = BEGIN state */
}; 


/* Search HMM structure:
 * We use a special structure to store HMM's in integer log odds
 * form, for use in the various Viterbi alignment procedures.
 * The format is very compact and a bit fragile because this
 * is optimized for speed & memory access patterns.
 */
struct shmm_s {
  int  M;                       /* length of the model           */
  int *m_emit[26];		/* 26 x M+1 emission scores      */
  int *i_emit[26];		/* 26 x M+1 insertion scores     */
  int *t;                       /* 9 x M+1 state transitions: 
				   order is dd,di,dm,id,ii,im,md,mi,mm */
  /* plus extra annotation:
   */
  int   flags;
  char *ref;
  char *cs;
  char *name;
};
/* Order of transition probabilities in shmm_s
 */
#define Tdd 0
#define Tdi 1
#define Tdm 2
#define Tid 3
#define Tii 4
#define Tim 5
#define Tmd 6
#define Tmi 7
#define Tmm 8


/* Prior information: expectations for state transition probabilities,
 * emission probabilities, and alphas (regularizers) describing the 
 * strengths of belief in the prior. alphas of 1 result in Laplace
 * small sample size corrections; very high alphas result in "hard-wiring"
 * of probabilities to the prior expectation.
 */
struct prior_s {
  int   strategy;		/* PRI_SIMPLE, etc.                          */

  int   tnum;			/* number of transition Dirichlet mixtures   */
  float tw[MAXDCHLET][NINPUTS];/* weight matrix for struct prior perceptron */
  float tq[MAXDCHLET];		/* probabilities of tnum components          */
  float tm[MAXDCHLET][3];	/* mat transition terms per mix component    */
  float td[MAXDCHLET][3];	/* del transition terms per mix component    */
  float ti[MAXDCHLET][3];	/* ins transition terms per mix component    */

  int   mnum;			/* number of mat emission Dirichlet mixtures */
  float mw[MAXDCHLET][NINPUTS]; /* weight matrix for struct prior perceptron */
  float mq[MAXDCHLET];		/* probabilities of mnum components          */
  float mat[MAXDCHLET][MAXABET];/* match emission terms per mix component    */

  int   inum;			/* number of ins emission Dirichlet mixtures */
  float iw[MAXDCHLET][NINPUTS]; /* weight matrix for struct prior perceptron */
  float iq[MAXDCHLET];		/* probabilities of inum components          */
  float ins[MAXDCHLET][MAXABET];/* insert emission m_x per mix component     */
};

#define PRI_SIMPLE 0		/* simple single-component Dirichlets (orig) */
#define PRI_PAM    1		/* ad hoc PAM-based mixture prior            */
#define PRI_MIX    2		/* mixture prior, a la Brown/Haussler        */
#define PRI_STRUCT 3		/* structure annotation based priors         */



/* Traceback structure for alignments of model to sequence.
 * Each array in a trace_s is 0..tlen-1.
 * Element 0 and tlen-1 are dummy "alignments" to BEGIN and END.
 */
struct trace_s {
  int   tlen;			/* length of traceback                       */
  int  *nodeidx;                /* index of aligned node, 0..hmm->M+1        */
  char *statetype;              /* state type used for alignment             */
  int  *rpos;                   /* position in raw sequence, 0..slen-1 or -1 */ 
};


/* Bookkeeping for individual cells of the dynamic programming calculations
 */
				/* used by most Viterbi procedures */
struct vit_s {
  int  score_m;			/* score to match here  */
  int  score_d;			/* score to delete here */
  int  score_i;			/* score to insert here */
};

				/* used by simulated annealing */
struct sa_s {
  double score_m;
  double score_d;
  double score_i;
};

				/* used by forward-backwards */
struct forback_s {
  float score_m;
  float score_d;
  float score_i;
};

				/* used in fragviterbi.c (hmmfs) */
struct fvit_s {
  int  score_m;			/* score to match here   */
  int  score_d;			/* score to delete here  */
  int  score_i;			/* score to insert here  */
  int  tback_m;			/* traceback from match  */
  int  tback_d;			/* traceback from delete */
  int  tback_i;			/* traceback from insert */
};


/* some define's for storing backtracking info
 */
#define FROM_NOWHERE 127
#define FROM_MATCH   0
#define FROM_INSERT  1
#define FROM_DELETE  2

/* same numbers, different names -- just for readability, makes
 * more semantic sense after the tracebacks
 */
#define MATCH  FROM_MATCH
#define INSERT FROM_INSERT
#define DELETE FROM_DELETE

#define BEGIN  MATCH
#define END    MATCH

#endif /* STRUCTSH_INCLUDED */
