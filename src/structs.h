/* structs.h
 * Data structures used in HMMER.
 * Also, a few miscellaneous macros and global variable declarations.
 * 
 * SVN $Id$
 */

#ifndef STRUCTSH_INCLUDED
#define STRUCTSH_INCLUDED

#include "config.h"

#include "squid.h"
#include "ssi.h"

#include "customstructs.h"
#include "plan7.h"


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
extern char  Alphabet[MAXCODE+1]; /* "ACDEFGHIKLMNPQRSTVWYBZX" for example */
extern int   Alphabet_type;       /* hmmNUCLEIC or hmmAMINO                */
extern int   Alphabet_size;       /* uniq alphabet size: 4 or 20           */
extern int   Alphabet_iupac;      /* total size of alphabet + IUPAC degen. */
extern char  Degenerate[MAXCODE][MAXABET];
extern int   DegenCount[MAXCODE];
#define hmmNOTSETYET 0
#define hmmNUCLEIC   2		/* compatibility with squid's kRNA   */
#define hmmAMINO     3		/* compatibility with squid's kAmino */

/* Declaration of default Plan7 dynamic programming matrix structure.
 */
struct dpmatrix_s {
  int **xmx;			/* special scores [0.1..N][BECJN]     */
  int **mmx;			/* match scores [0.1..N][0.1..M]      */
  int **imx;			/* insert scores [0.1..N][0.1..M-1.M] */
  int **dmx;			/* delete scores [0.1..N][0.1..M-1.M] */

  /* The other trick brought in w/ the Lindahl Altivec port; dp matrix
   * is retained and grown, rather than reallocated for every HMM or sequence.
   * Keep track of current allocated-for size in rows (sequence length N)
   * and columns (HMM length M). Also keep track of pad sizes: how much
   * we should overallocate rows or columns when we reallocate. If pad = 0,
   * then we're not growable in this dimension.
   */
  int maxN;			/* alloc'ed for seq of length N; N+1 rows */
  int maxM;			/* alloc'ed for HMM of length M; M+1 cols */

  int padN;			/* extra pad in sequence length/rows */
  int padM;			/* extra pad in HMM length/columns   */
};

/* Declaration of Plan7 shadow matrix structure.
 * In general, allowed values are STM, STI, etc.
 * However, E state has M possible sources, from 1..M match states;
 * hence the esrc array.
 */
struct dpshadow_s {
  char **xtb;			/* special state traces [0.1..N][BECJN]     */
  char **mtb;			/* match state traces [0.1..N][0.1..M]      */
  char **itb;			/* insert state traces [0.1..N][0.1..M-1.M] */
  char **dtb;			/* delete state traces [0.1..N][0.1..M-1.M] */
  int   *esrc;                  /* E trace is special; must store a M state number 1..M */
};

/* Structure: HMMFILE
 * 
 * Purpose:   An open HMM file or HMM library. See hmmio.c.
 */
struct hmmfile_s {
  FILE    *f;			/* pointer to file opened for reading           */
  SSIFILE *ssi;			/* pointer to open SSI index, or NULL           */
  int (*parser)(struct hmmfile_s *, struct plan7_s **);  /* parsing function    */
  int   is_binary;		/* TRUE if format is a binary one               */
  int   byteswap;               /* TRUE if binary and byteswapped               */

  /* Ewan (GeneWise) needs the input API to know the offset of each
   * HMM on the disk, as it's being read. This might be enough
   * support for him. hmmindex also uses this. Ewan, see
   * HMMFilePositionByIndex() for an example of how to use this
   * opaque offset type in the SSI API - the call you need 
   * is SSISetFilePosition().
   */
  int       is_seekable;	/* TRUE if we use offsets in this HMM file  */
  int       mode;		/* type of offset                           */
  SSIOFFSET offset;		/* Disk offset for beginning of current HMM */
};
typedef struct hmmfile_s HMMFILE; 


/* Plan 7 model state types
 * used in traceback structure
 */
#define STBOGUS 0
#define STM     1
#define STD     2
#define STI     3
#define STS     4
#define STN     5
#define STB     6
#define STE     7
#define STC     8
#define STT     9
#define STJ     10     

/* Structure: p7trace_s
 * 
 * Traceback structure for alignments of model to sequence.
 * Each array in a trace_s is 0..tlen-1.
 * Element 0 is always to STATE_S. Element tlen-1 is always to STATE_T. 
 */
struct p7trace_s {
  int   tlen;                   /* length of traceback                      */
  char *statetype;              /* state type used for alignment            */
  int  *nodeidx;                /* idx of aligned node, 1..M if M,D,I; or 0 */
  int  *pos;                    /* position in dsq, 1..L, or 0 if none      */ 
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

  double *expect;		/* expected counts of hits            */
  int    fit_type;		/* flag indicating distribution type  */
  double param[3];		/* parameters used for fits           */
  double chisq;			/* chi-squared val for goodness of fit*/
  double chip;			/* P value for chisquared             */
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
  char   *name;			/* name of the target               */
  char   *acc;			/* accession of the target          */
  char   *desc;			/* description of the target        */
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

/* struct threshold_s 
 * Contains score/evalue threshold settings.
 *
 * made first for hmmpfam:
 * Since we're going to loop over all HMMs in a Pfam (or pfam-like)
 * database in main_loop_{serial,pvm}, and we're going to
 * allow autocutoffs using Pfam GA, NC, TC lines, we will need
 * to reset those cutoffs with each HMM in turn. Therefore the
 * main loops need to know whether they're supposed to be
 * doing autocutoff. This amount of info was unwieldy enough
 * to pass through the argument list that I put it
 * in a structure.
 */
struct threshold_s {
  float  globT;			/* T parameter: keep only hits > globT bits */
  double globE;			/* E parameter: keep hits < globE E-value   */
  float  domT;			/* T parameter for individual domains       */
  double domE;			/* E parameter for individual domains       */
				/* autosetting of cutoffs using Pfam annot: */
  enum { CUT_NONE, CUT_GA, CUT_NC, CUT_TC } autocut;
  int   Z;			/* nseq to base E value calculation on      */
};

/**********************************************************
 * PVM parallelization
 **********************************************************/
#ifdef  HMMER_PVM

/* Message tags 
 */
#define HMMPVM_INIT         0	/* an initialization packet to all slaves */
#define HMMPVM_WORK         1	/* a work packet sent to a slave */
#define HMMPVM_RESULTS      2	/* a results packet sent back to master */
#define HMMPVM_TASK_TROUBLE 3	/* a notification of bad things in a slave task */
#define HMMPVM_HOST_TROUBLE 4	/* a notification of bad things in a PVM host */

/* error codes
 */
#define HMMPVM_OK         0
#define HMMPVM_NO_HMMFILE 1
#define HMMPVM_NO_INDEX   2	
#define HMMPVM_BAD_INIT   3	/* failed to initialize a slave somehow */

#endif


/**********************************************************
 * Plan 9: obsolete HMMER1.x code. We still need these structures
 * for reading old HMM files (e.g. backwards compatibility)
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
struct plan9_s {
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

#define MATCH  0
#define INSERT 1
#define DELETE 2
#define BEGIN  MATCH
#define END    MATCH

#endif /* STRUCTSH_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/

