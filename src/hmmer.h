/* The all-encompassing include file for HMMER.
 *
 *    1. P7_HMM:         a core model.
 *    2. P7_PROFILE:     a scoring profile, and its implicit model.
 *    3. P7_BG:          a null (background) model.
 *    4. P7_TRACE:       a traceback path (alignment of seq to profile).
 *    5. P7_HMMFILE:     an HMM save file or database, open for reading.
 *    6. P7_GMX:         a "generic" dynamic programming matrix
 *    7. P7_DPRIOR:      mixture Dirichlet prior for profile HMMs
 *    8. P7_SPENSEMBLE:  segment pair ensembles for domain locations
 *    9. P7_ALIDISPLAY:  an alignment formatted for printing
 *   10. P7_DOMAINDEF:   reusably managing workflow in annotating domains
 *   11. P7_TOPHITS:     ranking lists of top-scoring hits
 *   12. Inclusion of the architecture-specific optimized implementation.
 *   13. Declaration of functions in HMMER's exposed API.
 *   14. Copyright and license information.
 * 
 * SRE, Wed Jan  3 13:46:42 2007 [Janelia] [Philip Glass, The Fog of War]
 * SVN $Id$
 */
#ifndef P7_HMMERH_INCLUDED
#define P7_HMMERH_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"
#include "esl_alphabet.h"	/* ESL_DSQ, ESL_ALPHABET */
#include "esl_dmatrix.h"	/* ESL_DMATRIX           */
#include "esl_msa.h"		/* ESL_MSA               */
#include "esl_random.h"		/* ESL_RANDOMNESS        */
#include "esl_sq.h"		/* ESL_SQ                */
#include "esl_histogram.h"      /* ESL_HISTOGRAM         */
#include "esl_dirichlet.h"	/* ESL_MIXDCHLET         */

/* Search modes.
 */
#define p7_NO_MODE   0
#define p7_LOCAL     1		/* multihit local:  "fs" mode   */
#define p7_GLOCAL    2		/* multihit glocal: "ls" mode   */
#define p7_UNILOCAL  3		/* unihit local: "sw" mode      */
#define p7_UNIGLOCAL 4		/* unihit glocal: "s" mode      */

#define p7_IsLocal(mode)  (mode == p7_LOCAL || mode == p7_UNILOCAL)
#define p7_IsMulti(mode)  (mode == p7_LOCAL || mode == p7_GLOCAL)


/*****************************************************************
 * 1. P7_HMM: a core model.
 *****************************************************************/

/* Flag codes for hmm->flags.
 * Flags marked with ! may not be changed nor used for other meanings;
 * such flags were stored in old HMM files, and we must preserve their
 * meaning to preserve reverse compatibility.
 */
#define p7H_HASBITS (1<<0)    /* obsolete (was: model has log-odds scores)       !*/
#define p7H_DESC    (1<<1)    /* description exists                              !*/
#define p7H_RF      (1<<2)    /* #RF annotation available                        !*/
#define p7H_CS      (1<<3)    /* #CS annotation available                        !*/
#define p7H_XRAY    (1<<4)    /* obsolete (was: structural data available)       !*/
#define p7H_HASPROB (1<<5)    /* obsolete (was: model in probability form)       !*/
#define p7H_HASDNA  (1<<6)    /* obsolete (was: protein HMM->DNA seq params set) !*/
#define p7H_STATS   (1<<7)    /* model has EVD stats calibrated                  !*/
#define p7H_MAP     (1<<8)    /* alignment map is available                      !*/
#define p7H_ACC     (1<<9)    /* accession number is available                   !*/
#define p7H_GA      (1<<10)   /* gathering thresholds available                  !*/
#define p7H_TC      (1<<11)   /* trusted cutoffs available                       !*/
#define p7H_NC      (1<<12)   /* noise cutoffs available                         !*/
#define p7H_CA      (1<<13)   /* surface accessibilities available               !*/

/* Indices of Plan7 main model state transitions, hmm->t[k][]
 */
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

/* How the hmm->t[k] vector is interpreted as separate probability vectors.
 */
#define P7H_TMAT(hmm, k) ((hmm)->t[k])
#define P7H_TINS(hmm, k) ((hmm)->t[k]+3)
#define P7H_TDEL(hmm, k) ((hmm)->t[k]+5)
#define p7H_NTMAT 3
#define p7H_NTDEL 2
#define p7H_NTINS 2


/* Some notes:
 *   0. The model might be either in counts or probability form.
 *   1. t[0] is special: t[0][TMM,TMI,TMD] are the begin->M_1,I_0,D_1 entry probabilities,
 *      t[0][TIM,TII] are the I_0 transitions, and delete state 0 doesn't
 *      exist. Therefore D[0] transitions and mat[0] emissions are unused.
 *      To simplify some normalization code, we adopt a convention that these are set
 *      to valid probability distributions: 1.0 for t[0][TDM] and mat[0][0],
 *      and 0 for the rest.
 *   2. t[M] is also special: TMD and TDD are 0 because there is no next delete state;
 *      TDM is therefore 1.0 by definition. TMM and TDM are interpreted as the
 *      M->E and D->E end transitions. t[M][TDM] must be 1.0, therefore.
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
  char  *name;       /* name of the model                     (mandatory) */ /* String, \0-terminated */
  char  *acc;	     /* accession number of model (Pfam)      (p7_ACC)    */ /* String, \0-terminated */
  char  *desc;       /* brief (1-line) description of model   (p7_DESC)   */ /* String, \0-terminated */
  char  *rf;         /* reference line from alignment 1..M    (p7_RF)     */ /* String; 0=' ', M+1='\0' */
  char  *cs;         /* consensus structure line      1..M    (p7_CS)     */ /* String; 0=' ', M+1='\0' */
  char  *ca;	     /* consensus accessibility line  1..M    (p7_CA)     */ /* String; 0=' ', M+1='\0' */
  char  *comlog;     /* command line(s) that built model      (mandatory) */ /* String, \0-terminated */
  int    nseq;	     /* number of training sequences          (mandatory) */
  float  eff_nseq;   /* effective number of seqs (<= nseq)    (mandatory) */
  char  *ctime;	     /* creation date                         (mandatory) */
  int   *map;	     /* map of alignment cols onto model 1..M (p7_MAP)    */ /* Array; 0=0 */
  int    checksum;   /* checksum of training sequences        (mandatory) */

  float  evparam[p7_NEVPARAM]; 	/* parameters for determining E-values                        */
  float  cutoff[p7_NCUTOFFS]; 	/* per-seq/per-domain gathering, trusted, noise cutoffs       */

  off_t  offset;              /* HMM record offset on disk */
  int    flags;               /* status flags */
  const ESL_ALPHABET *abc;    /* ptr to alphabet info (hmm->abc->K is alphabet size) */
} P7_HMM;



/*****************************************************************
 * 2. P7_PROFILE: a scoring profile, and its implicit model.
 *****************************************************************/

/* Indices for special state types in the length model, gm->xsc[x][]
 */
enum p7p_xstates_e { 
  p7P_E = 0,
  p7P_N = 1,
  p7P_J = 2,
  p7P_C = 3
};
#define p7P_NXSTATES 4

/* Indices for transitions from the length modeling scores gm->xsc[][x]
 */
enum p7p_xtransitions_e {
  p7P_LOOP = 0,
  p7P_MOVE = 1
};
#define p7P_NXTRANS 2

/* Indices for transition scores gm->tsc[k][] */
/* order is optimized for dynamic programming */
enum p7p_tsc_e {
  p7P_MM = 0, 
  p7P_IM = 1, 
  p7P_DM = 2, 
  p7P_BM = 3, 
  p7P_MD = 4, 
  p7P_DD = 5, 
  p7P_MI = 6, 
  p7P_II = 7, 
};
#define p7P_NTRANS 8

/* Indices for residue emission score vectors
 */
enum p7p_rsc_e {
  p7P_MSC = 0, 
  p7P_ISC = 1
};
#define p7P_NR 2

/* Accessing transition, emission scores */
/* _BM is specially stored off-by-one: [k-1][p7P_BM] is score for entering at Mk */
#define p7P_TSC(gm, k, s) ((gm)->tsc[(k) * p7P_NTRANS + (s)])
#define p7P_MSC(gm, k, x) ((gm)->rsc[x][(k) * p7P_NR + p7P_MSC])
#define p7P_ISC(gm, k, x) ((gm)->rsc[x][(k) * p7P_NR + p7P_ISC])


typedef struct p7_profile_s {
  float  *tsc;          /* transitions  [0.1..M-1][0..p7P_NTRANS-1], hand-indexed  */
  float **rsc;          /* emissions [0..Kp-1][0.1..M][p7P_NR], hand-indexed       */
  float   xsc[p7P_NXSTATES][p7P_NXTRANS]; /* special transitions [NECJ][LOOP,MOVE] */

  int     mode;        	/* configured algorithm mode (e.g. p7_LOCAL)               */ 
  int     L;		/* current configured target seq length                    */
  int     allocM;	/* max # of nodes allocated in this structure              */
  int     M;		/* number of nodes in the model                            */
  float   nj;		/* expected # of uses of J; precalculated from loop config */

  /* Info, most of which is a copy from parent HMM:                                       */
  char  *name;			/* unique name of model                                   */
  char  *acc;			/* unique accession of model, or NULL                     */
  char  *desc;                  /* brief (1-line) description of model, or NULL           */
  char  *rf;                    /* reference line from alignment 1..M; *rf=0 means unused */
  char  *cs;                    /* consensus structure line      1..M, *cs=0 means unused */
  char  *consensus;		/* consensus residues to display in alignments, 1..M      */
  float  evparam[p7_NEVPARAM]; 	/* parameters for determining E-values                    */
  float  cutoff[p7_NCUTOFFS]; 	/* per-seq/per-domain gathering, trusted, noise cutoffs   */
  const ESL_ALPHABET *abc;	/* copy of pointer to appropriate alphabet                */
} P7_PROFILE;




/*****************************************************************
 * 3. P7_BG: a null (background) model.
 *****************************************************************/

typedef struct p7_bg_s {
  ESL_ALPHABET *abc;		/* reference to alphabet in use       */
  
  float  p1;			/* null model's self-loop probability */
  float *f;			/* residue frequencies [0..K-1] */
} P7_BG;

/*****************************************************************
 * 4. P7_TRACE:  a traceback (alignment of seq to profile).
 *****************************************************************/

/* Traceback structure for alignment of a model to a sequence.
 *
 * A traceback only makes sense in a triplet (tr, gm, dsq), for a
 * given profile or HMM (with nodes 1..M) and a given digital sequence
 * (with positions 1..L).
 * 
 * A traceback may be relative to a profile (usually) or to a core
 * model (as a special case in model construction; see build.c). You
 * can tell the difference by looking at the first statetype,
 * tr->st[0]; if it's a p7T_S, it's for a profile, and if it's p7T_B,
 * it's for a core model.
 * 
 * A profile's N,C,J states emit on transition, not on state, so a
 * path of N emits 0 residues, NN emits 1 residue, NNN emits 2
 * residues, and so on. By convention, the trace always associates an
 * emission-on-transition with the trailing (destination) state, so
 * the first N, C, or J is stored in a trace as a nonemitter (i=0).
 *
 * A i coords in a traceback are usually 1..L with respect to an
 * unaligned digital target sequence, but in the special case of
 * traces faked from existing MSAs (as in hmmbuild), the coords may
 * be 1..alen relative to an MSA's columns.
 */

/* State types */
enum p7t_statetype_e {
  p7T_BOGUS =  0,
  p7T_M     =  1,
  p7T_D     =  2,
  p7T_I     =  3,
  p7T_S     =  4,
  p7T_N     =  5,
  p7T_B     =  6, 
  p7T_E     =  7,
  p7T_C     =  8, 
  p7T_T     =  9, 
  p7T_J     = 10,
  p7T_X     = 11, 	/* missing data: used esp. for local entry/exits */
};
#define p7T_NSTATETYPES 12

typedef struct p7_trace_s {
  int   N;		/* length of traceback                       */
  int   nalloc;	        /* allocated length of traceback             */
  char  *st;		/* state type code                   [0..N-1]*/
  int   *k;		/* node index; 1..M if M,D,I; else 0 [0..N-1]*/
  int   *i;		/* position in dsq, 1..L; 0 if none  [0..N-1]*/

  /* The following section is data generated by "indexing" a trace's domains */
  int   ndom;		/* number of domains in trace (= # of B or E states) */
  int  *tfrom,   *tto;	/* locations of B/E states in trace (0..tr->N-1)     */
  int  *sqfrom,  *sqto;	/* first/last M-emitted residue on sequence (1..L)   */
  int  *hmmfrom, *hmmto;/* first/last M state on model (1..M)                */
  int   ndomalloc;	/* current allocated size of these stacks            */

} P7_TRACE;


/*****************************************************************
 * 5. P7_HMMFILE:  an HMM save file or database, open for reading.
 *****************************************************************/

typedef struct p7_hmmfile_s {
  FILE         *f;		 /* pointer to stream for reading                */
  char         *fname;	         /* name of the HMM file; [STDIN] if -           */

  int (*parser)(struct p7_hmmfile_s *, ESL_ALPHABET **, P7_HMM **);  /* parsing function */

  int           do_gzip;	/* TRUE if f is "gzip -dc |" (will pclose(f))    */ 
  int           do_stdin;       /* TRUE if f is stdin (won't close f)            */
  ESL_SSI      *ssi;		/* open SSI index file; or NULL if none.         */
} P7_HMMFILE;




/*****************************************************************
 * 6. P7_GMX: a "generic" dynamic programming matrix
 *****************************************************************/

enum p7g_scells_e {
  p7G_M = 0,
  p7G_I = 1,
  p7G_D = 2,
};
#define p7G_NSCELLS 3

enum p7g_xcells_e {
  p7G_E  = 0,
  p7G_N  = 1,
  p7G_J  = 2,
  p7G_B  = 3,
  p7G_C  = 4
};
#define p7G_NXCELLS 5


typedef struct p7_gmx_s {
  int  M;		/* actual model dimension (model 1..M)    */
  int  L;		/* actual sequence dimension (seq 1..L)   */
  
  int      allocR;      /* current allocated # of rows : L+1 <= validR <= allocR                */
  int      validR;	/* # of rows actually pointing at DP memory                             */
  int      allocW;	/* current set row width :  M+1 <= allocW                               */
  uint64_t ncells;	/* total # of allocated cells in 2D matrix : ncells >= (validR)(allocW) */

  float **dp;           /* logically [0.1..L][0.1..M][0..p7G_NSCELLS-1]; indexed [i][k*p7G_NSCELLS+s] */
  float  *xmx;          /* logically [0.1..L][0..p7G_NXCELLS-1]; indexed [i*p7G_NXCELLS+s]            */

  float  *dp_mem;
} P7_GMX;



/*****************************************************************
 * 7. P7_DPRIOR: mixture Dirichlet prior for profile HMMs
 *****************************************************************/

typedef struct p7_dprior_s {
  ESL_MIXDCHLET *tm;		/*  match transitions */
  ESL_MIXDCHLET *ti;		/* insert transitions */
  ESL_MIXDCHLET *td;		/* delete transitions */
  ESL_MIXDCHLET *em;		/*  match emissions   */
  ESL_MIXDCHLET *ei;		/* insert emissions   */
} P7_DPRIOR;



/*****************************************************************
 * 8. P7_SPENSEMBLE: segment pair ensembles for domain locations
 *****************************************************************/

/* struct p7_spcoord_s:
 *    a coord quad defining a segment pair. 
 */
struct p7_spcoord_s { 
  int idx; 	/* backreference index: which trace a seg came from, or which cluster a domain came from */
  int i, j;	/* start,end in a target sequence (1..L)  */
  int k, m;     /* start,end in a query model (1..M)      */
};

/* Structure: P7_SPENSEMBLE
 *
 * Collection and clustering of an ensemble of sampled segment pairs,
 * in order to define domain locations using their posterior
 * probability distribution (as opposed to Viterbi MAP tracebacks).
 */
typedef struct p7_spensemble_s {
  /* Section 1: a collected ensemble of segment pairs                                       */
  int                  nsamples;    /* number of sampled traces                             */
  struct p7_spcoord_s *sp;	    /* array of sampled seg pairs; [0..n-1]                 */
  int                  nalloc;	    /* allocated size of <sp>                               */
  int                  n;	    /* number of seg pairs in <sp>                          */

  /* Section 2: then the ensemble is clustered by single-linkage clustering                 */
  int *workspace;                   /* temp space for Easel SLC algorithm: 2*n              */
  int *assignment;                  /* each seg pair's cluster index: [0..n-1] = (0..nc-1)  */
  int  nc;	                    /* number of different clusters                         */

  /* Section 3: then endpoint distribution is examined within each large cluster            */
  int *epc;	                    /* array counting frequency of each endpoint            */
  int  epc_alloc;	            /* allocated width of <epc>                             */
  
  /* Section 4: finally each large cluster is resolved into domain coords                   */
  struct p7_spcoord_s *sigc;	    /* array of coords for each domain, [0..nsigc-1]        */
  float               *prob;	    /* posterior probability of each domain [0..nsigc-1]    */
  int                  nsigc;	    /* number of "significant" clusters, domains            */
  int                  nsigc_alloc; /* current allocated max for nsigc                      */
} P7_SPENSEMBLE;


/*****************************************************************
 * 9. P7_ALIDISPLAY: an alignment formatted for printing
 *****************************************************************/

/* Structure: P7_ALIDISPLAY
 * 
 * Alignment of a sequence domain to an HMM, formatted for printing.
 * 
 * For an alignment of L residues and names C chars long, requires
 * 5L + 2C + 30 bytes; for typical case of L=100,C=10, that's
 * 0.5 Kb.
 */
typedef struct p7_alidisplay_s {
  char *rfline;                 /* reference coord info                 */
  char *csline;                 /* consensus structure info             */
  char *model;                  /* aligned query consensus sequence     */
  char *mline;                  /* "identities", conservation +'s, etc. */
  char *aseq;                   /* aligned target sequence              */
  int   N;			/* length of strings                    */
  char *hmmname;		/* name of HMM                          */
  char *hmmacc;			/* accession of HMM                     */
  char *hmmdesc;		/* description of HMM                   */
  int   hmmfrom;		/* start position on HMM (1..M, or -1)  */
  int   hmmto;			/* end position on HMM (1..M, or -1)    */
  char *sqname;			/* name of target sequence              */
  char *sqacc;			/* accession of target sequence         */
  char *sqdesc;			/* description of target sequence       */
  long  sqfrom;			/* start position on sequence (1..L)    */
  long  sqto;		        /* end position on sequence   (1..L)    */

  char *mem;			/* memory used for the char data above  */
} P7_ALIDISPLAY;


/*****************************************************************
 * 10. P7_DOMAINDEF: reusably managing workflow in defining domains
 *****************************************************************/

typedef struct p7_dom_s { 
  int            ienv, jenv;
  int            iali, jali;
  float          envsc;  	/* Forward score in envelope ienv..jenv; without null2 correction */
  float          domcorrection;	/* null2 correction to add null score when calculating a per-domain score */
  float          bitscore;	/* overall score in bits, null corrected, if this were the only domain in seq */
  double         pvalue;	/* P-value of the bitscore */
  P7_ALIDISPLAY *ad; 
} P7_DOMAIN;

/* Structure: P7_DOMAINDEF
 * 
 * This is a container for all the necessary information for domain
 * definition procedures in <p7_domaindef.c>, including a bunch of
 * heuristic thresholds. The structure is reusable to minimize the
 * number of allocation/free cycles that need to be done when
 * processing a large number of sequences. You create the structure
 * with <p7_domaindef_Create()>; after you're done with defining
 * domains on a sequence, you call <p7_domaindef_Reuse()> before using
 * it on the next sequence; and when you're completely done, you free
 * it with <p7_domaindef_Destroy()>. All memory management is handled
 * internally; you don't need to reallocate anything yourself.
 */
typedef struct p7_domaindef_s {
  /* for posteriors of being in a domain, B, E */
  float *mocc;			/* mocc[i=1..L] = prob that i is emitted by core model (is in a domain)       */
  float *btot; 			/* btot[i=1..L] = cumulative expected times that domain starts at or before i */
  float *etot;			/* etot[i=1..L] = cumulative expected times that domain ends at or before i   */
  int    L;
  int    Lalloc;

  /* the ad hoc null2 model: 1..L nat scores for each residue, log f'(x_i) / f(x_i) */
  float *n2sc;

  /* rng and reusable memory for stochastic tracebacks */
  ESL_RANDOMNESS *r;
  P7_SPENSEMBLE  *sp;
  P7_TRACE       *tr;		/* reusable space for a trace of a domain                  */
  P7_TRACE       *gtr;		/* reusable space for a traceback of the entire target seq */

  /* Heuristic thresholds that control the region definition process */
  /* "rt" = "region threshold", for lack of better term  */
  float  rt1;   	/* controls when regions are called. mocc[i] post prob >= dt1 : triggers a region around i */
  float  rt2;		/* controls extent of regions. regions extended until mocc[i]-{b,e}occ[i] < dt2            */
  float  rt3;		/* controls when regions are flagged for split: if expected # of E preceding B is >= dt3   */
  
  /* Heuristic thresholds that control the stochastic traceback/clustering process */
  int    nsamples;	/* collect ensemble of this many stochastic traces */
  float  min_overlap;	/* 0.8 means >= 80% overlap of (smaller/larger) segment to link, both in seq and hmm            */
  int    of_smaller;	/* see above; TRUE means overlap denom is calc'ed wrt smaller segment; FALSE means larger       */
  int    max_diagdiff;	/* 4 means either start or endpoints of two segments must be within <=4 diagonals of each other */
  float  min_posterior;	/* 0.25 means a cluster must have >= 25% posterior prob in the sample to be reported            */
  float  min_endpointp;	/* 0.02 means choose widest endpoint with post prob of at least 2%                              */

  /* storage of the results; domain locations, scores, alignments          */
  P7_DOMAIN *dcl;
  int        ndom;	 /* number of domains defined, in the end.         */
  int        nalloc;     /* number of domain structures allocated in <dcl> */

  /* Additional results storage */
  float  nexpected;     /* posterior expected number of domains in the sequence (from posterior arrays) */
  int    nregions;	/* number of regions evaluated */
  int    nclustered;	/* number of regions evaluated by clustering ensemble of tracebacks */
  int    noverlaps;	/* number of envelopes defined in ensemble clustering that overlap w/ prev envelope */
  int    nenvelopes;	/* number of envelopes handed over for domain definition, null2, alignment, and scoring. */

} P7_DOMAINDEF;


/*****************************************************************
 * 11. P7_TOPHITS: ranking lists of top-scoring hits
 *****************************************************************/

/* Structure: P7_HIT
 * 
 * Info about a high-scoring database hit, kept so we can output a
 * sorted list of high hits at the end.
 *
 * sqfrom and sqto are the coordinates that will be shown
 * in the results, not coords in arrays... therefore, reverse
 * complements have sqfrom > sqto
 */
typedef struct p7_hit_s {
  char   *name;			/* name of the target               */
  char   *acc;			/* accession of the target          */
  char   *desc;			/* description of the target        */
  double sortkey;		/* number to sort by; big is better    */

  float  score;			/* bit score of the hit                        */
  float  pre_score;		/* bit score before null2 correction           */
  float  sum_score;		/* bit score reconstructed from sum of domains */

  double pvalue;		/* P-value of the score               */
  double pre_pvalue;		/* P-value of the pre_score           */
  double sum_pvalue;		/* P-value of the sum_score           */

  int    ndom;		/* total # of domains in this seq   */
  float  nexpected;     /* posterior expected number of domains in the sequence (from posterior arrays) */
  int    nregions;	/* number of regions evaluated */
  int    nclustered;	/* number of regions evaluated by clustering ensemble of tracebacks */
  int    noverlaps;	/* number of envelopes defined in ensemble clustering that overlap w/ prev envelope */
  int    nenvelopes;	/* number of envelopes handed over for domain definition, null2, alignment, and scoring. */

  P7_DOMAIN *dcl;		/* domain coordinate list and alignment display */
} P7_HIT;


/* Structure: P7_TOPHITS
 * 

 * merging when we prepare to output results. "hit" list is NULL and
 * unavailable until after we do a sort.  
 */
typedef struct p7_tophits_s {
  P7_HIT **hit;         /* sorted pointer array                     */
  P7_HIT  *unsrt;	/* unsorted data storage                    */
  int      Nalloc;	/* current allocation size                  */
  int      N;		/* number of hits in list now               */
  int      is_sorted;	/* TRUE when h->hit valid for all N hits    */
} P7_TOPHITS;




/*****************************************************************
 * 12. The optimized implementation.
 *****************************************************************/
#if   defined (p7_IMPL_SSE)
#include "impl_sse.h"
#elif defined (p7_IMPL_VMX)
#include "impl_vmx.h"
#else
#include "impl_dummy.h"
#endif

/*****************************************************************
 * 13. Other routines in HMMER's exposed API.
 *****************************************************************/

/* build.c */
extern int p7_Handmodelmaker(ESL_MSA *msa,                P7_HMM **ret_hmm, P7_TRACE ***ret_tr);
extern int p7_Fastmodelmaker(ESL_MSA *msa, float symfrac, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);

/* dp_generic.c */
extern int p7_GViterbi     (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);
extern int p7_GForward     (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);
extern int p7_GBackward    (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);
extern int p7_GHybrid      (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *opt_fwdscore, float *opt_hybscore);
extern int p7_GMSP         (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);
extern int p7_GTrace       (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr);
extern int p7_StochasticTrace(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr);
extern int p7_GViterbiScore(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,       P7_GMX *gx, float *ret_sc);

/* dp_optaccuracy.c */
extern int p7_PosteriorDecoding(int L, const P7_PROFILE *gm, const P7_GMX *fwd,      P7_GMX *bck, P7_GMX *pp);
extern int p7_OptimalAccuracyDP(int L, const P7_PROFILE *gm, const P7_GMX *pp,       P7_GMX *gx,  float *ret_e);
extern int p7_OATrace          (int L, const P7_PROFILE *gm, const P7_GMX *pp, const P7_GMX *gx,  P7_TRACE *tr);


/* emit.c */
extern int p7_CoreEmit   (ESL_RANDOMNESS *r, const P7_HMM *hmm,                                        ESL_SQ *sq, P7_TRACE *tr);
extern int p7_ProfileEmit(ESL_RANDOMNESS *r, const P7_HMM *hmm, const P7_PROFILE *gm, const P7_BG *bg, ESL_SQ *sq, P7_TRACE *tr);

/* errors.c */
extern void p7_Die (char *format, ...);
extern void p7_Fail(char *format, ...);

/* evalues.c */
extern int p7_Lambda(P7_HMM *hmm, P7_BG *bg, double *ret_lambda);
extern int p7_Mu (ESL_RANDOMNESS *r, P7_PROFILE *gm, P7_BG *bg, int L, int N, double lambda,               double *ret_mu);
extern int p7_Tau(ESL_RANDOMNESS *r, P7_PROFILE *gm, P7_BG *bg, int L, int N, double lambda, double tailp, double *ret_tau);

/* eweight.c */
extern int  p7_EntropyWeight(const P7_HMM *hmm, const P7_BG *bg, const P7_DPRIOR *pri, double infotarget, double *ret_Neff);

/* heatmap.c (evolving now, intend to move this to Easel in the future) */
extern double dmx_upper_max(ESL_DMATRIX *D);
extern double dmx_upper_min(ESL_DMATRIX *D);
extern double dmx_upper_element_sum(ESL_DMATRIX *D);
extern double dmx_upper_norm(ESL_DMATRIX *D);
extern int    dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max);

/* island.c */
extern int   p7_island_Viterbi(ESL_DSQ *dsq, int L, P7_PROFILE *gm, P7_GMX *mx, ESL_HISTOGRAM *h);

/* hmmer.c */
extern void  p7_banner(FILE *fp, char *progname, char *banner);
extern float p7_SILO2Lod(int silo);
extern int   p7_AminoFrequencies(float *f);

/* logsum.c */
extern int   p7_FLogsumInit(void);
extern float p7_FLogsum(float a, float b);
extern int   p7_ILogsumInit(void);
extern int   p7_ILogsum(int s1, int s2);


/* modelconfig.c */
extern int p7_ProfileConfig(const P7_HMM *hmm, const P7_BG *bg, P7_PROFILE *gm, int L, int mode);
extern int p7_ReconfigLength  (P7_PROFILE *gm, int L);
extern int p7_ReconfigMultihit(P7_PROFILE *gm, int L);
extern int p7_ReconfigUnihit  (P7_PROFILE *gm, int L);

/* modelstats.c */
extern double p7_MeanMatchInfo           (const P7_HMM *hmm, const P7_BG *bg);
extern double p7_MeanMatchEntropy        (const P7_HMM *hmm);
extern double p7_MeanMatchRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg);
extern double p7_MeanForwardScore        (const P7_HMM *hmm, const P7_BG *bg);
extern int    p7_MeanPositionRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg, double *ret_entropy);
extern int    p7_hmm_CompositionKLDist(P7_HMM *hmm, P7_BG *bg, float *ret_KL, float **opt_avp);

/* mpisupport.c */
#ifdef HAVE_MPI
extern int p7_hmm_MPISend(P7_HMM *hmm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int p7_hmm_MPIPackSize(P7_HMM *hmm, MPI_Comm comm, int *ret_n);
extern int p7_hmm_MPIPack(P7_HMM *hmm, char *buf, int n, int *position, MPI_Comm comm);
extern int p7_hmm_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ESL_ALPHABET **abc, P7_HMM **ret_hmm);
extern int p7_hmm_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, P7_HMM **ret_hmm);

extern int p7_profile_MPISend(P7_PROFILE *gm, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int p7_profile_MPIRecv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc, const P7_BG *bg,
			      char **buf, int *nalloc,  P7_PROFILE **ret_gm);
#endif /*HAVE_MPI*/


/* null2.c */
extern int p7_null2_ByExpectation(P7_DOMAINDEF *ddef, const P7_PROFILE *gm, const ESL_DSQ *dsq, int ienv, int jenv,
				  const P7_GMX *pp, P7_GMX *wrk, float *opt_null2);
extern int p7_null2_BySampling   (P7_DOMAINDEF *ddef, const P7_PROFILE *gm, const ESL_DSQ *dsq, int ireg, int jreg, P7_GMX *gx);


/* p7_alidisplay.c */
extern P7_ALIDISPLAY *p7_alidisplay_Create(const P7_TRACE *tr, int which, const P7_PROFILE *gm, const ESL_SQ *sq);
extern void           p7_alidisplay_Destroy(P7_ALIDISPLAY *ad);
extern int            p7_alidisplay_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth);

/* p7_bg.c */
extern P7_BG *p7_bg_Create(const ESL_ALPHABET *abc);
extern P7_BG *p7_bg_CreateUniform(const ESL_ALPHABET *abc);
extern int    p7_bg_Dump(FILE *ofp, P7_BG *bg);
extern void   p7_bg_Destroy(P7_BG *bg);
extern int    p7_bg_SetLength(P7_BG *bg, int L);
extern int    p7_bg_NullOne(const P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc);

/* p7_domaindef.c */
extern P7_DOMAINDEF *p7_domaindef_Create (ESL_RANDOMNESS *r);
extern int           p7_domaindef_Fetch  (P7_DOMAINDEF *ddef, int which, int *opt_i, int *opt_j, float *opt_sc, P7_ALIDISPLAY **opt_ad);
extern int           p7_domaindef_Reuse  (P7_DOMAINDEF *ddef);
extern int           p7_domaindef_DumpPosteriors(FILE *ofp, P7_DOMAINDEF *ddef);
extern void          p7_domaindef_Destroy(P7_DOMAINDEF *ddef);

extern int p7_domaindef_ByViterbi            (P7_PROFILE *gm, const ESL_SQ *sq, P7_GMX *gx1, P7_GMX *gx2, P7_DOMAINDEF *ddef);
extern int p7_domaindef_ByPosteriorHeuristics(P7_PROFILE *gm, const ESL_SQ *sq, P7_GMX *fwd, P7_GMX *bck, P7_DOMAINDEF *ddef);


/* p7_gmx.c */
extern P7_GMX *p7_gmx_Create(int allocM, int allocL);
extern int     p7_gmx_GrowTo(P7_GMX *gx, int allocM, int allocL);
extern void    p7_gmx_Destroy(P7_GMX *gx);
extern int     p7_gmx_Dump(FILE *fp, P7_GMX *gx);
extern int     p7_gmx_DumpWindow(FILE *fp, P7_GMX *gx, int istart, int iend, int kstart, int kend, int show_specials);


/* p7_hmm.c */
/*      1. The P7_HMM object: allocation, initialization, destruction. */
extern P7_HMM *p7_hmm_Create(int M, const ESL_ALPHABET *abc);
extern P7_HMM *p7_hmm_CreateShell(void);
extern int     p7_hmm_CreateBody(P7_HMM *hmm, int M, const ESL_ALPHABET *abc);
extern void    p7_hmm_Destroy(P7_HMM *hmm);
extern int     p7_hmm_CopyParameters(const P7_HMM *src, P7_HMM *dest);
extern P7_HMM *p7_hmm_Clone(const P7_HMM *hmm);
extern int     p7_hmm_Scale(P7_HMM *hmm, double scale);
extern int     p7_hmm_Zero(P7_HMM *hmm);
extern char   *p7_hmm_DescribeStatetype(char st);
/*      2. Convenience routines for setting fields in an HMM. */
extern int     p7_hmm_SetName(P7_HMM *hmm, char *name);
extern int     p7_hmm_SetAccession(P7_HMM *hmm, char *acc);
extern int     p7_hmm_SetDescription(P7_HMM *hmm, char *desc);
extern int     p7_hmm_AppendComlog(P7_HMM *hmm, int argc, char **argv);
extern int     p7_hmm_SetCtime(P7_HMM *hmm);
/*      3. Renormalization and rescaling counts in core HMMs. */
extern int     p7_hmm_Rescale(P7_HMM *hmm, float scale);
extern int     p7_hmm_Renormalize(P7_HMM *hmm);
/*      4. Debugging and development code. */
extern int     p7_hmm_Dump(FILE *fp, P7_HMM *hmm);
extern int     p7_hmm_Sample        (ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, P7_HMM **ret_hmm);
extern int     p7_hmm_SampleUngapped(ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, P7_HMM **ret_hmm);
extern int     p7_hmm_SampleEnumerable(ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, P7_HMM **ret_hmm);
extern int     p7_hmm_SampleUniform (ESL_RANDOMNESS *r, int M, ESL_ALPHABET *abc, 
				     float tmi, float tii, float tmd, float tdd,  P7_HMM **ret_hmm);
extern int     p7_hmm_Compare(P7_HMM *h1, P7_HMM *h2, float tol);
extern int     p7_hmm_Validate(P7_HMM *hmm, char *errbuf, float tol);
/*      5. Other routines in the API */
extern int     p7_hmm_CalculateOccupancy(const P7_HMM *hmm, float *occ);



/* p7_hmmfile.c */
extern int  p7_hmmfile_Open(char *filename, char *env, P7_HMMFILE **ret_hfp);
extern void p7_hmmfile_Close(P7_HMMFILE *hfp);
extern int  p7_hmmfile_Write(FILE *fp, P7_HMM *hmm);
extern int  p7_hmmfile_Read(P7_HMMFILE *hfp, ESL_ALPHABET **ret_abc,  P7_HMM **ret_hmm);
extern int  p7_hmmfile_PositionByKey(P7_HMMFILE *hfp, const char *key);

/* p7_prior.c */
extern P7_DPRIOR *p7_dprior_CreateAmino(void);
extern P7_DPRIOR *p7_dprior_CreateNucleic(void);
extern P7_DPRIOR *p7_dprior_CreateLaplace(ESL_ALPHABET *abc);
extern void       p7_dprior_Destroy(P7_DPRIOR *pri);
extern int        p7_ParameterEstimation(P7_HMM *hmm, const P7_DPRIOR *pri);

/* p7_profile.c */
extern P7_PROFILE *p7_profile_Create(int M, const ESL_ALPHABET *abc);
extern P7_PROFILE *p7_profile_Clone(const P7_PROFILE *gm);
extern int         p7_profile_Copy(const P7_PROFILE *src, P7_PROFILE *dst);
extern int         p7_profile_SetNullEmissions(P7_PROFILE *gm);
extern int         p7_profile_Reuse(P7_PROFILE *gm);
extern void        p7_profile_Destroy(P7_PROFILE *gm);
extern int         p7_profile_IsLocal(const P7_PROFILE *gm);
extern int         p7_profile_IsMultihit(const P7_PROFILE *gm);
extern int         p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1, 
				   char st2, int k2, float *ret_tsc);
extern int         p7_profile_Validate(const P7_PROFILE *gm, char *errbuf, float tol);
extern int         p7_profile_Compare(P7_PROFILE *gm1, P7_PROFILE *gm2, float tol);

/* p7_spensemble.c */
P7_SPENSEMBLE *p7_spensemble_Create(int init_n, int init_epc, int init_sigc);
extern int     p7_spensemble_Reuse(P7_SPENSEMBLE *sp);
extern int     p7_spensemble_Add(P7_SPENSEMBLE *sp, int sampleidx, int i, int j, int k, int m);
extern int     p7_spensemble_Cluster(P7_SPENSEMBLE *sp, 
				     float min_overlap, int of_smaller, int max_diagdiff, 
				     float min_posterior, float min_endpointp,
				     int *ret_nclusters);
extern int     p7_spensemble_GetClusterCoords(P7_SPENSEMBLE *sp, int which,
					      int *ret_i, int *ret_j, int *ret_k, int *ret_m, float *ret_p);
extern void    p7_spensemble_Destroy(P7_SPENSEMBLE *sp);

/* p7_tophits.c */
extern P7_TOPHITS *p7_tophits_Create(void);
extern int         p7_tophits_Grow(P7_TOPHITS *h);
extern int         p7_tophits_CreateNextHit(P7_TOPHITS *h, P7_HIT **ret_hit);
extern int         p7_tophits_Add(P7_TOPHITS *h,
				  char *name, char *acc, char *desc, 
				  double sortkey, 
				  float score,    double pvalue, 
				  float mothersc, double motherp,
				  int sqfrom, int sqto, int sqlen,
				  int hmmfrom, int hmmto, int hmmlen, 
				  int domidx, int ndom,
				  P7_ALIDISPLAY *ali);
extern int         p7_tophits_Sort(P7_TOPHITS *h);
extern int         p7_tophits_Merge(P7_TOPHITS *h1, P7_TOPHITS *h2);
extern int         p7_tophits_GetMaxNameLength(P7_TOPHITS *h);
extern void        p7_tophits_Destroy(P7_TOPHITS *h);


/* p7_trace.c */
extern P7_TRACE *p7_trace_Create(void);
extern int  p7_trace_Reuse(P7_TRACE *tr);
extern int  p7_trace_Grow(P7_TRACE *tr);
extern int  p7_trace_GrowIndex(P7_TRACE *tr);
extern int  p7_trace_GrowTo(P7_TRACE *tr, int N);
extern int  p7_trace_GrowIndexTo(P7_TRACE *tr, int ndom);
extern void p7_trace_Destroy(P7_TRACE *tr);
extern void p7_trace_DestroyArray(P7_TRACE **tr, int N);
extern int  p7_trace_Validate(P7_TRACE *tr, ESL_ALPHABET *abc, ESL_DSQ *sq, char *errbuf);
extern int  p7_trace_Dump(FILE *fp, P7_TRACE *tr, P7_PROFILE *gm, ESL_DSQ *dsq);

extern int  p7_trace_Append(P7_TRACE *tr, char st, int k, int i);
extern int  p7_trace_Reverse(P7_TRACE *tr);
extern int  p7_trace_Index(P7_TRACE *tr);
extern int  p7_trace_Count(P7_HMM *hmm, ESL_DSQ *dsq, float wt, P7_TRACE *tr);
extern int  p7_trace_Score(P7_TRACE *tr, ESL_DSQ *dsq, P7_PROFILE *gm, float *ret_sc);
extern int  p7_trace_GetDomainCount(P7_TRACE *tr, int *ret_ndom);
extern int  p7_trace_StateUseCounts(const P7_TRACE *tr, int *counts);
extern int  p7_trace_GetDomainCoords(P7_TRACE *tr, int which, int *ret_i1, int *ret_i2,
				     int *ret_k1, int *ret_k2);




/* seqmodel.c */
extern int p7_Seqmodel(ESL_ALPHABET *abc, ESL_DSQ *dsq, int M, char *name,
		       ESL_DMATRIX *P, float *f, double popen, double pextend,
		       P7_HMM **ret_hmm);

#endif /*P7_HMMERH_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
