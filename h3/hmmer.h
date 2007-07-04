/* The all-encompassing include file for HMMER.
 *
 *    1. P7_HMM:     a core model.
 *    2. P7_PROFILE: a scoring profile, and its implicit model.
 *    3. P7_BG:      a null (background) model.
 *    4. P7_TRACE:   a traceback path (alignment of seq to profile).
 *    5. P7_HMMFILE: an HMM save file or database, open for reading.
 *    6. P7_GMX:     a "generic" dynamic programming matrix
 *    7. P7_DPRIOR:  mixture Dirichlet prior for profile HMMs
 *    8. Other routines in HMMER's exposed API.
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
#include "esl_sqio.h"		/* ESL_SQ                */
#include "esl_histogram.h"      /* ESL_HISTOGRAM         */
#include "esl_dirichlet.h"	/* ESL_MIXDCHLET         */

/*****************************************************************
 * 1. P7_HMM: a core model.
 *****************************************************************/

/* Some notes:
 *   0. The model might be either in counts or probability form.
 *   1. t[0] is special: t[0][TMM,TMI,TMD] are the begin->M_1,I_0,D_1 entry probabilities,
 *      t[0][TIM,TII] are the I_0 transitions, and delete state 0 doesn't
 *      exist. Therefore D[0] transitions and mat[0] emissions are unused.
 *      To simplify some normalization code, we adopt a convention that these are set
 *      to valid probability distributions: 1.0 for t[0][TDM] and mat[0][0],
 *      and 0 for the rest.
 *   2. t[M] is also special: TMD and TDD are 0 because there is no next delete state;
 *      TDM is therefore 1.0 by definition; and TMM, TDM are interpreted as the
 *      M->E and D->E end transitions. 
 */
typedef struct p7_hmm_s {
  /*::cexcerpt::plan7_core::begin::*/
  int     M;                    /* length of the model (# nodes)          */
  float **t;                    /* transition prob's. t[(0),1..M][0..6]   */
  float **mat;                  /* match emissions.  mat[1..M][0..K-1]    */ 
  float **ins;                  /* insert emissions. ins[1..M][0..K-1]    */
  /*::cexcerpt::plan7_core::end::*/

  /* Annotation. Everything but <name> is optional. Flags are set when
   * optional values are set. All the char *'s are proper nul-terminated
   * strings, not just arrays. (hmm->map is an int array).
   */
  char  *name;                  /* name of the model                     (mandatory) */ /* String, \0-terminated */
  char  *acc;			/* accession number of model (Pfam)      (p7_ACC)    */ /* String, \0-terminated */
  char  *desc;                  /* brief (1-line) description of model   (p7_DESC)   */ /* String, \0-terminated */
  char  *rf;                    /* reference line from alignment 1..M    (p7_RF)     */ /* String; 0=' ', M+1='\0' */
  char  *cs;                    /* consensus structure line      1..M    (p7_CS)     */ /* String; 0=' ', M+1='\0' */
  char  *ca;			/* consensus accessibility line  1..M    (p7_CA)     */ /* String; 0=' ', M+1='\0' */
  char  *comlog;		/* command line(s) that built model      (mandatory) */ /* String, \0-terminated */
  int    nseq;			/* number of training sequences          (mandatory) */
  float  eff_nseq;		/* effective number of seqs (<= nseq)    (mandatory) */
  char  *ctime;			/* creation date                         (mandatory) */
  int   *map;			/* map of alignment cols onto model 1..M (p7_MAP)    */ /* Array; 0=0 */
  int    checksum;              /* checksum of training sequences        (mandatory) */

  /* Pfam-specific score cutoffs.
   * 
   * ga1, ga2 are valid if p7_GA is set in flags.
   * tc1, tc2 are valid if p7_TC is set in flags.
   * nc1, nc2 are valid if p7_NC is set in flags.
   */
  float  ga1, ga2;	/* per-seq/per-domain gathering thresholds (bits) (p7_GA) */
  float  tc1, tc2;	/* per-seq/per-domain trusted cutoff (bits)       (p7_TC) */
  float  nc1, nc2;	/* per-seq/per-domain noise cutoff (bits)         (p7_NC) */

  /* Things we keep references to.  */
  const ESL_ALPHABET        *abc; /* ptr to alphabet info (hmm->abc->K is alphabet size) */
  const struct p7_profile_s *gm;  /* generic search profile (incomplete type: P7_PROFILE declared below) */
  const struct p7_bg_s      *bg;  /* null background model  (incomplete type: P7_BG declared below)      */

  /* When an HMM is input from disk, remember its disk position. */
  off_t  offset;		/* HMM record offset on disk */

  int flags;
} P7_HMM;

/* Flag codes for hmm->flags.
 * Flags marked with ! may not be changed nor used for other meanings;
 * such flags were stored in old HMM files, and we must preserve their
 * meaning to preserve reverse compatibility.
 */
#define p7_HASBITS (1<<0)    /* obsolete (was: model has log-odds scores)       !*/
#define p7_DESC    (1<<1)    /* description exists                              !*/
#define p7_RF      (1<<2)    /* #RF annotation available                        !*/
#define p7_CS      (1<<3)    /* #CS annotation available                        !*/
#define p7_XRAY    (1<<4)    /* obsolete (was: structural data available)       !*/
#define p7_HASPROB (1<<5)    /* obsolete (was: model in probability form)       !*/
#define p7_HASDNA  (1<<6)    /* obsolete (was: protein HMM->DNA seq params set) !*/
#define p7_STATS   (1<<7)    /* obsolete (was: model has EVD stats calibrated)  !*/
#define p7_MAP     (1<<8)    /* alignment map is available                      !*/
#define p7_ACC     (1<<9)    /* accession number is available                   !*/
#define p7_GA      (1<<10)   /* gathering thresholds available                  !*/
#define p7_TC      (1<<11)   /* trusted cutoffs available                       !*/
#define p7_NC      (1<<12)   /* noise cutoffs available                         !*/
#define p7_CA      (1<<13)   /* surface accessibilities available               !*/

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
#define p7_TBM  7		/* B->Mk transition appears in optimized profiles  */

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
#define p7_STX   11 	/* missing data: used esp. for local entry/exits */
#define p7_NSTATETYPES 12


/*****************************************************************
 * 2. P7_PROFILE: a scoring profile, and its implicit model.
 *****************************************************************/

typedef struct p7_profile_s {
  int    mode;         	/* configured algorithm mode (e.g. p7_LOCAL)   */ 
  int    M;		/* number of nodes in the model                */
  int  **tsc;           /* transition scores     [0.6][1.M-1]          */
  int  **msc;           /* match emission scores [0.Kp-1][1.M]         */
  int  **isc;           /* ins emission scores   [0.Kp-1][1.M-1]       */
  int    xsc[4][2];     /* N,E,C,J transitions   [][LOOP,MOVE]         */
  int   *bsc;           /* begin transitions     [1.M]                 */
  int   *esc;		/* end transitions       [1.M]                 */

  /* We also have some probabilities relevant to the search profile but
   * not to the core model.
   */
  float  xt[4][2];	/* [NECJ][MOVE,LOOP] transitions               */
  float *begin;		/* 1..M begin "probabilities"                  */
  float *end;		/* 1..M end "probabilities"                    */

  /* Objects we keep references to */
  const ESL_ALPHABET    *abc;	/* copy of pointer to appropriate alphabet     */
  const struct p7_hmm_s *hmm;	/* who's your daddy                            */
  const struct p7_bg_s  *bg;	/* background null model                       */
  
  /* Numerical correction detail for long sequences */
  int           do_lcorrect;	/* TRUE to apply score correction      */
  float         lscore;		/* the correction to apply per residue */

  /* Flag(s) indicating test/debugging modes */
  int h2_mode;        		/* TRUE if model config is HMMER2 style */

} P7_PROFILE;

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

/* Indices for extra state transitions
 * Used for indexing hmm->xt[][].
 */
#define p7_MOVE 0          /* trNB, trEC, trCT, trJB */
#define p7_LOOP 1          /* trNN, trEJ, trCC, trJJ */

/* Search modes.
 */
#define p7_NO_MODE   0
#define p7_LOCAL     1		/* multihit local:  "fs" mode   */
#define p7_GLOCAL    2		/* multihit glocal: "ls" mode   */
#define p7_UNILOCAL  3		/* unihit local: "sw" mode      */
#define p7_UNIGLOCAL 4		/* unihit glocal: "s" mode      */


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
 * A traceback only makes sense in a triplet (tr, hmm, ax|dsq),
 * for a given HMM (with nodes 1..M) and a given digital sequence 
 * (with positions 1..L).
 * 
 * A traceback may be relative to a profile (usually) or to a core
 * model (in model construction; see build.c). You can tell the
 * difference by looking at the first statetype, tr->st[0]; if it's a
 * p7_STS, it's for a profile, and if it's p7_STB, it's for a core
 * model.
 * 
 * A profile's N,C,J states emit on transition, not on state, so a
 * path of N emits 0 residues, NN emits 1 residue, NNN emits 2
 * residues, and so on. By convention, the trace always associates an
 * emission-on-transition with the trailing (destination) state, so
 * the first N, C, or J is stored in a trace as a nonemitter (i=0).
 *
 * A i coords in a traceback may be relative to an aligned sequence or
 * an unaligned sequence in digital mode.
 */
typedef struct p7_trace_s {
  int   N;		/* length of traceback                       */
  int   nalloc;		/* allocated length of traceback             */
  char *st;		/* state type code                   [0..N-1]*/
  int  *k;		/* node index; 1..M if M,D,I; else 0 [0..N-1]*/
  int  *i;		/* position in dsq, 1..L; 0 if none  [0..N-1]*/
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

typedef struct p7_gmx_s {
  int  M;		/* actual model dimension (model 1..M)    */
  int  L;		/* actual sequence dimension (seq 1..L)   */
  
  size_t ncells;	/* current cell allocation limit: >= (M+1)*(L+1) */
  size_t nrows;    	/* current row allocation limit:  >= L+1  */

  int **xmx;		/* special scores [0.1..L][BECJN]      */
  int **mmx;		/* match scores   [0.1..L][0.1..M]     */
  int **imx;		/* insert scores  [0.1..L][0.1..M-1.M] */
  int **dmx;		/* delete scores  [0.1..L][0.1..M-1.M] */

  int  *xmx_mem;	
  int  *mmx_mem;
  int  *imx_mem;
  int  *dmx_mem;
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
 * 8. Other routines in HMMER's exposed API.
 *****************************************************************/

/* build.c */
extern int p7_Handmodelmaker(ESL_MSA *msa,                P7_HMM **ret_hmm, P7_TRACE ***ret_tr);
extern int p7_Fastmodelmaker(ESL_MSA *msa, float symfrac, P7_HMM **ret_hmm, P7_TRACE ***ret_tr);

/* dp_generic.c */
extern int p7_GViterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *mx, int *ret_sc);
extern int p7_GForward(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *mx, int *ret_sc);
extern int p7_GHybrid (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *mx, int *opt_fwdscore, int *opt_hybscore);
extern int p7_GTrace  (const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *mx, P7_TRACE *tr);

/* emit.c */
extern int p7_CoreEmit      (ESL_RANDOMNESS *r, P7_HMM     *hmm, ESL_SQ *sq, P7_TRACE *tr);
extern int p7_ProfileEmit   (ESL_RANDOMNESS *r, P7_PROFILE *gm,  ESL_SQ *sq, P7_TRACE *tr);
extern int p7_H2_ProfileEmit(ESL_RANDOMNESS *r, P7_PROFILE *gm,  ESL_SQ *sq, P7_TRACE *tr);

/* errors.c */
extern void p7_Die (char *format, ...);
extern void p7_Fail(char *format, ...);

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
extern int   p7_Prob2SILO(float p, float null);
extern int   p7_LL2SILO(float ll, float null);
extern float p7_SILO2Prob(int sc, float null);
extern float p7_SILO2Bitscore(int sc);
extern int   p7_AminoFrequencies(float *f);
extern int   p7_ILogsum(int s1, int s2);

/* modelconfig.c */
extern int p7_ProfileConfig(const P7_HMM *hmm, const P7_BG *bg, P7_PROFILE *gm, int mode);
extern int p7_ReconfigLength(P7_PROFILE *gm, int L);
extern int p7_H2_ProfileConfig(const P7_HMM *hmm, const P7_BG *bg, P7_PROFILE *gm, int mode);

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


/* p7_bg.c */
extern P7_BG *p7_bg_Create(const ESL_ALPHABET *abc);
extern P7_BG *p7_bg_CreateUniform(const ESL_ALPHABET *abc);
extern int    p7_bg_Dump(FILE *ofp, P7_BG *bg);
extern void   p7_bg_Destroy(P7_BG *bg);
extern int    p7_bg_SetLength(P7_BG *bg, int L);
extern int    p7_bg_NullOne(const P7_BG *bg, const ESL_DSQ *dsq, int L, int *ret_sc);

/* p7_gmx.c */
extern P7_GMX *p7_gmx_Create(int allocM, int allocL);
extern int     p7_gmx_GrowTo(P7_GMX *gx, int allocM, int allocL);
extern void    p7_gmx_Destroy(P7_GMX *gx);


/* p7_hmm.c */
/*      1. The P7_HMM object: allocation, initialization, destruction. */
extern P7_HMM *p7_hmm_Create(int M, const ESL_ALPHABET *abc);
extern P7_HMM *p7_hmm_CreateShell(void);
extern int     p7_hmm_CreateBody(P7_HMM *hmm, int M, const ESL_ALPHABET *abc);
extern void    p7_hmm_Destroy(P7_HMM *hmm);
extern int     p7_hmm_CopyParameters(const P7_HMM *src, P7_HMM *dest);
extern P7_HMM *p7_hmm_Duplicate(const P7_HMM *hmm);
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
extern int     p7_hmm_Validate(P7_HMM *hmm, float tol, char *errbuf);
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
extern int         p7_profile_SetNullEmissions(P7_PROFILE *gm);
extern void        p7_profile_Destroy(P7_PROFILE *gm);

extern int         p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1, 
				   char st2, int k2, int *ret_tsc);
extern int         p7_profile_Validate(const P7_PROFILE *gm, float tol);
extern int         p7_profile_Compare(P7_PROFILE *gm1, P7_PROFILE *gm2, float tol);

/* p7_trace.c */
extern P7_TRACE *p7_trace_Create(int N);
extern int  p7_trace_Reuse(P7_TRACE *tr);
extern int  p7_trace_Grow(P7_TRACE *tr);
extern int  p7_trace_GrowTo(P7_TRACE *tr, int N);
extern void p7_trace_Destroy(P7_TRACE *tr);
extern void p7_trace_DestroyArray(P7_TRACE **tr, int N);
extern int  p7_trace_Validate(P7_TRACE *tr, ESL_ALPHABET *abc, ESL_DSQ *sq, char *errbuf);
extern int  p7_trace_Dump(FILE *fp, P7_TRACE *tr, P7_PROFILE *gm, ESL_DSQ *dsq);

extern int  p7_trace_Append(P7_TRACE *tr, char st, int k, int i);
extern int  p7_trace_Reverse(P7_TRACE *tr);
extern int  p7_trace_Count(P7_HMM *hmm, ESL_DSQ *dsq, float wt, P7_TRACE *tr);
extern int  p7_trace_Score(P7_TRACE *tr, ESL_DSQ *dsq, P7_PROFILE *gm, int *ret_sc);
extern int  p7_trace_GetDomainCount(P7_TRACE *tr, int *ret_ndom);
extern int  p7_trace_StateUseCounts(const P7_TRACE *tr, int *counts);
extern int  p7_trace_GetDomainCoords(P7_TRACE *tr, int which, int *ret_i1, int *ret_i2,
				     int *ret_k1, int *ret_k2);




/* seqmodel.c */
extern int p7_Seqmodel(ESL_ALPHABET *abc, ESL_DSQ *dsq, int M, ESL_DMATRIX *P, 
		       float *f, double tmi, double tii, double tmd, double tdd,
		       P7_HMM **ret_hmm);

#endif /*P7_HMMERH_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
